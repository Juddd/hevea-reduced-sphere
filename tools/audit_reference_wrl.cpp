#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

Vec3 operator-(Vec3 a, Vec3 b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 cross(Vec3 a, Vec3 b) {
  return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}

double norm(Vec3 v) {
  return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

struct CapFit {
  double center = 0.0;
  double radius = 0.0;
  double eta = 0.0;
  double boundaryLatitude = 0.0;
  double maximumResidual = 0.0;
  std::size_t rows = 0;
};

struct RibbonFit {
  double firstNonAxisymmetricLatitude = 0.0;
  double firstMode21Latitude = 0.0;
  double firstMode142Latitude = 0.0;
  double maximumAxisymmetricDeviation = 0.0;
  double maximumMode21 = 0.0;
  double maximumMode142 = 0.0;
};

struct MetricAudit {
  double mean = 0.0;
  double maximum = 0.0;
  double maximumLongitude = 0.0;
  double maximumLatitude = 0.0;
  double maximumDeltaE = 0.0;
  double maximumDeltaF = 0.0;
  double maximumDeltaG = 0.0;
};

MetricAudit auditRoundMetric(const std::vector<Vec3>& points,
                             std::size_t longitudeCount,
                             std::size_t meridianRows, double scale) {
  const double pi = std::acos(-1.0);
  const double dx = 2.0*pi/longitudeCount;
  const double dy = pi/(meridianRows+1);
  auto at = [&](long longitude, std::size_t row) -> Vec3 {
    const long wrapped = (longitude%static_cast<long>(longitudeCount)+
        static_cast<long>(longitudeCount))%static_cast<long>(longitudeCount);
    return points[2+static_cast<std::size_t>(wrapped)*meridianRows+row];
  };
  MetricAudit result;
  double sum = 0.0, maximum = 0.0;
  std::size_t count = 0;
  for (std::size_t row = 1; row+1 < meridianRows; ++row) {
    const double latitude = pi/2.0-(row+1)*dy;
    const double roundE = std::cos(latitude)*std::cos(latitude);
    for (std::size_t longitude = 0; longitude < longitudeCount; ++longitude) {
      const Vec3 west = at(static_cast<long>(longitude)-1, row);
      const Vec3 east = at(static_cast<long>(longitude)+1, row);
      const Vec3 north = at(static_cast<long>(longitude), row-1);
      const Vec3 south = at(static_cast<long>(longitude), row+1);
      const Vec3 fx{(east.x-west.x)/(2.0*dx*scale),
                    (east.y-west.y)/(2.0*dx*scale),
                    (east.z-west.z)/(2.0*dx*scale)};
      const Vec3 fy{(north.x-south.x)/(2.0*dy*scale),
                    (north.y-south.y)/(2.0*dy*scale),
                    (north.z-south.z)/(2.0*dy*scale)};
      const double e = fx.x*fx.x+fx.y*fx.y+fx.z*fx.z;
      const double f = fx.x*fy.x+fx.y*fy.y+fx.z*fy.z;
      const double g = fy.x*fy.x+fy.y*fy.y+fy.z*fy.z;
      const double error = std::sqrt((e-roundE)*(e-roundE)+2.0*f*f+(g-1.0)*(g-1.0));
      sum += error;
      if (error > maximum) {
        maximum = error;
        result.maximumLongitude = longitude*dx;
        result.maximumLatitude = latitude;
        result.maximumDeltaE = e-roundE;
        result.maximumDeltaF = f;
        result.maximumDeltaG = g-1.0;
      }
      ++count;
    }
  }
  result.mean = sum/count;
  result.maximum = maximum;
  return result;
}

using Coefficients5 = std::array<long double, 5>;

Coefficients5 solveNormalEquations(
    std::array<std::array<long double, 5>, 5> matrix,
    Coefficients5 rightHandSide) {
  for (std::size_t column = 0; column < 5; ++column) {
    std::size_t pivot = column;
    for (std::size_t row = column+1; row < 5; ++row)
      if (std::abs(matrix[row][column]) > std::abs(matrix[pivot][column])) pivot = row;
    if (std::abs(matrix[pivot][column]) < 1e-24L)
      throw std::runtime_error("singular profile least-squares system");
    std::swap(matrix[pivot], matrix[column]);
    std::swap(rightHandSide[pivot], rightHandSide[column]);
    const long double diagonal = matrix[column][column];
    for (std::size_t k = column; k < 5; ++k) matrix[column][k] /= diagonal;
    rightHandSide[column] /= diagonal;
    for (std::size_t row = 0; row < 5; ++row) {
      if (row == column) continue;
      const long double factor = matrix[row][column];
      for (std::size_t k = column; k < 5; ++k)
        matrix[row][k] -= factor*matrix[column][k];
      rightHandSide[row] -= factor*rightHandSide[column];
    }
  }
  return rightHandSide;
}

std::array<long double, 3> solveThree(
    std::array<std::array<long double, 3>, 3> matrix,
    std::array<long double, 3> rightHandSide) {
  for (std::size_t column = 0; column < 3; ++column) {
    std::size_t pivot = column;
    for (std::size_t row = column+1; row < 3; ++row)
      if (std::abs(matrix[row][column]) > std::abs(matrix[pivot][column])) pivot = row;
    if (std::abs(matrix[pivot][column]) < 1e-28L)
      throw std::runtime_error("singular constrained profile system");
    std::swap(matrix[pivot], matrix[column]);
    std::swap(rightHandSide[pivot], rightHandSide[column]);
    const long double diagonal = matrix[column][column];
    for (std::size_t k = column; k < 3; ++k) matrix[column][k] /= diagonal;
    rightHandSide[column] /= diagonal;
    for (std::size_t row = 0; row < 3; ++row) {
      if (row == column) continue;
      const long double factor = matrix[row][column];
      for (std::size_t k = column; k < 3; ++k)
        matrix[row][k] -= factor*matrix[column][k];
      rightHandSide[row] -= factor*rightHandSide[column];
    }
  }
  return rightHandSide;
}

struct ProfileRecovery {
  double yInfinity = 0.0;
  double endpointGap = 0.0;
  double fitMaximum = 0.0;
  std::array<double, 6> freeCoefficients{};
};

ProfileRecovery recoverInitialProfile(const std::vector<Vec3>& points,
                                      std::size_t longitudeCount, std::size_t meridianRows,
                                      double scale, double eta) {
  // y_1,3 = N*_1,3/(2 N_1,3) = 997/(2*334.92) = 1.4884...
  // from the paper's own table.  Outside D_1,3 the final public mesh equals
  // f_1,2, and outside D_1,2 that in turn equals the untouched f0.  Therefore
  // this narrow band is authoritative initial-profile data.
  constexpr long double bandMinimum = 1.490L;
  constexpr long double bandMaximum = 1.525L;
  constexpr long double squaredCenter =
      (bandMinimum*bandMinimum+bandMaximum*bandMaximum)/2.0L;
  constexpr long double squaredScale =
      (bandMaximum*bandMaximum-bandMinimum*bandMinimum)/2.0L;
  std::array<std::array<long double, 5>, 5> normal{};
  Coefficients5 rhsX{}, rhsZ{};
  std::vector<std::array<long double, 3>> samples;
  const long double pi = std::acos(-1.0L);
  for (std::size_t row = 0; row < meridianRows; ++row) {
    const long double latitude = pi/2.0L-(row+1)*pi/(meridianRows+1);
    if (std::abs(latitude) < bandMinimum || std::abs(latitude) > bandMaximum) continue;
    long double radialSum = 0.0L, zSum = 0.0L;
    for (std::size_t longitude = 0; longitude < longitudeCount; ++longitude) {
      const long double phi = -2.0L*pi*longitude/longitudeCount;
      const auto& p = points[2+longitude*meridianRows+row];
      radialSum += p.x*std::cos(phi)+p.y*std::sin(phi);
      zSum += p.z;
    }
    const long double x = (radialSum/longitudeCount)/scale;
    const long double zOverY = ((zSum/longitudeCount)/scale)/latitude;
    const long double u = (latitude*latitude-squaredCenter)/squaredScale;
    std::array<long double, 5> powers{1.0L, u, u*u, u*u*u, u*u*u*u};
    for (std::size_t i = 0; i < 5; ++i) {
      rhsX[i] += powers[i]*x; rhsZ[i] += powers[i]*zOverY;
      for (std::size_t j = 0; j < 5; ++j) normal[i][j] += powers[i]*powers[j];
    }
    samples.push_back({latitude, x, zOverY});
  }
  long double constrainedY = 0.0L;
  long double constrainedScore = std::numeric_limits<long double>::infinity();
  long double constrainedMaximum = 0.0L;
  std::array<long double, 6> constrainedCoefficients{};
  for (int scan = 0; scan <= 20000; ++scan) {
    const long double yInfinity = 1.530L+scan*(0.010L/20000.0L);
    const long double cosine = std::cos(yInfinity), sine = std::sin(yInfinity);
    const long double a = yInfinity*sine/2.0L;
    const long double q0 = sine-eta;
    const long double q1 = (q0-yInfinity*cosine)/2.0L;
    std::array<std::array<long double, 3>, 3> matrix{};
    std::array<long double, 3> rightX{}, rightZ{};
    for (const auto& sample : samples) {
      const long double u = sample[0]/yInfinity, t = 1.0L-u*u;
      const std::array<long double, 3> basis{t*t, t*t*t, t*t*t*t};
      const long double residualX = sample[1]-(cosine+a*t);
      const long double residualZ = sample[2]*yInfinity-(q0+q1*t);
      for (std::size_t i = 0; i < 3; ++i) {
        rightX[i] += basis[i]*residualX; rightZ[i] += basis[i]*residualZ;
        for (std::size_t j = 0; j < 3; ++j) matrix[i][j] += basis[i]*basis[j];
      }
    }
    const auto coefficientsX = solveThree(matrix, rightX);
    const auto coefficientsZ = solveThree(matrix, rightZ);
    long double score = 0.0L, maximum = 0.0L;
    for (const auto& sample : samples) {
      const long double u = sample[0]/yInfinity, t = 1.0L-u*u;
      const long double predictedX = cosine+a*t+coefficientsX[0]*t*t+
          coefficientsX[1]*t*t*t+coefficientsX[2]*t*t*t*t;
      const long double predictedZOverY = (q0+q1*t+coefficientsZ[0]*t*t+
          coefficientsZ[1]*t*t*t+coefficientsZ[2]*t*t*t*t)/yInfinity;
      const long double errorX = predictedX-sample[1];
      const long double errorZ = predictedZOverY-sample[2];
      score += errorX*errorX+errorZ*errorZ;
      maximum = std::max(maximum, std::max(std::abs(errorX), std::abs(errorZ)));
    }
    if (score < constrainedScore) {
      constrainedScore = score; constrainedMaximum = maximum; constrainedY = yInfinity;
      constrainedCoefficients = {coefficientsX[0], coefficientsX[1], coefficientsX[2],
                                 coefficientsZ[0], coefficientsZ[1], coefficientsZ[2]};
    }
  }
  const Coefficients5 polynomialX = solveNormalEquations(normal, rhsX);
  const Coefficients5 polynomialZ = solveNormalEquations(normal, rhsZ);
  auto evaluate = [&](const Coefficients5& coefficients, long double y) {
    const long double u = (y*y-squaredCenter)/squaredScale;
    long double value = 0.0L, derivativeU = 0.0L;
    for (int k = 4; k >= 0; --k) value = value*u+coefficients[static_cast<std::size_t>(k)];
    for (int k = 4; k >= 1; --k)
      derivativeU = derivativeU*u+k*coefficients[static_cast<std::size_t>(k)];
    return std::array<long double, 2>{value, derivativeU*2.0L*y/squaredScale};
  };
  double fitMaximum = 0.0;
  for (const auto& sample : samples) {
    const auto x = evaluate(polynomialX, sample[0]);
    const auto q = evaluate(polynomialZ, sample[0]);
    fitMaximum = std::max(fitMaximum, static_cast<double>(
        std::max(std::abs(x[0]-sample[1]), std::abs(q[0]-sample[2]))));
  }
  long double bestY = 0.0L, bestError = std::numeric_limits<long double>::infinity();
  for (int index = 0; index <= 200000; ++index) {
    const long double y = 1.515L+index*(0.025L/200000.0L);
    const auto x = evaluate(polynomialX, y);
    const auto q = evaluate(polynomialZ, y);
    const long double z = y*q[0], dz = q[0]+y*q[1];
    const long double error = std::hypot(x[0]-std::cos(y), x[1]+std::sin(y))+
        std::hypot(z-(std::sin(y)-eta), dz-std::cos(y));
    if (error < bestError) { bestError = error; bestY = y; }
  }
  auto compose = [&](const Coefficients5& polynomial, long double multiplier) {
    Coefficients5 result{};
    const long double a = (bestY*bestY-squaredCenter)/squaredScale;
    const long double b = -bestY*bestY/squaredScale;
    constexpr long double choose[5][5]{{1,0,0,0,0},{1,1,0,0,0},{1,2,1,0,0},
      {1,3,3,1,0},{1,4,6,4,1}};
    for (std::size_t k = 0; k < 5; ++k)
      for (std::size_t j = 0; j <= k; ++j)
        result[j] += multiplier*polynomial[k]*choose[k][j]*
            std::pow(a, static_cast<int>(k-j))*std::pow(b, static_cast<int>(j));
    return result;
  };
  const auto tX = compose(polynomialX, 1.0L);
  const auto tZ = compose(polynomialZ, bestY);
  (void)bestError; (void)fitMaximum; (void)tX; (void)tZ;
  return {static_cast<double>(constrainedY),
          static_cast<double>(std::sqrt(constrainedScore/samples.size())),
          static_cast<double>(constrainedMaximum),
          {static_cast<double>(constrainedCoefficients[0]),
           static_cast<double>(constrainedCoefficients[1]),
           static_cast<double>(constrainedCoefficients[2]),
           static_cast<double>(constrainedCoefficients[3]),
           static_cast<double>(constrainedCoefficients[4]),
           static_cast<double>(constrainedCoefficients[5])}};
}

ProfileRecovery recoverGlobalMeanProfile(const std::vector<Vec3>& points,
                                         std::size_t longitudeCount,
                                         std::size_t meridianRows,
                                         double scale, double eta,
                                         double yInfinity) {
  const long double pi = std::acos(-1.0L);
  const long double y13 = 997.0L/(2.0L*334.92L);
  const long double cosine = std::cos(yInfinity), sine = std::sin(yInfinity);
  const long double a = yInfinity*sine/2.0L;
  const long double q0 = sine-eta;
  const long double q1 = (q0-yInfinity*cosine)/2.0L;
  std::array<std::array<long double, 3>, 3> matrixX{}, matrixZ{};
  std::array<long double, 3> rightX{}, rightZ{};
  struct Sample { long double y, x, z; };
  std::vector<Sample> samples;
  for (std::size_t row = 0; row < meridianRows; ++row) {
    const long double latitude = pi/2.0L-(row+1)*pi/(meridianRows+1);
    if (std::abs(latitude) > y13) continue;
    long double radialSum = 0.0L, zSum = 0.0L;
    for (std::size_t longitude = 0; longitude < longitudeCount; ++longitude) {
      const long double phi = -2.0L*pi*longitude/longitudeCount;
      const auto& p = points[2+longitude*meridianRows+row];
      radialSum += p.x*std::cos(phi)+p.y*std::sin(phi);
      zSum += p.z;
    }
    const long double x = (radialSum/longitudeCount)/scale;
    const long double z = (zSum/longitudeCount)/scale;
    const long double u = latitude/yInfinity, t = 1.0L-u*u;
    const std::array<long double, 3> basisX{t*t,t*t*t,t*t*t*t};
    const std::array<long double, 3> basisZ{u*t*t,u*t*t*t,u*t*t*t*t};
    const long double residualX = x-(cosine+a*t);
    const long double residualZ = z-u*(q0+q1*t);
    for (std::size_t i = 0; i < 3; ++i) {
      rightX[i] += basisX[i]*residualX; rightZ[i] += basisZ[i]*residualZ;
      for (std::size_t j = 0; j < 3; ++j) {
        matrixX[i][j] += basisX[i]*basisX[j];
        matrixZ[i][j] += basisZ[i]*basisZ[j];
      }
    }
    samples.push_back({latitude,x,z});
  }
  const auto coefficientsX = solveThree(matrixX,rightX);
  const auto coefficientsZ = solveThree(matrixZ,rightZ);
  long double sumSquares = 0.0L, maximum = 0.0L;
  for (const auto& sample : samples) {
    const long double u=sample.y/yInfinity,t=1.0L-u*u;
    const long double x=cosine+a*t+coefficientsX[0]*t*t+
        coefficientsX[1]*t*t*t+coefficientsX[2]*t*t*t*t;
    const long double z=u*(q0+q1*t+coefficientsZ[0]*t*t+
        coefficientsZ[1]*t*t*t+coefficientsZ[2]*t*t*t*t);
    const long double error=std::hypot(x-sample.x,z-sample.z);
    sumSquares+=error*error;maximum=std::max(maximum,error);
  }
  return {yInfinity,static_cast<double>(std::sqrt(sumSquares/samples.size())),
          static_cast<double>(maximum),
          {static_cast<double>(coefficientsX[0]),static_cast<double>(coefficientsX[1]),
           static_cast<double>(coefficientsX[2]),static_cast<double>(coefficientsZ[0]),
           static_cast<double>(coefficientsZ[1]),static_cast<double>(coefficientsZ[2])}};
}

RibbonFit inspectRibbonModes(const std::vector<Vec3>& points,
                             std::size_t longitudeCount,
                             std::size_t meridianRows) {
  std::vector<Vec3> sums(meridianRows), mode21Cos(meridianRows),
      mode21Sin(meridianRows), mode142Cos(meridianRows), mode142Sin(meridianRows);
  std::vector<double> sumSquares(meridianRows, 0.0);
  const double twoPi = 2.0*std::acos(-1.0);
  for (std::size_t longitude = 0; longitude < longitudeCount; ++longitude) {
    const double phi = -twoPi*longitude/longitudeCount;
    const double cosine = std::cos(phi), sine = std::sin(phi);
    const double c21 = std::cos(21.0*phi), s21 = std::sin(21.0*phi);
    const double c142 = std::cos(142.0*phi), s142 = std::sin(142.0*phi);
    for (std::size_t row = 0; row < meridianRows; ++row) {
      const auto& p = points[2+longitude*meridianRows+row];
      const Vec3 local{p.x*cosine+p.y*sine, -p.x*sine+p.y*cosine, p.z};
      auto accumulate = [&](Vec3& target, double weight) {
        target.x += weight*local.x;
        target.y += weight*local.y;
        target.z += weight*local.z;
      };
      accumulate(sums[row], 1.0);
      accumulate(mode21Cos[row], c21); accumulate(mode21Sin[row], s21);
      accumulate(mode142Cos[row], c142); accumulate(mode142Sin[row], s142);
      sumSquares[row] += local.x*local.x+local.y*local.y+local.z*local.z;
    }
  }
  RibbonFit result;
  constexpr double detectionThreshold = 2e-6;
  bool foundAxisymmetric = false, found21 = false, found142 = false;
  for (std::size_t row = 0; row < meridianRows; ++row) {
    const double inverse = 1.0/longitudeCount;
    const Vec3 mean{sums[row].x*inverse, sums[row].y*inverse, sums[row].z*inverse};
    const double variance = std::max(0.0, sumSquares[row]*inverse-
        (mean.x*mean.x+mean.y*mean.y+mean.z*mean.z));
    const double deviation = std::sqrt(variance);
    auto modeAmplitude = [&](Vec3 cosinePart, Vec3 sinePart) {
      const double scale = 2.0*inverse;
      return scale*std::sqrt(cosinePart.x*cosinePart.x+cosinePart.y*cosinePart.y+
          cosinePart.z*cosinePart.z+sinePart.x*sinePart.x+sinePart.y*sinePart.y+
          sinePart.z*sinePart.z);
    };
    const double amplitude21 = modeAmplitude(mode21Cos[row], mode21Sin[row]);
    const double amplitude142 = modeAmplitude(mode142Cos[row], mode142Sin[row]);
    result.maximumAxisymmetricDeviation = std::max(result.maximumAxisymmetricDeviation, deviation);
    result.maximumMode21 = std::max(result.maximumMode21, amplitude21);
    result.maximumMode142 = std::max(result.maximumMode142, amplitude142);
    const double latitude = std::acos(-1.0)/2.0-
        (row+1)*std::acos(-1.0)/(meridianRows+1);
    if (!foundAxisymmetric && deviation > detectionThreshold) {
      result.firstNonAxisymmetricLatitude = latitude; foundAxisymmetric = true;
    }
    if (!found21 && amplitude21 > detectionThreshold) {
      result.firstMode21Latitude = latitude; found21 = true;
    }
    if (!found142 && amplitude142 > detectionThreshold) {
      result.firstMode142Latitude = latitude; found142 = true;
    }
  }
  return result;
}

CapFit fitTranslatedCap(const std::vector<Vec3>& points, std::size_t meridianRows,
                        bool north) {
  if (meridianRows < 80) throw std::runtime_error("too few meridian rows for cap fit");
  const Vec3 pole = points[north ? 0 : 1];
  const std::size_t first = 2;
  const auto pointAt = [&](std::size_t offset) -> const Vec3& {
    return points[first + (north ? offset : meridianRows - 1 - offset)];
  };

  // The published 2001x2001 mesh samples latitude uniformly.  Near a pole,
  // r_k = R sin(k*pi/2000), which gives a well-conditioned scale estimate
  // even though z is rounded to four decimals.  This also avoids contaminating
  // a generic circle fit with the corrugated ribbon, whose boundary is not
  // known in advance.
  const double latitudeIntervals = static_cast<double>(meridianRows+1);
  double radialDotSine = 0.0, sineSquared = 0.0;
  for (std::size_t offset = 0; offset < 20; ++offset) {
    const auto& p = pointAt(offset);
    const double angle = (offset+1)*std::acos(-1.0)/latitudeIntervals;
    const double sine = std::sin(angle);
    radialDotSine += std::hypot(p.x, p.y)*sine;
    sineSquared += sine*sine;
  }
  const double radius = radialDotSine/sineSquared;
  const double center = pole.z + (north ? -radius : radius);
  if (!(radius > 0.0)) throw std::runtime_error("translated-cap radius is invalid");

  // Coordinates in the public WRL are rounded to 1e-4.  A 2e-3 radial
  // tolerance is deliberately much larger than quantisation, yet far below
  // the first visible corrugation.  Permit two isolated rounded outliers but
  // stop at three consecutive non-spherical rows.
  constexpr double residualTolerance = 2e-3;
  std::size_t lastGood = 0;
  int consecutiveBad = 0;
  double maximumResidual = 0.0;
  for (std::size_t offset = 0; offset < meridianRows; ++offset) {
    const auto& p = pointAt(offset);
    const double radial = std::hypot(p.x, p.y);
    const double residual = std::abs(std::hypot(radial, p.z-center)-radius);
    if (residual <= residualTolerance) {
      lastGood = offset;
      consecutiveBad = 0;
      maximumResidual = std::max(maximumResidual, residual);
    } else if (++consecutiveBad == 3) {
      break;
    }
  }
  const auto& boundary = pointAt(lastGood);
  const double sineLatitude = std::clamp((boundary.z-center)/radius, -1.0, 1.0);
  return {center, radius, std::abs(center)/radius, std::asin(std::abs(sineLatitude)),
          maximumResidual, lastGood+1};
}

} // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 2) {
      std::cerr << "usage: audit_reference_wrl mesh.wrl\n";
      return 2;
    }
    std::ifstream input(argv[1]);
    if (!input) throw std::runtime_error("cannot open WRL file");

    enum class Section { seekPoints, points, seekFaces, faces, done };
    Section section = Section::seekPoints;
    std::vector<Vec3> points;
    points.reserve(4'100'000);
    std::vector<long> face;
    face.reserve(8);
    Vec3 minimum{}, maximum{};
    double boundingRadius = 0.0;
    long long faceCount = 0;
    long long triangleCount = 0;
    long long invalidFaces = 0;
    double area = 0.0;
    std::string line;

    auto finishFace = [&] {
      if (face.size() < 3) {
        ++invalidFaces;
      } else {
        const auto base = static_cast<std::size_t>(face.front());
        if (base >= points.size()) {
          ++invalidFaces;
        } else {
          for (std::size_t k = 1; k + 1 < face.size(); ++k) {
            const auto b = static_cast<std::size_t>(face[k]);
            const auto c = static_cast<std::size_t>(face[k+1]);
            if (b >= points.size() || c >= points.size()) {
              ++invalidFaces;
              continue;
            }
            area += 0.5 * norm(cross(points[b]-points[base], points[c]-points[base]));
            ++triangleCount;
          }
        }
      }
      ++faceCount;
      face.clear();
    };

    while (section != Section::done && std::getline(input, line)) {
      if (section == Section::seekPoints) {
        if (line.find("point [") != std::string::npos) section = Section::points;
      } else if (section == Section::points) {
        if (line.find(']') != std::string::npos) {
          section = Section::seekFaces;
          continue;
        }
        Vec3 point;
        if (std::sscanf(line.c_str(), "%lf %lf %lf", &point.x, &point.y, &point.z) != 3)
          continue;
        if (points.empty()) minimum = maximum = point;
        minimum.x = std::min(minimum.x, point.x);
        minimum.y = std::min(minimum.y, point.y);
        minimum.z = std::min(minimum.z, point.z);
        maximum.x = std::max(maximum.x, point.x);
        maximum.y = std::max(maximum.y, point.y);
        maximum.z = std::max(maximum.z, point.z);
        boundingRadius = std::max(boundingRadius, norm(point));
        points.push_back(point);
      } else if (section == Section::seekFaces) {
        if (line.find("coordIndex [") != std::string::npos) section = Section::faces;
      } else if (section == Section::faces) {
        if (line.find(']') != std::string::npos) {
          if (!face.empty()) finishFace();
          section = Section::done;
          continue;
        }
        const char* cursor = line.c_str();
        while (*cursor != '\0') {
          char* end = nullptr;
          const long value = std::strtol(cursor, &end, 10);
          if (end == cursor) {
            ++cursor;
            continue;
          }
          cursor = end;
          if (value == -1) finishFace();
          else face.push_back(value);
        }
      }
    }
    if (section != Section::done || points.empty() || triangleCount == 0)
      throw std::runtime_error("first IndexedFaceSet is incomplete");

    const double pi = std::acos(-1.0);
    const double intrinsicRadius = std::sqrt(area/(4.0*pi));
    const auto longitudeCount = static_cast<std::size_t>(
        std::gcd(static_cast<long long>(points.size()-2), faceCount));
    if (longitudeCount == 0 || (points.size()-2)%longitudeCount != 0)
      throw std::runtime_error("cannot infer structured WRL dimensions");
    const std::size_t meridianRows = (points.size()-2)/longitudeCount;
    const CapFit northCap = fitTranslatedCap(points, meridianRows, true);
    const CapFit southCap = fitTranslatedCap(points, meridianRows, false);
    const RibbonFit ribbon = inspectRibbonModes(points, longitudeCount, meridianRows);
    const double referenceScale = (northCap.radius+southCap.radius)/2.0;
    const MetricAudit referenceMetric = auditRoundMetric(
        points, longitudeCount, meridianRows, referenceScale);
    const ProfileRecovery profile = recoverInitialProfile(
        points, longitudeCount, meridianRows, referenceScale,
        (northCap.eta+southCap.eta)/2.0);
    const ProfileRecovery globalProfile = recoverGlobalMeanProfile(
        points, longitudeCount, meridianRows, referenceScale,
        (northCap.eta+southCap.eta)/2.0, profile.yInfinity);
    std::cout << std::setprecision(12)
              << "reference_wrl=pass points=" << points.size()
              << " faces=" << faceCount
              << " triangles=" << triangleCount
              << " invalid_faces=" << invalidFaces
              << " area=" << area
              << " intrinsic_radius=" << intrinsicRadius
              << " bounding_radius=" << boundingRadius
              << " reduction_ratio=" << boundingRadius/intrinsicRadius
              << " structured_grid=" << longitudeCount << 'x' << meridianRows+2
              << " north_cap_radius=" << northCap.radius
              << " north_cap_center=" << northCap.center
              << " north_eta=" << northCap.eta
              << " north_y_infinity=" << northCap.boundaryLatitude
              << " north_cap_rows=" << northCap.rows
              << " north_cap_residual=" << northCap.maximumResidual
              << " south_cap_radius=" << southCap.radius
              << " south_cap_center=" << southCap.center
              << " south_eta=" << southCap.eta
              << " south_y_infinity=" << southCap.boundaryLatitude
              << " south_cap_rows=" << southCap.rows
              << " south_cap_residual=" << southCap.maximumResidual
              << " first_nonaxisymmetric_latitude=" << ribbon.firstNonAxisymmetricLatitude
              << " first_mode21_latitude=" << ribbon.firstMode21Latitude
              << " first_mode142_latitude=" << ribbon.firstMode142Latitude
              << " max_nonaxisymmetric_deviation=" << ribbon.maximumAxisymmetricDeviation
              << " max_mode21=" << ribbon.maximumMode21
              << " max_mode142=" << ribbon.maximumMode142
              << " reference_round_mean=" << referenceMetric.mean
              << " reference_round_max=" << referenceMetric.maximum
              << " reference_round_max_x=" << referenceMetric.maximumLongitude
              << " reference_round_max_y=" << referenceMetric.maximumLatitude
              << " reference_round_max_delta=" << referenceMetric.maximumDeltaE << ','
              << referenceMetric.maximumDeltaF << ',' << referenceMetric.maximumDeltaG
              << " recovered_y_infinity=" << profile.yInfinity
              << " recovered_endpoint_gap=" << profile.endpointGap
              << " recovered_fit_max=" << profile.fitMaximum
              << " recovered_profile_coefficients=" << profile.freeCoefficients[0] << ','
              << profile.freeCoefficients[1] << ',' << profile.freeCoefficients[2] << ','
              << profile.freeCoefficients[3] << ',' << profile.freeCoefficients[4] << ','
              << profile.freeCoefficients[5]
              << " global_profile_rms=" << globalProfile.endpointGap
              << " global_profile_max=" << globalProfile.fitMaximum
              << " global_profile_coefficients=" << globalProfile.freeCoefficients[0] << ','
              << globalProfile.freeCoefficients[1] << ',' << globalProfile.freeCoefficients[2] << ','
              << globalProfile.freeCoefficients[3] << ',' << globalProfile.freeCoefficients[4] << ','
              << globalProfile.freeCoefficients[5]
              << " bounds=" << minimum.x << ',' << maximum.x << ';'
              << minimum.y << ',' << maximum.y << ';'
              << minimum.z << ',' << maximum.z << '\n';
    return invalidFaces == 0 ? 0 : 1;
  } catch (std::exception const& error) {
    std::cerr << "audit_reference_wrl: " << error.what() << '\n';
    return 2;
  }
}
