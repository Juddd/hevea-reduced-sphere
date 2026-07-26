#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double pi = 3.1415926535897932384626433832795;

struct V3 {
  double x = 0.0, y = 0.0, z = 0.0;
  V3& operator+=(V3 b) { x += b.x; y += b.y; z += b.z; return *this; }
};

V3 operator*(double s, V3 a) { return {s*a.x,s*a.y,s*a.z}; }
double squaredNorm(V3 a) { return a.x*a.x+a.y*a.y+a.z*a.z; }

struct RowStats {
  double latitude = 0.0;
  V3 mean;
  double nonAxisymmetricRms = 0.0;
  double mode21 = 0.0;
  double mode142 = 0.0;
};

RowStats analyseRow(const std::vector<V3>& row, double latitude,
                    double longitudeSign, double coordinateScale,
                    int longitudeStride) {
  if (!(coordinateScale > 0.0) || longitudeStride < 1)
    throw std::runtime_error("invalid row-analysis scale or stride");
  V3 sum, mode21Cos, mode21Sin, mode142Cos, mode142Sin;
  double sumSquares = 0.0;
  std::size_t count = 0;
  for (std::size_t longitude = 0; longitude < row.size();
       longitude += static_cast<std::size_t>(longitudeStride)) {
    const double phi = longitudeSign*2.0*pi*longitude/row.size();
    const double cosine = std::cos(phi), sine = std::sin(phi);
    const V3 p = (1.0/coordinateScale)*row[longitude];
    const V3 local{p.x*cosine+p.y*sine,-p.x*sine+p.y*cosine,p.z};
    sum += local;
    sumSquares += squaredNorm(local);
    const double c21 = std::cos(21.0*phi), s21 = std::sin(21.0*phi);
    const double c142 = std::cos(142.0*phi), s142 = std::sin(142.0*phi);
    mode21Cos += c21*local; mode21Sin += s21*local;
    mode142Cos += c142*local; mode142Sin += s142*local;
    ++count;
  }
  if (count < 400) throw std::runtime_error("too few longitude samples for mode 142");
  const double inverse = 1.0/count;
  const V3 mean = inverse*sum;
  const double variance = std::max(0.0,sumSquares*inverse-squaredNorm(mean));
  auto amplitude = [&](V3 cosinePart,V3 sinePart) {
    return 2.0*inverse*std::sqrt(squaredNorm(cosinePart)+squaredNorm(sinePart));
  };
  return {latitude,mean,std::sqrt(variance),
          amplitude(mode21Cos,mode21Sin),amplitude(mode142Cos,mode142Sin)};
}

void writeCsv(const std::string& path,const std::vector<RowStats>& rows) {
  std::ofstream output(path);
  if (!output) throw std::runtime_error("cannot create envelope CSV: "+path);
  output << "latitude,meanRadial,meanTangential,meanZ,nonAxisymmetricRms,mode21,mode142\n"
         << std::setprecision(17);
  for (const auto& row : rows)
    output << row.latitude << ',' << row.mean.x << ',' << row.mean.y << ','
           << row.mean.z << ',' << row.nonAxisymmetricRms << ','
           << row.mode21 << ',' << row.mode142 << '\n';
}

std::vector<RowStats> analyseBinary(const std::string& path,int longitudeStride,
                                    int latitudeStride) {
  std::ifstream input(path,std::ios::binary);
  if (!input) throw std::runtime_error("cannot open binary grid: "+path);
  std::array<char,8> magic{}; int nx = 0, ny = 0; double ymax = 0.0;
  input.read(magic.data(),magic.size());
  input.read(reinterpret_cast<char*>(&nx),sizeof(nx));
  input.read(reinterpret_cast<char*>(&ny),sizeof(ny));
  input.read(reinterpret_cast<char*>(&ymax),sizeof(ymax));
  if (magic != std::array<char,8>{'H','E','V','S','P','H','1','\0'} ||
      nx < 400 || ny < 3 || !(ymax > 0.0) || longitudeStride < 1 ||
      latitudeStride < 1 || nx/longitudeStride < 400)
    throw std::runtime_error("invalid binary grid header or sampling stride");
  const std::streamoff header = static_cast<std::streamoff>(magic.size()+
      2*sizeof(int)+sizeof(double));
  const std::streamoff rowBytes = static_cast<std::streamoff>(nx)*sizeof(V3);
  std::vector<V3> row(static_cast<std::size_t>(nx));
  std::vector<RowStats> result;
  result.reserve(static_cast<std::size_t>((ny+latitudeStride-1)/latitudeStride));
  for (int j = 0; j < ny; j += latitudeStride) {
    input.seekg(header+static_cast<std::streamoff>(j)*rowBytes);
    input.read(reinterpret_cast<char*>(row.data()),rowBytes);
    if (!input) throw std::runtime_error("truncated binary grid row");
    const double latitude = -ymax+2.0*ymax*j/(ny-1);
    result.push_back(analyseRow(row,latitude,1.0,1.0,longitudeStride));
  }
  return result;
}

std::vector<V3> readWrlPoints(const std::string& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open WRL: "+path);
  bool points = false;
  std::vector<V3> result;
  result.reserve(4'100'000);
  std::string line;
  while (std::getline(input,line)) {
    if (!points) {
      if (line.find("point [") != std::string::npos) points = true;
      continue;
    }
    if (line.find(']') != std::string::npos) break;
    V3 point;
    if (std::sscanf(line.c_str(),"%lf %lf %lf",&point.x,&point.y,&point.z)==3)
      result.push_back(point);
  }
  if (!points || result.size() < 1000) throw std::runtime_error("WRL point section is incomplete");
  return result;
}

std::vector<RowStats> analyseWrl(const std::string& path,int nx,int latitudeStride,
                                double coordinateScale) {
  const auto points = readWrlPoints(path);
  if (nx < 400 || latitudeStride < 1 || points.size() < 2 ||
      (points.size()-2)%static_cast<std::size_t>(nx) != 0)
    throw std::runtime_error("WRL dimensions do not match the structured point array");
  const int ny = static_cast<int>((points.size()-2)/static_cast<std::size_t>(nx));
  std::vector<V3> row(static_cast<std::size_t>(nx));
  std::vector<RowStats> result;
  result.reserve(static_cast<std::size_t>((ny+latitudeStride-1)/latitudeStride));
  for (int j = 0; j < ny; j += latitudeStride) {
    for (int i = 0; i < nx; ++i)
      row[static_cast<std::size_t>(i)] = points[2+static_cast<std::size_t>(i)*ny+j];
    const double latitude = pi/2.0-(j+1)*pi/(ny+1);
    result.push_back(analyseRow(row,latitude,-1.0,coordinateScale,1));
  }
  return result;
}

std::vector<RowStats> readEnvelopeCsv(const std::string& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open envelope CSV: "+path);
  std::vector<RowStats> rows;
  std::string line;
  std::getline(input,line);
  while (std::getline(input,line)) {
    std::replace(line.begin(),line.end(),',',' ');
    std::istringstream values(line);
    RowStats row;
    if (!(values >> row.latitude >> row.mean.x >> row.mean.y >> row.mean.z
                 >> row.nonAxisymmetricRms >> row.mode21 >> row.mode142))
      throw std::runtime_error("invalid envelope CSV row");
    rows.push_back(row);
  }
  if (rows.size() < 2 || !std::is_sorted(rows.begin(),rows.end(),
      [](const RowStats& a,const RowStats& b){ return a.latitude < b.latitude; }))
    throw std::runtime_error("envelope CSV must contain sorted latitude rows");
  return rows;
}

double interpolateMode21(const std::vector<RowStats>& rows,double latitude) {
  const auto upper=std::lower_bound(rows.begin(),rows.end(),latitude,
    [](const RowStats& row,double value){ return row.latitude < value; });
  if (upper==rows.begin()) {
    if (std::abs(upper->latitude-latitude)>1e-12)
      throw std::runtime_error("reference latitude is outside envelope CSV");
    return upper->mode21;
  }
  if (upper==rows.end()) throw std::runtime_error("reference latitude is outside envelope CSV");
  if (std::abs(upper->latitude-latitude)<=1e-12) return upper->mode21;
  const auto lower=upper-1;
  const double weight=(latitude-lower->latitude)/(upper->latitude-lower->latitude);
  return lower->mode21+weight*(upper->mode21-lower->mode21);
}

int compareEnvelope(const std::string& actualPath,const std::string& referencePath,
                    double latitudeMaximum,double rmsLimit,double maximumLimit) {
  if (!(latitudeMaximum>0.0&&rmsLimit>0.0&&maximumLimit>0.0))
    throw std::runtime_error("envelope comparison limits must be positive");
  const auto actual=readEnvelopeCsv(actualPath);
  std::ifstream input(referencePath);
  if (!input) throw std::runtime_error("cannot open reference envelope: "+referencePath);
  std::string line;std::getline(input,line);
  double sumSquares=0.0,maximum=0.0;std::size_t count=0;
  while (std::getline(input,line)) {
    std::replace(line.begin(),line.end(),',',' ');
    std::istringstream values(line);
    double latitude=0.0,north=0.0,south=0.0;
    if (!(values>>latitude>>north>>south))
      throw std::runtime_error("invalid reference envelope row");
    if (latitude>latitudeMaximum+1e-12) continue;
    const double reference=(north+south)/2.0;
    if (!(reference>0.0)) throw std::runtime_error("reference envelope must be positive");
    const double measured=(interpolateMode21(actual,latitude)+
                           interpolateMode21(actual,-latitude))/2.0;
    const double relative=(measured-reference)/reference;
    sumSquares+=relative*relative;
    maximum=std::max(maximum,std::abs(relative));
    ++count;
  }
  if (count==0) throw std::runtime_error("no reference samples inside comparison band");
  const double rms=std::sqrt(sumSquares/count);
  const bool passed=rms<=rmsLimit&&maximum<=maximumLimit;
  std::cout<<std::setprecision(12)<<"envelope_compare="<<(passed?"pass":"fail")
           <<" samples="<<count<<" latitude_max="<<latitudeMaximum
           <<" relative_rms="<<rms<<" relative_max="<<maximum
           <<" rms_limit="<<rmsLimit<<" max_limit="<<maximumLimit<<'\n';
  return passed?0:3;
}

void printSummary(const std::vector<RowStats>& rows,const std::string& output) {
  double maxRms = 0.0, max21 = 0.0, max142 = 0.0;
  for (const auto& row : rows) {
    maxRms = std::max(maxRms,row.nonAxisymmetricRms);
    max21 = std::max(max21,row.mode21);
    max142 = std::max(max142,row.mode142);
  }
  std::cout << std::setprecision(12) << "corrugation_envelope=pass rows=" << rows.size()
            << " max_nonaxisymmetric_rms=" << maxRms << " max_mode21=" << max21
            << " max_mode142=" << max142 << " output=" << output << '\n';
}

} // namespace

int main(int argc,char** argv) {
  try {
    if (argc == 7 && std::string(argv[1]) == "--compare")
      return compareEnvelope(argv[2],argv[3],std::stod(argv[4]),
                             std::stod(argv[5]),std::stod(argv[6]));
    if (argc == 6 && std::string(argv[1]) == "--binary") {
      const auto rows = analyseBinary(argv[2],std::stoi(argv[4]),std::stoi(argv[5]));
      writeCsv(argv[3],rows); printSummary(rows,argv[3]); return 0;
    }
    if (argc == 7 && std::string(argv[1]) == "--wrl") {
      const auto rows = analyseWrl(argv[2],std::stoi(argv[4]),std::stoi(argv[5]),
                                   std::stod(argv[6]));
      writeCsv(argv[3],rows); printSummary(rows,argv[3]); return 0;
    }
    std::cerr << "usage:\n"
              << "  corrugation_envelope --binary GRID.bin OUTPUT.csv LONGITUDE_STRIDE LATITUDE_STRIDE\n"
              << "  corrugation_envelope --wrl MESH.wrl OUTPUT.csv LONGITUDE_COUNT LATITUDE_STRIDE SCALE\n"
              << "  corrugation_envelope --compare ACTUAL.csv REFERENCE.csv LATITUDE_MAX RMS_LIMIT MAX_LIMIT\n";
    return 2;
  } catch (const std::exception& error) {
    std::cerr << "corrugation_envelope: " << error.what() << '\n'; return 2;
  }
}
