#pragma once

#include <array>
#include <cmath>

namespace hevea::sphere {

// Deterministic symmetric profile.  The production degree-thirteen profile
// was obtained by the reproducible public-WRL envelope search documented in
// docs/profile-certificate.md; the endpoint Hermite data fix the other terms.
inline constexpr double yInfinity = 1.5358895;
inline constexpr double certifiedRadius = 0.52;
inline constexpr double certifiedEta = 0.5;
inline constexpr std::array<double, 6> freeCoefficients{
    -0.51864396687563019, 0.12700737424508893, 0.071711183124306049,
    -0.41551333272698654, 0.47383274947983156, -0.34817367193977911};
using ExtendedCoefficients = std::array<double, 10>;
inline constexpr ExtendedCoefficients extendedFreeCoefficients{
    -0.52164358694220492, -0.024827807054374828, 1.0452244736858278,
    -1.3656657598369506, 0.54300716429149642, -0.38685452062281694,
    0.31457715720002261, -0.61976108254554241, 1.1664678885754221,
    -0.73680735226233729};
inline constexpr char coefficientHash[] =
    "profile-20260725-author-envelope-paper-flow-certified-v5";

struct ProfilePoint {
  double x, z, dx, dz;
};

inline ProfilePoint evaluate(double y, const std::array<double, 6> &v) {
  const double u = y / yInfinity;
  const double t = 1.0 - u * u;
  const double c = std::cos(yInfinity), s = std::sin(yInfinity);
  const double a = yInfinity * s / 2.0;
  const double q0 = s - certifiedEta;
  const double q1 = (q0 - yInfinity * c) / 2.0;
  const double x = c + a*t + v[0]*t*t + v[1]*t*t*t + v[2]*t*t*t*t;
  const double q = q0 + q1*t + v[3]*t*t + v[4]*t*t*t + v[5]*t*t*t*t;
  const double dx = (-2.0*u/yInfinity) *
      (a + 2.0*v[0]*t + 3.0*v[1]*t*t + 4.0*v[2]*t*t*t);
  const double dz = (q - 2.0*u*u *
      (q1 + 2.0*v[3]*t + 3.0*v[4]*t*t + 4.0*v[5]*t*t*t)) / yInfinity;
  return {x, u*q, dx, dz};
}

inline ProfilePoint evaluate(double y, const ExtendedCoefficients &v) {
  const double u = y / yInfinity;
  const double t = 1.0 - u * u;
  const double t2 = t*t, t3 = t2*t, t4 = t3*t, t5 = t4*t, t6 = t5*t;
  const double c = std::cos(yInfinity), s = std::sin(yInfinity);
  const double a = yInfinity * s / 2.0;
  const double q0 = s - certifiedEta;
  const double q1 = (q0 - yInfinity * c) / 2.0;
  const double x = c + a*t + v[0]*t2 + v[1]*t3 + v[2]*t4 + v[3]*t5 + v[4]*t6;
  const double q = q0 + q1*t + v[5]*t2 + v[6]*t3 + v[7]*t4 + v[8]*t5 + v[9]*t6;
  const double dx = (-2.0*u/yInfinity) * (a + 2.0*v[0]*t + 3.0*v[1]*t2 +
      4.0*v[2]*t3 + 5.0*v[3]*t4 + 6.0*v[4]*t5);
  const double dz = (q - 2.0*u*u * (q1 + 2.0*v[5]*t + 3.0*v[6]*t2 +
      4.0*v[7]*t3 + 5.0*v[8]*t4 + 6.0*v[9]*t5)) / yInfinity;
  return {x, u*q, dx, dz};
}

inline ProfilePoint evaluate(double y) { return evaluate(y, extendedFreeCoefficients); }

} // namespace hevea::sphere
