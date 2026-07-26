#include "reduced_sphere_profile.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

int main() {
  using namespace hevea::sphere;
  const auto south = evaluate(-yInfinity), north = evaluate(yInfinity);
  const double c = std::cos(yInfinity), s = std::sin(yInfinity);
  double gap = 0.0;
  auto take = [&](double value) { gap = std::max(gap, std::abs(value)); };
  take(south.x-c); take(north.x-c);
  take(south.z-(certifiedEta-s)); take(north.z-(s-certifiedEta));
  take(south.dx-s); take(north.dx+s); take(south.dz-c); take(north.dz-c);

  double rho12 = 1e100, rho3 = 1e100, ball = 1e100;
  double mean = 0.0, maximum = 0.0, minWy = 1e100;
  constexpr int samples = 200001;
  long long metricCount=0;
  for (int k=0; k<samples; ++k) {
    const double u = -1.0 + 2.0*k/(samples-1.0), y = u*yInfinity;
    const auto p = evaluate(y);
    const auto extended = evaluate(y,extendedFreeCoefficients);
    take(extended.x-p.x); take(extended.z-p.z);
    take(extended.dx-p.dx); take(extended.dz-p.dz);
    const double e=p.x*p.x, g=p.dx*p.dx+p.dz*p.dz;
    const double r12=std::cos(y)*std::cos(y)-e;
    const double r3=std::sin(y)*std::sin(y)+e-g;
    if(std::abs(u)<=0.75)rho12=std::min(rho12,r12);
    if(std::abs(u)<=0.84)rho3=std::min(rho3,r3);
    ball=std::min(ball,certifiedRadius*certifiedRadius-p.x*p.x-p.z*p.z);
    const double err=std::hypot(r12,1.0-g);
    mean+=err; maximum=std::max(maximum,err);
    // For a diagonal initial metric the two diagonal primitive fields have
    // the same y component; direction three is identically one.
    if(std::abs(u)<=0.75){const double wy=(e+r12)/(e+g+r12);minWy=std::min(minWy,wy);}
    ++metricCount;
  }
  mean/=metricCount;
  std::cout << "boundary_gap=" << gap << " rho12_margin=" << rho12
            << " rho3_margin=" << rho3 << " ball_margin=" << ball
            << " min_wy=" << minWy << " metric_mean=" << mean
            << " metric_max=" << maximum << " hash=" << coefficientHash << '\n';
  return gap < 1e-10 && rho12 > 1e-6 && rho3 > 1e-6 && ball > 1e-6 &&
         minWy > 0.02 && std::abs(mean-0.90) <= 0.05 &&
         std::abs(maximum-1.17) <= 0.08 ? 0 : 1;
}
