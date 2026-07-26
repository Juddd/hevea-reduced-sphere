#include "reduced_sphere_profile.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char **argv) {
  long long samples=1000000;
  if (argc==3 && std::string(argv[1])=="--dense") samples=std::atoll(argv[2]);
  if (samples<1000) return 2;
  using namespace hevea::sphere;
  double r12=1e100,r3=1e100,ball=1e100,mean=0,meanSquare=0,
    weightedMean=0,weightSum=0,
    maximum=0,minWy=1e100;
  double at12=0,at3=0;
  for(long long k=0;k<samples;k++){
    const double u=-1.0+2.0*k/(samples-1.0),y=u*yInfinity;
    const auto p=evaluate(y);const double e=p.x*p.x,g=p.dx*p.dx+p.dz*p.dz;
    const double a=std::cos(y)*std::cos(y)-e,b=std::sin(y)*std::sin(y)+e-g;
    if(std::abs(u)<=0.75&&a<r12){r12=a;at12=y;}
    if(std::abs(u)<=0.84&&b<r3){r3=b;at3=y;}
    ball=std::min(ball,certifiedRadius*certifiedRadius-e-p.z*p.z);
    const double error=std::hypot(a,1-g),weight=std::cos(y);
    mean+=error;meanSquare+=error*error;weightedMean+=weight*error;weightSum+=weight;
    maximum=std::max(maximum,error);
    if(std::abs(u)<=0.75)minWy=std::min(minWy,(e+a)/(e+g+a));
  }
  const auto lo=evaluate(-yInfinity),hi=evaluate(yInfinity);
  const double c=std::cos(yInfinity),s=std::sin(yInfinity);
  const double gap=std::max({std::abs(lo.x-c),std::abs(hi.x-c),
    std::abs(lo.z-(certifiedEta-s)),std::abs(hi.z-(s-certifiedEta)),
    std::abs(lo.dx-s),std::abs(hi.dx+s),std::abs(lo.dz-c),std::abs(hi.dz-c)});
  std::cout<<std::setprecision(15)<<"samples="<<samples<<" boundary_gap="<<gap
    <<" rho12_margin="<<r12<<" rho12_at="<<at12<<" rho3_margin="<<r3
    <<" rho3_at="<<at3<<" ball_margin="<<ball<<" min_wy="<<minWy
    <<" predicted_min_jacobian="<<minWy<<" metric_mean="<<mean/samples
    <<" metric_rms="<<std::sqrt(meanSquare/samples)
    <<" area_weighted_metric_mean="<<weightedMean/weightSum
    <<" metric_max="<<maximum<<" coefficient_hash="<<coefficientHash<<'\n';
  return gap<1e-10&&r12>1e-6&&r3>1e-6&&ball>1e-6&&minWy>0.02?0:1;
}
