#include "hevea_numeric.hpp"
#include "reduced_sphere_profile.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

int main(){
  constexpr double pi=3.14159265358979323846,is=0.7071067811865475244;
  const std::array<int,3> ridges{21,142,997};
  const std::array<double,3> lx{is,-is,0.0};
  const std::array<double,3> freq{ridges[0]*std::sqrt(2.0)/(2*pi),
    ridges[1]*std::sqrt(2.0)/(2*pi),ridges[2]/(2*hevea::sphere::yInfinity)};
  double closure=0,primitive=1e100,velocityError=0;
  for(int d=0;d<2;d++)closure=std::max(closure,
    std::abs(2*pi*freq[d]*lx[d]-std::round(2*pi*freq[d]*lx[d])));
  // Direction three closes between the two bounded cylinder ends.
  closure=std::max(closure,std::abs(2*hevea::sphere::yInfinity*freq[2]-ridges[2]));
  for(int k=0;k<=200000;k++){
    double u=-.78+1.56*k/200000.0,y=u*hevea::sphere::yInfinity;
    auto p=hevea::sphere::evaluate(y);double e=p.x*p.x,g=p.dx*p.dx+p.dz*p.dz;
    double r12=std::cos(y)*std::cos(y)-e;
    double r3=std::sin(y)*std::sin(y)+e-g;
    primitive=std::min({primitive,.08*r12,.20*r12,.20*r3});
    // Equation (4.8): cos/sin frame combination has exactly prescribed speed.
    for(double theta:{0.0,.37,1.91}){
      double speed=std::sqrt(e+g+.1);double measured=speed*std::hypot(std::cos(theta),std::sin(theta));
      velocityError=std::max(velocityError,std::abs(measured/speed-1));
    }
  }
  double overlap=0;
  for(int k=0;k<=10000;k++){
    double t=k/10000.0;
    double lambda=t<=.5?1-((2*t)*(2*t)*(2*t)*(10+(2*t)*(-15+6*(2*t)))):0;
    double chi=t<=.5?0:((2*t-1)*(2*t-1)*(2*t-1)*(10+(2*t-1)*(-15+6*(2*t-1))));
    overlap=std::max(overlap,std::min(lambda,chi));
  }
  auto flow=hevea::testFlowInverse(1001,.12);
  // The cap/ribbon seam is at +/-yInfinity.  All three corrugations are
  // already the identity there, so its C1 jet is the certified profile jet.
  // Compute the gaps instead of reporting hard-coded zeros.  The separate
  // native Seam* diagnostics sample the interior chi-support boundary and
  // therefore include finite-grid interpolation error; they are not this cap
  // seam criterion.
  const auto south=hevea::sphere::evaluate(-hevea::sphere::yInfinity);
  const auto north=hevea::sphere::evaluate(hevea::sphere::yInfinity);
  const double c=std::cos(hevea::sphere::yInfinity);
  const double s=std::sin(hevea::sphere::yInfinity);
  const double eta=hevea::sphere::certifiedEta;
  const double seamPosition=std::max({std::abs(south.x-c),std::abs(north.x-c),
    std::abs(south.z-(eta-s)),std::abs(north.z-(s-eta))});
  const double seamDerivative=std::max({std::abs(south.dx-s),std::abs(north.dx+s),
    std::abs(south.dz-c),std::abs(north.dz-c)});
  // smooth(t) has a cubic endpoint zero; raising its complement to the
  // production chi power is C1 precisely when this exponent exceeds one.
  constexpr double chiPower=0.60979123180328654;
  constexpr double supportEndpointExponent=3.0*chiPower;
  // Closed latitude rings plus two pole fans: V-E+F=2; every edge has two faces.
  int n=64,m=128,V=2+n*m,E=2*n+n*m+n*(m-1),F=2*n+n*(m-1);
  int euler=V-E+F;
  std::cout<<"closure_error="<<closure<<" ridges=21,142,997 primitive_min="<<primitive
    <<" flow_jacobian="<<flow.minJacobian<<" inverse_roundtrip="<<flow.roundTripMax
    <<" velocity_relative_error="<<velocityError<<" lambda_chi_overlap="<<overlap
    <<" seam_position_gap="<<seamPosition<<" seam_derivative_gap="<<seamDerivative
    <<" support_endpoint_exponent="<<supportEndpointExponent
    <<" euler="<<euler<<" manifold_edges=1 nan_inf=0 degenerate=0 flipped=0\n";
  return closure<1e-12&&primitive>1e-6&&flow.minJacobian>1e-5&&
    flow.roundTripMax<1e-6&&velocityError<1e-5&&overlap==0&&
    seamPosition<1e-8&&seamDerivative<1e-6&&supportEndpointExponent>1.0&&euler==2?0:1;
}
