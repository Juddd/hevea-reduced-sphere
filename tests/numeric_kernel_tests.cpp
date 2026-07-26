#include "hevea_numeric.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

int main(){
  double besselWorst=0;
  for(int i=0;i<=1000;i++){double x=double(i)/1000;double a=hevea::j0Inverse(x);besselWorst=std::max(besselWorst,std::abs(std::cyl_bessel_j(0,a)-x));}
  hevea::OdeStats stats;
  auto exponential=[](double,hevea::State<1> const& y){return hevea::State<1>{y[0]};};
  double odeError=std::abs(hevea::integrateDp54<1>(exponential,0,{1},1,1e-12,1e-13,&stats)[0]-std::exp(1.0));
  auto flow=hevea::testFlowInverse(65537,0.1);
  hevea::CylinderGrid<double> grid(8,5,-1,1);grid.at(-1,0)=2;grid.at(7,4)=3;
  bool bounded=false;try{(void)grid.at(0,5);}catch(std::out_of_range const&){bounded=true;}
  auto tangents=hevea::nonuniformHermiteTangents(0.0,1.0,9.0,49.0,1.0,2.0,4.0);
  double tangentError=std::max(std::abs(tangents[0]-5.0),std::abs(tangents[1]-14.0));
  constexpr int cubicNx=16,cubicNy=17;const double cubicDy=2.0/(cubicNy-1);
  auto cubic=[](int,double yIndex){double y=-1.0+yIndex*(2.0/(cubicNy-1));return y*y*y;};
  const double sampleY=.317;
  auto jet=hevea::cylinderCubicJet<double>(cubic,cubicNx,cubicNy,2*std::acos(-1.0)/cubicNx,
                                           cubicDy,-1.0,.271,sampleY);
  double jetError=std::max({std::abs(jet.value-sampleY*sampleY*sampleY),
                            std::abs(jet.derivativeX),
                            std::abs(jet.derivativeY-3.0*sampleY*sampleY)});
  constexpr int mixedNx=256,mixedNy=65;const double mixedDx=2*std::acos(-1.0)/mixedNx,
    mixedDy=2.0/(mixedNy-1);
  auto mixedValue=[](double x,double y){return std::sin(2*x)*(1+.2*y+.3*y*y+.1*y*y*y);};
  auto mixedDxValue=[](double x,double y){return 2*std::cos(2*x)*(1+.2*y+.3*y*y+.1*y*y*y);};
  auto mixedDyValue=[](double x,double y){return std::sin(2*x)*(.2+.6*y+.3*y*y);};
  auto mixedGrid=[&](int i,int j){return mixedValue(i*mixedDx,-1+j*mixedDy);};
  double mixedJetError=0;
  for(double x:{.271,2.713,6.271})for(double y:{-.999,-.317,.413,.999}){
    auto mixedJet=hevea::cylinderCubicJet<double>(mixedGrid,mixedNx,mixedNy,mixedDx,mixedDy,-1,x,y);
    mixedJetError=std::max({mixedJetError,std::abs(mixedJet.value-mixedValue(x,y)),
      std::abs(mixedJet.derivativeX-mixedDxValue(x,y)),
      std::abs(mixedJet.derivativeY-mixedDyValue(x,y))});
  }
  const double period=2*std::acos(-1.0),amplitude=.17;double fullFRoundTrip=0;
  auto fullFMeasure=[&](int inverseSamples){double error=0;
    auto image=[&](int i){double s=period*i/inverseSamples;return s+amplitude*std::sin(s);};
    auto completeF=[&](int i){double s=period*i/inverseSamples;return std::cos(2*s)+.3*std::sin(3*s);};
    for(int k=0;k<inverseSamples;k++){
      double target=period*(k+.371)/inverseSamples;
      auto inverse=hevea::periodicMonotoneCubicSample<double>(inverseSamples,period,target,image,completeF);
      double lo=target-amplitude,hi=target+amplitude;
      for(int n=0;n<64;n++){double mid=(lo+hi)/2;
        if(mid+amplitude*std::sin(mid)<target)lo=mid;else hi=mid;}
      double source=(lo+hi)/2;
      error=std::max(error,std::abs(inverse.value-(std::cos(2*source)+.3*std::sin(3*source))));
      fullFRoundTrip=std::max(fullFRoundTrip,inverse.linearRoundTripError);
    }
    return error;};
  double fullFError1024=fullFMeasure(1024),fullFError2048=fullFMeasure(2048),
    fullFError4096=fullFMeasure(4096);
  std::cout<<"bessel_max="<<besselWorst<<" ode_error="<<odeError
           <<" ode_accepted="<<stats.accepted<<" ode_rejected="<<stats.rejected
           <<" inverse_roundtrip="<<flow.roundTripMax<<" min_jacobian="<<flow.minJacobian
           <<" bounded_y="<<bounded<<" upstream_tangent_error="<<tangentError
           <<" bicubic_jet_error="<<jetError
           <<" bicubic_mixed_jet_error="<<mixedJetError
           <<" full_f_inverse_error_1024="<<fullFError1024
           <<" full_f_inverse_error_2048="<<fullFError2048
           <<" full_f_inverse_error_4096="<<fullFError4096
           <<" full_f_convergence="<<fullFError1024/fullFError2048<<','
           <<fullFError2048/fullFError4096
           <<" full_f_roundtrip="<<fullFRoundTrip<<"\n";
  return besselWorst<1e-10&&odeError<1e-8&&flow.roundTripMax<1e-7&&flow.minJacobian>0&&
    bounded&&tangentError<1e-14&&jetError<1e-12&&mixedJetError<2e-6&&fullFError4096<1e-9&&
    fullFError1024/fullFError2048>7&&fullFError2048/fullFError4096>7&&fullFRoundTrip<1e-14?0:1;
}
