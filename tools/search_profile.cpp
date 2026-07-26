#include "reduced_sphere_profile.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

namespace {

using Coefficients = std::array<double, 6>;

struct Sample {
  double u, t, t2, t3, t4, sine2, cosine2;
};

struct Metrics {
  double score=std::numeric_limits<double>::infinity();
  double mean=0, maximum=0, xmax=0, zmax=0, radius=0;
  double xAtZmax=0;
  double minDz=std::numeric_limits<double>::infinity();
  std::array<double,3> minWy{std::numeric_limits<double>::infinity(),
    std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()};
  std::array<double,3> minPrimitive{std::numeric_limits<double>::infinity(),
    std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()};
  double rho12=std::numeric_limits<double>::infinity();
  double rho3=std::numeric_limits<double>::infinity();
  double coreRho12=std::numeric_limits<double>::infinity();
  double coreRho3=std::numeric_limits<double>::infinity();
  double middleRho3=std::numeric_limits<double>::infinity();
  double outerRms=0,outerMax=0;
  double globalRms=0,globalMax=0;
  double globalDerivativeRms=0,globalDerivativeMax=0;
  double envelopeRms=0,envelopeMax=0;
};

double square(double x){return x*x;}

double requiredPrimitive3(){
  static const double value=[] {
    if(const char* text=std::getenv("HEVEA_SEARCH_MIN_PRIMITIVE3"))return std::stod(text);
    return .033;
  }();
  return value;
}

double minimumMetricMean(){
  static const double value=[] {
    if(const char* text=std::getenv("HEVEA_SEARCH_MIN_METRIC_MEAN"))return std::stod(text);
    return .855;
  }();
  return value;
}

double targetFraction(){
  static const double value=[] {
    if(const char* text=std::getenv("HEVEA_SEARCH_TARGET_FRACTION"))return std::stod(text);
    return .215;
  }();
  return value;
}

struct M2 { double e=0,f=0,g=0; };

double smooth(double x){
  x=std::clamp(x,0.0,1.0);
  return x*x*x*(10+x*(-15+6*x));
}

double lambda(double y,double previous,double outer){
  const double midpoint=(previous+outer)/2;
  return y<=previous?1.0:(y>=midpoint?0.0:
    1-smooth((y-previous)/(midpoint-previous)));
}

double bilinear(M2 const&m,double ax,double ay,double bx,double by){
  return m.e*ax*bx+m.f*(ax*by+ay*bx)+m.g*ay*by;
}

double quadratic(M2 const&m,double x,double y){
  return m.e*x*x+2*m.f*x*y+m.g*y*y;
}

double inverseJ0(double value){
  value=std::clamp(value,0.0,1.0);
  if(value>=1.0-1e-15)return 0.0;
  double lo=0,hi=2.404825557695773,x=hi*std::sqrt(1-value);
  for(int iteration=0;iteration<10;iteration++){
    const double f=std::cyl_bessel_j(0,x)-value;
    if(f>0)lo=x;else hi=x;
    const double derivative=-std::cyl_bessel_j(1,x);
    const double next=std::abs(derivative)>1e-12?x-f/derivative:(lo+hi)/2;
    x=next>lo&&next<hi?next:(lo+hi)/2;
  }
  return x;
}

double firstModeAmplitude(Coefficients const&v,double y,double fraction){
  using namespace hevea::sphere;
  const auto p=evaluate(y,v);const double invSqrt2=1/std::sqrt(2.0);
  const M2 current{p.x*p.x,0,p.dx*p.dx+p.dz*p.dz};
  const double rho=fraction*(std::cos(y)*std::cos(y)-current.e);
  const double lx=invSqrt2,ly=invSqrt2,vx=-invSqrt2,vy=invSqrt2;
  const M2 mu{current.e+rho*lx*lx,rho*lx*ly,current.g+rho*ly*ly};
  const double ux=lx,uy=ly;
  const double z=-bilinear(mu,ux,uy,vx,vy)/quadratic(mu,vx,vy);
  const double wx=ux+z*vx,wy=uy+z*vy;
  const double oldSpeed=std::sqrt(quadratic(current,wx,wy));
  const double targetSpeed=std::sqrt(quadratic(mu,wx,wy));
  const double alpha=inverseJ0(oldSpeed/targetSpeed);
  return 2*targetSpeed*std::cyl_bessel_j(1,alpha)/(21*std::sqrt(2.0));
}

Metrics evaluate(Coefficients const&v,std::vector<Sample>const&samples){
  using namespace hevea::sphere;
  constexpr Coefficients officialOuterCoefficients{
    -0.356091328875,-0.00326357027649,0.44899529457,
    -0.311525848255,0.0997272648148,-0.21343506489};
  constexpr Coefficients officialGlobalMeanCoefficients{
    -0.358054775023,0.0509868430662,0.0000858502342762,
    -0.311467514785,0.0864614730285,-0.00861692529856};
  constexpr double profileRadius=.515;
  const double c=std::cos(yInfinity),s=std::sin(yInfinity);
  const double a=yInfinity*s/2,q0=s-certifiedEta,q1=(q0-yInfinity*c)/2;
  const double targetFractionValue=targetFraction();
  constexpr double y13=997.0/(2.0*334.92);
  const std::array<double,3> previous{.70,1.28,1.38};
  const std::array<double,3> outer{1.28,1.38,y13};
  const double invSqrt2=1/std::sqrt(2.0);
  const std::array<std::array<double,4>,3> basis{{
    {invSqrt2,invSqrt2,-invSqrt2,invSqrt2},
    {-invSqrt2,invSqrt2,-invSqrt2,-invSqrt2},
    {0,1,-1,0}}};
  Metrics m;double violation=0,sum=0,maxRadius2=0,outerSum=0,globalSum=0,globalDerivativeSum=0;
  std::size_t outerCount=0;
  for(std::size_t k=0;k<samples.size();k++){
    auto const&p=samples[k];
    double x=c+a*p.t+v[0]*p.t2+v[1]*p.t3+v[2]*p.t4;
    double q=q0+q1*p.t+v[3]*p.t2+v[4]*p.t3+v[5]*p.t4;
    double z=p.u*q;
    double dx=(-2*p.u/yInfinity)*(a+2*v[0]*p.t+3*v[1]*p.t2+4*v[2]*p.t3);
    double dz=(q-2*p.u*p.u*(q1+2*v[3]*p.t+3*v[4]*p.t2+4*v[5]*p.t3))/yInfinity;
    m.minDz=std::min(m.minDz,dz);
    double g=dx*dx+dz*dz,r12=p.cosine2-x*x,r3=p.sine2+x*x-g;
    double error=std::hypot(r12,1-g),radius2=x*x+z*z;
    sum+=error;m.maximum=std::max(m.maximum,error);m.xmax=std::max(m.xmax,x);
    if(z>m.zmax){m.zmax=z;m.xAtZmax=x;}maxRadius2=std::max(maxRadius2,radius2);
    if(k+1<samples.size()){
      m.rho12=std::min(m.rho12,r12);m.rho3=std::min(m.rho3,r3);
      if(x<0)violation+=square(x);
      if(k>0&&z<0)violation+=square(z);
      if(r12<0)violation+=square(r12);
      if(r3<0)violation+=square(r3);
    }
    if(yInfinity*p.u>=y13){
      const auto reference=hevea::sphere::evaluate(yInfinity*p.u,officialOuterCoefficients);
      const double outerError=std::hypot(x-reference.x,z-reference.z);
      outerSum+=outerError*outerError;m.outerMax=std::max(m.outerMax,outerError);outerCount++;
    }
    const auto globalReference=hevea::sphere::evaluate(yInfinity*p.u,officialGlobalMeanCoefficients);
    const double globalError=std::hypot(x-globalReference.x,z-globalReference.z);
    const double globalDerivativeError=std::hypot(dx-globalReference.dx,dz-globalReference.dz);
    globalSum+=globalError*globalError;m.globalMax=std::max(m.globalMax,globalError);
    globalDerivativeSum+=globalDerivativeError*globalDerivativeError;
    m.globalDerivativeMax=std::max(m.globalDerivativeMax,globalDerivativeError);
    if(p.u<=.75)m.coreRho12=std::min(m.coreRho12,r12);
    if(p.u<=.84)m.coreRho3=std::min(m.coreRho3,r3);
    if(p.u>=.45&&p.u<=.60)m.middleRho3=std::min(m.middleRho3,r3);
    M2 initial{x*x,0,g},current=initial;
    M2 round{p.cosine2,0,1};
    const double y=p.u*yInfinity;
    for(int direction=0;direction<3;direction++){
      const double blend=lambda(y,previous[direction],outer[direction]);
      M2 desired{initial.e+blend*targetFractionValue*(round.e-initial.e),
        initial.f+blend*targetFractionValue*(round.f-initial.f),
        initial.g+blend*targetFractionValue*(round.g-initial.g)};
      M2 defect{desired.e-current.e,desired.f-current.f,desired.g-current.g};
      const double rho=direction==0?defect.e+defect.f:
        (direction==1?defect.e-defect.f:defect.g-defect.e);
      if(y<=previous[direction])m.minPrimitive[direction]=std::min(m.minPrimitive[direction],rho);
      const auto&b=basis[direction];
      const double lx=b[0],ly=b[1],vx=b[2],vy=b[3];
      M2 mu{current.e+rho*lx*lx,current.f+rho*lx*ly,current.g+rho*ly*ly};
      const double ux=lx/(lx*lx+ly*ly),uy=ly/(lx*lx+ly*ly);
      const double z=-bilinear(mu,ux,uy,vx,vy)/std::max(1e-14,quadratic(mu,vx,vy));
      const double wy=uy+z*vy;
      if(y<=outer[direction])m.minWy[direction]=std::min(m.minWy[direction],wy);
      current={current.e+rho*lx*lx,current.f+rho*lx*ly,current.g+rho*ly*ly};
    }
  }
  m.mean=sum/samples.size();m.radius=std::sqrt(maxRadius2);
  m.outerRms=std::sqrt(outerSum/std::max<std::size_t>(1,outerCount));
  m.globalRms=std::sqrt(globalSum/samples.size());
  m.globalDerivativeRms=std::sqrt(globalDerivativeSum/samples.size());
  constexpr std::array<std::array<double,2>,14> referenceEnvelope{{
    {{.05,.02136070191256071}},{{.10,.021239845639787338}},
    {{.15,.021101155760065916}},{{.20,.02093983680283769}},
    {{.25,.020686933593500177}},{{.30,.020336700888866903}},
    {{.35,.01996043506531299}},{{.40,.01956974396646415}},
    {{.45,.019083986753773544}},{{.50,.01843256447224033}},
    {{.55,.017580914283186404}},{{.60,.01652427827340567}},
    {{.65,.015279632336890612}},{{.70,.013871590509470856}}
  }};
  double envelopeSquares=0;
  for(const auto& reference:referenceEnvelope){
    const double error=firstModeAmplitude(v,reference[0],targetFractionValue)-reference[1];
    envelopeSquares+=error*error;m.envelopeMax=std::max(m.envelopeMax,std::abs(error));
  }
  m.envelopeRms=std::sqrt(envelopeSquares/referenceEnvelope.size());
  violation+=square(std::max(0.0,1e-6-m.coreRho12));
  violation+=square(std::max(0.0,1e-6-m.coreRho3));
  violation+=square(std::max(0.0,1e-6-m.middleRho3));
  violation+=square(std::max(0.0,1e-4-m.minDz));
  violation+=1e6*square(std::max(0.0,-m.rho12));
  violation+=1e6*square(std::max(0.0,-m.rho3));
  violation+=square(std::max(0.0,.08-m.minPrimitive[0]));
  violation+=square(std::max(0.0,.005-m.minPrimitive[1]));
  violation+=square(std::max(0.0,requiredPrimitive3()-m.minPrimitive[2]));
  violation+=square(std::max(0.0,.10-m.coreRho3));
  for(int direction=0;direction<3;direction++)
    violation+=square(std::max(0.0,.01-m.minWy[direction]));
  const double fullZmax=std::max(m.zmax,1-certifiedEta);
  violation+=square(std::max(0.0,m.radius-profileRadius));
  violation+=square(std::max(0.0,minimumMetricMean()-m.mean));
  violation+=square(std::max(0.0,m.mean-.95));
  violation+=square(std::max(0.0,std::abs(m.maximum-1.17)-.08));
  violation+=square(std::max(0.0,std::abs(m.xmax-fullZmax)-.05));
  violation+=square(std::max(0.0,m.outerRms-3e-4));
  violation+=square(std::max(0.0,m.outerMax-1e-3));
  if(violation>0)m.score=1e12+1e15*violation;
  else m.score=1e8*(square(m.mean-.90)+square(m.maximum-1.17))+
    1e9*square(m.globalRms)+1e8*square(m.globalMax)+
    1e8*square(m.globalDerivativeRms)+1e7*square(m.globalDerivativeMax)+
    3e11*square(m.envelopeRms)+5e10*square(m.envelopeMax)+
    1e4*square(m.xmax-fullZmax)+1e8*square(m.outerRms)-1e3*m.minPrimitive[2];
  return m;
}

void print(Coefficients const&v,Metrics const&m){
  std::cout<<std::setprecision(17)<<"score="<<m.score<<" coefficients=";
  for(std::size_t i=0;i<v.size();i++)std::cout<<(i?",":"")<<v[i];
  std::cout<<" mean="<<m.mean<<" max="<<m.maximum<<" xmax="<<m.xmax
    <<" zmax="<<m.zmax<<" aspect="<<m.xmax/m.zmax<<" radius="<<m.radius
    <<" x_at_zmax="<<m.xAtZmax
    <<" min_dz="<<m.minDz
    <<" rho12="<<m.rho12<<" rho3="<<m.rho3
    <<" core_rho12="<<m.coreRho12<<" core_rho3="<<m.coreRho3
    <<" middle_rho3="<<m.middleRho3
    <<" min_primitive="<<m.minPrimitive[0]<<','<<m.minPrimitive[1]<<','<<m.minPrimitive[2]
    <<" min_wy="<<m.minWy[0]<<','<<m.minWy[1]<<','<<m.minWy[2]
    <<" outer_rms="<<m.outerRms<<" outer_max="<<m.outerMax
    <<" global_rms="<<m.globalRms<<" global_max="<<m.globalMax
    <<" global_derivative_rms="<<m.globalDerivativeRms
    <<" global_derivative_max="<<m.globalDerivativeMax
    <<" envelope_rms="<<m.envelopeRms<<" envelope_max="<<m.envelopeMax<<'\n';
}

} // namespace

int main(int argc,char**argv){
  constexpr int sampleCount=1201;
  using namespace hevea::sphere;
  std::cout<<std::setprecision(17)
    <<"CONFIG required_primitive3="<<requiredPrimitive3()
    <<" minimum_metric_mean="<<minimumMetricMean()
    <<" target_fraction="<<targetFraction()<<" seed=20260724\n";
  std::vector<Sample>samples; samples.reserve(sampleCount);
  for(int k=0;k<sampleCount;k++){
    double u=static_cast<double>(k)/(sampleCount-1),t=1-u*u,y=u*yInfinity;
    samples.push_back({u,t,t*t,t*t*t,t*t*t*t,square(std::sin(y)),square(std::cos(y))});
  }
  if(argc==8&&std::string(argv[1])=="--evaluate"){
    Coefficients coefficients{};
    for(int i=0;i<6;i++)coefficients[i]=std::stod(argv[i+2]);
    print(coefficients,evaluate(coefficients,samples));
    return 0;
  }
  int generations=argc>1?std::max(1,std::atoi(argv[1])):2000;
  int population=argc>2?std::max(32,std::atoi(argv[2])):384;
  std::array<double,6>lo{-2,-2,-2,-2,-2,-2},hi{2,2,2,2,2,2};
  std::vector<Coefficients>p(population),trial(population);
  std::vector<Metrics>pm(population),tm(population);
  std::mt19937_64 random(20260724);std::uniform_real_distribution<double>unit(0,1);
  Coefficients seedA=freeCoefficients;
  Coefficients seedB{-0.358054775023,0.0509868430662,0.0000858502342762,
    -0.311467514785,0.0864614730285,-0.00861692529856};
  Coefficients seedC{-0.49979018368893402,-0.64805485818266528,0.82703775697302617,
    -0.31963859438593079,-0.25853563532003854,0.36877538103474};
  p[0]=seedA;p[1]=seedB;p[2]=seedC;
  for(int i=3;i<population;i++)for(int j=0;j<6;j++){
    const auto&centerSeed=i%3==0?seedA:(i%3==1?seedB:seedC);
    double center=centerSeed[j],span=hi[j]-lo[j];
    double value=i<population*9/10?center+(unit(random)-.5)*.05*span:lo[j]+unit(random)*span;
    p[i][j]=std::clamp(value,lo[j],hi[j]);
  }
  #pragma omp parallel for
  for(int i=0;i<population;i++)pm[i]=evaluate(p[i],samples);
  for(int generation=0;generation<generations;generation++){
    for(int i=0;i<population;i++){
      int a,b,c;do a=static_cast<int>(unit(random)*population)%population;while(a==i);
      do b=static_cast<int>(unit(random)*population)%population;while(b==i||b==a);
      do c=static_cast<int>(unit(random)*population)%population;while(c==i||c==a||c==b);
      int forced=static_cast<int>(unit(random)*6)%6;
      for(int j=0;j<6;j++){
        double value=(unit(random)<.90||j==forced)?p[a][j]+.72*(p[b][j]-p[c][j]):p[i][j];
        if(value<lo[j]||value>hi[j])value=lo[j]+unit(random)*(hi[j]-lo[j]);
        trial[i][j]=value;
      }
    }
    #pragma omp parallel for
    for(int i=0;i<population;i++)tm[i]=evaluate(trial[i],samples);
    for(int i=0;i<population;i++)if(tm[i].score<pm[i].score){p[i]=trial[i];pm[i]=tm[i];}
    if(generation%100==0||generation+1==generations){
      int best=static_cast<int>(std::min_element(pm.begin(),pm.end(),[](auto const&a,auto const&b){return a.score<b.score;})-pm.begin());
      std::cout<<"generation="<<generation<<' ';print(p[best],pm[best]);
    }
  }
  std::vector<int>order(population);for(int i=0;i<population;i++)order[i]=i;
  std::sort(order.begin(),order.end(),[&](int a,int b){return pm[a].score<pm[b].score;});
  std::cout<<"FINAL_CANDIDATES\n";
  for(int k=0;k<std::min(128,population);k++)print(p[order[k]],pm[order[k]]);
}
