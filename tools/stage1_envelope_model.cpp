#include "reduced_sphere_profile.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Coefficients = std::array<double,10>;
using Complex = std::complex<double>;

struct Target {
  double latitude = 0.0;
  double north = 0.0;
  double south = 0.0;
};

struct State {
  double flowShift = 0.0;
  std::array<Complex,3> mode{};
};

hevea::sphere::ProfilePoint evaluateProfile(double y,const Coefficients& coefficients) {
  using namespace hevea::sphere;
  const double u=y/yInfinity,t=1.0-u*u;
  const double t2=t*t,t3=t2*t,t4=t3*t,t5=t4*t,t6=t5*t;
  const double cosine=std::cos(yInfinity),sine=std::sin(yInfinity);
  const double a=yInfinity*sine/2.0,q0=sine-certifiedEta;
  const double q1=(q0-yInfinity*cosine)/2.0;
  const double x=cosine+a*t+coefficients[0]*t2+coefficients[1]*t3+
    coefficients[2]*t4+coefficients[3]*t5+coefficients[4]*t6;
  const double q=q0+q1*t+coefficients[5]*t2+coefficients[6]*t3+
    coefficients[7]*t4+coefficients[8]*t5+coefficients[9]*t6;
  const double dx=(-2.0*u/yInfinity)*(a+2.0*coefficients[0]*t+
    3.0*coefficients[1]*t2+4.0*coefficients[2]*t3+
    5.0*coefficients[3]*t4+6.0*coefficients[4]*t5);
  const double dz=(q-2.0*u*u*(q1+2.0*coefficients[5]*t+
    3.0*coefficients[6]*t2+4.0*coefficients[7]*t3+
    5.0*coefficients[8]*t4+6.0*coefficients[9]*t5))/yInfinity;
  return {x,u*q,dx,dz};
}

State operator+(const State& a,const State& b) {
  State result{a.flowShift+b.flowShift,a.mode};
  for(std::size_t k=0;k<result.mode.size();++k)result.mode[k]+=b.mode[k];
  return result;
}

State operator*(double scale,const State& state) {
  State result{scale*state.flowShift,state.mode};
  for(auto& value:result.mode)value*=scale;
  return result;
}

double smooth(double value) {
  value=std::clamp(value,0.0,1.0);
  return value*value*value*(10.0+value*(-15.0+6.0*value));
}

double poweredStep(double value,double power) {
  value=std::clamp(value,0.0,1.0);
  if(power==0.0)return smooth(value);
  if(!(std::isfinite(power)&&power>1.0/3.0))
    throw std::runtime_error("transition shape must be zero or greater than one third");
  if(value==0.0||value==1.0)return value;
  const double plateau=std::clamp(smooth(value),0.0,1.0);
  const double left=std::pow(plateau,power),right=std::pow(1.0-plateau,power);
  return left/(left+right);
}

double inverseJ0(double value) {
  value=std::clamp(value,0.0,1.0);
  if(value>=1.0-1e-15)return 0.0;
  double lower=0.0,upper=2.404825557695773;
  double current=upper*std::sqrt(1.0-value);
  for(int iteration=0;iteration<20;++iteration) {
    const double residual=std::cyl_bessel_j(0,current)-value;
    if(residual>0.0)lower=current;else upper=current;
    const double slope=-std::cyl_bessel_j(1,current);
    const double next=std::abs(slope)>1e-12?current-residual/slope:(lower+upper)/2.0;
    current=next>lower&&next<upper?next:(lower+upper)/2.0;
  }
  return current;
}

State derivative(double y,const State& state,const Coefficients& coefficients,
                 double previous,double outer,double power,double split,double fraction) {
  constexpr double inverseSqrt2=0.70710678118654752440084436210485;
  constexpr int ridgeCount=21;
  const auto profile=evaluateProfile(y,coefficients);
  const double e=profile.x*profile.x;
  const double g=profile.dx*profile.dx+profile.dz*profile.dz;
  const double transition=previous+split*(outer-previous);
  const double absoluteY=std::abs(y);
  const double lambda=absoluteY<=previous?1.0:(absoluteY>=transition?0.0:
    1.0-poweredStep((absoluteY-previous)/(transition-previous),power));
  const double rho=lambda*fraction*(std::cos(y)*std::cos(y)-e);
  const double muE=e+rho/2.0,muF=rho/2.0,muG=g+rho/2.0;
  const double ux=inverseSqrt2,uy=inverseSqrt2;
  const double vx=-inverseSqrt2,vy=inverseSqrt2;
  const double muUv=muE*ux*vx+muF*(ux*vy+uy*vx)+muG*uy*vy;
  const double muVv=muE*vx*vx+2.0*muF*vx*vy+muG*vy*vy;
  const double zeta=-muUv/std::max(1e-15,muVv);
  double wx=ux+zeta*vx,wy=uy+zeta*vy;
  if(wy<0.0){wx=-wx;wy=-wy;}
  const double targetSpeed=std::sqrt(std::max(1e-15,
    muE*wx*wx+2.0*muF*wx*wy+muG*wy*wy));
  const double oldSpeed=std::sqrt(std::max(1e-15,e*wx*wx+g*wy*wy));
  const double alpha=inverseJ0(oldSpeed/targetSpeed);
  const double normalScale=std::sqrt(std::max(1e-15,g));
  const double normalRadial=profile.dz/normalScale;
  const double normalZ=-profile.dx/normalScale;
  const double phase=ridgeCount*(state.flowShift+y);
  const Complex oscillation=std::polar(
    2.0*targetSpeed*std::cyl_bessel_j(1,alpha)/wy,phase);
  State result;
  result.flowShift=wx/wy;
  result.mode={oscillation*(normalRadial*std::cos(state.flowShift)),
               oscillation*(normalRadial*std::sin(state.flowShift)),
               oscillation*normalZ};
  return result;
}

State rk4Step(double y,double step,const State& state,const Coefficients& coefficients,
              double previous,double outer,double power,double split,double fraction) {
  const State k1=derivative(y,state,coefficients,previous,outer,power,split,fraction);
  const State k2=derivative(y+step/2.0,state+(step/2.0)*k1,
    coefficients,previous,outer,power,split,fraction);
  const State k3=derivative(y+step/2.0,state+(step/2.0)*k2,
    coefficients,previous,outer,power,split,fraction);
  const State k4=derivative(y+step,state+step*k3,
    coefficients,previous,outer,power,split,fraction);
  return state+(step/6.0)*(k1+2.0*k2+2.0*k3+k4);
}

double amplitude(const State& state,double y,double previous,double outer,double power,double split) {
  const double transition=previous+split*(outer-previous);
  const double absoluteY=std::abs(y);
  const double keep=absoluteY<=transition?0.0:(absoluteY>=outer?1.0:
    poweredStep((absoluteY-transition)/(outer-transition),power));
  double squareSum=0.0;
  for(const auto& value:state.mode)squareSum+=std::norm(value);
  return (1.0-keep)*std::sqrt(squareSum);
}

std::vector<Target> readTargets(const std::string& path) {
  std::ifstream input(path);
  if(!input)throw std::runtime_error("cannot open reference envelope: "+path);
  std::string header;
  std::getline(input,header);
  std::vector<Target> targets;
  Target target;char comma1=0,comma2=0;
  while(input>>target.latitude>>comma1>>target.north>>comma2>>target.south) {
    if(comma1!=','||comma2!=','||!(target.latitude>0.0))
      throw std::runtime_error("invalid reference envelope row");
    targets.push_back(target);
  }
  if(targets.empty())throw std::runtime_error("reference envelope is empty");
  return targets;
}

struct Prediction {
  std::vector<double> north,south;
  double rms=0.0,maximum=0.0,symmetryMaximum=0.0;
};

Prediction predict(const Coefficients& coefficients,double previous,double outer,
                   double power,double split,double fraction,const std::vector<Target>& targets);

struct ProfileMetrics {
  double penalty=0.0,mean=0.0,maximum=0.0,radius=0.0,outerRms=0.0,outerMax=0.0;
  double globalRms=0.0,globalMax=0.0;
  double globalDerivativeRms=0.0,globalDerivativeMax=0.0;
  double maxSlope=0.0;
  std::array<double,3> minPrimitive{std::numeric_limits<double>::infinity(),
    std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()};
  std::array<double,3> minWy{std::numeric_limits<double>::infinity(),
    std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()};
};

double square(double value){return value*value;}

ProfileMetrics profileMetrics(const Coefficients& coefficients,double firstPrevious,
                              double firstOuter,double firstPower,double firstSplit) {
  constexpr Coefficients officialOuterCoefficients{
    -0.356091328875,-0.00326357027649,0.44899529457,0.0,0.0,
    -0.311525848255,0.0997272648148,-0.21343506489,0.0,0.0};
  constexpr Coefficients officialGlobalMeanCoefficients{
    -0.358054775023,0.0509868430662,0.0000858502342762,0.0,0.0,
    -0.311467514785,0.0864614730285,-0.00861692529856,0.0,0.0};
  constexpr int sampleCount=1201;
  constexpr double outerStart=997.0/(2.0*334.92);
  const std::array<double,3> previous{firstPrevious,firstOuter,1.38};
  const std::array<double,3> outer{firstOuter,1.38,outerStart};
  constexpr double inverseSqrt2=0.70710678118654752440084436210485;
  const std::array<std::array<double,4>,3> basis{{
    {inverseSqrt2,inverseSqrt2,-inverseSqrt2,inverseSqrt2},
    {-inverseSqrt2,inverseSqrt2,-inverseSqrt2,-inverseSqrt2},
    {0.0,1.0,-1.0,0.0}}};
  ProfileMetrics metrics;double errorSum=0.0,outerSquareSum=0.0;
  double globalSquareSum=0.0,globalDerivativeSquareSum=0.0;int outerCount=0;
  for(int index=0;index<sampleCount;++index) {
    const double y=hevea::sphere::yInfinity*index/(sampleCount-1.0);
    const auto point=evaluateProfile(y,coefficients);
    const double speedSquare=point.dx*point.dx+point.dz*point.dz;
    const double rho12=std::cos(y)*std::cos(y)-point.x*point.x;
    const double rho3=std::sin(y)*std::sin(y)+point.x*point.x-speedSquare;
    const double metricError=std::hypot(rho12,1.0-speedSquare);
    errorSum+=metricError;metrics.maximum=std::max(metrics.maximum,metricError);
    metrics.radius=std::max(metrics.radius,std::hypot(point.x,point.z));
    if(index+1<sampleCount) {
      metrics.penalty+=square(std::max(0.0,1e-6-point.x));
      if(index>0)metrics.penalty+=square(std::max(0.0,-point.z));
      metrics.penalty+=square(std::max(0.0,1e-4-point.dz));
      metrics.penalty+=1e6*square(std::max(0.0,-rho12));
      metrics.penalty+=1e6*square(std::max(0.0,-rho3));
      if(y<=0.75)metrics.penalty+=square(std::max(0.0,1e-6-rho12));
      if(y<=0.84)metrics.penalty+=square(std::max(0.0,0.10-rho3));
    }
    if(y>=outerStart) {
      const auto reference=evaluateProfile(y,officialOuterCoefficients);
      const double error=std::hypot(point.x-reference.x,point.z-reference.z);
      outerSquareSum+=error*error;metrics.outerMax=std::max(metrics.outerMax,error);++outerCount;
    }
    const auto globalReference=evaluateProfile(y,officialGlobalMeanCoefficients);
    const double globalError=std::hypot(point.x-globalReference.x,
                                        point.z-globalReference.z);
    const double globalDerivativeError=std::hypot(point.dx-globalReference.dx,
                                                  point.dz-globalReference.dz);
    globalSquareSum+=globalError*globalError;
    globalDerivativeSquareSum+=globalDerivativeError*globalDerivativeError;
    metrics.globalMax=std::max(metrics.globalMax,globalError);
    metrics.globalDerivativeMax=std::max(metrics.globalDerivativeMax,
                                         globalDerivativeError);
    double currentE=point.x*point.x,currentF=0.0,currentG=speedSquare;
    const double roundE=std::cos(y)*std::cos(y),roundG=1.0;
    for(int direction=0;direction<3;++direction) {
      const double split=direction==0?firstSplit:0.5;
      const double transition=previous[direction]+split*(outer[direction]-previous[direction]);
      const double power=direction==0?firstPower:0.0;
      const double lambda=y<=previous[direction]?1.0:(y>=transition?0.0:
        1.0-poweredStep((y-previous[direction])/(transition-previous[direction]),power));
      const double desiredE=point.x*point.x+lambda*0.237*(roundE-point.x*point.x);
      const double desiredF=0.0;
      const double desiredG=speedSquare+lambda*0.237*(roundG-speedSquare);
      const double defectE=desiredE-currentE,defectF=desiredF-currentF,defectG=desiredG-currentG;
      const double rho=direction==0?defectE+defectF:
        (direction==1?defectE-defectF:defectG-defectE);
      if(y<=previous[direction])metrics.minPrimitive[direction]=
        std::min(metrics.minPrimitive[direction],rho);
      const auto& vectors=basis[direction];
      const double lx=vectors[0],ly=vectors[1],vx=vectors[2],vy=vectors[3];
      const double muE=currentE+rho*lx*lx,muF=currentF+rho*lx*ly,muG=currentG+rho*ly*ly;
      const double muUv=muE*lx*vx+muF*(lx*vy+ly*vx)+muG*ly*vy;
      const double muVv=muE*vx*vx+2.0*muF*vx*vy+muG*vy*vy;
      const double zeta=-muUv/std::max(1e-15,muVv);
      const double wx=lx+zeta*vx,wy=ly+zeta*vy;
      if(y<=outer[direction])metrics.minWy[direction]=
        std::min(metrics.minWy[direction],wy);
      if(direction==0&&y<=outer[direction])metrics.maxSlope=std::max(
        metrics.maxSlope,std::abs(wx/std::max(1e-15,std::abs(wy))));
      currentE=muE;currentF=muF;currentG=muG;
    }
  }
  metrics.mean=errorSum/sampleCount;
  metrics.outerRms=std::sqrt(outerSquareSum/std::max(1,outerCount));
  metrics.globalRms=std::sqrt(globalSquareSum/sampleCount);
  metrics.globalDerivativeRms=std::sqrt(globalDerivativeSquareSum/sampleCount);
  metrics.penalty+=square(std::max(0.0,metrics.radius-0.515));
  metrics.penalty+=square(std::max(0.0,0.850-metrics.mean));
  metrics.penalty+=square(std::max(0.0,metrics.mean-0.95));
  metrics.penalty+=square(std::max(0.0,std::abs(metrics.maximum-1.17)-0.08));
  metrics.penalty+=square(std::max(0.0,metrics.outerRms-3e-4));
  metrics.penalty+=square(std::max(0.0,metrics.outerMax-1e-3));
  metrics.penalty+=square(std::max(0.0,0.08-metrics.minPrimitive[0]));
  metrics.penalty+=square(std::max(0.0,0.005-metrics.minPrimitive[1]));
  metrics.penalty+=square(std::max(0.0,0.033-metrics.minPrimitive[2]));
  metrics.penalty+=square(std::max(0.0,0.12-metrics.minWy[0]));
  metrics.penalty+=square(std::max(0.0,metrics.maxSlope-10.0));
  for(std::size_t direction=1;direction<metrics.minWy.size();++direction)
    metrics.penalty+=square(std::max(0.0,0.01-metrics.minWy[direction]));
  return metrics;
}

struct SearchMetrics {
  double score=std::numeric_limits<double>::infinity();
  Prediction envelope;
  ProfileMetrics profile;
  double relativeRms=0.0,weightedRms=0.0;
};

using Candidate=std::array<double,14>;

SearchMetrics evaluateCandidate(const Candidate& candidate,const std::vector<Target>& targets) {
  Coefficients coefficients{};
  std::copy_n(candidate.begin(),coefficients.size(),coefficients.begin());
  SearchMetrics metrics;metrics.profile=profileMetrics(
    coefficients,candidate[10],candidate[11],candidate[12],candidate[13]);
  if(metrics.profile.penalty>0.0) {
    metrics.score=1e12+1e15*metrics.profile.penalty;
    return metrics;
  }
  metrics.envelope=predict(coefficients,candidate[10],candidate[11],candidate[12],
    candidate[13],0.237,targets);
  double relativeSquares=0.0,weightedSquares=0.0,weightSum=0.0;
  std::size_t relativeCount=0;
  for(std::size_t k=0;k<targets.size();++k) {
    for(const auto& pair:std::array{std::pair{metrics.envelope.north[k],targets[k].north},
                                    std::pair{metrics.envelope.south[k],targets[k].south}}) {
      const double scale=std::max(0.001,pair.second);
      if(targets[k].latitude<=1.10) {
        relativeSquares+=square((pair.first-pair.second)/scale);++relativeCount;
      }
      const double weight=targets[k].latitude>=0.95&&targets[k].latitude<=1.10?8.0:
        (targets[k].latitude>=0.65&&targets[k].latitude<0.95?5.0:
        (targets[k].latitude>1.10?0.2:2.0));
      weightedSquares+=weight*square(pair.first-pair.second);weightSum+=weight;
    }
  }
  metrics.relativeRms=std::sqrt(relativeSquares/relativeCount);
  metrics.weightedRms=std::sqrt(weightedSquares/weightSum);
  metrics.score=2e8*square(metrics.weightedRms)+5e7*square(metrics.envelope.maximum)+
    1e3*square(metrics.relativeRms)+1e8*square(metrics.profile.outerRms)+
    2e9*square(metrics.profile.globalRms)+2e8*square(metrics.profile.globalMax)+
    2e8*square(metrics.profile.globalDerivativeRms)+
    2e7*square(metrics.profile.globalDerivativeMax);
  return metrics;
}

void printCandidate(const Candidate& candidate,const SearchMetrics& metrics) {
  std::cout<<std::setprecision(17)<<"score="<<metrics.score<<" candidate=";
  for(std::size_t k=0;k<candidate.size();++k)std::cout<<(k?",":"")<<candidate[k];
  std::cout<<" envelope_rms="<<metrics.envelope.rms
    <<" envelope_max="<<metrics.envelope.maximum
    <<" weighted_rms="<<metrics.weightedRms
    <<" relative_rms="<<metrics.relativeRms
    <<" symmetry_max="<<metrics.envelope.symmetryMaximum
    <<" profile_mean="<<metrics.profile.mean
    <<" profile_max="<<metrics.profile.maximum
    <<" radius="<<metrics.profile.radius
    <<" outer_rms="<<metrics.profile.outerRms
    <<" outer_max="<<metrics.profile.outerMax
    <<" global_rms="<<metrics.profile.globalRms
    <<" global_max="<<metrics.profile.globalMax
    <<" global_derivative_rms="<<metrics.profile.globalDerivativeRms
    <<" global_derivative_max="<<metrics.profile.globalDerivativeMax
    <<" max_slope="<<metrics.profile.maxSlope
    <<" min_primitive="<<metrics.profile.minPrimitive[0]<<','
    <<metrics.profile.minPrimitive[1]<<','<<metrics.profile.minPrimitive[2]
    <<" min_wy="<<metrics.profile.minWy[0]<<','<<metrics.profile.minWy[1]<<','
    <<metrics.profile.minWy[2]<<" profile_penalty="<<metrics.profile.penalty<<'\n';
}

void search(const std::vector<Target>& targets,int generations,int populationSize) {
  populationSize=std::max(32,populationSize);
  const std::array<double,14> lower{-1.5,-1.5,-1.5,-1.5,-1.5,
    -1.5,-1.5,-1.5,-1.5,-1.5,0.60,1.279,0.34,0.50};
  const std::array<double,14> upper{1.5,1.5,1.5,1.5,1.5,
    1.5,1.5,1.5,1.5,1.5,0.85,1.281,3.0,0.95};
  std::vector<Candidate> population(populationSize),trial(populationSize);
  std::vector<SearchMetrics> populationMetrics(populationSize),trialMetrics(populationSize);
  Candidate seed{-0.51864396687563019,0.12700737424508893,0.071711183124306049,0.0,0.0,
    -0.41551333272698654,0.47383274947983156,-0.34817367193977911,0.0,0.0,
    0.70,1.28,1.0,0.75};
  population[0]=seed;
  Candidate alternative{-0.4837821910492075,0.09560497747537516,0.06124524274815164,0.0,0.0,
    -0.45997384695584637,0.5803069930800483,-0.39048803622838424,0.0,0.0,
    0.70,1.28,1.0,0.75};
  population[1]=alternative;
  Candidate envelopeFit{-0.4387207637325457,-0.47031798438283273,0.45853381730287407,
    -0.3832106409205932,0.47516460967737623,-0.30504721604748886,
    -0.32507474330689717,-0.73498131027527891,-0.060570593449980392,
    1.069397303645947,0.70,1.28,1.0,0.75};
  population[2]=envelopeFit;
  Candidate fullFlowEnvelopeFit{-0.52164358694220492,-0.024827807054374828,
    1.0452244736858278,-1.3656657598369506,0.54300716429149642,
    -0.38685452062281694,0.31457715720002261,-0.61976108254554241,
    1.1664678885754221,-0.73680735226233729,0.63199298223094058,
    1.2792974481304009,0.60979123180328654,0.80193065649632189};
  population[3]=fullFlowEnvelopeFit;
  std::mt19937_64 random(20260724);std::uniform_real_distribution<double> uniform(0.0,1.0);
  for(int index=4;index<populationSize;++index)for(std::size_t parameter=0;parameter<seed.size();++parameter) {
    const double span=upper[parameter]-lower[parameter];
    const double center=index%3==0?seed[parameter]:
      (index%3==1?alternative[parameter]:envelopeFit[parameter]);
    const double radius=parameter<10?0.08*span:0.30*span;
    population[index][parameter]=std::clamp(center+(uniform(random)-0.5)*radius,
                                             lower[parameter],upper[parameter]);
  }
  #pragma omp parallel for
  for(int index=0;index<populationSize;++index)
    populationMetrics[index]=evaluateCandidate(population[index],targets);
  for(int generation=0;generation<generations;++generation) {
    for(int index=0;index<populationSize;++index) {
      int a=0,b=0,c=0;
      do a=static_cast<int>(uniform(random)*populationSize)%populationSize;while(a==index);
      do b=static_cast<int>(uniform(random)*populationSize)%populationSize;while(b==index||b==a);
      do c=static_cast<int>(uniform(random)*populationSize)%populationSize;while(c==index||c==a||c==b);
      const int forced=static_cast<int>(uniform(random)*seed.size())%seed.size();
      for(std::size_t parameter=0;parameter<seed.size();++parameter) {
        double value=(uniform(random)<0.90||static_cast<int>(parameter)==forced)?
          population[a][parameter]+0.65*(population[b][parameter]-population[c][parameter]):
          population[index][parameter];
        if(value<lower[parameter]||value>upper[parameter])
          value=lower[parameter]+uniform(random)*(upper[parameter]-lower[parameter]);
        trial[index][parameter]=value;
      }
    }
    #pragma omp parallel for
    for(int index=0;index<populationSize;++index)
      trialMetrics[index]=evaluateCandidate(trial[index],targets);
    for(int index=0;index<populationSize;++index)
      if(trialMetrics[index].score<populationMetrics[index].score) {
        population[index]=trial[index];populationMetrics[index]=trialMetrics[index];
      }
    if(generation%20==0||generation+1==generations) {
      const auto best=std::min_element(populationMetrics.begin(),populationMetrics.end(),
        [](const auto& a,const auto& b){return a.score<b.score;});
      const auto index=static_cast<std::size_t>(best-populationMetrics.begin());
      std::cout<<"generation="<<generation<<' ';printCandidate(population[index],*best);
    }
  }
  std::vector<int> order(populationSize);
  for(int index=0;index<populationSize;++index)order[index]=index;
  std::sort(order.begin(),order.end(),[&](int a,int b){
    return populationMetrics[a].score<populationMetrics[b].score;});
  std::cout<<"FINAL_CANDIDATES\n";
  for(int rank=0;rank<std::min(32,populationSize);++rank)
    printCandidate(population[order[rank]],populationMetrics[order[rank]]);
}

Prediction predict(const Coefficients& coefficients,double previous,double outer,
                   double power,double split,double fraction,const std::vector<Target>& targets) {
  if(!(previous>0.0&&previous<outer&&outer<hevea::sphere::yInfinity&&
       split>0.0&&split<1.0&&fraction>0.0&&fraction<1.0))
    throw std::runtime_error("invalid model parameters");
  struct Request {double y;std::size_t index;bool north;};
  std::vector<Request> requests;
  requests.reserve(2*targets.size());
  for(std::size_t k=0;k<targets.size();++k) {
    requests.push_back({-targets[k].latitude,k,false});
    requests.push_back({targets[k].latitude,k,true});
  }
  std::sort(requests.begin(),requests.end(),[](const auto& a,const auto& b){return a.y<b.y;});
  Prediction prediction;
  prediction.north.resize(targets.size());prediction.south.resize(targets.size());
  State state;double y=-outer;
  constexpr double maximumStep=0.00075;
  for(const auto& request:requests) {
    if(request.y<=-outer||request.y>=outer) {
      (request.north?prediction.north:prediction.south)[request.index]=0.0;
      continue;
    }
    while(y+maximumStep<request.y) {
      state=rk4Step(y,maximumStep,state,coefficients,previous,outer,power,split,fraction);
      y+=maximumStep;
    }
    if(y<request.y) {
      state=rk4Step(y,request.y-y,state,coefficients,previous,outer,power,split,fraction);
      y=request.y;
    }
    (request.north?prediction.north:prediction.south)[request.index]=
      amplitude(state,request.y,previous,outer,power,split);
  }
  double squareSum=0.0;std::size_t count=0;
  for(std::size_t k=0;k<targets.size();++k) {
    const std::array pairs{std::pair{prediction.north[k],targets[k].north},
                           std::pair{prediction.south[k],targets[k].south}};
    for(const auto& pair:pairs) {
      const double error=pair.first-pair.second;
      squareSum+=error*error;prediction.maximum=std::max(prediction.maximum,std::abs(error));++count;
    }
    prediction.symmetryMaximum=std::max(prediction.symmetryMaximum,
      std::abs(prediction.north[k]-prediction.south[k]));
  }
  prediction.rms=std::sqrt(squareSum/count);
  return prediction;
}

} // namespace

int main(int argc,char** argv) {
  try {
    if(argc==5&&std::string(argv[1])=="--search") {
      const auto targets=readTargets(argv[2]);
      search(targets,std::max(1,std::atoi(argv[3])),std::max(32,std::atoi(argv[4])));
      return 0;
    }
    if(argc!=17) {
      std::cerr<<"usage: stage1_envelope_model REFERENCE.csv C0 C1 C2 C3 C4 C5 C6 C7 C8 C9 "
                 "PREVIOUS OUTER POWER SPLIT FRACTION\n";
      return 2;
    }
    Coefficients coefficients{};
    for(std::size_t k=0;k<coefficients.size();++k)coefficients[k]=std::stod(argv[k+2]);
    const double previous=std::stod(argv[12]),outer=std::stod(argv[13]);
    const double power=std::stod(argv[14]),split=std::stod(argv[15]);
    const double fraction=std::stod(argv[16]);
    const auto targets=readTargets(argv[1]);
    const auto profile=profileMetrics(coefficients,previous,outer,power,split);
    const auto prediction=predict(coefficients,previous,outer,power,split,fraction,targets);
    std::cout<<std::setprecision(17)
      <<"stage1_envelope_model=pass rms="<<prediction.rms
      <<" max="<<prediction.maximum
      <<" symmetry_max="<<prediction.symmetryMaximum
      <<" previous="<<previous<<" outer="<<outer<<" power="<<power
      <<" split="<<split<<" fraction="<<fraction<<'\n';
    std::cout<<"profile mean="<<profile.mean<<" max="<<profile.maximum
      <<" radius="<<profile.radius<<" global_rms="<<profile.globalRms
      <<" global_max="<<profile.globalMax<<" global_derivative_rms="
      <<profile.globalDerivativeRms<<" global_derivative_max="
      <<profile.globalDerivativeMax<<" min_primitive="<<profile.minPrimitive[0]<<','
      <<profile.minPrimitive[1]<<','<<profile.minPrimitive[2]<<" min_wy="
      <<profile.minWy[0]<<','<<profile.minWy[1]<<','<<profile.minWy[2]
      <<" penalty="<<profile.penalty<<'\n';
    std::cout<<"latitude,referenceNorth,predictedNorth,referenceSouth,predictedSouth\n";
    for(std::size_t k=0;k<targets.size();++k)
      std::cout<<targets[k].latitude<<','<<targets[k].north<<','<<prediction.north[k]<<','
               <<targets[k].south<<','<<prediction.south[k]<<'\n';
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"stage1_envelope_model: "<<error.what()<<'\n';return 2;
  }
}
