/*
 * Reduced-sphere convex-integration experiment.
 * The generic corrugation and numerical-integration machinery adapts ideas
 * and code patterns from https://github.com/HeveaProject/Hevea (GPLv3).
 * Sphere-specific profile, characteristic-flow and validation code is this
 * project's implementation of Bartzos et al. (2017). See NOTICE and LICENSE.
 */
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "src/reduced_sphere_profile.hpp"
#include "src/hevea_numeric.hpp"

constexpr double pi = 3.1415926535897932384626433832795;

struct V3 {
  double x = 0, y = 0, z = 0;
  V3 operator+(V3 b) const { return {x+b.x,y+b.y,z+b.z}; }
  V3 operator-(V3 b) const { return {x-b.x,y-b.y,z-b.z}; }
  V3 operator*(double s) const { return {x*s,y*s,z*s}; }
  V3 operator/(double s) const { return *this*(1.0/s); }
  V3& operator+=(V3 b) { x+=b.x; y+=b.y; z+=b.z; return *this; }
};
V3 operator*(double s,V3 a){return a*s;}
double dot(V3 a,V3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
V3 cross(V3 a,V3 b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
double norm(V3 a){return std::sqrt(std::max(0.0,dot(a,a)));}
V3 unit(V3 a){double n=norm(a); return n>1e-14?a/n:V3{0,0,1};}

struct M2 { double e=0,f=0,g=0; };
double eval(M2 m,double x,double y){return m.e*x*x+2*m.f*x*y+m.g*y*y;}
double bilinear(M2 m,double ax,double ay,double bx,double by){
  return m.e*ax*bx+m.f*(ax*by+ay*bx)+m.g*ay*by;
}

struct Grid {
  int nx,ny; double ymax,dx,dy; std::vector<V3> p;
  Grid(int nx_,int ny_,double ymax_):nx(nx_),ny(ny_),ymax(ymax_),
    dx(2*pi/nx_),dy(2*ymax_/(ny_-1)),p(static_cast<size_t>(nx_)*ny_){}
  V3& at(int i,int j){i=(i%nx+nx)%nx; return p[static_cast<size_t>(j)*nx+i];}
  V3 at(int i,int j)const{i=(i%nx+nx)%nx; return p[static_cast<size_t>(j)*nx+i];}
  double y(int j)const{return -ymax+j*dy;}
};

hevea::CylinderCubicJet<V3> sampleJet(Grid const& q,double x,double y){
  return hevea::cylinderCubicJet<V3>(
    [&](int i,int j){return q.at(i,j);},q.nx,q.ny,q.dx,q.dy,-q.ymax,x,y);
}

double smooth(double t){t=std::clamp(t,0.0,1.0);return t*t*t*(10+t*(-15+6*t));}
double poweredStep(double t,bool lambda,int direction){
  t=std::clamp(t,0.0,1.0);
  static const std::array<double,3> lambdaPowers=[] {
    std::array<double,3> powers{0.60979123180328654,0.0,0.0};
    if(const char* value=std::getenv("HEVEA_LAMBDA_POWERS")){
      std::stringstream input(value);char comma=0;
      for(size_t k=0;k<powers.size();k++)if(!(input>>powers[k])||
          (k+1<powers.size()&&(!(input>>comma)||comma!=',')))
        throw std::runtime_error("HEVEA_LAMBDA_POWERS requires three comma-separated numbers");
    }else if(const char* value=std::getenv("HEVEA_LAMBDA_POWER"))powers.fill(std::stod(value));
    return powers;}();
  static const std::array<double,3> chiPowers=[] {
    std::array<double,3> powers{0.60979123180328654,0.0,0.0};
    if(const char* value=std::getenv("HEVEA_CHI_POWERS")){
      std::stringstream input(value);char comma=0;
      for(size_t k=0;k<powers.size();k++)if(!(input>>powers[k])||
          (k+1<powers.size()&&(!(input>>comma)||comma!=',')))
        throw std::runtime_error("HEVEA_CHI_POWERS requires three comma-separated numbers");
    }else if(const char* value=std::getenv("HEVEA_CHI_POWER"))powers.fill(std::stod(value));
    return powers;}();
  const double power=(lambda?lambdaPowers:chiPowers)[static_cast<size_t>(direction)];
  if(power==0.0)return smooth(t);
  if(!(std::isfinite(power)&&power>1.0/3.0))throw std::runtime_error(
    std::string(lambda?"HEVEA_LAMBDA_POWER":"HEVEA_CHI_POWER")+
    " must be finite and greater than one third");
  if(t==0.0||t==1.0)return t;
  const double plateau=std::clamp(smooth(t),0.0,1.0);
  const double left=std::pow(plateau,power),right=std::pow(1.0-plateau,power);
  return left/(left+right);
}
double lambdaStep(double t,int direction){return poweredStep(t,true,direction);}
double chiStep(double t,int direction){return poweredStep(t,false,direction);}
double transitionSplit(int direction){
  static const std::array<double,3> splits=[] {
    std::array<double,3> values{0.80193065649632189,0.5,0.5};
    if(const char* value=std::getenv("HEVEA_TRANSITION_SPLITS")){
      std::stringstream input(value);char comma=0;
      for(size_t k=0;k<values.size();k++)if(!(input>>values[k])||
          (k+1<values.size()&&(!(input>>comma)||comma!=',')))
        throw std::runtime_error("HEVEA_TRANSITION_SPLITS requires three comma-separated numbers");
    }else if(const char* value=std::getenv("HEVEA_TRANSITION_SPLIT"))values.fill(std::stod(value));
    if(!std::all_of(values.begin(),values.end(),[](double split){
         return std::isfinite(split)&&split>0.0&&split<1.0;}))
      throw std::runtime_error("HEVEA_TRANSITION_SPLITS must lie strictly between zero and one");
    return values;}();
  return splits[static_cast<size_t>(direction)];
}
double wrap(double x){x=std::fmod(x,2*pi);return x<0?x+2*pi:x;}

V3 sample(Grid const& q,double x,double y){
  x=wrap(x); double gx=x/q.dx, gy=(y+q.ymax)/q.dy;
  int i=static_cast<int>(std::floor(gx)); int j=std::clamp(static_cast<int>(std::floor(gy)),0,q.ny-2);
  double a=gx-i,b=std::clamp(gy-j,0.0,1.0);
  V3 p00=q.at(i,j),p10=q.at(i+1,j),p01=q.at(i,j+1),p11=q.at(i+1,j+1);
  return (1-a)*(1-b)*p00+a*(1-b)*p10+(1-a)*b*p01+a*b*p11;
}

void derivatives(Grid const& q,std::vector<V3>& fx,std::vector<V3>& fy){
  fx.resize(q.p.size());fy.resize(q.p.size());
  #pragma omp parallel for
  for(int j=0;j<q.ny;j++)for(int i=0;i<q.nx;i++){
    size_t k=static_cast<size_t>(j)*q.nx+i;
    fx[k]=(q.at(i-2,j)-8.0*q.at(i-1,j)+8.0*q.at(i+1,j)-q.at(i+2,j))/(12*q.dx);
    if(j==0)fy[k]=(-25.0*q.at(i,0)+48.0*q.at(i,1)-36.0*q.at(i,2)+
                    16.0*q.at(i,3)-3.0*q.at(i,4))/(12*q.dy);
    else if(j==1)fy[k]=(-3.0*q.at(i,0)-10.0*q.at(i,1)+18.0*q.at(i,2)-
                         6.0*q.at(i,3)+q.at(i,4))/(12*q.dy);
    else if(j==q.ny-2)fy[k]=(3.0*q.at(i,q.ny-1)+10.0*q.at(i,q.ny-2)-
                              18.0*q.at(i,q.ny-3)+6.0*q.at(i,q.ny-4)-
                              q.at(i,q.ny-5))/(12*q.dy);
    else if(j==q.ny-1)fy[k]=(25.0*q.at(i,q.ny-1)-48.0*q.at(i,q.ny-2)+
                              36.0*q.at(i,q.ny-3)-16.0*q.at(i,q.ny-4)+
                              3.0*q.at(i,q.ny-5))/(12*q.dy);
    else fy[k]=(q.at(i,j-2)-8.0*q.at(i,j-1)+8.0*q.at(i,j+1)-q.at(i,j+2))/(12*q.dy);
  }
}

template<class T>T hermite(double t,T const&p0,T const&p1,T const&v0,T const&v1){
  double h01=t*t*(3-2*t),h10=t*(t*(t-2)+1),h11=t*t*(t-1);
  return p0+h01*(p1-p0)+h10*v0+h11*v1;
}

template<class T,class Getter>T periodicCubicSample(Getter get,int nx,int ny,double dx,double dy,
                                                    double ymax,double x,double y){
  x=wrap(x);double gx=x/dx,gy=(y+ymax)/dy;int i=static_cast<int>(std::floor(gx));
  int j=std::clamp(static_cast<int>(std::floor(gy)),0,ny-2);
  double u=gx-i,v=std::clamp(gy-j,0.0,1.0);
  int jAux=std::clamp(j,2,ny-4);std::array<T,6> row{};
  auto periodic=[&](int ii,int jj){ii=(ii%nx+nx)%nx;return get(ii,jj);};
  for(int k=0;k<6;k++){
    int jj=jAux+k-2;
    T vi=(periodic(i-2,jj)-periodic(i+2,jj)+8.0*(periodic(i+1,jj)-periodic(i-1,jj)))/12.0;
    T vi1=(periodic(i-1,jj)-periodic(i+3,jj)+8.0*(periodic(i+2,jj)-periodic(i,jj)))/12.0;
    row[k]=hermite(u,periodic(i,jj),periodic(i+1,jj),vi,vi1);
  }
  int offset=j-jAux,index=2+offset;
  auto derivative=[&](int at){
    if(at==0&&offset==0)return (row[0]-row[4]+8.0*(row[3]-row[1]))/12.0;
    if(at==1&&offset==0)return (row[1]-row[5]+8.0*(row[4]-row[2]))/12.0;
    int base=std::clamp(index+at,1,4);
    return (row[base+1]-row[base-1])/2.0;
  };
  return hermite(v,row[index],row[index+1],derivative(0),derivative(1));
}

template<class T>T periodicSample(std::vector<T> const& a,int nx,int ny,double dx,double dy,double ymax,double x,double y){
  return periodicCubicSample<T>([&](int i,int j){return a[static_cast<size_t>(j)*nx+i];},
                                nx,ny,dx,dy,ymax,x,y);
}

std::array<double,2> periodicSampleW(std::vector<std::array<double,2>> const& a,int nx,int ny,
                                     double dx,double dy,double ymax,double x,double y){
  std::array<double,2> out{};
  for(int c=0;c<2;c++)out[c]=periodicCubicSample<double>(
    [&](int i,int j){return a[static_cast<size_t>(j)*nx+i][c];},nx,ny,dx,dy,ymax,x,y);
  return out;
}

double invJ0(double value){
  value=std::clamp(value,0.0,1.0); double lo=0,hi=2.404825557695773;
  for(int k=0;k<48;k++){double m=(lo+hi)/2;if(std::cyl_bessel_j(0,m)>value)lo=m;else hi=m;}
  return (lo+hi)/2;
}

Grid initialSurface(int nx,int ny,double ymax,double ballRadius,double eta){
  if(std::abs(ymax-hevea::sphere::yInfinity)>1e-12 ||
     std::abs(ballRadius-hevea::sphere::certifiedRadius)>1e-12 ||
     std::abs(eta-hevea::sphere::certifiedEta)>1e-12)
    throw std::runtime_error("uncertified initial profile parameters");
  Grid q(nx,ny,ymax);auto coefficients=hevea::sphere::freeCoefficients;
  auto extendedCoefficients=hevea::sphere::extendedFreeCoefficients;
  const bool useExtended=std::getenv("HEVEA_PROFILE_COEFFICIENTS")==nullptr||
    std::getenv("HEVEA_EXTENDED_PROFILE_COEFFICIENTS")!=nullptr;
  if(const char* value=std::getenv("HEVEA_EXTENDED_PROFILE_COEFFICIENTS")){
    std::stringstream input(value);char comma=0;
    for(size_t k=0;k<extendedCoefficients.size();k++){
      if(!(input>>extendedCoefficients[k])||(k+1<extendedCoefficients.size()&&
          (!(input>>comma)||comma!=',')))
        throw std::runtime_error(
          "HEVEA_EXTENDED_PROFILE_COEFFICIENTS requires ten comma-separated numbers");
    }
  }
  if(const char* value=std::getenv("HEVEA_PROFILE_COEFFICIENTS")){
    std::stringstream input(value);char comma=0;
    for(size_t k=0;k<coefficients.size();k++){
      if(!(input>>coefficients[k])||(k+1<coefficients.size()&&(!(input>>comma)||comma!=',')))
        throw std::runtime_error("HEVEA_PROFILE_COEFFICIENTS requires six comma-separated numbers");
    }
  }
  #pragma omp parallel for
  for(int j=0;j<ny;j++)for(int i=0;i<nx;i++){
    double y=q.y(j),x=i*q.dx;auto p=useExtended?
      hevea::sphere::evaluate(y,extendedCoefficients):hevea::sphere::evaluate(y,coefficients);
    q.at(i,j)={p.x*std::cos(x),p.x*std::sin(x),p.z};
  }
  return q;
}

struct Fields {std::vector<V3> fx,fy,n;std::vector<M2> metric;};
Fields makeFields(Grid const& q){
  Fields z;derivatives(q,z.fx,z.fy);z.n.resize(q.p.size());z.metric.resize(q.p.size());
  #pragma omp parallel for
  for(size_t k=0;k<q.p.size();k++){
    z.n[k]=unit(cross(z.fx[k],z.fy[k]));
    z.metric[k]={dot(z.fx[k],z.fx[k]),dot(z.fx[k],z.fy[k]),dot(z.fy[k],z.fy[k])};
  }
  return z;
}

struct StepStats {double minRho=1e100,minWy=1e100,maxSlope=0,minSpacing=1e100;
  double minRhoX=0,minRhoY=0;
  double minJacobian=1e100,inverseRoundtrip=0,velocityRelativeError=0,phaseClosureError=0;
  double seamPositionGap=0,seamDerivativeGap=0,targetMean=0,targetMax=0;};

Grid corrugate(Grid const& old,Grid const& initial,int direction,int ridgeCount,double previous,double outer,
               double targetFraction,double correctionScale,double phaseOffset,StepStats& stats){
  Fields fld=makeFields(old),initialFields=makeFields(initial);int nx=old.nx,ny=old.ny;
  double invSqrt2=1/std::sqrt(2.0),lx,ly,vx,vy,N;
  if(direction==0){lx=invSqrt2;ly=invSqrt2;vx=-invSqrt2;vy=invSqrt2;N=ridgeCount*std::sqrt(2.0)/(2*pi);}
  else if(direction==1){lx=-invSqrt2;ly=invSqrt2;vx=-invSqrt2;vy=-invSqrt2;N=ridgeCount*std::sqrt(2.0)/(2*pi);}
  else{lx=0;ly=1;vx=-1;vy=0;N=ridgeCount/(2*outer);}

  double minRho=1e100,minWy=1e100,maxSlope=0;int invalidPrimitive=0;
  #pragma omp parallel for reduction(min:minRho,minWy) reduction(max:maxSlope,invalidPrimitive)
  for(int j=0;j<ny;j++)for(int i=0;i<nx;i++){
    size_t k=static_cast<size_t>(j)*nx+i;double y=old.y(j),ay=std::abs(y);
    double transition=previous+transitionSplit(direction)*(outer-previous);
    double lambda=ay<=previous?1.0:(ay>=transition?0.0:
      1-lambdaStep((ay-previous)/(transition-previous),direction));
    M2 m=fld.metric[k],m0=initialFields.metric[k],round{std::cos(y)*std::cos(y),0,1};
    M2 targetFull{m0.e+targetFraction*(round.e-m0.e),
                  m0.f+targetFraction*(round.f-m0.f),
                  m0.g+targetFraction*(round.g-m0.g)};
    M2 desired{m0.e+lambda*(targetFull.e-m0.e),m0.f+lambda*(targetFull.f-m0.f),
               m0.g+lambda*(targetFull.g-m0.g)};
    M2 defect{desired.e-m.e,desired.f-m.f,desired.g-m.g};
    double rho=direction==0?defect.e+defect.f:(direction==1?defect.e-defect.f:defect.g-defect.e);
    if(ay<=previous)minRho=std::min(minRho,rho);
    if(correctionScale>0&&ay<=previous&&rho<=1e-6)invalidPrimitive=1;
    rho*=correctionScale;
    M2 mu{m.e+rho*lx*lx,m.f+rho*lx*ly,m.g+rho*ly*ly};
    double ux=lx/(lx*lx+ly*ly),uy=ly/(lx*lx+ly*ly);
    double z=-bilinear(mu,ux,uy,vx,vy)/std::max(1e-14,eval(mu,vx,vy));
    double wx=ux+z*vx,wy=uy+z*vy;
    if(ay<=outer){minWy=std::min(minWy,wy);
      maxSlope=std::max(maxSlope,std::abs(wx/std::max(std::abs(wy),1e-14)));}
  }
  for(int j=0;j<ny;j++)for(int i=0;i<nx;i++){
    size_t k=static_cast<size_t>(j)*nx+i;double y=old.y(j),ay=std::abs(y);
    if(ay>previous)continue;
    M2 m=fld.metric[k],m0=initialFields.metric[k],round{std::cos(y)*std::cos(y),0,1};
    M2 targetFull{m0.e+targetFraction*(round.e-m0.e),m0.f+targetFraction*(round.f-m0.f),m0.g+targetFraction*(round.g-m0.g)};
    M2 defect{targetFull.e-m.e,targetFull.f-m.f,targetFull.g-m.g};
    double rho=direction==0?defect.e+defect.f:(direction==1?defect.e-defect.f:defect.g-defect.e);
    if(rho<=stats.minRho){stats.minRho=rho;stats.minRhoX=i*old.dx;stats.minRhoY=y;}
  }
  if(invalidPrimitive)throw std::runtime_error("non-positive primitive coordinate before corrugation: min="+
    std::to_string(minRho)+" x="+std::to_string(stats.minRhoX)+" y="+std::to_string(stats.minRhoY));
  stats.minWy=minWy;stats.maxSlope=maxSlope;
  stats.phaseClosureError=direction<2?
    std::abs(2*pi*N*lx-std::round(2*pi*N*lx)):
    std::abs(2*outer*N-ridgeCount);
  if(const char* value=std::getenv("HEVEA_STOP_BEFORE_FLOW_DIRECTION"))
    if(std::atoi(value)==direction)return old;
  if(correctionScale==0.0){
    stats.minSpacing=old.dx;stats.minJacobian=1.0;
    return old;
  }

  int j0=std::max(0,static_cast<int>(std::ceil((old.ymax-outer)/old.dy)));
  int j1=std::min(ny-1,static_cast<int>(std::floor((old.ymax+outer)/old.dy)));
  Grid raw=old,flowIncrement(nx,ny,old.ymax);
  std::vector<double> flowX(old.p.size());
  struct RhsValue {double flowDerivative;V3 incrementDerivative;double targetSpeed;double measuredSpeed;};
  auto rhs=[&](double x,double y){
    auto currentJet=sampleJet(old,x,y),initialJet=sampleJet(initial,x,y);
    V3 fx=currentJet.derivativeX,fy=currentJet.derivativeY;
    M2 m{dot(fx,fx),dot(fx,fy),dot(fy,fy)};
    M2 m0{dot(initialJet.derivativeX,initialJet.derivativeX),
          dot(initialJet.derivativeX,initialJet.derivativeY),
          dot(initialJet.derivativeY,initialJet.derivativeY)};
    double ay=std::abs(y),transition=previous+transitionSplit(direction)*(outer-previous);
    double lambda=ay<=previous?1.0:(ay>=transition?0.0:
      1-lambdaStep((ay-previous)/(transition-previous),direction));
    M2 round{std::cos(y)*std::cos(y),0,1};
    M2 desired{m0.e+lambda*targetFraction*(round.e-m0.e),
               m0.f+lambda*targetFraction*(round.f-m0.f),
               m0.g+lambda*targetFraction*(round.g-m0.g)};
    M2 defect{desired.e-m.e,desired.f-m.f,desired.g-m.g};
    double rho=(direction==0?defect.e+defect.f:
      (direction==1?defect.e-defect.f:defect.g-defect.e))*correctionScale;
    M2 mu{m.e+rho*lx*lx,m.f+rho*lx*ly,m.g+rho*ly*ly};
    double ux=lx/(lx*lx+ly*ly),uy=ly/(lx*lx+ly*ly);
    double z=-bilinear(mu,ux,uy,vx,vy)/std::max(1e-14,eval(mu,vx,vy));
    double wx=ux+z*vx,wy=uy+z*vy;
    if(wy<1e-8){wx=-wx;wy=-wy;}
    double rr=std::sqrt(std::max(1e-14,eval(mu,wx,wy)));
    V3 oldDerivative=fx*wx+fy*wy;
    double aa=invJ0(std::clamp(norm(oldDerivative)/rr,0.0,1.0));
    V3 nn=unit(cross(fx,fy));
    V3 tangent=unit(fx*wx+fy*wy);double phase=lx*x+ly*y;
    double theta=aa*std::cos(2*pi*N*phase+phaseOffset);
    V3 newDerivative=rr*(std::cos(theta)*tangent+std::sin(theta)*nn);
    // Integrating the increment is algebraically equivalent to integrating the
    // complete corrugated map, but it avoids re-interpolating all high-frequency
    // geometry inherited from earlier stages.  In particular, rho=0 gives an
    // exactly zero increment and therefore an exact discrete identity step.
    return RhsValue{wx/wy,(newDerivative-oldDerivative)/wy,rr,norm(newDerivative)};
  };
  double velocityRelativeError=0;
  const bool equatorAnchor=(direction==2&&
    std::getenv("HEVEA_SOUTH_STAGE3_ANCHOR")==nullptr)||
    std::getenv("HEVEA_EQUATOR_ANCHOR_ALL")!=nullptr;
  const bool centerIncrement=!equatorAnchor&&direction<2&&
    std::getenv("HEVEA_CENTER_INCREMENT")!=nullptr;
  #pragma omp parallel for
  for(int i=0;i<nx;i++){
    double localVelocityRelativeError=0;
    auto derivative=[&](double x,double y,double sign){
      auto value=rhs(x,y);
      localVelocityRelativeError=std::max(localVelocityRelativeError,
        std::abs(value.measuredSpeed-value.targetSpeed)/
        std::max(1e-14,value.targetSpeed));
      return hevea::State<4>{sign*value.flowDerivative,sign*value.incrementDerivative.x,
        sign*value.incrementDerivative.y,sign*value.incrementDerivative.z};
    };
    if(equatorAnchor){
      const int anchor=std::clamp(static_cast<int>(std::llround(old.ymax/old.dy)),j0,j1);
      hevea::State<4> state{i*old.dx,0,0,0};
      flowX[static_cast<size_t>(anchor)*nx+i]=state[0];flowIncrement.at(i,anchor)={};
      for(int j=anchor;j<j1;j++){
        const double y=old.y(j),h=old.dy;
        auto forward=[&](double tau,hevea::State<4> const& z){
          return derivative(z[0],y+tau,1.0);};
        state=hevea::integrateDp54<4>(forward,0,state,h,1e-9,1e-11);
        const size_t k=static_cast<size_t>(j+1)*nx+i;flowX[k]=state[0];
        flowIncrement.at(i,j+1)={state[1],state[2],state[3]};
      }
      state={i*old.dx,0,0,0};
      for(int j=anchor;j>j0;j--){
        const double y=old.y(j),h=old.dy;
        auto backward=[&](double tau,hevea::State<4> const& z){
          return derivative(z[0],y-tau,-1.0);};
        state=hevea::integrateDp54<4>(backward,0,state,h,1e-9,1e-11);
        const size_t k=static_cast<size_t>(j-1)*nx+i;flowX[k]=state[0];
        flowIncrement.at(i,j-1)={state[1],state[2],state[3]};
      }
    }else{
      hevea::State<4> state{i*old.dx,0,0,0};
      flowX[static_cast<size_t>(j0)*nx+i]=state[0];flowIncrement.at(i,j0)={};
      for(int j=j0;j<j1;j++){
        const double y=old.y(j),h=old.dy;
        auto forward=[&](double tau,hevea::State<4> const& z){
          auto value=rhs(z[0],y+tau);
        localVelocityRelativeError=std::max(localVelocityRelativeError,
          std::abs(value.measuredSpeed-value.targetSpeed)/
          std::max(1e-14,value.targetSpeed));
        return hevea::State<4>{value.flowDerivative,value.incrementDerivative.x,
          value.incrementDerivative.y,value.incrementDerivative.z};};
        state=hevea::integrateDp54<4>(forward,0,state,h,1e-9,1e-11);
        const size_t k=static_cast<size_t>(j+1)*nx+i;flowX[k]=state[0];
        flowIncrement.at(i,j+1)={state[1],state[2],state[3]};
      }
      if(centerIncrement){
        const V3 halfTotal{state[1]/2,state[2]/2,state[3]/2};
        for(int j=j0;j<=j1;j++)flowIncrement.at(i,j)=flowIncrement.at(i,j)-halfTotal;
      }
    }
    #pragma omp critical
    velocityRelativeError=std::max(velocityRelativeError,localVelocityRelativeError);
  }
  stats.velocityRelativeError=velocityRelativeError;
  for(int j=j0;j<=j1;j++)for(int i=0;i<nx;i++){
    const size_t k=static_cast<size_t>(j)*nx+i;
    const V3 increment=flowIncrement.at(i,j);
    if(!std::isfinite(flowX[k])||!std::isfinite(increment.x)||
       !std::isfinite(increment.y)||!std::isfinite(increment.z))
      throw std::runtime_error("non-finite forward flow at x="+
        std::to_string(i*old.dx)+" y="+std::to_string(old.y(j)));
  }
  for(int j=j0;j<=j1;j++){
    size_t row=static_cast<size_t>(j)*nx;double first=flowX[row];
    double last=first+2*pi;
    for(int i=0;i<nx;i++){
      double target=i*old.dx;
      while(target<first)target+=2*pi;
      while(target>=last)target-=2*pi;
      auto inverse=hevea::periodicMonotoneCubicSample<V3>(nx,2*pi,target,
        [&](int index){return flowX[row+index];},
        [&](int index){return flowIncrement.at(index,j);});
      stats.minSpacing=std::min(stats.minSpacing,inverse.cellSpacing);
      stats.inverseRoundtrip=std::max(stats.inverseRoundtrip,inverse.linearRoundTripError);
      raw.at(i,j)=old.at(i,j)+inverse.value;
    }
  }
  stats.minJacobian=stats.minSpacing/old.dx;
  if(!(stats.minSpacing>0))throw std::runtime_error("folded forward flow map");

  Grid out=old;
  #pragma omp parallel for
  for(int j=0;j<ny;j++)for(int i=0;i<nx;i++){
    double ay=std::abs(old.y(j));
    double transition=previous+transitionSplit(direction)*(outer-previous);
    double keep=ay<=transition?0.0:(ay>=outer?1.0:
      chiStep((ay-transition)/(outer-transition),direction));
    out.at(i,j)=(1-keep)*raw.at(i,j)+keep*old.at(i,j);
  }
  for(double seamY:{-outer,outer})for(int i=0;i<nx;i++){
    const double x=i*old.dx;
    const auto outJet=sampleJet(out,x,seamY),oldJet=sampleJet(old,x,seamY);
    stats.seamPositionGap=std::max(stats.seamPositionGap,norm(outJet.value-oldJet.value));
    stats.seamDerivativeGap=std::max({stats.seamDerivativeGap,
      norm(outJet.derivativeX-oldJet.derivativeX),
      norm(outJet.derivativeY-oldJet.derivativeY)});
  }
  double targetSum=0,targetMax=0;long long targetCount=0;
  #pragma omp parallel for reduction(+:targetSum,targetCount) reduction(max:targetMax)
  for(int j=1;j<ny-1;j++)for(int i=0;i<nx;i++){
    size_t k=static_cast<size_t>(j)*nx+i;double y=old.y(j),ay=std::abs(y);
    M2 m0=initialFields.metric[k],round{std::cos(y)*std::cos(y),0,1};
    double transition=previous+transitionSplit(direction)*(outer-previous);
    double lambda=ay<=previous?1.0:(ay>=transition?0.0:
      1-lambdaStep((ay-previous)/(transition-previous),direction));
    M2 desired{m0.e+lambda*targetFraction*(round.e-m0.e),
      m0.f+lambda*targetFraction*(round.f-m0.f),
      m0.g+lambda*targetFraction*(round.g-m0.g)};
    V3 ox=(out.at(i-2,j)-8.0*out.at(i-1,j)+8.0*out.at(i+1,j)-out.at(i+2,j))/(12*out.dx);
    V3 oy;
    if(j==1)oy=(-3.0*out.at(i,0)-10.0*out.at(i,1)+18.0*out.at(i,2)-
                 6.0*out.at(i,3)+out.at(i,4))/(12*out.dy);
    else if(j==out.ny-2)oy=(3.0*out.at(i,out.ny-1)+10.0*out.at(i,out.ny-2)-
                              18.0*out.at(i,out.ny-3)+6.0*out.at(i,out.ny-4)-
                              out.at(i,out.ny-5))/(12*out.dy);
    else oy=(out.at(i,j-2)-8.0*out.at(i,j-1)+8.0*out.at(i,j+1)-out.at(i,j+2))/(12*out.dy);
    M2 actual{dot(ox,ox),dot(ox,oy),dot(oy,oy)};
    double error=std::sqrt((actual.e-desired.e)*(actual.e-desired.e)+
      2*(actual.f-desired.f)*(actual.f-desired.f)+(actual.g-desired.g)*(actual.g-desired.g));
    targetSum+=error;targetMax=std::max(targetMax,error);targetCount++;
  }
  stats.targetMean=targetSum/targetCount;stats.targetMax=targetMax;
  return out;
}

struct Report {double roundMean=0,roundMax=0,targetMean=0,targetMax=0,bound=0;
  double roundMaxX=0,roundMaxY=0,roundMaxE=0,roundMaxF=0,roundMaxG=0;
  double minRho1=1e100,minRho2=1e100,minRho3=1e100,minRho3Y=0;};
Report report(Grid const&q,Grid const&initial,double targetFraction){
  Fields f=makeFields(q),f0=makeFields(initial);Report r;double sum=0,tsum=0;size_t count=0;
  for(int j=1;j<q.ny-1;j++)for(int i=0;i<q.nx;i++){
    size_t k=static_cast<size_t>(j)*q.nx+i;double y=q.y(j);M2 m=f.metric[k];
    double re=std::sqrt((m.e-std::cos(y)*std::cos(y))*(m.e-std::cos(y)*std::cos(y))+2*m.f*m.f+(m.g-1)*(m.g-1));
    double de=std::cos(y)*std::cos(y)-m.e,df=-m.f,dg=1-m.g;
    r.minRho1=std::min(r.minRho1,de+df);r.minRho2=std::min(r.minRho2,de-df);
    if(dg-de<r.minRho3){r.minRho3=dg-de;r.minRho3Y=y;}
    M2 m0=f0.metric[k];M2 target{m0.e+targetFraction*(std::cos(y)*std::cos(y)-m0.e),
      (1-targetFraction)*m0.f,m0.g+targetFraction*(1-m0.g)};
    double te=std::sqrt((m.e-target.e)*(m.e-target.e)+2*(m.f-target.f)*(m.f-target.f)+(m.g-target.g)*(m.g-target.g));
    sum+=re;tsum+=te;if(re>r.roundMax){r.roundMax=re;r.roundMaxX=i*q.dx;r.roundMaxY=y;
      r.roundMaxE=m.e-std::cos(y)*std::cos(y);r.roundMaxF=m.f;r.roundMaxG=m.g-1;}
    r.targetMax=std::max(r.targetMax,te);count++;
  }
  for(auto p:q.p)r.bound=std::max(r.bound,norm(p));
  r.roundMean=sum/count;r.targetMean=tsum/count;return r;
}

void writeVTK(std::string const& name,Grid const&q,double eta){
  int capRows=std::max(12,q.ny/16),rings=2*(capRows-1)+q.ny;
  std::vector<V3> pts;pts.reserve(static_cast<size_t>(rings)*q.nx+2);
  pts.push_back({0,0,-1+eta});
  for(int k=1;k<capRows;k++){
    double y=-pi/2+(pi/2-q.ymax)*k/capRows;
    for(int i=0;i<q.nx;i++){double x=i*q.dx;pts.push_back({std::cos(y)*std::cos(x),std::cos(y)*std::sin(x),std::sin(y)+eta});}
  }
  pts.insert(pts.end(),q.p.begin(),q.p.end());
  for(int k=1;k<capRows;k++){
    double y=q.ymax+(pi/2-q.ymax)*k/capRows;
    for(int i=0;i<q.nx;i++){double x=i*q.dx;pts.push_back({std::cos(y)*std::cos(x),std::cos(y)*std::sin(x),std::sin(y)-eta});}
  }
  int north=static_cast<int>(pts.size());pts.push_back({0,0,1-eta});
  std::ofstream o(name);if(!o)throw std::runtime_error("cannot open VTK output");
  o<<"# vtk DataFile Version 3.0\nHevea reduced sphere\nASCII\nDATASET POLYDATA\nPOINTS "<<pts.size()<<" double\n"<<std::setprecision(12);
  for(auto p:pts)o<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
  size_t cells=static_cast<size_t>(rings-1)*q.nx+2*q.nx;
  o<<"POLYGONS "<<cells<<' '<<(static_cast<size_t>(rings-1)*q.nx*5+2*q.nx*4)<<"\n";
  for(int i=0;i<q.nx;i++)o<<"3 0 "<<1+(i+1)%q.nx<<' '<<1+i<<'\n';
  for(int j=0;j<rings-1;j++)for(int i=0;i<q.nx;i++){
    int k=1+j*q.nx+i,n=1+j*q.nx+(i+1)%q.nx;
    o<<"4 "<<k<<' '<<n<<' '<<n+q.nx<<' '<<k+q.nx<<'\n';
  }
  int last=1+(rings-1)*q.nx;
  for(int i=0;i<q.nx;i++)o<<"3 "<<last+i<<' '<<last+(i+1)%q.nx<<' '<<north<<'\n';
}

Grid sampled(Grid const&q,int maxNx=800,int maxNy=6000){
  int nx=std::min(q.nx,maxNx),ny=std::min(q.ny,maxNy);Grid out(nx,ny,q.ymax);
  #pragma omp parallel for
  for(int j=0;j<ny;j++)for(int i=0;i<nx;i++){
    int i0=i*q.nx/nx,i1=std::max(i0+1,(i+1)*q.nx/nx);
    int j0=j*q.ny/ny,j1=std::max(j0+1,(j+1)*q.ny/ny);V3 sum;int count=0;
    for(int jj=j0;jj<std::min(j1,q.ny);jj++)for(int ii=i0;ii<std::min(i1,q.nx);ii++){
      sum+=q.at(ii,jj);count++;
    }
    out.at(i,j)=sum/static_cast<double>(count);
  }
  return out;
}

void writeOutputs(std::string const& stem,Grid const&q,double eta){
  writeVTK(stem+".vtk",q,eta);
  int previewNx=800,previewNy=6000;
  if(const char* value=std::getenv("HEVEA_PREVIEW_NX"))previewNx=std::max(32,std::atoi(value));
  if(const char* value=std::getenv("HEVEA_PREVIEW_NY"))previewNy=std::max(64,std::atoi(value));
  Grid preview=sampled(q,previewNx,previewNy);
  writeVTK("sampled_"+stem+".vtk",preview,eta);
}

void writePreview(std::string const& stem,Grid const&q,double eta){
  int previewNx=800,previewNy=6000;
  if(const char* value=std::getenv("HEVEA_PREVIEW_NX"))previewNx=std::max(32,std::atoi(value));
  if(const char* value=std::getenv("HEVEA_PREVIEW_NY"))previewNy=std::max(64,std::atoi(value));
  writeVTK("sampled_"+stem+".vtk",sampled(q,previewNx,previewNy),eta);
}

void writeBinary(std::string const& name,Grid const& q){
  std::ofstream o(name,std::ios::binary);if(!o)throw std::runtime_error("cannot open binary output");
  const std::array<char,8> magic{'H','E','V','S','P','H','1','\0'};
  o.write(magic.data(),magic.size());o.write(reinterpret_cast<char const*>(&q.nx),sizeof(q.nx));
  o.write(reinterpret_cast<char const*>(&q.ny),sizeof(q.ny));
  o.write(reinterpret_cast<char const*>(&q.ymax),sizeof(q.ymax));
  o.write(reinterpret_cast<char const*>(q.p.data()),static_cast<std::streamsize>(q.p.size()*sizeof(V3)));
  if(!o)throw std::runtime_error("failed writing binary output");
}

Grid readBinary(std::string const& name,int xStride=1){
  std::ifstream input(name,std::ios::binary);if(!input)throw std::runtime_error("cannot read binary input");
  std::array<char,8> magic{};int nx=0,ny=0;double ymax=0;
  input.read(magic.data(),magic.size());input.read(reinterpret_cast<char*>(&nx),sizeof(nx));
  input.read(reinterpret_cast<char*>(&ny),sizeof(ny));
  input.read(reinterpret_cast<char*>(&ymax),sizeof(ymax));
  if(magic!=std::array<char,8>{'H','E','V','S','P','H','1','\0'}||nx<5||ny<5||
     xStride<1||nx%xStride!=0||nx/xStride<5)
    throw std::runtime_error("invalid binary input header");
  Grid q(nx/xStride,ny,ymax);
  if(xStride==1){
    input.read(reinterpret_cast<char*>(q.p.data()),static_cast<std::streamsize>(q.p.size()*sizeof(V3)));
  }else{
    std::vector<V3> row(nx);
    for(int j=0;j<ny;j++){
      input.read(reinterpret_cast<char*>(row.data()),static_cast<std::streamsize>(row.size()*sizeof(V3)));
      if(!input)break;
      for(int i=0;i<q.nx;i++)q.at(i,j)=row[i*xStride];
    }
  }
  if(!input)throw std::runtime_error("truncated binary input");
  return q;
}

int main(int argc,char**argv){
  try{
    if(argc!=9)throw std::runtime_error("expected 8 arguments");
    int nx=std::stoi(argv[1]),ny=std::stoi(argv[2]);
    std::array<int,3> ridges{std::stoi(argv[3]),std::stoi(argv[4]),std::stoi(argv[5])};
    double ballRadius=std::stod(argv[6]),eta=std::stod(argv[7]),targetFraction=std::stod(argv[8]);
    double ymax=hevea::sphere::yInfinity;
    Grid q=initialSurface(nx,ny,ymax,ballRadius,eta),initial=q;
    const bool diagnosticsOnly=std::getenv("HEVEA_DIAGNOSTICS_ONLY")!=nullptr;
    const bool previewOnly=std::getenv("HEVEA_PREVIEW_ONLY")!=nullptr;
    const int stopAfter=std::getenv("HEVEA_STOP_AFTER_STAGE")?
      std::atoi(std::getenv("HEVEA_STOP_AFTER_STAGE")):3;
    const bool binaryOutput=std::getenv("HEVEA_BINARY_OUTPUT")!=nullptr;
    int startDirection=0;
    if(const char* directory=std::getenv("HEVEA_DIAGNOSTIC_RESUME_STAGE2_DIRECTORY")){
      if(!diagnosticsOnly||previewOnly||binaryOutput)
        throw std::runtime_error("Stage-2 resume is restricted to diagnostics-only runs");
      const std::string root=directory;
      int xStride=1;
      if(const char* value=std::getenv("HEVEA_DIAGNOSTIC_RESUME_X_STRIDE"))xStride=std::atoi(value);
      initial=readBinary(root+"/reduced_sphere_stage=0.bin",xStride);
      q=readBinary(root+"/reduced_sphere_stage=2_dir=1_ridges=142.bin",xStride);
      if(initial.nx!=nx||initial.ny!=ny||q.nx!=nx||q.ny!=ny||
         std::abs(initial.ymax-ymax)>1e-12||std::abs(q.ymax-ymax)>1e-12)
        throw std::runtime_error("resumed grids do not match requested dimensions");
      startDirection=2;
      std::cout<<"RESUME Stage=2 Directory="<<root<<"\n";
    }
    Report r;
    if(startDirection==0){
      if(previewOnly)writePreview("reduced_sphere_stage=0",q,eta);
      else if(!diagnosticsOnly)writeOutputs("reduced_sphere_stage=0",q,eta);
      if(binaryOutput)writeBinary("reduced_sphere_stage=0.bin",q);
      r=report(q,initial,targetFraction);
      std::cout<<"METRIC Stage=0 RoundMax="<<r.roundMax<<" RoundMean="<<r.roundMean<<" BoundingRadius="<<r.bound
               <<" RoundMaxX="<<r.roundMaxX<<" RoundMaxY="<<r.roundMaxY
               <<" RoundMaxComponents="<<r.roundMaxE<<','<<r.roundMaxF<<','<<r.roundMaxG<<"\n";
    }
    std::array<double,4> ribbonFraction{0.41148336662952684,0.83293586428607069,
                                        1.38/ymax,(997.0/(2.0*334.92))/ymax};
    if(const char* value=std::getenv("HEVEA_RIBBON_FRACTIONS")){
      std::stringstream input(value);char comma=0;
      for(size_t k=0;k<ribbonFraction.size();k++){
        if(!(input>>ribbonFraction[k])||(k+1<ribbonFraction.size()&&
           (!(input>>comma)||comma!=',')))
          throw std::runtime_error("HEVEA_RIBBON_FRACTIONS requires four comma-separated numbers");
      }
      if(!(0<ribbonFraction[0]&&ribbonFraction[0]<ribbonFraction[1]&&
           ribbonFraction[1]<ribbonFraction[2]&&ribbonFraction[2]<ribbonFraction[3]&&
           ribbonFraction[3]<1))
        throw std::runtime_error("HEVEA_RIBBON_FRACTIONS must be strictly increasing in (0,1)");
    }
    std::array<double,3> previous{ribbonFraction[0]*ymax,ribbonFraction[1]*ymax,
                                  ribbonFraction[2]*ymax};
    std::array<double,3> outer{ribbonFraction[1]*ymax,ribbonFraction[2]*ymax,
                               ribbonFraction[3]*ymax};
    std::array<double,3> correctionScale{1.0,1.0,1.0};
    if(const char* value=std::getenv("HEVEA_CORRECTION_SCALES")){
      std::stringstream input(value);char comma=0;
      for(size_t k=0;k<correctionScale.size();k++){
        if(!(input>>correctionScale[k])||(k+1<correctionScale.size()&&
           (!(input>>comma)||comma!=',')))
          throw std::runtime_error("HEVEA_CORRECTION_SCALES requires three comma-separated numbers");
      }
      if(!std::all_of(correctionScale.begin(),correctionScale.end(),
          [](double scale){return std::isfinite(scale)&&scale>=0;}))
        throw std::runtime_error("HEVEA_CORRECTION_SCALES must be finite and non-negative");
    }
    std::array<double,3> phaseOffset{0.0,0.0,0.0};
    if(const char* value=std::getenv("HEVEA_PHASE_OFFSETS")){
      std::stringstream input(value);char comma=0;
      for(size_t k=0;k<phaseOffset.size();k++){
        if(!(input>>phaseOffset[k])||(k+1<phaseOffset.size()&&
           (!(input>>comma)||comma!=',')))
          throw std::runtime_error("HEVEA_PHASE_OFFSETS requires three comma-separated numbers");
      }
      if(!std::all_of(phaseOffset.begin(),phaseOffset.end(),
          [](double offset){return std::isfinite(offset);}))
        throw std::runtime_error("HEVEA_PHASE_OFFSETS must be finite");
    }
    double stageFraction=targetFraction;
    if(const char* value=std::getenv("HEVEA_STAGE_FRACTION"))stageFraction=std::stod(value);
    if(!(std::isfinite(stageFraction)&&stageFraction>0&&stageFraction<1))
      throw std::runtime_error("stage target fraction must be finite and in (0,1)");
    std::array<double,3> targetSchedule{stageFraction,stageFraction,stageFraction};
    if(const char* value=std::getenv("HEVEA_STAGE_FRACTIONS")){
      std::stringstream input(value);char comma=0;
      for(size_t k=0;k<targetSchedule.size();k++)if(!(input>>targetSchedule[k])||
          (k+1<targetSchedule.size()&&(!(input>>comma)||comma!=',')))
        throw std::runtime_error("HEVEA_STAGE_FRACTIONS requires three comma-separated numbers");
      if(!std::all_of(targetSchedule.begin(),targetSchedule.end(),[](double fraction){
           return std::isfinite(fraction)&&fraction>0.0&&fraction<1.0;}))
        throw std::runtime_error("HEVEA_STAGE_FRACTIONS must lie strictly between zero and one");
    }
    std::cout<<"CONFIG RibbonFractions="<<ribbonFraction[0]<<','<<ribbonFraction[1]<<','
             <<ribbonFraction[2]<<','<<ribbonFraction[3]
             <<" CorrectionScales="<<correctionScale[0]<<','<<correctionScale[1]<<','
             <<correctionScale[2]<<" PhaseOffsets="<<phaseOffset[0]<<','<<phaseOffset[1]<<','
             <<phaseOffset[2]<<" TargetFractions="<<targetSchedule[0]<<','<<targetSchedule[1]<<','
             <<targetSchedule[2]
             <<" LambdaPowers="<<(std::getenv("HEVEA_LAMBDA_POWERS")?std::getenv("HEVEA_LAMBDA_POWERS"):
                (std::getenv("HEVEA_LAMBDA_POWER")?std::getenv("HEVEA_LAMBDA_POWER"):
                 "0.60979123180328654,0,0"))
             <<" ChiPowers="<<(std::getenv("HEVEA_CHI_POWERS")?std::getenv("HEVEA_CHI_POWERS"):
                (std::getenv("HEVEA_CHI_POWER")?std::getenv("HEVEA_CHI_POWER"):
                 "0.60979123180328654,0,0"))
             <<" TransitionSplits="<<(std::getenv("HEVEA_TRANSITION_SPLITS")?
                std::getenv("HEVEA_TRANSITION_SPLITS"):
                (std::getenv("HEVEA_TRANSITION_SPLIT")?std::getenv("HEVEA_TRANSITION_SPLIT"):
                 "0.80193065649632189,0.5,0.5"))
             <<" CenterIncrement="<<(std::getenv("HEVEA_CENTER_INCREMENT")?"true":"false")
             <<" EquatorAnchorAll="<<(std::getenv("HEVEA_EQUATOR_ANCHOR_ALL")?"true":"false")
             <<" Stage3Anchor="<<(std::getenv("HEVEA_SOUTH_STAGE3_ANCHOR")?"south":"equator")
             <<" ExtendedProfile="<<((std::getenv("HEVEA_PROFILE_COEFFICIENTS")==nullptr||
                std::getenv("HEVEA_EXTENDED_PROFILE_COEFFICIENTS")!=nullptr)?"true":"false")<<"\n";
    for(int d=startDirection;d<3;d++){
      StepStats s;q=corrugate(q,initial,d,ridges[d],previous[d],outer[d],
                              targetSchedule[d],correctionScale[d],phaseOffset[d],s);
      std::string file="reduced_sphere_stage="+std::to_string(d+1)+"_dir="+std::to_string(d)+"_ridges="+std::to_string(ridges[d])+".vtk";
      if(previewOnly)writePreview(file.substr(0,file.size()-4),q,eta);
      else if(!diagnosticsOnly)writeOutputs(file.substr(0,file.size()-4),q,eta);
      r=report(q,initial,targetSchedule[d]);
      if(binaryOutput)writeBinary(file.substr(0,file.size()-4)+".bin",q);
      std::cout<<"METRIC Stage="<<d+1<<" RoundMax="<<r.roundMax<<" RoundMean="<<r.roundMean
               <<" TargetMax="<<s.targetMax<<" TargetMean="<<s.targetMean
               <<" BoundingRadius="<<r.bound<<" MinPrimitive="<<s.minRho
               <<" MinWy="<<s.minWy<<" MaxSlope="<<s.maxSlope
               <<" MinFlowSpacing="<<s.minSpacing<<" RoundMaxX="<<r.roundMaxX
               <<" MinFlowJacobian="<<s.minJacobian
               <<" InverseFlowRoundtrip="<<s.inverseRoundtrip
               <<" VelocityRelativeError="<<s.velocityRelativeError
               <<" PhaseClosureError="<<s.phaseClosureError
               <<" SeamPositionGap="<<s.seamPositionGap
               <<" SeamDerivativeGap="<<s.seamDerivativeGap
               <<" RoundMaxY="<<r.roundMaxY<<" RoundMaxComponents="<<r.roundMaxE<<','
               <<r.roundMaxF<<','<<r.roundMaxG<<" RoundRho1="<<r.minRho1
               <<" RoundRho2="<<r.minRho2<<" RoundRho3="<<r.minRho3
               <<" RoundRho3Y="<<r.minRho3Y<<"\n";
      if(d+1>=stopAfter)return 0;
    }
    return 0;
  }catch(std::exception const&e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 2;}
}
