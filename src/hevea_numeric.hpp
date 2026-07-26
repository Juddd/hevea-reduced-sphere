#pragma once

/* Numerical convex-integration utilities adapted from the ideas and code
 * structure of https://github.com/HeveaProject/Hevea (GPLv3). */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace hevea {

inline double j0Inverse(double value) {
  if (!std::isfinite(value) || value < 0.0 || value > 1.0)
    throw std::domain_error("J0 inverse argument must lie in [0,1]");
  double lo = 0.0, hi = 2.4048255576957727686;
  for (int i = 0; i < 64; ++i) {
    double mid = (lo + hi) / 2.0;
    if (std::cyl_bessel_j(0, mid) > value) lo = mid;
    else hi = mid;
  }
  return (lo + hi) / 2.0;
}

struct OdeStats {
  std::size_t accepted = 0;
  std::size_t rejected = 0;
};

template<std::size_t N> using State = std::array<double, N>;

template<class T>
std::array<T,2> nonuniformHermiteTangents(
    T const& previous,T const& left,T const& right,T const& next,
    double previousSpacing,double spacing,double nextSpacing) {
  if(!(previousSpacing>0.0&&spacing>0.0&&nextSpacing>0.0))
    throw std::invalid_argument("non-positive nonuniform cell spacing");
  return {0.5*((right-left)+(spacing/previousSpacing)*(left-previous)),
          0.5*((right-left)+(spacing/nextSpacing)*(next-right))};
}

template<class T>
T cubicHermite(double t,T const& left,T const& right,
               T const& leftTangent,T const& rightTangent) {
  const double h01=t*t*(3.0-2.0*t);
  const double h10=t*(t*(t-2.0)+1.0);
  const double h11=t*t*(t-1.0);
  return left+h01*(right-left)+h10*leftTangent+h11*rightTangent;
}

template<class T>
T cubicHermiteDerivative(double t,T const& left,T const& right,
                         T const& leftTangent,T const& rightTangent) {
  const double dh01=6.0*t*(1.0-t);
  const double dh10=(t-1.0)*(3.0*t-1.0);
  const double dh11=t*(3.0*t-2.0);
  return dh01*(right-left)+dh10*leftTangent+dh11*rightTangent;
}

template<class T> struct PeriodicMonotoneSample {
  T value{};
  double cellSpacing = 0.0;
  double linearRoundTripError = 0.0;
};

// This is the x-periodic, y-bounded inverse-flow resampling used by
// CYL_embedding::cyl_to_torus in the original Hevea code.  The samples are
// stored in their (strictly monotone) forward-flow order; values are then
// interpolated as a function of the non-uniform image coordinate.  The value
// type is generic: callers may resample a complete map or a numerically more
// stable map increment.
template<class T,class PositionGetter,class ValueGetter>
PeriodicMonotoneSample<T> periodicMonotoneCubicSample(
    int count,double period,double target,
    PositionGetter&& position,ValueGetter&& value) {
  if(count<4||!(period>0.0)||!std::isfinite(target))
    throw std::invalid_argument("invalid periodic monotone interpolation input");
  const double first=position(0);
  if(!std::isfinite(first))
    throw std::invalid_argument("non-finite periodic monotone node");
  target=first+std::fmod(target-first,period);
  if(target<first)target+=period;

  auto node=[&](int index){
    int cycle=index/count,base=index%count;
    if(base<0){base+=count;--cycle;}
    const double coordinate=position(base)+cycle*period;
    if(!std::isfinite(coordinate))
      throw std::invalid_argument("non-finite periodic monotone node");
    return std::pair<double,T>{coordinate,value(base)};
  };
  int low=0,high=count;
  while(low<high){
    const int middle=low+(high-low)/2;
    if(node(middle).first<=target)low=middle+1;else high=middle;
  }
  const int right=low,left=right-1;
  const auto previous=node(left-1),p0=node(left),p1=node(right),next=node(right+1);
  const double spacing=p1.first-p0.first;
  if(!(previous.first<p0.first&&p0.first<p1.first&&p1.first<next.first))
    throw std::invalid_argument("non-positive periodic monotone cell spacing");
  const double t=(target-p0.first)/spacing;
  const auto tangents=nonuniformHermiteTangents(
      previous.second,p0.second,p1.second,next.second,
      p0.first-previous.first,spacing,next.first-p1.first);
  return {cubicHermite(t,p0.second,p1.second,tangents[0],tangents[1]),
          spacing,std::abs((p0.first+t*spacing)-target)};
}

template<class T>
T shiftedFourthOrderDerivative(std::array<T,6> const& values,int shift,int node) {
  if(node!=0&&node!=1)throw std::invalid_argument("Hermite node must be 0 or 1");
  T value{};
  if(shift==0)
    value=(values[node]-values[node+4])+8.0*(values[node+3]-values[node+1]);
  else if(shift==2&&node==0)
    value=-25.0*values[0]+48.0*values[1]-36.0*values[2]+16.0*values[3]-3.0*values[4];
  else if((shift==2&&node==1)||(shift==1&&node==0))
    value=-3.0*values[0]-10.0*values[1]+18.0*values[2]-6.0*values[3]+values[4];
  else if((shift==1&&node==1)||(shift==-1&&node==0))
    value=(values[1-node]-values[5-node])+8.0*(values[4-node]-values[2-node]);
  else if((shift==-1&&node==1)||(shift==-2&&node==0))
    value=-1.0*values[1]+6.0*values[2]-18.0*values[3]+10.0*values[4]+3.0*values[5];
  else if(shift==-2&&node==1)
    value=3.0*values[1]-16.0*values[2]+36.0*values[3]-48.0*values[4]+25.0*values[5];
  else
    throw std::invalid_argument("invalid shifted derivative stencil");
  return value/12.0;
}

template<class T> struct CylinderCubicJet {
  T value{};
  T derivativeX{};
  T derivativeY{};
};

template<class T,class Getter>
CylinderCubicJet<T> cylinderCubicJet(
    Getter&& get,int nx,int ny,double dx,double dy,double ymin,double x,double y) {
  if(nx<6||ny<6||!(dx>0.0&&dy>0.0))
    throw std::invalid_argument("invalid bicubic cylinder grid");
  const double period=nx*dx;
  x=std::fmod(x,period);
  if(x<0.0)x+=period;
  const double gx=x/dx;
  const double gy=std::clamp((y-ymin)/dy,0.0,static_cast<double>(ny-1));
  const int i=static_cast<int>(std::floor(gx));
  const int j=std::clamp(static_cast<int>(std::floor(gy)),0,ny-2);
  const double u=gx-i;
  const double v=std::clamp(gy-j,0.0,1.0);
  const int jAux=std::clamp(j,2,ny-4);
  const int shift=jAux-j;
  const int index=2-shift;
  auto periodic=[&](int ix,int iy){ix%=nx;if(ix<0)ix+=nx;return get(ix,iy);};
  std::array<T,6> rows{},rowsDx{};
  for(int k=0;k<6;k++){
    const int iy=jAux+k-2;
    const T left=periodic(i,iy),right=periodic(i+1,iy);
    const T leftTangent=((periodic(i-2,iy)-periodic(i+2,iy))+
                         8.0*(periodic(i+1,iy)-periodic(i-1,iy)))/12.0;
    const T rightTangent=((periodic(i-1,iy)-periodic(i+3,iy))+
                          8.0*(periodic(i+2,iy)-periodic(i,iy)))/12.0;
    rows[k]=cubicHermite(u,left,right,leftTangent,rightTangent);
    rowsDx[k]=cubicHermiteDerivative(u,left,right,leftTangent,rightTangent);
  }
  const T yTangentLeft=shiftedFourthOrderDerivative(rows,shift,0);
  const T yTangentRight=shiftedFourthOrderDerivative(rows,shift,1);
  const T dxTangentLeft=shiftedFourthOrderDerivative(rowsDx,shift,0);
  const T dxTangentRight=shiftedFourthOrderDerivative(rowsDx,shift,1);
  return {
    cubicHermite(v,rows[index],rows[index+1],yTangentLeft,yTangentRight),
    cubicHermite(v,rowsDx[index],rowsDx[index+1],dxTangentLeft,dxTangentRight)/dx,
    cubicHermiteDerivative(v,rows[index],rows[index+1],yTangentLeft,yTangentRight)/dy
  };
}

template<std::size_t N, class Function>
State<N> integrateDp54(Function&& f, double t0, State<N> state, double t1,
                       double relativeTolerance = 1e-11,
                       double absoluteTolerance = 1e-12,
                       OdeStats* statistics = nullptr) {
  if (!(t1 >= t0) || relativeTolerance <= 0 || absoluteTolerance <= 0)
    throw std::invalid_argument("invalid ODE interval or tolerance");
  OdeStats local;
  double t = t0, h = std::max(t1 - t0, 1e-8);
  auto add = [](State<N> const& a, double scale, State<N> const& b) {
    State<N> out{}; for (std::size_t i=0;i<N;i++) out[i]=a[i]+scale*b[i]; return out;
  };
  while (t < t1) {
    h = std::min(h, t1 - t);
    auto k1=f(t,state);
    auto s2=add(state,h*(1.0/5.0),k1); auto k2=f(t+h/5.0,s2);
    State<N> s3{},s4{},s5{},s6{},s7{};
    for(std::size_t i=0;i<N;i++)s3[i]=state[i]+h*(3.0*k1[i]/40.0+9.0*k2[i]/40.0);
    auto k3=f(t+3.0*h/10.0,s3);
    for(std::size_t i=0;i<N;i++)s4[i]=state[i]+h*(44.0*k1[i]/45.0-56.0*k2[i]/15.0+32.0*k3[i]/9.0);
    auto k4=f(t+4.0*h/5.0,s4);
    for(std::size_t i=0;i<N;i++)s5[i]=state[i]+h*(19372.0*k1[i]/6561.0-25360.0*k2[i]/2187.0+64448.0*k3[i]/6561.0-212.0*k4[i]/729.0);
    auto k5=f(t+8.0*h/9.0,s5);
    for(std::size_t i=0;i<N;i++)s6[i]=state[i]+h*(9017.0*k1[i]/3168.0-355.0*k2[i]/33.0+46732.0*k3[i]/5247.0+49.0*k4[i]/176.0-5103.0*k5[i]/18656.0);
    auto k6=f(t+h,s6);
    for(std::size_t i=0;i<N;i++)s7[i]=state[i]+h*(35.0*k1[i]/384.0+500.0*k3[i]/1113.0+125.0*k4[i]/192.0-2187.0*k5[i]/6784.0+11.0*k6[i]/84.0);
    auto k7=f(t+h,s7);
    State<N> fourth{}; double error=0;
    for(std::size_t i=0;i<N;i++){
      fourth[i]=state[i]+h*(5179.0*k1[i]/57600.0+7571.0*k3[i]/16695.0+393.0*k4[i]/640.0-92097.0*k5[i]/339200.0+187.0*k6[i]/2100.0+k7[i]/40.0);
      double scale=absoluteTolerance+relativeTolerance*std::max(std::abs(state[i]),std::abs(s7[i]));
      error=std::max(error,std::abs(s7[i]-fourth[i])/scale);
    }
    if(error<=1.0){state=s7;t+=h;local.accepted++;}
    else local.rejected++;
    double factor=error==0?5.0:std::clamp(0.9*std::pow(error,-0.2),0.2,5.0);
    h*=factor;
    if(h < std::numeric_limits<double>::epsilon()*std::max(1.0,std::abs(t)))
      throw std::runtime_error("ODE step underflow");
  }
  if(statistics)*statistics=local;
  return state;
}

template<class T> class CylinderGrid {
 public:
  CylinderGrid(std::size_t nx,std::size_t ny,double ymin,double ymax)
      : nx_(nx),ny_(ny),ymin_(ymin),ymax_(ymax),data_(nx*ny) {
    if(nx<2||ny<2||!(ymax>ymin))throw std::invalid_argument("invalid cylinder grid");
  }
  T& at(std::ptrdiff_t ix,std::size_t iy){return data_[index(ix,iy)];}
  T const& at(std::ptrdiff_t ix,std::size_t iy)const{return data_[index(ix,iy)];}
  std::size_t nx()const{return nx_;} std::size_t ny()const{return ny_;}
 private:
  std::size_t index(std::ptrdiff_t ix,std::size_t iy)const{
    if(iy>=ny_)throw std::out_of_range("cylinder y is bounded, not periodic");
    auto n=static_cast<std::ptrdiff_t>(nx_);ix%=n;if(ix<0)ix+=n;
    return iy*nx_+static_cast<std::size_t>(ix);
  }
  std::size_t nx_,ny_;double ymin_,ymax_;std::vector<T> data_;
};

struct InverseFlowReport { double roundTripMax=0; double minJacobian=1; };

inline InverseFlowReport testFlowInverse(std::size_t samples,double amplitude) {
  if(samples<3||std::abs(amplitude)>=0.5)throw std::invalid_argument("invalid flow test");
  std::vector<double> image(samples);
  OdeStats stats;
  auto field=[=](double,State<1> const& x){return State<1>{amplitude*std::sin(2.0*3.14159265358979323846*x[0])};};
  for(std::size_t i=0;i<samples;i++)image[i]=integrateDp54<1>(field,0,{double(i)/(samples-1)},1,1e-12,1e-13,&stats)[0];
  InverseFlowReport report;
  double dx=1.0/(samples-1);
  for(std::size_t i=0;i+1<samples;i++)report.minJacobian=std::min(report.minJacobian,(image[i+1]-image[i])/dx);
  for(std::size_t i=0;i<samples;i++){
    double target=image[i];auto it=std::lower_bound(image.begin(),image.end(),target);
    std::size_t hi=std::min<std::size_t>(samples-1,std::max<std::size_t>(1,it-image.begin()));
    std::size_t lo=hi-1;double u=(target-image[lo])/(image[hi]-image[lo]);
    double inverse=(lo+u)*dx;report.roundTripMax=std::max(report.roundTripMax,std::abs(inverse-i*dx));
  }
  return report;
}

} // namespace hevea
