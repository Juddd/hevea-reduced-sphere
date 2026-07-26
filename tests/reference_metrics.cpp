#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

namespace {
constexpr double pi = 3.14159265358979323846;
struct V { double x,y,z; };
V operator-(V a,V b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
V operator/(V a,double s){return {a.x/s,a.y/s,a.z/s};}
double dot(V a,V b){return a.x*b.x+a.y*b.y+a.z*b.z;}
struct M {double e,f,g;};
M metric(V fx,V fy){return {dot(fx,fx),dot(fx,fy),dot(fy,fy)};}
double maxdiff(M a,M b){return std::max({std::abs(a.e-b.e),std::abs(a.f-b.f),std::abs(a.g-b.g)});}
std::array<double,3> primitive(M m){return {m.e+m.f,m.e-m.f,m.g-m.e};}
M reconstruct(std::array<double,3> r){return {(r[0]+r[1])/2,(r[0]-r[1])/2,(r[0]+r[1])/2+r[2]};}
V sphere(double x,double y,double s=1){return {s*std::cos(y)*std::cos(x),s*std::cos(y)*std::sin(x),s*std::sin(y)};}
V torus(double x,double y,double R,double r){double q=R+r*std::cos(y);return {q*std::cos(x),q*std::sin(x),r*std::sin(y)};}
template<class F>M finiteMetric(F f,double x,double y){constexpr double h=1e-5;return metric((f(x+h,y)-f(x-h,y))/(2*h),(f(x,y+h)-f(x,y-h))/(2*h));}
bool near(double x,double tol){return std::abs(x)<=tol;}
}

int main(){
  double worst=0;
  for(int j=0;j<25;j++)for(int i=0;i<31;i++){
    double x=2*pi*i/31,y=-1.3+2.6*j/24;
    M round{std::cos(y)*std::cos(y),0,1};
    worst=std::max(worst,maxdiff(finiteMetric([](double a,double b){return sphere(a,b);},x,y),round));
    M half{.25*round.e,0,.25};
    worst=std::max(worst,maxdiff(finiteMetric([](double a,double b){return sphere(a,b,.5);},x,y),half));
    auto rho=primitive(round);worst=std::max(worst,maxdiff(reconstruct(rho),round));
  }
  double torusWorst=0,R=.5,r=.2;
  for(int j=0;j<25;j++)for(int i=0;i<31;i++){
    double x=2*pi*i/31,y=2*pi*j/25;
    M exact{std::pow(R+r*std::cos(y),2),0,r*r};
    torusWorst=std::max(torusWorst,maxdiff(finiteMetric([&](double a,double b){return torus(a,b,R,r);},x,y),exact));
  }
  std::cout<<"standard_scaled_primitive_max="<<worst<<"\nflat_torus_metric_max="<<torusWorst<<"\n";
  if(!near(worst,1e-6)||!near(torusWorst,1e-6))return 1;
  return 0;
}
