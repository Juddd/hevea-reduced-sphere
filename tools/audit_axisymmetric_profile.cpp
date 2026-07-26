#include "reduced_sphere_profile.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct V3 { double x=0,y=0,z=0; };

double square(double value){return value*value;}

} // namespace

int main(int argc,char**argv){
  try{
    if(argc!=2)throw std::runtime_error("usage: audit_axisymmetric_profile GRID.bin");
    std::ifstream input(argv[1],std::ios::binary);
    if(!input)throw std::runtime_error("cannot open grid");
    std::array<char,8>magic{};std::int32_t nx=0,ny=0;double ymax=0;
    input.read(magic.data(),magic.size());input.read(reinterpret_cast<char*>(&nx),sizeof(nx));
    input.read(reinterpret_cast<char*>(&ny),sizeof(ny));
    input.read(reinterpret_cast<char*>(&ymax),sizeof(ymax));
    if(std::string(magic.data(),7)!="HEVSPH1"||nx<2||ny<3||!(ymax>0))
      throw std::runtime_error("invalid grid header");

    // Degree-nine fit of the axisymmetric mean of the authors' public f_1,3
    // WRL, independently recovered by audit_reference_wrl.
    constexpr std::array<double,6> reference{
      -0.358054775023,0.0509868430662,0.0000858502342762,
      -0.311467514785,0.0864614730285,-0.00861692529856};
    const double twoPi=2*std::acos(-1.0),dy=2*ymax/(ny-1);
    std::vector<double> cosine(nx),sine(nx);std::vector<V3>row(nx);
    for(int i=0;i<nx;i++){double angle=twoPi*i/nx;cosine[i]=std::cos(angle);sine[i]=std::sin(angle);}
    long double sumSquared=0;double maximum=0,maximumY=0,maximumTangential=0;
    for(int j=0;j<ny;j++){
      input.read(reinterpret_cast<char*>(row.data()),static_cast<std::streamsize>(row.size()*sizeof(V3)));
      if(!input)throw std::runtime_error("truncated grid data");
      long double radial=0,tangential=0,z=0;
      for(int i=0;i<nx;i++){
        radial+=row[i].x*cosine[i]+row[i].y*sine[i];
        tangential+=-row[i].x*sine[i]+row[i].y*cosine[i];z+=row[i].z;
      }
      radial/=nx;tangential/=nx;z/=nx;
      const double y=-ymax+j*dy;const auto target=hevea::sphere::evaluate(y,reference);
      const double error=std::sqrt(square(static_cast<double>(radial)-target.x)+
                                   square(static_cast<double>(z)-target.z));
      sumSquared+=error*error;maximumTangential=std::max(maximumTangential,std::abs(static_cast<double>(tangential)));
      if(error>maximum){maximum=error;maximumY=y;}
    }
    std::cout<<std::setprecision(12)<<"axisymmetric_profile=pass grid="<<nx<<'x'<<ny
      <<" reference=public-wrl-global-degree9 rms="<<std::sqrt(static_cast<double>(sumSquared/ny))
      <<" max="<<maximum<<" max_y="<<maximumY
      <<" max_tangential_mean="<<maximumTangential<<'\n';
    return 0;
  }catch(std::exception const&error){
    std::cerr<<"audit_axisymmetric_profile: "<<error.what()<<'\n';return 2;
  }
}
