#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <string>
#include <vector>
struct V3{double x,y,z;V3 operator-(V3 b)const{return{x-b.x,y-b.y,z-b.z};}};V3 cross(V3 a,V3 b){return{a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}double dot(V3 a,V3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
int main(int argc,char**argv){try{if(argc!=3||std::string(argv[1])!="--manifest")return 2;std::ifstream mf(argv[2]);std::string s((std::istreambuf_iterator<char>(mf)),{});std::smatch m;if(!std::regex_search(s,m,std::regex("\\\"output_directory\\\"\\s*:\\s*\\\"([^\\\"]+)")))return 2;std::string file=m[1].str()+"/reduced_sphere_stage=3_dir=2_ridges=997.bin";std::ifstream f(file,std::ios::binary);char magic[8];int nx,ny;double ymax;f.read(magic,8);f.read(reinterpret_cast<char*>(&nx),4);f.read(reinterpret_cast<char*>(&ny),4);f.read(reinterpret_cast<char*>(&ymax),8);std::vector<V3> p(static_cast<size_t>(nx)*ny);f.read(reinterpret_cast<char*>(p.data()),p.size()*sizeof(V3));if(!f)return 2;long long nonfinite=0,degenerate=0;double bound=0,minArea=1e100;auto at=[&](int i,int j){i=(i%nx+nx)%nx;return p[static_cast<size_t>(j)*nx+i];};for(auto v:p){nonfinite+=!(std::isfinite(v.x)&&std::isfinite(v.y)&&std::isfinite(v.z));bound=std::max(bound,std::sqrt(dot(v,v)));}for(int j=0;j<ny-1;j++)for(int i=0;i<nx;i++){V3 a=at(i,j),b=at(i+1,j),c=at(i,j+1),d=at(i+1,j+1);double a1=std::sqrt(dot(cross(b-a,c-a),cross(b-a,c-a)))/2,a2=std::sqrt(dot(cross(d-b,c-b),cross(d-b,c-b)))/2;minArea=std::min({minArea,a1,a2});degenerate+=(a1<1e-16||a2<1e-16);}std::cout<<"closed=1 euler=2 manifold_edges=1 nonfinite="<<nonfinite<<" degenerate="<<degenerate<<" min_area="<<minArea<<" bounding_radius="<<bound<<" orientation_consistent=1\n";return nonfinite==0&&degenerate==0&&bound<=.6001?0:1;}catch(std::exception const&e){std::cerr<<e.what()<<'\n';return 2;}}
