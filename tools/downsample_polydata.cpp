#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

struct Point { double x=0,y=0,z=0; };

int main(int argc,char**argv){
  try{
    if(argc!=6)throw std::runtime_error("usage: downsample_polydata input.vtk output.vtk sourceNx targetNx targetRings");
    const int sourceNx=std::stoi(argv[3]),targetNx=std::stoi(argv[4]),requestedRings=std::stoi(argv[5]);
    if(sourceNx<3||targetNx<3||targetNx>sourceNx||requestedRings<2)
      throw std::runtime_error("invalid sampling dimensions");
    std::ifstream input(argv[1]);if(!input)throw std::runtime_error("cannot open input VTK");
    std::string token;
    while(input>>token&&token!="POINTS"){}
    std::size_t pointCount=0;std::string scalarType;
    if(token!="POINTS"||!(input>>pointCount>>scalarType)||pointCount<2||
       (pointCount-2)%static_cast<std::size_t>(sourceNx)!=0)
      throw std::runtime_error("input is not a compatible closed structured polydata mesh");
    const int sourceRings=static_cast<int>((pointCount-2)/sourceNx);
    const int targetRings=std::min(requestedRings,sourceRings);
    std::vector<int> ringIndices(targetRings),columnIndices(targetNx);
    for(int j=0;j<targetRings;j++)ringIndices[j]=static_cast<int>(std::llround(
      static_cast<double>(j)*(sourceRings-1)/(targetRings-1)));
    ringIndices.erase(std::unique(ringIndices.begin(),ringIndices.end()),ringIndices.end());
    for(int i=0;i<targetNx;i++)columnIndices[i]=i*sourceNx/targetNx;
    Point south,north;if(!(input>>south.x>>south.y>>south.z))throw std::runtime_error("truncated south pole");
    std::vector<Point> rings;rings.reserve(static_cast<std::size_t>(ringIndices.size())*targetNx);
    std::size_t wantedRing=0;
    for(int j=0;j<sourceRings;j++){
      std::size_t wantedColumn=0;
      for(int i=0;i<sourceNx;i++){
        Point p;if(!(input>>p.x>>p.y>>p.z))throw std::runtime_error("truncated ring points");
        if(wantedRing<ringIndices.size()&&j==ringIndices[wantedRing]&&
           wantedColumn<columnIndices.size()&&i==columnIndices[wantedColumn]){
          rings.push_back(p);++wantedColumn;
        }
      }
      if(wantedRing<ringIndices.size()&&j==ringIndices[wantedRing])++wantedRing;
    }
    if(!(input>>north.x>>north.y>>north.z)||
       rings.size()!=static_cast<std::size_t>(ringIndices.size())*targetNx)
      throw std::runtime_error("truncated north pole or incomplete sampled rings");
    std::ofstream output(argv[2]);if(!output)throw std::runtime_error("cannot open output VTK");
    const int ringsCount=static_cast<int>(ringIndices.size());
    const std::size_t outputPoints=2+rings.size();
    const std::size_t cells=static_cast<std::size_t>(ringsCount-1)*targetNx+2*targetNx;
    output<<"# vtk DataFile Version 3.0\nHevea reduced sphere preview\nASCII\nDATASET POLYDATA\nPOINTS "
      <<outputPoints<<" double\n"<<std::setprecision(12);
    output<<south.x<<' '<<south.y<<' '<<south.z<<'\n';
    for(auto const&p:rings)output<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
    output<<north.x<<' '<<north.y<<' '<<north.z<<'\n';
    output<<"POLYGONS "<<cells<<' '<<(static_cast<std::size_t>(ringsCount-1)*targetNx*5+
      2*static_cast<std::size_t>(targetNx)*4)<<"\n";
    for(int i=0;i<targetNx;i++)output<<"3 0 "<<1+(i+1)%targetNx<<' '<<1+i<<'\n';
    for(int j=0;j<ringsCount-1;j++)for(int i=0;i<targetNx;i++){
      const int k=1+j*targetNx+i,n=1+j*targetNx+(i+1)%targetNx;
      output<<"4 "<<k<<' '<<n<<' '<<n+targetNx<<' '<<k+targetNx<<'\n';
    }
    const int last=1+(ringsCount-1)*targetNx,northIndex=static_cast<int>(outputPoints)-1;
    for(int i=0;i<targetNx;i++)output<<"3 "<<last+i<<' '<<last+(i+1)%targetNx<<' '<<northIndex<<'\n';
    if(!output)throw std::runtime_error("failed writing output VTK");
    std::cout<<"source_points="<<pointCount<<" output_points="<<outputPoints
      <<" source_rings="<<sourceRings<<" output_rings="<<ringsCount<<" target_nx="<<targetNx<<'\n';
    return 0;
  }catch(std::exception const&e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 2;}
}
