#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct V3 {
  double x, y, z;
  V3 operator-(V3 other) const { return {x-other.x,y-other.y,z-other.z}; }
  V3 operator+(V3 other) const { return {x+other.x,y+other.y,z+other.z}; }
  V3 operator*(double scale) const { return {x*scale,y*scale,z*scale}; }
  V3 operator/(double scale) const { return {x/scale,y/scale,z/scale}; }
};

V3 operator*(double scale,V3 value){return value*scale;}
double dot(V3 a,V3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}

struct Grid {
  int nx=0,ny=0;
  double ymax=0,dx=0,dy=0;
  std::vector<V3> points;
  V3 at(int i,int j)const{
    i=(i%nx+nx)%nx;
    return points[static_cast<std::size_t>(j)*nx+i];
  }
};

struct Metric {double e=0,f=0,g=0;};
struct Values {
  double roundMean=0,roundMax=0,targetMean=0,targetMax=0,bound=0;
  double roundMaxX=0,roundMaxY=0,roundMaxE=0,roundMaxF=0,roundMaxG=0;
  double targetMaxX=0,targetMaxY=0,targetMaxE=0,targetMaxF=0,targetMaxG=0;
  double transverseDerivativeMean=0,transverseDerivativeMax=0;
};

std::string readText(std::string const&file){
  std::ifstream input(file);
  if(!input)throw std::runtime_error("cannot read "+file);
  return {(std::istreambuf_iterator<char>(input)),{}};
}

std::string stringField(std::string const&json,std::string const&key){
  std::smatch match;
  std::regex pattern("\\\""+key+"\\\"\\s*:\\s*\\\"([^\\\"]+)\\\"");
  if(!std::regex_search(json,match,pattern))throw std::runtime_error("manifest field missing: "+key);
  return match[1];
}

double realField(std::string const&json,std::string const&key){
  std::smatch match;
  std::regex pattern("\\\""+key+"\\\"\\s*:\\s*([-+0-9.eE]+)");
  if(!std::regex_search(json,match,pattern))throw std::runtime_error("manifest field missing: "+key);
  return std::stod(match[1]);
}

template<std::size_t N>
std::array<double,N> realArrayField(std::string const&json,std::string const&key){
  std::smatch match;
  std::regex pattern("\\\""+key+"\\\"\\s*:\\s*\\[([^\\]]+)\\]");
  if(!std::regex_search(json,match,pattern))throw std::runtime_error("manifest field missing: "+key);
  std::string values=match[1];std::replace(values.begin(),values.end(),',',' ');
  std::istringstream input(values);std::array<double,N> result{};
  for(double& value:result)if(!(input>>value))throw std::runtime_error("bad manifest array: "+key);
  double extra=0;if(input>>extra)throw std::runtime_error("bad manifest array length: "+key);
  return result;
}

Grid load(std::string const&file){
  std::ifstream input(file,std::ios::binary);
  if(!input)throw std::runtime_error("cannot read "+file);
  char magic[8];Grid grid;
  input.read(magic,8);
  input.read(reinterpret_cast<char*>(&grid.nx),sizeof(grid.nx));
  input.read(reinterpret_cast<char*>(&grid.ny),sizeof(grid.ny));
  input.read(reinterpret_cast<char*>(&grid.ymax),sizeof(grid.ymax));
  if(std::string(magic,7)!="HEVSPH1"||grid.nx<5||grid.ny<5)
    throw std::runtime_error("bad binary grid: "+file);
  grid.dx=2*std::acos(-1.0)/grid.nx;
  grid.dy=2*grid.ymax/(grid.ny-1);
  grid.points.resize(static_cast<std::size_t>(grid.nx)*grid.ny);
  input.read(reinterpret_cast<char*>(grid.points.data()),
             static_cast<std::streamsize>(grid.points.size()*sizeof(V3)));
  if(!input)throw std::runtime_error("truncated binary grid: "+file);
  return grid;
}

std::array<V3,2> derivatives(Grid const&grid,int i,int j){
  V3 x=(grid.at(i-2,j)-8.0*grid.at(i-1,j)+8.0*grid.at(i+1,j)-grid.at(i+2,j))/(12*grid.dx);
  V3 y;
  if(j==1)y=(-3.0*grid.at(i,0)-10.0*grid.at(i,1)+18.0*grid.at(i,2)-
              6.0*grid.at(i,3)+grid.at(i,4))/(12*grid.dy);
  else if(j==grid.ny-2)y=(3.0*grid.at(i,grid.ny-1)+10.0*grid.at(i,grid.ny-2)-
                           18.0*grid.at(i,grid.ny-3)+6.0*grid.at(i,grid.ny-4)-
                           grid.at(i,grid.ny-5))/(12*grid.dy);
  else y=(grid.at(i,j-2)-8.0*grid.at(i,j-1)+8.0*grid.at(i,j+1)-grid.at(i,j+2))/(12*grid.dy);
  return {x,y};
}

Metric metric(Grid const&grid,int i,int j){
  auto d=derivatives(grid,i,j);
  return {dot(d[0],d[0]),dot(d[0],d[1]),dot(d[1],d[1])};
}

double smooth(double value){
  value=std::clamp(value,0.0,1.0);
  return value*value*value*(10+value*(-15+6*value));
}

double poweredStep(double value,double power){
  value=std::clamp(value,0.0,1.0);
  if(power==0.0)return smooth(value);
  if(!(std::isfinite(power)&&power>1.0/3.0))
    throw std::runtime_error("manifest lambda power must be zero or greater than one third");
  if(value==0.0||value==1.0)return value;
  const double plateau=std::clamp(smooth(value),0.0,1.0);
  const double left=std::pow(plateau,power),right=std::pow(1.0-plateau,power);
  return left/(left+right);
}

double blend(double y,double previous,double outer,double power,double split){
  const double transition=previous+split*(outer-previous);
  return y<=previous?1.0:(y>=transition?0.0:
    1.0-poweredStep((y-previous)/(transition-previous),power));
}

Values measure(Grid const&grid,Grid const&initial,Grid const&previousGrid,int stage,double fraction,
               std::array<double,3> const&previous,std::array<double,3> const&outer,
               std::array<double,3> const&lambdaPowers,std::array<double,3> const&transitionSplits){
  double roundSum=0,targetSum=0,transverseSum=0;long long count=0;
  constexpr double inverseSqrtTwo=.7071067811865475244;
  const std::array<std::array<double,2>,3> transverse{{
    {-inverseSqrtTwo,inverseSqrtTwo},{-inverseSqrtTwo,-inverseSqrtTwo},{-1,0}}};
  Values values;
  for(int j=1;j<grid.ny-1;j++){
    double y=-grid.ymax+j*grid.dy,lambda=blend(std::abs(y),previous[stage-1],outer[stage-1],
      lambdaPowers[stage-1],transitionSplits[stage-1]);
    double roundE=std::cos(y)*std::cos(y);
    for(int i=0;i<grid.nx;i++){
      Metric actual=metric(grid,i,j),base=metric(initial,i,j);
      auto currentDerivative=derivatives(grid,i,j),previousDerivative=derivatives(previousGrid,i,j);
      V3 currentV=currentDerivative[0]*transverse[stage-1][0]+currentDerivative[1]*transverse[stage-1][1];
      V3 previousV=previousDerivative[0]*transverse[stage-1][0]+previousDerivative[1]*transverse[stage-1][1];
      double transverseError=std::sqrt(dot(currentV-previousV,currentV-previousV));
      double roundError=std::sqrt((actual.e-roundE)*(actual.e-roundE)+
        2*actual.f*actual.f+(actual.g-1)*(actual.g-1));
      Metric target{base.e+lambda*fraction*(roundE-base.e),
        base.f+lambda*fraction*(0-base.f),base.g+lambda*fraction*(1-base.g)};
      double targetError=std::sqrt((actual.e-target.e)*(actual.e-target.e)+
        2*(actual.f-target.f)*(actual.f-target.f)+(actual.g-target.g)*(actual.g-target.g));
      roundSum+=roundError;targetSum+=targetError;transverseSum+=transverseError;
      values.transverseDerivativeMax=std::max(values.transverseDerivativeMax,transverseError);
      if(roundError>values.roundMax){
        values.roundMax=roundError;values.roundMaxX=i*grid.dx;values.roundMaxY=y;
        values.roundMaxE=actual.e-roundE;values.roundMaxF=actual.f;
        values.roundMaxG=actual.g-1;
      }
      if(targetError>values.targetMax){
        values.targetMax=targetError;values.targetMaxX=i*grid.dx;values.targetMaxY=y;
        values.targetMaxE=actual.e-target.e;values.targetMaxF=actual.f-target.f;
        values.targetMaxG=actual.g-target.g;
      }
      count++;
    }
  }
  for(V3 point:grid.points)values.bound=std::max(values.bound,std::sqrt(dot(point,point)));
  values.roundMean=roundSum/count;values.targetMean=targetSum/count;
  values.transverseDerivativeMean=transverseSum/count;
  return values;
}

std::array<Values,3> productionValues(std::string const&log){
  std::array<Values,3> values{};
  std::regex pattern("METRIC Stage=([123]) RoundMax=([-+0-9.eE]+) RoundMean=([-+0-9.eE]+) "
    "TargetMax=([-+0-9.eE]+) TargetMean=([-+0-9.eE]+) BoundingRadius=([-+0-9.eE]+)");
  int found=0;
  for(std::sregex_iterator it(log.begin(),log.end(),pattern),end;it!=end;++it){
    int stage=std::stoi((*it)[1]);
    values[stage-1].roundMean=std::stod((*it)[3]);
    values[stage-1].roundMax=std::stod((*it)[2]);
    values[stage-1].targetMean=std::stod((*it)[5]);
    values[stage-1].targetMax=std::stod((*it)[4]);
    values[stage-1].bound=std::stod((*it)[6]);
    found++;
  }
  if(found!=3)throw std::runtime_error("native log does not contain three stage reports");
  return values;
}

double largestDifference(Values const&a,Values const&b){
  return std::max({std::abs(a.roundMean-b.roundMean),std::abs(a.roundMax-b.roundMax),
    std::abs(a.targetMean-b.targetMean),std::abs(a.targetMax-b.targetMax),std::abs(a.bound-b.bound)});
}

} // namespace

int main(int argc,char**argv){
  try{
    if(argc!=3||std::string(argv[1])!="--manifest")return 2;
    std::string manifestText=readText(argv[2]);
    std::string directory=stringField(manifestText,"output_directory");
    std::string nativeLog=stringField(manifestText,"native_log");
    double fraction=realField(manifestText,"target_fraction");
    const auto ribbonFractions=realArrayField<4>(manifestText,"ribbon_fractions");
    const auto lambdaPowers=realArrayField<3>(manifestText,"lambda_powers");
    const auto transitionSplits=realArrayField<3>(manifestText,"transition_splits");
    Grid initial=load(directory+"/reduced_sphere_stage=0.bin");
    const std::array<double,3> previous{ribbonFractions[0]*initial.ymax,
      ribbonFractions[1]*initial.ymax,ribbonFractions[2]*initial.ymax};
    const std::array<double,3> outer{ribbonFractions[1]*initial.ymax,
      ribbonFractions[2]*initial.ymax,ribbonFractions[3]*initial.ymax};
    Grid previousGrid;
    std::array<Values,3> production=productionValues(readText(nativeLog));
    constexpr std::array<double,3> paperMean{.83,.73,.66},paperMax{1.03,.95,.94};
    constexpr std::array<double,3> paperTargetMean{.14,.07,.03},paperTargetMax{.24,.16,.18};
    constexpr std::array<double,3> meanTolerance{.06,.06,.06},maxTolerance{.12,.12,.12};
    constexpr std::array<double,3> targetMeanTolerance{.04,.03,.02},targetMaxTolerance{.08,.06,.07};
    bool ok=true;
    std::cout<<std::setprecision(10);
    for(int stage=1;stage<=3;stage++){
      std::string suffix=stage==1?"_dir=0_ridges=21":stage==2?"_dir=1_ridges=142":"_dir=2_ridges=997";
      Grid grid=load(directory+"/reduced_sphere_stage="+std::to_string(stage)+suffix+".bin");
      Values measured=measure(grid,initial,stage==1?initial:previousGrid,stage,fraction,
        previous,outer,lambdaPowers,transitionSplits);
      double difference=largestDifference(measured,production[stage-1]);
      bool stageOk=std::abs(measured.roundMean-paperMean[stage-1])<=meanTolerance[stage-1]&&
        std::abs(measured.roundMax-paperMax[stage-1])<=maxTolerance[stage-1]&&
        std::abs(measured.targetMean-paperTargetMean[stage-1])<=targetMeanTolerance[stage-1]&&
        std::abs(measured.targetMax-paperTargetMax[stage-1])<=targetMaxTolerance[stage-1]&&
        measured.bound<=.5201&&difference<1e-4;
      std::cout<<"stage="<<stage
        <<" round_mean="<<measured.roundMean<<" paper_round_mean="<<paperMean[stage-1]
        <<" round_max="<<measured.roundMax<<" paper_round_max="<<paperMax[stage-1]
        <<" target_mean="<<measured.targetMean<<" paper_target_mean="<<paperTargetMean[stage-1]
        <<" target_max="<<measured.targetMax<<" paper_target_max="<<paperTargetMax[stage-1]
        <<" bound="<<measured.bound<<" production_delta="<<difference
        <<" round_max_at="<<measured.roundMaxX<<','<<measured.roundMaxY
        <<" round_max_components="<<measured.roundMaxE<<','<<measured.roundMaxF<<','<<measured.roundMaxG
        <<" target_max_at="<<measured.targetMaxX<<','<<measured.targetMaxY
        <<" target_max_components="<<measured.targetMaxE<<','<<measured.targetMaxF<<','<<measured.targetMaxG
        <<" transverse_derivative_mean="<<measured.transverseDerivativeMean
        <<" transverse_derivative_max="<<measured.transverseDerivativeMax
        <<" status="<<(stageOk?"pass":"fail")<<'\n';
      if(stage==3){
        int i=static_cast<int>(std::llround(measured.targetMaxX/grid.dx))%grid.nx;
        int j=std::clamp(static_cast<int>(std::llround(
          (measured.targetMaxY+grid.ymax)/grid.dy)),1,grid.ny-2);
        double y=-grid.ymax+j*grid.dy;
        Metric actual=metric(grid,i,j),before=metric(previousGrid,i,j),base=metric(initial,i,j);
        double lambda=blend(std::abs(y),previous[2],outer[2],lambdaPowers[2],transitionSplits[2]);
        Metric target{base.e+lambda*fraction*(std::cos(y)*std::cos(y)-base.e),
          base.f+lambda*fraction*(0-base.f),base.g+lambda*fraction*(1-base.g)};
        std::cout<<"stage3_target_point_before_components="
          <<before.e-target.e<<','<<before.f-target.f<<','<<before.g-target.g
          <<" stage3_target_point_after_components="
          <<actual.e-target.e<<','<<actual.f-target.f<<','<<actual.g-target.g
          <<" stage3_target_point_metric_delta="
          <<actual.e-before.e<<','<<actual.f-before.f<<','<<actual.g-before.g
          <<" stage3_target_point_ideal_rho3="
          <<(target.g-before.g)-(target.e-before.e)<<'\n';
      }
      if(stage==2){
        Values beforeThird=measure(grid,initial,grid,3,fraction,
          previous,outer,lambdaPowers,transitionSplits);
        std::cout<<"pre_stage3_target_mean="<<beforeThird.targetMean
          <<" pre_stage3_target_max="<<beforeThird.targetMax
          <<" pre_stage3_target_max_at="<<beforeThird.targetMaxX<<','<<beforeThird.targetMaxY
          <<" pre_stage3_target_max_components="<<beforeThird.targetMaxE<<','
          <<beforeThird.targetMaxF<<','<<beforeThird.targetMaxG<<'\n';
      }
      ok&=stageOk;
      previousGrid=std::move(grid);
    }
    return ok?0:1;
  }catch(std::exception const&error){
    std::cerr<<"verify error: "<<error.what()<<'\n';return 2;
  }
}
