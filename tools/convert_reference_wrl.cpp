#include <algorithm>
#include <array>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct V3 { double x = 0.0, y = 0.0, z = 0.0; };

std::vector<V3> readPoints(const std::string& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open WRL: " + path);
  bool inPoints = false;
  std::vector<V3> points;
  points.reserve(4'100'000);
  std::string line;
  while (std::getline(input, line)) {
    if (!inPoints) {
      if (line.find("point [") != std::string::npos) inPoints = true;
      continue;
    }
    if (line.find(']') != std::string::npos) break;
    V3 point;
    if (std::sscanf(line.c_str(), "%lf %lf %lf", &point.x, &point.y,
                    &point.z) == 3)
      points.push_back(point);
  }
  if (!inPoints || points.size() < 1000)
    throw std::runtime_error("WRL point section is incomplete");
  return points;
}

void writeVtk(const std::string& path, const std::vector<V3>& source,
              int sourceNx, int outputNx, double scale) {
  if (sourceNx < 4 || outputNx < 4 || outputNx > sourceNx ||
      !(scale > 0.0) || source.size() < 2 ||
      (source.size() - 2) % static_cast<std::size_t>(sourceNx) != 0)
    throw std::runtime_error("invalid structured WRL dimensions or scale");
  const int rows = static_cast<int>((source.size() - 2) /
                                    static_cast<std::size_t>(sourceNx));
  const std::size_t pointCount = 2 + static_cast<std::size_t>(outputNx) * rows;
  const std::size_t cellCount = 2 * static_cast<std::size_t>(outputNx) +
      static_cast<std::size_t>(rows - 1) * outputNx;
  std::ofstream output(path);
  if (!output) throw std::runtime_error("cannot create VTK: " + path);
  output << "# vtk DataFile Version 3.0\nAuthor public reduced sphere\nASCII\n"
         << "DATASET POLYDATA\nPOINTS " << pointCount << " double\n"
         << std::setprecision(12);
  std::array<double, 3> low{
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity()};
  std::array<double, 3> high{-low[0], -low[1], -low[2]};
  auto emit = [&](V3 point) {
    point.x /= scale; point.y /= scale; point.z /= scale;
    output << point.x << ' ' << point.y << ' ' << point.z << '\n';
    const std::array<double, 3> values{point.x, point.y, point.z};
    for (int component = 0; component < 3; ++component) {
      low[component] = std::min(low[component], values[component]);
      high[component] = std::max(high[component], values[component]);
    }
  };
  emit(source[0]);
  emit(source[1]);
  for (int row = 0; row < rows; ++row)
    for (int outputLongitude = 0; outputLongitude < outputNx;
         ++outputLongitude) {
      const int sourceLongitude = outputLongitude * sourceNx / outputNx;
      emit(source[2 + static_cast<std::size_t>(sourceLongitude) * rows + row]);
    }

  output << "POLYGONS " << cellCount << ' '
         << (2 * static_cast<std::size_t>(outputNx) * 4 +
             static_cast<std::size_t>(rows - 1) * outputNx * 5) << '\n';
  auto index = [=](int longitude, int row) {
    longitude = (longitude % outputNx + outputNx) % outputNx;
    return 2 + row * outputNx + longitude;
  };
  for (int longitude = 0; longitude < outputNx; ++longitude)
    output << "3 0 " << index(longitude, 0) << ' '
           << index(longitude + 1, 0) << '\n';
  for (int row = 0; row + 1 < rows; ++row)
    for (int longitude = 0; longitude < outputNx; ++longitude)
      output << "4 " << index(longitude, row) << ' '
             << index(longitude, row + 1) << ' '
             << index(longitude + 1, row + 1) << ' '
             << index(longitude + 1, row) << '\n';
  for (int longitude = 0; longitude < outputNx; ++longitude)
    output << "3 1 " << index(longitude + 1, rows - 1) << ' '
           << index(longitude, rows - 1) << '\n';
  if (!output) throw std::runtime_error("failed writing VTK: " + path);
  std::cout << std::setprecision(12)
            << "reference_wrl_conversion=pass source_grid=" << sourceNx << 'x'
            << rows + 2 << " output_grid=" << outputNx << 'x' << rows + 2
            << " scale=" << scale << " bounds=" << low[0] << ',' << high[0]
            << ';' << low[1] << ',' << high[1] << ';' << low[2] << ','
            << high[2] << " output=" << path << '\n';
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 6) {
      std::cerr << "usage: convert_reference_wrl INPUT.wrl OUTPUT.vtk "
                   "SOURCE_LONGITUDES OUTPUT_LONGITUDES SCALE\n";
      return 2;
    }
    const auto points = readPoints(argv[1]);
    writeVtk(argv[2], points, std::stoi(argv[3]), std::stoi(argv[4]),
             std::stod(argv[5]));
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "convert_reference_wrl: " << error.what() << '\n';
    return 2;
  }
}
