#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {
constexpr double pi = 3.1415926535897932384626433832795;

struct V3 {
  double x = 0, y = 0, z = 0;
  V3 operator-(V3 const& other) const {
    return {x - other.x, y - other.y, z - other.z};
  }
  V3& operator+=(V3 const& other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }
};

V3 operator/(V3 value, double divisor) {
  return {value.x / divisor, value.y / divisor, value.z / divisor};
}

struct Grid {
  int nx = 0, ny = 0;
  double ymax = 0;
  std::vector<V3> points;

  V3 const& at(int i, int j) const {
    i = (i % nx + nx) % nx;
    return points[static_cast<std::size_t>(j) * nx + i];
  }
};

struct Bounds {
  std::array<double, 3> low{
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity()};
  std::array<double, 3> high{
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity()};

  void add(V3 const& point) {
    std::array<double, 3> const value{point.x, point.y, point.z};
    for (int component = 0; component < 3; ++component) {
      low[component] = std::min(low[component], value[component]);
      high[component] = std::max(high[component], value[component]);
    }
  }
};

struct PaperManifest {
  fs::path sourceFile;
  fs::path outputDirectory;
  std::array<int, 2> grid{};
  std::array<int, 3> ridges{};
  double ballRadius = 0;
  double eta = 0;
  double targetFraction = 0;
  std::string stage3Anchor;
  std::string configHash;
};

std::string readText(fs::path const& file) {
  std::ifstream input(file);
  if (!input) throw std::runtime_error("cannot read manifest: " + file.string());
  return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

std::string jsonString(std::string const& json, std::string const& key) {
  std::smatch match;
  std::regex const expression("\\\"" + key +
                              "\\\"\\s*:\\s*\\\"([^\\\"]*)\\\"");
  if (!std::regex_search(json, match, expression))
    throw std::runtime_error("manifest is missing string key: " + key);
  return match[1].str();
}

double jsonNumber(std::string const& json, std::string const& key) {
  std::smatch match;
  std::regex const expression("\\\"" + key +
                              "\\\"\\s*:\\s*([-+0-9.eE]+)");
  if (!std::regex_search(json, match, expression))
    throw std::runtime_error("manifest is missing number key: " + key);
  return std::stod(match[1].str());
}

template <std::size_t Size>
std::array<int, Size> jsonIntegerArray(std::string const& json,
                                      std::string const& key) {
  std::smatch match;
  std::regex const expression("\\\"" + key +
                              "\\\"\\s*:\\s*\\[([^\\]]+)\\]");
  if (!std::regex_search(json, match, expression))
    throw std::runtime_error("manifest is missing integer array key: " + key);
  std::string values = match[1].str();
  std::replace(values.begin(), values.end(), ',', ' ');
  std::istringstream input(values);
  std::array<int, Size> result{};
  for (int& value : result)
    if (!(input >> value))
      throw std::runtime_error("malformed integer array in manifest: " + key);
  int extra = 0;
  if (input >> extra)
    throw std::runtime_error("wrong integer array length in manifest: " + key);
  return result;
}

PaperManifest readManifest(fs::path source, bool requireCertified) {
  if (fs::is_directory(source)) source /= "manifest.json";
  source = fs::weakly_canonical(source);
  std::string const json = readText(source);
  PaperManifest manifest;
  manifest.sourceFile = source;
  manifest.outputDirectory = jsonString(json, "output_directory");
  manifest.grid = jsonIntegerArray<2>(json, "grid");
  manifest.ridges = jsonIntegerArray<3>(json, "ridges");
  manifest.ballRadius = jsonNumber(json, "ball_radius");
  manifest.eta = jsonNumber(json, "eta");
  manifest.targetFraction = jsonNumber(json, "target_fraction");
  manifest.stage3Anchor = jsonString(json, "stage3_anchor");
  manifest.configHash = jsonString(json, "config_hash");
  if (fs::weakly_canonical(manifest.outputDirectory) != source.parent_path())
    throw std::runtime_error("manifest output_directory does not name its own directory");
  if (manifest.grid[0] < 16 || manifest.grid[1] < 16 ||
      !std::all_of(manifest.ridges.begin(), manifest.ridges.end(),
                   [](int ridges) { return ridges > 0; }) ||
      !(manifest.ballRadius > 0) || !std::isfinite(manifest.eta) ||
      !(manifest.targetFraction > 0 && manifest.targetFraction < 1) ||
      manifest.stage3Anchor.empty())
    throw std::runtime_error("manifest has invalid reduced-sphere parameters");
  if (requireCertified &&
      (manifest.grid != std::array<int, 2>{4000, 20000} ||
       manifest.ridges != std::array<int, 3>{21, 142, 997} ||
       std::abs(manifest.ballRadius - 0.52) > 1e-12 ||
       std::abs(manifest.eta - 0.5) > 1e-12 ||
       std::abs(manifest.targetFraction - 0.237) > 1e-12 ||
       manifest.stage3Anchor != "equator"))
    throw std::runtime_error("manifest is not the certified Phase-5 paper configuration");
  return manifest;
}

Grid readGrid(fs::path const& file, std::array<int, 2> const& expectedGrid) {
  std::ifstream input(file, std::ios::binary);
  std::array<char, 8> magic{};
  Grid grid;
  input.read(magic.data(), 8);
  input.read(reinterpret_cast<char*>(&grid.nx), sizeof(grid.nx));
  input.read(reinterpret_cast<char*>(&grid.ny), sizeof(grid.ny));
  input.read(reinterpret_cast<char*>(&grid.ymax), sizeof(grid.ymax));
  if (std::string(magic.data(), 7) != "HEVSPH1" ||
      grid.nx != expectedGrid[0] || grid.ny != expectedGrid[1] ||
      !(grid.ymax > 0))
    throw std::runtime_error("bad or mismatched paper grid: " + file.string());
  grid.points.resize(static_cast<std::size_t>(grid.nx) * grid.ny);
  input.read(reinterpret_cast<char*>(grid.points.data()),
             static_cast<std::streamsize>(grid.points.size() * sizeof(V3)));
  if (!input) throw std::runtime_error("truncated paper grid: " + file.string());
  return grid;
}

Grid blockAverage(Grid const& source, int nx, int ny) {
  nx = std::min(nx, source.nx);
  ny = std::min(ny, source.ny);
  Grid result{nx, ny, source.ymax,
              std::vector<V3>(static_cast<std::size_t>(nx) * ny)};
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      int const i0 = i * source.nx / nx;
      int const i1 = std::max(i0 + 1, (i + 1) * source.nx / nx);
      int const j0 = j * source.ny / ny;
      int const j1 = std::max(j0 + 1, (j + 1) * source.ny / ny);
      V3 sum;
      int count = 0;
      for (int sourceJ = j0; sourceJ < std::min(j1, source.ny); ++sourceJ)
        for (int sourceI = i0; sourceI < std::min(i1, source.nx); ++sourceI) {
          sum += source.at(sourceI, sourceJ);
          ++count;
        }
      result.points[static_cast<std::size_t>(j) * nx + i] = sum / count;
    }
  }
  return result;
}

Bounds writeSurface(fs::path const& file, Grid const& grid, double eta) {
  int const capRows = std::max(12, grid.ny / 16);
  int const rings = 2 * (capRows - 1) + grid.ny;
  std::vector<V3> points;
  points.reserve(static_cast<std::size_t>(rings) * grid.nx + 2);
  points.push_back({0, 0, -1 + eta});
  for (int k = 1; k < capRows; ++k) {
    double const y = -pi / 2 + (pi / 2 - grid.ymax) * k / capRows;
    for (int i = 0; i < grid.nx; ++i) {
      double const x = 2 * pi * i / grid.nx;
      points.push_back({std::cos(y) * std::cos(x), std::cos(y) * std::sin(x),
                        std::sin(y) + eta});
    }
  }
  points.insert(points.end(), grid.points.begin(), grid.points.end());
  for (int k = 1; k < capRows; ++k) {
    double const y = grid.ymax + (pi / 2 - grid.ymax) * k / capRows;
    for (int i = 0; i < grid.nx; ++i) {
      double const x = 2 * pi * i / grid.nx;
      points.push_back({std::cos(y) * std::cos(x), std::cos(y) * std::sin(x),
                        std::sin(y) - eta});
    }
  }
  int const north = static_cast<int>(points.size());
  points.push_back({0, 0, 1 - eta});

  std::ofstream output(file);
  if (!output) throw std::runtime_error("cannot write surface: " + file.string());
  std::size_t const cells = static_cast<std::size_t>(rings - 1) * grid.nx +
                            2 * static_cast<std::size_t>(grid.nx);
  output << "# vtk DataFile Version 3.0\nHevea reduced sphere grid\nASCII\n"
         << "DATASET POLYDATA\nPOINTS " << points.size() << " double\n"
         << std::setprecision(12);
  Bounds bounds;
  for (V3 const& point : points) {
    output << point.x << ' ' << point.y << ' ' << point.z << '\n';
    bounds.add(point);
  }
  output << "POLYGONS " << cells << ' '
         << (static_cast<std::size_t>(rings - 1) * grid.nx * 5 +
             2 * static_cast<std::size_t>(grid.nx) * 4)
         << '\n';
  for (int i = 0; i < grid.nx; ++i)
    output << "3 0 " << 1 + (i + 1) % grid.nx << ' ' << 1 + i << '\n';
  for (int j = 0; j < rings - 1; ++j) {
    for (int i = 0; i < grid.nx; ++i) {
      int const a = 1 + j * grid.nx + i;
      int const b = 1 + j * grid.nx + (i + 1) % grid.nx;
      output << "4 " << a << ' ' << b << ' ' << b + grid.nx << ' '
             << a + grid.nx << '\n';
    }
  }
  int const last = 1 + (rings - 1) * grid.nx;
  for (int i = 0; i < grid.nx; ++i)
    output << "3 " << last + i << ' ' << last + (i + 1) % grid.nx << ' '
           << north << '\n';
  if (!output) throw std::runtime_error("failed writing surface: " + file.string());
  return bounds;
}

void writePatch(fs::path const& file, Grid const& grid, int stage) {
  int const width = std::min(grid.nx, 360);
  int const height = std::min(grid.ny, stage < 3 ? 1200 : 800);
  int const i0 = (grid.nx - width) / 2;
  int const j0 = (grid.ny - height) / 2;
  std::ofstream output(file);
  if (!output) throw std::runtime_error("cannot write detail patch: " + file.string());
  std::size_t const cells = static_cast<std::size_t>(height - 1) * (width - 1);
  output << "# vtk DataFile Version 3.0\nHevea unfiltered stage " << stage
         << " detail\nASCII\nDATASET POLYDATA\nPOINTS "
         << static_cast<std::size_t>(width) * height << " double\n"
         << std::setprecision(12);
  for (int j = 0; j < height; ++j)
    for (int i = 0; i < width; ++i) {
      V3 const point = grid.at(i0 + i, j0 + j);
      output << point.x << ' ' << point.y << ' ' << point.z << '\n';
    }
  output << "POLYGONS " << cells << ' ' << cells * 5 << '\n';
  for (int j = 0; j < height - 1; ++j)
    for (int i = 0; i < width - 1; ++i) {
      int const a = j * width + i;
      output << "4 " << a << ' ' << a + 1 << ' ' << a + 1 + width << ' '
             << a + width << '\n';
    }
  if (!output) throw std::runtime_error("failed writing detail patch: " + file.string());
}

double spectralEnergy(std::vector<V3> const& signal, int frequency) {
  std::array<std::complex<double>, 3> coefficient{};
  for (std::size_t index = 0; index < signal.size(); ++index) {
    double const angle = -2 * pi * frequency * index / signal.size();
    std::complex<double> const phase{std::cos(angle), std::sin(angle)};
    coefficient[0] += signal[index].x * phase;
    coefficient[1] += signal[index].y * phase;
    coefficient[2] += signal[index].z * phase;
  }
  return std::norm(coefficient[0]) + std::norm(coefficient[1]) +
         std::norm(coefficient[2]);
}

int localSpectralPeak(std::vector<V3> const& signal, int expected) {
  int best = 0;
  double bestEnergy = -1;
  int const radius = std::max(8, expected / 20);
  for (int frequency = std::max(1, expected - radius);
       frequency <= expected + radius; ++frequency) {
    double const energy = spectralEnergy(signal, frequency);
    if (energy > bestEnergy) {
      bestEnergy = energy;
      best = frequency;
    }
  }
  return best;
}

std::vector<V3> longitudeIncrement(Grid const& before, Grid const& after) {
  int const row = before.ny / 2;
  std::vector<V3> signal(before.nx);
  for (int i = 0; i < before.nx; ++i) {
    V3 const value = after.at(i, row) - before.at(i, row);
    double const x = 2 * pi * i / before.nx;
    signal[i] = {value.x * std::cos(x) + value.y * std::sin(x),
                 -value.x * std::sin(x) + value.y * std::cos(x), value.z};
  }
  return signal;
}

V3 interpolateY(Grid const& grid, int column, double y) {
  double const coordinate = (y + grid.ymax) * (grid.ny - 1) / (2 * grid.ymax);
  int const lower = std::clamp(static_cast<int>(std::floor(coordinate)), 0, grid.ny - 2);
  double const fraction = std::clamp(coordinate - lower, 0.0, 1.0);
  V3 const a = grid.at(column, lower);
  V3 const b = grid.at(column, lower + 1);
  return {a.x + fraction * (b.x - a.x),
          a.y + fraction * (b.y - a.y),
          a.z + fraction * (b.z - a.z)};
}

std::vector<V3> latitudeIncrement(Grid const& before, Grid const& after,
                                  int ridgeCount) {
  double const outer = ridgeCount / (2 * 334.92);
  int const samples = 4 * ridgeCount;
  int const column = before.nx / 2;
  std::vector<V3> signal(samples);
  for (int index = 0; index < samples; ++index) {
    double const y = -outer + 2 * outer * index / samples;
    signal[index] = interpolateY(after, column, y) -
                    interpolateY(before, column, y);
  }
  return signal;
}

double directionEnergy(Grid const& before, Grid const& after, int ridgeCount,
                       int ySign) {
  std::array<std::complex<double>, 3> coefficient{};
  int constexpr rows = 160;
  int constexpr columns = 500;
  for (int jSample = 0; jSample < rows; ++jSample) {
    int const j = before.ny / 2 - before.ny / 10 +
                  (2 * before.ny / 10) * jSample / (rows - 1);
    double const y = -before.ymax + 2 * before.ymax * j / (before.ny - 1);
    for (int iSample = 0; iSample < columns; ++iSample) {
      int const i = iSample * before.nx / columns;
      double const x = 2 * pi * i / before.nx;
      V3 const global = after.at(i, j) - before.at(i, j);
      V3 const value{global.x * std::cos(x) + global.y * std::sin(x),
                     -global.x * std::sin(x) + global.y * std::cos(x), global.z};
      double const angle = -ridgeCount * (x + ySign * y);
      std::complex<double> const phase{std::cos(angle), std::sin(angle)};
      coefficient[0] += value.x * phase;
      coefficient[1] += value.y * phase;
      coefficient[2] += value.z * phase;
    }
  }
  return std::norm(coefficient[0]) + std::norm(coefficient[1]) +
         std::norm(coefficient[2]);
}

std::string jsonEscape(std::string const& text) {
  std::ostringstream output;
  for (char character : text) {
    if (character == '\\' || character == '"') output << '\\';
    output << character;
  }
  return output.str();
}

void writeExportManifest(fs::path const& file, PaperManifest const& source,
                         std::vector<fs::path> const& files,
                         fs::path const& fullResolutionFinal,
                         std::array<int, 3> const& counted,
                         std::array<double, 2> const& directionRatios,
                         Bounds const& globalBounds, double sourceYmax) {
  int constexpr globalNx = 400;
  int constexpr globalNy = 2000;
  int constexpr capRows = globalNy / 16;
  int constexpr rings = 2 * (capRows - 1) + globalNy;
  std::size_t constexpr pointCount = static_cast<std::size_t>(rings) * globalNx + 2;
  std::size_t constexpr cellCount =
      static_cast<std::size_t>(rings - 1) * globalNx + 2 * globalNx;
  std::ofstream output(file);
  if (!output) throw std::runtime_error("cannot write export manifest");
  output << std::setprecision(17)
         << "{\n  \"source_manifest\": \""
         << jsonEscape(source.sourceFile.string()) << "\",\n"
         << "  \"source_config_hash\": \"" << source.configHash << "\",\n"
         << "  \"full_resolution_final\": \""
         << jsonEscape(fullResolutionFinal.string()) << "\",\n"
         << "  \"eta\": " << source.eta << ",\n"
         << "  \"source_ymax\": " << sourceYmax << ",\n"
         << "  \"global_grid\": [" << globalNx << ", " << globalNy << "],\n"
         << "  \"cap_rows\": " << capRows << ",\n"
         << "  \"topology\": {\"connected_components\": 1, \"closed\": true, "
         << "\"euler_characteristic\": 2, \"point_count\": " << pointCount
         << ", \"cell_count\": " << cellCount << "},\n"
         << "  \"ridges_expected\": [" << source.ridges[0] << ", "
         << source.ridges[1] << ", " << source.ridges[2] << "],\n"
         << "  \"ridges_counted\": [" << counted[0] << ", " << counted[1]
         << ", " << counted[2] << "],\n"
         << "  \"direction_energy_ratios\": [" << directionRatios[0] << ", "
         << directionRatios[1] << "],\n"
         << "  \"global_bounds\": [[" << globalBounds.low[0] << ", "
         << globalBounds.high[0] << "], [" << globalBounds.low[1] << ", "
         << globalBounds.high[1] << "], [" << globalBounds.low[2] << ", "
         << globalBounds.high[2] << "]],\n  \"files\": [\n";
  for (std::size_t index = 0; index < files.size(); ++index)
    output << "    \"" << jsonEscape(files[index].string()) << "\""
           << (index + 1 == files.size() ? "\n" : ",\n");
  output << "  ]\n}\n";
  if (!output) throw std::runtime_error("failed writing export manifest");
}
}  // namespace

int main(int argc, char** argv) {
  try {
    bool const diagnostic = argc == 4 && std::string(argv[1]) == "--diagnostic";
    if ((!diagnostic && argc != 3) || (diagnostic && argc != 4)) {
      std::cerr << "usage: export_paper_meshes [--diagnostic] "
                   "MANIFEST_OR_INPUT_DIR OUTPUT_DIR\n";
      return 2;
    }
    int const argumentOffset = diagnostic ? 1 : 0;
    PaperManifest const manifest = readManifest(argv[1 + argumentOffset],
                                                !diagnostic);
    fs::path const outputDirectory = fs::absolute(argv[2 + argumentOffset]);
    fs::create_directories(outputDirectory);
    auto stageFile = [&](int stage) {
      if (stage == 0)
        return manifest.outputDirectory / "reduced_sphere_stage=0.bin";
      return manifest.outputDirectory /
             ("reduced_sphere_stage=" + std::to_string(stage) + "_dir=" +
              std::to_string(stage - 1) + "_ridges=" +
              std::to_string(manifest.ridges[stage - 1]) + ".bin");
    };

    std::vector<fs::path> exportedFiles;
    Grid previous = readGrid(stageFile(0), manifest.grid);
    fs::path const stage0Global = outputDirectory / "global-stage0.vtk";
    writeSurface(stage0Global, blockAverage(previous, 300, 1000), manifest.eta);
    exportedFiles.push_back(stage0Global);

    std::array<int, 3> counted{};
    std::array<double, 2> directionRatios{};
    Grid final;
    for (int stage = 1; stage <= 3; ++stage) {
      Grid current = readGrid(stageFile(stage), manifest.grid);
      std::vector<V3> const signal = stage < 3
          ? longitudeIncrement(previous, current)
          : latitudeIncrement(previous, current, manifest.ridges[2]);
      counted[stage - 1] = localSpectralPeak(signal, manifest.ridges[stage - 1]);
      if (stage < 3) {
        int const expectedSign = stage == 1 ? 1 : -1;
        double const expected = directionEnergy(previous, current,
                                                manifest.ridges[stage - 1],
                                                expectedSign);
        double const opposite = directionEnergy(previous, current,
                                                manifest.ridges[stage - 1],
                                                -expectedSign);
        directionRatios[stage - 1] = expected / std::max(opposite, 1e-300);
      }

      fs::path const global = outputDirectory /
                              ("global-stage" + std::to_string(stage) + ".vtk");
      writeSurface(global, blockAverage(current, 300, 1000), manifest.eta);
      exportedFiles.push_back(global);
      fs::path const detail = outputDirectory /
          ("detail-stage" + std::to_string(stage) + "-ridges" +
           std::to_string(manifest.ridges[stage - 1]) + ".vtk");
      writePatch(detail, current, stage);
      exportedFiles.push_back(detail);
      if (stage == 3) final = std::move(current);
      else previous = std::move(current);
    }

    if (counted != manifest.ridges || directionRatios[0] <= 1 ||
        directionRatios[1] <= 1) {
      std::ostringstream diagnostic;
      diagnostic << "geometric ridge audit failed: counted=" << counted[0] << ','
                 << counted[1] << ',' << counted[2] << " expected="
                 << manifest.ridges[0] << ',' << manifest.ridges[1] << ','
                 << manifest.ridges[2] << " direction_ratios="
                 << directionRatios[0] << ',' << directionRatios[1];
      throw std::runtime_error(diagnostic.str());
    }

    fs::path const finalPreview =
        outputDirectory / "sampled-reduced-sphere-stage3.vtk";
    Bounds const bounds =
        writeSurface(finalPreview, blockAverage(final, 400, 2000), manifest.eta);
    exportedFiles.push_back(finalPreview);
    fs::path const exportManifest = outputDirectory / "figure9-mesh-manifest.json";
    writeExportManifest(exportManifest, manifest, exportedFiles, stageFile(3), counted,
                        directionRatios, bounds, final.ymax);
    std::cout << std::setprecision(8)
              << "export=pass mode=" << (diagnostic ? "diagnostic" : "certified")
              << " source_config_hash=" << manifest.configHash
              << " final=400x2000 detail=unfiltered eta=" << manifest.eta
              << " ridges_counted=" << counted[0] << ',' << counted[1] << ','
              << counted[2] << " direction_ratios=" << directionRatios[0] << ','
              << directionRatios[1] << " manifest=" << exportManifest << '\n';
    return 0;
  } catch (std::exception const& error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 1;
  }
}
