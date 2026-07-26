#include "WolframLibrary.h"

#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <filesystem>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

namespace {
std::string executablePath() {
  if (char const* overridePath = std::getenv("HEVEA_REDUCED_SPHERE_EXECUTABLE")) {
    if (*overridePath != '\0') return std::filesystem::absolute(overridePath).string();
  }
  Dl_info info{};
  if (dladdr(reinterpret_cast<void*>(&executablePath), &info) == 0 ||
      info.dli_fname == nullptr)
    return "hevea_reduced_sphere";
  return (std::filesystem::absolute(info.dli_fname).parent_path() /
          "hevea_reduced_sphere").string();
}

std::string jsonEscape(std::string const& text) {
  std::ostringstream out;
  for (unsigned char ch : text) {
    switch (ch) {
      case '\\': out << "\\\\"; break;
      case '"': out << "\\\""; break;
      case '\n': out << "\\n"; break;
      case '\r': out << "\\r"; break;
      case '\t': out << "\\t"; break;
      default:
        if (ch < 0x20) out << "\\u00" << "0123456789abcdef"[ch >> 4]
                           << "0123456789abcdef"[ch & 15];
        else out << static_cast<char>(ch);
    }
  }
  return out.str();
}

std::string resultJson(int exitCode, std::string const& outputDirectory,
                       std::string const& log) {
  std::ostringstream json;
  json << "{\"ExitCode\":" << exitCode << ",\"OutputDirectory\":\""
       << jsonEscape(outputDirectory) << "\",\"Log\":\"" << jsonEscape(log) << "\"}";
  return json.str();
}

int setResultString(MArgument result, std::string const& value) {
  char* copy = static_cast<char*>(std::malloc(value.size() + 1));
  if (!copy) return LIBRARY_MEMORY_ERROR;
  std::memcpy(copy, value.c_str(), value.size() + 1);
  MArgument_setUTF8String(result, copy);
  return LIBRARY_NO_ERROR;
}
}

extern "C" {
DLLEXPORT mint WolframLibrary_getVersion(void) { return WolframLibraryVersion; }
DLLEXPORT int WolframLibrary_initialize(WolframLibraryData) { return LIBRARY_NO_ERROR; }
DLLEXPORT void WolframLibrary_uninitialize(WolframLibraryData) {}

DLLEXPORT int heveaReducedSphereRun_mma(WolframLibraryData libData, mint argc,
                                        MArgument* args, MArgument result) {
  if (argc != 12) return LIBRARY_FUNCTION_ERROR;
  char* rawOutput = MArgument_getUTF8String(args[0]);
  std::string outputDirectory = rawOutput ? rawOutput : "";
  libData->UTF8String_disown(rawOutput);
  mint nx = MArgument_getInteger(args[1]), ny = MArgument_getInteger(args[2]);
  mint r0 = MArgument_getInteger(args[3]), r1 = MArgument_getInteger(args[4]);
  mint r2 = MArgument_getInteger(args[5]);
  double ballRadius = MArgument_getReal(args[6]);
  double eta = MArgument_getReal(args[7]);
  double targetFraction = MArgument_getReal(args[8]);
  mint timeLimit = MArgument_getInteger(args[9]);
  mint outputMode = MArgument_getInteger(args[10]);
  mint requiredFreeBytes = MArgument_getInteger(args[11]);
  if (outputDirectory.empty() || nx < 32 || ny < 64 || r0 < 1 || r1 < 1 || r2 < 1 ||
      !std::isfinite(ballRadius) || ballRadius <= 0.1 || ballRadius >= 1.0 ||
      !std::isfinite(eta) || eta <= 0.0 || eta >= 1.5 ||
      !std::isfinite(targetFraction) || targetFraction <= 0.0 || targetFraction > 1.0 ||
      timeLimit < 1 || timeLimit > 86400 || outputMode < 0 || outputMode > 1 ||
      requiredFreeBytes < 0) return LIBRARY_TYPE_ERROR;
  std::string const executable = executablePath();
  try {
    std::filesystem::create_directories(outputDirectory);
    const auto space = std::filesystem::space(outputDirectory);
    if (space.available < static_cast<std::uintmax_t>(requiredFreeBytes))
      return setResultString(result, resultJson(75, outputDirectory,
        "ERROR: insufficient disk space: available="+std::to_string(space.available)+
        " required="+std::to_string(requiredFreeBytes)+"\n"));
  }
  catch (...) { return LIBRARY_FUNCTION_ERROR; }

  int pipefd[2];
  if (pipe(pipefd) != 0) return LIBRARY_FUNCTION_ERROR;
  pid_t pid = fork();
  if (pid < 0) { close(pipefd[0]); close(pipefd[1]); return LIBRARY_FUNCTION_ERROR; }
  if (pid == 0) {
    close(pipefd[0]); dup2(pipefd[1], STDOUT_FILENO); dup2(pipefd[1], STDERR_FILENO);
    close(pipefd[1]);
    if (chdir(outputDirectory.c_str()) != 0) _exit(126);
    setenv("OMP_NUM_THREADS", "16", 1); setenv("OMP_DYNAMIC", "FALSE", 1);
    if (outputMode == 1) setenv("HEVEA_PREVIEW_ONLY", "1", 1);
    else unsetenv("HEVEA_PREVIEW_ONLY");
    alarm(static_cast<unsigned int>(timeLimit));
    std::string sx = std::to_string(nx), sy = std::to_string(ny);
    std::string sr0 = std::to_string(r0), sr1 = std::to_string(r1), sr2 = std::to_string(r2);
    std::string sb = std::to_string(ballRadius), se = std::to_string(eta);
    std::string st = std::to_string(targetFraction);
    execl(executable.c_str(), executable.c_str(), sx.c_str(), sy.c_str(), sr0.c_str(), sr1.c_str(),
          sr2.c_str(), sb.c_str(), se.c_str(), st.c_str(), static_cast<char*>(nullptr));
    _exit(errno == ENOENT ? 127 : 126);
  }
  close(pipefd[1]);
  std::string log; char buffer[8192]; ssize_t count;
  while ((count = read(pipefd[0], buffer, sizeof(buffer))) > 0)
    log.append(buffer, static_cast<std::size_t>(count));
  close(pipefd[0]);
  int status = 0; waitpid(pid, &status, 0);
  int exitCode = WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
  return setResultString(result, resultJson(exitCode, outputDirectory, log));
}
}
