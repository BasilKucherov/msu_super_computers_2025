#include "arg_parser.hpp"

#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>

namespace utils {
WaveParams ArgParser::parse(int argc, char *argv[]) {
  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "--help") == 0 ||
        std::strcmp(argv[i], "-h") == 0) {
      std::cout << getHelpMessage() << std::endl;
      throw ParserException("Help requested");
    }
  }

  WaveParams params = {};

  bool hasLx = false, hasLy = false, hasLz = false;
  bool hasNx = false, hasNy = false, hasNz = false;
  bool hasT = false;
  bool hasTau = false, hasK = false;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];

    if (i + 1 >= argc) {
      throw ParserException("Missing value for argument: " + arg);
    }

    std::string value = argv[i + 1];

    try {
      if (arg == "--Lx") {
        params.Lx = std::stod(value);
        hasLx = true;
        ++i;
      } else if (arg == "--Ly") {
        params.Ly = std::stod(value);
        hasLy = true;
        ++i;
      } else if (arg == "--Lz") {
        params.Lz = std::stod(value);
        hasLz = true;
        ++i;
      } else if (arg == "--Nx") {
        params.Nx = std::stoi(value);
        hasNx = true;
        ++i;
      } else if (arg == "--Ny") {
        params.Ny = std::stoi(value);
        hasNy = true;
        ++i;
      } else if (arg == "--Nz") {
        params.Nz = std::stoi(value);
        hasNz = true;
        ++i;
      } else if (arg == "--T") {
        params.T = std::stod(value);
        hasT = true;
        ++i;
      } else if (arg == "--tau") {
        params.tau = std::stod(value);
        hasTau = true;
        ++i;
      } else if (arg == "--K") {
        params.K = std::stoi(value);
        hasK = true;
        ++i;
      } else {
        throw ParserException("Unknown argument: " + arg);
      }
    } catch (const std::invalid_argument &e) {
      throw ParserException("Invalid value for " + arg + ": " + value);
    } catch (const std::out_of_range &e) {
      throw ParserException("Value out of range for " + arg + ": " + value);
    }
  }

  std::ostringstream missing;
  if (!hasLx)
    missing << " --Lx";
  if (!hasLy)
    missing << " --Ly";
  if (!hasLz)
    missing << " --Lz";
  if (!hasNx)
    missing << " --Nx";
  if (!hasNy)
    missing << " --Ny";
  if (!hasNz)
    missing << " --Nz";
  if (!hasT)
    missing << " --T";

  if (!missing.str().empty()) {
    throw ParserException("Missing required arguments:" + missing.str());
  }
  validateParams(params, hasTau, hasK);
  computeMissingTimeParam(params, hasTau, hasK);

  return params;
}

std::string ArgParser::getHelpMessage() {
  std::ostringstream help;
  help
      << "Wave Equation Solver - Command Line Arguments\n"
      << "==============================================\n\n"
      << "REQUIRED ARGUMENTS:\n"
      << "  --Lx <value>    Domain length in X direction (double, > 0)\n"
      << "  --Ly <value>    Domain length in Y direction (double, > 0)\n"
      << "  --Lz <value>    Domain length in Z direction (double, > 0)\n"
      << "  --Nx <value>    Number of grid points in X direction (int32, > 0)\n"
      << "  --Ny <value>    Number of grid points in Y direction (int32, > 0)\n"
      << "  --Nz <value>    Number of grid points in Z direction (int32, > 0)\n"
      << "  --T <value>     Total simulation time (double, > 0)\n\n"
      << "TIME STEPPING (exactly one required):\n"
      << "  --tau <value>   Time step size (double, > 0)\n"
      << "                  If provided, K will be computed as K = ceil(T / "
         "tau)\n"
      << "  --K <value>     Number of time steps (int32, > 0)\n"
      << "                  If provided, tau will be computed as tau = T / "
         "K\n\n"
      << "OTHER:\n"
      << "  --help, -h      Display this help message\n\n"
      << "EXAMPLE USAGE:\n"
      << "  ./solver --Lx 1.0 --Ly 1.0 --Lz 1.0 --Nx 100 --Ny 100 --Nz 100 --T "
         "10.0 --tau 0.01\n"
      << "  ./solver --Lx 2.0 --Ly 2.0 --Lz 2.0 --Nx 200 --Ny 200 --Nz 200 --T "
         "5.0 --K 1000\n";

  return help.str();
}

void ArgParser::validateParams(const WaveParams &params, bool hasTau,
                               bool hasK) {
  if (params.Lx <= 0.0) {
    throw ParserException("Lx must be positive (got " +
                          std::to_string(params.Lx) + ")");
  }
  if (params.Ly <= 0.0) {
    throw ParserException("Ly must be positive (got " +
                          std::to_string(params.Ly) + ")");
  }
  if (params.Lz <= 0.0) {
    throw ParserException("Lz must be positive (got " +
                          std::to_string(params.Lz) + ")");
  }

  if (params.Nx <= 0) {
    throw ParserException("Nx must be positive (got " +
                          std::to_string(params.Nx) + ")");
  }
  if (params.Ny <= 0) {
    throw ParserException("Ny must be positive (got " +
                          std::to_string(params.Ny) + ")");
  }
  if (params.Nz <= 0) {
    throw ParserException("Nz must be positive (got " +
                          std::to_string(params.Nz) + ")");
  }

  if (params.T <= 0.0) {
    throw ParserException("T must be positive (got " +
                          std::to_string(params.T) + ")");
  }

  if (hasTau && params.tau <= 0.0) {
    throw ParserException("tau must be positive (got " +
                          std::to_string(params.tau) + ")");
  }
  if (hasK && params.K <= 0) {
    throw ParserException("K must be positive (got " +
                          std::to_string(params.K) + ")");
  }

  if (!hasTau && !hasK) {
    throw ParserException("Exactly one of --tau or --K must be provided");
  }
  if (hasTau && hasK) {
    throw ParserException(
        "Cannot provide both --tau and --K. Provide exactly one.");
  }
}

void ArgParser::computeMissingTimeParam(WaveParams &params, bool hasTau,
                                        bool hasK) {
  if (hasTau) {
    params.K = static_cast<int32_t>(std::ceil(params.T / params.tau));
  } else if (hasK) {
    params.tau = params.T / static_cast<double>(params.K);
  }
}
} // namespace utils
