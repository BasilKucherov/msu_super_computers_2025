#include "arg_parser.hpp"

#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>

namespace utils {
CmdArgs ArgParser::parse(int argc, char *argv[], bool parse_mpi_grid) {
  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "--help") == 0 ||
        std::strcmp(argv[i], "-h") == 0) {
      std::cout << getHelpMessage(parse_mpi_grid) << std::endl;
      throw ParserException("Help requested");
    }
  }

  CmdArgs params = {};

  bool hasLx = false, hasLy = false, hasLz = false;
  bool hasHx = false, hasHy = false, hasHz = false;
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
      } else if (arg == "--hx") {
        params.hx = std::stod(value);
        hasHx = true;
        ++i;
      } else if (arg == "--hy") {
        params.hy = std::stod(value);
        hasHy = true;
        ++i;
      } else if (arg == "--hz") {
        params.hz = std::stod(value);
        hasHz = true;
        ++i;
      } else if (arg == "--Nx") {
        if (!parse_mpi_grid) {
          throw ParserException("--Nx is only valid in MPI mode");
        }
        params.Nx = std::stoi(value);
        hasNx = true;
        ++i;
      } else if (arg == "--Ny") {
        if (!parse_mpi_grid) {
          throw ParserException("--Ny is only valid in MPI mode");
        }
        params.Ny = std::stoi(value);
        hasNy = true;
        ++i;
      } else if (arg == "--Nz") {
        if (!parse_mpi_grid) {
          throw ParserException("--Nz is only valid in MPI mode");
        }
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
  if (!hasHx)
    missing << " --hx";
  if (!hasHy)
    missing << " --hy";
  if (!hasHz)
    missing << " --hz";
  if (!hasT)
    missing << " --T";

  if (parse_mpi_grid) {
    if (!hasNx)
      missing << " --Nx";
    if (!hasNy)
      missing << " --Ny";
    if (!hasNz)
      missing << " --Nz";
  }

  if (!missing.str().empty()) {
    throw ParserException("Missing required arguments:" + missing.str());
  }
  validateParams(params, hasTau, hasK, parse_mpi_grid);
  computeMissingTimeParam(params, hasTau, hasK);

  return params;
}

std::string ArgParser::getHelpMessage(bool mpi_mode) {
  std::ostringstream help;
  help << "Wave Equation Solver - Command Line Arguments\n"
       << "==============================================\n\n"
       << "REQUIRED ARGUMENTS:\n"
       << "  --Lx <value>    Domain length in X direction (double, > 0)\n"
       << "  --Ly <value>    Domain length in Y direction (double, > 0)\n"
       << "  --Lz <value>    Domain length in Z direction (double, > 0)\n"
       << "  --hx <value>    Spatial step size in X direction (double, > 0)\n"
       << "  --hy <value>    Spatial step size in Y direction (double, > 0)\n"
       << "  --hz <value>    Spatial step size in Z direction (double, > 0)\n"
       << "  --T <value>     Total simulation time (double, > 0)\n\n";

  if (mpi_mode) {
    help
        << "MPI PROCESS GRID (required in MPI mode):\n"
        << "  --Nx <value>    Number of MPI processes in X direction (int32, > "
           "0)\n"
        << "  --Ny <value>    Number of MPI processes in Y direction (int32, > "
           "0)\n"
        << "  --Nz <value>    Number of MPI processes in Z direction (int32, > "
           "0)\n\n";
  }

  help << "TIME STEPPING (exactly one required):\n"
       << "  --tau <value>   Time step size (double, > 0)\n"
       << "                  If provided, K will be computed as K = ceil(T / "
          "tau)\n"
       << "  --K <value>     Number of time steps (int32, > 0)\n"
       << "                  If provided, tau will be computed as tau = T / "
          "K\n\n"
       << "OTHER:\n"
       << "  --help, -h      Display this help message\n\n"
       << "EXAMPLE USAGE:\n";

  if (mpi_mode) {
    help << "  mpirun -np 8 ./solver --Lx 1.0 --Ly 1.0 --Lz 1.0 --hx 0.01 "
            "--hy 0.01 --hz 0.01 --Nx 2 --Ny 2 --Nz 2 --T 10.0 --tau 0.01\n";
  } else {
    help << "  ./solver --Lx 1.0 --Ly 1.0 --Lz 1.0 --hx 0.01 --hy 0.01 --hz "
            "0.01 --T 10.0 --tau 0.01\n"
         << "  ./solver --Lx 2.0 --Ly 2.0 --Lz 2.0 --hx 0.02 --hy 0.02 --hz "
            "0.02 --T 5.0 --K 1000\n";
  }

  return help.str();
}

void ArgParser::validateParams(const CmdArgs &params, bool hasTau, bool hasK,
                               bool parse_mpi_grid) {
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

  if (params.hx <= 0.0) {
    throw ParserException("hx must be positive (got " +
                          std::to_string(params.hx) + ")");
  }
  if (params.hy <= 0.0) {
    throw ParserException("hy must be positive (got " +
                          std::to_string(params.hy) + ")");
  }
  if (params.hz <= 0.0) {
    throw ParserException("hz must be positive (got " +
                          std::to_string(params.hz) + ")");
  }

  if (parse_mpi_grid) {
    if (params.Nx.value() <= 0) {
      throw ParserException("Nx must be positive (got " +
                            std::to_string(params.Nx.value()) + ")");
    }
    if (params.Ny.value() <= 0) {
      throw ParserException("Ny must be positive (got " +
                            std::to_string(params.Ny.value()) + ")");
    }
    if (params.Nz.value() <= 0) {
      throw ParserException("Nz must be positive (got " +
                            std::to_string(params.Nz.value()) + ")");
    }
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

void ArgParser::computeMissingTimeParam(CmdArgs &params, bool hasTau,
                                        bool hasK) {
  if (hasTau) {
    params.K = static_cast<int32_t>(std::ceil(params.T / params.tau));
  } else if (hasK) {
    params.tau = params.T / static_cast<double>(params.K);
  }
}
} // namespace utils
