#ifndef ARG_PARSER_HPP
#define ARG_PARSER_HPP

#include <optional>
#include <stdexcept>
#include <string>

namespace utils {
struct CmdArgs {
  double Lx;
  double Ly;
  double Lz;

  std::optional<int> Nx;
  std::optional<int> Ny;
  std::optional<int> Nz;

  double T;
  double tau;
  int K;
  int N;
  bool debug;
};

class ParserException : public std::runtime_error {
public:
  explicit ParserException(const std::string &message)
      : std::runtime_error(message) {}
};

class ArgParser {
public:
  static CmdArgs parse(int argc, char *argv[], bool parse_mpi_grid = false);

private:
  static std::string getHelpMessage(bool mpi_mode);
  static void validateParams(const CmdArgs &params, bool hasTau, bool hasK,
                             bool parse_mpi_grid);
  static void computeMissingTimeParam(CmdArgs &params, bool hasTau, bool hasK);
};
} // namespace utils

#endif
