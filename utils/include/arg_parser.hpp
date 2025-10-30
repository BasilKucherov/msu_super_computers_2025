#ifndef ARG_PARSER_HPP
#define ARG_PARSER_HPP

#include <stdexcept>
#include <string>

namespace utils {
struct CmdArgs {
  double Lx;
  double Ly;
  double Lz;

  int Nx;
  int Ny;
  int Nz;

  double T;
  double tau;
  int K;
};

class ParserException : public std::runtime_error {
public:
  explicit ParserException(const std::string &message)
      : std::runtime_error(message) {}
};

class ArgParser {
public:
  static CmdArgs parse(int argc, char *argv[]);

private:
  static std::string getHelpMessage();
  static void validateParams(const CmdArgs &params, bool hasTau, bool hasK);
  static void computeMissingTimeParam(CmdArgs &params, bool hasTau, bool hasK);
};
} // namespace utils

#endif
