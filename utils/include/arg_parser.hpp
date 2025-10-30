#ifndef ARG_PARSER_HPP
#define ARG_PARSER_HPP

#include <stdexcept>
#include <string>

struct WaveParams {
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
  static WaveParams parse(int argc, char *argv[]);

private:
  static std::string getHelpMessage();
  static void validateParams(const WaveParams &params, bool hasTau, bool hasK);
  static void computeMissingTimeParam(WaveParams &params, bool hasTau,
                                      bool hasK);
};

#endif
