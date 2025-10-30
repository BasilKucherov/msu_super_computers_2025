#ifndef SOLVER_PARAMS_HPP
#define SOLVER_PARAMS_HPP

#include "arg_parser.hpp"
#include "variant_config.hpp"

namespace sequential {
struct SolverParams {
  double Lx, Ly, Lz;
  double hx, hy, hz;

  double T;
  double tau;
  int K;

  double a_squared;
  double a_t;

  std::function<double(double, double, double, double)> u_analytical;

  SolverParams(const utils::CmdArgs &args, const utils::VariantConfig &config)
      : Lx(args.Lx), Ly(args.Ly), Lz(args.Lz), hx(args.hx), hy(args.hy),
        hz(args.hz), T(args.T), tau(args.tau), K(args.K),
        a_squared(config.a_squared), a_t(config.a_t),
        u_analytical(config.u_analytical) {}
};
} // namespace sequential

#endif
