#ifndef MPI_SOLVER_PARAMS_HPP
#define MPI_SOLVER_PARAMS_HPP

#include "arg_parser.hpp"
#include "variant_config.hpp"

namespace mpi_solver {
struct SolverParams {
  double Lx, Ly, Lz;
  int N;

  double T;
  double tau;
  int K;

  double a_squared;
  double a_t;
  bool debug;

  std::function<double(double, double, double, double)> u_analytical;

  SolverParams(const utils::CmdArgs &args, const utils::VariantConfig &config)
      : Lx(args.Lx), Ly(args.Ly), Lz(args.Lz), N(args.N), T(args.T),
        tau(args.tau), K(args.K), a_squared(config.a_squared), a_t(config.a_t),
        debug(args.debug), u_analytical(config.u_analytical) {}
};
} // namespace mpi_solver

#endif
