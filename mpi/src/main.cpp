#include "arg_parser.hpp"
#include "mpi_solver.hpp"
#include "mpi_solver_params.hpp"
#include "variant_config.hpp"
#include <iostream>
#include <mpi.h>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int retval = 0;

  try {
    utils::CmdArgs args = utils::ArgParser::parse(argc, argv);
    utils::VariantConfig config =
        utils::VariantConfig::variant6(args.Lx, args.Ly, args.Lz);
    mpi_solver::SolverParams params(args, config);
    mpi_solver::MPISolver solver(params);
    solver.solve();
  } catch (const std::exception &e) {
    if (rank == 0) {
      std::cerr << "Error: " << e.what() << std::endl;
    }
    retval = 1;
  }

  MPI_Finalize();
  return retval;
}
