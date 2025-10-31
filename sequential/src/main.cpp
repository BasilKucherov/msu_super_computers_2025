#include "arg_parser.hpp"
#include "solver.hpp"
#include "solver_params.hpp"
#include "variant_config.hpp"
#include <iostream>

int main(int argc, char *argv[]) {
  try {
    utils::CmdArgs args = utils::ArgParser::parse(argc, argv);
    utils::VariantConfig config =
        utils::VariantConfig::variant6(args.Lx, args.Ly, args.Lz);
    sequential::SolverParams params(args, config);
    sequential::SequentialSolver solver(params);
    solver.solve();
  } catch (const utils::ParserException &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}