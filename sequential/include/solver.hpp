#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "solver_params.hpp"
#include <vector>

namespace sequential {

class SequentialSolver {
public:
  explicit SequentialSolver(const SolverParams &params);

  void solve();

private:
  const SolverParams &params_;

  int N_;
  size_t total_points_;

  std::vector<double> u_prev_;
  std::vector<double> u_curr_;
  std::vector<double> u_next_;
  std::vector<double> u_analytical_;

  void initializeU0();
  void computeU1();
  void computeTimeStep();
  void applyBoundaryConditions(std::vector<double> &u);
};

} // namespace sequential

#endif
