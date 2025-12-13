#ifndef MPI_SOLVER_HPP
#define MPI_SOLVER_HPP

#include "mpi_solver_params.hpp"
#include <mpi.h>
#include <vector>

namespace mpi_solver {

class MPISolver {
public:
  explicit MPISolver(const SolverParams &params);
  ~MPISolver();

  void solve();

private:
  const SolverParams &params_;

  MPI_Comm cart_comm_;
  int rank_;
  int nprocs_;
  int dims_[3];
  int coords_[3];
  int neighbors_[6];

  int Nx_, Ny_, Nz_;
  double hx_, hy_, hz_;
  double invhx2_, invhy2_, invhz2_;
  double coeff_;

  int local_nx_, local_ny_, local_nz_;
  int global_i_start_, global_j_start_, global_k_start_;

  int local_nx_ghost_, local_ny_ghost_, local_nz_ghost_;
  size_t local_total_size_;

  std::vector<double> u_prev_;
  std::vector<double> u_curr_;
  std::vector<double> u_next_;
  std::vector<double> u_analytical_;

  std::vector<double> send_x_left_, send_x_right_;
  std::vector<double> recv_x_left_, recv_x_right_;
  std::vector<double> send_y_down_, send_y_up_;
  std::vector<double> recv_y_down_, recv_y_up_;
  std::vector<double> send_z_back_, send_z_front_;
  std::vector<double> recv_z_back_, recv_z_front_;

  void setupTopology();
  void computeLocalBounds();
  void allocateBuffers();

  void initializeU0();
  void computeU1();
  void computeTimeStep();
  void applyBoundaryConditions(std::vector<double> &u);
  double computeLaplacian(const std::vector<double> &u, int i, int j,
                          int k) const;

  void exchangeGhosts(std::vector<double> &u);

  void packSendBuffers(const std::vector<double> &u);
  void startAsyncExchange(MPI_Request *reqs, int &req_count);
  void unpackRecvBuffers(std::vector<double> &u);
  void computeInterior(std::vector<double> &u_next,
                       const std::vector<double> &u_curr,
                       const std::vector<double> &u_prev, double coeff);
  void computeBoundaryShell(std::vector<double> &u_next,
                            const std::vector<double> &u_curr,
                            const std::vector<double> &u_prev, double coeff);

  inline size_t globalToLocalIdx(int i, int j, int k) const {
    return k + local_nz_ghost_ * (j + local_ny_ghost_ * i);
  }

  inline int periodicX(int i) const { return (i % Nx_ + Nx_) % Nx_; }
  inline int periodicZ(int k) const { return (k % Nz_ + Nz_) % Nz_; }

  void computeLocalAnalyticalSolution(std::vector<double> &u_analytical,
                                      double t);
  double computeLocalMaxError(const std::vector<double> &u_numerical,
                              const std::vector<double> &u_analytical);
  double computeLocalSumSquaredError(const std::vector<double> &u_numerical,
                                     const std::vector<double> &u_analytical);
};

} // namespace mpi_solver

#endif
