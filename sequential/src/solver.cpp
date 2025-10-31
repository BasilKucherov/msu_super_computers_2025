#include "solver.hpp"
#include "error_metrics.hpp"
#include "grid_utils.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>

namespace sequential {

SequentialSolver::SequentialSolver(const SolverParams &params)
    : params_(params) {
  N_ = params_.N;
  Nx_ = N_;
  Ny_ = N_ + 1;
  Nz_ = N_;
  hx_ = params_.Lx / Nx_;
  hy_ = params_.Ly / (Ny_ - 1);
  hz_ = params_.Lz / Nz_;
  invhx2_ = 1.0 / (hx_ * hx_);
  invhy2_ = 1.0 / (hy_ * hy_);
  invhz2_ = 1.0 / (hz_ * hz_);

  coeff_ = params_.tau * params_.tau * params_.a_squared;
  total_points_ = Nx_ * Ny_ * Nz_;

  u_prev_.resize(total_points_);
  u_curr_.resize(total_points_);
  u_next_.resize(total_points_);
  u_analytical_.resize(total_points_);
}

void SequentialSolver::initializeU0() {
  utils::computeAnalyticalSolution(u_prev_, params_.u_analytical, params_.Lx,
                                   params_.Ly, params_.Lz, Nx_, Ny_, Nz_, 0.0);
  applyBoundaryConditions(u_prev_);
}

void SequentialSolver::computeU1() {
  const double coeff = 0.5 * coeff_;

  for (int i = 0; i < Nx_; ++i) {
    for (int j = 1; j < Ny_ - 1; ++j) {
      for (int k = 0; k < Nz_; ++k) {
        const size_t idx = utils::idx3d(i, j, k, Ny_, Nz_);
        const double laplacian = computeLaplacian(u_prev_, i, j, k);
        u_curr_[idx] = u_prev_[idx] + coeff * laplacian;
      }
    }
  }

  applyBoundaryConditions(u_curr_);
}

void SequentialSolver::computeTimeStep() {
  for (int i = 0; i < Nx_; ++i) {
    for (int j = 1; j < Ny_ - 1; ++j) {
      for (int k = 0; k < Nz_; ++k) {
        const size_t idx = utils::idx3d(i, j, k, Ny_, Nz_);
        const double laplacian = computeLaplacian(u_curr_, i, j, k);
        u_next_[idx] = 2.0 * u_curr_[idx] - u_prev_[idx] + coeff_ * laplacian;
      }
    }
  }

  applyBoundaryConditions(u_next_);
  u_prev_.swap(u_curr_);
  u_curr_.swap(u_next_);
}

void SequentialSolver::applyBoundaryConditions(std::vector<double> &u) {
  for (int i = 0; i < Nx_; ++i) {
    for (int k = 0; k < Nz_; ++k) {
      u[utils::idx3d(i, 0, k, Ny_, Nz_)] = 0.0;
      u[utils::idx3d(i, Ny_ - 1, k, Ny_, Nz_)] = 0.0;
    }
  }
}

double SequentialSolver::computeLaplacian(const std::vector<double> &u, int i,
                                          int j, int k) const {
  const size_t idx = utils::idx3d(i, j, k, Ny_, Nz_);

  const double d2u_dx2 =
      (u[utils::idx3d(utils::ip(i, Nx_), j, k, Ny_, Nz_)] - 2.0 * u[idx] +
       u[utils::idx3d(utils::im(i, Nx_), j, k, Ny_, Nz_)]) *
      invhx2_;

  const double d2u_dy2 =
      (u[utils::idx3d(i, j + 1, k, Ny_, Nz_)] - 2.0 * u[idx] +
       u[utils::idx3d(i, j - 1, k, Ny_, Nz_)]) *
      invhy2_;

  const double d2u_dz2 =
      (u[utils::idx3d(i, j, utils::ip(k, Nz_), Ny_, Nz_)] - 2.0 * u[idx] +
       u[utils::idx3d(i, j, utils::im(k, Nz_), Ny_, Nz_)]) *
      invhz2_;

  return d2u_dx2 + d2u_dy2 + d2u_dz2;
}

void SequentialSolver::solve() {
  std::cout << "Starting Sequential Wave Equation Solver" << std::endl;
  std::cout << "Mode: " << (params_.debug ? "DEBUG" : "TEST") << std::endl;
  std::cout << "Grid size: " << Nx_ << "*" << Ny_ << "*" << Nz_ << " = "
            << total_points_ << " points" << std::endl;
  std::cout << "Time steps: " << params_.K << ", τ = " << params_.tau
            << ", T = " << params_.T << std::endl;
  std::cout << "Domain: [0, " << params_.Lx << "] * [0, " << params_.Ly
            << "] * [0, " << params_.Lz << "]" << std::endl;
  std::cout << std::endl;

  std::cout << std::fixed << std::setprecision(24);

  if (params_.debug) {
    initializeU0();
    utils::computeAnalyticalSolution(u_analytical_, params_.u_analytical,
                                     params_.Lx, params_.Ly, params_.Lz, Nx_,
                                     Ny_, Nz_, 0.0);
    auto metrics = utils::computeErrorMetrics(u_prev_, u_analytical_);
    std::cout << "Step 0, t=0.000000: max_error=" << metrics.max_abs_error
              << ", rmse=" << metrics.rmse << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    computeU1();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    utils::computeAnalyticalSolution(u_analytical_, params_.u_analytical,
                                     params_.Lx, params_.Ly, params_.Lz, Nx_,
                                     Ny_, Nz_, params_.tau);
    metrics = utils::computeErrorMetrics(u_curr_, u_analytical_);
    std::cout << "Step 1, t=" << params_.tau
              << ": max_error=" << metrics.max_abs_error
              << ", rmse=" << metrics.rmse << ", time=" << elapsed.count()
              << "s" << std::endl;

    double total_time = elapsed.count();

    for (int n = 2; n <= params_.K; ++n) {
      double t_current = n * params_.tau;

      start = std::chrono::high_resolution_clock::now();
      computeTimeStep();
      end = std::chrono::high_resolution_clock::now();
      elapsed = end - start;
      total_time += elapsed.count();

      utils::computeAnalyticalSolution(u_analytical_, params_.u_analytical,
                                       params_.Lx, params_.Ly, params_.Lz, Nx_,
                                       Ny_, Nz_, t_current);
      metrics = utils::computeErrorMetrics(u_curr_, u_analytical_);
      std::cout << "Step " << n << ", t=" << t_current
                << ": max_error=" << metrics.max_abs_error
                << ", rmse=" << metrics.rmse << ", time=" << elapsed.count()
                << "s" << std::endl;
    }

    std::cout << std::endl;
    std::cout << "Total computation time: " << total_time << "s" << std::endl;
    std::cout << "Final errors - max_error: " << metrics.max_abs_error
              << ", rmse: " << metrics.rmse << std::endl;
  } else {
    auto start = std::chrono::high_resolution_clock::now();
    initializeU0();
    computeU1();
    for (int n = 2; n <= params_.K; ++n) {
      computeTimeStep();
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_time = end - start;

    double t_final = params_.K * params_.tau;
    utils::computeAnalyticalSolution(u_analytical_, params_.u_analytical,
                                     params_.Lx, params_.Ly, params_.Lz, Nx_,
                                     Ny_, Nz_, t_final);
    auto metrics = utils::computeErrorMetrics(u_curr_, u_analytical_);

    std::cout << std::endl;
    std::cout << "Total computation time: " << total_time.count() << "s"
              << std::endl;
    std::cout << "Time per step: " << total_time.count() / params_.K << "s"
              << std::endl;
    std::cout << "Final errors at t=" << t_final
              << " - max_error: " << metrics.max_abs_error
              << ", rmse: " << metrics.rmse << std::endl;
  }
}

} // namespace sequential
