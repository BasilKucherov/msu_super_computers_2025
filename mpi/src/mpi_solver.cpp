#include "mpi_solver.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace mpi_solver {

MPISolver::MPISolver(const SolverParams &params) : params_(params) {
  Nx_ = params_.N;
  Ny_ = params_.N + 1;
  Nz_ = params_.N;

  hx_ = params_.Lx / Nx_;
  hy_ = params_.Ly / (Ny_ - 1);
  hz_ = params_.Lz / Nz_;

  invhx2_ = 1.0 / (hx_ * hx_);
  invhy2_ = 1.0 / (hy_ * hy_);
  invhz2_ = 1.0 / (hz_ * hz_);

  coeff_ = params_.tau * params_.tau * params_.a_squared;

  setupTopology();
  computeLocalBounds();
  allocateBuffers();
}

MPISolver::~MPISolver() {
  if (cart_comm_ != MPI_COMM_NULL && cart_comm_ != MPI_COMM_WORLD) {
    MPI_Comm_free(&cart_comm_);
  }
}

void MPISolver::setupTopology() {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs_);

  std::fill(std::begin(dims_), std::end(dims_), 0);
  MPI_Dims_create(nprocs_, 3, dims_);

  int periods[3] = {1, 0, 1};
  int reorder = 1;

  MPI_Cart_create(MPI_COMM_WORLD, 3, dims_, periods, reorder, &cart_comm_);
  MPI_Cart_coords(cart_comm_, rank_, 3, coords_);

  // left, right
  MPI_Cart_shift(cart_comm_, 0, 1, &neighbors_[0], &neighbors_[1]);
  // down, up
  MPI_Cart_shift(cart_comm_, 1, 1, &neighbors_[2], &neighbors_[3]);
  // back, front
  MPI_Cart_shift(cart_comm_, 2, 1, &neighbors_[4], &neighbors_[5]);
}

void MPISolver::computeLocalBounds() {
  int nx_per_process_ = Nx_ / dims_[0];
  int nx_leftover = Nx_ % dims_[0];
  local_nx_ = nx_per_process_ + (coords_[0] < nx_leftover ? 1 : 0);
  global_i_start_ = std::min(nx_leftover, coords_[0]) * (nx_per_process_ + 1) +
                    std::max(coords_[0] - nx_leftover, 0) * nx_per_process_;

  int ny_per_process_ = Ny_ / dims_[1];
  int ny_leftover = Ny_ % dims_[1];
  local_ny_ = ny_per_process_ + (coords_[1] < ny_leftover ? 1 : 0);
  global_j_start_ = std::min(ny_leftover, coords_[1]) * (ny_per_process_ + 1) +
                    std::max(coords_[1] - ny_leftover, 0) * ny_per_process_;

  int nz_per_process_ = Nz_ / dims_[2];
  int nz_leftover = Nz_ % dims_[2];
  local_nz_ = nz_per_process_ + (coords_[2] < nz_leftover ? 1 : 0);
  global_k_start_ = std::min(nz_leftover, coords_[2]) * (nz_per_process_ + 1) +
                    std::max(coords_[2] - nz_leftover, 0) * nz_per_process_;

  local_nx_ghost_ = local_nx_ + 2;
  local_ny_ghost_ = local_ny_ + 2;
  local_nz_ghost_ = local_nz_ + 2;
  local_total_size_ = local_nx_ghost_ * local_ny_ghost_ * local_nz_ghost_;
}

void MPISolver::allocateBuffers() {
  u_prev_.resize(local_total_size_);
  u_curr_.resize(local_total_size_);
  u_next_.resize(local_total_size_);
  u_analytical_.resize(local_total_size_);

  size_t x_face_size = local_ny_ * local_nz_;
  send_x_left_.resize(x_face_size);
  send_x_right_.resize(x_face_size);
  recv_x_left_.resize(x_face_size);
  recv_x_right_.resize(x_face_size);

  size_t y_face_size = local_nx_ * local_nz_;
  send_y_down_.resize(y_face_size);
  send_y_up_.resize(y_face_size);
  recv_y_down_.resize(y_face_size);
  recv_y_up_.resize(y_face_size);

  size_t z_face_size = local_nx_ * local_ny_;
  send_z_back_.resize(z_face_size);
  send_z_front_.resize(z_face_size);
  recv_z_back_.resize(z_face_size);
  recv_z_front_.resize(z_face_size);
}

void MPISolver::exchangeGhosts(std::vector<double> &u) {
  MPI_Request reqs[12];
  int req_count = 0;

  for (int j = 1; j <= local_ny_; ++j) {
    for (int k = 1; k <= local_nz_; ++k) {
      size_t buf_idx = (k - 1) + local_nz_ * (j - 1);
      send_x_left_[buf_idx] = u[globalToLocalIdx(1, j, k)];
      send_x_right_[buf_idx] = u[globalToLocalIdx(local_nx_, j, k)];
    }
  }

  MPI_Isend(send_x_left_.data(), send_x_left_.size(), MPI_DOUBLE, neighbors_[0],
            0, cart_comm_, &reqs[req_count++]);
  MPI_Irecv(recv_x_right_.data(), recv_x_right_.size(), MPI_DOUBLE,
            neighbors_[1], 0, cart_comm_, &reqs[req_count++]);

  MPI_Isend(send_x_right_.data(), send_x_right_.size(), MPI_DOUBLE,
            neighbors_[1], 1, cart_comm_, &reqs[req_count++]);
  MPI_Irecv(recv_x_left_.data(), recv_x_left_.size(), MPI_DOUBLE, neighbors_[0],
            1, cart_comm_, &reqs[req_count++]);

  for (int i = 1; i <= local_nx_; ++i) {
    for (int k = 1; k <= local_nz_; ++k) {
      size_t buf_idx = (k - 1) + local_nz_ * (i - 1);
      send_y_down_[buf_idx] = u[globalToLocalIdx(i, 1, k)];
      send_y_up_[buf_idx] = u[globalToLocalIdx(i, local_ny_, k)];
    }
  }

  MPI_Isend(send_y_down_.data(), send_y_down_.size(), MPI_DOUBLE, neighbors_[2],
            2, cart_comm_, &reqs[req_count++]);
  MPI_Irecv(recv_y_up_.data(), recv_y_up_.size(), MPI_DOUBLE, neighbors_[3], 2,
            cart_comm_, &reqs[req_count++]);

  MPI_Isend(send_y_up_.data(), send_y_up_.size(), MPI_DOUBLE, neighbors_[3], 3,
            cart_comm_, &reqs[req_count++]);
  MPI_Irecv(recv_y_down_.data(), recv_y_down_.size(), MPI_DOUBLE, neighbors_[2],
            3, cart_comm_, &reqs[req_count++]);

  for (int i = 1; i <= local_nx_; ++i) {
    for (int j = 1; j <= local_ny_; ++j) {
      size_t buf_idx = (j - 1) + local_ny_ * (i - 1);
      send_z_back_[buf_idx] = u[globalToLocalIdx(i, j, 1)];
      send_z_front_[buf_idx] = u[globalToLocalIdx(i, j, local_nz_)];
    }
  }

  MPI_Isend(send_z_back_.data(), send_z_back_.size(), MPI_DOUBLE, neighbors_[4],
            4, cart_comm_, &reqs[req_count++]);
  MPI_Irecv(recv_z_front_.data(), recv_z_front_.size(), MPI_DOUBLE,
            neighbors_[5], 4, cart_comm_, &reqs[req_count++]);

  MPI_Isend(send_z_front_.data(), send_z_front_.size(), MPI_DOUBLE,
            neighbors_[5], 5, cart_comm_, &reqs[req_count++]);
  MPI_Irecv(recv_z_back_.data(), recv_z_back_.size(), MPI_DOUBLE, neighbors_[4],
            5, cart_comm_, &reqs[req_count++]);

  MPI_Waitall(req_count, reqs, MPI_STATUSES_IGNORE);

  for (int j = 1; j <= local_ny_; ++j) {
    for (int k = 1; k <= local_nz_; ++k) {
      size_t buf_idx = (k - 1) + local_nz_ * (j - 1);
      u[globalToLocalIdx(0, j, k)] = recv_x_left_[buf_idx];
      u[globalToLocalIdx(local_nx_ + 1, j, k)] = recv_x_right_[buf_idx];
    }
  }

  if (neighbors_[2] != MPI_PROC_NULL) {
    for (int i = 1; i <= local_nx_; ++i) {
      for (int k = 1; k <= local_nz_; ++k) {
        size_t buf_idx = (k - 1) + local_nz_ * (i - 1);
        u[globalToLocalIdx(i, 0, k)] = recv_y_down_[buf_idx];
      }
    }
  }

  if (neighbors_[3] != MPI_PROC_NULL) {
    for (int i = 1; i <= local_nx_; ++i) {
      for (int k = 1; k <= local_nz_; ++k) {
        size_t buf_idx = (k - 1) + local_nz_ * (i - 1);
        u[globalToLocalIdx(i, local_ny_ + 1, k)] = recv_y_up_[buf_idx];
      }
    }
  }

  for (int i = 1; i <= local_nx_; ++i) {
    for (int j = 1; j <= local_ny_; ++j) {
      size_t buf_idx = (j - 1) + local_ny_ * (i - 1);
      u[globalToLocalIdx(i, j, 0)] = recv_z_back_[buf_idx];
      u[globalToLocalIdx(i, j, local_nz_ + 1)] = recv_z_front_[buf_idx];
    }
  }
}

void MPISolver::computeLocalAnalyticalSolution(
    std::vector<double> &u_analytical, double t) {
  for (int i = 1; i <= local_nx_; ++i) {
    int global_i = global_i_start_ + (i - 1);
    double x = global_i * hx_;

    for (int j = 1; j <= local_ny_; ++j) {
      int global_j = global_j_start_ + (j - 1);
      double y = global_j * hy_;

      for (int k = 1; k <= local_nz_; ++k) {
        int global_k = global_k_start_ + (k - 1);
        double z = global_k * hz_;

        u_analytical[globalToLocalIdx(i, j, k)] =
            params_.u_analytical(x, y, z, t);
      }
    }
  }
}

void MPISolver::initializeU0() {
  computeLocalAnalyticalSolution(u_prev_, 0.0);
  exchangeGhosts(u_prev_);
  applyBoundaryConditions(u_prev_);
}

void MPISolver::computeU1() {
  const double coeff = 0.5 * coeff_;

  for (int i = 1; i <= local_nx_; ++i) {
    for (int j = 1; j <= local_ny_; ++j) {
      int global_j = global_j_start_ + (j - 1);
      if (global_j == 0 || global_j == Ny_ - 1)
        continue;

      for (int k = 1; k <= local_nz_; ++k) {
        const double laplacian = computeLaplacian(u_prev_, i, j, k);
        u_curr_[globalToLocalIdx(i, j, k)] =
            u_prev_[globalToLocalIdx(i, j, k)] + coeff * laplacian;
      }
    }
  }

  exchangeGhosts(u_curr_);
  applyBoundaryConditions(u_curr_);
}

void MPISolver::computeTimeStep() {
  for (int i = 1; i <= local_nx_; ++i) {
    for (int j = 1; j <= local_ny_; ++j) {
      int global_j = global_j_start_ + (j - 1);
      if (global_j == 0 || global_j == Ny_ - 1)
        continue;

      for (int k = 1; k <= local_nz_; ++k) {
        size_t idx = globalToLocalIdx(i, j, k);
        const double laplacian = computeLaplacian(u_curr_, i, j, k);
        u_next_[idx] = 2.0 * u_curr_[idx] - u_prev_[idx] + coeff_ * laplacian;
      }
    }
  }

  exchangeGhosts(u_next_);
  applyBoundaryConditions(u_next_);

  u_prev_.swap(u_curr_);
  u_curr_.swap(u_next_);
}

void MPISolver::applyBoundaryConditions(std::vector<double> &u) {
  if (global_j_start_ == 0) {
    int local_j_boundary = 1;
    for (int i = 1; i <= local_nx_; ++i) {
      for (int k = 1; k <= local_nz_; ++k) {
        u[globalToLocalIdx(i, local_j_boundary, k)] = 0.0;
      }
    }
  }

  int global_j_end = global_j_start_ + local_ny_ - 1;
  if (global_j_end == Ny_ - 1) {
    int local_j_boundary = local_ny_;
    for (int i = 1; i <= local_nx_; ++i) {
      for (int k = 1; k <= local_nz_; ++k) {
        u[globalToLocalIdx(i, local_j_boundary, k)] = 0.0;
      }
    }
  }
}

double MPISolver::computeLaplacian(const std::vector<double> &u, int i, int j,
                                   int k) const {
  const size_t idx = globalToLocalIdx(i, j, k);

  const double d2u_dx2 = (u[globalToLocalIdx(i + 1, j, k)] - 2.0 * u[idx] +
                          u[globalToLocalIdx(i - 1, j, k)]) *
                         invhx2_;

  const double d2u_dy2 = (u[globalToLocalIdx(i, j + 1, k)] - 2.0 * u[idx] +
                          u[globalToLocalIdx(i, j - 1, k)]) *
                         invhy2_;

  const double d2u_dz2 = (u[globalToLocalIdx(i, j, k + 1)] - 2.0 * u[idx] +
                          u[globalToLocalIdx(i, j, k - 1)]) *
                         invhz2_;

  return d2u_dx2 + d2u_dy2 + d2u_dz2;
}

double
MPISolver::computeLocalMaxError(const std::vector<double> &u_numerical,
                                const std::vector<double> &u_analytical) {
  double local_max = 0.0;
  for (int i = 1; i <= local_nx_; ++i) {
    for (int j = 1; j <= local_ny_; ++j) {
      for (int k = 1; k <= local_nz_; ++k) {
        size_t idx = globalToLocalIdx(i, j, k);
        double error = std::abs(u_numerical[idx] - u_analytical[idx]);
        if (error > local_max)
          local_max = error;
      }
    }
  }
  return local_max;
}

double MPISolver::computeLocalSumSquaredError(
    const std::vector<double> &u_numerical,
    const std::vector<double> &u_analytical) {
  double sum_sq = 0.0;
  for (int i = 1; i <= local_nx_; ++i) {
    for (int j = 1; j <= local_ny_; ++j) {
      for (int k = 1; k <= local_nz_; ++k) {
        size_t idx = globalToLocalIdx(i, j, k);
        double error = u_numerical[idx] - u_analytical[idx];
        sum_sq += error * error;
      }
    }
  }
  return sum_sq;
}

void MPISolver::solve() {
  size_t total_global_points = Nx_ * Ny_ * Nz_;

  if (rank_ == 0) {
    std::cout << "Starting MPI Wave Equation Solver" << std::endl;
    std::cout << "Mode: " << (params_.debug ? "DEBUG" : "TEST") << std::endl;
    std::cout << "Processes: " << nprocs_ << " (" << dims_[0] << "x" << dims_[1]
              << "x" << dims_[2] << ")" << std::endl;
    std::cout << "Global grid size: " << Nx_ << "*" << Ny_ << "*" << Nz_
              << " = " << total_global_points << " points" << std::endl;
    std::cout << "Time steps: " << params_.K << ", τ = " << params_.tau
              << ", T = " << params_.T << std::endl;
    std::cout << "Domain: [0, " << params_.Lx << "] * [0, " << params_.Ly
              << "] * [0, " << params_.Lz << "]" << std::endl;
    std::cout << std::endl;
    std::cout << std::scientific << std::setprecision(16);
  }

  if (params_.debug) {
    initializeU0();
    computeLocalAnalyticalSolution(u_analytical_, 0.0);

    double local_max = computeLocalMaxError(u_prev_, u_analytical_);
    double local_sum_sq = computeLocalSumSquaredError(u_prev_, u_analytical_);

    double global_max, global_sum_sq;
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, cart_comm_);
    MPI_Allreduce(&local_sum_sq, &global_sum_sq, 1, MPI_DOUBLE, MPI_SUM,
                  cart_comm_);
    double global_rmse = std::sqrt(global_sum_sq / total_global_points);

    if (rank_ == 0) {
      std::cout << "Step 0, t=0.000000: max_error=" << global_max
                << ", rmse=" << global_rmse << std::endl;
    }

    auto start = std::chrono::high_resolution_clock::now();
    computeU1();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    computeLocalAnalyticalSolution(u_analytical_, params_.tau);
    local_max = computeLocalMaxError(u_curr_, u_analytical_);
    local_sum_sq = computeLocalSumSquaredError(u_curr_, u_analytical_);

    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, cart_comm_);
    MPI_Allreduce(&local_sum_sq, &global_sum_sq, 1, MPI_DOUBLE, MPI_SUM,
                  cart_comm_);
    global_rmse = std::sqrt(global_sum_sq / total_global_points);

    if (rank_ == 0) {
      std::cout << "Step 1, t=" << params_.tau << ": max_error=" << global_max
                << ", rmse=" << global_rmse << ", time=" << elapsed.count()
                << "s" << std::endl;
    }

    double total_time = elapsed.count();

    for (int n = 2; n <= params_.K; ++n) {
      double t_current = n * params_.tau;

      start = std::chrono::high_resolution_clock::now();
      computeTimeStep();
      end = std::chrono::high_resolution_clock::now();
      elapsed = end - start;
      total_time += elapsed.count();

      computeLocalAnalyticalSolution(u_analytical_, t_current);
      local_max = computeLocalMaxError(u_curr_, u_analytical_);
      local_sum_sq = computeLocalSumSquaredError(u_curr_, u_analytical_);

      MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX,
                    cart_comm_);
      MPI_Allreduce(&local_sum_sq, &global_sum_sq, 1, MPI_DOUBLE, MPI_SUM,
                    cart_comm_);
      global_rmse = std::sqrt(global_sum_sq / total_global_points);

      if (rank_ == 0) {
        std::cout << "Step " << n << ", t=" << t_current
                  << ": max_error=" << global_max << ", rmse=" << global_rmse
                  << ", time=" << elapsed.count() << "s" << std::endl;
      }
    }

    if (rank_ == 0) {
      std::cout << std::endl;
      std::cout << "Total computation time: " << total_time << "s" << std::endl;
      std::cout << "Final errors - max_error: " << global_max
                << ", rmse: " << global_rmse << std::endl;
    }
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
    computeLocalAnalyticalSolution(u_analytical_, t_final);

    double local_max = computeLocalMaxError(u_curr_, u_analytical_);
    double local_sum_sq = computeLocalSumSquaredError(u_curr_, u_analytical_);

    double global_max, global_sum_sq;
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, cart_comm_);
    MPI_Allreduce(&local_sum_sq, &global_sum_sq, 1, MPI_DOUBLE, MPI_SUM,
                  cart_comm_);
    double global_rmse = std::sqrt(global_sum_sq / total_global_points);

    if (rank_ == 0) {
      std::cout << std::endl;
      std::cout << "Total computation time: " << total_time.count() << "s"
                << std::endl;
      std::cout << "Time per step: " << total_time.count() / params_.K << "s"
                << std::endl;
      std::cout << "Final errors at t=" << t_final
                << " - max_error: " << global_max << ", rmse: " << global_rmse
                << std::endl;
    }
  }
}

} // namespace mpi_solver
