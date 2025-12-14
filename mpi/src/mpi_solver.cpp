#include "mpi_solver.hpp"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

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

  for (int i = 1; i <= local_nx_; ++i) {
    for (int k = 1; k <= local_nz_; ++k) {
      size_t buf_idx = (k - 1) + local_nz_ * (i - 1);
      send_y_down_[buf_idx] = u[globalToLocalIdx(i, 1, k)];
      send_y_up_[buf_idx] = u[globalToLocalIdx(i, local_ny_, k)];
    }
  }

  for (int i = 1; i <= local_nx_; ++i) {
    for (int j = 1; j <= local_ny_; ++j) {
      size_t buf_idx = (j - 1) + local_ny_ * (i - 1);
      send_z_back_[buf_idx] = u[globalToLocalIdx(i, j, 1)];
      send_z_front_[buf_idx] = u[globalToLocalIdx(i, j, local_nz_)];
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

  MPI_Isend(send_y_down_.data(), send_y_down_.size(), MPI_DOUBLE, neighbors_[2],
            2, cart_comm_, &reqs[req_count++]);
  MPI_Irecv(recv_y_up_.data(), recv_y_up_.size(), MPI_DOUBLE, neighbors_[3], 2,
            cart_comm_, &reqs[req_count++]);

  MPI_Isend(send_y_up_.data(), send_y_up_.size(), MPI_DOUBLE, neighbors_[3], 3,
            cart_comm_, &reqs[req_count++]);
  MPI_Irecv(recv_y_down_.data(), recv_y_down_.size(), MPI_DOUBLE, neighbors_[2],
            3, cart_comm_, &reqs[req_count++]);

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

void MPISolver::packSendBuffers(const std::vector<double> &u) {
  for (int j = 1; j <= local_ny_; ++j) {
    for (int k = 1; k <= local_nz_; ++k) {
      size_t buf_idx = (k - 1) + local_nz_ * (j - 1);
      send_x_left_[buf_idx] = u[globalToLocalIdx(1, j, k)];
      send_x_right_[buf_idx] = u[globalToLocalIdx(local_nx_, j, k)];
    }
  }

  for (int i = 1; i <= local_nx_; ++i) {
    for (int k = 1; k <= local_nz_; ++k) {
      size_t buf_idx = (k - 1) + local_nz_ * (i - 1);
      send_y_down_[buf_idx] = u[globalToLocalIdx(i, 1, k)];
      send_y_up_[buf_idx] = u[globalToLocalIdx(i, local_ny_, k)];
    }
  }

  for (int i = 1; i <= local_nx_; ++i) {
    for (int j = 1; j <= local_ny_; ++j) {
      size_t buf_idx = (j - 1) + local_ny_ * (i - 1);
      send_z_back_[buf_idx] = u[globalToLocalIdx(i, j, 1)];
      send_z_front_[buf_idx] = u[globalToLocalIdx(i, j, local_nz_)];
    }
  }
}

void MPISolver::startAsyncExchange(MPI_Request *reqs, int &req_count) {
  req_count = 0;

  MPI_Isend(send_x_left_.data(), send_x_left_.size(), MPI_DOUBLE, neighbors_[0],
            0, cart_comm_, &reqs[req_count++]);
  MPI_Irecv(recv_x_right_.data(), recv_x_right_.size(), MPI_DOUBLE,
            neighbors_[1], 0, cart_comm_, &reqs[req_count++]);

  MPI_Isend(send_x_right_.data(), send_x_right_.size(), MPI_DOUBLE,
            neighbors_[1], 1, cart_comm_, &reqs[req_count++]);
  MPI_Irecv(recv_x_left_.data(), recv_x_left_.size(), MPI_DOUBLE, neighbors_[0],
            1, cart_comm_, &reqs[req_count++]);

  MPI_Isend(send_y_down_.data(), send_y_down_.size(), MPI_DOUBLE, neighbors_[2],
            2, cart_comm_, &reqs[req_count++]);
  MPI_Irecv(recv_y_up_.data(), recv_y_up_.size(), MPI_DOUBLE, neighbors_[3], 2,
            cart_comm_, &reqs[req_count++]);

  MPI_Isend(send_y_up_.data(), send_y_up_.size(), MPI_DOUBLE, neighbors_[3], 3,
            cart_comm_, &reqs[req_count++]);
  MPI_Irecv(recv_y_down_.data(), recv_y_down_.size(), MPI_DOUBLE, neighbors_[2],
            3, cart_comm_, &reqs[req_count++]);

  MPI_Isend(send_z_back_.data(), send_z_back_.size(), MPI_DOUBLE, neighbors_[4],
            4, cart_comm_, &reqs[req_count++]);
  MPI_Irecv(recv_z_front_.data(), recv_z_front_.size(), MPI_DOUBLE,
            neighbors_[5], 4, cart_comm_, &reqs[req_count++]);

  MPI_Isend(send_z_front_.data(), send_z_front_.size(), MPI_DOUBLE,
            neighbors_[5], 5, cart_comm_, &reqs[req_count++]);
  MPI_Irecv(recv_z_back_.data(), recv_z_back_.size(), MPI_DOUBLE, neighbors_[4],
            5, cart_comm_, &reqs[req_count++]);
}

void MPISolver::unpackRecvBuffers(std::vector<double> &u) {
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

void MPISolver::computeInterior(std::vector<double> &u_next,
                                const std::vector<double> &u_curr,
                                const std::vector<double> &u_prev,
                                double coeff) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 2; i <= local_nx_ - 1; ++i) {
    for (int j = 2; j <= local_ny_ - 1; ++j) {
      int global_j = global_j_start_ + (j - 1);
      if (global_j == 0 || global_j == Ny_ - 1)
        continue;

      for (int k = 2; k <= local_nz_ - 1; ++k) {
        size_t idx = globalToLocalIdx(i, j, k);
        const double laplacian = computeLaplacian(u_curr, i, j, k);
        u_next[idx] = 2.0 * u_curr[idx] - u_prev[idx] + coeff * laplacian;
      }
    }
  }
}

void MPISolver::computeBoundaryShell(std::vector<double> &u_next,
                                     const std::vector<double> &u_curr,
                                     const std::vector<double> &u_prev,
                                     double coeff) {
  auto compute_point = [&](int i, int j, int k) {
    int global_j = global_j_start_ + (j - 1);
    if (global_j == 0 || global_j == Ny_ - 1)
      return;
    size_t idx = globalToLocalIdx(i, j, k);
    const double laplacian = computeLaplacian(u_curr, i, j, k);
    u_next[idx] = 2.0 * u_curr[idx] - u_prev[idx] + coeff * laplacian;
  };

  for (int j = 1; j <= local_ny_; ++j) {
    for (int k = 1; k <= local_nz_; ++k) {
      compute_point(1, j, k);
      if (local_nx_ > 1) {
        compute_point(local_nx_, j, k);
      }
    }
  }

  for (int i = 2; i <= local_nx_ - 1; ++i) {
    for (int k = 1; k <= local_nz_; ++k) {
      compute_point(i, 1, k);
      if (local_ny_ > 1) {
        compute_point(i, local_ny_, k);
      }
    }
  }

  for (int i = 2; i <= local_nx_ - 1; ++i) {
    for (int j = 2; j <= local_ny_ - 1; ++j) {
      compute_point(i, j, 1);
      if (local_nz_ > 1) {
        compute_point(i, j, local_nz_);
      }
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

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
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

  applyBoundaryConditions(u_curr_);
}

void MPISolver::computeTimeStep() {
  auto t0 = std::chrono::high_resolution_clock::now();

  MPI_Request reqs[12];
  int req_count;

  auto t1 = std::chrono::high_resolution_clock::now();
  packSendBuffers(u_curr_);
  auto t2 = std::chrono::high_resolution_clock::now();

  startAsyncExchange(reqs, req_count);
  auto t3 = std::chrono::high_resolution_clock::now();

  computeInterior(u_next_, u_curr_, u_prev_, coeff_);
  auto t4 = std::chrono::high_resolution_clock::now();

  MPI_Waitall(req_count, reqs, MPI_STATUSES_IGNORE);
  auto t5 = std::chrono::high_resolution_clock::now();

  unpackRecvBuffers(u_curr_);
  auto t6 = std::chrono::high_resolution_clock::now();

  computeBoundaryShell(u_next_, u_curr_, u_prev_, coeff_);
  auto t7 = std::chrono::high_resolution_clock::now();

  applyBoundaryConditions(u_next_);
  auto t8 = std::chrono::high_resolution_clock::now();

  u_prev_.swap(u_curr_);
  u_curr_.swap(u_next_);

  timings_.pack_buffers += std::chrono::duration<double>(t2 - t1).count();
  timings_.start_async += std::chrono::duration<double>(t3 - t2).count();
  timings_.compute_interior += std::chrono::duration<double>(t4 - t3).count();
  timings_.wait_comm += std::chrono::duration<double>(t5 - t4).count();
  timings_.unpack_buffers += std::chrono::duration<double>(t6 - t5).count();
  timings_.compute_boundary += std::chrono::duration<double>(t7 - t6).count();
  timings_.apply_bc += std::chrono::duration<double>(t8 - t7).count();
  timings_.total_step += std::chrono::duration<double>(t8 - t0).count();
  timings_.num_steps++;
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

    printProfilingResults();
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

    printProfilingResults();
  }
}

void MPISolver::printProfilingResults() {
  struct TimingData {
    double pack_buffers;
    double start_async;
    double compute_interior;
    double wait_comm;
    double unpack_buffers;
    double compute_boundary;
    double apply_bc;
    double total_step;
    int num_steps;
    int rank;
    int local_nx;
    int local_ny;
    int local_nz;
  };

  TimingData local_data;
  local_data.pack_buffers = timings_.pack_buffers;
  local_data.start_async = timings_.start_async;
  local_data.compute_interior = timings_.compute_interior;
  local_data.wait_comm = timings_.wait_comm;
  local_data.unpack_buffers = timings_.unpack_buffers;
  local_data.compute_boundary = timings_.compute_boundary;
  local_data.apply_bc = timings_.apply_bc;
  local_data.total_step = timings_.total_step;
  local_data.num_steps = timings_.num_steps;
  local_data.rank = rank_;
  local_data.local_nx = local_nx_;
  local_data.local_ny = local_ny_;
  local_data.local_nz = local_nz_;

  std::vector<TimingData> all_data(nprocs_);

  MPI_Gather(&local_data, sizeof(TimingData), MPI_BYTE, all_data.data(),
             sizeof(TimingData), MPI_BYTE, 0, cart_comm_);

  if (rank_ == 0) {
    std::cout << std::endl;
    std::cout << "============================================================="
                 "===================\n";
    std::cout << "PROFILING RESULTS (per-process breakdown)\n";
    std::cout << "============================================================="
                 "===================\n";
    std::cout << std::fixed << std::setprecision(4);

    std::cout << std::setw(4) << "Rank" << std::setw(12) << "LocalGrid"
              << std::setw(10) << "PackBuf" << std::setw(10) << "StartAsync"
              << std::setw(10) << "Interior" << std::setw(10) << "WaitComm"
              << std::setw(10) << "Unpack" << std::setw(10) << "Boundary"
              << std::setw(10) << "ApplyBC" << std::setw(10) << "Total"
              << std::endl;
    std::cout << std::string(106, '-') << std::endl;

    double sum_pack = 0, sum_async = 0, sum_interior = 0, sum_wait = 0;
    double sum_unpack = 0, sum_boundary = 0, sum_bc = 0, sum_total = 0;
    double max_wait = 0, min_wait = 1e9;
    double max_interior = 0, min_interior = 1e9;
    int slowest_rank = -1, fastest_rank = -1;

    for (int i = 0; i < nprocs_; ++i) {
      const TimingData &d = all_data[i];

      char grid_str[32];
      snprintf(grid_str, sizeof(grid_str), "%dx%dx%d", d.local_nx, d.local_ny,
               d.local_nz);

      std::cout << std::setw(4) << d.rank << std::setw(12) << grid_str
                << std::setw(10) << d.pack_buffers << std::setw(10)
                << d.start_async << std::setw(10) << d.compute_interior
                << std::setw(10) << d.wait_comm << std::setw(10)
                << d.unpack_buffers << std::setw(10) << d.compute_boundary
                << std::setw(10) << d.apply_bc << std::setw(10) << d.total_step
                << std::endl;

      sum_pack += d.pack_buffers;
      sum_async += d.start_async;
      sum_interior += d.compute_interior;
      sum_wait += d.wait_comm;
      sum_unpack += d.unpack_buffers;
      sum_boundary += d.compute_boundary;
      sum_bc += d.apply_bc;
      sum_total += d.total_step;

      if (d.wait_comm > max_wait) {
        max_wait = d.wait_comm;
        slowest_rank = d.rank;
      }
      if (d.wait_comm < min_wait) {
        min_wait = d.wait_comm;
        fastest_rank = d.rank;
      }
      if (d.compute_interior > max_interior)
        max_interior = d.compute_interior;
      if (d.compute_interior < min_interior)
        min_interior = d.compute_interior;
    }

    std::cout << std::string(106, '-') << std::endl;

    std::cout << std::setw(4) << "AVG" << std::setw(12) << "" << std::setw(10)
              << sum_pack / nprocs_ << std::setw(10) << sum_async / nprocs_
              << std::setw(10) << sum_interior / nprocs_ << std::setw(10)
              << sum_wait / nprocs_ << std::setw(10) << sum_unpack / nprocs_
              << std::setw(10) << sum_boundary / nprocs_ << std::setw(10)
              << sum_bc / nprocs_ << std::setw(10) << sum_total / nprocs_
              << std::endl;

    std::cout << std::endl;
    std::cout << "SUMMARY:\n";
    std::cout << "  Interior compute: min=" << min_interior
              << "s, max=" << max_interior
              << "s, imbalance=" << (max_interior / min_interior - 1.0) * 100
              << "%\n";
    std::cout << "  Wait for comm:    min=" << min_wait << "s (rank "
              << fastest_rank << "), max=" << max_wait << "s (rank "
              << slowest_rank << ")\n";

    double avg_interior = sum_interior / nprocs_;
    double avg_wait = sum_wait / nprocs_;
    double overlap_efficiency = 0.0;
    if (avg_interior > 0) {
      overlap_efficiency =
          std::min(1.0, avg_interior / (avg_interior + avg_wait)) * 100;
    }
    std::cout << "  Comm/Compute overlap efficiency: " << overlap_efficiency
              << "%\n";
    std::cout << "  (100% = communication fully hidden, lower = waiting for "
                 "network)\n";
    std::cout << "============================================================="
                 "===================\n";
  }
}

} // namespace mpi_solver
