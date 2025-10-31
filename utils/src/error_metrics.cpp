#include "error_metrics.hpp"
#include "grid_utils.hpp"
#include <cmath>

namespace utils {
void computeAnalyticalSolution(
    std::vector<double> &u_analytical_out,
    const std::function<double(double, double, double, double)>
        &u_analytical_func,
    double Lx, double Ly, double Lz, int Nx, int Ny, int Nz, double t) {

  const double hx = Lx / Nx;
  const double hy = Ly / (Ny - 1);
  const double hz = Lz / Nz;

  for (int i = 0; i < Nx; ++i) {
    double x = i * hx;
    for (int j = 0; j < Ny; ++j) {
      double y = j * hy;
      for (int k = 0; k < Nz; ++k) {
        double z = k * hz;
        size_t idx = idx3d(i, j, k, Nx, Ny, Nz);
        u_analytical_out[idx] = u_analytical_func(x, y, z, t);
      }
    }
  }
}

ErrorMetrics computeErrorMetrics(const std::vector<double> &u_numerical,
                                 const std::vector<double> &u_analytical) {
  ErrorMetrics metrics;
  metrics.max_abs_error = 0.0;
  double sum_squared_error = 0.0;

  size_t N = u_numerical.size();

  for (size_t i = 0; i < N; ++i) {
    double error = std::abs(u_numerical[i] - u_analytical[i]);
    metrics.max_abs_error = std::max(metrics.max_abs_error, error);
    sum_squared_error += error * error;
  }

  metrics.rmse = std::sqrt(sum_squared_error / static_cast<double>(N));

  return metrics;
}
} // namespace utils
