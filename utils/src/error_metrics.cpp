#include "error_metrics.hpp"
#include <cmath>

namespace utils {
void computeAnalyticalSolution(
    std::vector<double> &u_analytical_out,
    const std::function<double(double, double, double, double)>
        &u_analytical_func,
    double Lx, double Ly, double Lz, int N, double t) {

  double hx = Lx / N;
  double hy = Ly / N;
  double hz = Lz / N;

  for (int i = 0; i < N; ++i) {
    double x = i * hx;
    for (int j = 0; j < N; ++j) {
      double y = j * hy;
      for (int k = 0; k < N; ++k) {
        double z = k * hz;
        size_t idx = k + N * (j + N * i);
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

  metrics.mse = std::sqrt(sum_squared_error / static_cast<double>(N));

  return metrics;
}
} // namespace utils
