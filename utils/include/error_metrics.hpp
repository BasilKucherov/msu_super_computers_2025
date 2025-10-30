#ifndef ERROR_METRICS_HPP
#define ERROR_METRICS_HPP

#include <functional>
#include <vector>

namespace utils {

struct ErrorMetrics {
  double max_abs_error;
  double mse;
};

void computeAnalyticalSolution(
    std::vector<double> &u_analytical_out,
    const std::function<double(double, double, double, double)>
        &u_analytical_func,
    double Lx, double Ly, double Lz, int N, double t);

ErrorMetrics computeErrorMetrics(const std::vector<double> &u_numerical,
                                 const std::vector<double> &u_analytical);

} // namespace utils

#endif
