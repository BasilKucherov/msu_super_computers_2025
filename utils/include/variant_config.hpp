#ifndef VARIANT_CONFIG_HPP
#define VARIANT_CONFIG_HPP

#include <cmath>
#include <functional>

namespace utils {
struct VariantConfig {
  double a_squared;
  double a_t;
  std::function<double(double, double, double, double)> u_analytical;

  static VariantConfig variant6(double Lx, double Ly, double Lz) {
    VariantConfig config;
    config.a_squared = 1.0 / 9.0;
    config.a_t = (M_PI / 3.0) *
                 std::sqrt(4.0 / (Lx * Lx) + 1.0 / (Ly * Ly) + 4.0 / (Lz * Lz));
    config.u_analytical = [Lx, Ly, Lz, a_t = config.a_t](double x, double y,
                                                         double z, double t) {
      return std::sin(2.0 * M_PI * x / Lx) * std::sin(M_PI * y / Ly + M_PI) *
             std::sin(2.0 * M_PI * z / Lz + 2.0 * M_PI) *
             std::cos(a_t * t + M_PI);
    };

    return config;
  }
};
} // namespace utils

#endif
