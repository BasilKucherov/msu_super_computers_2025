#ifndef GRID_UTILS_HPP
#define GRID_UTILS_HPP

#include <cstddef>

namespace utils {
inline size_t idx3d(int i, int j, int k, int Ny, int Nz) {
  return k + Nz * (j + Ny * i);
}

inline int ip(int i, int N) { return (i + 1 == N) ? 0 : i + 1; }
inline int im(int i, int N) { return (i == 0) ? N - 1 : i - 1; }
} // namespace utils

#endif
