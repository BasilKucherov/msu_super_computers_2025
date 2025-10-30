#ifndef GRID_UTILS_HPP
#define GRID_UTILS_HPP

#include <cstddef>

namespace utils {
inline size_t idx3d(int i, int j, int k, int N) { return k + N * (j + N * i); }
} // namespace utils

#endif
