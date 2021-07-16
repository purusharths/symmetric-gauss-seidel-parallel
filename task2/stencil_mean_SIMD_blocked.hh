#pragma once
#include <algorithm>
#include <iostream>

void stencil_mean_SIMD_blocked(int n, int k, double *grid, double *local_mean, int blocksize = 4) {
  for (int i = 0; i < n; i = i + (blocksize)) {
    for (int j = 0; j < n; j = j + (blocksize)) {

    for (int blockRow = i; blockRow < i + blocksize; blockRow++) {
        for (int blockCol = j; blockCol < j + blocksize; blockCol++) {
          double sum = 0, points = 0;
          for (int mm = std::max(blockCol - k, 0);
               mm <= std::min(blockCol + k, n - 1); mm++) {
            for (int ll = std::max(blockRow - k, 0);
                 ll <= std::min(blockRow + k, n - 1); ll++) {
              sum += grid[ll * n + mm];
              points++;
            }
          }
          local_mean[blockCol * n + blockRow] =
              static_cast<double>(sum / points); // /
          sum = 0;
        }
      }
    }
  }
}
