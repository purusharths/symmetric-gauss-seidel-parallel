#pragma once
#include <iostream>

void stencil_mean_vanilla(int n, int k, double *grid, double *local_mean) {
  //   double max_points = (2 * k + 1) * (2 * k + 1);
  // std::cout << "vanilla" << std::endl;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double sum = 0;
      double points = 0; // max_points;
      int col_max = std::min(j + k, n - 1);
      int col_min = std::max(j - k, 0);
      int row_max = std::min(i + k, n - 1);
      int row_min = std::max(i - k, 0);
      int m = col_max - col_min;
      int o = row_max - row_min;
      int mn = (m + 1) * (o + 1);
      // for (int mm = std::max(j - k, 0); mm <= std::min(j + k, n - 1); mm++) {
      // for (int ll = std::max(i - k, 0); ll <= std::min(i + k, n - 1); ll++) {

      for (int mm = row_min; mm <= row_max; mm++) {
        for (int ll = col_min; ll <= col_max; ll++) {
          sum += grid[ll * n + mm];
          // points += 1;
        }
        // local_mean[j * n + i] = sum; previous location
      }
      local_mean[j * n + i] =
          static_cast<double>(sum / mn); // / change to multiply
    }
  }
}