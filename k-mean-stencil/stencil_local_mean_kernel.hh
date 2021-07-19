#pragma once
#include <iostream>

void stencil_mean_vanilla(int n, int k, double *grid, double *local_mean) {
  //   double max_points = (2 * k + 1) * (2 * k + 1);
  // std::cout << "vanilla" << std::endl;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double sum = 0;
      double points = 0; // max_points;
      for (int mm = std::max(j - k, 0); mm <= std::min(j + k, n - 1); mm++) {
        for (int ll = std::max(i - k, 0); ll <= std::min(i + k, n - 1); ll++) {
          sum += grid[ll * n + mm];
          points += 1;
        }
        // local_mean[j * n + i] = sum; previous location
      }
      local_mean[j * n + i] = static_cast<double>(sum / points); // / change to multiply
    }
  }
}