#pragma once
#include <cmath>
#include <iostream>

void backward_gauss_sidel(int n, int k, double *grid, double *local_mean) {
  for (int i = n - 1; i >= 0; i--) {
    for (int j = n - 1; j >= 0; j--) {
      double sum = 0;
      double points = 1; // max_points;
      //   int row_min = std::max(i - k, 0);
      //   int row_max = std::min(i + k, n - 1);
      std::cout << grid[j * n + i] << std::endl;
      for (int mm = std::max(j - k, 0); mm <= std::min(j + k, n - 1); mm++) {
        for (int ll = std::max(i - k, 0); ll <= std::min(i + k, n - 1); ll++) {
          // sum += grid[ll *  n + mm];
          if (mm > i || (mm == i && ll > j)) {
            // B -
          } else {
            // A -
          }
        }
      }

      local_mean[j * n + i] =
          static_cast<double>(sum / points); // / change to multiply
      std::cout << sum << std::endl;
    }
  }
}
