#pragma once
#include <cmath>
#include <iostream>

void backward_gauss_sidel(int n, int k, double *grid, double *local_mean) {
  for (int i = n - 1; i >= 0; i--) {
    for (int j = n - 1; j >= 0; j--) {
      double sum = 0;
      double points = 1; // max_points;
      int row_min = std::max(i - k, 0);
      int row_max = std::min(i + k, n - 1);
      int col_min = std::max(j - k, 0);
      int col_max = std::min(j + k, n - 1);
      std::cout << "\t rowrange: (" << col_min << ", " << col_max << ")\n\t"
                << " colrange: (" << row_min << ", " << row_max
                << ") :    " << grid[i * n + j] << std::endl;
      for (int mm = row_min; mm <= row_max; mm++) {
        for (int ll = col_min; ll <= col_max; ll++) {
          std::cout << "\n->" << mm << ", " << ll << ": " << std::endl;
          // sum += grid[ll *  n + mm];
          if (mm > i || (mm == i && ll > j)) {
            // B -
            std::cout << grid[ll * n + mm] << " B- ";
          } else {
            // A -
            std::cout << grid[ll * n + mm] << " A- ";
          }
        }
      }

      local_mean[j * n + i] =
          static_cast<double>(sum / points); // / change to multiply
      // std::cout << sum << std::endl;
    }
  }
}
