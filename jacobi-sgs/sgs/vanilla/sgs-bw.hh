#pragma once
#include <cmath>
#include <iostream>

void backward_gauss_sidel(int n, int k, double *grid, double *local_mean,
                          double *alpha) {
  // std::cout << "Backward Gauss sidel" << std::endl;
  int K = (2 * k + 1);
  for (int i = n - 1; i >= 0; i--) {
    for (int j = n - 1; j >= 0; j--) {

      double points = 1; // max_points;
      int row_min = std::max(i - k, 0);
      int row_max = std::min(i + k, n - 1);
      int col_min = std::max(j - k, 0);
      int col_max = std::min(j + k, n - 1);

      double sum = 0;
      for (int mm = row_min; mm <= row_max; mm++) {
        for (int ll = col_min; ll <= col_max; ll++) {
          double alph = alpha[(mm - (i - k)) * K + (ll - (j - k))];
          if (mm > i || (mm == i && ll > j)) {
            // B -
            sum += alph * local_mean[ll * n + mm];
          } else {
            // A -
            sum += alph * grid[ll * n + mm];
          }
        }
      }

      local_mean[j * n + i] =
          static_cast<double>(sum); // / change to multiply
      // std::cout << sum << std::endl;
    }
  }
}
