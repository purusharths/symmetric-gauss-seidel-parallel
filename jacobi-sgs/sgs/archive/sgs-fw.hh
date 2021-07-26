#pragma once
#include <cmath>
#include <iostream>

void forward_gauss_sidel(int n, int k, double *grid, double *local_mean) {
  std::cout << "Vanilla GS" << std::endl;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double sum = 0;
      double points = 1; // max_points;
      int row_min = std::max(i - k, 0);
      int row_max = std::min(i + k, n - 1);
      int col_min = std::max(j - k, 0);
      int col_max = std::min(j + k, n - 1);

      for (int mm = row_min; mm <= row_max; mm++) {
        for (int ll = col_min; ll <= col_max; ll++) {

          if (mm < i || (mm == i && ll < j)) { // mm >= k && ll > i) {
            // think here: which array
            // B+
            std::cout << local_mean[ll * n + mm] << " " << ll << ", " << mm ;
            sum += local_mean[ll * n + mm];
          } else {
            // A+
            sum += grid[ll * n + mm];
          }
        }
      }std::cout << std::endl;

      local_mean[j * n + i] =
          static_cast<double>(sum / points); // / change to multiply
    }
  }
}