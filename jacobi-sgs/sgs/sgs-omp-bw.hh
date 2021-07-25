#pragma once
#include <cmath>
#include <iostream>
#include <omp.h>

void backward_gauss_sidel_omp(int n, int k, double *grid, double *local_mean) {

  for (int i = n - 1; i >= 0; i--) {           // i loop
    for (int c = n - 1; c >= n - k - 1; c--) { // c loop

      int loop_count = 1 + std::min(n - 1 - i, (c / (k + 1)));

#pragma omp parallel for
      for (int l = 0; l < loop_count; l++) {
        int p = i + l;
        int q = c - (l * (k + 1));

        double sum = 0;

        int row_min = std::max(p - k, 0);
        int row_max = std::min(p + k, n - 1);
        int col_min = std::max(q - k, 0);
        int col_max = std::min(q + k, n - 1);

        for (int mm = row_min; mm <= row_max; mm++) {
          for (int ll = col_min; ll <= col_max; ll++) {

            if (mm > p || (mm == p && ll > q)) {
              // B -
              sum += local_mean[ll * n + mm];
            } else {
              // A -
              sum += grid[ll * n + mm];
            }
          }
        }

        local_mean[q * n + p] = static_cast<double>(sum); //
      }
    }
  }

  for (int c = n - k - 1; c >= 0; c--) {
    int loop_count = 1 + std::min(c, (k + 1));
    int i = 0;
    #pragma omp parallel for
    for (int l = 0; l < loop_count; l++) {
      int p = i + l;
      int q = c - (l * (k + 1));

      double sum = 0;

      int row_min = std::max(p - k, 0);
      int row_max = std::min(p + k, n - 1);
      int col_min = std::max(q - k, 0);
      int col_max = std::min(q + k, n - 1);

      for (int mm = row_min; mm <= row_max; mm++) {
        for (int ll = col_min; ll <= col_max; ll++) {

          if (mm > p || (mm == p && ll > q)) {
            // B -
            sum += local_mean[ll * n + mm];
          } else {
            // A -
            sum += grid[ll * n + mm];
          }
        }
      }

      local_mean[q * n + p] = static_cast<double>(sum); //
    }
  }
}