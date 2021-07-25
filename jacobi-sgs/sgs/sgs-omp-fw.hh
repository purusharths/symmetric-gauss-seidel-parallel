#pragma once
#include <cmath>
#include <iostream>
#include <omp.h>

void forward_gauss_sidel_omp(int n, int k, double *grid, double *local_mean) {
  // int j = 0;
  for (int i = 0; i < n; i++) {       // i loop
    for (int c = 0; c < k + 1; c++) { // c loop

      double points = 1; // max_points;

      // int p = i;
      // int q = c, l;
      int loop_count = std::min(i, (n - 1 - c) / (k + 1)) + 1;

#pragma omp parallel for
      for (int l = 0; l < loop_count; l++) {
        int p = i - (1 * l);
        int q = c + ((k + 1) * l);
        double sum = 0;

        int row_min = std::max(p - k, 0);
        int row_max = std::min(p + k, n - 1);
        int col_min = std::max(q - k, 0);
        int col_max = std::min(q + k, n - 1);
        // int tid = omp_get_thread_num();

        // loop for k - box
        // #pragma omp critical
        //         std::cout << "(" << p << ", " << q << ") with thread: " <<
        //         std::endl;
        #pragma omp parallel for
        for (int mm = row_min; mm <= row_max; mm++) {
          #pragma omp simd reduction(+:sum)
          for (int ll = col_min; ll <= col_max; ll++) {
            if (mm < p || (mm == p && ll < q)) { // mm >= k && ll > i) {
                                                 // B+

              sum += local_mean[ll * n + mm];
            } else {
              // A+
              sum += grid[ll * n + mm];
            }
          }
        }

        local_mean[q * n + p] =
            static_cast<double>(sum / points); // / change to multiply
                                               // p = p - 1;
                                               // q = q + (k + 1);
      }
    }
  }

  for (int c = k + 1; c < n; c++) {
    int i = n - 1;
    int loop_count = 1 + std::floor((n - 1 - c) / (k + 1));
#pragma omp parallel for
    for (int l = 0; l < loop_count; l++) {
      int p = i - l;
      int q = c + (l * (k + 1));
        double sum = 0;

        int row_min = std::max(p - k, 0);
        int row_max = std::min(p + k, n - 1);
        int col_min = std::max(q - k, 0);
        int col_max = std::min(q + k, n - 1);

        for (int mm = row_min; mm <= row_max; mm++) {
          #pragma omp simd reduction(+:sum)
          for (int ll = col_min; ll <= col_max; ll++) {
            if (mm < p || (mm == p && ll < q)) { 
              // B+
              sum += local_mean[ll * n + mm];
            } else {
              // A+
              sum += grid[ll * n + mm];
            }
          }
        }
        local_mean[q * n + p] =
            static_cast<double>(sum);

    }
  }

}