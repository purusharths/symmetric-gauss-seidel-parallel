#pragma once
#include <cmath>
#include <iostream>
#include <omp.h>

void forward_gauss_sidel_omp(int n, int k, double *grid, double *local_mean) {
  int j = 0;
  for (int i = 0; i < n; i++) {       // i loop
    for (int c = 0; c < k + 1; c++) { // c loop

      double points = 1; // max_points;

      int p = i;
      int q = c;
      int loop_count = std::min(p, (n - 1 - q) / (k + 1)) + 1;

#pragma omp parallel for shared(p, q) 
      for (int l = 1; l <= loop_count; l++) {
        double sum = 0;

        int row_min = std::max(p - k, 0);
        int row_max = std::min(p + k, n - 1);
        int col_min = std::max(q - k, 0);
        int col_max = std::min(q + k, n - 1);
        // loop for k - box
        for (int mm = row_min; mm <= row_max; mm++) {
          for (int ll = col_min; ll <= col_max; ll++) {

            if (mm < p || (mm == p && ll < q)) { // mm >= k && ll > i) {
              // B+
              // #pragma omp critical
              sum += local_mean[ll * n + mm];
            } else {
              // A+
              // #pragma omp critical
              sum += grid[ll * n + mm];
            }
          }
        }
        // std::cout << p << ", " << q << ": " << sum << "|\t Rowmin/max:: " <<
        // row_min
        // << ", " << row_max <<"\t Col Min/Max: " << col_min << ", "<<col_max
        // << std::endl;
        local_mean[q * n + p] =
            static_cast<double>(sum / points); // / change to multiply
                                               // #pragma omp barrier
#pragma omp critical
        {
          p = p - 1;
          q = q + (k + 1);
        }
      }
    }
  }
}
