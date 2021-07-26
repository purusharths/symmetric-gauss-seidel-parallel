#pragma once
#include <cmath>
#include <iostream>
#include <omp.h>
#include <tuple>
#include <vector>

void forward_gauss_sidel_omp(int n, int k, double *grid, double *local_mean) {
  std::vector<std::tuple<int, int>> wave_coords(n);
  for (int i = 0; i < n; i++) {       // i loop
    for (int c = 0; c < k + 1; c++) { // c loop

      double points = 1; // max_points;

      int p = i;
      int q = c, l;
      int loop_count = std::min(p, (n - 1 - q) / (k + 1)) + 1;

      for (int l = 1; l <= loop_count; l++) {
        wave_coords.emplace_back(p, q);
        p = p - 1;
        q = q + (k + 1);
      }

    #pragma omp parallel for
      for (int v = 0; v < loop_count; v++) {
        int p = std::get<0>(wave_coords[v]);
        int q = std::get<1>(wave_coords[v]);

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
        #pragma omp simd reduction(+:sum)
        for (int mm = row_min; mm <= row_max; mm++) {
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
      }
      wave_coords.clear();
    }
  }
}