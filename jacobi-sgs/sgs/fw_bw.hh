#pragma once
#include <cmath>
#include <iostream>
#include <omp.h>

void gauss_sidel_omp(int n, int k, int iterations, double *grid, double *local_mean,
                     double *alpha) {
  int K = (2 * k + 1);
  for (int i = 0; i < iterations; i++) {
    for (int i = 0; i < n; i++) {       // i loop
      for (int c = 0; c < k + 1; c++) { // c loop

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

          for (int mm = row_min; mm <= row_max; mm++) {
#pragma omp simd reduction(+ : sum)
            for (int ll = col_min; ll <= col_max; ll++) {
              double alph = alpha[(mm - (p - k)) * K + (ll - (q - k))];
              // double alph = 1;
              if (mm < p || (mm == p && ll < q)) { // B+
                sum += alph * local_mean[ll * n + mm];
              } else { // A+
                sum += alph * grid[ll * n + mm];
              }
            }
          }
          // #pragma omp critical
          local_mean[q * n + p] = static_cast<double>(sum);
        }
      }
    }
    // loop for last few points
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
#pragma omp simd reduction(+ : sum)
          for (int ll = col_min; ll <= col_max; ll++) {
            double alph = alpha[(mm - (p - k)) * K + (ll - (q - k))];
            // double alph = 1;
            if (mm < p || (mm == p && ll < q)) {
              // B+
              sum += alph * local_mean[ll * n + mm];
            } else {
              // A+
              sum += alph * grid[ll * n + mm];
            }
          }
        }
#pragma omp critical
        local_mean[q * n + p] = static_cast<double>(sum);
      }
      // #pragma omp barrier
    }
// end fw
// end fw
// end fw
// -------------------
#pragma omp barrier
    // at the end of forward iteration
    // grid -> u^0
    // local_mean -> u^(1/2)
    std::swap(grid, local_mean);
    // Now:
    // grid -> u^(1/2)
    // local_mean -> u^0 (will be overwritten by u^1)

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
              double alph = alpha[(mm - (p - k)) * K + (ll - (q - k))];
              if (mm > p || (mm == p && ll > q)) {
                // B -
                sum += alph * local_mean[ll * n + mm];
              } else {
                // A -
                sum += alph * grid[ll * n + mm];
              }
            }
          }
          // #pragma omp critical
          local_mean[q * n + p] = static_cast<double>(sum); //
        }
      }
    }

    for (int c = n - k - 1; c >= 0; c--) {
      int loop_count = 1 + std::min(c, (k + 1));
      int i = 0;
#pragma omp for
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
            double alph = alpha[(mm - (p - k)) * K + (ll - (q - k))];
            if (mm > p || (mm == p && ll > q)) {
              // B -
              sum += alph * local_mean[ll * n + mm];
            } else {
              // A -
              sum += alph * grid[ll * n + mm];
            }
          }
        }
        // #pragma omp critical
        local_mean[q * n + p] = static_cast<double>(sum); //
      }
    }

    // grid-> u^1/2
    // local_mean -> u^1
    std::swap(local_mean, grid);
    //   grid->u^1
    // local_mean -> u^1/2
  }

    if (local_mean != grid) {
    std::swap(grid, local_mean);
  }

}