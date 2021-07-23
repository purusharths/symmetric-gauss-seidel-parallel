#pragma once
#include <cmath>
#include <iostream>
#include <type_traits>
#include <utility>
#include <omp.h>

void jacobi_vanilla(int n, int k, int iterations, double *grid,
                    double *computed_vals, double *alpha) {
  double *uold = grid;
  double *unew = computed_vals;
  for (int iter = 0; iter < iterations; iter++) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double sum = 0;
        double points = 0; // max_points;
        int row_min = std::max(i - k, 0);
        int row_max = std::min(i + k, n - 1);
        int col_min = std::max(j - k, 0);
        int col_max = std::min(j + k, n - 1);
        
        int c = col_max - col_min;
        int r = row_max - row_min;
        int size = (c + 1) * (r + 1);
        int K = (2 * k + 1);

        // std::max(j - k, 0); mm <= std::min(j + k, n - 1); mm++) WRONG
        for (int mm = row_min; mm <= row_max; mm++) {
          // std::max(i - k, 0); ll <= std::min(i + k, n - 1); WRONG
          for (int ll = col_min; ll <= col_max; ll++) {
            double alph = alpha[(mm - (i - k)) * K + (ll - (j - k))];
            alph = 1;
            sum += alph * uold[ll * n + mm];
          }
        }
        unew[j * n + i] = static_cast<double>(sum / size);
      }
    }
    std::swap(uold, unew);
  }
  if (computed_vals != uold) {
    std::swap(grid, computed_vals);
  }
}


void jacobi_omp(int n, int k, int iterations, double *grid,
                    double *computed_vals, double *alpha) {
  double *uold = grid;
  double *unew = computed_vals;
  for (int iter = 0; iter < iterations; iter++) {
    #pragma omp parallel for num_threads(8)
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double sum = 0;
        double points = 0; // max_points;
        int row_min = std::max(i - k, 0);
        int row_max = std::min(i + k, n - 1);
        int col_min = std::max(j - k, 0);
        int col_max = std::min(j + k, n - 1);
        
        int c = col_max - col_min;
        int r = row_max - row_min;
        int size = (c + 1) * (r + 1);
        int K = (2 * k + 1);

        // std::max(j - k, 0); mm <= std::min(j + k, n - 1); mm++) WRONG
        for (int mm = row_min; mm <= row_max; mm++) {
          // std::max(i - k, 0); ll <= std::min(i + k, n - 1); WRONG
          #pragma omp simd reduction(+ : sum)
          for (int ll = col_min; ll <= col_max; ll++) {
            // double alph = alpha[(mm - (i - k)) * K + (ll - (j - k))];
            // alph = 1;
            sum += uold[ll * n + mm];
          }
        }
        unew[j * n + i] = static_cast<double>(sum / size);
      }
    }
    #pragma omp critical
    std::swap(uold, unew);
  }
  if (computed_vals != uold) {
    std::swap(grid, computed_vals);
  }
}