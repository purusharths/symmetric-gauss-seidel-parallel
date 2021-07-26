#pragma once
#include <cmath>
#include <iostream>
#include <omp.h>
#include <type_traits>
#include <utility>

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

        for (int mm = row_min; mm <= row_max; mm++) {
#pragma omp simd reduction(+ : sum)
          for (int ll = col_min; ll <= col_max; ll++) {
            double alph = alpha[(mm - (i - k)) * K + (ll - (j - k))];
            sum += alph * uold[ll * n + mm];
          }
        }
        unew[j * n + i] = static_cast<double>(sum);
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
#pragma omp parallel for
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

        for (int mm = row_min; mm <= row_max; mm++) {
#pragma omp simd reduction(+ : sum)
          for (int ll = col_min; ll <= col_max; ll++) {
            double alph = alpha[(mm - (i - k)) * K + (ll - (j - k))];
            // alph = 1;
            sum += alph * uold[ll * n + mm];
          }
        }
        unew[j * n + i] = static_cast<double>(sum);
      }
    }
#pragma omp barrier // implicit with opennp. Don't be paranoid.
    std::swap(uold, unew);
  }
  if (computed_vals != uold) {
    std::swap(grid, computed_vals);
  }
}