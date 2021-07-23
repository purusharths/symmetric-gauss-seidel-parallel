#pragma once
#include <cmath>
#include <iostream>

void forward_gauss_sidel(int n, int k, double *grid, double *local_mean) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      std::cout << "\t" << i << ", " << j << std::endl;
      double sum = 0;
      double points = 1; // max_points;
      int row_min = std::max(i - k, 0);
      int row_max = std::min(i + k, n - 1);
      int col_min = std::max(j - k, 0);
      int col_max = std::min(j + k, n - 1);
      std::cout << "\t rowrange: (" << col_min << ", " << col_max << ")\n\t"
                << " colrange: (" << row_min << ", " << row_max << ")"
                << std::endl;
      for (int mm = row_min; mm <= row_max;
           mm++) { // std::max(j - k, 0); mm <= std::min(j + k, n - 1); mm++) {
        for (int ll = col_min; ll <= col_max;
             ll++) { // std::max(i - k, 0); ll <= std::min(i + k, n - 1); ll++)
                     // {
          //   std::cout << ll << " ";
          double alpha =
              (ll == mm) ? 50
                         : std::pow(std::max(std::abs(ll), std::abs(mm)), -2);
          std::cout << "\n->" << mm << ", " << ll << ": " << std::endl;
          if (mm < i || (mm == i && ll < j)) { // mm >= k && ll > i) {
            // think here: which array
            // B+
            std::cout << grid[ll * n + mm] << " B ";
            // sum += alpha * local_mean[ll * n + mm];
          } else {
            // A+
            // std::cout << std::endl;
            std::cout << grid[ll * n + mm] << " A ";
            // sum += alpha * grid[ll * n + mm];
          }

          sum += alpha * grid[ll * n + mm];
        }
        std::cout << std::endl;
      }
      std::cout << " -*- " << std::endl;
      local_mean[j * n + i] =
          static_cast<double>(sum / points); // / change to multiply
    }
  }
}
