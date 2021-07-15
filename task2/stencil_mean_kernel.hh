#pragma once
#include <iostream>

void stencil_mean_vanilla(int n, int k, double *grid, double *local_mean) {
  double max_points = (2 * k + 1) * (2 * k + 1);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double sum = 0;
      int m = j - k;
      int l = i - k;
      int max_m = j + k;
      int max_l = i + k;
      double points = 0;//max_points;
      // if i,j are boundry values, then k-point stencil will not be a perfect
      // square and will change its shape
      if (i < k || j < k || i > n - k - 1 || j > n - k - 1) {
        if (j - k < 0)
          m = 0;
        if (i - k < 0)
          l = 0;
        if (j + k > n)
          max_m = n;
        if (j + k > n)
          max_l = n;
        // points = (max_m - m) * (max_l - l);
      }
      for (int mm = m; mm <= max_m; mm++) {
        for (int ll = l; ll <= max_l; ll++) {
          sum += grid[ll * n + mm];
          points+=1;
        }
        // local_mean[j * n + i] = sum; previous location
      }
      local_mean[j * n + i] = static_cast<double>(sum/points); // / 
    }
    std::cout << std::endl;
  }
}