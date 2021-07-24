#include <future>
#include <iostream>
#include <thread>
#include <vector>
#include <algorithm>

void jacobi_vanilla(int n, int k, int i_min, int i_max, int j_min, int j_max,
                    double *grid, double *computed_vals, double *alpha) {
  double *uold = grid;
  double *unew = computed_vals;
  for (int iter = 0; iter < 1; iter++) {
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

int main() {
  //   int num_threads = 8;
  //   for (int i = 0; i < num_threads; i++) {
  //       std::async(jacobi, n, k, )
  //   }

  //   std::vector<> threads;
  //   auto f0 = std::async(jacobi_vanilla, n, k, __, __, __, __, grid,
  //                        computed_vals, alpha);
  // distribute rows
  const int SIZE = 100;
  int p = 8;
  int p_size = SIZE > p ? p : SIZE;

  std::vector<int> iter(p_size);

  for (int i = 0; i < SIZE; i++) {
    iter[i % p_size] += 1;
  }

  std::reverse(iter.begin(), iter.end());
  std::vector<std::tuple<int, int>> loop_domains;


  int val_start = 0;
  int val_end = 0;
  for (unsigned int j = 0; j < iter.size(); ++j) {
      val_end+=iter[j];
    printf("[%d] = %d | %d, %d\n", j, (int)iter[j], val_start, val_end);
    val_start += iter[j];
  }
  return 0;

  // swap pointers
}