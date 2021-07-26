#include <iostream>
#include <omp.h>
// Pointer Bug!
#include "utilities.hh"

class SGS {
public:
  SGS(int n, int k, int iterations) : n(n), k(k), iterations(iterations) {
    coeff_stencil(k, alpha);
    fill_array(n, grid);
    for (int i1 = 0; i1 < n; i1++) {
      for (int i0 = 0; i0 < n; i0++) {
        local_mean[i1 * n + i0] = 0.0;
      }
    }
  }
  void calculate(int l, int i, int c, int K, int k, int n, double *alpha,
                 double *grid, double *local_mean);
  void symmetric_gauss_seidel(int n, int k, int iterations, double *grid,
                              double *local_mean, double *alpha);
  void start();
  ~SGS() {
    print_grid(local_mean, n);
    delete[] grid;
    delete[] local_mean;
    delete[] alpha;
  }

private:
  int n, k, iterations;
  double *grid = new (std::align_val_t(64)) double[n * n];
  double *local_mean = new (std::align_val_t(64)) double[n * n];
  int stencil_size = 2 * k + 1;
  double *alpha = new double[stencil_size * stencil_size];
};

void SGS::start() {
  symmetric_gauss_seidel(n, k, iterations, grid, local_mean, alpha);
}

void SGS::calculate(int l, int i, int c, int K, int k, int n, double *alpha,
                    double *grid, double *local_mean) {

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

void SGS::symmetric_gauss_seidel(int n, int k, int iterations, double *uold,
                                 double *unew, double *alpha) {
  double *grid = uold;
  double *local_mean = unew;
  double *a = alpha;
  int K = (2 * k + 1);
  for (int i = 0; i < iterations; i++) {
    std::cout << "jey";
    for (int i = 0; i < n; i++) {       // i loop
      for (int c = 0; c < k + 1; c++) { // c loop

        int loop_count = std::min(i, (n - 1 - c) / (k + 1)) + 1;

#pragma omp parallel for shared(alpha, grid, local_mean, K)
        for (int l = 0; l < loop_count; l++) {
          calculate(l, i, c, K, k, n, a, grid, local_mean);
        }
      }
      // loop for last few points
      for (int c = k + 1; c < n; c++) {
        int i = n - 1;
        int loop_count = 1 + std::floor((n - 1 - c) / (k + 1));
#pragma omp parallel for shared(alpha, grid, local_mean, K)
        for (int l = 0; l < loop_count; l++) {
          calculate(l, i, c, K, k, n, a, grid, local_mean);
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
            calculate(l, i, c, K, k, n, a, grid, local_mean);
          }
        }
      }

      for (int c = n - k - 1; c >= 0; c--) {
        int loop_count = 1 + std::min(c, (k + 1));
        int i = 0;
#pragma omp for
        for (int l = 0; l < loop_count; l++) {
          calculate(l, i, c, K, k, n, a, grid, local_mean);
        }
      }

      // grid-> u^1/2
      // local_mean -> u^1
      std::swap(local_mean, grid);
      //   grid->u^1
      // local_mean -> u^1/2
    }
  }

  if (local_mean != grid) {
    std::swap(grid, local_mean);
  }
}

int main(int argc, char **argv) {
  int n; // 0 -> n = n+1
  int k;
  int blocksize;
  int iterations;
  if (argc < 2 || argc > 4) {
    std::cout << "usage: ./sgs_omp <grid size> <stencil size (k)> <iterations>"
              << std::endl;
    exit(1);
  }
  sscanf(argv[1], "%d", &n);
  sscanf(argv[2], "%d", &k);
  // sscanf(argv[3], "%d", &blocksize);
  if (argc > 3) {
    sscanf(argv[3], "%d", &iterations);
  } else {
    iterations = 100;
  }
  std::cout << "Iterations: " << iterations;

  if (k > n) {
    std::cout << "Stencil Size must be smaller than Matrix Size" << std::endl;
    exit(1);
  }

  SGS sgs(n, k, iterations);
  sgs.start();
}