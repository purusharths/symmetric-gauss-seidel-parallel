#include <cmath>
#include <iostream>
#include <vector>

#include "stencil_vanilla.hh"
#include "stencil_blocked.hh"
#include "stencil_simd.hh"
#include "stencil_omp_simd.hh"
#include "time_experiment.hh"
#include "utilities.hh"

using NUMBER = double;
const int N = 32 * 1024 * 1024; // problem size

//  stencil_mean_vanilla(n, k, grid, local_mean);
// stencil_mean_blocked(n, k, grid, local_mean, blocksize);
// stencil_simd(n, k, grid, local_mean, blocksize);

class StencilVanilla { //: public Experiment {
  int n, k;
  double *grid = new (std::align_val_t(64)) double[n * n];
  double *local_mean = new (std::align_val_t(64)) double[n * n];

public:
  // construct an experiment
  StencilVanilla(int n_, int k_) : n(n_), k(k_) {
    fill_array(n, grid);
    fill_array(n, local_mean);
  }
  ~StencilVanilla(){
    delete [] grid;
    delete [] local_mean;
  }
  // run an experiment; can be called several times
  void run() const { stencil_mean_vanilla(n, k, grid, local_mean); }
  // report number of operations
  double operations() const {
    return (n * n) * std::pow(2 * k + 1, 2) - 1 + 1;
  } // no operations
};
// ---------------------------------------------------------------------------------

class StencilBlocked { //: public Experiment {
  int n, k, blocked;
  double *grid = new (std::align_val_t(64)) double[n * n];
  double *local_mean = new (std::align_val_t(64)) double[n * n];

public:
  // construct an experiment
  StencilBlocked(int n_, int k_, int b) : n(n_), k(k_), blocked(b) {
    fill_array(n, grid);
    fill_array(n, local_mean);
  }
  ~StencilBlocked(){
    delete [] grid;
    delete [] local_mean;
  }
  // run an experiment; can be called several times
  void run() const { stencil_mean_blocked(n, k, grid, local_mean, blocked); }
  // report number of operations
  double operations() const {
    return (n * n) * std::pow(2 * k + 1, 2) - 1 + 1;
  } // no operations
};
// ---------------------------------------------------------------------------------

class StencilSIMD_OMP { //: public Experiment {
  int n, k, blocked;
  double *grid = new (std::align_val_t(64)) double[n * n];
  double *local_mean = new (std::align_val_t(64)) double[n * n];

public:
  // construct an experiment
  StencilSIMD_OMP(int n_, int k_, int b) : n(n_), k(k_), blocked(b) {
    fill_array(n, grid);
    fill_array(n, local_mean);
  }
  ~StencilSIMD_OMP(){
    delete [] grid;
    delete [] local_mean;
  }
  // run an experiment; can be called several times
  void run() const { stencil_simd(n, k, grid, local_mean, blocked); }
  // report number of operations
  double operations() const {
    return (n * n) * std::pow(2 * k + 1, 2) - 1 + 1;
  } // no operations
};

class StencilSIMD { //: public Experiment {
  int n, k, blocked;
  double *grid = new (std::align_val_t(64)) double[n * n];
  double *local_mean = new (std::align_val_t(64)) double[n * n];

public:
  // construct an experiment
  StencilSIMD(int n_, int k_, int b) : n(n_), k(k_), blocked(b) {
    fill_array(n, grid);
    fill_array(n, local_mean);
  }
  ~StencilSIMD(){
    delete [] grid;
    delete [] local_mean;
  }
  // run an experiment; can be called several times
  void run() const { stencil_simd(n, k, grid, local_mean, blocked); }
  // report number of operations
  double operations() const {
    return (n * n) * std::pow(2 * k + 1, 2) - 1 + 1;
  } // no operations
};

class StencilSIMD_LargeVec { //: public Experiment {
  int n, k, blocked;
  double *grid = new (std::align_val_t(64)) double[n * n];
  double *local_mean = new (std::align_val_t(64)) double[n * n];

public:
  // construct an experiment
  StencilSIMD_LargeVec(int n_, int k_, int b) : n(n_), k(k_), blocked(b) {
    fill_array(n, grid);
    fill_array(n, local_mean);
  }
  ~StencilSIMD_LargeVec(){
    delete [] grid;
    delete [] local_mean;
  }
  // run an experiment; can be called several times
  void run() const { stencil_simd(n, k, grid, local_mean, blocked); }
  // report number of operations
  double operations() const {
    return (n * n) * std::pow(2 * k + 1, 2) - 1 + 1;
  } // no operations
};
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
void benchmark(int k, int blocksize) {
  // std::cout << N * sizeof(NUMBER) / 1024 / 1024 << " MByte per vector"
  //           << std::endl;
  double time_factor = 1e6;
  std::vector<int> sizes = {200, 300, 400, 500, 600, 700, 800};

  std::cout << "Stencil Size: " << k << ", Block Size:" << blocksize << std::endl;
  // std::cout << "**** Vanilla ****\n";
  std::cout << "experiment,n,time,repetitions,Gflops/s,GByte/s" << std::endl;
  for (auto n : sizes) {
    if ((n) % blocksize != 0) {
      std::cout << "Blocksize must be a multiple of n" << std::endl;
      exit(1);
    }

    StencilVanilla e(n, k);
    auto d = time_experiment(e);
    double flops = d.first * e.operations() / d.second * 1e6 / 1e9;

    std::cout << "vanilla, " << n << ", " << d.second / time_factor << ", " << d.first << ", "
              << flops << ", " << flops * sizeof(NUMBER) << std::endl;
  }
  // std::cout << "\n\n----Experiment Over---- \n\n";

  // ---------------------------------------------------------------------------------
  // std::cout << "**** Blocked ****\n";
  // std::cout << "experiment, n, time (us), repetitions, Gflops/s, GByte/s"
            // << std::endl;
  for (auto n : sizes) {
    if ((n) % blocksize != 0) {
      std::cout << "Blocksize must be a multiple of n" << std::endl;
      exit(1);
    }

    StencilBlocked e(n, k, blocksize);
    auto d = time_experiment(e);
    double flops = d.first * e.operations() / d.second * 1e6 / 1e9;

    std::cout << "blocked, " << n << ", " << d.second / time_factor << ", " << d.first << ", "
              << flops << ", " << flops * sizeof(NUMBER) << std::endl;
  }
  // std::cout << "\n\n----Experiment Over---- \n\n";
  // ---------------------------------------------------------------------------------

  // std::cout << "**** SIMD OMP ****\n";
  // std::cout << "experiment, n, time (us), repetitions, Gflops/s, GByte/s"
            // << std::endl;
  for (auto n : sizes) {
    if ((n) % blocksize != 0) {
      std::cout << "Blocksize must be a multiple of n" << std::endl;
      exit(1);
    }

    StencilSIMD_OMP e(n, k, blocksize);
    auto d = time_experiment(e);
    double flops = d.first * e.operations() / d.second * 1e6 / 1e9;

    std::cout << "simd omp, " << n << ", " << d.second / time_factor<< ", " << d.first
              << ", " << flops << ", " << flops * sizeof(NUMBER) << std::endl;
  }
  // std::cout << "\n\n----Experiment Over---- \n\n";

  // ---------------------------------------------------------------------------------

  // std::cout << "**** SIMD ****\n";

  // std::cout << "experiment, n, time (us), repetitions, Gflops/s, GByte/s"
            // << std::endl;

  for (auto n : sizes) {
    if ((n) % blocksize != 0) {
      std::cout << "Blocksize must be a multiple of n" << std::endl;
      exit(1);
    }

    StencilSIMD e(n, k, blocksize);
    auto d = time_experiment(e);
    double flops = d.first * e.operations() / d.second * 1e6 / 1e9;
    std::cout << "simd 4d, " << n << ", " << d.second / time_factor<< ", " << d.first << ", "
              << flops << ", " << flops * sizeof(NUMBER) << std::endl;
  }
  // std::cout << "\n\n----Experiment Over---- \n\n";

  // ---------------------------------------------------------------------------------

  // std::cout << "**** SIMD Large Vec ****\n";
  // std::cout << "experiment, n, time (us), repetitions, Gflops/s, GByte/s"
            // << std::endl;
  for (auto n : sizes) {
    if ((n) % blocksize != 0) {
      std::cout << "Blocksize must be a multiple of n" << std::endl;
      exit(1);
    }
    if (k > 6) {
      StencilSIMD_LargeVec e(n, k, blocksize);
      auto d = time_experiment(e);
      double flops = d.first * e.operations() / d.second * 1e6 / 1e9;
      std::cout << "simd 8d, " << n << ", " << d.second / time_factor<< ", " << d.first
                << ", " << flops << ", " << flops * sizeof(NUMBER) << std::endl;
    } else {
      std::cout << "Stencil smaller than vector size. Skipping..." << std::endl;
    }
    // std::cout << "n=" << n << " took " << d.second << " us for " << d.first
    //           << " repetitions"
    //           << " " << flops << " Gflops/s"
    //           << " " << flops * sizeof(NUMBER) << " GByte/s" << std::endl;
  }
  // std::cout << "\n\n----Experiment Over---- \n\n";
}

int main(int argc, char **argv) {
  int k, blocksize;
  sscanf(argv[1], "%d", &k);
  sscanf(argv[2], "%d", &blocksize);
  benchmark(k, blocksize);
}