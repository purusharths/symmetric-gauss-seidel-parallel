Install python dependency for plotting:
`pip install -r requirements.txt`
CPP dependency: VCL, OpenMP

The folder structure is as follows:
```
.
├── code
│   ├── jacobi-sgs
│   │   ├── jacobi
│   │   │   ├── ...
│   │   └── sgs
│   │       ├── ...
│   │       └── vanilla
│   │           ├── ...
│   └── k-mean-stencil
│       ├── ...
├── README.md
├── requirements.txt
```

## K-Stencil Mean
For local mean kernel.
To compile: `make` inside the `k-mean-stencil` to get the `bm_stencil_mean`.
Command line args are: `./bm_stencil_mean <k> <block size>`

## Jacobi
For jacobi iteration.
To compile: `make` inside `jaccobi` folder. `jacobi`, `jacobi_seq` and `bm_jacobi` files are created.
`jacobi` is the paralel version and `jacobi_seq` is the sequential version. <br>
Command line args are:
(main) `./ (jacobi/jacobi_seq) <grid size> <stencil size (k)> <iterations>` if iteration flag is not given, it's taken as 100. 
(benchmark)`./bm_jacobi <k>`. Default iterations are 100, that can be changed [here](https://github.com/purusharths/hasc-proj/blob/master/code/jacobi-sgs/jacobi/bm_jacobi.cc#L76)


## Symmetric Gauss Sidel
For Symmetric Gauss Sidel
To compile: `make` inside `sgs` folder. `sgs`, `sgs_seq` and `bm_sgs` files are created.
`sgs` is the paralel version and `sgs_seq` is the sequential version. <br>
Command line args are: 
(main) `./ (sgs/sgs_seq) <grid size> <stencil size (k)> <iterations>` if iteration flag is not given, it's taken as 100. 
(bechmark) `./bm_sgs <k> <iterations>`

`make test` to create unit_test.


