
Uses guide for JacobiG_nonsym:

The codes are used to tensor diagonalization for nonsymmetric (real or complex) tensors. 

Please refer to the following paper:

@article{sheng2021Jacobi,
author={Z. Sheng and J. Li and Q. Ni}, 
title={Jacobi-type algorithms for homogeneous polynomial optimization on Stiefel manifolds with applications to tensor approximations}, 
journal={arXiv:2110.01886},
year={2021}
}

Uses guide for JacobiG:

The codes are used to tensor diagonalization for symmetric (complex) tensors,
and are downloaded from https://github.com/kdu/jacobi-G-unitary-matlab.

The main contributors to this code are Usevich, Li and Comon.

Please refer to their paper as follows:

@article{Usevich2020,
   title={Approximate matrix and tensor diagonalization by unitary transformations: convergence of {J}acobi-type algorithms},
   author={Usevich, K. and Li, J. and Comon, P.},
   journal={SIAM Journal on Optimization},
   volume={30},
   number={4},
   pages={2998--3028},
   year={2020},
}

In fact, some codes of JacobiG are modified based on jacobi-G-unitary-matlab.

The main functions are JacobiG_nonsym_complex.m, JacobiG_nonsym_complex_r.m, JacobiG_2.m and JacobiG_2P.m. 

Please see the examplary files Jacobi_test_example.m and main_test.m for how to use them. 
