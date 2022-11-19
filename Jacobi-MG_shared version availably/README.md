
Uses guide for JacobiG_nonsym: Written by sz and jz.

The codes are used to tensor diagonalization for nonsymmetric (real or complex) tensors. 

Please refer to the following paper:

@article{sheng2021Jacobi,
author={Z. Sheng and J. Li and Q. Ni}, 
title={Jacobi-type algorithms for homogeneous polynomial optimization on Stiefel manifolds with applications to tensor approximations}, 
journal={arXiv:2110.01886},
year={2021}
}


Installation
 
Download the zip file Jacobi-MG_shared_version_availably.zip and extract it to a directory where you like to save. 

Make sure Matlab paths are added correctly. 
To compare with some Riemannian algorithms, you need to install the Manopt 
  - Manopt_6.0 https://www.manopt.org/downloads.html 


How to use the code? 

The main functions are JacobiG_nonsym_complex.m and JacobiG_nonsym_complex_r.m. 

Please see the examplary files Jacobi_test_example.m and test_tensor.m for how to use it. 
