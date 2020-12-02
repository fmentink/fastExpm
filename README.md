# fastExpm.m 

This function efficiently implements matrix exponential matrix for sparse and full matrices. 
This code is based on scaling, taylor series and scaling.
Works for both CPU and GPU (with appropriate toolbox). 
Use: 
- P = fastExpm(H) --> default settings
- P = fastExpm(H,convergenceCriteria)
- P = fastExpm(H,convergenceCriteria, nonZeroTol)

Two criteria are used to speed the computation and preserve the sparsity.
- [1] convergenceCriteria determine the threshold for the Taylor series (default 1e-6)
- [2] nonZeroTol strips elements smaller than nonZeroTol at each computation step to preserve sparsity, default is 1e-14.
The code automatically switches from sparse to full if sparsity is below 15 % to maintain speed.
If H is a gpuArray, P will be computed on the GPU.

This code has been originally developed by Ilya Kuprov (http://spindynamics.org/) and has been adapted by F. Mentink-Vigier (fmentink@magnet.fsu.edu).
If you use this code, please cite 
- H. J. Hogben, M. Krzystyniak, G. T. P. Charnock, P. J. Hore and I. Kuprov, Spinach – A software library for simulation of spin dynamics in large spin systems, J. Magn. Reson., 2011, 208, 179–194.
- I. Kuprov, Diagonalization-free implementation of spin relaxation theory for large spin systems., J. Magn. Reson., 2011, 209, 31–38.
