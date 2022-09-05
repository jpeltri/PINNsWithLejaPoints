# PINNsWithLejaPoints

Some numerical experiments on solving the heat equation with a physics-informed neural network with a focus on using quadrature based on leja-points.
This serves as supplementary material to my masters thesis.

The jupyter notebook contains all the actual experiments. The matlab code is for constructing sparse grids and uses the Sparse Grids Matlab Kit (available at https://sites.google.com/view/sparse-grids-kit).
Note that the API used to interface with these bits of code has a relatively narrow window of compatible Matlab and Python versions, see https://www.mathworks.com/content/dam/mathworks/mathworks-dot-com/support/sysreq/files/python-compatibility.pdf for the full list.
