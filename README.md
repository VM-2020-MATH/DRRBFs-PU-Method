There are two scripts in MATLAB that implement the DRRBF-PU (Direct Rational RBF Partition of Unity) method, for estimating the derivatives.

These scripts are

"D_Rational_RBFs_PU_2D.m" and 

"D_Rational_RBFs_PU_2D.m". 

To reproduce the results,

the user can directly run

"D_Rational_RBFs_PU_2D.m" (see also the comments type proposed in the file).

Authors: Vahid Mohammadi 1, Stefano De Marchi 2

1. Department of Mathematics, Faculty of Science, Shahid Rajaee Teacher Training University, Tehran, 16785-163, Iran

2. Department of Mathematics "Tullio Levi-Civita",
   University of Padua, Italy

   
   Note that:
   
a. To construct the local approximations in each patch a Matérn radial kernel is utilized.

b. A compactly supported Wendland's function is used as the PU weights.

c. The user can change the radial kernel and PU weights (it also is possible to apply a discontinuous weight function).

Reference Paper:

A Note on the Direct Approximation of Derivatives in Rational Radial Basis Functions Partition of Unity Method, Submitted.
