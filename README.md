There are two scripts in MATLAB that implement the D-RRBF-PU (Direct Rational RBF Partition of Unity) method, for estimating the derivatives.

These scripts are:

"D_Rational_RBFs_PU_2D.m" 

and the function 

"Weight_Rational_RBFs_PU_2D.m". 

To reproduce the results,

the user can directly run

"D_Rational_RBFs_PU_2D.m" (see also the comments type proposed in the file).

Also, two other scripts in MATLAB are given showing the implementation of the D-RRBF-PU method combined with a fourth-order Runge-Kutta method for solving the Convection-Diffusion equation in two dimensions (see Subsection 5.1, A rotating Gaussian pulse).

These scripts are:

"Convection_Diffusion_Ex1.m" 

and the function

"Weight_Rational_RBF_PU_2D_CD.m"

Authors: Vahid Mohammadi 1, Stefano De Marchi 2

1. Department of Mathematics, Faculty of Science, Shahid Rajaee Teacher Training University, Tehran, 16785-163, Iran

2. Department of Medicine and Padova Neuroscience Center, University of Padua, Italy

   
   Note that:
   
a. To construct the local approximations in each patch a Matérn radial kernel is utilized.

b. A compactly supported Wendland's function is used as the PU weights.

c. The user can change the radial kernel and PU weights (it also is possible to apply a discontinuous weight function).

d. If the user wants to reproduce the second example of Convection-diffusion equation (See Subsection 5.2, the inviscid Burgers equation), he/she needs to change accordingly. 

Reference Paper:

A Note on the Direct Approximation of Derivatives in Rational Radial Basis Functions Partition of Unity Method and Its
Application to the Convection-Diffusion Equations, Submitted.
