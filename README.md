# CFD_Training_2
Programming exercises for different methods.  
The purpose of this repository is to give ideas for building a complete numerical simulation program. The purpose of this repository is to give ideas for building a complete numerical simulation program. Starting from the FDM, each example uses a different method to solve the corresponding partial differential equation. The source code is stored in the path CFD\_Training\_2/ CFD\_FDM_2. The spectral method(high order method) in FEM was added later.   
The methods used in each source file and the partial differential equations solved are given below:  

-  FDM/ HLLC\_Riemann_Solver.cpp  
	-  Solving the 1d Euler equation
	-  Using WENO-5 Scheme to reconstruct the left and right side fluxes at the interface
	-  Using Runge-Kutta-3 Scheme for time integration
	-  Using HLLC scheme to approximate Riemann solver
	-  ![](https://github.com/adin888/CFD_Training_2/blob/main/CFD_FDM_2/FDM/HLLC_Riemann_Solver.png)
-  FDM/ Roe\_Riemann_Solver.cpp
	-  Solving the 1d Euler equation
	-  Using WENO-5 Scheme to reconstruct the left and right side fluxes at the interface
	-  Using Runge-Kutta-3 Scheme for time integration
	-  Using Reo approximate Riemann solver 
	-  ![](https://github.com/adin888/CFD_Training_2/blob/main/CFD_FDM_2/FDM/Roe_Riemann_Solver.png)
-  FDM/ Rusanov\_Riemann_Solver.cpp
	-  Solving the 1d Euler equation
	-  Using WENO-5 Scheme to reconstruct the left and right side fluxes at the interface
	-  Using Runge-Kutta-3 Scheme for time integration
	-  Using Rusanov scheme to approximate Riemann solver
	-  ![](https://github.com/adin888/CFD_Training_2/blob/main/CFD_FDM_2/FDM/Rusanov_Riemann_Solver.png)
-  FDM/ FFT.cpp
	-  Using Fast Fourier Transform to slove Poisson equation for the periodic domain
	-  Using additional library FFTW3
	-  ![](https://github.com/adin888/CFD_Training_2/blob/main/CFD_FDM_2/FDM/FFT.png)
-  FDM/ FST.cpp
	-  -Using Fast Sine Transform to slove Poisson equation for the periodic domain
	-  ![](https://github.com/adin888/CFD_Training_2/blob/main/CFD_FDM_2/FDM/FST.png)
-  FDM/ Gauss\_seidel.cpp
	-  Using Gauss-Seidel iterative methods to slove Poisson equation with Dirichlet boundary condition
	-  Using the ratio of L2 to its initial value during the iteration as the criterion for convergence
	-  ![](https://github.com/adin888/CFD_Training_2/blob/main/CFD_FDM_2/FDM/Gauss_seidel.png)
-  FDM/ Conjugate\_Gradient.cpp
	-  Using Conjugate Gradient methods to slove Poisson equation with Dirichlet boundary condition
	-  Using the ratio of L2 to its initial value during the iteration as the criterion for convergence
	-  ![](https://github.com/adin888/CFD_Training_2/blob/main/CFD_FDM_2/FDM/Conjugate_Gradient.png)
-  FDM/ Conjugate\_Gradient.cpp
	-  viscous incompressible flow
	-  a square cavity consisting of three rigid walls with no-slip conditions and a lid moving with a tangential unit velocity
	-  Using Runge-Kutta-3 Scheme for time integration
	-  Using second-order Arakawa scheme for the nonlinear terms
	-  ![](https://github.com/adin888/CFD_Training_2/blob/main/CFD_FDM_2/FDM/Lid_Driven_Cavity_w.png)
	-  ![](https://github.com/adin888/CFD_Training_2/blob/main/CFD_FDM_2/FDM/Lid_Driven_Cavity_psi.png)
-   HOM(SEM)/ HOM1D\_Helmholtz_v2
	-   Using Galerkin method to solve 1d-Helmholtz equation
	-   Using Gauss–Lobatto–Legendre points to a
	-   Using Conjugate gradient(CG) method oder Gauss-Seidel(GS) method to solve system of equations
	-   ![](https://github.com/adin888/CFD_Training_2/blob/main/CFD_FDM_2/HOM(SEM)/HOM1D_HelmholtzCG.png)
	-   ![](https://github.com/adin888/CFD_Training_2/blob/main/CFD_FDM_2/HOM(SEM)/HOM1D_HelmholtzGS.png)  
  
###Note：
You need to configure the Eigen library before using the code in this program. 
