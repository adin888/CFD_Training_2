#include"../Header/cfd_head.h"
#include"../Header/Jacobi.h"
/*
* - Using Galerkin method to solve 1d-Helmholtz equation
* - Using GLL points
* - Using CG(Conjugate gradient method) to solve system of equations
* - Using vectors to save variables 
*/
ArrayXd Test_f2(double lambda, ArrayXd x)
{
	ArrayXd y;

	y = (lambda + PI * PI) * sin(PI * x);
	return y;
}

ArrayXd Sol_u2(ArrayXd x)
{
	ArrayXd y;

	y = sin(PI * x);
	return y;
}

void AssembleMatrix2(int pNode, int np, ArrayXd& v)
{
	for (int i = 1; i < np; i++)
	{
		if ((i - (pNode + 1)) % pNode == 0)
		{
			v(i) += v(i + 1);
			v(i + 1) = v(i);
		}
	}
}


/* Conjugate gradient method */
ArrayXd CGHelmholtz1D2(int pNode, int ne, int np, int maxIter, double de, double lambda, double tolerance,
	double bcn, ArrayXXd massMatrix, ArrayXXd stiffnessMatrix, ArrayXd fh, ArrayXd uh)
{
	int i, j;

	double alpha;
	double beta;

	double rms;
	double initRms;

	ArrayXd r = ArrayXd::Zero(np - ne + 1);
	ArrayXd p = ArrayXd::Zero(np - ne + 1);
	ArrayXd Ap = ArrayXd::Zero(np - ne + 1);
	ArrayXd rOld = ArrayXd::Zero(np - ne + 1);

	VectorXd u = VectorXd::Zero(np - ne + 1);
	//VectorXd fe(fh);
	VectorXd f = VectorXd::Zero(np - ne + 1);

	MatrixXd Q = MatrixXd::Zero(np, np - ne + 1);          //distribution matrix
	MatrixXd Ae = MatrixXd::Zero(pNode + 1, pNode + 1);
	MatrixXd diagA = MatrixXd::Zero(np, np);
	MatrixXd A = MatrixXd::Zero(np - ne + 1, np - ne + 1);       //(lambda*M+L)

	// Defining the distribution matrix
	j = 0;
	for (i = 0; i < np; i++)
	{
		j = (int)(i / (pNode + 1));
		Q(i, i - j) = 1.0;
	}

	// Get the global coefficient matrix (lambda*M+L) 

	Ae = (lambda * de / 2.0 * massMatrix + 2.0 / de * stiffnessMatrix).matrix();
	for (i = 0; i < ne; i++)
	{
		diagA.block(i * (pNode + 1), i * (pNode + 1), pNode + 1, pNode + 1) = Ae;
	}
	A = Q.transpose() * diagA * Q;

	f = Q.transpose() * (VectorXd)fh;

	/* Calculate the initial residual */
	r = f - A * u;

	/* Dirichlet BC at x = -1 and Neumann-BC at x = 1 boundary */
	r(0) = 0;
	r(np - ne) += bcn;

	rms = pow(r, 2).sum() / (pNode + 1) / ne;
	initRms = rms;

	p = r;

	for (int iterNum = 1; iterNum < maxIter; iterNum++)
	{

		Ap = A * p.matrix();

		/* Dirichlet BC at x = -1 */
		Ap(0) = 0;

		alpha = ((r * r).sum()) /
			((Ap * p).sum() + 1e-16);

		u.array() += alpha * p;

		rOld = r;

		r -= alpha * Ap;
		//rms = pow((f - A * u).array(), 2).sum() / (pNode + 1) / ne;
		rms = pow(r, 2).sum() / (pNode + 1) / ne;

		cout << "The count of iteration is " << iterNum << endl;
		cout << "rms = " << setprecision(6) << rms / initRms << endl;

		if (rms / initRms <= tolerance)
			break;

		beta = ((r * r).sum()) /
			((rOld * rOld).sum() + 1e-16);

		p = r + beta * p;

	}
	r = f - A * u; //0 = -PI£¬ N = PI
	uh = Q * u;
	return uh;
}

/* Matrix-based Gauss¨CSeidel method  */

ArrayXd GSHelmholtz1D2(int pNode, int ne, int np, int maxIter, double de, double lambda, double tolerance,
	double bcn, ArrayXXd massMatrix, ArrayXXd stiffnessMatrix, ArrayXd fh, ArrayXd uh)
{
	int i, j;                     //global points

	double rms;
	double initRms;
	double dbc = uh(0, 0);									//Dirichlet BC

	VectorXd u = VectorXd::Zero(np - ne + 1);
	VectorXd f = VectorXd::Zero(np - ne + 1);
	VectorXd r = VectorXd::Zero(np - ne + 1);

	MatrixXd Q = MatrixXd::Zero(np, np - ne + 1);          //distribution matrix
	MatrixXd Ae = MatrixXd::Zero(pNode + 1, pNode + 1);
	MatrixXd diagA = MatrixXd::Zero(np, np);
	MatrixXd A = MatrixXd::Zero(np - ne + 1, np - ne + 1);       //(lambda*M+L)
	MatrixXd Linv = MatrixXd::Zero(np - ne, np - ne);     //left hand side of the expression for uh, lower truangular
	MatrixXd Us = MatrixXd::Zero(np - ne, np - ne);   //strictly upper triangular

	// Defining the distribution matrix
	j = 0;
	for (i = 0; i < np; i++)
	{
		j = (int)(i / (pNode + 1));
		Q(i, i - j) = 1.0;
	}

	// Get the global coefficient matrix (lambda*M+L) 

	Ae = (lambda * de / 2.0 * massMatrix + 2.0 / de * stiffnessMatrix).matrix();
	for (i = 0; i < ne; i++)
	{
		diagA.block(i * (pNode + 1), i * (pNode + 1), pNode + 1, pNode + 1) = Ae;
	}
	A = Q.transpose() * diagA * Q;

	// Get the global f
	f = Q.transpose() * (VectorXd)fh;

	// Calculate the inverse lower triangular and strictly upper triangular
	Linv = A.block(1, 1, np - ne, np - ne).triangularView<Lower>();
	if (Linv.determinant() != 0.0)
		Linv = Linv.inverse();
	else
		cout << "Error: the matrix L is not invertible." << endl;
	Us = A.block(1, 1, np - ne, np - ne).triangularView<StrictlyUpper>();

	// Dirichlet BC at x = -1 
	u(0) = dbc;
	f(np - ne) += bcn;
	r = f - A * u;
	r(0) = 0;


	rms = pow(r.array(), 2).sum() / (pNode + 1) / ne;
	initRms = rms;

	for (int iterNum = 1; iterNum < maxIter; iterNum++)
	{
		u.block(1, 0, np - ne, 1) = Linv * (f.block(1, 0, np - ne, 1) - Us * u.block(1, 0, np - ne, 1));
		r.block(1, 0, np - ne, 1) = f.block(1, 0, np - ne, 1) - A.block(1, 1, np - ne, np - ne) * u.block(1, 0, np - ne, 1);
		//r = f - A * u;
		// Dirichlet BC at x = -1 
		//r(0) = 0.0;

		rms = pow(r.array(), 2).sum() / (pNode + 1) / ne;

		cout << "The count of iteration is " << iterNum << endl;
		cout << "rms = " << setprecision(6) << rms / initRms << endl;

		if (rms / initRms <= tolerance)
			break;

	}

	uh = Q * u;
	return uh;
}

/* Element-based Gauss¨CSeidel method */

ArrayXd GSHelmholtz1D3(int pNode, int ne, int np, int maxIter, double de, double lambda, double tolerance,
	double bcn, ArrayXXd massMatrix, ArrayXXd stiffnessMatrix, ArrayXd fh, ArrayXd uh)
{
	int i, j;                     //global points

	double rms;
	double initRms;
	double dbc = uh(0, 0);									//Dirichlet BC
	double temp1, temp2;

	VectorXd u = VectorXd::Zero(np - ne + 1);
	VectorXd f = VectorXd::Zero(np - ne + 1);
	VectorXd r = VectorXd::Zero(np - ne + 1);

	MatrixXd Q = MatrixXd::Zero(np, np - ne + 1);          //distribution matrix
	MatrixXd Ae = MatrixXd::Zero(pNode + 1, pNode + 1);
	MatrixXd diagA = MatrixXd::Zero(np, np);
	MatrixXd A = MatrixXd::Zero(np - ne + 1, np - ne + 1);       //(lambda*M+L)

	// Defining the distribution matrix
	j = 0;
	for (i = 0; i < np; i++)
	{
		j = (int)(i / (pNode + 1));
		Q(i, i - j) = 1.0;
	}

	// Get the global coefficient matrix (lambda*M+L) 

	Ae = (lambda * de / 2.0 * massMatrix + 2.0 / de * stiffnessMatrix).matrix();
	for (i = 0; i < ne; i++)
	{
		diagA.block(i * (pNode + 1), i * (pNode + 1), pNode + 1, pNode + 1) = Ae;
	}
	A = Q.transpose() * diagA * Q;

	// Get the global f
	f = Q.transpose() * (VectorXd)fh;

	// Dirichlet BC at x = -1 
	u(0) = dbc;
	f(0) -= bcn;
	f(np - ne) += bcn;
	/*for (i = 1; i < (np - ne + 1); i++)
	{
		f(i) -= A(i, 0) * dbc;
	}*/
	r = f - A * u;
	r(0) = 0;

	rms = pow(r.array(), 2).sum() / (pNode + 1) / ne;
	initRms = rms;

	for (int iterNum = 1; iterNum < maxIter; iterNum++)
	{   
		for (i = 0; i < (np - ne + 1); i++)
		{
			temp1 = 0.0;
			temp2 = 0.0;
			for ( j = 0; j < i; j++)
			{
				temp1 += A(i, j) * u(j);
			}
			for (j = i + 1; j < (np - ne + 1); j++)
			{
				temp2 += A(i, j) * u(j);
			}
			u(i) = (f(i) - temp1 - temp2) / A(i, i);
		}
		
		//u.block(1, 0, np - ne, 1) = Linv * (f.block(1, 0, np - ne, 1) - Us * u.block(1, 0, np - ne, 1));
		//r.block(1, 0, np - ne, 1) = f.block(1, 0, np - ne, 1) - A.block(1, 1, np - ne, np - ne) * u.block(1, 0, np - ne, 1);
		r = f - A * u;
		// Dirichlet BC at x = -1 
		u(0) = dbc;
		r(0) = 0.0;

		rms = pow(r.array(), 2).sum() / (pNode + 1) / ne;

		cout << "The count of iteration is " << iterNum << endl;
		cout << "rms = " << setprecision(6) << rms / initRms << endl;

		if (rms / initRms <= tolerance)
			break;
	}

	uh = Q * u;
	return uh;
}

void HOM1DHelmholtz2()
{
	int pNode = 7;			 //polynomial degree
	int nx = 200;            //number of plotting points per element
	int ne = 20;              //number of elements
	int np = (pNode + 1) * ne;
	int i, j, q;

	double dx = 2.0 / nx;
	double de = 2.0 / ne;    //width of element
	double lambda = 1.0;
	double bcn = -PI;        //Neumann BC
	double errMax, errL2;
	double maxIter = 1e4;

	ArrayXd xi(pNode + 1);                       //locations of the nodal points inside the element
	ArrayXd wi = ArrayXd::Zero(pNode + 1);       //the integration weights
	ArrayXXd diff = ArrayXd::Zero(pNode + 1);    //the differentiation matrix

	ArrayXd xh(np);                  //the node locations
	ArrayXd fh(np);
	ArrayXd uex(np);                 //exact value of function
	ArrayXd uh = ArrayXd::Zero(np);

	ArrayXd massMatrixDiag(pNode + 1);
	ArrayXXd massMatrix = ArrayXXd::Zero(pNode + 1, pNode + 1);
	ArrayXXd stiffnessMatrix = ArrayXXd::Zero(pNode + 1, pNode + 1);

	Jacobi J;

	/* Get the GLL points, weights, diffMatrix in standard element */

	xi = J.GLLPoints(pNode);
	wi = J.GLLWeights(xi);
	diff = J.GLLDiffMatrix(xi);

	/* Obtain the node locations in global system */

	for (i = 0; i < np; i++)
	{
		if (i % (pNode + 1) == 0)
			j = 0;
		xh(i) = -1.0 + de * int(i / (pNode + 1)) + (xi(j) + 1.0) / 2.0 * de;
		j++;
	}

	/* Calculate the system matrix */

	massMatrixDiag = wi;
	for (i = 0; i < pNode + 1; i++)
	{
		massMatrix(i, i) = massMatrixDiag(i);
	}

	for (i = 0; i < pNode + 1; i++)
	{
		for (j = 0; j < pNode + 1; j++)
		{
			for (q = 0; q < pNode + 1; q++)
			{
				stiffnessMatrix(i, j) += wi(q) * diff(q, i) * diff(q, j);
			}
		}
	}

	/* Calculate F in the system (lambda*M + L)*u = F */

	fh = Test_f2(lambda, xh);

	for (i = 0; i < np; i++)
	{
		if (i % (pNode + 1) == 0)
			j = 0;
		fh(i) *= de / 2.0 * massMatrixDiag(j);
		j++;
	}

	/* Set the Dirichlet-BC at the left boundary */

	uh(0) = Sol_u2(xh)(0);
	//fh(np - 1) += bcn;

	//uh = CGHelmholtz1D2(pNode, ne, np, maxIter, de, lambda, 1e-5,
	//	bcn, massMatrix, stiffnessMatrix, fh, uh);

	//uh = GSHelmholtz1D2(pNode, ne, np, maxIter, de, lambda, 1e-5,
	//	bcn, massMatrix, stiffnessMatrix, fh, uh);

	uh = GSHelmholtz1D3(pNode, ne, np, maxIter, de, lambda, 1e-5,
		bcn, massMatrix, stiffnessMatrix, fh, uh);

	/* calculate the exact result and error */

	uex = Sol_u2(xh);
	errMax = (uex - uh).abs().maxCoeff();
	errL2 = sqrt(pow(uex - uh, 2).sum() / (pNode + 1) / ne);

	/* write the result */
	ofstream dataOut;

	dataOut.open("u_final.csv", ios::out | ios::trunc);

	dataOut << setprecision(16) << "errMax" << "," << errMax << "," << "errL2" << "," << errL2 << endl;
	dataOut << "x" << "," << "uex" << "," << "uh" << endl;
	for (i = 0; i < np; i++)
	{
		dataOut << setprecision(16) << xh(i) << "," << uex(i) << "," << uh(i) << endl;
	}

	dataOut.close();
}