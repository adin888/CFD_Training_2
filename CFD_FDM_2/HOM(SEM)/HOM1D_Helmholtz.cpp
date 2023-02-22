#include"../Header/cfd_head.h"
#include"../Header/Jacobi.h"
/*
* -Using Galerkin method to solve 1d-Helmholtz equation
* -Using GLL points
* -Using CG(Conjugate gradient method) to solve system of equations
*/
ArrayXXd Test_f(double lambda, ArrayXXd x)
{
	ArrayXXd y;

	y = (lambda + PI * PI) * sin(PI * x);
	return y;
}

ArrayXXd Sol_u(ArrayXXd x)
{
	ArrayXXd y;

	y = sin(PI * x);
	return y;
}

void AssembleMatrix(int pNode, int ne, ArrayXXd &m)
{
	for (int i = 1; i < ne; i++)
	{
		m(0, i) += m(pNode, i - 1);
		m(pNode, i - 1) = m(0, i);
	}
}

void AssembleMatrix(int pNode, int ne, MatrixXd& m)
{
	for (int i = 1; i < ne; i++)
	{
		m(0, i) += m(pNode, i - 1);
		m(pNode, i - 1) = m(0, i);
	}
}

/* Conjugate gradient method */
ArrayXXd CGHelmholtz1D(int pNode, int ne, int maxIter, double de, double lambda, double tolerance,
	ArrayXXd massMatrix, ArrayXXd stiffnessMatrix, ArrayXXd fh, ArrayXXd uh)
{
	double alpha;
	double beta;

	double rms;
	double initRms;

	ArrayXXd r = ArrayXXd::Zero(pNode + 1, ne);
	ArrayXXd p = ArrayXXd::Zero(pNode + 1, ne);
	ArrayXXd Ap = ArrayXXd::Zero(pNode + 1, ne);
	ArrayXXd rOld = ArrayXXd::Zero(pNode + 1, ne);

	/* Calculate the initial residual */
	r = fh.matrix() - 
		(lambda * de / 2.0 * massMatrix + 2.0 / de * stiffnessMatrix).matrix() * uh.matrix();
	/* Assemble the elementwise matrix in global matrix */
	AssembleMatrix(pNode, ne, r);

	/* Dirichlet BC at x = -1 */
	r(0, 0) = 0;

	rms = pow(r, 2).sum() / (pNode + 1) / ne;
	initRms = rms;

	p = r;

	for (int iterNum = 1; iterNum < maxIter; iterNum++)
	{

		Ap = (lambda * de / 2.0 * massMatrix + 2.0 / de * stiffnessMatrix).matrix() * p.matrix();
		AssembleMatrix(pNode, ne, Ap);

		/* Dirichlet BC at x = -1 */
		Ap(0, 0) = 0;

		alpha = ((r * r).sum()) /
			((Ap * p).sum() + 1e-16);

		uh += alpha * p;

		rOld = r;

		r -= alpha * Ap;

		rms = pow(r, 2).sum() / (pNode + 1) / ne;

		cout << "The count of iteration is " << iterNum << endl;
		cout << "rms = " << setprecision(6) << rms / initRms << endl;

		if (rms / initRms <= tolerance)
			break;

		beta = ((r * r).sum()) /
			((rOld * rOld).sum() + 1e-16);

		p = r + beta * p;

	}
	return uh;
}

/* Gauss¨CSeidel method (Not applicable to current variable writing!) */

ArrayXXd GSHelmholtz1D(int pNode, int ne, int maxIter, double de, double lambda, double tolerance,
	ArrayXXd massMatrix, ArrayXXd stiffnessMatrix, ArrayXXd fh, ArrayXXd uh)
{
	int i, j;
	int gNx = (pNode + 1) * ne;                          //global points

	double rms;
	double initRms;
	double dbc = uh(0, 0);									//Dirichlet BC

	VectorXd ue = VectorXd::Zero(gNx);
	VectorXd fe = VectorXd::Zero(gNx);
	VectorXd u = VectorXd::Zero(gNx - ne + 1);
	VectorXd f = VectorXd::Zero(gNx - ne + 1);
	VectorXd r = VectorXd::Zero(gNx - ne + 1);

	MatrixXd Q = MatrixXd::Zero(gNx, gNx - ne + 1);          //distribution matrix
	MatrixXd Ae = MatrixXd::Zero(pNode + 1, pNode + 1);
	MatrixXd diagA = MatrixXd::Zero(gNx, gNx);
	MatrixXd A = MatrixXd::Zero(gNx - ne + 1, gNx - ne + 1);       //(lambda*M+L)
	MatrixXd Linv = MatrixXd::Zero(gNx - ne + 1, gNx - ne + 1);     //left hand side of the expression for uh, lower truangular
	MatrixXd Us = MatrixXd::Zero(gNx - ne + 1, gNx - ne + 1);   //strictly upper triangular

	// Defining the distribution matrix
	j = 0;
	for (i = 0; i < gNx; i++)
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
	//cout << diagA << endl;
	//cout << A << endl;
	// Calculate the inverse lower triangular and strictly upper triangular
	Linv = A.triangularView<Lower>();
	if (Linv.determinant() != 0.0)
		Linv = Linv.inverse();
	else
		cout << "Error: the matrix L is not invertible." << endl;
	Us = A.triangularView<StrictlyUpper>();
	// Get the global f
	//AssembleMatrix(pNode, ne, fh);
	for (i = 0; i < gNx; i++)
	{
		if (i % (pNode + 1) == 0)
			j = 0;
		fe(i) = fh(j, (int)(i / (pNode + 1)));
		j++;
	}
	//cout << fh << endl;
	//cout << fe << endl;

	f = Q.transpose() * fe;
	//cout<<f  << endl;
	r = f - A * u;
	// Dirichlet BC at x = -1 
	u(0) = dbc;
	r(0) = 0;

	rms = pow(r.array(), 2).sum() / (pNode + 1) / ne;
	initRms = rms;

	for (int iterNum = 1; iterNum < maxIter; iterNum++)
	{
		u = Linv * (f - Us * u);
		r = f - A * u;

		// Dirichlet BC at x = -1 
		//u(0) = dbc;
		r(0) = 0.0;

		rms = pow(r.array(), 2).sum() / (pNode + 1) / ne;

		cout << "The count of iteration is " << iterNum << endl;
		cout << "rms = " << setprecision(6) << rms / initRms << endl;

		if (rms / initRms <= tolerance)
			break;

	}

	ue = Q * u;
	for (i = 0; i < gNx; i++)
	{
		if (i % (pNode + 1) == 0)
			j = 0;
		uh(j, (int)(i / (pNode + 1))) = ue(i);
		j++;
	}

	return uh;
}


void HOM1DHelmholtz()
{
	int pNode = 7;			 //polynomial degree
	int nx = 200;            //number of plotting points per element
	int ne = 20;              //number of elements
	int i, j, e, q;

	double dx = 2.0 / nx;    
	double de = 2.0 / ne;    //width of element
	double lambda = 1.0;
	double bcn = -PI;        //Neumann BC
	double errMax, errL2;
	double maxIter = 1e4;

	ArrayXd xi(pNode + 1);                       //locations of the nodal points inside the element
	ArrayXd wi = ArrayXd::Zero(pNode + 1);       //the integration weights
	ArrayXXd diff = ArrayXd::Zero(pNode + 1);    //the differentiation matrix

	ArrayXXd xh(pNode + 1, ne);                  //the node locations
	ArrayXXd fh(pNode + 1, ne);
	ArrayXXd uex(pNode + 1, ne);                 //exact value of function
	ArrayXXd uh = ArrayXXd::Zero(pNode + 1, ne);   

	ArrayXd massMatrixDiag(pNode + 1);
	ArrayXXd massMatrix = ArrayXXd::Zero(pNode + 1, pNode + 1);
	ArrayXXd stiffnessMatrix = ArrayXXd::Zero(pNode + 1, pNode + 1);

	Jacobi J;

	/* Get the GLL points, weights, diffMatrix in standard element */

	xi = J.GLLPoints(pNode);
	wi = J.GLLWeights(xi);
	diff = J.GLLDiffMatrix(xi);

	/* Obtain the node locations in global system */

	for (e = 0; e < ne; e++)
	{
		xh.col(e) = -1.0 + de * e + (xi + 1.0) / 2.0 * de;
	}
	
	/* Calculate the system matrix */

	massMatrixDiag = wi;
	for ( i = 0; i < pNode + 1; i++)
	{
		massMatrix(i, i) = massMatrixDiag(i);
	}

	for (i = 0; i < pNode + 1; i++)
	{
		for (j = 0; j < pNode + 1; j++)
		{
			for ( q = 0; q < pNode + 1; q++)
			{
				stiffnessMatrix(i, j) += wi(q) * diff(q, i) * diff(q, j);
			}
		}
	}

	/* Calculate F in the system (lambda*M + L)*u = F */

	fh = Test_f(lambda, xh);

	for (e = 0; e < ne; e++)
	{
		fh.col(e) *= de/2.0 * massMatrixDiag;
	}

	/* Set the Dirichlet-BC at the left boundary and Neumann-BC at the right boundary */

	uh(0, 0) = Sol_u(xh)(0, 0);
	fh(pNode, ne - 1) += bcn;

	uh = CGHelmholtz1D(pNode, ne, maxIter, de, lambda, 1e-8,
		massMatrix, stiffnessMatrix, fh, uh);

	//uh = GSHelmholtz1D(pNode, ne, maxIter, de, lambda, 1e-4,
	//	massMatrix, stiffnessMatrix, fh, uh);

	/* calculate the exact result and error */

	uex = Sol_u(xh);
	errMax = (uex - uh).abs().maxCoeff();
	errL2 = ComputeL2norm(pNode + 1, ne, uex - uh);

	/* write the result */
	ofstream dataOut;

	dataOut.open("u_final.csv", ios::out | ios::trunc);

	dataOut << setprecision(16) << "errMax" << "," << errMax << "," << "errL2" << "," << errL2 << endl;
	dataOut << "x" << "," << "uex" << "," << "uh" << endl;
	for (e = 0; e < ne; e++)
	{
		for ( i = 0; i < pNode + 1; i++)
		{
			dataOut << setprecision(16) << xh(i, e) << "," << uex(i, e) << ","<< uh(i, e) << endl;
		}

	}

	dataOut.close();
}