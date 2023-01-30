#include"cfd_head.h"
#include"Jacobi.h"

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

ArrayXXd CGHelmholtz1D(int pNode, int ne, int maxIter, double de, double lambda, double tolerance,
	ArrayXXd massMatrix, ArrayXXd stiffnessMatrix, ArrayXXd fh, ArrayXXd uh)
{
	int i; 

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

	for (int iterNum = 1; iterNum < maxIter*10; iterNum++)
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
		//
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
	double maxIter = 1000;

	ArrayXd xi(pNode + 1);                       //locations of the nodal points inside the element
	ArrayXd wi = ArrayXd::Zero(pNode + 1);       //the integration weights
	ArrayXXd diff = ArrayXd::Zero(pNode + 1);    //the differentiation matrix

	ArrayXXd xh(pNode + 1, ne);                  //the node locations
	ArrayXXd fh(pNode + 1, ne);
	ArrayXXd uex(pNode + 1, ne);                 //exact value of function
	ArrayXXd uh = ArrayXXd::Zero(pNode + 1, ne);

	ArrayXd xiPlot(pNode + 1);
	ArrayXXd xPlot(nx + 1, ne);
	ArrayXXd uhPlote(nx + 1, ne);            
	ArrayXXd uexPlot(nx + 1, ne);         

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

	uh = CGHelmholtz1D(pNode, ne, maxIter, de, lambda, 1e-4,
		massMatrix, stiffnessMatrix, fh, uh);

	/* calculate the exact result and error */

	uex = Sol_u(xh);
	errMax = (uex - uh).abs().maxCoeff();
	errL2 = ComputeL2norm(pNode + 1, ne, uex - uh);

	/* write the result */
	ofstream dataOut;

	dataOut.open("test_data.csv", ios::out | ios::trunc);

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