#include"../Header/cfd_head.h"
#include"../Header/Jacobi.h"

/*
* -This is about the basic use of the Galerkin method in SEM£¨spectral method£©/HOM£¨high order method£©. 
* -It includes approximating the original function and finding integrals and differentiators.
* -Approximation of different functions using GL points or GLL points.
*/
ArrayXd Function(ArrayXd x)
{
	ArrayXd y;

	//y = sin(PI * x);
	y = 1.0 / (1.0 + 25.0 * x * x);
	//y = abs(x);

	return y;
}

double IFunction()
{
	double y;

	//y = 0.0;
	y = 2.0 / 5.0 * atan(5.0);
	//y = 1.0;

	return y;
}

ArrayXd DiffFunction(ArrayXd x)
{
	ArrayXd y;

	//y = PI * cos(PI * x);
	y = -1.0 / (pow(1.0 + 25.0 * x * x, 2.0)) * (50.0 * x);
	//y = x.cwiseSign();

	return y;
}
/* Fitting functions using higher order methods(spectral methods) */

void FittingFunction()
{
	int pNode = 15;
	int nx = 500;
	int poly = 0;            // 0=GL, 1=GLL

	double dx = 2.0 / nx;
	double integral = 0.0;
	double integralExact;

	ArrayXd xi(pNode + 1);                       //locations of the nodal points inside the element
	ArrayXd wi = ArrayXd::Zero(pNode + 1);      //the integration weights
	ArrayXd ui = ArrayXd::Zero(pNode + 1);      //the coeff. of the basisfunction
	ArrayXd dui = ArrayXd::Zero(pNode + 1);     //the coeff. of the derivative of BF
	MatrixXd diff = ArrayXd::Zero(pNode + 1);    //the differentiation matrix

	ArrayXd x(nx + 1);
	ArrayXd ue(nx + 1), due(nx + 1);            //exact value of function
	ArrayXd up = ArrayXd::Zero(nx + 1), dup = ArrayXd::Zero(nx + 1); //approximate value of function

	Jacobi J;

	for (int i = 0; i < nx + 1; i++)
	{
		x(i) = -1.0 + dx * i;
	}

	/* Obtain the points and the integration and differentiation parameters */
	switch (poly)
	{
	case 0:
	{
		xi = J.GLPoints(pNode);
		wi = J.GLWeights(xi);
		diff = J.GLDiffMatrix(xi);
		break;
	}
	case 1:
	{
		xi = J.GLLPoints(pNode);
		wi = J.GLLWeights(xi);
		diff = J.GLLDiffMatrix(xi);
		break;
	}
	default:
		break;
	}

	/* Obtain the coefficients of polynomials */
	ui = Function(xi);

	/* Calculate and compare integrals */
	integralExact = IFunction();

	for (int i = 0; i <= pNode; i++)
	{
		integral += wi(i) * ui(i);
	}

	cout << setprecision(16) << "The numerical integral is " << integral << endl;
	cout << setprecision(16) << "The exact integral is " << integralExact << endl;
	cout << setprecision(16) << "The Error of integration is " << abs(integral - integralExact) << endl;

	/* Calculate the derivative */
	for (int i = 0; i <= pNode; i++)
	{
		for (int j = 0; j <= pNode; j++)
		{
			dui(i) += diff(i, j) * ui(j);
		}
	}

	/* Approximation of the function */
	switch (poly)
	{
	case 0:
	{
		for (int i = 0; i < nx + 1; i++)
		{
			for (int j = 0; j < pNode + 1; j++)
			{
				up(i) += J.GLPolynomial(j, x(i), xi) * ui(j);
				dup(i) += J.GLPolynomial(j, x(i), xi) * dui(j);
			}
		}
		break;
	}
	case 1:
	{
		for (int i = 0; i < nx + 1; i++)
		{
			for (int j = 0; j < pNode + 1; j++)
			{
				up(i) += J.GLLPolynomial(j, x(i), xi) * ui(j);
				dup(i) += J.GLLPolynomial(j, x(i), xi) * dui(j);
			}
		}
		break;
	}
	default:
		break;
	}

	/* Exact Value of the function */
	ue = Function(x);
	due = DiffFunction(x);

	/* write the result */
	ofstream dataOut;

	dataOut.open("test_data.csv", ios::out | ios::trunc);

	dataOut << "x" << "," << "ue" << "," << "due" << "," << "up" << "," << "dup" << endl;
	for (int i = 0; i < nx + 1; i++)
	{
		dataOut << setprecision(16) << x(i) << ","<<ue(i) << "," 
			<<due(i) << ","<<up(i) << ","<<dup(i)<<endl;
	}

	dataOut.close();
}