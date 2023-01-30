#include "Jacobi.h"

#pragma region Jacobi Polynomial
double Jacobi::Polynomial(double alpha, double beta, int nOrder, double x)
{
	double a1, a2, a3, a4;            //Intermediate variables
	double P0, P1, Pn;                //the value of Jacobi polynomial P^{a,b}_{n}(x)

	P0 = 1.0;
	P1 = (alpha - beta + (alpha + beta + 2.0) * x) / 2.0;

	if (nOrder < 2)
		Pn = (1.0 - nOrder) * P0 + nOrder * P1;
	else if (alpha > -1.0 && beta > -1.0)
	{
		for (int i = 1; i < nOrder; i++)
		{
			a1 = 2.0 * (i + 1.0) * (i + alpha + beta + 1.0) * (2.0 * i + alpha + beta);
			a2 = (2.0 * i + alpha + beta + 1.0) * (alpha * alpha - beta * beta);
			a3 = (2.0 * i + alpha + beta) * (2.0 * i + alpha + beta + 1.0) * (2.0 * i + alpha + beta + 2.0);
			a4 = 2.0 * (i + alpha) * (i + beta) * (2.0 * i + alpha + beta + 2.0);
			Pn = ((a2 + a3 * x) * P1 - a4 * P0) / a1;
			P0 = P1;
			P1 = Pn;
		}
	}
	else
	{
		Pn = 1e10;
		cout << "erro: alpha or beta in JacobiPolynomial is out of range" << endl;
	}

	return Pn;
}
double Jacobi::PolynomialDerivative(double alpha, double beta, int nOrder, double x)
{
	double dPn;  //first derivative of P^ { a,b }_{ n }

	return dPn = 0.5 * (alpha + beta + nOrder + 1.0)
		* Polynomial(alpha + 1.0, beta + 1.0, nOrder - 1.0, x);
}
ArrayXd Jacobi::PolynomialZeros(double alpha, double beta, int nOrder)
{
	int maxIter = 100;
	double r, s, d;
	double Pn, dPn;
	ArrayXd xZeros(nOrder);      //Zeros of Jacobi polynomials with n Order
	ArrayXd tem(nOrder);

	for (int i = 0; i < nOrder; i++)
	{
		r = -(double)(cos(PI * (2.0 * i + 1.0) / (2.0 * nOrder)));
		if (i > 0)
			r = 0.5 * (r + xZeros(i - 1));
		for (int j = 1; j <= maxIter; j++)
		{
			if (i > 0)
			{
				s = 0.0;
				for (int k = 0; k < nOrder; k++)
					s += 1.0 / (r - xZeros(k));
			}
			else
				s = 1.0;
			Pn = Polynomial(alpha, beta, nOrder, r);
			dPn = PolynomialDerivative(alpha, beta, nOrder, r);
			d = Pn / (dPn - s * Pn);
			r = r - d;
			if (abs(d) < eps)
				break;
		}
		xZeros(i) = r;
	}
	return xZeros;
}

#pragma endregion

#pragma region Gauss-Legendre Lagrangre quadrature points and polynomials
/* p = ID of polynomial */
double Jacobi::GLPolynomial(int p, double x, ArrayXd xi)
{
	int pOrder = xi.size() - 1;
	double Pn;

	if (abs(x - xi(p)) < eps)
		Pn = 1.0;
	else
		Pn = Polynomial(0.0, 0.0, pOrder + 1, x)
		/ (PolynomialDerivative(0.0, 0.0, pOrder + 1, xi(p)) * (x - xi(p)));
	return Pn;
}
/* n = ID of point */
double Jacobi::GLPolynomialDerivative(int p, int n, ArrayXd xi)
{
	int pOrder = xi.size() - 1;
	double dPn;

	if (p == n)
		dPn = xi(n) / (1.0 - xi(n) * xi(n));
	else
		dPn = PolynomialDerivative(0.0, 0.0, pOrder + 1, xi(n))
		/ (PolynomialDerivative(0.0, 0.0, pOrder + 1, xi(p)) * (xi(n) - xi(p)));
	return dPn;
}

ArrayXd Jacobi::GLPoints(int pOrder)
{
	ArrayXd x(pOrder + 1);           //Gauss-Legendre points

	x = PolynomialZeros(0.0, 0.0, pOrder + 1);
	return x;
}

ArrayXd Jacobi::GLWeights(ArrayXd x)
{
	ArrayXd w(x.size());

	int pOrder = x.size() - 1;
	for (int i = 0; i <= pOrder; i++)
	{
		w(i) = 2.0 / ((1.0 - x(i) * x(i))
			* pow(PolynomialDerivative(0.0, 0.0, pOrder + 1, x(i)), 2));
	}
	return w;
}

MatrixXd Jacobi::GLDiffMatrix(ArrayXd xi)
{
	int pOrder = xi.size() - 1;
	MatrixXd diff(pOrder + 1, pOrder + 1);

	for (int i = 0; i < pOrder + 1; i++)
	{
		for (int j = 0; j < pOrder + 1; j++)
		{
			diff(i, j) = GLPolynomialDerivative(j, i, xi);
		}
	}
	return diff;
}
#pragma endregion

#pragma region Gauss-Lobatto-Legendre Lagrangre quadrature points and polynomials

/* p = ID of polynomial */
double Jacobi::GLLPolynomial(int p, double x, ArrayXd xi)
{
	int pOrder = xi.size() - 1;
	double Pn;

	if (abs(x - xi(p)) < eps)
		Pn = 1;
	else
		Pn = (x - 1.0) * (x + 1.0) * PolynomialDerivative(0.0, 0.0, pOrder, x)
		/ (pOrder * (pOrder + 1.0) * Polynomial(0.0, 0.0, pOrder, xi(p)) * (x - xi(p)));
	return Pn;
}
/* n = ID of point */
double Jacobi::GLLPolynomialDerivative(int p, int n, ArrayXd xi)
{
	int pOrder = xi.size() - 1;
	double dPn;

	if (p == 0 && n == 0)
		dPn = -pOrder * (pOrder + 1.0) / 4.0;
	else if (p == pOrder && n == pOrder)
		dPn = pOrder * (pOrder + 1.0) / 4.0;
	else if (p == n)
		dPn = 0.0;
	else
		dPn = Polynomial(0.0, 0.0, pOrder, xi(n))
		/ (Polynomial(0.0, 0.0, pOrder, xi(p)) * (xi(n) - xi(p)));
	return dPn;
}

ArrayXd Jacobi::GLLPoints(int pOrder)
{
	ArrayXd x(pOrder + 1);           //Gauss-Lobatto-Legendre points

	x(0) = -1;
	x.segment(1, pOrder - 1) = PolynomialZeros(1.0, 1.0, pOrder - 1);
	x(pOrder) = 1;
	return x;
}

ArrayXd Jacobi::GLLWeights(ArrayXd x)
{
	ArrayXd w(x.size());

	int pOrder = x.size() - 1;
	for (int i = 0; i <= pOrder; i++)
	{
		w(i) = 2.0 / ((pOrder*(pOrder+1.0))
			* pow(Polynomial(0.0, 0.0, pOrder, x(i)), 2));
	}
	return w;
}

/* The same as the diffMatrix of GL */
MatrixXd Jacobi::GLLDiffMatrix(ArrayXd xi)
{
	int pOrder = xi.size() - 1;
	MatrixXd diff(pOrder + 1, pOrder + 1);

	for (int i = 0; i < pOrder + 1; i++)
	{
		for (int j = 0; j < pOrder + 1; j++)
		{
			diff(i, j) = GLLPolynomialDerivative(j, i, xi);
		}
	}
	return diff;
}
#pragma endregion

