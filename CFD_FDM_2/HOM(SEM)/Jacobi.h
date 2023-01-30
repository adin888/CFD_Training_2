#pragma once
#include"cfd_head.h"
class Jacobi
{ /* 可以把下列重复的变量通过构造函数一次性获得 ，就不用每个函数都输入了 */
public:
	double Polynomial(double, double, int, double);
	double PolynomialDerivative(double, double, int, double);
	ArrayXd PolynomialZeros(double, double, int);

	double GLPolynomial(int, double, ArrayXd);
	double GLPolynomialDerivative(int, int, ArrayXd);
	ArrayXd GLPoints(int);
	ArrayXd GLWeights(ArrayXd);
	MatrixXd GLDiffMatrix(ArrayXd);

	double GLLPolynomial(int, double, ArrayXd);
	double GLLPolynomialDerivative(int, int, ArrayXd);
	ArrayXd GLLPoints(int);
	ArrayXd GLLWeights(ArrayXd);
	MatrixXd GLLDiffMatrix(ArrayXd);
private:
	double eps = 1.0e-20;
};

