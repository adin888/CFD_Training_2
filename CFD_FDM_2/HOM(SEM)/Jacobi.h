#pragma once
#include"cfd_head.h"
class Jacobi
{ /* ���԰������ظ��ı���ͨ�����캯��һ���Ի�� ���Ͳ���ÿ�������������� */
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

