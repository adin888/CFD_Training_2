#include"../Header/cfd_head.h"

/*
* -Using Conjugate Gradient methods to slove Poisson equation with Dirichlet boundary condition
* -Using the ratio of L2 to its initial value during the iteration as the criterion for convergence
* -Exact result:
*           ue(i, j) = (x(i) * x(i) - 1.0) * (y(i) * y(i) - 1.0);
*           f(i, j) = -2.0 * (2.0 - x(i) * x(i) - y(i) * y(i));
*/

MatrixXd ConjugateGradient(int nx, int ny, int maxIter, double dx, double dy, double rms, double initRms, double tolerance,
	MatrixXd r, MatrixXd f, MatrixXd un)
{
	double alpha;
	double beta;

	MatrixXd p = MatrixXd::Zero(nx + 1, ny + 1);
	MatrixXd Ap = MatrixXd::Zero(nx + 1, ny + 1);
	MatrixXd rOld = MatrixXd::Zero(nx + 1, ny + 1);
    MatrixXd rr = MatrixXd::Zero(nx + 1, ny + 1);

	CalculateResidual(nx, ny, dx, dy, f, un, r);

	rms = ComputeL2norm(nx, ny, r);
	initRms = rms;

	p = r;

	for (int iterNum = 1; iterNum < maxIter; iterNum++)
	{
		for (int j = 1; j < ny; j++)
		{
			for (int i = 1; i < nx; i++)
			{
				Ap(i, j) = (p(i + 1, j) - 2.0 * p(i, j) + p(i - 1, j)) / dx / dx 
						+ (p(i, j + 1) - 2.0 * p(i, j) + p(i, j - 1)) / dy / dy;
			}
		}

		alpha = ((r.block(1, 1, nx - 1, ny - 1).array() * r.block(1, 1, nx - 1, ny - 1).array()).sum()) /
			((Ap.block(1, 1, nx - 1, ny - 1).array() * p.block(1, 1, nx - 1, ny - 1).array()).sum() + 1e-16);

		un.block(1, 1, nx - 1, ny - 1) += alpha * p.block(1, 1, nx - 1, ny - 1);

		rOld = r;

		r.block(1, 1, nx - 1, ny - 1) -= alpha * Ap.block(1, 1, nx - 1, ny - 1);
        //
        CalculateResidual(nx, ny, dx, dy, f, un, rr);
		rms = ComputeL2norm(nx, ny, rr);

		cout << "The count of iteration is " << iterNum << endl;
		cout << "rms = " << setprecision(6) << rms / initRms << endl;

		if (rms / initRms <= tolerance)
			break;

		beta = ((r.block(1, 1, nx - 1, ny - 1).array() * r.block(1, 1, nx - 1, ny - 1).array()).sum()) /
			((rOld.block(1, 1, nx - 1, ny - 1).array() * rOld.block(1, 1, nx - 1, ny - 1).array()).sum() + 1e-16);

		p = r + beta * p;

	}
	return un;
}

void ConjugateGradientSolver()
{
    double x_l = 0.0;
    double x_r = 1.0;
    double y_l = 0.0;
    double y_r = 1.0;
    int nx = 128;
    int ny = 128;
    double dx = (x_r - x_l) / nx;
    double dy = (y_r - y_l) / ny;

    double tolerance = 1.0e-5;
    int maxIter = 1e5;
    double rms = 0.0;
    double initRms = 0.0;

    VectorXd x(nx + 1);
    VectorXd y(ny + 1);

    MatrixXd ue(nx + 1, ny + 1);
    MatrixXd delu(nx + 1, ny + 1);
    MatrixXd f(nx + 1, ny + 1);
    MatrixXd un = MatrixXd::Zero(nx + 1, ny + 1);
    MatrixXd r = MatrixXd::Zero(nx + 1, ny + 1);

    for (int i = 0; i < nx + 1; i++)
    {
        x[i] = x_l + dx * i;  //Assign node locations
    }
    for (int i = 0; i < ny + 1; i++)
    {
        y[i] = y_l + dy * i;  //Assign node locations
    }
    for (int i = 0; i < nx + 1; i++)
    {
        for (int j = 0; j < ny + 1; j++)
        {
            ue(i, j) = (x(i) * x(i) - 1.0) * (y(j) * y(j) - 1.0);
            f(i, j) = -2.0 * (2.0 - x(i) * x(i) - y(j) * y(j));
        }
    }

    /* Dirichlet boundary condition */
    un.row(0) = ue.row(0);
    un.row(nx) = ue.row(nx);
    un.col(0) = ue.col(0);
    un.col(ny) = ue.col(ny);

    un = ConjugateGradient(nx, ny, maxIter, dx, dy, rms, initRms, tolerance, r, f, un);

    delu = ue - un;
    ofstream u;
    ofstream du;

    u.open("u_final.csv", ios::out | ios::trunc);
    du.open("u_delta.csv", ios::out | ios::trunc);

    for (int i = 0; i < nx + 1; i++)
    {
        u << setprecision(16) << x[i] << ",";
        du << setprecision(16) << x[i] << ",";
        for (int j = 0; j < ny + 1; j++)
        {
            u << setprecision(16) << un(i, j) << ",";
            du << setprecision(16) << delu(i, j) << ",";
        }
        u << endl;
        du << endl;
    }

    u.close();
    du.close();
}