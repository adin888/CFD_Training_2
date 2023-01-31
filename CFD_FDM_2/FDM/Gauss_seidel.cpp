#include"../Header/cfd_head.h"

/*
* -Using Gauss-Seidel iterative methods to slove Poisson equation with Dirichlet boundary condition
* -Using the ratio of L2 to its initial value during the iteration as the criterion for convergence
* -Exact result: 
*           ue(i, j) = (x(i) * x(i) - 1.0) * (y(i) * y(i) - 1.0);
*           f(i, j) = -2.0 * (2.0 - x(i) * x(i) - y(i) * y(i));
*/


MatrixXd GaussSeidel(int nx, int ny, int maxIter, double dx, double dy, double rms, double initRms, double tolerance,
	MatrixXd r, MatrixXd f, MatrixXd un)
{
	CalculateResidual(nx, ny, dx, dy, f, un, r);

	rms = ComputeL2norm(nx, ny, r);
	initRms = rms;

	double dxy = -2.0 / dx / dx - 2.0 / dy / dy;

	for (int iterNum = 1; iterNum < maxIter; iterNum++)
	{
		CalculateResidual(nx, ny, dx, dy, f, un, r);

		un = un + r / dxy;

		CalculateResidual(nx, ny, dx, dy, f, un, r);

		rms = ComputeL2norm(nx, ny, r);

		cout << "The count of iteration is " << iterNum << endl;
		cout << "rms = " << setprecision(6) << rms / initRms << endl;

		if (rms / initRms <= tolerance)
			break;
	}
	return un;
}

void GaussSeidelSolver()
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

    MatrixXd un = MatrixXd::Zero(nx + 1, ny + 1);
    MatrixXd r = MatrixXd::Zero(nx + 1, ny + 1);
    MatrixXd ue(nx + 1, ny + 1);
    MatrixXd delu(nx + 1, ny + 1);
    MatrixXd f(nx + 1, ny + 1);

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

    un = GaussSeidel(nx, ny, maxIter, dx, dy, rms, initRms, tolerance, r, f, un);

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