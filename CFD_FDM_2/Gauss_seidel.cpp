#include"cfd_head.h"

void CalculateResidual(int nx, int ny, double dx, double dy, MatrixXd f, MatrixXd un, MatrixXd &r)
{
	double dxx, dyy;
	for (int j = 1; j < nx; j++)
	{
		for (int i = 1; i < ny; i++)
		{
			dxx = (un(i + 1, j) - 2 * un(i, j) + un(i - 1, j)) / dx / dx;
			dyy = (un(i, j + 1) - 2 * un(i, j) + un(i, j - 1)) / dy / dy;
			r(i, j) = f(i, j) - dxx - dyy;
		}
	}
}

MatrixXd GaussSeidel(int nx, int ny, int maxIter, double dx, double dy, double rms, double initRms, double tolerance,
	MatrixXd r, MatrixXd f, MatrixXd un)
{
	CalculateResidual(nx, ny, dx, dy, f, un, r);

	rms = ComputeL2norm(nx, ny, r);
	initRms = rms;

	double dxy = -2.0 / dx / dx - 2.0 / dy / dy;

	for (int iterNum = 1; iterNum < 5*maxIter; iterNum++)
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
    int nx = 512;
    int ny = 512;
    double dx = (x_r - x_l) / nx;
    double dy = (y_r - y_l) / ny;

    double tolerance = 1.0e-10;
    int maxIter = 2e6;
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
    for (int i = 0; i < ny; i++)
    {
        y[i] = y_l + dy * i;  //Assign node locations
    }
    for (int i = 0; i < nx + 1; i++)
    {
        for (int j = 0; j < ny + 1; j++)
        {
            ue(i, j) = (x(i) * x(i) - 1.0) * (y(i) * y(i) - 1.0);
            f(i, j) = -2.0 * (2.0 - x(i) * x(i) - y(i) * y(i));
        }
    }

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