#include"../Header/cfd_head.h"

/* 
* -viscous incompressible flow 
* -A square cavity consisting of three rigid walls with no-slip conditions and a lid moving with a tangential unit velocity
* -Using Runge - Kutta - 3 Scheme for time integration
* -Using second-order Arakawa scheme for the nonlinear terms dphi/dy*dw/dx-dphi/dx*dw/dy
* -Total time t=0.1s
*/


void JBC(int nx, int ny, double dx, double dy, MatrixXd& w, MatrixXd psi)
{
	VectorXd one = VectorXd::Ones(nx + 1);

	/* second order approximation for Jensen-vorticity boundary condition */
	/* left and right */
	w.col(0) = (-4.0 * psi.col(1) + 0.5 * psi.col(2)) / (dx * dx);
	w.col(nx) = (-4.0 * psi.col(nx - 1) + 0.5 * psi.col(nx - 2)) / (dx * dx);

	/* bottom and top */
	w.row(0) = (-4.0 * psi.row(1) + 0.5 * psi.row(2)) / (dy * dy);
	w.row(ny) = (-4.0 * psi.row(ny - 1) + 0.5 * psi.row(ny - 2)) / (dy * dy) - one.transpose() * 3.0 / dy;
}

void LDCRHS(int nx, int ny, double dx, double dy, double re, MatrixXd w, MatrixXd psi, MatrixXd &r)
{
	double xx = 1.0 / (re * dx * dx);
	double yy = 1.0 / (re * dy * dy);
	double xy = 1.0 / (4.0 * dx * dy);

	double jac, j1, j2, j3;

	for (int j = 1; j < ny; j++)
	{
		for (int i = 1; i < nx; i++)
		{
			j1 = xy * ((w(i + 1, j) - w(i - 1, j)) * (psi(i, j + 1) - psi(i, j - 1))
				- (w(i, j + 1) - w(i, j - 1)) * (psi(i + 1, j) - psi(i - 1, j)));

			j2 = xy * (w(i + 1, j) * (psi(i + 1, j + 1) - psi(i + 1, j - 1))
				- w(i - 1, j) * (psi(i - 1, j + 1) - psi(i - 1, j - 1))
				- w(i, j + 1) * (psi(i + 1, j + 1) - psi(i - 1, j + 1))
				+ w(i, j - 1) * (psi(i + 1, j - 1) - psi(i - 1, j - 1)));

			j3 = xy * (w(i + 1, j + 1) * (psi(i, j + 1) - psi(i + 1, j))
				- w(i - 1, j - 1) * (psi(i - 1, j) - psi(i, j - 1))
				- w(i - 1, j + 1) * (psi(i, j + 1) - psi(i - 1, j))
				+ w(i + 1, j - 1) * (psi(i + 1, j) - psi(i, j - 1)));

			jac = (j1 + j2 + j3) / 3.0;

			r(i, j) = -jac + xx * (w(i + 1, j) - 2.0 * w(i, j) + w(i - 1, j))
				+ yy * (w(i, j + 1) - 2.0 * w(i, j) + w(i, j - 1));
		}
	}
}

void LDCNumerical(int nx, int ny, int nt, double dx, double dy, double dt, double re, MatrixXd &w, MatrixXd &psi)
{
	MatrixXd r = MatrixXd::Zero(nx + 1, ny + 1);
	MatrixXd wTem = MatrixXd::Zero(nx + 1, ny + 1);
	MatrixXd psiOld = MatrixXd::Zero(nx + 1, ny + 1);

	ArrayXXd rmsArray = ArrayXXd::Zero(nx + 1, ny + 1);

	double rms;

	for (int t = 1; t < nt; t++)
	{
		psiOld = psi;
		/* 1st */
		LDCRHS(nx, ny, dx, dy, re, w, psi, r);

		wTem.block(1, 1, nx - 1, ny - 1) = w.block(1, 1, nx - 1, ny - 1) + dt * r.block(1, 1, nx - 1, ny - 1);

		JBC(nx, ny, dx, dy, wTem, psi);

		psi = Numerical_FST(nx, ny, dx, dy, -wTem, psi);
		/* 2nd */
		LDCRHS(nx, ny, dx, dy, re, wTem, psi, r);

		wTem.block(1, 1, nx - 1, ny - 1) = 0.75 * w.block(1, 1, nx - 1, ny - 1)
										+ 0.25 * wTem.block(1, 1, nx - 1, ny - 1)
										+ 0.25 * dt * r.block(1, 1, nx - 1, ny - 1);

		JBC(nx, ny, dx, dy, wTem, psi);

		psi = Numerical_FST(nx, ny, dx, dy, -wTem, psi);

		/* 3rd */

		LDCRHS(nx, ny, dx, dy, re, wTem, psi, r);

		w.block(1, 1, nx - 1, ny - 1) = (1.0/3.0) * w.block(1, 1, nx - 1, ny - 1)
			+ (2.0/3.0) * wTem.block(1, 1, nx - 1, ny - 1)
			+ (2.0 / 3.0) * dt * r.block(1, 1, nx - 1, ny - 1);

		JBC(nx, ny, dx, dy, w, psi);

		psi = Numerical_FST(nx, ny, dx, dy, -w, psi);

		rmsArray = psi - psiOld;
		rms = (rmsArray * rmsArray).sum() / ((nx + 1) * (ny + 1));

		cout << "The count of iteration is " << t << endl;
		cout << "rms = " << setprecision(6) << rms << endl;
	}	
}

void LDC_Solver()
{
	double x_l = 0.0;
	double x_r = 1.0;
	double y_l = 0.0;
	double y_r = 1.0;
	int nx = 64;
	int ny = 64;
	double dx = (x_r - x_l) / nx;
	double dy = (y_r - y_l) / ny;

	double dt = 0.001;
	double tEnd = 0.1;
	int nt = ceil(tEnd / dt);

	double re = 100.0;

	VectorXd x(nx + 1);
	VectorXd y(ny + 1);

	MatrixXd w = MatrixXd::Zero(nx + 1, ny + 1);                       //Vorticity
	MatrixXd psi = MatrixXd::Zero(nx + 1, ny + 1);                     //Stream

	for (int i = 0; i < nx + 1; i++)
	{
		x[i] = x_l + dx * i;  //Assign node locations
	}
	for (int i = 0; i < ny; i++)
	{
		y[i] = y_l + dy * i;  //Assign node locations
	}

	LDCNumerical(nx, ny, nt, dx, dy, dt, re, w, psi);

	ofstream wStream;
	ofstream psiStream;

	wStream.open("w_final.csv", ios::out | ios::trunc);
	psiStream.open("psi_final.csv", ios::out | ios::trunc);

	for (int i = 0; i < nx + 1; i++)
	{
		wStream << setprecision(16) << x[i] << ",";
		psiStream << setprecision(16) << x[i] << ",";
		for (int j = 0; j < ny + 1; j++)
		{
			wStream << setprecision(16) << w(j, i) << ",";
			psiStream << setprecision(16) << psi(j, i) << ",";
		}
		wStream << endl;
		psiStream << endl;
	}

	wStream.close();
	psiStream.close();
}