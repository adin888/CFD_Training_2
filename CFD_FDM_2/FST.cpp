#include"cfd_head.h"

MatrixXd Numerical_FST(int nx, int ny, double dx, double dy, MatrixXd f, MatrixXd u)
{
    MatrixXd fTemp(nx - 1, ny - 1);
    MatrixXd fF(nx - 1, ny - 1);
    MatrixXd uF(nx - 1, ny - 1);

    double xx = 2.0 / (dx * dx);
    double yy = 2.0 / (dy * dy);
    double xy = -xx - yy;

    double eps = 1.0e-6;
    double sumTemp = 0;

    fTemp = f.block(1, 1, nx - 1, ny - 1);

    for (int n = 0; n < ny - 1; n++)
    {
        for (int m = 0; m < nx - 1; m++)
        {
            for (int j = 0; j < ny - 1; j++)
            {
                for (int i = 0; i < nx - 1; i++)
                {
                    sumTemp += fTemp(i, j) * sin(PI * m * i / nx) * sin(PI * n * j / ny);
                }
            }
            fF(m, n) = sumTemp;
            sumTemp = 0;
        }
    }

    for (int i = 0; i < nx - 1; i++)
    {
        for (int j = 0; j < ny - 1; j++)
        {
            uF(i, j) = fF(i, j) / (xy + xx * cos(PI * i / nx) + yy * cos(PI * j / ny) + 1e-16);
        }
    }

    for (int j = 0; j < ny - 1; j++)
    {
        for (int i = 0; i < nx - 1; i++)
        {
            for (int n = 0; n < ny - 1; n++)
            {
                for (int m = 0; m < nx - 1; m++)
                {
                    sumTemp += uF(m, n) * sin(PI * m * i / nx) * sin(PI * n * j / ny);
                }
            }
            u(i + 1, j + 1) = sumTemp * 4 / nx / ny;
            sumTemp = 0;
        }
    }

    return u;
}

void FST_Solver()
{
    double x_l = 0.0;
    double x_r = 1.0;
    double y_l = 0.0;
    double y_r = 1.0;
    int nx = 128;
    int ny = 128;
    double dx = (x_r - x_l) / nx;
    double dy = (y_r - y_l) / ny;

    VectorXd x(nx + 1);
    VectorXd y(ny + 1);

    MatrixXd un = MatrixXd::Zero(nx + 1, ny + 1);
    MatrixXd ue = MatrixXd::Zero(nx + 1, ny + 1);
    MatrixXd f = MatrixXd::Zero(nx + 1, ny + 1);
    ArrayXXd delu = ArrayXXd::Zero(nx + 1, ny + 1);

    double km = 16.0;
    double c1 = pow((1.0 / km), 2);
    double c2 = -8.0 * PI * PI;

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
            ue(i, j) = sin(2.0 * PI * x(i)) * sin(2.0 * PI * y(j))
                + c1 * sin(km * 2.0 * PI * x(i)) * sin(km * 2.0 * PI * y(j));
            f(i, j) = c2 * sin(2.0 * PI * x(i)) * sin(2.0 * PI * y(j))
                + c2 * sin(km * 2.0 * PI * x(i)) * sin(km * 2.0 * PI * y(j));
        }
    }

    un = Numerical_FST(nx, ny, dx, dy, f, un);

    delu = (ue - un).array()/ue.array();
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