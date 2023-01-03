#include"cfd_head.h"

MatrixXd Numerical_FFT(int nx, int ny, double dx, double dy, MatrixXd f, MatrixXd u)
{
    VectorXd wx(nx);
    VectorXd wy(ny);
    
    int rounded = floor((nx + 1) / 2 + 1);
    MatrixXcd fF(rounded, ny + 1);
    MatrixXcd uF(rounded, ny + 1);

    MatrixXd udebug(u.rows(), u.cols());

    
    double xx = 2.0 / (dx * dx);
    double yy = 2.0 / (dy * dy);
    double xy = -xx - yy;

    double dw = 2.0 * PI / nx;                     // nx = ny in this example
    double eps = 1.0e-6;
    
    int half = ceil(nx / 2);
    for (int i = 0; i < half; i++)
    {
        wx(i) = dw * i;
        wx(i + half) = dw * (i - half);
    }

    wx(0) = eps;

    wy = wx;

    //fft(nx, ny, f, ff);
    fftw_plan plan;
    plan = fftw_plan_dft_r2c_2d(nx + 1, ny + 1, f.data(), (fftw_complex*)fF.data(), FFTW_ESTIMATE);

    fftw_execute(plan);
    
    fftw_destroy_plan(plan);

    for (int i = 0; i < rounded; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            //uF(i, j) = fF(i, j) / (xy + xx * cos(wx(i)) + yy * cos(wy(j)));
            uF(i, j) = fF(i, j) / (-pow(wx(i),2)- pow(wy(i), 2));
        }
    }

    plan = fftw_plan_dft_c2r_2d(nx + 1, ny + 1, (fftw_complex*)fF.data(), u.data(), FFTW_ESTIMATE);

    fftw_execute(plan);

    fftw_destroy_plan(plan);

    u.row(nx) = u.row(0);
    u.col(ny) = u.col(0);
    
    //ifft(nx, ny, u, uf);
   
    return u / ((nx + 1) * (ny + 1));
}

void FFT_Solver()
{
    double x_l = 0.0;
    double x_r = 1.0;
    double y_l = 0.0;
    double y_r = 1.0;
    int nx = 512;
    int ny = 512;
    double dx = (x_r - x_l) / nx;
    double dy = (y_r - y_l) / ny;

    VectorXd x(nx + 1);
    VectorXd y(ny + 1);

    MatrixXd un(nx + 1, ny + 1);
    MatrixXd ue(nx + 1, ny + 1);
    MatrixXd delu(nx + 1, ny + 1);
    MatrixXd f(nx + 1, ny + 1);

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
    for (int i = 0; i < nx + 1 ; i++)
    {
        for (int j = 0; j < ny + 1; j++)
        {
            ue(i,j)= sin(2.0 * PI * x(i)) * sin(2.0 * PI * y(j))
                + c1 * sin(km * 2.0 * PI * x(i)) * sin(km * 2.0 * PI * y(j));
            f(i, j) = c2 * sin(2.0 * PI * x(i)) * sin(2.0 * PI * y(j)) 
                + c2 * sin(km * 2.0 * PI * x(i)) * sin(km * 2.0 * PI * y(j));
        }
    }

    un = Numerical_FFT(nx, ny, dx, dy, f, un);

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