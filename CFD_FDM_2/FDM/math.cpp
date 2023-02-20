#include"cfd_head.h"


void fft(int nx, int ny, vector< vector<double> > f, fftw_complex *ff)
{
    
    fftw_plan plan = fftw_plan_dft_r2c_2d(nx, ny, (double*)f.data(), (fftw_complex*)ff, FFTW_ESTIMATE); 

    fftw_execute(plan);

    fftw_destroy_plan(plan);

}

void ifft(int nx, int ny, vector< vector<double> > &u, vector< vector<fftw_complex> > uf)
{

    fftw_plan plan = fftw_plan_dft_c2r_2d(nx, ny, (fftw_complex*)uf.data(), (double*)u.data(), FFTW_ESTIMATE); 

    fftw_execute(plan);

    fftw_destroy_plan(plan);

}

/*
* Calculate L - 2 norm of a vector
*/
double ComputeL2norm(int nx, int ny, MatrixXd erro)
{
    double rms = 0.0;
    for (int i = 1; i < nx; i++)
    {
        for (int j = 1; j < ny; j++)
        {
            rms = rms + pow(erro(i, j), 2);
        } 
    }
    rms = sqrt(rms / ((nx - 1) * (ny - 1)));
    return rms;
}

void CalculateResidual(int nx, int ny, double dx, double dy, MatrixXd f, MatrixXd un, MatrixXd& r)
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

