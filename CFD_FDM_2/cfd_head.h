#pragma once
#include <iostream>
#include<fstream>
#include<iomanip>
#include<math.h>
#include<algorithm>
#include<vector>
#include <cstring>
#include <fftw3.h>
#include <Eigen\Dense>

using namespace std;
using namespace Eigen;

#define PI acos(-1)

#pragma region Roe_Riemann

void Roe_Riemann_Solver();
vector< vector<double> > roe(int nx, double gamma, vector< vector<double> > qL, vector< vector<double> > qR,
    vector< vector<double> > f, vector< vector<double> > fL, vector< vector<double> > fR);
vector< vector<double> > RRweno5L(int nx, vector< vector<double> > u, vector< vector<double> > fL);
vector< vector<double> > RRweno5R(int nx, vector< vector<double> > u, vector< vector<double> > fR);

#pragma endregion

#pragma region HLLC_Riemann

void HLLC_Riemann_Solver();
vector< vector<double> > HLLC(int nx, double gamma, vector< vector<double> > qL, vector< vector<double> > qR,
    vector< vector<double> > f, vector< vector<double> > fL, vector< vector<double> > fR);

#pragma endregion

#pragma region Rusanov_Riemann

void Rusanov_Riemann_Solver();
vector< vector<double> > Rusanov(int nx, double gamma, vector< vector<double> > qL, vector< vector<double> > qR,
    vector< vector<double> > f, vector< vector<double> > fL, vector< vector<double> > fR);
vector<double> Wavespeed(int nx, double gamma, vector< vector<double> > qL, vector< vector<double> > qR, vector<double> ps);

#pragma endregion

#pragma region FFTSolver
void FFT_Solver();

#pragma endregion

#pragma region GaussSeidelSolver
void GaussSeidelSolver();

#pragma endregion

#pragma region Function

vector< vector<double> > fluxes(int nx, double gamma, vector< vector<double> > q, vector< vector<double> > f);
double wL(double u1, double u2, double u3, double u4, double u5);
double wR(double u1, double u2, double u3, double u4, double u5);
void fft(int nx, int ny, vector< vector<double> > f, fftw_complex *ff);
void ifft(int nx, int ny, vector< vector<double> > &u, vector< vector<fftw_complex> > uf);
double ComputeL2norm(int nx, int ny, MatrixXd erro);
#pragma endregion

