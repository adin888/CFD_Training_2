#pragma once
#include <iostream>
#include<fstream>
#include<iomanip>
#include<math.h>
#include<algorithm>
#include<vector>

using namespace std;

#define PI acos(-1)

void Roe_Riemann_Solver();
vector< vector<double> > fluxes(int nx, double gamma, vector< vector<double> > q, vector< vector<double> > f);
vector< vector<double> > roe(int nx, double gamma, vector< vector<double> > qL, vector< vector<double> > qR,
    vector< vector<double> > f, vector< vector<double> > fL, vector< vector<double> > fR);

double wL(double u1, double u2, double u3, double u4, double u5);
double wR(double u1, double u2, double u3, double u4, double u5);
vector< vector<double> > RRweno5L(int nx, vector< vector<double> > u, vector< vector<double> > fL);
vector< vector<double> > RRweno5R(int nx, vector< vector<double> > u, vector< vector<double> > fR);
