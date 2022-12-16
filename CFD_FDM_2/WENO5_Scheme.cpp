#include"cfd_head.h"

/*
* WENO reconstruction for left states at the interface in Roe Riemann Solver
*/
vector< vector<double> > RRweno5L(int nx, vector< vector<double> > u, vector< vector<double> > fL)
{
    double u1, u2, u3, u4, u5;              //Reconstruction with 5 points

    for (int m = 0; m < 3; m++)
    {
        int i = -1;
        u1 = u[i + 3][m];
        u2 = u[i + 2][m];
        u3 = u[i + 1][m];
        u4 = u[i + 1][m];
        u5 = u[i + 2][m];
        fL[i + 1][m] = wL(u1, u2, u3, u4, u5);

        i = 0;
        u1 = u[i + 1][m];
        u2 = u[i][m];
        u3 = u[i][m];
        u4 = u[i + 1][m];
        u5 = u[i + 2][m];
        fL[i + 1][m] = wL(u1, u2, u3, u4, u5);

        i = 1;
        u1 = u[i - 1][m];
        u2 = u[i - 1][m];
        u3 = u[i][m];
        u4 = u[i + 1][m];
        u5 = u[i + 2][m];
        fL[i + 1][m] = wL(u1, u2, u3, u4, u5);

        for (int i = 2; i < nx - 2; i++)
        {
            u1 = u[i - 2][m];
            u2 = u[i - 1][m];
            u3 = u[i][m];
            u4 = u[i + 1][m];
            u5 = u[i + 2][m];
            fL[i + 1][m] = wL(u1, u2, u3, u4, u5);
        }

        i = nx - 2;
        u1 = u[i - 2][m];
        u2 = u[i - 1][m];
        u3 = u[i][m];
        u4 = u[i + 1][m];
        u5 = u[i + 1][m];
        fL[i + 1][m] = wL(u1, u2, u3, u4, u5);

        i = nx - 1;
        u1 = u[i - 2][m];
        u2 = u[i - 1][m];
        u3 = u[i][m];
        u4 = u[i][m];
        u5 = u[i - 1][m];
        fL[i + 1][m] = wL(u1, u2, u3, u4, u5);
    }

    return fL;
}

/*
* WENO reconstruction for right states at the interface in Roe Riemann Solver
*/
vector< vector<double> > RRweno5R(int nx, vector< vector<double> > u, vector< vector<double> > fR)
{
    double u1, u2, u3, u4, u5;              //Reconstruction with 5 points

    for (int m = 0; m < 3; m++)
    {
        int i = 0;
        u1 = u[i + 1][m];
        u2 = u[i][m];
        u3 = u[i][m];
        u4 = u[i + 1][m];
        u5 = u[i + 2][m];
        fR[i][m] = wR(u1, u2, u3, u4, u5);

        i = 1;
        u1 = u[i - 1][m];
        u2 = u[i - 1][m];
        u3 = u[i][m];
        u4 = u[i + 1][m];
        u5 = u[i + 2][m];
        fR[i][m] = wR(u1, u2, u3, u4, u5);

        for (int i = 2; i < nx - 2; i++)
        {
            u1 = u[i - 2][m];
            u2 = u[i - 1][m];
            u3 = u[i][m];
            u4 = u[i + 1][m];
            u5 = u[i + 2][m];
            fR[i][m] = wR(u1, u2, u3, u4, u5);
        }

        i = nx - 2;
        u1 = u[i - 2][m];
        u2 = u[i - 1][m];
        u3 = u[i][m];
        u4 = u[i + 1][m];
        u5 = u[i + 1][m];
        fR[i][m] = wR(u1, u2, u3, u4, u5);

        i = nx - 1;
        u1 = u[i - 2][m];
        u2 = u[i - 1][m];
        u3 = u[i][m];
        u4 = u[i][m];
        u5 = u[i - 1][m];
        fR[i][m] = wR(u1, u2, u3, u4, u5);

        i = nx;
        u1 = u[i - 2][m];
        u2 = u[i - 1][m];
        u3 = u[i - 1][m];
        u4 = u[i - 2][m];
        u5 = u[i - 3][m];
        fR[i][m] = wR(u1, u2, u3, u4, u5);
    }

    return fR;
}

/*
* nonlinear weights for upwind direction
*/
double wL(double u1, double u2, double u3, double u4, double u5)
{
    double eps = 1.0e-6;
    double u_L;

    /* smoothness indicators */
    double s1 = pow((13.0 / 12.0) * (u1 - 2.0 * u2 + u3), 2) + pow(0.25 * (u1 - 4.0 * u2 + 3.0 * u3), 2);
    double s2 = pow((13.0 / 12.0) * (u2 - 2.0 * u3 + u4), 2) + pow(0.25 * (u2 - u4), 2);
    double s3 = pow((13.0 / 12.0) * (u3 - 2.0 * u4 + u5), 2) + pow(0.25 * (3.0 * u3 - 4.0 * u4 + u5), 2);

    double d1 = 1.0 / 10.0;
    double d2 = 3.0 / 5.0;
    double d3 = 3.0 / 10.0;

    double c1 = d1 / pow(eps + s1, 2);
    double c2 = d2 / pow(eps + s2, 2);
    double c3 = d3 / pow(eps + s3, 2);

    double w1 = c1 / (c1 + c2 + c3);
    double w2 = c2 / (c1 + c2 + c3);
    double w3 = c3 / (c1 + c2 + c3);

    double b1 = u1 / 3.0 - 7.0 / 6.0 * u2 + 11.0 / 6.0 * u3;
    double b2 = -u2 / 6.0 + 5.0 / 6.0 * u3 + u4 / 3.0;
    double b3 = u3 / 3.0 + 5.0 / 6.0 * u4 - u5 / 6.0;

    return u_L = w1 * b1 + w2 * b2 + w3 * b3;
}

/*
* nonlinear weights for downwind direction
*/
double wR(double u1, double u2, double u3, double u4, double u5)
{
    double eps = 1.0e-6;
    double u_R;

    /* smoothness indicators */
    double s1 = pow((13.0 / 12.0) * (u1 - 2.0 * u2 + u3), 2) + pow(0.25 * (u1 - 4.0 * u2 + 3.0 * u3), 2);
    double s2 = pow((13.0 / 12.0) * (u2 - 2.0 * u3 + u4), 2) + pow(0.25 * (u2 - u4), 2);
    double s3 = pow((13.0 / 12.0) * (u3 - 2.0 * u4 + u5), 2) + pow(0.25 * (3.0 * u3 - 4.0 * u4 + u5), 2);

    double d1 = 3.0 / 10.0;
    double d2 = 3.0 / 5.0;
    double d3 = 1.0 / 10.0;

    double c1 = d1 / pow(eps + s1, 2);
    double c2 = d2 / pow(eps + s2, 2);
    double c3 = d3 / pow(eps + s3, 2);

    double w1 = c1 / (c1 + c2 + c3);
    double w2 = c2 / (c1 + c2 + c3);
    double w3 = c3 / (c1 + c2 + c3);

    double b1 = -u1 / 6.0 + 5.0 / 6.0 * u2 + u3 / 3.0;
    double b2 = u2 / 3.0 + 5.0 / 6.0 * u3 - u4 / 6.0;
    double b3 = 11.0 / 6.0 * u3 - 7.0 / 6.0 * u4 + u5 / 3.0;

    return u_R = w1 * b1 + w2 * b2 + w3 * b3;
}