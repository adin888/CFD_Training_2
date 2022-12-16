#include"cfd_head.h"

/*
* -Using Reo approximate Riemann solver 
* -Using WENO-5 Scheme to reconstruct the left and right side fluxes at the interface
* -Using Runge-Kutta-3 Scheme for time integration
* -Data is saved every 0.025s
*/

/*
* Calculate right hand term of the inviscid Burgers equation
* r = -udu/dx
*/
vector< vector<double> > RRrhs(int nx, double dx, double gamma, vector< vector<double> > q, vector< vector<double> > r)
{
    vector< vector<double> > qL(nx + 1, vector<double>(3));          //left and right states of field at discret nodal points
    vector< vector<double> > qR(nx + 1, vector<double>(3));
    vector< vector<double> > fL(nx + 1, vector<double>(3));          //left and right states of fluxes at the interface   
    vector< vector<double> > fR(nx + 1, vector<double>(3));          

    vector< vector<double> > f(nx + 1, vector<double>(3));           //flux at the interface computed using Rusanov scheme


    qL = RRweno5L(nx, q, qL);
    qR = RRweno5R(nx, q, qR);

    fL = fluxes(nx, gamma, qL, fL);
    fR = fluxes(nx, gamma, qR, fR);

    f = roe(nx, gamma, qL, qR, f, fL, fR);

    for (int i = 0; i < nx; i++)
    {
        for (int m = 0; m < 3; m++)
        {
            r[i][m] = -(f[i + 1][m] - f[i][m]) / dx;
        }
        
    }

    return r;
}

vector< vector<double> > roe(int nx, double gamma, vector< vector<double> > qL, vector< vector<double> > qR, 
    vector< vector<double> > f, vector< vector<double> > fL, vector< vector<double> > fR)
{
    vector<double> dq(3);                              //Temporary variable used to compute f
    vector<double> ddq(3);
    vector<double> dddq(3);
    double rL, uL, eL, pL, hL, rR, uR, eR, pR, hR;
    double alpha, beta, phi2;
    double uu, hh, aa;
    double J1, J2, J3;
    double L11, L12, L13, L21, L22, L23, L31, L32, L33;
    double R11, R12, R13, R21, R22, R23, R31, R32, R33;

    for (int i = 0; i < nx + 1; i++)
    {
        /* Compute Left and right states of rho(r), u, p, h, e */
        rL = qL[i][0];
        uL = qL[i][1] / rL;
        eL = qL[i][2] / rL;
        pL = (gamma - 1.0) * (eL * rL - 0.5 * rL * (uL * uL));
        hL = eL + pL / rL;

        rR = qR[i][0];
        uR = qR[i][1] / rR;
        eR = qR[i][2] / rR;
        pR = (gamma - 1.0) * (eR * rR - 0.5 * rR * (uR * uR));
        hR = eR + pR / rR;

        /* The Roe averaging formulas */
        alpha = 1.0 / (sqrt(abs(rL)) + sqrt(abs(rR)));

        uu = (sqrt(abs(rL)) * uL + sqrt(abs(rR)) * uR) * alpha;
        hh = (sqrt(abs(rL)) * hL + sqrt(abs(rR)) * hR) * alpha;
        aa = sqrt(abs((gamma - 1.0) * (hh - 0.5 * uu * uu)));

        /* Eigenvalues of the Jacobian matrix */
        J1 = abs(uu);
        J2 = abs(uu + aa);
        J3 = abs(uu - aa);

        /* Left and Right eigenvector matrix*/
        beta = 0.5 / (aa * aa);
        phi2 = 0.5 * (gamma - 1.0) * uu * uu;

        L11 = 1.0 - phi2 / (aa * aa);
        L12 = (gamma - 1.0) * uu / (aa * aa);
        L13 = -(gamma - 1.0) / (aa * aa);
        L21 = phi2 - uu * aa;
        L22 = aa - (gamma - 1.0) * uu;
        L23 = gamma - 1.0;
        L31 = phi2 + uu * aa;
        L32 = -aa - (gamma - 1.0) * uu;
        L33 = gamma - 1.0;

        R11 = 1.0;
        R12 = uu;
        R13 = phi2 / (gamma - 1.0);
        R21 = beta;
        R22 = beta * (uu + aa);
        R23 = beta * (hh + uu * aa);
        R31 = beta;
        R32 = beta * (uu - aa);
        R33 = beta * (hh - uu * aa);
         /* Compute f using Roe Solver */
        for (int m = 0; m < 3; m++)
        {
            dq[m] = 0.5 * (qR[i][m] - qL[i][m]);
        }

        ddq[0] = J1 * (L11 * dq[0] + L12 * dq[1] + L13 * dq[2]);
        ddq[1] = J2 * (L21 * dq[0] + L22 * dq[1] + L23 * dq[2]);
        ddq[2] = J3 * (L31 * dq[0] + L32 * dq[1] + L33 * dq[2]);

        dddq[0] = R11 * ddq[0] + R12 * ddq[1] + R13 * ddq[2];
        dddq[1] = R21 * ddq[0] + R22 * ddq[1] + R23 * ddq[2];
        dddq[2] = R31 * ddq[0] + R32 * ddq[1] + R33 * ddq[2];

        for (int m = 0; m < 3; m++)
        {
            f[i][m] = 0.5 * ((fR[i][m]) + fL[i][m]) - dddq[m];
        }
    }
    return f;
}
vector < vector< vector<double> > > numerical_Roe(int nx, int ns, int nt, double dx, double dt, vector<double> x, vector < vector< vector<double> > > q)
{
    vector< vector<double> > qn(nx, vector<double>(3));
    vector< vector<double> > qt(nx, vector<double>(3));
    vector< vector<double> > r(nx, vector<double>(3));
    double rho, u, p, e;

    int freq = ceil(nt / ns);
    
    double gamma = 1.4;            //specific gas ratio

    /* Sod's Riemann problem */
    double rhoL = 1.0;               
    double uL = 0.0;
    double pL = 1.0;

    double rhoR = 0.125;
    double uR = 0.0;
    double pR = 0.1;

    double xc = 0.5;               //seperator location
    
    for (int i = 0; i < nx; i++)
    {
        if (x[i] > xc)
        {
            rho = rhoR;
            u = uR;
            p = pR;
        }
        else
        {
            rho = rhoL;
            u = uL;
            p = pL;
        }
        e = p / (rho * (gamma - 1.0)) + 0.5 * u * u;

        qn[i][0] = rho;
        qn[i][1] = rho * u;
        qn[i][2] = rho * e;
    }

    for (int i = 0; i < nx; i++)
    {
        for (int m = 0; m < 3; m++)
        {
            q[i][m][0] = qn[i][m];
        }
    }

    for (int j = 1; j < nt + 1; j++)
    {
        r = RRrhs(nx, dx, gamma, qn, r);
        for (int i = 0; i < nx; i++)
        {
            for (int m = 0; m < 3; m++)
            {
                qt[i][m] = qn[i][m] + dt * r[i][m];
            }
        }

        r = RRrhs(nx, dx, gamma, qt, r);
        for (int i = 0; i < nx; i++)
        {
            for (int m = 0; m < 3; m++)
            {
            qt[i][m] = 0.75 * qn[i][m] + 0.25 * qt[i][m] + 0.25 * dt * r[i][m];
            }
        }

        r = RRrhs(nx, dx, gamma, qt, r);
        for (int i = 0; i < nx; i++)
        {
            for (int m = 0; m < 3; m++)
            {
                qn[i][m] = (1.0 / 3.0) * qn[i][m] + (2.0 / 3.0) * qt[i][m] + (2.0 / 3.0) * dt * r[i][m];
            } 
        }

        if (j % freq == 0)
        {
            for (int i = 0; i < nx; i++)
            {
                for (int m = 0; m < 3; m++)
                {
                    q[i][m][j / freq] = qn[i][m];
                }
            }
        }

    }
    return q;
}

void Roe_Riemann_Solver()
{
    double x_l = 0.0;
    double x_r = 1.0;
    int nx = 256;
    double dx = (x_r - x_l) / nx;

    double t = 0.2;
    double dt = 0.0001;
    int nt = ceil(t / dt);

    int ns = 20;   //Save ten sets of data results
    double ds = t / ns;          //The time interval for saving data

    vector<double> x(nx);
    vector < vector< vector<double> > > q(nx, vector<vector<double>>(3, vector<double>(ns + 1)));

    for (int i = 0; i < nx; i++)
    {
        x[i] = x_l + dx * i + 0.5 * dx;  //Assign node locations
    }

    q = numerical_Roe(nx, ns, nt, dx, dt, x, q);

    ofstream rho;
    ofstream u;
    ofstream e;
    rho.open("rho_final.csv", ios::out | ios::trunc);
    u.open("u_final.csv", ios::out | ios::trunc);
    e.open("e_final.csv", ios::out | ios::trunc);

    rho << "x" << ",";
    u << "x" << ",";
    e << "x" << ",";
    for (int i = 0; i < ns + 1; i++)
    {
        rho << i * ds << ",";
        u << i * ds << ",";
        e << i * ds << ",";
    }
    rho << endl;
    u << endl;
    e << endl;

    for (int i = 0; i < nx; i++)
    {
        rho << setprecision(16) << x[i] << ",";
        u << setprecision(16) << x[i] << ",";
        e << setprecision(16) << x[i] << ",";
        for (int j = 0; j < ns + 1; j++)
        {
            rho << setprecision(16) << q[i][0][j] << ",";
            u << setprecision(16) << q[i][1][j] / q[i][0][j] << ",";
            e << setprecision(16) << q[i][2][j] / q[i][0][j] << ",";
        }
        rho << endl;
        u << endl;
        e << endl;
    }

    rho.close();
    u.close();
    e.close();
}