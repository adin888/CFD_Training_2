#include"../Header/cfd_head.h"

vector< vector<double> > fluxes(int nx, double gamma, vector< vector<double> > q, vector< vector<double> > f)
{
	double p;
	for (int i = 0; i < nx + 1; i++)
	{
		p = (gamma - 1.0) * (q[i][2] - 0.5 * q[i][1] * q[i][1] / q[i][0]);
		f[i][0] = q[i][1];
		f[i][1] = q[i][1] * q[i][1] / q[i][0] + p;
		f[i][2] = q[i][1] * q[i][2] / q[i][0] + p * q[i][1] / q[i][0];
	}
	return f;
}