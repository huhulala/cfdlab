#include "uvp.h"
#include <math.h>

void calculate_fg(double Re, double GX, double GY, double alpha, double dt,
		double dx, double dy, int imax, int jmax, double **U, double **V,
		double **F, double **G)
{
	/* See formular 9 and 10 in combination with formulars 4 and 5 */
	int i;
	int j;

	for(i=1; i<=imax; i++)
	{
		for(j=1; j<=jmax; j++)
		{
			if(i<=imax-1)
				F[i][j] = U[i][j] + dt * (
						/* 1/Re * (d²u/dx² + d²u/dy²) */
						1/Re * (1/dx * (U[i-1][j] - 2*U[i][j] + U[i+1][j])/(dx*dx) +
								(U[i][j-1] - 2*U[i][j] + U[i][j+1])/(dy*dy))
						/* - du²/dx */
						- 1/dx*(((U[i][j] + U[i+1][j])/2)*((U[i][j] + U[i+1][j])/2) -
								((U[i-1][j] + U[i][j])/2)*((U[i-1][j] + U[i][j])/2)) +
								alpha/dx * (fabs(U[i][j] +U[i+1][j])/2*(U[i][j]-U[i+1][j])/2 -
										fabs(U[i-1][j]+U[i][j])/2*(U[i-1][j]-U[i][j])/2)
						/* - duv/dy */
						- 1/dy*((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])/4 -
								(V[i][j-1] + V[i+1][j-1])*(U[i][j-1]+U[i][j])/4) +
								alpha/dy * (fabs(V[i][j] +V[i+1][j])/2*(U[i][j]-U[i][j+1])/2 -
										fabs(V[i][j-1]+V[i+1][j-1])/2*(U[i][j-1]-U[i][j])/2)
						+ GX);
			if(j<=jmax-1)
				G[i][j] = V[i][j] + dt * (
						/* 1/Re * (d²v/dx² + d²v/dy²) */
						1/Re *(1/dx * (V[i-1][j] - 2*V[i][j] + V[i+1][j])/(dx*dx) +
								(V[i][j-1] - 2*V[i][j] + V[i][j+1])/(dy*dy))
						/* - duv/dx */
						- 1/dx*((U[i][j] + U[i][j+1])*(V[i][j] + V[i+1][j])/4 -
								(U[i-1][j] + U[i-1][j+1])*(V[i-1][j] + V[i][j])/4) +
								alpha/dx * (fabs(U[i][j] + U[i][j+1])*(V[i][j] - V[i+1][j])/4 -
										fabs(U[i-1][j] + U[i-1][j+1])*(V[i-1][j] - V[i][j])/4)
						/* - dv²/dy */
						- 1/dy*(((V[i][j] + V[i][j+1])/2)*((V[i][j] + V[i][j+1])/2) -
								((V[i][j-1] + V[i][j])/2)*((V[i][j-1] + V[i][j])/2)) +
								alpha/dy * (fabs(V[i][j] + V[i][j+1])*(V[i][j] - V[i][j+1])/4 -
										fabs(V[i][j-1] + V[i][j])*(V[i][j-1] - V[i][j])/4)
						+ GY);
		}
	}
}


void calculate_rs(double dt, double dx, double dy, int imax, int jmax,
		double **F, double **G, double **RS)
{
	/* Right hand side of formular 11 */
	int i;
	int j;

	for(i=1; i<=imax; i++)
	{
		for(j=1; j<=jmax; j++)
		{
			RS[i][j] = 1/dt * ((F[i][j]-F[i-1][j])/dx + (G[i][j]-G[i][j-1])/dy);
		}
	}
}


void calculate_dt(double Re, double tau, double *dt, double dx, double dy,
		int imax, int jmax, double **U, double **V)
{
	/* See formular 13 */
	double umax = U[0][0];
	double vmax = V[0][0];
	double dt1;
	double dt2;
	double dt3;
	int i;
	int j;
	for(i=1; i<=imax; i++)
	{
		for(j=1; j<=jmax; j++)
		{
			if(U[i][j] > umax)
				umax = U[i][j];
			if(V[i][j] > vmax)
				vmax = V[i][j];
		}
	}

	dt1 = Re/2 / (1/(dx*dx) + 1/(dy*dy));
	dt2 = dx/fabs(umax);
	dt3 = dy/fabs(vmax);

	if(dt1 > dt2)
		if(dt3 > dt1)
			*dt = tau * dt3;
		else
			*dt = tau * dt1;
	else
		if(dt3 > dt2)
			*dt = tau * dt3;
		else
			*dt = tau * dt2;
}


void calculate_uv(double dt, double dx, double dy, int imax, int jmax,
		double **U, double **V, double **F, double **G, double **P)
{
	/* See formular 7 and 8 */
	int i;
	int j;

	for(i=1; i<=imax; i++)
	{
		for(j=1; j<=jmax; j++)
		{
			if(i<=imax-1)
				U[i][j] = F[i][j] - dt/dx*(P[i+1][j] - P[i][j]);
			if(j<=jmax-1)
				V[i][j] = G[i][j] - dt/dy*(P[i][j+1] - P[i][j]);
		}
	}
}
