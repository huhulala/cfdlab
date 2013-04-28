#include "uvp.h"
#include <math.h>
#include "fd.h"

void calculate_fg(double Re, double GX, double GY, double alpha, double dt,
		double dx, double dy, int imax, int jmax, double **U, double **V,
		double **F, double **G) {

	/* see formula 9 and 10 in combination with formulas 4 and 5 */
	int i;
	int j;

	for (j = 1; j <= jmax; j++)
	for (i = 1; i <= imax; i++)
	{
		/********** calculate F **********/
		if (i <= imax - 1)
			F[i][j] = U[i][j] + dt * (
			/* 1/Re * (d²u/dx² + d²u/dy²) */
			1 / Re * (d2udx2(i, j, U, dx) + d2udy2(i, j, U, dy))
			/* - du²/dx */
			- du2dx(i, j, U, dx, alpha)
			/* - duv/dy */
			- duvdy(i, j, U, V, dy, alpha) + GX);

		/********** calculate G **********/
		if (j <= jmax - 1)
				G[i][j] = V[i][j] + dt * (
				/* 1/Re * (d²v/dx² + d²v/dy²) */
				1 / Re * (d2vdx2(i, j, V, dx) + d2vdy2(i, j, V, dy))
				/* - duv/dx */
				- duvdx(i, j, U, V, dx, dy, alpha)
				/* - dv²/dy */
				- dv2dy(i, j, V, dy, alpha) + GY);
	}

	/* calculate boundary values -  see formula 17 */
	for (j = 1; j <= jmax; j++)
	{
		F[0][j] = U[0][j];
		F[imax][j] = U[imax][j];
	}
	for (i = 1; i <= imax; i++)
	{
		G[i][0] = V[i][0];
		G[i][jmax] = V[i][jmax];
	}
}

void calculate_rs(double dt, double dx, double dy, int imax, int jmax,
		double **F, double **G, double **RS) {

	/* right hand side of formula 11 */
	int i;
	int j;

	for (j = 1; j <= jmax; j++)
	for (i = 1; i <= imax; i++)
	{
		RS[i][j] = 1 / dt * ((F[i][j] - F[i - 1][j]) / dx + (G[i][j]
					- G[i][j - 1]) / dy);
	}
}

void calculate_dt(double Re, double tau, double *dt, double dx, double dy,
		int imax, int jmax, double **U, double **V) {
	/* See formula 13 */
	double umax = fabs(U[1][1]);
	double vmax = fabs(V[1][1]);
    double dtcon, dxcon, dycon;
    double min;
	int i;
	int j;

	for (j = 1; j <= jmax; j++)
	for (i = 1; i <= imax; i++)
	{
		if (fabs(U[i][j]) > umax)umax = fabs(U[i][j]);
		if (fabs(V[i][j]) > vmax)vmax = fabs(V[i][j]);
	}

	/* conditions */
	dtcon = Re/(2*(1/(dx*dx) + 1/(dy*dy)));
	dxcon = dx/fabs(umax);
	dycon = dy/fabs(vmax);

	/* determine smalles condition */
    min = dtcon;
	if(min > dxcon)min = dxcon;
	if(min > dycon) min = dycon;

	/* calculate dt */
	*dt = tau * min;
}

void calculate_uv(double dt, double dx, double dy, int imax, int jmax,
		double **U, double **V, double **F, double **G, double **P) {
	/* see formula 7 and 8 */
	int i,j;
	double dtdx = dt / dx;
	double dtdy = dt / dy;

	for (j = 1; j <= jmax; j++)
	for (i = 1; i <= imax; i++)
	{
		if (i <= imax - 1)U[i][j] = F[i][j] - dtdx * (P[i + 1][j] - P[i][j]);
		if (j <= jmax - 1)V[i][j] = G[i][j] - dtdy * (P[i][j + 1] - P[i][j]);
	}
}
