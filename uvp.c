#include "uvp.h"
#include <math.h>
#include "fd.h"
#include <mpi.h>

void calculate_fg(double Re, double GX, double GY, double alpha, double dt,
		double dx, double dy, int imax, int jmax, double **U, double **V,
		double **F, double **G, int rank_l, int rank_r, int rank_b, int rank_t) {

	/* see formula 9 and 10 in combination with formulas 4 and 5 */
	int i;
	int j;
	int istart, jstart;

	istart = 1;
	jstart = 1;

	if(rank_l == MPI_PROC_NULL)
		istart = 2;
	if(rank_b == MPI_PROC_NULL)
		jstart = 2;

	for (j = jstart; j <= jmax; j++)
	for (i = istart; i <= imax; i++)
	{
		/********** calculate F **********/
		if (!(rank_r == MPI_PROC_NULL && i > imax-1)) /* not a bounder */
		{
				F[i][j] = U[i][j] + dt * (
				/* 1/Re * (d²u/dx² + d²u/dy²) */
				1 / Re * (d2udx2(i, j, U, dx) + d2udy2(i, j, U, dy))
				/* - du²/dx */
				- du2dx(i, j, U, dx, alpha)
				/* - duv/dy */
				- duvdy(i, j, U, V, dy, alpha) + GX);
		}
		/********** calculate G **********/
		if (!(rank_t == MPI_PROC_NULL && j > jmax-1))
		{
				G[i][j] = V[i][j] + dt * (
				/* 1/Re * (d²v/dx² + d²v/dy²) */
				1 / Re * (d2vdx2(i, j, V, dx) + d2vdy2(i, j, V, dy))
				/* - duv/dx */
				- duvdx(i, j, U, V, dx, dy, alpha)
				/* - dv²/dy */
				- dv2dy(i, j, V, dy, alpha) + GY);
		}
	}

	/* calculate boundary values -  see formula 17 */
	if(rank_l == MPI_PROC_NULL)
	{
		for (j = 1; j <= jmax; j++)
		{
			F[1][j] = U[1][j];
		}
	}
	if(rank_r == MPI_PROC_NULL)
	{
		for(j = 1; j <= jmax; j++)
			F[imax][j] = U[imax][j];
	}
	if(rank_b == MPI_PROC_NULL)
	{
		for (i = 1; i <= imax; i++)
			G[i][1] = V[i][1];
	}
	if(rank_t == MPI_PROC_NULL)
	{
		for(i = 1; i <= imax; i++)
			G[i][jmax] = V[i][jmax];
	}
}

void calculate_rs(double dt, double dx, double dy, int imax, int jmax,
		double **F, double **G, double **RS) {

	/* right hand side of formula 11 */
	int i;
	int j;

	for (j = 2; j <= jmax; j++)
	for (i = 2; i <= imax; i++)
	{
		RS[i][j] = 1 / dt * ((F[i][j] - F[i - 1][j]) / dx + (G[i][j]
					- G[i][j - 1]) / dy);
	}
}

void calculate_dt(double Re, double tau, double *dt, double dx, double dy,
		int imax, int jmax, double **U, double **V) {

	/* See formula 13 */
	double umax = fabs(U[1][2]);
	double vmax = fabs(V[2][1]);
    double global_Umax = 0;
    double global_Vmax = 0;
    double dtcon, dxcon, dycon;
    double min;
	int i;
	int j;

	for (j = 1; j <= jmax; j++)
	for (i = 1; i <= imax; i++)
	{
		/* Check all computed values of the domain */
		if (fabs(U[i][j]) > umax && j>1) umax = fabs(U[i][j]);
		if (fabs(V[i][j]) > vmax && i>1) vmax = fabs(V[i][j]);
	}

    /* determines global umax and vmax with MPI_REDUCE and send it to the master threa */
    MPI_Reduce(&umax, &global_Umax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&vmax, &global_Vmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    /* broadcasts umax and vmax */
    MPI_Bcast(&global_Umax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&global_Vmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	/* now use global max conditions */
	dtcon = Re/(2*(1/(dx*dx) + 1/(dy*dy)));
	dxcon = dx/fabs(global_Umax);
	dycon = dy/fabs(global_Vmax);

	/* determine smalles condition */
    min = dtcon;
	if(min > dxcon)min = dxcon;
	if(min > dycon) min = dycon;

	/* calculate dt */
	*dt = tau * min;
}

void calculate_uv(double dt, double dx, double dy, int imax, int jmax,
		double **U, double **V, double **F, double **G, double **P,
		int rank_l, int rank_r, int rank_b, int rank_t) {
	/* see formula 7 and 8 */
	int i,j;
	int istart, jstart;
	double dtdx = dt / dx;
	double dtdy = dt / dy;


	istart = 1;
	jstart = 1;

	if(rank_l == MPI_PROC_NULL)
		istart = 2;
	if(rank_b == MPI_PROC_NULL)
		jstart = 2;

	/* Calculate like formulas 7 - 10 on sheet 1 with shifted indices */
	for (j = jstart; j <= jmax; j++)
	for (i = istart; i <= imax; i++)
	{
		if (j > 1 && !(rank_r == MPI_PROC_NULL && i > imax-1)) U[i][j] = F[i][j] - dtdx * (P[i + 1][j] - P[i][j]);
		if (i > 1 && !(rank_t == MPI_PROC_NULL && j > jmax-1)) V[i][j] = G[i][j] - dtdy * (P[i][j + 1] - P[i][j]);
	}
}
