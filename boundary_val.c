#include "boundary_val.h"
#include <mpi.h>

void boundaryvalues(int imax, int jmax, int rank_l, int rank_r, int rank_b, int rank_t, double **U, double **V)
{
	/* See formula 14 and 15, set bounders only if domain has bounders */
	int i;
	int j;
	for(i=1; i<=imax+1; i++)
	{
		/* upper and lower bounder without corners */
		if(rank_b == MPI_PROC_NULL)
		{
			U[i][1] = -U[i][2];
			V[i][1] = 0.0;
		}
		/* U[i][jmax+1] = -U[i][jmax]; Moving Wall instead of no-slip */
		if(rank_t == MPI_PROC_NULL)
		{
			U[i][jmax+1] = -U[i][jmax] + 2;
			V[i][jmax] = 0.0;
		}
	}
	for(j=1; j<=jmax+1; j++)
	{
		/* left and right border without corners */
		if(rank_l == MPI_PROC_NULL)
		{
			U[1][j] = 0.0;
			V[1][j] = -V[2][j];
		}
		if(rank_r == MPI_PROC_NULL)
		{
			U[imax][j] = 0.0;
			V[imax+1][j] = -V[imax][j];
		}
	}
}
