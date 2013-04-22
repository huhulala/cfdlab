#include "boundary_val.h"

void boundaryvalues(int imax, int jmax, double **U, double **V)
{
	int i;
	int j;
	for(i=0; i<=imax; i++)
	{
		/* upper and lower bounder with corners */
		U[i][0] = 1;
		U[i][jmax] = 0.0;
		V[i][0] = 0.0;
		V[i][jmax] = 0.0;
	}
	for(j=0; j<=jmax; j++)
	{
		/* left and right border without corners */
		U[0][j] = 0.0;
		U[imax][j] = 0.0;
		V[0][j] = 0.0;
		V[imax][j] = 0.0;
	}
}

