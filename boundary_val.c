#include "boundary_val.h"

void boundaryvalues(int imax, int jmax, double **U, double **V)
{
	/* See formula 14 and 15 */
	int i;
	int j;
	for(i=1; i<=imax; i++)
	{
		/* upper and lower bounder without corners */
		U[i][0] = -U[i][1];
		/* U[i][jmax+1] = -U[i][jmax]; Moving Wall instead of no-slip TODO: Correct? */
		U[i][jmax+1] = -U[i][jmax] + 2;
		V[i][0] = 0.0;
		V[i][jmax] = 0.0;
	}
	for(j=1; j<=jmax; j++)
	{
		/* left and right border without corners */
		U[0][j] = 0.0;
		U[imax][j] = 0.0;
		V[0][j] = -V[1][j];
		V[imax+1][j] = -V[imax][j];
	}
}
