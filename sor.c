#include "sor.h"
#include "parallel.h"
#include "helper.h"
#include <math.h>
#include <stdio.h>
#include <mpi.h>

void sor(
  double omg,
  double dx,
  double dy,
  int    imax,
  int    jmax,
  int imax_local,
  int jmax_local,
  int rank_l,
  int rank_r,
  int rank_b,
  int rank_t,
  double **P,
  double **RS,
  double *res,
  double *bufSend,
  double *bufRecv,
  MPI_Status *status,
  int chunk
) {
  int i,j;
  double rloc, rloc2;
  char message[80];
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
  /* SOR iteration */


  pressure_comm1(P, imax_local, jmax_local, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, status, chunk);
  for(i = 2; i <= imax_local; i++) {
	for(j = 2; j<=jmax_local; j++) {
	  P[i][j] = (1.0-omg)*P[i][j]
			  + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
	}
  }

  /* exchanges pressure between processes */
  //pressure_comm(P, imax_local, jmax_local, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, status, chunk);
  pressure_comm2(P, imax_local, jmax_local, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, status, chunk);


  /* compute the residual */
  rloc = 0;
  rloc2 = 0;
  *res = 0;
  for(i = 2; i <= imax_local; i++) {
    for(j = 2; j <= jmax_local; j++) {
      rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
              ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }
  /*rloc = rloc/(imax*jmax);
  rloc = sqrt(rloc);*/

  /* set residual */
  /* sums all local residuals and sends result to master */
  MPI_Reduce(&rloc, &rloc2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  rloc2 = sqrt(rloc2/(double)(imax*jmax));
  *res = rloc2;
  MPI_Bcast(res, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  /* set boundary values */
	/* bottom row */
	if(rank_b == MPI_PROC_NULL)
	{
		for(i = 0; i <= imax_local+1; i++)
		{
			P[i][1] = P[i][2];
		}
	}
	/* top row */
	if(rank_t == MPI_PROC_NULL)
	{
		for(i = 0; i <= imax_local+1; i++)
		{
			P[i][jmax_local+1] = P[i][jmax_local];
		}
	}
	/* left column */
	if(rank_l == MPI_PROC_NULL)
	{
		for(j = 0; j <= jmax_local+1; j++)
		{
			P[1][j] = P[2][j];
		}
	}
	/* right column */
	if(rank_r == MPI_PROC_NULL)
	{
		for(j = 0; j <= jmax_local+1; j++)
		{
			P[imax_local+1][j] = P[imax_local][j];
		}
	}
}

