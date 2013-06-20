#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

/*  MPI message flags */
/* processes indices */
#define OMG_I   1
#define OMG_J   2
/* neighborhood flags */
#define RANK_L    3
#define RANK_R    4
#define RANK_T    5
#define RANK_B    6

#define BOUNDL 7
#define BOUNDR 8
#define BOUNDT 9
#define BOUNDB 10

/* produces a stderr text output  */
void Program_Message(char *txt);

/* produces a stderr textoutput and synchronize all processes */
void Programm_Sync(char *txt);

/* all processes will produce a text output, be synchronized and finished */
void Programm_Stop(char *txt);



/* split the domain  */
void init_parallel(
	int iproc,
    int jproc,
    int imax,
    int jmax,
    int *myrank,
    int *il,
    int *ir,
    int *jb,
    int *jt,
    int *rank_l,
    int *rank_r,
    int *rank_b,
    int *rank_t,
    int *omg_i,
    int *omg_j,
    int num_proc);

/* Exchanges pressure values after sor interation between neighboring processes */
/* Exchanges pressure between nighbouring processes */
void pressure_comm(double **P,
				   int imax,
				   int jmax,
                   int rank_l,
                   int rank_r,
                   int rank_b,
                   int rank_t,
                   double *bufSend,
                   double *bufRecv,
                   MPI_Status *status,
                   int chunk);
/* Before computing */
void pressure_comm1(double **P,
				   int imax,
				   int jmax,
                   int rank_l,
                   int rank_r,
                   int rank_b,
                   int rank_t,
                   double *bufSend,
                   double *bufRecv,
                   MPI_Status *status,
                   int chunk);
/* After calculating P */
void pressure_comm2(double **P,
				   int imax,
				   int jmax,
                   int rank_l,
                   int rank_r,
                   int rank_b,
                   int rank_t,
                   double *bufSend,
                   double *bufRecv,
                   MPI_Status *status,
                   int chunk);
/* Exchanges velocity values after sor interation between neighboring processes */
void uv_comm(double **U,
				   double ** V,
				   int imax,
				   int jmax,
                   int rank_l,
                   int rank_r,
                   int rank_b,
                   int rank_t,
                   double *bufSend,
                   double *bufRecv,
                   MPI_Status *status,
                   int chunk);

/* Exchanges velocities between nighbouring processes */
void uv_comm1(double **U,
              double **V,
              int il,
              int ir,
              int jb,
              int jt,
              int rank_l,
              int rank_r,
              int rank_b,
              int rank_t,
              double *bufSend,
              double *bufRecv,
              MPI_Status *status,
              int chunk);
