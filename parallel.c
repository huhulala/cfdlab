#include "parallel.h"


void Program_Message(char *txt)
/* produces a stderr text output  */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
}


void Programm_Sync(char *txt)
/* produces a stderr textoutput and synchronize all processes */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */  
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
}


void Programm_Stop(char *txt)
/* all processes will produce a text output, be synchronized and finished */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
   fprintf(stderr,"-STOP- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(1);
}

void init_parallel(
		int iproc,
        int jproc,
        int imax,
        int jmax,
        int *rank,
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
        int num_proc)
{
	int i = 0;
	int j = 0;

	int processCounter = 1;
	int bound = 0;
	int bottomDomainIdx = MPI_PROC_NULL;
	int leftDomainIdx   = MPI_PROC_NULL;
	int topDomainIdx    = MPI_PROC_NULL;
	int rightDomainIdx   = MPI_PROC_NULL;
	MPI_Status status;
/*	MPI_Request request;*/


printf("init parallel setup: \n");

	/* master process initialize everything here */
    if (*rank == 0)
    {
    	/*set process 0 values */
    	*omg_i=0;
    	*omg_j=0;

    	*rank_l = MPI_PROC_NULL ;
    	*rank_b = MPI_PROC_NULL ;

        /* top neighbor */
        if (j == jproc-1) *rank_t = MPI_PROC_NULL;    /*  Top row */
        else *rank_t = iproc;                        /*  Add one block width  */

        /* right neighbor */
        if (i == iproc-1) *rank_r = MPI_PROC_NULL;  /* Right column */
        else *rank_r =  1;


        /* bottom bound */
        *jb = (jmax / jproc) * j + 1;

        /* top bound */
        if (j == jproc-1) *jt = jmax;
        else *jt = (j + 1) * (jmax / jproc);

        /* left bound */
        *il = (imax / iproc) * i + 1;

        /* right bound */
        *ir =  (i + 1) * (imax / iproc);


    	/* loop through the subdomains and count the processes:
    	 * bottom left cell is zero */
        for (j = 0; j < jproc; ++j)
        for (i = 0; i < iproc; ++i)
        {
        	  if(i==0 && j ==0) continue;

        	  MPI_Send(&i, 1, MPI_INT, processCounter, OMG_I, MPI_COMM_WORLD);
        	  MPI_Send(&j, 1, MPI_INT, processCounter, OMG_J, MPI_COMM_WORLD);

        	  /* Determine neighborhood relations */
        	  /* bottom neighbor */
              if (j == 0) bottomDomainIdx = MPI_PROC_NULL;      /* Bottom row */
              else bottomDomainIdx = processCounter - iproc;    /* Subtract one block width */
              MPI_Send(&bottomDomainIdx, 1, MPI_INT, processCounter, RANK_B, MPI_COMM_WORLD);

              /* left neighbor */
              if (i == 0) leftDomainIdx = MPI_PROC_NULL;        /* Left column */
              else leftDomainIdx = processCounter - 1;          /* Subtract one  */
              MPI_Send(&leftDomainIdx, 1, MPI_INT, processCounter, RANK_L, MPI_COMM_WORLD);

              /* top neighbor */
              if (j == jproc-1) topDomainIdx = MPI_PROC_NULL;    /*  Top row */
              else topDomainIdx = processCounter + iproc;        /*  Add one block width  */
              MPI_Send(&topDomainIdx, 1, MPI_INT, processCounter, RANK_T, MPI_COMM_WORLD);

              /* right neighbor */
              if (i == iproc-1) rightDomainIdx = MPI_PROC_NULL;  /* Right column */
              else rightDomainIdx = processCounter + 1;         /* Add one  */
             MPI_Send(&rightDomainIdx, 1, MPI_INT, processCounter, RANK_R, MPI_COMM_WORLD);


              /* Determine bounds  */
              /* (bikini) bottom bound */
              bound = (jmax / jproc) * j + 1;
              MPI_Send(&bound, 1, MPI_INT, processCounter, BOUNDB, MPI_COMM_WORLD);

              /* top bound */
              if (j == jproc-1) bound = jmax;
              else bound = (j + 1) * (jmax / jproc);
              MPI_Send(&bound, 1, MPI_INT, processCounter, BOUNDT, MPI_COMM_WORLD);

              /* left bound */
              bound = (imax / iproc) * i + 1;
              MPI_Send(&bound, 1, MPI_INT, processCounter, BOUNDL, MPI_COMM_WORLD);

              /* right bound */
              if (i == iproc-1) bound = imax;
              else bound = (i + 1) * (imax / iproc);
              MPI_Send(&bound, 1, MPI_INT, processCounter, BOUNDR, MPI_COMM_WORLD);
              ++processCounter;
        }
    }
    else
    {


    	/* Receive all the signals */
    	MPI_Recv(omg_i, 1, MPI_INT, 0, OMG_I, MPI_COMM_WORLD, &status);
    	MPI_Recv(omg_j, 1, MPI_INT, 0, OMG_J, MPI_COMM_WORLD, &status);

    	MPI_Recv(rank_b, 1, MPI_INT, 0, RANK_B, MPI_COMM_WORLD, &status);
    	MPI_Recv(rank_r, 1, MPI_INT, 0, RANK_R, MPI_COMM_WORLD, &status);
    	MPI_Recv(rank_t, 1, MPI_INT, 0, RANK_T, MPI_COMM_WORLD, &status);
    	MPI_Recv(rank_l, 1, MPI_INT, 0, RANK_L, MPI_COMM_WORLD, &status);

    	MPI_Recv(jb, 1, MPI_INT, 0, BOUNDB, MPI_COMM_WORLD, &status);
    	MPI_Recv(jt, 1, MPI_INT, 0, BOUNDT, MPI_COMM_WORLD, &status);
    	MPI_Recv(il, 1, MPI_INT, 0, BOUNDL, MPI_COMM_WORLD, &status);
    	MPI_Recv(ir, 1, MPI_INT, 0, BOUNDR, MPI_COMM_WORLD, &status);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

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
                   int chunk)
{
    int i;
    int j;
    /*** SEND LEFT, receive from right ***/
    /*copying data to be send prom pressure matrix to send buffer*/
    /*if there is no left neighbour there is not need to copy data, since it will not be sent anyway*/
    if(rank_l != MPI_PROC_NULL)
    {
        for(j = 1; j <= jmax; j++)
        {
            bufSend[j-1] = P[2][j];
        }
    }
    /*trasfering data*/
    MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_l, chunk, bufRecv, jmax, MPI_DOUBLE, rank_r, chunk, MPI_COMM_WORLD, status);
    /*copying received data from receive buffer to pressure matric*/
    if(rank_r != MPI_PROC_NULL) /*if there is no right neighbour there is no need to copy data. Since this is a system boundary it will not be changed anyway*/
    {
        for(j = 1; j <= jmax; j++)
        {
            P[imax+1][j] = bufRecv[j-1];
        }
    }

    /*** SEND RIGHT, receive from left ***/
    if(rank_r != MPI_PROC_NULL)
    {
        for(j = 1; j <= jmax; j++)
        {
            bufSend[j-1] = P[imax][j];
        }
    }
    MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_r, chunk, bufRecv, jmax, MPI_DOUBLE, rank_l, chunk, MPI_COMM_WORLD, status);
    if(rank_l != MPI_PROC_NULL){
        for(j = 1; j <= jmax; j++)
        {
            P[1][j] = bufRecv[j-1];
        }
    }

    /*** SEND UP, receive from bottom ***/
    if(rank_t != MPI_PROC_NULL)
    {
        for(i = 1; i <= imax; i++)
        {
            bufSend[i-1] = P[i][jmax];
        }
    }
    MPI_Sendrecv(bufSend, imax, MPI_DOUBLE, rank_t, chunk, bufRecv, imax, MPI_DOUBLE, rank_b, chunk, MPI_COMM_WORLD, status);
    if(rank_b != MPI_PROC_NULL)
    {
        for(i = 1; i <= imax; i++)
        {
            P[i][1] = bufRecv[i-1];
        }
    }

    /*** SEND DOWN, receive from up ***/
    if(rank_b != MPI_PROC_NULL)
    {
        for(i = 1; i <= imax; i++)
        {
            bufSend[i-1] = P[i][2];
        }
    }
    MPI_Sendrecv(bufSend, imax, MPI_DOUBLE, rank_b, chunk, bufRecv, imax, MPI_DOUBLE, rank_t, chunk, MPI_COMM_WORLD, status);
    if(rank_t != MPI_PROC_NULL)
    {
        for(i = 1; i <= imax; i++)
        {
            P[i][jmax+1] = bufRecv[i-1];
        }
    }
}


void pressure_comm1(double **P, int imax, int jmax, int rank_l, int rank_r, int rank_b, int rank_t,
                   double *bufSend, double *bufRecv, MPI_Status *status, int chunk)
{
	/* Computes Exchange into the right and top direction */
    int i;
    int j;

    if(rank_l != MPI_PROC_NULL){
    	MPI_Recv(bufRecv, jmax-1, MPI_DOUBLE, rank_l, chunk, MPI_COMM_WORLD, status);
        for(j = 2; j <= jmax; j++)
        {
            P[1][j] = bufRecv[j-2];
        }
    }
    if(rank_b != MPI_PROC_NULL)
    {
    	MPI_Recv(bufRecv, imax-1, MPI_DOUBLE, rank_b, chunk, MPI_COMM_WORLD, status);
    	for(i = 2; i <= imax; i++)
        {
            P[i][1] = bufRecv[i-2];
        }
    }

}


void pressure_comm2(double **P, int imax, int jmax, int rank_l, int rank_r, int rank_b, int rank_t,
                   double *bufSend, double *bufRecv, MPI_Status *status, int chunk)
{
	int i,j;
	/* Send all other values not sended in comm1 */

    if(rank_r != MPI_PROC_NULL)
    {
        for(j = 2; j <= jmax; j++)
        {
            bufSend[j-2] = P[imax][j];
        }
        MPI_Send(bufSend, jmax-1, MPI_DOUBLE, rank_r, chunk, MPI_COMM_WORLD);
    }

    if(rank_t != MPI_PROC_NULL)
    {
        for(i = 2; i <= imax; i++)
        {
            bufSend[i-2] = P[i][jmax];
        }
        MPI_Send(bufSend, imax-1, MPI_DOUBLE, rank_t, chunk, MPI_COMM_WORLD);
    }

    if(rank_r != MPI_PROC_NULL)
    {
    	MPI_Recv(bufRecv, jmax-1, MPI_DOUBLE, rank_r, chunk, MPI_COMM_WORLD, status);
		for(j = 2; j <= jmax; j++)
		{
			P[imax+1][j] = bufRecv[j-2];
		}
    }

    if(rank_l != MPI_PROC_NULL)
    {
        for(j = 2; j <= jmax; j++)
        {
            bufSend[j-2] = P[2][j];
        }
        MPI_Send(bufSend, jmax-1, MPI_DOUBLE, rank_l, chunk, MPI_COMM_WORLD);
    }

    if(rank_t != MPI_PROC_NULL)
    {
    	MPI_Recv(bufRecv, imax-1, MPI_DOUBLE, rank_t, chunk, MPI_COMM_WORLD, status);
        for(i = 2; i <= imax; i++)
        {
            P[i][jmax+1] = bufRecv[i-2];
        }
    }

    if(rank_b != MPI_PROC_NULL)
    {
        for(i = 2; i <= imax; i++)
        {
            bufSend[i-2] = P[i][2];
        }
        MPI_Send(bufSend, imax-1, MPI_DOUBLE, rank_b, chunk, MPI_COMM_WORLD);
    }
}

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
                   int chunk)
{
    int i, j;
    /******************** Velocity U ********************/

    /* Send left, receive right */

    if (rank_l != MPI_PROC_NULL){
        for (j = 2; j <= jmax; j++)
        {
            bufSend[j-2] = U[2][j];
        }
    }
    MPI_Sendrecv(bufSend, jmax-1, MPI_DOUBLE, rank_l, chunk, bufRecv, jmax-1, MPI_DOUBLE, rank_r, chunk, MPI_COMM_WORLD, status);
    if (rank_r != MPI_PROC_NULL){
        for (j = 2; j <= jmax; j++){
             U[imax+1][j]=bufRecv[j-2];
        }
    }

     /* Send right, receive left */

    if (rank_r != MPI_PROC_NULL){
        for (j = 2; j <= jmax; j++)
        {
            bufSend[j-2] = U[imax-1][j];
        }
    }
    MPI_Sendrecv(bufSend, jmax-1, MPI_DOUBLE, rank_r, chunk, bufRecv, jmax-1, MPI_DOUBLE, rank_l, chunk, MPI_COMM_WORLD, status);
    if (rank_l !=MPI_PROC_NULL){
        for (j= 2;j <= jmax;j++){
             U[0][j]=bufRecv[j-2];
        }
    }

     /* Send up, receive bottom */
    if (rank_t != MPI_PROC_NULL){
        for (i = 1; i <= imax; i++)
        {
            bufSend[i-1] = U[i][jmax];
        }
    }
       MPI_Sendrecv(bufSend, imax, MPI_DOUBLE, rank_t, chunk, bufRecv, imax, MPI_DOUBLE, rank_b, chunk, MPI_COMM_WORLD, status);
    if (rank_b !=MPI_PROC_NULL){
        for (i = 1; i <= imax; i++){
             U[i][1]=bufRecv[i-1];
        }
    }

    /* Send to bottom, receive from top */
    if (rank_b != MPI_PROC_NULL){
        for (i = 1; i <= imax; i++)
        {
            bufSend[i-1] = U[i][2];
        }
    }

        MPI_Sendrecv(bufSend, imax, MPI_DOUBLE, rank_b, chunk, bufRecv, imax, MPI_DOUBLE, rank_t, chunk, MPI_COMM_WORLD, status);

    if (rank_t !=MPI_PROC_NULL){
        for (i = 1; i <= imax;i++){
             U[i][jmax+1]=bufRecv[i-1];
        }
    }


     /********************* Velocity V *************************/

    /* Send-Recv from Right to Left  */

    if (rank_l != MPI_PROC_NULL){
        for (j = 1; j <= jmax; j++)
        {
            bufSend[j-1] = V[2][j];
        }
    }

        MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_l, chunk, bufRecv, jmax, MPI_DOUBLE, rank_r, chunk, MPI_COMM_WORLD, status);

    if (rank_r !=MPI_PROC_NULL){
        for (j = 1; j <= jmax; j++){
             V[imax+1][j]=bufRecv[j-1];
        }
    }

     /* Send-Recv from left to right */

    if (rank_r != MPI_PROC_NULL){
        for (j = 1; j <= jmax; j++)
        {
            bufSend[j-1] = V[imax][j];
        }
    }

        MPI_Sendrecv(bufSend, jmax, MPI_DOUBLE, rank_r, chunk, bufRecv, jmax, MPI_DOUBLE, rank_l, chunk, MPI_COMM_WORLD, status);

    if (rank_l !=MPI_PROC_NULL){
        for (j = 1 ; j <= jmax; j++){
             V[1][j] = bufRecv[j-1];
        }
    }


     /* Send-Recv from bottom to top */

    if (rank_t != MPI_PROC_NULL){
        for (i = 2; i <= imax; i++)
        {
            bufSend[i-2] = V[i][jmax-1];
        }
    }

        MPI_Sendrecv(bufSend, imax-1, MPI_DOUBLE, rank_t, chunk, bufRecv, imax-1, MPI_DOUBLE, rank_b, chunk, MPI_COMM_WORLD, status);

    if (rank_b != MPI_PROC_NULL){
        for (i= 2;i <= imax;i++){
             V[i][0] = bufRecv[i-2];
        }
    }

    /* Send-Recv Bot-Top */

    if (rank_b != MPI_PROC_NULL){
        for (i= 2; i <= imax; i++)
        {
            bufSend[i-2] = V[i][2];
        }
    }

        MPI_Sendrecv(bufSend, imax-1, MPI_DOUBLE, rank_b, chunk, bufRecv, imax-1, MPI_DOUBLE, rank_t, chunk, MPI_COMM_WORLD, status);

    if (rank_t !=MPI_PROC_NULL){
        for (i = 2; i <= imax; i++){
             V[i][jmax+1] = bufRecv[i-2];
        }
    }
}

