#include "helper.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include "visual.h"
#include <stdio.h>

/* include MPI */
#include <mpi.h>
#include "parallel.h"

/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
int main(int argn, char** args) {

	char* szFileName;
	double Re,UI,VI,PI,GX,GY;
	double t_end,xlength,ylength;
	double dt,dx,dy;
	double alpha,omg,tau;
	double eps, dt_value,res, t;
	int itermax,it,imax,jmax, n, imax_local, jmax_local;
	int n_div;

	/* Output filename */
	char* szProblem;

	/* arrays */
	double **U;
	double **V;
	double **P;
	double **RS;
	double **F;
	double **G;


	/**** MPI parameters ****/
    MPI_Status status;
    int nbrOfProcesses; /* the number of MPI units */
    int rankOfMainProcess; /* the rank of the main unit */
  	char main_processor_name[MPI_MAX_PROCESSOR_NAME];  /* the main unit name */
  	int name_len; /* length of the main unit name */

    int iproc; /* x-direction domain decomposition */
    int jproc; /* y-direction domain decomposition */

    int omg_i; /* i-index of this subdomain*/
    int omg_j; /* j-index of this subdomain*/

    int il; /* x-left bound */
    int ir; /* x-right bound */
    int jb; /* j-bottom bound */
    int jt; /* j-top bound */

    int rank_l; /* left  domain ID */
    int rank_r; /* left  domain ID */
    int rank_b; /* right domain ID */
    int rank_t; /* top   domain ID */

    int isFollowingSorIteration; /* flag indicating if there is a next sor interation */
    char message[80];

    int bufferLength; /* maximum length of a subdomain side */
    double *sendBuffer; /* buffer for sending P and RS data */
    double *reicBuffer; /* buffer for receiving P and RS data */

    /**** Initialize MPI  ****/
    MPI_Init(&argn, &args);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankOfMainProcess);
    MPI_Comm_size(MPI_COMM_WORLD, &nbrOfProcesses);

    // Get the name of the processor
  	MPI_Get_processor_name(main_processor_name, &name_len);
	szProblem = "./cavity_output/Output";
	szFileName = "./cavity100.dat";
	read_parameters(szFileName, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength,
			&ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau,
			&itermax, &eps, &dt_value, &iproc, &jproc);

  	printf("NumProc: %d, iproc: %d, jproc: %d\n", nbrOfProcesses, iproc, jproc);
    /**** is the number of processes equal to the number of subdomains ? ****/
    if (nbrOfProcesses % (iproc*jproc) != 0) {
        if (rankOfMainProcess == 0)
        {
            printf("number of processes and subdomains must be equal!! \n");
        }
        MPI_Finalize();
        return -1;
    }

    /* MPI setup ok! */
	printf("MPI setup: %s, rank %d out of %d processors\n",
	         main_processor_name, rankOfMainProcess, nbrOfProcesses);

    /**** INIT PARALLEL METHOD: subdomain initialization ****/
    init_parallel(iproc,jproc,imax,jmax,&rankOfMainProcess,&il,&ir,&jb,&jt,&rank_l,
                  &rank_r,&rank_b,&rank_t,&omg_i,&omg_j,nbrOfProcesses);

    printf("process rank: %i matrix setup: il:%i ir:%i jb:%i jt:%i, rank_l: %d, rank_r: %d, rank_b: %d, rank_t: %d\n",
    		rankOfMainProcess,il, ir,jb,jt, rank_l, rank_r, rank_b, rank_t);


    /* set imax_local, jmax_local according to the subdomain (= last index of domain cell) */
    /* The first domain cell has index 2 */
    imax_local = ir - il + 2;
    jmax_local = jt - jb + 2;

	/*************** the algorithm (see section 5) *************************/

	/* allocate memory for the arrays */
	U = matrix(0, imax_local+1, 0, jmax_local+1);
	V = matrix(0, imax_local+1, 0, jmax_local+1);
	P = matrix(0, imax_local+1, 0, jmax_local+1);
	F = matrix(0, imax_local+1, 0, jmax_local+1);
	G = matrix(0, imax_local+1, 0, jmax_local+1);
	RS = matrix(0, imax_local+1, 0, jmax_local+1);

	/* Determine buffer length */
	bufferLength = imax_local + 3;
    if(bufferLength < jmax_local)
    	bufferLength = jmax_local + 3;

    /* allocate memory for the send and receive buffers */
    sendBuffer = (double *) malloc((size_t)(bufferLength * sizeof(double)));
    reicBuffer = (double *) malloc((size_t)(bufferLength * sizeof(double)));


	/* init uvp */
	init_uvp(UI, VI, PI, imax_local, jmax_local, U, V, P);

	t = 0.0;
	n = 0;

	while(t < t_end)
	{
		/*calculate the delta values*/
		calculate_dt(Re, tau, &dt, dx, dy, imax_local, jmax_local, U, V);
		/*calculate the boundary values*/
		boundaryvalues(imax_local, jmax_local, rank_l, rank_r, rank_b, rank_t, U, V);
		/*calculate F&G*/
		calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax_local, jmax_local, U, V, F, G, rank_l, rank_r, rank_b, rank_t);
		/*calculate righthand site*/
		calculate_rs(dt, dx, dy, imax_local, jmax_local, F, G, RS);

		/* set initial residual*/
		it = 0;
		res = eps + 1;

		isFollowingSorIteration=1;
		while (isFollowingSorIteration > 0)
		{
			sor(omg, dx, dy, imax, jmax, imax_local, jmax_local, rank_l, rank_r, rank_b, rank_t, P, RS, &res,
						sendBuffer, reicBuffer, &status, it);
			it++;

			/* master process determines the next iteration flag*/
		    if(rankOfMainProcess == 0)
		    {
		    	if(it < itermax && res > eps) isFollowingSorIteration = 1;
		    	else isFollowingSorIteration = 0;
		    }
		    MPI_Bcast(&isFollowingSorIteration, 1, MPI_INT, 0, MPI_COMM_WORLD);
/*		    sprintf(message, "isFollowingSor: %d, res: %f, eps: %f, t: %f", isFollowingSorIteration, res, eps, t);
		    Programm_Sync(message);*/
		}
		/* calculate uv */
		calculate_uv(dt, dx, dy, imax_local, jmax_local, U, V, F, G, P, rank_l, rank_r, rank_b, rank_t);
		uv_comm(U, V, imax_local, jmax_local, rank_l, rank_r,  rank_b, rank_t, sendBuffer, reicBuffer, &status, n);

		n_div = (dt_value / dt);
		if (n_div!=0 && n % n_div == 0)
		{
			output_vtk(U, V, P, il, ir, jb, jt, imax, jmax, omg_i, omg_j, iproc, jproc, dx, dy,
			                n, szProblem);
			Programm_Sync("Wrote output");

		}		t = t + dt;
		n++;
		if(rankOfMainProcess == 0)
		{

			printf("F\n");
			print_matrix(F, 1, imax_local+1, 1, jmax_local+1);
			/*printf("RS\n");
			print_matrix(RS, 1, imax_local+1, 1, jmax_local+1);
			printf("P (rank %d) (Sor-Its: %d, res: %f)", rankOfMainProcess, it, res);
			print_matrix(P, 1, imax_local+1, 1, jmax_local+1);*/
			printf("U\n");
			print_matrix(U, 1, imax_local+1, 1, jmax_local+1);
		}
		Programm_Sync("Dummy");
		if(rankOfMainProcess == 1)
		{
			printf("F\n");
			print_matrix(F, 1, imax_local+1, 1, jmax_local+1);
			/*printf("RS\n");
			print_matrix(RS, 1, imax_local+1, 1, jmax_local+1);
			printf("P (rank %d) (Sor-Its: %d, res: %f)", rankOfMainProcess, it, res);
			print_matrix(P, 1, imax_local+1, 1, jmax_local+1);*/
			printf("U\n");
			print_matrix(U, 1, imax_local+1, 1, jmax_local+1);
		}
		Programm_Sync("Dummy2");
		if(rankOfMainProcess == 2)
		{
			printf("F\n");
			print_matrix(F, 1, imax_local+1, 1, jmax_local+1);
			/*printf("RS\n");
			print_matrix(RS, 1, imax_local+1, 1, jmax_local+1);
			printf("P (rank %d) (Sor-Its: %d, res: %f)", rankOfMainProcess, it, res);
			print_matrix(P, 1, imax_local+1, 1, jmax_local+1);*/
			printf("U\n");
			print_matrix(U, 1, imax_local+1, 1, jmax_local+1);
		}

		sprintf(message, "n = %d, it = %d", n, it);
		Programm_Sync(message);
		if(n>1)
			break;
	}

	/* free arrays */
	free_matrix(U, 0, imax_local + 1, 0, jmax_local + 1);
	free_matrix(V, 0, imax_local + 1, 0, jmax_local + 1);
	free_matrix(P, 0, imax_local + 1, 0, jmax_local + 1);
	free_matrix(F, 0, imax_local + 1, 0, jmax_local + 1);
	free_matrix(G, 0, imax_local + 1, 0, jmax_local + 1);
	free_matrix(RS, 0, imax_local + 1, 0, jmax_local + 1);

	Programm_Stop("Finally...");
	MPI_Finalize();
	return -1;
}
