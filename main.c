#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include <stdio.h>

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
	double Re;
	double UI;
	double VI;
	double PI;
	double GX;
	double GY;
	double t_end;
	double xlength;
	double ylength;
	double dt;
	double dx;
	double dy;
	int imax;
	int jmax;
	double alpha;
	double omg;
	double tau;
	int itermax;
	double eps;
	double dt_value;

	int it;
	double res;
	double t;
	int n;
	int i,j;

	/* Output file */
	char* szProblem;

	/* Arrays */
	double **U;
	double **V;
	double **P;
	double **RS;
	double **F;
	double **G;

	szProblem = "./Out/Output";
	/* Algorithm (see section 5) */
	szFileName = "./cavity100.dat";
	read_parameters(szFileName, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength,
			&ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau,
			&itermax, &eps, &dt_value);
	t = 0.0;
	n = 0;

	/* allocate memory for the arrays */
	U = matrix(0, imax + 1, 0, jmax + 1);
	V = matrix(0, imax + 1, 0, jmax + 1);
	P = matrix(0, imax + 1, 0, jmax + 1);
	F = matrix(0, imax + 1, 0, jmax + 1);
	G = matrix(0, imax + 1, 0, jmax + 1);
	RS = matrix(0, imax + 1, 0, jmax + 1);

	init_uvp(UI, VI, PI, imax, jmax, U, V, P);

	while (t < t_end && n < 300) /* TODO n-Schranke wieder entfernen! */
	{
		/*calculate the delta values*/
		calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V);
		/*calculate the boundary values*/
		boundaryvalues(imax, jmax, U, V);
		/*calculate F&G*/
		calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G);
		/*calculate righthand site*/
		calculate_rs(dt, dx, dy, imax, jmax, F, G, RS);

		/* set initial residual*/
		it = 0;
		res = eps + 1;

		while (it < itermax && res > eps) {

			for (i = 1; i <= imax; i++) {
				P[i][0] = P[i][1];
				P[i][jmax + 1] = P[i][jmax];
			}

			for (j = 1; j <= jmax; j++) {
				P[0][j] = P[1][j];
				P[imax + 1][j] = P[imax][j];
			}

			sor(omg, dx, dy, imax, jmax, P, RS, &res);
			it++;
		}
		calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P);
		write_vtkFile(szProblem, n, xlength, ylength, imax, jmax, dx, dy, U, V,
				P);
		t = t + dt;
		n++;
	}

	write_vtkFile(szProblem, n, xlength, ylength, imax, jmax, dx, dy, U, V,
				P);

	free_matrix(U, 0, imax + 1, 0, jmax + 1);
	free_matrix(V, 0, imax + 1, 0, jmax + 1);
	free_matrix(P, 0, imax + 1, 0, jmax + 1);
	free_matrix(F, 0, imax + 1, 0, jmax + 1);
	free_matrix(G, 0, imax + 1, 0, jmax + 1);
	free_matrix(RS, 0, imax + 1, 0, jmax + 1);

	return -1;
}
