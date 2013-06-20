#ifndef __HELPER_H__
#define __HELPER_H__

/* includefiles */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>

#ifdef PI			   /* to avoid problems with Linux */
#undef PI
#endif

#define FREE_ARG char*
#define MAX_LINE_LENGTH 1024	   /* maximum length of input lines */

extern clock_t last_timer_reset;   /* stores the last timer value */
                                   /* definition in helper.c      */

/* ----------------------------------------------------------------------- */
/*                         auxiliary functions                             */
/* ----------------------------------------------------------------------- */

int min( int a, int b);	       
int max( int a, int b);
double fmin( double a, double b);
double fmax( double a, double b);


/* ----------------------------------------------------------------------- */
/*                         error handling                                  */
/* ----------------------------------------------------------------------- */
/* ERROR(s) writes an error message and terminates the program             */
/*                                                                         */
/* Example:                                                                */
/* ERROR("File not found !");                                              */

#define ERROR(s)    errhandler( __LINE__, __FILE__, s)
#define ERROUT stdout
void  errhandler( int nLine, const char *szFile, const char *szString );


/* ----------------------------------------------------------------------- */
/*                   reading from a datafile                               */
/* ----------------------------------------------------------------------- */
/* The foloowing three macros help reading values from the parameter file. */
/* If a variable cannot be found, the program stops with an error message. */
/*                                                                         */
/* Example:                                                                */
/* READ_INT( "MyFile.dat", imax );                                         */
/* READ_STRING( szFile, szProblem );                                       */

#define READ_INT( szFileName, VarName)    read_int   ( szFileName, #VarName, &(VarName) ) 
#define READ_DOUBLE( szFileName, VarName) read_double( szFileName, #VarName, &(VarName) )
#define READ_STRING( szFileName, VarName) read_string( szFileName, #VarName,  (VarName) )

void read_string( const char* szFilename, const char* szName, char*  sValue);
void read_int   ( const char* szFilename, const char* szName, int*    nValue);
void read_double( const char* szFilename, const char* szName, double*  Value);


/* ----------------------------------------------------------------------- */
/*                   writing matrices to a file                            */
/* ----------------------------------------------------------------------- */
/* write_matrix(...) wites a matrice to a file                               */
/* the file has the following format                                       */
/*                                                                         */
/*    -----------------------------------------                            */
/*    |  xlength          |  float  |  ASCII  |                            */
/*    ----------------------------------------|                            */
/*    |  ylength          |  float  |  ASCII  |                            */
/*    ----------------------------------------|                            */
/*    |  nrl              |  int    |  ASCII  |                            */
/*    ----------------------------------------|                            */
/*    |  nrh              |  int    |  ASCII  |    1. call of the          */
/*    ----------------------------------------|                            */
/*    |  ncl              |  int    |  ASCII  |                            */
/*    ----------------------------------------|                            */
/*    |  nch              |  int    |  ASCII  |    1. call of the          */
/*    ----------------------------------------|    function with            */
/*    |  m[nrl][ncl]      |  float  |  binaer |    bFirst == 1             */
/*    ----------------------------------------|                            */
/*    |  m[nrl][ncl+1]    |  float  |  binaer |                            */
/*    ----------------------------------------|                            */
/*    |                  .                    |                            */
/*                       .                                                 */
/*    |                  .                    |                            */
/*    -----------------------------------------                            */
/*    |  m[nrh][nch]      |  float  |  binary |                            */
/*    -----------------------------------------------------------------    */
/*    |  m[nrl][ncl]      |  float  |  binary |                            */
/*    ----------------------------------------|                            */
/*    |  m[nrl][ncl+1]    |  float  |  binary |     2. call with           */
/*    ----------------------------------------|     bFirst == 0            */
/*    |                  .                    |                            */
/*                       .                                                 */
/*    |                  .                    |                            */
/*    -----------------------------------------                            */
/*    |  m[nrh][nch]      |  float  |  binary |                            */
/*    ------------------------------------------------------------------   */
/*                                                                         */
/*                                                                         */

void write_matrix( const char* szFileName,       /* name of the file */
		   double **m,		       /* matrix */
		   int nrl,		       /* first column */
		   int nrh,		       /* last column */
		   int ncl,		       /* first row */
		   int nch,		       /* last row */
                   double xlength,	       /* size of the geometry in 
                                                  x-direction              */
                   double ylength,	       /* size of the geometry in
                                                  y-direction              */
		   int fFirst, 	       /* 0 == append, else overwrite*/
		   int bPart );                 /* 0 for files corresponding to subdomains, 
                                                  in this case xlength, ylenght, nrl, nrh, ncl, nch are written
                                                  at the beginning of the file,   
                                                  1 for whole region, no additional range information is written */         

void print_matrix(
		   double **m,		       /* matrix */
		   int nrl,		       /* first column */
		   int nrh,		       /* last column */
		   int ncl,		       /* first row */
		   int nch		       /* last row */
	 );

size_t read_matrix( FILE *fh,                    /* file handle */
		    double **m,		         /* matrix */
		    int nrl,		         /* first row */
		    int nrh,		         /* last row */
		    int ncl,		         /* first line */
		    int nch );


/* ----------------------------------------------------------------------- */
/*                                matrices                                 */
/* ----------------------------------------------------------------------- */
/* matrix(...)        storage allocation for a matrix (nrl..nrh, ncl..nch)   */
/* free_matrix(...)   storage deallocation                                   */
/* init_matrix(...)   initialization of all matrix entries with a fixed      */
/*                  (floating point) value                                 */
/* imatrix(...)       analog for matrices with integer-entries               */
/*                                                                         */
/*                                                                         */
/* Example:                                                                */
/*    U     = matrix ( 0 , imax+1 , 0 , jmax+1 );                            */
/*    init_matrix( U , 0, imax+1, 0, jmax+1, 0 );                            */
/*    free_matrix( U,  0, imax+1, 0, jmax+1 );                               */
 
double **matrix( int nrl, int nrh, int ncl, int nch );
void free_matrix( double **m, int nrl, int nrh, int ncl, int nch );
void init_matrix( double **m, int nrl, int nrh, int ncl, int nch, double a);

int  **imatrix( int nrl, int nrh, int ncl, int nch );
void free_imatrix( int **m, int nrl, int nrh, int ncl, int nch );
void init_imatrix( int **m, int nrl, int nrh, int ncl, int nch, int a);


/* ----------------------------------------------------------------------- */
/*                         useful macros                                   */
/* ----------------------------------------------------------------------- */
/*  The following macros can be helpful to display variables during the    */
/*  runtime of the program.                                                */
/*  If you start the program in a shell from xemacs, you can jump to the   */
/*  respectove rows by switching to the compilation-minor-mode.            */
/*                                                                         */
/*  DUMP_POSITION()           dumps the actual position within the program */
/*  DUMP_MESSAGE( .)          dump a message in addition                   */
/*  DUMP_INT(..)              dump an integer variable                     */
/*                                                                         */
/*  DUMP_MATRIX_TO_FILE(..)                                                */
/*  DUMP_INT_TO_FILE(..)      writes the value of the variable in          */
/*  DUMP_DOUBLE_TO_FILE(..)   a tracefile                                  */
/*  DUMP_STRING_TO_FILE(..)                                                */
/*                                                                         */
/*  RESET_TIMER()     set timer to zero                                    */
/*  DUMP_TIMER()      dump time that has passed since the last             */
/*                    RESET_TIMER()                                        */
/*                                                                         */

#define DUMPOUT stdout

#define DUMP_POSITION() fprintf( DUMPOUT, "%s:%d Dumpposition \n", __FILE__, __LINE__ )
#define DUMP_MESSAGE(s) fprintf( DUMPOUT, "%s:%d %s\n",            __FILE__, __LINE__, s  )

#define DUMP_INT(n)     fprintf( DUMPOUT, "%s:%d %s = %d\n", __FILE__, __LINE__, #n, n )
#define DUMP_DOUBLE(d)  fprintf( DUMPOUT, "%s:%d %s = %f\n", __FILE__, __LINE__, #d, d )
#define DUMP_STRING(s)  fprintf( DUMPOUT, "%s:%d %s = %s\n", __FILE__, __LINE__, #s, s )

#define RESET_TIMER()   last_timer_reset = clock()
#define DUMP_TIMER()    fprintf( DUMPOUT, "%s:%d Timer: %f\n", __FILE__, __LINE__, (float)(clock()-last_timer_reset)/(float)CLOCKS_PER_SEC )

#define DUMP_MATRIX_TO_FILE( m, nrl, nrh, ncl, nch, xlength, ylength) \
        {  \
           static nCount = 0; \
	   char szFileName[100];  \
	   sprintf( szFileName, "%s__%d__%s.out", __FILE__, __LINE__, #m); \
           write_matrix( szFileName, m, nrl, nrh, ncl, nch, xlength, ylength, nCount == 0); \
	   ++nCount; \
        }

#define DUMP_INT_TO_FILE(n) \
        {  \
           static nCount = 0; \
           FILE *fh = 0; \
	   char szFileName[100];  \
	   sprintf( szFileName, "%s__%d__%s.out", __FILE__, __LINE__, #n); \
	   if( nCount == 0) \
              fh = fopen( szFileName, "w"); \
           else  \
              fh = fopen( szFileName, "a"); \
           if( fh )  \
              fprintf( fh, "%d:%d\n", nCount, n ); \
           else  \
              ERROR("Fehler beim Dumpen");  \
           fclose(fh);  \
	   ++nCount; \
        }

#define DUMP_DOUBLE_TO_FILE(d) \
        {  \
           static nCount = 0; \
           FILE *fh = 0; \
	   char szFileName[100];  \
	   sprintf( szFileName, "%s__%d__%s.out", __FILE__, __LINE__, #d); \
	   if( nCount == 0) \
              fh = fopen( szFileName, "w"); \
           else  \
              fh = fopen( szFileName, "a"); \
           if( fh )  \
              fprintf( fh, "%d:%f\n", nCount, d ); \
           else  \
              ERROR("Fehler beim Dumpen");  \
           fclose(fh);  \
	   ++nCount; \
        }

#define DUMP_STRING_TO_FILE(s) \
        {  \
           static nCount = 0; \
           FILE *fh = 0; \
	   char szFileName[100];  \
	   sprintf( szFileName, "%s__%d__%s.out", __FILE__, __LINE__, #s); \
	   if( nCount == 0) \
              fh = fopen( szFileName, "w"); \
           else  \
              fh = fopen( szFileName, "a"); \
           if( fh )  \
              fprintf( fh, "%d:%s\n", nCount, s ); \
           else  \
              ERROR("Fehler beim Dumpen");  \
           fclose(fh);  \
	   ++nCount; \
        }

#endif     

