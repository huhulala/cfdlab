#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int imax,
  int jmax,
  int rank_l,
  int rank_r,
  int rank_b,
  int rank_t,
  double **U,
  double **V
);

#endif
