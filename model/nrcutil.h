/*
 * Programmed by Hsiao-Te Su
 *
 * Numerical Recipe in C utilities
 */
#ifndef NRCUTIL_H
#define NRCUTIL_H

void nrerror( const char *error_text );

float *vector( int nl, int nh );

void free_vector( float *v, int nl, int nh );

int *ivector( int nl, int nh );

void free_ivector( int *v, int nl, int nh );

float **matrix( int nrl, int nrh, int ncl, int nch );

void free_matrix( float **m, int nrl, int nrh, int ncl, int nch );

/*
 * Allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
 * a declared in the standard C manner as a[nrow][ncol], where 
 * nrow = nrh - nrl + 1 and ncol = nch - ncl + 1.  The routine should
 * be called with the address &a[0][0] as the first argument.
 */
float **convert_matrix( float *a, int nrl, int nrh, int ncl, int nch );

void free_convert_matrix( float **b, int nrl, int nrh, int ncl, int nch );
#endif // NRCUTIL_H
