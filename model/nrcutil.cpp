/*
 * Programmed by Hsiao-Te Su
 *
 * Numerical Recipe in C utilities
 */
//#include "stdafx.h"
//#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrcutil.h"

void nrerror( const char *error_text )
{
  fprintf( stderr, "Numerical Recipes run-time error...\n" );
  fprintf( stderr, "%s\n", error_text );  
  fprintf( stderr, "...now exiting to system...\n" );
  exit( 1 );
}

float *vector( int nl, int nh )
{
  float *v;

  v = (float *) malloc( (unsigned) (nh - nl + 1)*sizeof(float) );
  if( !v )
    nrerror( "allocation failure in vector()" );

  return( v - nl );
}

void free_vector( float *v, int nl, int nh )
{
  free( (char *) (v + nl) );
}

int *ivector( int nl, int nh )
{
  int *v;

  v = (int *) malloc( (unsigned) (nh - nl + 1)*sizeof(int) );
  if( !v )
    nrerror( "allocation failure in vector()" );

  return( v - nl );
}

void free_ivector( int *v, int nl, int nh )
{
  free( (char *) (v + nl) );
}

float **matrix( int nrl, int nrh, int ncl, int nch )
{
  int i;
  float **m;

  m = (float **) malloc((unsigned) (nrh - nrl + 1)*sizeof(float*));
  if( !m ) nrerror( "allocation failure 1 in matrix()" );
  m -= nrl;

  for( i = nrl; i <= nrh; i++ )
  {
    m[i] = (float *) malloc((unsigned) (nch - ncl + 1)*sizeof(float));
    if( !m[i] ) nrerror( "allocation failure 2 in matrix()" );
    m[i] -= ncl;
  }
  return( m );
}

void free_matrix( float **m, int nrl, int nrh, int ncl, int nch )
{
  int i;

  for( i = nrh; i >= nrl; i-- ) free( (char *) (m[i] + ncl) );
  free( (char*) (m + nrl) );
}


/*
 * Allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
 * a declared in the standard C manner as a[nrow][ncol], where 
 * nrow = nrh - nrl + 1 and ncol = nch - ncl + 1.  The routine should
 * be called with the address &a[0][0] as the first argument.
 */
float **convert_matrix( float *a, int nrl, int nrh, int ncl, int nch )
{
  int i, j, nrow, ncol;
  float **m;

  nrow = nrh - nrl + 1;
  ncol = nch - ncl + 1;

  /* allocate pointers to rows */
  m = (float **) malloc( (unsigned) (nrow)*sizeof(float*) );
  if( !m )
    nrerror( "allocation failure in convert_matrix()" );
  m -= nrl;
  for( i = 0, j = nrl; i <= nrow - 1; i++, j++ )
    m[j] = a + ncol*i - ncl;    /* set pointers to rows */

  return( m );
}

void free_convert_matrix( float **b, int nrl, int nrh, int ncl, int nch )
{
  free( (char *) (b + nrl) );
}
