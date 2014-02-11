/*
 * Programmed by Hsiao-Te Su
 * Mon 09-12-1994
 *
 * A generic matrix class using floats only 
 */
//#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <memory.h>
#include <iostream>
#include <math.h>
#include "nrcutil.h"
#include "matrix.h"

#include "iia_err.h"

FMatrix::FMatrix()                        // default constructor
{
  xdim = 0;
  ydim = 0;
  elem = NULL;
}

FMatrix::FMatrix( int y, int x )  
{
  assert( x >= 0 && y >= 0 );

  xdim = x;
  ydim = y;

  if( x > 0 && y > 0 )
  {
		elem = new float[ xdim*ydim ];
		if( elem == NULL )
			crash( "Out of memory in FMatrix( y, x )" );
	}
	else
		elem = NULL;
}

FMatrix::FMatrix( const FMatrix &M )
{
  xdim = M.XDim();
  ydim = M.YDim();

  if( xdim > 0 && ydim > 0 )
  {
    elem = new float[ xdim*ydim ];
    if( elem == NULL )
    	crash( "Out of memory in FMatrix's copy constructor" );

    for( int iy = 0; iy < ydim; iy++ )
      for( int ix = 0; ix < xdim; ix++ )
        elem[ iy*xdim + ix ] = M.rElem( iy, ix );
  }
  else
    elem = NULL;
}

/*
 * Reset the dimension of the matrix, previous information is lost
 * Returns 0 on memory allocation failure, 1 on success
 */
int FMatrix::SetDim( int y, int x )
{
	int ret;

  if( elem != NULL )
    delete []elem;

  xdim = x;
  ydim = y;

	if( xdim == 0 || ydim == 0 )
	{
		elem = NULL;
		ret = 1;
	}
	else
	{
		elem = new float[ xdim * ydim ];
		if( elem == NULL )
			ret = 0;
		else
			ret = 1;
  }

  return ret;
}

/* 
 * Print out the matrix 
 */
void FMatrix::Print() const
{
#if 0
  for( int iy = 0; iy < ydim; iy++ )
  {
    for( int ix = 0; ix < xdim; ix++ )
      cout << elem[ iy*xdim + ix ] << " ";
    cout << endl;
  }
#else
/*  TRACE( "\n" );
  for( int iy = 0; iy < ydim; iy++ )
  {
    for( int ix = 0; ix < xdim; ix++ )
      TRACE( "%g ", elem[ iy*xdim + ix ] );
    TRACE( "\n" );
  }
  TRACE( "\n" );
*/
#endif
}


/*
 * Assignment operator for integers
 * 0 - zeros the matrix
 * 1 - sets it to identy
 * anything else is an error
 */
int FMatrix::operator=( int x )              // sets matrix to 0 or 1
{
  if( x == 0 )
    zeroAll();
  else if( x == 1 )
  {
    assert( xdim == ydim );
    for( int ix = 0; ix < xdim; ix++ )
      for( int iy = 0; iy < ydim; iy++ )
        elem[ iy*xdim + ix ] = (float) (ix == iy); 
  }
  else
    assert( 0 );

  return( x );
}

/*
 * copy matrix
 */
FMatrix& FMatrix::operator=( const FMatrix &M )
{
  if( this == &M )
    return( *this ); 

  if( M.XDim() > 0 && M.YDim() > 0 )
  {
    if( SetDim( M.XDim(), M.YDim() ) == 0 )
    	crash( "Out of memory in FMatrix assignment operator" );

    for( int iy = 0; iy < ydim; iy++ )
      for( int ix = 0; ix < xdim; ix++ )
        elem[ iy*xdim + ix ] = M.rElem( iy, ix );
  }
  else
  {
    xdim = M.XDim();
    ydim = M.YDim();
    elem = NULL;
  }

  return( *this );
}


/*
 * Increment 
 */
FMatrix& FMatrix::operator+=( const FMatrix &M )
{
  // check dimension
  assert( XDim() == M.XDim() && YDim() == M.YDim() );

  for( int iy = 0; iy < ydim; iy++ )
    for( int ix = 0; ix < xdim; ix++ )
      elem[ iy*xdim + ix ] += M.rElem( iy, ix );

  return( *this );
}



/*
 * Multiply the matrix
 * nxm matrix multiply by mxr matrix gives nxr matrix
 */
FMatrix FMatrix::operator*( const FMatrix &M ) const
{
  // check dimension
  assert( XDim() == M.YDim() );

  // assign convience variables
  int n = YDim();
  int m = XDim();
  int r = M.XDim();

  FMatrix ret( n, r );

  ret = 0;

  for( int in = 0; in < n; in++ )
    for( int ir = 0; ir < r; ir++ )
      for( int im = 0; im < m; im++ )
        ret.Elem( in, ir ) += rElem( in, im ) * M.rElem( im, ir );

  return( ret );
}

/*
 * Multiply the M matrix by a scalar x
 */
FMatrix operator*( float x, const FMatrix &M ) 
{
  // assign convience variables
  int n = M.YDim();
  int m = M.XDim();

  FMatrix ret( n, m );

  for( int in = 0; in < n; in++ )
    for( int im = 0; im < m; im++ )
      ret.Elem( in, im ) = x * M.rElem( in, im );

  return( ret );
}


/*
 * Multiply the M matrix by a scalar x
 */
FMatrix operator*( const FMatrix &M, float x ) 
{
  return( x*M );
}


/*
 * Invert the matrix if x = -1 
 * use gaussj
 */
FMatrix FMatrix::operator^( int x ) const
{
  assert( x == -1 );      /* may want to extend to 2... in the future */

  // check dimension
  assert( XDim() == YDim() );

  // convience variables
  int n = XDim();
  int size = n*n;

  // make a copy of the matrix for nrc routine to work on
  float *copy = new float[ size ];
  memcpy( copy, elem, size*sizeof(float) );

  float **a = convert_matrix( copy, 1, n, 1, n );
  float **b = matrix( 1, n, 1, 1 );
  b[1][1] = 1.0F;
  for( int i = 2; i <= n; i++ )
    b[i][1] = 0.0F;

  gaussj( a, n, b, 1 );

  FMatrix ret( n, n );
  for( int iy = 0; iy < n; iy++ )
    for( int ix = 0; ix < n; ix++ )
      ret.Elem( iy, ix ) = copy[ iy*n + ix ];

  delete [] copy;
  free_matrix( b, 1, n, 1, 1 );
  free_convert_matrix( a, 1, n, 1, n );

  return( ret );
}


/*
 * Transpose
 */
FMatrix FMatrix::T()
{
  FMatrix ret( xdim, ydim );

  for( int iy = 0; iy < ydim; iy++ )
    for( int ix = 0; ix < xdim; ix++ )
      ret.Elem( ix, iy ) = rElem( iy, ix );

  return( ret );
}



void FMatrix::zeroAll()
{
  int size = xdim*ydim;

  for( int i = 0; i < size; i++ )
    elem[i] = 0.0F;
}




FMatrix::~FMatrix()                        // default destructor
{
  delete [] elem;
  elem = 0;
}


/*
 * The following is taken from Numerical Recipe in C
 */


#define TINY 0; 

/*
 * Given an nxn matrix a[1..n][1..n], this routine replaces it by the LU 
 * decomposition of a rowwise permutation of itself.  a and n are input.
 * a is output; indx[1..n] is an output vector which records the row
 * permutation effected by the partial pivoting; d is output as +-1
 * depending on whether the number of row interchanges was even
 * or odd, respectively.  This routine is used in combination with
 * lubksb to solve linear equations or invert a matrix.
 */
static void ludcmp( float **a, int n, int *indx, float *d )
{
  int i, imax, j, k;
  float big, dum, sum, temp;
  float *vv;                  /* vv stroes the implicit scaling for each row */

  vv = vector( 1, n );
  *d = 1.0;                    /* No row interchanges yet */

  for( i = 1; i <= n; i++ )    /* Loop over rows to get the implicit scaling information */
  {
    big = 0.0;
    for( j = 1; j <= n; j++ )
      if( (temp = (float) fabs( a[i][j] )) > big ) big = temp;
    if( big == 0.0 ) nrerror( "Singular matrix in routine LUDCMP" );
    /* No nonzero largest element */
    vv[i] = 1.0F/big;          /* save the scaling factor */
  }

  for( j = 1; j <= n; j++ )
  {                            /* This is the loop over columns of Crout's method */
    for( i = 1; i < j; i++ )
    {
      sum = a[i][j];
      for( k = 1; k < i; k++ ) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;                /* Initialize for the search for largest pivot element */
    for( i = j; i <= n; i++ )
    {                        
      sum = a[i][j];
      for( k = 1; k < j; k++ )
        sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      if( (dum = vv[i]*((float) fabs( sum ))) >= big )
      {                       /* Is the figure of merit for the pivot better than the best so far */  
        big = dum;
        imax = i;
      }
    }
    if( j != imax )
    {                        /* Do we need to interchange rows? */
      for( k = 1; k <= n; k++ )
      {                      /* Yes we do */
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      *d = -(*d);            /* and change the parity of d */
      vv[imax] = vv[j];      /* Also interchange the scale factor */
    }
    indx[j] = imax;
    if( a[j][j] == 0.0 ) a[j][j] = TINY;
    /* 
     * If the pviot element is zero the matrix is singular (at least to the
     * precision of the algorithm).  For some applications on singular matrices
     * it is desirable to substitute TINY for zero.
     */
    if( j != n )
    {                        /* Now, finally, divide by the pivot element */
      dum = 1.0F/(a[j][j]);
      for( i = j + 1; i <= n; i++ )
        a[i][j] *= dum;
    }
  }
  free_vector( vv, 1, n );
}

/*
 * Solves the set of n linear equations AX = B.  Here a[1..n][1..n] is input,
 * not as the matrix A, but rather its LU decomposition, determined by the
 * routine ludcmp.  indx[1..n] is input as the permutation vector returned
 * by ludcmp. b[1..n] is input as the right-hand side vector B, and returns
 * with the solution vector X. a, n, and indx are not modified by this
 * routine and can be left in place for successive calls with different
 * right-hand sides b.  This routine takes into account the possibility that 
 * b will begin with many zero elements, so it is efficient for use in 
 * matrix inversion
 */
static void lubksb( float **a, int n, int *indx, float b[] )
{
  int i, ii = 0, ip, j;
  float sum;

  for( i = 1; i <= n; i++ )
  {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if( ii )
      for( j = ii; j <= i - 1; j++ )
        sum -= a[i][j]*b[j];
    else if( sum ) 
      ii = i;
    b[i] = sum;
  }

  for( i = n; i >= 1; i-- )
  {
    sum = b[i];
    for( j = i + 1; j <= n; j++ )
      sum -= a[i][j]*b[j];
    b[i] = sum/a[i][i];
  }
}

#define SWAP(a,b) {float temp = (a); (a) = (b); (b) = temp;}

/*
 * Linear equation solution by Gauss-Jordan elimination, a[1..n][1..n]
 * is an input matrix of n by n elements.  b[1..n][1..m] is an input
 * matrix of size n by m containing the m right-hand side vectors.  On
 * output, a is replaced by its matrix inverse, and b is replaced by
 * the corresponding set of solution vectors.
 */
void gaussj(float **a, int n, float **b, int m)
{
  int *indxc,*indxr,*ipiv;
  int i,icol,irow,j,k,l,ll;
  float big,dum,pivinv;
  
  indxc=ivector(1,n);
  indxr=ivector(1,n);
  ipiv=ivector(1,n);
  for(j=1;j<=n;j++) ipiv[j]=0;
  for(i=1;i<=n;i++) {
    big=0.0;
    for(j=1;j<=n;j++)
      if(ipiv[j] != 1)
        for(k=1;k<=n;k++) {
          if(ipiv[k] == 0) {
            if(fabs(a[j][k]) >= big) {
              big=(float)fabs(a[j][k]);
              irow=j;
              icol=k;
            }
          } else if(ipiv[k] > 1) {
            fprintf(stderr,"GAUSSJ: Singular matrix-1");
            free_ivector(ipiv,1,n);
            free_ivector(indxr,1,n);
            free_ivector(indxc,1,n);
            return;
          }
        }
    ++(ipiv[icol]);
    if(irow != icol) {
      for(l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l]);
      for(l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l]);
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if(a[icol][icol] == 0.0) {
      fprintf(stderr, "GAUSSJ: singular matrix-2");
      free_ivector(ipiv,1,n);
      free_ivector(indxr,1,n);
      free_ivector(indxc,1,n);
      return;
    }
    pivinv=1.0F/a[icol][icol];
    a[icol][icol]=1.0;
    for(l=1;l<=n;l++) a[icol][l] *= pivinv;
    for(l=1;l<=m;l++) b[icol][l] *= pivinv;
    for(ll=1;ll<=n;ll++)
      if(ll != icol) {
        dum=a[ll][icol];
        a[ll][icol]=0.0;
        for(l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
        for(l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
      }
  }
  for(l=n;l>=1;l--) {
    if(indxr[l] != indxc[l])
      for(k=1;k<=n;k++)
        SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  }
  free_ivector(ipiv,1,n);
  free_ivector(indxr,1,n);
  free_ivector(indxc,1,n);
}
