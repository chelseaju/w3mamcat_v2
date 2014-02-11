/*
 * Programmed by Hsiao-Te Su
 * Mon 09-12-1994
 *
 * A generic matrix class using floats only 
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <assert.h>

class FMatrix
{
public:
  FMatrix();                        // default constructor
  FMatrix( int ydim, int xdim );  
  FMatrix( const FMatrix &M );      // memberwise init constructor

  ~FMatrix();                        // default destructor

public:
  int SetDim( int ydim, int xdim );    // reset the dimension of the matrix
  int XDim() const { return xdim; };  // read only access to xdim
  int YDim() const { return ydim; };  // read only access to ydim
  float& Elem( int iy, int ix )       // direct access to the (ix,iy)th element
  { 
    assert( iy >= 0 && iy < ydim && ix >= 0 && ix < xdim );   // so I won't shoot myself in the foot
    return elem[iy*xdim+ix];         
  };
  float  rElem( int iy, int ix ) const // read only access to the (ix,iy)th element
  { 
    assert( iy >= 0 && iy < ydim && ix >= 0 && ix < xdim );
    return elem[iy*xdim+ix]; 
  }; 
    // same as above, except offset by 1
  float& Elem1( int iy, int ix ) { return Elem( iy-1, ix-1 ); };        // direct access to the (ix,iy)th element
  float  rElem1( int iy, int ix ) const { return rElem( iy-1, ix-1 ); }; // read only access to the (ix,iy)th element
  void Print() const;                    // print out the matrix

  // operators
  int operator=( int x );                      // sets matrix to 0 or 1
  FMatrix& operator=( const FMatrix &M );      // copy matrix
  FMatrix& operator+=( const FMatrix &M );    // increment
  FMatrix operator*( const FMatrix &M ) const;// multiply matrices
  FMatrix operator^( int x ) const;          // x = -1 for inversion

  FMatrix T();                          // Transpose;  

private:
  void zeroAll();                        // zeros the element array

  // Data
private:
  int xdim;        // horizontal dimension of the matrix
  int ydim;        // vertical dimension of the matrix
  float *elem;    // 2D array of floats, xdim is the smallest index
};

FMatrix operator*( const FMatrix &M, float x );        // scalar multipliation
FMatrix operator*( float x, const FMatrix &M );        // scalar multipliation

static void ludcmp( float **a, int n, int *indx, float *d );
static void lubksb( float **a, int n, int *indx, float b[] );
static void gaussj( float **a, int n, float **b, int m );

#endif MATRIX_H
