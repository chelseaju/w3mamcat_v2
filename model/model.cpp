// Programmed by Hsiao-Te Su
// 12/26/93
 
//#include "stdafx.h"
#include <stdio.h>  
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
//#include <malloc.h>
#include <math.h>
#include <assert.h>
#include "iia_num.h"
#include "model.h"
#include "iia_mem.h" 
#include "iia_err.h"
#include "iia_key.h"
#include "nrcutil.h"
                             
#define LINDELL 1				// Robert Lindell's way of V1

#ifdef DEBUG
#undef DEBUG
#endif

#define DEBUG  0

#define BUFSIZE  256


char* Model::timestr[] = { "sec", "min", "hr", "day" };
char* Model::volstr[] = { "ml", "L", "ul" };
char* Model::massstr[] = { "%Dose", "gm", "kg", "lb" };
char* Model::endomustr[] = { "gm", "kg", "lb" };



// Constructors

Model::Model()      // default constructor
{
  Type = Mam;    // default type
  N = 0;
  C1 = 1;
  outunit = Conc;
  bw = 1;
  dose = 100;			// need to be consistent with the default massunit = %dose 
  for( int i = 0; i < MAXPARAM; i++ )
    constraint[i] = UNCONSTRAINED;
  errmod = CCV; 
  timen = 0;
  ccv = (float) 0.0;
  csd = (float) 0.0;
  cvar = (float) 0.0;
  varb = (float) 0.0;
  varc = (float) 0.0;
  vard = (float) 0.0;
  timeunit = T_SEC;
  volunit = V_ML;
	massunit = M_PDOSE;
	endomu = ENDO_M_G;
}



// Destructors
Model::~Model()
{
}


/*
 * Check for internal consistencies
 */
bool Model::Check()
{
  bool result = true;
  char errMsg[256];
  int i;

  // Check model parameters
  if( PoolNum() <= 0 || PoolNum() > 10 )
  {
    sprintf( errMsg, "The number of pools is out of bounds\n" );
    warn( errMsg );
    result = false;
  }
  else
  {
    for( i = 1; i <= PoolNum(); i++ )
    {
      if( A(i) <= 0 )
      {
        sprintf( errMsg, "A(%d) needs to be greater than 0\n", i );
        warn( errMsg );
        result = false;
      }
      else if( L(i) >= 0 )
      {
        sprintf( errMsg, "L(%d) needs to be less than 0\n", i );
        warn( errMsg );
        result = false;
      }
      if( result == false )
        break;
    }
  }

  // Check sample inputs
  if( result )
  {
    if( SampleN() < 0 )
    {
      warn( "Number of Samples needs to be greater than 0" );
      result = false;
    }
    else
      for( i = 1; i <= SampleN(); i++ )
        if( SampleTime( i ) < 0 )
        {
          sprintf( errMsg, "Sample Time(%d) needs to be greater than 0\n", i );
          warn( errMsg );
          result = false;
          break;
        }
  }

  // Check Error model
  if( result && SampleN() > 0 )
  {
    switch( ErrMod() )
    {
      case CCV:
        if( Ccv() < 0 )
        {
          warn( "CV needs to be greater than or equal to 0" );
          result = false;
        }
        break;
      case CSD:
        if( Csd() < 0 )
        {
          warn( "Standard deviation needs to be greater than or equal to 0" );
          result = false;
        }
        break;
      case CVAR:
        if( Cvar() < 0 )
        {
          warn( "Variance needs to be greater than or equal to 0" );
          result = false;
        }
        break;
      case VVAR:
        break;
      case VCV:
        for( i = 1; i <= SampleN(); i++ )
        {
          if( Vcv( i ) < 0 )
          {
            sprintf( errMsg, "CV(%d) needs to be greater than or equal to 0\n", i );
            warn( errMsg );
            result = false;
            break;
          }
        }
        break;
      default:
        assert( 0 );
        break;
    }
  }

  // Check Dose
  if( MassUnit() == M_PDOSE && Dose() != 100 )
  {
  	result = false;
  	warn( "Dose needs to be 100 since A(i) units are %Dose" );
  }

  return result;
}

/*
 * Calculate the model
 */
int Model::Calculate()
{
  float    *org_a = new float[N+1];
  int			 i;

	// sort A and L
	SortAL();

	// store the original A 
  for( i = 1; i <= N; i++ )
  	org_a[i] = A(i);

	// scale the A and L by BodyWeight and Dose
	// Probably should change any constant factors here in the future
  for( i = 1; i <= N; i++ )
  	a_exp[i-1] *= (BodyWeight()/Dose());
	
	// do the computation
	ALtoAB();
	ABtoKG();  
	UnID();  
	ConID();
	VQID();     
	MfluxID();
	Global();
	GetCorr();

	// restore the A and L
  for( i = 1; i <= N; i++ )
  	A(i) = org_a[i]; 

	delete [] org_a;
  return( 0 );
}



// helping functions
             
// calculations 

/*
 * Evaluates the sum of exponential function at x via the formula
 * Sum( Ai*e^Li )
 */
float  Model::SumOfExp( float x )
{
  float  val = 0.0;
  int i;

  for( i = 1; i <= N; i++ )
    val += A(i)*(float) exp( (double) x*L(i) );
    
  return( val );    
}

/*
 * Evaluates the transfer function at x via the formula
 * Sum( Ai/s-Li )
 */
float  Model::Transfer( float x )
{
    int i;
  for(i = 1; i <= N; i++ )
    if( x == L(i) )
      crash( "attempting to evaluate the transfer function at lambda[i]" );
    
  float  val = 0.0;
  for( i = 1; i <= N; i++ )
    val += A(i)/(x - L(i));
    
  return( val );    
}

                   
// Finds the root of the transfer function between x1 and x2 
// with in tolerance level( +- tol )
// This uses the Brent Method as presented in Numerical Recipes
// in C on page 268.

#define ITMAX 100
#define EPS  3.0e-8F

float  Model::FindRoot( float x1, float x2, float tol )
{                         
  int   iter;
  float  a = x1, b = x2, c, d, e, min1, min2;
  float  fa = Transfer(a), fb = Transfer(b), fc, p, q, r, s, tol1, xm;
                                                           
  if( fb*fa > 0.0) 
    crash("Root must be bracketed in ZBRENT");
  fc = fb;
  for( iter = 1; iter <= ITMAX; iter++ )
  {
    if( fb*fc > 0.0 )
    {
      c = a;          // Rename a, b, c and adjust bounding interval d
      fc = fa;
      e = d = b-a;
    }   
    if( fabs(fc) < fabs(fb) )
    {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }           
    tol1 = 2.0F*EPS*(float) fabs(b) + 0.5F*tol;  // Convergence check
    xm = 0.5F*(c-b);
    if( fabs(xm) <= tol1 || fb == 0.0 )
      return b;
    if( fabs(e) >= tol1 && fabs(fa) > fabs(fb) )
    {
      s = fb/fa;      // Attempt inverse quadratic interpolation
      if( a == c )
      {
        p = 2.0F*xm*s;
        q = 1.0F-s;
      }             
      else
      {
        q = fa/fc;
        r = fb/fc;
        p = s*(2.0F*xm*q*(q-r)-(b-a)*(r-1.0F));
        q = (q-1.0F)*(r-1.0F)*(s-1.0F);
      }                               
      if( p > 0.0 ) q = -q;  // Check whether in bounds
      p = (float) fabs(p);
      min1 = 3.0F*xm*q - (float) fabs(tol1*q);
      min2 = (float) fabs(e*q);
      if( 2.0*p < (min1 < min2 ? min1 : min2) )
      {
        e = d;        // Accept interpolation
        d = p/q;
      }           
      else
      {
        d = xm;        // interpolation failed, use bisection
        e = d;
      }         
    }
    else                // Bounds decreasing too slowly, use bisection
    {
      d = xm;
      e = d;
    }
    
    a = b;          // Move last best guess to a
    fa = fb;
    if( fabs(d) > tol1 )  // Evaluate new trial root
      b += d;
    else
      b += (xm > 0.0 ? (float) fabs(tol1) : -(float) fabs(tol1) );
    fb = Transfer(b);  
  }                                                           
  crash("Maximum number of iterations exceeded in ZBRENT");
  return( 0 );        // to satisfy the compiler
}                   

                   
// This function converts A and L to Alpha and Beta
// A's are the coefficients in the sum of exponential equation 
// L's are the exponents in the sum of exponential equation
// Alpha's are the coeficients of the demonator of the transfer ftn
// Beta's are the coeficients of the numerator of the transfer ftn
//
// The detail algorithm is documented in: 
// Landaw, Chen & DiStefano III
// "An Algorithm for the Identifiable Parameter Combinations of the
// General Mammillary Compartmental Model" Mathematical Biosciences
// 72:199-212(1984) 
// 
// It is a simple recursive relation, just need to be a little careful
// with boundary conditions
int Model::ALtoAB()
{      
  int    i,k;

  // initial condition, i = 1 and k = 1
  Beta(1) = A(1);  
  Alpha(1) = -L(1);
  
  for( k = 2; k <= N; k++ )
  {
    // i = k boundary condition
    Beta(k) = Beta(k-1) + A(k);    // recursive formula
    Alpha(k) = Alpha(k-1) - L(k);

    // normal recursive formula
    for( i = k - 1; i >= 2; i-- )
    {             
      Beta(i) = Beta(i-1) - L(k)*Beta(i) + A(k)*Alpha(i);  // recursive formula
      Alpha(i) = Alpha(i-1) - L(k)*Alpha(i);
    }  

    // i = 1 boundary condition
    Beta(1) = - L(k)*Beta(1) + A(k)*Alpha(1);  // recursive formula
    Alpha(1) = - L(k)*Alpha(1);
  }

#if 0
	// This is original mampool scaling

  float scale;

	scale = 1.0/100;			

  for( i = 1; i <= N; i++ )
    Beta(i) *= (bw*scale);           
#endif
    
  return( 0 );
}

// converts alpha and beta to kappa and gamma
// kappa and gamma are the identifible parameter combinations
// The algorithm is documented in the same article as above
//                                          
// For mamillary models:
// Basically, kappas are the roots of the transfer ftn and
// gammas are the weighted harmonic mean of the (kappa-L)^2    
//
// For catenary models:
// Kii and gamma are evaluated via the coefficients of the eigenvalue determinant
// See Algorithms for the Identifiable Parameter Combinations and
// Parameter Bounds of Unidentifiable Catenary Compartmental Models
int Model::ABtoKG()
{   
  if( Type == Mam )            // Mamillary Algorithm
  {                         
    float  delta;           
    float  lower, upper;
    //+======================================================
    //|
    //| MODIFIEDCODE
    //|
    //| KB 5/10/2001
    //|
    //| Force float representation with an "F".
    //|
    //+======================================================
    //float  tol = 5.0e-8;
    float  tol = 5.0e-8F;
    
      // computing kappa
    if( N > 1 )
      K(1,1) = -Alpha(N) + Beta(N-1)/Beta(N);
    else
      K(1,1) = -Alpha(N);
      int i; // SLR
    for( int i=2; i <= N; i++ )
    {                        
      // find 2 pts that bound the root
      // assumes that L's are sorted
      delta = (L(i) - L(i-1))/10.0F;
      lower = L(i-1) + delta;
      upper = L(i) - delta;
      while( Transfer( lower )*Transfer( upper ) > 0.0 )
      {
        delta *= 0.5F;    
        lower = L(i-1) + delta;
        upper = L(i) - delta;
      }
    
      K(i,i) = FindRoot( lower, upper, tol );
    }  

    // computing gamma
    float  sumN = 0.0F;      // numerator sum
    float  sumD = 0.0F;      // denominator sum
  
    for( i = 1; i <= N; i++ )
      sumN += A(i);
    
    for( int j = 2; j <= N; j++ )
    {       
      sumD = 0.0F;
      for( i = 1; i <= N; i++ )
        sumD += A(i)/(( K(j,j)-L(i) )*( K(j,j)-L(i) ));
      
      Gamma(j) = sumN/sumD;
    }             
  
    gamma[0] = 0.0F;        // not used
  }
  else                // Catenary algorithm
  {               
    int   n = PoolNum();
    float  **p = matrix( 1, n, 1, n );
    int    i, j, m;

    for( m = 1; m <= n; m++ )        // boundary condition
      for( j = n-m+2; j <= n; j++ )
        p[m][j] = 0.0F;

    for( j = 1; j <= n; j++ )        // initialization
      p[1][j] = Alpha(n-j+1);
    for( j = 1; j <= n-1; j++ )
      p[2][j] = Beta(n-j)/Beta(n);
      
    for( i = 1; i <= n-2; i++ )
    {  
      // evaluate the Ks
      K(i,i) = p[i+1][1] - p[i][1];
                         
      // evaluate the Gammas
      Gamma(i+1) = p[i+1][2] - p[i][2] - K(i,i)*p[i+1][1];
      
      // advance the p's
      for( j = 1; j <= n-i-1; j++ )
        p[i+2][j] = ( p[i+1][j+2] - p[i][j+2] - K(i,i)*p[i+1][j+1] )/Gamma(i+1);
    }       

    K(N-1,N-1) = p[N][1] - p[N-1][1];
    Gamma(N) = -p[N-1][2] - K(N-1,N-1)*p[N][1];

    K(N,N) = -p[N][1];
      
    free_matrix( p, 1, n, 1, n );
  }

  return( 0 );  
}   
 
// Compute the max and min bounds of the unidentifible parameters
// This version does not consider constraints
// The Mammillary algorithm is new, derived as in catpool.
int  Model::UnID()
{       
  int i,j;
    
  if( Type == Mam )          // Mamillary algorithm
  {                                                        
    for( i = 2; i <= N; i++ )  
    {
      K(1,i,Upper) = -K(i,i);
      K(i,1,Lower) = Gamma(i)/K(1,i,Upper);
    }   
    
    for( i = 2; i <= N; i++ )
    {
      K(i,1,Upper) = -K(1,1);
      for( j = 2; j <= N; j++ )
        if( j != i )
          K(i,1,Upper) -= K(j,1,Lower);
      K(1,i,Lower) = Gamma(i)/K(i,1,Upper);
    }   
    
    K(0,1,Upper) = -K(1,1);
    for( i = 2; i <= N; i++ )
      K(0,1,Upper) -= K(i,1,Lower);                                     
      
    for( i = 2; i <= N; i++ )
      K(0,i,Upper) = -K(i,i) - K(1,i,Lower);
  
    if( N > 1 )    
      for( i = 1; i <= N; i++ )
        K(0,i,Lower) = 0.0;      
    else
      K(0,1,Lower) = -K(1,1);
  }  
  else                // Catenary algorithm
  {                              
    if( N > 1 )
    {   
      K(2,1,Upper) = -K(1,1);
      for( i = 2; i <= N-1; i++ )    
        K(i+1,i,Upper) = -K(i,i) - Gamma(i)/K(i,i-1,Upper);
      
      for( i = 1; i <= N-1; i++ )
        K(i,i+1,Lower) = Gamma(i+1)/K(i+1,i,Upper);
      
      K(N-1,N,Upper) = -K(N,N);
      for( i = N-1; i >= 2; i-- )                              
        K(i-1,i,Upper) = -K(i,i) - Gamma(i+1)/K(i,i+1,Upper);
      
      for( i = 2; i <= N; i++ )
        K(i,i-1,Lower) = Gamma(i)/K(i-1,i,Upper);  
    }
    
    if( N > 1 )  
      K(0,1,Upper) = -K(1,1) - K(2,1,Lower);
    else
      K(0,1,Upper) = -K(1,1);
      
    for( i = 2; i <= N-1; i++ )
      K(0,i,Upper) = -K(i,i) - K(i-1,i,Lower) - K(i+1,i,Lower);

    if( N > 1 )
      K(0,N,Upper) = -K(N,N) - K(N-1,N,Lower);

    if( N > 1 )    
      for( i = 1; i <= N; i++ )
        K(0,i,Lower) = 0.0;
    else
      K(0,1,Lower) = -K(1,1);
  }

  return( 0 );
}


// Compute the max & min bounds of the constrained case
// There's no documentation on this right now, but Paolo
// does have a paper pending for catpool
int  Model::ConID()
{                
  int i;

  if( N == 1 )
  {
    // only 1 pool, so identifible  
    ConK(0,1,Lower) = K(0,1,Lower);
    ConK(0,1,Upper) = K(0,1,Upper);
  }
  else
  {
    // check to see if the constraints are within
    // the unconstrained bounds.  no solution would
    // exist if not.
    if( CheckConstr() != 0 )
      crash( "Constraints are out of bounds" );
        
    // copy the unconstrained bounds to the constrained bounds
    for( i = 0; i < 3*N-2; i++ )
    {
      upperC[i] = upperK[i];
      lowerC[i] = lowerK[i];
    }
    
    if( Type == Mam )    // mamillary
    {                        
      // check if there are over constraints         
      // this part would probably be pulled out later
      
      for( i = 2; i <= N; i++ )
        if( (Constr(0,i) != UNCONSTRAINED && Constr(1,i) != UNCONSTRAINED) ||
          (Constr(1,i) != UNCONSTRAINED && Constr(i,1) != UNCONSTRAINED) )
          crash( "This model is over constrained" );
          
      int over = (Constr(0,1) != UNCONSTRAINED);
      for( i = 2; i <= N; i++ )
        over = over && (Constr(i,1) != UNCONSTRAINED);
      if( over )
        crash( "This model is over constrained" );      
              
      // First check for constraints on the peripheral pools
      for( i = 2; i <= N; i++ )
      {
        int count = 0;  // count the number of constrains
        if( Constr(0,i) != UNCONSTRAINED )
          count++;
        if( Constr(1,i) != UNCONSTRAINED )
          count++;
        if( Constr(i,1) != UNCONSTRAINED )
          count++;
          
        switch( count )
        {
          case 0:    // no constrain on pool i                
                      // don't need to adjust the bounds
            break;
          case 1:    // one constrain on pool i
            if( Constr(0,i) != UNCONSTRAINED )
            {
              // calculate the identifible ones
              ConK(0,i,Upper) = Constr(0,i);
              ConK(0,i,Lower) = ConK(0,i,Upper);
              
              ConK(1,i,Upper) = -K(i,i) - ConK(0,i,Upper);
              ConK(1,i,Lower) = ConK(1,i,Upper); 
              
              ConK(i,1,Upper) = Gamma(i)/ConK(1,i,Lower);
              ConK(i,1,Lower) = ConK(i,1,Upper);
            }
            else if( Constr(1,i) != UNCONSTRAINED )
            {
              // calculate the identifible ones
              ConK(1,i,Upper) = Constr(1,i);
              ConK(1,i,Lower) = ConK(1,i,Upper); 
              
              ConK(i,1,Upper) = Gamma(i)/ConK(1,i,Lower);
              ConK(i,1,Lower) = ConK(i,1,Upper);
                                               
              ConK(0,i,Upper) = -K(i,i) - ConK(1,i,Upper);
              ConK(0,i,Lower) = ConK(0,i,Upper);
            }
            else  // Constr(i,1)
            {
              // calculate the identifible ones
              ConK(i,1,Upper) = Constr(i,1);
              ConK(i,1,Lower) = ConK(i,1,Upper);
                                               
              ConK(1,i,Upper) = Gamma(i)/ConK(i,1,Lower);
              ConK(1,i,Lower) = ConK(1,i,Upper); 
              
              ConK(0,i,Upper) = -K(i,i) - ConK(1,i,Upper);
              ConK(0,i,Lower) = ConK(0,i,Upper);
            }
            break;
          default:  // more than one constrain on pool i
            crash( "Conflicting constrains have been applied" );
            break;
        }        
      }
      
      // now check for K01 constraint
      if( Constr(0,1) != UNCONSTRAINED )
      {
        ConK(0,1,Upper) = Constr(0,1);
        ConK(0,1,Lower) = ConK(0,1,Upper);        
      }
      
      // fix up others
      MamFix();
    }
    else          // catenary
    { 
      // first check if there are over constraints
      // this part would probably be pulled out later
      for( i = 1; i <= N-1; i++ )
        if( Constr(i,i+1) != UNCONSTRAINED && Constr(i+1,i) != UNCONSTRAINED )
          crash( "This model is over constrained" );
      if( Constr(0,1) != UNCONSTRAINED && Constr(2,1) != UNCONSTRAINED )
        crash( "This model is over constrained" );
      for( i = 2; i <= N-1; i++ )
        if( Constr(0,i) != UNCONSTRAINED && 
          Constr(i-1,i) != UNCONSTRAINED &&
          Constr(i+1,i) != UNCONSTRAINED )
          crash( "This model is over constrained" );
      if( Constr(0,N) != UNCONSTRAINED && Constr(N-1,N) != UNCONSTRAINED )
        crash( "This model is over constrained" );    
    
      // check for constraints on koi,
      // since it's used for the rest of the recursion
      for( i = 1; i <= N; i++ )
        if( Constr(0,i) != UNCONSTRAINED )
        {
          ConK(0,i,Upper) = Constr(0,i);
          ConK(0,i,Lower) = ConK(0,i,Upper);
        }                                     
      
      // Now check for constraints on ki-1,i & ki,i-1
      for( i = 2; i <= N; i++ )
        if( Constr(i-1,i) != UNCONSTRAINED )
          if( Constr(i,i-1) != UNCONSTRAINED )
            crash( "This model is over-constrained." );
          else
          {       
            // only Ki-1,i is constrained
            ConK(i-1,i,Upper) = Constr(i-1,i);
            ConK(i-1,i,Lower) = ConK(i-1,i,Upper);
            
            ConK(i,i-1,Upper) = Gamma(i)/ConK(i-1,i,Lower);
            ConK(i,i-1,Lower) = ConK(i,i-1,Upper);
          }
        else
          if( Constr(i,i-1) != UNCONSTRAINED )
          {        
            // only Ki,i-1 is constrained
            ConK(i,i-1,Upper) = Constr(i,i-1);
            ConK(i,i-1,Lower) = ConK(i,i-1,Upper);
            
            ConK(i-1,i,Upper) = Gamma(i)/ConK(i,i-1,Lower);
            ConK(i-1,i,Lower) = ConK(i-1,i,Upper);
          }
          else {}

      // now check if any parameters associated with the first
      // pool or the last pool are constrained.  This allows
      // all other parameters associated with that pool to be solvable
      
      // checking pool 1
      int count = 0;
      if( Constr(0,1) != UNCONSTRAINED )
        count++;
      if( Constr(2,1) != UNCONSTRAINED )
        count++;
      if( Constr(1,2) != UNCONSTRAINED )
        count++;
      switch( count )
      {
        case 0:    // no constrain on pool 1
                    // don't need to adjust the bounds
          break;
        case 1:    // one constrain on pool 1
          if( Constr(0,1) != UNCONSTRAINED )
          {
            // calculate the identifible ones
            ConK(2,1,Upper) = -K(1,1) - ConK(0,1,Upper);
            ConK(2,1,Lower) = ConK(2,1,Upper);
            
            ConK(1,2,Lower) = Gamma(2)/ConK(2,1,Upper);
            ConK(1,2,Upper) = ConK(1,2,Lower);
          }
          else if( Constr(1,2) != UNCONSTRAINED )
          {
            // calculate the identifible ones
            ConK(0,1,Upper) = -K(1,1) - ConK(2,1,Lower);
            ConK(0,1,Lower) = ConK(0,1,Upper);
          }
          else  //Constr(2,1)
          {
            // calculate the identifible ones
            ConK(0,1,Upper) = -K(1,1) - ConK(2,1,Lower);
            ConK(0,1,Lower) = ConK(0,1,Upper);
          }
          break;
        default:  // more than one constrain on pool i
          crash( "Conflicting constrains have been applied to catenary pool #1" );
          break;
      }
      
      // checking pool n
      count = 0;
      if( Constr(0,N) != UNCONSTRAINED )
        count++;
      if( Constr(N-1,N) != UNCONSTRAINED )
        count++;
      if( Constr(N,N-1) != UNCONSTRAINED )
        count++;
      switch( count )
      {
        case 0:    // no constrain on pool 1
                    // don't need to adjust the bounds
          break;
        case 1:    // one constrain on pool 1
          if( Constr(0,N) != UNCONSTRAINED )
          {
            // calculate the identifible ones
            ConK(N-1,N,Upper) = -K(N,N) - ConK(0,N,Upper);
            ConK(N-1,N,Lower) = ConK(N-1,N,Upper);
            
            ConK(N,N-1,Lower) = Gamma(N)/ConK(N-1,N,Upper);
            ConK(N,N-1,Upper) = ConK(N,N-1,Lower);
          }
          else if( Constr(N,N-1) != UNCONSTRAINED )
          {
            // calculate the identifible ones
            ConK(0,N,Upper) = -K(N,N) - ConK(N-1,N,Lower);
            ConK(0,N,Lower) = ConK(0,N,Upper);
          }
          else  //Constr(N-1,N)
          {
            // calculate the identifible ones
            ConK(0,N,Upper) = -K(N,N) - ConK(N-1,N,Lower);
            ConK(0,N,Lower) = ConK(0,N,Upper);
          }
          break;
        default:  // more than one constrain on pool i
          crash( "Conflicting constrains have been applied to the last catenary pool" );
          break;
      }                  
      
      // Now apply the recursive relations to find bounds
      if( ConK(2,1,Upper) != ConK(2,1,Lower) )      // if not identified
      {
        ConK(2,1,Upper) = -K(1,1) - ConK(0,1,Lower);    
        ConK(1,2,Lower) = Gamma(2)/ConK(2,1,Upper);
      }
      for( i = 2; i <= N-1; i++ )                
        if( ConK(i+1,i,Lower) != ConK(i+1,i,Upper) )  // if not identified
        {
          ConK(i+1,i,Upper) = -K(i,i) - Gamma(i)/ConK(i,i-1,Upper) - ConK(0,i,Lower);
          ConK(i,i+1,Lower) = Gamma(i+1)/ConK(i+1,i,Upper);
          
          if( ConK(i,i-1,Upper) == ConK(i,i-1,Lower) && 
            ConK(0,i,Upper) == ConK(0,i,Lower) )
          {                                         
            // if the parameters above are identifible, then
            // so are these parameters
            ConK(i+1,i,Lower) = ConK(i+1,i,Upper);
            ConK(i,i+1,Upper) = ConK(i,i+1,Lower);
          }
        }                           
        
      if( ConK(N-1,N,Upper) != ConK(N-1,N,Lower) )    // if not identified
      {
        ConK(N-1,N,Upper) = -K(N,N) - ConK(0,N,Lower);  
        ConK(N,N-1,Lower) = Gamma(N)/ConK(N-1,N,Upper);
      }                                                   
      for( i = N-1; i >= 2; i-- )
        if( ConK(i-1,i,Lower) != ConK(i-1,i,Upper) )  // if not identified
        {
          ConK(i-1,i,Upper) = -K(i,i) - Gamma(i+1)/ConK(i,i+1,Upper) - ConK(0,i,Lower);
          ConK(i,i-1,Lower) = Gamma(i)/ConK(i-1,i,Upper);
          
          if( ConK(i,i+1,Upper) == ConK(i,i+1,Lower) && 
            ConK(0,i,Upper) == ConK(0,i,Lower) )
          {                                         
            // if the parameters above are identifible, then
            // so are these parameters
            ConK(i-1,i,Lower) = ConK(i-1,i,Upper);
            ConK(i,i-1,Upper) = ConK(i,i-1,Lower);
          }
        }        
                
      // compute k0i max
      if( ConK(0,1,Upper) != ConK(0,1,Lower) ) 
      {
        ConK(0,1,Upper) = -K(1,1) - ConK(2,1,Lower);
        if( ConK(2,1,Lower) == ConK(2,1,Upper) )
          ConK(0,1,Lower) = ConK(0,1,Upper);
      }
        
      for( i = 2; i <= N-1; i++ )
        if( ConK(0,i,Upper) != ConK(0,i,Lower) )    // if not identified
        {
          ConK(0,i,Upper) = -K(i,i) - ConK(i-1,i,Lower) - ConK(i+1,i,Lower);
          
          // if the neighboring pools are identifible, then this pool
          // must be identifible
          if( ConK(i-1,i,Upper) == ConK(i-1,i,Lower) &&
            ConK(i+1,i,Upper) == ConK(i+1,i,Lower) )
            ConK(0,i,Lower) = ConK(0,i,Upper);
        }
          
      if( ConK(0,N,Upper) != ConK(0,N,Lower) )      
      {
        ConK(0,N,Upper) = -K(N,N) - ConK(N-1,N,Lower);
        if( ConK(N-1,N,Upper) == ConK(N-1,N,Lower) )
          ConK(0,N,Lower) = ConK(0,N,Upper);
      }
    }    
  }             
  return( 0 );
}


// compute the bounds of compartment volumn
// algorithm is straight forward except the total volumn
int    Model::VQID()
{   
  // V1 is always identifible, and its formula is the same
  // for both mampool and catpool            
  V(1,Upper) = 1/Beta(N);
  V(1,Lower) = V(1,Upper);       
  ConV(1,Upper) = V(1,Upper);
  ConV(1,Lower) = V(1,Upper);
  Q(1,Upper) = C1*V(1,Upper);
  Q(1,Lower) = Q(1,Upper);
  ConQ(1,Upper) = Q(1,Upper);
  ConQ(1,Lower) = Q(1,Lower);
  
  if( N > 1 )
  {
    if( Type == Mam )
    {
        int i; // SLR
      for(i = 2; i <= N; i++ )
      {
        // unconstrained
        V(i,Lower) = -V(1,Upper)*K(i,1,Lower)/K(i,i);
        V(i,Upper) = -V(1,Upper)*K(i,1,Upper)/K(i,i);
        Q(i,Lower) = C1*V(i,Lower);
        Q(i,Upper) = C1*V(i,Upper);
                          
          // constrained
        ConV(i,Lower) = -V(1,Upper)*ConK(i,1,Lower)/K(i,i);
        ConV(i,Upper) = -V(1,Upper)*ConK(i,1,Upper)/K(i,i);
        ConQ(i,Lower) = C1*ConV(i,Lower);
        ConQ(i,Upper) = C1*ConV(i,Upper);
      }                                 
      
      // total                  
      // unconstrained
      V(N+1,Upper) = V(1,Upper) + V(N,Upper);
      for( i = 2; i <= N-1; i++ )
        V(N+1,Upper) += V(i,Lower);    
      if( K(0,1,Lower) == K(0,1,Upper) && K(0,1,Lower) == 0.0 )
      {                                  
        V(N+1,Upper) = 0.0;
        for( i = 1; i <= N; i++ )
          V(N+1,Upper) += V(i,Lower);
      }
      else
      {        
        V(N+1,Lower) = V(1,Upper) + V(2,Upper);
        for( i = 2; i <= N; i++ )
          V(N+1,Lower) += V(i,Lower);
      }      
      Q(N+1,Upper) = C1*V(N+1,Upper);
      Q(N+1,Lower) = C1*V(N+1,Lower);   
                      
      // constrained
      ConV(N+1,Upper) = ConV(1,Upper) + ConV(N,Upper);
      for( i = 2; i <= N-1; i++ )
        ConV(N+1,Upper) += ConV(i,Lower);    
      if( ConK(0,1,Lower) == ConK(0,1,Upper) && ConK(0,1,Lower) == 0.0 )
      {                                  
        ConV(N+1,Upper) = 0.0;
        for( i = 1; i <= N; i++ )
          ConV(N+1,Upper) += ConV(i,Lower);
      }
      else
      {        
        ConV(N+1,Lower) = ConV(1,Upper) + ConV(2,Upper);
        for( i = 2; i <= N; i++ )
          ConV(N+1,Lower) += ConV(i,Lower);
      }      
      ConQ(N+1,Upper) = C1*ConV(N+1,Upper);
      ConQ(N+1,Lower) = C1*ConV(N+1,Lower);
    }
    else
    {
      // Catenary algorithm

      float  detK2;
      float  sum1, sum2;
                  
      sum1 = 0.0F;
      sum2 = 0.0F;
        int i; // SLR
      // compute detK2
      for(  i = 1; i <= N; i++ )
      {
        sum1 += A(i)/L(i);
        sum2 += A(i);
      }          
      detK2 = sum1/sum2;
      for( i = 1; i <= N; i++ )
        detK2 *= L(i);                      
        
      Q(N,Upper) = Q(1,Upper)/detK2;
      Q(N,Lower) = Q(N,Upper);
      ConQ(N,Upper) = Q(N,Upper);
      ConQ(N,Lower) = Q(N,Upper);
      for( i = 1; i <= N-1; i++ )    
      {
        Q(N,Upper) *= K(i+1,i,Upper);
        Q(N,Lower) *= K(i+1,i,Lower);
        ConQ(N,Upper) *= ConK(i+1,i,Upper);
        ConQ(N,Lower) *= ConK(i+1,i,Lower);
      }
      
      Q(N-1,Upper) = Q(1,Upper)*K(N,N)/detK2; 
      
      Q(N-1,Lower) = Q(N-1,Upper);
      ConQ(N-1,Upper) = Q(N-1,Upper);
      ConQ(N-1,Lower) = Q(N-1,Upper);
      for( i = 1; i <= N-2; i++ )
      {
        Q(N-1,Upper) *= K(i+1,i,Upper);
        Q(N-1,Lower) *= K(i+1,i,Lower);
        ConQ(N-1,Upper) *= ConK(i+1,i,Upper);
        ConQ(N-1,Lower) *= ConK(i+1,i,Lower);
      }      
                                              
      if( N%2 == 0 )
      {
        Q(N,Upper) *= -1.0F;
        Q(N,Lower) *= -1.0F;
        ConQ(N,Upper) *= -1.0F;
        ConQ(N,Lower) *= -1.0F;
      }
      else
      {
        Q(N-1,Upper) *= -1.0F;
        Q(N-1,Lower) *= -1.0F;
        ConQ(N-1,Upper) *= -1.0F;
        ConQ(N-1,Lower) *= -1.0F;
      }              

      for( i = N-2; i >= 2; i-- )
      {
        Q(i,Upper) = (K(i+1,i+1)*Q(i+1,Upper) + K(i+1,i+2,Lower)*Q(i+2,Upper))/(-K(i+1,i,Upper));
        Q(i,Lower) = (K(i+1,i+1)*Q(i+1,Lower) + K(i+1,i+2,Upper)*Q(i+2,Lower))/(-K(i+1,i,Lower));
        ConQ(i,Upper) = (K(i+1,i+1)*ConQ(i+1,Upper) + ConK(i+1,i+2,Lower)*ConQ(i+2,Upper))/(-ConK(i+1,i,Upper));
        ConQ(i,Lower) = (K(i+1,i+1)*ConQ(i+1,Lower) + ConK(i+1,i+2,Upper)*ConQ(i+2,Lower))/(-ConK(i+1,i,Lower));
      }
                
      Q(N+1,Lower) = 0.0F;
      Q(N+1,Upper) = 0.0F;
      ConQ(N+1,Lower) = 0.0F;
      ConQ(N+1,Upper) = 0.0F;
      for( i = 1; i <= N; i++ )
      {                      
        Q(N+1,Lower) += Q(i,Lower);
        Q(N+1,Upper) += Q(i,Upper);
        ConQ(N+1,Lower) += ConQ(i,Lower);
        ConQ(N+1,Upper) += ConQ(i,Upper);
      }                               
      
      for( i = 1; i <= N+1; i++ )
      {
        V(i,Upper) = Q(i,Upper)/C1;
        V(i,Lower) = Q(i,Lower)/C1;
        ConV(i,Upper) = ConQ(i,Upper)/C1;
        ConV(i,Lower) = ConQ(i,Lower)/C1;
      }   
    }
  }
  
  return( 0 );
}


// compute bounds on mass flux
int  Model::MfluxID()
{  
  if( Type == Mam )
  { 
    float  factor = -K(1,1);
  
      int i; // SLR
    for( i = 2; i <= N; i++ )
      factor += Gamma(i)/K(i,i);
    
    for( i = 1; i <= N; i++ )
    {
      if( K(0,i,Upper) == K(0,i,Lower) && 
        Q(i,Upper) == Q(i,Lower) )           
      {
        Mflux(0,i,Upper) = K(0,i,Upper)*Q(i,Upper);
        Mflux(0,i,Lower) = Mflux(0,i,Upper);
      }                                       
      else  
      {
        Mflux(0,i,Upper) = Q(1,Upper)*factor;
        Mflux(0,i,Lower) = 0.0F;
      }                           
      if( ConK(0,i,Upper) == ConK(0,i,Lower) &&
        ConQ(i,Upper) == ConQ(i,Lower) )
      {
        ConMflux(0,i,Upper) = ConK(0,i,Upper)*ConQ(i,Upper);
        ConMflux(0,i,Lower) = ConMflux(0,i,Upper);
      }                                             
      else
      {
        ConMflux(0,i,Upper) = Mflux(0,i,Upper);
        ConMflux(0,i,Lower) = Mflux(0,i,Lower);
      }                                          
    }
    
    for( i = 2; i <= N; i++ )
    {
      Mflux(1,i,Lower) = -Q(1,Upper)*Gamma(i)/K(i,i);
      Mflux(1,i,Upper) = Mflux(1,i,Lower);
      ConMflux(1,i,Lower) = Mflux(1,i,Lower);
      ConMflux(1,i,Upper) = Mflux(1,i,Lower);
    }                                          
  
    for( i = 2; i <= N; i++ )
    {
      Mflux(i,1,Lower) = -K(i,i)*Q(i,Lower);
      Mflux(i,1,Upper) = -K(i,i)*Q(i,Upper);
      ConMflux(i,1,Lower) = -K(i,i)*ConQ(i,Lower);
      ConMflux(i,1,Upper) = -K(i,i)*ConQ(i,Upper);
    }
  }
  else    // catenary algorithm
  {
      int i;  // SLR
    for(  i = 1; i <= N-2; i++ )
    {
      Mflux(i+1,i,Upper) = K(i+1,i,Upper)*Q(i,Upper);
      Mflux(i+1,i,Lower) = K(i+1,i,Lower)*Q(i,Lower);
      ConMflux(i+1,i,Upper) = ConK(i+1,i,Upper)*ConQ(i,Upper);
      ConMflux(i+1,i,Lower) = ConK(i+1,i,Lower)*ConQ(i,Lower);
    }   

    if( N > 1 )
    {
       Mflux(N,N-1,Upper) = -K(N,N)*Q(N,Upper);
       Mflux(N,N-1,Lower) = -K(N,N)*Q(N,Lower);
       ConMflux(N,N-1,Upper) = -K(N,N)*ConQ(N,Upper);
       ConMflux(N,N-1,Lower) = -K(N,N)*ConQ(N,Lower);       
    }                                                 
    
    for( i = 2; i <= N; i++ )
    {
      Mflux(i-1,i,Upper) = Q(i,Upper)*K(i-1,i,Lower);
      Mflux(i-1,i,Lower) = Q(i,Lower)*K(i-1,i,Upper);
      ConMflux(i-1,i,Upper) = ConQ(i,Upper)*ConK(i-1,i,Lower);
      ConMflux(i-1,i,Lower) = ConQ(i,Lower)*ConK(i-1,i,Upper);
    }
    
    for( i = 1; i <= N; i++ )
    {
      Mflux(0,i,Lower) = 0.0F;
      Mflux(0,i,Upper) = K(0,i,Upper)*Q(i,Upper);
      ConMflux(0,i,Lower) = ConK(0,i,Lower)*ConQ(i,Lower);
      ConMflux(0,i,Upper) = ConK(0,i,Upper)*ConQ(i,Upper);
    }
  }

  return( 0 );
}


// Compute some global parameters, These include:
// pcr, pr, MRT's, and VD's
int Model::Global()
{
  float   AUC = 0.0F;    // area under the curve
  float   AUMC = 0.0F;   // area under the moment curve

    int i; // SLR
  for(  i = 1; i <= N; i++ )
  {
    AUC += -A(i)/L(i);
    AUMC += A(i)/(L(i) * L(i));
  }

  pcr = 1.0F/AUC;
  pr = pcr * C1;

  lowerMRT = AUMC/AUC;
  upperMRT = 0.0F;
  for( i = 1; i <= N; i++ )
    upperMRT += -1.0F/L(i);

  lowerVD = lowerMRT * pcr;
  upperVD = upperMRT * pcr;
  return( 0 );
}


/*
 * Compute the correlation matrix 
 */
int Model::GetCorr()
{
  // update the variance array
  if( SampleN() > 0 )
  {
    UpdateVar();

    FMatrix M = GetCovId();  
    FMatrix N = GetdFdP(1);  
    FMatrix ConN = GetdFdP(0);
    
    M.Print();
   
    cov = N*M*N.T();                //computes the unconstrained covariance matrix
    concov = ConN*M*ConN.T();       //computes the constrained covariance matrix
    
    corr = cov;                  // first set it to the covariance matrix   
		// compute the correlation matrix
      int ix; // SLR
    for( ix = 1; ix <= corr.XDim(); ix++ )
    	for( int iy = ix+1; iy <= corr.YDim(); iy++ )
    	{
				corr.Elem1( iy, ix ) /= (float) sqrt( corr.Elem1( iy, iy )*corr.Elem1( ix, ix ) ); 
				corr.Elem1( ix, iy ) = corr.Elem1( iy, ix );
		}

    for( ix = 1; ix <= corr.XDim(); ix++ )
			corr.Elem1( ix, ix ) = 1;
    	
  }
  else
  	corr.SetDim( 0, 0 );
  return( 0 );
}

/*
 * Return the sensitive function matrix 
 */
 
FMatrix Model::GetdFdP( int uorc )
{
  int n = PoolNum();
  int iy, ix;
  int i, j;
  FMatrix ret;
  
  if( OutputUnit() == Mass )
  	ret.SetDim( 5*n-4, 2*n-1 );
  else 
		ret.SetDim( 5*n-4, 2*n );

  if( MType() == Mam )
  {
    // Mammillary model
    // Matrix layout:
    //         K11, ..., Knn, G2, ..., Gn, 1/V1
    // K01max
    // ...
    // K0nmax
    // K12max
    // K12min
    // ...          MATRIX CONTENT
    // K1nmax
    // K1nmin
    // K21max
    // K21min
    // ...
    // Kn1max
    // Kn1min

    // K_0i^max
    for( ix = 1; ix <= 2*n-1; ix++ )
    {
      if( ix <= n )
      {
        // Kjj
        j = ix;

        iy = 1;

        // K0i max
        if( j == 1 )
          ret.Elem1( iy++, ix ) = -1;
        else
        {                                                                               
          if( uorc==1 )                                                                 //added clause
            ret.Elem1( iy++, ix ) = -Gamma(j)/square(K(j,j));
          else
            ret.Elem1( iy++, ix ) = -Gamma(j)/square(K(j,j));
        }

        for( i = 2; i <= n; i++ )
        {
          if( j == 1 )
          {
            if( uorc==1 )                                                                 //added clause
              ret.Elem1( iy++, ix ) = -Gamma(i)/square(K(i,1,Upper));
            else
              ret.Elem1( iy++, ix ) = -Gamma(i)/square(ConK(i,1,Upper));
          }
          else if( j == i )
            ret.Elem1( iy++, ix ) = -1;
          else
          {
            if ( uorc==1 )                                                                //added clause
              ret.Elem1( iy++, ix ) = -Gamma(i)*Gamma(j)/square(K(j,j)*K(i,1,Upper)); 
            else
              ret.Elem1( iy++, ix ) = -Gamma(i)*Gamma(j)/square(K(j,j)*ConK(i,1,Upper));
          }
        }

        // K1i max and K1i min
        for( i = 2; i <= n; i++ )
        {
          // max
          if( j == i )
            ret.Elem1( iy++, ix ) = -1;
          else
            ret.Elem1( iy++, ix ) = 0;

          // min
          if( j == 1 )
          {
            if ( uorc==1)                                                                    //added clause
              ret.Elem1( iy++, ix ) = Gamma(i)/square(K(i,1,Upper)); 
            else
              ret.Elem1( iy++, ix ) = Gamma(i)/square(ConK(i,1,Upper));
          }
          else if( j == i )
            ret.Elem1( iy++, ix ) = 0;
          else
          {
            if (uorc == 1)                                                                    //added clause
              ret.Elem1( iy++, ix ) = Gamma(i)*Gamma(j)/square(K(j,j)*K(i,1,Upper));
            else
              ret.Elem1( iy++, ix ) = Gamma(i)*Gamma(j)/square(K(j,j)*ConK(i,1,Upper));
          }
        }

        // Ki1max and Ki1min 
        for( i = 2; i <= n; i++ )
        {
          // max
          if( j == 1 )
            ret.Elem1( iy++, ix ) = -1;
          else if( j == i )
            ret.Elem1( iy++, ix ) = 0;
          else
          {                                                                                  //added clause
            if (uorc == 1)
              ret.Elem1( iy++, ix ) = -Gamma(j)/square(K(j,j)); 
            else
              ret.Elem1( iy++, ix ) = -Gamma(j)/square(K(j,j));
          }  

          // min
          if( j == i )
          {
            if  (uorc == 1 )                                                                 //added clause
              ret.Elem1( iy++, ix ) = Gamma(j)/square(K(j,j));
            else
              ret.Elem1( iy++, ix ) = Gamma(j)/square(K(j,j));
          }
          else
            ret.Elem1( iy++, ix ) = 0;
        }
        assert( iy == (5*n-4+1) );  // self-check
      }
      else
      {
        // Gj
        j = ix - n + 1;
        iy = 1;

        // K0i max
        if ( uorc == 1)                                                              //added clause
          ret.Elem1( iy++, ix ) = 1.0F/K(j,j); 
        else
          ret.Elem1( iy++, ix ) = 1.0F/K(j,j);

        for( i = 2; i <= n; i++ )
        {
          if( j == i )                                                              //added clause
          {
            if ( uorc ==1 )
              ret.Elem1( iy++, ix ) = -1/K(i,1,Upper); 
            else
              ret.Elem1( iy++, ix ) = -1/ConK(i,1,Upper);
          }
          else
          {
            if ( uorc ==1)                                                           //added clause
              ret.Elem1( iy++, ix ) = Gamma(i)/(K(j,j)*square(K(i,1,Upper))); 
            else
              ret.Elem1( iy++, ix ) = Gamma(i)/(K(j,j)*square(ConK(i,1,Upper)));
          }
        }

        // K1i max and K1i min
        for( i = 2; i <= n; i++ )
        {
          // max
          ret.Elem1( iy++, ix ) = 0;

          // min
          if( j == i )
          {                                                                           //added clause
            if ( uorc == 1 )
              ret.Elem1( iy++, ix ) = 1.0F/K(i,1,Upper); 
            else
              ret.Elem1( iy++, ix ) = 1.0F/ConK(i,1,Upper);
          }
          else
          {                                                                           //added clause
            if ( uorc == 1 )
              ret.Elem1( iy++, ix ) = -Gamma(i)/(K(j,j)*square(K(i,1,Upper))); 
            else
              ret.Elem1( iy++, ix ) = -Gamma(i)/(K(j,j)*square(ConK(i,1,Upper)));
          }
        }

        // Ki1max and Ki1min 
        for( i = 2; i <= n; i++ )
        {
          // max
          if( j == i )
            ret.Elem1( iy++, ix ) = 0;
          else                                                                        //added clause
          {
            if ( uorc == 1 )
              ret.Elem1( iy++, ix ) = 1.0F/K(j,j); 
            else
              ret.Elem1( iy++, ix ) = 1.0F/K(j,j);
          }

          // min
          if( j == i )
          {                                                                           //added clause
            if ( uorc == 1 )
              ret.Elem1( iy++, ix ) = -1.0F/K(j,j); 
            else
              ret.Elem1( iy++, ix ) = -1.0F/K(j,j);
          }
          else
            ret.Elem1( iy++, ix ) = 0;
        }
        assert( iy == (5*n-4+1) );
      }
    }

		if( OutputUnit() == Conc )
		{
			// 1/V1
			for( iy = 1; iy <= 5*n-4; iy++ )
				ret.Elem1( iy, 2*n ) = 0; 
   	}
   	
   	if (uorc == 0)                                                                //added clause
   	{
   	  for (int jx = 0; jx <= 3*n-2; jx++) 
   	  {
   	    if (constraint[jx] != UNCONSTRAINED)
   	    {
   	      for (int jy = 0; jy <= 2*n-2; jy++)
   	      {
   	        if ( jx <= n - 1 )
   	          ret.Elem(jx,jy) = 0;
   	        else 
   	        {
   	          int jz = 2*jx - n;
   	          ret.Elem(jz++,jy) = 0;
   	          ret.Elem(jz,jy) = 0;
   	        }
   	      }
   	    }
      }  
    }
    
        
  }
  else
  {
    // Catenary model
    // Matrix layout:
    //         K11, ..., Knn, G2, ..., Gn, 1/V1
    // K01max
    // ...
    // K0nmax
    // K12max
    // K12min
    // ...          MATRIX CONTENT
    // Kn-1,nmax
    // Kn-1,nmin
    // K21max
    // K21min
    // ...
    // Kn,n-1max
    // Kn,n-1min

    // These are more difficult, need temporary arrays
    // confused? so am I!

    // i always belong to the index that runs from 1 to n-1

    float **kkmax0 = matrix( 1, n, 1, n );    // dk0imax/dkjj
    float **kgmax0 = matrix( 1, n, 2, n );    // dk0imax/dgj
    float **kkmax1 = matrix( 1, n-1, 1, n );  // dki,i+1max/dkjj
    float **kkmin1 = matrix( 1, n-1, 1, n );   // dki,i+1min/dkjj
    float **kkmax2 = matrix( 1, n-1, 1, n );  // dki+1,imax/dkjj
    float **kkmin2 = matrix( 1, n-1, 1, n );  // dki+1,imin/dkjj
    float **kgmax1 = matrix( 1, n-1, 2, n );  // dki,i+1max/dgj
    float **kgmin1 = matrix( 1, n-1, 2, n );  // dki,i+1min/dgj
    float **kgmax2 = matrix( 1, n-1, 2, n );  // dki+1,imax/dgj
    float **kgmin2 = matrix( 1, n-1, 2, n );  // dki+1,imin/dgj

    // dki+1,imax/dkjj
    for( i = 1; i <= n-1; i++ )
      for( j = 1; j <= n; j++ )
      {
        if( j == i )
          kkmax2[i][j] = -1;
        else if( j < i )
        {
          float prod = 1.0;
          for( int r = j+1; r <= i; r++ )
            prod *= Gamma(r)/square(K(r,r-1,Upper));
          kkmax2[i][j] = -1*prod;
        }
        else
          kkmax2[i][j] = 0;
      }

    // dki+1,imax/dgj
    for( i = 1; i <= n-1; i++ )
      for( j = 2; j <= n; j++ )
      {
        if( j == i )
          kgmax2[i][j] = -1/K(i,i-1,Upper);
        else if( j < i )
        {
          float prod = 1.0F;
          for( int r = j+1; r <= i; r++ )
            prod *= Gamma(r)/square(K(r,r-1,Upper));
          kgmax2[i][j] = (-1/K(j,j-1,Upper))*prod;
        }
        else
          kgmax2[i][j] = 0;
      }

    // dki,i+1min/dkjj
    for( i = 1; i <= n-1; i++ )
      for( j = 1; j <= n; j++ )
      {
        if( j <= i )
          kkmin1[i][j] = -(Gamma(i+1)/square(K(i+1,i,Upper)))*kkmax2[i][j];
        else
          kkmin1[i][j] = 0;
      }

    // dki,i+1min/dgj
    for( i = 1; i <= n-1; i++ )
      for( j = 2; j <= n; j++ )
      {
        if( j > i+1 )
          kgmin1[i][j] = 0;
        else if( j == i+1 )
          kgmin1[i][j] = 1/K(i+1,i,Upper);
        else
          kgmin1[i][j] = -(Gamma(i+1)/square(K(i+1,i,Upper)))*kgmax2[i][j];
      }

    // dki,i+1max/dkjj
    for( i = 1; i <= n-1; i++ )
      for( j = 1; j <= n; j++ )
      {
        if( j == i+1 )
          kkmax1[i][j] = -1;
        else if( j > i+1 )
        {
          float prod = 1.0F;
          for( int r = i+1; r <= j-1; r++ )
            prod *= Gamma(r+1)/square(K(r,r+1,Upper));
          kkmax1[i][j] = -prod;
        }
        else
          kkmax1[i][j] = 0;
      }

    // dki,i+1max/dgj
    for( i = 1; i <= n-1; i++ )
      for( j = 2; j <= n; j++ )
      {
        if( j == i+2 )
          kgmax1[i][j] = -1/K(i,i+1,Upper);
        else if( j > i+2 )
        {
          float prod = 1.0F;
          for( int r = i+1; r <= j-1; r++ )
            prod *= Gamma(r+1)/square(K(r,r+1,Upper));
          kgmax1[i][j] = -prod/K(j-1,j,Upper);
        }
        else
          kgmax1[i][j] = 0;
      }

    // dki+1,imin/dkjj
    for( i = 1; i <= n-1; i++ )
      for( j = 1; j <= n; j++ )
        kkmin2[i][j] = -Gamma(i+1)*kkmax1[i][j]/square(K(i,i+1,Upper));

    // dki+1,i,min/dgj
    for( i = 1; i <= n-1; i++ )
      for( j = 2; j <= n; j++ )
        if( j == i+1 )
          kgmin2[i][j] = 1/K(i,i+1,Upper);
        else
          kgmin2[i][j] = -Gamma(i+1)*kgmax1[i][j]/square(K(i,i+1,Upper));

    // dk0imax/dkjj
    for( i = 1; i <= n; i++ )
      for( j = 1; j <= n; j++ )
      {
        if( i == 1 )
          kkmax0[i][j] = -delta(1,j) - kkmin2[1][j];
        else if( i < n )

          kkmax0[i][j] = -delta(i,j) - kkmin1[i-1][j] - kkmin2[i][j];
        else
          kkmax0[i][j] = -delta(n,j) - kkmin1[n-1][j];
      }

    // dk0imax/dgj
    for( i = 1; i <= n; i++ )
      for( j = 2; j <= n; j++ )
      {
        if( i == 1 )
          kgmax0[i][j] = -kgmin2[1][j];
        else if( i < n )
          kgmax0[i][j] = -kgmin1[i-1][j] - kgmin2[i][j];
        else
          kgmax0[i][j] = -kgmin1[n-1][n];
      } 

    // assign results to the matrix
    for( ix = 1; ix <= 2*n-1; ix++ )
    {
      if( ix <= n )
      {
        // Kii
        j = ix;

        iy = 1;

        // k0imax
        for( i = 1; i <= n; i++ )
          ret.Elem1( iy++, ix ) = kkmax0[i][j];

        // ki,i+1max, ki,i+1min
        for( i = 1; i <= n-1; i++ )
        {
          ret.Elem1( iy++, ix ) = kkmax1[i][j];
          ret.Elem1( iy++, ix ) = kkmin1[i][j];
        }

        // ki+1,imax, ki+1,imin
        for( i = 1; i <= n-1; i++ )
        {
          ret.Elem1( iy++, ix ) = kkmax2[i][j];
          ret.Elem1( iy++, ix ) = kkmin2[i][j];
        }

        assert( iy == (5*n-4+1) );
      }
      else
      {
        // Gj
        j = ix - n + 1;

        iy = 1;

        // k0imax
        for( i = 1; i <= n; i++ )
          ret.Elem1( iy++, ix ) = kgmax0[i][j];

        // ki,i+1max, ki,i+1min
        for( i = 1; i <= n-1; i++ )
        {
          ret.Elem1( iy++, ix ) = kgmax1[i][j];
          ret.Elem1( iy++, ix ) = kgmin1[i][j];
        }

        // ki+1,imax, ki+1,imin
        for( i = 1; i <= n-1; i++ )
        {
          ret.Elem1( iy++, ix ) = kgmax2[i][j];
          ret.Elem1( iy++, ix ) = kgmin2[i][j];
        }

        assert( iy == (5*n-4+1) );
      }
    }

		if( OutputUnit() == Conc )
		{
			// 1/V1
			for( iy = 1; iy <= 5*n-4; iy++ )
				ret.Elem1( iy, 2*n ) = 0; 
		}

    free_matrix( kkmax0, 1, n, 1, n );
    free_matrix( kgmax0, 1, n, 2, n );
    free_matrix( kkmax1, 1, n-1, 1, n );
    free_matrix( kkmin1, 1, n-1, 1, n );
    free_matrix( kkmax2, 1, n-1, 1, n );
    free_matrix( kkmin2, 1, n-1, 1, n );
    free_matrix( kgmax1, 1, n-1, 2, n );
    free_matrix( kgmin1, 1, n-1, 2, n );
    free_matrix( kgmax2, 1, n-1, 2, n );
    free_matrix( kgmin2, 1, n-1, 2, n );
  }

  return( ret );
}




/*
 * Return the sensitive function matrix 

FMatrix Model::GetdFdP()
{
  int n = PoolNum();
  int iy, ix;
  int i, j;
  FMatrix ret;

  if( OutputUnit() == Mass )
  	ret.SetDim( 5*n-4, 2*n-1 );
  else 
		ret.SetDim( 5*n-4, 2*n ); 
/*
 Something wrong HERE get a revised copy of this from 
	David Schiebel at some point in time.  Now everything is 
	just going to calculated as a catenary model correlation
	matrix

  - S. Russell



  if( MType() == Mam )
  {
	  goto CATENARY_PART;  // ***** GOTO FUNCTION TO CAT ******
    // Mammillary model
    // Matrix layout:
    //         K11, ..., Knn, G2, ..., Gn, 1/V1
    // K01max
    // ...
    // K0nmax
    // K12max
    // K12min
    // ...          MATRIX CONTENT
    // K1nmax
    // K1nmin
    // K21max
    // K21min
    // ...
    // Kn1max
    // Kn1min

    // K_0i^max
    for( ix = 1; ix <= 2*n-1; ix++ )
    {
      if( ix <= n )
      {
        // Kjj
        j = ix;

        iy = 1;

        // K0i max

		// Added by Tuong
		//====================================
	
		if(Constr(0,j) != UNCONSTRAINED)	
		{
			ret.Elem1( iy++, ix ) = 0;		//constrained 
		}
		else
		{
		//====================================
			if( j == 1 )
			{
				ret.Elem1( iy++, ix ) = -1;
			}
			else
			{
				ret.Elem1( iy++, ix ) = -Gamma(j)/square(K(j,j));
			}
		}



        for( i = 2; i <= n; i++ )
        {
			// Added by Tuong
			//====================================
			if(Constr(0,i) != UNCONSTRAINED)	
			{
				ret.Elem1( iy++, ix ) = 0;		//constrained 
			}
			else
			{
			//====================================

        //  if( j == 1 )
        //    ret.Elem1( iy++, ix ) = -Gamma(i)/square(K(i,1,Upper));
        //  else if( j == i )
		  if(j == i)
            ret.Elem1( iy++, ix ) = -1;
          else
            ret.Elem1( iy++, ix ) = -Gamma(i)*Gamma(j)/square(K(j,j)*K(i,1,Upper));

			}
        }

        // K1i max and K1i min
        for( i = 2; i <= n; i++ )
        {
			// Added by Tuong
			//====================================
			if(Constr(1,i) != UNCONSTRAINED)	
			{
				//max
				ret.Elem1( iy++, ix ) = 0;		//constrained 
				//min
				ret.Elem1( iy++, ix ) = 0;		//constrained 
			}
			else
			{
			//====================================

				// max
				if( j == i )
				   ret.Elem1( iy++, ix ) = -1;
				 else
				   ret.Elem1( iy++, ix ) = 0;

				 // min
				 // Added by Tuong
				//====================================
				if(Constr(i, 1) != UNCONSTRAINED)	//if Ki1 max is constrained
				{
					//constrained 
					ret.Elem1( iy++, ix ) = 0;		// = 0 for ix<=n	
				}
				else
				{
				//====================================

					if( j == 1 )
						ret.Elem1( iy++, ix ) = Gamma(i)/square(K(i,1,Upper));
					else if( j == i )
						ret.Elem1( iy++, ix ) = 0;
					else
						ret.Elem1( iy++, ix ) = Gamma(i)*Gamma(j)/square(K(j,j)*K(i,1,Upper));
				}
			}
        }

        // Ki1max and Ki1min 
        for( i = 2; i <= n; i++ )
        {
			// Added by Tuong
			//====================================
			if(Constr(i,1) != UNCONSTRAINED)	
			{
				//max
				ret.Elem1( iy++, ix ) = 0;		//constrained 
				//min
				ret.Elem1( iy++, ix ) = 0;		//constrained 
			}
			else
			{
			//====================================
				// max
				if( j == 1 )
					ret.Elem1( iy++, ix ) = -1;
				else if( j == i )
					ret.Elem1( iy++, ix ) = 0;
				else
					ret.Elem1( iy++, ix ) = -Gamma(j)/square(K(j,j));

				// min
				if( j == i )
					ret.Elem1( iy++, ix ) = Gamma(j)/square(K(j,j));
				else
					ret.Elem1( iy++, ix ) = 0;
			}
        }
        assert( iy == (5*n-4+1) );  // self-check
      }
      else
      {
        // Gj
        j = ix - n + 1;
        iy = 1;

        // K0i max

		// Added by Tuong
		//====================================
		if(Constr(0,i) != UNCONSTRAINED)	
		{
			ret.Elem1( iy++, ix ) = 0;		//constrained 
		}
		else
		{
		//====================================
			ret.Elem1( iy++, ix ) = 1.0F/K(j,j);
		}

        for( i = 2; i <= n; i++ )
        {
			// Added by Tuong
			//====================================
			if(Constr(0,i) != UNCONSTRAINED)	
			{
				ret.Elem1( iy++, ix ) = 0;		//constrained 
			}
			else
			{
			//====================================
				if( j == i )
					ret.Elem1( iy++, ix ) = -1/K(i,1,Upper);
				else
					ret.Elem1( iy++, ix ) = Gamma(i)/(K(j,j)*square(K(i,1,Upper)));
			}
        }

        // K1i max and K1i min
        for( i = 2; i <= n; i++ )
        {
			// Added by Tuong
			//====================================
			if(Constr(1,i) != UNCONSTRAINED)	
			{
				//max
				ret.Elem1( iy++, ix ) = 0;		//constrained 
				//min
				ret.Elem1( iy++, ix ) = 0;		//constrained 
			}
			else
			{
			//====================================

				// max
		        ret.Elem1( iy++, ix ) = 0;

				// min
				if( j == i )
					ret.Elem1( iy++, ix ) = 1.0F/K(i,1,Upper);
				else
				{
					// Added by Tuong
					//====================================
					if(Constr(i,1) != UNCONSTRAINED)	
					{
						ret.Elem1( iy++, ix ) = 0;
					}
					else
					{
					//====================================
						ret.Elem1( iy++, ix ) = -Gamma(i)/(K(j,j)*square(K(i,1,Upper)));
					}
				}
			}
        }

        // Ki1max and Ki1min 
        for( i = 2; i <= n; i++ )
        {
			// Added by Tuong
			//====================================
			if(Constr(i,1) != UNCONSTRAINED)	
			{
				//max
				ret.Elem1( iy++, ix ) = 0;		//constrained 
				//min
				ret.Elem1( iy++, ix ) = 0;		//constrained 
			}
			else
			{
			//====================================

				// max
				if( j == i )
					ret.Elem1( iy++, ix ) = 0;
				else
					ret.Elem1( iy++, ix ) = 1.0F/K(j,j);

				// min
				if( j == i )
					ret.Elem1( iy++, ix ) = -1.0F/K(j,j);
				else
					ret.Elem1( iy++, ix ) = 0;
			}
        }
        assert( iy == (5*n-4+1) );
      }
    }

		if( OutputUnit() == Conc )
		{
			// 1/V1
			for( iy = 1; iy <= 5*n-4; iy++ )
				ret.Elem1( iy, 2*n ) = 0; 
   	}
  }
  else
  {
CATENARY_PART:
    // Catenary model
    // Matrix layout:
    //         K11, ..., Knn, G2, ..., Gn, 1/V1
    // K01max
    // ...
    // K0nmax
    // K12max
    // K12min
    // ...          MATRIX CONTENT
    // Kn-1,nmax
    // Kn-1,nmin
    // K21max
    // K21min
    // ...
    // Kn,n-1max
    // Kn,n-1min

    // These are more difficult, need temporary arrays
    // confused? so am I!

    // i always belong to the index that runs from 1 to n-1

    float **kkmax0 = matrix( 1, n, 1, n );    // dk0imax/dkjj
    float **kgmax0 = matrix( 1, n, 2, n );    // dk0imax/dgj
    float **kkmax1 = matrix( 1, n-1, 1, n );  // dki,i+1max/dkjj
    float **kkmin1 = matrix( 1, n-1, 1, n );   // dki,i+1min/dkjj
    float **kkmax2 = matrix( 1, n-1, 1, n );  // dki+1,imax/dkjj
    float **kkmin2 = matrix( 1, n-1, 1, n );  // dki+1,imin/dkjj
    float **kgmax1 = matrix( 1, n-1, 2, n );  // dki,i+1max/dgj
    float **kgmin1 = matrix( 1, n-1, 2, n );  // dki,i+1min/dgj
    float **kgmax2 = matrix( 1, n-1, 2, n );  // dki+1,imax/dgj
    float **kgmin2 = matrix( 1, n-1, 2, n );  // dki+1,imin/dgj

    // dki+1,imax/dkjj
    for( i = 1; i <= n-1; i++ )
      for( j = 1; j <= n; j++ )
      {
		// Added by Tuong
		//====================================
		if(Constr(i+1,i) != UNCONSTRAINED)
		{
			kkmax2[i][j] = 0;
		}
		else
		{
		//====================================
			if( j == i )
				kkmax2[i][j] = -1;
			else if( j < i )
			{
				// Added by Tuong
				//====================================
				if(Constr(i,i-1) != UNCONSTRAINED)
				{
					kkmax2[i][j] = 0;
				}
				else
				{
				//====================================
					float prod = 1.0;
					for( int r = j+1; r <= i; r++ )
						prod *= Gamma(r)/square(K(r,r-1,Upper));
					kkmax2[i][j] = -1*prod;
				}
			}
			else
				kkmax2[i][j] = 0;
		}
      }

    // dki+1,imax/dgj
    for( i = 1; i <= n-1; i++ )
      for( j = 2; j <= n; j++ )
      {
	    // Added by Tuong
		//====================================
		if(Constr(i+1,i) != UNCONSTRAINED)
		{
			kgmax2[i][j] = 0;
		}
		else
		{
		//====================================
			if( j == i )
				kgmax2[i][j] = -1/K(i,i-1,Upper);
			else if( j < i )
			{
				// Added by Tuong
				//====================================
				if(Constr(i,i-1) != UNCONSTRAINED)
				{
					kgmax2[i][j]  = 0;
				}
				else
				{
				//====================================
					float prod = 1.0F;
					for( int r = j+1; r <= i; r++ )
						prod *= Gamma(r)/square(K(r,r-1,Upper));
					kgmax2[i][j] = (-1/K(j,j-1,Upper))*prod;
				}
			}
			else
				kgmax2[i][j] = 0;
		}
      }

    // dki,i+1min/dkjj
    for( i = 1; i <= n-1; i++ )
      for( j = 1; j <= n; j++ )
      {
		// Added by Tuong
		//====================================
		if(Constr(i,i+1) != UNCONSTRAINED)
		{
			kkmin1[i][j] = 0;
		}
		else
		{
		//====================================
			if( j <= i )
			{
				// Added by Tuong
				//====================================
				if(Constr(i+1,i) != UNCONSTRAINED)
				{
					kkmin1[i][j] = 0;
				}
				else
				{
				//====================================
					kkmin1[i][j] = -(Gamma(i+1)/square(K(i+1,i,Upper)))*kkmax2[i][j];
				}
			}
			else
				kkmin1[i][j] = 0;
		}
      }

    // dki,i+1min/dgj
    for( i = 1; i <= n-1; i++ )
      for( j = 2; j <= n; j++ )
      {
		// Added by Tuong
		//====================================
		if(Constr(i,i+1) != UNCONSTRAINED)
		{
			kgmin1[i][j] = 0;
		}
		else
		{
		//====================================
			if( j > i+1 )
				kgmin1[i][j] = 0;
			else if( j == i+1 )
				kgmin1[i][j] = 1/K(i+1,i,Upper);
			else
			{
				// Added by Tuong
				//====================================
				if(Constr(i+1,i) != UNCONSTRAINED)
				{
					kgmin1[i][j] = 0;
				}
				else
				{
				//====================================
					kgmin1[i][j] = -(Gamma(i+1)/square(K(i+1,i,Upper)))*kgmax2[i][j];
				}
			}
		}
      }

    // dki,i+1max/dkjj
    for( i = 1; i <= n-1; i++ )
      for( j = 1; j <= n; j++ )
      {
		// Added by Tuong
		//====================================
		if(Constr(i,i+1) != UNCONSTRAINED)
		{
			kkmax1[i][j] = 0;
		}
		else
		{
		//====================================
			if( j == i+1 )
				kkmax1[i][j] = -1;
			else if( j > i+1 )
			{
				// Added by Tuong
				//====================================
				if(Constr(i,i+1) != UNCONSTRAINED)
				{
					kkmax1[i][j] = 0;
				}
				else
				{
				//====================================
					float prod = 1.0F;
					for( int r = i+1; r <= j-1; r++ )
						prod *= Gamma(r+1)/square(K(r,r+1,Upper));
					kkmax1[i][j] = -prod;
				}
			}
			else
				kkmax1[i][j] = 0;
		}
      }

    // dki,i+1max/dgj
    for( i = 1; i <= n-1; i++ )
      for( j = 2; j <= n; j++ )
      {
		// Added by Tuong
		//====================================
		if(Constr(i,i+1) != UNCONSTRAINED)
		{
			kgmax1[i][j] = 0;
		}
		else
		{
		//====================================
			if( j == i+2 )
				kgmax1[i][j] = -1/K(i,i+1,Upper);
			else if( j > i+2 )
			{
				// Added by Tuong
				//====================================
				if(Constr(i,i+1) != UNCONSTRAINED)
				{
					kgmax1[i][j] = 0;
				}
				else
				{
				//====================================
					float prod = 1.0F;
					for( int r = i+1; r <= j-1; r++ )
						prod *= Gamma(r+1)/square(K(r,r+1,Upper));
					kgmax1[i][j] = -prod/K(j-1,j,Upper);
				}
			}
			else
				kgmax1[i][j] = 0;
		}
      }

    // dki+1,imin/dkjj
    for( i = 1; i <= n-1; i++ )
      for( j = 1; j <= n; j++ )
	  {
		// Added by Tuong
		//====================================
		if(Constr(i+1,i) != UNCONSTRAINED)
		{
			kkmin2[i][j] = 0;
		}
		else
		{
			if(Constr(i,i+1) != UNCONSTRAINED)
			{
				kkmin2[i][j] = 0;
			}
			else
			{
		//====================================
				kkmin2[i][j] = -Gamma(i+1)*kkmax1[i][j]/square(K(i,i+1,Upper));
			}
		}
	  }

    // dki+1,i,min/dgj
    for( i = 1; i <= n-1; i++ )
      for( j = 2; j <= n; j++ )
	  {
		// Added by Tuong
		//====================================
		if(Constr(i+1,i) != UNCONSTRAINED)
		{
			kgmin2[i][j] = 0;
		}
		else
		{
		//====================================
			if( j == i+1 )
				kgmin2[i][j] = 1/K(i,i+1,Upper);
			else
			{
				// Added by Tuong
				//====================================
				if(Constr(i,i+1) != UNCONSTRAINED)
				{
					kgmin2[i][j] = 0;
				}
				else
				{
				//====================================
					kgmin2[i][j] = -Gamma(i+1)*kgmax1[i][j]/square(K(i,i+1,Upper));
				}
			}
		}
	  }

    // dk0imax/dkjj
    for( i = 1; i <= n; i++ )
      for( j = 1; j <= n; j++ )
      {
		// Added by Tuong
		//====================================
		if(Constr(0,i) != UNCONSTRAINED)
		{
			kkmax0[i][j] = 0;
		}
		else
		{
		//====================================
			if( i == 1 )
			{
				// Added by Tuong
				//====================================
				if(Constr(2,1) != UNCONSTRAINED)
				{
					kkmax0[i][j] = -(float)delta(1,j);
				}
				else
				{
				//====================================
					kkmax0[i][j] = -delta(1,j) - kkmin2[1][j];
				}
			}
			else if( i < n )
			{
				// Added by Tuong
				//====================================
				if((Constr(i-1,i) != UNCONSTRAINED) && (Constr(i+1,i) != UNCONSTRAINED))
				{
					kkmax0[i][j] = -(float)delta(i,j);
				}
				else if(Constr(i-1,i) != UNCONSTRAINED)
				{
					kkmax0[i][j] = -delta(i,j) - kkmin2[i][j];
				}
				else if(Constr(i+1,i) != UNCONSTRAINED)
				{
					kkmax0[i][j] = -delta(i,j) - kkmin1[i-1][j];
				}
				else
				{
				//====================================
					kkmax0[i][j] = -delta(i,j) - kkmin1[i-1][j] - kkmin2[i][j];
				}
			}
			else
			{
				// Added by Tuong
				//====================================
				if(Constr(n-1,n) != UNCONSTRAINED)
				{
					kkmax0[i][j] = -(float)delta(n,j);
				}
				else
				{
				//====================================
					kkmax0[i][j] = -delta(n,j) - kkmin1[n-1][j];
				}
			}
		}
      }

    // dk0imax/dgj
    for( i = 1; i <= n; i++ )
      for( j = 2; j <= n; j++ )
      {
		// Added by Tuong
		//====================================
		if(Constr(0,i) != UNCONSTRAINED)
		{
			kgmax0[i][j] = 0;
		}
		else
		{
		//====================================
			if( i == 1 )
			{
				// Added by Tuong
				//====================================
				if(Constr(2,1) != UNCONSTRAINED)
				{
					kgmax0[i][j] = 0;
				}
				else
				{
				//====================================
					kgmax0[i][j] = -kgmin2[1][j];
				}
			}
			else if( i < n )
			{
				// Added by Tuong
				//====================================
				if((Constr(i-1,1) != UNCONSTRAINED) && (Constr(i+1,1) != UNCONSTRAINED))
				{
					kgmax0[i][j] = 0;
				}
				else if(Constr(i-1,1) != UNCONSTRAINED)
				{
					kgmax0[i][j] = -kgmin2[i][j];
				}
				else if(Constr(i+1,1) != UNCONSTRAINED)
				{
					kgmax0[i][j] = -kgmin1[i-1][j];
				}
				else
				{
				//====================================
					kgmax0[i][j] = -kgmin1[i-1][j] - kgmin2[i][j];
				}
			}
			else
			{
				// Added by Tuong
				//====================================
				if(Constr(n-1,n) != UNCONSTRAINED)
				{
					kgmax0[i][j] = 0;
				}
				else
				{
				//====================================
					kgmax0[i][j] = -kgmin1[n-1][n];
				}
			}
		}
      } 

    // assign results to the matrix
    for( ix = 1; ix <= 2*n-1; ix++ )
    {
      if( ix <= n )
      {
        // Kii
        j = ix;

        iy = 1;

        // k0imax
        for( i = 1; i <= n; i++ )
          ret.Elem1( iy++, ix ) = kkmax0[i][j];

        // ki,i+1max, ki,i+1min
        for( i = 1; i <= n-1; i++ )
        {
          ret.Elem1( iy++, ix ) = kkmax1[i][j];
          ret.Elem1( iy++, ix ) = kkmin1[i][j];
        }

        // ki+1,imax, ki+1,imin
        for( i = 1; i <= n-1; i++ )
        {
          ret.Elem1( iy++, ix ) = kkmax2[i][j];
          ret.Elem1( iy++, ix ) = kkmin2[i][j];
        }

        assert( iy == (5*n-4+1) );
      }
      else
      {
        // Gj
        j = ix - n + 1;

        iy = 1;

        // k0imax
        for( i = 1; i <= n; i++ )
          ret.Elem1( iy++, ix ) = kgmax0[i][j];

        // ki,i+1max, ki,i+1min
        for( i = 1; i <= n-1; i++ )
        {
          ret.Elem1( iy++, ix ) = kgmax1[i][j];
          ret.Elem1( iy++, ix ) = kgmin1[i][j];
        }

        // ki+1,imax, ki+1,imin
        for( i = 1; i <= n-1; i++ )
        {
          ret.Elem1( iy++, ix ) = kgmax2[i][j];
          ret.Elem1( iy++, ix ) = kgmin2[i][j];
        }

        assert( iy == (5*n-4+1) );
      }
    }

		if( OutputUnit() == Conc )
		{
			// 1/V1
			for( iy = 1; iy <= 5*n-4; iy++ )
				ret.Elem1( iy, 2*n ) = 0; 
		}

    free_matrix( kkmax0, 1, n, 1, n );
    free_matrix( kgmax0, 1, n, 2, n );
    free_matrix( kkmax1, 1, n-1, 1, n );
    free_matrix( kkmin1, 1, n-1, 1, n );
    free_matrix( kkmax2, 1, n-1, 1, n );
    free_matrix( kkmin2, 1, n-1, 1, n );
    free_matrix( kgmax1, 1, n-1, 2, n );
    free_matrix( kgmin1, 1, n-1, 2, n );
    free_matrix( kgmax2, 1, n-1, 2, n );
    free_matrix( kgmin2, 1, n-1, 2, n );
  }

  return( ret );
}
*/

/*
 * Returns the covariance matrix of the identifiable parameters
 */
FMatrix Model::GetCovId()
{
  FMatrix M;        // the Fisher information matrix
  int i;

  M = GetM( 1 );
  for( i = 2; i <= SampleN(); i++ )
    M += GetM( i );

  return( M^(-1) );
}


/*
 * Compute Mi
 */
FMatrix Model::GetM( int i )
{
  FMatrix dYdP = GetdYdP( SampleTime(i) );

//	TRACE( "time = %d\n", i );
  dYdP.Print();

  return( (1/Var(i))*dYdP.T()*dYdP );
}


/*
 * Compute the dYdP row vector
 * parameters are ordered as:
 * 
 * K11, ..., Knn, G2, ..., Gn, 1/V1
 */
FMatrix Model::GetdYdP( float time )
{
	int ParaNum;

	if( OutputUnit() == Mass )
		ParaNum = 2*PoolNum() - 1;
	else
		ParaNum = 2*PoolNum();    		// number of parameters

  FMatrix ret( 1, ParaNum );        // the row vector returned

  for( int i = 1; i <= ParaNum; i++ )
    ret.Elem1( 1, i ) = GetdYdPij( time, i ); 

  return( ret );
}


/*
 * Compute an element of the dYdP row vector
 * time = sample time
 * paramI = parameter index
 * Parameters are ordered as:
 *   1   2 ...   n n+1 ... 2n-1 2n
 * K11 K22 ... Knn  G2 ... Gn   1/V1
 *
 * Code layout is made to be as close to the formula as possible
 * Formula is based on the index of Robert Lindell's paper.
 */
float Model::GetdYdPij( float time, int paramI )
{
  float ret = 0.0F;       // returned value
  int n = PoolNum();      // number of pools
  float t = time;         // makes things close to the formula
  float V1 = V(1, Upper);	// V1 is identifiable anyway

  // range check
  if( OutputUnit() == Mass )
		assert( paramI >= 1 && paramI <= 2*n-1 );
  else
		assert( paramI >= 1 && paramI <= 2*n );

  if( MType() == Mam )
  {
    // Mammillary model
    if( paramI == 1 )
    {
      // parameter is K11
      int i, j;
      float sum_i;
      float sum_j;

      sum_i = 0.0F;
      for( i = 1; i <= n; i++ )
      {
        sum_j = 0.0F;
        for( j = i + 1; j <= n; j++ )
          sum_j += (float) (A(i)*A(j)*(exp(L(i)*t) - exp( L(j)*t ))/(L(i) - L(j)));

        sum_i += (float) (square(A(i))*t*exp(L(i)*t) + 2*sum_j);
      }

      ret = sum_i;

#if LINDELL
			ret *= V1;
#else
			if( OutputUnit() == Conc )
				ret /= V1;
#endif
    }
    else if( paramI <= n )
    {
      // parameter is Kii, i > 1
      int i = paramI;    // making these indexes the same as the formula    
      int j;
      int k;
      float sum_j;
      float sum_k;

      sum_j = 0;
      for( j = 1; j <= n; j++ )
      {
        sum_k = 0;
        for( k = j + 1; k <= n; k++ )
          sum_k += (float) (A(j)*A(k)*( 
            (exp(L(j)*t)/square(L(j) - K(i,i)) - 
             exp(L(k)*t)/square(L(k) - K(i,i)))/(L(j) - L(k)) +
             (exp(K(i,i)*t)/((L(j) - K(i,i))*(L(k) - K(i,i))))*
             (t + (L(j) + L(k) - 2*K(i,i))/((L(k) - K(i,i))*(L(j) - K(i,i))))));

        sum_j += (float) (square(A(j)/(L(j) - K(i,i)))*
                 (t*exp(L(j)*t) + t*exp(K(i,i)*t) + 
                  2*(exp(K(i,i)*t) - exp(L(j)*t))/(L(j) - K(i,i))) + 2*sum_k);
      }

			ret = Gamma(i)*sum_j;

#if LINDELL
			ret *= V1;
#else 
			if( OutputUnit() == Conc )
				ret /= V1;
#endif
    }
    else if( paramI != 2*n )
    {
      // parameter is Gamma
      int i = paramI - n + 1;
      int j;
      int k;
      float sum_j;
      float sum_k;

      sum_j = 0;
      for( j = 1; j <= n; j++ )
      {
        sum_k = 0;
        for( k = j + 1; k <= n; k++ )
          sum_k += (float) (A(j)*A(k)*((exp(L(j)*t)/(L(j) - K(i,i)) - 
                                       exp(L(k)*t)/(L(k) - K(i,i)))/(L(j) - L(k)) +
                                      exp(K(i,i)*t)/((L(j) - K(i,i))*(L(k) - K(i,i)))));

        sum_j += (float) (square(A(j))*(t*exp(L(j)*t) + 
                              (exp(K(i,i)*t) - exp(L(j)*t))/(L(j) - K(i,i)))/
                (L(j) - K(i,i)) + 2*sum_k);
      }

     	ret = sum_j; 

#if LINDELL
			ret *= V1;
#else
			if( OutputUnit() == Conc )
				ret /= V1;
#endif
    }
    else
    {
    	assert( OutputUnit() == Conc );

    	// parameter is 1/V1
    	ret = SumOfExp( t );
    }
  }
  else
  {
    // Catenary model
    if( paramI <= n )
    {
      // parameter is Kii
      int i = paramI;    // making these indexes the same as the formula    
      int j;
      float prod_j = 1.0F;

      if( i > 1 )
        for( j = 2; j <= i; j++ )
          prod_j *= Gamma(j);

      ret = prod_j*GetDiDj( i+1, i+1, t );

#if LINDELL
			ret *= V1;
#else
      if( OutputUnit() == Conc )
				ret /= V1;
#endif
    }
    else if( paramI != 2*n ) 
    {
      // parameter is Gamma
      int i = paramI - n + 1;
      int j;
      float prod_j = 1.0F;

      if( i > 2 )
        for( j = 2; j <= i - 1; j++ )
          prod_j *= Gamma(j);

      ret = prod_j*GetDiDj( i, i+1, t );

#if LINDELL
			ret *= V1;
#else
      if( OutputUnit() == Conc )
      	ret /= V1;
#endif
    }
    else
    {
			assert( OutputUnit() == Conc );

    	// parameter is 1/V1
    	ret = SumOfExp( t );
    }
  }

  return( ret );
}


/*
 *   returns Lap^-1( D_i*D_j/D_1^2 )
 */
float Model::GetDiDj( int i, int j, float time )
{
  int n = PoolNum();
  float **d = matrix( 1, n+1, 1, n );    // the dij array 
  float t = time;
  assert( d != NULL );
  float ret;
  int k, r;
  float sum_k, sum_r;

  // compute the entire dij array first
  GetDij( d, 1, n+1, 1, n );

  sum_k = 0.0F;
  for( k = 1; k <= n; k++ )
  {
    sum_r = 0.0F;
    for( r = k+1; r <= n; r++ )
      sum_r += (float) ((d[j][r]*d[i][k] + d[i][r]*d[j][k])*
                        (exp(L(k)*t) - exp(L(r)*t))/
                        (L(k) - L(r)));
    sum_k += (float) (d[i][k]*d[j][k]*t*exp(L(k)*t) + sum_r);
  }
  ret = sum_k;

  // clean up
  free_matrix( d, 1, n+1, 1, n );

  return( ret );
}



/*
 * Get the dij array, these are temporary values used Lindell's formula
 *
 * d      - the array
 * You have to read the paper to understand this!
 */
void Model::GetDij( float **d, int nrl, int nrh, int ncl, int nch ) 
{
  int n = PoolNum();
  assert( nrl == 1 && nrh == n+1 && ncl ==1 && nch == n );
  int i, j;

  // recursive formula for computing dij
  for( j = 1; j <= n; j++ )
  {
    d[1][j] = 0.0F;
    d[2][j] = A(j);
    for( i = 3; i <= n+1; i++ )
      d[i][j] = ((L(j) - K(i-2,i-2))*d[i-1][j] - d[i-2][j])/Gamma(i-1);
  }
}




/*
 * Compute the variance array
 */
int Model::UpdateVar()
{
  int i;

  if( SampleN() > 0 )
  {
    switch( ErrMod() )
    {
      case CCV:
        for( i = 1; i <= SampleN(); i++ )
          Var(i) = square(Ccv()*SumOfExp( SampleTime(i) )/100);
        break;
      case CSD:
        for( i = 1; i <= SampleN(); i++ )
          Var(i) = square(Csd());
        break;
      case CVAR:
        for( i = 1; i <= SampleN(); i++ )
          Var(i) = Cvar();
        break;
      case VVAR:
        for( i = 1; i <= SampleN(); i++ )
          Var(i) = Varb() + Varc() * (float) pow( (double) SumOfExp( SampleTime(i) ), (double) Vard() );
        break;
      case VCV:
        for( i = 1; i <= SampleN(); i++ )
          Var(i) = square(Vcv(i)*SumOfExp( SampleTime(i) )/100);
        break;
      default:
        assert( 0 );
        break;
    }
  }

//  for( i = 1; i <= SampleN(); i++ )
//  	TRACE( "Var(%d) = %g\n", i, Var(i) );

  return( 0 );
}

                   
                   
// member access functions 
ModelType&  Model::MType()
{
  return( Type );
}

int&    Model::PoolNum()
{
  return( N );
}  

float&  Model::SSConc1()
{
  return( C1 );
}

OutputType& Model::OutputUnit()
{
  return( outunit );
}

float& Model::BodyWeight()
{
  return( bw );
}

float& Model::Dose()
{
  return( dose );
}

float&  Model::A( int i )
{
  assert( i >= 1 && i <= N );
  return( a_exp[i-1] );
}                           

float&  Model::L( int i )
{
  assert( i >= 1 && i <= N );
  return( lambda[i-1] );
}
    
float&  Model::Alpha( int i )
{
  assert( i >= 1 && i <= N );
  return( alpha[i-1] );
}

float&  Model::Beta( int i )
{
  assert( i >= 1 && i <= N );
  return( beta[i-1] );
}      

float&  Model::Gamma( int i )
{                  
  assert( i >= 1 && i <= N );
  return( gamma[i-1] );
}


float&  Model::K( int i, int j, BoundType bt )
{      
  assert( i >= 0 && i <= N && j >= 0 && j <= N ); 
  assert( bt == Exact || bt == Lower || bt == Upper );
  
  if( i == j && bt == Exact )
    return( kappa[i-1] );  
  else if( bt == Lower )  
    return( lowerK[ IJtoI(i,j) ] );
  else 
    return( upperK[ IJtoI(i,j) ] );
}      


float&  Model::Constr( int i, int j )
{
  return( constraint[ IJtoI(i,j) ] );
}         


float&  Model::ConK( int i, int j, BoundType bt )
{      
  assert( i >= 0 && i <= N && j >= 0 && j <= N && i != j ); 
  assert( bt == Lower || bt == Upper );
  
  if( bt == Lower )  
    return( lowerC[ IJtoI(i,j) ] );
  else 
    return( upperC[ IJtoI(i,j) ] );
}            

// V(N+1) is the total volumn bounds
float&  Model::V( int i, BoundType bt )
{
  assert( i >= 1 && i <= N+1 );
  assert( bt == Upper || bt == Lower );

  if( bt == Upper )
    return( upperV[i-1] );
  else 
    return( lowerV[i-1] );
}                                   

float&  Model::ConV( int i, BoundType bt )
{             
  assert( i >= 1 && i <= N+1 );
  assert( bt == Upper || bt == Lower );
  
  if( bt == Upper )
    return( upperCV[i-1] );
  else 
    return( lowerCV[i-1] );
}

// Q(N+1) is the total volumn bounds
float&  Model::Q( int i, BoundType bt )
{
  assert( i >= 1 && i <= N+1 );
  assert( bt == Upper || bt == Lower );
  
  if( bt == Upper )
    return( upperQ[i-1] );
  else 
    return( lowerQ[i-1] );
}                                   

float&  Model::ConQ( int i, BoundType bt )
{    
  assert( i >= 1 && i <= N+1 );
  assert( bt == Upper || bt == Lower );
  
  if( bt == Upper )
    return( upperCQ[i-1] );
  else 
    return( lowerCQ[i-1] );
}

float&  Model::Mflux( int i, int j, BoundType bt )
{      
  assert( i >= 0 && i <= N && j >= 0 && j <= N && i != j );
  assert( bt == Lower || bt == Upper );
  
  if( bt == Lower )  
    return( lowerMflux[ IJtoI(i,j) ] );
  else 
    return( upperMflux[ IJtoI(i,j) ] );
}            

float&  Model::ConMflux( int i, int j, BoundType bt )
{      
  assert( i >= 0 && i <= N && j >= 0 && j <= N && i != j );
  assert( bt == Lower || bt == Upper );
  
  if( bt == Lower )  
    return( lowerCMflux[ IJtoI(i,j) ] );
  else 
    return( upperCMflux[ IJtoI(i,j) ] );
}            

float& Model::PCR()
{
  return( pcr );
}

float& Model::PR()
{
  return( pr );
}

float& Model::MRT( BoundType bt )
{
  assert( bt == Upper || bt == Lower );

  if( bt == Upper )
    return( upperMRT );
  else 
    return( lowerMRT );
}

float& Model::VD( BoundType bt )
{
  assert( bt == Upper || bt == Lower );

  if( bt == Upper )
    return( upperVD );
  else 
    return( lowerVD );
}

int& Model::SampleN()
{
  return( timen );
}

ErrModType& Model::ErrMod()
{
  return( errmod );
}

float& Model::SampleTime( int i )
{
  assert( i >= 1 && i <= timen );
  return( time[i-1] );
}

float& Model::Ccv()
{
  return( ccv );
}

float& Model::Csd()
{
  return( csd );
}

float& Model::Cvar()
{
  return( cvar );
}

float& Model::Varb()
{
  return( varb );
}

float& Model::Varc()
{
  return( varc );
}

float& Model::Vard()
{
  return( vard );
}

float& Model::Vcv( int i )
{
  assert( i >= 1 && i <= timen );
  return( vcv[i-1] );
}

float& Model::Var( int i )
{
  assert( i >= 1 && i <= timen );
  return( var[i-1] );
}

FMatrix Model::Cov()
{
  return( cov );
}

FMatrix Model::ConCov()
{
  return( concov );
}

FMatrix Model::Corr()
{
  return( corr );
}


TimeUnitType& Model::TimeUnit()
{
	return( timeunit );
}

VolUnitType&		Model::VolUnit()
{
	return( volunit );
}

MassUnitType& Model::MassUnit()
{
	return( massunit );
}

EndoMUType& Model::EndoMU()
{
	return( endomu );
}

char* Model::GetTimeStr()	
{
	return( timestr[(int) timeunit] );
}

char* Model::GetVolStr()	
{
	return( volstr[(int) volunit] );
}

char* Model::GetMassStr()
{
	return( massstr[(int) massunit] );
}

char* Model::GetEndoMUStr()
{
	return( endomustr[(int) endomu] );
}

// misc

void Model::SortAL()
{
  // Sort A and L so that L is in increasing order
  // insertion sort
  
  int   i,j;
  float  Atmp, Ltmp;
  
  for( j = 2; j <= N; j++ )
  {
    Atmp = A(j);
    Ltmp = L(j);
    i = j - 1;
    while( i > 0 && L(i) > Ltmp )
    {
      L(i+1) = L(i);
      A(i+1) = A(i);
      i--;
    }    
    L(i+1) = Ltmp;
    A(i+1) = Atmp;
  }
}        

// Maps ordinary indexes to array index
int    Model::IJtoI( int i, int j )          
{            
  if( i < 0 || j < 0 || i > N || j > N )
    crash( "IJtoI:  index out of bounds" );

  if( Type == Mam )
  {
    if( i == 0 )
      return( j-1 );        // map k0i to i-1
    else if( i == 1 )
      return( j-2+N );    // map k1i to i-2+N
    else if( j == 1 )
      return( i+2*N-3 );    // map ki1 to i-2 + 2N-1
    else                                  
    {
      crash( "invalid Mam index in IJtoI" );  // should never occur
      return( -1 );              // invalid
    }
  }
  else if( Type == Cat )
  {                         
    if( i == 0 )
      return( j-1 );      // map k0i to i-1
    else if( i == j-1 )
      return( j-2+N );    // same map as above
    else if( i-1 == j )
      return( i+2*N-3 );    // same map as above
    else
    {
      crash( "invalid Cat index in IJtoI" );  // should never occur
      return( -1 );              // invalid
    }
  }
  else                        
  {
    crash( "invalid model type in IJtoI" );    // should never occur
    return( -1 );
  }
}                                                       

/*
 * map array index to ordinary index
 * array is zero based
 */
void   Model::ItoIJ( int k, int* one, int* two )
{       
  int    i,j;

  if( Type == Mam )
  {       
    if( k > 3*N-3 )
      crash( "invalid index number in ItoIJ" );  // should never occur
    else if( k <= N-1 )
    {       
      i = 0;
      j = k+1;
    }
    else if( k <= 2*N-2 )
    {    
      i = 1;
      j = k-N+2;
    }
    else
    {             
      i = k-(2*N-1)+2;
      j = 1;
    }   
    *one = i;
    *two = j;
  }
  else if( Type == Cat )
  {    
     if( k > 3*N-3 )
      crash( "invalid index number in ItoIJ" );  // should never occur
    else if( k <= N-1 )
    {       
      i = 0;
      j = k+1;
    }
    else if( k <= 2*N-2 )
    {
      j = k-N+2;
      i = j-1;
    }
    else
    {             
      i = k-(2*N-1)+2;
      j = i-1;
    }   
    *one = i;
    *two = j;
  }
  else
    crash( "invalid model type in ItoIJ" );    // should never occur
}

/*void  Model::Serialize(CArchive& ar)
{
#define FBS 256
  CStdioFile  *pFile = (CStdioFile *) ar.GetFile();
  char  tmp[FBS];
  int   i, j;

  if(  ar.IsStoring() )
  {
    // Model Type
    if( MType() == Mam )
      sprintf( tmp, "Mammillary Model\n" );
    else
      sprintf( tmp, "Catenary Model\n" );
    pFile->WriteString( tmp );

    // Number of pools
    sprintf( tmp, "Number of Pools = %d\n", PoolNum() );
    pFile->WriteString( tmp );

    // Output Type 
    if( OutputUnit() == Conc )
      sprintf( tmp, "Concentration Output\n", SSConc1() );
    else
      sprintf( tmp, "Mass Output\n", SSConc1() );
    pFile->WriteString( tmp );

    // Concentration in pool 1
    sprintf( tmp, "Steady State Concentration in Pool1 = %g\n", SSConc1() );
    pFile->WriteString( tmp );

    // Bodyweight
		sprintf( tmp, "Body weight = %g\n", BodyWeight() );
    pFile->WriteString( tmp );

    // Units
    sprintf( tmp, "Time Unit = %d\n", (int) TimeUnit() );
    pFile->WriteString( tmp );
    sprintf( tmp, "Volume Unit = %d\n", (int) VolUnit() );
    pFile->WriteString( tmp );
    sprintf( tmp, "Mass Unit = %d\n", (int) MassUnit() );
    pFile->WriteString( tmp );
    sprintf( tmp, "Endogenous Mass Unit = %d\n", (int) EndoMU() );
    pFile->WriteString( tmp );

    // Coefficient & Exponent
    for( i = 1; i <= PoolNum(); i++ )
    {
      sprintf( tmp, "A(%d) = %g\tL(%d) = %g\n", i, A(i), i, L(i) );
      pFile->WriteString( tmp );
    }

    // Alpha and Beta
    for( i = 1; i <= PoolNum(); i++ )
    {
      sprintf( tmp, "Alpha(%d) = %g\tBeta(%d) = %g\n", i, Alpha(i), i, Beta(i) );
      pFile->WriteString( tmp );
    }

    // Kii
    for( i = 1; i <= PoolNum(); i++ )
    {
      sprintf( tmp, "K%d,%d = %g\n", i, i, K( i, i ) );
      pFile->WriteString( tmp );
    }

    // Gamma
    for( i = 2; i <= PoolNum(); i++ )
    {
      sprintf( tmp, "K%d,%d K%d,%d = %g\n", i-1, i, i, i-1, Gamma(i) );
      pFile->WriteString( tmp );
    }

    // Kij
    for( i = 0; i < 3*PoolNum() - 2; i++ )
    {
      int one, two;
      ItoIJ( i, &one, &two );
      sprintf( tmp, "%g < %g < K%d,%d < %g < %g\n", 
               K( one, two, Lower ), ConK( one, two, Lower ),
               one, two,
               ConK( one, two, Upper ), K( one, two, Upper ) );
      pFile->WriteString( tmp );
    }

    // Constraints
    for( i = 0; i < 3*PoolNum() - 2; i++ )
    {
      int one, two;
      ItoIJ( i, &one, &two );
      sprintf( tmp, "Constraint on K%d,%d = %g\n", 
               one, two, Constr( one, two ) );
      pFile->WriteString( tmp );
    }

    // V
    for( i = 1; i <= PoolNum(); i++ )
    {
      sprintf( tmp, "%g < %g < V%d < %g < %g\n", 
               V( i, Lower ), ConV( i, Lower ),
               i,
               ConV( i, Upper ), V( i, Upper ) );
      pFile->WriteString( tmp );
    }

    // Q
    for( i = 1; i <= PoolNum(); i++ )
    {
      sprintf( tmp, "%g < %g < Q%d < %g < %g\n", 
               Q( i, Lower ), ConQ( i, Lower ),
               i,
               ConQ( i, Upper ), Q( i, Upper ) );
      pFile->WriteString( tmp );
    }

    // Mass Flux
    for( i = 0; i < 3*PoolNum() - 2; i++ )
    {
      int one, two;
      ItoIJ( i, &one, &two );
      sprintf( tmp, "%g < %g < K%d,%d Q%d < %g < %g\n", 
               Mflux( one, two, Lower ), ConMflux( one, two, Lower ),
               one, two, two,
               ConMflux( one, two, Upper ), Mflux( one, two, Upper ) );
      pFile->WriteString( tmp );
    }

    // Plasma Clearance Rate
    sprintf( tmp, "Plasma Clearance Rate = %g\n", PCR() );
    pFile->WriteString( tmp );

    // Plasma Production Rate
    sprintf( tmp, "Plasma Production Rate = %g\n", PR() );
    pFile->WriteString( tmp );

    // Mean Residence Time
    sprintf( tmp, "%g < MRT < %g\n", MRT( Lower ), MRT( Upper ) );
    pFile->WriteString( tmp );

    // VD
    sprintf( tmp, "%g < VD < %g\n", VD( Lower ), VD( Upper ) );
    pFile->WriteString( tmp );

    // Number of samples 
    sprintf( tmp, "Number of Samples = %d\n", SampleN() ); 
    pFile->WriteString( tmp );

    // Sampling time
    for( i = 1; i <= SampleN(); i++ )
    {
      sprintf( tmp, "T(%d) = %g\n", i, SampleTime( i ) );
      pFile->WriteString( tmp );
    }

    // Error Model
    sprintf( tmp, "Error Model Code %d\n", (int) ErrMod() );
    pFile->WriteString( tmp );
    switch( ErrMod() )
    {
      case CCV:
        sprintf( tmp, "Constant CV = %g\n", Ccv() ); 
        pFile->WriteString( tmp );
        break;
      case CSD:
        sprintf( tmp, "Constant SD = %g\n", Csd() ); 
        pFile->WriteString( tmp );
        break;
      case CVAR:
        sprintf( tmp, "Constant VAR = %g\n", Cvar() ); 
        pFile->WriteString( tmp );
        break;
      case VVAR:
        sprintf( tmp, "B = %g\n", Varb() ); 
        pFile->WriteString( tmp );
        sprintf( tmp, "C = %g\n", Varc() ); 
        pFile->WriteString( tmp );
        sprintf( tmp, "D = %g\n", Vard() ); 
        pFile->WriteString( tmp );
        break;
      case VCV:
        for( i = 1; i <= SampleN(); i++ )
        {
          sprintf( tmp, "CV(%d) = %g\n", i, Vcv( i ) );
          pFile->WriteString( tmp );
        }
        break;
      default:
        assert( 0 );
        break;
    }

    // Correlation Matrix
    if( SampleN() > 0 )
			for( i = 1; i <= corr.XDim(); i++ )
			{
				CString row;

				for( j = 1; j <= corr.YDim(); j++ )
				{
					sprintf( tmp, "%g ", corr.rElem1( j, i ) );
					row += tmp;
				}

				row += "\n";
				pFile->WriteString( row );
			}

    // Covariance Matrix
    if( SampleN() > 0 )
			for( i = 1; i <= cov.XDim(); i++ )
			{
				CString row;

				for( j = 1; j <= cov.YDim(); j++ )
				{
					sprintf( tmp, "%g ", cov.rElem1( j, i ) );
					row += tmp;
				}

				row += "\n";
				pFile->WriteString( row );
			}

		// Dose
		sprintf( tmp, "Dose = %g\n", Dose() );
    pFile->WriteString( tmp );
  }
  else
  {
    // Model Type
    pFile->ReadString( tmp, FBS );
    if( strstr( tmp, "Mam" ) )
      MType() = Mam;
    else
      MType() = Cat;

    // Number of pools
    pFile->ReadString( tmp, FBS );
    sscanf( tmp, "Number of Pools = %d\n", &PoolNum() );

    // Output Type 
    pFile->ReadString( tmp, FBS );
    if( strstr( tmp, "Conc" ) ) 
    	OutputUnit() = Conc;
    else
    	OutputUnit() = Mass;

    // Concentration in pool 1
    pFile->ReadString( tmp, FBS );
    sscanf( tmp, "Steady State Concentration in Pool1 = %g\n", &SSConc1() );

    // Bodyweight
    pFile->ReadString( tmp, FBS );
		sscanf( tmp, "Body weight = %g\n", &BodyWeight() );

    // Units
    int ucode;
    pFile->ReadString( tmp, FBS );
    sscanf( tmp, "Time Unit = %d\n", &ucode );
    TimeUnit() = (TimeUnitType) ucode;
    pFile->ReadString( tmp, FBS );
    sscanf( tmp, "Volume Unit = %d\n", &ucode );
    VolUnit() = (VolUnitType) ucode;
    pFile->ReadString( tmp, FBS );
    sscanf( tmp, "Mass Unit = %d\n", &ucode );
    MassUnit() = (MassUnitType) ucode;
    pFile->ReadString( tmp, FBS );
    sscanf( tmp, "Endogenous Mass Unit = %d\n", &ucode );
    EndoMU() = (EndoMUType) ucode;

    // Coefficient & Exponent
    for( i = 1; i <= PoolNum(); i++ )
    {
      pFile->ReadString( tmp, FBS );
      sscanf( tmp, "A(%d) = %g\tL(%d) = %g\n", &j, &A(i), &j, &L(i) );
    }

    // Alpha and Beta
    for( i = 1; i <= PoolNum(); i++ )
    {
      pFile->ReadString( tmp, FBS );
      sscanf( tmp, "Alpha(%d) = %g\tBeta(%d) = %g\n", &j, &Alpha(i), &j, &Beta(i) );
    }

    // Kii
    for( i = 1; i <= PoolNum(); i++ )
    {
      pFile->ReadString( tmp, FBS );
      sscanf( tmp, "K%d,%d = %g\n", &j, &j, &K( i, i ) );
    }

    // Gamma
    for( i = 2; i <= PoolNum(); i++ )
    {
      pFile->ReadString( tmp, FBS );
      sscanf( tmp, "K%d,%d K%d,%d = %g\n", &j, &j, &j, &j, &Gamma(i) );
    }

    // Kij
    for( i = 0; i < 3*PoolNum() - 2; i++ )
    {
      int one, two;
      ItoIJ( i, &one, &two );
      pFile->ReadString( tmp, FBS );
      sscanf( tmp, "%g < %g < K%d,%d < %g < %g\n", 
               &K( one, two, Lower ), &ConK( one, two, Lower ),
               &j, &j, 
               &ConK( one, two, Upper ), &K( one, two, Upper ) );
    }

    // Constraints
    for( i = 0; i < 3*PoolNum() - 2; i++ )
    {
      int one, two;
      ItoIJ( i, &one, &two );
      pFile->ReadString( tmp, FBS );
      sscanf( tmp, "Constraint on K%d,%d = %g\n", 
               &j, &j, &Constr( one, two ) );
    }

    // V
    for( i = 1; i <= PoolNum(); i++ )
    {
      pFile->ReadString( tmp, FBS );
      sscanf( tmp, "%g < %g < V%d < %g < %g\n", 
               &V( i, Lower ), &ConV( i, Lower ),
               &j,
               &ConV( i, Upper ), &V( i, Upper ) );
    }

    // Q
    for( i = 1; i <= PoolNum(); i++ )
    {
      pFile->ReadString( tmp, FBS );
      sscanf( tmp, "%g < %g < Q%d < %g < %g\n", 
               &Q( i, Lower ), &ConQ( i, Lower ),
               &j,
               &ConQ( i, Upper ), &Q( i, Upper ) );
    }

    // Mass Flux
    for( i = 0; i < 3*PoolNum() - 2; i++ )
    {
      int one, two;
      ItoIJ( i, &one, &two );
      pFile->ReadString( tmp, FBS );
      sscanf( tmp, "%g < %g < K%d,%d Q%d < %g < %g\n", 
               &Mflux( one, two, Lower ), &ConMflux( one, two, Lower ),
               &j, &j, &j,
               &ConMflux( one, two, Upper ), &Mflux( one, two, Upper ) );
    }

    // Plasma Clearance Rate
    pFile->ReadString( tmp, FBS );
    sscanf( tmp, "Plasma Clearance Rate = %g\n", &PCR() );

    // Plasma Production Rate
    pFile->ReadString( tmp, FBS );
    sscanf( tmp, "Plasma Production Rate = %g\n", &PR() );

    // Mean Residence Time
    pFile->ReadString( tmp, FBS );
    sscanf( tmp, "%g < MRT < %g\n", &MRT( Lower ), &MRT( Upper ) );

    // VD
    pFile->ReadString( tmp, FBS );
    sscanf( tmp, "%g < VD < %g\n", &VD( Lower ), &VD( Upper ) );

    // Number of samples 
    pFile->ReadString( tmp, FBS );
    sscanf( tmp, "Number of Samples = %d\n", &SampleN() ); 

    // Sampling time
    for( i = 1; i <= SampleN(); i++ )
    {
      int foo;
      pFile->ReadString( tmp, FBS );
      sscanf( tmp, "T(%d) = %g\n", &foo, &SampleTime( i ) );
    }

    // Error Model
    pFile->ReadString( tmp, FBS );
    sscanf( tmp, "Error Model Code %d\n", &ErrMod() );
    switch( ErrMod() )
    {
      case CCV:
        pFile->ReadString( tmp, FBS );
        sscanf( tmp, "Constant CV = %g\n", &Ccv() ); 
        break;
      case CSD:
        pFile->ReadString( tmp, FBS );
        sscanf( tmp, "Constant SD = %g\n", &Csd() ); 
        break;
      case CVAR:
        pFile->ReadString( tmp, FBS );
        sscanf( tmp, "Constant VAR = %g\n", &Cvar() ); 
        break;
      case VVAR:
        pFile->ReadString( tmp, FBS );
        sscanf( tmp, "B = %g\n", &Varb() ); 
        pFile->ReadString( tmp, FBS );
        sscanf( tmp, "C = %g\n", &Varc() ); 
        pFile->ReadString( tmp, FBS );
        sscanf( tmp, "D = %g\n", &Vard() ); 
        break;
      case VCV:
        for( i = 1; i <= SampleN(); i++ )
        {
          int foo;
          pFile->ReadString( tmp, FBS );
          sscanf( tmp, "CV(%d) = %g\n", &foo, &Vcv( i ) );
        }
        break;
      default:
        assert( 0 );
        break;
    }

    // Correlation Matrix
    if( SampleN() > 0 )
    {
			corr.SetDim( 5*PoolNum()-4, 5*PoolNum()-4 );
			for( i = 1; i <= 5*PoolNum()-4; i++ )
			{
				pFile->ReadString( tmp, FBS );
				char *p = tmp;

				for( j = 1; j <= 5*PoolNum()-4; j++ )
				{
					sscanf( p, "%g ", &(corr.Elem1( j, i )) );
					p = strstr( p, " " ) + 1;
				}
			}
    }

    // Covariance Matrix
    if( SampleN() > 0 )
    {
			cov.SetDim( 5*PoolNum()-4, 5*PoolNum()-4 );
			for( i = 1; i <= 5*PoolNum()-4; i++ )
			{
				pFile->ReadString( tmp, FBS );
				char *p = tmp;

				for( j = 1; j <= 5*PoolNum()-4; j++ )
				{
					sscanf( p, "%g ", &(cov.Elem1( j, i )) );
					p = strstr( p, " " ) + 1;
				}
			}
    }

		// Dose
    pFile->ReadString( tmp, FBS );
		sscanf( tmp, "Dose = %g\n", &Dose() );
  }
#undef FBS
} */

// Checks to see if the constraints are within bounds
// This should work for either mamillary model or catenary mode
int  Model::CheckConstr()
{    
  for( int i = 0; i < 3*N-2; i++ )
    if( constraint[i] != UNCONSTRAINED )
      if( constraint[i] < lowerK[i] || 
        constraint[i] > upperK[i] )
      {
#if  DEBUG
        cerr << i << " " << lowerK[i] << " ";
        cerr << constraint[i] << " " << upperK[i] << endl;
#endif 
        return( -1 );
      }
  return( 0 );
}      


// Adjust the other bounds 
void  Model::MamFix()
{
  int    j,m;
  float  sum;
  int    count = 0;    // count the number of unconstrained parameters
  int    which;      // the latest constrained pool number
  
  for( j = 2; j <= N; j++ )
    if( Constr(j,1) == UNCONSTRAINED &&
      Constr(1,j) == UNCONSTRAINED &&
      Constr(0,j) == UNCONSTRAINED )  // do not modify the already constrained ones
    {        
      count++;        // some book keeping
      which = j;
      
      sum = 0.0F;
      for( m = 2; m <= N; m++ )
        if( m != j )
          sum += ConK(m,1,Lower);
      sum += ConK(0,1,Lower);
      ConK(j,1,Upper) = -K(1,1) - sum;
      ConK(1,j,Lower) = Gamma(j)/ConK(j,1,Upper);
      ConK(0,j,Upper) = -K(j,j) - ConK(1,j,Lower);
    }                    
  
  if( Constr(0,1) == UNCONSTRAINED )
  {                                                  
    count++;        // some book keeping 
    which = 1;

    ConK(0,1,Upper) = -K(1,1);
    for( j = 2; j <= N; j++ )
      ConK(0,1,Upper) -= ConK(j,1,Lower);
  }                                                
  
  if( count == 1 )      
  {
    // This means that there are n-1 constraints
    // So the only unconstrained pool should also
    // be identifible
    
    if( which != 1 )
    {                    
      ConK(which,1,Lower) = ConK(which,1,Upper);
      ConK(1,which,Upper) = ConK(1,which,Lower);
      ConK(0,which,Lower) = ConK(0,which,Upper);
    }
    else
    {           
      ConK(0,1,Lower) = ConK(0,1,Upper);
    }
  }
}

float  Model::Corr_rElem( int iy, int ix ) const { 
    return corr.rElem(iy, ix); 
}; 

float  Model::Cov_rElem( int iy, int ix ) const { 
    return cov.rElem(iy, ix); 
}; 

float  Model::ConCov_rElem( int iy, int ix ) const { 
    return concov.rElem(iy, ix); 
};

float  Model::Corr_rElem1( int iy, int ix ) const { 
    return corr.rElem1(iy, ix); 
}; 

float  Model::Cov_rElem1( int iy, int ix ) const { 
    return cov.rElem1(iy, ix); 
}; 

float  Model::ConCov_rElem1( int iy, int ix ) const { 
    return concov.rElem1(iy, ix); 
};

// read only access to xdim
int Model::Corr_XDim() const { 
	
	return corr.XDim(); 
};  

// read only access to ydim
int Model::Corr_YDim() const { 
	
	return corr.YDim(); 
};  


void Model::Reset() {

	int i;
	
	for (i=0; i < MAXPOOL; i++) {
		a_exp[i] = 0;      // A's and L's in the exponential
		lambda[i] = 0;    // equation  
        alpha[i] = 0;      // alphas and betas in the 
		beta[i] = 0;      // transfer function
		kappa[i] = 0;      // identifible parameter combinations
		gamma[i] = 0;              
	}

	for (i=0; i < (MAXPOOL+1); i++) {
		upperV[i] = 0;    // unconstrained volumn, use N+1 to store total      
		lowerV[i] = 0;           
		upperCV[i] = 0;    // constrained volumn, use N+1 to store total
		lowerCV[i] = 0;               
		upperQ[i] = 0;    // unconstrained mass, use N+1 to store total
		lowerQ[i] = 0;                    
		upperCQ[i] = 0;    // constrained mass, use N+1 to store total
		lowerCQ[i] = 0;                    
	}

	for (i=0; i < MAXPARAM; i++) {
		constraint[i] = -1;// -1 = unconstrained, others = constrained to that value
		upperK[i] = 0;    // unidentifible parameters
		lowerK[i] = 0;             
		upperC[i] = 0;    // constrained unidentifible parameters
		lowerC[i] = 0;                                                    
		upperMflux[i] = 0;    // unconstrained mass flux
		lowerMflux[i] = 0;
		upperCMflux[i] = 0;    // constrained mass flux
		lowerCMflux[i] = 0;    
	}

	for (i=0; i < MAXNTIME; i++) {
		time[i] =0;      // sample time
		vcv[i] = 0;      // variable variances
		var[i] = 0;      // variances
	}
	  // input variables
	Type = Mam;              // Mam or Cat
	N = 0;                  // number of pools 
	C1 = 0;                  // steady state value in pool1
	outunit = Conc;          // steady state unit type
	bw = 0;                  // bodyweight
	dose= 0;								// tracer dose
	errmod = CCV;    // Error Model 
    timen = 0;              // number of samples
    ccv = 0;                // constant cv
    csd = 0;                // constant sd
    cvar = 0;                // constant var
    varb = 0;                // variance formula var = b + c * z^d
    varc = 0;
    vard = 0;
  
  // output variables
	pcr = 0;  // pool 1 plasma clearnace rate
	pr = 0;    // pool 1 production rate
	upperMRT = 0;    // mean residence time
    lowerMRT = 0;
	upperVD = 0;    // equivalent distribution volume
    lowerVD = 0;

}

