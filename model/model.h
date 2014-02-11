// Programmed By Hsiao-Te Su
// 12/26/93  
#ifndef INC_MODEL_H
#define INC_MODEL_H

#include "matrix.h"

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define  UNCONSTRAINED  -1.0F
#define MAXPOOL    10
#define MAXPARAM  3*MAXPOOL-2
#define MAXNTIME    100
                             
typedef enum {Mam, Cat} ModelType;           
typedef enum {Exact, Lower, Upper} BoundType;
typedef enum {Mass, Conc} OutputType;
typedef enum {Gram, Kg} BWType;
typedef enum {CCV, CSD, CVAR, VVAR, VCV} ErrModType;
typedef enum {T_SEC, T_MIN, T_HOUR, T_DAY} TimeUnitType;
typedef enum {V_ML, V_L, V_UL} VolUnitType;
typedef enum {M_PDOSE, M_G, M_KG, M_LB} MassUnitType;
typedef enum {ENDO_M_G, ENDO_M_KG, ENDO_M_LB} EndoMUType;
           
class Model                       
{        
public:
  Model();            // default constructor
  ~Model();            // default destructor 
  
  bool Check();            // check internal consistencies
  int Calculate();        // performs internal calculations

  // helping functions
public:             

private:  
  // calculation 
  float  SumOfExp( float x );      // evaluates the sum of exponential ftn
  float  Transfer( float x );    // evaluates the transfer ftn
  float  FindRoot( float xl, float x2, float tol );  // finds the root of the transfer ftn
  int   ALtoAB();        // converts A and L to alpha and beta 
  int   ABtoKG();        // alpha & beta to identifible parameters
  int    UnID();          // compute the max & min bounds of the Un-ID's
  int    ConID();        // compute the max & min bounds of the constrained IDs
  int    VQID();          // compute the bounds of compartment volumn
  int    MfluxID();      // compute the bounds of mass fluxes
  int    Global();        // compute some global parameters
  int    GetCorr();        // compute the correlation matrix
  
  // misc
  float square( float x ) { return( x*x ); };  
  int    delta( int i, int j ) { return( (i==j?1:0) ); };  // discrete delta function
  void   SortAL();
  int    CheckConstr();  // checks to see if the constraints are in bounds
  void  MamFix();        // fix up other parameters when some constraints are applied
  int    UpdateVar();    // compute the variances
  FMatrix GetCovId();    // compute the convariance matrix of the identifiable paramteres
  FMatrix GetM( int i );  // compute Mi
  FMatrix GetdFdP(int);      // compute the sensitive functions
  FMatrix GetdYdP( float time );  // compute the dYdP row vector
  float    GetdYdPij( float time, int paramI );  // compute an element of dYdP
  float   GetDiDj( int i, int j, float time );  // returns Lap^-1( D_i*D_j/D_1^2 )
  void    GetDij( float **d, int ncl, int nch, int nrl, int nrh );  // get the di array

public:
  int    IJtoI( int i, int j );  // map ordinary index to array index
  void  ItoIJ( int k, int* one, int* two );     // map array index to ordinary index
  //void  Serialize(CArchive& ar); // file i/o


  // member access functions    
  // Be careful!  These access functions return references. 
public:
  ModelType&  MType();         // access ftn for model type
  int&    PoolNum();           // access ftn for pool number
  float&  SSConc1();           // access ftn for steady state value in pool1
  OutputType& OutputUnit();        // access ftn for steady state units
  float&  BodyWeight();        // access ftn for bodyweight value in pool1
  float&	Dose();							 // access ftn for dose  
  float&  A( int i );          // access ftn for a_exp
  float&  L( int i );          // access ftn for lambda
  float&  Alpha( int i );      // access ftn for alpha
  float&  Beta( int i );       // access ftn for beta 
  float&  Gamma( int i );      // access ftn for gamma
  float&  K( int i, int k,         // access ftn for all of the K's
           BoundType bt = Exact);       
  float&  Constr( int i, int k );  // access ftn for constraint
  float&  ConK( int i, int k,      // access ftn for constrained unidentifible parameters
          BoundType bt = Exact );          
  float&  V( int i, BoundType bt = Exact );  // access ftn for unconstrained volumn
  float&  ConV( int i, BoundType bt = Exact );  // access ftn for constrained volumn
  float&  Q( int i, BoundType bt = Exact );     // access tfn for unconstrained mass
  float&  ConQ( int i, BoundType bt = Exact );  // access ftn for constrained mass
  float&  Mflux( int i, int j, BoundType bt = Exact );// access ftn for unconstrained mass flux
  float&  ConMflux( int i, int j, BoundType bt = Exact );// access ftn for constrained mass flux
  float&  PCR();                         // access ftn for plasma clearance rate
  float&  PR();                          // access ftn for plasma production rate
  float&  MRT( BoundType bt = Exact );   // access ftn for mean residence time
  float&  VD( BoundType bt = Exact );    // access ftn for equivalent distribution volumn
  int&    SampleN();                    // access ftn for number of samples
  float&  SampleTime( int i );          // access ftn for sample times
  ErrModType& ErrMod();                // access ftn for error model
  float&  Ccv();                        // constant cv
  float&  Csd();                        // constant sd
  float&  Cvar();                      // constant var
  float&  Varb();                      // variance formula var = b + c * z^d
  float&  Varc();
  float&  Vard();
  float&  Vcv( int i );                // access ftn for variable CV
  float&  Var( int i );                // access ftn for variances
  FMatrix Cov();                      // read-only access to the covariance matrix
  FMatrix ConCov();                    // read-only access to the constrained covariance matrix
  FMatrix Corr();                     // read-only access to the correlation matrix
  TimeUnitType& TimeUnit();						// access ftn for time unit
  VolUnitType&		VolUnit();					// access ftn for volumn/weight unit
  MassUnitType& MassUnit();						// access ftn for mass unit
  EndoMUType& EndoMU();								// access ftn for endogenous mass unit
  char *GetTimeStr();									// access ftn for time unit string
  char *GetVolStr();									// access ftn for volume unit string
  char *GetMassStr();									// access ftn for mass unit string
  char *GetEndoMUStr();								// access ftn for endogenous mass unit string
    
  // members
private:
  // input variables
  ModelType  Type;              // Mam or Cat
  int      N;                  // number of pools 
  float    C1;                  // steady state value in pool1
  OutputType  outunit;          // steady state unit type
  float    bw;                  // bodyweight
  float		 dose;								// tracer dose
  float    a_exp[MAXPOOL];      // A's and L's in the exponential
  float    lambda[MAXPOOL];    // equation    
  float    constraint[MAXPARAM];// -1 = unconstrained, others = constrained to that value
  ErrModType  errmod;    // Error Model 
  int      timen;              // number of samples
  float    time[MAXNTIME];      // sample time
  float    ccv;                // constant cv
  float    csd;                // constant sd
  float    cvar;                // constant var
  float    varb;                // variance formula var = b + c * z^d
  float    varc;
  float    vard;
  float    vcv[MAXNTIME];      // variable variances
  float    var[MAXNTIME];      // variances
  
  // output variables
  float    alpha[MAXPOOL];      // alphas and betas in the 
  float    beta[MAXPOOL];      // transfer function
  float    kappa[MAXPOOL];      // identifible parameter combinations
  float    gamma[MAXPOOL];              
  float    upperK[MAXPARAM];    // unidentifible parameters
  float    lowerK[MAXPARAM];             
  float    upperC[MAXPARAM];    // constrained unidentifible parameters
  float    lowerC[MAXPARAM];                                                    
  float    upperV[MAXPOOL+1];    // unconstrained volumn, use N+1 to store total      
  float    lowerV[MAXPOOL+1];           
  float    upperCV[MAXPOOL+1];    // constrained volumn, use N+1 to store total
  float    lowerCV[MAXPOOL+1];               
  float    upperQ[MAXPOOL+1];    // unconstrained mass, use N+1 to store total
  float    lowerQ[MAXPOOL+1];                    
  float    upperCQ[MAXPOOL+1];    // constrained mass, use N+1 to store total
  float    lowerCQ[MAXPOOL+1];                    
  float    upperMflux[MAXPARAM];    // unconstrained mass flux
  float    lowerMflux[MAXPARAM];
  float    upperCMflux[MAXPARAM];    // constrained mass flux
  float    lowerCMflux[MAXPARAM];    
  float    pcr;  // pool 1 plasma clearnace rate
  float    pr;    // pool 1 production rate
  float    upperMRT,    // mean residence time
           lowerMRT;
  float    upperVD,    // equivalent distribution volume
           lowerVD;
  FMatrix cov;
  FMatrix concov;
  FMatrix corr;      // the correlation matrix
  TimeUnitType timeunit;
  VolUnitType	volunit;
  MassUnitType massunit;
  EndoMUType endomu;
	static char *timestr[];
	static char *volstr[];
	static char *massstr[];
	static char *endomustr[];
public:
    float  Corr_rElem( int iy, int ix ) const; // read only access to the (ix,iy)th correlation element
    float  Cov_rElem( int iy, int ix ) const; // read only access to the (ix,iy)th covariance element
	float  ConCov_rElem( int iy, int ix ) const; // read only access to the (ix,iy)th constrained covariance element
	// same as above, except offset by 1, ie., first element at (1,1)
    float  Corr_rElem1( int iy, int ix ) const; // read only access to the (ix,iy)th element
    float  Cov_rElem1( int iy, int ix ) const; // read only access to the (ix,iy)th element
	float  ConCov_rElem1( int iy, int ix ) const; // read only access to the (ix,iy)th constrained element
    int	Corr_XDim() const;  // read only access to xdim
    int	Corr_YDim() const;  // read only access to ydim
	void Reset(); // resets all the private values of the model

};
#endif // INC_MODEL_H
