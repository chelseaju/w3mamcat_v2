// Programmed By Hsiao-Te Su
// 12/26/93  
#ifndef INC_MAMCATMODEL_H
#define INC_MAMCATMODEL_H

#include "matrix.h"
#include "model.h"

//	Model p;
 

class MamcatModel                       
{  
//private:
  
public:
	MamcatModel() {};            // default constructor
	~MamcatModel() {};            // default destructor 
  
  bool Check();             // check internal consistencies
  int Calculate();         // performs internal calculations

public:
  // map ordinary index to array index
	int		IJtoI( int i, int j );

  // map array index to ordinary index
	void	ItoIJ( int k, int* one, int* two );

public:
  // access ftn for model type
	ModelType  MType();     
	void	MType( ModelType x );

  // access ftn for pool number
	int     PoolNum();      
	void	PoolNum( int x );
  
  // access ftn for steady state value in pool1
	float   SSConc1();
	void	SSConc1( float x );

  // access ftn for steady state units
	OutputType OutputUnit();
	void    OutputUnit( OutputType x );

  // access ftn for bodyweight value in pool1
	float   BodyWeight();
	void    BodyWeight( float x );

  // access ftn for dose  
	float	Dose();					 
	void	Dose( float x );

  // access ftn for a_exp
	float   A( int i );
	void	A( int i, float x );

  // access ftn for lambda
	float   L( int i ); 
	void	L( int i, float x );

  // access ftn for alpha
	float   Alpha( int i );
	void    Alpha( int i, float x );

  // access ftn for beta 
	float   Beta( int i);
	void    Beta( int i, float x );

  // access ftn for gamma
	float   Gamma( int i );
	void	Gamma( int i, float x );

  // access ftn for all of the K's
  float		K( int i, int k, BoundType bt = Exact);
  void		K( int i, int k, float x, BoundType bt = Exact);

  // access ftn for constraint
  float		Constr( int i, int k );
  void		Constr( int i, int k , float x );

  // access ftn for constrained unidentifible parameters
  float	ConK( int i, int k, BoundType bt = Exact );
  void	ConK( int i, int k, float x, BoundType bt = Exact );

  // access ftn for unconstrained volumn
  float  V( int i, BoundType bt = Exact );
  void  V( int i, float x, BoundType bt = Exact );

  // access ftn for constrained volumn
  float  ConV( int i, BoundType bt = Exact );
  void  ConV( int i, float x, BoundType bt = Exact );

  // access tfn for unconstrained mass
  float Q( int i, BoundType bt = Exact );
  void  Q( int i, float x, BoundType bt = Exact );

  // access ftn for constrained mass
  float  ConQ( int i, BoundType bt = Exact );
  void ConQ( int i, float x, BoundType bt = Exact );

  // access ftn for unconstrained mass flux
  float  Mflux( int i, int j, BoundType bt = Exact );
  void  Mflux( int i, int j, float x, BoundType bt = Exact );

  // access ftn for constrained mass flux
  float  ConMflux( int i, int j, BoundType bt = Exact );
  void ConMflux( int i, int j, float x, BoundType bt = Exact );

  // access ftn for plasma clearance rate
  float		PCR();
  void		PCR( float x );

  // access ftn for plasma production rate
  float		PR();
  void		PR( float x );

  // access ftn for mean residence time
  float		MRT( BoundType bt = Exact );
  void		MRT( float x, BoundType bt = Exact);

  // access ftn for equivalent distribution volumn
  float		VD( BoundType bt = Exact );
  void		VD( float x, BoundType bt = Exact);

  // access ftn for number of samples
  int		SampleN();
  void		SampleN( int x );

  // access ftn for sample times
  float		SampleTime( int i );
  void		SampleTime( int i, float x);

  // access ftn for error model
  ErrModType ErrMod();
  void		ErrMod( ErrModType x );

  // constant cv
  float		Ccv();
  void		Ccv( float x );

  // constant sd
  float		Csd();
  void		Csd( float x );

  // constant var
  float		Cvar();
  void		Cvar( float x );

  // variance formula var = b + c * z^d
  float		Varb();
  void		Varb( float x );

  float		Varc();
  void		Varc( float x );

  float		Vard();
  void		Vard( float x );

  // access ftn for variable CV
  float		Vcv( int i );
  void		Vcv ( int i, float x );

  // access ftn for variances
  float		Var( int i );
  void		Var( int i, float x);

  // read-only access to the covariance matrix
  //FMatrix Cov();               
  // read-only access to the correlation matrix
  //FMatrix Corr();          
  
  // access ftn for time unit
  TimeUnitType TimeUnit();
  void		TimeUnit( TimeUnitType x );

  // access ftn for volumn/weight unit
  VolUnitType VolUnit();
  void		VolUnit( VolUnitType x );

  // access ftn for mass unit
  MassUnitType MassUnit();
  void MassUnit( MassUnitType x );

  // access ftn for endogenous mass unit
  EndoMUType EndoMU();
  void EndoMU( EndoMUType x );

  // access ftn for time unit string
  char *GetTimeStr();
  
  // access ftn for volume unit string
  char *GetVolStr();								

  // access ftn for mass unit string
  char *GetMassStr();
  
  // access ftn for endogenous mass unit string
  char *GetEndoMUStr();
    
  // members
public:
	// read only access to the (ix,iy)th correlation element
    float  Corr_rElem( int iy, int ix ) const ;

	// read only access to the (ix,iy)th covariance element
    float  Cov_rElem( int iy, int ix ) const ;

	// read only access to the (ix,iy)th constrained covariance element
    float  ConCov_rElem( int iy, int ix ) const ;


	// same as above, except offset by 1, ie., first element at (1,1)
	// read only access to the (ix,iy)th element
    float  Corr_rElem1( int iy, int ix ) const;

	// read only access to the (ix,iy)th element
    float  Cov_rElem1( int iy, int ix ) const ;

	// read only access to the (ix,iy)th constrained covariance element
    float  ConCov_rElem1( int iy, int ix ) const ;

	// read only access to xdim
    int	Corr_XDim() const ;

	// read only access to ydim
    int	Corr_YDim() const ;

	// resets all private model values
	void Reset(); 

};
#endif // INC_MAMCATMODEL_H
