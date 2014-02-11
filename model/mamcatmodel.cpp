// Programmed By Hsiao-Te Su
// 12/26/93  

#include "matrix.h"
#include "model.h"
#include "mamcatmodel.h"

Model p;

bool MamcatModel::Check() { return p.Check(); }          // check internal consistencies
int MamcatModel::Calculate() { return p.Calculate(); }        // performs internal calculations


// map ordinary index to array index
int		MamcatModel::IJtoI( int i, int j ) { return p.IJtoI(i, j); }

  // map array index to ordinary index
void	MamcatModel::ItoIJ( int k, int* one, int* two ) { p.ItoIJ(k, one, two); }


  // access ftn for model type
ModelType  MamcatModel::MType() { return (p.MType()); }
void	MamcatModel::MType( ModelType x ) { p.MType() = x;}

  // access ftn for pool number
int     MamcatModel::PoolNum() { return p.PoolNum(); }     
void	MamcatModel::PoolNum( int x ) { p.PoolNum() = x; }
  
 // access ftn for steady state value in pool1
float   MamcatModel::SSConc1() { return p.SSConc1(); }
void	MamcatModel::SSConc1( float x ) { p.SSConc1() = x; }

  // access ftn for steady state units
OutputType MamcatModel::OutputUnit() { return p.OutputUnit(); }
void    MamcatModel::OutputUnit( OutputType x ) { p.OutputUnit() = x; }

  // access ftn for bodyweight value in pool1
float   MamcatModel::BodyWeight() { return p.BodyWeight(); }
void    MamcatModel::BodyWeight( float x ) { p.BodyWeight() = x; }

  // access ftn for dose  
float	MamcatModel::Dose() { return p.Dose(); }		 
void	MamcatModel::Dose( float x ) { p.Dose() = x; }

  // access ftn for a_exp
float   MamcatModel::A( int i ) { return p.A(i); }
void	MamcatModel::A( int i, float x ) { p.A(i) = x; }

  // access ftn for lambda
float   MamcatModel::L( int i ) { return p.L(i); }  
void	MamcatModel::L( int i, float x ) { p.L(i) = x; }

  // access ftn for alpha
float   MamcatModel::Alpha( int i ) { return p.Alpha(i); }
void    MamcatModel::Alpha( int i, float x ) { p.Alpha(i) = x; }

  // access ftn for beta 
float   MamcatModel::Beta( int i) { return p.Beta(i); }
void    MamcatModel::Beta( int i, float x ) { p.Beta(i) = x; }
  // access ftn for gamma

float   MamcatModel::Gamma( int i ) { return p.Gamma(i); }
void	MamcatModel::Gamma( int i, float x ) { p.Gamma(i) = x; }

  // access ftn for all of the K's
float	MamcatModel::K( int i, int k, BoundType bt ) { return p.K(i, k, bt); }
void	MamcatModel::K( int i, int k, float x, BoundType bt ) { p.K(i, k, bt) = x; }

  // access ftn for constraint
float	MamcatModel::Constr( int i, int k ) { return p.Constr(i,k); }
void	MamcatModel::Constr( int i, int k , float x ) { p.Constr(i,k) = x; }

  // access ftn for constrained unidentifible parameters
float	MamcatModel::ConK( int i, int k, BoundType bt  ) { return p.ConK(i, k, bt); }
void	MamcatModel::ConK( int i, int k, float x, BoundType bt ) { p.ConK(i, k, bt) = x; }

  // access ftn for unconstrained volumn
float  MamcatModel::V( int i, BoundType bt ) { return p.V(i, bt); }
void  MamcatModel::V( int i, float x, BoundType bt ) { p.V(i, bt) = x; }

  // access ftn for constrained volumn
float  MamcatModel::ConV( int i, BoundType bt  ) { return p.ConV(i, bt); }
void  MamcatModel::ConV( int i, float x, BoundType bt ) { p.ConV(i, bt) = x; }

  // access tfn for unconstrained mass
float MamcatModel::Q( int i, BoundType bt ) { return p.Q(i, bt); }
void  MamcatModel::Q( int i, float x, BoundType bt ) { p.Q(i, bt) = x; }

  // access ftn for constrained mass
float  MamcatModel::ConQ( int i, BoundType bt  ) { return p.ConQ(i, bt); }
void MamcatModel::ConQ( int i, float x, BoundType bt ) { p.ConQ(i, bt) = x; }

  // access ftn for unconstrained mass flux
float  MamcatModel::Mflux( int i, int j, BoundType bt  ) { return p.Mflux(i, j, bt); }
void  MamcatModel::Mflux( int i, int j, float x, BoundType bt  ) { p.Mflux(i, j, bt) = x; }

  // access ftn for constrained mass flux
float  MamcatModel::ConMflux( int i, int j, BoundType bt  ) { return p.ConMflux(i, j, bt); }
void MamcatModel::ConMflux( int i, int j, float x, BoundType bt  ) { p.ConMflux(i, j, bt) = x; }

  // access ftn for plasma clearance rate
float	MamcatModel::PCR() { return p.PCR(); }
void	MamcatModel::PCR( float x ) { p.PCR() = x; }

  // access ftn for plasma production rate
float	MamcatModel::PR() { return p.PR(); }
void	MamcatModel::PR( float x ) { p.PR() = x; }

  // access ftn for mean residence time
float	MamcatModel::MRT( BoundType bt ) { return p.MRT(bt); }
void	MamcatModel::MRT( float x, BoundType bt ) { p.MRT(bt) = x; }

  // access ftn for equivalent distribution volumn
float	MamcatModel::VD( BoundType bt ) { return p.VD(bt); }
void	MamcatModel::VD( float x, BoundType bt ) { p.VD(bt) = x; }

  // access ftn for number of samples
int		MamcatModel::SampleN() { return p.SampleN(); }
void	MamcatModel::SampleN( int x ) { p.SampleN() = x; }

  // access ftn for sample times
  float		MamcatModel::SampleTime( int i ) { return p.SampleTime(i); }
  void		MamcatModel::SampleTime( int i, float x) { p.SampleTime(i) = x; }

  // access ftn for error model
  ErrModType MamcatModel::ErrMod() { return p.ErrMod(); }
  void		MamcatModel::ErrMod( ErrModType x ) { p.ErrMod() = x; }

  // constant cv
  float		MamcatModel::Ccv() { return p.Ccv(); }
  void		MamcatModel::Ccv( float x ) { p.Ccv() = x; }

  // constant sd
  float		MamcatModel::Csd() { return p.Csd(); }
  void		MamcatModel::Csd( float x ) { p.Csd() = x; }

  // constant var
  float		MamcatModel::Cvar() { return p.Cvar(); };
  void		MamcatModel::Cvar( float x ) { p.Cvar() = x; }

  // variance formula var = b + c * z^d
  float		MamcatModel::Varb() { return p.Varb(); }
  void		MamcatModel::Varb( float x ) { p.Varb() = x; };

  float		MamcatModel::Varc() { return p.Varc(); }
  void		MamcatModel::Varc( float x ) { p.Varc() = x; }

  float		MamcatModel::Vard() { return p.Vard(); };
  void		MamcatModel::Vard( float x ) { p.Vard() = x; }

  // access ftn for variable CV
  float		MamcatModel::Vcv( int i ) { return p.Vcv(i); }
  void		MamcatModel::Vcv ( int i, float x ) { p.Vcv(i) = x; }

  // access ftn for variances
  float		MamcatModel::Var( int i ) { return p.Var(i); }
  void		MamcatModel::Var( int i, float x) { p.Var(i) = x; }

  // read-only access to the covariance matrix
  //FMatrix Cov();               
  // read-only access to the correlation matrix
  //FMatrix Corr();          
  
  // access ftn for time unit
  TimeUnitType MamcatModel::TimeUnit() { return p.TimeUnit(); }
  void		MamcatModel::TimeUnit( TimeUnitType x ) { p.TimeUnit() = x; }

  // access ftn for volumn/weight unit
  VolUnitType MamcatModel::VolUnit() { return p.VolUnit(); };
  void		MamcatModel::VolUnit( VolUnitType x ) { p.VolUnit() = x; }

  // access ftn for mass unit
  MassUnitType MamcatModel::MassUnit() { return p.MassUnit(); }
  void MamcatModel::MassUnit( MassUnitType x ) { p.MassUnit() = x; }

  // access ftn for endogenous mass unit
  EndoMUType MamcatModel::EndoMU() { return p.EndoMU(); }
  void MamcatModel::EndoMU( EndoMUType x ) { p.EndoMU() = x; }

  // access ftn for time unit string
  char * MamcatModel::GetTimeStr() { return p.GetTimeStr(); }
  
  // access ftn for volume unit string
  char * MamcatModel::GetVolStr() { return p.GetVolStr(); }							

  // access ftn for mass unit string
  char * MamcatModel::GetMassStr() { return p.GetMassStr(); }
  
  // access ftn for endogenous mass unit string
  char * MamcatModel::GetEndoMUStr() { return p.GetEndoMUStr(); }
    
  // members
// public memebers
	// read only access to the (ix,iy)th correlation element
    float  MamcatModel::Corr_rElem( int iy, int ix ) const { 
		return p.Corr_rElem(iy, ix);
	}

	// read only access to the (ix,iy)th covariance element
    float  MamcatModel::Cov_rElem( int iy, int ix ) const {
		return p.Cov_rElem(iy, ix);
	}

	// read only access to the (ix,iy)th constrained covariance element
    float  MamcatModel::ConCov_rElem( int iy, int ix ) const {
		return p.ConCov_rElem(iy, ix);
	}

	// same as above, except offset by 1, ie., first element at (1,1)
	// read only access to the (ix,iy)th element
    float  MamcatModel::Corr_rElem1( int iy, int ix ) const {
		return p.Corr_rElem1(iy, ix);
	}

	// same as above, except offset by 1, ie., first element at (1,1)
	// read only access to the (ix,iy)th element
    float  MamcatModel::Cov_rElem1( int iy, int ix ) const {
		return p.Cov_rElem1(iy, ix);
	}

	// same as above, except offset by 1, ie., first element at (1,1)
	// read only access to the (ix,iy)th constrained covariance element
    float  MamcatModel::ConCov_rElem1( int iy, int ix ) const {
		return p.ConCov_rElem1(iy, ix);
	}

	// read only access to xdim
    int	MamcatModel::Corr_XDim() const { return p.Corr_XDim(); }

	// read only access to ydim
    int	MamcatModel::Corr_YDim() const { return p.Corr_YDim(); }

	void MamcatModel::Reset() { 
		p.Reset();
	} 

