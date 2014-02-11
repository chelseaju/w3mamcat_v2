// Programmed by Hsiao-Te su
// 1/14/94  
//
// define the keys used for file storage
#ifndef IIA_KEY_H
#define IIA_KEY_H
     
     // input variables      
#define  title      "iia file, do NOT modify!"
#define beginModel    "Model: "
#define beginPoolNum  "Number of Pools: "
#define beginA      "Begin Coefficients"
#define endA      "End Coefficients"
#define beginL      "Begin Exponents"
#define endL      "End Exponents"
#define beginC1      "Begin steady state concentration in pool 1"
#define endC1      "End steady state concentration in pool 1"
#define beginConstr    "Begin Constraints"
#define endConstr    "End Constraints"

  // output variables
#define beginAlpha    "Begin Alpha"
#define endAlpha    "End Alpha"
#define beginBeta    "Begin Beta"
#define endBeta      "End Beta"
#define beginKappa    "Begin Kappa"
#define endKappa    "End Kappa"
#define beginGamma    "Begin Gamma"
#define endGamma    "End Gamma"

#define beginUK0i    "Begin K0i"
#define endUK0i      "End Koi"
#define beginUK1i    "Begin K1i"
#define endUK1i      "End K1i"
#define beginUKi1    "Begin Ki1"
#define endUKi1      "End Ki1"              

#define beginCK0i    "Begin Constrained K0i"
#define endCK0i      "End Constrained K0i"
#define beginCK1i    "Begin Constrained K1i"
#define endCK1i      "End Constrained K1i"
#define beginCKi1    "Begin Constrained Ki1"
#define endCKi1      "End Constrained Ki1"

#define beginUKii1    "Begin Ki,i+1"
#define endUKii1    "End Ki,i+1"
#define beginUKi1i    "Begin Ki+1,i"
#define endUKi1i    "End Ki+1,i"

#define beginCKii1    "Begin Constrained Ki,i+1"
#define endCKii1    "End Constrained Ki,i+1"
#define beginCKi1i    "Begin Constrained Ki+1,i"
#define endCKi1i    "End Constrained Ki+1,i"

#define beginUVi      "Begin Vi"
#define endUVi        "End Vi"
#define beginUQi      "Begin Qi"
#define endUQi        "End Qi"

#define beginCVi      "Begin Constrained Vi"
#define endCVi        "End Constrained Vi"
#define beginCQi      "Begin Constrained Qi"
#define endCQi        "End Constrained Qi"

#define beginUMF0i      "Begin Mass Flux 0i"
#define endUMF0i      "End Mass Flux oi"
#define beginUMF1i      "Begin Mass Flux 1i"
#define endUMF1i      "End Mass Flux 1i"
#define beginUMFi1      "Begin Mass Flux i1"
#define endUMFi1      "End Mass Flux i1"              

#define beginCMF0i      "Begin Constrained Mass Flux 0i"
#define endCMF0i      "End Constrained Mass Flux 0i"
#define beginCMF1i      "Begin Constrained Mass Flux 1i"
#define endCMF1i      "End Constrained Mass Flux 1i"
#define beginCMFi1      "Begin Constrained Mass Flux i1"
#define endCMFi1      "End Constrained Mass Flux i1"

#define beginUMFii1      "Begin Mass Flux i,i+1"
#define endUMFii1      "End Mass Flux i,i+1"
#define beginUMFi1i      "Begin Mass Flux i+1,i"
#define endUMFi1i      "End Mass Flux i+1,i"

#define beginCMFii1      "Begin Constrained Mass Flux i,i+1"
#define endCMFii1      "End Constrained Mass Flux i,i+1"
#define beginCMFi1i      "Begin Constrained Mass Flux i+1,i"
#define endCMFi1i      "End Constrained Mass Flux i+1,i"

#endif  // IIA_KEY_H
