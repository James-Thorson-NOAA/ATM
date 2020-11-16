#include <TMB.hpp>

// Function to subtract diagonal by colSum as mass-balance
template<class Type>
matrix<Type> subtract_colsum_from_diagonal( int n_g, matrix<Type> mat_gg ){

  // Get columnwise sum
  vector<Type> colsum( n_g );
  colsum = mat_gg.colwise().sum();

  // Loop through diagonal
  for( int g=0; g<n_g; g++ ){
    mat_gg(g,g) = -1.0 * colsum(g);
  }

  return mat_gg;
}

// Duplicate with Eigen::SparseMatrix
template<class Type>
Eigen::SparseMatrix<Type> subtract_colsum_from_diagonal( int n_g, Eigen::SparseMatrix<Type> mat_gg, matrix<int> mat_ij ){

  // Get columnwise sum
  Eigen::SparseMatrix<Type> newmat_gg( n_g, n_g );
  newmat_gg = mat_gg;

  // Loop through diagonal
  for( int r=0; r<mat_ij.rows(); r++ ){
    newmat_gg.coeffRef( mat_ij(r,1), mat_ij(r,1) ) -= mat_gg.coeffRef( mat_ij(r,0), mat_ij(r,1) );
  }

  return newmat_gg;
}

// Function to calculate matrix exponential, or approximate it using Euler method
template<class Type>
matrix<Type> matexp( int n_g, int log2steps, matrix<Type> mat_gg ){

  if( log2steps > 100 ){
    return expm(mat_gg);
  }else{
    matrix<Type> I_gg( n_g, n_g );
    I_gg.setIdentity();
    mat_gg = I_gg + mat_gg / pow(2, log2steps);
    for( int step=0; step<log2steps; step++ ){
      mat_gg = mat_gg * mat_gg;
    }
    return( mat_gg );
  }
}

// Duplicate with Eigen::SparseMatrix
template<class Type>
Eigen::SparseMatrix<Type> matexp( int n_g, int log2steps, Eigen::SparseMatrix<Type> mat_gg ){

  //if( log2steps > 100 ){
  //  return expm(mat_gg);
  //}else{
    Eigen::SparseMatrix<Type> I_gg( n_g, n_g );
    I_gg.setIdentity();
    mat_gg = I_gg + mat_gg / pow(2, log2steps);
    for( int step=0; step<log2steps; step++ ){
      mat_gg = mat_gg * mat_gg;
    }
    return( mat_gg );
  //}
}

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Global namespaces
  using namespace Eigen;

  // Data
  DATA_INTEGER( log2steps );

  DATA_ARRAY( X_guyk );
  DATA_IMATRIX( uy_tz );
  DATA_IMATRIX( satellite_iz );
  DATA_MATRIX( survey_jz );
  DATA_VECTOR( duration_u );
  DATA_SPARSE_MATRIX( A_gg );
  DATA_IMATRIX( A_ij );

  // Parameters
  PARAMETER( ln_sigma );
  PARAMETER_VECTOR( alpha_logit_ratio_k );
  
  // Dimensions
  int n_g = X_guyk.col(0).col(0).col(0).size();
  int n_u = X_guyk.col(0).col(0).size() / n_g;
  int n_y = X_guyk.col(0).size() / n_g / n_u;
  int n_k = X_guyk.size() / n_g / n_u / n_y;
  int n_t = uy_tz.rows();
  int n_i = satellite_iz.rows();

  // Transform inputs
  Type sigma2 = exp( 2.0 * ln_sigma );
  vector<Type> alpha_k( n_k );
  alpha_k = sigma2 * (2.0 * invlogit(alpha_logit_ratio_k) - 1.0);

  // Global variables
  vector<Type> nll_i( n_i );
  Type jnll = 0;
  array<Type> prob_satellite_igt( n_i, n_g, n_t );
  //array<Type> logprob_satellite_igt( n_i, n_g, n_t );
  vector<Type> Preference_g( n_g );
  vector<Type> init_g( n_g );

  Eigen::SparseMatrix<Type> Diffusion_gg( n_g, n_g );
  Eigen::SparseMatrix<Type> Taxis_gg( n_g, n_g );
  Eigen::SparseMatrix<Type> Mprime_gg( n_g, n_g );
  Eigen::SparseMatrix<Type> Mprimesum_gg( n_g, n_g );
  Eigen::SparseMatrix<Type> Movement_gg( n_g, n_g );
  Mprimesum_gg.setZero();

  // Loop through times
  for( int t=0; t<n_t; t++ ){
    // Diffusion-rate matrix
    Diffusion_gg = sigma2 * A_gg;
    Diffusion_gg = subtract_colsum_from_diagonal( n_g, Diffusion_gg, A_ij );

    // Advection-rate matrix
    Preference_g.setZero();
    for( int g=0; g<n_g; g++ ){
    for( int k=0; k<n_k; k++ ){
      Preference_g(g) += X_guyk(g,uy_tz(t,0),uy_tz(t,1),k) * alpha_k(k);
    }}
    for( int r=0; r<A_ij.rows(); r++ ){
      Taxis_gg.coeffRef(A_ij(r,0),A_ij(r,1)) = Preference_g(A_ij(r,0)) - Preference_g(A_ij(r,1));
    }
    Taxis_gg = subtract_colsum_from_diagonal( n_g, Taxis_gg, A_ij );

    // Movement probability matrix
    Mprime_gg = Diffusion_gg + Taxis_gg;
    Mprimesum_gg += Mprime_gg;
    Movement_gg = matexp( n_g, log2steps, Mprime_gg );

    // Apply to satellite tags
    // TODO:  Explore logspace_add for numerical stability
    for( int i=0; i<n_i; i++ ){
      if( satellite_iz(i,2) == t ){
        for( int g=0; g<n_g; g++ ){
          prob_satellite_igt(i,g,t) = Movement_gg.coeffRef( g, satellite_iz(i,0) );
          //logprob_satellite_igt(i,g,t) = log(Movement_gg( g, satellite_iz(i,0) ));
        }
      }
      if( (satellite_iz(i,2)<t) & (satellite_iz(i,3)>=t) ){
        for( int g1=0; g1<n_g; g1++ ){
        for( int g2=0; g2<n_g; g2++ ){
          prob_satellite_igt(i,g2,t) += Movement_gg.coeffRef(g2,g1) * prob_satellite_igt(i,g1,t-1);
          //logprob_satellite_igt(i,g2,t) = logspace_add( logprob_satellite_igt(i,g2,t), log(Movement_gg(g2,g1))+logprob_satellite_igt(i,g1,t-1) );
        }}
      }
    }
  }

  // Calculate log-likelihood for sat-tags
  for( int i=0; i<n_i; i++ ){
    nll_i(i) = -1 * log(prob_satellite_igt( i, satellite_iz(i,1), satellite_iz(i,3) ));
    //nll_i(i) = -1 * logprob_satellite_igt( i, satellite_iz(i,1), satellite_iz(i,3) );
  }
  jnll = sum(nll_i);

  REPORT( sigma2 );
  REPORT( alpha_k );
  //REPORT( Diffusion_gg );
  //REPORT( Taxis_gg );
  //REPORT( Mprime_gg );
  //REPORT( Movement_gg );
  //REPORT( Mprimesum_gg );
  //REPORT( logprob_satellite_igt );
  //REPORT( prob_satellite_igt );
  //REPORT( jnll );
  //REPORT( Preference_g );
  //REPORT( nll_i );

  return jnll;
}
