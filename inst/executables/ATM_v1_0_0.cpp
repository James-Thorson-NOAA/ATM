#include <TMB.hpp>

// Function to subtract diagonal by colSum as mass-balance
template<class Type>
//Eigen::SparseMatrix<Type> subtract_colsum( int n_g, Eigen::SparseMatrix<Type> mat_gg ){
matrix<Type> subtract_colsum( int n_g, matrix<Type> mat_gg ){

  // New matrix
  //Eigen::SparseMatrix<Type> newmat_gg;
  matrix<Type> newmat_gg( n_g, n_g );

  // Get columnwise sum
  vector<Type> colsum( n_g );
  colsum = mat_gg.colwise().sum();

  // Loop through diagonal
  for( int g=0; g<n_g; g++ ){
    //newmat_gg.coeffRef(g,g) = -1.0 * colsum(g);
    newmat_gg(g,g) = -1.0 * colsum(g);
  }

  return newmat_gg;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Global namespaces
  using namespace Eigen;

  // Data
  DATA_ARRAY( X_guyk );
  DATA_IMATRIX( uy_tz );
  DATA_IMATRIX( satellite_iz );
  DATA_MATRIX( survey_jz );
  DATA_VECTOR( duration_u );
  //DATA_SPARSE_MATRIX( A_gg );
  DATA_MATRIX( A_gg );

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
  Type jnll = 0;
  //array<Type> prob_satellite_igt( n_i, n_g, n_t );
  array<Type> logprob_satellite_igt( n_i, n_g, n_t );
  //Eigen::SparseMatrix<Type> Diffusion_gg;
  //Eigen::SparseMatrix<Type> Taxis_gg;
  //Eigen::SparseMatrix<Type> Mprime_gg;
  //Eigen::SparseMatrix<Type> I_gg( n_g, n_g );
  matrix<Type> Diffusion_gg( n_g, n_g );
  matrix<Type> Taxis_gg( n_g, n_g );
  matrix<Type> Mprime_gg( n_g, n_g );
  matrix<Type> Movement_gg( n_g, n_g );
  matrix<Type> I_gg( n_g, n_g );
  I_gg.setIdentity();
  vector<Type> Preference_g( n_g );
  vector<Type> init_g( n_g );

  // Loop through times
  for( int t=0; t<n_t; t++ ){
    // Diffusion-rate matrix
    Diffusion_gg = sigma2 * A_gg;
    Diffusion_gg = subtract_colsum( n_g, Diffusion_gg );

    // Advection-rate matrix
    //Preference_g = as.vector( X_guyk[,data_list$uy_tz[tI,'u'],data_list$uy_tz[tI,'y'],] %*% alpha )
    Preference_g.setZero();
    for( int g=0; g<n_g; g++ ){
    for( int u=0; u<n_u; u++ ){
    for( int y=0; y<n_y; y++ ){
    for( int k=0; k<n_k; k++ ){
      Preference_g(g) += X_guyk(g,u,y,k) * alpha_k(k);
    }}}}
    for( int g1=0; g1<n_g; g1++ ){
    for( int g2=0; g2<n_g; g2++ ){
      Taxis_gg(g1,g2) = A_gg(g1,g2) * Preference_g(g1) * Preference_g(g2);
    }}
    Taxis_gg = subtract_colsum( n_g, Taxis_gg );

    // Movement probability matrix
    Mprime_gg = Diffusion_gg + Taxis_gg;
    Movement_gg = expm( Mprime_gg );

    // Apply to satellite tags
    // TODO:  Explore logspace_sum for numerical stability
    for( int i=0; i<n_i; i++ ){
      if( satellite_iz(i,2) == t ){
        for( int g=0; g<n_g; g++ ){
          // prob_satellite_igt(i,g,t) = Movement_gg( g, satellite_iz(i,0) );
          logprob_satellite_igt(i,g,t) = log(Movement_gg( g, satellite_iz(i,0) ));
        }
      }
      if( (satellite_iz(i,2)<t) & (satellite_iz(i,3)>=t) ){
        for( int g1=0; g1<n_g; g1++ ){
        for( int g2=0; g2<n_g; g2++ ){
          //prob_satellite_igt(i,g2,t) += Movement_gg(g2,g1) * prob_satellite_igt(i,g1,t-1);
          logprob_satellite_igt(i,g2,t) = logspace_add( logprob_satellite_igt(i,g2,t), log(Movement_gg(g2,g1))+logprob_satellite_igt(i,g1,t-1) );
        }}
      }
    }
  }

  // Calculate log-likelihood for sat-tags
  for( int i=0; i<n_i; i++ ){
    //jnll -= log(prob_satellite_igt( i, satellite_iz(i,1), satellite_iz(i,3) ));
    jnll -= logprob_satellite_igt( i, satellite_iz(i,1), satellite_iz(i,3) );
  }

  REPORT( sigma2 );
  REPORT( alpha_k );
  REPORT( Diffusion_gg );
  REPORT( Taxis_gg );
  REPORT( Mprime_gg );
  REPORT( Movement_gg );
  REPORT( logprob_satellite_igt );
  REPORT( jnll );

  return jnll;
}
