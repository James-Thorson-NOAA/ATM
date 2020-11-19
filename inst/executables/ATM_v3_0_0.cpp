#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

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

// Function to calculate matrix exponential, or approximate it using Euler method
template<class Type>
matrix<Type> matrix_exponential( int n_g, int log2steps, matrix<Type> mat_gg ){

  if( (log2steps <=0 ) || (log2steps > 100) ){
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

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Global namespaces
  using namespace R_inla;
  using namespace Eigen;
  using namespace density;

  // Data
  DATA_INTEGER( log2steps );
  DATA_ARRAY( X_guyk );
  DATA_IMATRIX( uy_tz );
  DATA_IMATRIX( satellite_iz );
  DATA_VECTOR( duration_u );
  DATA_MATRIX( A_gg );

  // Survey data objects
  DATA_STRUCT( spde_aniso, spde_aniso_t );
  DATA_VECTOR( b_j );
  DATA_IVECTOR( t_j );
  DATA_IVECTOR( g_j );

  // Parameters
  PARAMETER( ln_sigma );
  PARAMETER_VECTOR( alpha_logit_ratio_k );

  // SPDE parameters
  PARAMETER_VECTOR(ln_H_input); // Anisotropy parameters
  PARAMETER(ln_kappa);
  PARAMETER(ln_sigma_epsilon0);
  PARAMETER(ln_sigma_epsilon);
  PARAMETER(ln_phi);
  PARAMETER(power_prime);
  PARAMETER_VECTOR(Beta_t);
  PARAMETER_ARRAY(ln_d_st);     // Annual log-density

  // Dimensions
  int n_g = X_guyk.col(0).col(0).col(0).size();
  int n_u = X_guyk.col(0).col(0).size() / n_g;
  int n_y = X_guyk.col(0).size() / n_g / n_u;
  int n_k = X_guyk.size() / n_g / n_u / n_y;
  int n_t = uy_tz.rows();
  int n_i = satellite_iz.rows();
  int n_j = b_j.size();
  int n_s = ln_d_st.col(0).size();

  // Transform inputs
  Type sigma2 = exp( 2.0 * ln_sigma );
  vector<Type> alpha_k( n_k );
  alpha_k = sigma2 * (2.0 * invlogit(alpha_logit_ratio_k) - 1.0);

  // Global variables
  Type jnll = 0;
  vector<Type> nll_i( n_i );
  vector<Type> nll_j( n_j );
  vector<Type> nll_t( n_t );

  // Global variables
  array<Type> prob_satellite_igt( n_i, n_g, n_t );
  //array<Type> logprob_satellite_igt( n_i, n_g, n_t );
  matrix<Type> Diffusion_gg( n_g, n_g );
  matrix<Type> Taxis_gg( n_g, n_g );
  matrix<Type> Mprime_gg( n_g, n_g );
  array<Type> Mprime_ggt( n_g, n_g, n_t );
  matrix<Type> Mprimesum_gg( n_g, n_g );
  matrix<Type> Movement_gg( n_g, n_g );
  matrix<Type> Preference_gt( n_g, n_t );
  vector<Type> init_g( n_g );
  matrix<Type> Msum_gg( n_g, n_g );
  Mprimesum_gg.setZero();
  Preference_gt.setZero();

  // Loop through times
  for( int t=0; t<n_t; t++ ){
    // Diffusion-rate matrix
    Diffusion_gg = sigma2 * A_gg;
    Diffusion_gg = subtract_colsum_from_diagonal( n_g, Diffusion_gg );

    // Advection-rate matrix
    for( int g=0; g<n_g; g++ ){
    for( int k=0; k<n_k; k++ ){
      Preference_gt(g,t) += X_guyk(g,uy_tz(t,0),uy_tz(t,1),k) * alpha_k(k);
    }}
    for( int g1=0; g1<n_g; g1++ ){
    for( int g2=0; g2<n_g; g2++ ){
      Taxis_gg(g1,g2) = A_gg(g1,g2) * (Preference_gt(g1,t) - Preference_gt(g2,t));
    }}
    Taxis_gg = subtract_colsum_from_diagonal( n_g, Taxis_gg );

    // Movement probability matrix
    Mprime_gg = Diffusion_gg + Taxis_gg;
    Mprimesum_gg += Mprime_gg;
    Movement_gg = matrix_exponential( n_g, log2steps, Mprime_gg );
    for( int g1=0; g1<n_g; g1++ ){
    for( int g2=0; g2<n_g; g2++ ){
      Mprime_ggt(g1,g2,t) = Mprime_gg(g1,g2);
    }}

    // Apply to satellite tags
    // TODO:  Explore logspace_sum for numerical stability
    for( int i=0; i<n_i; i++ ){
      if( satellite_iz(i,2) == t ){
        for( int g=0; g<n_g; g++ ){
          prob_satellite_igt(i,g,t) = Movement_gg( g, satellite_iz(i,0) );
          //logprob_satellite_igt(i,g,t) = log(Movement_gg( g, satellite_iz(i,0) ));
        }
      }
      if( (satellite_iz(i,2)<t) & (satellite_iz(i,3)>=t) ){
        for( int g1=0; g1<n_g; g1++ ){
        for( int g2=0; g2<n_g; g2++ ){
          prob_satellite_igt(i,g2,t) += Movement_gg(g2,g1) * prob_satellite_igt(i,g1,t-1);
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

  // Assemble predicted density
  array<Type> dhat_st( n_s, n_t );
  dhat_st.setZero();
  for( int t=0; t<n_t; t++ ){
    if( t==0 ){
      for( int g=0; g<n_g; g++ ){
        dhat_st(g,t) = exp(Beta_t(t));
      }
    }else{
      for( int g1=0; g1<n_g; g1++ ){
      for( int g2=0; g2<n_g; g2++ ){
        dhat_st(g2,t) += exp(Beta_t(t)) * Movement_gg(g2,g1) * exp(ln_d_st(g1,t-1));
      }}
    }
    for( int s_extra=n_g; s_extra<n_s; s_extra++ ){
      dhat_st(s_extra,t) = exp(Beta_t(t));
    }
  }

  // Anisotropy elements
  matrix<Type> H(2,2);
  H(0,0) = exp(ln_H_input(0));
  H(1,0) = ln_H_input(1);
  H(0,1) = ln_H_input(1);
  H(1,1) = (1.0+ln_H_input(1)*ln_H_input(1)) / exp(ln_H_input(0));

  // Survey data
  Eigen::SparseMatrix<Type> Q( n_s, n_s );
  GMRF_t<Type> gmrf_Q;
  Q = Q_spde(spde_aniso, exp(ln_kappa), H);
  gmrf_Q = GMRF( Q );
  Type logtau = log( 1.0 / (exp(ln_kappa) * sqrt(4.0*M_PI)) );

  // Log-likelihood from GMRF
  nll_t(0) = SCALE( gmrf_Q, exp(-logtau) * exp(ln_sigma_epsilon0) )( ln_d_st.col(0) - log(dhat_st.col(0)) );
  for( int t=1; t<n_t; t++ ){
    nll_t(t) = SCALE( gmrf_Q, exp(-logtau) * exp(ln_sigma_epsilon) )( ln_d_st.col(t) - log(dhat_st.col(t)) );
  }

  // Log-likelihood for survey data
  Type phi = exp(ln_phi);
  Type power = 1.0 + invlogit( power_prime );
  for( int j=0; j<n_j; j++ ){
    nll_j(j) = -1 * dtweedie( b_j(j), exp(ln_d_st(g_j(j),t_j(j))), phi, power, true );
  }

  // Calculate annualized movement
  array<Type> Mannual_ggt( n_g, n_g, n_t );
  for( int t1=0; t1<n_t; t1++ ){
    Mprime_gg.setZero();
    for( int t2=t1; (t2<n_t) & (t2<(t1+n_u)); t2++ ){
    for( int g1=0; g1<n_g; g1++ ){
    for( int g2=0; g2<n_g; g2++ ){
      Mprime_gg(g1,g2) = Mprime_ggt(g1,g2,t2);
    }}}
    Mprime_gg = matrix_exponential( n_g, log2steps, Mprime_gg );
    for( int g1=0; g1<n_g; g1++ ){
    for( int g2=0; g2<n_g; g2++ ){
      Mannual_ggt(g1,g2,t1) = Mprime_gg(g1,g2);
    }}
  }
  Msum_gg = matrix_exponential( n_g, log2steps, Mprimesum_gg );

  // Cumulate jnll
  jnll += sum(nll_i);
  jnll += sum(nll_j);
  jnll += sum(nll_t);

  // Report stuff out
  REPORT( sigma2 );
  REPORT( alpha_k );
  REPORT( Diffusion_gg );
  REPORT( Taxis_gg );
  REPORT( Mprime_ggt );
  REPORT( Movement_gg );
  REPORT( Mprimesum_gg );
  REPORT( Mannual_ggt );
  REPORT( Msum_gg );
  //REPORT( logprob_satellite_igt );
  REPORT( prob_satellite_igt );
  REPORT( Preference_gt );
  REPORT( dhat_st );
  REPORT( jnll );
  REPORT( nll_i );
  REPORT( nll_j );
  REPORT( nll_t);

  return jnll;
}
