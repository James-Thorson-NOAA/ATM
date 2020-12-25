#include <TMB.hpp>

// Function to import R list for user-defined Options_vec and Options, packaged as list Options_list in TmbData
// Try asScalarType
//template<class Type>
//struct options_list {
//  vector<Type> constant_tail_probability;
//  vector<int> log2steps;
//  vector<Type> alpha_ratio_bounds;
//  options_list(SEXP x){ // Constructor
//    constant_tail_probability = asVector<Type>(getListElement(x,"constant_tail_probability"));;
//    log2steps = asVector<int>(getListElement(x,"log2steps"));;
//    alpha_ratio_bounds = asVector<Type>(getListElement(x,"alpha_ratio_bounds"));;
//  }
//};

// Posfun -- see https://github.com/kaskr/adcomp/wiki/Code--snippets
// See more: https://www.coin-or.org/CppAD/Doc/doxydoc/html/cond__exp_8hpp_source.html
template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0));   // left, right, if-true, if-false
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2)-x/eps));
}

// square
template<class Type>
Type square(Type x){
  return pow(x,2);
}

// sqrt
template<class Type>
Type sqrt(Type x){
  return pow(x,0.5);
}

// new version
//Type posfun_hinge( Type x, Type rate ){
//  return log( exp(Type(0)) + exp(x*rate) ) / rate;
//  return logspace_add(Type(0), x*rate) / rate;
//}

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Simulate from tweedie
// Adapted from tweedie::rtweedie function in R
template<class Type>
Type rtweedie( Type mu, Type phi, Type power){
  Type lambda = pow(mu, Type(2.0) - power) / (phi * (Type(2.0) - power));
  Type alpha = (Type(2.0) - power) / (Type(1.0) - power);
  Type gam = phi * (power - Type(1.0)) * pow(mu, power - Type(1.0));
  Type N = rpois(lambda);
  Type B = rgamma(-N * alpha, gam);   /// Using Shape-Scale parameterization
  return B;
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

  // Options
  //DATA_STRUCT( Options_list, options_list );

  // Data
  DATA_INTEGER( log2steps );
  DATA_SCALAR( constant_tail_probability );
  DATA_SCALAR( alpha_ratio_bounds );
  DATA_INTEGER( diffusion_bounds );
  DATA_INTEGER( report_early );
  DATA_ARRAY( X_guyk );
  DATA_ARRAY( Z_guyl );
  DATA_IMATRIX( uy_tz );
  DATA_IMATRIX( satellite_iz );
  DATA_IMATRIX( conventional_hz );
  DATA_ARRAY( E_guy );
  DATA_VECTOR( duration_u );
  DATA_MATRIX( A_gg );

  // Survey data objects
  DATA_STRUCT( spde_aniso, spde_aniso_t );
  DATA_VECTOR( b_j );
  DATA_IVECTOR( t_j );
  DATA_IVECTOR( g_j );

  // Fishery data objects
  DATA_VECTOR( b_f );
  DATA_IVECTOR( t_f );
  DATA_IVECTOR( g_f );

  // Parameters
  PARAMETER_VECTOR( ln_sigma_l );
  PARAMETER_VECTOR( alpha_logit_ratio_k );

  // SPDE parameters
  PARAMETER_VECTOR(ln_H_input); // Anisotropy parameters
  PARAMETER(ln_kappa);
  PARAMETER(ln_sigma_epsilon0);
  PARAMETER(ln_sigma_epsilon);
  PARAMETER(ln_phi);
  PARAMETER(power_prime);
  PARAMETER(ln_CV);
  PARAMETER(lambda);           // log-ratio of fishery gear performance to survey gear performance
  PARAMETER_VECTOR(Beta_t);
  PARAMETER_ARRAY(ln_d_st);     // Annual log-density

  // Dimensions
  int n_g = X_guyk.col(0).col(0).col(0).size();
  int n_u = X_guyk.col(0).col(0).size() / n_g;
  int n_y = X_guyk.col(0).size() / n_g / n_u;
  int n_k = X_guyk.size() / n_g / n_u / n_y;
  int n_l = Z_guyl.size() / n_g / n_u / n_y;
  int n_t = uy_tz.rows();
  int n_h = conventional_hz.rows();
  int n_i = satellite_iz.rows();
  int n_j = b_j.size();
  int n_f = b_f.size();
  int n_s = ln_d_st.col(0).size();

  // Transform inputs
  //vector<Type> sigma2_l = exp( 2.0 * ln_sigma_l );
  vector<Type> alpha_k( n_k );
  if( alpha_ratio_bounds > 0 ){
    alpha_k = alpha_ratio_bounds * exp(2.0*ln_sigma_l(0)) * (2.0 * invlogit(alpha_logit_ratio_k) - 1.0);
  }else{
    alpha_k = alpha_logit_ratio_k;
  }

  // Global variables
  Type jnll = 0;
  vector<Type> nll_h( n_h );
  vector<Type> nll_i( n_i );
  vector<Type> nll_j( n_j );
  vector<Type> nll_t( n_t );
  vector<Type> nll_f( n_f );
  nll_h.setZero();
  nll_i.setZero();
  nll_j.setZero();
  nll_t.setZero();
  nll_f.setZero();

  if( report_early == 1 ){
    REPORT( alpha_k );
    REPORT( ln_sigma_l );
    return( jnll );
  }

  // Global variables
  array<Type> prob_satellite_igt( n_i, n_g, n_t );
  array<Type> prob_conventional_hgt( n_h, n_g, n_t );
  array<Type> prob_satellite_ig( n_i, n_g );
  array<Type> prob_conventional_hg( n_h, n_g );
  //array<Type> logprob_satellite_igt( n_i, n_g, n_t );
  matrix<Type> Diffusion_gg( n_g, n_g );
  matrix<Type> Taxis_gg( n_g, n_g );
  matrix<Type> Mprime_gg( n_g, n_g );
  array<Type> Mprime_ggt( n_g, n_g, n_t );
  matrix<Type> Mprimesum_gg( n_g, n_g );
  matrix<Type> Movement_gg( n_g, n_g );
  matrix<Type> Preference_gt( n_g, n_t );
  matrix<Type> ln_sigma_gt( n_g, n_t );
  vector<Type> init_g( n_g );
  matrix<Type> Msum_gg( n_g, n_g );
  Mprimesum_gg.setZero();
  Preference_gt.setZero();
  ln_sigma_gt.setZero();

  // Loop through times
  //array<Type> logspaceadd_itgg( n_i, n_t, n_g, n_g );
  for( int t=0; t<n_t; t++ ){
    // Advection-rate matrix
    for( int g=0; g<n_g; g++ ){
    for( int k=0; k<n_k; k++ ){
      Preference_gt(g,t) += X_guyk(g,uy_tz(t,0),uy_tz(t,1),k) * alpha_k(k);
    }}
    for( int g1=0; g1<n_g; g1++ ){
    for( int g2=0; g2<n_g; g2++ ){
      Taxis_gg(g1,g2) = A_gg(g1,g2) * (Preference_gt(g1,t) - Preference_gt(g2,t));
    }}

    // Diffusion-rate matrix
    for( int g=0; g<n_g; g++ ){
    for( int l=0; l<n_l; l++ ){
      ln_sigma_gt(g,t) += Z_guyl(g,uy_tz(t,0),uy_tz(t,1),l) * ln_sigma_l(l);
    }}
    for( int g1=0; g1<n_g; g1++ ){
    for( int g2=0; g2<n_g; g2++ ){
      Diffusion_gg(g1,g2) = A_gg(g1,g2) * exp(2.0 * ln_sigma_gt(g2,t));
    }}
    if( diffusion_bounds == 1 ){
      // Ensure that Diffusion_gg(g1,g2)+Taxis_gg(g1,g2) > 0 for all g1 != g2
      Type minval;
      for( int g2=0; g2<n_g; g2++ ){
        minval = min(vector<Type>(Taxis_gg.col(g2)));
        for( int g1=0; g1<n_g; g1++ ){
          Diffusion_gg(g1,g2) += -1 * A_gg(g1,g2) * minval;
        }
      }
    }

    // Do mass-balance after min for Taxis_gg
    Diffusion_gg = subtract_colsum_from_diagonal( n_g, Diffusion_gg );
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
    // TODO:  Explore logspace_add for numerical stability
    for( int i=0; i<n_i; i++ ){
      if( satellite_iz(i,2) == t ){
        for( int g=0; g<n_g; g++ ){
          prob_satellite_igt(i,g,t) = Movement_gg( g, satellite_iz(i,0) );
          // Low tail probability inflation
          prob_satellite_igt(i,g,t) = constant_tail_probability + (1.0 - n_g*constant_tail_probability) * prob_satellite_igt(i,g,t);
          // Log-space calculations
          //logprob_satellite_igt(i,g,t) = log(Movement_gg( g, satellite_iz(i,0) ));
       }
      }
      if( (satellite_iz(i,2)<t) & (satellite_iz(i,3)>=t) ){
        for( int g2=0; g2<n_g; g2++ ){
          for( int g1=0; g1<n_g; g1++ ){
            prob_satellite_igt(i,g2,t) += Movement_gg(g2,g1) * prob_satellite_igt(i,g1,t-1);
            // Log-space calculationss
            //logspaceadd_itgg(i,t,g1,g2) = logspace_add( logprob_satellite_igt(i,g2,t), log(Movement_gg(g2,g1))+logprob_satellite_igt(i,g1,t-1) );
            //logprob_satellite_igt(i,g2,t) = logspaceadd_itgg(i,t,g1,g2);
          }
          // Low tail probability inflation
          prob_satellite_igt(i,g2,t) = constant_tail_probability + (1.0 - n_g*constant_tail_probability) * prob_satellite_igt(i,g2,t);
          // Log-space calculationss
          //logprob_satellite_igt(i,g2,t) = log(logprob_satellite_igt(i,g2,t));    // Very weird that this is needed, but appears to be !
        }
      }
    }

    // Apply to conventional tags
    for( int h=0; h<n_h; h++ ){
      if( conventional_hz(h,2) == t ){
        for( int g=0; g<n_g; g++ ){
          prob_conventional_hgt(h,g,t) = Movement_gg( g, conventional_hz(h,0) );
          // Low tail probability inflation
          prob_conventional_hgt(h,g,t) = constant_tail_probability + (1.0 - n_g*constant_tail_probability) * prob_conventional_hgt(h,g,t);
        }
      }
      if( (conventional_hz(h,2)<t) & (conventional_hz(h,3)>=t) ){
        for( int g2=0; g2<n_g; g2++ ){
          for( int g1=0; g1<n_g; g1++ ){
            prob_conventional_hgt(h,g2,t) += Movement_gg(g2,g1) * prob_conventional_hgt(h,g1,t-1);
          }
          // Low tail probability inflation
          prob_conventional_hgt(h,g2,t) = constant_tail_probability + (1.0 - n_g*constant_tail_probability) * prob_conventional_hgt(h,g2,t);
        }
      }
    }
  }

  // Calculate log-likelihood for sat-tags
  for( int i=0; i<n_i; i++ ){
    nll_i(i) = -1 * log( prob_satellite_igt(i,satellite_iz(i,1),satellite_iz(i,3)) );
    //nll_i(i) = -1 * logprob_satellite_igt( i, satellite_iz(i,1), satellite_iz(i,3) );
    for( int g=0; g<n_g; g++ ){
      prob_satellite_ig(i,g) = prob_satellite_igt(i,g,satellite_iz(i,3));
    }
  }

  // Calculate log-likelihood for conventional-tags
  for( int h=0; h<n_h; h++ ){
    nll_h(h) = -1 * log( E_guy(conventional_hz(h,1),uy_tz(conventional_hz(h,3),0),uy_tz(conventional_hz(h,3),1)) * prob_conventional_hgt(h,conventional_hz(h,1),conventional_hz(h,3)) );
    for( int g=0; g<n_g; g++ ){
      prob_conventional_hg(h,g) = prob_conventional_hgt(h,g,conventional_hz(h,3));
    }
  }

  // Survey and fishery data -- Skip if survey AND fishery data not present
  if( (n_j>0) | (n_f>0) ){
    // Anisotropy elements
    matrix<Type> H(2,2);
    H(0,0) = exp(ln_H_input(0));
    H(1,0) = ln_H_input(1);
    H(0,1) = ln_H_input(1);
    H(1,1) = (1.0+ln_H_input(1)*ln_H_input(1)) / exp(ln_H_input(0));

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
    REPORT( dhat_st );

    // GMRF precision
    Eigen::SparseMatrix<Type> Q( n_s, n_s );
    GMRF_t<Type> gmrf_Q;
    Q = Q_spde(spde_aniso, exp(ln_kappa), H);
    gmrf_Q = GMRF( Q );
    Type logtau = log( 1.0 / (exp(ln_kappa) * sqrt(4.0*M_PI)) );

    // Log-likelihood from GMRF
    array<Type> dtilda_st( n_s, n_t );
    dtilda_st.col(0) = (ln_d_st.col(0) - log(dhat_st.col(0))) / exp(ln_sigma_epsilon0);
    nll_t(0) = SCALE( gmrf_Q, exp(-logtau) * exp(ln_sigma_epsilon0) )( ln_d_st.col(0) - log(dhat_st.col(0)) );
    for( int t=1; t<n_t; t++ ){
      dtilda_st.col(t) = (ln_d_st.col(t) - log(dhat_st.col(t))) / exp(ln_sigma_epsilon);
      nll_t(t) = SCALE( gmrf_Q, exp(-logtau) * exp(ln_sigma_epsilon) )( ln_d_st.col(t) - log(dhat_st.col(t)) );
    }
    REPORT( dtilda_st );

    // Log-likelihood for survey data
    Type phi = exp(ln_phi);
    Type power = 1.0 + invlogit( power_prime );
    vector<Type> bhat_j( n_j );
    for( int j=0; j<n_j; j++ ){
      bhat_j(j) = exp(ln_d_st(g_j(j),t_j(j)));
      nll_j(j) = -1 * dtweedie( b_j(j), bhat_j(j), phi, power, true );
      SIMULATE{
        b_j(j) = rtweedie( bhat_j(j), phi, power );   // Defined above
      }
      REPORT( bhat_j );
    }

    // Log-likelihood for fishery data
    Type CV = exp(ln_CV);
    vector<Type> bhat_f( n_f );
    for( int f=0; f<n_f; f++ ){
      bhat_f(f) = exp( ln_d_st(g_f(f),t_f(f)) + lambda );
      // shape = 1/CV^2;   scale = mean*CV^2
      nll_f(f) = -1 * dgamma( b_f(f), 1/square(CV), bhat_f(f)*square(CV), true );
      SIMULATE{
        b_f(f) = rgamma( 1/square(CV), bhat_f(f)*square(CV) );   // Defined above
      }
      REPORT( bhat_f );
    }
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
  jnll += sum(nll_h);
  jnll += sum(nll_i);
  jnll += sum(nll_j);
  jnll += sum(nll_t);
  jnll += sum(nll_f);

  // Report stuff out
  //REPORT( sigma2 );
  REPORT( alpha_k );
  REPORT( ln_sigma_l );
  REPORT( Diffusion_gg );
  REPORT( Taxis_gg );
  REPORT( Mprime_ggt );
  REPORT( Movement_gg );
  REPORT( Mprimesum_gg );
  REPORT( Mannual_ggt );
  REPORT( Msum_gg );
  //REPORT( logprob_satellite_igt );
  //REPORT( logspaceadd_itgg );
  REPORT( prob_satellite_igt );
  REPORT( prob_conventional_hgt );
  REPORT( prob_satellite_ig );
  REPORT( prob_conventional_hg );
  REPORT( Preference_gt );
  REPORT( ln_sigma_gt );
  REPORT( jnll );
  REPORT( nll_h );
  REPORT( nll_i );
  REPORT( nll_j );
  REPORT( nll_t);
  REPORT( nll_f);
  ADREPORT( alpha_k );

  SIMULATE{
    REPORT( b_j );
  }

  return jnll;
}
