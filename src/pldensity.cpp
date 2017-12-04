#define ARMA_64BIT_WORD 1
#define _USE_MATH_DEFINES

#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 


// // This inline function creates one sample of multivariate t-distribution
// // [[Rcpp::export]]
// inline arma::vec rSt(
//     const arma::vec& mu,
//     const arma::mat& Sigma,
//     const double df) {
//   
//   int d = mu.size();
//   vec normal =  chol(Sigma).t() * randn(d);
//   double chi2 = rchisq(1, df)[0];
//   return mu + normal / sqrt(chi2 / df);
// }

// This inline function evaluates the density of t distribution
inline double dst(
    const arma::vec& x,
    const arma::vec& mu,
    const arma::mat& Sigma, // Sigma^{-1}
    const double df) {
  
  int d = mu.n_elem;
  vec xcentered = x - mu;
  double innerterm =  - 0.5 * (df + d) * 
    log(1.0 + as_scalar(xcentered.t() * inv_sympd(Sigma) * xcentered) / df);
  double ldet;
  double sign;
  log_det(ldet, sign, Sigma); // compute and store the logarithm of determinant
  double extterm = log(tgamma(0.5 * (df + d))) - log(tgamma(0.5 * df)) - 0.5 * d * log(df * M_PI) - 0.5 * ldet;
  double density = exp(extterm + innerterm);
  return !isnan(density) ? density : 0.0;
}

// Samples from a multivariate categorical variable
inline arma::uvec resample(int N, vec prob) {
  vec probsum = cumsum(prob) / sum(prob);
  uvec out(N);
  for (int i = 0; i < N; i++) {
    double u = unif_rand();
    int j = 0;
    while (u > probsum[j]) {j++;}
    out[i] = j;
  }
  return out;
}

// Defines a particle, which carries a essensial state
struct Particle {
  int m; // Total clusters
  vec n; // Observations per cluster
  mat mu; // Mean (dim1) per cluster
  cube S; // SS (dim1, dim2) per cluster
  Particle (int m_, arma::vec& n_, arma::mat& mu_, arma::cube& S_) 
    : m(m_), n(n_), mu(mu_), S(S_) {};
  Particle (const arma::vec& x) {
    m = 1;
    n.resize(1);
    n[0] = 1;
    mu.resize(x.n_elem, 1);
    mu.col(0) = x;
    S.resize(x.n_elem, x.n_elem, 1);
    S.slice(0).fill(0.0);
  }
  Particle (const Particle& z) { // copy constructor
    m = z.m;
    n = z.n;
    mu = z.mu;
    S = z.S;
  }; 
  Particle () {};
};

// Defines the prior, to avoid passing along all the parameters
struct DPNormalPrior { // mu ~ N(lambda, S / kappa), S^-1 ~ W(Omega, nu)
  const double alpha;
  const arma::vec& lambda; 
  const double kappa; 
  const double nu;
  const arma::mat& Omega;
  DPNormalPrior(
    const double a_,
    const arma::vec& l_,
    const double k_, 
    double n_, 
    const arma::mat& O_) 
    : alpha(a_), lambda(l_), kappa(k_), nu(n_), Omega(O_) {}
};

// Evaluate contricution to predictive of cluster, j == z.m is the prior
inline arma::vec predictive(
    const arma::vec& x,
    const Particle& z,
    const DPNormalPrior& p
) {
  // Initialize
  double d = x.n_elem;
  vec dpred(z.m + 1);
  
  // Cluster contribution
  for (int j = 0; j < z.m; j++) {
    vec aj = (p.kappa * p.lambda + z.n[j] * z.mu.col(j)) / (p.kappa + z.n[j]);
    mat Dj = z.S.slice(j) + p.kappa * z.n[j] / 
      (p.kappa +  z.n[j]) * (p.lambda - z.mu.col(j)) * (p.lambda - z.mu.col(j)).t();
    double cj = 2 * p.nu + z.n[j] - d + 1.0;
    mat Bj = 2.0 * (p.kappa + z.n[j] + 1.0) / (p.kappa + z.n[j]) / 
      cj * (p.Omega  + 0.5 * Dj);
    dpred[j] = z.n[j] * dst(x, aj, Bj, cj);
  }
  
  // Prior contribution
  vec a0 = p.lambda;
  double c0 = 2.0 * p.nu - d + 1.0;
  const mat& B0 = 2.0 * (p.kappa + 1.0) / p.kappa / c0 * p.Omega;
  dpred[z.m] = p.alpha * dst(x, a0, B0, c0);
  
  return dpred;
}


// Evaluates the density of a new point given a particle
inline Particle propagate(
    const arma::vec& xnew, 
    const Particle& z,
    const DPNormalPrior& p // iteration timestamp
) {
  // Initialize output
  Particle out(z);
  
  // Propagate allocation
  vec cp = predictive(xnew, z, p);
  int k = resample(1, cp)[0]; 
  
  // If new cluster
  if (k == z.m) {
    out.n.insert_rows(out.m, 1);
    out.mu.insert_cols(out.m, xnew);
    out.S.insert_slices(out.m, 1);
    out.m += 1;
    out.n[k] += 1;
  } else {
    out.n[k] += 1;
    out.mu.col(k) =  (z.n[k] * z.mu.col(k) + xnew) / out.n[k];
    out.S.slice(k) += (xnew * xnew.t()) + z.n[k] * (z.mu.col(k) * z.mu.col(k).t()) -  
      out.n[k] * (out.mu.col(k) * out.mu.col(k).t());  
  }
  
  return out;
}

//' @title Dirichlet Process Normal Mixture Kernel Density Estimation
//' @description Bla Bla
//' @param x a matrix of observations
//' @param alpha concentration parameter of Dirichlett process
//' @export
// [[Rcpp::export]]
Rcpp::List dp_normal_mix(
    const arma::mat& x,
    const int N,
    const double alpha,
    const arma::vec& lambda,
    const double kappa,
    const double nu,
    const arma::mat& Omega
) {
  // Total observations
  const int T = x.n_rows;
  const int d = x.n_cols;
  
  // Save prior in structure
  const DPNormalPrior prior(alpha, lambda, kappa, nu, Omega);
  
  // Initialize N particles
  std::vector<Particle> particle(N);
  for (int i = 0; i < N; i++) {
    vec randstart(d);
    for (int j = 0; j < d; j++) {
      randstart[j] = unif_rand();
    }
    particle[i] = Particle(randstart);
  }
  
  // Update every particle
  for (int t = 1; t < T; t++) {
    // Resample 
    vec weight(N);
    for (int i = 0; i < N; i++) {
      weight[i] = sum(predictive(x.row(t).t(), particle[i], prior));
    }
    uvec new_idx = resample(N, weight);
    
    // Propagate
    std::vector<Particle> temp_particle(particle);
    for (int i = 0; i < N; i++) {
      particle[i] = propagate(x.row(t).t(), temp_particle[new_idx[i]], prior);
    }
  }
  
  // Wrap all the output in lists for R
  Rcpp::List param = Rcpp::List::create(
    Named("alpha") = alpha,
    Named("lambda") = wrap(lambda),
    Named("kappa") = wrap(kappa),
    Named("nu") = nu,
    Named("Omega") = wrap(Omega)
  );
  
  Rcpp::List particle_list;
  for (int i = 0; i < N; i++) {
    Rcpp::List particlei = Rcpp::List::create(
      Named("m") = particle[i].m,
      Named("n") = wrap(particle[i].n),
      Named("mu") = wrap(particle[i].mu),
      Named("S") = wrap(particle[i].S));
    particle_list.push_back(particlei);
  }
  
  Rcpp::List out = Rcpp::List::create(
    Named("N") = N,
    Named("param") = param,
    Named("particle_list") = particle_list
  );
  
  out.attr("class") = "PL";
  return out;
}

//' @title Eval Point Density
//' @description Bla Bla
//' @export
// [[Rcpp::export]]
arma::vec dp_normal_deval(
    const Rcpp::List& model,
    const arma::mat& xnew,
    const int nparticles = 50
  ) {
  // Check list is of appropriate type
  if (!Rf_inherits(model, "PL")) throw "List must be a PL object";

  // Allocate space
  const uword K = xnew.n_rows;
  vec density(K, fill::zeros);

  // Prior structure
  const int N = model["N"];
  Rcpp::List param = model["param"];
  vec lambda = param["lambda"];
  mat Omega = param["Omega"];
  double alpha = param["alpha"], kappa = param["kappa"], nu = param["nu"];
  const DPNormalPrior prior(alpha, lambda, kappa, nu, Omega);
  const int N0 = min(N, nparticles);
  
  Rcpp::List particle_list = model["particle_list"];
  
  for (int i = 0; i < N0; i++) {
    Rcpp::List particle_param = particle_list[i];
    int m = particle_param["m"];
    vec n = particle_param["n"];
    mat mu = particle_param["mu"];
    cube S = particle_param["S"];
    Particle z(m, n, mu, S);
    for (uword k = 0; k < K; k++) {
      density(k) += sum(predictive(xnew.row(k).t(), z, prior)) / N0;
    }
  }
  
  return density;
}

//' @title Take marginals
//' @description Bla Bla
//' @export
// [[Rcpp::export]]
Rcpp::List dp_normal_marginal(
    const Rcpp::List& model,
    const arma::uvec& dims
) {
  // Check list is of appropriate type
  if (!Rf_inherits(model, "PL")) throw "List must be a PL object";
  const uvec dims_ = dims - 1;
  
  // Parameters
  const int N = model["N"];
  Rcpp::List param_full = model["param"];
  const double alpha = param_full["alpha"], kappa =  param_full["kappa"], nu = param_full["nu"];
  vec lambda = param_full["lambda"];
  vec lambdanew = lambda.elem(dims_);
  mat Omega = param_full["Omega"];
  mat Omeganew = Omega.submat(dims_, dims_);
  Rcpp::List param_reduced = Rcpp::List::create(
    Named("alpha") = alpha, 
    Named("lambda") = wrap(lambdanew),
    Named("kappa") = kappa,
    Named("nu") = nu,
    Named("Omega") = wrap(Omeganew)
  );

  // Copy only selected dimensions
  Rcpp::List particle_list_reduced;
  Rcpp::List particle_list_full = model["particle_list"];
  for (int i = 0; i < N; i++) {
    Rcpp::List particle = particle_list_full[i];
    int m = particle["m"];
    vec n = particle["n"];
    mat mu = particle["mu"];
    mat munew = mu.rows(dims_);
    cube S = particle["S"];
    cube Snew(dims_.n_elem, dims_.n_elem, m);
    for (int k = 0; k < m; k++) {
      Snew.slice(k) = S.slice(k).submat(dims_, dims_);
    }
    Rcpp::List particlei = Rcpp::List::create(
      Named("m") = m,
      Named("n") = n,
      Named("mu") = wrap(munew),
      Named("S") = wrap(Snew));
    particle_list_reduced.push_back(particlei);
  }

  Rcpp::List model_reduced = Rcpp::List::create(
    Named("N") = N,
    Named("param") = param_reduced,
    Named("particle_list") = particle_list_reduced
  );

  model_reduced.attr("class") = "PL";
  return model_reduced;
}


//' @title Eval Point Conditional density
//' @description Bla Bla
//' @export
// [[Rcpp::export]]
arma::mat dp_normal_deval_conditional(
    const Rcpp::List& model,
    const arma::mat& xnew,
    const arma::uvec& eval_dims,
    const arma::uvec& condition_dims,
    const arma::mat& condition_values,
    const int nparticles = 50
) {
  // Check list is of appropriate type
  if (!Rf_inherits(model, "PL"))
    throw "List must be a PL object";

  // Allocate space
  const uvec edims = eval_dims - 1;
  const uvec cdims = condition_dims - 1;
  const int total_dims = edims.n_elem + cdims.n_elem;
  const uword K = xnew.n_rows, L = condition_values.n_rows;
  mat density(K, L, fill::zeros);

  // Marginal model
  Rcpp::List marginal = dp_normal_marginal(model, cdims + 1);

  // Priors
  Rcpp::List param = model["param"];
  vec lambda = param["lambda"];
  vec mlambda = lambda(cdims);
  mat Omega = param["Omega"];
  mat mOmega = Omega.submat(cdims, cdims);
  const DPNormalPrior prior((double) param["alpha"], lambda, (double) param["kappa"], (double) param["nu"], Omega);
  const DPNormalPrior mprior((double) param["alpha"], mlambda, (double) param["kappa"], (double) param["nu"], mOmega);

  //  Eval conditional densities
  const int N0 = min((int) model["N"], nparticles);
  Rcpp::List particle_list = model["particle_list"];
  Rcpp::List mparticle_list = marginal["particle_list"];
  for (int i = 0; i < N0; i++) {
    Rcpp::List particle_param = particle_list[i];
    Rcpp::List mparticle_param = mparticle_list[i];
    int m = particle_param["m"];
    vec n = particle_param["n"];
    mat mu = particle_param["mu"];
    mat mmu = mparticle_param["mu"];
    cube S = particle_param["S"];
    cube mS = mparticle_param["S"];
    Particle zfull(m, n, mu, S);
    Particle zmarginal(m, n, mmu, mS);
    for (uword l = 0; l < L; l++) {
      for (uword k = 0; k < K; k++) {
        vec xfull(total_dims, fill::zeros);
        xfull(edims) = xnew.row(k).t();
        xfull(cdims) = condition_values.row(l).t();
        double djoint = sum(predictive(xfull, zfull, prior));
        double dmarginal = sum(predictive(condition_values.row(l).t(), zmarginal, mprior));
        density(k, l) += djoint / dmarginal / N0;
      }
    }
  }
  
  return density;
}