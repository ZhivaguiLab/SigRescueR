functions {
  // Cosine similarity function
  real cosine_similarity(vector x, vector y) {
    return dot_product(x, y) / (sqrt(dot_product(x, x)) * sqrt(dot_product(y, y)));
  }

  // approximate inversion: mu -> lambda
  real mu_to_lambda(real mu, real nu) {
    return pow(mu + (nu - 1) / (2 * nu), nu);
  }

  // Approximate log normalizing constant Z(lambda, nu)
  real log_Z_compute(real lambda, real nu, int k_max) {
    vector[k_max+1] terms;
    for (k in 0:k_max) {
      terms[k+1] = k * log(lambda) - nu * lgamma(k + 1);
    }
    return log_sum_exp(terms);
  }

  // Vectorized COM-Poisson log PMF for arrays
  real com_poisson(int[] y, vector lambda, real nu, int[] k_max) {
    int N = size(y); // how many length
    vector[N] logp; // log probabilty by mutation context
    for (n in 1:N) {
      real lambda_n = mu_to_lambda(lambda[n], nu);
      logp[n] = y[n] * log(lambda_n)
                - nu * lgamma(y[n] + 1)
                - log_Z_compute(lambda_n, nu, k_max[n]);
    }
    return sum(logp);   // sum over all observations
  }
}

data {
  int<lower=1> N;                          // number of mutation channels (96)
  int<lower=1> K;                          // number of background signatures
  int<lower=0> treatment[N];               // observed mutation counts for 1 sample
  matrix<lower=0>[N, K] s_b;               // normalized background signature
  real<lower=0> total_mut;                 // total mutation in treatment
  int<lower=0> k_max[N];                   // k number for normalization constant... Z(lambda,v)
}

parameters {
  vector<lower=0>[K] theta_b;              // exposure to background signature
  real<lower=0> theta_t;                   // exposure to treatment signature (latent)
  simplex[N] s_t;                          // treatment only signature (latent variable) - normalized b/c simplex (sums to 1)
  //simplex[N] s_b_latent[K];                // K rows, each sums to 1
  real<lower=1.1, upper=10> nu;            // lower bound 1 ensures underdispersion
}

transformed parameters {
  vector<lower=0>[N] lambda;               // profile count based
  vector<lower=0>[N] reconstructed;        // profile not count based
  vector[K] norm_theta_b;
  real norm_theta_t;

  real total_theta = sum(theta_b) + theta_t;

  norm_theta_b = theta_b / total_theta;  // normalize
  norm_theta_t = theta_t / total_theta;  // normalize

  vector[N] total = rep_vector(0, N);      // initialize to zero to total all background
  for (k in 1:K) {
    total += norm_theta_b[k] * s_b[,k];//s_b_latent[k];         // for each bg sig, sum of all activity * bg signature
  }

  reconstructed = total + (norm_theta_t * s_t);   // bg sig + tx sig
  lambda = total_mut * reconstructed;             // count based of reconstructed using observed total mut
}

model {
  // Priors
  theta_b ~ gamma(2, 2);                   // mean = 1,
  theta_t ~ gamma(1, 5);                   // mean = 0.2, less importance in reconstruction
  nu ~ normal(2,1);                        // mean = 2, sd 1

  // Model Count
  //treatment ~ poisson(lambda);
  target += com_poisson(treatment, lambda, nu, k_max);

  //Model Shape
  //Observed normalized
  vector[N] treatment_pseudo = fmax(to_vector(treatment), 1e-3); // Add pseudocouont, dirichlet doesn't like negative
  vector[N] treatment_prop = treatment_pseudo / sum(treatment_pseudo);
  //Reconstructed normalized
  vector[N] reconstructed_pseudo = fmax(reconstructed, 1e-3);
  vector[N] reconstructed_prop = reconstructed_pseudo / sum(reconstructed_pseudo);

  treatment_prop ~ dirichlet(reconstructed_prop * 100);  // 100 = shape precision

  // Reward if cos_sim is greater than 0.95
  real cos_sim = cosine_similarity((to_vector(treatment)/sum(to_vector(treatment))), (reconstructed/sum(reconstructed)));
  target += normal_lpdf(cos_sim | 1, 0.0255);

}

generated quantities {
  real cos_sim = cosine_similarity((to_vector(treatment)/sum(to_vector(treatment))), (reconstructed/sum(reconstructed)));
}
