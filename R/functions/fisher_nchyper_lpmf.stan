// fisher_hyper_lpmf.stan
real fisher_hyper_lpmf(int k, real mu, int m1, int m2, int N) {
  if (mu <= 0) return not_a_number();
  
  // Bounds for hypergeometric distribution
  int k_min = max({0, m1 + m2 - N});
  int k_max = min({m1, m2});
  
  // Safety check for support bounds
  if (k < k_min || k > k_max) return negative_infinity();
  
  // Log-numerator for the observed value k
  real log_num = lchoose(m1, k) + lchoose(N - m1, m2 - k) + k * log(mu);
  
  // Normalization constant calculation
  int range_size = k_max - k_min + 1;
  vector[range_size] terms;
  
  // FIX: Explicitly declare 'int' for the loop iterator
  for (i in 0:(range_size - 1)) {
    int curr_i = k_min + i;
    terms[i + 1] = lchoose(m1, curr_i) + lchoose(N - m1, m2 - curr_i) + curr_i * log(mu);
  }
  
  return log_num - log_sum_exp(terms);
}
