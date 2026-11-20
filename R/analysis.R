# R/analysis.R
# -----------------------------------------------------------------------------
# Runs the 10-factor analysis once and returns every result the pages quote,
# so all pages share one set of settings, seeds and numbers.
# Requires R/generate_data.R, R/metropolis_sampler.R, R/risk_functions.R and
# R/diagnostics.R to be sourced first. Runs in about a second.
# -----------------------------------------------------------------------------

analysis_settings <- function() {
  list(
    n_iter         = 3000,
    burn_in        = 500,
    proposal_sd    = 0.00032,
    start_sd       = 3,               # overdispersed starts, in posterior SDs
    seed_main      = 42,              # chain 1, started at zero
    seeds_extra    = c(99, 123, 456), # chains 2-4: all high, all low, mixed
    seed_signs     = 2026,            # random signs for the "mixed" start
    alpha          = 0.95,            # VaR / ES confidence level
    rhat_threshold = 1.01,            # Vehtari et al. (2021)
    ess_target     = 400,             # Vehtari et al. (2021)
    audit_factors  = c("Tech_Growth", "US_Treas_10Y", "US_LargeCap")
  )
}

run_analysis <- function(csv_path = "data/synthetic_returns.csv", settings = analysis_settings()) {
  s <- settings
  fixture     <- build_risk_target(csv_path)
  d           <- fixture$d
  asset_names <- fixture$asset_names
  R_mat       <- fixture$returns_mat
  S           <- fixture$cov_emp
  n_obs       <- nrow(R_mat)
  weights     <- rep(1 / d, d)
  n_keep      <- s$n_iter - s$burn_in

  # Chains: chain 1 from zero, chains 2-4 from overdispersed starts
  run_chain <- function(start, seed) {
    names(start) <- asset_names
    run_metropolis(fixture$target_fn, start,
                   iterations = s$n_iter, proposal_sd = s$proposal_sd, seed = seed)
  }
  post_center <- colMeans(R_mat)
  post_sd     <- sqrt(diag(S) / n_obs)   # flat-prior approximation of the posterior SD

  set.seed(s$seed_signs)
  signs  <- sample(c(-1, 1), d, replace = TRUE)
  starts <- list(
    post_center + s$start_sd * post_sd,          # all high
    post_center - s$start_sd * post_sd,          # all low
    post_center + s$start_sd * post_sd * signs   # mixed
  )
  chains <- c(
    list(run_chain(rep(0, d), s$seed_main)),
    purrr::map2(starts, s$seeds_extra, run_chain)
  )
  draws_clean <- purrr::map(chains, \(ch) dplyr::slice_tail(ch$draws, n = n_keep))
  post_burn   <- draws_clean[[1]]

  # Portfolio tail risk from chain 1
  risk_table <- evaluate_portfolio_tail_risk(
    draws_df = post_burn, weights = weights, cov_mat = S, alpha = s$alpha
  )

  # Exact conjugate posterior of mu and of portfolio VaR / ES
  S_inv  <- solve(S)
  V_post <- solve(n_obs * S_inv + fixture$inv_prior_var * diag(d))
  m_post <- setNames(as.numeric(V_post %*% (S_inv %*% colSums(R_mat))), asset_names)
  sd_post <- sqrt(diag(V_post))
  sig_p  <- sqrt(as.numeric(t(weights) %*% S %*% weights))
  mu_p   <- sum(weights * m_post)
  sd_mup <- sqrt(as.numeric(t(weights) %*% V_post %*% weights))
  k_var  <- qnorm(s$alpha)
  k_es   <- dnorm(qnorm(s$alpha)) / (1 - s$alpha)
  z_ci   <- qnorm(0.95)   # 90% interval, matching the 5th-95th percentiles in evaluate_portfolio_tail_risk()
  exact_var <- -mu_p + k_var * sig_p
  exact_es  <- -mu_p + k_es * sig_p
  exact <- list(
    m_post = m_post, sd_post = sd_post, V_post = V_post,
    var = exact_var, es = exact_es, sd_mup = sd_mup, z_ci = z_ci,
    k_var = k_var, k_es = k_es,
    var_ci = exact_var + c(-1, 1) * z_ci * sd_mup,
    es_ci  = exact_es + c(-1, 1) * z_ci * sd_mup
  )
  width_ratio <- (risk_table$ci_95[1] - risk_table$ci_05[1]) / (2 * z_ci * sd_mup)

  draws_mat <- as.matrix(post_burn[, asset_names])
  factor_bench <- tibble::tibble(
    Factor     = asset_names,
    MCMC_Mean  = colMeans(draws_mat),
    Exact_Mean = m_post,
    MCMC_SD    = apply(draws_mat, 2, sd),
    Exact_SD   = sd_post
  ) |>
    dplyr::mutate(
      SD_Ratio      = MCMC_SD / Exact_SD,
      Mean_Error_SD = (MCMC_Mean - Exact_Mean) / Exact_SD
    )

  # Convergence diagnostics
  rhat <- purrr::map_dbl(setNames(asset_names, asset_names), \(a)
    compute_split_rhat(purrr::map(draws_clean, \(dc) dc[[a]])))
  chain_matrix <- function(f) sapply(draws_clean, \(dc) dc[[f]])
  af <- setNames(s$audit_factors, s$audit_factors)
  ess_chain1  <- purrr::map_dbl(af, \(f) compute_ess_robust(post_burn[[f]]))
  ess_basic   <- purrr::map_dbl(af, \(f) round(posterior::ess_basic(post_burn[[f]]), 1))
  ess_bulk    <- purrr::map_dbl(af, \(f) round(posterior::ess_bulk(chain_matrix(f))))
  ess_tail    <- purrr::map_dbl(af, \(f) round(posterior::ess_tail(chain_matrix(f))))

  list(
    settings = s, fixture = fixture, d = d, asset_names = asset_names,
    n_obs = n_obs, n_keep = n_keep, n_chains = length(chains), weights = weights,
    post_center = post_center, post_sd = post_sd,
    chains = chains, draws_clean = draws_clean, post_burn = post_burn,
    acceptance_rate = chains[[1]]$acceptance_rate,
    risk_table = risk_table, exact = exact, width_ratio = width_ratio,
    factor_bench = factor_bench,
    rhat = rhat, n_rhat_fail = sum(rhat >= s$rhat_threshold),
    ess_chain1 = ess_chain1, ess_basic = ess_basic,
    ess_bulk = ess_bulk, ess_tail = ess_tail
  )
}
