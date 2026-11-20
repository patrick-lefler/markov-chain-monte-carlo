### The Brilliance of Markov Chain Monte Carlo
#### Taming High-Dimensional Risk Models Through Random Walks
Author: Patrick Lefler 

Original Publish Date: 2026-09-24

Rendered link: https://patrick-lefler.com/markov-chain-monte-carlo/

### Project Introduction

An R Metropolis sampler replaces intractable grid integration with a random walk, pricing 10-factor portfolio VaR and Expected Shortfall in seconds.

### Overview

Correlated multi-asset portfolio risk requires a joint posterior over every 
factor exposure at once, and grid integration scales exponentially with 
dimension: about ten billion evaluations at ten factors. This project builds 
a Random-Walk Metropolis sampler from first principles in pure R and applies 
it to a ten-factor portfolio spanning equity, credit, duration, and commodity 
exposures.

The model is deliberately conjugate. With covariance fixed at its sample 
estimate, the exact posterior of expected returns is known in closed form, so 
every MCMC output (means, standard deviations, credible intervals, VaR, and 
Expected Shortfall) can be checked against the true answer. The project is 
therefore a model-validation case study: it shows how a sampler with a 
textbook acceptance rate can still fail to converge and understate risk, and 
how multi-chain diagnostics and an exact benchmark expose that failure before 
it reaches a risk report.

### Tech Stack

Language: R

Framework: Quarto multi-page website (Sandstone theme, floating sidebar)

Primary Libraries: dplyr, ggplot2, tidyr, tibble, readr, scales, gt, reactable, posterior, MASS, sessioninfo

Environment: renv (package versions pinned in `renv.lock`)

Deployment/Output: Five-page HTML site rendered to `docs/` (index.qmd, 01-mechanics.qmd, 02-risk-application.qmd, 03-diagnostics.qmd, 04-reproducibility.qmd)

### Repository Structure

```
markov-chain-monte-carlo/
├── _quarto.yml                  # Quarto configuration (pages, sidebar, theme, freeze)
├── _brand.yml                   # Colors and fonts
├── index.qmd                    # Intuition & 2D random walk visualization
├── 01-mechanics.qmd             # Metropolis-Hastings engine & step-size tuning
├── 02-risk-application.qmd      # 10-factor portfolio risk model & tail estimates
├── 03-diagnostics.qmd           # Convergence diagnostics (split-R-hat, ACF, ESS)
├── 04-reproducibility.qmd       # Rebuild steps, environment, seeds, limitations
├── license.qmd                  # Renders the LICENSE file
├── LICENSE                      # MIT License
├── renv.lock                    # Pinned package versions
├── R/
│   ├── diagnostics.R            # ESS (Geyer) and split-R-hat
│   ├── generate_data.R          # Generates the synthetic 10-asset returns file
│   ├── metropolis_sampler.R     # Pure R Metropolis-Hastings sampler
│   ├── plots.R                  # ggplot2 theme and plotting functions
│   ├── repro_env.R              # Environment table for the Reproducibility page
│   ├── risk_functions.R         # Joint log-posterior & VaR / Expected Shortfall
│   ├── verify_diagnostics.R     # Manual sanity checks for the ESS estimator
│   ├── verify_risk_functions.R  # Manual sanity checks for the risk functions
│   └── verify_sampler.R         # Manual sanity checks for the sampler
├── assets/images/               # Logo and favicon
├── data/
│   └── synthetic_returns.csv    # Frozen 250-day x 10-factor returns matrix
└── docs/                        # Rendered site
```

To rebuild: `Rscript -e 'renv::restore()'`, then `quarto render`. See the Reproducibility page for details.
 
### Key Findings
- **Point estimates are right.** The sampler reproduces the exact posterior 
  means of portfolio VaR (95%) at 1.05% and Expected Shortfall (95%) at 1.33%, 
  agreeing with the closed-form answer to within about half a basis point.
- **Uncertainty is understated.** The MCMC 90% credible interval for VaR 
  (0.99%–1.11%) is about 17% narrower than the exact interval (0.98%–1.12%), 
  and nine of ten factors show posterior standard deviations 3–18% below 
  the true values. On a $100M portfolio, the exact interval represents about 
  $140,000 of VaR uncertainty from mean-estimation error alone; the sampler 
  reports roughly $120,000.
- **The chains have not converged.** Across four chains started at ±3 
  posterior standard deviations, all ten factors fail the split-R̂ < 1.01 
  test (maximum 1.057, Tech Growth), and the 500-iteration burn-in is too 
  short for overdispersed starts.
- **Effective sample size is far too low.** Pooled across four chains 
  (10,000 draws), bulk ESS is only 65–316 for the audited factors, below the 
  recommended 400. On tail ESS only Treasuries clear 400 (640); Tech Growth is 132 and US Large Cap 244. 
  Estimates were cross-checked against the `posterior` package.
- **A healthy acceptance rate is not evidence of convergence.** The chain 
  accepts 28.9% of proposals, within the conventional target, yet mixes 
  poorly because a single isotropic step size cannot match factors with very 
  different scales and strong correlations. Covariance-preconditioned or 
  gradient-based samplers are the indicated fix.
- **Parameter recovery is a data problem, not a sampler problem.** With 250 
  daily observations, Tech Growth's posterior sits 2.8 standard deviations 
  above its true generating value. One year of data cannot reliably separate 
  high-drift and low-drift assets.

### License

This project is licensed under the MIT License. See the LICENSE file for details.

### Contact

Patrick Lefler [https://www.linkedin.com/in/patricklefler/] | [patrick-lefler.github.io] | [https://substack.com/@pflefler]
