# mrlocus 0.0.26

* Note that the homepage for v0.0.26 and onward has moved to:
  https://thelovelab.github.io/mrlocus

# mrlocus 0.0.25

* Add convenience function `normalizedAllelicSpread`, a useful
  statistic for comparing sigma across different GWAS datasets.

# mrlocus 0.0.23

* Fixing small bug: trimClusters didn't trim `r2` as data.frame.

# mrlocus 0.0.22

* Added `trimClusters` helper function, with default
  to trim down from most significant signal cluster to
  the least.
* Warning if < 2 signal clusters provided (not recommended).

# mrlocus 0.0.17

* `collapseHighCorSNPs` now outputs the SNP ids that were collapsed.

# mrlocus 0.0.16

* Adding prior predictive checks at the suggestion of Rob Moccia.

# mrlocus 0.0.15

* Adapting `sd_sigma` argument in `fitSlope` based on scale of the
  posterior coefficients for B, `beta_hat_b`.
* Also output priors from `fitSlope`.

# mrlocus 0.0.14

* This is the first major release.
