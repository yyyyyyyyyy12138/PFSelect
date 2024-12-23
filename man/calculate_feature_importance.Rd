% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PFselect.R
\name{calculate_feature_importance}
\alias{calculate_feature_importance}
\title{Calculate Local Feature Importance}
\usage{
calculate_feature_importance(
  coefs,
  genetic_variants,
  genetic_variants_knockoff,
  interaction_terms,
  interaction_terms_knockoff,
  unpenalized_covariates,
  Z
)
}
\arguments{
\item{coefs}{Coefficients extracted from Lasso model.}

\item{genetic_variants}{Genetic variants matrix.}

\item{genetic_variants_knockoff}{Knockoff genetic variants matrix.}

\item{interaction_terms}{Interaction terms for original genetic variants.}

\item{interaction_terms_knockoff}{Interaction terms for knockoff genetic variants.}

\item{unpenalized_covariates}{Matrix of covariates (not penalized).}

\item{Z}{Heterogeneity variable (e.g., EUR or PCs).}
}
\value{
List with feature importance matrices for original and knockoff variants.
}
\description{
Calculate Local Feature Importance
}
