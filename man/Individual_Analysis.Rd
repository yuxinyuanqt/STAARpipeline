% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Individual_Analysis.R
\name{Individual_Analysis}
\alias{Individual_Analysis}
\title{Individual-variant analysis using score test}
\usage{
Individual_Analysis(
  chr,
  start_loc = NULL,
  end_loc = NULL,
  individual_results = NULL,
  genofile,
  obj_nullmodel,
  mac_cutoff = 20,
  subset_variants_num = 5000,
  QC_label = "annotation/filter",
  variant_type = c("variant", "SNV", "Indel"),
  geno_missing_imputation = c("mean", "minor"),
  tol = .Machine$double.eps^0.25,
  max_iter = 1000,
  SPA_p_filter = TRUE,
  p_filter_cutoff = 0.05,
  use_ancestry_informed = FALSE,
  find_weight = FALSE
)
}
\arguments{
\item{chr}{chromosome.}

\item{start_loc}{starting location (position) of the genetic region for each individual variant to be analyzed using score test.}

\item{end_loc}{ending location (position) of the genetic region for each individual variant to be analyzed using score test.}

\item{individual_results}{the data frame of (significant) individual variants of interest for ancestry-informed analysis.
The first 4 columns should correspond to chromosome (CHR), position (POS), reference allele (REF), and alternative allele (ALT).}

\item{genofile}{an object of opened annotated GDS (aGDS) file.}

\item{obj_nullmodel}{an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function,
or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{\link{genesis2staar_nullmodel}} function.}

\item{mac_cutoff}{the cutoff of minimum minor allele count in
defining individual variants (default = 20).}

\item{subset_variants_num}{the number of variants to run per subset for each time (default = 5e3).}

\item{QC_label}{channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").}

\item{variant_type}{type of variant included in the analysis. Choices include "variant", "SNV", or "Indel" (default = "variant").}

\item{geno_missing_imputation}{method of handling missing genotypes. Either "mean" or "minor" (default = "mean").}

\item{tol}{a positive number specifying tolerance, the difference threshold for parameter
estimates in saddlepoint approximation algorithm below which iterations should be stopped (default = ".Machine$double.eps^0.25").}

\item{max_iter}{a positive integer specifying the maximum number of iterations for applying the saddlepoint approximation algorithm (default = "1000").}

\item{SPA_p_filter}{logical: are only the variants with a score-test-based p-value smaller than a pre-specified threshold use the SPA method to recalculate the p-value, only used for imbalanced case-control setting (default = TRUE).}

\item{p_filter_cutoff}{threshold for the p-value recalculation using the SPA method, only used for imbalanced case-control setting (default = 0.05)}

\item{use_ancestry_informed}{logical: is ancestry-informed association analysis used to estimate p-values (default = FALSE).}

\item{find_weight}{logical: should the ancestry group-specific weights and weighting scenario-specific p-values for each base test be saved as output (default = FALSE).}
}
\value{
A data frame containing the score test p-value and the estimated effect size of the minor allele for each individual variant in the given genetic region, or as provided in \code{individual_results}
for ancestry-informed variant analysis. The first 4 columns correspond to chromosome (CHR), position (POS), reference allele (REF), and alternative allele (ALT).
If \code{find_weight} is TRUE, returns a list containing the ancestry-informed score test p-values and the estimated effect size of the minor allele for each individual variant provided in \code{individual_results}.
The ensemble weights under two sampling scenarios and p-values under scenarios 1, 2, and combined for each base test are saved as well.
}
\description{
The \code{Individual_Analysis} function takes in chromosome, starting location, ending location, an user-defined variant list for
ancestry-informed analyses, the object of opened annotated GDS file, and the object from fitting the null model to analyze the association between a
quantitative/dichotomous phenotype (including imbalanced case-control design) and each individual variant in a genetic region by using score test.
For multiple phenotype analysis (\code{obj_nullmodel$n.pheno > 1}),
the results correspond to multi-trait score test p-values by leveraging
the correlation structure between multiple phenotypes.
For ancestry-informed analysis, the results correspond to ensemble p-values across base tests,
with the option to return a list of base weights and p-values for each base test.
}
\references{
Chen, H., et al. (2016). Control for population structure and relatedness for binary traits
in genetic association studies via logistic mixed models. \emph{The American Journal of Human Genetics}, \emph{98}(4), 653-666.
(\href{https://doi.org/10.1016/j.ajhg.2016.02.012}{pub})

Li, Z., Li, X., et al. (2022). A framework for detecting
noncoding rare-variant associations of large-scale whole-genome sequencing
studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
(\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
}
