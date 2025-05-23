% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MultiSTAAR_sp.R
\name{MultiSTAAR_sp}
\alias{MultiSTAAR_sp}
\title{Multi-trait STAAR (MultiSTAAR_sp) procedure using omnibus test (dgCMatrix version)}
\usage{
MultiSTAAR_sp(
  genotype_sp,
  MAF = NULL,
  obj_nullmodel,
  annotation_phred = NULL,
  rare_maf_cutoff = 0.01,
  rv_num_cutoff = 2,
  rv_num_cutoff_max = 1e+09
)
}
\arguments{
\item{genotype_sp}{a sparse genotype matrix (dgCMatrix) of dimension n*p,
where n is the sample size and p is the number of genetic variants.
The matrix should be extracted using \code{STAARpipeline::Genotype_sp_extraction},
and it has been flipped: only the minor allele is stored (coded as 1), while the major allele is not stored (coded as 0).}

\item{MAF}{a numeric vector of minor allele frequencies for the variants in \code{genotype_sp}.
It can be computed using \code{Matrix::colMeans(genotype_sp,na.rm = TRUE)/2}.}

\item{obj_nullmodel}{an object from fitting the null model, which is the
output from \code{\link{fit_null_glmmkin_multi}} function for unrelated or
related samples. Note that \code{\link{fit_null_glmmkin_multi}}
is a wrapper of the \code{\link{glmmkin}} function from the \code{\link{GMMAT}} package.}

\item{annotation_phred}{a data frame or matrix of functional annotation data
of dimension p*q (or a vector of a single annotation score with length p).
Continuous scores should be given in PHRED score scale, where the PHRED score
of j-th variant is defined to be -10*log10(rank(-score_j)/total) across the genome. (Binary)
categorical scores should be taking values 0 or 1, where 1 is functional and 0 is
non-functional. If not provided, MultiSTAAR will perform the multi-trait
SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1), ACAT-V(1,25), ACAT-V(1,1)
and ACAT-O tests (default = NULL).}

\item{rare_maf_cutoff}{the cutoff of maximum minor allele frequency in
defining rare variants (default = 0.01).}

\item{rv_num_cutoff}{the cutoff of minimum number of variants of analyzing
a given variant-set (default = 2).}

\item{rv_num_cutoff_max}{the cutoff of maximum number of variants of analyzing
a given variant-set (default = 1e+09).}
}
\value{
a list with the following members:

\code{num_variant}: the number of variants with minor allele frequency > 0 and less than
\code{rare_maf_cutoff} in the given variant-set that are used for performing the
variant-set using MultiSTAAR.

\code{cMAC}: the cumulative minor allele count of variants with
minor allele frequency > 0 and less than \code{rare_maf_cutoff} in the given variant-set.

\code{RV_label}: the boolean vector indicating whether each variant in the given
variant-set has minor allele frequency > 0 and less than \code{rare_maf_cutoff}.

\code{results_STAAR_O}: the multi-trait STAAR-O (MultiSTAAR-O) p-value
that aggregated multi-trait SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1),
ACAT-V(1,25), and ACAT-V(1,1) together with p-values of each test weighted by
each annotation using Cauchy method.

\code{results_ACAT_O}: the multi-trait ACAT-O p-value that aggregated multi-trait SKAT(1,25),
SKAT(1,1), Burden(1,25), Burden(1,1), ACAT-V(1,25), and ACAT-V(1,1) using Cauchy method.

\code{results_STAAR_S_1_25}: a vector of multi-trait STAAR-S(1,25) p-values,
including multi-trait SKAT(1,25) p-value weighted by MAF, the multi-trait SKAT(1,25)
p-values weighted by each annotation, and a multi-trait STAAR-S(1,25)
p-value by aggregating these p-values using Cauchy method.

\code{results_STAAR_S_1_1}: a vector of multi-trait STAAR-S(1,1) p-values,
including multi-trait SKAT(1,1) p-value weighted by MAF, the multi-trait SKAT(1,1)
p-values weighted by each annotation, and a multi-trait STAAR-S(1,1)
p-value by aggregating these p-values using Cauchy method.

\code{results_STAAR_B_1_25}: a vector of multi-trait STAAR-B(1,25) p-values,
including multi-trait Burden(1,25) p-value weighted by MAF, the multi-trait Burden(1,25)
p-values weighted by each annotation, and a multi-trait STAAR-B(1,25)
p-value by aggregating these p-values using Cauchy method.

\code{results_STAAR_B_1_1}: a vector of multi-trait STAAR-B(1,1) p-values,
including multi-trait Burden(1,1) p-value weighted by MAF, the multi-trait Burden(1,1)
p-values weighted by each annotation, and a multi-trait STAAR-B(1,1)
p-value by aggregating these p-values using Cauchy method.

\code{results_STAAR_A_1_25}: a vector of multi-trait STAAR-A(1,25) p-values,
including multi-trait ACAT-V(1,25) p-value weighted by MAF, the multi-trait ACAT-V(1,25)
p-values weighted by each annotation, and a multi-trait STAAR-A(1,25)
p-value by aggregating these p-values using Cauchy method.

\code{results_STAAR_A_1_1}: a vector of multi-trait STAAR-A(1,1) p-values,
including multi-trait ACAT-V(1,1) p-value weighted by MAF, the multi-trait ACAT-V(1,1)
p-values weighted by each annotation, and a multi-trait STAAR-A(1,1)
p-value by aggregating these p-values using Cauchy method.
}
\description{
The \code{MultiSTAAR_sp} function is a sparse matrix (dgcMatrix) version of \code{MultiSTAAR}.
It takes in a sparse genotype matrix, the object from fitting the null
model, and functional annotation data to analyze the association between multiple
(quantitative) phenotypes and a variant-set by using MultiSTAAR procedure.
For each variant-set, the multi-trait STAAR-O (MultiSTAAR-O) p-value is a p-value
from an omnibus test that aggregated multi-trait SKAT(1,25), SKAT(1,1),
Burden(1,25), Burden(1,1), ACAT-V(1,25), and ACAT-V(1,1) together with p-values
of each test weighted by each annotation using Cauchy method.
}
\references{
Li, X., Chen, H., et al. (2023). A statistical framework for powerful
multi-trait rare variant analysis in large-scale whole-genome sequencing studies.
\emph{bioRxiv}. (\href{https://doi.org/10.1101/2023.10.30.564764}{preprint})

Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
in silico functional annotations empowers rare variant association analysis of
large whole-genome sequencing studies at scale. \emph{Nature Genetics}, \emph{52}(9), 969-983.
(\href{https://doi.org/10.1038/s41588-020-0676-4}{pub})

Li, Z., Li, X., et al. (2022). A framework for detecting
noncoding rare-variant associations of large-scale whole-genome sequencing
studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
(\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})

Liu, Y., et al. (2019). Acat: A fast and powerful p value combination
method for rare-variant analysis in sequencing studies.
\emph{The American Journal of Human Genetics}, \emph{104}(3), 410-421.
(\href{https://doi.org/10.1016/j.ajhg.2019.01.002}{pub})

Li, Z., Li, X., et al. (2020). Dynamic scan procedure for
detecting rare-variant association regions in whole-genome sequencing studies.
\emph{The American Journal of Human Genetics}, \emph{104}(5), 802-814.
(\href{https://doi.org/10.1016/j.ajhg.2019.03.002}{pub})
}
