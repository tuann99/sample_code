\name{runSAM}
\alias{runSAM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ SAM analysis with siggenes package }
\description{
     Performs a Significance Analysis of Microarrays (SAM). It is
     possible to perform one and two class analyses using either a
     modified t-statistic or a (standardized) Wilcoxon rank statistic,
     and a multiclass analysis using a modified F-statistic. Moreover,
     this function provides a SAM procedure for categorical data such
     as SNP data and the possibility to employ an user-written score
     function.
}
\usage{
runSAM(data, labels, nbpermut = 500, q = 0.05, plot = TRUE, method
="d.stat",var.equal = TRUE, include.zero = FALSE, paired = FALSE,
seed=123)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A matrix, a data frame, or an ExpressionSet object. Each
    row of 'data' (or 'exprs(data)', respectively) must correspond to a gene, and each column to a sample.}
  \item{labels}{A vector of length 'ncol(data)' containing the class labels
    of the samples.
    In the two class unpaired case, 'labels' should be a vector
    containing 0's (specifying the samples of, e.g., the control
    group) and 1's (specifying, e.g., the case group). 
    
    In the two class paired case, 'labels' can be either a numeric
    vector or a numeric matrix.  If it is a vector, then 'labels' has
    to consist of the integers between -1 and  -n/2 (e.g., before
    treatment group) and between 1 and n/2 (e.g., after treatment
    group), where n is the length of 'labels' and k is paired with
    -k, k=1,...,n/2. If 'labels' is a matrix, one column should
    contain -1's and 1's specifying, e.g., the before and the
    after treatment samples, respectively, and the other column
    should contain integer between 1 and n/2 specifying the n/2
    pairs of observations.
  }
  \item{nbpermut}{ A numeric value specifying the number of permutation. }
  \item{q}{ A numeric value specifying the FDR threshold. see details. }
  \item{plot}{ A logical value specifying if drawing plots or not. }
  \item{method}{A character string or a name specifying the method/function
          that should be used in the computation of the expression
          scores d. 

          If 'method = d.stat', a modified t-statistic or F-statistic,
          respectively, will be computed as proposed by Tusher et al.
          (2001). 

          If 'method = wilc.stat', a Wilcoxon rank sum statistic or
          Wilcoxon signed rank statistic will be used as expression
          score. 

          For an analysis of categorical data such as SNP data, 
          'method' can be set to 'chisq.stat'. In this case Pearson's
          ChiSquare statistic is computed for each row. 

          If the variables are ordinal and a trend test should be
          applied  (e.g., in the two-class case, the Cochran-Armitage
          trend test), 'method = trend.stat' can be employed.
}
  \item{var.equal}{ A logical value. If 'method=d.stat', TRUE for
    student test , FALSE for Welch test. }
  \item{include.zero}{ A numeric value specifying if sO=0 is possible. }
  \item{paired}{ A logical value specifying if paired test or not. }
  \item{seed}{ Seed initialization for results reproducibility.}
}
\details{
  SAM has it own FDR procedure which allows to find significant genes for
  a fixed threshold 'q'.
  The genes' signicance found by SAM is not based on the adjusted
          pvalues (qvalues). That's why we do not report them.
  %If we look for the significant genes using the adjusted pvalues (<q),the size of the genes list is very different (bigger in most cases).
  }

\value{
  A matrix with the probes ID, the statistics, the raw p-values, and the significance (according to SAM FDR procedure).
}
\references{
  Schwender, H., Krause, A. and Ickstadt, K. (2003). Comparison of
  the Empirical Bayes and the Significance Analysis of Microarrays.
  _Technical Report_, SFB 475, University of Dortmund, Germany.
  
  Schwender, H. (2004). Modifying Microarray Analysis Methods for 
  Categorical Data - SAM and PAM for SNPs. To appear in:
  _Proceedings of the the 28th Annual Conference of the GfKl_.
  
  Tusher, V.G., Tibshirani, R., and Chu, G. (2001). Significance
  analysis of microarrays applied to the ionizing radiation
  response. _PNAS_, 98, 5116-5121.
}
\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{ \code{\link{sam}} }

\examples{
## load data
data(marty)

\dontrun{
## filtering data
marty <- expFilter(marty, threshold=3.5, graph=FALSE)
}

##Class label 0/1
marty.type.num <- ifelse(marty.type.cl=="Her2+",0,1)

## run sam analysis on example set
example.subset <- marty[1:100,]
samOUT <- runSAM(example.subset, marty.type.num, nbpermut=50, q=0.05, plot=TRUE)
samSIGN <- samOUT[which(samOUT[,"Significant"]),]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
\keyword{htest}
