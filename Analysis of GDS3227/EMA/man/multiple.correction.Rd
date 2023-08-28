\name{multiple.correction}
\alias{multiple.correction}
\title{Multiple testing correction }
\description{
  Given a set of p-values, returns p-values adjusted using one of several methods.
}
\usage{
multiple.correction(pval, typeFDR, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pval}{Vector of pvalues.}
  \item{typeFDR}{The correction method.}
  \item{q}{The error rate to use for the Two-stages procedure (FDR-TST).}
}
\details{
  The multiple correction methods include Bonferroni correction
  ("FWER"), Benjamini-Hochberg standard false discovery rate correction
  ("FDR-BH"), Benjamini-Hochberg Adaptive Procedure ("FDR-TST") and the
  Qvalue procedure (Storey).
}
\value{
 A vector of adjusted pvalues.
}
\references{
- Benjamini Y and Hochberg Y. (1995) Controlling the false discovery rate: A practical
and powerful approach to multiple testing. . Journal of the Royal Statistical Society,
Series B, 57: 289-300.
- Storey JD. (2002) A direct approach to false discovery rates. Journal of the Royal
Statistical Society, Series B, 64: 479-498.
- Benjamini Y., Kenigsberg E., Reiner A., Yekutieli D (2005). FDR
adjustments of Microarray Experiments.
}
\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{\code{\link{p.adjust}}}
\examples{
\dontrun{
## load data
data(marty)

## filtering data
marty <- expFilter(marty, threshold=3.5, graph=FALSE)


##Class label 0/1
marty.type.num <- ifelse(marty.type.cl=="Her2+",0,1)

##Example dataset
example.subset<-marty[1:100,]

## run differential analysis Basal vs HER2+
out <- runTtest(example.subset, labels=marty.type.num, typeFDR="FDR-BH")
## OR :
out2 <- multiple.correction(out$RawpValue, typeFDR="FDR-BH")
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
