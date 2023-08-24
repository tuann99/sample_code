\name{runWilcox}
\alias{runWilcox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computing Multiple Wilcoxon Tests }
\description{
    This function provides a convenient way to compute the wilcoxon
    statistics, for each row of a data frame using the multtest
  package. It returns the raw and adjusted pvalues for each genes as
  well as the significance of the genes and a quantile-quantile plot.
 
}
\usage{
runWilcox(data, labels, typeFDR = "FDR-BH", q = 0.05, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A matrix, a data frame, or an ExpressionSet object. Each
    row of 'data' (or 'exprs(data)', respectively) must correspond to a gene, and each column to a sample.}
  \item{labels}{A vector of integers corresponding to observation (column) class labels. For 2 classes, the labels must be 0 and 1.}
  \item{typeFDR}{The method to apply fo the multiple testing correction.}
  \item{q}{A numeric value specifying the pvalue threshold.}
  \item{plot}{A logical value specifying if drawing plots or not.}
}

\value{
    A matrix with the probes ID, the statistics, the raw p-value and the adjust p-value
}
\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{\code{\link[multtest]{mt.teststat}},\code{\link{multiple.correction}}}
\examples{
## load data
data(marty)

## filtering data
##marty <- expFilter(marty, threshold=3.5, graph=FALSE)

##Class label 0/1
marty.type.num <- ifelse(marty.type.cl=="Her2+",0,1)

## run differential analysis
example.subset<-marty[1:100,]
out <- runWilcox(example.subset, labels=marty.type.num, typeFDR="FDR-BH", plot=FALSE)

}
