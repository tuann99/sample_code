\name{runIndTest}
\alias{runIndTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Computing Differential Analysis for each gene}
\description{
  This function computes test statistics, e.g., two-sample Welch
  t-statistics, t-statistics, or wilcoxon, independently for each row of
  a data frame.
}
\usage{
runIndTest(data, labels, gene.names = NULL, plot = TRUE, dirname= NULL,
grp.name=c("Group1","Group2")) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a matrix, a data frame, or an ExpressionSet object. Each row of 'data' (or 'exprs(data)', respectively) must correspond to a gene, and each column to a sample.}
  \item{labels}{A vector of integers corresponding to observation (column) class labels. For 2 classes, the labels must be 0 and 1.}
  \item{gene.names}{A vector of description or name for each gene.}
  \item{plot}{A logical value specifying if drawing plots or not.}
  \item{dirname}{If specified, the .png plots are created in the
    directory.}
  \item{grp.name}{Vector with the name of the two groups}
}
\details{
 For each gene independently, the function tests for the normality
 (Shapiro test) and the variance equality (F test) of each
 groups. According to the results, a welch test, a student test or a
 wilcoxon test is performed.
}
\value{
  A matrix with the gene names, the statistics, and the p-values.

}
\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{\code{\link{shapiro.test}}, \code{\link{var.test}},\code{\link{t.test}},\code{\link{wilcox.test}}}
\examples{
## load data
data(marty)

##random choice of genes - in practice genes of interest
geneOfInterest<-sample(1:ncol(marty),5)

##Class label 0/1
marty.type.num <- ifelse(marty.type.cl=="Her2+",0,1)

## run differential analysis
out <- runIndTest(marty[geneOfInterest,], labels=marty.type.num)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
