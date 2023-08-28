\name{foldchange}
\alias{foldchange}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute foldchange}
\description{
 Caculates the foldchange for each genes.
}
\usage{
foldchange(data, labels, unlog = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A matrix or a data frame of expression data. Each row of 'data' must correspond to a gene, and each column to a sample.}
  \item{labels}{A vector of length 'ncol(data)' containing the class
    labels of the samples.
    'labels' should be a numeric vector 0 / 1.
    0 specifying the samples of, e.g., the control group
    and 1 specifying, e.g., the case group.}
  \item{unlog}{'TRUE' if the data have to be unlog before the foldchange calculation.}
}

\value{
 A vector of length 'n genes' with the foldchange for each genes.
}
\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\examples{
## load data
data(marty)

\dontrun{
## filtering data
marty <- expFilter(marty, threshold=3.5, graph=FALSE)
}

##Class label 0/1
marty.type.num <- ifelse(marty.type.cl=="Her2+",0,1)

## run folchange on 50 genes
fcOUT <- foldchange(marty[1:50,], marty.type.num)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
