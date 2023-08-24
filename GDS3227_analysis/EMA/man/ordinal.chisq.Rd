\name{ordinal.chisq}
\alias{ordinal.chisq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Chisq test for ordinal values}
\description{
  This function computes a chisq test for ordinal values
}
\usage{
ordinal.chisq(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a contingency table with ordinal values in column}
}

\details{
This function applies a ordinal chisq test as described in
http://www.uvm.edu/~dhowell/StatPages/More_Stuff/OrdinalChisq/OrdinalChiSq.html.
The results are identical to those returned by SPSS

}
\value{
  Chisq statistics and pvalues
}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}
\examples{
##General hypothesis : was that participants who had experienced more traumatic events
## during childhood would be more likely to drop out of treatment.
trau=matrix(c(25,13,9,10,6,31,21,6,2,3), byrow=TRUE, nrow=2)
colnames(trau)<-c("0","1","2","3","4+")
rownames(trau)<-c("dropout","remain")
ordinal.chisq(trau)

##Association between grade and local breast cancer relapse
grade <- matrix(c(16,42,71,4,27,49), ncol=3, byrow=TRUE)
colnames(grade)<-c("low", "intermediate","high")
rownames(grade)<-c("0", "1")

ordinal.chisq(grade)

}
\keyword{utilities}

