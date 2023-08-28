\name{test.nested.model}
\alias{test.nested.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test for nested ANOVA models}
\description{
The function compute the F-statistic for nested ANOVA model.
}
\usage{
test.nested.model(X, X0, Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{this is the design matrix corresponding to the full model. X
    must be of full rank.}
  \item{X0}{this is the design matrix corresponding to the sub-model. X0
    is such that span(X0) is a vectorial subspace of span(X). X0
    must be of full rank.}
  \item{Y}{ this is a matrix corresponding to the dependent variable
    (samples in row, variables (e.g. genes) in column.}
}
\details{
  The F-statistic  is computed and the p-value is returned. The
  F-statistic has r-r0 and n-r degrees of freedom where r is the rank of
  matrix X, r0 is the rand of matrix X0 and n is the number of observations.
}
\value{
 A list with the following items:
  \item{theta }{The estimation of the paratemer of the full model.}
  \item{F }{Value of the F-statistic.}
  \item{pvalue}{The p-value.}
  \item{residuals}{The residuals for the full model.}
  \item{sigma2}{The estimation of the variance for the full model.}
  \item{X}{The matrix design X}
  ...
}
% \references{ ~put references to the literature/web site here ~ }
\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}
% \note{ ~~further notes~~ 

%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
Y <- matrix(rnorm(21), 21, 1)
ef <- gl(3,7)
X <- lm(Y ~ ef, x = TRUE)$x
X0 <- lm(Y ~ 1, x = TRUE)$x

res.test <- test.nested.model(X, X0, Y)

}
