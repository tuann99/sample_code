\name{PLS}
\Rdversion{1.1}
\alias{PLS}
\title{
Partial Least Squares
}
\description{
Partial Least Squares
}
\usage{
PLS(E, F, n = 1, scale = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{E}{
    Matrix which correspond to the predictor variables.
%%     ~~Describe \code{E} here~~
}
\item{F}{
  Matrix which corresponds to the outcome variables.
%%     ~~Describe \code{F} here~~
}
  \item{n}{
Number of components to compute.
}
  \item{scale}{
    If TRUE, each variable is sacled to have a unit variance.
    %%     ~~Describe \code{scale} here~~
}
}
\details{
The function implements the PLS1 algorithm
(Wold, 1966; Wold et al., 1982) (details about this algorithm and other
implementations can be found in Tenenhaus (1998)).



%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list whith the following elements:
\item{T}{PLS components. T = EW(P'W)^(-1)}
\item{E0}{The input E matrix. If scale=TRUE, then E is scaled.}
\item{F0}{The input F matrix. If scale=TRUE, then F is scaled.}
\item{C}{Covariance between the PLS componants and the outcome variables
  F.}
\item{Wstar}{EW(P'W)^(-1)}
\item{P}{Regression vectors.}
\item{param.center}{Means of E columns.}
\item{param.scaled}{Standard-deviation of E colums.}
%%  \item{comp2 }{Description of 'comp2'}

}
\references{
                        
Tenenhaus, M. (1998).   La regression PLS - Theorie et
pratique. Editions TECHNIP.

Wold, H. (1966). Estimation of principal component and related models by
 iterative least squares. In Krishnaiah, P. R., editor, Multivariate
 Analysis, pages 391-420. New-York: Academic Press.

                                             
Wold, S., Martens, H., and Wold, H. (1982). The multivariate calibration problem in chemistry solved by the PLS method. In Ruhe, A. and Kastrom, B., editors, Matrix Pencils, Lecture Notes in Mathematics, pages 286-293. Springer Berlin / Heidelberg.
}

\examples{
X <- matrix(rnorm(100), 20, 5)
Y <- matrix(c(rep(0,10), rep(1,10)))
res <- PLS(X,Y)
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
