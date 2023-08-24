\name{runMFA}
\alias{runMFA}
\title{Function to perform a Multiple Factor Analysis.}
\description{
  This function allows to perform a Multiple Factor Analysis and to
  build a report with the main statistics and the main graphics.
}
\usage{
runMFA(Data, group = NULL, ncp = 5, name.group = NULL, type = NULL,
ind.sup = NULL, num.group.sup = NULL, graph = TRUE,
report.file = NULL, report.pdf = NULL)
}
\arguments{
  \item{Data}{A data.frame or a list, see details}
  \item{group}{A vector indicating the number of variables in each group}
  \item{ncp}{Number of components to keep}
  \item{name.group}{Names of the groups of variables, default is NULL}
  \item{type}{the type of variables in each group; three possibilities:
    "c" or "s" for quantitative variables (the difference is that for
    "s" variables are scaled to unit variance), "n" for qualitative
    variables; by default, all variables are quantitative and scaled to
    unit variance}
  \item{ind.sup}{A vector indicating the indexes of the supplementary individuals}
  \item{num.group.sup}{The indexes of the illustrative groups (by
    default, NULL and no group are illustrative)}
  \item{graph}{boolean, if TRUE a graph is displayed}
  \item{report.file}{Name of the txt file for the text report}
  \item{report.pdf}{Name of the pdf file for the graphical report}
}
\details{
  Data can either be a \code{data.frame} with all the groups binded or a
  list containing all the groups such \code{list(group1=group1,
    group2=group2, ...)}. In the case of a \code{data.frame},
  \code{group} must be indicated.

  The analysis report can also be build later by using functun \code{\link{MFAreport}}.
}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{\code{\link{MFA}}, \code{\link{plot.MFA}}, \code{\link{MFAreport}}}
\examples{
data(marty)
## Perform a MFA on splitted data
resMFA <- runMFA(Data=list(group1=t(marty[1:100,]),
group2=t(marty[101:200,])), type=c("c", "c"), graph=FALSE)
\dontrun{
## plot global analysis with partial individuals
plot(resMFA, choix="ind", partial="all")

## plot groups link
plot(resMFA, choix="group")
}

}
\keyword{multivariate}
