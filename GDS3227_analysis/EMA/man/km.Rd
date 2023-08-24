\name{km}
\alias{km}
\title{
Compute survival curves and test difference between the curves
}
\description{
Computes and draw estimates of survival curves for censored data using Kaplan-Meier's method.
In case of >=2 survival curves, use logrank test to assess the difference between the survival curves.
Missing values are removed
}
\usage{
km(time, status, group = NULL, xlab="Time (years)", ylab="", ...)
}
\arguments{
  \item{time}{
numeric, this is the follow up time (used with right censored data)
}
  \item{status}{
The status indicator, normally 0=alive, 1=dead (numeric). Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death)
}
  \item{group}{
indicates the group to which is assigned each observation (factor).
For one groupe only (no comparison), group is set to null (default)
}
  \item{xlab}{
(optional): a character string, xlabel of the Kaplan Meier's plot
}
  \item{ylab}{
(optional): a character string, ylabel of the Kaplan Meier's plot
}
\item{...}{
(optional): Additional graphical parameters
}
}

\value{
A list with
\item{fit.km}{results provided by Kaplan Meier analysis. See the R help on survfit for details}
\item{lr}{results provided by logrank analysis. See the R help on survdiff for details}
\item{p.lr}{pvalue of the logrank test}
}
\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{
\code{\link{survfit}}, \code{\link{survdiff}}
}
\examples{
require(survival)
data(leukemia)
time<-leukemia$time
status<-leukemia$status
# One group
res<-km(time,status)
res$fit.km

# Two groups
group<-leukemia$x
res<-km(time,status,group)
res<-km(time,status,group,title="Kaplan Meier curve")
res<-km(time,status,group,title="Kaplan Meier curve",pdfname="My survival curve")
res<-km(time,status,group,pdfname="My survival curve",pdfwidth=11.69,pdfheight=8.27)
res$fit.km
names(res$fit.km)
res$lr
names(res$lr)
res$p.lr
}
\keyword{survival}


