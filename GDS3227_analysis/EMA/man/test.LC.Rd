\name{test.LC}
\alias{test.LC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test linear combinations of parameters of a linear model}
\description{
This function performs either a global significance F-test of one or several linear combinations of the parameters (contrast) of a (generalized) linear model, or tests for the significance of each contrast simultaneously, with or without multiple test correction.}
\usage{
test.LC(C, X, Y, global=FALSE, cor.multtest=TRUE,typeFDR="FDR-BH")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{C}{Vector or matrix of the linear combinations of factors we want to test.}
  \item{X}{Design matrix of the (generalized) linear model adjusted to the data.}
  \item{Y}{Vector or response matrix (e.g. a gene expression matrix) to which the (generalized) linear model is applied. If a matrix is specified the test(s) will be done on each column of the matrix.}
  \item{global}{logical indicating if one want to perform a F test of the global hypothesis H_0=0, or to test each contrast simultaneously }
  \item{cor.multtest}{logical.If global=FALSE indicates if one want to apply a multiple test correction for the computation of the p-values}
  \item{typeFDR}{If global=FALSE and cor.multtest=TRUE this argument is passed to the function multiple.correction(), it specifies which correction method to apply. See '?multiple.correction' for more details.}
}
\details{
The design matrix X can be extracted from a (generalized) linear model with the function \code{\link{model.matrix}}.
If Y is a response matrix each column shall be a gene/individual response on which the linear combinations will be tested.
If one wish to test several linear combinations at the same time, C must be a matrix p*n with n the  number of columns in the design matrix X (i.e. the number of parameters in the (generalized) linear model) and p the number of linear combinations, so that each row stands for one linear combination. Else C shall be a vector of length n.
In case where both C and Y are matrices, each linear combination of the parameters (i.e. each row of C) will be tested for nullity on each gene/individual (i.e. each column of Y). 
}
\value{
A list of matrices :
  \item{Estimate }{A matrix containing the estimated values: F[i,j] is the estimated value of the linear combination i (on gene j). }
  \item{F }{A matrix containing the F values: F[i,j] is the F value of the test of linear combination i on gene j. }
  \item{pvalue }{A matrix containing the p-values: pvalue[i,j] is the p-value of the test of linear combination i on gene j. }
  \item{Y.pred }{A matrix containing the predicted response vectors}
  \item{resid }{A matrix containing the residuals of the linear model}
  \item{sigma2 }{A matrix with one row containing the residual variance for each gene}
  \item{theta }{A matrix containing the estimates of the effects of the original model. These estimates are the ones also obtained with summary.lm()}
}
\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}
\seealso{\code{\link{model.matrix}},\code{\link{lm}}}
\examples{
data(marty)

##Class label 0/1
marty.type.num <- ifelse(marty.type.cl=="Her2+",0,1)

#Annotation of the grade of tumor
grade=factor(sample(c(1:3),23,rep=TRUE),labels=c("I","II","III"))
typeXfactor=paste(marty.type.num,grade,sep=":")

##dataset fro example
normData<-marty[1:100,]

rt<-runTtest(normData, labels=marty.type.num, plot=FALSE)

normData.DE<-normData[which(rt["AdjpValue"]<0.05),]

marty.lm=lm(t(normData.DE)~as.factor(typeXfactor))
X=model.matrix(marty.lm)

#We want to test Basal vs Her2+ within each grade 
LC<-matrix(c(0,0,0,-1,0,0,0,1,0,0,-1,0,0,0,1,0,0,-1),ncol=6,byrow=TRUE)
#We also want to test grade II vs grade III 
LC2=c(0,1,-1,0,1,-1)
LC=rbind(LC,LC2)
row.names(LC)=c("B:I-H:I","B:II-H:II","B:III-H:III","II-III")
marty.LC=test.LC(LC,X,t(normData.DE))
marty.LC$pvalue

#List of the probesets differentially expressed for each of the four tests :
ll=list()
for(i in 1:nrow(marty.LC$pvalue)){
	ll[[i]]=as.matrix(marty.LC$pvalue[i,which(marty.LC$pvalue[i,]<0.05)])
	rownames(ll[[i]])=names(which(marty.LC$pvalue[i,]<0.05))
}
names(ll)=rownames(marty.LC$pvalue)
print(ll)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{models}
