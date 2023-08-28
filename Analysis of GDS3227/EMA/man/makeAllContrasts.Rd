\name{makeAllContrasts}
\alias{makeAllContrasts}
\title{Create all pairwise contrasts}
\description{
  Create the matrix of all pairwise contrasts between parameters for the
  test.LC function
}
\usage{
makeAllContrasts(X, annot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{the design matrix}
  \item{annot}{the annotation table of the individuals with the columns
    corresponding to the model variables}
  
}
\details{
  This function create a contrasts matrix that can be used with the
  \code{test.LC} function. All pairwise comparisons for each
  variable stratified into all other variables are returned.
}


\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\examples{
data(marty)

##Class label 0/1
marty.type.num <- ifelse(marty.type.cl=="Her2+",0,1)

#Annotation of the grade of tumor
grade <- factor(sample(c(1:3),23,rep=TRUE),labels=c("I","II","III"))

annot <- data.frame(type=marty.type.num, grade=grade)
rownames(annot) <- colnames(marty)

marty.lm <- lm(marty[1,]~ annot$grade * annot$type)
X <- model.matrix(marty.lm)

LC <- makeAllContrasts(X, annot)

}

\keyword{models}

