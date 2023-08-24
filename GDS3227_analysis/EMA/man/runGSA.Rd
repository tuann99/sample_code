\name{runGSA}
\alias{runGSA}


\title{GSA analysis}

\description{
     Determines the significance of  pre-defined sets of genes  with
     respect to an outcome variable, such as a group indicator, a
     quantitative variable or a survival time
}

\usage{
runGSA(nData, labels, gmtfile, chip = "hgu133plus2", np = 1000,
 minsize = 10, maxsize = 800, resp.type = "Two class unpaired",
fdr = 0.25)
}


\arguments{
  \item{nData}{a matrix or a data frame of expression data.
    Each row of 'data' must correspond to a gene, and each column to a
    sample.
  }
  \item{labels}{ a vector of length 'ncol(data)' containing the class labels
    of the samples.
    In "Two class unpaired", 'labels' should be a vector containing 0's
    (specifying the samples of, e.g., the control group) and 1's
    (specifying, e.g., the case group) or more if multiclass (0,1,2...)
    In "Two class paired" for paired outcomes, coded -1,1 (first pair),
    -2,2 (second pair), etc.. 
  }
  \item{gmtfile}{ a character string corresponding to the path file of a
    gmt file, corresponding to a gene set collection (a list).
  }
  \item{chip}{ a character string corresponding to the chip type of the
    data.
  }
  \item{np}{ a numerical value corresponding to the number of
    permutations.
  }
  \item{minsize}{ a numerical value corresponding to the minimum number
    of genes in  genesets to be considered.
  }
  \item{maxsize}{ a numerical value corresponding to the minimum number
    of genes in  genesets to be considered.
  }
  \item{resp.type}{
    Problem type: "quantitative" for a continuous parameter;
    "Two class unpaired" ; "Survival" for censored survival
    outcome; "Multiclass" : more than 2 groups.
    "Two class paired" for paired outcomes.
  }
  \item{fdr}{ a numerical value corresponding to the fdr threshold.
  }

}


\details{
The GSA package is presented as an improvement of the GSEA approach.
It differs from a GSEA in its use of the "maxmean" statistic: this is the mean of the positive or negative part of gene scores in the gene set, whichever is large in absolute values.\\
Efron and Tibshirani shows that this is often more powerful than the modified KS statistic used in GSEA. GSA also does "restandardization" of the genes (rows), on top of the
permutation of columns (done in GSEA).
  }

\value{
  A list of geneset found
  If it is a LIST, use
  \item{FDRcut}{ a numerical value corresponding to the threshold FDR.
  }
  \item{negative}{ a character matrix corresponding to the downexpressed
  gene sets found.
  }
  \item{positive}{ a character matrix corresponding to the upexpressed
    gene sets found.
  }
  \item{nsets.neg}{ a numerical value corresponding to the number of
    downexpressed gene sets found.
  }
  \item{nsets.pos}{ a numerical value corresponding to the number of
    upexpressed gene sets found.
  }
}

\references{
  Efron, B. and   Tibshirani, R.   On testing the significance of
  sets of genes. Stanford tech report rep 2006.
  http://www-stat.stanford.edu/~tibs/ftp/GSA.pdf
  
  Subramanian, A. and Tamayo, P. Mootha, V. K. and Mukherjee, S. and
  Ebert, B. L. and Gillette, M. A. and Paulovich, A. and Pomeroy, S.
  L. and Golub, T. R. and Lander, E. S. and Mesirov, J. P. (2005)  A
  knowledge-based approach for interpreting genome-wide expression
  profiles. PNAS. 102, pg 15545-15550.
}
\author{ Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe }

\seealso{ \code{\link{GSA}}}

\examples{
\dontrun{
require(hgu133plus2.db)

## Two class unpaired comparison
## load data
data(marty)

## filtering data
marty <- expFilter(marty, threshold=3.5, graph=FALSE)

##Class label 0/1
marty.type.num <- ifelse(marty.type.cl=="Her2+",0,1)

## run sam analysis
gsaOUT <- runGSA(marty, marty.type.num ,
   gmtfile="./c2.kegg.v2.5.symbols.gmt", chip="hgu133plus2")
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
