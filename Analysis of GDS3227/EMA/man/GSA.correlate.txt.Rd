\name{GSA.correlate.txt}
\alias{GSA.correlate.txt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Correlation between Genes collection and Genes Array }
\description{
  Function to summarize the gene set used and with the genes of interest.
}

\usage{
GSA.correlate.txt(GSA.genesets.obj, genenames)
}


\arguments{
  \item{GSA.genesets.obj}{ a gmt object obtain by 'GSA.read.gmt'
    function.
  }
  \item{genenames}{ a character with genenames corresponding to
    probesets.
  }
}

\value{
  only print output
}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{
  \code{\link{GSA.read.gmt}} 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print}
\keyword{internal}
