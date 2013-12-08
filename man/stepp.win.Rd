\name{stepp.win}
\alias{stepp.win}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  a method to create the stepp window object
}

\description{
  This is the constructor function to create a stepp window object.
}

\usage{
  stepp.win(type, r1, r2)
}

\arguments{
  \item{type}{type of stepp window; only "sliding" is supported}
  \item{r1}{minimum number of patients allowed in overlapping windows}
  \item{r2}{size of subpopulation in each window}
}

\details{
  This is the functional interface to construct a stepp window object besides
  using new("stwin", ...).    

}


\author{
 Wai-ki Yip
}


\seealso{
  \code{\linkS4class{stwin}},     \code{\linkS4class{stsubpop}},  \code{\linkS4class{stmodelKM}}, 
  \code{\linkS4class{stmodelCOX}},\code{\linkS4class{stmodelCI}}, \code{\linkS4class{stmodelGLM}}, 
  \code{\linkS4class{steppes}},   \code{\linkS4class{stmodel}}, 
                                  \code{\link{stepp.subpop}},     \code{\link{stepp.KM}},
  \code{\link{stepp.CI}},         \code{\link{stepp.COX}},        \code{\link{stepp.GLM}},
  \code{\link{stepp.test}},       \code{\link{estimate}},         \code{\link{generate}}

}

\examples{
  # create a stepp window object of type "sliding", 
  # subpopulation size is 200 and allows only 50 patients
  # between overlapping windows
  mywin <- stepp.win(type="sliding", r1=50, r2=200)

  # print a summary of the stepp window object
  summary(mywin) 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
