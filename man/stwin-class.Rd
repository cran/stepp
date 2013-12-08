\name{stwin-class}
\Rdversion{1.1}
\docType{class}
\alias{stwin-class}
\alias{summary,stwin-method}

\title{Class \code{"stwin"}}
\description{
  This the S4 class for the stepp window object.  The stepp window class describes the way
  to set up the subpopulation for a stepp analysis.  Each window represents a subpopulation 
  along the covariate of interest.  The only supported window type is "sliding" which means
  subpopulations overlap one another as the analysis is done along the covariate
  of interest.  Two parameters control how the subpopulations are constructed: r2 
  specifies the size of each subpopulation and r1 specifies the minimum number of 
  patients allowed to overlap between two adjacent subpopulations.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("stwin", type="sliding", r1=5, r2=20 ) or
  the constructor function stepp.win}.
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"} \cr
		stepp window type; currently only support "sliding" }
    \item{\code{r1}:}{Object of class \code{"numeric"} \cr 
		minimum number of patients allowed to overlap between subpopulations}
    \item{\code{r2}:}{Object of class \code{"numeric"} \cr
		size of each subpopulation }
    \item{\code{nevent}:}{Object of class \code{"numeric"} \cr 
		reserve for future use }
  }
}
\section{Methods}{
  \describe{
    \item{summary}{\code{signature(object = "stwin")}: \cr
		print a summary of the stepp windows object }
	 }
}

\author{
  Wai-Ki Yip
}

\seealso{
                                  \code{\linkS4class{stsubpop}},  \code{\linkS4class{stmodelKM}}, 
  \code{\linkS4class{stmodelCOX}},\code{\linkS4class{stmodelCI}}, \code{\linkS4class{stmodelGLM}}, 
  \code{\linkS4class{steppes}},   \code{\linkS4class{stmodel}}, 
  \code{\link{stepp.win}},        \code{\link{stepp.subpop}},     \code{\link{stepp.KM}},
  \code{\link{stepp.CI}},         \code{\link{stepp.COX}},        \code{\link{stepp.GLM}},
  \code{\link{stepp.test}},       \code{\link{estimate}},         \code{\link{generate}}

}
\examples{
showClass("stwin")

# create a stepp window of type "sliding" with (r2) size of subpopulation
# in each window to be 200 and (r1) allows only 50 patients in the
# overlapping windows
mywin <- new("stwin", type="sliding", r1=50, r2=200)

# print a summary of the stepp window object created
summary(mywin)
}
\keyword{classes}
