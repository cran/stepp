\name{stmodelCOX-class}
\Rdversion{1.1}
\docType{class}
\alias{stmodelCOX-class}
\alias{estimate,stmodelCOX-method}
\alias{print,stmodelCOX-method}
\alias{test,stmodelCOX-method}

\title{Class \code{"stmodelCOX"}}
\description{
	This is the stepp model for survival data using COX proportional hazard
}
\section{Objects from the Class}{
	Objects can be created by calls of the form \code{new("stmodelCOX", ...)} or by \cr
  	the construction function stmodel.COX.
}
\section{Slots}{
  \describe{
    \item{\code{MM}:}{Object of class \code{"ANY"} \cr
	a model matrix created using \code{model.matrix} in R to specify additional covariates for adjustment;
	default is NULL }
    \item{\code{coltrt}:}{Object of class \code{"numeric"} \cr
	the treatment variable }
    \item{\code{survTime}:}{Object of class \code{"numeric"} \cr
	the time to event variable }
    \item{\code{censor}:}{Object of class \code{"numeric"} \cr
	the censor variable }
    \item{\code{trts}:}{Object of class \code{"numeric"} \cr
	a vector containing the codes for the 2 treatment arms, 1st and 2nd treatment arms, respectively}
    \item{\code{timePoint}:}{Object of class \code{"numeric"} \cr
	timepoint to estimate survival}
  }
}
\section{Extends}{
Class \code{"\linkS4class{stmodelKM}"}, directly.
Class \code{"\linkS4class{stmodel}"}, directly.
}
\section{Methods}{
  \describe{
    \item{estimate}{\code{signature(.Object = "stmodelCOX")}: \cr
	estimate the effect in absolute and relative scale of the overall and each subpopulation }
    \item{print}{\code{signature(x = "stmodelCOX")}: \cr
	print the estimate, covariance matrices and statistics }
    \item{test}{\code{signature(.Object = "stmodelCOX")}: \cr
	perform the permutation tests and obtain various statistics }
  }
}

\value{
  The new method returns the stmodelCOX object.

  The estimate method returns a list with the following fields:

  \item{model}{the stepp model - "COXe"}
  \item{sObs1}{a vector of effect estimates of all subpopulations based on the 1st treatment}
  \item{sSE1}{a vector of standard errors of effect estimates of all subpopulations based on the 1st treatment}
  \item{oObs1}{effect estimate of the entire population based on the 1st treatment}
  \item{oSE1}{the standard error of the effect estimate of the entire population based on the 1st treatment}
  \item{sObs2}{a vector of effect estimates of all subpopulations based on the 1st treatment}
  \item{sSE2}{a vector of standard errors of effect estimates of all subpopulations based on the 1st treatment}
  \item{oObs2}{effect estimate of the entire population based on the 1st treatment}
  \item{oSE2}{the standard error of the effect estimate of the entire population based on the 1st treatment}
  \item{skmw}{Wald's statistics for the effect estimate differences between the two treatments}
  \item{logHR}{a vector of log hazard ratio estimate of the subpopulations comparing 1st and 2nd treatments}
  \item{logHRSE}{a vector of standard error of the log hazard ratio estimate of the subpopulations comparing 1st and 2nd treatments}
  \item{ologHR}{the log hazard ratio estimate of the entire population comparing 1st and 2nd treatments}
  \item{ologHRSE}{the standard error of the log hazard ratio estimate of the entire population comparing 1st and 2nd treatments}
  \item{logHRw}{Wald's statistics for the log hazard ratio between the two treatments}

  The test method returns a list with the following fields:

  \item{model}{the stepp model - "COXt"}
  \item{sigma}{the covariance matrix for subpopulations based on effect differences}
  \item{hasigma}{the homogeneous association covariance matrix for subpopulations based on effect differences}
  \item{HRsigma}{the covariance matrix for the subpopulations based on hazard ratio}
  \item{haHRsigma}{the homogeneous association covariance matrix for subpopulations based on hazard ratio}
  \item{pvalue}{the supremum pvalue based on effect difference}
  \item{chi2pvalue}{the chisquare pvalue based on effect difference}
  \item{hapvalue}{the homogeneous association pvalue based on effect difference}
  \item{HRpvalue}{the supremum pvalue based on hazard ratio}
  \item{haHRpvalue}{the homogeneous association pvalue based on hazard ratio}  
  \item{GEEpvalue}{the GEE based pvalue}
}

\author{
  Wai-Ki Yip
}

\seealso{
  \code{\linkS4class{stwin}},     \code{\linkS4class{stsubpop}},  \code{\linkS4class{stmodelKM}}, 
  \code{\linkS4class{stmodelCOX}},\code{\linkS4class{stmodelCI}}, \code{\linkS4class{stmodelGLM}}, 
  \code{\linkS4class{steppes}},   \code{\linkS4class{stmodel}}, 
  \code{\link{stepp.win}},        \code{\link{stepp.subpop}},     \code{\link{stepp.KM}},
  \code{\link{stepp.CI}},                                         \code{\link{stepp.GLM}},
  \code{\link{stepp.test}},       \code{\link{estimate}},         \code{\link{generate}}

}
\examples{
showClass("stmodelCOX")
#GENERATE TREATMENT VARIABLE:
N <- 1000
Txassign <- sample(c(1,2), N, replace=TRUE, prob=c(1/2, 1/2))
n1 <- length(Txassign[Txassign==1])
n2 <- N - n1
#GENERATE A COVARIATE:
covariate <- rnorm(N, 55, 7)
#GENERATE SURVIVAL AND CENSORING VARIABLES ASSUMING A TREATMENT COVARIATE INTERACTION:
Entry <- sort( runif(N, 0, 5) )
SurvT1 <- .5
beta0 <-  -65 / 75
beta1 <- 2 / 75
Surv <- rep(0, N)
lambda1 <- -log(SurvT1) / 4
Surv[Txassign==1] <- rexp(n1, lambda1)
Surv[Txassign==2] <- rexp(n2, (lambda1*(beta0+beta1*covariate[Txassign==2])))
EventTimes <- rep(0, N)
EventTimes <- Entry + Surv
censor <- rep(0, N)
time <- rep(0,N)
for ( i in 1:N )
   {
     censor[i] <- ifelse( EventTimes[i] <= 7, 1, 0 )
     time[i] <- ifelse( EventTimes[i] < 7, Surv[i], 7 - Entry[i] ) 
   }

c1 <- rnorm(1000)

# create the model matrix to represent extra covariate for adjustment
mm0 <- model.matrix(~c1)

# create the stepp model for survival data using COX PH model for analysis
modCOX  <- new("stmodelCOX", coltrt=Txassign, survTime=time, censor=censor, trts=c(1,2), 
               timePoint=4, MM=mm0)

}
\keyword{classes}
