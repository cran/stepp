\name{stepp.COX}
\alias{stepp.COX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	a method to create the stmodelCOX object
}
\description{
  	This is the constructor function for the stmodelCOX object.  This object sets up the data with
	a stepp model using Cox method for analysis.
}
\usage{
	stepp.COX(coltrt, survTime, censor, trts, timePoint, MM)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coltrt}{the treatment variable}
  \item{survTime}{the time to event variable}
  \item{censor}{the censor variable}
  \item{trts}{a vector containing the codes for the 2 treatment arms, 1st and 2nd treatment arms, respectively}
  \item{timePoint}{timepoint to estimate survival}
  \item{MM}{a model matrix for any additional covariates; default is NULL}
}
\details{
  One can specify additional covariates for adjustment when using the COX PH model using model.matrix in R.
  \code{mm <- model.matrix(~cov1+cov2+...+covn)}. Since there is no intercept term in the COX model, the 
  intercept column will be automatically stripped if it is used.
}
\value{
  It returns a stmodelCOX object.
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
modCOX  <- stepp.COX(coltrt=Txassign, survTime=time, censor=censor, 
                     trts=c(1,2), timePoint=4, MM=mm0)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
