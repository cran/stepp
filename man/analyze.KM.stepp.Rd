\name{analyze.KM.stepp}
\alias{analyze.KM.stepp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Analyze survival data using Kaplan-Meier method
}

\description{
  This method will be deprecated in the future.  Please use the constructor function stepp.KM to
  create a STEPP Kaplan-Meier model for future development.

  A method to explore the treatment-covariate interactions in survival data arising 
  from two treatment arms of a clinical trial. The treatment effects are measured using survival 
  functions at a specified time point estimated from the Kaplan-Meier method and the hazard ratio 
  based on observed-minus-expected estimation. A permutation distribution approach to inference 
  is implemented, based on permuting the covariate values within each treatment group. 
  The statistical significance of observed heterogeneity of treatment effects is calculated using 
  permutation tests:\cr\cr
  1) for the maximum difference between each subpopulation effect and the overall population 
     treatment effect or supremum based test statistic; \cr
  2) for the difference between each subpopulation effect and the overall population treatment 
     effect, which resembles the chi-square statistic. 
}

\usage{
  analyze.KM.stepp(coltrt, coltime, colcens, covar, trts, patspop, minpatspop, 
    timest, noperm=2500, ncex = 0.7, legendy = 30, pline = -2.5, 
    color = c("red", "black"),
    xlabel = "Subpopulations by Median Covariate", 
    ylabel = "?-year Disease-Free Survival", 
    tlegend = c("1st Treatment", "2nd Treatment"), 
    nlas = 3, pointwise = FALSE)
}

\arguments{
  \item{coltrt}{the treatment variable}
  \item{coltime}{the time to event variable}
  \item{colcens}{the censoring variable}
  \item{covar}{the covariate of interest}
  \item{trts}{a vector containing the codes for the 2 treatment arms, 1st and 2nd treatment arms, respectively}
  \item{patspop}{larger parameter(r2) for subpopulation construction that determines how many patients are in each subpopulation}
  \item{minpatspop}{smaller parameter(r1) for subpopulation construction that determines the largest number of patients in common among consecutive subpopulations}
  \item{timest}{timepoint to estimate survival}
  \item{noperm}{the desired number of permutations; must be 0 or above}
  \item{ncex}{optional - specify the size of the text for the sample size annotation, that is, the character expansion factor}
  \item{legendy}{optional - the vertical location of the legend according to the units on the y-axis}
  \item{pline}{optional - specify the vertical location of the p-value, that is, on which margin line, starting at 0 counting outwards}
  \item{color}{optional - a vector containing the line colors for the 1st and 2nd treatment, respectively}
  \item{ylabel}{optional - specify the label for the y-axis}
  \item{xlabel}{optional - specify the label for the x-axis}

  \item{tlegend}{optional - a vector containing the treatment labels, 1st and 2nd treatment, respectively}
  \item{nlas}{optional - specify the las parameter (0,1,2,3) to determine the orientation of the sample size annotation}
  \item{pointwise}{optional -specify pointwise confidence intervals (pointwise=TRUE), or confidence bands (pointwise=FALSE, default) to be displayed}
}

\details{
  A statistical method to explore treatment-covariate interactions in survival data arising from two 
  treatment arms of a clinical trial. The method is based on constructing overlapping subpopulations 
  of patients with respect to one covariate of interest, and in observing the pattern of the treatment 
  effects estimated across subpopulations. A plot of these treatment effects is called STEPP, or 
  Subpopulation Treatment Effect Pattern Plot. STEPP uses the permutation distribution based approach
  for inference.

  One can explore the window parameters without invoking the permutation analysis by setting noperm 
  to 0.  In that case, pvalue and the covarance matrix will not be available.

}

\value{
  analyze.KM.stepp generates a Subpopulation Treatment Effect Pattern Plot (STEPP) displaying 
  the p-value from the test for interaction. Descriptive summaries of the dataset and the estimated 
  variance-covariance matrix are returned through the steppes object. See documentation on steppes 
  object for details on how you can use it. 

}

\references{
  Bonetti M, Gelber RD. Patterns of treatment effects in subsets of patients in clinical trials. Biostatistics 2004; 5(3):465-481.

  Bonetti M, Zahrieh D, Cole BF, Gelber RD. A small sample study of the STEPP approach to assessing treatment-covariate interactions in survival data. Statistics in Medicine 2009; 28(8):1255-68. 

}

\author{
 Wai-ki Yip, David Zahrieh, Marco Bonetti, Bernard Cole, Ann Lazar, Richard Gelber
}

\note{
  STEPP is an exploratory tool, with graphical features that make it easy for clinicians to interpret 
  the results of the analysis. The method provides an opportunity to detect interactions other than 
  those that may be apparent for example by using Cox models. Positive results should prompt the need 
  for confirmation from other datasets investigating similar treatment comparisons. It should also be 
  clear that STEPP is not meant to estimate specific cutpoints in the range of values of the covariate 
  of interest, but rather to provide some indication on ranges of values where treatment effect might 
  have a particular behavior.

  STEPP considers the case in which the subpopulations are constructed according to a sliding window 
  pattern. The larger parameter (patspop) determines how many patients are in each subpopulation, and 
  the smaller parameter (minpatspop) determines the largest number of patients in common among 
  consecutive subpopulations. A minimum of 80-100 patients should probably be in each subpopulation, 
  but that is not strictly necessary. The difference (patspop-minpatspop) is the approximate number of 
  patients replaced between any two subsequent subpopulations, and can be used to determine the number 
  of subpopulations once patspop is fixed. The choice of the values of the parameters patspop and 
  minpatspop to be used does change the appearance of the plot and the corresponding p-value. It is 
  probably reasonable to experiment with a few combinations to ensure that the significance (or lack 
  of significance) is stable with respect to that choice.

  For best results, consider implementing 2500 permutations of the covariate (vector of subpopulations)
  to obtain a rich distribution in which to draw inference. 
}

\section{Warning }{
  This function together with other old functions will be depreciated in the future.  A new set of S4 classes are
  implemented to replace old interfaces.  Please use them for future development.

  A few tips to keep in mind:

  The variables coltrt, coltime, colcens, and covar must be numeric. No formatting allowed.

  If you receive the error "Error in solve.default(sigma): system is computationally singular; reciprocal condition number = 0" then we recommend changing the seed by re-running analyze.KM.stepp. If this error persists after several runs, then the program cannot provide reliable results. Please try modifying your choices of the two parameters minpatspop(r1) and patspop(r2) that define the subpopulation.

  The number of permutations specified in noperm, the sample size, and the number of subpopulations generated will affect how long analyze.KM.stepp takes to execute. The results are stable if 2500 or more permutations are specified. Furthermore, varying the number of subpopulations will affect inference.

  The time point selected to estimate survival in timest must be in the same units (e.g., months) as the coltime variable.

  The order of the treatments in the vector trts must be in the same order in the vector tlegend. 
}

\seealso{
  \code{\linkS4class{stwin}},     \code{\linkS4class{stsubpop}},  \code{\linkS4class{stmodelKM}}, 
  \code{\linkS4class{steppes}},   \code{\linkS4class{stmodel}}, 
  \code{\link{stepp.win}},        \code{\link{stepp.subpop}},     \code{\link{stepp.KM}},
  \code{\link{stepp.test}},       \code{\link{estimate}},         \code{\link{generate}}

  Old functions to be deprecated: stepp, stepp_summary, stepp_print, stepp_plot, and analyze.CumInc.stepp. 

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

#CALL analyze.KM.stepp to analyze the data
# Warning: In this example, the permutations have been set to 0 to allow 
# the stepp function to finish in a short amount of time.  IT IS RECOMMEND
# TO USE AT LEAST 2500 PERMUTATIONS TO PROVIDE STABLE RESULTS.
output <- analyze.KM.stepp ( coltrt=Txassign, coltime=time, 
  colcens=censor, covar=covariate, trts=c(1,2), patspop=300,
  minpatspop=200, timest=4, noperm=0, ncex=0.70, legendy=30,
  pline=-2.5, color=c("red", "black"),
  xlabel="Subpopulations by Median Age", ylabel="4-year Disease-Free Survival", 
  tlegend=c("Treatment A", "Treatment B"), nlas=3, pointwise=FALSE)
}
