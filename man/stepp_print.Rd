% --- Source file: stepp_print.Rd ---
\name{stepp_print}
\alias{stepp_print}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Subpopulation Treatment Effect Pattern Plot (STEPP) Accessory Function - stepp_print
}
\description{
A method to print the computational results of a stepp analysis that is captured in the stepp object.   
}
\usage{
stepp_print(x, estimate=TRUE, cov=TRUE, test=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the stepp object created using stepp}
  \item{estimate}{print the estimate, default is TRUE}
  \item{cov}{print the covariance matrix, default is TRUE}
  \item{test}{print the permutation test p-values, default is TRUE}
}
\details{
  stepp_print produces the result of the estimates, the covariance matrix and the pvalues generated from the permutation approach.  One has the option of printing any of the three outputs by setting the corresponding options to FALSE.  The default is TRUE for all three outputs.

  See other related stepp functions:
  1. stepp
  2. stepp_summary
  3. stepp_print
  4. stepp_plot
  5. analyze.KM.stepp
  6. analyze.CumInc.stepp

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.


