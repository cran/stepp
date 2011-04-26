% --- Source file: stepp_plot.Rd ---
\name{stepp_plot}
\alias{stepp_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Subpopulation Treatment Effect Pattern Plot (STEPP) Accessory Function - stepp_plot
}
\description{
A method to generate stepp plots from the results of a stepp analysis that is captured in the stepp object.   
}
\usage{
stepp_plot(x, legendy = 30, pline = -2.5, color = c("red", "black"), ylabel= "Specify Timepoint & Endpoint", xlabel="Subpopulations by Median Covariate", ncex = 0.7, tlegend=c("Specify 1st Treatment", "Specify 2nd Treatment"), nlas = 0, alpha = 0.05, pointwise = FALSE, diff = TRUE, ci = TRUE, pv = TRUE, showss = TRUE, ylimit=c(0,100,-100,100,0,3), dev="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the stepp object created using stepp.}
  \item{legendy}{optional - the vertical location of the legend according to the units on the y-axis; default is 30.}
  \item{pline}{optional - specify the vertical location of the p-value, that is, on which MARgin line, starting at 0 counting outwards; default is -2.5.}
  \item{color}{optional - a vector containing the line colors for the 1st and 2nd treatment, respectively; default is "red" and "black".}
  \item{ylabel}{optional - specify the label for the y-axis; default is "Specify Timepoint & Endpoint".}
  \item{ylimit}{optional - specify the limits of the y values on the y-axis; default is c(0,100,-100,100,0,3) for the 3 plots.}
  \item{xlabel}{optional - specify the label for the x-axis; default is "Subpopulations by Median Covariate".}
  \item{ncex}{optional - specify the size of the text for the sample size annotation, that is, the character expansion factor; default is 0.7.}
  \item{tlegend}{optional - a vector containing the treatment labels, 1st and 2nd treatment, respectively; default is c("Specify 1st Treatment", "Specify 2nd Treatment").}
  \item{nlas}{optional - specify the las paramter (0,1,2,3) to determine the orientation of the sample size annotation; default is 0.}
  \item{alpha}{optional - specify the significance level; default is 0.05.}
  \item{pointwise}{optional - specify pointwise confidence intervals (pointwise=TRUE), or confidence bands; default is FALSE.}
  \item{diff}{optional - specify if you want just the first plot; default is TRUE (all three plots).}
  \item{ci}{optional - specify if you want to display the conf. interval or band; default is TRUE.}
  \item{pv}{optional - pvalue will be displayed in the plots; FALSE: pvalue will not be displayed in the plots; default is TRUE.}
  \item{showss}{optional - show the sample size; FALSE: sample size will not be shown in the plots; default is TRUE.}
  \item{dev}{optional - output device control.  Output the plots to the specified output format.
	postscript for .ps; eps (encapsulated postscript) for .eps; pdf for .pdf; png for .png; bmp for .bmp; tiff for .tif and jpeg for .jpeg.  
      The filenames are always SteppPlot1, SteppPlot2 and SteppPlot3 in the R working directory; default is "" }
}
\details{
  stepp_plot produces the graphics (STEPP plots) for analysis.  It produces three plots based on the result captured in the stepp object.

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


