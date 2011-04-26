#h.packageName <- "stepp"


#####################################################################
#
# Plot the result of the stepp analysis
#
#####################################################################
stepp_plot <- function (x,			# stepp obj 
							# graphics parameters - optional
			legendy = 30,		#   the vertical location of the legend according to the units on the y-axis
			pline = -2.5,		#   the vertical location of the p-value, starting at 0 counting outwards
			color = c("red", "black"),
							#   a vector containing the line colors for the 1st and 2nd trt, respectively
			ylabel = "Specify Timepoint & Endpoint",   	# label for the y-axis
			xlabel = "Subpopulations by Median Covariate", 	# label for the x-axis
    			ncex = 0.7,			#   the size of the text for the sample size annotation
			tlegend = c("Specify 1st Treatment", "Specify 2nd Treatment"),
							#   a vector containing the treatment labels, 1st and 2nd trt, respectively
			nlas = 0,			#   the las paramter (0,1,2,3) - the orientation of the sample size annotation
							# plot options - optional
			alpha = 0.05,		#   sig. level
			pointwise = FALSE,	#   pointwise confidence intervals (pointwise=TRUE),
							#   or confidence bands (pointwise=FALSE, default) to be displayed
			diff = TRUE,		#   generate 2 additional plots comparing the diff between measures
			ci = TRUE,			#   display the conf. interval or band
			pv = TRUE,			#   display the supremum pvalue
			showss = TRUE,		#   display sample size on the x-axis
			ylimit = c(0,100,-100,100,0,3),			# y limits for the 3 graphs
			dev = ""			#   graphics device for output; default to Null Device
		       )
{
    skmObs1 <- x$est$skmObs1
    skmObs2 <- x$est$skmObs2
    skmSE1  <- x$est$skmSE1
    skmSE2  <- x$est$skmSE2
    nsubpop <- x$sp$nsubpop
    npatsub <- x$sp$npatsub
    medians <- x$sp$medians
    r1      <- x$arg$minpatspop
    r2      <- x$arg$patspop
    pvalue  <- x$res$pvalue
    HRpvalue<- x$res$HRpvalue
    logHR   <- x$est$logHR
    logHRSE <- x$est$logHRSE

#   for(i in 1:3) win.graph()
    for(i in 1:3) {
	if (dev == "") dev.new()
	else
	if (dev == "postscript") {
	  fname = paste("SteppPlot",as.character(i),".ps",sep="")
	  postscript(file=fname)
	} 
	else  
	if (dev == "eps") {
	  fname = paste("SteppPlot",as.character(i),".eps",sep="")
	  postscript(file=fname)
	} 
	else 
	if (dev == "pdf") {
	  fname = paste("SteppPlot",as.character(i),".pdf",sep="")
	  pdf(file=fname)
	} 
	else 
	if (dev == "png") {
	  fname = paste("SteppPlot",as.character(i),".png",sep="")
	  png(filename=fname)
	} 
	else 
	if (dev == "bmp") {
	  fname = paste("SteppPlot",as.character(i),".bmp",sep="")
	  bmp(filename=fname)
	} 
	else 
	if (dev == "tiff") {
	  fname = paste("SteppPlot",as.character(i),".tif",sep="")
	  tiff(filename=fname)
	} 
	else 
	if (dev == "jpeg") {
	  fname = paste("SteppPlot",as.character(i),".jpeg",sep="")
	  jpeg(filename=fname)
	} 
    }
    devlst <- dev.list()

#   generate the first stepp plot
#     STEPP analysis of treatment effect as measured
#     by KM/HR or cumulative incidence.
#
    dev.set(devlst[length(devlst)-2])
    skmObs <- rep(0, nsubpop * 2)
    xvalues <- rep(0, nsubpop * 2)
    skmObs <- c(skmObs1, skmObs2)
    xvalues <- c(medians, medians)
    skmObs <- skmObs * 100
    group <- rep(1:0, each = nsubpop)
    lbls <- rep(" ", nsubpop)
    ssize <- rep(" ", nsubpop)
    for (i in 1:nsubpop) ssize[i] <- paste(c("(n=", npatsub[i], ")"), collapse = "")
    p <- paste(c("Supremum p-value = ", pvalue), collapse = "")
    par(mfrow = c(1, 1), omi = c(0.5, 0.5, 0.5, 0.5))

    if (length(ylimit) < 2) yl <- c(0,100)
    else yl <- ylimit[1:2]
    plot(xvalues, skmObs, axes = FALSE, ylim = yl, ylab = ylabel, xlab = "")
    points(xvalues[group == 1], skmObs[group == 1], lty = 1, lwd = 2, pch = 19, type = "o", col = color[1], bg = color[1])
    points(xvalues[group == 0], skmObs[group == 0], lty = 1, lwd = 2, pch = 24, type = "o", col = color[2], bg = color[2])
    axis(1, at = xvalues, font = 1)
    axis(2, at = c(0, (0:5) * 20), font = 1)
    if (nlas != 3 & nlas != 2) {
        if (showss) {
	    mtext(ssize, side = 1, at = xvalues, line = 2, font = 1, cex = ncex, adj = 0.5, las = nlas)
          mtext(xlabel, side = 1, line = 3.5)
	  } else
	  {
	    mtext(xlabel, side = 1, line = 2)
	  }
    }
    if (nlas == 3 | nlas == 2) {
        if (showss) {
	    mtext(ssize, side = 1, at = xvalues, line = 3.5, font = 1, cex = ncex, adj = 0.5, las = nlas)
          mtext(xlabel, side = 1, line = 0.8, outer = TRUE)
	  } else
	  {
	    mtext(xlabel, side = 1, line = 2)
	  }
    }
    legend(min(xvalues), legendy, pch = c(19, 24), lty = 1, lwd = 2, col = color, pt.bg = color, legend = tlegend, bty = "n")
    
    if (pv) mtext(p, side = 1, at = (min(xvalues) + 0.2 * (max(xvalues) - min(xvalues))), line = pline)
    
#
    if(diff){
#   	pointwise is specified, generate two additional plots
#
       dev.set(devlst[length(devlst)-1])
       if(pointwise) zcrit <- qnorm(1-alpha/2)
       else zcrit <- qnorm(1-alpha/(2*nsubpop))
       skmObs <- rep(0, nsubpop)
       xvalues <- rep(0, nsubpop)
       skmObs <- skmObs1-skmObs2
       se <- sqrt(skmSE1^2+skmSE2^2)
       xvalues <- medians
       skmObs <- skmObs * 100
       lbls <- rep(" ", nsubpop)
       ssize <- rep(" ", nsubpop)
       for (i in 1:nsubpop) ssize[i] <- paste(c("(n=", npatsub[i], ")"), collapse = "")
       p <- paste(c("Supremum p-value = ", pvalue), collapse = "")
       par(mfrow = c(1, 1), omi = c(0.5, 0.5, 0.5, 0.5))
	 if (length(ylimit) < 4) yl = c(-100,100)
	 else yl = ylimit[3:4]
       plot(xvalues, skmObs, axes = FALSE, ylim = yl, ylab = paste("Difference in",ylabel), xlab = "")
       points(xvalues, skmObs, lty = 1, lwd = 2, pch = 19, type = "o", col = color[1], bg = color[1])
       if(ci){
          cilim <- skmObs - zcrit*se*100
          cilim <- pmax(cilim,-100)
          lines(xvalues, cilim, lty = 2, lwd = 2, col = color[1], bg = color[1])
          cilim <- skmObs + zcrit*se*100
          cilim <- pmin(cilim,100)
          lines(xvalues, cilim, lty = 2, lwd = 2, col = color[1], bg = color[1])
       }
       lines(c(min(xvalues),max(xvalues)),c(0,0),lty=1)
       axis(1, at = xvalues, font = 1)
       axis(2, at = c(0, (-5:5) * 20), font = 1)
       if (nlas != 3 & nlas != 2) {
           if (showss) {
	       mtext(ssize, side = 1, at = xvalues, line = 2, font = 1, cex = ncex, adj = 0.5, las = nlas)
             mtext(xlabel, side = 1, line = 3.5)
	     } else {
	       mtext(xlabel, side = 1, line = 2)
	     }
       }
       if (nlas == 3 | nlas == 2) {
           if (showss) {
	       mtext(ssize, side = 1, at = xvalues, line = 3.5, font = 1, cex = ncex, adj = 0.5, las = nlas)
             mtext(xlabel, side = 1, line = 0.8, outer = TRUE)
	     } else {
	       mtext(xlabel, side = 1, line = 2)
	     }	
       }
       if (pv) mtext(p, side = 1, at = (min(xvalues) + 0.2 * (max(xvalues) - min(xvalues))), line = pline)
      
       dev.set(devlst[length(devlst)])
       p <- paste(c("Supremum HR p-value = ", HRpvalue), collapse = "")
       par(mfrow = c(1, 1), omi = c(0.5, 0.5, 0.5, 0.5))
       if (length(ylimit) < 6) yl = c(0,3)
	 else yl = ylimit[5:6]
       plot(xvalues, exp(logHR), axes = FALSE, ylim = yl, ylab = "Hazard Ratio", xlab = "")
       points(xvalues, exp(logHR), lty = 1, lwd = 2, pch = 19, type = "o", col = color[1], bg = color[1])
       if(ci){
          cilim <- logHR - zcrit*logHRSE
          cilim <- exp(cilim)
          lines(xvalues, cilim, lty = 2, lwd = 2, col = color[1], bg = color[1])
          cilim <- logHR + zcrit*logHRSE
          cilim <- exp(cilim)
          lines(xvalues, cilim, lty = 2, lwd = 2, col = color[1], bg = color[1])
       }
       lines(c(min(xvalues),max(xvalues)),c(1,1),lty=1)
       axis(1, at = xvalues, font = 1)
       axis(2, at = c(0, (0:15)*.2), font = 1)
       if (nlas != 3 & nlas != 2) {
           if (showss) {
	       mtext(ssize, side = 1, at = xvalues, line = 2, font = 1, cex = ncex, adj = 0.5, las = nlas)
             mtext(xlabel, side = 1, line = 3.5)
	     } else {
		 mtext(xlabel, side = 1, line = 2)
	     }
       }
       if (nlas == 3 | nlas == 2) {
           if (showss) {
		 mtext(ssize, side = 1, at = xvalues, line = 3.5, font = 1, cex = ncex, adj = 0.5, las = nlas)
             mtext(xlabel, side = 1, line = 0.8, outer = TRUE)
	     } else
	     {
		 mtext(xlabel, side = 1, line = 2)
	     }
       }
       if (pv) mtext(p, side = 1, at = (min(xvalues) + 0.2 * (max(xvalues) - min(xvalues))), line = pline)
     }

    for(i in 1:3) {
	if (dev != "") dev.off(devlst[i])
    }
}




#####################################################################
#
# Generate a summary of the stepp analysis
#
#####################################################################
stepp_summary<-function(x)
{
#
#   print version and time
#
    cat("\n")
    write(paste("STEPP Version ", x$version),file="")
    cat("\n")
    write(paste("Creation date and time: ", x$cdate),file="")
    cat("\n")
#
#   The time reported is system dependent.  Usually, proc.time is 
#   accurate to 10 ms as in Windows.
    write(paste("Time elapsed: ", x$etime["elapsed"]/60, " secs"),file="")
    cat("\n")
#
    trttype <- x$arg$trttype
    timePoint <- x$arg$timest
    r1 <- x$arg$minpatspop
    r2 <- x$arg$patspop
    noPerms <- x$arg$noperm
    survTime <- x$arg$coltime
    censor <- x$arg$colcens
    type <- x$arg$coltype
    covariate <- x$arg$colvar
#
#   create 0/1 treatment assignment 
#
    txassign <- rep(NA, length(x$arg$coltrt))
    txassign[which(x$arg$coltrt == x$arg$trts[1])]<-1
    txassign[which(x$arg$coltrt == x$arg$trts[2])]<-0
#
#   set up the data
#
    if (trttype == "KM") indata <- cbind(txassign, survTime, censor, covariate)
    else indata <- cbind(txassign, survTime, type, covariate)
    indata <- indata[apply(indata, 1, function(x) !any(is.na(x))),, drop = FALSE]
    txassign <- indata[,1]
    survTime <- indata[,2]
    if (trttype == "KM") censor <- indata[, 3]
    else type <- indata[,3]
    covariate <- indata[,4]
#
#   treatment group summary information
#
    if (trttype == "KM") write ("Using KM estimates", file="")
    else write ("Using Cumulative Incidences", file="")
    cat("\n") 
    write(paste("STEPP is comparing treatment",x$arg$trts[1],"versus treatment",x$arg$trts[2]),file="")
    cat("\n")
    N <- length(indata[,1])
    n1 <- length(txassign[txassign == 1])
    n2 <- N-n1
    write(paste("            Sample size in treatment",x$arg$trts[1],":", n1),file="")
    write(paste("            Sample size in treatment",x$arg$trts[2],":", n2),file="")
    write(paste(" Total sample size (excluding missing) :", N),file="")
    write(paste("                Number of permutations :", noPerms),file="")
    write(paste("        Number of patients per subpopulation (patspop r2):", r2),file="")
    write(paste("Largest number of patients in common among consecutive subpopulations(minpatspop r1):", r1),file="")
# 
#   Subpopulation summary
#
    write(paste("      Number Of Subpopulations Created :", x$sp$nsubpop),file="")
    cat("\n")
    write("Subpopulation Summary Information",file="")
    nper <- apply(x$sp$subpop,2,sum)
    temp <- matrix(c(1:x$sp$nsubpop,x$sp$medians,x$sp$minz,x$sp$maxz,nper),ncol=5)
    write("                                  Covariate Summary                  Sample",file="")
    write("     Subpopulation        Median       Minimum       Maximum          Size",file="")
    for (i in 1:x$sp$nsubpop) {
       write(paste(format(temp[i,1],width=12),format(temp[i,2],width=19,nsmall=2),
          format(temp[i,3],width=13,nsmall=4),format(temp[i,4],width=13,nsmall=4),
          format(temp[i,5],width=13)),file="")
    }
}

#################################################################
#
# Generate a print output of the stepp analysis
#
#################################################################

stepp_print<-function(x, estimate=TRUE, cov=TRUE, test=TRUE){
    timePoint <- x$arg$timest
#
#
#   1. survival or cumulative incidence estimates
#

    if (x$arg$noperm == 0)
	print (paste("WARNING: The number of permutations must be greater than 0; no p value can be generated"))

    if (estimate){
      if (x$arg$trttype == "KM"){
        cat("\n")
        write(paste("Survival estimates for treatment group", x$arg$trts[1], 
            "at time point", timePoint), file = "")
        temp <- matrix(c(1:x$sp$nsubpop, round(x$est$skmObs1, digits = 4), round(x$est$skmSE1, 
            digits = 4)), ncol = 3)
        write("                         Survival", file = "")
        write("     Subpopulation     Probability      Std. Err.", 
            file = "")
        for (i in 1:x$sp$nsubpop) {
            write(paste(format(temp[i, 1], width = 12), format(temp[i, 
                2], width = 19, nsmall = 4), format(temp[i, 3], width = 15, 
                nsmall = 4)), file = "")
        }
        write(paste("        Overall", format(round(x$est$overallSkmObs1, 
            digits = 4), nsmall = 4, width = 16), format(round(x$est$overallSkmSE1, 
            digits = 4), nsmall = 4, width = 15)), file = "")
        cat("\n")
        write(paste("Survival estimates for treatment group", x$arg$trts[2], 
            "at time point", timePoint), file = "")
        temp <- matrix(c(1:x$sp$nsubpop, round(x$est$skmObs2, digits = 4), round(x$est$skmSE2, 
            digits = 4)), ncol = 3)
        write("                         Survival", file = "")
        write("     Subpopulation     Probability      Std. Err.", 
            file = "")
        for (i in 1:x$sp$nsubpop) {
            write(paste(format(temp[i, 1], width = 12), format(temp[i, 
                2], width = 19, nsmall = 4), format(temp[i, 3], width = 15, 
                nsmall = 4)), file = "")
        }
        write(paste("        Overall", format(round(x$est$overallSkmObs2, 
            digits = 4), nsmall = 4, width = 16), format(round(x$est$overallSkmSE2, 
            digits = 4), nsmall = 4, width = 15)), file = "")
        cat("\n")
        write(paste("Survival differences at time point", timePoint), 
            file = "")
        temp <- matrix(c(1:x$sp$nsubpop, round(x$est$skmObs1 - x$est$skmObs2, digits = 4), 
          round(sqrt(x$est$skmSE1^2 + x$est$skmSE2^2), digits = 4)), ncol = 3)
        write("                         Survival", file = "")
        write("     Subpopulation      Difference      Std. Err.", 
            file = "")
        for (i in 1:x$sp$nsubpop) {
            write(paste(format(temp[i, 1], width = 12), format(temp[i, 
                2], width = 19, nsmall = 4), format(temp[i, 3], width = 15, 
                nsmall = 4)), file = "")
        }
        write(paste("        Overall", format(round(x$est$overallSkmObs1 - 
            x$est$overallSkmObs2, digits = 4), nsmall = 4, width = 16), 
            format(round(sqrt(x$est$overallSkmSE1^2 + x$est$overallSkmSE2^2), 
                digits = 4), nsmall = 4, width = 15)), file = "")
        cat("\n")
        write("Hazard ratio estimates", file = "")
        temp <- matrix(c(1:x$sp$nsubpop, round(x$est$logHR, digits = 6), round(x$est$logHRSE, 
            digits = 6), round(exp(x$est$logHR), digits = 2)), ncol = 4)
        write("     Subpopulation        Log HR       Std. Err.       Hazard Ratio", 
            file = "")
        for (i in 1:x$sp$nsubpop) {
            write(paste(format(temp[i, 1], width = 12), format(temp[i, 
                2], width = 19, nsmall = 6), format(temp[i, 3], width = 14, 
                nsmall = 6), format(temp[i, 4], width = 15, nsmall = 2)), 
                file = "")
        }
        write(paste("        Overall", format(round(x$est$overallLogHR, 
            digits = 6), nsmall = 6, width = 16), format(round(x$est$overallLogHRSE, 
            digits = 6), nsmall = 6, width = 14), format(round(exp(x$est$overallLogHR), 
            digits = 2), nsmall = 2, width = 15)), file = "")
        cat("\n")
      }
      else {
      #
      #   print out the cumulative incidence results 
      #
      cat("\n")
      write(paste("Cumulative incidence estimates for treatment group",x$arg$trts[1],"at time point",timePoint),file="")
      temp <- matrix(c(1:x$sp$nsubpop,round(x$est$skmObs1,digits=4),round(x$est$skmSE1,digits=4)),ncol=3)
      write("                        Cumulative",file="")
      write("     Subpopulation      Incidence        Std. Err.",file="")
      for (i in 1:x$sp$nsubpop) {
          write(paste(format(temp[i,1],width=12),format(temp[i,2],width=19,nsmall=4),
              format(temp[i,3],width=15,nsmall=4)),file="")
      }
      write(paste("        Overall",format(round(x$est$overallSkmObs1,digits=4),nsmall=4,width=16),
          format(round(x$est$overallSkmSE1,digits=4),nsmall=4,width=15)),file="")
      cat("\n")
      write(paste("Cumulative incidence estimates for treatment group",x$arg$trts[2],"at time point",timePoint),file="")
      temp <- matrix(c(1:x$sp$nsubpop,round(x$est$skmObs2,digits=4),round(x$est$skmSE2,digits=4)),ncol=3)
      write("                        Cumulative",file="")
      write("     Subpopulation      Incidence        Std. Err.",file="")
      for (i in 1:x$sp$nsubpop) {
          write(paste(format(temp[i,1],width=12),format(temp[i,2],width=19,nsmall=4),
              format(temp[i,3],width=15,nsmall=4)),file="")
      }
      write(paste("        Overall",format(round(x$est$overallSkmObs2,digits=4),nsmall=4,width=16),
          format(round(x$est$overallSkmSE2,digits=4),nsmall=4,width=15)),file="")
      cat("\n")
      write(paste("Cumulative incidence differences at time point",timePoint),file="")
      temp <- matrix(c(1:x$sp$nsubpop,round(x$est$skmObs1-x$est$skmObs2,digits=4),round(sqrt(x$est$skmSE1^2+x$est$skmSE2^2),digits=4)),ncol=3)
      write("                        Cumulative",file="")
      write("                        Incidence",file="")
      write("     Subpopulation      Difference       Std. Err.",file="")
      for (i in 1:x$sp$nsubpop) {
          write(paste(format(temp[i,1],width=12),format(temp[i,2],width=19,nsmall=4),
              format(temp[i,3],width=15,nsmall=4)),file="")
      }
      write(paste("        Overall",
          format(round(x$est$overallSkmObs1-x$est$overallSkmObs2,digits=4),nsmall=4,width=16),
          format(round(sqrt(x$est$overallSkmSE1^2+x$est$overallSkmSE2^2),digits=4),nsmall=4,width=15)),file="")
      cat("\n")
      write("Hazard ratio estimates",file="")
      temp <- matrix(c(1:x$sp$nsubpop,round(x$est$logHR,digits=6),round(x$est$logHRSE,digits=6),round(exp(x$est$logHR),digits=2)),ncol=4)
      write("     Subpopulation        Log HR       Std. Err.       Hazard Ratio",file="")
      for (i in 1:x$sp$nsubpop) {
          write(paste(format(temp[i,1],width=12),format(temp[i,2],width=19,nsmall=6),
              format(temp[i,3],width=14,nsmall=6),format(temp[i,4],width=15,nsmall=2)),file="")
      }
      write(paste("        Overall",format(round(x$est$overallLogHR,digits=6),nsmall=6,width=16),
          format(round(x$est$overallLogHRSE,digits=6),nsmall=6,width=14),
          format(round(exp(x$est$overallLogHR),digits=2),nsmall=2,width=15)),file="")
      cat("\n")
      }	
    }

#
#   2. covariance matrices
#
    if (cov & x$arg$noperm > 0){
      cat("\n")
      if (x$arg$trttype == "KM")
        write(paste("The covariance matrix of the Kaplan-Meier differences at", 
              x$arg$timest, "time units for the", x$sp$nsubpop, "subpopulations is:"), 
              file = "")
      else
        write(paste("The covariance matrix of the cumulative incidence differences at", 
              x$arg$timest, "time units for the", x$sp$nsubpop, "subpopulations is:"), 
              file = "")

      print(x$res$sigma)

      cat("\n")
      write(paste("The covariance matrix of the log hazard ratios for the", 
            x$sp$nsubpop, "subpopulations is:"), file = "")
      print(x$res$HRsigma)
      cat("\n")
    }
#
#   3. Supremum test and Chi-square test results
#
    if (test & x$arg$noperm > 0){
      cat("\n")
      write(paste("Supremum test results"), file = "")
      if (x$arg$trttype == "KM")
        write(paste("Interaction P-value based on Kaplan-Meier estimates :", x$res$pvalue), file = "")
      else
        write(paste("Interaction P-value based on cumulative incidence estimates :", x$res$pvalue), file = "")
      write(paste("Interaction P-value based on hazard ratio estimates :", x$res$HRpvalue), file = "")
      cat("\n")
      write(paste("Chi-square test results"), file = "")
      if (x$arg$trttype == "KM")
        write(paste("Interaction P-value based on Kaplan-Meier estimates :", 
              x$res$chi2pvalue), file = "")
      else
        write(paste("Interaction P-value based on cumulative incidence estimates :", 
              x$res$chi2pvalue), file = "")
      write(paste("Interaction P-value based on hazard ratio estimates :", 
            x$res$HRchi2pvalue), file = "")
      cat("\n")
    }
}



#####################################################################
#
# Perform a stepp analysis with KM/HR and create a stepp object
#
#####################################################################
stepp <- function(
		  trttype,	# type of treatment effect
				#   "KM" or "CI"
    		  coltrt,	# treatment variable
    		  coltime,	# time to event variable (must be zero or positive)
    		  colcens=0,	# censoring variable (for KM only)
		  coltype=0,	# codes for different causes of failure (for CI only)
				#   0 for censored observations
				#   1 for event of interest
				#   2 for other causes of failure
		  colvar,	# covariate of interest
		  trts,		# a vector containing the codes for the 2 treatment arms
				#   1st and 2nd arms respectively
		  patspop,	# larger parameter (r2) for subpopulation construction
				#   that determines how many patients are in each subpopulation
		  minpatspop,	# smaller parameter (r1) for subpopulation construction
				#   that determines the largest number of patients in common
				#   among consecutive subpopulations
		  timest,	# timepoint to estimate survival
		  noperm	# the desired number of permutations; must be greater than zero
	)
{
#   include library of survival functions
    library(survival)

# internal functions
# +++++
version.stepp <- function()
{
    # return the current version number of STEPP
    # May 17, 2011
    "2.3-2"
}

cuminc.HR <- function (ftime, fstatus, group, strata, rho = 0, cencode = 0, 
    subset, na.action = na.omit) 
{
#
    library(cmprsk)
    d <- data.frame(time = ftime, cause = fstatus, group = as.factor(if (missing(group)) 
        rep(1, length(ftime))
    else group), strata = as.factor(if (missing(strata)) 
        rep(1, length(ftime))
    else strata))
    if (!missing(subset)) 
        d <- d[subset, ]
    tmp <- nrow(d)
    d <- na.action(d)
    if (nrow(d) != tmp) 
        cat(format(tmp - nrow(d)), "cases omitted due to missing values\n")
    no <- nrow(d)
    cg <- "  "
    nst <- length(levels(d$strata))
    d <- d[order(d$time), ]
    ugg <- table(d$group)
    d$group <- factor(d$group, names(ugg)[ugg > 0])
    ugg <- levels(d$group)
    censind <- ifelse(d$cause == cencode, 0, 1)
    uc <- table(d$cause[censind == 1])
    if (is.factor(d$cause)) 
        uclab <- names(uc)[uc > 0]
    else uclab <- as.numeric(names(uc)[uc > 0])
    nc <- length(uclab)
    ng <- length(ugg)
    if (ng > 1) {
        ng1 <- ng - 1
        ng2 <- ng * ng1/2
        v <- matrix(0, nrow = ng1, ncol = ng1)
        storage.mode(v) <- "double"
        vt <- double(ng2)
        s <- double(ng1)
    }
    pf <- vector("list", ng * nc)
    stat <- double(nc)
    l <- 0
    for (ii in 1:nc) {
        causeind <- ifelse(d$cause == uclab[ii], 1, 0)
        for (jj in 1:length(ugg)) {
            cg <- c(cg, paste(ugg[jj], uclab[ii]))
            l <- l + 1
            cgind <- d$group == ugg[jj]
            ncg <- length(cgind[cgind])
            n2 <- length(unique(d$time[cgind & causeind == 1]))
            n2 <- 2 * n2 + 2
            tmp <- double(n2)
            z <- .Fortran("cinc", as.double(d$time[cgind]), as.integer(censind[cgind]), 
                as.integer(causeind[cgind]), as.integer(ncg), 
                x = tmp, f = tmp, v = tmp, PACKAGE = "cmprsk")
            pf[[l]] <- list(time = z$x, est = z$f, var = z$v)
        }
        if (ng > 1) {
            causeind <- 2 * censind - causeind
            z2 <- .Fortran("crstm", as.double(d$time), as.integer(causeind), 
                as.integer(d$group), as.integer(d$strata), as.integer(no), 
                as.double(rho), as.integer(nst), as.integer(ng), 
                s, v, as.double(d$time), as.integer(causeind), 
                as.integer(d$group), vt, s, vt, double((4 + 3 * 
                  ng) * ng), integer(4 * ng), PACKAGE = "cmprsk")
            stat[ii] <- -1
            a <- qr(z2[[10]])
            if (a$rank == ncol(a$qr)) {
                b <- diag(dim(a$qr)[1])
                stat[ii] <- z2[[9]] %*% qr.coef(a, b) %*% z2[[9]]
            }
#
            if (ii == 1) {
               ome <- (-1)*z2[[9]]
               omevar <- z2[[10]]
            }
        }
    }
    names(pf) <- cg[2:length(cg)]
    if (ng > 1) {
        names(stat) <- uclab
        stat <- list(Tests = cbind(stat = stat, pv = 1 - pchisq(stat, 
            ng - 1), df = rep(ng - 1, length(stat))))
        pf <- c(pf, stat)
#
        omeres <- list(ome=ome,omevar=omevar)
        pf <- c(pf,omeres)
    }
    pf
}


kmest1<-function(y,m,n,ndf,t,s,v,ntpt,tpt,nrr,ndd)
{
  # initialization 
  f<-1.0
  kr<-n
  nrr[1]<-n
  ndd[1]<-0

  for (i in 2:ntpt){
    nrr[i]<-0
    ndd[i]<-0
  }

  ltp<-1
  var<-0.0
  l<-1
  t[1]<-0
  s[1]<-1
  v[1]<-0
  i<-1

  # main loop 
  while (i<= n) {
    k<-i+1
    k2<-0

    while (k2<=0) {
      if (k > n) k2<-1
      else {
	if (y[k] != y[i]) k2<-1 else k<-k+1
      }
	
    } # end while 

    k<-k-1
    nd<-0

    for (j in i:k) nd<-nd+m[j]

    while (ltp<=ntpt && y[i]>tpt[ltp+1]) {
      ltp<-ltp+1
      nrr[ltp]<-kr
    }  # end while 

    ndd[ltp]<-ndd[ltp]+nd
    if (nd>0) {
      t1<- nd / kr  
      f<-f*(1-t1)
      if (nd<kr) var<-var+t1/(kr-nd)
      t[l+1]<-y[i]
      s[l+1]<-s[l]
      v[l+1]<-v[l]
      l<-l+2
      t[l]<-y[i]
      s[l]<-f
      v[l]<-var*f*f
    } # end if 

    i<-k+1
    kr<-n-k     # kr<-n-k-1
    k<-i
  } # end while - main loop 

  l<-l+1
  t[l]<-y[n]
  s[l]<-s[l-1]
  v[l]<-v[l-1]

  # return a list of all the arguments
  kmtest1<-list(y=y,m=m,n=n,ndf=ndf,t=t,s=s,v=v,ntpt=ntpt,tpt=tpt,nrr=as.integer(nrr),ndd=as.integer(ndd))
  kmtest1
  
}


kmest <-
function (time,status,group,tpt,pv=TRUE,pv.strat,pv.sub,rho=0,subset,na.action=na.omit) 
{
  d <- data.frame(time=time,status=status,
    group=as.factor(if (missing(group)) rep(1,length(time)) else group),
    pv.strat=as.factor(if (missing(pv.strat)) rep(1,length(time)) else pv.strat),
    pv.sub = as.factor(if (missing(pv.sub)) rep(TRUE, length(time)) else pv.sub))
  if (!missing(subset)) d <- d[subset,]
  tmp <- nrow(d)
  d <- na.action(d)
  if (nrow(d) != tmp) cat(format(tmp-nrow(d)),'cases omitted due to missing values\n')
  no <- nrow(d)
  if (any(d$status != 0 & d$status != 1)) stop("invalid status values")
  T1 <- table(d$group)
  subgl <- T1>0
  T8 <- d$group[d$status==1]
  if (length(T8) > 0) {
    T8 <- cbind(T1,table(d$group[d$status==1]))
  } else {
    T8 <- cbind(T1,rep(0,length(T8)))
  }
  T8 <- T8[subgl,,drop=FALSE]
  T1 <- names(T1)[subgl]
  time <- ifelse(d$time<=0,.00001,d$time)
  if (missing(tpt)) {
    tpt <- pretty(time)}
  else {
    ym <- round(max(time),2)
    tpt <- c(0,tpt,ym)
  }
  ntpt <- length(tpt)-1
  lev <- vector("character",ntpt)
  for (i in 1:ntpt) lev[i] <- paste(format(tpt[i]),format(tpt[i+1]),sep="-")
  nrd <- rep(0, ntpt)
  o <- order(time)
  Tl <- length(T1)
  z <- as.list(1:Tl)
  for(i in 1:Tl) {
    Ty <- (time[o])[d$group[o]==T1[i]]
    Tm <- (d$status[o])[d$group[o]==T1[i]]
    ndf <- length(unique(Ty[Tm==1]))
    t <- double(2*ndf+2)

    a <- kmest1(Ty, Tm, length(Ty), ndf, t,t,t,ntpt, tpt, nrd, nrd)

    tt <- paste(a[[11]],a[[10]],sep="/")
    names(tt) <- lev
    z[[i]] <- list(time=a[[5]],est=a[[6]],var=a[[7]],tint=tt,nnd=T8[i,])
    }
  names(z) <- T1
  if (pv & Tl > 1) {
    stop("internal error")
    # pv <- logrank(time=d$time,status=d$status,group=d$group,strata=d$pv.strat,rho=rho,subset=(d$pv.sub=='TRUE'))$pv
    # attr(z,"pv") <- pv
  }

  class(z) <- "kmest"
  z
}




#
tpest1<-function(x,n,ind,tp5,ntp)
{
  l <- ntp
  for (i in rev(1:ntp)){
    if (x[n] >= tp5[i]) break;
    ind[l]<-0
    l <- l - 1
  }

  if (l <=0) return (ind);

  if (x[n] == tp5[l]){
    ind[l] <- n
    l <- l - 1
  }

  # assuming unique values sorted in ascending order
  k <- n-1
  loop <- TRUE
  while(loop){
    if (l <= 0) return (ind);

   loop <- FALSE
    for (i in rev(1:k)) {
      if (x[k] <= tp5[l]){
        ind[l] <- k+1
        l <- l - 1
        loop <- TRUE
        break #out of the for loop
      } else 
      k <- k-1
    } # end for loop
  } # end while loop

  # error in the following loop corrected 9-28-04
  for (i in 1:l) ind[i] <- 0
  return (ind);	
}


tpest <-
function(w,times) 
{
  if (!is.null(w$Tests)) w <- w[names(w) != 'Tests']
  ng <- length(w)
  times <- sort(unique(times))
  nt <- length(times)
  storage.mode(times) <- "double"
  storage.mode(nt) <- "integer"
  ind <- matrix(0,ncol=nt,nrow=ng)
  oute <- matrix(NA,ncol=nt,nrow=ng)
  outv <- oute
  storage.mode(ind) <- "integer"
  slct <- rep(TRUE,ng)
  for (i in 1:ng) {
    if (is.null((w[[i]])$est)) { slct[i] <- FALSE} else { 
       z1 <- as.integer(tpest1(w[[i]][[i]], length(w[[i]][[i]]), ind[i,], times, nt))
       ind[i,] <- z1
       oute[i, ind[i,]>0] <- w[[i]][[2]][z1]
       if (length(w[[i]])>2) outv[i,ind[i,]>0] <- w[[i]][[3]][z1]
	}
    }
  
  dimnames(oute) <- list(names(w)[1:ng],as.character(times))
  dimnames(outv) <- dimnames(oute)
  list(est=oute[slct,,drop=FALSE],var=outv[slct,,drop=FALSE])
}

#
# internal functions end
#
# Main Program starts here
#
    print(paste("STEPP Version ", version.stepp()))

#   initialize parameters
#   ----------------------
    t0 <- proc.time()
    stime <- Sys.time()
#
#   create input object - to be save later
#
    input.obj <- list(trttype	 = trttype,
		      coltrt     = coltrt,
                      coltime    = coltime,
                      colcens    = colcens,
		      coltype    = coltype,
                      colvar     = colvar,	
                      trts       = trts,
                      patspop    = patspop,
		      minpatspop = minpatspop,
                      timest     = timest,
                      noperm     = noperm)

    timePoint <- timest
    r1 <- minpatspop
    r2 <- patspop
    noPerms <- noperm
    survTime <- coltime
    censor <- colcens
    covariate <- colvar
    type <- coltype
#
#   Perform error check on input to make sure it is reasonable:
#   -----------------------------------------------------------
#   1. check that there are two treatments
#
    if (length(trts) != 2) stop("Two Treatment Arms Are Needed !")
#
#   2. check r1 and r2
#
    if (r1 >= r2) stop("minpatspop MUST be less than patspop")
#
#   3. check noperms
    if (noperm < 0) stop("no. of permutation must be greater than or equal to zero")
#
#   4. check trttype
#
    if (trttype != "KM" & trttype !="CI") stop ("Invalid Treatment Type !")
#
#   5. check if time is negative
#
    for (ii in 1:length(survTime)){
	if (!is.na(survTime[ii])){
	  if (survTime[ii] < 0) stop ("Time to event must be >= 0 !")
	}
    }    

#
#   Error checking done - real work start here
#   ------------------------------------------
#
#   create 0/1 treatment assignment 
#
    txassign <- rep(NA, length(coltrt))
    txassign[which(coltrt == trts[1])] <- 1
    txassign[which(coltrt == trts[2])] <- 0
#
#   set up the data
#
    if (trttype == "KM") indata <- cbind(txassign, survTime, censor, covariate)
    else indata <- cbind(txassign, survTime, type, covariate)
    indata <- indata[apply(indata, 1, function(x) !any(is.na(x))),, drop = FALSE]
    txassign <- indata[,1]
    survTime <- indata[,2]
    if (trttype == "KM") censor <- indata[, 3]
    else type <- indata[,3]
    covariate <- indata[,4]

#   get the overlapping subgroups
#
    zvals <- rep(0, length(covariate))
    absfreq <- rep(0, length(covariate))
    sortedz <- sort(covariate)
    zvals[1] <- sortedz[1]
    j <- 1
    absfreq[1] <- 1
    for (i in 2:length(covariate)) {
        if (sortedz[i] != zvals[j]) {
            j <- j + 1
            zvals[j] <- sortedz[i]
            absfreq[j] <- 1
        }
        else {
            absfreq[j] <- absfreq[j] + 1
        }
    }
    zvals <- zvals[1:j]
    absfreq <- absfreq[1:j]
    cumfreq <- absfreq
    for (i in 2:length(cumfreq)) cumfreq[i] <- cumfreq[i] + cumfreq[i - 1]
    I0 <- rep(0, 15)
    I1 <- rep(0, 15)
    I0[1] <- 1
    I1[1] <- sum(cumfreq < r2) + 1
    stopflag <- 0
    nsubpop <- 2
    while (stopflag == 0) {
        indinf <- I0[nsubpop - 1] + 1
        while ((cumfreq[I1[nsubpop - 1]] - cumfreq[indinf - 1]) > r1) {
            indinf <- indinf + 1
        }
        I0[nsubpop] <- indinf
        indsup <- I1[nsubpop - 1]
        while (((cumfreq[indsup] - cumfreq[I0[nsubpop]] + absfreq[I0[nsubpop]]) < r2) && (stopflag == 0)) {
            indsup <- indsup + 1
            stopflag <- 1 * (indsup == length(zvals))
        }
        I1[nsubpop] <- indsup
        nsubpop <- nsubpop + 1
    }
    nsubpop <- nsubpop - 1
    npatsub <- rep(0, nsubpop)
    minz <- rep(NA, nsubpop)
    maxz <- rep(NA, nsubpop)
    npatsub[1] <- cumfreq[I1[1]]
    for (i in 2:nsubpop) npatsub[i] <- cumfreq[I1[i]] - cumfreq[I0[i] - 1]
    I0 <- I0[1:nsubpop]
    I1 <- I1[1:nsubpop]
    npats <- length(covariate)
    subpop <- matrix(rep(0, (npats * nsubpop)), ncol = nsubpop)
    medians <- rep(0, nsubpop)
    for (i in 1:nsubpop) {
        subpop[, i] <- (covariate >= zvals[I0[i]]) * (covariate <= zvals[I1[i]])
        medians[i] <- round((median(covariate[subpop[, i] == 1])), digits = 2)
        minz[i] <- round(zvals[I0[i]],digits=4)
        maxz[i] <- round(zvals[I1[i]],digits=4)
    }
#
#   create subpopulation object - to be saved
#
    sp.obj<-list(nsubpop = nsubpop,
                 subpop  = subpop,
		 npatsub = npatsub,
                 medians = medians,
                 minz    = minz,
                 maxz    = maxz)

    if (trttype == "KM"){
	#
	#   compute survival estimates with KM/HR
	#
    	skmObs1 <- rep(0, nsubpop)
    	skmObs2 <- rep(0, nsubpop)
    	skmSE1 <- rep(0, nsubpop)
    	skmSE2 <- rep(0, nsubpop)
    	logHR <- rep(0,nsubpop)
    	logHRSE <- rep(0,nsubpop)
    	for (i in 1:nsubpop) {
          skmObs1[i] <- max(tpest(kmest(survTime[txassign==1 & subpop[,i]==1],
              censor[txassign==1 & subpop[,i]==1]),timePoint)$est,0)
          skmObs2[i] <- max(tpest(kmest(survTime[txassign==0 & subpop[,i]==1],
              censor[txassign==0 & subpop[,i]==1]),timePoint)$est,0)
          skmSE1[i] <- sqrt(tpest(kmest(survTime[txassign==1 & subpop[,i]==1],
              censor[txassign==1 & subpop[,i]==1]),timePoint)$var)
          skmSE2[i] <- sqrt(tpest(kmest(survTime[txassign==0 & subpop[,i]==1],
              censor[txassign==0 & subpop[,i]==1]),timePoint)$var)
          LogRank <- survdiff(Surv(survTime[subpop[,i]==1],censor[subpop[,i]==1])~txassign[subpop[,i]==1])
          logHR[i] <- -(LogRank$obs[1]-LogRank$exp[1])/LogRank$var[1,1]
          logHRSE[i] <- sqrt(1/LogRank$var[1,1])
    	}
    	overallSkmObs1 <- max(tpest(kmest(survTime[txassign==1],censor[txassign==1]),timePoint)$est, 0)
    	overallSkmObs2 <- max(tpest(kmest(survTime[txassign==0],censor[txassign==0]),timePoint)$est, 0)
    	overallSkmSE1 <- sqrt(tpest(kmest(survTime[txassign==1],censor[txassign==1]),timePoint)$var)
    	overallSkmSE2 <- sqrt(tpest(kmest(survTime[txassign==0],censor[txassign==0]),timePoint)$var)
    	LogRank <- survdiff(Surv(survTime,censor) ~ txassign)
    	overallLogHR <- -(LogRank$obs[1]-LogRank$exp[1])/LogRank$var[1,1]
    	overallLogHRSE <- sqrt(1/LogRank$var[1,1])
    	if (sum(is.na(skmObs1)) != 0 | sum(is.na(skmObs2)) != 0 | 
            is.na(overallSkmObs1) != FALSE | is.na(overallSkmObs2) != FALSE) {
          	cat("\n")
        	print(paste("Unable to estimate the survival time at ", 
            		timePoint, " time-unit(s) within a subpopulation(s)."))
        	print(paste("Consider creating larger subpopulations or selecting a different timepoint for estimation."))
        	stop()
    	}
	#
	#   create the estimates (KM & HR) object - to be saved
	#
    	est.obj <- list(skmObs1        = skmObs1,
        	        skmSE1         = skmSE1,
                        overallSkmObs1 = overallSkmObs1,
                        overallSkmSE1  = overallSkmSE1,
                        skmObs2        = skmObs2,
                        skmSE2         = skmSE2,
                        overallSkmObs2 = overallSkmObs2,
                        overallSkmSE2  = overallSkmSE2,
                        logHR          = logHR,
                        logHRSE        = logHRSE,
                        overallLogHR   = overallLogHR,
                        overallLogHRSE = overallLogHRSE)

	#
	#   do the permutations
	#
    	pvalue <- NA
    	differences <- matrix(rep(0, (noPerms * nsubpop)), ncol = nsubpop)
    	logHRs <- matrix(rep(0, (noPerms * nsubpop)), ncol = nsubpop)
    	tPerm <- rep(0, noPerms)
    	no <- 0
    	p <- 0
 
   	terminate <- 0
    	Ntemp <- nrow(subpop)
    	IndexSet1 <- (1:Ntemp)[txassign == 1]
    	IndexSet2 <- (1:Ntemp)[txassign == 0]
    	cat("\nComputing the pvalue\n")
    	while (no < noPerms) {
      	  ## PERMUTE THE VECTOR OF SUBPOPULATIONS WITHIN TREATMENTS ##
      	  if (no%%100 == 0) cat(no*100/noPerms, "% completed\n")

      	  Subpop <- as.matrix(subpop)
      	  permuteSubpop <- matrix(0, nrow = Ntemp, ncol = nsubpop)
      	  permuteSubpop[txassign == 1] <- Subpop[sample(IndexSet1),]
      	  permuteSubpop[txassign == 0] <- Subpop[sample(IndexSet2),]
      	  subpop <- permuteSubpop
      	  skm1 <- rep(0, nsubpop)
      	  skm2 <- rep(0, nsubpop)
      	  slogHRs <- rep(0, nsubpop)
      	  for (i in 1:nsubpop) {
        	skm1[i] <- max(tpest(kmest(survTime[txassign == 1 & subpop[, i] == 1], censor[txassign == 1 & subpop[, i] == 1]), timePoint)$est, 0)
        	skm2[i] <- max(tpest(kmest(survTime[txassign == 0 & subpop[, i] == 1], censor[txassign == 0 & subpop[, i] == 1]), timePoint)$est, 0)
        	LogRank <- survdiff(Surv(survTime[subpop[,i] == 1],censor[subpop[,i] == 1]) ~ txassign[subpop[,i] == 1])
        	slogHRs[i] <- -(LogRank$obs[1]-LogRank$exp[1])/LogRank$var[1,1]
      	  }
      	  overallSkm1 <- max(tpest(kmest(survTime[txassign == 1], censor[txassign == 1]), timePoint)$est, 0)
      	  overallSkm2 <- max(tpest(kmest(survTime[txassign == 0], censor[txassign == 0]), timePoint)$est, 0)
      	  LogRank <- survdiff(Surv(survTime,censor) ~ txassign)
      	  overallSlogHR <- -(LogRank$obs[1]-LogRank$exp[1])/LogRank$var[1,1]
      	  if (sum(is.na(skm1)) == 0 & sum(is.na(skm2)) == 0 & is.na(overallSkm1) == FALSE & is.na(overallSkm2) == FALSE) {
        	no <- no + 1
        	p <- p + 1
        	for (s in 1:nsubpop) {
          	  differences[p, s] <- (skm1[s] - skm2[s]) - (overallSkm1 - overallSkm2)
          	  logHRs[p,s] <- slogHRs[s]-overallSlogHR
        	}
      	  }
      	  terminate <- terminate + 1
      	  if (terminate >= noPerms + 10000) {
        	print(paste("After permuting ", noPerms, "plus 10000, or ", 
        	      noPerms + 10000, " times, the program is unable to generate the permutation distribution based on ", 
        	      noPerms, "permutations of the data"))
        	print(paste("Consider creating larger subpopulations or selecting a different timepoint for estimation"))
        	stop()
      	  }
    	}
    	cat("100% completed\n")
    }
    else {
	#
	#   compute cumulative incidence estimates
	#
    	ObsCI1 <- rep(NA,nsubpop)
    	ObsCI2 <- rep(NA,nsubpop)
    	ObsCISE1 <- rep(NA,nsubpop)
    	ObsCISE2 <- rep(NA,nsubpop)
    	overallObsCI1 <- NA
    	overallObsCISE1 <- NA
    	overallObsCI2 <- NA
    	overallObsCISE2 <- NA
    	logHR <- rep(0,nsubpop)
    	logHRSE <- rep(0,nsubpop)
    	for (i in 1:nsubpop) {
            result <- cuminc.HR(survTime[subpop[,i]==1],type[subpop[,i]==1],txassign[subpop[,i]==1])
            if (max(result$"1 1"$time) >= timePoint) {
                index <- sum(result$"1 1"$time <= timePoint)
                ObsCI1[i] <- result$"1 1"$est[index]
                ObsCISE1[i] <- sqrt(result$"1 1"$var[index])
            }
            if (max(result$"0 1"$time) >= timePoint) {
                index <- sum(result$"0 1"$time <= timePoint)
                ObsCI2[i] <- result$"0 1"$est[index]
                ObsCISE2[i] <- sqrt(result$"0 1"$var[index])
            }
            logHR[i] <- result$ome/result$omevar
            logHRSE[i] <- sqrt(1/result$omevar)
    	}
    	result <- cuminc.HR(survTime,type,txassign)
    	if (max(result$"1 1"$time) >= timePoint) {
            index <- length(result$"1 1"$time[result$"1 1"$time <= timePoint])
            overallObsCI1 <- result$"1 1"$est[index]
            overallObsCISE1 <- sqrt(result$"1 1"$var[index])
    	}
    	if (max(result$"0 1"$time) >= timePoint) {
            index <- sum(result$"0 1"$time <= timePoint)
            overallObsCI2 <- result$"0 1"$est[index]
            overallObsCISE2 <- sqrt(result$"0 1"$var[index])
        }
    	overallLogHR <- result$ome/result$omevar
    	overallLogHRSE <- sqrt(1/result$omevar)
    	if (sum(is.na(ObsCI1)) != 0 | sum(is.na(ObsCI2)) != 0 | 
            is.na(overallObsCI1) != FALSE | is.na(overallObsCI2) != FALSE) {
              cat("\n")
              print(paste("Unable to estimate the survival time at ", 
                    timePoint, " time-unit(s) within a subpopulation(s)."))
              print(paste("Consider creating larger subpopulations or selecting a different timepoint for estimation."))
              stop()
    	}
	#
	#   create the estimates (Cumulative Incidence) object - to be saved
	#
    	est.obj <- list(skmObs1          = ObsCI1,
        	        skmSE1           = ObsCISE1,
                        overallSkmObs1   = overallObsCI1,
                        overallSkmSE1    = overallObsCISE1,
                        skmObs2          = ObsCI2,
                        skmSE2           = ObsCISE2,
                        overallSkmObs2   = overallObsCI2,
                        overallSkmSE2    = overallObsCISE2,
                        logHR            = logHR,
                        logHRSE          = logHRSE,
                        overallLogHR     = overallLogHR,
                        overallLogHRSE   = overallLogHRSE)


	#
	#   do the permutations
	#
    	pvalue <- NA
       	differences <- matrix(rep(0, (noPerms * nsubpop)), ncol = nsubpop)
       	logHRs <- matrix(rep(0, (noPerms * nsubpop)), ncol = nsubpop)
       	tPerm <- rep(0, noPerms)
       	no <- 0
       	p <- 0
       	terminate <- 0
        Ntemp <- nrow(subpop)
        IndexSet1 <- (1:Ntemp)[txassign == 1]
        IndexSet2 <- (1:Ntemp)[txassign == 0]
	cat("\nComputing the pvalue\n")
        while (no < noPerms) {
	  ## PERMUTE THE VECTOR OF SUBPOPULATIONS WITHIN TREATMENTS ##
	  if (no%%100 == 0) cat(no*100/noPerms, "% completed\n")

          Subpop <- as.matrix(subpop)
          permuteSubpop <- matrix(0, nrow = Ntemp, ncol = nsubpop)
          permuteSubpop[txassign == 1] <- Subpop[sample(IndexSet1),]
          permuteSubpop[txassign == 0] <- Subpop[sample(IndexSet2),]
          subpop <- permuteSubpop
          sCI1 <- rep(NA,nsubpop)
          sCI2 <- rep(NA,nsubpop)
          overallsCI1 <- NA
          overallsCI2 <- NA
          slogHRs <- rep(0, nsubpop)
          for (i in 1:nsubpop) {
              result <- cuminc.HR(survTime[subpop[,i]==1],type[subpop[,i]==1],txassign[subpop[,i]==1])
              if (max(result$"1 1"$time) >= timePoint) {
                  index <- sum(result$"1 1"$time <= timePoint)
                  sCI1[i] <- result$"1 1"$est[index]
              }
              if (max(result$"0 1"$time) >= timePoint) {
                  index <- sum(result$"0 1"$time <= timePoint)
                  sCI2[i] <- result$"0 1"$est[index]
              }
              slogHRs[i] <- result$ome/result$omevar
          }
          result <- cuminc.HR(survTime,type,txassign)
          if (max(result$"1 1"$time) >= timePoint) {
              index <- length(result$"1 1"$time[result$"1 1"$time <= timePoint])
              overallsCI1 <- result$"1 1"$est[index]
          }
          if (max(result$"0 1"$time) >= timePoint) {
              index <- sum(result$"0 1"$time <= timePoint)
              overallsCI2 <- result$"0 1"$est[index]
          }
          overallSlogHR <- result$ome/result$omevar
          if (sum(is.na(sCI1)) == 0 & sum(is.na(sCI2)) == 0 & is.na(overallsCI1) == FALSE & is.na(overallsCI2) == FALSE) {
              no <- no + 1
              p <- p + 1
              for (s in 1:nsubpop) {
                  differences[p,s] <- (sCI1[s]-sCI2[s])-(overallsCI1-overallsCI2)
                  logHRs[p,s] <- slogHRs[s]-overallSlogHR
              }
          }
          terminate <- terminate + 1
          if (terminate >= noPerms + 10000) {
              print(paste("After permuting ", noPerms, "plus 10000, or ", 
                noPerms + 10000, " times, the program is unable to generate the permutation distribution based on ", 
                noPerms, "permutations of the data"))
              print(paste("Consider creating larger subpopulations or selecting a different timepoint for estimation"))
              stop()
            }
      }
    }

    if (noperm != 0)
    {	
    #
    #   generate covariance matrix, p-value based on KM
    #

      sigma <- var(differences)

    #   catch the error if we encounter a singularity problem when inverting the matrix sigma
    #
      tryCatch({
        sigmainv <- solve(sigma)
        }, error = function(ex){
		cat(geterrmessage())
		cat("\n")
		set.seed<-4593432
		stop("STEPP has encountered an unexpected problem. Please retry. If problem persists, try different values of minpatspop(r1) and patspop(r2).")
	        }
      )
      sigma1<-sigma

      chi2perm <- rep(NA,noPerms)
      for (i in 1:noPerms){
        res <- matrix(differences[i,],ncol=1)
        chi2perm[i] <- t(res)%*%sigmainv%*%res
      }
      obsDifferences <- matrix(rep(NA, nsubpop),ncol=1)

      if (trttype == "KM"){
        for (o in 1:nsubpop)
	  obsDifferences[o,1] <- (skmObs1[o] - skmObs2[o]) - (overallSkmObs1 - overallSkmObs2)
      }
      else{
        for (o in 1:nsubpop)
	  obsDifferences[o,1] <- (ObsCI1[o] - ObsCI2[o]) - (overallObsCI1-overallObsCI2)
      }
      chi2obs <- c(t(obsDifferences)%*%sigmainv%*%obsDifferences)
      chi2pvalue <- sum(ifelse(chi2perm > chi2obs, 1, 0))/noPerms
      signif(chi2pvalue, 6)
#
      sigma <- sqrt(diag(var(differences)))
      stdDifferences <- t(apply(differences, 1, function(x) x/sigma))
      tPerm <- apply(abs(stdDifferences), 1, max)
      obsDifferences <- matrix(rep(0, nsubpop))

      if (trttype == "KM"){
        for (o in 1:nsubpop)
       	obsDifferences[o] <- (skmObs1[o] - skmObs2[o]) - (overallSkmObs1 - overallSkmObs2)
      }
      else {
        for (o in 1:nsubpop)
	  obsDifferences[o] <- (ObsCI1[o]-ObsCI2[o])-(overallObsCI1-overallObsCI2)
      }
      obsDifferences <- t(obsDifferences)
      stdObsDifferences <- apply(obsDifferences, 1, function(x) x/sigma)
      tObs <- apply(abs(stdObsDifferences), 2, max)
      pvalue <- sum(ifelse(tPerm > tObs, 1, 0))/noPerms
      signif(pvalue, 6)
#
#     generate covariance matrix, p-value based on hazard ratios
#
      sigma <- var(logHRs)

#     catch the error if we encounter a singularity problem when inverting the matrix sigma
#
      tryCatch({
      	sigmainv <- solve(sigma)
	      }, error = function(ex){
		cat(geterrmessage())
		cat("\n")
		set.seed<-4593432
		stop("STEPP has encountered an unexpected problem. Please retry. If problem persists, try different values of minpatspop(r1) and patspop(r2).")
	  	}
        )
      sigma2 <-sigma

      chi2perm <- rep(NA,noPerms)
      for (i in 1:noPerms){
        res <- matrix(logHRs[i,],ncol=1)
        chi2perm[i] <- t(res)%*%sigmainv%*%res
      }
      obsLogHR <- matrix(rep(NA, nsubpop),ncol=1)
      for (o in 1:nsubpop) obsLogHR[o] <- logHR[o] - overallLogHR
        chi2obs <- c(t(obsLogHR)%*%sigmainv%*%obsLogHR)
        HRchi2pvalue <- sum(ifelse(chi2perm > chi2obs, 1, 0))/noPerms
#
      sigma <- sqrt(diag(var(logHRs)))
      stdLogHRs <- t(apply(logHRs, 1, function(x) x/sigma))
      tPerm <- apply(abs(stdLogHRs), 1, max)
      obsLogHR <- matrix(rep(0, nsubpop))
      for (o in 1:nsubpop) obsLogHR[o] <- logHR[o]-overallLogHR
      obsLogHR <- t(obsLogHR)
      stdObsLogHR <- apply(obsLogHR, 1, function(x) x/sigma)
      tObs <- apply(abs(stdObsLogHR), 2, max)
      HRpvalue <- sum(ifelse(tPerm > tObs, 1, 0))/noPerms
#
#     prepare stepp object to be return
#
      res.obj <- list(sigma      = sigma1,
                    HRsigma      = sigma2,
                    pvalue       = pvalue,
                    HRpvalue     = HRpvalue,
                    chi2pvalue   = chi2pvalue,
                    HRchi2pvalue = HRchi2pvalue)
    }
    else
    {
	# create a res.obj if no permutation is done
	res.obj <- list(sigma      = NA,
                    HRsigma      = NA,
                    pvalue       = NA,
                    HRpvalue     = NA,
                    chi2pvalue   = NA,
                    HRchi2pvalue = NA)
    }
#
    stepp.obj <- list(version = version.stepp(),
                      cdate   = stime,
                      etime   = proc.time()-t0,
                      arg     = input.obj,
                      sp      = sp.obj,
                      est     = est.obj,
                      res     = res.obj)

    attr(stepp.obj, "class") <- "stepp"	
    stepp.obj
}




######################################################################
#
# Old interface - keep for compatibility reason
#
#######################################################################
analyze.KM.stepp<-function(coltrt, coltime, colcens, colvar, trts, patspop, minpatspop, 
    timest, noperm, legendy = 30, pline = -2.5, color = c("red", "black"),
    ylabel = "Specify Timepoint & Endpoint", xlabel = "Subpopulations by Median Covariate", 
    ncex = 0.7, tlegend = c("Specify 1st Treatment", "Specify 2nd Treatment"), 
    nlas = 0, alpha = 0.05, pointwise = FALSE)
{
   x<-stepp("KM", coltrt, coltime, colcens, , colvar, trts, patspop, minpatspop, timest, noperm)
   stepp_summary(x)
   stepp_print(x)
   stepp_plot(x, legendy, pline, color, ylabel, xlabel, ncex, tlegend, nlas, alpha, pointwise)
   x
}


analyze.CumInc.stepp<-function(coltrt, coltime, coltype, colvar, trts, patspop, minpatspop, 
    timest, noperm, legendy = 30, pline = -2.5, color = c("red", "black"),
    ylabel = "Specify Timepoint & Endpoint", xlabel = "Subpopulations by Median Covariate", 
    ncex = 0.7, tlegend = c("Specify 1st Treatment", "Specify 2nd Treatment"), 
    nlas = 0, alpha = 0.05, pointwise = FALSE)
{
   x<-stepp("CI", coltrt, coltime, , coltype, colvar, trts, patspop, minpatspop, timest, noperm)
   stepp_summary(x)
   stepp_print(x)
   stepp_plot(x, legendy, pline, color, ylabel, xlabel, ncex, tlegend, nlas, alpha, pointwise)
   x
}


