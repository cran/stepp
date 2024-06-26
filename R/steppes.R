#################################################################
#
# steppes.R
#
#############################
# stepp estimate/statistics #
#############################
setClass("steppes",
  representation(
    subpop   = "stsubpop",  # stepp subpopulation
    model    = "stmodel", # stepp model
    effect   = "ANY",   # list of absolute effect est
    testresults   = "ANY",   # test statistics
    nperm    = "numeric"  # number of permutations 0-n
  ),
  prototype = c(NULL, NULL, NULL, NULL, NULL, NULL)
)

setMethod("initialize", "steppes",
  function(.Object, subpop, model, effect, testresults, nperm, ...) {
    if (missing(subpop)) {
      subpop <- new("stsubpop")
    }
    if (missing(model)) {
      model <- new("stmodelKM")
    }
    if (missing(effect)) {
      effect <- NULL
    }
    if (missing(testresults)) {
      testresults <- NULL
    }
    if (missing(nperm)) {
      nperm <- 0
    }
    .Object@subpop <- subpop
    .Object@model <- model
    .Object@effect <- effect
    .Object@testresults <- testresults
    .Object@nperm <- nperm
    if (!validObject(.Object)) stop("")
    callNextMethod(.Object, ...)
  }
)

setValidity("steppes", 
  function(object) {
    status <- TRUE

    # add conditions to verify

    return(status)
  }
)

setMethod("estimate",
  signature = "steppes",
  definition = function(.Object, sp, model) {
    .Object@subpop <- sp
    .Object@model  <- model
    .Object@effect <- estimate(model, sp)
    return(.Object)
  }
)

setMethod("test",
  signature = "steppes",
  definition = function(.Object, nperm=100, showstatus=TRUE) {
    if (is.null(.Object@subpop) | is.null(.Object@model)) {
      print("You need to estimate the effects first before testing for interaction.")
    }
    else {
      if (length(.Object@model@trts) == 1) nperm <- 0
      .Object@nperm  <- nperm
      .Object@testresults <- test(.Object@model, nperm, .Object@subpop, .Object@effect, showstatus = showstatus)
    }
    return(.Object)
  }
)

setGeneric("group", function(.Object, criteria)
  standardGeneric("group")) 

setMethod("group",
  signature = "steppes",
  definition = function(.Object, criteria) {
    j <- which(criteria$trtid == .Object@model@trts)
    if (criteria$scale == "A") {
      # absolute scale
      TrtEff <- .Object@effect$TrtEff
      diff <- TrtEff[[j]]$sObs-TrtEff[[1]]$sObs
      if (criteria$threshold >=0) {
        mv <- diff >= criteria$threshold
      }
      else {
        mv <- diff <= criteria$threshold
      }
    } else if (criteria$scale == "R") {
      # relative scale
      Ratios <- .Object@effect$Ratios
      lratio <- Ratios[[j - 1]]$logHR
      lth <- log(criteria$threshold)
      if (lratio >=0){
        mv <- (lratio >= lth)
      }
      else {
        mv <- (lratio <= lth)
      }
    } else {
      stop("Unknown scale for merge criteria.")
    }

    return(list(mv=mv, left.side=.Object@subpop@minc, right.side=.Object@subpop@maxc))
  }
)

setGeneric("reestimate", function(.Object, sp)
  standardGeneric("reestimate"))  

setMethod("reestimate",
  signature = "steppes",
  definition = function(.Object, sp) {
    obj.new <- new("steppes")
    obj.new <- estimate(obj.new, sp, .Object@model)
    return(obj.new)
  }
)

setGeneric("cutpoint", function(.Object, criteria, side, debug)
  standardGeneric("cutpoint"))  

setMethod("cutpoint",
  signature = "steppes",
  definition = function(.Object, criteria, side, debug = 0) {
    grp <- group(.Object, criteria)

    mv <- grp$mv
    cv <- which(mv)
    if (debug){
      print(mv)
      print(length(cv))
      print(cv)
    }
    if (sum(mv) == 0){
      if (side == "L") cp <- round(.Object@subpop@maxc[length(mv)], digits = 4)
      else cp <- round(.Object@subpop@minc[1], digits = 4)
    } 
    else {
      if (side == "L") cutpoint.index <- which(mv)[1]
      else cutpoint.index <- cv[length(cv)]

      if (debug) print(.Object@subpop@medianz)
      cp <- .Object@subpop@medianz[cutpoint.index]
    }   
  
    return(cp)
  }
)

setGeneric("find.edge", function(.Object, criteria, debug)
  standardGeneric("find.edge")) 

setMethod("find.edge",
  signature = "steppes",
  definition = function(.Object, criteria, debug = 0) {
    grp   <- group(.Object, criteria)
    #print(paste("ogrp=",grp))
    if (sum(grp$mv) == 0) {
      rv <- NULL  # no edge if no subpopulation satisfies the criteria
    } else {
      s.new <- stepp.merge.subpop(.Object@subpop, grp$mv)
      #print(paste("grp$mv=", grp$mv))
      est1  <- reestimate(.Object, s.new)

      grp.new  <- group(est1, criteria)
      if (debug) print(paste("ngrp=",grp.new))
     
      if (sum(grp.new$mv) == 0){
        left  <- c(.Object@subpop@minc[1])
        right <- c(.Object@subpop@maxc[.Object@subpop@nsubpop])
        if (debug){
            print(paste("left.1=", left))
            print(paste("right.1=", right))
        }
      } else {
        grp.list <- which(grp.new$mv == 1)
        #print(paste("grp.list=", grp.list))
        ngrp <- length(grp.list)
        left <- rep(NA, ngrp)
        right <- rep(NA, ngrp)
        for (i in 1:ngrp){
          # edge analysis
          # Left edge
          if (debug) print("Do Left")
          if (grp.list[i] == 1){
            left[i] <- est1@subpop@minc[1]
            if (debug) print(paste("left.2i=", left[i]))
          } else {
            left.edge  <- stepp.edge.subpop(s.new, grp.list[i], "L")
            if (debug <= -2) summary(left.edge)
            est.left   <- stepp.test(left.edge, .Object@model, nperm=0)
            if (debug <= -35) print(est.left,estimate=TRUE, cov=FALSE, test=FALSE)
            if (debug == -30) plot(est.left, ncex=0.70,legendy=0, pv=FALSE,
              color=c("red", "black"), xlabel="Subpopulations by Median Value",
              ylabel="Outcome", tlegend=c("trt 0", "trt 1"), nlas=3, 
              pointwise=TRUE, noyscale=TRUE, rug=FALSE)

            left[i]    <- cutpoint(est.left,  criteria, "L", debug=debug)
            if (debug) print(paste("left.2i=", left[i]))
          }

          # Right edge
          if (debug) print("Do Right")
          if (grp.list[i] == length(grp.new$mv)){
            right[i] <- est1@subpop@maxc[grp.list[i]]
            if (debug) print(paste("right.2i=", right[i]))
          } else {
            right.edge <- stepp.edge.subpop(s.new, grp.list[i], "R")
            if (debug <= -2) summary(right.edge)
            est.right  <- stepp.test(right.edge, .Object@model, nperm=0)
            if (debug <= -35) print(est.right,estimate=TRUE, cov=FALSE, test=FALSE)
            if (debug <= -40) plot(est.right, ncex=0.70,legendy=0, pv=FALSE,
              color=c("red", "black"), xlabel="Subpopulations by Median Value",
              ylabel="Outcome", tlegend=c("trt 0", "trt 1"), nlas=3, 
              pointwise=TRUE, noyscale=TRUE, rug=FALSE)

            right[i] <- cutpoint(est.right, criteria, "R", debug=debug)
            if (debug) print(paste("right.3i=", right[i]))
          }
        }
      }
      rv = list(left = left, right = right)
    }
    return(rv)
  }
)

setGeneric("edge.boot", function(.Object, criteria, edge, j, seed, nsample, showstatus, debug)
  standardGeneric("edge.boot")) 

setMethod("edge.boot",
  signature = "steppes",
  definition = function(.Object, criteria, edge, j=2, seed=17, nsample=1000, showstatus=showstatus, debug=0) {
    ori.subpop <- .Object@subpop
    npatients  <- dim(ori.subpop@subpop)[1]
    bmatrix    <- matrix(0, nrow=npatients, ncol=nsample)

    set.seed (seed)
    # generate bootstrap samples 
    for (i in 1:nsample){
        bmatrix[,i] <- sample(seq(1,npatients), replace=TRUE)
    }

    nsub        <- length(edge$left)
    left.edge.boot  <- matrix(NA, nrow=nsample, ncol=nsub)
    right.edge.boot <- matrix(NA, nrow=nsample, ncol=nsub)

    ncount <- 0   
    # do this for each bootstrap sample
    if (showstatus) {
      title <- paste("\nEdge analysis with ", nsample)
      title <- paste(title, "bootstrap samples\n")
      cat(title)
      pb <- txtProgressBar(min=0, max=nsample-1, style=3)
    }

    for (i in 1:nsample) {
      if (showstatus) setTxtProgressBar(pb, i)

      # create the subpopulation
      subpop.i <- stepp.subpop(ori.subpop@win, ori.subpop@covar[bmatrix[,i]])

      # call to create a new model
      model.i  <- subgroup(.Object@model, bmatrix[,i])

      # estimate
      est.i    <- stepp.test(subpop.i, model.i, nperm=0)

      # find edge
      edge.i <- find.edge(est.i, criteria, debug = (debug < 0 | debug == i))
      if (is.null(edge.i)) {
        # if no edge is found; just set the edge to be the population edge
        # i.e. the largest variability allowed in the situation.
        edge.i <- list(left = c(.Object@subpop@minc[1]),
          right = c(.Object@subpop@maxc[.Object@subpop@nsubpop]))
      }
      if (debug < 0 | debug == i) {
        print(paste("edge.i.left",edge.i$left))
        print(paste("edge.i.right",edge.i$right))
      }

      # map ref edges to the boot sample edges
      for (k in 1:nsub) {
        if (length(edge.i$left) == 0) {
          left.edge.boot[i,k] <- min(cov)
          right.edge.boot[i,k] <- max(cov)
        } else {
          l <- edge$left[k]
          mdl <- min(abs(edge.i$left-l))
          lindex  <- which(mdl == abs(edge.i$left-l))[1]
          left.edge.boot[i,k] <- edge.i$left[lindex]
          r <- edge$right[k]
          mdr <- min(abs(edge.i$right-r))
          rindex <- which(mdr == abs(edge.i$right-r))[1]
          right.edge.boot[i,k]<- edge.i$right[rindex]
        }
        if (debug < 0 | debug == i) print(paste("final=",i,k,left.edge.boot[i,k],right.edge.boot[i,k]))
      }
    }
    if (showstatus) close(pb)

    # generate 95% confidence interval and other bootstrap statistics
    left.ci.hi  <- rep(0, nsub)
    left.ci.lo  <- rep(0, nsub)
    right.ci.hi <- rep(0, nsub)
    right.ci.lo <- rep(0, nsub)
    for (k in 1:nsub){
      left.temp <- sort(left.edge.boot[,k])
      left.ci.hi[k] <- left.temp[round(nsample*0.975)]
      left.ci.lo[k] <- left.temp[round(nsample*0.025)]
      right.temp <- sort(right.edge.boot[,k])
      right.ci.hi[k] <- right.temp[round(nsample*0.975)]
      right.ci.lo[k] <- right.temp[round(nsample*0.025)]
    }
  
    return(list(left = edge$left, right = edge$right, left.ci.hi = left.ci.hi, left.ci.lo = left.ci.lo, 
      right.ci.hi = right.ci.hi, right.ci.lo = right.ci.lo,
      bsample.left.edge = left.edge.boot, bsample.right.edge = right.edge.boot))
  }
)

setMethod("summary",
  signature = "steppes",
  definition = function(object) {
    summary(object@subpop)

    # print number of patients in each subpopulation for each treatment
    if (object@subpop@init) {
      subpop   <- object@subpop
      nsubpop  <- subpop@nsubpop
      covar   <- subpop@covar
      minc     <- subpop@minc
      maxc     <- subpop@maxc
    
      model    <- object@model
      trts     <- model@trts
      coltrt   <- model@coltrt
      ntrts    <- length(trts)
      txassign <- rep(-1, length(model@coltrt))
   
      for (j in 1:ntrts) txassign[which(coltrt == trts[j])] <- j
    
      cat("\n")
      write("Treatments Sample Size Information (with only specified treatments)", file="")
      trttitle <- paste(paste("trt",trts),"  ", collapse="")
      write(paste("     Subpopulation    ", trttitle, "Total"), file="")
      for (i in 1:nsubpop) {
        trtj <- rep(0, ntrts)
        for (j in 1:ntrts) {
          subpopj <- covar[covar >= minc[i] & covar <= maxc[i] & txassign == j]
          trtj[j] <- trtj[j] + length(subpopj)
        }
        subtrt <- paste(format(trtj,width=8),sep="",collapse="")
        if (subpop@win@type=="tail-oriented" & i==(length(subpop@win@r1)+1)){
          write(paste(format(i,width=12), "     ", subtrt, "     ", sum(trtj), "Entire Cohort"),file="")
        } else {
          write(paste(format(i,width=12), "     ", subtrt, "     ", sum(trtj)),file="")
        }
      }
      cat("\n")
      write("(To display subpopulation distribution together with the first STEPP Plot, use the 'subplot' option in the plot function.)",file="")
    }
    # cat("\n")
  }
)

setMethod("print",
  signature = "steppes",
  definition = function(x, estimate=TRUE, cov=TRUE, test=TRUE, ...) {
  ntrts <- x@effect$ntrts
  n <- 0
  # cat("\n")
  for (j in 1:ntrts) {
    nj <- sum(x@model@coltrt==x@model@trts[j])
    n  <- n + nj
    if (ntrts > 1) {
      write(paste0("Sample size in treatment ", x@model@trts[j], ": ", nj), file="")
    }
  }
  write(paste0("Total sample size (excluding missing data): ", n), file="")
  cat("\n")
  print(x@model, x, estimate, cov, test, ...)

  # if (test & ntrts > 1)
  #   write("Note: The p-values are not adjusted for multiple testing.",file="")
  }
)

#
# Internal worker routine for plot
.Stepp.plot <- function(x, y, legendy, pline, at, color, ylabel, xlabel, ncex, tlegend,
  nlas, alpha, pointwise, ci, pv, showss, ylimit, which, noyscale, rug,
  lsty, marker, subset, subplot, legend_diff, ...) {

  Result        <- x@testresults$Res
  nperm         <- x@nperm
  TrtEff        <- x@effect$TrtEff
  Ratios        <- x@effect$Ratios
  nsubpop       <- x@subpop@nsubpop
  npatsub       <- x@subpop@npatsub
  medians       <- x@subpop@medianz
  minc          <- x@subpop@minc
  maxc          <- x@subpop@maxc
  covar         <- x@subpop@covar
  r1            <- x@subpop@win@r1
  r2            <- x@subpop@win@r2
  coltrt        <- x@model@coltrt
  trts          <- x@model@trts
  issinglegroup <- (x@effect$ntrts == 1)
  
  # checks
  if (!is.numeric(which) || any(which < 1) || any(which > 3)) {
    stop("'which' must be in 1:3")
  }

  # set up the defaults
  ntrts <- length(trts)
  txassign <- rep(NA, length(coltrt))
  for (j in 1:ntrts) txassign[which(coltrt == trts[j])] <- j

  if (is.null(subset)) subset <- rep(1, ntrts)
  if (length(nperm) == 0) {
    pv <- FALSE # no pvalue displayed if permutation test not done
  } else if (nperm == 0) {
    pv <- FALSE # no pvalue displayed if permutation test not done
  }

  which <- sort(which)
  show <- rep(FALSE, 3)
  show[which] <- TRUE

  # Apply graphics defaults for lsty, marker
  # Check to make sure it is consistent
  if (is.null(lsty)) lsty = c(1,2,3,4,5,6,1,2)[1:ntrts]
  if (is.null(marker)) marker = seq(1:25)[1:ntrts]
  if (is.null(color)) color = c("black", "red", "blue", "orange", "green", "brown", "yellow", "magenta")[1:ntrts]

  if (length(color) > 8 | length(lsty) >  8 | length(marker) > 8 |
      length(color) != ntrts | length(lsty) != ntrts | length(marker) != ntrts) {
    print("Graphics options: color, lsty and marker must be NULL or equal to the number of treatments and no more than 8.")
    stop()
  }
  if (sum(subset) < 1) {
    print("Subset must have at least one treatment specified.")
    stop()
  }
  n <- length(which)
  if (is(x@model, "stmodelGLM")) noyscale <- TRUE
  if (subset[1] == 1 & sum(subset) == 1) n <- 1

  ssize <- rep(" ", nsubpop)
  tpatsub <- rep(0, nsubpop)

  for (i in 1:nsubpop) {
    trtj <- rep(0, ntrts)
    for (j in 1:ntrts) {
      subpopj <- covar[covar >= minc[i] & covar <= maxc[i] & txassign == j]
      trtj[j] <- trtj[j] + sum(!is.na(subpopj))
    }
    tpatsub[i] <- sum(trtj)
  }
  for (i in 1:nsubpop) ssize[i] <- paste(c("(n=", tpatsub[i], ")"), collapse = "")

  # parold <- par(no.readonly = TRUE)

  if (show[1L]) {
    # first plot (estimates)
    ncurves <- sum(subset)
    sell <- which(subset == 1)

    skmObs  <- NULL
    for (j in 1:ntrts) {
      if (subset[j] == 1) skmObs <- c(skmObs, TrtEff[[j]]$sObs)
    }
    xvalues <- rep(medians, ncurves)
    if (!noyscale) skmObs <- skmObs * 100
    group <- rep(1:ncurves, each = nsubpop)
    lbls <- rep(" ", nsubpop)
    
    if (!issinglegroup) {
      p <- "supremum pv = "
      for (j in 2:ntrts) {
        if (subset[j] == 1) {
          p <- paste(c(p, Result[[j - 1]]$pvalue), collapse = "")
        }
      }
    } else {
      p <- ""
    }

    # generate subpopulation distribution subplot if needed
    if (subplot) {
      par(fig = c(0, 1, 0, 0.8), omi = c(0.4, 0.4, 0.1, 0.1))
      # par(fig = c(0, 1, 0.3, 1), omi = c(0.4, 0.4, 0.1, 0.1))
    } else {
      par(mfrow = c(1, 1), omi = c(0.4, 0.4, 0.1, 0.1))
    }

    if (noyscale) {
      yl <- c(min(skmObs), max(skmObs))
      # legendy <- max(skmObs) - legendy
    } else {
      if (length(ylimit) < 2) {
        yl <- c(0, 100)
      } else {
        yl <- ylimit[1:2]
      }
    }

    if (noyscale) {
      plot(xvalues, skmObs, axes = TRUE, xaxt="n", ylim = yl, ylab = ylabel, xlab = "", type = "n")
    } else {
      plot(xvalues, skmObs, axes = FALSE, ylim = yl, ylab = ylabel, xlab = "", type = "n")
    }

    for (j in 1:ncurves) {
      points(xvalues[group == j], skmObs[group == j], lty = lsty[sell][j], lwd = 2,
        pch = marker[sell][j], type = "o", col = color[sell][j], bg = color[sell][j])
    }

    axis(1, at = xvalues, font = 1)

    if (!noyscale) axis(2, at = c(0, (0:10) * 10), font = 1)

    if (nlas != 3 & nlas != 2) {
      if (showss) {
        mtext(ssize, side = 1, at = xvalues, line = 2, font = 1, cex = ncex, adj = 0.5, las = nlas)
        mtext(xlabel, side = 1, line = 3.5)
      } else {
        mtext(xlabel, side = 1, line = 2)
      }
    } else {
      if (showss) {
        mtext(ssize, side = 1, at = xvalues, line = 3.5, font = 1, cex = ncex, adj = 0.5, las = nlas)
        mtext(xlabel, side = 1, line = 0.8, outer = TRUE)
      } else {
        mtext(xlabel, side = 1, line = 2)
      }
    }

    if (!issinglegroup) {
      legend(min(xvalues), legendy, pch = marker[sell], lty = lsty[sell], lwd = 2, 
        col = color[sell], pt.bg = color[sell], legend = tlegend[sell], bty = "n")
    }

    if (pv) {
      if (is.na(at)) {
        mtext(p, side = 1, at = (min(xvalues) + 0.2 * (max(xvalues) - min(xvalues))), line = pline)
      } else {
        mtext(p, side = 1, at = at, line = pline)
      }
    }
    if (rug) rug(xvalues)

    # WKY 
    # generate subpopulation distribution subplot
    if (subplot) {
      par(fig = c(0, 1, 0.58, 1), new = TRUE)
      # par(fig = c(0, 1, 0, 0.45), new = TRUE)
      n.y <- length(xvalues[group == 1])
      y <- seq(1, n.y)
      x.min <- rep(0, n.y)
      x.max <- rep(0, n.y)
      for (j in 1:n.y){
        x.min[j] <- minc[j]
        x.max[j] <- maxc[j]
      }
      # plot(xvalues[group == 1], y, axes = 0, cex = 0.5, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      plot(xvalues[group == 1], y, cex = 0.5, xaxt = "n", yaxt = "n", xlab = "", ylab = "")

      for (j in 1:n.y) {
        lines(x = c(x.min[j], x.max[j]), y = c(y[j], y[j]))
      }
    }
  }

  if (!issinglegroup) {
    # pointwise is specified, generate two additional plots
    subset[1] <- 1  # make sure the baseline trt is included
    ncurves <- sum(subset)
    sell <- which(subset == 1)     

    ndcurves <- ncurves - 1
    if (pointwise) {
      zcrit <- qnorm(1 - alpha/2)
    } else {
      zcrit <- qnorm(1 - alpha/(2*nsubpop))
    }
    xvalues <- rep(medians, ndcurves)
    group <- rep(1:ndcurves, each = nsubpop)

    if (show[2L]) {
      # second plot (treatment effect estimates)
      if (!rstudioapi::isAvailable()) {
        if (show[1L]) {
          dev.new()
        }
      }

      if (ndcurves != 0) {
        skmObs <- NULL
        se <- NULL
        for (j in 2:ntrts) {
          if (subset[j] == 1) {
            # skmObs  <- c(skmObs, TrtEff[[j]]$sObs - TrtEff[[1]]$sObs)
            skmObs <- c(skmObs, TrtEff[[1]]$sObs - TrtEff[[j]]$sObs)
            se <- c(se, sqrt(TrtEff[[j]]$sSE^2 + TrtEff[[1]]$sSE^2))
          }
        }
        if (!noyscale) skmObs <- skmObs * 100

        lbls <- rep(" ", nsubpop)

        # ssize <- rep(" ", nsubpop)
        # for (i in 1:nsubpop) ssize[i] <- paste(c("(n=", tpatsub[i], ")"), collapse = "")

        # already done, no need to repeat for pvalue
        par(mfrow = c(1, 1), omi = c(0.4, 0.4, 0.1, 0.1))

        if (noyscale) {
          ext <- (max(skmObs) - min(skmObs))*0.5
          yl <- c(min(skmObs) - ext, max(skmObs) + ext)
        } else {
          if (length(ylimit) < 4) {
            yl <- c(-100, 100)
          } else {
            yl <- ylimit[3:4]
          }
        }
      
        if (noyscale) {
          plot(xvalues, skmObs, axes = TRUE, xaxt="n", ylim = yl,
            ylab = paste("Difference in",ylabel), xlab = "", type = "n")
        } else {
          plot(xvalues, skmObs, axes = FALSE, ylim = yl,
            ylab = paste("Difference in",ylabel), xlab = "", type = "n")
        }

        for (j in 1:ndcurves) {
          points(xvalues[group == j], skmObs[group == j], lty = lsty[sell][j + 1],
            lwd = 2, pch = marker[sell][j + 1], type = "o", col = color[sell][j + 1],
            bg = color[sell][j + 1])
        }
        if (rug) rug(xvalues)

        if (ci) {
          if (!noyscale) {
            cilim <- skmObs - zcrit*se*100
            cilim <- pmax(cilim,-100)
          } else {
            cilim <- skmObs - zcrit*se
          }

          for (j in 1:ndcurves) {
            lines(xvalues[group == j], cilim[group == j], lty = 2, lwd = 2,
              col = color[sell][j + 1], bg = color[sell][j + 1])
          }
          if (!noyscale) {
            cilim <- skmObs + zcrit*se*100
            cilim <- pmin(cilim,100)
          } else { 
            cilim <- skmObs + zcrit*se
          }

          for (j in 1:ndcurves) {
            lines(xvalues[group == j], cilim[group == j], lty = 2, lwd = 2,
              col = color[sell][j + 1], bg = color[sell][j + 1])
          }
        }

        abline(h = 0, lty = 1)

        axis(1, at = xvalues, font = 1)
        if (!noyscale) axis(2, at = c(0, (-10:10) * 10), font = 1)
        if (nlas != 3 & nlas != 2) {
          if (showss) {
            mtext(ssize, side = 1, at = xvalues, line = 2, font = 1, cex = ncex, adj = 0.5, las = nlas)
            mtext(xlabel, side = 1, line = 3.5)
          } else {
            mtext(xlabel, side = 1, line = 2)
          }
        } else {
          if (showss) {
            mtext(ssize, side = 1, at = xvalues, line = 3.5, font = 1, cex = ncex, adj = 0.5, las = nlas)
            mtext(xlabel, side = 1, line = 0.8, outer = TRUE)
          } else {
            mtext(xlabel, side = 1, line = 2)
          }
        }

        if (pv) {
          if (is.na(at)) {
            mtext(p, side = 1, at = (min(xvalues) + 0.2 * (max(xvalues) - min(xvalues))), line = pline)
          } else {
            mtext(p, side = 1, at = at, line = pline)
          }
        }

        if (!is.null(legend_diff) & length(legend_diff) == 2) {
          plc <- switch(legend_diff[1], "topleft", "topright", "bottomright", "bottomleft")
          legend(x = plc, legend = paste(tlegend, collapse = " vs. "), bty = "n")
        }
      }
    }

    if (show[3L]) {
      # third plot (ratio estimates)
      if (!rstudioapi::isAvailable()) {
        if (any(show[1:2])) {
          dev.new()
        }
      }

      if (!is(x@model, "stmodelGLM")) {
        # Rpvalue <- x@testresults$HRpvalue
        # logR <- x@effect$logHR
        # logRSE <- x@effect$logHRSE
        text1 <- "Supremum HR p-value = "
        if (is(x@model, "stmodelCI")) {
          text2 <- "Subdistribution Hazard Ratio"
        } else {
          text2  <- "Hazard Ratio"
        }
      } else {
        # Rpvalue <- x@testresults$logRpvalue
        # logR <- x@effect$logR
        # logRSE <- x@effect$logRSE
        if (x@model@glm == "gaussian") {
          text1 <- "Supremum Effect Ratio p-value = "
          text2 <- paste(ylabel, "Ratio")
        } else if (x@model@glm == "binomial") {
          text1 <- "Supremum Odds Ratio p-value = "
          text2 <- "Odds Ratio"
        } else if (x@model@glm == "poisson") {
          text1 <- "Supremum Risks Ratio p-value = "
          text2 <- "Risks Ratio"
        }
      }
      p <- "supremum pv = "
      for (j in 2:ntrts) {
        if (subset[j] == 1) {
          p <- paste(c(p, Result[[j - 1]]$HRpvalue), collapse = "")
        }
      }
      par(mfrow = c(1, 1), omi = c(0.5, 0.5, 0.5, 0.5))
      if (length(ylimit) < 6) {
        yl <- c(0,3)
      } else {
        yl <- ylimit[5:6]
      }

      logHR <- NULL
      logHRSE <- NULL
      for (j in 2:ntrts) {
        if (subset[j] == 1) {
          logHR <- c(logHR, Ratios[[j - 1]]$logHR)
          logHRSE <- c(logHRSE, Ratios[[j - 1]]$logHRSE)
        }
      }

      sObs <- exp(logHR)
      plot(xvalues, exp(logHR), axes = FALSE, ylim = yl, ylab = text2, xlab = "", type = "n")
      for (j in 1:ndcurves) {
        points(xvalues[group == j], sObs[group == j], lty = lsty[sell][j + 1],
          lwd = 2, pch = marker[sell][j + 1], type = "o", col = color[sell][j + 1],
          bg = color[sell][j + 1])
      }

      if (ci) {
        cilim <- logHR - zcrit*logHRSE
        cilim <- exp(cilim)
        for (j in 1:ndcurves) {
          lines(xvalues[group == j], cilim[group == j], lty = 2, lwd = 2,
            col = color[sell][j + 1], bg = color[sell][j + 1])
        }
        cilim <- logHR + zcrit*logHRSE
        cilim <- exp(cilim)
        for (j in 1:ndcurves) {
          lines(xvalues[group == j], cilim[group == j], lty = 2, lwd = 2,
            col = color[sell][j + 1], bg = color[sell][j + 1])
        }
      }

      abline(h = 1, lty = 1)
      axis(1, at = xvalues, font = 1)
      axis(2, at = c(0, (0:15)*.2), font = 1)
      if (nlas != 3 & nlas != 2) {
        if (showss) {
        mtext(ssize, side = 1, at = xvalues, line = 2, font = 1, cex = ncex, adj = 0.5, las = nlas)
        mtext(xlabel, side = 1, line = 3.5)
        } else {
          mtext(xlabel, side = 1, line = 2)
        }
      } else {
        if (showss) {
          mtext(ssize, side = 1, at = xvalues, line = 3.5, font = 1, cex = ncex, adj = 0.5, las = nlas)
          mtext(xlabel, side = 1, line = 0.8, outer = TRUE)
        } else {
          mtext(xlabel, side = 1, line = 2)
        }
      }

      if (pv) {
        if (is.na(at)) {
          mtext(p, side = 1, at = (min(xvalues) + 0.2 * (max(xvalues) - min(xvalues))), line = pline)
        } else {
          mtext(p, side = 1, at = at, line = pline)
        }
      }

      if (!is.null(legend_diff) & length(legend_diff) == 2) {
        plc <- switch(legend_diff[2], "topleft", "topright", "bottomright", "bottomleft")
        legend(x = plc, legend = paste(tlegend, collapse = " vs. "), bty = "n")
      }
    }
  }

  # par(parold)
} # end of worker function

# S3 method
plot.steppes <- function(x, y, legendy = 30, pline = -2.5, color = c("red", "black"),
  ylabel = "Specify Timepoint & Endpoint",
  xlabel = "Subpopulations by Median Covariate",
  ncex = 0.7, tlegend = c("Specify 1st Treatment", "Specify 2nd Treatment"),
  nlas = 0, alpha = 0.05, pointwise = FALSE, ci = TRUE,
  pv = TRUE, showss = TRUE, ylimit = c(0,100,-100,100,0,3), which = c(1, 2, 3),
  noyscale = FALSE, rug = FALSE, at = NA, subplot = FALSE, legend_diff = 1, ...) {
  return(.Stepp.plot(x, y, legendy, pline, at, color, ylabel, xlabel, ncex, tlegend,
    nlas, alpha, pointwise, ci, pv, showss, ylimit, which,
    noyscale, rug, subplot, legend_diff, ...))
}

setMethod("plot",
  signature = "steppes",
  definition = function(x, y, 
    # graphics parameters - optional
    legendy = 30,   #   the vertical location of the legend according to the units on the y-axis
    pline = -2.5,   #   the vertical location of the p-value, starting at 0 counting outwards
    color = NULL,
            #   a vector containing the line colors for trts
    ylabel = "Specify Timepoint & Endpoint",    # label for the y-axis
    xlabel = "Subpopulations by Median Covariate",  # label for the x-axis
    ncex = 0.7,     #   the size of the text for the sample size annotation
    tlegend = c("Specify 1st Treatment", "Specify 2nd Treatment"),
            #   a vector containing the treatment labels, 1st and 2nd trt, respectively
    nlas = 0,     #   the las paramter (0,1,2,3) - the orientation of the sample size annotation
            # plot options - optional
    alpha = 0.05,   #   sig. level
    pointwise = FALSE,  #   pointwise confidence intervals (pointwise=TRUE),
            #   or confidence bands (pointwise=FALSE, default) to be displayed
    ci = TRUE,      #   display the conf. interval or band
    pv = TRUE,      #   display the supremum pvalue
    showss = TRUE,    #   display sample size on the x-axis
    ylimit = c(0,100,-100,100,0,3),
            #   y limits for the 3 graphs
    which = c(1, 2, 3),     #   which plots to produce
    noyscale = FALSE,   #   do not scale y axis to %
    rug = FALSE,    #   put a rug plot for each
    at = NA,      #   centering position for pline
    lsty = NULL,    #   line style used in the plots (no. 1-6 as in R)
    marker = NULL,    #   marker style used in the plots (correspond to marker style in R)
    subset = NULL,    #   subsets of treatment effects to be plotted; default is all
    subplot = FALSE,    #   show subpopulation distribution plot together with the plot
    legend_diff = rep(1, 2),    #   placements of the legend in the 2nd and 3rd graph
                                #   (1 = "topleft", 2 = "topright", etc., NULL = no legend)
    ... ) {
  return(.Stepp.plot(x, y, legendy, pline, at, color, ylabel, xlabel, ncex, tlegend,
        nlas, alpha, pointwise, ci, pv, showss, ylimit, which,
        noyscale, rug, lsty, marker, subset, subplot, legend_diff, ...))
  }
)

# constructor and access functions for steppes

# 1. estimate treatment effects and perform permutation test
stepp.test <- function(subpop, model, nperm, showstatus = TRUE) {
  result <- new("steppes")
  result <- estimate(result, subpop, model)
  result <- test(result, nperm, showstatus = showstatus)
  return(result)
}

# 2. edge analysis for stepp
stepp.edge <- function(est, criteria, j = 2, boot = 0, seed = 17, showstatus = TRUE, debug = 0) {
  bsample <- NULL
  ci    <- NULL
  edge.ref <- find.edge(est, criteria, debug)
  result   <- edge.ref
  if ((!is.null(edge.ref)) & (boot > 0)) {
    result <- edge.boot(est, criteria, edge.ref, j, seed, boot, showstatus, debug)
  }

  return (result)
}

stepp.group <- function(est, criteria) {
  return(group(est, criteria))
}

stepp.reestimate <- function(est, sp) {
  return(reestimate(est, sp))
}

stepp.cutpoint <- function(est, criteria, side, debug = 0) {
  return(cutpoint(est, criteria, side, debug))
}
