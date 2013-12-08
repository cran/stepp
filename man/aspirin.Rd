\name{aspirin}
\alias{aspirin}
\docType{data}
\title{
  the aspirin data set.
}

\description{
  The "aspirin" data set is the data from the Aspirin/Folate Polyp Prevention Study Group collaborative investigating the effect of
  aspirin on the outcome of any occurrence of a colorectal adenoma and any occurrence of an advanced leson.  There
  are 1,121 subjects randomized to three doses of aspirin (placebo, 81 mg/day and 325 mg/day).  The subjects were
  followed for three years and then underwent colonoscopy.  Follow-up information is available for 1,084 subjects.
  The remainder did not undergo a follow-up exam.

  The data set contained here is a "simulated" data set based on the actual data from the Aspirin/Folate Polyp Prevention Study.  
}

\usage{data(aspirin)}
\format{
  The data set has the following columns:
    ID   - Patient ID number
    DOSE - Randomized aspirin dose (0=placebo)
    AGE  - Age of the patient
    HT   - Height (in meters)
    WT   - Weight (in Kg)
    G    - Gender (M or F)
    AD   - Any adenoma (0=no, 1=yes, .=no exam)
    LE   - Any advanced lesion (0=no, 1=yes, .=noexam)
}

\source{
  The Aspirin/Folate Polyp Prevention Study Group collaborative is acknowledged for permission to use the data from 
  the Aspirin/Folate Polyp Prevention study.
}

\references{
  Baron J et al (2003), "A Randomized Trial of Aspirin to prevent Colorectal Adenomas." n engl j med 2003;234:891-9.
}

\examples{
data(aspirin)

# remove cases with missing data
aspirinc <- aspirin[complete.cases(aspirin),]

# model matrix can only handle numeric data
# convert to 0, 1 encoding for gender
GENDER    <- as.numeric( aspirinc[,"G"] == "M")
aspirinc  <- cbind(aspirinc, GENDER)
remove("GENDER")

# make a subset of patients with placebo and 81 mg
attach(aspirinc)
subset1  <- DOSE == 0 | DOSE == 81
aspirin1 <- aspirinc[subset1,]
detach(aspirinc)

# set up treatment assignment
trtA     <- rep(0, dim(aspirin1)[1])
trtA[aspirin1[,"DOSE"] == 81] <- 1

# STEPP analysis A: placebo vs 81 mg aspirin
attach(aspirin1)
inc_win     <- stepp.win(type="sliding", r1=30, r2=100)
inc_sp      <- stepp.subpop(swin=inc_win, cov=AGE)

ADorLE      <- as.numeric(AD | LE)
modelA      <- stepp.GLM(coltrt=trtA, trts=c(0,1), colY=ADorLE, glm="binomial", debug=0)
# Warning: In this example, the permutations have been set to 50 to allow the function
# to finish in a short amount of time.  IT IS RECOMMEND TO USE AT LEAST 2500 PERMUTATIONS TO 
# PROVIDE STABLE RESULTS.
steppGLMA   <- stepp.test(inc_sp, modelA, nperm=50)
summary(steppGLMA)
print(steppGLMA)
plot(steppGLMA, ncex=0.70, legendy=30,
	pline=-4.5, color=c("red","black"),
	xlabel="Subpopulations by Median Age", ylabel="Risk",
	tlegend=c("Placebo", "81 mg aspirin"), nlas=3, pointwise=FALSE, noyscale=TRUE, rug=FALSE)

detach(aspirin1)
}

\keyword{datasets}

