\name{RobSAOD-package}
\alias{RobSAOD-package}
\alias{RobSAOD}
\docType{package}
\title{
Robust Significance Analysis of Omics Data: A Weighted Approach for Biomarkers and Drug Discovery in Omics Data 
\packageTitle{RobSAOD}
}
\description{
This package implements a weight based significance analysis for identifying differentially expressed (DE) genes from noisy gene expression dataset (e.g.,microarray). This method is outlier-robust because a weight function is used to produce very smaller weight for outlying cell of the data matrix.
}

\details{
  \tabular{ll}{
        Package: \tab RobSAOD\cr
        Type: \tab Package\cr
        Version: \tab 1.0\cr
        Date: \tab 2022-08-06\cr
        License: \tab GPL \cr
        Depends: \tab R (>=2.10)\cr
    }
	

	Package RobSAOD has the following functions:
    \tabular{ll}{
         wmean(): \tab This function calculates the weighted mean \cr                                  \tab function.
\cr
         wVarSs(): \tab This function calculates the weighted sum of squares\cr                              \tab function.
\cr
         RobSAOD(): \tab This function calculates the weighted significance analysis of Omics Data  \cr

	   
    }
}
\author{
Nishith Kumar <nk.bru09@gmail.com>  and 
Md Ashad Alam <malam@tulane.edu>
%%  ~~who you are~~
}


\keyword{ package }





