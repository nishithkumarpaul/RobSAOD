\name{RobSAOD}
\alias{RobSAOD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weighted significance analysis for identifying differentially expressed(DE) genes.
%%  ~~function to do ... ~~
}
\description{
This is a function of weighted significance analysis for identifying differentially expressed (DE) genes.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RobSAOD(x, ngroup = 3, noGrpMem = c(50,40,30))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a data matrix of  gene expression (e.g.,microarray), whose row contains genes and column contains subjects
%%     ~~Describe \code{x} here~~
}
  \item{ngroup}{
number of group in the dataset
%%     ~~Describe \code{ngroup} here~~
}
  \item{noGrpMem}{
a vector that contains the number of subjects of each group.
%%     ~~Describe \code{noGrpMem} here~~
}
}

\value{
calculate the p-value
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Nishith Kumar <nk.bru09@gmail.com> and 
Md Ashad Alam <malam@tulane.edu>
%%  ~~who you are~~
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
data(dummyData)
pvalue<-RobSAOD(dummyData,ngroup = 3, noGrpMem = c(50,40,30))
}
