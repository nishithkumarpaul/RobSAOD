\name{wmean}
\alias{wmean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of Weighted Mean
%%  ~~function to do ... ~~
}
\description{
This function calculates the weighted mean. The user must put the weight in this function; otherwise, it will calculate the classical mean. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wmean(x, w = c(rep(1, length(x))))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a numerical vector
%%     ~~Describe \code{x} here~~
}
  \item{w}{
value of the weight of the vector. The default weight is 1 for each observation of the vector.
%%     ~~Describe \code{w} here~~
}
}
\value{
Computed the weighted mean.
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
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
set.seed(1234)
x<-c(rnorm(10,3,0.25)) 
x1<-c(rnorm(7,3,0.25),rnorm(3,20,2)) #contain 20 outliers#
wmean(x,w= c(rep(1,10)))
mean(x)
wmean(x1,w= c(c(rep(1,7),c(rep(0,3)))))
mean(x1)
}
