\name{wVarSs}
\alias{wVarSs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weighted sum of square calculation
%%  ~~function to do ... ~~
}
\description{
This function calculates the weighted sum of the square. The user must put the weight in this function, otherwise, it will calculate the classical sum of the square. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wVarSs(x, w = c(rep(1, length(x))), center = 0)
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
  \item{center}{
Put the center value.Default is 0. 
%%     ~~Describe \code{center} here~~
}
}
\value{
Computed the weighted sum of square.
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
wVarSs(x,w= c(rep(1,10)),center=3)
wVarSs(x1,w= c(c(rep(1,7),c(rep(0.15,3)))),center=3)

}


