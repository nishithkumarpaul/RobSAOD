\name{dummyData}
\alias{dummyData}
\docType{data}
\title{
Artificially generated  gene expression (e.g.,microarray) data
}
\description{
This dataset is a matrix containing measurements of gene expression of different subjects.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("dummyData")}
\format{
  A data matrix (e.g.,microarray) with 2000 rows and 120 columns. Each row is a gene. The columns are the measurements of gene expression from different subjects.
}
\details{
This data set contains three groups. The first group contains 50 subjects, the second group includes 40 subjects, and the third group contains 30 subjects. of 2000 genes, the first 300 are differentially expressed (DE) genes, and the rest are equally expressed (DE) genes. This data also contain outliers. 
%%  ~~ If necessary, more details than the __description__ above ~~
}
%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(dummyData)
## str(dummyData) 
heatmap(dummyData)
}
\keyword{datasets}
