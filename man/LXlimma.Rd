\name{LXlimma}
\alias{LXlimma}
\title{LXlimma::LXlimma}
\usage{
LXlimma(gene_matrix)
}
\description{
LXlimma can be used to analysis gene expression matrix.
}
\examples{

if(!requireNamespace("devtools"))
  install.packages("devtools")

library(devtools)

install_github("gluck4668/LXlimma")

library(LXlimma)

??LXlimma

#-------------------
data(gene_matrix)

#------------------

rm(list=ls())


gene_matrix = "gene_matrix.xlsx"


LXlimma(gene_matrix)

}
