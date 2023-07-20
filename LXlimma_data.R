
library(openxlsx)

gene_matrix <- read.xlsx("gene_matrix.xlsx")

usethis::use_data(gene_matrix,overwrite = T)

rm(list=ls())

data(gene_matrix)


