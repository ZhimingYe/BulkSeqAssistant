% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\name{wilDEGs}
\alias{wilDEGs}
\title{A DEG analysis method by wilcox test}
\usage{
wilDEGs(expr_matrix, sample_df, pval, fdr, FC, adj.method = "BH")
}
\arguments{
\item{adj.method}{}
}
\description{
input1:data_matrix，rowname is gene，colname is sample; input2:sample_group(df),column must contains Sample_ID and Type
}
\author{
Zhiming Ye
}
