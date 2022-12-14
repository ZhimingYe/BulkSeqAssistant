% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\name{PlotNES}
\alias{PlotNES}
\title{A useful plot with NES and pvalue to display the result of GSEA analysis}
\usage{
PlotNES(
  gseRes,
  UpperLim,
  LowerLim,
  Col1 = "#B39CD0",
  Col2 = "#FBEAFF",
  SepCol = F,
  HiExprGroupName = "A",
  LoExprGroupName = "B",
  WarpWidth = 50
)
}
\arguments{
\item{gseRes}{GSEA result generated by clusterProfiler.}

\item{UpperLim}{Display the 1 to <UpperLim> terms.}

\item{LowerLim}{Display the <LowerLim> to the final result. It is usually the total number minus the required number and then adding one.}

\item{Col1}{Color 1}

\item{Col2}{Color 2}

\item{SepCol}{Method of display result. Display according groups or P values. if True display by groups, false displayed by pval.default is false.}

\item{HiExprGroupName}{Char vector, when SepCol=T, define the hi-expr group name.}

\item{LoExprGroupName}{Char vector, when SepCol=T, define the low-expr group name.}

\item{WarpWidth}{Wrap length}
}
\value{
a ggplot2 object.
}
\description{
A useful plot with NES and pvalue to display the result of GSEA analysis
}
\author{
Zhiming Ye
}
