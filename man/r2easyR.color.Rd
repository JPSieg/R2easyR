% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R2ReasyR_functions.R
\name{r2easyR.color}
\alias{r2easyR.color}
\title{Assigns colors to reactivity values}
\usage{
r2easyR.color(
  data_frame,
  palette,
  no_data = "dimgrey",
  abs_reactivity_threshold = 0,
  manual.scale = FALSE,
  write_legend = TRUE
)
}
\arguments{
\item{data_frame}{A data_frame containing reactivity values in a column labeled "Reactivity"}

\item{palette}{A palette contisting of a vector containing 35 R colors. Compatible palettes are easily generated with r2easyR.palettes or r2easyR.custom.palette.}

\item{no_data}{The color you want the nucleotide to have when there is no data. Default = "dimgrey"}

\item{abs_reactivity_threshold}{Minimum threshold for reactivity data to be mapped to a color. Default = 0. Using a higher threshold prevents numerous, but low reactivity values from cluttering up the finished picture.}

\item{manual.scale}{Manually specify a min and max. Useful if you want to color reactivity from multiple transcripts using 1 scale. Default is "FALSE". To set the scale, use manual.scale = c(min, max)}

\item{write_legend}{Option to not write the legend or the Rxn plot. Default TRUE. If set to FALSE it will skip time consuming graphics writing steps. Good for running in bulk.}
}
\value{
A data frame containing a R2R label and Color column
}
\description{
Assigns a color to each reactivity value in a dataframe. Also, saves a
PDF legend and a PDF depicting the reactivity values in your current working directory.
}
