% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ry.tRNA.R
\name{ry.tRNA}
\alias{ry.tRNA}
\title{Custom R2easyR script that writes R2R inputs based on Ryota's tRNA structure seq data}
\usage{
ry.tRNA(
  ct_file = "CT files",
  shape_file = "shape files",
  output_file = "R2R_files"
)
}
\arguments{
\item{ct_file}{Path to the directory containing .ct files}

\item{shape_file}{Path to the directory containing .shape files}

\item{output_file}{Path to the directory you want the R2R .sto and metafiles to print in}
}
\description{
This function reads a .ct file, reformats it to include dot-bracket secondary structure info,
adds shape reactivities from a .shape file, maps them to a R2easyR color pallet, and
and prints .sto and a R2R meta files, and rewrites the .sto files to include label line
and place explicit information to get tRNA that are drawn with the traditional layout.
Run R2R on the master.r2r_meta file to generate all tRNA at once.
}
