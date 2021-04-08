#'A function that reformats Peter's double tabdelimited files for reading in to the R2easyR pipeline
#'
#'This function edits a double tab delimited ct (conectivity file) to make it compatible with the read.ct function in R2easyR.
#'
#'@param file.ct Path to the ct file that you want to reformat
#' @export
peter.ct.reformat = function(file.ct){

con = file("FSE.ct")

Lines = readLines(con)

Lines = Lines[-which(Lines == "")]

for (i in 1:length(Lines)){
  Lines[i] <- gsub("\t\t", "   ", Lines[i])
  #print(gsub("\t\t", "   ", Lines[i]))
}

Lines[1] <- "N   Nucleotide  N-1   N+1   BP   N"

Lines

writeLines(Lines, con)

close(con)
}
