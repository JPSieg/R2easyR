#'This function edits a R2R stockholm file to recolor nucleotides in a secondary structure.
#'
#'@param R2R.sto Path to the R2R stockholm file that contains the drawing information you want to optomize
#'@param Nucleotides A vector listing the nucleotide numbers you want to color
#'@param Color The color you want the nucleotides to be. Use R color formats. Default is Dim grey
#' @export
r2easyR.grey_letters_editor = function(R2R.sto,
                               Nucleotides,
                               Color = "dimgrey"){

  con <- file(R2R.sto)
  lines <- readLines(con)

  lines <- lines[1:(length(lines)-1)]

  new.lines <- paste("#=GF R2R nuc_color #:", Nucleotides - 1, " rgb:", col2rgb(Color)[1,1], ",", col2rgb(Color)[2,1], ",",col2rgb(Color)[3,1], sep ="")

  lines <- c(lines, new.lines, "//")

  print(lines)

  writeLines(lines, con)
  close(con)
}
