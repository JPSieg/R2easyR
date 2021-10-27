#'Read a .connect file file
#'
#'This function reads a .ct file, from the CoFold RNA folding algorithm, into a R data frame.
#'
#'@param data_file Path to the data file
#'@param T.to.U Default = TRUE then read connect will convert all "T"s to "U"s. If false read.connect will leave "T"s.
#'@return A data frame containing the information in the .ct file
#' @export
read.connect = function(data_file, T.to.U = TRUE){
  d = read.delim(data_file, sep = "\t", skip = 1, header = FALSE)
  colnames(d) <- c("N", "Nucleotide", "N-1", "N+1", "BP", "N")
  print(head(d))
  if (T.to.U){
    d$Nucleotide[which(d$Nucleotide == "T")] = "U"
  }
  output <- d
}
