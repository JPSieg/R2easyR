#'Read a Stockholm file
#'
#'This function reads a Stockholm formated allignment from the RFAM database into a R
#'data.frame. It trims the long header at the front of the Stockholm file and keeps the
#'data. It also separates the Accension number and the sequence location into
#'separate column in the data frame. Data in data.frames are stored as strings.
#'
#'@param data_file Path to the data file
#'@return A data frame containing the information in the Stockholm file
#' @export
read.stock = function(data_file){
  stodata <- read.table(data_file, sep = "\t", stringsAsFactors = FALSE)
  accnum <- c()
  res <- c()
  RNA <- c()
  for (i in c(1:(length(stodata$V1)-1))){
    a <- c()
    b <- c()
    a <- strsplit(stodata$V1[i], "/")[[1]]
    accnum[i] <- a[1]
    b <- strsplit(a[2], " ")[[1]]
    res[i] <- b[1]
    RNA[i] <- gsub("-", ".", b[length(b)])
  }
  output <- data.frame("Accension" = accnum,
                       "Residues" = res,
                       "RNAseq" = RNA)
}
