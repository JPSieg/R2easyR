#'Makes a list of accension numbers that are the same
#'
#'This function cross references the accenssion numbers in two data
#'frames to determine which accenssion numbers are in both data frames. It
#'requires that both data frames have a column titled "Accension", with the
#'Acension numbers that you want to cross reference.
#'
#'@param data1 The first data frame you want to cross reference
#'@param data2 The second data frame you want to cross reference
#'@return A data frame containing rows in data1 that contain accension numbers that are also in data2
#' @export
check.accn = function(data1, data2){
  both <- c()
  for (i in c(1:length(data1$Accension))){
    for (j in c(1:length(data2$Accension))){
      print(paste(i, j, sep = "/"))
      print(toString(data1$Accension[i]))
      print(toString(data2$Accension[j]))
      if (toString(data1$Accension[i]) == toString(data2$Accension[j])){
        print("True")
        both[i] <- i
      }
      else(print("False"))
    }
  }
  output <- data1$Accension[both]
}
