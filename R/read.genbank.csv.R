#'Read a csv file from genbank
#'
#'This function reads a csv formated allignment from the GENBANK database into a R
#'data.frame. It also separates the Accension number from other data into
#'separate column in the data frame. Data in data.frames are stored as strings.
#'
#'@param csv_file Path to the csv formated data file in GENBANK
#'@return A data frame containing the information in the csv formated data file
#' @export
read.genbank.csv = function(csv_file){
  csvdata <- readr::read.csv(csv_file)
  chromosome <- c()
  accn <- c()
  for (i in c(1:length(csvdata$X.Organism.Name))){
    a <- c()
    b <- c()
    c <- c()
    if (csvdata$Replicons[i] == ""){
      chromosome[i] <- ""
      accn[i] <- toString(csvdata$WGS[i])
    }
    else{
      a <- strsplit(toString(csvdata$Replicons[i]), "/")[[1]]
      if (length(a) > 1){
        b <- strsplit(a, ":")[[1]]
        if (length(b) == 3){
          chromosome[i] <- b[1]
          accn[i] <- b[3]
        }
        if (length(b == 2)){
          c <- strsplit(b[2], "_")[[1]]
          chromosome[i] <- b[1]
          accn[i] <- a[2]
        }
      }
      else{
        b <- strsplit(a, ":")[[1]]
        if (length(b) == 3){
          chromosome[i] <- b[1]
          accn[i] <- b[3]
        }
        else{
          chromosome[i] <- b[1]
          accn[i] <- b[2]
        }
      }
    }
  }
  for (i in c(1:length(csvdata$X.Organism.Name))){
    if (length(strsplit(accn[i], ";")[[1]]) > 1){
      accn[i] <- strsplit(accn[i], ";")[[1]][1]
    }
  }
  output <- data.frame("Species" <- csvdata$X.Organism.Name,
                       "Group" <- csvdata$Organism.Groups,
                       "Chromosome" <- chromosome,
                       "Accension" <- accn)
  colnames(output) <- c("Species", "Group", "Chromosome", "Accension")
  output <- output
}
