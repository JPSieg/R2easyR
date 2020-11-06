#'Loads an RNABOB output file into a R data frame
#'
#'This function reads a RNABOBoutput file into a R data frame
#'
#'@param data_file Path to the RNABOB output file
#'@return A data frame containing the data from the RNABOB output file
#' @export
read.rnabob = function(data_file){
  con = file(data_file, "r")

  database <- c()
  descriptor <- c()
  seq.f <- c()
  seq.t <- c()
  accen.number <- c()
  species <- c()
  sequence <- c()
  n <- 1

  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) != 0 ) {
      if(2 == length(strsplit(line, split = ":")[[1]])){
        if(strsplit(line, split = ":")[[1]][1] == "Database file"){ database <-  gsub(" ", "", strsplit(line, split = ":")[[1]][2]) }
        if(strsplit(line, split = ":")[[1]][1] == "Descriptor file"){ descriptor <-  gsub(" ", "", strsplit(line, split = ":")[[1]][2]) }
      }
      if(line == " seq-f  seq-t     name     description"){
        while ( TRUE ){
          a <- readLines(con, n = 1)
          if (length(strsplit(a, split = " ")[[1]]) > 4){
            if (length(which(strsplit(a, split = " ")[[1]] == "")) != 0){
              b <- strsplit(a, split = " ")[[1]][-which(strsplit(a, split = " ")[[1]] == "")]
            }
            if (length(which(strsplit(a, split = " ")[[1]] == "")) == 0){
              b <- strsplit(a, split = " ")[[1]]
            }
            seq.f <- c(seq.f, b[1])
            seq.t <- c(seq.t, b[2])
            accen.number <- c(accen.number, b[3])
            species <- c(species, toString(b[c(4:5)]))
          }
          if (length(strsplit(a, split = " ")[[1]]) == 1){
            c <- strsplit(a, "")[[1]]
            c <- c[-which(c == "|")]
            sequence <- c(sequence, gsub(", ", "", toString(c)))
          }
          if ( a == "" ) {
            break
          }
        }
      }
    }
    if ( length(line) == 0 ) {
      break
    }
  }

  close(con)

  if (length(sequence) > 0){
    output <- data.frame(database,
                         descriptor,
                         seq.f,
                         seq.t,
                         accen.number,
                         species,
                         sequence)
  }
  if (length(sequence) == 0){
    output <- data.frame(database,
                         descriptor,
                         seq.f = NA,
                         seq.t = NA,
                         accen.number = NA,
                         species = NA,
                         sequence = NA)
  }
  print(head(output))
  output <- output
}
