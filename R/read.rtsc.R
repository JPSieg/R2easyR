#'Read a .rtsc file
#'
#'This function reads a .rtsc file from StructureFold2 into R
#'
#'@param data_file Path to the data file
#'@return A list of vectors containing stop count data from the different RNA in the .rtsc file
#' @export
read.rtsc = function(data_file){
  con = file(data_file)
  
  lines = readLines(con)
  
  lines = lines[-which(lines == "")]
  
  names.vector = lines[seq(1, length(lines) - 1, by = 2)]
  stops.vector = lines[seq(2, length(lines), by = 2)]
  
  parse.stops = function(x){as.numeric(as.character(strsplit(x, split = "\t")[[1]]))}
  
  Reactivity.list = lapply(stops.vector, FUN = parse.stops)
  
  names(Reactivity.list) = names.vector
  
  Output <- Reactivity.list
}
