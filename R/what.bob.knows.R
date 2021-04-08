#'Identifies what known RNA were identified by a motif seach using RNAbob
#'
#'@param known.Accenssion Vector containing known accenssion numbers
#'@param known.Sequence Vector containing known RNA sequences
#'@param hit.Accenssion Vector containing hit accension numbers
#'@param hit.Sequence Vector containing hit sequences
#'@param hit.Descriptor Vector containing hit descriptors 
#'@return A R data frame containing the number of times RNAbob found a known RNA, the descriptors that found each RNA, and the RNA sequences that RNAbob found
#' @export
What.bob.knows = function(known.Accenssion,
                         known.Sequence,
                         hit.Accenssion,
                         hit.Sequence,
                         hit.Descriptor){
  found <- c()
  n.hits <- c()
  descriptors <- c()
  
  for (i in 1:length(known.Accenssion)){
    found[i] <- "hits"
    descriptors[i] <- "descriptors"
    known.Sequence[i] <- gsub("U", "T", gsub("[.]", "", x = known.Sequence[i]))
    known <- strsplit(known.Sequence[i], "")[[1]]
    if (length(which(hit.Accenssion == known.Accenssion[i])) > 0){
      hit <- hit.Sequence[which(hit.Accenssion == known.Accenssion[i])]
      hit.desc <- hit.Descriptor[which(hit.Accenssion == known.Accenssion[i])]
      same.seq <- 0
      length.same <- 0
      for(j in 1:length(hit)){
        if (is.na(hit[j]) == FALSE){
          hit.split <- strsplit(hit[j], split = "")[[1]]
          for (k in 1:length(known)){
            if (known[k] == hit.split[1]){
              if (length(known[k:length(known)]) >= length(hit.split)){
                if(gsub(", ", "", toString(known[k:(k + length(hit.split) - 1)])) == hit[j]){
                  found[i] <- paste(found[i], hit[j], sep = ";")
                  descriptors[i] <- paste(descriptors[i], hit.desc[j], sep = ";")
                }
              }
            } 
          }
        }
      }
    }
    n.hits[i] <- length(strsplit(found[i], split = ";")[[1]]) - 1
  }
  output <- data.frame("Accension" = known.Accenssion,
                       "Sequence" = known.Sequence,
                       "n.hits" = n.hits,
                       "Descriptors" = descriptors,
                       "Hits" = found)
  print(head(output))
  output <- output
}
