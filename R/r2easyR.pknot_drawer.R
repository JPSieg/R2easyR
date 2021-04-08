#'Edits a R2R stockholm file for drawing psuedoknot labels on RNA secondary structures
#'
#'Useful for adding psuedoknot labels to RNA secondary structures drawn with R2R. Counts the labels already specified by the stockholm
#'and adjusts accordingly.
#'
#'@param R2R.sto Path to the R2R stockholm file that contains the drawing information you want to add a pknot to
#'@param pknot A vector containing pknot nucleotide numbers
#' @export
r2easyR.pknot_drawer = function(R2R.sto,
                               pknot){
  if (class(pknot) != "list"){
    pknot = list(pknot)
  }

  for (l in 1:length(pknot)){

    con <- file(R2R.sto)
    lines <- readLines(con)

    lines <- lines[1:(length(lines)-1)]

    #Check to see if the pknot region is removed from the dotbracket secondary structure

    dot.bracket.line <- strsplit(strsplit(lines[3], split = "\t")[[1]][2], split = "")[[1]]

    for (m in 1:length(pknot[[l]])){
      #print(m)
      if (dot.bracket.line[pknot[[l]][m]] != "."){
        print(paste("Warning: The residue at position", m, "already has a secondary structure specified, which could result in innaccurate secondary structure depictions. Make sure you have removed pknots from the dotbracket column of the data frame before writing the R2R stockholm files"))
      }
    }

    #index of labels to use in the R2R label line

    label.index <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "x", "t", "u", "v", "w", "x", "y", "z")

    label.line <- strsplit(strsplit(lines[4], split = "\t")[[1]][2], split = "")[[1]]

    #Filters out labels that are already present in the label line (in case another pknot is already specified)

    if (length(which(label.index == intersect(label.index, label.line))) != 0){
      label.index <- label.index[-which(label.index == intersect(label.index, label.line))]
    }

    #add labels to the label line while keeping track of what labels have been added and not using the same labels for non contiguous nucleotides

    firstN <- TRUE
    label.index.pos <- 1
    unique.labels <- c()
    gaps <- c()

    for (i in 1:length(pknot[[l]])){
      if (firstN){
        label.line[pknot[[l]][i]] <- label.index[label.index.pos]
        firstN <- FALSE
        unique.labels <- c(unique.labels, label.index[label.index.pos])
      }else{
        if (pknot[[l]][i]-1 == pknot[[l]][i -1]){
          label.line[pknot[[l]][i]] <- label.index[label.index.pos]
        }else{
          label.index.pos <- label.index.pos + 1
          label.line[pknot[[l]][i]] <- label.index[label.index.pos]
          unique.labels <- c(unique.labels, label.index[label.index.pos])
          gaps <- c(gaps, i)
        }
      }
    }

    #Add an edited label line

    lines[4] <- paste("#=GC R2R_LABEL\t", gsub(", ", "", toString(label.line)), sep = "")

    #Counts how many pknots are already in the label line

    n.pknot = (length(c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "x", "t", "u", "v", "w", "x", "y", "z")) - length(label.index))/2 + 1

    #Find middle of contiguous pknot regions for label. Use gaps, eg where there are breaks in contiguous pknot sequences

    pknot.label.index <- c()

    if (length(gaps) == 1){
      pknot.label.index <- c(pknot.label.index, ceiling(mean(pknot[[l]][1:(gaps-1)])), ceiling(mean(pknot[[l]][(gaps):length(pknot[[l]])])))
    }else{
      #Don't know how to do this yet. I think I need to use a for loop to recursively specify mean pknot regions.
      print("You are recieving this error message because the creator did not have the time to figure this out when he wrote this code. please email jpsieg225@gmail.com or jus841@psu.edu if you encounter this error message")
    }

    #Generate Peter's label commands

    pknot.label.line <- rep(".", length(label.line))

    new.lines <- c()

    for (i in 1:length(unique.labels)){
      new.lines <- c(new.lines, paste('#=GF R2R outline_nuc', unique.labels[i]))
      new.lines <- c(new.lines, paste('#=GF R2R outline_nuc', unique.labels[i]))
      new.lines <- c(new.lines, paste('#=GF R2R tick_label pk:', unique.labels[i], ' pk', n.pknot, sep = ""))
      pknot.label.line[pknot.label.index[i]] <- unique.labels[i]
    }

    new.lines <- c(paste('#=GC R2R_XLABEL_pk\t', gsub(", ", "", toString(pknot.label.line)), sep = ""), new.lines)

    lines <- c(lines, new.lines, "//")

    print(lines[4])
    print(new.lines)

    writeLines(lines, con)
    close(con)
  }

}
