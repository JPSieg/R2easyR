#'A that converts a RFAM formatted sto file to a connectivity (CT) table file
#'
#'Parses a RFAM formatted Stockholm alignment file, pulls out the sequence you want, and maps the consensus secondary structure onto a
#'connectivity table (CT).
#'
#'@param Path.to.sto Path to the sto file that you want to reformat
#'@param RNA.name Name of the RNA that you want to pull out
#'@param Output.name Prefix you want to put on the output CT file. Default = "RNA". Can be a path to a CT file storage directory.
#'@param Remove.non.canonical Do you want to remove non canonical and non GU wobble pairs produced from the consensus allignment. Default = TRUE. Set to Remove.non.canonical = FALSE.
#'@param Keep.pknots Do you want to keep pknots? Default = TRUE. Set to Keep.pknots = FALSE to remove pknots.
#' @export
RFAM_sto_to_ct = function(Path.to.sto,
                          RNA.name,
                          Output.name = "RNA",
                          Remove.non.canonical = TRUE,
                          Keep.pknots = TRUE){
  ####1. Read sto file####

  con = file(Path.to.sto)
  Lines = readLines(con)
  close(con)

  Lines <- Lines[-length(Lines)]
  Lines <- Lines[-1]
  Lines <- Lines[-which(Lines == "")]

  ####2.) Parse sto file lines using a for loop####

  #Identify our RNA and identify #=GC SS_cons

  for (i in 1:length(Lines)){
    if (strsplit(Lines[i], split = " ")[[1]][1] == "#=GF"){}else{
      if (strsplit(Lines[i], split = " ")[[1]][1] == "#=GC"){
        if (strsplit(Lines[i], split = " ")[[1]][2] == "SS_cons"){
          ss_cons_line = Lines[i]
        }else{}
      }else{
        if (strsplit(Lines[i], split = " ")[[1]][1] == RNA.name){RNA.seq.line = Lines[i]}
      }
    }
  }

  #Make a vector from the seq and ss_cons line

  ss_cons.vector <- strsplit(strsplit(ss_cons_line, split = " ")[[1]][-c(1, 2, which(strsplit(ss_cons_line, split = " ")[[1]] == ""))], split = "")[[1]]
  RNA.seq.vector <- strsplit(strsplit(RNA.seq.line, split = " ")[[1]][-c(1, which(strsplit(RNA.seq.line, split = " ")[[1]] == ""))], split = "")[[1]]

  #Remove characters that are missing in the seq vector and from both the seq and ss_cons line

  remove.index.vector <- which(RNA.seq.vector == "-")
  ss_cons.vector <- ss_cons.vector[-remove.index.vector]
  RNA.seq.vector <- RNA.seq.vector[-remove.index.vector]

  #ss_cons.vector
  #RNA.seq.vector

  ####3.) Replace all characters that mean single stranded with "." in the ss_cons vector####

  if (length(which(ss_cons.vector == ",")) != 0){ss_cons.vector[which(ss_cons.vector == ",")] <- "."}
  if (length(which(ss_cons.vector == "-")) != 0){ss_cons.vector[which(ss_cons.vector == "-")] <- "."}
  if (length(which(ss_cons.vector == "_")) != 0){ss_cons.vector[which(ss_cons.vector == "_")] <- "."}
  if (length(which(ss_cons.vector == ";")) != 0){ss_cons.vector[which(ss_cons.vector == ";")] <- "."}

  ss_cons.vector

  ####4.) Make an index of characters to pull out into their own vectors####

  pknot.character.index <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z")
  ss.character.index <- c("<", "(", "[", "{")
  ss.character.rev.index <- c(">", ")", "]", "}")

  ####5.) Split the GC_sscons line into a list of vectors. Have only one type of character in each vector. Exclude pknots here.####

  characters <- intersect(ss.character.index, ss_cons.vector)

  ss.cons.list <- {}

  for (i in 1:length(characters)){
    forward <- ss.character.index[which(ss.character.index == characters[i])]
    reverse <- ss.character.rev.index[which(ss.character.index == characters[i])]
    ss.vector <- rep(".", length(ss_cons.vector))
    ss.vector[which(ss_cons.vector == forward)] <- "<"
    ss.vector[which(ss_cons.vector == reverse)] <- ">"
    if (i == 1){
      ss.cons.list[[1]] <- ss.vector
    }else{
      ss.cons.list[[length(ss.cons.list) + 1]] <- ss.vector
    }
  }

  ####6.) Split the GC_sscons line into a list of vectors. Have only one type of character in each vector. Include pknots here.####

  if(Keep.pknots){
    characters <- intersect(pknot.character.index, ss_cons.vector)

    for (i in 1:length(characters)){
      forward <- pknot.character.index[which(pknot.character.index == characters[i])]
      reverse <-  tolower(forward)
      ss.vector <- rep(".", length(ss_cons.vector))
      ss.vector[which(ss_cons.vector == forward)] <- "<"
      ss.vector[which(ss_cons.vector == reverse)] <- ">"
      ss.cons.list[[length(ss.cons.list) + 1]] <- ss.vector
    }
  }

  #ss.cons.list

  ####7.) Recursively determine BP####

  list.BP <- {}

  for (i in 1:length(ss.cons.list)){
    dot.bracket <- ss.cons.list[[i]]
    BP <- rep(0, length(ss.cons.list[[i]]))
    while(length(which(dot.bracket != ".")) != 0){
      kill = TRUE
      if (length(which(dot.bracket == "<")) != 0 & length(which(dot.bracket == ">")) != 0){
        kill = FALSE
      }
      if (kill){
        dot.bracket[which(dot.bracket != ".")] <- "."
      }else{
        for (j in 1:length(dot.bracket)){
          if (dot.bracket[j] == "<"){
            if (dot.bracket[j + 1] == "."){
              if (dot.bracket[j + min(which(dot.bracket[(j+1):length(dot.bracket)] != "."))] == ">"){
                BP[j] <- j + min(which(dot.bracket[(j+1):length(dot.bracket)] != "."))
                BP[j + min(which(dot.bracket[(j+1):length(dot.bracket)] != "."))] <- j
                dot.bracket[j] <- "."
                dot.bracket[j + min(which(dot.bracket[(j+1):length(dot.bracket)] != "."))] <- "."
              }
            }
          }
        }
      }
      #print(dot.bracket)
      #print(BP)
    }
    list.BP[[i]] <- BP
  }

  ####8.) Make a vector for each easy ct column (N, Nucleotide, N+1, and N-1)####

  N = 1:length(RNA.seq.vector)
  Nucleotide = RNA.seq.vector
  N.minus.1 = N - 1
  N.plus.1 = N + 1
  BP = rep(0, length(N))

  ####9.) Consolidate the list of BP vectors into a single vector for the BP vector in the CT file####

  for (i in 1:length(list.BP)){
    replace.index <- which(list.BP[[i]] != 0)
    BP[replace.index] <- list.BP[[i]][replace.index]
  }

  ####10.) Make a data frame that contains the vectors for the ct file####

  output = data.frame(N, Nucleotide, N.minus.1, N.plus.1, BP, N)

  ####11.) Filter out non cannonical BPs that are not GC wobbles####

  keep.vector = c("AU", "UA", "GC", "CG", "GU", "UG")

  if (Remove.non.canonical){
    for (i in 1:length(output$N)){
      BP.ID = paste(output$Nucleotide[i], output$Nucleotide[output$BP[i]], sep = "")
      #print(BP.ID)
      if (length(which(keep.vector == BP.ID)) == 0){
        output$BP[i] <- 0
      }
    }
  }

  print(head(output))

  ####12.) Use write.delim to write a CT file with space delimeters compatible with RNAstructure####

  ?write

  ?write.table

  write.table(output, file = gsub("/", "_", paste(Output.name, "_", RNA.name, ".ct", sep = "")), sep = "    ")

  ####Make output####

  output = output
}
