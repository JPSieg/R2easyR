#'A function that optomizes the stem layout for a RNA 2-structure drawn with R2R
#'
#'This function edits a R2R stockholm file to provide a better stem layout for most RNA. By default, R2R will place
#'adjacent stems on the same side of the line of the RNA backbone. This will cause the stems to clash in the figure
#'if the stems are to close together in the 1-sequence. The r2easyR.stem_editor prevents this by flipping every other
#'stem to the oposite side of the backbone using the R2R label line and the R2R place_explict command. The r2easyR.stem_editor
#'currently supports optimizing the layout of up a structure with up to 26 stems.
#'
#'@param R2R.sto Path to the R2R stockholm file that contains the drawing information you want to optomize
#'@param source Where you got the stockholm. Did you make it with "R2easyR" or did you make it with "R2R"? Used to parse the secondary structure line.
#' @export
r2easyR.stem_editor = function(R2R.sto,
                               source = "R2easyR"){

  con <- file(R2R.sto)
  lines <- readLines(con)

  if (source == "R2easyR"){
    SS_cons <- strsplit(strsplit(lines[3], split = "\t")[[1]][2], split = "")[[1]]
  }
  if (source == "R2R"){
    SS_cons_line <- lines[which(grepl("#=GC SS_cons", lines))]
    SS_cons_line = gsub("#=GC SS_cons", "", SS_cons_line)
    SS_cons_line = gsub(" ", "", SS_cons_line)
    SS_cons = strsplit(SS_cons_line, split = "")[[1]]
    SS_cons[which(SS_cons == "(")] = "<"
    SS_cons[which(SS_cons == ")")] = ">"
  }

  open_pairs <- 0
  stems <- 0
  R2R_LABEL <- 0
  ends_in_a_helix <- FALSE
  next.N.start.stem <- FALSE
  zoro_n_juction <- c()

  open_labels <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M")
  close_labels <- c("N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z")

  length(SS_cons)

  for (i in 1:length(SS_cons)){
    if(next.N.start.stem) {
      open_pairs = open_pairs + 1
      stems = stems + 1
      R2R_LABEL[i] = open_labels[(stems - 1)]
      next.N.start.stem <- FALSE
    } else {
      if (stems > 1){
        if (open_pairs != 0){
          if (SS_cons[i] == "<"){
            open_pairs = open_pairs + 1
            R2R_LABEL[i] = "."
          }
          if (SS_cons[i] == ">"){
            open_pairs = open_pairs - 1
            R2R_LABEL[i] = "."
            if (open_pairs == 0){
              if (stems %% 2 == 0){
                if (i != length(SS_cons)){
                  R2R_LABEL[i] = close_labels[(stems - 1)]
                }
                if (i != length(SS_cons)){
                  if (SS_cons[(i+1)] == "<"){
                    next.N.start.stem <- TRUE
                    zoro_n_juction[stems] <- TRUE
                  }else{
                    zoro_n_juction[stems] <- FALSE
                  }
                }
                if (i == length(SS_cons)){
                  ends_in_a_helix <- TRUE
                  next.N.start.stem <- FALSE
                }
              }else{
                if (i != length(SS_cons)){
                  if (SS_cons[(i+1)] == "<"){
                    next.N.start.stem <- TRUE
                    zoro_n_juction[stems] <- TRUE
                  }else{
                    zoro_n_juction[stems] <- FALSE
                  }
                }
                if (i == length(SS_cons)){
                  ends_in_a_helix <- TRUE
                  next.N.start.stem <- FALSE
                }
              }
            }
          }
          if (SS_cons[i] == "."){
            R2R_LABEL[i] = "."
          }
        }
        if (open_pairs == 0){
          if (SS_cons[i] == "<"){
            stems = stems + 1
            open_pairs = open_pairs + 1
            if (stems %% 2 == 0){
              R2R_LABEL[i] = open_labels[(stems - 1)]
            }else{
              R2R_LABEL[i] = "."
            }
          }
          if (SS_cons[i] == "."){
            R2R_LABEL[i] = "."
          }
        }
      }
      if (stems == 1){
        if (open_pairs != 0){
          if (SS_cons[i] == "<"){
            open_pairs = open_pairs + 1
            R2R_LABEL[i] = "."
          }
          if (SS_cons[i] == ">"){
            open_pairs = open_pairs - 1
            R2R_LABEL[i] = "."
            if (open_pairs == 0){
              if (i != length(SS_cons)){
                if (SS_cons[(i+1)] == "<"){
                  next.N.start.stem <- TRUE
                  zoro_n_juction[stems] <- TRUE
                }else{
                  zoro_n_juction[stems] <- FALSE
                }
              }
            }
          }
          if (SS_cons[i] == "."){
            R2R_LABEL[i] = "."
          }
        }
        if (open_pairs == 0){
          if (SS_cons[i] == "<"){
            stems = stems + 1
            open_pairs = open_pairs + 1
            R2R_LABEL[i] = open_labels[(stems - 1)]
          }
          if (SS_cons[i] == "."){
            R2R_LABEL[i] = "."
          }
        }
      }
      if (stems == 0){
        if (open_pairs == 0){
          if (SS_cons[i] == "<"){
            stems = stems + 1
            open_pairs = open_pairs + 1
            R2R_LABEL[i] = "."
          }
          if (SS_cons[i] != "<"){
            R2R_LABEL[i] = "."
          }
        }
      }
    }

    #print(paste("i =", i))
    #print(paste("open pairs", open_pairs))
    #print(paste("stems", stems))
    #print(next.N.start.stem)
    #print(R2R_LABEL)

  }

  lines[4] <- paste("#=GC R2R_LABEL", gsub(", ", "", toString(R2R_LABEL)), sep = "\t")

  if (SS_cons[length(SS_cons)] == ">"){ends_in_a_helix == TRUE}else{ends_in_a_helix == FALSE}

  if (length(zoro_n_juction) <= 1){zoro_n_juction <- c(zoro_n_juction, c(FALSE, FALSE, FALSE, FALSE))}

  if (length(SS_cons) != length(SS_cons[which(SS_cons == ".")])){
    if (ends_in_a_helix == FALSE){
      for ( i in 1:stems){
        if (i %% 2 == 0){
          if (zoro_n_juction[i-1]){
            lines <- c(lines[1:(length(lines)-1)],
                       paste("#=GF R2R place_explicit ", open_labels[(i - 1)], " ", open_labels[(i - 1)], "-- -45 1 0 0 0 0 f", sep = ""),
                       paste("#=GF R2R place_explicit ", close_labels[(i - 1)], "++ ", close_labels[(i - 1)], " 45 1 0 0 0 90 f", sep = ""),
                       "//")
          }else{
            lines <- c(lines[1:(length(lines)-1)],
                       paste("#=GF R2R place_explicit ", open_labels[(i - 1)], " ", open_labels[(i - 1)], "-- 45 1 0 0 0 90 f", sep = ""),
                       paste("#=GF R2R place_explicit ", close_labels[(i - 1)], "++ ", close_labels[(i - 1)], " 45 1 0 0 0 90 f", sep = ""),
                       "//")
          }
        }
      }
    }
    if (ends_in_a_helix == TRUE){
      for ( i in 1:(stems - 1)){
        if (i %% 2 == 0){
          if (zoro_n_juction[i-1]){
            lines <- c(lines[1:(length(lines)-1)],
                       paste("#=GF R2R place_explicit ", open_labels[(i - 1)], " ", open_labels[(i - 1)], "-- -45 1 0 0 0 0 f", sep = ""),
                       paste("#=GF R2R place_explicit ", close_labels[(i - 1)], "++ ", close_labels[(i - 1)], " 45 1 0 0 0 90 f", sep = ""),
                       "//")
          }else{
            lines <- c(lines[1:(length(lines)-1)],
                       paste("#=GF R2R place_explicit ", open_labels[(i - 1)], " ", open_labels[(i - 1)], "-- 45 1 0 0 0 90 f", sep = ""),
                       paste("#=GF R2R place_explicit ", close_labels[(i - 1)], "++ ", close_labels[(i - 1)], " 45 1 0 0 0 90 f", sep = ""),
                       "//")
          }
        }
      }
      if (zoro_n_juction[(stems-1)]){
        lines <- c(lines[1:(length(lines)-1)],
                   paste("#=GF R2R place_explicit ", open_labels[(stems - 1)], " ", open_labels[(stems - 1)], "-- -45 1 0 0 0 0 f", sep = ""),
                   "//")
      }else{
        lines <- c(lines[1:(length(lines)-1)],
                   paste("#=GF R2R place_explicit ", open_labels[(stems - 1)], " ", open_labels[(stems - 1)], "-- 45 1 0 0 0 90 f", sep = ""),
                   "//")
      }
    }
  }

  print(lines)

  writeLines(lines, con)
  close(con)
}
