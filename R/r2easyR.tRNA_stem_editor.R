#'A function that optomizes the stem layout for a tRNA 2-structure drawn with R2R
#'
#'This function edits a R2R stockholm file to provide a better stem layout for tRNA.
#'
#'@param R2R.sto Path to the R2R stockholm file that contains the drawing information you want to optomize
#' @export
r2easyR.tRNA_stem_editor = function(R2R.sto){
  conn <- file(R2R.sto, open = "r")
  lines <- readLines(conn)
  close(conn)
  for (j in 1:length(lines)){ #adds a label "a" to properly format the 3' end
    a <- lines[j]
    a <- strsplit(a, split = "\t")[[1]]
    if (a[1] == "#=GC SS_cons"){
      label.line <- c()
      ss <- c()
      b <- strsplit(a[2], split = "")[[1]]
      for (k in length(b):1){
        if (b[k] != ">"){
          label.line <- c(".", label.line)
        }
        if (b[k] == ">"){
          ss <- c(ss, 1)
          if (length(ss) == 1){
            label.line <- c("a", label.line)
          }
          if (length(ss) != 1){
            label.line <- c(".", label.line)
          }
        }
      }
      d <- FALSE
      if (b[1] == "."){
        label.line[min(which(b == "<"))] <- "b"
        d <- TRUE
      }
    }
    if (a[1] == "#=GC R2R_LABEL"){
      lines[j] <- paste("#=GC R2R_LABEL", "\t", gsub(", ", "", toString(label.line)), sep = "")
    }
  }
  lines <- lines[-length(lines)]
  if (d == FALSE){
    lines <- c(lines, c("#=GF R2R tick_label_regular_numbering 0 20 firstNucNum 1",
                        "#=GF R2R set_dir pos0 90 f",
                        "#=GF R2R place_explicit a++ a 45 1 0 0 0 90 f",
                        "//"))
  }
  if (d == TRUE){
    lines <- c(lines, c("#=GF R2R tick_label_regular_numbering 0 20 firstNucNum 1",
                        "#=GF R2R place_explicit b b-- 45 1 0 0 0 90 f",
                        "#=GF R2R place_explicit a++ a 45 1 0 0 0 90 f",
                        "//"))
  }
  fileConn <- file(R2R.sto)
  writeLines(lines,
             fileConn)
  close(fileConn)
}
