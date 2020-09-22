#'Custom R2easyR script that writes R2R inputs based on Ryota's tRNA structure seq data
#'
#'This function reads a .ct file, reformats it to include dot-bracket secondary structure info,
#'adds shape reactivities from a .shape file, maps them to a R2easyR color pallet, and
#'and prints .sto and a R2R meta files, and rewrites the .sto files to include label line
#'and place explicit information to get tRNA that are drawn with the traditional layout.
#'Run R2R on the master.r2r_meta file to generate all tRNA at once.
#'
#'@param ct_file Path to the directory containing .ct files
#'@param shape_file Path to the directory containing .shape files
#'@param output_file Path to the directory you want the R2R .sto and metafiles to print in
#' @export
ry.tRNA = function(ct_file = "CT files",
                   shape_file = "shape files",
                   output_file = "R2R_files"){
  ####Load in CT files####

  df <- lapply(paste(ct_file, list.files(ct_file), sep = "/"), R2easyR::read.ct)

  ####Add dot bracket####

  df <- lapply(df, R2easyR::add.dot.bracket)

  ####Add reactvivity data####

  for (i in 1:length(df)){
    df[[i]]$Reactivity <- R2easyR::read.shape(paste(shape_file, list.files(shape_file)[i], sep = "/"))
  }

  df

  ####Generate color palettes####

  palettes <- R2easyR::r2easyR.palettes()

  ####Make a vector of file names####

  file.names <- list.files(ct_file)

  for (i in 1:length(file.names)){
    file.names[i] <- strsplit(toString(file.names[i]), split = "-maxexpect")[[1]][1]
  }

  file.names

  ####Map reactivity to a palette and print the file####

  for (i in 1:length(df)){
    df[[i]]$Reactivity[which(df[[i]]$Reactivity <= 0)] <- 0
  }

  df.color <- {}

  for (i in 1:length(df)){
    print(i)
    df.color[[i]] <- R2easyR::r2easyR.color(df[[i]],
                                            palette = palettes$YlOrRd.c,
                                            abs_reactivity_threshold = 0.2,
                                            manual.scale = c(0, 5))
    R2easyR::r2easyR(paste(output_file, file.names[i], sep = "/"), df.color[[i]], RNA_name = file.names[i], colors = "circles")
  }

  ####Insert drawing information into the R2R stockholm file####

  for (i in 1:length(df.color)){
    path <- paste(output_file, "/", file.names[i], ".sto", sep = "")
    conn <- file(path, open = "r")
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
    fileConn <- file(paste(output_file, "/", file.names[i], ".sto", sep = ""))
    writeLines(lines,
               fileConn)
    close(fileConn)
  }

  ####Write the master R2R meta file####

  lines <- paste(output_file, "/" , file.names, ".sto", "\t", "oneseq", "\t", file.names, sep = "")

  fileConn <- file("master.r2r_meta")
  writeLines(lines,
             fileConn)
  close(fileConn)

}
