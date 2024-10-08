#'Read a .ct file
#'
#'This function reads a .ct file, from RNA folding algorithms, into a R data frame.
#'
#'@param data_file Path to the data file
#'@return A data frame containing the information in the .ct file
#' @export
read.ct = function(data_file){
  path <- data_file
  conn <- file(path,open="r")
  lines <- readLines(conn)
  close(conn)
  a <- lines[-1]
  N <- c()
  Nucleotide <- c()
  N.minus.1 <- c()
  N.plus.1 <- c()
  BP <- c()
  N.1 <- c()
  for (i in 1:length(a)){
    b <- strsplit(toString(a[i]), " ")
    b <- b[[1]][-which(b[[1]] == "")]
    N <- c(N, b[1])
    Nucleotide <- c(Nucleotide, b[2])
    N.minus.1 <- c(N.minus.1, b[3])
    N.plus.1 <- c(N.plus.1, b[4])
    BP <- c(BP, b[5])
    N.1 <- c(N.1, b[6])
  }
  d <- data.frame("N" = N,
                  "Nucleotide" = Nucleotide,
                  "N-1" = N.minus.1,
                  "N+1" = N.plus.1,
                  "BP" = BP,
                  "N" = N.1)
  colnames(d) <- c("N", "Nucleotide", "N-1", "N+1", "BP", "N")
  print(head(d))
  output <- d
}

#'Converts RNA secondary structure information in a .ct file into dot-bracket notation
#'
#'This function reformats information in a dataframe, made from reading a .ct file into R,
#'into the dot-bracket notation.
#'
#'@param ctdata_frame A data frame made by reading a .ct file into R with read.ct
#'@return A data frame containing a collumn with dot-bracket formated RNA secondary structure
#' @export
add.dot.bracket = function(ctdata_frame){
  dotbracket <- c()
  for (i in c(1:length(ctdata_frame$N))){
   if (ctdata_frame$BP[i] == 0){
     dotbracket[i] <- "."
   }
    else{
      if (as.numeric(as.character(ctdata_frame$BP[i])) > as.numeric(as.character(ctdata_frame$N[i]))){ dotbracket[i] <- "<" }
      if (as.numeric(as.character(ctdata_frame$BP[i])) < as.numeric(as.character(ctdata_frame$N[i]))){ dotbracket[i] <- ">" }
    }
  }
  ctdata_frame$Dotbracket <- dotbracket
  output <- ctdata_frame
}

#'Read a .react file
#'
#'This function reads a .react file from StructureFold2 into R
#'
#'@param data_file Path to the data file
#'@return A list of vectors containing reactivity data from the different RNA in the .react file
#' @export
read.react = function(data_file){
  path <- data_file
  conn <- file(path,open="r")
  lines <- readLines(conn)
  close(conn)
  odds <- seq(1, length(lines), 2)
  evens <- seq(2, length(lines), 2)
  reactivity <- {}
  for (i in c(1:length(evens))){
    reactivity[[i]] <- strsplit(lines[evens[i]], "\t")[[1]]
  }
  names(reactivity) <- lines[odds]
  for (i in c(1:length(reactivity))){
    for(j in c(1:length(reactivity[[i]]))){
      if (toString(reactivity[[i]][j]) == "NA"){reactivity[[i]][j] <- NA}
      if (toString(reactivity[[i]][j]) != "NA"){reactivity[[i]][j] <- as.numeric(toString(reactivity[[i]][j]))}
    }
  }
  for (i in c(1:length(reactivity))){
    reactivity[[i]] <- as.numeric(as.character(reactivity[[i]]))
  }
  output <- reactivity
}

#'Read a .shape file
#'
#'This function reads a .react file from StructureFold2 into R
#'
#'@param data_file Path to the data file
#'@return A list of vectors containing reactivity data from the different RNA in the .react file
#' @export
read.shape = function(data_file){
  path <- data_file
  conn <- file(path,open="r")
  lines <- readLines(conn)
  close(conn)
  reactivity <- c()
  for (i in 1:length(lines)){
    a <- lines[i]
    reactivity <- c(reactivity, as.numeric(strsplit(toString(a), "\t")[[1]][2]))
  }
  reactivity[which(reactivity == -999)] <- NA
  print(reactivity)
  output <- reactivity
}

#'Generates a custom palette
#'
#'Makes a custom R2Reasy palette from a vector containing less than 35 colors
#'
#'@param colors A vector containing less than 35 R colors
#'@return A custom palette that can be used by r2easyR.color
#' @export
r2easyR.custom.palette = function(colors){
  a <- rep(colors[1], (floor(35/length(colors)) + (35 - (floor(35/length(colors))*length(colors)))))
  for (i in c(2:length(colors))){
    a <- c(a, rep(colors[i], floor(35/length(colors))))
  }
  palette <- a
}

#'Generates a list of palettes
#'
#'Makes a list of 59 color palettes for r2easyR.colors. Also saves a PDF depicting
#'all of the color palettes in your current working directory. Color palettes are generated using RColorBrewer
#'or Viridis. ".l" palettes are recommended for coloring letters because they exclude light shades, which would be
#'hard to see against a white background. ".c" palettes are recommended for coloring circles behind letters
#'because they exclude dark colors, which obscure the letter inside the circle. Viridis palettes, Viridis, Magma,
#'Plasma, Inferno, and cividis are not recommended because they result in cluttered secondary structures.
#'
#'@return A list of 59 vectors containing Viridis and Colorbrewer palettes
#' @export
r2easyR.palettes = function(){
  palettes_names <- c("Viridis",
                     "Magma",
                     "Plasma",
                     "Inferno",
                     "Cividis",
                     "YlOrRd",
                     "YlOrBr",
                     "YlGnBu",
                     "YlGn",
                     "Reds",
                     "RdPu",
                     "Purples",
                     "PuRd",
                     "PuBuGn",
                     "PuBu",
                     "OrRd",
                     "Oranges",
                     "Greys",
                     "Greens",
                     "GnBu",
                     "BuPu",
                     "BuGn",
                     "Blues",
                     "Spectral",
                     "RdYlGn",
                     "RdYlBu",
                     "RdGy",
                     "RdBu",
                     "PuOr",
                     "PRGn",
                     "PiYG",
                     "BrBG",
                     "YlOrRd",
                     "YlOrBr",
                     "YlGnBu",
                     "YlGn",
                     "Reds",
                     "RdPu",
                     "Purples",
                     "PuRd",
                     "PuBuGn",
                     "PuBu",
                     "OrRd",
                     "Oranges",
                     "Greys",
                     "Greens",
                     "GnBu",
                     "BuPu",
                     "BuGn",
                     "Blues",
                     "Spectral",
                     "RdYlGn",
                     "RdYlBu",
                     "RdGy",
                     "RdBu",
                     "PuOr",
                     "PRGn",
                     "PiYG",
                     "BrBG")
  palettes <- {}
  palettes[[1]] <- viridis::viridis(35, end = 0.9) #Viridis palettes
  palettes[[2]] <- viridis::magma(35, end = 0.85)
  palettes[[3]] <- viridis::plasma(35, end = 0.9)
  palettes[[4]] <- viridis::inferno(35, end = 0.9)
  palettes[[5]] <- viridis::cividis(35, end = 1)
  for (i in c(6:23)){ #Pos or neg palettes from color brewer palettes for letters
    palettes[[i]] <- c(rep(RColorBrewer::brewer.pal(9, name = palettes_names[i])[5], 7),
                      rep(RColorBrewer::brewer.pal(9, name = palettes_names[i])[6], 7),
                      rep(RColorBrewer::brewer.pal(9, name = palettes_names[i])[7], 7),
                      rep(RColorBrewer::brewer.pal(9, name = palettes_names[i])[8], 7),
                      rep(RColorBrewer::brewer.pal(9, name = palettes_names[i])[9], 7))
  }
  for (i in c(24:32)){#Neg to pos palettes from color brewer palettes for letters
    palettes[[i]] <- c(rep(RColorBrewer::brewer.pal(11, name = palettes_names[i])[11], 4),
                      rep(RColorBrewer::brewer.pal(11, name = palettes_names[i])[10], 4),
                      rep(RColorBrewer::brewer.pal(11, name = palettes_names[i])[9], 4),
                      rep(RColorBrewer::brewer.pal(11, name = palettes_names[i])[8], 5),
                      rep(RColorBrewer::brewer.pal(11, name = palettes_names[i])[4], 5),
                      rep(RColorBrewer::brewer.pal(11, name = palettes_names[i])[3], 5),
                      rep(RColorBrewer::brewer.pal(11, name = palettes_names[i])[2], 4),
                      rep(RColorBrewer::brewer.pal(11, name = palettes_names[i])[1], 4))
  }
  for (i in c(33:50)){ #Pos or neg palettes from color brewer palettes for circles
    palettes[[i]] <- c(rep(RColorBrewer::brewer.pal(9, name = palettes_names[i])[2], 5),
                      rep(RColorBrewer::brewer.pal(9, name = palettes_names[i])[3], 6),
                      rep(RColorBrewer::brewer.pal(9, name = palettes_names[i])[4], 6),
                      rep(RColorBrewer::brewer.pal(9, name = palettes_names[i])[5], 6),
                      rep(RColorBrewer::brewer.pal(9, name = palettes_names[i])[6], 6),
                      rep(RColorBrewer::brewer.pal(9, name = palettes_names[i])[7], 6))
  }
  for (i in c(51:59)){#Neg to pos palettes from color brewer palettes for circles
    palettes[[i]] <- c(rep(RColorBrewer::brewer.pal(11, name = palettes_names[i])[10], 5),
                      rep(RColorBrewer::brewer.pal(11, name = palettes_names[i])[9], 4),
                      rep(RColorBrewer::brewer.pal(11, name = palettes_names[i])[8], 4),
                      rep(RColorBrewer::brewer.pal(11, name = palettes_names[i])[7], 4),
                      rep(RColorBrewer::brewer.pal(11, name = palettes_names[i])[6], 4),
                      rep(RColorBrewer::brewer.pal(11, name = palettes_names[i])[4], 4),
                      rep(RColorBrewer::brewer.pal(11, name = palettes_names[i])[3], 5),
                      rep(RColorBrewer::brewer.pal(11, name = palettes_names[i])[2], 5))
  }
  names(palettes) <- c("Viridis",
                      "Magma",
                      "Plasma",
                      "Inferno",
                      "Cividis",
                      "YlOrRd.l",
                      "YlOrBr.l",
                      "YlGnBu.l",
                      "YlGn.l",
                      "Reds.l",
                      "RdPu.l",
                      "Purples.l",
                      "PuRd.l",
                      "PuBuGn.l",
                      "PuBu.l",
                      "OrRd.l",
                      "Oranges.l",
                      "Greys.l",
                      "Greens.l",
                      "GnBu.l",
                      "BuPu.l",
                      "BuGn.l",
                      "Blues.l",
                      "Spectral.l",
                      "RdYlGn.l",
                      "RdYlBu.l",
                      "RdGy.l",
                      "RdBu.l",
                      "PuOr.l",
                      "PRGn.l",
                      "PiYG.l",
                      "BrBG.l",
                      "YlOrRd.c",
                      "YlOrBr.c",
                      "YlGnBu.c",
                      "YlGn.c",
                      "Reds.c",
                      "RdPu.c",
                      "Purples.c",
                      "PuRd.c",
                      "PuBuGn.c",
                      "PuBu.c",
                      "OrRd.c",
                      "Oranges.c",
                      "Greys.c",
                      "Greens.c",
                      "GnBu.c",
                      "BuPu.c",
                      "BuGn.c",
                      "Blues.c",
                      "Spectral.c",
                      "RdYlGn.c",
                      "RdYlBu.c",
                      "RdGy.c",
                      "RdBu.c",
                      "PuOr.c",
                      "PRGn.c",
                      "PiYG.c",
                      "BrBG.c")
  yvalues <- c(59:1)
  df <- data.frame("xaxis" = c(1:35)/35,
                   "yaxis" = yvalues[1],
                   "Palettes" = names(palettes)[1],
                   "Colour" = palettes[[1]])
  for (i in c(2:length(palettes_names))){
    df <- rbind(df, data.frame("xaxis" = c(1:35)/35,
                               "yaxis" = yvalues[i],
                               "Palettes" = names(palettes)[i],
                               "Colour" = palettes[[i]]))
  }
  graph <- ggplot2::ggplot(data = df, ggplot2::aes(x = xaxis, y = yaxis)) +
    ggplot2::geom_point(show.legend = FALSE, size = 3, shape = 15, colour = df$Colour) +
    ggplot2::geom_text(data = df[df$xaxis == 1, ], ggplot2::aes(label = names(palettes)), nudge_x = 0.15, nudge_y = 0.08) +
    ggplot2::xlab("Scale") +
    ggplot2::ylab("") +
    ggplot2::scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1.2)) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.line = ggplot2::element_line(colour = 'black', size = 1.5),
          axis.ticks = ggplot2::element_line(colour = "black", size = 1.5),
          axis.text.x = ggplot2::element_text(color = "Black", size = 16),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.line.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_text(color = "Black", size = 18),
          axis.title.y = ggplot2::element_text(color = "Black", size = 18),
          legend.text = ggplot2::element_text(color = "Black", size = 16),
          legend.title = ggplot2::element_text(color = "Black", size = 16))
  ggplot2::ggsave(filename = "palettes.pdf", path = getwd(), plot = graph, scale = 2.5, width = 5, height = 12, units = "cm", dpi = 300)
  output <- palettes
}

#'Assigns colors to reactivity values
#'
#'Assigns a color to each reactivity value in a dataframe. Also, saves a
#'PDF legend and a PDF depicting the reactivity values in your current working directory.
#'
#'@param data_frame A data_frame containing reactivity values in a column labeled "Reactivity"
#'@param palette A palette contisting of a vector containing 35 R colors. Compatible palettes are easily generated with r2easyR.palettes or r2easyR.custom.palette.
#'@param abs_reactivity_threshold Minimum threshold for reactivity data to be mapped to a color. Default = 0. Using a higher threshold prevents numerous, but low reactivity values from cluttering up the finished picture.
#'@param no_data The color you want the nucleotide to have when there is no data. Default = "dimgrey"
#'@param manual.scale Manually specify a min and max. Useful if you want to color reactivity from multiple transcripts using 1 scale. Default is "FALSE". To set the scale, use manual.scale = c(min, max)
#'@param write_legend Option to not write the legend or the Rxn plot. Default TRUE. If set to FALSE it will skip time consuming graphics writing steps. Good for running in bulk.
#'@return A data frame containing a R2R label and Color column
#' @export
r2easyR.color = function(data_frame,
                         palette,
                         no_data = "dimgrey",
                         abs_reactivity_threshold = 0,
                         manual.scale = FALSE,
                         write_legend = TRUE){
  if (min(data_frame$Reactivity, na.rm = TRUE) >= 0){ #For reactivity data where all values are greater than 0
    a <- c()
    for (i in c(1:length(data_frame$Reactivity))){
      if (is.na(data_frame$Reactivity[i])){a[i] <- NA}
      else{
        if (abs(as.numeric(toString(data_frame$Reactivity[i]))) < abs_reactivity_threshold){a[i] <- NA}
        if (abs(as.numeric(toString(data_frame$Reactivity[i]))) >= abs_reactivity_threshold){a[i] <- as.numeric(toString(data_frame$Reactivity[i]))}
      }
    }
    if (length(manual.scale) == 1){
      values <- max(a, na.rm = TRUE)*(c(1:length(palette))/length(palette))
    }
    if (length(manual.scale) != 1){
      values <- manual.scale[2]*(c(1:length(palette))/length(palette))
    }
    b <- c()
    c <- c()
    for (i in c(1:length(a))){
      if(is.na(a[i])){
        b[i] <- no_data
        c[i] <- "0"
      }
      else{
        for (j in c(length(values):1)){
          if (a[i] <= values[j]){
            b[i] <- palette[j]
            c[i] <- "1"
          }
          if (a[i] > max(values)){
            b[i] <- palette[35]
            c[i] <- "1"
          }
        }
      }
    }
  df <- data.frame("xaxis" = c(1:length(a)),
                   "yaxis" = a,
                   "Colour" = b)
  df2 <- data.frame("xaxis" = c(1),
                    "yaxis" = values,
                    "Colour" = palette)
  levels(df2$Colour) <- c(levels(df2$Colour), c("white"))
  for (i in c(1:length(df2$yaxis))){
    if (df2$yaxis[i] < abs_reactivity_threshold){df2$Colour[i] <- "white"}
  }
  legend <- ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= 0,ymax= df2$yaxis[1]), fill= df2$Colour[1]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[1], ymax= df2$yaxis[2]), fill= df2$Colour[2]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[2], ymax= df2$yaxis[3]), fill= df2$Colour[3]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[3], ymax= df2$yaxis[4]), fill= df2$Colour[4]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[4], ymax= df2$yaxis[5]), fill= df2$Colour[5]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[5], ymax= df2$yaxis[6]), fill= df2$Colour[6]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[6], ymax= df2$yaxis[7]), fill= df2$Colour[7]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[7], ymax= df2$yaxis[8]), fill= df2$Colour[8]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[8], ymax= df2$yaxis[9]), fill= df2$Colour[9]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[9], ymax= df2$yaxis[10]), fill= df2$Colour[10]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[10], ymax= df2$yaxis[11]), fill= df2$Colour[11]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[11], ymax= df2$yaxis[12]), fill= df2$Colour[12]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[12], ymax= df2$yaxis[13]), fill= df2$Colour[13]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[13], ymax= df2$yaxis[14]), fill= df2$Colour[14]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[14], ymax= df2$yaxis[15]), fill= df2$Colour[15]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[15], ymax= df2$yaxis[16]), fill= df2$Colour[16]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[16], ymax= df2$yaxis[17]), fill= df2$Colour[17]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[17], ymax= df2$yaxis[18]), fill= df2$Colour[18]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[18], ymax= df2$yaxis[19]), fill= df2$Colour[19]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[19], ymax= df2$yaxis[20]), fill= df2$Colour[20]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[20], ymax= df2$yaxis[21]), fill= df2$Colour[21]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[21], ymax= df2$yaxis[22]), fill= df2$Colour[22]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[22], ymax= df2$yaxis[23]), fill= df2$Colour[23]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[23], ymax= df2$yaxis[24]), fill= df2$Colour[24]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[24], ymax= df2$yaxis[25]), fill= df2$Colour[25]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[25], ymax= df2$yaxis[26]), fill= df2$Colour[26]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[26], ymax= df2$yaxis[27]), fill= df2$Colour[27]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[27], ymax= df2$yaxis[28]), fill= df2$Colour[28]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[28], ymax= df2$yaxis[29]), fill= df2$Colour[29]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[29], ymax= df2$yaxis[30]), fill= df2$Colour[30]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[30], ymax= df2$yaxis[31]), fill= df2$Colour[31]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[31], ymax= df2$yaxis[32]), fill= df2$Colour[32]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[32], ymax= df2$yaxis[33]), fill= df2$Colour[33]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[33], ymax= df2$yaxis[34]), fill= df2$Colour[34]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[34], ymax= df2$yaxis[35]), fill= df2$Colour[35]) +
    ggplot2::xlab("") +
    ggplot2::xlim(0, 1) +
    ggplot2::ylab("Reactivity") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste("N = (no data)\nor\nReactivity < ", abs_reactivity_threshold)) +
    ggplot2::theme(axis.line = ggplot2::element_line(colour = 'black', size = 1.5),
                   axis.ticks = ggplot2::element_line(colour = "black", size = 1.5),
                   title = ggplot2::element_text(color = no_data, size =16),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(color = "Black", size = 16),
                   axis.title.x = ggplot2::element_text(color = "Black", size = 18),
                   axis.title.y = ggplot2::element_text(color = "Black", size = 18),
                   legend.text = ggplot2::element_text(color = "Black", size = 16),
                   legend.title = ggplot2::element_text(color = "Black", size = 16))
  if (write_legend){
    ggplot2::ggsave(filename = "Legend.pdf", path = getwd(), plot = legend, scale = 2.5, width = 2.85, height = 5, units = "cm", dpi = 300)
  }
  rxnplot <- ggplot2::ggplot(data = df, ggplot2::aes(x = xaxis, y = yaxis)) +
    ggplot2::geom_point(show.legend = FALSE, size =2, colour = df$Colour) +
    ggplot2::xlab("Residue") +
    ggplot2::ylab("Reactivity") +
    ggplot2::theme_classic()+
    ggplot2::theme(axis.line = ggplot2::element_line(colour = 'black', size = 1.5),
                   axis.ticks = ggplot2::element_line(colour = "black", size = 1.5),
                   axis.text.x = ggplot2::element_text(color = "Black", size = 16),
                   axis.text.y = ggplot2::element_text(color = "Black", size = 16),
                   axis.title.x = ggplot2::element_text(color = "Black", size = 18),
                   axis.title.y = ggplot2::element_text(color = "Black", size = 18),
                   legend.text = ggplot2::element_text(color = "Black", size = 16),
                   legend.title = ggplot2::element_text(color = "Black", size = 16))
  if (write_legend){
    ggplot2::ggsave(filename = "Rxn_plot.pdf", path = getwd(), plot = rxnplot, scale = 2.5, width = 7, height = 5, units = "cm", dpi = 300)
  }
  }
  if (min(data_frame$Reactivity, na.rm = TRUE) < 0){ #For reactivity data where some values are less than 0
    if (length(manual.scale) == 1){
      a <- c()
      for (i in c(1:length(data_frame$Reactivity))){
        if (is.na(data_frame$Reactivity[i])){a[i] <- NA}
        else{
          if (abs(as.numeric(toString(data_frame$Reactivity[i]))) < abs_reactivity_threshold){a[i] <- NA}
          if (abs(as.numeric(toString(data_frame$Reactivity[i]))) >= abs_reactivity_threshold){a[i] <- as.numeric(toString(data_frame$Reactivity[i]))}
        }
      }
      values1 <- -max(abs(a), na.rm = TRUE)*(c(17:1)/17)
      values2 <- 0
      values3 <- max(abs(a), na.rm = TRUE)*(c(1:17)/17)
      values <- c(values1, values2, values3)
    }
    if (length(manual.scale) != 1){
      a <- c()
      for (i in c(1:length(data_frame$Reactivity))){
        if (is.na(data_frame$Reactivity[i])){a[i] <- NA}
        else{
          if (abs(as.numeric(toString(data_frame$Reactivity[i]))) < abs_reactivity_threshold){a[i] <- NA}
          if (abs(as.numeric(toString(data_frame$Reactivity[i]))) >= abs_reactivity_threshold){a[i] <- as.numeric(toString(data_frame$Reactivity[i]))}
        }
      }
      values1 <- manual.scale[1]*(c(17:1)/17)
      values2 <- 0
      values3 <- manual.scale[2]*(c(1:17)/17)
      values <- c(values1, values2, values3)
    }
    b <- c()
    c <- c()
    for (i in c(1:length(a))){
      if(is.na(a[i])){
        b[i] <- no_data
        c[i] <- "0"
      }
      else{
        for (j in c((length(values)):1)){
          if (a[i] <= values[j]){
            b[i] <- palette[c(1:35)[j]]
            c[i] <- "1"
          }
          if (a[i] > max(values)){
            b[i] <- palette[35]
            c[i] <- "1"
          }
        }
      }
    }
    df <- data.frame("xaxis" = c(1:length(a)),
                     "yaxis" = a,
                     "Colour" = b)
    df2 <- data.frame("xaxis" = c(1),
                      "yaxis" = values,
                      "Colour" = palette)
    levels(df2$Colour) <- c(levels(df2$Colour), c("white"))
    for (i in c(1:length(df2$yaxis))){
      if (abs(df2$yaxis[i]) < abs_reactivity_threshold){df2$Colour[i] <- "white"}
    }
    legend <- ggplot2::ggplot() +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[1], ymax= df2$yaxis[2]), fill= df2$Colour[1]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[2], ymax= df2$yaxis[3]), fill= df2$Colour[2]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[3], ymax= df2$yaxis[4]), fill= df2$Colour[3]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[4], ymax= df2$yaxis[5]), fill= df2$Colour[4]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[5], ymax= df2$yaxis[6]), fill= df2$Colour[5]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[6], ymax= df2$yaxis[7]), fill= df2$Colour[6]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[7], ymax= df2$yaxis[8]), fill= df2$Colour[7]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[8], ymax= df2$yaxis[9]), fill= df2$Colour[8]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[9], ymax= df2$yaxis[10]), fill= df2$Colour[9]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[10], ymax= df2$yaxis[11]), fill= df2$Colour[10]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[11], ymax= df2$yaxis[12]), fill= df2$Colour[11]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[12], ymax= df2$yaxis[13]), fill= df2$Colour[12]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[13], ymax= df2$yaxis[14]), fill= df2$Colour[13]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[14], ymax= df2$yaxis[15]), fill= df2$Colour[14]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[15], ymax= df2$yaxis[16]), fill= df2$Colour[15]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[16], ymax= df2$yaxis[17]), fill= df2$Colour[16]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[17], ymax= df2$yaxis[18]), fill= df2$Colour[17]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[17], ymax= df2$yaxis[18]), fill= df2$Colour[18]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[18], ymax= df2$yaxis[19]), fill= df2$Colour[19]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[19], ymax= df2$yaxis[20]), fill= df2$Colour[20]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[20], ymax= df2$yaxis[21]), fill= df2$Colour[21]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[21], ymax= df2$yaxis[22]), fill= df2$Colour[22]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[22], ymax= df2$yaxis[23]), fill= df2$Colour[23]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[23], ymax= df2$yaxis[24]), fill= df2$Colour[24]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[24], ymax= df2$yaxis[25]), fill= df2$Colour[25]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[25], ymax= df2$yaxis[26]), fill= df2$Colour[26]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[26], ymax= df2$yaxis[27]), fill= df2$Colour[27]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[27], ymax= df2$yaxis[28]), fill= df2$Colour[28]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[28], ymax= df2$yaxis[29]), fill= df2$Colour[29]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[29], ymax= df2$yaxis[30]), fill= df2$Colour[30]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[30], ymax= df2$yaxis[31]), fill= df2$Colour[31]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[31], ymax= df2$yaxis[32]), fill= df2$Colour[32]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[32], ymax= df2$yaxis[33]), fill= df2$Colour[33]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[33], ymax= df2$yaxis[34]), fill= df2$Colour[34]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = 0.15, ymin= df2$yaxis[34], ymax= df2$yaxis[35]), fill= df2$Colour[35]) +
      ggplot2::xlab("") +
      ggplot2::xlim(0, 1) +
      ggplot2::ylab("Reactivity") +
      ggplot2::theme_classic() +
      ggplot2::labs(title = paste("N = (no data)\nor\nReactivity < ", abs_reactivity_threshold)) +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = 'black', size = 1.5),
                     axis.ticks = ggplot2::element_line(colour = "black", size = 1.5),
                     title = ggplot2::element_text(color = no_data, size =16),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.line.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(color = "Black", size = 16),
                     axis.title.x = ggplot2::element_text(color = "Black", size = 18),
                     axis.title.y = ggplot2::element_text(color = "Black", size = 18),
                     legend.text = ggplot2::element_text(color = "Black", size = 16),
                     legend.title = ggplot2::element_text(color = "Black", size = 16))
    if (write_legend){
      ggplot2::ggsave(filename = "Legend.pdf", path = getwd(), plot = legend, scale = 2.5, width = 2.85, height = 5, units = "cm", dpi = 300)
    }
    rxnplot <- ggplot2::ggplot(data = df, ggplot2::aes(x = xaxis, y = yaxis)) +
      ggplot2::geom_point(show.legend = FALSE, size =2, colour = df$Colour) +
      ggplot2::xlab("Residue") +
      ggplot2::ylab("Reactivity") +
      ggplot2::theme_classic()+
      ggplot2::theme(axis.line = ggplot2::element_line(colour = 'black', size = 1.5),
                     axis.ticks = ggplot2::element_line(colour = "black", size = 1.5),
                     axis.text.x = ggplot2::element_text(color = "Black", size = 16),
                     axis.text.y = ggplot2::element_text(color = "Black", size = 16),
                     axis.title.x = ggplot2::element_text(color = "Black", size = 18),
                     axis.title.y = ggplot2::element_text(color = "Black", size = 18),
                     legend.text = ggplot2::element_text(color = "Black", size = 16),
                     legend.title = ggplot2::element_text(color = "Black", size = 16))
    if (write_legend){
      ggplot2::ggsave(filename = "Rxn_plot.pdf", path = getwd(), plot = rxnplot, scale = 2.5, width = 7, height = 5, units = "cm", dpi = 300)
    }
  }
  data_frame$Labels <- c
  data_frame$Colors <- b
  output <- data_frame
  }

#'Generates a Stockholm file and a R2R meta file that can be read by R2R to draw a secondary structure.
#'
#'@param output Path to the output Stockholm file
#'@param data_frame A data_frame containing columns labeled "Nucleotide", "Dotbracket", "Labels", and "Colors".
#'@param RNA_name The name of the RNA that you are drawing
#'@param colors How R2R will draw reactivity colors. No colors = "NA". Colored letters = "letters". Colored circles = "circles".
#'@return A Stockholm formated file.
#' @export
r2easyR.write = function(output,
                   data_frame,
                   RNA_name = "default",
                   colors = "NA"){
  if(colors == "NA"){
    fileConn<-file(paste(output, ".sto", sep = ""))
    writeLines(c("# STOCKHOLM 1.0",
                 paste(RNA_name, gsub("T", "U", gsub(" ", "", gsub(",", "", toString(R.utils::capitalize(data_frame$Nucleotide))))), sep = "\t"),
                 paste("#=GC SS_cons", gsub(", ", "", toString(data_frame$Dotbracket)), sep = "\t"),
                 paste("#=GC R2R_LABEL", gsub(", ", "", toString(rep(".", length(data_frame$N)))), sep = "\t"),
                 paste("#=GC cons", gsub("T", "U", gsub(" ", "", gsub(",", "", toString(R.utils::capitalize(data_frame$Nucleotide))))), sep = "\t"),
                 paste("#=GC conss", gsub(", ", "", toString(rep("2", length(data_frame$Nucleotide)))), sep = "\t"),
                 paste("#=GC cov_SS_cons", gsub(", ", "", toString(rep("3", length(data_frame$Nucleotide)))), sep = "\t"),
                 paste("#=GF R2R tick_label_regular_numbering ", data_frame$N[1], " 10 firstNucNum ", data_frame$N[1], sep = ""),
                 "//"),
               fileConn)
    close(fileConn)
  }
  if(colors == "circles"){
    a <- c()
    n <- which(data_frame$Labels != "0")
    if (length(n) != 0){
      for (i in c(1:length(n))){
        a[i] <- paste("#=GF R2R circle_nuc #:", n[i]-1, " rgb:", col2rgb(data_frame$Colors[n[i]])[1,1], ",", col2rgb(data_frame$Colors[n[i]])[2,1], ",", col2rgb(data_frame$Colors[n[i]])[3,1], sep ="")
      }
      fileConn<-file(paste(output, ".sto", sep = ""))
      writeLines(c("# STOCKHOLM 1.0",
                   paste(RNA_name, gsub("T", "U", gsub(" ", "", gsub(",", "", toString(R.utils::capitalize(data_frame$Nucleotide))))), sep = "\t"),
                   paste("#=GC SS_cons", gsub(", ", "", toString(data_frame$Dotbracket)), sep = "\t"),
                   paste("#=GC R2R_LABEL", gsub(", ", "", toString(rep(".", length(data_frame$Labels)))), sep = "\t"),
                   paste("#=GC cons", gsub("T", "U", gsub(" ", "", gsub(",", "", toString(R.utils::capitalize(data_frame$Nucleotide))))), sep = "\t"),
                   paste("#=GC conss", gsub(", ", "", toString(rep("2", length(data_frame$Nucleotide)))), sep = "\t"),
                   paste("#=GC cov_SS_cons", gsub(", ", "", toString(rep("3", length(data_frame$Nucleotide)))), sep = "\t"),
                   gsub(", ", "\n",toString(unique(a))),
                   "#=GF R2R SetDrawingParam nucShrinkWithCircleNuc 1 pairBondScaleWithCircleNuc 1",
                   paste("#=GF R2R tick_label_regular_numbering ", data_frame$N[1], " 10 firstNucNum ", data_frame$N[1], sep = ""),
                   "//"),
                 fileConn)
      close(fileConn)
    }
    if (length(n) == 0){
      fileConn<-file(paste(output, ".sto", sep = ""))
      writeLines(c("# STOCKHOLM 1.0",
                   paste(RNA_name, gsub("T", "U", gsub(" ", "", gsub(",", "", toString(R.utils::capitalize(data_frame$Nucleotide))))), sep = "\t"),
                   paste("#=GC SS_cons", gsub(", ", "", toString(data_frame$Dotbracket)), sep = "\t"),
                   paste("#=GC R2R_LABEL", gsub(", ", "", toString(rep(".", length(data_frame$Labels)))), sep = "\t"),
                   paste("#=GC cons", gsub("T", "U", gsub(" ", "", gsub(",", "", toString(R.utils::capitalize(data_frame$Nucleotide))))), sep = "\t"),
                   paste("#=GC conss", gsub(", ", "", toString(rep("2", length(data_frame$Nucleotide)))), sep = "\t"),
                   paste("#=GC cov_SS_cons", gsub(", ", "", toString(rep("3", length(data_frame$Nucleotide)))), sep = "\t"),
                   "#=GF R2R SetDrawingParam nucShrinkWithCircleNuc 1 pairBondScaleWithCircleNuc 1",
                   paste("#=GF R2R tick_label_regular_numbering ", data_frame$N[1], " 10 firstNucNum ", data_frame$N[1], sep = ""),
                   "//"),
                 fileConn)
      close(fileConn)
    }
  }
  if(colors == "letters"){
    a <- c()
    for (i in c(1:length(data_frame$Colors))){
      a[i] <- paste("#=GF R2R nuc_color #:", i-1, " rgb:", col2rgb(data_frame$Colors[i])[1,1], ",", col2rgb(data_frame$Colors[i])[2,1], ",",col2rgb(data_frame$Colors[i])[3,1], sep ="")
    }
    fileConn<-file(paste(output, ".sto", sep = ""))
    writeLines(c("# STOCKHOLM 1.0",
                 paste(RNA_name, gsub("T", "U", gsub(" ", "", gsub(",", "", toString(R.utils::capitalize(data_frame$Nucleotide))))), sep = "\t"),
                 paste("#=GC SS_cons", gsub(", ", "", toString(data_frame$Dotbracket)), sep = "\t"),
                 paste("#=GC R2R_LABEL", gsub(", ", "", toString(rep(".", length(data_frame$Labels)))), sep = "\t"),
                 paste("#=GC cons", gsub("T", "U", gsub(" ", "", gsub(",", "", toString(R.utils::capitalize(data_frame$Nucleotide))))), sep = "\t"),
                 paste("#=GC conss", gsub(", ", "", toString(rep("2", length(data_frame$Nucleotide)))), sep = "\t"),
                 paste("#=GC cov_SS_cons", gsub(", ", "", toString(rep("3", length(data_frame$Nucleotide)))), sep = "\t"),
                 gsub(", ", "\n",toString(unique(a))),
                 paste("#=GF R2R tick_label_regular_numbering ", data_frame$N[1], " 10 firstNucNum ", data_frame$N[1], sep = ""),
                 "//"),
               fileConn)
    close(fileConn)
  }
  meta <- file(paste(output, ".r2r_meta", sep = ""))
  writeLines(c(paste(output, ".sto", "\t", "oneseq", "\t", RNA_name, sep = "")),
             meta)
  close(meta)
}
