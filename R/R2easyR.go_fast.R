#'Prodeces R2R inputs in bulk with optimized layouts
#'
#'This function reads a folder of RNAstructure CT files and a folder of ".shape" formated reactivity data and produces
#'R2R inputs with optimized stem layouts in bulk. Enables you to make lots of RNA secondary structures with two lines of code.
#'
#'@param CT_folder Path to the directory containing .ct files
#'@param Reactivity_folder Path to the directory containing .shape files
#'@param  R2R_folder Path to the directory you want the R2R .sto and metafiles to print in
#'@param abs_reactivity_threshold Minimum threshold for reactivity data to be mapped to a color. Default = 0. Using a higher threshold prevents numerous, but low reactivity values from cluttering up the finished picture.
#'@param manual.scale Manually specify a min and max. Useful if you want to color reactivity from multiple transcripts using 1 scale. Default is "FALSE". To set the scale, use manual.scale = c(min, max)
#'@param RNA_name The name of the RNA that you are drawing
#'@param RNA_to_make_legend_on You will need a legend to interperet the results. This program will print one for you. It will also print a reaction plot for you to see what the pallets look like on data. Set this to determine which ".ct" file you want to print the reaction plot and legend on (Legends will be the same for every structure because all structures use the sane manual scale)
#'@export
r2easyR.go_fast = function(CT_folder = "CT",
                                  Reactivity_folder = "Reactivity",
                                  R2R_folder = "R2R_files",
                                  abs_reactivity_threshold = 0.01,
                                  manual.scale = c(0, 2),
                                  RNA_name = "RNA",
                                  RNA_to_make_legend_on = 1){
  print("Reading CT files")
  
  df <- lapply(paste(CT_folder, list.files(CT_folder), sep = "/"), R2easyR::read.ct)
  
  print("Done reading CT files")
  
  print("Reading reactivity data")
  
  Rxn <- lapply(paste(Reactivity_folder, list.files(Reactivity_folder), sep = "/"), R2easyR::read.shape)
  
  print("Done reading reactivity data")
  
  print("Adding dot bracket")
  
  df <- lapply(df, R2easyR::add.dot.bracket)
  
  print("Done adding dot bracket")
  
  print("Adding in reactivity data")
  
  for (i in 1:length(df)){
    df[[i]]$Reactivity <- Rxn[[i]]
  }
  
  print("Done adding in reactivity data")
  
  print("Generating palettes")
  
  palettes <- R2easyR::r2easyR.palettes()
  
  print("Done generating palettes")
  
  print("Writing R2R files")
  
  names <- c()
  
  for (i in 1:length(df)){
    print(paste(i, "out of", length(df), sep = " "))
    df[[i]] <- R2easyR::r2easyR.color(df[[i]],
                                      palette = palettes$YlOrRd.c,
                                      abs_reactivity_threshold = abs_reactivity_threshold,
                                      manual.scale = manual.scale,
                                      write_legend = FALSE)
    names[i] <- strsplit(list.files("CT")[i], split = ".ct")[[1]][1]
    R2easyR::r2easyR.write(paste(R2R_folder, names[i], sep = "/"),
                           df[[i]],
                           RNA_name = "mRNA",
                           colors = "circles")
  }
  
  print("Done writing R2R files")
  
  print("Generating legend and rxn plot")
  
  df[[RNA_to_make_legend_on]] <- R2easyR::r2easyR.color(df[[RNA_to_make_legend_on]],
                                                        palette = palettes$YlOrRd.c,
                                                        abs_reactivity_threshold = abs_reactivity_threshold,
                                                        manual.scale = manual.scale,
                                                        write_legend = TRUE)
  
  print("Done generating legend and rxn plot")
  
  print("Editing default stem layout")
  
  for (i in 1:length(df)){
    print(paste(i, "out of", length(df), sep = " "))
    R2easyR::r2easyR.stem_editor(paste(R2R_folder, "/", names[i], ".sto", sep = ""))
  }
  
  print("Done editing default stem layout")
  
  print("Writing the master R2R meta file")
  
  lines <- paste(paste(R2R_folder, "/", names, ".sto", sep = ""), "oneseq", "mRNA", sep = "\t")
  
  fileConn <- file("master.r2r_meta")
  writeLines(lines,
             fileConn)
  
  close(fileConn)
  
  print("Done writing the master R2R meta file")
  
  print("Done")
}

