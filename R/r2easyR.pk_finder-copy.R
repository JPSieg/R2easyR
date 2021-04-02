#'Identifies non-nested base pairs (pseudoknots) from an RNA secondary structure specified
#'by a connectivity table (.ct) formatted data frame.
#'
#'Capable of identifying multiple pseudoknots and pseudoknots within pseudoknots.
#'Returns a modified data frame with pseudoknotted base pairs deleted and a list
#'of vectors of the pseudoknot base pairs. Pseudoknots can be passed to the r2easyR.pknot_drawer.
#'Pseudoknots are troublesome for structure drawing programs because their non-nested nature
#'fools the drawing program into pairing the incorrect bases.
#'
#'
#'
#'@param ctdata_frame A data frame made by reading a .ct file into R with read.ct
#'@return A list of two elements, the first is the modified data frame with the pseudoknots deleted from the dotbracket column and the second is a list of vectors identifying pseudonotted bases.
#' @export
r2easyR.pk_finder = function(ctdata_frame){

  #df <- read.ct(ctdata_frame)

  df = ctdata_frame

  list_stems <- list()

  current_stem <-c()

  already_counted <- c()
  
  df$N <- as.character(df$N)
  df$BP <- as.character(df$BP)

  #this serves to gather all stems
  for (i in 1:length(df$N)){
    N_i <- as.integer(df$N[i])
    BP_i <- as.integer(df$BP[i])

    n <- min(as.integer(df$N[i]),as.integer(df$BP[i]))
    bp <- max(as.integer(df$N[i]),as.integer(df$BP[i]))

    if (i == length(df$N)){
      list_stems <- c(list_stems, list(current_stem))
      current_stem <- c()
    }
    if (BP_i != 0){
      if (n %in% already_counted == TRUE){
        next
      } else {
        already_counted <- c(already_counted, n)
      }

      pair <- c(as.integer(n), as.integer(bp))

      if (length(current_stem)==0){
        current_stem <- c(current_stem, list(pair))
      }else if (as.integer(pair[1])-1 == current_stem[[length(current_stem)]][1] && as.integer(pair[2])+1 == current_stem[[length(current_stem)]][2]){
        current_stem <- c(current_stem, list(pair))
      }else {
        list_stems <- c(list_stems, list(current_stem))
        current_stem <- c()
        current_stem <-list(pair)
      }

    }
  }
  #print(list_stems)

  #this serves to remove all fully nested hairpins
  indecies_to_keep <- c()
  for (k in 1:length(list_stems)) {
    stem = list_stems[k]
    #print(stem)
    stem[[1]][[1]]
    pk = FALSE
    already_counted <-c()
    stem_most_i = stem[[1]][[1]][1]
    stem_most_j = stem[[1]][[1]][2]
    #print('scanning stem...')
    #print(stem_most_i)
    #print(stem_most_j)
    for (i in 1:length(df$N)){
      n <- min(as.integer(df$N[i]),as.integer(df$BP[i]))
      bp <- max(as.integer(df$N[i]),as.integer(df$BP[i]))
      if (n != 0){
        if (n %in% already_counted == TRUE){
          next
        } else {
          already_counted <- c(already_counted, n)
        }
        if (n %in% stem_most_i:stem_most_j==TRUE && bp %in% stem_most_i:stem_most_j==FALSE) {
          #print('pseudoknot detected')
          indecies_to_keep <- c(indecies_to_keep, k)
          break
        }
        if (bp %in% stem_most_i:stem_most_j==TRUE && n %in% stem_most_i:stem_most_j==FALSE) {
          #print('pseudoknot detected')
          indecies_to_keep <- c(indecies_to_keep, k)
          break
        }
      }
    }
  }

  #print(indecies_to_keep)

  list_pk_stems <- list()

  for (i in indecies_to_keep){
    #print(i)
    list_pk_stems <- c(list_pk_stems, list_stems[i])
  }


  #print(list_pk_stems)

  list_competitors <- list()

  for (ii in 1:length(list_pk_stems)) {
    already_counted <-c()
    stem=list_pk_stems[[ii]]
    stem_most_i = stem[[1]][1]
    stem_most_j = stem[[1]][2]
    conflicting_bases <- c()
    for (i in 1:length(df$N)){
      n <- min(as.integer(df$N[i]),as.integer(df$BP[i]))
      bp <- max(as.integer(df$N[i]),as.integer(df$BP[i]))
      if (n != 0){
        if (n %in% already_counted == TRUE){
          next
        } else {
          already_counted <- c(already_counted, n)
        }
        if (n %in% stem_most_i:stem_most_j==TRUE && bp %in% stem_most_i:stem_most_j==FALSE) {
          #print('pseudoknot detected')
          conflicting_bases <- c(conflicting_bases, n, bp)
        }
        if (bp %in% stem_most_i:stem_most_j==TRUE && n %in% stem_most_i:stem_most_j==FALSE) {
          #print('pseudoknot detected')
          conflicting_bases <- c(conflicting_bases, n, bp)
        }
      }
    }
    #print('conflict')
    #print(conflicting_bases)
    conflicting_bases=list(conflicting_bases)
    list_competitors <- c(list_competitors, conflicting_bases)
  }

  #print('list')
  #print(list_competitors)

  #this serves to combine like-nested stems
  indecies_to_remove <- c()

  for (i in 1:length(list_pk_stems)) {
    one_struc_pair = list_pk_stems[[i]][[1]]
    one_struc_conflict = list_competitors[[i]]
    for (j in i:length(list_pk_stems)){
      two_struc_pair = list_pk_stems[[j]][[1]]
      two_struc_conflict = list_competitors[[j]]
      smae_conflict = TRUE
      for (iii in one_struc_conflict){
        if (iii %in% two_struc_conflict == FALSE){
          smae_conflict = FALSE
          break
        }
      }
      for (iii in two_struc_conflict){
        if (iii %in% one_struc_conflict == FALSE){

          smae_conflict = FALSE
          break
        }
      }

      #if (two_struc_pair[1] %in% one_struc_pair[1]:one_struc_pair[2] == TRUE && two_struc_pair[2] %in% one_struc_pair[1]:one_struc_pair[2] == TRUE && one_struc_pair != two_struc_pair && smae_conflict == TRUE){
      if (smae_conflict == TRUE && one_struc_pair != two_struc_pair){
        list_pk_stems[[i]] <- c(list_pk_stems[[i]], list_pk_stems[[j]])
        if (j %in% indecies_to_remove== FALSE) {
          indecies_to_remove <- c(indecies_to_remove, j)
          #print('how')
          #print(one_struc_conflict)
          #print(two_struc_conflict)
        }
      }
    }
  }

  #print(indecies_to_remove)
  #print('index')
  if (length(indecies_to_remove) != 0){
    for (i in 1:length(indecies_to_remove)) {
      list_pk_stems[[indecies_to_remove[[length(indecies_to_remove)+1-i]]]] <- NULL
      list_competitors[[indecies_to_remove[[length(indecies_to_remove)+1-i]]]] <- NULL
    }
  }

  #print(list_pk_stems)

  #this serves to decide which competing nests are pseudoknots

  pknots <- c()

  for (i in 1:length(list_pk_stems)) {
    one_struc_pair = list_pk_stems[[i]][[1]]
    for (j in i:length(list_pk_stems)){
      two_struc_pair = list_pk_stems[[j]][[1]]
      if (two_struc_pair[1] %in% one_struc_pair[1]:one_struc_pair[2] == TRUE && two_struc_pair[2] %in% one_struc_pair[1]:one_struc_pair[2] == FALSE && one_struc_pair != two_struc_pair){
        if (length(list_pk_stems[[i]]) >= length(list_pk_stems[[j]])){
          for (pk_pair in list_pk_stems[[j]]){
            pknots <- c(pknots, pk_pair[1])
          }
          for (pk_pair in list_pk_stems[[j]]){
            pknots <- c(pknots, pk_pair[2])
          }
        } else{
          for (pk_pair in list_pk_stems[[i]]){
            pknots <- c(pknots, pk_pair[1])
          }
          for (pk_pair in list_pk_stems[[i]]){
            pknots <- c(pknots, pk_pair[2])
          }
        }
      } else if (two_struc_pair[1] %in% one_struc_pair[1]:one_struc_pair[2] == FALSE && two_struc_pair[2] %in% one_struc_pair[1]:one_struc_pair[2] == TRUE && one_struc_pair != two_struc_pair){
        if (length(list_pk_stems[[i]]) >= length(list_pk_stems[[j]])){
          for (pk_pair in list_pk_stems[[j]]){
            pknots <- c(pknots, pk_pair[1])
          }
          for (pk_pair in list_pk_stems[[j]]){
            pknots <- c(pknots, pk_pair[2])
          }
        } else{
          for (pk_pair in list_pk_stems[[i]]){
            pknots <- c(pknots, pk_pair[1])
          }
          for (pk_pair in list_pk_stems[[i]]){
            pknots <- c(pknots, pk_pair[2])
          }
        }
      }
    }
  }
  #print(pknots)
  wait_pknots <- pknots
  pknots <- c()
  for (i in wait_pknots){
    if (i  %in% pknots == FALSE){
      pknots <- c(pknots, i)
    }
  }
  poly_pk = FALSE

  for (i in 1:(length(pknots)-1)){
    as.integer(pknots[i])
    if (as.integer(pknots[i]) - as.integer(pknots[i+1]) >0 && as.integer(pknots[i]) - as.integer(pknots[i+1]) != 1) {
      poly_pk=TRUE
    }
  }

  already_counted <- c()

  pknot_pairs <- c()
  pknot2 <-c()
  for (i in pknots){
    for (ii in 1:length(df$N)){
      if (i == df$N[ii]){
        if (df$N[ii] %in% already_counted == FALSE && df$BP[ii] != 0) {
          already_counted <- c(already_counted,df$N[ii], df$BP[ii])
          pair <- list(df$N[ii], df$BP[ii])
          pknot_pairs <- c(pknot_pairs,pair)
        }
      }
    }
    if (poly_pk == TRUE ){
    }
    for (i in  1:length(pknot_pairs)){
      if (i%%2 == 1){
        for (ii in 1:length(pknot_pairs)){
          if (ii%%2==1){
            pknot_pairs[i]
            pknot_pairs[i+1]
            if (as.integer(pknot_pairs[ii]) %in% as.integer(pknot_pairs[i]):as.integer(pknot_pairs[i+1]) == TRUE && as.integer(pknot_pairs[ii+1]) %in% as.integer(pknot_pairs[i]):as.integer(pknot_pairs[i+1]) == FALSE) {
              if (as.integer(pknot_pairs[ii]) %in% pknot2 == FALSE) {
                pknot2 <- c(pknot2, as.integer(pknot_pairs[ii]), as.integer(pknot_pairs[ii+1]))
              }
            }
          }
        }
      }
    }
  }


  list_pseudoknot_stems <- list()
  stem <- c()
  for (i in  1:(length(pknot_pairs))){
    if (i == length(pknot_pairs)){
      list_pseudoknot_stems <- c(list_pseudoknot_stems, list(stem))
      stem <- c()
    }
    if (i%%2 == 1){
      pair_i = as.integer(pknot_pairs[i])
      pair_j = as.integer(pknot_pairs[i+1])

      pair <- c(pair_i, pair_j)
      if (length(stem)==0){
        stem <- c(stem, c(pair))
      }else if (as.integer(pair[1])-1 == stem[length(stem)-1] && as.integer(pair[2])+1 == stem[length(stem)]){
        stem <- c(stem, c(pair))
        #print(stem)
      }else {
        list_pseudoknot_stems <- c(list_pseudoknot_stems, list(stem))

        stem <- c()
        stem <-c(pair)
      }
    }

  }
  #print(stem)

  #below serves to delete pseudoknotted base pairs from ct data frame

  df=R2easyR::add.dot.bracket(df)

  for (i in 1:length(list_pseudoknot_stems)){
    df$Dotbracket[list_pseudoknot_stems[[i]]] <- "."
    list_pseudoknot_stems[[i]] = sort(list_pseudoknot_stems[[i]])
  }
  
  print(list_pseudoknot_stems)

  output = list(df, list_pseudoknot_stems)

  names(output) = c("r2easyR.dataframe", "pknot.list")

  output=output

}
