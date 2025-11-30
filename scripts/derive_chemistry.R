################################################################################
#### Project: Metabolism traits
#### Title:   Function | Derive annotation chemistry
#### Author:  Tom Walker (thomas.walker@unine.ch)
#### Date:    19 April 2022
#### ---------------------------------------------------------------------------

derive_chemistry <- function(mtbs){
  ## Format input data ----
  # filter compound data to remove any missing smiles
  isPresent <- !is.na(mtbs$compounds$smiles) & mtbs$compounds$smiles != ""
  compoundsOK <- mtbs$compounds[isPresent, ]
  # collapse compounds to smiles level
  compounds <- compoundsOK %>%
    split(., .$smiles) %>%
    lapply(function(x){
      summarise(x, across(class:my_class, take_most_common))
    }) %>%
    # bind to data frame
    do.call(rbind, .) %>%
    as.data.frame %>%
    rownames_to_column("smiles")
  # collapse metabolite data to smiles level
  metabolites <- mtbs$presAbs[, isPresent] %>%
    # transpose and make data frame
    t %>% 
    as.data.frame %>%
    # add smiles information, group, take mean and ungroup
    mutate(smiles = compoundsOK$smiles) %>%
    group_by(smiles) %>%
    summarise(across(everything(), mean)) %>%
    ungroup
  
  ## Prepare data tables ----
  # subset and reorder compound data
  compReady <- metabolites %>%
    # select smiles code only
    select(smiles) %>%
    # join back to original data (effectively subset to sorted-unique)
    left_join(., compounds) %>%
    # generate smiles ID (easier formatting)
    mutate(SID = paste0("S", 1:nrow(.))) %>%
    # reorder columns and make data frame
    select(SID, smiles, class:my_class) %>%
    as.data.frame
  # re-transpose metabolite data
  mtbsReady <- metabolites %>%
    select(-smiles) %>%
    t %>%
    as.data.frame
  # change mtbs column names
  colnames(mtbsReady) <- compReady$SID
  # reprocess to presence absence
  mtbsReady[mtbsReady > 0] <- 1

  ## Derive chemistry ----
  # extract smiles and parse to correct format
  smilesOnly <- as.vector(compReady$smiles)
  names(smilesOnly) <- compReady$SID
  smilesParsed <- rcdk::parse.smiles(smilesOnly, omit.nulls = T)
  # extract parameters (all; need to uncomment descriptors above if subset)
  redundant <- c(2, 7, 8, 11, 15, 17, 18, 20, 21, 24, 29, 33:38, 41, 43:45)
  descriptors <- rcdk::get.desc.names()[-redundant]
  rawChem <- rcdk::eval.desc(smilesParsed, descriptors)
  # format output
  compOut <- compReady %>%
    left_join(., rownames_to_column(rawChem, "SID")) %>%
    select(SID, smiles:my_class, Fsp3:nAcid)
  # make rownames
  rownames(compOut) <- compOut$SID
  
  ## Populate chemistry across presence-absence matrix ----
  chemSpecies <- compOut %>%
    # select chemistry
    select(Fsp3:nAcid) %>%
    # apply populate chemistry function to all columns
    apply(2, function(x){
      populate_chemistry(
        mtbs = mtbsReady,
        trait = x,
        rows = nrow(mtbsReady),
        missing2na = T
      )
    }) %>%
    # make data frame
    as.data.frame
  # delete species where all absent (shouldn't exist)
  chemSpecies <- chemSpecies[, colSums(chemSpecies) > 0]

  ## Return ----
  out <- list(
    presAbs = mtbsReady,
    chemistry = chemSpecies,
    compounds = compOut
  )
  return(out)
}





