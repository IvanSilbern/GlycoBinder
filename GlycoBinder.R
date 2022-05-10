################################# Copyright ####################################
#          
#<GlycoBinder: an R script to streamline glycoproteomics data processing>   
#    Copyright (C) <2019>  <Ivan Silbern [isilber@mpibpc.mpg.de]
#                           Kuan-Ting Pan [kapn@mpibpc.mpg.de]>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
################################################################################

ptm <- proc.time()

version <- "1.0.0"
message("\n***** GlycoBinder version ", version, " *****\n")

##### Load required packages #####

# function for checking / installing packages
install_packages <- function(packages, install = TRUE, verbose = verbose){
  # function checks if the required packages are installed and
  # tries to install missing packages if install = TRUE
  # if packages cannot be installed it stops the execution
  
  # check if all packages are installed
  not_installed_packages <- packages[!packages %in% installed.packages()]
  
  if(length(not_installed_packages) > 0){
    
    message("Following packages are not found\n", paste(not_installed_packages, collapse = " "))
    
    if(install){
      
      for(i in seq_along(not_installed_packages)){
        
        message("\nInstalling ",  not_installed_packages[i])
        try(install.packages(not_installed_packages[i]))
        
      }
      
    }
    
  }
  
  not_installed_packages <- packages[!packages %in% installed.packages()]
  if(length(not_installed_packages) > 0) {
    
    stop(paste("Following packages could not be installed, please conduct installation manually\n", not_installed_packages, collapse = " " ))
    
  }
  
}

# packages needed
packages_to_load <- c("future.apply", "data.table", "dplyr", "stringr")

# check / install packages
install_packages(packages_to_load, install = TRUE)

# attach packages
library("data.table")
library("future.apply")

# set future options and parameters
options(future.globals.maxSize = +Inf)

##### custom functions ######

# write logs
writeLog <- function(log_gb, new1 = "", new2 = ""){
  
  if(file.exists("Glycobinder.log.txt")) log_gb <- readLines("Glycobinder.log.txt") else log_gb <- character()
  writeLines(append(log_gb, c(paste0("[", Sys.time(), "] ", new1),
                              new2)), "Glycobinder.log.txt")
}

# collect arguments
collectArgs <- function(arg, collection, default = NA,
                        transform_fun = NULL, check_fun = NULL,
                        verbose = TRUE){
  
  if(any(collection == arg)) x <- collection[which(collection == arg) + 1] else {
    
    if(verbose) message("set ", arg, " to default: ", default)
    return(default)
    
  }
  if(!is.null(transform_fun)) x <- suppressWarnings(transform_fun(x))
  if(!is.null(check_fun)) check <- check_fun(x) else check <- TRUE
  if(check == TRUE) return(x) else {
    
    if(verbose) message("set ", arg, " to default: ", default)
    return(default)
    
  }
  
}

# read MGF files
readMgf    <- function(path){
  
  # function reads mgf file into a data table of scans
  # scans are separated  based on "BEGIN IONS" and "END IONS" lines
  # Data table contains:
  # scan_start = line number of "BEGIN IONS" string
  # scan_end = line number of "END IONS" string
  # ions_start = line number from which ion information starts
  # scan_number = number of the scan
  # scan_header = entire scan header as a character verctor
  # ions = matrix of ions m/z and corresponding intensities
  
  mgf             <- fread(file = path, header = FALSE, sep = NULL, colClasses = "character")
  start.positions <- stringr::str_which(mgf[[1]], "^BEGIN IONS$")
  end.positions   <- stringr::str_which(mgf[[1]], "^END IONS$")
  mgf_scan_nr     <- as.integer(stringr::str_match(stringr::str_subset(mgf[[1]], "^TITLE="), "\\.([0-9]+)\\.")[, 2])
  
  ionsBegin <- function(mgf, start, end){
    
    stringr::str_which(mgf[start:end], "^[0-9.]+\\s[0-9.]+$")[1] + (start - 1)
    
  }
  
  # determine ion positions
  ions.positions <- unlist(Map(ionsBegin, start = start.positions, end = end.positions, MoreArgs = list(mgf = mgf[[1]])))
  
  # remove empty scans
  if(any(is.na(ions.positions))){
    
    start.positions <- start.positions[!is.na(ions.positions)]
    end.positions   <- end.positions[!is.na(ions.positions)]
    mgf_scan_nr     <- mgf_scan_nr[!is.na(ions.positions)]
    ions.positions  <- ions.positions[!is.na(ions.positions)]
    
  }
  
  scan_headers   <- Map(function(mgf, start, end) mgf[start:end], start = start.positions, end = (ions.positions-1), MoreArgs = list(mgf = mgf[[1]]))
  ions           <- Map(function(mgf, start, end){
    
    temp <- matrix(as.numeric(stringr::str_split_fixed(mgf[start:end], " ", n = 2)), ncol = 2, byrow = FALSE)
    temp <- subset(temp, temp[, 2] != 0)
    return(temp)
    
    
  }, start = ions.positions,  end = (end.positions - 1), MoreArgs = list(mgf = mgf[[1]]))
  
  data.table(scan_start  = start.positions,
             scan_end    = end.positions,
             ions_start  = ions.positions,
             scan_number = mgf_scan_nr,
             scan_header = scan_headers,
             ions        = ions,
             key         = "scan_number")
  
}
readMgfHeader    <- function(path){
  
  # function reads mgf file into a data table of scans
  # it preserves only SCAN NUMBERS and SCAN HEADERS without ion information
  # scans are separated  based on "BEGIN IONS" and "END IONS" lines
  # Data table contains:
  # scan_number = number of the scan
  # scan_header = entire scan header as a character verctor
  
  mgf             <- fread(file = path, header = FALSE, sep = NULL, colClasses = "character")
  start.positions <- stringr::str_which(mgf[[1]], "^BEGIN IONS$")
  end.positions   <- stringr::str_which(mgf[[1]], "^END IONS$")
  mgf_scan_nr     <- as.integer(stringr::str_match(stringr::str_subset(mgf[[1]], "^TITLE="), "\\.([0-9]+)\\.")[, 2])
  
  ionsBegin <- function(mgf, start, end){
    
    stringr::str_which(mgf[start:end], "^[0-9.]+\\s[0-9.]+$")[1] + (start - 1)
    
  }
  
  ions.positions <- unlist(Map(ionsBegin, start = start.positions, end = end.positions, MoreArgs = list(mgf = mgf[[1]])))
  scan_headers   <- Map(function(mgf, start, end) mgf[start:end], start = start.positions, end = (ions.positions-1), MoreArgs = list(mgf = mgf[[1]]))
  
  data.table(scan_number = mgf_scan_nr,
             scan_header = scan_headers,
             key         = "scan_number")
  
}

# find marker ions
ion_match <- function(marker_ions,
                      mgf,
                      ion_match_tolerance,
                      tolerance_unit){
  
  dt <- data.table()
  for(i in 1:mgf[, .N]){
    
    ions      <- mgf$ions[[i]][, 1]
    intensity <- mgf$ions[[i]][, 2]
    delta_m <- abs(outer(ions, marker_ions$Mass, FUN = "-"))
    
    # provide tolerance level for each ms3 ion
    if(tolerance_unit == "ppm"){
      
      tolerances <- ion_match_tolerance*ions/10^6
      
    } else if(tolerance_unit == "Th"){
      
      tolerances <- rep(ion_match_tolerance, rep = length(ions))
      
    } else {
      
      message("ion matching tolerance is 1 ppm")
      tolerances <- ion_match_tolerance*ions/10^6
      
    }
    
    # which ions passed the tolerance cutoff
    passed <- delta_m < tolerances | dplyr::near(delta_m, tolerances)
    
    temp   <- data.table()
    temp   <- temp[, lapply(data.frame(passed), function(x) sum(intensity[x]))]
    dt <- rbind(dt, temp)
    
  }
  
  names(dt) <- paste0(marker_ions$Marker, "_", trunc(marker_ions$Mass))
  dt[, Scan := mgf$scan_number]
  return(dt)
  
}

# Form peptide groups based on Sequence Window
aggregateIds <- function(list_ids, start_id){
  
  # Function forms a closed group of indices around starting integer number.
  # It returns a vector of TRUE or FALSE values, where each TRUE/FALSE value corresponds to an element in the list.
  # TRUE means that the list element contains any indicies from the formed group of indicies.
  # Function takes a list of integer vectors (list_ids) and an integer (start_id) as a starting point.
  # It finds which list elements contain the starting integer number.
  # Other indices that are present within those list elements together with the starting number form a group of indices.
  # The algorithm continues to search for list elements that contain any of the numbers from the newly formed group of indices. 
  # If new list elements contain other indicies that are not a part of the group yet, they are added to the group.
  # Algorithm continues as long as the group of indices grows.
  # Once no new indices are entering the group, a closed group of indices within the data set is formed.
  # TRUE is returned for elements containing any of the indicies from the closed group of indices.
  # FALSE is returned for all other elements
  # Then it looks for the list elements that contain the starting integer number and other numbers
  # that were present in the same list elements
  # Arguments:
  # list_ids = list of integer vectors (e.g. modified peptide ids belonging to each of the sequence windows)
  # start_id = integer number (id) to start with
  # Value:
  # Logical vector, TRUE corresponds to elements in the list_ids that contain any of ids from the formed group of ids
  
  temp_id <- start_id
  length_before <- 0
  length_after  <- 1
  
  while(length_before < length_after){
    
    length_before <- length(temp_id)
    
    temp_sw <- unlist(lapply(list_ids, function(x) any(x %in% temp_id)))
    temp_id <- unique(unlist(list_ids[temp_sw]))
    
    length_after <- length(temp_id)
    
  }
  
  return(temp_sw)
  
}

# run external tools through the interface to the command line
run_pglyco    <- function(file_name, out_dir, pglyco_config, pglyco_path){
  
  # Function runs pGlyco through a system call
  # Arguments:
  # file_name = name of the file
  # out_idr = output directory
  # pglyco_config = path to pglyco config file
  # pglyco_path = path to pglyco executable
  # Value:
  # NULL
  
  ### change working directory to pglyco location
  if(!dir.exists(pglyco_path)) stop("Cannot find folder location of pGlyco")
  
  wd <- getwd()
  setwd(pglyco_path)
  
  ### change pglyco_config template
  
  # number of processes
  pglyco_config[grepl("^process=", pglyco_config)]    <- paste0("process=1")
  
  # total number of files to analyze
  pglyco_config[grepl("^spectrum_total=", pglyco_config)] <- paste0("spectrum_total=1")
  
  #output_dir
  pglyco_config[grepl("^output_dir=", pglyco_config)] <- paste0("output_dir=", normalizePath(out_dir, winslash = "\\"))
  
  # file to process
  pglyco_config <- pglyco_config[1:which(grepl("^spectrum_total=", pglyco_config))] # should be the last line before file paths
  process_file  <- paste0(wd, "/", file_name)
  if(!file.exists(process_file)) { # check if the file exists 
    
    writeLines("Not Finished", paste0(out_dir, "done.txt"))
    stop(paste0("Cannot find file to process ", process_file))
    
  }
  pglyco_config <- c(pglyco_config, paste0("file1", "=", process_file))
  
  # path to pglyco config file
  pglyco_config_file <- paste0(out_dir, "/pGlyco_task.pglyco")
  
  # write pGlyco config file
  writeLines(pglyco_config, pglyco_config_file)
  
  # prepare system calls
  cmd_pglycodb  <- paste0('pGlycodb.exe ',         '"', pglyco_config_file, '"')
  cmd_pglycofdr <- paste0('pGlycoFDR.exe', " -p ", '"', pglyco_config_file, '"', " -r ", '"', paste0(out_dir, "/pGlycoDB-GP.txt"), '"')
  cmd_pglycopro <- paste0('pGlycoProInfer.exe ',   '"', pglyco_config_file, '"')
  
  shell(cmd = paste0(cmd_pglycodb, ' && ', cmd_pglycofdr, ' && ', cmd_pglycopro, ' & echo "done" >> ', '"', out_dir, '/done.txt', '"'), wait = FALSE, ignore.stdout = TRUE)
  
  # switch back to the original directory
  setwd(wd)
  
}

# add core fucosylation information
add_corefuc <- function(df, scans){
  
  ids <- lapply(stringr::str_split(df$pGlyco_ids, pattern = "[;/]"), as.integer)
  CoreFuc <- unlist(lapply(ids, function(x){
    
    paste0(scans$CoreFuc[scans$id %in% x], collapse = ";")
    
  }))
  df[, CoreFuc := CoreFuc]
  df
  
}

# check if any structure contain fucose
check_f_struct <- function(df, scans){
  
  ids <- lapply(stringr::str_split(df$pGlyco_ids, pattern = "[;/]"), as.integer)
  f_struct <- unlist(lapply(ids, function(x){
    
    any(grepl("F", scans$PlausibleStruct[scans$id %in% x]))
    
  }))
  
  df[, FucoseStruct := f_struct]
  df
  
}

# check glycans regarding the core fucose type
check_core_fucose <- function(df, scans){
  
  ids <- lapply(stringr::str_split(df$pGlyco_ids, pattern = "[;/]"), as.integer)
  core_fucose <- unlist(lapply(ids, function(x){
    
    if(all(scans$CoreFuc[scans$id %in% x] %in% c("11", "12"))) return("Yes")
    else if(all(scans$CoreFuc[scans$id %in% x] %in% c("0"))) return("No")
    else return("Ambiguous")
    
  }))
  
  df[, CoreFucoseOnly := core_fucose]
  df
  
}

# aggregate parent peak areas
addPeakArea <- function(df, scans, FUN = "sum"){
  
  ids <- lapply(stringr::str_split(df$pGlyco_ids, pattern = "[;/]"), as.integer)
  area <- unlist(lapply(ids, function(x){
    
    ppa <- scans[scans$id %in% x, c("RawName", "Peptide", "Glycan(H,N,A,G,F)",
                                    "PrecursorMZ", "PrecursorCharge", "ParentPeakArea")]
    ppa[, PrecursorMZ := signif(PrecursorMZ, 5)]
    ppa <- ppa[order(-ParentPeakArea)]
    ppa <- ppa[!duplicated(ppa[, -c("ParentPeakArea")])]$ParentPeakArea
    ParentPeakArea  <- eval(call(get("FUN"), ppa, na.rm = TRUE))
    return(ParentPeakArea)
    
  }))
  
  df[, ParentPeakArea := area]
  return(df)
  
}

# convert reporter ion intensities into %
percentIntensity <- function(df, int_names){
  
  df[, TotalRepIonIntensity := sum(unlist(.SD), na.rm = TRUE),   .SD = int_names, by = "id"]
  df[, paste0(int_names, "Percent") := .SD/TotalRepIonIntensity, .SD = int_names, by = "id"]
  df[, paste0(int_names, "Percent") := 100*.SD/TotalRepIonIntensity, .SD = int_names, by = "id"]
  df[, paste0(int_names, "_ParentPeakArea") := .SD*ParentPeakArea*0.01, .SD = paste0(int_names, "Percent"), by = "id"]
  return(df)
  
}

# Define Glycan/Glycan Antenna types
defineGlycoType <- function(df){
  
  df[, composition := lapply(stringr::str_split(df$`Glycan(H,N,A,G,F)`, " "), as.integer)]
  df[, Hex    := unlist(lapply(composition, "[", 1))]
  df[, HexNAc := unlist(lapply(composition, "[", 2))]
  df[, NeuAc  := unlist(lapply(composition, "[", 3))]
  df[, NeuGc  := unlist(lapply(composition, "[", 4))]
  df[, Fuc    := unlist(lapply(composition, "[", 5))]
  
  df[, GlycanType := "Other"]
  df[, GlycanAntennaType := "Other"]
  df[(HexNAc == 2 &
        Hex >= 1 &
        Hex <= 4 &
        NeuAc == 0), GlycanAntennaType := "Paucimannose"]
  df[(HexNAc == 2 &
        Hex >= 5 &
        Hex <= 9 &
        NeuAc == 0 &
        Fuc == 0), GlycanAntennaType := "High-mannose"]
  df[(HexNAc == 2 &
        Hex >= 10 &
        Hex <= 12 &
        NeuAc == 0 &
        Fuc == 0), GlycanAntennaType := "Initiation"]
  df[(HexNAc == 3 &
        Hex >= 3), GlycanAntennaType := "Hybrid/A1"]
  df[(HexNAc == 4 &
        Hex >= 3), GlycanAntennaType := "A2/A1B"]
  df[(HexNAc == 5 &
        Hex >= 3), GlycanAntennaType := "A3/A2B"]
  df[(HexNAc >= 6 &
        Hex >= 3), GlycanAntennaType := "A4/A3B"]
  df[, GlycanType := GlycanAntennaType]
  df[GlycanAntennaType %in% c("A2/A1B", "A3/A2B", "A4/A3B"), GlycanType := "Complex"]
  
  df <- df[, .SD, .SDcols = c(names(df)[1:which(names(df) == "Glycan(H,N,A,G,F)")],
                              "GlycanType", "GlycanAntennaType",
                              setdiff(names(df)[(which(names(df) == "Glycan(H,N,A,G,F)")+1):length(df)], c("GlycanType", "GlycanAntennaType"))
  )]
  
  return(df[, -c("composition", "Hex", "HexNAc", "NeuAc", "NeuGc", "Fuc")])
  
}

# Define Glycotope for "complex" or "hybrid" glycan types
defineGlycoTope <- function(df, limit_intensity = 0.1){
  
  if(!all(c("NeuAc_H2O_274",
            "NeuAc_292",
            "NeuAc_Hex_454",
            "NeuAc_HexNAc_495",
            "Hex_HexNAc_Fuc_512",
            "NeuAc_NeuAc_583",
            "NeuAc_Hex_HexNAc_657",
            "Hex_HexNAc_2Fuc_658",
            "NeuAc_Hex_HexNAc_Fuc_803",
            "2NeuAc_Hex_HexNAc_948",
            "3NeuAc_Hex_NexNAc_1239") %in% names(df))){
    
    message("Skip Glycotope checking")
    df[, Glycotope_A := ""]
    df[, Glycotope_F := ""]
    return(df)
    
  }
  
  df[, Glycotope_F := ""]
  df[, Glycotope_A := ""]
  df[GlycanType %in% c("Complex", "Hybrid/A1") &
       NeuAc_Hex_HexNAc_Fuc_803 > limit_intensity &
       (NeuAc_H2O_274 > limit_intensity | NeuAc_292 > limit_intensity),
     Glycotope_F := "AHN(F)"]
  df[GlycanType %in% c("Complex", "Hybrid/A1") &
       NeuAc_Hex_HexNAc_Fuc_803 < limit_intensity &
       Hex_HexNAc_2Fuc_658 > limit_intensity,
     Glycotope_F := "H(F)N(F)"]
  df[GlycanType %in% c("Complex", "Hybrid/A1") &
       NeuAc_Hex_HexNAc_Fuc_803 < limit_intensity &
       Hex_HexNAc_2Fuc_658 < limit_intensity &
       Hex_HexNAc_Fuc_512 > limit_intensity,
     Glycotope_F := "HN(F)"]
  df[Glycotope_F != "" &
       (NeuAc_H2O_274 > limit_intensity | NeuAc_292 > limit_intensity) &
       `3NeuAc_Hex_NexNAc_1239` > limit_intensity,
     Glycotope_A := "AAAHNorAAHN(A)"]
  df[Glycotope_A == "AAAHNorAAHN(A)" &
       NeuAc_HexNAc_495 > limit_intensity,
     Glycotope_A := "AAHN(A)"]
  df[Glycotope_F != "" &
       (NeuAc_H2O_274 > limit_intensity | NeuAc_292 > limit_intensity) &
       `3NeuAc_Hex_NexNAc_1239` < limit_intensity &
       `2NeuAc_Hex_HexNAc_948` > limit_intensity,
     Glycotope_A := "AHN(A)orAAHN"]
  df[Glycotope_A == "AHN(A)orAAHN" &
       NeuAc_HexNAc_495 > limit_intensity &
       NeuAc_NeuAc_583 < limit_intensity,
     Glycotope_A := "AHN(A)"]
  df[Glycotope_A == "AHN(A)orAAHN" &
       NeuAc_HexNAc_495 < limit_intensity &
       NeuAc_NeuAc_583 > limit_intensity,
     Glycotope_A := "AAHN"]
  df[Glycotope_F != "" &
       (NeuAc_H2O_274 > limit_intensity | NeuAc_292 > limit_intensity) &
       `3NeuAc_Hex_NexNAc_1239` < limit_intensity &
       `2NeuAc_Hex_HexNAc_948` < limit_intensity &
       NeuAc_Hex_HexNAc_657 > limit_intensity,
     Glycotope_A := "AHNorHN(A)"]
  df[Glycotope_A == "AHNorHN(A)" &
       NeuAc_HexNAc_495 > limit_intensity &
       NeuAc_Hex_454 < limit_intensity,
     Glycotope_A := "HN(A)"]
  df[Glycotope_A == "AHNorHN(A)" &
       NeuAc_HexNAc_495 < limit_intensity &
       NeuAc_Hex_454 > limit_intensity,
     Glycotope_A := "AHN"]
  
  df[, Glycotope := Glycotope_F]
  df[Glycotope_F != "" & Glycotope_A != "", Glycotope := paste0(Glycotope_F, " & ", Glycotope_A)]
  
  df <- df[, .SD, .SDcols = c(names(df)[1:which(names(df) == "GlycanAntennaType")],
                              "Glycotope", "Glycotope_F", "Glycotope_A",
                              setdiff(names(df)[(which(names(df) == "GlycanAntennaType")+1):length(df)],
                                      c("Glycotope", "Glycotope_F", "Glycotope_A"))
  )]
  
  return(df[, -c("Glycotope")])
  
}

#check for glycotope conflict
glycotopeConflict <- function(df){
  
  df[, composition := stringr::str_split(`Glycan(H,N,A,G,F)`, " ")]
  df[, n_neuac := as.integer(unlist(lapply(composition, "[", 3)))]
  df[, n_fuc   := as.integer(unlist(lapply(composition, "[", 5)))]
  # df[, Glycotope_F := unlist(lapply(stringr::str_split(Glycotope, " & "), "[", 1))]
  # df[, Glycotope_A := unlist(lapply(stringr::str_split(Glycotope, " & "), "[", 2))]
  
  df[, Glycotope_F_unmatched := ""]
  df[, Glycotope_A_unmatched := ""]
  df[Glycotope_F %in% c("HN(F)", "AHN(F)") &
       !(CoreFuc %in% c(11, 12)) &
       n_fuc < 1, Glycotope_F_unmatched := "+"]
  df[Glycotope_F %in% c("H(F)N(F)") &
       !(CoreFuc %in% c(11, 12)) &
       n_fuc < 2, Glycotope_F_unmatched := "+"]
  df[Glycotope_F %in% c("HN(F)", "AHN(F)") &
       CoreFuc %in% c(11, 12) &
       n_fuc < 2, Glycotope_F_unmatched := "+"]
  df[Glycotope_F %in% c("H(F)N(F)") &
       CoreFuc %in% c(11, 12) &
       n_fuc < 3, Glycotope_F_unmatched := "+"]
  df[grepl("A", Glycotope_F) &
       n_neuac < 1, Glycotope_F_unmatched := "+"]
  df[grepl("A", Glycotope_A) &
       n_neuac < 1, Glycotope_A_unmatched := "+"]
  df[Glycotope_A %in% c("AHN(A)", "AAHN", "AHN(A)orAAHN") &
       n_neuac < 2, Glycotope_A_unmatched := "+"]
  df[Glycotope_A %in% c("AAHN(A)", "AAAHN", "AAAHNorAAHN(A)") &
       n_neuac < 3, Glycotope_A_unmatched := "+"]
  
  df <- df[, -c("composition", "n_neuac", "n_fuc")]
  df <- df[, .SD, .SDcols = c(names(df)[1:which(names(df) == "Glycotope_A")],
                              "Glycotope_F_unmatched", "Glycotope_A_unmatched",
                              setdiff(names(df)[(which(names(df) == "Glycotope_A")+1):length(df)],
                                      c("Glycotope_F_unmatched", "Glycotope_A_unmatched"))
  )]
  
  return(df)
  
}

##### define some variables #####

raw_file_extension      <- "raw"
supported_reporter_ions <- c("TMT0",
                             "TMT2",
                             "TMT6",
                             "TMT10",
                             "TMT11",
                             "TMT16",
                             "iTRAQ4",
                             "iTRAQ8",
                             "not_labeled")

# parallelizing strategy for the "future" package
plan_strategy <- "multisession"

# pglyco
reversed_regexpr <- "^REV_"   # decoy sequences 
pglyco_separator <- "/"       # separator in pGlyco

# default config file for pGlyco
pglyco_config_default <- c("[version]",
                           "pGlyco_version=pGlyco_Tag20190101",
                           "[flow]",
                           "glyco_type=N-Glyco",
                           "glycandb=pGlyco.gdb",
                           "process=1",
                           "output_dir=C:\\",
                           "[protein]",
                           "fasta=",
                           "enzyme=Trypsin_KR-C",
                           "digestion=not_support_in_cur_version",
                           "max_miss_cleave=2",
                           "max_peptide_len=40",
                           "min_peptide_len=6",
                           "max_peptide_weight=4000",
                           "min_peptide_weight=600",
                           "[modification]",
                           "fix_total=3",
                           "fix1=Carbamidomethyl[C]",
                           "fix2=",
                           "fix3=",
                           "max_var_modify_num=3",
                           "var_total=1",
                           "var1=Oxidation[M]",
                           "[search]",
                           "search_precursor_tolerance=10",
                           "search_precursor_tolerance_type=ppm",
                           "search_fragment_tolerance=20",
                           "search_fragment_tolerance_type=ppm",
                           "spec_file_type=mgf",
                           "spectrum_total=")

# pglyco3_config_default <- c("[version]",
#                             "pGlyco_version=pGlyco3.0.rc2",
#                             "pGlyco_type=pGlycoDB",
#                             "[ini]",
#                             "glycoini=glyco.ini",
#                             "modini=modification.ini",
#                             "[glycan]",
#                             "glycan_type=N-Glycan",
#                             "glycan_db=pGlyco-N-Human.gdb",
#                             "glycan_fix_mod=",
#                             "glycan_var_mod=",
#                             "max_var_mod_on_glycan=1",
#                             "max_glycan_db_size=100000",
#                             "[protein]",
#                             "fasta=", 
#                             "enzyme=Trypsin KR _ C",
#                             "digestion=specific",
#                             "max_miss_cleave=2",
#                             "min_peptide_len=6",
#                             "max_peptide_len=40",
#                             "min_peptide_weight=600",
#                             "max_peptide_weight=4000",
#                             "protein_fix_mod=Carbamidomethyl[C]",
#                             "protein_var_mod=Oxidation[M],Acetyl[ProteinN-term]",
#                             "max_var_mod_on_peptide=2",
#                             "[search]",
#                             "precursor_tolerance=10",
#                             "precursor_tolerance_type=ppm",
#                             "fragment_tolerance=20",
#                             "fragment_tolerance_type=ppm",
#                             "fragmentation_type=HCD",
#                             "spec_file_type=raw",
#                             "spectrum_total=1",
#                             "file1=",
#                             "top_n_peaks=300",
#                             "output_top_n=5",
#                             "percolator=0",
#                             "FDR=0.01",
#                             "FMM_for_peptide_FDR=0",
#                             "pGlycoSite_check_db=1",
#                             "output_dir=")

marker_ions <- c("Mass;Marker;Formula
                  109.028;Hex_fragm;C6H4O2
                  115.039;Hex_fragm;C5H6O3
                  126.055;HexNAc_fragm;C6H7O2N1
                  127.039;Hex_2xH2O;C6H6O3
                  138.055;HexNAc_fragm;C7H7O2N1
                  144.066;HexNAc_fragm;C6H9O3N1
                  163.060;Hex;C6H10O5
                  168.066;HexNAc_2xH2O;C8H9O3N1
                  186.076;HexNAc_2xH2O;C8H11O4N1
                  204.087;HexNAc;C8H13O5N1
                  274.092;NeuAc_H2O;C11H15O7N1
                  290.087;NeuGc_H2O;C11H15O8N1
                  292.103;NeuAc;C11H17O8N1
                  308.098;NeuGc;C11H17O9N1
                  366.140;Hex_HexNAc;C14H23O10N1
                  454.156;NeuAc_Hex;;
                  495.182;NeuAc_HexNAc;;
                  512.197;Hex_HexNAc_Fuc;;
                  583.198;NeuAc_NeuAc;;
                  657.235;NeuAc_Hex_HexNAc;C25H40O18N2;
                  673.230;Hex_HexNAc_NeuGc;C25H40O19N2
                  658.255;Hex_HexNAc_2Fuc;;
                  803.293;NeuAc_Hex_HexNAc_Fuc;;
                  948.330;2NeuAc_Hex_HexNAc;;
                  1239.426;3NeuAc_Hex_NexNAc;;
                 ")

# base glycotope types
base_glycotopes <- c("HN(F)", "AHN", "HN(A)", "H(F)N(F)", "AHN(F)", "AHN(A)", "AAHN", "AAHN(A)", "AAAHN",
                     "AHNorHN(A)", "AHN(A)orAAHN", "AAAHNorAAHN(A)")

# variables that will be set through arguments
verbose                   <- TRUE   # verbosity
nr_threads                <- 2      # number of threads to use 
ion_match_tolerance       <- 1      # width of the tolerance window for ion matching
tolerance_unit            <- "ppm"  # unit for the tolerance window
reporter_ion              <- "TMT6" # type of reporter ions
seq_wind_size             <- 7      # number of aa from each side of the modified residue
pglyco_fdr_threshold      <- 0.02  # total FDR cutoff
report_intermediate_files <- FALSE
skip_marker_ions          <- FALSE
parent_area_fun           <- "sum"

##### collect arguments and set working directory #####
args <- (commandArgs(TRUE))

# log information
message("Provided arguments:\n",
        paste(as.character(args), collapse = "\n"), "\n")
writeLog(log_gb, "Glycobinder started with arguments:", paste(as.character(args), collapse = " "))

# working directory
wd <- collectArgs("--wd", args, default = getwd(),
                  check_fun = function(x) all(length(x) == 1 & dir.exists(x)))
setwd(wd)

##### find raw files in the working directory #####
raw_file_names <- list.files(pattern = paste0("\\.", raw_file_extension, "$"))

if(length(raw_file_names) == 0) stop(paste0("Cannot find .", raw_file_extension, " files in the specified directory"))
message(paste0("Found raw files:\n", paste(raw_file_names, collapse = "\n")))

if(any(grepl("\\.", gsub(paste0(".", raw_file_extension), "", raw_file_names)))) stop("Dots . are not allowed in raw file names (except .raw file extension)")


##### collect further arguments #####

# verbose
if(any(grepl("--verbose", args))) verbose <- TRUE else verbose <- FALSE

# tolerance unit for merging MS2 and MS3 spectra
tolerance_unit <- collectArgs("--tol_unit", args, default = "ppm", 
                              check_fun = function(x) all(length(x) == 1 & x %in% c("ppm", "Th")))

# tolerance for merging MS2 and MS3 spectra
ion_match_tolerance <- collectArgs("--match_tol", args, default = 1,
                                   transform_fun = function(x) as.integer(x),
                                   check_fun = function(x) all(length(x) == 1 & !is.na(x) & x != 0))

# reporter ion type parameter for RawTools
reporter_ion <- collectArgs("--reporter_ion", args, default = NA,
                            check_fun = function(x) all(length(x) == 1 & x %in% supported_reporter_ions))
if(is.na(reporter_ion)) stop("Reporter ion type is not provided or not correct.")


# pglyco fdr threshold
pglyco_fdr_threshold <- collectArgs("--pglyco_fdr_threshold", args, default = 0.02,
                                    transform_fun = function(x) as.numeric(x),
                                    check_fun = function(x) all(length(x) == 1 & x & !is.na(x) & x > 0 & x <= 1))

# define sequence window size, +/- n amino acids around modified residue
seq_wind_size <- collectArgs("--seq_wind_size", args, default = 7,
                             transform_fun = function(x) as.numeric(x),
                             check_fun = function(x) all(length(x) == 1 & x & !is.na(x) & x >= 1 & x < 100))

# function for aggregating parent peak areas
parent_area_fun <-  collectArgs("--parent_area_fun", args, default = "sum",
                                check_fun = function(x) all(length(x) == 1 & x %in% c("sum", "median", "mean", "max", "min")))

# number of available threads
nr_of_processors <- shell("set NUMBER_OF_PROCESSORS", intern = TRUE)
nr_of_processors <- gsub("NUMBER_OF_PROCESSORS=", "", nr_of_processors)
nr_of_processors <- as.integer(nr_of_processors)

nr_threads <- collectArgs("--nr_threads", args, default = max(nr_of_processors - 2, 1), verbose = verbose,
                          transform_fun = function(x) as.integer(x),
                          check_fun = function(x) all(length(x) == 1 & !is.na(x) & x >= 1 & x <= nr_of_processors))
if(nr_threads > length(raw_file_names)) nr_threads <- length(raw_file_names)

# perform second pglyco search on the reduced fasta file
if(any(grepl("--no_second_search", args))) second_search <- FALSE else second_search <- TRUE

#report intermediate results
if(any(grepl("--report_intermediate_files", args))) report_intermediate_files <- TRUE else report_intermediate_files <- FALSE

# skip search for marker ions
if(any(grepl("--skip_marker_ions", args))) skip_marker_ions <- TRUE


##### read fasta file #####
local({
  
  # find fasta file
  fasta_files <- list.files(pattern = "\\.[Ff][Aa][Ss][Tt][Aa]$")
  if(length(fasta_files) == 0) stop("Cannot find any .fasta files")
  if(verbose) message("Read fasta file: ", paste(fasta_files, collapse = "\n"))
  if(length(fasta_files) > 1){
    
    fasta.name <<- "combined.fasta"
    
  } else {
    
    fasta.name <<- fasta_files
    
  }  
  
  readFASTA <- function(fasta_files){

    fasta <- character()
    for(i in seq_along(fasta_files)){

      fasta <- c(fasta, fread(fasta_files[i], sep = NULL, header = FALSE)[[1]])

    }

    # extract header positions
    # Lines containing headers
    ff_header_pos <- which(grepl("^>", fasta))
    ff_headers    <- fasta[ff_header_pos]

    fasta[ff_header_pos] <- "___"
    fasta <- paste0(fasta, collapse = "")
    fasta <- unlist(stringr::str_split(fasta, "___"))
    fasta <- fasta[-1]
    fasta2 <- vector("character", length(ff_header_pos) * 2)
    fasta2[seq(1, length(fasta2), 2)] <- ff_headers
    fasta2[seq(2, length(fasta2), 2)] <- fasta
    
    return(fasta2)

  }
  
  fasta         <<- readFASTA(fasta_files)
  ff_header_pos <<- seq(1, length(fasta), 2)
  ff_headers    <<- stringr::str_split_fixed(fasta[ff_header_pos],
                                            pattern = " ", n = 2)[, 1]
  ff_headers    <<- gsub("^>", "", ff_headers)
  fasta.N2J     <- fasta
  fasta.N2J[-ff_header_pos] <- gsub("N(.[STC])", "J\\1", fasta.N2J[-ff_header_pos])
  
  # write down the modified fasta file
  writeLines(fasta.N2J, paste0(fasta.name, ".N2J")) # don't use 'fwrite' function: it causes problems with pGlyco

})
##### run RawTools #####

if(!all(file.exists(paste0("rawtools_output\\",
                           raw_file_names,
                           "_Matrix.txt")))) {
  
  local({
    
    # output folder
    if(!dir.exists("rawtools_output")) dir.create("rawtools_output")
      
    # prepare a string of arguments
    RawTools_args    <- c('-parse',
                          '-d', paste0('"', wd, '"'),
                          '-out', paste0('"', wd, '/rawtools_output', '"'),
                          '-R', '-u')
    
    if(reporter_ion != "not_labeled"){
      
      RawTools_args <- c( RawTools_args, '-q', '-r', reporter_ion)
      
    }
    
    if(verbose) message(paste0("Running Rawtools with arguments ",
                        paste(RawTools_args, collapse = " ")))
    
    # run
    system2(command = "RawTools", args = RawTools_args, wait = TRUE)
    
  })
  
} else {
  
  message("Found all _Matrix.txt files. Skip RawTools processing")
  
}

# check which files contain ms3 spectra
contain_ms3 <- logical()
for(i in seq_along(raw_file_names)){
  
  contain_ms3[i] <- any(names(fread(paste0("rawtools_output\\",
                                           raw_file_names[i],
                                           "_Matrix.txt"),
                                    nrows = 1)) == "MS3ScanNumber")
  
}

##### check if output already exists #####

# check should be performed after the raw tools processing is finished
# because it finds out which raw files contain MS3 scans
local({
  
  msconv_output <- gsub("\\.raw$", "\\.mgf",raw_file_names[contain_ms3])
  msconv_output <- paste0("msconvert_output/", msconv_output)
  output_msconvert_exists <<- all(file.exists(msconv_output))
  
  pparse_mod <- gsub("\\.raw$", "_pParse_mod.mgf", raw_file_names)
  pparse_mod <- paste0("pparse_output/", pparse_mod)
  output_pparse_exists <<- all(file.exists(pparse_mod))
  
  pparse_out <- list.files(path = "pparse_output", pattern = "_[A-Z]+FT\\.mgf")
  pparse_out <- gsub("_[A-Z]+FT\\.mgf", "", pparse_out)
  output_pparse_exists <<- all(gsub("\\.raw$", "", raw_file_names) %in%
                                 pparse_out)
  
  pparsemod_out <- gsub("\\.raw$", "_pParse_mod.mgf", raw_file_names)
  pparsemod_out <- paste0("pparse_output/", pparsemod_out)
  output_pparsemod_exists <<- all(file.exists(pparsemod_out))
  
  output_pglyco1_exists   <<- any(grepl("-Pro.txt$",
                                        list.files(path = "pglyco_output")))
  output_pglyco2_exists   <<- any(grepl("-Pro2.txt",
                                        list.files(path = "pglyco_output")))

})

##### run MSConvert #####

if(any(contain_ms3) && # there are files with ms3 spectra 
   
   # mgf files for ms3 spectra containing raw files do not exist
   !output_msconvert_exists &&
   
   # pParse_mod.mgf do not exist
   !output_pparsemod_exists &&
   
   # no second search and no pglyco 1 ouput or
   # second search and no pglyco 2 output
   ((!second_search && !output_pglyco1_exists) ||                  
    (second_search  && !output_pglyco2_exists))
   
){
  
  # any(contain_ms3) == TRUE
  # output_msconvert_exists == FALSE
  # output_pparsemod_exists == FALSE
  # either second_search == FALSE & output_pglyco1_exists == FALSE or
  # or     second_search == TRUE  & output_pglyco2_exists == FALSE
  
  local({
    
    if(!dir.exists("msconvert_output"))  dir.create("msconvert_output")
    
    message("Run msconvert\n")
    
    run_msconvert <- function(file_name, out_dir){
      
      # Function runs msconvert through a system call
      # Arguments:
      # file_name = name of the file
      # out_idr = output directory
      # Value:
      # NULL
      
      wd <- getwd()
      
      # prepare system calls
      msConvert_args <- paste(paste0('"', file_name, '"'), '--outdir',
                              paste0('"', wd, '/msconvert_output', '"'),
                              '--mgf', '--ignoreUnknownInstrumentError',
                              '--singleThreaded', '--filter',
                              paste0('"', 'peakPicking vendor','"'), '--filter',
                              paste0('"', 'defaultArrayLength 1-', '"'),
                              '--filter',
                              paste0('"', 'titleMaker <RunId>.<ScanNumber>.<ScanNumber>.<ChargeState>', '"'),
                              collapse = " ")
      
      shell(cmd = paste("msconvert", msConvert_args, '& echo "done" >>',
                        paste0('"', out_dir, '/msconvert_done.txt', '"'),
                        collapse = " "), wait = FALSE, ignore.stdout = TRUE)
      
      
    }
    
    # file names to process
    files_to_process <- raw_file_names[contain_ms3]
    
    # number of files to process
    nr_files <- length(files_to_process)
    
    # maximum iterations given the number of available threads and number of files to process
    max_iterations <- ceiling(nr_files/nr_threads)
    if(is.na(max_iterations) || max_iterations < 1) stop("Wrong number of max_iterations")
    
    processes <- 0
    for(i in 1:max_iterations){
      
      processes <- (max(processes) + 1) : min(i*nr_threads, nr_files)
      
      if(verbose) message("msConvert processes\n", paste(files_to_process[processes], collapse = "\n"))
      
      for(j in processes){
        
        out_dir <- paste0(wd, "/msconvert_output/msconvert_process", j)
        dir.create(out_dir, showWarnings = verbose)
        
        run_msconvert(file_name = files_to_process[j], out_dir = out_dir)
        
      }
      
      while(any(!file.exists(paste0(wd,"/msconvert_output/msconvert_process", processes, "/msconvert_done.txt")))){
        
        #if(verbose) {
        
        cat(".")
        
        #} 
        
        Sys.sleep(2)
        
      }
      cat("\n")
      
    }
    
    # remove created folders
    list_dirs <- list.dirs(path = "msconvert_output", recursive = FALSE)
    to_delete <- list_dirs[grepl("msconvert_process[0-9]+", list_dirs)]
    unlink(to_delete, recursive = TRUE)
    
  })
  
} else {
  
  message("Skip MSConvert processing")
  
}

##### run pParse #####

if(!output_pparse_exists && !output_pparsemod_exists){  
  
  # if files pParse_mod.mgf already exist, there is no need to run pParse
  # pParse processing will be initiated only if neither of pParse files exist
  # output_pparse_exists == FALSE &
  # output_pparsemod_exists == FALSE
  
  local({
    
    if(!dir.exists("pparse_output"))  dir.create("pparse_output")
    
    message("\nrun pParse")
    
    run_pparse    <- function(file_name, out_dir, pparse_path){
      
      # Function runs pParse through a system call
      # Arguments:
      # file_name = name of the file
      # out_idr = output directory
      # pparse_path = path to pparse executable
      # Value:
      # NULL
      
      ### change working directory to pglyco location
      if(!dir.exists(pparse_path)) stop("Cannot find folder location of pparse")
      
      wd <- getwd()
      setwd(pparse_path)
      
      # prepare system calls
      cmd_pparse  <- paste('pParse.exe', '-D', paste0('"', normalizePath(wd), '\\', file_name, '"'), '-O', paste0('"', normalizePath(wd), '\\pparse_output', '"'), '-p', '0')
      
      shell(cmd = paste0(cmd_pparse, ' & echo "done" >> ', '"', out_dir, '/pparse_done.txt', '"'), wait = FALSE, ignore.stdout = TRUE)
      
      # switch back to the original directory
      setwd(wd)
      
    }
    
    # for pParse we need to find a path where the pParse is located and then set wd to its location
    # afterwards we can change the wd back to the location of the raw files
    
    pparse_path <- system2("where", args = "pparse", stdout = TRUE)
    pparse_path <- file.path(gsub("p[Pp]arse\\.exe", "", pparse_path))[[1]]
    
    if(!dir.exists(pparse_path)) stop("Cannot find location folder of pParse.exe
                                      Please check if the path is set
                                      in the environmental variables.")
    
    # number of files to process
    nr_files <- length(raw_file_names)
    
    # file names to process
    files_to_process <- raw_file_names
    
    # maximum iterations given the number of available threads and number of files to process
    max_iterations <- ceiling(nr_files/nr_threads)
    if(is.na(max_iterations) || max_iterations < 1) stop("Wrong number of max_iterations")
    
    processes <- 0
    for(i in 1:max_iterations){
      
      processes <- (max(processes) + 1) : min(i*nr_threads, nr_files)
      
      if(verbose) message("pParse processes\n",
                          paste(files_to_process[processes], collapse = "\n"))
      
      for(j in processes){
        
        out_dir <- paste0(wd, "/pparse_output/pparse_process", j)
        
        if(!dir.exists(out_dir)) dir.create(out_dir, showWarnings = verbose)
        
        run_pparse(file_name = files_to_process[j], out_dir = out_dir, pparse_path = pparse_path)
        
      }
      
      while(any(!file.exists(paste0(wd,"/pparse_output/pparse_process", processes, "/pparse_done.txt")))){
        
        #if(verbose) {
        
        cat(".")
        
        #} 
        
        Sys.sleep(2)
        
      }
      cat("\n")
      
    }
    
  })
  
  if(!report_intermediate_files){
    
    # remove created folders
    list_dirs <- list.dirs("pparse_output", recursive = FALSE)
    to_delete <- list_dirs[grepl("pparse_process[0-9]+", list_dirs)]
    unlink(to_delete, recursive = TRUE)
    
  }
  
  if(!report_intermediate_files){
    
    suppressWarnings(
      
      invisible(
        
        file.remove(paste0("pparse_output/",
                           c(gsub("\\.raw$", ".xtract", raw_file_names),      # pParse xtract file
                             gsub("\\.raw$", ".ms1", raw_file_names),         # pParse ms1 
                             gsub("\\.raw$", ".ms2", raw_file_names),         # pParse ms2  
                             gsub("\\.raw$", ".csv", raw_file_names)          # pparse csv  
                           )))
        
      )
    )
    
  }
  
} else {
  
  message("Found all pParse output files. Skip pParse processing")
  
}


# check pParse output

local({
    
  pparse_out <- list.files(path = "pparse_output", pattern = "_[A-Z]+FT\\.mgf")
  pparse_out <- gsub("_[A-Z]+FT\\.mgf", "", pparse_out)
  output_pparse_exists <<- all(gsub("\\.raw$", "", raw_file_names) %in%
                                 pparse_out)
  if(!output_pparse_exists &&
     !output_pparsemod_exists) stop("Not all pParse output files found.
                                    Check pParse log file")
  
  if(verbose) message("Found FT.mgf files\n",
                      paste(list.files(pattern = "_[A-Z]+FT\\.mgf"),
                            collapse = "\n"))
  
  # check RawTools output
  rawtools_out <- paste0("rawtools_output\\", raw_file_names, "_Matrix.txt")
  output_rawtools_exists  <<- all(file.exists(rawtools_out))
  if(!output_rawtools_exists) stop("Not all _Matrix.txt files found.
                                   Check RawTools log file.")
  
  if(verbose) message("Found _Matrix.txt files\n",
                      paste(list.files(pattern = "_Matrix\\.txt"),
                            collapse = "\n"))
  
  # check MSconvert output
  msconv_output <- gsub("\\.raw$", "\\.mgf",raw_file_names[contain_ms3])
  msconv_output <- paste0("msconvert_output/", msconv_output)
  output_msconvert_exists <<- all(file.exists(msconv_output))
  if(any(contain_ms3) && !output_msconvert_exists && !output_pparsemod_exists &&
     ((!second_search && !output_pglyco1_exists) ||
      (second_search  && !output_pglyco2_exists))) stop("Not all .mgf files found.
                                                        Check MSconvert parameters")
  
  if(verbose) message("Found .mgf files\n",
                      paste(list.files(pattern = "\\.mgf$"), collapse = "\n"))

})

##### merge MS2 and MS3 scans #####

# raw files without ms3 scans are renamed to "_pParse_mod.mgf" for consistency
if(any(!contain_ms3)) {
  
  local({
    
    suppressWarnings({ 
      
      # find pparse output files
      pparse_out_files <- list.files(path = "pparse_output",
                                     pattern = "_[A-Z]+FT\\.mgf")
      
      # pparse mgf files from raw files that do not contain MS3 scans
      pparse_out_files <- pparse_out_files[gsub("_[A-Z]+FT\\.mgf", "",
                                                pparse_out_files) %in%
                                           gsub("\\.raw", "",
                                                raw_file_names[!contain_ms3])] 
      
      # rename
      file.rename(paste0("pparse_output/", pparse_out_files),
                  paste0("pparse_output/", gsub("(.+)_[A-Z]+FT\\.mgf$",
                                                "\\1_pParse_mod.mgf",
                                                pparse_out_files)))
      
    })
    
  })
  
}

# names of pParse files
pparse_out_files <- list.files(path = "pparse_output",
                               pattern = "_[A-Z]+FT\\.mgf")
pparse_mod_files <- paste0("pparse_output/",
                           gsub("\\.raw$", "_pParse_mod.mgf",
                                raw_file_names[contain_ms3]))

if(any(contain_ms3) && !all(file.exists(pparse_mod_files))){
  
  # there are raw files that contain MS3 ions
  # Merged MS2 and MS3 spectra are not available
  # for all raw files containing MS3 spectra

  # initiate workers
  plan(strategy = plan_strategy, workers = nr_threads, gc = TRUE)
  
  local({
    
    ptm_merge_ms <- proc.time() 
    message("Merging MS2 and MS3 scans")
    
    # merge MS2 and MS3 spectra
    merge_ms2_ms3 <- function(matrix_data, mgf_tab, 
                              ion_match_tolerance,
                              tolerance_unit){
      
      # Function combines ms2 and ms3 ion intensities
      # Arguments:
      # matrix_data = _Matrix.txt file from RawTools loaded as a data table
      # mgf_tab = msconvert mgf file read in with readMgf function
      # ion_match_tolerance = tolerance window for matching ions
      # tolerance_unit = units of the tolerance window
      #
      # Value:
      # !Function returns NULL, but it modifies mgf_tab in place!
      
      invisible(
        
        lapply(1:matrix_data[, .N], function(row_nr){
          
          scan_ms2 <- which(mgf_tab$scan_number == matrix_data$MS2ScanNumber[row_nr])[1]
          scan_ms3 <- which(mgf_tab$scan_number == matrix_data$MS3ScanNumber[row_nr])[1]
          
          ms2 <- mgf_tab[scan_ms2, ions][[1]]
          ms3 <- mgf_tab[scan_ms3, ions][[1]]
          
          if(length(ms2) == 0 || length(ms3) == 0) return(NULL)
          
          # first do a rough match by integer ion mz
          not_matching <- which(!floor(ms3[, 1]) %in% c(floor(ms2[, 1]),
                                                        (floor(ms2[, 1]) - 1),
                                                        (floor(ms2[, 1]) + 1)))
          
          # combine spectra and stop if there are no matches
          if(length(not_matching) == length(ms3[, 1])){
            
            ms2 <- rbind(ms2, ms3)
            ms2 <- ms2[order(ms2[, 1]), ]
            set(mgf_tab, i = scan_ms2, j = "ions", value = list(list(ms2)))
            return(NULL)
            
          }
          
          # check remaining for overlap precisely
          # indices of ms3 ions that might have a match to ms2 ions
          probable_match <- setdiff(seq_along(ms3[, 1]), not_matching)
          # mz values of ms3 ions that might have a match to ms2 ions
          test_ms3_ions  <- ms3[probable_match, 1]
          # mz differences between ms3 and ms2 ions
          differences    <- abs(outer(test_ms3_ions, ms2[, 1], FUN = "-")) 
          
          # provide tolerance level for each ms3 ion
          if(tolerance_unit == "ppm"){
            
            tolerances <- ion_match_tolerance*ms3[probable_match, 1]/10^6
            
          } else if(tolerance_unit == "Th"){
            
            tolerances <- rep(ion_match_tolerance,
                              times = length(probable_match))
            
          } else {
            
            message("ion matching tolerance is 1Th")
            tolerances <- rep(1, times = length(probable_match))
            
          }
          
          # which ions passed the tolerance cutoff
          passed <- data.table(which(differences < tolerances |
                                     dplyr::near(differences, tolerances),
                                     arr.ind = TRUE))
          
          # combine spectra and stop if there are no ions
          # that passed tolerance cutoff
          if(passed[, .N] == 0) {
            
            ms2 <- rbind(ms2, ms3)
            ms2 <- ms2[order(ms2[, 1]), ]
            set(mgf_tab, i = scan_ms2, j = "ions", value = list(list(ms2)))
            return(NULL)
            
          }
          
          # if one ms3 ion (row) corresponds to several ms2 ions,
          # their indices will be saved as a list 
          passed <- passed[, .(ms2_ind = list(c(col))), by = row]
          
          # one ms3 ion (row) -> one ms2 ion (ms2)
          # decide which ions to merge
          # ms2 column is a list, ms2[[1]] is required
          
          if(all(lengths(passed$ms2_ind) == 1)){
            
            passed[, best_ms2 := unlist(ms2_ind)]
            
          } else {
            
            passed[, best_ms2 := ms2_ind[[1]][1], by = "row"]
            passed[lengths(ms2_ind) > 1,
                   best_ms2 := ms2_ind[which.min(
                     
                     abs(ms2[ms2_ind[[1]]] - test_ms3_ions[row])
                     
                     )],
                   by = "row"]
            
          }
          
          # one ms2 spectrum (row) -> one ms3 spectrum
          # decide which ions to merge
          # group by ms2 spectra
          # note: ms3 column is a list, ms3[[1]] is required
          passed2 <- passed[, .(ms3_ind = list(c(row))), by = best_ms2]
          
          if(all(lengths(passed2$ms3_ind) == 1)){
            
            passed2[, best_ms3 := unlist(ms3_ind)]
            
          } else {
            
            passed2[, best_ms3 := ms3_ind[[1]][1], by = "best_ms2"]
            passed2[lengths(ms3_ind) > 1,
                    best_ms3 := ms3_ind[which.min(
                      
                      abs(ms2[best_ms2] - test_ms3_ions[ms3_ind[[1]]])
                      
                      )], by = "best_ms2"]
            
          }
          
          
          # get an actual index in ms3 (ms3_ind was pointing to probable_match)
          passed2[, best_ms3 := probable_match[best_ms3], by = "best_ms2"]
          
          ms2[passed2[["best_ms2"]], 2] <- rowSums(
            
            matrix(c(ms2[passed2[["best_ms2"]], 2],
                     ms3[passed2[["best_ms3"]], 2]), 
                   byrow = FALSE, ncol = 2)
            
            )
          # remove matched ms3 ions
          ms3 <- ms3[-passed2[["best_ms3"]], ]
          
          # combine ms2 and not-matching ms3 ions
          ms2 <- rbind(ms2, ms3) 
          ms2 <- ms2[order(ms2[, 1]), ]
          
          # update mgf
          set(mgf_tab, i = scan_ms2, j = "ions", value = list(list(ms2)))
          return(NULL)
          
        })
      )
      
    }
    
    # use future_lapply function for parallelization
    future_lapply(raw_file_names[contain_ms3], FUN = function(file_name){
      
      file_name_core <- gsub("\\.raw", "", file_name)
      
      # read matrix file
      mtx <- fread(paste0("rawtools_output/", file_name, "_Matrix.txt"))
      
      # read mgf file (msconvert output)
      mgf_tab <- readMgf(paste0("msconvert_output/", file_name_core, ".mgf"))
      
      invisible(  
        merge_ms2_ms3(mgf_tab = mgf_tab,
                      matrix_data = mtx,
                      ion_match_tolerance = ion_match_tolerance,
                      tolerance_unit = tolerance_unit)
        
      )
      
      save(mgf_tab, file = paste0("pparse_output/",
                                  file_name_core, "_merged.Rdata"))
      
    })
    
    if(verbose) print(proc.time() - ptm_merge_ms)
    
  })
  
  local({
    
    ptm_merge_ms <- proc.time() 
    message("Update pParse mgf with merged spectra")
    
    # replace MS2 spectra in pParse by merged MS2/MS3 spectra
    change_ions_pParse <- function(mgf_tab, mgf_pParse, file_name){
      
      # function keeps scan headers in the mgf from pParse
      # it replaces the MS2 scans with merged MS2/MS3 scans
      # save file as a new pdf
      # Arguments:
      # mgf_tab = mgf file in form of data table with merged MS2/MS3 scans
      # (after merge_ms2_ms3 function)
      # mgf_pParse = path to pParse mgf file
      # file_name = file name for the output file
      # step = how many spectra to process at once, needed for optimization
      # Value
      # The function returns NULL; mgf file is saved under file_name
      
      if(file.exists(file_name)) suppressWarnings(file.remove(file_name)) 
      
      mgf_pParse <- readMgfHeader(mgf_pParse)
      mgf_pParse <- merge(mgf_pParse, mgf_tab[, c("scan_number", "ions")],
                          by = "scan_number", all.x = TRUE)
      
      if(nrow(mgf_pParse) < 1) stop()
      
      mgf_pParse[, scan_id := 1:.N]
      fwrite(list(unlist(
        
          mgf_pParse[, .(scan = list(
          
            c(unlist(scan_header),
              paste(ions[[1]][, 1], ions[[1]][, 2], sep = " "),
              "END IONS")
          
        )), by = scan_id]$scan
      
      )), file_name)
      
    }
 
    future_lapply(raw_file_names[contain_ms3], FUN = function(file_name){
    
      file_name_core <- gsub("\\.raw", "", file_name)
    
      # replace ion intensities in mgf_pParse by merged intensities mgf_tab
      pparse_files <- gsub("_[A-Z]+FT\\.mgf", "", pparse_out_files)
      pparse_out <- pparse_out_files[pparse_files %in% file_name_core]
      pparse_out <- paste0("pparse_output/", pparse_out)
      pparse_mod <- paste0("pparse_output/", file_name_core, "_pParse_mod.mgf")
      
      # load mgf_tab
      load(paste0("pparse_output/", file_name_core, "_merged.Rdata"))
      
      change_ions_pParse(mgf_tab = mgf_tab,
                         mgf_pParse = pparse_out,
                         file_name = pparse_mod)
      
    })
    
    if(verbose) print(proc.time() - ptm_merge_ms)
    
  })
  
  plan(strategy = "sequential")
  
} else {
  
  message("Found all pParse_mod.mgf files. Skip merging MS2 and MS3 spectra")
  
}

cat("\n")

##### run pGlyco #####

if(!output_pglyco1_exists &&                 # output from the first pGlyco search does not exist
   !(second_search && output_pglyco2_exists) # no need to process if second search is needed and
   # the output from the second pGlyco search already exists
){ 
  
  local({
    
    ptm_pglyco <- proc.time() 
    
    # create output folder
    suppressWarnings(
      
      dir.create("pglyco_output\\")
      
    )
    
    message("\n Run pGlyco ")
    
    # find location of pglyco executable
    
    pglyco_path <- system2("where", args = "pglyco", stdout = TRUE)
    pglyco_path <- file.path(gsub("pGlyco\\.exe", "", pglyco_path))[[1]]
    
    if(!dir.exists(pglyco_path)){ # another way to find pGlyco
      
      pglyco_path <- shell("set path", intern = TRUE)
      pglyco_path <- unlist(strsplit(pglyco_path, ";"))
      pglyco_path <- pglyco_path[grepl("p[Gg]lyco", pglyco_path)][[1]]
      
    }
    
    if(!dir.exists(pglyco_path)) stop("Cannot find location folder of pGlyco.exe. Please check if the path is set in the environmental variables.")
    
    # check if configuration file already exists. If not, modify the template accordingly
    
    if(any(grepl("^pGlyco_task.*\\.pglyco$", list.files()))){ #check if configuration file already exists
      
      pglyco_config_file <- list.files(pattern = "^[Pp][Gg]lyco_task\\.[Pp][Gg][Ll][Yy][Cc][Oo]$")
      pglyco_config <- readLines(pglyco_config_file)
      
      #fasta
      pglyco_config[grepl("^fasta=", pglyco_config)]      <- paste0("fasta=", wd,  "/", fasta.name, ".N2J")
      
      writeLines(pglyco_config, "pGlyco_task.pglyco")
      
      
    } else {
      
      message("pGlyco configuration was not provided. Use default settings")
      
      pglyco_config <- pglyco_config_default
      
      ### modify configuration file
      
      # pglyco version
      setwd(pglyco_path)
      pglyco_version <- system2("pglycodb", args = "--version", stdout = TRUE)
      pglyco_config[grepl("^pGlyco_version", pglyco_config)] <- paste0("pGlyco_version=", pglyco_version)
      setwd(wd)
      
      # number of processes
      pglyco_config[grepl("^process=", pglyco_config)]    <- paste0("process=", nr_threads)
      
      #output_dir
      pglyco_config[grepl("^output_dir=", pglyco_config)] <- paste0("output_dir=", normalizePath(paste0(wd, "/pglyco_output"), winslash = "\\"))
      
      #fasta
      pglyco_config[grepl("^fasta=", pglyco_config)]      <- paste0("fasta=", wd,  "/", fasta.name, ".N2J")
      
      #fixed modifications
      
      # fixed modification should match the reporter_ion_type parameter for RawTools
      if(reporter_ion != "not_labeled"){
        
        if(reporter_ion == "TMT0"){
          
          fixed_mod_reporter <- "TMT"
          
        } else {
          
          fixed_mod_reporter <- paste0(reporter_ion, "plex")
          
        }
        
        pglyco_config[grepl("^fix2=", pglyco_config)]       <- paste0("fix2=", fixed_mod_reporter, "[K]")
        pglyco_config[grepl("^fix3=", pglyco_config)]       <- paste0("fix3=", fixed_mod_reporter, "[AnyN-term]")
        
      }
      # total number of files to analyze
      
      pglyco_config[grepl("^spectrum_total=", pglyco_config)] <- paste0("spectrum_total=", length(raw_file_names))
      
      # add paths to files to analyze
      
      pglyco_config <- pglyco_config[1:which(grepl("^spectrum_total=", pglyco_config))] # should be the last line before file paths
      for(i in seq_along(raw_file_names)){
        
        pglyco_config <- c(pglyco_config, paste0("file", i, "=", wd, "/pparse_output/", gsub("\\.raw$", "_pParse_mod.mgf", raw_file_names[i])))  
        
      }
      
      writeLines(pglyco_config, "pGlyco_task.pglyco")
      
    }
    
    # path to pglyco config file
    pglyco_config_file <- paste0(wd, "/", list.files(pattern = "pGlyco_task\\.pglyco")[1])
    
    # number of files to process
    nr_files <- length(raw_file_names)
    
    # file names to process
    files_to_process <- paste0("pparse_output/", gsub("\\.raw$", "_pParse_mod.mgf", raw_file_names))
    
    # maximum iterations given the number of available threads and number of files to process
    max_iterations <- ceiling(nr_files/nr_threads)
    if(is.na(max_iterations) || max_iterations < 1) stop("Wrong number of max_iterations")
    
    pglyco_config <- readLines(pglyco_config_file)
    
    processes <- 0
    for(i in 1:max_iterations){
      
      processes <- (max(processes) + 1) : min(i*nr_threads, nr_files)
      
      if(verbose) message("pGlyco processes\n", paste(files_to_process[processes], collapse = "\n"))
      
      for(j in processes){
        
        out_dir <- paste0(wd, "/pglyco_output/pglyco_process", j)
        if(!dir.exists(out_dir)) dir.create(out_dir, showWarnings = verbose)
        
        run_pglyco(file_name = files_to_process[j], out_dir = out_dir, pglyco_config = pglyco_config, pglyco_path = pglyco_path)
        
      }
      
      while(any(!file.exists(paste0(wd,"/pglyco_output/pglyco_process", processes, "/done.txt")))){
        
        #if(verbose) {
        
        cat(".")
        
        #} 
        
        Sys.sleep(2)
        
      }
      cat("\n")
      
    }
    
    df_pglyco <- data.table()
    for(i in seq_along(raw_file_names)){
      
      if(!file.exists(paste0("pglyco_output/pglyco_process", i, "/pGlycoDB-GP-FDR-Pro.txt"))){
        
        message("No pGlyco result for ", raw_file_names[i])
        next
        
      }
      
      temp <- fread(paste0("pglyco_output/pglyco_process", i, "/pGlycoDB-GP-FDR-Pro.txt"))
      df_pglyco <- rbind(df_pglyco, temp)
      
      #remove the label "done"
      file.remove(paste0("pglyco_output/pglyco_process", i, "/done.txt"))
      
    }
    
    # check the length of df_pglyco
    if(nrow(df_pglyco) == 0) stop("No pGlyco results")
    
    fwrite(df_pglyco, "pglyco_output/pGlycoDB-GP-FDR-Pro.txt", sep = "\t", quote = FALSE, row.names = FALSE)
    
    if(!report_intermediate_files && file.exists("pGlycoDB-GP-FDR-Pro.txt")){
      
      list_dirs <- list.dirs(path = "pglyco_output/", recursive = FALSE)
      to_delete <- list_dirs[grepl("pglyco_process[0-9]+", list_dirs)]
      unlink(to_delete, recursive = TRUE)
      
      
    }
    
    
    if(verbose) proc.time() - ptm_pglyco
    
  })
  
  
} else {
  
  message("Found -Pro.txt file. Skip pGlyco processing")
  
}

##### second pGlyco search #####

# check that second pglyco output is present
if(second_search){
  
  if(!output_pglyco2_exists){ # check if pGlyco result already exists
    
    local({
      
      ptm_pglyco <- proc.time()  
      
      # find pGlyco result file from the first search
      pglyco_out   <- list.files(path = "pglyco_output", pattern = "-Pro.txt$", full.names = TRUE)
      
      if(length(pglyco_out) == 0) stop(paste0("Cannot find pGlyco output files in the specified directory"))
      if(verbose) message(paste0("Found pGlyco files:\n", paste(pglyco_out, collapse = "\n")))
      if(length(pglyco_out) > 1 && verbose) message("Found several pGlyco output files - only the first one will be used ", pglyco_out[1])
      pglyco_out <- pglyco_out[1]
      
      pglyco_file <- fread(pglyco_out, header = T)
      pglyco_file <- pglyco_file[TotalFDR < (pglyco_fdr_threshold)]
      proteins    <- unique(unlist(stringr::str_split(pglyco_file[["Proteins"]], "/")))
      proteins    <- proteins[!grepl("^REV_", proteins)]
      
      if(!file.exists(gsub("(.+)\\.[Ff][Aa][Ss][Tt][Aa]$", "\\1_sub.fasta.N2J", fasta.name))){
        
        ff_n2j <- fread(paste0(fasta.name, ".N2J"), sep = NULL, header = FALSE)[[1]]
        take   <- ff_headers %in% proteins
        start  <- ff_header_pos[which(take)]
        end    <- ff_header_pos[which(take) + 1] - 1 #line where starts the next header - 1
        end[is.na(end)] <- length(ff_n2j)
        
        writeLines(unlist(Map(function(start, end) ff_n2j[start:end], start = start, end = end)), gsub("(.+)\\.fasta$", "\\1_sub.fasta.N2J", fasta.name))
        
      }
      
      message("\n Run second pGlyco search ")
      
      # find location of pglyco
      
      pglyco_path <- system2("where", args = "pglyco", stdout = TRUE)
      pglyco_path <- file.path(gsub("pGlyco\\.exe", "", pglyco_path))[[1]]
      
      if(!dir.exists(pglyco_path)){ # another way to find pGlyco
        
        pglyco_path <- shell("set path", intern = TRUE)
        pglyco_path <- unlist(strsplit(pglyco_path, ";"))
        pglyco_path <- pglyco_path[grepl("p[Gg]lyco", pglyco_path)][[1]]
        
      }
      
      if(!dir.exists(pglyco_path)) stop("Cannot find location folder of pGlyco.exe. Please check if the path is set in the environmental variables.")               
      
      # check if configuration file already exists. If not, modify the template accordingly
      
      if(any(grepl("^pGlyco_task.*\\.pglyco$", list.files()))){ #check if configuration file already exists
        
        pglyco_config_file <- list.files(pattern = "^[Pp][Gg]lyco_task\\.[Pp][Gg][Ll][Yy][Cc][Oo]$")
        
        pglyco_config <- readLines(pglyco_config_file)
        
        # change fasta file name to the name of the reduced fasta file
        pglyco_config[grepl("^fasta=", pglyco_config)]      <- paste0("fasta=", wd,  "/", gsub("(.+)\\.fasta$", "\\1_sub.fasta.N2J", fasta.name))
        
        writeLines(pglyco_config, "pGlyco_task.pglyco")
        
      } else {
        
        if(verbose) message("pGlyco configuration was not provided. Use default settings")
        
        pglyco_config <- pglyco_config_default
        
        ### modify configuration file
        
        # pglyco version
        setwd(pglyco_path)
        pglyco_version <- system2("pglycodb", args = "--version", stdout = TRUE)
        pglyco_config[grepl("^pGlyco_version", pglyco_config)] <- paste0("pGlyco_version=", pglyco_version)
        setwd(wd)
        
        # number of processes
        pglyco_config[grepl("^process=", pglyco_config)]   <- paste0("process=", nr_threads)
        
        #output_dir
        pglyco_config[grepl("^output_dir=", pglyco_config)] <- paste0("output_dir=", normalizePath(wd, winslash = "\\"))
        
        #fasta
        pglyco_config[grepl("^fasta=", pglyco_config)]      <- paste0("fasta=", wd,  "/", gsub("(.+)\\.fasta$", "\\1_sub.fasta.N2J", fasta.name))
        
        #fixed modifications
        
        # fixed modification should match the reporter_ion_type parameter for RawTools
        if(reporter_ion != "not_labeled"){
          
          if(reporter_ion == "TMT0"){
            
            fixed_mod_reporter <- "TMT"
            
          } else {
            
            fixed_mod_reporter <- paste0(reporter_ion, "plex")
            
          }
          
          pglyco_config[grepl("^fix2=", pglyco_config)]       <- paste0("fix2=", fixed_mod_reporter, "[K]")
          pglyco_config[grepl("^fix3=", pglyco_config)]       <- paste0("fix3=", fixed_mod_reporter, "[AnyN-term]")
          
        }
        
        # total number of files to analyze
        
        pglyco_config[grepl("^spectrum_total=", pglyco_config)] <- paste0("spectrum_total=", length(raw_file_names))
        
        # add paths to files to analyze
        
        pglyco_config <- pglyco_config[1:which(grepl("^spectrum_total=", pglyco_config))] # should be the last line before file paths
        for(i in seq_along(raw_file_names)){
          
          pglyco_config <- c(pglyco_config, paste0("file", i, "=", wd, "/", gsub("\\.raw$", "_pParse_mod.mgf", raw_file_names[i])))  
          
        }
        
        writeLines(pglyco_config, "pGlyco_task.pglyco")
        
      }
      
      # path to pglyco config file
      pglyco_config_file <- paste0(wd, "/", list.files(pattern = "pGlyco_task.*\\.pglyco")[1])
      
      # number of files to process
      nr_files <- length(raw_file_names)
      
      # file names to process
      files_to_process <- paste0("pparse_output/", gsub("\\.raw$", "_pParse_mod.mgf", raw_file_names))
      
      # maximum iterations given the number of available threads and number of files to process
      max_iterations <- ceiling(nr_files/nr_threads)
      if(is.na(max_iterations) || max_iterations < 1) stop("Wrong number of max_iterations")
      
      pglyco_config <- readLines(pglyco_config_file)
      
      processes <- 0
      for(i in 1:max_iterations){
        
        processes <- (max(processes) + 1) : min(i*nr_threads, nr_files)
        
        if(verbose) message("pGlyco processes\n", paste(files_to_process[processes], collapse = "\n"))
        
        for(j in processes){
          
          out_dir <- paste0(wd, "/pglyco_output/pglyco_process", j)
          dir.create(out_dir, showWarnings = verbose)
          
          run_pglyco(file_name = files_to_process[j], out_dir = out_dir, pglyco_config = pglyco_config, pglyco_path = pglyco_path)
          
        }
        
        while(any(!file.exists(paste0(wd,"/pglyco_output/pglyco_process", processes, "/done.txt")))){
          
          #if(verbose) {
          
          cat(".")
          
          #} 
          
          Sys.sleep(2)
          
        }
        cat("\n")
        
      }
      
      df_pglyco <- data.table()
      for(i in seq_along(raw_file_names)){
        
        if(!file.exists(paste0("pglyco_output/pglyco_process", i, "/pGlycoDB-GP-FDR-Pro.txt"))){
          
          message("No pGlyco result for ", raw_file_names[i])
          next
          
        }
        
        temp <- fread(paste0("pglyco_output/pglyco_process", i, "/pGlycoDB-GP-FDR-Pro.txt"))
        df_pglyco <- rbind(df_pglyco, temp)
        
        #remove the label "done"
        file.remove(paste0("pglyco_output/pglyco_process", i, "/done.txt"))
        
      }
      
      # check the length of df_pglyco
      if(nrow(df_pglyco) == 0) stop("No pGlyco results")
      
      # write combined pGlyco results
      fwrite(df_pglyco, "pglyco_output/pGlycoDB-GP-FDR-Pro2.txt", sep = "\t", quote = FALSE, row.names = FALSE)
      
      if(!report_intermediate_files && file.exists("pglyco_output/pGlycoDB-GP-FDR-Pro2.txt")){
        
        list_dirs <- list.dirs("pglyco_output", recursive = FALSE)
        to_delete <- list_dirs[grepl("pglyco_process[0-9]+", list_dirs)]
        unlink(to_delete, recursive = TRUE)
        
        
      }
      
      
      if(verbose) proc.time() - ptm_pglyco
      
      
      
    })
    
  } else {
    
    message("Found -Pro2.txt file. Skip pGlyco processing")
    
  }  
  
}  

##### combine pglyco output and RawTools output #####
local({
  
  # find pGlyco output file
  if(second_search) pattern <- "-Pro2.txt$" else pattern <- "-Pro.txt$"
  pglyco_out <- list.files(path = "pglyco_output", pattern = pattern, full.names = TRUE)
  
  if(length(pglyco_out) == 0) stop(paste0("Cannot find pGlyco output files"))
  if(verbose) message(paste0("Found pGlyco files:\n", paste(pglyco_out, collapse = "\n")))
  
  # check that matrix files are present
  matrix_files  <- list.files(path = "rawtools_output", pattern = "_Matrix.txt", full.names = FALSE)
  if(length(matrix_files) == 0) stop(paste0("Cannot find RawTools _Matrix files in the specified directory"))
  if(verbose) message(paste0("Found _Matrix.txt files:\n", paste(matrix_files, collapse = "\n")))
  if(verbose) message("Combining pglyco and RawTools output")
  
  addRawtoolsData <- function(search_files,
                              rawtools_files){
    
    # read the pglyco output file
    search_data <- data.table()
    for(i in seq_along(pglyco_out)){
      
      search_data <- rbind(search_data,
                           fread(search_files[i], header = T, na.strings = "NA",
                                 stringsAsFactors = FALSE, key = "Scan"))  
      
    }
    
    # format the raw file names
    search_data[, RawName := gsub("\\..+$", "", search_data$GlySpec)]
    
    # read rawtools output
    matrix_names <- stringr::str_split(rawtools_files, "\\.", simplify = TRUE)[, 1]
    rawtools_data <- data.table()
    for(i in seq_along(rawtools_files)){
      
      temp <- fread(paste0("rawtools_output/", rawtools_files[i]),
                        stringsAsFactors = FALSE, header = TRUE, key = "MS2ScanNumber")
      #remove an empty column if it is appended
      if(sum(is.na(temp[, c(ncol(..temp))])) == nrow(temp)) temp <- temp[, -ncol(temp), with = FALSE]
      temp[, RawName := matrix_names[i]]
      rawtools_data <- rbind(rawtools_data, temp)
      
    }
    search_data <- merge(search_data, rawtools_data,
                         by.x = c("RawName", "Scan"), by.y = c("RawName", "MS2ScanNumber"))
    
    return(search_data)
    
  }
  
  search_data <- addRawtoolsData(search_files = pglyco_out, rawtools_files = matrix_files)

  fwrite(search_data, "pglyco_output/pglyco_quant_results.txt",
         na = "NA", row.names = FALSE, quote = FALSE, sep = "\t")
  
})

##### identify oxonium ions #####
if(!skip_marker_ions){
  
  local({
    
    message("\nidentify marker ions")
    
    if(!dir.exists("pglyco_output")) dir.create("pglyco_output")
    
    if(!file.exists("marker_ions.txt")) {
      
      marker_ions <- fread(marker_ions, sep = ";")
      fwrite(marker_ions, "marker_ions.txt", sep = "\t")
      
    } else {
      
      marker_ions_temp <- fread("marker_ions.txt")
      
      marker_ions_temp[, Mass := as.numeric(Mass)]
      marker_ions_temp <- marker_ions_temp[!is.na(Mass)]
      marker_ions_temp <- marker_ions_temp[!duplicated(Mass)]
      
      
      if(nrow(marker_ions_temp) == 0){
        
        message("A problem has occured with supplied oxonium ion table. Use default oxonium ions")
        marker_ions <- fread(marker_ions, sep = ";")
        
      } else {
        
        marker_ions <- marker_ions_temp
        
      }
      
    }
    
    plan(strategy = plan_strategy, workers = nr_threads, gc = TRUE)
    pglyco_output <- fread("pglyco_output/pglyco_quant_results.txt")
    
    # use future_lapply function for parallelization
    list_dt <- future_lapply(raw_file_names, FUN = function(file_name){
      
      # read mgf file (msconvert output)
      if(file.exists(paste0("pparse_output/", gsub(".raw", "_merged.Rdata", file_name)))){
        
        load(paste0("pparse_output/", gsub(".raw", "_merged.Rdata", file_name)))
        
      } else {
        
        mgf_tab <- readMgf(paste0("pparse_output/", gsub(".raw", "_pParse_mod.mgf", file_name)))
        # remove duplicated scan numbers (pParse assigns multiple precursor masses per scan)
        mgf_tab <- mgf_tab[!duplicated(mgf_tab$scan_number)] 
        
      }
      
      # keep only identified scans 
      mgf_tab <- mgf_tab[scan_number %in% pglyco_output[RawName ==  gsub(".raw$", "", file_name)]$Scan]
      
      ion_match(marker_ions, mgf_tab, 0.01, "Th")
      
    })
    
    list_dt <- Map(function(dt, RawName) return(dt[, RawName := RawName]),
                   dt = list_dt, RawName = gsub(".raw", "", raw_file_names))
    marker_dt <- rbindlist(list_dt)
    marker_dt <- marker_dt[!duplicated(marker_dt[, c("Scan", "RawName")])]
    
    # add marker ion intensities to pglyco output
    pglyco_output <- merge(pglyco_output, marker_dt, by = c("Scan", "RawName"), all.x = TRUE)
    fwrite(marker_dt, "pglyco_output/marker_ions_identified.txt", sep = "\t")
    fwrite(pglyco_output, "pglyco_output/pglyco_quant_results.txt", sep = "\t")
    
    plan(strategy = "sequential")
    
  })
  
}

##### define Glycotypes and glycotopes #####

if(!file.exists("pglyco_output/pglyco_quant_results.txt")) stop("Cannot find combined result file from pGlyco and RawTools output")
local({
  
  df <- fread("pglyco_output/pglyco_quant_results.txt", sep = "\t")

  fwrite(
    
    glycotopeConflict(
    
      defineGlycoTope(
      
        defineGlycoType(df)
      
      )
    
    ), "pglyco_output/pglyco_quant_results.txt", sep = "\t")
  
})

##### Combine reporter ion intensities of each Glycoform on a particular site #####

# Combine reporter ion intensities based on Peptide + GlyID
if(!file.exists("pglyco_output/pglyco_quant_results.txt")) stop("Cannot find combined result file from pGlyco and RawTools output")
if(verbose) message("Transforming output tables")

local({
  
  ptm_transform_out <- proc.time()
  
  df_scans <- fread("pglyco_output/pglyco_quant_results.txt", sep = "\t")
  df_scans <- df_scans[order(Peptide)]
  df_scans <- df_scans[TotalFDR < pglyco_fdr_threshold]
  df_scans[, id := 1:df_scans[, .N]]
  df_scans[, Unique_site := !grepl("/", Proteins)] # site is unique if only one protein id is reported
  
  # intensity names
  int_names <- grep("\\d\\d+[NC]?Intensity$", names(df_scans), value = TRUE)
  
  # marker ion names
  if(!skip_marker_ions){
    
    mion_names <- unlist(fread("pglyco_output/marker_ions_identified.txt", nrow = 1, header = FALSE))
    mion_names <- mion_names[!mion_names %in% c("Scan", "RawName")]
    
  } else {
    
    mion_names <- character()
    
  }
  # combine intensities for modified peptides
  modpept_int_sum        <- df_scans[, lapply(.SD, sum), by = .(Peptide, GlySite, `Glycan(H,N,A,G,F)`),
                                     .SDcols = c(int_names, mion_names)]
  modpept_precursor_data <- df_scans[, lapply(.SD, paste, collapse = pglyco_separator), by = .(Peptide, GlySite, `Glycan(H,N,A,G,F)`),
                                     .SDcols = c("id", "RawName", "Scan", "PrecursorMZ", "Charge", "Mod", "ParentPeakArea")]
  modpept_pept_data      <- df_scans[, head(.SD, 1), by = .(Peptide, GlySite, `Glycan(H,N,A,G,F)`),
                                     .SDcols = c("Proteins", "ProSite", "Unique_site")]
  modpept_glycan_data    <- df_scans[, head(.SD, 1), by = .(Peptide, GlySite, `Glycan(H,N,A,G,F)`),
                                     .SDcols = c("GlyID", "PlausibleStruct", "GlycanType", "GlycanAntennaType", "GlyFrag", "GlyMass")]
  modpept_glycan_data2   <- df_scans[, lapply(.SD, paste, collapse = pglyco_separator), by = .(Peptide, GlySite, `Glycan(H,N,A,G,F)`),
                                     .SDcols = c("Glycotope_F", "Glycotope_A", "Glycotope_F_unmatched", "Glycotope_A_unmatched")]
  
  df_modpept <- cbind(modpept_precursor_data,
                      modpept_pept_data[, -c("Peptide", "GlySite", "Glycan(H,N,A,G,F)")],
                      modpept_glycan_data[, -c("Peptide", "GlySite", "Glycan(H,N,A,G,F)")],
                      modpept_glycan_data2[, -c("Peptide", "GlySite", "Glycan(H,N,A,G,F)")],
                      modpept_int_sum    [, -c("Peptide", "GlySite", "Glycan(H,N,A,G,F)")])
  
  # pGlyco_ids correspond to the row number of df_scans, df_modpept gets its own ids
  names(df_modpept)[names(df_modpept) == "id"] <- "pGlyco_ids" 
  df_modpept[, id := 1:df_modpept[, .N]]
  
  # Proteins and sequence windows
  
  # one protein - site per row
  df_prot <- df_modpept[, list(Protein_single = unlist(strsplit(Proteins, pglyco_separator)),
                               ProSite_single = unlist(strsplit(as.character(ProSite),  pglyco_separator))),
                        by = id]
  
  # remove reversed sequences
  df_prot <- df_prot[!grepl(reversed_regexpr, df_prot$Protein_single), ] 
  
  # find header position in the fasta file
  df_prot[, ff_header_pos :=  ..ff_header_pos[match(df_prot$Protein_single, ..ff_headers)]]
  df_prot[, ff_seq_stop   := (..ff_header_pos[match(df_prot$Protein_single, ..ff_headers) + 1] - 1)]
  df_prot[is.na(ff_seq_stop), ff_seq_stop := length(fasta)]
  
  if(verbose && sum(is.na(df_prot$ff_header_pos)) > 0) message("\nProteins not found in fasta file: ", sum(is.na(df_prot$ff_header_pos)))
  
  # extract sequence windows
  df_prot[, seq_wind := unlist(lapply(seq_along(df_prot$ff_header_pos), function(i){
    
    temp_header_pos <- df_prot$ff_header_pos[i]
    temp_seq_stop   <- df_prot$ff_seq_stop[i]
    temp_seq_pos    <- as.integer(df_prot$ProSite_single[i])
    
    temp_sequence <- paste0(fasta[(temp_header_pos + 1) : temp_seq_stop], collapse = "")
    temp_seq_wind <- substr(temp_sequence, start = temp_seq_pos - seq_wind_size, stop = temp_seq_pos + seq_wind_size)
    
    # if the site is located at the end of the sequence: append "___"
    
    missing_start <- temp_seq_pos - (seq_wind_size + 1)
    missing_end   <- nchar(temp_sequence) - (temp_seq_pos + seq_wind_size)
    
    if(missing_start < 0) temp_seq_wind <- paste0(paste0(rep("_", times = abs(missing_start)), collapse = ""), temp_seq_wind)
    if(missing_end < 0)   temp_seq_wind <- paste0(temp_seq_wind, paste0(rep("_", times = abs(missing_end)),   collapse = ""))
    
    temp_seq_wind
    
  }))]
  
  # sequence windows without space holders
  df_prot[, seq_wind_noSpace := gsub("_", "", seq_wind)]
  
  # merge with df_modpept
  df_prot <- merge(df_prot[, -c("ff_header_pos", "ff_seq_stop")], df_modpept, by = "id", all = TRUE)
  
  # Protein ranking
  
  # Extract information about proteins in order to rank them
  # one protein per row
  df_prot_single <- df_prot[, list(Nr_Unique     = sum(Unique_site[!duplicated(ProSite_single)]),
                                   Nr_Sites      = length(unique(ProSite_single)),
                                   Nr_Glycoforms = length(id),
                                   Sites         = paste(unique(ProSite_single), collapse = "/"),
                                   Ids           = paste(unique(id), collapse = "/")),
                            by = Protein_single]
  
  # From SwissProt?
  df_prot_single[, SwissProt := grepl("^sp\\|", df_prot_single$Protein_single)]
  
  # is protein an Isoform?
  df_prot_single[, Protein_Isoform := grepl("[A-Za-z0-9]-\\d", df_prot_single$Protein_single)]
  
  # sort proteins = rank
  df_prot_single <- df_prot_single[order(-Nr_Unique, -Nr_Sites, -Nr_Glycoforms, -SwissProt, Protein_Isoform, Protein_single)]
  
  # assign Protein_id that can be also used to rank proteins
  df_prot_single[, Protein_id := 1:df_prot_single[, .N]]
  
  # add protein rank to df_prot table
  df_prot <- merge(df_prot, df_prot_single[, c("Protein_single", "Protein_id")], by = "Protein_single", all.x = TRUE)
  
  # Select sequence windows that explain modified peptdies
  
  # aggregate modpept_id for each sequence window
  sw_modpept <- df_prot[, lapply(.SD, function(x) paste0(unique(x), collapse = ";")),
                        by = .(seq_wind), .SDcols = c("id")]
  
  # aggregate Protein ids for each sequence window
  sw_protein <- df_prot[, lapply(.SD, function(x) paste0(x, collapse = ";")),
                        by = .(seq_wind), .SDcols = c("Protein_single", "ProSite_single")]
  
  # extract maximal (minimal number) protein rank for a sequence window
  sw_max_rank <- df_prot[, lapply(.SD, min), by = .(seq_wind), .SDcols = "Protein_id"]
  
  # df_sw table combines sequence windows and modified peptides that explain those sequence windows
  ids      <- unique(df_prot$id)
  list_ids <- strsplit(sw_modpept$id, split = ";")
  list_ids <- lapply(list_ids, as.numeric)
  
  df_sw <- data.table()
  while(length(ids) > 0){
    
    # find all peptides that form a group
    aggregated_sw  <- aggregateIds(list_ids, start_id = ids[1])
    sum(aggregated_sw)
    
    # extract modpept ids belonging to the same group
    unique_ids     <- unique(unlist(list_ids[aggregated_sw]))
    
    # how many peptides from the group fit into each sequence window
    ids_per_sw     <- lengths(list_ids[aggregated_sw])
    
    # order sequence windows based on 1. the number of peptides in the group it can accomodate and then if equal
    # 2. the maximal rank of the proteins it belongs to
    sw_ordered     <- sw_modpept$seq_wind[aggregated_sw][order(-ids_per_sw, sw_max_rank$Protein_id[aggregated_sw])]
    id_ordered     <- list_ids[aggregated_sw][order(-ids_per_sw, sw_max_rank$Protein_id[aggregated_sw])]
    
    # check, which sequence window can be attributed to the same peptides
    same <- lapply(id_ordered, function(x) {
      
      which(unlist(lapply(id_ordered, function(y) all.equal(x, y))) == TRUE)
      
    })
    
    # combine sequence windows with the similar peptide sets
    sw_list <-lapply(seq_along(same), function(i){
      
      sw_ordered[same[[i]]]
      
    })
    
    # remove duplicated peptide sets
    id_ordered <- id_ordered[!duplicated(sw_list)]
    sw_list    <- sw_list[!duplicated(sw_list)]
    
    # if there are more than one sequence window per peptide group, 
    # attribute shared peptides to the top sequence window (accomodates the majority of peptides)
    if(length(id_ordered) > 1){
      
      for(i in 2:length(id_ordered)){
        
        id_ordered[[i]] <- id_ordered[[i]][!id_ordered[[i]] %in% unlist(id_ordered[1:(i-1)])]
        
      }
      
    }
    
    # prepare a data table
    temp_df <- rbindlist(
      
      lapply(seq_along(id_ordered), function(i){
        
        proteins <- unlist(stringr::str_split(sw_protein$Protein_single[sw_protein$seq_wind %in% sw_list[[i]]], ";"))
        sites    <- unlist(stringr::str_split(sw_protein$ProSite_single[sw_protein$seq_wind %in% sw_list[[i]]], ";"))
        
        sites          <- paste0(sites[!duplicated(proteins)],    collapse = ";")
        proteins       <- paste0(proteins[!duplicated(proteins)], collapse = ";")
        
        data.table(seq_wind = paste(sw_list[[i]], collapse = ";"),
                   Proteins = proteins,
                   ProSites = sites,
                   modpept_ids = paste(id_ordered[[i]], collapse = ";"))
        
        
        
      })
      
    )
    
    # extend previous data table
    df_sw          <- rbind(df_sw, temp_df)
    
    # remove ids that were used from the list
    ids            <- ids[!ids %in% unique_ids] 
    
  }
  
  # convert factors to characters
  df_sw[] <- lapply(df_sw, as.character)
  
  # reorder proteins based on their rank
  temp <- lapply(seq_along(df_sw$seq_wind), function(i){
    
    temp     <- df_sw[i, ]
    temp[, Leading_Protein := NA]
    temp[, Leading_ProSite := NA]
    
    proteins <- unlist(strsplit(temp$Proteins, split = ";"))
    sites    <- unlist(strsplit(temp$ProSites, split = ";"))
    
    take     <- match(proteins, df_prot_single$Protein_single)
    ranks    <- df_prot_single$Protein_id[take]
    
    proteins <- proteins[order(ranks)]
    sites    <- sites[order(ranks)]
    
    sites    <- sites[!is.na(proteins)]
    proteins <- proteins[!is.na(proteins)]
    
    temp[, Leading_Protein := ..proteins[1]]
    temp[, Leading_ProSite := ..sites[1]]
    
    temp[, Proteins := paste0(..proteins, collapse = ";")]
    temp[, ProSites := paste0(..sites,    collapse = ";")]
    
    return(temp)
    
  })
  
  df_sw <- rbindlist(temp)
  
  # add sequence window information and re-ordered proteins to df_modpept
  df_modpept2 <- df_sw[, list(modpept_id = unlist(strsplit(modpept_ids, split = ";"))), by = seq_wind]
  df_modpept2 <- merge(df_modpept2, df_sw, by = "seq_wind") 
  df_modpept2[, modpept_id := as.integer(df_modpept2$modpept_id)]
  
  df_modpept  <- merge(df_modpept[, -c("Proteins", "ProSite")], df_modpept2[, -c("modpept_ids")], by.x = "id", by.y = "modpept_id", all.x = TRUE)
  
  
  # combine intensities for glycoform
  # glycoform = distinct sequence window + specific Glycan structure
  # distinct from df_modpept, because peptides containing missed cleavage sites will be combined 
  
  df_glycof_data  <- df_modpept[, lapply(.SD, paste, collapse = ";"), by = .(seq_wind, `Glycan(H,N,A,G,F)`),
                                .SDcols = c("id", "Scan", "pGlyco_ids", "Peptide", "GlySite", "GlyID")]
  # df_glycof_data2 <- df_modpept[, lapply(.SD, function(x) {
  #   
  #   glycotopes <- unlist(stringr::str_split(x, "[/;&]"))
  #   glycotopes <- gsub(" ", "", glycotopes)
  #   glycotopes <- glycotopes[!is.na(glycotopes) & glycotopes != ""]
  #   glycotopes <- unique(glycotopes)
  #   return(paste(glycotopes, collapse = ";"))
  #   
  # }), by = .(seq_wind, `Glycan(H,N,A,G,F)`), .SDcols = c("Glycotope_F", "Glycotope_A")]
  df_glycof_data2 <- df_modpept[, .(Glycotope = unlist(Map(function(gf, ga, unf, una){
    
    gf <- unlist(stringr::str_split(gf, "[/;]"))
    ga <- unlist(stringr::str_split(ga, "[/;]"))
    unf <- unlist(stringr::str_split(unf, "[/;]"))
    una <- unlist(stringr::str_split(una, "[/;]"))
    
    glycotopes <- gf[unf != "+"]
    glycotopes <- c(glycotopes, ga[una != "+"])
    glycotopes <- unique(glycotopes)
    glycotopes <- glycotopes[!is.na(glycotopes) & glycotopes != "" & glycotopes != "NA"]
    
    return(paste(glycotopes, collapse = ";"))
    
  }, gf = Glycotope_F, ga = Glycotope_A, unf = Glycotope_F_unmatched, una = Glycotope_A_unmatched))),
  by = .(seq_wind, `Glycan(H,N,A,G,F)`)]
  
  df_glycof_data2 <- df_glycof_data2[!duplicated(df_glycof_data2[, c("seq_wind", "Glycan(H,N,A,G,F)")])]
  df_glycof_data3 <- df_modpept[, head(.SD, 1),                       by = .(seq_wind, `Glycan(H,N,A,G,F)`), .SDcols = c("GlycanType", "GlycanAntennaType", "GlyMass", "Proteins", "ProSites", "Leading_Protein", "Leading_ProSite")]
  df_glycof_int   <- df_modpept[, lapply(.SD, sum, na.rm = TRUE),     by = .(seq_wind, `Glycan(H,N,A,G,F)`), .SDcols = c(int_names, mion_names)]
  
  df_modpept[, .(Glycotope_F = paste(Glycotope_F, collapse = ";"),
                 Glycotope_A = paste(Glycotope_A, collapse = ";"),
                 Glycotope_F_unmatched = paste(Glycotope_F_unmatched, collapse = ";"),
                 Glycotope_A_unmatched = paste(Glycotope_A_unmatched, collapse = ";")),
             by = .(seq_wind, `Glycan(H,N,A,G,F)`)]
  
  df_glycof <- cbind(df_glycof_data,
                     df_glycof_data2[, -c("seq_wind", "Glycan(H,N,A,G,F)")],
                     df_glycof_data3[, -c("seq_wind", "Glycan(H,N,A,G,F)")],
                     df_glycof_int[,   -c("seq_wind", "Glycan(H,N,A,G,F)")])
  names(df_glycof)[names(df_glycof) == "id"] <- "modpept_ids"
  df_glycof[, id := 1:df_glycof[, .N]]
  df_glycof <- df_glycof[, .SD, .SDcols = c("id", names(df_glycof)[names(df_glycof) != "id"])]
  
  rm(list = c("df_glycof_data",
              "df_glycof_data2",
              "df_glycof_int"))
  
  # combine intensities for glycosites
  df_gsite_data  <- df_modpept[, lapply(.SD, paste, collapse = ";"), by = .(seq_wind), .SDcols = c("id", "Scan", "pGlyco_ids", "Peptide", "GlySite", "GlyID", "Glycan(H,N,A,G,F)", "GlycanType", "GlycanAntennaType", "GlyMass")]
  df_gsite_data2 <- df_modpept[, .(Glycotope = unlist(Map(function(gf, ga, unf, una){
    
    gf <- unlist(stringr::str_split(gf, "[/;]"))
    ga <- unlist(stringr::str_split(ga, "[/;]"))
    unf <- unlist(stringr::str_split(unf, "[/;]"))
    una <- unlist(stringr::str_split(una, "[/;]"))
    
    glycotopes <- gf[unf != "+"]
    glycotopes <- c(glycotopes, ga[una != "+"])
    glycotopes <- glycotopes[!is.na(glycotopes) & glycotopes != "" & glycotopes != "NA"]
    glycotopes <- unique(glycotopes)
    
    return(paste(glycotopes, collapse = ";"))
    
  }, gf = Glycotope_F, ga = Glycotope_A, unf = Glycotope_F_unmatched, una = Glycotope_A_unmatched))),
  by = .(seq_wind)]
  df_gsite_data2 <- df_gsite_data2[!duplicated(df_gsite_data2[, c("seq_wind")])]
  
  df_gsite_data3 <- df_modpept[, head(.SD, 1),                       by = .(seq_wind), .SDcols = c("Proteins", "ProSites", "Leading_Protein", "Leading_ProSite")]
  df_gsite_int   <- df_modpept[, lapply(.SD, sum, na.rm = TRUE),     by = .(seq_wind), .SDcols = c(int_names, mion_names)]
  
  df_gsite <- cbind(df_gsite_data,
                    df_gsite_data2[, -c("seq_wind")],
                    df_gsite_data3[, -c("seq_wind")],
                    df_gsite_int[, -c("seq_wind")])
  names(df_gsite)[names(df_gsite) == "id"] <- "modpept_ids"
  df_gsite[, id := 1:df_gsite[, .N]]
  df_gsite <- df_gsite[, .SD, .SDcols = c("id", names(df_gsite)[names(df_gsite) != "id"])]
  
  # make a combined column of proteins and sites:
  df_gsite[, Protein_Site := lapply(Map(function(x, y){
    
    temp_x <- unlist(stringr::str_split(x, ";"))
    temp_y <- unlist(stringr::str_split(y, ";"))
    
    paste(temp_x, temp_y, sep = "_")
    
  }, x = df_gsite$Proteins, y = df_gsite$ProSites), paste, collapse = ";")]
  
  
  rm(list = c("df_gsite_data",
              "df_gsite_data2",
              "df_gsite_int"))
  
  # combine intensities for glycans
  
  df_glycan_data  <- df_modpept[, lapply(.SD, paste, collapse = ";"),
                                by = .(`Glycan(H,N,A,G,F)`),
                                .SDcols = c("id",
                                            "pGlyco_ids",
                                            "Scan",
                                            "Leading_Protein",
                                            "Leading_ProSite")]
  
  df_glycan_data2 <- df_modpept[, .(Glycotope = unlist(Map(function(gf, ga, unf, una){
    
    gf <- unlist(stringr::str_split(gf, "[/;]"))
    ga <- unlist(stringr::str_split(ga, "[/;]"))
    unf <- unlist(stringr::str_split(unf, "[/;]"))
    una <- unlist(stringr::str_split(una, "[/;]"))
    
    glycotopes <- gf[unf != "+"]
    glycotopes <- c(glycotopes, ga[una != "+"])
    glycotopes <- glycotopes[!is.na(glycotopes) & glycotopes != "" & glycotopes != "NA"]
    glycotopes <- unique(glycotopes)
    
    return(paste(glycotopes, collapse = ";"))
    
  }, gf = Glycotope_F, ga = Glycotope_A, unf = Glycotope_F_unmatched, una = Glycotope_A_unmatched))),
  by = .(`Glycan(H,N,A,G,F)`)]
  df_glycan_data2 <- df_glycan_data2[!duplicated(df_glycan_data2[, c("Glycan(H,N,A,G,F)")])]
  
  df_glycan_data3 <- df_modpept[, head(.SD, 1), 
                                by = .(`Glycan(H,N,A,G,F)`),
                                .SDcols = c("GlyID",
                                            "GlycanType",
                                            "GlycanAntennaType",
                                            "PlausibleStruct",
                                            "GlyFrag",
                                            "GlyMass")]
  
  df_glycan_int   <- df_modpept[, lapply(.SD, sum, na.rm = TRUE),
                                by = .(`Glycan(H,N,A,G,F)`),
                                .SDcols = c(int_names, mion_names)]
  df_glycan <- cbind(df_glycan_data,
                     df_glycan_data2[, -c("Glycan(H,N,A,G,F)")],
                     df_glycan_data3[, -c("Glycan(H,N,A,G,F)")],
                     df_glycan_int  [, -c("Glycan(H,N,A,G,F)")])
  names(df_glycan)[names(df_glycan) == "id"] <- "modpept_ids"
  df_glycan$id <- 1:df_glycan[, .N]
  df_glycan <- df_glycan[, .SD, .SDcols = c("id", names(df_glycan)[names(df_glycan) != "id"])]
  
  rm(list = c("df_glycan_data",
              "df_glycan_data2",
              "df_glycan_int"))
  
  # combine intensities for glycotopes
  # temp <- df_scans[, list(Glycotope = unlist(stringr::str_split(Glycotope, "&"))), by = "id"]
  # temp[, Glycotope := gsub(" ", "", Glycotope)]
  
  temp <- rbind(df_scans[, .(Glycotope = Glycotope_F,
                             Glycotope_unmatched = Glycotope_F_unmatched), by = "id"],
                df_scans[, .(Glycotope = Glycotope_A,
                             Glycotope_unmatched = Glycotope_A_unmatched), by = "id"])
  temp <- temp[!is.na(Glycotope) & Glycotope != ""]
  temp <- temp[Glycotope_unmatched != "+"]
  temp <- merge(temp, df_scans[, -c("Glycotope_F", "Glycotope_A", "Glycotope_F_unmatched", "Glycotope_A_unmatched")], by = "id", all.x = TRUE)
  
  df_glycotope_data  <- temp[, lapply(.SD, paste, collapse = ";"),
                             by = .(Glycotope),
                             .SDcols = c("id")]
  
  df_glycotope_int   <- temp[, lapply(.SD, sum, na.rm = TRUE),
                             by = .(Glycotope),
                             .SDcols = c(int_names, mion_names)]
  df_glycotope <- cbind(df_glycotope_data,
                        df_glycotope_int  [, -c("Glycotope")])
  
  df_glycotope <- rbind(df_glycotope,
                        data.table(Glycotope = base_glycotopes[!base_glycotopes %in% df_glycotope$Glycotope]),
                        fill = TRUE)
  df_glycotope[, Glycotope := factor(Glycotope, levels = ..base_glycotopes)]
  df_glycotope <- df_glycotope[order(Glycotope)]
  #df_glycotope <- df_glycotope[!is.na(Glycotope) & Glycotope != ""]
  
  names(df_glycotope)[names(df_glycotope) == "id"] <- "pGlyco_ids"
  df_glycotope[, id := 1:.N]
  df_glycotope <- df_glycotope[, .SD, .SDcols = c("id", names(df_glycotope)[names(df_glycotope) != "id"])]
  
  rm(list = c("df_glycotope_data",
              "df_glycotope_int"))
  
  dat <- list(
    
    scans   = df_scans,
    modpept = df_modpept,
    sites   = df_gsite,
    gforms  = df_glycof,
    glycans = df_glycan,
    glycotopes = df_glycotope
    
  )
  
  ##### check fucosylation #####
  use_tables <- c("sites", "gforms", "glycans")
  for(i in use_tables){
    
    dat[[i]] <- add_corefuc(dat[[i]], df_scans)
    dat[[i]] <- check_f_struct(dat[[i]], df_scans)
    dat[[i]] <- check_core_fucose(dat[[i]], df_scans)
    
  }
  
  ##### Add parent peak area #####
  scans <- df_scans[ParentPeakFound == TRUE] 
  use_tables <- c("modpept", "sites", "gforms", "glycans", "glycotopes")
  for(i in use_tables){
    
    if(nrow(dat[[i]]) == 0) next
    dat[[i]] <- addPeakArea(df = dat[[i]], scans = scans, FUN = parent_area_fun)
    
  }
  
  # convert reporter ion intensities into %
  if(reporter_ion != "not_labeled"){
    
    use_tables <- c("scans", "modpept", "sites", "gforms", "glycans", "glycotopes")
    for(i in use_tables){
      
      if(nrow(dat[[i]]) == 0) next
      dat[[i]] <- percentIntensity(dat[[i]], int_names = int_names)
      
    }
    
  }
  
  # write the tables
  if(verbose) message("Write output tables")
  fwrite(dat[["scans"]],   "pglyco_output\\pGlyco_Scans.txt", sep = "\t")
  fwrite(dat[["modpept"]], "pglyco_output\\pGlyco_modified_peptides.txt", sep = "\t")
  fwrite(dat[["sites"]],   "pglyco_output\\pGlyco_glycosites.txt", sep = "\t")
  fwrite(dat[["gforms"]],  "pglyco_output\\pGlyco_glycoforms.txt", sep = "\t")
  fwrite(dat[["glycans"]], "pglyco_output\\pGlyco_glycans.txt", sep = "\t")
  fwrite(dat[["glycotopes"]], "pglyco_output\\pGlyco_glycotopes.txt", sep = "\t")
  
  if(verbose) proc.time() - ptm_transform_out
  
})

message("\nDONE!\n")

if(verbose) proc.time() - ptm
writeLog(log_gb, new1 = paste0("Processing time: "), new2 = c(paste0("user: ", round((proc.time() - ptm)[1], 4)),
                                                              paste0("system: ", round((proc.time() - ptm)[2], 4)),
                                                              paste0("elapsed: ", round((proc.time() - ptm)[3], 4))))

