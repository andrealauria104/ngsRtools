#!/usr/bin/env Rscript
# Command line/Env variables ---
suppressPackageStartupMessages(library(docopt))
'Usage:
   get-illumina-ssheet2.r [-f <sheetfile> -s <sstart> -o <outdir>]

Options:
   -f, --sheetfile Path to excel workbook 
   -s, --sstart Sheet start row [default: 1]
   -o, --outdir outdir directory [default: .]
   -d, --date Report date in file name [default: FALSE]

Author: Andrea Lauria' -> doc

opts <- docopt(doc)
if(is.null(opts$sheetfile)) {
  message(doc)
  quit(save = "no", status = 0, runLast = TRUE)
}

sheetfile <- as.character(opts$sheetfile)
sstart    <- as.integer(opts$sstart)
outdir    <- as.character(opts$outdir)
rdate      <- as.logical(opts$date)
# 0. Resources ----
# packages
library(openxlsx)
options(stringsAsFactors=F)
# functions
read_sheetfile <- function(sheetfile, sstart)
{
  # Read Excel workbook
  rsheetfile <- tryCatch(
    { 
      wb      <- loadWorkbook(sheetfile) 
      sheets  <- names(wb) # retrive sheet names
      sheets  <- grep('dual|double|single',sheets, value = T)
      rsheets <- lapply(sheets, function(sheetname) 
        {
        y <- read.xlsx(sheetfile
          , sheet = sheetname
          , startRow = sstart
          , skipEmptyRows = T
          , skipEmptyCols = T
          , colNames = F)
        if(length(y[is.na(y)])!=0) y[is.na(y)] <- ""
	return(y)
      })
      names(rsheets) <- gsub("[[:space:]]","",sheets)
      names(rsheets) <- gsub("double","dual",sheets)
      return(rsheets)
    }, 
    error = function(e) {
      message(e)
      return(NA)
    } ) # read sheets
  return(rsheetfile)
}

write_to_file <- function(rsheetfile, sheetfile, outdir, rdate) 
{
  # Write outdir to file
  idate <- unlist(regmatches(sheetfile, gregexpr("\\d+", sheetfile)))
  outnames <- c("SampleSheet_dualindex.csv","SampleSheet_singleindex.csv")
  if(rdate) outnames <- gsub('.csv',paste0('_',idate,'.csv'),outnames)
  names(outnames) <- c('dual','single')
  
  for(i in names(rsheetfile)) {
    outfile <- paste0(outdir,"/",outnames[i])
    message(" -- ", outfile)
    write.table(rsheetfile[[i]], file = outfile, col.names = F, row.names = F, quote = F, sep = ',')
  }
}

# 1. Read Excel workbook ----
message("")
message("======================================================")
message("== [*] Illumina SampleSheet from Excel workbook [*] ==")
message("======================================================")
message("")

message(" -- reading data from: ", sheetfile)
rsheetfile <- read_sheetfile(sheetfile = sheetfile, sstart = sstart)

# 2. Write outdir to file ----
message(" -- writing data to: ", outdir)

write_to_file(rsheetfile, sheetfile, outdir, rdate)
