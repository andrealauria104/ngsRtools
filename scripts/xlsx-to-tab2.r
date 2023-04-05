#!/usr/bin/env Rscript
# Command line/Env variables ---
suppressPackageStartupMessages(library(docopt))
'Usage:
   xlsx-to-tab2.r [-f <sheetfile> -s <sstart> -o <output>]

Options:
   -f, --sheetfile Path to excel workbook 
   -s, --sstart Sheet start row [default: 1]
   -o, --output Output file [default: workbook.txt]

Author: Andrea Lauria' -> doc

opts <- docopt(doc)
if(is.null(opts$sheetfile)) {
  message(doc)
  quit(save = "no", status = 0, runLast = TRUE)
}

sheetfile <- as.character(opts$sheetfile)
sstart    <- as.integer(opts$sstart)
output    <- as.character(opts$output)
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
      rsheets <- lapply(sheets, function(sheetname)
        {
        y <- read.xlsx(sheetfile
          , sheet = sheetname
          , startRow = sstart
          , skipEmptyRows = T
          , skipEmptyCols = T
          , colNames = T
          )
        return(y)
      })
      names(rsheets) <- gsub("[[:space:]]","",sheets)
      return(rsheets)
    },
    error = function(e) {
      message(e)
      return(NA)
    } ) # read sheets
  return(rsheetfile)
}

write_to_file <- function(rsheetfile, output) 
{
  # Write output to file
  for(i in 1:length(rsheetfile)) {
    
    write(paste0(">",names(rsheetfile)[i]), file = output, append = ifelse(i==1, F,T))
    write.table(rsheetfile[[i]]
		, file = output
                , append = T
                , col.names = T
                , row.names = F
                , quote = F
                , sep = "\t")
  }
}

# 1. Read Excel workbook ----
message("")
message("===================================================")
message("== [*] Excel to tab-delimited text converter [*] ==")
message("===================================================")
message("")

message(" -- reading data from: ", sheetfile)
rsheetfile <- read_sheetfile(sheetfile = sheetfile, sstart = sstart)

# 2. Write output to file ----
if(output=="workbook.txt") output <- gsub("\\.xls(x)?$","\\.txt",sheetfile)
message(" -- writing data to: ", output)

write_to_file(rsheetfile = rsheetfile, output = output)
