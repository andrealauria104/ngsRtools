#!/usr/bin/env Rscript
# Command line/Env variables ---
suppressWarnings(suppressMessages(require(docopt)))
'Usage:
   get-illumina-ssheet.r [-f <sheetfile> -s <sstart> -o <outdir>]

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
suppressWarnings(suppressMessages(require(xlsx)))

# functions
read_sheetfile <- function(sheetfile, sstart)
{
  # Read Excel workbook
  rsheetfile <- tryCatch(
    { 
      wb      <- loadWorkbook(sheetfile) 
      sheets  <- names(getSheets(wb)) # retrive sheet names
      sheets  <- grep('dual|double|single',sheets, value = T)
      rsheets <- lapply(sheets, function(sheetname) 
        {
        y <- suppressWarnings(read.xlsx(sheetfile, sheetName = sheetname, startRow = sstart, stringsAsFactors = F, header = F))
	      naidx <- which(apply(y, 1, function(x) all(is.na(x))))
	      if(length(naidx)>0) y <- y[-naidx,]
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
    df <- sapply(rsheetfile[[i]], as.character)
    df[is.na(df)] <- ""
    write.table(df, file = outfile, col.names = F, row.names = F, quote = F, sep = ',')
  }
}

# 1. Read Excel workbook ----
message("")
message("======================================================")
message("== [*] Illumina SampleSheet from Excek workbook [*] ==")
message("======================================================")
message("")

message(" -- reading data from: ", sheetfile)
rsheetfile <- read_sheetfile(sheetfile = sheetfile, sstart = sstart)

# 2. Write outdir to file ----
message(" -- writing data to: ", outdir)

write_to_file(rsheetfile, sheetfile, outdir, rdate)
