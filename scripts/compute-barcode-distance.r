#!/usr/bin/env Rscript
# Command line/Env variables ---
suppressWarnings(suppressMessages(require(docopt)))
'Usage:
   compute-barcode-distance.r [-f <sheetfile> -s <sstart> -d <dstart> -m <maxdist> -v <stdout> -o <output>]

Options:
   -f, --sheetfile Path to sheetfile (excel workbook) or 2-column table with sample-barcode list (.txt) 
   -s, --sstart Single index sheet start row [default: 1]
   -d, --dstart Double index sheet start row [default: 1]
   -m, --maxdist Check collisions up to max distance [default: 2]
   -v, --stdout Print results to std output [default: TRUE]
   -o, --output Output file [default: none]

Author: Andrea Lauria' -> doc

opts <- docopt(doc)
if(is.null(opts$sheetfile)) {
  message(doc)
  quit(save = "no", status = 0, runLast = TRUE)
}

sheetfile   <- as.character(opts$sheetfile)
singlestart <- as.integer(opts$sstart)
doublestart <- as.integer(opts$dstart)
mdist       <- as.integer(opts$maxdist)
stdout      <- as.logical(opts$stdout)
output      <- as.character(opts$output)
# 0. Resources ----
# packages
suppressWarnings(suppressMessages(require(sequenceAnalysisRtools)))

# functions
read_sheetfile <- function(sheetfile, sheetName, startRow)
{
  rsheetfile <- tryCatch(
    { 
      if(grepl("\\.xlsx", sheetfile)) {
        suppressWarnings(read.xlsx(sheetfile, sheetName = sheetName, startRow = startRow, stringsAsFactors = F))
      } else if(grepl("\\.txt", sheetfile)) {
        read.delim2(sheetfile, header = T, stringsAsFactors = F)
      } else if(grepl("\\.csv", sheetfile)) {
        read.delim2(sheetfile, header = T, stringsAsFactors = F, sep = ",")
      } 
    }, 
    error = function(e) {
      message(e)
      return(NA)
    } )
  
  if(class(rsheetfile)!="data.frame" || !any(grepl("Sample|index$|ID|id",colnames(rsheetfile)))) {
    stop(message("[!] Error reading file: ", sheetfile)) 
  } else {
    return(rsheetfile) 
  }
}
calc_allowed_mismatches <- function(x)
{
  if(x>0 && x<3) {
    return(0)
  } else if(3<=x && x<5) {
    return(1)
  } else if(x>=5) {
    return(2) 
  } else if(x==0) {
    return("demultiplexing error, exact barcode collision!") 
  }
}
write_to_file <- function(x, output, append = T) 
{
  if(is.character(x)) {
    write(paste0(x, collapse = "\t"), file = output, append = append)
  } else if(is.data.frame(x)) {
    suppressWarnings(write.table(x, file = output, append = append, col.names = T, row.names = F, quote = F, sep = "\t"))
  }
}
message("")
message("=========================================")
message("== [*] Barcode distance calculator [*] ==")
message("=========================================")
message("")
# 1. Read Sample Sheets ----

message(" -- reading data from: ", sheetfile)

singlesheet <- read_sheetfile(sheetfile, sheetName = "single", startRow = singlestart)

if(grepl("\\.xlsx", sheetfile)){
  doublesheet <- read_sheetfile(sheetfile, sheetName = "double", startRow = doublestart)
} 

# 2. Prepare single/double index ----
if(exists("doublesheet")) {
  # doublesheet$index <- substr(doublesheet$index, 1,6)
  index_comb <- c(unique(as.character(doublesheet$index)), as.character(singlesheet$index))
} else {
  index_comb <- as.character(singlesheet$index)
}

# 3. Calculate distances ----
cat("\n")
d <- calculate_string_distance(index_comb, sub.method = "trim")
d <- unlist(d)

x <- d[which(d<=mdist)]

message("\n -- minimum index distance = ", min(d))
message(" -- maximum allowed mismatches = ", calc_allowed_mismatches(min(d)))

# check collapsing indexes

if(min(d)<=mdist && stdout){
  
  message("\n -- check collapsing index combinations, total = ",length(x))
  
  for(i in names(x)) {
    cat("\n")
    cat(paste0(i, ", distance = ", x[i]),"\n")
    index_1 <- strsplit(i,"\\-")[[1]][1]
    index_2 <- strsplit(i,"\\-")[[1]][2]
    
    if(exists("doublesheet")) {
      iid <- grep("index_[1,2]",ls(), value = T)
      for(id in iid) {
        if(length(grep(get(id),doublesheet$index))!=0) {
          print(doublesheet[grep(get(id),doublesheet$index),])  
        }
        if(length(grep(get(id),singlesheet$index))!=0) {
          print(singlesheet[grep(get(id),singlesheet$index),])  
        }
      }
      
    } else {
      
      if(index_1!=index_2) {
        print(singlesheet[match(index_1,singlesheet$index),])
        print(singlesheet[match(index_2,singlesheet$index),])
      } else if(index_1 == index_2) {
        print(singlesheet[grep(index_1, singlesheet$index),])
      }
    }
  }
  cat("\n")
}

# 4. [Optional] Write output to file ----
if(output!='none' && min(d)<=mdist) {
  
  message("\n -- writing to file = ", output,"\n")
  
  for(i in names(x)) {
    
    index_1 <- strsplit(i,"\\-")[[1]][1]
    index_2 <- strsplit(i,"\\-")[[1]][2]
    
    write_to_file(paste0(i, ", distance = ", x[i]), output, append = ifelse(i == names(x)[1], F, T))
    if(exists("doublesheet")) {
      iid <- grep("index_[1,2]",ls(), value = T)
      for(id in iid) {
        
        if(length(grep(get(id),doublesheet$index))!=0) {
          # double --
          write_to_file(doublesheet[grep(get(id),doublesheet$index),], output = output)  
        }
        if(length(grep(get(id),singlesheet$index))!=0) {
          # single --
          write_to_file(singlesheet[grep(get(id),singlesheet$index),], output = output)  
        }
      }
      
    } else {
      # single
      write_to_file(colnames(singlesheet), output)
      
      index_1 <- strsplit(i,"\\-")[[1]][1]
      index_2 <- strsplit(i,"\\-")[[1]][2]
      if(index_1!=index_2) {
        write_to_file(as.character(singlesheet[match(index_1,singlesheet$index),]), output)
        write_to_file(as.character(singlesheet[match(index_2,singlesheet$index),]), output)
      } else if(index_1 == index_2) {
        write_to_file(as.character(singlesheet[grep(index_1, singlesheet$index),][1,]), output)
        write_to_file(as.character(singlesheet[grep(index_1, singlesheet$index),][2,]), output)
      }
      
    }
    write_to_file("", output = output)
  }
  rm(index_1, index_2)
}

message("\n ... done! \n")
