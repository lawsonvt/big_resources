#library(openxlsx)
library(readxl)
library(dplyr)

root_dir <- "/Users/mjl3p/Documents/projects/harrislab/ZBP1/ImmuneCellandToxoData/"

#xlsx_files <- list.files(root_dir, ".xlsx", full.names = T, recursive = T)
#xls_files <- list.files(root_dir, ".xls$", full.names = T, recursive = T)

excel_files <- list.files(root_dir, "(\\.xls$)|(\\.xlsx$)", full.names = T, recursive = T)

init_csv_files <- list.files(root_dir, ".csv", full.names = T, recursive = T)

# filter out open files
excel_files <- excel_files[!grepl("\\~\\$", excel_files)]

# first test "can we open every file"

for (file in excel_files) {
  
  print(file)
  
  print(excel_sheets(file))
  
}
# no issues

# are there duplicate files
base_files <- basename(excel_files)
names(base_files) <- excel_files

dupe_files <- base_files[duplicated(base_files)]
# does this matter?

# elisa files
elisa_files <- excel_files[grepl("(ELISA)|(elisa)", excel_files)]

# simple conversion

convert_files <- lapply(excel_files, function(file) {
  
  print(file)
  
  sheet_names <- excel_sheets(file)
  
  # if only one sheet, read that in
  if (length(sheet_names) == 1) {
    
    data <- read_excel(file)
    
  } else {
    # if multiple sheets, first see if any are labeled "Sheet"
    sheets_named_sheet <- grep("Sheet", sheet_names, value=T)
    
    # if none are, just read in first sheet
    if (length(sheets_named_sheet) == 0) {
      
      data <- read_excel(file)
      
    } else {
      # if at least one is labeled as Sheet, read in first one
      data <- read_excel(file, sheet = sheets_named_sheet[1])
      
    }
    
  }
  
  # create output file
  csv_file <- gsub("(\\.xls$)|(\\.xlsx$)", ".csv", file)
  
  write.csv(data, csv_file, row.names = F)
  
  return(csv_file)
})

# delete excel files
for (file in excel_files) {
  print(paste0("Deleting ", file))
  unlink(file)
}

# read in what was converted and pull out column names
col_names <- lapply(convert_files, function(file) {
  
  
  data <- read.csv(file)
  
  file_simple <- gsub(root_dir, "", file)
  
  data_columns <- colnames(data)
  
  # remove columns added in by read in function
  data_columns <- data_columns[!grepl("^\\.\\.\\.[0-9]+$", data_columns)]
  
  data.frame(filename=file_simple,
             raw_columns=paste0(data_columns, collapse="|"),
             output_columns=paste0("\n\t", paste0(data_columns, collapse="\n\t")))
  
})
col_names <- bind_rows(col_names)

# output a readme
write.table(col_names[,c(1,3)], file=paste0(root_dir, "raw_readme.txt"), quote=F, 
            row.names = F, col.names = F)

# eliminate all non-csv files

all_files <- list.files(root_dir, recursive = T, full.names = T)

csv_files <- list.files(root_dir, ".csv", recursive = T, full.names = T)

non_csv_files <- all_files[!all_files %in% csv_files]

unlink(non_csv_files)


