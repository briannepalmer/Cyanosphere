library(tidyr)
library(purrr)
library(readr)

path <- "/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/Cyanosphere/data/DRAM"


# in one pipeline:
DRAM  <- path %>% 
  
  # get csvs full paths. (?i) is for case insentitive
  list.files(pattern = "(?i)\\.tsv$", full.names = TRUE) %>% 
  
  # create a named vector: you need it to assign ids in the next step.
  # and remove file extection to get clean colnames
  set_names(tools::file_path_sans_ext(basename(.))) %>% 
  
  # read file one by one, bind them in one df and create id column 
  map_dfr(read_tsv) 


write.csv(DRAM, file = "DRAM.csv")

DRAM <- read.csv("DRAM.csv")
DRAM <- DRAM[,-1]

taxonomy <- read.csv("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/Cyanosphere/data/taxonomy.csv")
metadata <- read.csv("data/metdata_Feb2024.csv")

DRAM.Merged <- merge(taxonomy, DRAM, by = "genome", fill = TRUE)
DRAM.Merged <- merge(DRAM.Merged, metadata, by = "Host", fill = TRUE)

write.csv(DRAM.Merged, file = "DRAM_w_Tax.csv" )
