###########################################################################
## Lisa T. Haber                                                         ##
## November 2021                                                         ##
## FoRTE Subcanopy Leaf Functional Traits Manuscript Prep                ##
## 2018-2021 data integration & cleaning                                 ##
###########################################################################

## This code will extract and utilize subcanopy physiology (LICOR 6400XT) files from the shared FoRTE data drive.
## There are four separate data folders for the 2018-2021 field seasons' LICOR data on Google Drive. 

## load required packages
library(dplyr)
library(readr)
library(googledrive)
library(ggplot2)
library(tidyr)
library(lubridate)

#### I. Download successive years' data into data directories (only need to do this once)

###########NOTE: 12/6/2021 ALL data drives have changed location in migration of Google Drive to a university maintained share drive. Need to updat all links before running those download scripts again.

##################################################
### 2018:
##################################################
### 1. Bring in data
# Direct Google Drive link to "FoRTE/data/subcanopy_leaf_physiology/2018"
as_id("https://drive.google.com/drive/u/1/folders/1N4YvllZd-9zui9rQcJM_FKaM7C1DXf8X") %>% 
  drive_ls ->
  gdfiles

# Create a new data directory for files, if necessary
data_dir <- "data/LICOR_2018_RAW"
if(!dir.exists(data_dir)) dir.create(data_dir)

# Download data
# for(f in seq_len(nrow(gdfiles))) {
#   cat(f, "/", nrow(gdfiles), " Downloading ", gdfiles$name[f], "...\n", sep = "")
#   drive_download(gdfiles[f,], overwrite = TRUE, path = file.path(data_dir, gdfiles$name[f]))
# }

# Get a (fresh) list of the downloaded data we're working with
# Filenames we want end with eight digits and no file extension
files <- list.files(data_dir, pattern = "[0-9]{8}$", full.names = TRUE)
HEADER_PATTERN <- "\"OPEN \\d\\.\\d\\.\\d"
DATA_PATTERN <- "\\$STARTOFDATA\\$"

# Scan through all the data files and read data into list structure
filedata <- list()
for(f in files) {
  cat(" Reading ", f, "...\n", sep = "")
  text_raw <- readLines(f, skipNul = TRUE)
  data_start <- grep(DATA_PATTERN, text_raw)
  first_comment <- text_raw[data_start - 1] # there's always a comment on this line
  
  if(length(data_start)) {
    # What makes this tricky is that there can be additional comments WITHIN the data frame
    # Who on earth thought that was a good idea?!?
    data_raw <- text_raw[data_start+1:length(text_raw)] %>% na.omit
    line_lengths <- lapply(strsplit(data_raw, "\t"), length) %>% unlist
    data_rows <- line_lengths == line_lengths[1]
    comments <- paste(which(!data_rows), data_raw[!data_rows], sep = ". ") %>%
      paste(first_comment, ., sep = "; ") %>% 
      gsub('\"', "", .)
    
    # OK, now read the data into a data frame and add the 'comments'
    con <- textConnection(data_raw[data_rows])
    read.table(con, header = TRUE, stringsAsFactors = FALSE) %>% 
      mutate(Filename = basename(f),
             Timestamp = text_raw[grep(HEADER_PATTERN, text_raw) + 1],
             Comments = paste(comments, collapse = "; ")) ->
      filedata[[f]]
    close(con)
  }
}

# Combine data into a single data frame for analysis
filedata %>% 
  bind_rows %>% 
  as_tibble %>% 
  mutate(Timestamp = mdy_hms(Timestamp)) %>%  # change to a POSIXct object
  separate(Filename, into = c("Plot", "Species", "Sample", "Filename_date"), remove = FALSE) ->
  licordata2018

unique(licordata2018$Filename_date)

### 2: exclude excessive observations
# figure out which dates/files have more than 5 observations
excessobs <- licordata2018 %>%
  subset(Obs >= 6)

dates <- excessobs$Filename

tibble::view(excessobs)
print(unique(dates))

# remove excess observations
# figure out which rows in dataframe need removal (easiest to just do this by hand)
# these are observations ("logged" measurements) on the LICOR that I accidentally made before a leaf had fully stabilized, or just from hitting the log button too many times by accident
tibble::view(licordata2018)

# now, exclude them. Note that this means the affected files will have observations beginning at "2" or "6", not "1"
drop <- c(191, 567:571, 587:591, 662:666, 1117:1121, 1697:1701)

licordata2018_cleaned <- licordata2018[-drop, ] #this subsets the original dataframe to exclude unwanted measurements


### 3: attach unique subcanopy tree IDs to all trees from 2018
# use lookup table

ID_lookup <- read.csv("data/forte_subcanopy_tree_id_master_lookup.csv") %>%
  as_tibble()

licordata2018_cleaned <- left_join(licordata2018_cleaned, ID_lookup, by = "Filename", keep = TRUE) 
tibble::view(licordata2018_cleaned)

# drop variables we don't want retained in final data set
# this includes erroneous subplot information in column "Plot" deriving from the Licor filename itself, which was sometimes incorrectly entered in 2018

head(licordata2018_cleaned)

licordata2018_cleaned <- licordata2018_cleaned %>%
  select(-Plot, -Species.y, -Sample, -Filename.y, -ID_tag_2018, -ID_tag_added_2019_or_after, -Notes_2019, -Notes_2020, -Notes_2021) %>%
  rename(Species = Species.x) %>%
  rename(Filename = Filename.x) %>%
  mutate(Year = 2018)
  
# check to see what overall non-negative sample size is
licordata2018_means <- licordata2018_cleaned %>%
  group_by(Filename) %>%
  summarize(MeanPhoto = mean(Photo)) %>%
  filter(MeanPhoto >= 0) # I will leave the negative values in the "raw" data for fortedata; my analyses will ultimately exclude them

# we exclude negative photosynthetic values because they are not usable measurements

nrow(licordata2018_means)#390 trees represented in the data set; 2 excluded by negative measurements

### 4: time to add the covered IRGA area for leaves. LI-6400 assumes 2x3 cm2 coverage at IRGA mouth but some needle-leaf specimens were measured, and their measurements will have to be scaled accordingly

# bring in the 2018 leaf morphology data set
morph2018 <- read.csv("data/leaf_morphology/FoRTE_subcanopy_leaf_morphology_2018.csv") %>%
  as_tibble()

# join to LICOR data
morph2018 <- morph2018 %>%
  rename(final_ID = tag_number) %>%
  select(final_ID, IRGA_covered_area, Asat_multiplier)
  
scaled_licordata_2018 <- left_join(licordata2018_cleaned, morph2018, by = "final_ID") %>%
  select(-Filename_date) %>%
  mutate(girdled = "n") %>%
  select(-girdled.) %>%
  rename(girdled. = girdled)

# write.csv(scaled_licordata_2018, "data/cleaned_data_for_analyses/LICOR_leaf_phys_2018.csv", row.names = FALSE)

#######################################################
### 2019:
#######################################################
# Direct Google Drive link to "FoRTE/data/subcanopy_leaf_physiology/2018"
as_id("https://drive.google.com/drive/u/1/folders/1h8FKmyZCoLET9XZjFmt_1-7bqHshFwdV") %>% 
  drive_ls ->
  gdfiles

# Create a new data directory for files, if necessary
data_dir <- "data/LICOR_2019_RAW"
if(!dir.exists(data_dir)) dir.create(data_dir)

# Download data
# for(f in seq_len(nrow(gdfiles))) {
#   cat(f, "/", nrow(gdfiles), " Downloading ", gdfiles$name[f], "...\n", sep = "")
#   drive_download(gdfiles[f,], overwrite = TRUE, path = file.path(data_dir, gdfiles$name[f]))
# }

# Get a (fresh) list of the downloaded data we're working with
# Filenames we want end with four digits and no file extension
files <- list.files(data_dir, pattern = "[0-9]{4}$", full.names = TRUE)
HEADER_PATTERN <- "\"OPEN \\d\\.\\d\\.\\d"
DATA_PATTERN <- "\\$STARTOFDATA\\$"

# Scan through all the data files and read data into list structure
filedata <- list()
for(f in files) {
  cat(" Reading ", f, "...\n", sep = "")
  text_raw <- readLines(f, skipNul = TRUE)
  data_start <- grep(DATA_PATTERN, text_raw)
  first_comment <- text_raw[data_start - 1] # there's always a comment on this line
  
  if(length(data_start)) {
    # What makes this tricky is that there can be additional comments WITHIN the data frame
    # Who on earth thought that was a good idea?!?
    data_raw <- text_raw[data_start+1:length(text_raw)] %>% na.omit
    line_lengths <- lapply(strsplit(data_raw, "\t"), length) %>% unlist
    data_rows <- line_lengths == line_lengths[1]
    comments <- paste(which(!data_rows), data_raw[!data_rows], sep = ". ") %>%
      paste(first_comment, ., sep = "; ") %>%
      gsub('\"', "", .)
    
    # OK, now read the data into a data frame and add the 'comments'
    con <- textConnection(data_raw[data_rows])
    read.table(con, header = TRUE, stringsAsFactors = FALSE) %>% 
      mutate(Filename = basename(f),
             Timestamp = text_raw[grep(HEADER_PATTERN, text_raw) + 1],
             Comments = paste(comments, collapse = "; ")) ->
      filedata[[f]]
    close(con)
  }
}

# Combine data into a single data frame for analysis
filedata %>% 
  bind_rows %>% 
  as_tibble %>% 
  mutate(Timestamp = mdy_hms(Timestamp)) ->  # change to a POSIXct object
  # separate(Filename, into = c("Plot", "Species", "Sample", "Filename_date"), remove = FALSE) ->
  licordata2019 

## code to generate the clean 2019 phys data 
excessobs <- licordata2019 %>%
  subset(Obs >= 6)

dates <- excessobs$Filename

tibble::view(excessobs)
print(unique(dates))

# no excess observations in the 2019 data set! Crazy! Must have cleaned these at some earlier date....

# check to see what overall non-negative sample size is
licordata_2019_means <- licordata2019 %>%
  group_by(Filename) %>%
  summarize(MeanPhoto = mean(Photo)) %>%
  filter(MeanPhoto >= 0)

nrow(licordata_2019_means)

# join to lookup table to get species, subplot, etc.

licordata2019 <- licordata2019 %>%
  mutate(final_ID = as.numeric(licordata2019$Filename))

licordata2019_cleaned <- left_join(licordata2019, ID_lookup, by = "final_ID", keep = TRUE) 

# remove unwanted columns
licordata2019_cleaned <- licordata2019_cleaned %>%
  select(-Filename.y, -ID_tag_2018, -ID_tag_added_2019_or_after, -final_ID.x, -Notes_2019, -Notes_2020, -Notes_2021) %>%
  rename(final_ID = final_ID.y) %>%
  rename(Filename = Filename.x) %>%
  mutate(Year = 2019)

# join to leaf area data for scaling pine samples
# bring in morphology data from 2019
morph2019 <- read.csv("data/leaf_morphology/FoRTE_subcanopy_leaf_morphology_2019.csv")

morph2019 <- morph2019 %>%
  mutate(final_ID = as.numeric(tag_number)) %>%
  select(final_ID, IRGA_covered_area, Asat_multiplier) 

scaled_licordata_2019 <- left_join(licordata2019_cleaned, morph2019, by = "final_ID")

# write.csv(scaled_licordata_2019, "data/cleaned_data_for_analyses/LICOR_leaf_phys_2019.csv", row.names = FALSE)

# check2018 <- read.csv("data/cleaned_data_for_analyses/LICOR_leaf_phys_2018.csv")
# check2019 <- read.csv("data/cleaned_data_for_analyses/LICOR_leaf_phys_2019.csv")
# 
# columnssame <- merge(check2018, check2019)
# head(columnssame)

#######################################################
### 2020:
#######################################################
# Direct Google Drive link to "FoRTE/data/subcanopy_leaf_physiology/2020"
as_id("https://drive.google.com/drive/u/1/folders/1qUe1SwSO7v6scWOIUYghyOTOKebAbrWM") %>% 
  drive_ls ->
  gdfiles

# Create a new data directory for files, if necessary
data_dir <- "data/LICOR_2020_RAW"
if(!dir.exists(data_dir)) dir.create(data_dir)

# Download data
# for(f in seq_len(nrow(gdfiles))) {
#   cat(f, "/", nrow(gdfiles), " Downloading ", gdfiles$name[f], "...\n", sep = "")
#   drive_download(gdfiles[f,], overwrite = TRUE, path = file.path(data_dir, gdfiles$name[f]))
# }

# Get a (fresh) list of the downloaded data we're working with
# Filenames we want end with four digits and no file extension
files <- list.files(data_dir, pattern = "[0-9]{8}$", full.names = TRUE)
HEADER_PATTERN <- "\"OPEN \\d\\.\\d\\.\\d"
DATA_PATTERN <- "\\$STARTOFDATA\\$"

# Scan through all the data files and read data into list structure
filedata <- list()
for(f in files) {
  cat(" Reading ", f, "...\n", sep = "")
  text_raw <- readLines(f, skipNul = TRUE)
  data_start <- grep(DATA_PATTERN, text_raw)
  first_comment <- text_raw[data_start - 1] # there's always a comment on this line
  
  if(length(data_start)) {
    # What makes this tricky is that there can be additional comments WITHIN the data frame
    # Who on earth thought that was a good idea?!?
    data_raw <- text_raw[data_start+1:length(text_raw)] %>% na.omit
    line_lengths <- lapply(strsplit(data_raw, "\t"), length) %>% unlist
    data_rows <- line_lengths == line_lengths[1]
    comments <- paste(which(!data_rows), data_raw[!data_rows], sep = ". ") %>%
      paste(first_comment, ., sep = "; ") %>% 
      gsub('\"', "", .)
    
    # OK, now read the data into a data frame and add the 'comments'
    con <- textConnection(data_raw[data_rows])
    read.table(con, header = TRUE, stringsAsFactors = FALSE) %>% 
      mutate(Filename = basename(f),
             Timestamp = text_raw[grep(HEADER_PATTERN, text_raw) + 1],
             Comments = paste(comments, collapse = "; ")) ->
      filedata[[f]]
    close(con)
  }
}

# Combine data into a single data frame for analysis
filedata %>% 
  bind_rows %>% 
  as_tibble %>% 
  mutate(Timestamp = mdy_hms(Timestamp)) %>%  # change to a POSIXct object
  separate(Filename, into = c("Sample", "Filename_date"), remove = FALSE) ->
  licordata2020 

##############################################################
# ## code to generate the clean 2020 phys data 
excessobs <- licordata2020 %>%
  subset(Obs >= 6)

dates <- excessobs$Filename

tibble::view(excessobs)
print(unique(dates))

# remove excess observations
# figure out which rows in dataframe need removal (easiest to just do this by hand at this point)
tibble::view(licordata2020)

# now, exclude them. First two files had accidental 6th measurement taken, last one I logged way too early on first log and got a negative value
drop <- c(1150, 1401, 1842)

licordata2020_cleaned <- licordata2020[-drop, ]

tibble::view(licordata2020_cleaned)  

nrow(licordata2020_cleaned)  

# check to see what overall non-negative sample size is
licordata_2020_means <- licordata2020_cleaned %>%
  group_by(Filename) %>%
  summarize(MeanPhoto = mean(Photo)) %>%
  filter(MeanPhoto >= 0)

nrow(licordata_2020_means)

# join to lookup table to get species, subplot, etc.

licordata2020_cleaned <- licordata2020_cleaned %>%
  mutate(final_ID = as.numeric(licordata2020_cleaned$Sample))

licordata2020_cleaned <- left_join(licordata2020_cleaned, ID_lookup, by = "final_ID", keep = TRUE) 

# remove unwanted columns
licordata2020_cleaned <- licordata2020_cleaned %>%
  select(-Filename.y, -ID_tag_2018, -ID_tag_added_2019_or_after, -Sample, -Filename_date, -final_ID.x, -Notes_2019, -Notes_2020, -Notes_2021) %>%
  rename(final_ID = final_ID.y) %>%
  rename(Filename = Filename.x) %>%
  mutate(Year = 2020)

# join to leaf area data for scaling pine samples
# bring in morphology data from 2020
morph2020 <- read.csv("data/leaf_morphology/FoRTE_subcanopy_leaf_morphology_2020.csv") %>%
  as_tibble()

morph2020 <- morph2020 %>%
  mutate(final_ID = as.numeric(tag_number)) %>%
  select(final_ID, IRGA_covered_area, Asat_multiplier) 

scaled_licordata_2020 <- left_join(licordata2020_cleaned, morph2020, by = "final_ID")

# write.csv(scaled_licordata_2020, "data/cleaned_data_for_analyses/LICOR_leaf_phys_2020.csv", row.names = FALSE)


##################################################################
### 2021
##################################################################
# Direct Google Drive link to "FoRTE/data/subcanopy_leaf_physiology/2021"
as_id("https://drive.google.com/drive/u/1/folders/1DKCJwYsYpCR7VHp3twSc5APFRiDZVAT7") %>% 
  drive_ls ->
  gdfiles

# Create a new data directory for files, if necessary
data_dir <- "data/LICOR_2021_RAW"
if(!dir.exists(data_dir)) dir.create(data_dir)

# Download data
for(f in seq_len(nrow(gdfiles))) {
  cat(f, "/", nrow(gdfiles), " Downloading ", gdfiles$name[f], "...\n", sep = "")
  drive_download(gdfiles[f,], overwrite = TRUE, path = file.path(data_dir, gdfiles$name[f]))
}

# Get a (fresh) list of the downloaded data we're working with
# Filenames we want end with four digits and no file extension
files <- list.files(data_dir, pattern = "[0-9]{8}$", full.names = TRUE)
HEADER_PATTERN <- "\"OPEN \\d\\.\\d\\.\\d"
DATA_PATTERN <- "\\$STARTOFDATA\\$"

# Scan through all the data files and read data into list structure
filedata <- list()
for(f in files) {
  cat(" Reading ", f, "...\n", sep = "")
  text_raw <- readLines(f, skipNul = TRUE)
  data_start <- grep(DATA_PATTERN, text_raw)
  first_comment <- text_raw[data_start - 1] # there's always a comment on this line
  
  if(length(data_start)) {
    # What makes this tricky is that there can be additional comments WITHIN the data frame
    # Who on earth thought that was a good idea?!?
    data_raw <- text_raw[data_start+1:length(text_raw)] %>% na.omit
    line_lengths <- lapply(strsplit(data_raw, "\t"), length) %>% unlist
    data_rows <- line_lengths == line_lengths[1]
    comments <- paste(which(!data_rows), data_raw[!data_rows], sep = ". ") %>%
      paste(first_comment, ., sep = "; ") %>% 
      gsub('\"', "", .)
    
    # OK, now read the data into a data frame and add the 'comments'
    con <- textConnection(data_raw[data_rows])
    read.table(con, header = TRUE, stringsAsFactors = FALSE) %>% 
      mutate(Filename = basename(f),
             Timestamp = text_raw[grep(HEADER_PATTERN, text_raw) + 1],
             Comments = paste(comments, collapse = "; ")) ->
      filedata[[f]]
    close(con)
  }
}

# Combine data into a single data frame for analysis
filedata %>% 
  bind_rows %>% 
  as_tibble %>% 
  mutate(Timestamp = mdy_hms(Timestamp)) %>%  # change to a POSIXct object
  separate(Filename, into = c("Sample", "Filename_date"), remove = FALSE) ->
  licordata2021 

##############################################################
# ## code to generate the clean 2021 phys data 
excessobs <- licordata2021 %>%
  subset(Obs >= 6)

dates <- excessobs$Filename

tibble::view(excessobs)
print(unique(dates))

# remove excess observations
# figure out which rows in dataframe need removal (easiest to just do this by hand at this point)
tibble::view(licordata2021)

# now, exclude them. Different files had accidental logs for different reasons, either failure to stabilize on initial logs or excessive logs by accident. Cleaned according to field notebook notes from field data collection
drop <- c(206:208, 334:336, 797, 1148:1152, 1528:1531, 1822, 1883:1885)

licordata2021_cleaned <- licordata2021[-drop, ]

tibble::view(licordata2021_cleaned)  

nrow(licordata2021_cleaned)  

# check to see what overall non-negative sample size is
licordata_2021_means <- licordata2021_cleaned %>%
  group_by(Filename) %>%
  summarize(MeanPhoto = mean(Photo)) %>%
  filter(MeanPhoto >= 0)

nrow(licordata_2021_means)

# join to lookup table to get species, subplot, etc.

licordata2021_cleaned <- licordata2021_cleaned %>%
  mutate(final_ID = as.numeric(licordata2021_cleaned$Sample))

licordata2021_cleaned <- left_join(licordata2021_cleaned, ID_lookup, by = "final_ID", keep = TRUE) 

# remove unwanted columns
licordata2021_cleaned <- licordata2021_cleaned %>%
  select(-Filename.y, -ID_tag_2018, -ID_tag_added_2019_or_after, -Sample, -Filename_date, -final_ID.x, -Notes_2019, -Notes_2020, -Notes_2021) %>%
  rename(final_ID = final_ID.y) %>%
  rename(Filename = Filename.x) %>%
  mutate(Year = 2021)

# join to leaf area data for scaling pine samples
# bring in morphology data from 2021
morph2021 <- read.csv("data/leaf_morphology/FoRTE_subcanopy_leaf_morphology_2021.csv")

morph2021 <- morph2021 %>%
  mutate(final_ID = as.numeric(tag_number)) %>%
  select(final_ID, IRGA_covered_area, Asat_multiplier) 

scaled_licordata_2021 <- left_join(licordata2021_cleaned, morph2021, by = "final_ID")

# write.csv(scaled_licordata_2021, "data/cleaned_data_for_analyses/LICOR_leaf_phys_2021.csv", row.names = FALSE)

