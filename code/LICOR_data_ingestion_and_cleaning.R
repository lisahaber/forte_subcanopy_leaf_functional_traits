################################################
## Lisa T. Haber                              ##
## 11/8/2021                                  ##
## FoRTE Subcanopy Physiology                 ##
## 2018-2021 data integration & cleaning      ##
################################################

# This code will extract and utilize subcanopy physiology (LICOR 6400XT) files from the shared FoRTE data drive.
# There are four separate data folders for the 2018-2021 field seasons' LICOR data on Google Drive. 

# load required packages
library(dplyr)
library(readr)
library(googledrive)
library(ggplot2)
library(tidyr)
library(lubridate)

