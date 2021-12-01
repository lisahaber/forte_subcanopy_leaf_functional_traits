###########################################################################
## Lisa T. Haber                                                         ##
## November 2021                                                         ##
## FoRTE Subcanopy Leaf Functional Traits Manuscript Prep                ##
## 2018-2021 data analysis                                               ##
###########################################################################

## This code will use cleaned FoRTE subcanopy physiology data for preliminary analysis.
## The immediate goal is preparation for AGU 2021.
## The principle analyses here will be ANOVA/mixed models..... refer to proposal!!!

## load required packages
library(tidyverse)
library(lubridate)

## bring in the cleaned phys data
d2018 <- read.csv("data/cleaned_data_for_analyses/LICOR_leaf_phys_2018.csv") %>%
  as_tibble()
d2019 <- read.csv("data/cleaned_data_for_analyses/LICOR_leaf_phys_2019.csv") %>%
  as_tibble() 
d2020 <- read.csv("data/cleaned_data_for_analyses/LICOR_leaf_phys_2020.csv") %>%
  as_tibble()
d2021 <- read.csv("data/cleaned_data_for_analyses/LICOR_leaf_phys_2021.csv") %>%
  as_tibble()

## get rid of negative Asat values and create means for individual leaves from 5 observations per leaf
Asat2018 <- d2018 %>%
  group_by(Filename) %>%
  mutate(MeanPhoto = mean(Photo)) %>%
  filter(MeanPhoto >= 0)

# Filename column in 2019 data is a numeric vector rather than character, so need to convert
d2019$Filename <- as.character(d2019$Filename)

Asat2019 <- d2019 %>%
  group_by(Filename) %>%
  mutate(MeanPhoto = mean(Photo)) %>%
  filter(MeanPhoto >= 0)

Asat2020 <- d2020 %>%
  group_by(Filename) %>%
  mutate(MeanPhoto = mean(Photo)) %>%
  filter(MeanPhoto >= 0)

Asat2021 <- d2021 %>%
  group_by(Filename) %>%
  mutate(MeanPhoto = mean(Photo)) %>%
  filter(MeanPhoto >= 0)

## bring in the FoRTE subplots treatment assignments
library(readr)
urlfile="https://raw.githubusercontent.com/lisahaber/fortedata/master/inst/extdata/fd_subplots.csv" #github link for raw fortedata file with treatment and severity assignments
treatments <- read_csv(url(urlfile))

head(treatments)

## create a new column for the complete subplot identifier, so you can join to cleaned phys data
require(tidyverse)

treatments <- treatments %>%
  mutate(Zero = "0")

treatments <- treatments %>% 
  unite(Subplot, c("replicate", "Zero", "plot", "subplot"), sep="")

## merge all years' data
alldata <- rbind(Asat2018, Asat2019, Asat2020, Asat2021)

## now join to treatments data
treatments <- treatments %>%
  select(Subplot, disturbance_severity, treatment)

df <- left_join(alldata, treatments, by = "Subplot")

## REMEMBER that you have to SCALE for your Asat values for pines. This means you need a new, final_Asat column
df <- df %>%
  mutate(Scaled_Asat = MeanPhoto*Asat_multiplier)

# check the range of values -- some are definitely going to be outliers
range(df$Scaled_Asat, na.rm=TRUE)

na <- which(is.na(df$Scaled_Asat))
tibble::view(na)
tibble::view(df)

# check for number of values in excess of 25 umol/m2/sec. This would be a very high value for a full-sun aspen, maybe a good benchmark for maximum possible rate

length(df$Scaled_Asat[df$Scaled_Asat > 25])

outliers <- df %>%
  filter(Scaled_Asat > 18) %>%
  select(Obs, Photo, Filename, Species, girdled., Year, IRGA_covered_area, Asat_multiplier, MeanPhoto, Scaled_Asat)

unique(outliers$Species) # it is only pire and pist that have extremely high values in the data set
#amel also has high values but I have reason to believe those are fine 
#find out what amel values are

amel <- df %>%
  filter(Species == "amel") %>%
  group_by(Filename) %>%
  summarize(MeanPhoto = MeanPhoto)

tibble::view(amel)

# highest value for amel was in 2020, 17.34 umol/m2/sec. So, 18 seems reasonable threshold.

# note that there is one unknown tree (mistaken ID in 2019 filename, very likely) with NA's that needs to be removed
# ID 61, 2019 data set

## visualize initial trends
# check for normality
library(ggplot2)
# Basic density
p1 <- ggplot(df, aes(x=Scaled_Asat)) + 
  geom_density()

p2 <- ggplot(df, aes(x=Scaled_Asat)) +
  geom_histogram()



shapiro.test(df$MeanPhoto)


