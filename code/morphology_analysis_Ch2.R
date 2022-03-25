###########################################################################
## Lisa T. Haber                                                         ##
## December 2021                                                         ##
## FoRTE Subcanopy Leaf Functional Traits Manuscript Prep                ##
## 2018-2021 LMA data analysis                                           ##
###########################################################################

## Working up LMA data for analysis

## load required packages
library(tidyverse)

## load data, 2018-2021
morph2018 <- read.csv("data/leaf_morphology/FoRTE_subcanopy_leaf_morphology_2018.csv", header = TRUE) %>%
  as_tibble()

morph2018 <- morph2018 %>%
  mutate(Year = 2018) %>%
  rename(final_ID = tag_number) %>%
  select(final_ID, Subplot, Species, leaf_area_cm2, leaf_mass_g, LMA_gm2, Year)

morph2019 <- read.csv("data/leaf_morphology/FoRTE_subcanopy_leaf_morphology_2019.csv", header = TRUE) %>%
  as_tibble()

morph2019 <- morph2019 %>%
  mutate(Year = 2019) %>%
  rename(final_ID = tag_number)

morph2020 <- read.csv("data/leaf_morphology/FoRTE_subcanopy_leaf_morphology_2020.csv", header = TRUE) %>%
  as_tibble()

morph2020 <- morph2020 %>%
  mutate(Year = 2020) %>%
  rename(Subplot = subplot_id,
         final_ID = tag_number)

morph2021 <- read.csv("data/leaf_morphology/FoRTE_subcanopy_leaf_morphology_2021.csv", header = TRUE) %>%
  as_tibble()

morph2021 <- morph2021 %>%
  mutate(Year = 2021) %>%
  rename(Subplot = subplot_id,
         final_ID = tag_number)

## add columns for Species, Subplot to years where these variables are missing; clean and tidy the data for binding
ID_lookup <- read.csv("data/forte_subcanopy_tree_id_master_lookup.csv") %>%
  as_tibble()

#2019
morph2019 <- morph2019 %>%
  left_join(ID_lookup, by = "final_ID") %>%
  select(final_ID, Subplot, Species, leaf_area_cm2, leaf_mass_g, LMA_gm2, Year)

#2020
morph2020 <- morph2020 %>%
  left_join(ID_lookup, by = "final_ID") %>%
  select(final_ID, Subplot.y, Species, leaf_area_cm2, leaf_mass_g, LMA_gm2, Year) %>%
  rename(Subplot = Subplot.y)

#2021
morph2021 <- morph2021 %>%
  left_join(ID_lookup, by = "final_ID") %>%
  select(final_ID, Subplot.y, Species, leaf_area_cm2, leaf_mass_g, LMA_gm2, Year) %>%
  rename(Subplot = Subplot.y)

allmorph <- rbind(morph2018, morph2019, morph2020, morph2021) %>%
  as_tibble() %>%
  na.omit()

## visualize overall samples
# note, this is not very helpful; need to get to subplot scale first as experimental unit
library(ggplot2)
# Basic density
p1 <- ggplot(allmorph, aes(x=LMA_gm2)) + 
  geom_density()
p1

p2 <- ggplot(allmorph, aes(x=LMA_gm2)) +
  geom_histogram(aes(y=..density..)) +
  stat_function(fun = dnorm,
                args = list(mean = mean(allmorph$LMA_gm2),
                            sd = sd(allmorph$LMA_gm2)),
                col = "#1b98e0",
                size = 3)
p2

## find NA values in dataset
which(is.na(allmorph$LMA_gm2))

range(allmorph$LMA_gm2) # lots of very large values for pines, need to weed out crazy pine values

length(allmorph$LMA_gm2[allmorph$LMA_gm2 > 150])
# 59 samples with LMA over 150; 19 with LMA over 200

allmorph_cleaned <- allmorph %>%
  filter(LMA_gm2 <= 150)

p3 <- ggplot(allmorph_cleaned, aes(x=LMA_gm2)) + 
  geom_density()
p3

counts <- allmorph_cleaned %>% group_by(Year, Subplot) %>% summarize(n = n())
tibble::view(counts)

write.csv(allmorph_cleaned, "data/cleaned_data_for_analyses/lma_sc_all_years.csv", row.names=FALSE)
