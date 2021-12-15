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
## all dataframes from 2019, 2020, 2021 need to be cleaned of any trees that were girdled in early 2019
d2019$Filename <- as.character(d2019$Filename)

Asat2019 <- d2019 %>%
  group_by(Filename) %>%
  mutate(MeanPhoto = mean(Photo)) %>%
  filter(MeanPhoto >= 0,
         girdled. == "n")

Asat2020 <- d2020 %>%
  group_by(Filename) %>%
  mutate(MeanPhoto = mean(Photo)) %>%
  filter(MeanPhoto >= 0,
         girdled. == "n")

Asat2021 <- d2021 %>%
  group_by(Filename) %>%
  mutate(MeanPhoto = mean(Photo)) %>%
  filter(MeanPhoto >= 0,
         girdled. == "n")

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
treatments1 <- treatments %>%
  select(Subplot, disturbance_severity, treatment)

df <- left_join(alldata, treatments1, by = "Subplot")

## REMEMBER that you have to SCALE for your Asat values for pines. This means you need a new, final_Asat column
df <- df %>%
  mutate(Scaled_Asat = MeanPhoto*Asat_multiplier)

# check the range of values -- some are definitely going to be outliers
# this is mostly due to inaccuracies in pine needle area measurement and subsequent scaling for gas exchange values measured by IRGA
range(df$Scaled_Asat, na.rm=TRUE)

na <- which(is.na(df$Scaled_Asat))
tibble::view(na)
tibble::view(df)

# check for number of values in excess of 25 umol/m2/sec. This would be a very high value for a full-sun aspen, maybe a good benchmark for maximum possible rate

length(df$Scaled_Asat[df$Scaled_Asat > 25]) # 40 measurements (8 trees over time) with values in excess of 25

outliers <- df %>%
  filter(Scaled_Asat > 18) %>%
  select(Obs, Photo, Filename, Species, girdled., Year, IRGA_covered_area, Asat_multiplier, MeanPhoto, Scaled_Asat)

unique(outliers$Species) # it is only pire and pist that have extremely high values in the data set
#amel also has high values but I have reason to believe those are ecologically valid/meaningful
#find out what amel values are

amel <- df %>%
  filter(Species == "amel") %>%
  group_by(Filename) %>%
  summarize(MeanPhoto = MeanPhoto)

tibble::view(amel)

# highest value for amel was in 2020, 17.34 umol/m2/sec. So, 18 seems reasonable threshold.

# note that there is one unknown tree (mistaken ID in 2019 filename, very likely) with NA's that needs to be removed
# Tag # 61, 2019 data set

#### Last stage of data cleaning:
# 1. remove values of Scaled_Asat > 18
# 2. remove NAs


df1 <- df %>%
  na.omit() %>%
  filter(Scaled_Asat < 18)


## visualize overall samples
# note, this is not very helpful; need to get to subplot scale first as experimental unit
library(ggplot2)
# Basic density
p1 <- ggplot(df1, aes(x=Scaled_Asat)) + 
  geom_density()

p2 <- ggplot(df1, aes(x=Scaled_Asat)) +
  geom_histogram(aes(y=..density..)) +
  stat_function(fun = dnorm,
                args = list(mean = mean(df1$Scaled_Asat),
                            sd = sd(df1$Scaled_Asat)),
                col = "#1b98e0",
                size = 3)

###################################################################
## getting data summarized by subplot, year

subplot_mean_phys <- df1 %>%
  select(Cond, Filename, Species, Subplot, final_ID, Year, disturbance_severity, treatment, Asat_multiplier, Scaled_Asat) %>%
  group_by(Year, Subplot) %>%
  summarize(Subplot_mean_Asat = mean(Scaled_Asat))

subplot_mean_phys <- left_join(subplot_mean_phys, treatments1, by = "Subplot")

# add replicate to the dataset
subplot_mean_phys$Replicate <-
  ifelse(grepl("A...", subplot_mean_phys$Subplot),
         "A",
         ifelse(grepl("B...", subplot_mean_phys$Subplot),
                "B",
                ifelse(grepl("C...", subplot_mean_phys$Subplot),
                       "C",
                       ifelse(grepl("D...", subplot_mean_phys$Subplot),
                              "D", "Other"
                       ))))


## making boxplots
## make severity a factor!!!!
subplot_mean_phys$disturbance_severity <- as.factor(subplot_mean_phys$disturbance_severity)

bp1 <- subplot_mean_phys %>%
  ggplot(aes(x = disturbance_severity, y = Subplot_mean_Asat, group_by(disturbance_severity), fill = disturbance_severity)) + 
  theme_classic(base_size = 13) + 
  geom_boxplot(show.legend = FALSE) +
  xlab("Severity (% Gross Defoliation)") +
  ylab('Subcanopy CWM leaf photosynthetic rate (' *mu~ 'mol' ~CO[2]~ m^-2~s^-1*')') + 
  scale_fill_manual(values = c("#000000", "#009E73", "#0072B2", "#D55E00")) + 
  facet_grid(~Year)

bp1

################ Analysis 02.1 #################################################
##### Split-split plot ANOVA

#### First of two model runs

## using code from Max Grigri; same model run by him and also by Kayla Mathes on their summer 2019 data sets (note that this has been superseded by a new model version for Kayla's work)
library(agricolae)

## check normality of experimental unit-level Asat data
shapiro.test(subplot_mean_phys$Subplot_mean_Asat) # p = 0.18, looks normal!!

# clean up my all_data df for model run
subplot_mean_phys$disturbance_severity <- as.factor(subplot_mean_phys$disturbance_severity)
subplot_mean_phys$treatment <- as.factor(subplot_mean_phys$treatment)
subplot_mean_phys$Replicate <- as.factor(subplot_mean_phys$Replicate)
subplot_mean_phys$Year <- as.factor(subplot_mean_phys$Year)
# recode 0 to 0.00 to help me when with the post-hoc output
subplot_mean_phys$disturbance_severity <- recode_factor(subplot_mean_phys$disturbance_severity, "0" = "0.00")

### from here, need to work on it and haven't run it yet
model_run <- subplot_mean_phys %>%
  ungroup() %>%
  select(Replicate, disturbance_severity, treatment, Year, Subplot_mean_Asat) %>%
  filter(!is.nan(Subplot_mean_Asat))
# overall anova split-split plot model run
Asatmodel1 <- with(model_run, ssp.plot(Replicate, disturbance_severity, treatment, Year,
                                      Subplot_mean_Asat))
# post hoc analysis
# creating an object for each part of the model
gla <- Asatmodel1$gl.a
glb <- Asatmodel1$gl.b
glc <- Asatmodel1$gl.c
# doing the same as above for each error term
Ea <- Asatmodel1$Ea
Eb <- Asatmodel1$Eb
Ec <- Asatmodel1$Ec
# running an LSD post hoc test: Year
# LSD_output <- with(model_run, LSD.test(Subplot_mean_Asat, Year, glb, Eb,
#                                        console = TRUE))
# # LSD post hoc: Severity-Year
# LSD_output <- with(model_run, LSD.test(mean_Asat, Severity:Year, glc, Ec,
#                                        console = TRUE))
# 

##### Second of two model runs: use AOV, follow same format Kayla used in later model run

Asatmodel2 <- aov(Subplot_mean_Asat ~ disturbance_severity*treatment*Year + Error(Replicate/disturbance_severity/treatment/Year), data = subplot_mean_phys)

summary(Asatmodel2)

## post-hoc analysis
with(subplot_mean_phys, LSD.test(Subplot_mean_Asat, disturbance_severity:Year,72,1.47, console = TRUE))
# with(data, LSD.test(y, x, dferror, MSerror, console = TRUE))

## total agreement between both models! That's a good sign...

################ Analysis 02.2 #################################################
## Part 1: same analysis as above, but using VAI as outcome variable. Need VAI in all years
# 2018 - 2020 in fortedata, 2021 in preliminary data set from Kerstin N.
library(fortedata)

VAI_subplot_3yrs <- fortedata::fd_canopy_structure() %>%
  select(subplot_id, replicate, year, plot, subplot, vai_mean) %>%
  group_by(year, replicate, plot, subplot) %>%
  summarize(mean_vai = mean(vai_mean)) %>%
  mutate(Zero = "0")

unique(VAI_subplot_3yrs$year) #just checking -- yes, this data set is currently missing 2021 data

VAI_subplot_3yrs <- VAI_subplot_3yrs %>% 
  unite(Subplot, c("replicate", "Zero", "plot", "subplot"), sep="")

# now bring in 2021 data from Kerstin N
VAI_subplot_2021 <- read.csv("data/umbs_pcl_2021.csv") %>%
  select(Subplot, vai.mean) %>%
  mutate(year = 2021) %>%
  group_by(Subplot, year) %>%
  summarize(mean_vai = mean(vai.mean))

# now join tables
VAI_all_years <- rbind(VAI_subplot_3yrs, VAI_subplot_2021) %>%
  rename(Year = year) 

# now join to Asat data
# first refactor Year in VAI data
VAI_all_years$Year <- as.factor(VAI_all_years$Year)
Asat_VAI <- left_join(subplot_mean_phys, VAI_all_years, by = c("Subplot", "Year")) 

# look at trends in VAI over time
bp2 <- Asat_VAI %>%
  ggplot(aes(x = disturbance_severity, y = mean_vai, group_by(disturbance_severity), fill = disturbance_severity)) + 
  theme_classic(base_size = 13) + 
  geom_boxplot(show.legend = FALSE) +
  xlab("% Gross Defoliation") +
  ylab("Canopy VAI (unitless)") + 
  scale_fill_manual(values = c("#000000", "#009E73", "#0072B2", "#D55E00")) + 
  facet_grid(~Year)

bp2 

## run model for VAI, same form as 02.1 model
VAImodel1 <- aov(mean_vai ~ disturbance_severity*treatment*Year + Error(Replicate/disturbance_severity/treatment/Year), data = Asat_VAI)

summary(VAImodel1) ## This will not work because we are missing observations in some subplots on some years. We need to instead use a linear mixed model
#what we have is a split-split plot design with unbalanced repeated measures
library(lme4)
library(lmerTest) # this library will enable generation of p-values for effects, among other things

VAImodel2 <- lmer(mean_vai ~ disturbance_severity*Year + (1|Replicate), data = Asat_VAI)
summary(VAImodel2)

## running into issues...going to post question on Stack Exchange. Determine some info about data to share:
VAI_means <- Asat_VAI %>%
  na.omit() %>%
  group_by(Year) %>%
  summarize(mean(mean_vai),
            sd = sd(mean_vai))

######### code for Stack Overflow
# collected_vai <- rnorm(125, mean = 6, sd = 1)
# missing <- rep("NA", times = 3)
# all_vai <- c(collected_vai, missing)
# year1 <- rep(2018, times = 32)
# year2 <- rep(2019, times = 32)
# year3 <- rep(2020, times = 32)
# year4 <- rep(2021, times = 32)
# year <- c(year1, year2, year3, year4)
# disturbance_severity <- rep(c(0,45,65,85), each = 32)
# treatment <- rep(c("B" , "T"), each = 64)
# replicate <- rep(c("A", "B", "C", "D"), each = 32)
# data = cbind(all_vai, year, disturbance_severity, treatment, replicate)
# data <- as.data.frame(data)
# data$year <- as.factor(data$year)
# data$disturbance_severity <- as.factor(data$disturbance_severity)
# data$treatment <- as.factor(data$treatment)
# data$replicate <- as.factor(data$replicate)

####### Part 2 of 02.2
lmeAsat <- lmer(Subplot_mean_Asat ~ mean_vai*Year + (1|Year), data = Asat_VAI)
summary(lmeAsat) ## there is a significant interaction between disturbance 85% and year 2021

Asat_VAI_no2018 <- Asat_VAI %>%
  filter(Year != 2018)

lmeAsat_no2018 <- lmer(Subplot_mean_Asat ~ mean_vai*Year + (1|Year), data = Asat_VAI_no2018)
summary(lmeAsat_no2018) ## there is a significant interaction between disturbance 85% and year 2021


## graph, Asat vs. VAI
Asat_VAI_only2018 <- Asat_VAI %>%
  filter(Year == 2018)

p1 <- Asat_VAI %>%
  ggplot(aes(x = mean_vai, y = Subplot_mean_Asat)) + 
  theme_classic(base_size = 15) + 
  geom_point(show.legend = FALSE, size=3) +
  xlab("Canopy VAI (unitless)") +
  ylab('Subcanopy CWM leaf photosynthetic rate (' *mu~ 'mol' ~CO[2]~ m^-2~s^-1*')') +
  geom_point(data=Asat_VAI_only2018, 
             aes(x=mean_vai,y=Subplot_mean_Asat), 
             color='red',
             size=3)

p1


p11 <- Asat_VAI_no2018 %>%
  ggplot(aes(x = mean_vai, y = Subplot_mean_Asat)) + 
  theme_classic(base_size = 15) + 
  geom_point(show.legend = FALSE, size=3) +
  xlab("Canopy VAI (unitless)") +
  ylab('Subcanopy CWM leaf photosynthetic rate (' *mu~ 'mol' ~CO[2]~ m^-2~s^-1*')') +
  geom_smooth(method = "lm")

p11

fit1 <- lm(Subplot_mean_Asat ~ mean_vai, data = Asat_VAI_no2018)
summary(fit1)

## plotting relationships from LME is a bit more challenging; call new library and try to add effect estimates to a data frame
library(effects)
ef <- effect("mean_vai:Year", lmeAsat)
x <- as.data.frame(ef)
x

# another way to add to data frame
Asat_VAI_noNAs <- Asat_VAI %>%
  na.omit()

Asat_VAI_noNAs$lmeAsatfit <- predict(lmeAsat)

# try to visualize
p2 <- ggplot(Asat_VAI_noNAs,aes(mean_vai, Subplot_mean_Asat, group=Year, col=Year)) + 
  geom_line(aes(y=lmeAsatfit, lty=Year), size=1) +
  geom_point(size = 2.5) + 
  theme_classic(base_size=13) +
  xlab("Canopy VAI (unitless)") +
  ylab('Subcanopy CWM leaf photosynthetic rate (' *mu~ 'mol' ~CO[2]~ m^-2~s^-1*')') 
 
p2

###### Try Asat vs. VAI analysis using no-NA's data and an anova


################ Analysis 02.3 #################################################
## bring in subcanopy NPP data

Asat_npp <- read.csv("data/cleaned_data_for_analyses/subcan_NPP_phys.csv")
head(Asat_npp)

qqnorm(Asat_npp$Subcanopy_NPP)
qqline(Asat_npp$Subcanopy_NPP) ## npp data are very non-normal
hist(Asat_npp$Subcanopy_NPP)

Asat_npp$Year <- as.factor(Asat_npp$Year)
Asat_npp$disturbance_severity <- as.factor(Asat_npp$disturbance_severity)

as_tibble(Asat_npp)

## GLMM for NPP vs Asat 
npp_model <- glmer(Subcanopy_NPP ~ Subplot_mean_Asat*disturbance_severity + disturbance_severity*Year + (1|Year), data = Asat_npp, family = Gamma)
summary(npp_model)
plot(allEffects(npp_model)) ## right now, this is not making sense to me...

## logged NPP is response variable!
lmeNPP <- lmer(log(Subcanopy_NPP) ~ Subplot_mean_Asat*disturbance_severity*Year + (1|Year), data = Asat_npp)
summary(lmeNPP)

# try to visualize
p3 <- ggplot(Asat_npp,aes(Subplot_mean_Asat, log(Subcanopy_NPP), group=Year, col=Year)) + 
  geom_point(size = 4) + 
  theme_classic(base_size=18) +
  xlab('Subcanopy CWM leaf photosynthetic rate (' *mu~ 'mol' ~CO[2]~ m^-2~s^-1*')') +
  ylab("log(Subcanopy NPP) kg C *m^-2* yr") 

p3

p4 <- ggplot(Asat_npp,aes(Subplot_mean_Asat, Subcanopy_NPP, group=Year, col=Year)) + 
  geom_point(size = 4) + 
  theme_classic(base_size=18) +
  xlab('Subcanopy CWM leaf photosynthetic rate (' *mu~ 'mol' ~CO[2]~ m^-2~s^-1*')') +
  ylab("Subcanopy NPP kg C *m^-2* yr") 

p4

## checking NPP: running ANOVA on 3 year NPP
npp <- read.csv("data/cleaned_data_for_analyses/three_yr_sc_NPP.csv")
head(npp)

as_tibble(npp)

## we know NPP data are not-normal. test degree of skewness using "moments" package
library(moments)
skewness(npp$NPP)

## now transform vector of NPP data

npp$severity <- as.factor(npp$severity)
npp$treatment <- as.factor(npp$treatment)
npp$year <- as.factor(npp$year)
npp$replicate <- as.factor(npp$replicate)
npp$log_NPP <- log(npp$NPP) # log transforming NPP data to deal with skewness

## Try the original ANOVA model on log-transformed NPP data
NPPmodel_aov <- aov(log_NPP ~ severity*treatment*year + Error(replicate/severity/treatment/year), data = npp)

summary(NPPmodel_aov)


## run posthoc, severity:year
with(npp, LSD.test(log_NPP, severity:year,48,0.286, console = TRUE))
# with(data, LSD.test(y, x, dferror, MSerror, console = TRUE))

## visualize the subcanopy NPP data
bp3 <- npp %>%
  ggplot(aes(x = severity, y = NPP, group_by(severity), fill = severity)) + 
  theme_classic(base_size = 20) + 
  geom_boxplot(show.legend = FALSE) +
  xlab("Severity (% Gross Defoliation)") +
  ylab(expression(paste("ANPP" [w], " ( ",kgC," ",ha^-1," ",year^-1,")"))) + 
  scale_fill_manual(values = c("#000000", "#009E73", "#0072B2", "#D55E00")) + 
  facet_grid(~year)

bp3 
