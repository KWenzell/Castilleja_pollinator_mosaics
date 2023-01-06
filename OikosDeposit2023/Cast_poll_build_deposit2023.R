#Pollinator visitation to Castilleja sessiliflora and the C. purpurea complex
#Observations recorded 2017-2019

#BUILD- data import and cleaning
library(tidyverse)
library(vegan)
library(glmmTMB)
library(lme4)
library(performance)
library(reshape2)

library(dplyr)
library(ggplot2)
library(ggpubr)

#Focal dataset 2017-2019
setwd("~/CASE_Rdocs/Cast_pollinators/OikosDeposit2023/build/")
pollObs2 <- read.csv("PollObs_rawFocal.csv")

str(pollObs2)

#summarise raw data 
#####

#fill in number of open flowers for all focal patches per population-year
poOpenflrs <- group_by(pollObs2, Year, Pop) %>%
  mutate(OpenFlrs_patches = sum(unique(OpenFlrs_patches)) )

#calculate total number of flowers visited by each pollinator group per population-year
po <- group_by(poOpenflrs, Year, Pop, PollFxGroup) %>%
  summarise(pgVisit = sum(NoFlowersVisited), openFlrs = mean(OpenFlrs_patches) ) #mean here just gives same value

#calculate visitation rate (number of visits relative to number open flowers) and add as column 
po <- po %>% 
  mutate(visRate = pgVisit / openFlrs, popYr = paste(Pop, Year, sep = "-"))

#calculate rate per hour
po <- po %>% 
  mutate(visRateH = visRate / 7.33) #440 minutes per site/ 60= 7.33 hours per site

#repeat for nonfocal

#calculate proportion of visits 
#first, calculate total visits per pop-year
po <- group_by(po, popYr) %>%
  mutate(visTotal = sum(pgVisit))
#calculate proportion visit for poll groups 
po <- po %>%
  mutate(visProp = pgVisit / visTotal)

#replace NA and NaN with 0 (product of no_visitation sites visTotal=0, so visProp= NaN; also openFlrs= NA when no_visitation, so visRate= NA)
po[is.na(po)] <- 0

#Nonfocal dataset- 2019
pollObsNonfocal <- read.csv("PollObs_rawNonfocal.csv")


#fill in number of estimated available flowers for nonfocal observations per population-year
nfAvailflrs <- group_by(pollObsNonfocal, year, pop) %>%
  mutate(estTotalFlrs = sum(unique(estimatedTotalFlrs)) )

#calculate total number of estimated flowers visited by each pollinator group per population-year
nf <- group_by(nfAvailflrs, year, pop, pollFxGroup) %>%
  summarise(pgVisit = sum(totalNumberFlrsVisited), estOpenFlrs = mean(estTotalFlrs), dataset = "nf" )

#calculate visitation rate (number of visits relative to number estimated open flowers) and add as column 
nf <- nf %>% 
  mutate(visRate = pgVisit / estOpenFlrs, popYr = paste(pop, year, dataset, sep = "-"))

#calculate rate per hour
nf <- nf %>% 
  mutate(visRateH = visRate / 7.33) #440 minutes per site/ 60= 7.33 hours per site

#calculate NF proportion of visits 
#first, calculate NF total visits per pop-year
nf <- group_by(nf, popYr) %>%
  mutate(visTotal = sum(pgVisit))
#calculate NF proportion visit for poll groups 
nf <- nf %>%
  mutate(visProp = pgVisit / visTotal)

#replace NAs from no_visitation pops with 0's
nf[is.na(nf)] <- 0


#in dataset, there are 15 unique pollinator groups, need to combine some
#compare # visits by poll group
aggregate(pgVisit ~ PollFxGroup, po, sum)

#see # entries per poll group (# observed per popYr)
fct_count(po$PollFxGroup)


#combine functionally similar groups (where one group is infrequent (<15 visits)):
po$pollFxGrp2 <- fct_collapse(po$PollFxGroup,
                              bumble_lg_bee = c("bumblebee", "large_bee"),
                              butterfly_all = c("butterfly", "small_butterfly"),
                              small_med_bee = c("small_bee", "med_bee")
                              
)

#compare # visits by poll group
aggregate(pgVisit ~ pollFxGrp2, po, sum)

#see # entries per poll group (# observed per popYr)
fct_count(po$pollFxGrp2)

#now, combine all groups with <50 visits OR at <3 site-years* into "other" category:
#both carpenter bee and honeybee were recorded but 0 visits (carpenter bees at PTMS were nectar robbing)
#include 'other' with low frequency (<50 visits) groups: hoverfly, wasp_other, small_moth, and carpenter bees, honeybee, and no_visitation
#include in "other" bee-fly, which was only present at 1 population (3pop-Yrs in dataset, but were only observed at SIC)
po$pollFxGrp2 <- fct_collapse(po$pollFxGrp2,
                              other = c("carpenter_bee", "honeybee", "hoverfly", "no_visitation", 
                                        "small_moth", "wasp_other", "bee-fly")
                              
)

#compare # visits by poll group
aggregate(pgVisit ~ pollFxGrp2, po, sum)

#see # entries per poll group (# observed per popYr)
fct_count(po$pollFxGrp2)

#now pollFxGrp2 has 6 levels!



#repeat with nonfocal dataset, using same combinations:
#compare # visits by poll group
aggregate(pgVisit ~ pollFxGroup, nf, sum)
#see # entries per poll group (# observed per popYr)
fct_count(nf$pollFxGroup)


#combine functionally similar groups and "other" for infrequent categories:
nf$pollFxGrp2 <- fct_collapse(nf$pollFxGroup,
                              #bumble_lg_bee = c("bumblebee", "large_bee"), #no "large_bee" in nf 
                              butterfly_all = c("butterfly", "small_butterfly"),
                              #small_med_bee = c("small_bee", "med_bee") #no "med_bee" in nf
                              other = c("carpenter_bee", "no_visitation", "small_moth")
                              
)

#compare # visits by poll group
aggregate(pgVisit ~ pollFxGrp2, nf, sum)
#see # entries per poll group (# observed per popYr)
fct_count(nf$pollFxGrp2)



#pivot from long to wide dataframe, separate count, rate, and proportion data
#remove duplicates and sum values from collapsed factors
#FOCAL COUNT DATA
poCt <- group_by(po, popYr, Year, Pop, pollFxGrp2) %>%
  summarise(pgVisit = sum(pgVisit))

fCt <- poCt %>% 
  pivot_wider(names_from = pollFxGrp2, names_prefix = "ct_", values_from = pgVisit, values_fill=list(pgVisit=0))


#FOCAL PROPORTION DATA
poPr <- group_by(po, popYr, Year, Pop, pollFxGrp2) %>%
  summarise(visProp = sum(visProp))

fPr <- poPr %>% 
  pivot_wider(names_from = pollFxGrp2, names_prefix = "pr_", values_from = visProp, values_fill=list(visProp=0))


#FOCAL RATE DATA
poRt <- group_by(po, popYr, Year, Pop, pollFxGrp2) %>%
  summarise(visRate = sum(visRate))

fRt <- poRt %>% 
  pivot_wider(names_from = pollFxGrp2, names_prefix = "rt_", values_from = visRate, values_fill=list(visRate=0))


#REPEAT WITH NONFOCAL DATA
#NONFOCAL COUNT DATA
nfCt_l <- group_by(nf, popYr, year, pop, pollFxGrp2) %>%
  summarise(pgVisit = sum(pgVisit))

nfCt <- nfCt_l %>% 
  pivot_wider(names_from = pollFxGrp2, names_prefix = "ct_", values_from = pgVisit, values_fill=list(pgVisit=0))


#NONFOCAL PROPORTION DATA
nfPr_l <- group_by(nf, popYr, year, pop, pollFxGrp2) %>%
  summarise(visProp = sum(visProp))

nfPr <- nfPr_l %>% 
  pivot_wider(names_from = pollFxGrp2, names_prefix = "pr_", values_from = visProp, values_fill=list(visProp=0))

#NONFOCAL RATE DATA
nfRt_l <- group_by(nf, popYr, year, pop, pollFxGrp2) %>%
  summarise(visRate = sum(visRate))

nfRt <- nfRt_l %>% 
  pivot_wider(names_from = pollFxGrp2, names_prefix = "rt_", values_from = visRate, values_fill=list(visRate=0))

#link pollinator data with floral data
ca <- read.csv("CastFloral_v2-coreTraits_corrected_forPoll.csv")
#summarise by population
caPop <- group_by(ca, pop, taxa, sp, reg1, reg2) %>%
  summarise(corL = mean(corL), corW = mean(corW2), stigma = mean(stigma), lipL = mean(lipL), brLbW = mean(brLbW),
            red = mean(red), green = mean(green), blue = mean(blue), latN = mean(LatN), lonW = mean(LonW))

#join dataframes

#combine count, prop, rate back into one dataframe
fw <- left_join(fCt, fPr, by = c("popYr", "Year", "Pop"))
fw <- left_join(fw, fRt, by = c("popYr", "Year", "Pop"))

#join floral data
fwf <- left_join(fw, caPop, by = c("Pop"="pop"))

#add total number of visits back in, to weight binomial regression
fwf <- fwf  %>% mutate(visTotal = sum(c(ct_bumble_lg_bee, ct_butterfly_all, ct_other, ct_hawkmoth, ct_hummingbird, ct_small_med_bee)))


#check to see values repeated for repeated pops:
select(fwf, popYr, Pop, corL, latN)

#repeat above for nonfocal
#NONFOCAL: combine count, prop, rate back into one dataframe
nfw <- left_join(nfCt, nfPr, by = c("popYr", "year", "pop"))
nfw <- left_join(nfw, nfRt, by = c("popYr", "year", "pop"))

#join floral data to nonfocal
nfwf <- left_join(nfw, caPop, by = "pop")
#add total number of visits back in, to weight binomial regression
nfwf <- nfwf  %>% mutate(visTotal = sum(c(ct_bumblebee, ct_butterfly_all, ct_other, ct_hawkmoth, ct_hummingbird, ct_small_bee)))


#rename dataframes
fp <- fwf
nfp <- nfwf

#combine focal and nonfocal datasets

#first, add variable to each specifying how data was collected
#for focal> "focal_narrow"
fp <- fp %>% 
  mutate(dataset = paste("focal_narrow"))
# #remove junk columns (col#s added by R)
# fp <- fp[,2:41]
#rename column names for fp to match corresponding nfp columns (and remove annoying capitals)
fp <- fp %>%
  rename(year= Year,
         pop= Pop)

fp <- fp %>%
  relocate(dataset, .after = pop)

#for nonfocal/ wide-view: "nonfocal_wide"
nfp <- nfp %>% 
  mutate(dataset = paste("nonfocal_wide"))

# #remove junk columns (col#s added by R)
# nfp <- nfp[,2:43]
nfp <- nfp %>%
  relocate(dataset, .after = pop)

#rename column names for nfp to match corresponding fp columns
nfp <- nfp %>%
  rename(ct_bumble_lg_bee= ct_bumblebee,
         ct_small_med_bee= ct_small_bee,
         pr_bumble_lg_bee= pr_bumblebee,
         pr_small_med_bee= pr_small_bee,
         rt_bumble_lg_bee= rt_bumblebee,
         rt_small_med_bee= rt_small_bee)
         

#now, append nfp onto fp dataframe
tp <- rbind(fp,nfp)
tp$year <- as.factor(tp$year)
str(tp)

#export
write.csv(tp, "Cast_totalPollFlrFrt_dep.csv")

