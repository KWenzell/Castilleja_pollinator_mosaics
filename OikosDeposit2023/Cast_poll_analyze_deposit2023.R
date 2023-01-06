#20220103
#Code to deposit for Wenzell et al., 2023, Oikos
#Compiling code for Castilleja Pollinator Mosaic paper
#from Cast_poll_submission_20220330_cleaningv1.R
library(tidyverse)
library(vegan)
library(glmmTMB)
library(lme4)
library(performance)
library(betareg)
library(car)
library(reshape2)
library(dplyr)
library(ggpubr)
library(lsmeans)
library(mvnormtest)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(DescTools)


#master dataframe for analysis of pollinator observations: 
# #multiyear data fp and nfp combined (with variable for dataset and variable for year):
setwd("~/CASE_Rdocs/Cast_pollinators/OikosDeposit2023/analyze/")
myp1 <- read.csv("Cast_totalPollFlrFrt_dep.csv")
#remove junk columns (col#s added by R)
myp1 <- myp1[,2:38]
#change year to factor from integer
myp1$year <- as.factor(myp1$year)
str(myp1)
#remove SNRS rows w/ NAs in floral trait data:
myp2 <- myp1[-c(38, 39), ]


#FLORAL TRAITS
#####

#read in individual-level floral trait data:
ca <- read.csv("CastFloral_v2-coreTraits_corrected.csv")
#subset for only populations with pollinator data
cap <- ca[-c(330:359, 536:565, 626:655), ]
unique(cap$pop) #successfully removed SBK, SMAR, SOZ
#use this dataframe for individual-level floral traits

#NMDS for colors- dimension reduction
mds <- metaMDS(cap[ 11:13], distance= "gower", zerodist="add")

mds
# Stress:     0.03418663 
#add NMDS values to dataframe
cap$nmds1_rgb <- mds$points[,1]
cap$nmds2_rgb <- mds$points[,2]
str(cap)

#sp symbols for ordination
cap$pch_ind[cap$sp == "P"] <- 22 #square
cap$pch_ind[cap$sp == "S"] <- 21 #circle
cap$pch_ind[cap$sp == "C"] <- 23 #diamond
cap$pch_ind[cap$sp == "L"] <- 24 #triangle

#use hexadecimals for color values

#library(stringr)
#columns with RGB values
rgb <- cap[,11:13]
#function rgb_to_hex: https://rdrr.io/github/stuart-morrison/schemr/man/rgb_to_hex.html
rgb_to_hex <- function(rgb) {

  # Extract colour channels from first three colummns
  red <- unlist(rgb[ , 1])
  green <- unlist(rgb[ , 2])
  blue <- unlist(rgb[ , 3])

  # Test if any RGB values are outside of [0, 255]
  if (any(red < 0 | red > 255) |
      any(green < 0 | green > 255) |
      any(blue < 0 | blue > 255)) {
    stop("Colour channels should be between 0 and 255.")
  }

  # Convert colour channels to hex
  hex_value <- paste0("#",
                      # Convert red
                      str_pad(string = as.hexmode(round(red, 0)),
                              width = 2, pad = "0", side = "left"),
                      # Convert green
                      str_pad(string = as.hexmode(round(green, 0)),
                              width = 2, pad = "0", side = "left"),
                      # Convert blue
                      str_pad(string = as.hexmode(round(blue, 0)),
                              width = 2, pad = "0", side = "left")
  )

  return(hex_value)
}
cap$hex <- rgb_to_hex(rgb)

#Visualize NMDS_colors only with hex colors
fig <- ordiplot(mds, type = "none")
points(fig, "sites", pch= cap$pch_ind, bg= cap$hex,
       cex=1.3)
legend("topleft",  inset= c(0.05,0),    # location and inset- still not great
       bty="n", cex= 1,
       title="Species", xpd= TRUE,
       c("C. sessiliflora", "C. purpurea", "C. lindheimeri", "C. citrina"),
       pch= c(16, 15, 17, 18),  pt.cex= 1.5,
       pt.bg= "black")
fit <- envfit(mds, cap[, 11:13], perm = 1000)
fit
plot(fit, p.max = 0.05, col = "black", cex=1.2)
#NMDS1: higher=blue (cool colors= purple, pink, green), lower= red (warm colors= red, orange, yellow)
#NMDS2: higher= red (value and color); red and purple; lower= green (value and color); green and yellow
#quadrants: (-, +) top left= L/red;; (+,+) top right= P/purple;; (-,-) bottom left= C(S)/ yellow;; (-,+) bottom right= S/green-pink
### =FIGURE S1


#####
#Test for multivariate variation in all floral traits among species

#MANOVA w/ morphological floral traits and NMDS color values (not RGB)

#species and pop nested w/in sp
#double checck columns include 5 morph traits and 2 NMDS color values
traits2 <- as.matrix(cbind(cap[,c(6,7,8,9,10,16,17)]))
sp <- cap[,4]
pop <- cap[,2]
fit2 <- manova(traits2 ~ sp/pop)
summary(fit2)

#FLORAL TRAIT X SPECIES VIOLIN PLOTS

cl <- ggplot(data = cap, aes(x= factor(sp, level= c("C", "P", "L", "S")), y= corL))+
  geom_violin(aes(fill= sp), trim= T, alpha= 0.8)+
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="black")+
  scale_fill_manual(values=c('S'='#b9d06a','P'='#b5266e', 'L'='#f13918', 'C'='#f8e432'))+
  labs(x= "Species", y= "Mean corolla length (mm)")+
  theme_classic()

cl 
#anova

#corolla length- variation among species and population nested within species
model <- aov(corL~ sp/pop, data=cap)
summary(model)
TukeyHSD(model, conf.level=.95)
#ignore $sp:pop comparisons (too many pairwise comps), only look at $sp comparisons

cw <- ggplot(data = cap, aes(x= factor(sp, level= c("C", "P", "L", "S")), y= corW2))+
  geom_violin(aes(fill= sp), trim= T, alpha= 0.8)+
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="black")+
  scale_fill_manual(values=c('S'='#b9d06a','P'='#b5266e', 'L'='#f13918', 'C'='#f8e432'))+
  labs(x= "Species", y= "Mean corolla width (mm)")+
  theme_classic()

cw 

#ANOVA
model <- aov(corW2~sp/pop, data=cap)
summary(model)
TukeyHSD(model, conf.level=.95)
#report $sp comparisons


ll <- ggplot(data = cap, aes(x= factor(sp, level= c("C", "P", "L", "S")), y= lipL))+
  geom_violin(aes(fill= sp), trim= T, alpha= 0.8)+
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="black")+
  scale_fill_manual(values=c('S'='#b9d06a','P'='#b5266e', 'L'='#f13918', 'C'='#f8e432'))+
  labs(x= "Species", y= "Mean lip length (mm)")+
  theme_classic()

ll 
#anova
model <- aov(lipL~sp/pop, data=cap)
summary(model)
TukeyHSD(model, conf.level=.95)
#report $sp comparisons

st <- ggplot(data = cap, aes(x= factor(sp, level= c("C", "P", "L", "S")), y= stigma))+
  geom_violin(aes(fill= sp), trim= T, alpha= 0.8)+
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="black")+
  scale_fill_manual(values=c('S'='#b9d06a','P'='#b5266e', 'L'='#f13918', 'C'='#f8e432'))+
  labs(x= "Species", y= "Mean stigma exsertion (mm)")+
  theme_classic()

st 

#anova
model <- aov(stigma~sp/pop, data=cap)
summary(model)
TukeyHSD(model, conf.level=.95)
#report $sp comparisons

bl <- ggplot(data = cap, aes(x= factor(sp, level= c("C", "P", "L", "S")), y= brLbW))+
  geom_violin(aes(fill= sp), trim= T, alpha= 0.8)+
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="black")+
  scale_fill_manual(values=c('S'='#b9d06a','P'='#b5266e', 'L'='#f13918', 'C'='#f8e432'))+
  labs(x= "Species", y= "Mean bract lobe width (mm)")+
  theme_classic()

bl 

#anova
model <- aov(brLbW~sp/pop, data=cap)
summary(model)
TukeyHSD(model, conf.level=.95)
#report $sp comparisons


#with NMDS floral color axes
n1 <- ggplot(data = cap, aes(x= factor(sp, level= c("C", "P", "L", "S")), y= nmds1_rgb))+
  geom_violin(aes(fill= sp), trim= T, alpha= 0.8)+
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="black")+
  scale_fill_manual(values=c('S'='#b9d06a','P'='#b5266e', 'L'='#f13918', 'C'='#f8e432'))+
  labs(x= "Species", y= "Color (NMDS1)")+
  theme_classic()
n1 

#anova
model <- aov(nmds1_rgb~sp/pop, data=cap)
summary(model)
TukeyHSD(model, conf.level=.95)
#report $sp comparisons

n2 <- ggplot(data = cap, aes(x= factor(sp, level= c("C", "P", "L", "S")), y= nmds2_rgb))+
  geom_violin(aes(fill= sp), trim= T, alpha= 0.8)+
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="black")+
  #scale_fill_manual(values=c('S'='olivedrab','P'='mediumorchid', 'L'='tomato', 'C'='darkgoldenrod1'))+
  scale_fill_manual(values=c('S'='#b9d06a','P'='#b5266e', 'L'='#f13918', 'C'='#f8e432'))+  
  labs(x= "Species", y= "Color (NMDS2)")+
  #coord_flip()+
  #stat_compare_means(label= "p.signif", ref.group = ".all.")+
  theme_classic()
n2 

#anova
model <- aov(nmds2_rgb~sp/pop, data=cap)
summary(model)
TukeyHSD(model, conf.level=.95)
#report $sp comparisons

#violin plots
ggarrange(cl, ll, st, bl, cw, n1, n2, nrow= 2, ncol= 4, common.legend = T)
#FIG S2
#contrast letters added manually (Inkscape)
#####

#add population average NMDS color values to pollinator dataframe
avgNMDS <- group_by(cap, pop) %>%
  summarise(avgNMDS1 = mean(nmds1_rgb), avgNMDS2 = mean(nmds2_rgb))

#join dataframes
#combine count, prop, rate back into one dataframe
myp2 <- left_join(myp2, avgNMDS, by = c("pop"))



#Figure 2- MAP WITH MEDIAN HEX COLORS AND PIE CHARTS
#####
#Summarize average floral colors by population:
#export myp2 with population summaries for pollinator vis, floral traits, etc
#add hexadecimal population mean color values

#try hex colors
#library(stringr)
rgb <- myp2[,32:34]
#function rgb_to_hex: https://rdrr.io/github/stuart-morrison/schemr/man/rgb_to_hex.html
rgb_to_hex <- function(rgb) {
  
  # Extract colour channels from first three colummns
  red <- unlist(rgb[ , 1])
  green <- unlist(rgb[ , 2])
  blue <- unlist(rgb[ , 3])
  
  # Test if any RGB values are outside of [0, 255]
  if (any(red < 0 | red > 255) |
      any(green < 0 | green > 255) |
      any(blue < 0 | blue > 255)) {
    stop("Colour channels should be between 0 and 255.")
  }
  
  # Convert colour channels to hex
  hex_value <- paste0("#",
                      # Convert red
                      str_pad(string = as.hexmode(round(red, 0)),
                              width = 2, pad = "0", side = "left"),
                      # Convert green
                      str_pad(string = as.hexmode(round(green, 0)),
                              width = 2, pad = "0", side = "left"),
                      # Convert blue
                      str_pad(string = as.hexmode(round(blue, 0)),
                              width = 2, pad = "0", side = "left")
  )
  
  return(hex_value)
}
myp2$hex <- rgb_to_hex(rgb)

# BUT mean RGB color is ugly and less accurate (ie, brown) for red/green polymorphic pops (ie SBL, SRS). Try median/ mode?
#summarize RGB color values by population
mode <- function(codes){
  which.max(tabulate(codes))
}

PopColor <- group_by(cap, pop) %>%
  summarise(meanRed = mean(red),
            meanGreen = mean(green),
            meanBlue = mean(blue),
            
            medRed = median(red),
            medGreen = median(green),
            medBlue = median(blue),
            
            modeRed = mode(red),
            modeGreen = mode(green),
            modeBlue = mode(blue)
  )

#convert to hexadecimal

#mean value
MeanRgb <- PopColor[,2:4]
#function rgb_to_hex: https://rdrr.io/github/stuart-morrison/schemr/man/rgb_to_hex.html
rgb_to_hex <- function(MeanRgb) {
  
  # Extract colour channels from first three colummns
  red <- unlist(MeanRgb[ , 1])
  green <- unlist(MeanRgb[ , 2])
  blue <- unlist(MeanRgb[ , 3])
  
  # Test if any RGB values are outside of [0, 255]
  if (any(red < 0 | red > 255) |
      any(green < 0 | green > 255) |
      any(blue < 0 | blue > 255)) {
    stop("Colour channels should be between 0 and 255.")
  }
  
  # Convert colour channels to hex
  hex_value <- paste0("#",
                      # Convert red
                      str_pad(string = as.hexmode(round(red, 0)),
                              width = 2, pad = "0", side = "left"),
                      # Convert green
                      str_pad(string = as.hexmode(round(green, 0)),
                              width = 2, pad = "0", side = "left"),
                      # Convert blue
                      str_pad(string = as.hexmode(round(blue, 0)),
                              width = 2, pad = "0", side = "left")
  )
  
  return(hex_value)
}
PopColor$MeanRgbHex <- rgb_to_hex(MeanRgb)

#median value
MedRgb <- PopColor[,5:7]
#function rgb_to_hex: https://rdrr.io/github/stuart-morrison/schemr/man/rgb_to_hex.html
rgb_to_hex <- function(MedRgb) {
  
  # Extract colour channels from first three colummns
  red <- unlist(MedRgb[ , 1])
  green <- unlist(MedRgb[ , 2])
  blue <- unlist(MedRgb[ , 3])
  
  # Test if any RGB values are outside of [0, 255]
  if (any(red < 0 | red > 255) |
      any(green < 0 | green > 255) |
      any(blue < 0 | blue > 255)) {
    stop("Colour channels should be between 0 and 255.")
  }
  
  # Convert colour channels to hex
  hex_value <- paste0("#",
                      # Convert red
                      str_pad(string = as.hexmode(round(red, 0)),
                              width = 2, pad = "0", side = "left"),
                      # Convert green
                      str_pad(string = as.hexmode(round(green, 0)),
                              width = 2, pad = "0", side = "left"),
                      # Convert blue
                      str_pad(string = as.hexmode(round(blue, 0)),
                              width = 2, pad = "0", side = "left")
  )
  
  return(hex_value)
}
PopColor$MedRgbHex <- rgb_to_hex(MedRgb)
medianHexPop <- PopColor$MedRgbHex

#mode value
ModeRgb <- PopColor[,8:10]
#function rgb_to_hex: https://rdrr.io/github/stuart-morrison/schemr/man/rgb_to_hex.html
rgb_to_hex <- function(ModeRgb) {
  
  # Extract colour channels from first three colummns
  red <- unlist(ModeRgb[ , 1])
  green <- unlist(ModeRgb[ , 2])
  blue <- unlist(ModeRgb[ , 3])
  
  # Test if any RGB values are outside of [0, 255]
  if (any(red < 0 | red > 255) |
      any(green < 0 | green > 255) |
      any(blue < 0 | blue > 255)) {
    stop("Colour channels should be between 0 and 255.")
  }
  
  # Convert colour channels to hex
  hex_value <- paste0("#",
                      # Convert red
                      str_pad(string = as.hexmode(round(red, 0)),
                              width = 2, pad = "0", side = "left"),
                      # Convert green
                      str_pad(string = as.hexmode(round(green, 0)),
                              width = 2, pad = "0", side = "left"),
                      # Convert blue
                      str_pad(string = as.hexmode(round(blue, 0)),
                              width = 2, pad = "0", side = "left")
  )
  
  return(hex_value)
}
PopColor$ModeRgbHex <- rgb_to_hex(ModeRgb)

#export PopColor with population summary values for color/ RGB and hex
#write.csv(PopColor, "CastPopulationColor_summary.csv")
#####

#Map and pie charts
#####

# #edited/compiled in excel, for KS map editing (added Lat/Long, taxa/species codes)
# #plot populations on map
pmap <- read.csv("Cast_mypMedHex_forMap.csv")

# library(sf)
# library("rnaturalearth")
# library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
#class(world)

myHex <- pmap$MedianRgbHex
names(myHex) <- pmap$pop
myHex


ggplot(data = world) +
  geom_sf() +
  geom_point(data = pmap, aes(x = lonW, y = latN, fill= pop, shape= sp), size= 4) +
  geom_text(data = pmap, aes(x = lonW, y = latN, label= pop), hjust= -.2, vjust= 0.1, size= 3)+
  scale_fill_manual(values = myHex)+
  scale_shape_manual(values=c('S'=21,'P'=22, 'L'=24, 'C'=23))+
  coord_sf(xlim = c(-115, -83), ylim = c(25, 50), expand = FALSE)

#exported this map with correct locations for study populations and median RGB floral colors
#as base for building map figure (Fig. 2) in Adobe Illustrator (KAS) and Inkscape (KEW)


#pool counts across years for pie charts 
#better for overall proportion (Crawley 2015, p.257)
#####
popCt <- group_by(myp2, pop) %>%
  summarise(sum_ct_hawkmoth = sum(ct_hawkmoth), 
            sum_ct_hummingbird = sum(ct_hummingbird),
            sum_ct_small_med_bee = sum(ct_small_med_bee),
            sum_ct_bumble_lg_bee = sum(ct_bumble_lg_bee),
            sum_ct_butterfly_all = sum(ct_butterfly_all),
            sum_ct_other = sum(ct_other),
            sum_visTotal = sum(visTotal))


#make individual dataframes for each population
CHCL <- popCt[1, c(1, 2:7)]
#make long
CHCLl <- melt(CHCL, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

#repeat for all pops
#make individual dataframes for each population
CLA <- popCt[2, c(1, 2:7)]
CLAl <- melt(CLA, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

CMSC <- popCt[3, c(1, 2:7)]
CMSCl <- melt(CMSC, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

CQL <- popCt[4, c(1, 2:7)]
CQLl <- melt(CQL, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

LMN <- popCt[5, c(1, 2:7)]
LMNl <- melt(LMN, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

LRR <- popCt[6, c(1, 2:7)]
LRRl <- melt(LRR, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

LVH <- popCt[7, c(1, 2:7)]
LVHl <- melt(LVH, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

PCM <- popCt[8, c(1, 2:7)]
PCMl <- melt(PCM, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

PMT <- popCt[9, c(1, 2:7)]
PMTl <- melt(PMT, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

PTH <- popCt[10, c(1, 2:7)]
PTHl <- melt(PTH, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

PTMS <- popCt[11, c(1, 2:7)]
PTMSl <- melt(PTMS, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

SBL <- popCt[12, c(1, 2:7)]
SBLl <- melt(SBL, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

SCL <- popCt[13, c(1, 2:7)]
SCLl <- melt(SCL, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

SDC <- popCt[14, c(1, 2:7)]
SDCl <- melt(SDC, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

SFP <- popCt[15, c(1, 2:7)]
SFPl <- melt(SFP, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

SIC <- popCt[16, c(1, 2:7)]
SICl <- melt(SIC, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

SILB <- popCt[17, c(1, 2:7)]
SILBl <- melt(SILB, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

SMP <- popCt[18, c(1, 2:7)]
SMPl <- melt(SMP, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

SNG <- popCt[19, c(1, 2:7)]
SNGl <- melt(SNG, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

SPB <- popCt[20, c(1, 2:7)]
SPBl <- melt(SPB, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

SRS <- popCt[21, c(1, 2:7)]
SRSl <- melt(SRS, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

SSC <- popCt[22, c(1, 2:7)]
SSCl <- melt(SSC, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

SYH <- popCt[23, c(1, 2:7)]
SYHl <- melt(SYH, id.vars= c("pop"), variable.name = "pollGroup", value.name= "sum_ct_poll") 

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )


CHCL <- ggplot(CHCLl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  ggtitle("CHCL")+
  theme(axis.text.x=element_blank())
CHCL

CLA <- ggplot(CLAl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  ggtitle("CLA")+
  theme(axis.text.x=element_blank())
CLA

CMSC <- ggplot(CMSCl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  ggtitle("CMSC")+
  theme(axis.text.x=element_blank())
CMSC

CQL <- ggplot(CQLl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  ggtitle("CQL")+
  theme(axis.text.x=element_blank())
CQL

LMN <- ggplot(LMNl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  ggtitle("LMN")+
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  theme(axis.text.x=element_blank())
LMN

LRR <- ggplot(LRRl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  ggtitle("LRR")+
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  theme(axis.text.x=element_blank())
LRR

LVH <- ggplot(LVHl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  ggtitle("LVH")+
  theme(axis.text.x=element_blank())
LVH

PCM <- ggplot(PCMl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  ggtitle("PCM")+
  theme(axis.text.x=element_blank())
PCM

PMT <- ggplot(PMTl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  ggtitle("PMT")+
  theme(axis.text.x=element_blank())
PMT

PTH <- ggplot(PTHl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  ggtitle("PTH")+
  theme(axis.text.x=element_blank())
PTH

PTMS <- ggplot(PTMSl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  ggtitle("PTMS")+
  theme(axis.text.x=element_blank())
PTMS

SBL <- ggplot(SBLl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  ggtitle("SBL")+
  theme(axis.text.x=element_blank())
SBL

SCL <- ggplot(SCLl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  ggtitle("SCL")+
  theme(axis.text.x=element_blank())
SCL

SDC <- ggplot(SDCl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  ggtitle("SDC")+
  theme(axis.text.x=element_blank())
SDC

SFP <- ggplot(SFPl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  ggtitle("SFP")+
  theme(axis.text.x=element_blank())
SFP

SIC <- ggplot(SICl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  ggtitle("SIC")+
  theme(axis.text.x=element_blank())
SIC

SILB <- ggplot(SILBl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  ggtitle("SILB")+
  theme(axis.text.x=element_blank())
SILB

SMP <- ggplot(SMPl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  ggtitle("SMP")+
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  theme(axis.text.x=element_blank())
SMP

SNG <- ggplot(SNGl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  ggtitle("SNG")+
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  theme(axis.text.x=element_blank())
SNG

SPB <- ggplot(SPBl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  ggtitle("SPB")+
  theme(axis.text.x=element_blank())
SPB

SRS <- ggplot(SRSl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  ggtitle("SRS")+
  theme(axis.text.x=element_blank())
SRS

SSC <- ggplot(SSCl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  ggtitle("SSC")+
  theme(axis.text.x=element_blank())
SSC

SYH <- ggplot(SYHl, aes(x="", y=sum_ct_poll, fill=pollGroup))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  scale_fill_manual(values=c('sum_ct_hawkmoth'="#0074A2",'sum_ct_hummingbird'="#ED145B", 
                             'sum_ct_small_med_bee'="#EADF6E", 'sum_ct_bumble_lg_bee'="#FFA500", 
                             "sum_ct_butterfly_all" = "#80B7C3", "sum_ct_other" = "#6F944E"))+
  ggtitle("SYH")+
  theme(axis.text.x=element_blank())
SYH

ggarrange(CHCL, 	CLA, 	CMSC, 	CQL, 	LMN, 	LRR, 	LVH, 	PCM, 	PMT, 	PTH, 	PTMS, 	
          SBL, 	SCL, 	SDC, 	SFP, 	SIC, 	SILB, 	SMP, 	SNG, 	SPB, 	SRS, 	SSC, 	SYH, 
          nrow=6, ncol=4, common.legend=T)
#Exported these pie charts and add them to map in Illustrator (KAS) for FIG2


#####

#MULTIVARIATE VARIATION ACROSS ALL POLLINATOR GROUPS AMONG SPECIES AND POPULATIONS
#use adonis2 for PERMANOVA test, permutational nonparametric manova
pr <- vegdist(myp2[,11:16], method= "gower")
p.div <- adonis2(pr ~ sp/pop + year+ dataset, data= myp2, permutations= 1000, method= "gower")
p.div

#FIGURE3- BOX PLOTS OF POLLINATOR VISIT COUNTS AND PROPORTIONS-
#####

#Proportion
mypPpr <- select(myp2, popYr, year, pop, dataset, taxa, sp, latN, lonW, 
                 pr_bumble_lg_bee, pr_butterfly_all, pr_hummingbird, pr_other,
                 pr_hawkmoth, pr_small_med_bee)

#make long for boxplots
mplp <- melt(mypPpr, id.vars= c("popYr", "pop", "year", "dataset", "taxa", "sp", "latN", "lonW"), 
             variable.name = "pollGroup", value.name= "proportion_poll") 

#average by pop for visualization
mplpAv <- group_by(mplp, pop, taxa, sp, latN, lonW, pollGroup) %>%
  summarise(proportion_poll = mean(proportion_poll))

PxsAvP <-  ggplot(data= mplpAv, aes(x= sp, y= proportion_poll)) +
  geom_boxplot(aes(fill=sp))+
  scale_fill_manual(values=c('S'='#b9d06a','P'='#b5266e', 'L'='#f13918', 'C'='#f8e432'))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  facet_wrap(vars(factor(pollGroup, level= c("pr_hummingbird", "pr_hawkmoth", "pr_bumble_lg_bee", 
                                             "pr_small_med_bee", "pr_butterfly_all", "pr_other"),
                         labels= c("pr_hummingbird"= "Hummingbirds", "pr_hawkmoth"= "Hawkmoths",
                                   "pr_bumble_lg_bee"="Lg Bees", "pr_small_med_bee"="Sm/Med Bees",
                                   "pr_butterfly_all"= "Butterflies", "pr_other"= "Other"))), 
             nrow=6)+
  labs(title= "A", y= "Average proportion of visits")+
  theme_classic()

PxsAvP
#FIGURE 3A

#Count
#select count poll data
mypPct <- select(myp2, popYr, year, pop, dataset, taxa, sp, latN, lonW, 
                 ct_bumble_lg_bee, ct_butterfly_all, ct_hummingbird, ct_other,
                 ct_hawkmoth, ct_small_med_bee)

#make long for boxplots
mplc <- melt(mypPct, id.vars= c("popYr", "pop", "year", "dataset", "taxa", "sp", "latN", "lonW"), 
             variable.name = "pollGroup", value.name= "count_poll") 

#average by pop for visualization
mplcAv <- group_by(mplc, pop, taxa, sp, latN, lonW, pollGroup) %>%
  summarise(avg_count_poll = mean(count_poll))

PxsAvC <-   ggplot(data= mplcAv, aes(x= sp, y= avg_count_poll)) +
  geom_boxplot(aes(fill=sp))+
  scale_fill_manual(values=c('S'='#b9d06a','P'='#b5266e', 'L'='#f13918', 'C'='#f8e432'))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  facet_wrap(vars(factor(pollGroup, level= c("ct_hummingbird", "ct_hawkmoth", "ct_bumble_lg_bee", 
                                             "ct_small_med_bee", "ct_butterfly_all", "ct_other"),
                         labels= c("ct_hummingbird"= "Hummingbirds", "ct_hawkmoth"= "Hawkmoths",
                                   "ct_bumble_lg_bee"="Lg Bees", "ct_small_med_bee"="Sm/Med Bees",
                                   "ct_butterfly_all"= "Butterflies", "ct_other"= "Other"))), 
             nrow=6)+
  labs(title= "Count data", y= "Average number of visits")+
  theme_classic()

PxsAvC
#FIG S3 


#STATS
#GLMM- Do species vary in proportion of visits from each floral visitor/functional group?

#PROPORTION DATA
#for proportion, compare binomial and betabinomial-- using DHARMa to assess overdispersion

#hummingbirds
m <- glmmTMB(pr_hummingbird ~ sp + year + dataset, weights= visTotal, family = "betabinomial", data= myp2)
summary(m)
Anova(m)
#pairwise comparisons
lsmeans(m, pairwise ~ sp, adjust= "tukey")

#assess overdispersion in model
sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#ns

#binomial
m <- glmmTMB(pr_hummingbird ~ sp + year + dataset, weights= visTotal, family = "binomial", data= myp2)
summary(m)
Anova(m)
#pairwise comparisons
lsmeans(m, pairwise ~ sp, adjust= "tukey")

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#sig deviations, worse than betabinomial

#hawkmoths
m <- glmmTMB(pr_hawkmoth ~ sp + year + dataset, weights= visTotal, family = "betabinomial", data= myp2)
summary(m)
Anova(m)

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#NS

#binomial
m <- glmmTMB(pr_hawkmoth ~ sp + year + dataset, weights= visTotal, family = "binomial", data= myp2)
summary(m)
Anova(m)

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#sig deviations, worse than betabinomial

#small bee
m <- glmmTMB(pr_small_med_bee ~ sp + year + dataset, weights= visTotal, family = "betabinomial", data= myp2)
summary(m)
Anova(m)

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#ns

#binomial
m <- glmmTMB(pr_small_med_bee ~ sp + year + dataset, weights= visTotal, family = "binomial", data= myp2)
summary(m)
Anova(m)

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#sig deviations, worse than betabinomial

#lg bee
m <- glmmTMB(pr_bumble_lg_bee ~ sp + year + dataset, weights= visTotal, family = "betabinomial", data= myp2)
summary(m)
Anova(m)

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#ns

#binomial
m <- glmmTMB(pr_bumble_lg_bee ~ sp + year + dataset, weights= visTotal, family = "binomial", data= myp2)
summary(m)
Anova(m)

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#sig deviations, worse than betabinomial

#butterfly
m <- glmmTMB(pr_butterfly_all ~ sp + year + dataset, weights= visTotal, family = "betabinomial", data= myp2)
summary(m)
Anova(m)
lsmeans(m, pairwise ~ sp, adjust= "tukey")

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#NS

#binomial
m <- glmmTMB(pr_butterfly_all ~ sp + year + dataset, weights= visTotal, family = "binomial", data= myp2)
summary(m)
Anova(m)
lsmeans(m, pairwise ~ sp, adjust= "tukey")

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#sig deviations, worse than betabinomial

#other
m <- glmmTMB(pr_other ~ sp + year + dataset, weights= visTotal, family = "betabinomial", data= myp2)
summary(m)
Anova(m)

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB) 
#ns

#binomial
m <- glmmTMB(pr_other ~ sp + year + dataset, weights= visTotal, family = "binomial", data= myp2)
summary(m)
Anova(m)

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB) 
#sig deviations, worse than betabinomial

#COUNT DATA
#for counts, compared models with poisson distribution to negative binomial distribution (nbinom1) using DHARMa package

#hummingbirds
#model won't converge if sessiliflora included (due to all zeroes), so run with only purpurea species (taxa= "PC")
m <- glmmTMB(ct_hummingbird ~ sp + year + dataset, family = "nbinom1", data= subset(myp2, taxa== "PC"))
summary(m)
Anova(m)
#pairwise comparisons
lsmeans(m, pairwise ~ sp, adjust= "tukey")

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#ns

#poisson
m <- glmmTMB(ct_hummingbird ~ sp + year + dataset, family = "poisson", data= subset(myp2, taxa== "PC"))
summary(m)
Anova(m)

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#sig deviations, worse than nbinom1

#hawkmoths
m <- glmmTMB(ct_hawkmoth ~ sp + year + dataset, family = "nbinom1", data= myp2)
summary(m)
Anova(m)
sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#NS

#poisson
m <- glmmTMB(ct_hawkmoth ~ sp + year + dataset, family = "poisson", data= myp2)
summary(m)
Anova(m)
sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#sig deviations, worse than nbinom1

#small bee
m <- glmmTMB(ct_small_med_bee ~ sp + year + dataset, family = "nbinom1", data= myp2)
summary(m)
Anova(m)
lsmeans(m, pairwise ~ sp, adjust= "tukey")

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#NS

#poisson
m <- glmmTMB(ct_small_med_bee ~ sp + year + dataset, family = "poisson", data= myp2)
summary(m)
Anova(m)
sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
##sig deviations, worse than nbinom1

#lg bee
m <- glmmTMB(ct_bumble_lg_bee ~ sp + year + dataset, family = "nbinom1", data= myp2)
summary(m)
Anova(m)

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#ns

#poisson
m <- glmmTMB(ct_bumble_lg_bee ~ sp + year + dataset, family = "poisson", data= myp2)
summary(m)
Anova(m)

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#sig deviations, worse than nbinom1

#butterfly
m <- glmmTMB(ct_butterfly_all ~ sp + year + dataset, family = "nbinom1", data= myp2)
summary(m)
Anova(m)
lsmeans(m, pairwise ~ sp, adjust= "tukey")

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#ns

#poisson
m <- glmmTMB(ct_butterfly_all ~ sp + year + dataset, family = "poisson", data= myp2)
summary(m)
Anova(m)
lsmeans(m, pairwise ~ sp, adjust= "tukey")

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#sig deviations, worse than nbinom1

#other
m <- glmmTMB(ct_other ~ sp + year + dataset, family = "nbinom1", data= myp2)
summary(m)
Anova(m)

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#NS

#poisson
m <- glmmTMB(ct_other ~ sp + year + dataset, family = "poisson", data= myp2)
summary(m)
Anova(m)

sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#sig deviations, worse than nbinom1

#####


#FIGURE4- CASE POLL PROP BY LAT
#####
#set up/fix population median floral colors
myHex <- pmap$MedianRgbHex
names(myHex) <- pmap$pop
myHex

#need to change abbreviations of a few pops (SDC>SDCR, SNG>SNRM)
myHex1 <- myHex
names(myHex1)[14] <- "SDCR"
names(myHex1)[19] <- "SNRM"
myHex1
#####

#hawkmoth

hmls <- ggplot(data=subset(myp2, sp== "S"), aes(x= latN, y=pr_hawkmoth, fill= pop))+
  geom_point(size=3, pch= 21)+
  geom_smooth(inherit.aes = F, aes(x= latN, y=pr_hawkmoth, succ= ct_hawkmoth, fail= visTotal-ct_hawkmoth),
              linetype= "solid", color= "black",
              method="glm", method.args=list(family="betabinomial"),
              formula = cbind(succ, fail) ~ x,
              alpha=0.2, fullrange=FALSE)+
  scale_fill_manual(values=c(myHex1))+
  #labs(x= "Population Latitude (N)", y= "Hawkmoths")+
  labs(title= "A", x= NULL, y= "Hawkmoths")+
  theme_classic()

hmls

#hm
m <- glmmTMB(pr_hawkmoth ~ latN
             + year + dataset, weights=visTotal, family= betabinomial, data= subset(myp2, sp== "S"))
summary(m)
Anova(m)

# sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
# plot(sim_resid_glmmTMB)

#sm bee
sbls <- ggplot(data=subset(myp2, sp== "S"), aes(x= latN, y= pr_small_med_bee, fill= pop))+
  geom_point(size=3, pch= 21)+
  geom_smooth(inherit.aes = F, aes(x= latN, y=pr_small_med_bee, succ= ct_small_med_bee, fail= visTotal-ct_small_med_bee), 
              linetype= "solid", color= "black",
              method="glm", method.args=list(family="betabinomial"),
              formula = cbind(succ, fail) ~ x,
              alpha=0.2, fullrange=FALSE)+
  scale_fill_manual(values=c(myHex1))+
  #labs(x= "Population Latitude (N)", y= "Small/med bees")+
  labs(title= "B",x= NULL, y= "Small/med bees")+
  theme_classic()

sbls

#sb
m <- glmmTMB(pr_small_med_bee ~ latN
             + year + dataset, weights=visTotal, family= betabinomial, data= subset(myp2, sp== "S"))
summary(m)
Anova(m)

# sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
# plot(sim_resid_glmmTMB)


#bb
bbls <- ggplot(data=subset(myp2, sp== "S"), aes(x= latN, y= pr_bumble_lg_bee, fill= pop))+
  geom_point(size=3, pch= 21)+
  # geom_smooth(inherit.aes = F, aes(x= latN, y=pr_bumble_lg_bee, succ= ct_bumble_lg_bee, fail= visTotal-ct_bumble_lg_bee), 
  #             linetype= "dashed", color= "black",
  #             method="glm", method.args=list(family="betabinomial"),
  #             formula = cbind(succ, fail) ~ x,
  #             alpha=0.2, fullrange=FALSE)+
  scale_fill_manual(values=c(myHex1))+
  #labs(x= "Population Latitude (N)", y= "Large bees")+
  labs(title= "C", x= NULL, y= "Large bees")+
  theme_classic()

bbls

#bb
m <- glmmTMB(pr_bumble_lg_bee ~ latN
             + year + dataset, weights=visTotal, family= betabinomial, data= subset(myp2, sp== "S"))
summary(m)
Anova(m)



#bf
bfls <- ggplot(data=subset(myp2, sp== "S"), aes(x= latN, y= pr_butterfly_all, fill= pop))+
  geom_point(size=3, pch= 21)+
  # geom_smooth(inherit.aes = F, aes(x= latN, y=pr_butterfly_all, succ= ct_butterfly_all, fail= visTotal-ct_butterfly_all), 
  #             linetype= "dashed", color= "black",
  #             method="glm", method.args=list(family="betabinomial"),
  #             formula = cbind(succ, fail) ~ x,
  #             alpha=0.2, fullrange=FALSE)+
  scale_fill_manual(values=c(myHex1))+
  ylim(0,1)+
  #labs(x= "Geographic position across range (lat)", y= "Butterflies")+
  labs(title= "D",x= NULL, y= "Butterflies")+
  theme_classic()
#no trend line here: Warning message:
#glm.fit: fitted probabilities numerically 0 or 1 occurred 

bfls

m <- glmmTMB(pr_butterfly_all ~ latN
             + year + dataset, weights=visTotal, family= betabinomial, data= subset(myp2, sp== "S"))
summary(m)
#NA- can't run (because only one population/ one Lat value recorded butterfly visits)

#hummingbirds- no visits to sessiliflora

#other

ols <- ggplot(data=subset(myp2, sp== "S"), aes(x= latN, y= pr_other, fill= pop))+
  geom_point(size=3, pch= 21)+
  # geom_smooth(inherit.aes = F, aes(x= latN, y=pr_other, succ= ct_other, fail= visTotal-ct_other),
  #             linetype= "dashed", color= "black",
  #             method="glm", method.args=list(family="betabinomial"),
  #             formula = cbind(succ, fail) ~ x,
  #             alpha=0.2, fullrange=FALSE)+
  scale_fill_manual(values=c(myHex1))+
  #geom_smooth(method="loess")+
  ylim(0,1)+
  labs(title= "E", x= NULL, y= "Other visitors")+
  #labs(x= "Geographic position across range (lat)", y= "Other visitors")+
  theme_classic()

ols

m <- glmmTMB(pr_other ~ latN
             + year + dataset, weights=visTotal, family= betabinomial, data= subset(myp2, sp== "S"))
summary(m)
Anova(m)

# sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
# plot(sim_resid_glmmTMB)




#####

#FRUIT SET BY SPECIES AND BY POP LATITUDE
#####
#Fruit set by species
fs1 <- read.csv("Cast_fruitSet_raw.csv")
fs1$year <- as.factor(fs1$year)
str(fs1)
#calculate fruit-flower ratio
fs1 <- fs1 %>% 
  mutate(FrtSet = nFrts / nFlrs)
#DATAFRAME FOR FRUIT SET DATA
#check pop-years are unique
unique(fs1$popYear)

#FRT SET BY SPECIES
fsS <- ggplot(data = fs1, aes(x= sp, y= FrtSet))+
  geom_violin(aes(fill= sp), trim= T, alpha= 0.9)+
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="black")+
  scale_fill_manual(values=c('S'='#b9d06a','P'='#b5266e', 'L'='#f13918', 'C'='#f8e432'))+
  labs(title= "B", x= "Species", y= "Fruit set")+
  theme_classic()

fsS 
#FIG 3B

#STATS
#species- when pop included as random effect (nested or unnested), sig effect disappears
m <- glmmTMB(FrtSet ~ sp + year + (1|pop), weights=nFlrs, family = "betabinomial", data= fs1)
summary(m)
Anova(m)
# sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
# plot(sim_resid_glmmTMB)
# #not ideal, but better than binomial (below)
# 
# m <- glmmTMB(FrtSet ~ sp + year + (1|pop), weights=nFlrs, family = "binomial", data= fs1)
# summary(m)
# Anova(m)
# sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
# plot(sim_resid_glmmTMB)


#FRUIT SET BY GEOGRAPHY/LATITUDE (SESSILIFLORA ONLY)
####
fsl <- ggplot(data = subset(fs1, sp== "S"), aes(x= latN, y= FrtSet, shape= sp, fill= pop))+
  #geom_smooth(method="glm")+ #right now, plotting by pop
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", size= .6)+
  geom_smooth(inherit.aes = F, aes(x= latN, y=FrtSet, succ= nFrts, fail= nFlrs-nFrts),
              linetype= "solid", color= "black",
              method="glm", method.args=list(family="betabinomial"),
              formula = cbind(succ, fail) ~ x,
              alpha=0.2, fullrange=FALSE)+
  scale_fill_manual(values=c(myHex1))+
  scale_shape_manual(values=c('S'=21))+
  ylim(0,1)+
  labs(title= "F", x= NULL, y= "Fruit set")+
  theme_classic()

fsl

m <- glmmTMB(FrtSet ~ latN + year, weights=nFlrs, family = "betabinomial", data= subset(fs1, sp=="S"))
summary(m)
Anova(m)
# sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
# plot(sim_resid_glmmTMB)
# 
# m <- glmmTMB(FrtSet ~ latN + year, weights=nFlrs, family = "binomial", data= subset(fs1, sp=="S"))
# summary(m)
# Anova(m)
# sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
# plot(sim_resid_glmmTMB)
# #betabinomial is better

#FIGURE 4
ggarrange(hmls, sbls, bbls, bfls, ols, fsl, nrow=2, ncol=3, legend= F)
#Add shared x axis label: "Geographic position across range (N lat)"
#add in inkscape


#FLORAL TRAITS X POLLINATOR VISITATION

#scale floral traits
#in separate dataframe
myp3_st <- 
  myp2 %>% 
  mutate(corL_s = scale(corL),
         corW_s = scale(corW),
         stigma_s = scale(stigma),
         lipL_s = scale(lipL),
         brLbW_s = scale(brLbW),
         avgNMDS1_s = scale(avgNMDS1),
         avgNMDS2_s = scale(avgNMDS2)
  )
str(myp3_st)


#TABLE2- MULTIPLE REGRESSION
#####
#hummingbird
#betabinomial
m <- glmmTMB(pr_hummingbird ~ corL_s + corW_s + lipL_s + stigma_s + brLbW_s + avgNMDS1_s + avgNMDS2_s + 
               year + dataset, weights=visTotal, family = betabinomial, data= myp3_st)
summary(m)
Anova(m)

#check overdispersion
sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)

# #binomial model
# m <- glmmTMB(pr_hummingbird ~ corL_s + corW_s + lipL_s + stigma_s + brLbW_s + avgNMDS1_s + avgNMDS2_s + 
#                year + dataset, weights=visTotal, family = binomial, data= myp3_st)
# summary(m)
# Anova(m)
# 
# sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
# plot(sim_resid_glmmTMB)
# #looks worse than betabinomial

#hawkmoth
m <- glmmTMB(pr_hawkmoth~ corL_s + corW_s + lipL_s + stigma_s + brLbW_s + avgNMDS1_s + avgNMDS2_s +
               year+ dataset, weights=visTotal, family = betabinomial, data= myp3_st)
summary(m)
Anova(m)

# #check overdispersion
# sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
# plot(sim_resid_glmmTMB)
# #NS
# 
# #binomial error dist
# m <- glmmTMB(pr_hawkmoth~ corL_s + corW_s + lipL_s + stigma_s + brLbW_s + avgNMDS1_s + avgNMDS2_s +
#                year+ dataset, weights=visTotal, family = binomial, data= myp3_st)
# summary(m)
# Anova(m)
# 
# #check overdispersion
# sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
# plot(sim_resid_glmmTMB)
# #sig deviations, worse than betabinomial


#bumblebees
m <- glmmTMB(pr_bumble_lg_bee~ corL_s + corW_s + lipL_s + stigma_s + brLbW_s + avgNMDS1_s + avgNMDS2_s +
               year+ dataset, weights=visTotal, family = betabinomial, data= myp3_st)
summary(m)
Anova(m)

# sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
# plot(sim_resid_glmmTMB)
# #ns
# 
# #binomial
# m <- glmmTMB(pr_bumble_lg_bee~ corL_s + corW_s + lipL_s + stigma_s + brLbW_s + avgNMDS1_s + avgNMDS2_s +
#                year+ dataset, weights=visTotal, family = binomial, data= myp3_st)
# summary(m)
# Anova(m)
# 
# sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
# plot(sim_resid_glmmTMB)
# #sig deviations, worse than betabinomial

#small bees
m <- glmmTMB(pr_small_med_bee~ corL_s + corW_s + lipL_s + stigma_s + brLbW_s + avgNMDS1_s + avgNMDS2_s +
               year + dataset, weights=visTotal, family = betabinomial, data= myp3_st)
summary(m)
Anova(m)

# sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
# plot(sim_resid_glmmTMB)
# 
# #binomial
# m <- glmmTMB(pr_small_med_bee~ corL_s + corW_s + lipL_s + stigma_s + brLbW_s + avgNMDS1_s + avgNMDS2_s +
#                year + dataset, weights=visTotal, family = binomial, data= myp3_st)
# summary(m)
# Anova(m)
# 
# sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
# plot(sim_resid_glmmTMB)
# #looks worse (more deviations) than betabinomial


#butterfly
m <- glmmTMB(pr_butterfly_all ~ corL_s + corW_s + lipL_s + stigma_s + brLbW_s + avgNMDS1_s + avgNMDS2_s +
               year+ dataset, weights=visTotal, family = betabinomial, data= myp3_st)
summary(m)
Anova(m)

# sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
# plot(sim_resid_glmmTMB)
# 
# #binomial
# m <- glmmTMB(pr_butterfly_all ~ corL_s + corW_s + lipL_s + stigma_s + brLbW_s + avgNMDS1_s + avgNMDS2_s +
#                year+ dataset, weights=visTotal, family = binomial, data= myp3_st)
# summary(m)
# Anova(m)
# 
# sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
# plot(sim_resid_glmmTMB)
# #more deviations than betabinomial

#other
m <- glmmTMB(pr_other ~ corL_s + corW_s + lipL_s + stigma_s + brLbW_s + avgNMDS1_s + avgNMDS2_s +
               year+ dataset, weights=visTotal, family = betabinomial, data= myp3_st)
summary(m)
Anova(m)

# sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
# plot(sim_resid_glmmTMB)
# 
# #binomial
# m <- glmmTMB(pr_other ~ corL_s + corW_s + lipL_s + stigma_s + brLbW_s + avgNMDS1_s + avgNMDS2_s +
#                year+ dataset, weights=visTotal, family = binomial, data= myp3_st)
# summary(m)
# Anova(m)
# 
# sim_resid_glmmTMB<-DHARMa::simulateResiduals(m, 1000)
# plot(sim_resid_glmmTMB)
# #more deviations than betabinomial


#plot some relationships between floral traits and pollinators
#NMDS2 and lgbee
n2lb <- ggplot(data= myp3_st, aes(x= avgNMDS2_s, y= pr_bumble_lg_bee, fill= pop, shape= sp))+
  geom_point(size= 3)+
  scale_shape_manual(values=c('S'=21,'P'=22, 'L'=24, 'C'=23))+
  scale_fill_manual(values=myHex1)+
  labs(title= "A", x= "Floral color (Avg NMDS2 value)", y= "Large Bee visitation")+
  geom_smooth(inherit.aes = F, aes(x= avgNMDS2_s, y=pr_bumble_lg_bee, succ= ct_bumble_lg_bee, fail= visTotal-ct_bumble_lg_bee), 
              linetype= "solid", color= "black",
              method="glm", method.args=list(family="betabinomial"),
              formula = cbind(succ, fail) ~ x,
              alpha=0.2, fullrange=FALSE)+#scale_linetype_manual(values= c("dashed", "solid"))+
  theme_classic()

n2lb

#NMDS2 and butterflies
n2bf <- ggplot(data= myp3_st, aes(x= avgNMDS2_s, y= pr_butterfly_all, fill= pop, shape= sp))+
  geom_point(size= 3)+
  scale_shape_manual(values=c('S'=21,'P'=22, 'L'=24, 'C'=23))+
  scale_fill_manual(values=myHex1)+
  labs(title= "B", x= "Floral color (Avg NMDS2 value)", y= "Butterfly visitation")+
  ylim(0,1)+
  geom_smooth(inherit.aes = F, aes(x= avgNMDS2_s, y=pr_butterfly_all, succ= ct_butterfly_all, fail= visTotal-ct_butterfly_all), 
              linetype= "solid", color= "black",
              method="glm", method.args=list(family="betabinomial"),
              formula = cbind(succ, fail) ~ x,
              alpha=0.2, fullrange=FALSE)+#scale_linetype_manual(values= c("dashed", "solid"))+
  theme_classic()

n2bf

#stigma x small bee
stsb <- ggplot(data= myp3_st, aes(x= stigma, y= pr_small_med_bee, fill= pop, shape= sp))+
  geom_point(size= 3)+
  scale_shape_manual(values=c('S'=21,'P'=22, 'L'=24, 'C'=23))+
  scale_fill_manual(values=myHex1)+
  labs(title= "C",x= "Stigma exsertion (mm)", y= "Small/med bee visitation")+
  ylim(0,1)+
  geom_smooth(inherit.aes = F, aes(x= stigma, y=pr_small_med_bee, succ= ct_small_med_bee, fail= visTotal-ct_small_med_bee), 
              linetype= "solid", color= "black",
              method="glm", method.args=list(family="betabinomial"),
              formula = cbind(succ, fail) ~ x,
              alpha=0.2, fullrange=FALSE)+#scale_linetype_manual(values= c("dashed", "solid"))+
  theme_classic()
stsb

#FIGURE 5A-C
ggarrange(n2lb, n2bf, stsb, nrow=2, ncol=2, legend = F)

#####

#FIGURE 5D

#nmds with poll vis
#add shapes to code species on ordination
myp2$pch_ind[myp2$sp == "P"] <- 22 #square
myp2$pch_ind[myp2$sp == "S"] <- 21 #circle
myp2$pch_ind[myp2$sp == "C"] <- 23 #diamond
myp2$pch_ind[myp2$sp == "L"] <- 24 #triangle

#pop mean floral traits- values in columns 27-34
p.mds <- metaMDS(myp2[ ,27:34], distance= "gower", zerodist="add")
p.mds
stressplot(p.mds)
#stress= 0.107
myp2$nmds1 <- p.mds$points[,1]
myp2$nmds2 <- p.mds$points[,2]

#plotted with hex inflor colors
fig <- ordiplot(p.mds, type = "none")
points(fig, "sites", pch= myp2$pch_ind, bg=myp2$hex,
       cex=1.3)
#add species labels to centroids and draw st dev
for(i in unique(myp2$sp)) {
  ordiellipse(p.mds$point[grep(i,myp2$sp),],draw="polygon", alpha = 90, kind= "sd",
              groups=myp2$sp[myp2$sp==i], label= T)
}

#add trait loadings:
#for prop vis from pollinator groups
fit <- envfit(p.mds ~ pr_hummingbird + pr_hawkmoth + pr_other + pr_small_med_bee 
               + pr_bumble_lg_bee + pr_butterfly_all 
               ,data= myp2, perm = 1000)
fit
plot(fit, p.max = 0.05, col = "black", cex=1)

#final dataframe
#write.csv(myp2, "Castilleja_pollMosaics_analysisDatav1.csv")
