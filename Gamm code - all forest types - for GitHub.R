### Analysis Code – GAMMs for all forest types combined###

# Install / load packages 
library(gam)
library(dplyr)
library(plyr)
library(ggplot2)
library(tidyr)
library(mgcv)
library(nlme)
library(gamm4)
library(lme4)

# Set the working directory to access data files 
setwd ("location") 

##Data file columns were: Dataset name, biotic group (with subgroup and sub-subgroup for arthropods only),
##forest type (conifer, mixed, broadleaf), years-post-harvest (YPH), 
##code for the reference replicate within dataset, code for the treatment replicate within dataset,
##ChaoJaccardSimilarity value


##Data files  used in the analyses are available at:
##https://doi.org/10.17632/yt42chsskj.1

###########################################################################################################

# Load in data for the given biotic group
Alltrt <- read.csv("Biota_AllForest_TRTvsREF.csv")

#change "Forest.Type" to a factor variable
fForest.Type <- as.factor(Alltrt$Forest.Type)

#run GAMM for all forest types combined
gammAll <- gamm4(ChaoJaccardSimilarity ~ s(YPH, fx = FALSE, bs = 'cr'), random = ~(1|Dataset.Name/TRT.being.compared.with.REF), data = Alltrt, REML = TRUE) 
                                                                                          
#model testing for differences among forest types
gammAllFT <- gamm4(ChaoJaccardSimilarity ~ fForest.Type + s(YPH, by=fForest.Type, fx = FALSE, bs = 'cr'), random = ~(1|Dataset.Name/TRT.being.compared.with.REF), data = Alltrt, REML = TRUE)

#compare AIC of the model without versus with "by=fForest.type"
AIC(gammAll$mer, gammAllFT$mer)

#results from the model with all forest types combined
summary(gammAll$gam)
summary(gammAll$mer)
gam.check(gammAll$gam)

#results from the model with "by fForest.type"
summary(gammAllFT$gam)
anova(gammAllFT$gam)
gam.check(gammAllFT$gam)

##Evidence for differences in response curve among forest types determined by:
##Comparing AIC values between the two models
##Comparing significance and estimated d.f. of smoothers among forest types
##Checking significance of "forest type" as a predictor




