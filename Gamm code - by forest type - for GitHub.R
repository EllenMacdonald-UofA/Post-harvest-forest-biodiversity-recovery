###Code for running GAMMs by forest type separately

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

#install package Rmisc to make summarySE work
install.packages("Rmisc", dependencies = TRUE)

# Set the working directory 
setwd ("location")

# Load in data
#this is the code for one biotic group (or subgroup, or sub-subgroup for Arthropods) for one forest type
Trt <- read.csv("Biota_ForestType_TRTvsREF.csv")
Ref <- read.csv("Biota_ForestType_REFvsREF.csv")

##Data file columns were: Dataset name, biotic group (with subgroup and sub-subgroup for arthropods only),
##forest type (conifer, mixed, broadleaf), years-post-harvest (YPH, for treatment vs reference only), 
##Two columns with codes for the reference or treatment replicates within dataset,
##ChaoJaccardSimilarity value


##Data files  used in the analyses are available at:
##https://doi.org/10.17632/yt42chsskj.1

###########################################################################################################

#GAMM (k undefined) for post-harvest recovery of community composition similarity
##between harvested (clearcut) and reference (mature forest) replicates
##nested to account for the data structure with several treated versus reference replicates within each dataset
gamm <- gamm4(ChaoJaccardSimilarity ~ s(YPH, fx = FALSE, bs = 'cr'), random = ~(1|Dataset.Code/TRT.being.compared.with.REF), data = Trt, REML = TRUE) 

#get the full results from the model 
#note significance of smoother, edf, F, k'
#note n (number of similarity value comparisons), # of datasets, # of replicate comparisons within datasets
summary(gamm$gam)
summary(gamm$mer)
gam.check(gamm$gam)


##########################################################################################################################

#create predicted values from the model, for plotting
pre.gamm <- predict(gamm$gam, se.fit = TRUE)

##########################################################################################################################

## calculating the mean similarity between reference replicates within datasets
##to avoid excessive influence from studies with more replicates we
##calculated mean reference versus reference similarity within each dataset and
##then averaged that across datasets
Code <- list(Ref$Dataset.Code)
AVERAGED <- aggregate(Ref$ChaoJaccardSimilarity, 
                      by = Code, FUN = mean)

#code calculating the 95% confidence interval for the ref vs ref similarity (for plot error bars)
#This creates the function "summarySE" from the Rmisc package
#you just have to run this once in an R session
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) 
  
{
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
#End of code

# Run summarySE" to Output 95% CI for ref vs ref similarity values 

summarySE(Ref, measurevar = "ChaoJaccardSimilarity")
#This outputs the mean, s.d., s.e., and ci of the reference vs reference similarity values
#enter the value of C.I. that is produced here
#and run this next line to define ci, this is an example value for CI
ci <- 0.024

###################################################################################################################

##plot the GAMM curve, with standard error around it, for clearcut versus reference
##similarity over time post-harvest.
##Add the mean similarity between reference replicates, with 95% confidence, to the right of the graph

##geom_ribbon creates a shaded area on the plot, alpha determines the transparency
PLOT <- ggplot(data = Trt, aes(x = YPH, y = ChaoJaccardSimilarity)) +
  geom_point() +
  geom_ribbon(aes(ymin = c(pre.gamm$fit - pre.gamm$se.fit),
                  ymax = pre.gamm$fit), alpha=0.3) +
  geom_ribbon(aes(ymin = c(pre.gamm$fit + pre.gamm$se.fit),
                  ymax = pre.gamm$fit), alpha=0.3) +
  geom_line(colour="black", aes(y = pre.gamm$fit)) +             
##plotting the ref vs ref mean similarity (of the means by dataset.code)
##change position on the x axis (x =  ) based on the max YPH to place this to the right of the end of the GAMM curve
##don't need "data = " because you define x and y in the aes
   geom_point(aes(x = 85, y = mean(AVERAGED$x)), colour="black") +
##plotting the 95% confidence interval error bars for the mean similarity among reference replicates 
##change position on the x axis (x = ) to match that for mean reference versus reference similarity
##change the x scale and breaks depending on the maximum Years post-harvest
##don't need "data = " because you define x and y in the aes
  geom_errorbar(aes(ymin = mean(AVERAGED$x) - ci, ymax = mean(AVERAGED$x) + ci, x = 85), colour="black") +
  labs(
    x="Years Post Harvest",
    y="Chao Jaccard Similarity",
    title="Biota - Forest type") +
  scale_x_continuous(breaks = c(10,20,30, 40, 60, 80)) +
  scale_y_continuous(breaks = c(0.2,0.4,0.6,0.8,1)) + 
  theme_minimal()
PLOT

#####################################################################################################################
## Save the plot
## Figure file name = biota_foresttype.jpeg
ggsave("Biota_ForestType.jpeg", height = 3, width = 6, units = "in")

####################################################################################################################
