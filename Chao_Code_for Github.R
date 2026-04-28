# Chao-Jaccard Index Code #
###########################
#Code for calculating ChaoJaccard similarity of community composition between 
#treated (clearcut) replicates and Reference (mature forest) replicates
#and between Reference and Reference replicates within each dataset.
#Calculated for community composition of different biotic groups:
#arthropods (with subgroups: beetles, spiders, other arthropods; and sub-subgroups:
#carabids, saproxylic beetles, other beetles, spiders, other arthropods), birds,
#small mammals, bryophytes, lichens, vascular plants

#################################################################################
### Load in required packages
install.packages("CommEcol")
install.packages("vegan")
install.packages("tidyverse")
library(CommEcol)
library(vegan)
library(tidyverse) # this is for mutating the data frame to make blank cells = 0

### Set working directory
setwd("drive_location")

### Load in data
##data arranged as one row for each replicate within a dataset
##columns were: code specific to dataset, replicate within dataset and treatment (treated or reference),
##list of species
##values within cells were abundance or relative abundance or cover for each species

DF <- read_csv("Data_File_Name.csv")

### Create a new data frame that contains only the code column and species columns
## row.names has to be TRUE or it will try to run dis.chao() on the row names too
## the second entry into the new_DF will need to be changed based on how many species there are
new_DF <- data.frame(DF[,c(1, 8:317)], row.names = TRUE)

### Convert NA values to 0s and make sure code is set to be row names
noNA_DF <- new_DF %>%
  mutate_all(~replace(., is.na(.), 0))

### Run Chao index for each data set
disChao <- as.matrix(dis.chao(noNA_DF, index = "jaccard", version = "rare"))

### Create another data frame of similarity indices for each data set
simChao <- 1-disChao

### Create and save CSV files containing similarity values for the dataset
write.csv(simChao, 
          "/dataset_name_similarity.csv")

### Format Chao outputs to only include the pairs we want
##i.e. we don't need treatment versus treatment similarity
##in some cases we needed to constrain reference versus reference similarity 
##to data collected in the same sampling year.
##Specifics varied with the number of treatment and reference replicates

# extract treatment vs reference 
trtVSref <- data.frame(simChao[13:24, 1:12])

# extract reference vs reference values
refVSref <- data.frame(simChao[c(13:24), c(13:24)])

### save trtVSref and refVSref data as csv files

write.csv(trtVSref, 
          "/dataset_name_trtvsref.csv")
write.csv(refVSref_na, 
          "/dataset_name_refVSref.csv")

##the "treatment versus reference" and "reference versus reference"
##similarity values for each dataset were then manually combined into 
##one "treatment versus reference" and one "reference versus reference" file 
##for each biotic group
##each row was for a replicate "treatment versus reference" or "reference versus reference"
##similarity value 
##columns were: dataset_code, biotic_group, subgroup and sub-subgroup (for arthropods only)
##forest_type (conifer, mixed, broadleaf), years_post_harvest (for treatment versus reference only),
##two columns with codes for reference or treatment replicates within dataset that were compared with one another
##Chao-Jaccard similarity value.

##Those files were used into the analyses for all forest types combined.
##They were also split into separate files by forest type for the analyses
##done separately by forest type.
##Arthropod data were also split into separate files by sub-group and sub-sub-group
##for separate analyses (All arthropods, All Beetles, Carabid beetles, Saproxylic beetles, Spiders).

##Data files created and used in subsequent analyses are available at:
##https://doi.org/10.17632/yt42chsskj.1





