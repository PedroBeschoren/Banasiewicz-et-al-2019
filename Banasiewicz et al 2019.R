# This code was written by Pedro Beschoren da Costa in mid 2019
# it employs some classical microbial ecology methods (NMDS, PErmanova) to a nifD metagenomic sequence dataset
# It also utulizes a mosaic plot from vcd package


library("vcd")
library("dplyr")
library("tidyr")
library("ggplot2")
library("tibble")
library(vegan)
library("ggfortify")

#sload NifD sequence data

nifD_Metagenomic_CaseForm <- read.csv("Input/nifD_Metagenomic sequences table.csv") 

####Calcualte frequencies of sequences #####
# calculate frequency of each clade in a new dataframe 
CladeFreq<-as.data.frame(summary(nifD_Metagenomic_CaseForm$Clade))%>%
  rownames_to_column("Clade")%>% # remove rownames
  rename(Clade_Frequency="summary(nifD_Metagenomic_CaseForm$Clade)") # rename variable
nifD_Metagenomic_CladeFreq<- full_join(CladeFreq, nifD_Metagenomic_CaseForm, by= "Clade" ) # Join with original dataset
nifD_Metagenomic_CladeFreq$Clade<-as.factor(nifD_Metagenomic_CladeFreq$Clade)


# calculate frequency of each legume, remove rownames, change variable name, join with original dataset
LegumeFreq<-as.data.frame(summary(nifD_Metagenomic_CaseForm$Legume.species))%>%
  rownames_to_column("Legume.species")%>% # remove rownames
  rename(Legume_Frequency="summary(nifD_Metagenomic_CaseForm$Legume.species)") # rename variable
nifD_Metagenomic_CladeFreq<- full_join(LegumeFreq, nifD_Metagenomic_CladeFreq, by= "Legume.species" ) # Join with original dataset
nifD_Metagenomic_CladeFreq$Legume.species<-as.factor(nifD_Metagenomic_CladeFreq$Legume.species)

# calculate frequency of soil location, remove rownames, change variable name, join with original dataset
SoilFreq<-as.data.frame(summary(nifD_Metagenomic_CaseForm$Soil.data))%>%
  rownames_to_column("Soil.data")%>% # remove rownames
  rename(Soil_Frequency="summary(nifD_Metagenomic_CaseForm$Soil.data)") # rename variable
nifD_Metagenomic_CladeFreq<- full_join(SoilFreq, nifD_Metagenomic_CladeFreq, by= "Soil.data" ) # Join with original dataset
nifD_Metagenomic_CladeFreq$Soil.data<-as.factor(nifD_Metagenomic_CladeFreq$Soil.data)

# calculate frequency of primers, remove rownames, change variable name, join with original dataset
nifD_Metagenomic_CaseForm$Internal.external<-as.factor(nifD_Metagenomic_CaseForm$Internal.external)
PrimerFreq<-as.data.frame(summary(nifD_Metagenomic_CaseForm$Internal.external))%>%
  rownames_to_column("Internal.external")%>% # remove rownames
  rename(Primer_Frequency="summary(nifD_Metagenomic_CaseForm$Internal.external)") # rename variable
nifD_Metagenomic_CladeFreq<- full_join(PrimerFreq, nifD_Metagenomic_CladeFreq, by= "Internal.external" ) # Join with original dataset
nifD_Metagenomic_CladeFreq$Internal.external<-as.factor(nifD_Metagenomic_CladeFreq$Internal.external)




##### moisaic plots for primers and clades #####
#preparing dataframe
nifD_CladeFreq_1cladeFreq<-filter(nifD_Metagenomic_CladeFreq, Clade_Frequency>1) # remove single-seqeunce clades
str(nifD_CladeFreq_1cladeFreq)
summary(nifD_CladeFreq_1cladeFreq$Clade)
nifD_CladeFreq_1cladeFreq$Clade<-droplevels(nifD_CladeFreq_1cladeFreq$Clade) # drop unsed clades as factors
nifD_Metagenomic_xtabs<-xtabs(~Internal.external + Clade + Soil.data, data=nifD_CladeFreq_1cladeFreq) # chage table from case form to frequency form

#redable mosaics
mosaic(~ Internal.external + Clade,  data=nifD_CladeFreq_1cladeFreq, shade=TRUE, legend=TRUE, residuals_type="deviance", 
       labeling_args = list(rot_labels = c(top = 90, left = 0)),
       offset_varnames = c(top = 2), offset_labels = c(top = 0.8),
       margins = c(right = 0, bottom =0)) # Mosaic Plot, refer to "vcd_tutorial" for more information
      # OBS: the output of this plot was saved as a .svg file and edited in Cytoscape to generate a more clear figure






#############Permanova ################

#load soil chemistry 
SoilChem <- read.csv("Input/SoilChem.csv") 

# attach soil chemistry to NifD seqeunce data
FUllWithSoil <- full_join(SoilChem, nifD_Metagenomic_CaseForm, by= "Soil.data")


#adjusting data to species table format, with species per legume + soil
FullXtab <- as.data.frame(xtabs(~Clade+ Legume.species + Soil.data, data=FUllWithSoil)) # Each combination of legume and soil, for each species, forms a line in this table
FUllWithSoil_spread <- spread(FullXtab,Clade, Freq ) # clades forms columns, creating the species table
FUllWithSoil_spread_united <-unite(FUllWithSoil_spread, LegumePlusSoil, Legume.species, Soil.data, sep = " ") # Join soil and legume data on the same column

# set factors for PERMANOVA
OrdiantionFactos <- select(FUllWithSoil_spread,Legume.species,Soil.data) 
OrdiantionRownames <- unite(OrdiantionFactos, LegumePlusSoil, Legume.species, Soil.data, sep = " ")
OrdiantionFactos <- bind_cols(OrdiantionFactos,OrdiantionRownames)%>%
  column_to_rownames(var = "LegumePlusSoil") # now rownames ahve soil+legume; legume.species shows legumes, soil.data shows soil location

FUllWithSoil_spread_united<- column_to_rownames(FUllWithSoil_spread_united, var = "LegumePlusSoil") # Final species table format with metadata
str(FUllWithSoil_spread_united)
summary(FUllWithSoil_spread_united)

#remove rows with zero occurences (non-inoformative), and then remove unsed factors based on this removal
FUllWithSoil_spread_united<- FUllWithSoil_spread_united[which(rowSums(FUllWithSoil_spread_united) >= 1),]
OrdiantionFactos <-merge(OrdiantionFactos,FUllWithSoil_spread_united, by="row.names")%>%
  column_to_rownames(var = "Row.names")%>%
  select(Legume.species, Soil.data)

#check range to decide a transofrmation
range(FUllWithSoil_spread_united)
range(FUllWithSoil_spread_united^0.5) # use square root transformation to keep range between 0 and 10
range(FUllWithSoil_spread_united^0.25)

#creates distance matrix
DistaMatrix<-vegdist(FUllWithSoil_spread_united^0.5, method = "bray")


#permutation tests for LEGUME+SOIL##############
# legume*soil interactions could not be tested as there were soils with only one legume and legumes with only one soil

#check ebta dispersion
bd<-betadisper(DistaMatrix, OrdiantionFactos$Legume.species,Soil.data) #quantify beta dispersion, akin to multivariate leneve test for homogeniety of variances
anova(bd)# test homogeniety of variances, akin to a multivariate leneve test for homogeniety of variances
permutest(bd)# test homogeniety of variances, akin to a multivariate leneve test for homogeniety of variances

#run and check PERMANOVA
pmv<-adonis(FUllWithSoil_spread_united^0.5 ~ Soil.data + Legume.species, data=OrdiantionFactos, permutations=9999, method="bray")
pmv
##RESULT:###### valid and Significant soil + LEGUME  tests# THERE ARE SOIL EFFECTS, INDEPENDENT OF LEGUME EFFECTS


############## NMDS ################
#CalculateNMDS NMDS NifD sequences slit by soil and legume
NifDNMDS<-metaMDS (FUllWithSoil_spread_united) #metaMDS automatically square roots transforms and uses bray-curtis
stressplot(NifDNMDS) # check stress
mean(goodness(NifDNMDS))# check stress
plot(NifDNMDS) #circle = site, red cross = sp
ordiplot(NifDNMDS,type="n")
orditorp(NifDNMDS,display="species",col="red",air=0.01, cex=0.7, pcex=0)
orditorp(NifDNMDS,display="sites",cex=0.7,air=0.01,pcex=0)
#RESULT: valid stress, satisfactory solution for the data ; greatest differences are in the clades of soils E and F

################PCA OF SOIL DATA#############
#PCA of soil data
soil_PCA<-column_to_rownames(SoilChem,var="Soil.data") # switchgr rownames so that PCA object is only numeric

autoplot(prcomp(soil_PCA, scale= TRUE, center=TRUE), # calculates principal compenents and pltos with ggplot2
         data= soil_PCA, label = TRUE, label.size = 5, shape = FALSE, # add metadate, labels of objects
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 5) # add metadate, labels variables



