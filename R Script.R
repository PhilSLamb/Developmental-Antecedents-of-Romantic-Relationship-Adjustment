#Clear workspace
rm(list = ls())

#Load packages
library(mitools)
library(miceadds)
library(labelled)
library(mice)
library(mitml)
library(haven)
library(psych)
library(performance)
library(skimr)
library(dplyr)
library(apaTables)
library(semTools)

data <- read_sav(file.choose())
RomRel_temp = subset(data, SHINE_505_in_relshp == 1)

#Checking coding of demographic variables
table(RomRel_temp$CHWHITENH)
summary(RomRel_temp$Income_Needs_composite)
summary(RomRel_temp$MEDUCM01)
table(RomRel_temp$CSEX_M01)

#Recoding sex assigned at birth to 0's and 1's and as factor
RomRel_temp$sex <- RomRel_temp$CSEX_M01 - 1 #subtracting 1 from each value
table(RomRel_temp$sex)

#Renaming other demographic variables (race, income-to-needs, maternal education) for clarity/ease of use and recode race as factor
RomRel_temp$race = RomRel_temp$CHWHITENH
RomRel_temp$income = RomRel_temp$Income_Needs_composite
RomRel_temp$mat_edu = RomRel_temp$MEDUCM01

#Creating committed versus non-committed relationship variable
head(RomRel_temp$dasr1a)
table(RomRel_temp$dasr1a)

RomRel_temp$commit = ifelse(RomRel_temp$dasr1a == "a", 0,
                            ifelse(RomRel_temp$dasr1a != "a", 1, 
                                   RomRel_temp$dasr1a))

table(RomRel_temp$commit)
RomRel_temp$commit = as.numeric(RomRel_temp$commit)


#Subset dataset to what we need for analyses
names(RomRel_temp)
RomRel_temp1 = RomRel_temp[,c(1,3:5,10:13,97:101,63:94)]
names(RomRel_temp1)

#Removing labels (due to error message when running mice) and ID for multiple imputation
RomRel_Numeric=RomRel_temp1[,c(2:45)]
RomRel_Numeric=sapply(RomRel_Numeric,haven::zap_labels)
RomRel_Numeric=as.data.frame(RomRel_Numeric)
RomRel = RomRel_Numeric

#Rename to create copy for multiple imputation
romrel2=RomRel
names(romrel2)

  # number of X columns / X variables
  n_Xcols = 12
  
  # number of Y columns / Y variables
  n_Ycols = 32
  
  # Assumption: they come in blocks ... 
  # in particular, X vars comes first, then we have Ycols  last
  
  # Creating a list of WHICH columns should be imputed (only y variables) and which one should not (only X variables)
  imp_method_list = c(rep("", n_Xcols),  rep("pmm", n_Ycols))
  
#imputation
  RomRel_Imputed<-mice(romrel2, m=5, maxit = 50, method = imp_method_list, seed = 500)  
  

#Check data from the imputed datasets
summary(complete(RomRel_Imputed,1))
summary(complete(RomRel_Imputed,2)) 
summary(complete(RomRel_Imputed,3)) 
summary(complete(RomRel_Imputed,4)) 
summary(complete(RomRel_Imputed,5)) 


#Sum items to create a new column "DAS_Total"
long <- complete(RomRel_Imputed, action='long', include=TRUE)

# Generate new variable
names(long)
long$DAS_Total <- rowSums(long[,15:46])
long$DAS_Total1 <- long$DAS_Total/10

# Convert back to Mids
RomRel_Imputed2 <- as.mids(long)

#Descriptive stats for calculating M and SD of DAS_Total

descriptive <- with(RomRel_Imputed2, expr=c("M"=mean(DAS_Total, na.rm=T), "SD"=sd(DAS_Total, na.rm=T) ) )
descriptive1 <- with(RomRel_Imputed2, expr=c("M1"=mean(DAS_Total1, na.rm=T), "SD1"=sd(DAS_Total1, na.rm=T) ) )

# pool estimates
withPool_MI(descriptive)
withPool_MI(descriptive1)

################
##CORRELATIONS##
################


##Correlation Matrix (only for DAS_Total, the rest is done in SPSS since no MI was done on other variables)##

names(RomRel_Imputed2$data) #Column numbers

#Pooling Correlations using the "miceadds" package
micombine.cor(RomRel_Imputed2, variables=c(1:12,45), conf.level=0.95,
              method="pearson", nested=FALSE, partial=NULL)


################
##REGRESSIONS###
################

##REGRESSION ANALYSES##

#Block 1: Covariates --> DAS

reg.model1 = 'DAS_Total ~ sex + race + income + mat_edu
             '

out1 <- runMI(reg.model1, 
              data=RomRel_Imputed2,
              fun="sem", missing = "fiml.x")
summary(out1, standardized = TRUE, rsquare = TRUE)


#Block 2: Interpersonal composites + covariates --> DAS

reg.model2 = 'DAS_Total ~ MotherEIComposite + MatSensComposite + FQQ_Composite + sex + race + income + mat_edu
             '

out2 <- runMI(reg.model2, 
              data=RomRel_Imputed2,
              fun="sem", missing = "fiml.x")
summary(out2, standardized = TRUE, rsquare = TRUE)


#######################
##MODERATION ANALYSES##
#######################

mod.model1 = 'DAS_Total ~ MotherEIComposite + MatSensComposite + FQQ_Composite + sex + race + income + mat_edu + MotherEIComposite:sex + MatSensComposite:sex + FQQ_Composite:sex
             '

outmod1 <- runMI(mod.model1, 
              data=RomRel_Imputed2,
              fun="sem", missing = "fiml.x")
summary(outmod1, standardized = TRUE, rsquare = TRUE)



mod.model2 = 'DAS_Total ~ MotherEIComposite + MatSensComposite + FQQ_Composite + sex + race + income + mat_edu + MotherEIComposite:race + MatSensComposite:race + FQQ_Composite:race
             '

outmod2 <- runMI(mod.model2, 
              data=RomRel_Imputed2,
              fun="sem", missing = "fiml.x")
summary(outmod2, standardized = TRUE, rsquare = TRUE)


mod.model3 = 'DAS_Total ~ MotherEIComposite + MatSensComposite + FQQ_Composite + sex + race + income + mat_edu + MotherEIComposite:income + MatSensComposite:income + FQQ_Composite:income
             '

outmod3 <- runMI(mod.model3, 
                 data=RomRel_Imputed2,
                 fun="sem", missing = "fiml.x")
summary(outmod3, standardized = TRUE, rsquare = TRUE)


mod.model4 = 'DAS_Total ~ MotherEIComposite + MatSensComposite + FQQ_Composite + sex + race + income + mat_edu + MotherEIComposite:mat_edu + MatSensComposite:mat_edu + FQQ_Composite:mat_edu
             '

outmod4 <- runMI(mod.model4, 
                 data=RomRel_Imputed2,
                 fun="sem", missing = "fiml.x")
summary(outmod4, standardized = TRUE, rsquare = TRUE)


mod.model5 = 'DAS_Total ~ MotherEIComposite + MatSensComposite + FQQ_Composite + sex + race + income + mat_edu + commit + MotherEIComposite:commit + MatSensComposite:commit + FQQ_Composite:commit
             '

outmod5 <- runMI(mod.model5, 
                 data=RomRel_Imputed2,
                 fun="sem", missing = "fiml.x")
summary(outmod5, standardized = TRUE, rsquare = TRUE)



##MEDIATION ANALYSES## 

#examine early SEN predicting DAS via all three later measures (later SEN, later PAIR, and FQQ)

set.seed(500)

med.model1 = 'DAS_Total1 ~ c*earlymatsencomposite 
             latermatsencomposite ~ a1*earlymatsencomposite 
             latereicomposite ~ a2*earlymatsencomposite 
             FQQ_Composite ~ a3*earlymatsencomposite 
             DAS_Total1 ~ b1*latermatsencomposite + b2*latereicomposite + b3*FQQ_Composite 
             a1b1 := a1*b1
             a2b2 := a2*b2
             a3b3 := a3*b3
             '
outmed1 <- runMI(med.model1, 
              data=RomRel_Imputed2,
              fun="sem", missing = "fiml.x")

summary(outmed1, standardized = TRUE)
monteCarloCI(outmed1, append.samples = T)
  

#examine early PAIR predicting DAS via all three later measures (later SEN, later PAIR, and FQQ)

med.model2 = 'DAS_Total1 ~ c*earlyeicomposite 
             latermatsencomposite ~ a1*earlyeicomposite 
             latereicomposite ~ a2*earlyeicomposite 
             FQQ_Composite ~ a3*earlyeicomposite 
             DAS_Total1 ~ b1*latermatsencomposite + b2*latereicomposite + b3*FQQ_Composite 
             a1b1 := a1*b1
             a2b2 := a2*b2
             a3b3 := a3*b3
             '
outmed2 <- runMI(med.model2, 
                 data=RomRel_Imputed2,
                 fun="sem", missing = "fiml.x")

summary(outmed2, standardized = TRUE)
monteCarloCI(outmed2, append.samples = T)




