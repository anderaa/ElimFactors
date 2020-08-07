#install.packages("MASS")
library(MASS)#stepAIC

##
##   Make sure you have binData loaded into memory (from Script 2) 
##
binData <- read.csv('binData.csv')
head(binData)

#set up a basic first order logistic regression model
mdl.log <-glm(Observed~Population+Fecundity+Infectivity+
                Exposure_survival + Initial_immunity+
                In_migration+Incidence,
              family =binomial(link = 'logit'),data=binData)
summary(mdl.log)

#include quadratic and cubic terms without interactions (interactions are added by the step-model)
mdl.log.2 <-glm(Observed~Population+I(Population^2)+I(Population^3)+
                  Fecundity+I(Fecundity^2)+I(Fecundity^3)+
                  Infectivity+I(Infectivity^2)+I(Infectivity^3)+
                  Exposure_survival+I(Exposure_survival^2)+I(Exposure_survival^3)+
                  Initial_immunity+I(Initial_immunity^2)+I(Initial_immunity^3)+
                  Incidence+I(Incidence^2)+I(Incidence^3)+
                  In_migration+I(In_migration^2)+I(In_migration^3),
                family =binomial(link = 'logit'),data=binData)
summary(mdl.log.2)


####### step model #################################

####################################################
#                                                  #
#                  WARNING!!!                      #
#                                                  #
# Running the next line will take weeks to compute #                                                 
#                                                  #
#                                                  #    
####################################################
 
mdl.log.step <- stepAIC(mdl.log.2,scope=.~.^2)
binData$Predicted <- predict(mdl.log.step,binData,type="response")

write.csv(binData, 'binDataObsPred.csv')
