# Analysis Script

# Import data from simulations.
dependentVars   <- read.csv('simulationData.csv')[, -1]
explanatoryVars <- read.csv('paramCombos.csv')[, -1]

head(dependentVars)
head(dependentVars)

########################################################################################################################
# These are the basic models with p(elim given dogs remain) as the dependent variable.

# after 1 year
model1 <- glm(cbind(dependentVars[, 'qElimWithDogs1'], dependentVars[, 'qNoElimWithDogs1']) ~
                explanatoryVars[, 'initialPopSize'] +
                explanatoryVars[, 'varianceFactor'] +
                explanatoryVars[, 'fecundityFactor'] +
                explanatoryVars[, 'survivalProb'] +
                explanatoryVars[, 'initialFracImmune'] + 
                explanatoryVars[, 'inMigrationFactor'],
              family=binomial (link='logit'))
summary(model1)

# after 2 years
model2 <- glm(cbind(dependentVars[, 'qElimWithDogs2'], dependentVars[, 'qNoElimWithDogs2']) ~
                explanatoryVars[, 'initialPopSize'] +
                explanatoryVars[, 'varianceFactor'] +
                explanatoryVars[, 'fecundityFactor'] +
                explanatoryVars[, 'survivalProb'] +
                explanatoryVars[, 'initialFracImmune'] + 
                explanatoryVars[, 'inMigrationFactor'],
              family=binomial (link='logit'))
summary(model2)

# after 3 years
model3 <- glm(cbind(dependentVars[, 'qElimWithDogs3'], dependentVars[, 'qNoElimWithDogs3']) ~
                explanatoryVars[, 'initialPopSize'] +
                explanatoryVars[, 'varianceFactor'] +
                explanatoryVars[, 'fecundityFactor'] +
                explanatoryVars[, 'survivalProb'] +
                explanatoryVars[, 'initialFracImmune'] + 
                explanatoryVars[, 'inMigrationFactor'],
              family=binomial (link='logit'))
summary(model3)

# after 4 years
model4 <- glm(cbind(dependentVars[, 'qElimWithDogs4'], dependentVars[, 'qNoElimWithDogs4']) ~
                explanatoryVars[, 'initialPopSize'] +
                explanatoryVars[, 'varianceFactor'] +
                explanatoryVars[, 'fecundityFactor'] +
                explanatoryVars[, 'survivalProb'] +
                explanatoryVars[, 'initialFracImmune'] + 
                explanatoryVars[, 'inMigrationFactor'],
              family=binomial (link='logit'))
summary(model4)

# after 5 years
model5 <- glm(cbind(dependentVars[, 'qElimWithDogs5'], dependentVars[, 'qNoElimWithDogs5']) ~
                explanatoryVars[, 'initialPopSize'] +
                explanatoryVars[, 'varianceFactor'] +
                explanatoryVars[, 'fecundityFactor'] +
                explanatoryVars[, 'survivalProb'] +
                explanatoryVars[, 'initialFracImmune'] + 
                explanatoryVars[, 'inMigrationFactor'],
              family=binomial (link='logit'))
summary(model5)

########################################################################################################################
# This function returns the fitted value (i.e. the predicited probability of elimination). The estimtated
# model and a set of data points must be supplied as arguments. See use example below the function.
getFitted1to5 <- function(model, initialPopSize, varianceFactor, fecundityFactor,
                          survivalProb, initialFracImmune, inMigrationFactor) {
  # get the estimated coefficients
  coef <- coef(model)
  # get the linear combo (x'b)
  linearCombo <- coef[1] + coef[2] * initialPopSize + coef[3] * varianceFactor + 
                           coef[4] * fecundityFactor + coef[5] * survivalProb + 
                           coef[6] * initialFracImmune + coef[7] * inMigrationFactor
  # remove the name that we don't need
  linearCombo <- unname(linearCombo)
  # feed linear combo thru logistic cdf to get fitted probability
  prob = plogis(linearCombo)
  return(prob)
}

# use example - change model and values to desired
pElimination <- getFitted1to5(model=model1, initialPopSize=500, varianceFactor=1.5, fecundityFactor=0.5,
                              survivalProb=0, initialFracImmune=0, inMigrationFactor=0.05)
pElimination
########################################################################################################################


########################################################################################################################
# This function returns the marginal effects of each variable in the model. The marginal effect of a variable is
# dP/dx where P is the fitted probability and x is the variable. So it is the change in the fitted probability in
# response to a small change in the variable. We calculate it by getting the fitted probability when all variables
# are at there means. Then we increase a specific variable by a small amount and get the resulting change in P.
# Dividing the change in P by the amount we increased the variable yields the marginal effect. Note that our model is
# non-linear, so the marginal effects are not constant. We are calculating their values near the mean of the data, but
# they would be different at other locations. There is also an alterative caclulation that provides slightly different
# answeres that we can talk about.

getMarginalEffects1to5 <- function(model) {
  # get mean of each explantory variable
  means <- apply(explanatoryVars, 2, mean)
  # get fitted values at means...use previous function
  fittedAtMeans <- getFitted1to5(model, initialPopSize=means[1], varianceFactor=means[2], 
                                 fecundityFactor=means[3], survivalProb=means[4], 
                                 initialFracImmune=means[5], inMigrationFactor=means[6]) 
  # make a place to store the marginal effects we will calculate
  marginalEffects <- means
  for(i in 1:length(means)) {
    # increase the ith variable by a small amount
    newValues <- means
    increment <- 0.01
    newValues[i] <- newValues[i] + increment
    # get the new fitted value
    newFitted <- getFitted1to5(model, initialPopSize=newValues[1], varianceFactor=newValues[2], 
                               fecundityFactor=newValues[3], survivalProb=newValues[4], 
                               initialFracImmune=newValues[5], inMigrationFactor=newValues[6]) 
    # get the change in fitted values as a result of small change in the ith variable
    changeInFitted <- newFitted - fittedAtMeans
    # calc the ith marginal effect
    marginalEffects[i] <- changeInFitted / increment
  }
  return(marginalEffects)
}

# use example - change model to desired
getMarginalEffects1to5(model1)
########################################################################################################################


########################################################################################################################
# These models are same as above except they include an interaction variable: population size * variance of bites

# after 1 year
model6 <- glm(cbind(dependentVars[, 'qElimWithDogs1'], dependentVars[, 'qNoElimWithDogs1']) ~
                explanatoryVars[, 'initialPopSize'] +
                explanatoryVars[, 'varianceFactor'] +
                I(explanatoryVars[, 'initialPopSize'] * explanatoryVars[, 'varianceFactor']) +
                explanatoryVars[, 'fecundityFactor'] +
                explanatoryVars[, 'survivalProb'] +
                explanatoryVars[, 'initialFracImmune'] + 
                explanatoryVars[, 'inMigrationFactor'],
              family=binomial (link='logit'))
summary(model6)

# after 2 years
model7 <- glm(cbind(dependentVars[, 'qElimWithDogs2'], dependentVars[, 'qNoElimWithDogs2']) ~
                explanatoryVars[, 'initialPopSize'] +
                explanatoryVars[, 'varianceFactor'] +
                I(explanatoryVars[, 'initialPopSize'] * explanatoryVars[, 'varianceFactor']) +
                explanatoryVars[, 'fecundityFactor'] +
                explanatoryVars[, 'survivalProb'] +
                explanatoryVars[, 'initialFracImmune'] + 
                explanatoryVars[, 'inMigrationFactor'],
              family=binomial (link='logit'))
summary(model7)

# after 3 years
model8 <- glm(cbind(dependentVars[, 'qElimWithDogs3'], dependentVars[, 'qNoElimWithDogs3']) ~
                explanatoryVars[, 'initialPopSize'] +
                explanatoryVars[, 'varianceFactor'] +
                I(explanatoryVars[, 'initialPopSize'] * explanatoryVars[, 'varianceFactor']) +
                explanatoryVars[, 'fecundityFactor'] +
                explanatoryVars[, 'survivalProb'] +
                explanatoryVars[, 'initialFracImmune'] + 
                explanatoryVars[, 'inMigrationFactor'],
              family=binomial (link='logit'))
summary(model8)

# after 4 years
model9 <- glm(cbind(dependentVars[, 'qElimWithDogs4'], dependentVars[, 'qNoElimWithDogs4']) ~
                explanatoryVars[, 'initialPopSize'] +
                explanatoryVars[, 'varianceFactor'] +
                I(explanatoryVars[, 'initialPopSize'] * explanatoryVars[, 'varianceFactor']) +
                explanatoryVars[, 'fecundityFactor'] +
                explanatoryVars[, 'survivalProb'] +
                explanatoryVars[, 'initialFracImmune'] + 
                explanatoryVars[, 'inMigrationFactor'],
              family=binomial (link='logit'))
summary(model9)

# after 5 years
model10 <- glm(cbind(dependentVars[, 'qElimWithDogs5'], dependentVars[, 'qNoElimWithDogs5']) ~
                explanatoryVars[, 'initialPopSize'] +
                explanatoryVars[, 'varianceFactor'] +
                I(explanatoryVars[, 'initialPopSize'] * explanatoryVars[, 'varianceFactor']) +
                explanatoryVars[, 'fecundityFactor'] +
                explanatoryVars[, 'survivalProb'] +
                explanatoryVars[, 'initialFracImmune'] + 
                explanatoryVars[, 'inMigrationFactor'],
              family=binomial (link='logit'))
summary(model10)

########################################################################################################################
# This function returns the fitted value (i.e. the predicited probability of elimination). The estimtated
# model and a set of data points must be supplied as arguments. See use example below the function.
getFitted6to10 <- function(model, initialPopSize, varianceFactor, fecundityFactor,
                          survivalProb, initialFracImmune, inMigrationFactor) {
  # get the estimated coefficients
  coef <- coef(model)
  # get the linear combo (x'b)
  linearCombo <- coef[1] + coef[2] * initialPopSize + coef[3] * varianceFactor + 
    coef[4] * (initialPopSize * varianceFactor) + coef[5] * fecundityFactor + coef[6] * survivalProb + 
    coef[7] * initialFracImmune + coef[8] * inMigrationFactor
  # remove the name that we don't need
  linearCombo <- unname(linearCombo)
  # feed linear combo thru logistic cdf to get fitted probability
  prob = plogis(linearCombo)
  return(prob)
}

# use example - change model and values to desired
pElimination <- getFitted6to10(model=model6, initialPopSize=500, varianceFactor=1.5, fecundityFactor=0.5,
                              survivalProb=0, initialFracImmune=0, inMigrationFactor=0.05)
pElimination
########################################################################################################################


########################################################################################################################
# This function returns the marginal effects of each variable in the model. The marginal effect of a variable is
# dP/dx where P is the fitted probability and x is the variable. So it is the change in the fitted probability in
# response to a small change in the variable. We calculate it by getting the fitted probability when all variables
# are at there means. Then we increase a specific variable by a small amount and get the resulting change in P.
# Dividing the change in P by the amount we increased the variable yields the marginal effect. Note that our model is
# non-linear, so the marginal effects are not constant. We are calculating their values near the mean of the data, but
# they would be different at other locations. There is also an alterative caclulation that provides slightly different
# answeres that we can talk about.

getMarginalEffects6to10 <- function(model) {
  # get mean of each explantory variable
  means <- apply(explanatoryVars, 2, mean)
  # get fitted values at means...use previous function
  fittedAtMeans <- getFitted6to10(model, initialPopSize=means[1], varianceFactor=means[2], 
                                 fecundityFactor=means[3], survivalProb=means[4], 
                                 initialFracImmune=means[5], inMigrationFactor=means[6]) 
  # make a place to store the marginal effects we will calculate
  marginalEffects <- means
  for(i in 1:length(means)) {
    # increase the ith variable by a small amount
    newValues <- means
    increment <- 0.01
    newValues[i] <- newValues[i] + increment
    # get the new fitted value
    newFitted <- getFitted6to10(model, initialPopSize=newValues[1], varianceFactor=newValues[2], 
                               fecundityFactor=newValues[3], survivalProb=newValues[4], 
                               initialFracImmune=newValues[5], inMigrationFactor=newValues[6]) 
    # get the change in fitted values as a result of small change in the ith variable
    changeInFitted <- newFitted - fittedAtMeans
    # calc the ith marginal effect
    marginalEffects[i] <- changeInFitted / increment
  }
  return(marginalEffects)
}

# use example - change model to desired
getMarginalEffects6to10(model6)
########################################################################################################################






########################################################################################################################
########################################################################################################################
########################################################################################################################





# SAME AS ABOVE EXCEPT "SUCCESS" = FADEOUT _OR_ EXTINCTION
########################################################################################################################
# These are the basic models with p(elim) as the dependent variable.

# after 1 year
model11 <- glm(cbind(dependentVars[, 'qElim1'], dependentVars[, 'qNoElimWithDogs1']) ~
                explanatoryVars[, 'initialPopSize'] +
                explanatoryVars[, 'varianceFactor'] +
                explanatoryVars[, 'fecundityFactor'] +
                explanatoryVars[, 'survivalProb'] +
                explanatoryVars[, 'initialFracImmune'] + 
                explanatoryVars[, 'inMigrationFactor'],
              family=binomial (link='logit'))
summary(model11)

# after 2 years
model12 <- glm(cbind(dependentVars[, 'qElim2'], dependentVars[, 'qNoElimWithDogs2']) ~
                explanatoryVars[, 'initialPopSize'] +
                explanatoryVars[, 'varianceFactor'] +
                explanatoryVars[, 'fecundityFactor'] +
                explanatoryVars[, 'survivalProb'] +
                explanatoryVars[, 'initialFracImmune'] + 
                explanatoryVars[, 'inMigrationFactor'],
              family=binomial (link='logit'))
summary(model12)

# after 3 years
model13 <- glm(cbind(dependentVars[, 'qElim3'], dependentVars[, 'qNoElimWithDogs3']) ~
                explanatoryVars[, 'initialPopSize'] +
                explanatoryVars[, 'varianceFactor'] +
                explanatoryVars[, 'fecundityFactor'] +
                explanatoryVars[, 'survivalProb'] +
                explanatoryVars[, 'initialFracImmune'] + 
                explanatoryVars[, 'inMigrationFactor'],
              family=binomial (link='logit'))
summary(model13)

# after 4 years
model14 <- glm(cbind(dependentVars[, 'qElim4'], dependentVars[, 'qNoElimWithDogs4']) ~
                explanatoryVars[, 'initialPopSize'] +
                explanatoryVars[, 'varianceFactor'] +
                explanatoryVars[, 'fecundityFactor'] +
                explanatoryVars[, 'survivalProb'] +
                explanatoryVars[, 'initialFracImmune'] + 
                explanatoryVars[, 'inMigrationFactor'],
              family=binomial (link='logit'))
summary(model14)

# after 5 years
model15 <- glm(cbind(dependentVars[, 'qElim5'], dependentVars[, 'qNoElimWithDogs5']) ~
                explanatoryVars[, 'initialPopSize'] +
                explanatoryVars[, 'varianceFactor'] +
                explanatoryVars[, 'fecundityFactor'] +
                explanatoryVars[, 'survivalProb'] +
                explanatoryVars[, 'initialFracImmune'] + 
                explanatoryVars[, 'inMigrationFactor'],
              family=binomial (link='logit'))
summary(model15)

########################################################################################################################
# This function returns the fitted value (i.e. the predicited probability of elimination). The estimtated
# model and a set of data points must be supplied as arguments. See use example below the function.
getFitted11to15 <- function(model, initialPopSize, varianceFactor, fecundityFactor,
                          survivalProb, initialFracImmune, inMigrationFactor) {
  # get the estimated coefficients
  coef <- coef(model)
  # get the linear combo (x'b)
  linearCombo <- coef[1] + coef[2] * initialPopSize + coef[3] * varianceFactor + 
    coef[4] * fecundityFactor + coef[5] * survivalProb + 
    coef[6] * initialFracImmune + coef[7] * inMigrationFactor
  # remove the name that we don't need
  linearCombo <- unname(linearCombo)
  # feed linear combo thru logistic cdf to get fitted probability
  prob = plogis(linearCombo)
  return(prob)
}

# use example - change model and values to desired
pElimination <- getFitted11to15(model=model11, initialPopSize=500, varianceFactor=1.5, fecundityFactor=0.5,
                              survivalProb=0, initialFracImmune=0, inMigrationFactor=0.05)
pElimination
########################################################################################################################


########################################################################################################################
# This function returns the marginal effects of each variable in the model. The marginal effect of a variable is
# dP/dx where P is the fitted probability and x is the variable. So it is the change in the fitted probability in
# response to a small change in the variable. We calculate it by getting the fitted probability when all variables
# are at there means. Then we increase a specific variable by a small amount and get the resulting change in P.
# Dividing the change in P by the amount we increased the variable yields the marginal effect. Note that our model is
# non-linear, so the marginal effects are not constant. We are calculating their values near the mean of the data, but
# they would be different at other locations. There is also an alterative caclulation that provides slightly different
# answeres that we can talk about.

getMarginalEffects11to15 <- function(model) {
  # get mean of each explantory variable
  means <- apply(explanatoryVars, 2, mean)
  # get fitted values at means...use previous function
  fittedAtMeans <- getFitted11to15(model, initialPopSize=means[1], varianceFactor=means[2], 
                                 fecundityFactor=means[3], survivalProb=means[4], 
                                 initialFracImmune=means[5], inMigrationFactor=means[6]) 
  # make a place to store the marginal effects we will calculate
  marginalEffects <- means
  for(i in 1:length(means)) {
    # increase the ith variable by a small amount
    newValues <- means
    increment <- 0.01
    newValues[i] <- newValues[i] + increment
    # get the new fitted value
    newFitted <- getFitted11to15(model, initialPopSize=newValues[1], varianceFactor=newValues[2], 
                                 fecundityFactor=newValues[3], survivalProb=newValues[4], 
                                 initialFracImmune=newValues[5], inMigrationFactor=newValues[6]) 
    # get the change in fitted values as a result of small change in the ith variable
    changeInFitted <- newFitted - fittedAtMeans
    # calc the ith marginal effect
    marginalEffects[i] <- changeInFitted / increment
  }
  return(marginalEffects)
}

# use example - change model to desired
getMarginalEffects11to15(model11)
########################################################################################################################


########################################################################################################################
# These models are same as above except they include an interaction variable: population size * variance of bites

# after 1 year
model16 <- glm(cbind(dependentVars[, 'qElim1'], dependentVars[, 'qNoElimWithDogs1']) ~
                explanatoryVars[, 'initialPopSize'] +
                explanatoryVars[, 'varianceFactor'] +
                I(explanatoryVars[, 'initialPopSize'] * explanatoryVars[, 'varianceFactor']) +
                explanatoryVars[, 'fecundityFactor'] +
                explanatoryVars[, 'survivalProb'] +
                explanatoryVars[, 'initialFracImmune'] + 
                explanatoryVars[, 'inMigrationFactor'],
              family=binomial (link='logit'))
summary(model16)

# after 2 years
model17 <- glm(cbind(dependentVars[, 'qElim2'], dependentVars[, 'qNoElimWithDogs2']) ~
                explanatoryVars[, 'initialPopSize'] +
                explanatoryVars[, 'varianceFactor'] +
                I(explanatoryVars[, 'initialPopSize'] * explanatoryVars[, 'varianceFactor']) +
                explanatoryVars[, 'fecundityFactor'] +
                explanatoryVars[, 'survivalProb'] +
                explanatoryVars[, 'initialFracImmune'] + 
                explanatoryVars[, 'inMigrationFactor'],
              family=binomial (link='logit'))
summary(model17)

# after 3 years
model18 <- glm(cbind(dependentVars[, 'qElim3'], dependentVars[, 'qNoElimWithDogs3']) ~
                explanatoryVars[, 'initialPopSize'] +
                explanatoryVars[, 'varianceFactor'] +
                I(explanatoryVars[, 'initialPopSize'] * explanatoryVars[, 'varianceFactor']) +
                explanatoryVars[, 'fecundityFactor'] +
                explanatoryVars[, 'survivalProb'] +
                explanatoryVars[, 'initialFracImmune'] + 
                explanatoryVars[, 'inMigrationFactor'],
              family=binomial (link='logit'))
summary(model18)

# after 4 years
model19 <- glm(cbind(dependentVars[, 'qElim4'], dependentVars[, 'qNoElimWithDogs4']) ~
                explanatoryVars[, 'initialPopSize'] +
                explanatoryVars[, 'varianceFactor'] +
                I(explanatoryVars[, 'initialPopSize'] * explanatoryVars[, 'varianceFactor']) +
                explanatoryVars[, 'fecundityFactor'] +
                explanatoryVars[, 'survivalProb'] +
                explanatoryVars[, 'initialFracImmune'] + 
                explanatoryVars[, 'inMigrationFactor'],
              family=binomial (link='logit'))
summary(model19)

# after 5 years
model20 <- glm(cbind(dependentVars[, 'qElim5'], dependentVars[, 'qNoElimWithDogs5']) ~
                 explanatoryVars[, 'initialPopSize'] +
                 explanatoryVars[, 'varianceFactor'] +
                 I(explanatoryVars[, 'initialPopSize'] * explanatoryVars[, 'varianceFactor']) +
                 explanatoryVars[, 'fecundityFactor'] +
                 explanatoryVars[, 'survivalProb'] +
                 explanatoryVars[, 'initialFracImmune'] + 
                 explanatoryVars[, 'inMigrationFactor'],
               family=binomial (link='logit'))
summary(model20)

########################################################################################################################
# This function returns the fitted value (i.e. the predicited probability of elimination). The estimtated
# model and a set of data points must be supplied as arguments. See use example below the function.
getFitted16to20 <- function(model, initialPopSize, varianceFactor, fecundityFactor,
                           survivalProb, initialFracImmune, inMigrationFactor) {
  # get the estimated coefficients
  coef <- coef(model)
  # get the linear combo (x'b)
  linearCombo <- coef[1] + coef[2] * initialPopSize + coef[3] * varianceFactor + 
    coef[4] * (initialPopSize * varianceFactor) + coef[5] * fecundityFactor + coef[6] * survivalProb + 
    coef[7] * initialFracImmune + coef[8] * inMigrationFactor
  # remove the name that we don't need
  linearCombo <- unname(linearCombo)
  # feed linear combo thru logistic cdf to get fitted probability
  prob = plogis(linearCombo)
  return(prob)
}

# use example - change model and values to desired
pElimination <- getFitted16to20(model=model16, initialPopSize=500, varianceFactor=1.5, fecundityFactor=0.5,
                               survivalProb=0, initialFracImmune=0, inMigrationFactor=0.05)
pElimination
########################################################################################################################


########################################################################################################################
# This function returns the marginal effects of each variable in the model. The marginal effect of a variable is
# dP/dx where P is the fitted probability and x is the variable. So it is the change in the fitted probability in
# response to a small change in the variable. We calculate it by getting the fitted probability when all variables
# are at there means. Then we increase a specific variable by a small amount and get the resulting change in P.
# Dividing the change in P by the amount we increased the variable yields the marginal effect. Note that our model is
# non-linear, so the marginal effects are not constant. We are calculating their values near the mean of the data, but
# they would be different at other locations. There is also an alterative caclulation that provides slightly different
# answeres that we can talk about.

getMarginalEffects16to20 <- function(model) {
  # get mean of each explantory variable
  means <- apply(explanatoryVars, 2, mean)
  # get fitted values at means...use previous function
  fittedAtMeans <- getFitted16to20(model, initialPopSize=means[1], varianceFactor=means[2], 
                                  fecundityFactor=means[3], survivalProb=means[4], 
                                  initialFracImmune=means[5], inMigrationFactor=means[6]) 
  # make a place to store the marginal effects we will calculate
  marginalEffects <- means
  for(i in 1:length(means)) {
    # increase the ith variable by a small amount
    newValues <- means
    increment <- 0.01
    newValues[i] <- newValues[i] + increment
    # get the new fitted value
    newFitted <- getFitted16to20(model, initialPopSize=newValues[1], varianceFactor=newValues[2], 
                                 fecundityFactor=newValues[3], survivalProb=newValues[4], 
                                 initialFracImmune=newValues[5], inMigrationFactor=newValues[6]) 
    # get the change in fitted values as a result of small change in the ith variable
    changeInFitted <- newFitted - fittedAtMeans
    # calc the ith marginal effect
    marginalEffects[i] <- changeInFitted / increment
  }
  return(marginalEffects)
}

# use example - change model to desired
getMarginalEffects16to20(model16)
########################################################################################################################


# Make a graph
fitted <- rep(NA, nrow(dependentVars))
for(i in 1:nrow(dependentVars)) {
  fitted[i] <- getFitted16to20(model=model18, 
                               initialPopSize=explanatoryVars[i, 'initialPopSize'], 
                               varianceFactor=explanatoryVars[i, 'varianceFactor'], 
                               fecundityFactor=explanatoryVars[i, 'fecundityFactor'],
                               survivalProb=explanatoryVars[i, 'survivalProb'], 
                               initialFracImmune=explanatoryVars[i, 'initialFracImmune'], 
                               inMigrationFactor=explanatoryVars[i, 'inMigrationFactor'])
}
plot(explanatoryVars[explanatoryVars[, 'initialPopSize'] == 20814, 'varianceFactor'], 
     fitted[explanatoryVars[, 'initialPopSize'] == 20814],
     xaxt = 'n', col='green', ylim=c(0, 1),
     xlab = 'varianceFactor', ylab='fitted p(fade-out)')
axis(side=1, at=c(0.5, 0.75, 1, 1.25, 1.5))
points(explanatoryVars[explanatoryVars[, 'initialPopSize'] == 500, 'varianceFactor'], 
       fitted[explanatoryVars[, 'initialPopSize'] == 500], col='red')
legend( "bottomright", legend = c("small pop  ", "large pop"), bty='y', col = c("red", "green"), pch=1, horiz=TRUE)

     

