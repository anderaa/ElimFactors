library('foreach')
library('doParallel')

rm(list=ls())

########################################
# Set Parameter Values

# inputs for simulation
simulationYears <- 1  
iterations      <- 5  

# Build matrix of all combinations of parameters that are varied:
initialPopSizeList    <- 5400
fecundityFactorList   <- seq(0.5, 1.5, length.out=5)
varianceFactorList    <- c(1)  # variability of infectivity (not used) 
meanFactorList        <- seq(0.5, 1.5, length.out=5) # infectivity
survivalProbList      <- seq(0, 0.001, length.out=5)
initialFracImmuneList <- seq(0, 1, length.out=11)
inMigrationFactor     <- seq(0, 0.1, length.out=5)
maxAnnualIncidence    <- 0.03
# representing a 300% increase at each subsequent level
epidemicSizeList <- c(maxAnnualIncidence / 9, maxAnnualIncidence / 3, maxAnnualIncidence) 


paramCombos <- expand.grid(initialPopSizeList, fecundityFactorList, 
                           varianceFactorList, meanFactorList,
                           survivalProbList, initialFracImmuneList, inMigrationFactor, epidemicSizeList)
colnames(paramCombos) <- c('initialPopSize', 'fecundityFactor', 
                           'varianceFactor', 'meanFactor',
                           'survivalProb', 'initialFracImmune', 'inMigrationFactor','epidemicSize')

# Store results here.
qElimination        <- rep(NA, simulationYears)
qNoDogsLeft         <- rep(NA, simulationYears)
qElimWithDogsLeft   <- rep(NA, simulationYears)
qNoElimWithDogsLeft <- rep(NA, simulationYears)

startTime <- proc.time()

registerDoParallel(20)

simulationData <- foreach(c = 1:nrow(paramCombos), .combine=rbind, .inorder=FALSE) %dopar% {
#simulationData <- foreach(c = 68751:68752, .combine=rbind, .inorder=FALSE) %dopar% {   # for testing
  
  ########################################
  getDailyBudget <- function(j) {
    # Arguments: The year of the simulation (j)
    # Return:    The amount that can be spent each day given the management months specified and the annual budget
    # Purpose:   Called at the beginning of each year to allocate the budget across management campaigns during year
    dailyBudget <- mgtDayVector * (annualBudget[j] / sum(mgtDayVector))
    return(dailyBudget)
  }
  ########################################
  
  ########################################
  InitialPopulation <- function() {
    # Agruments: None.
    # Return:    The population matrix.
    # Purpose:   Construct the initial population matrix.
    
    popMatrix <- matrix(0, nrow=initialPopSize, ncol=length(traitList))
    colnames(popMatrix) <- traitList                 
    popMatrix[, 'age'] <- c(sample(seq(1, maxPuppyAge), initialPuppies, replace=TRUE),
                            sample(seq(maxPuppyAge + 1, maxJuvAge), initialJuveniles, replace=TRUE),
                            sample(seq(maxJuvAge + 1, maxAge), initialAdults, replace=TRUE))
    popMatrix[, 'female'] <- sample(c(0, 1), initialPopSize, replace=TRUE, 
                                    prob=c(1-initialFracFemale, initialFracFemale))
    popMatrix[, 'contracepted'] <- sample(c(0, 1), initialPopSize, replace=TRUE, 
                                          prob=c(1-initialFracContra, initialFracContra))
    popMatrix[, 'sterilized'] <- sample(c(0, 1), initialPopSize, replace=TRUE, 
                                        prob=c(1-initialFracSter, initialFracSter))
    popMatrix[, 'vaccinated'] <- sample(c(0, 1), initialPopSize, replace=TRUE, 
                                        prob=c(1-initialFracVacc, initialFracVacc))
    popMatrix[, 'immune'] <- sample(c(0, 1), initialPopSize, replace=TRUE, 
                                    prob=c(1-initialFracImmune, initialFracImmune))
    popMatrix[popMatrix[, 'contracepted']==1, 'timeContra'] <- 
      vapply(popMatrix[popMatrix[, 'contracepted']==1, 'age'], sample, size=1, FUN.VALUE=0) 
    popMatrix[popMatrix[, 'vaccinated']==1, 'timeVacc'] <- 
      vapply(popMatrix[popMatrix[, 'vaccinated']==1, 'age'], sample, size=1, FUN.VALUE=0)
    popMatrix[popMatrix[, 'age'] > maxJuvAge, 'adult'] <- 1
    popMatrix[popMatrix[, 'age'] <= maxPuppyAge, 'puppy'] <- 1
    popMatrix[, 'month'] <- 1
    popMatrix[, 'contactCost'] <- sample(marginalCost, nrow(popMatrix), replace=TRUE)
    popMatrix[, 'timeLimitExposed'] <- rgamma(nrow(popMatrix), shape=exposedTimeShape, rate=exposedTimeRate)
    popMatrix[, 'timeLimitInfective'] <- rgamma(nrow(popMatrix), shape=infectiveTimeShape, rate=infectiveTimeRate)
    
    return(popMatrix)
  }
  ########################################
  
  ########################################
  TotalCost <- function(k, dogDemoForStrategies) {
    # Arguments: Money allocated to contact or capture, number of dogs available for each strategy.
    # Return:    The total cost including contact or capture and treatment costs.
    # Purpose:   Calculates the total cost of management for a given year given the amount devoted to contact or 
    #            capture and current demographics.
    
    currentAbundance <- nrow(popMatrix)
    contactAllocation <- costSequence[k]
    
    # adjust contact cost for current population relative to k
    dogsContacted = contactMapping[k] * currentAbundance
    adjAllocation = contactAllocation * (currentAbundance / carryingCap)
    
    demoOfContacted <- dogsContacted * pmax((dogDemoForStrategies / nrow(popMatrix)),  0, na.rm=TRUE)
    totalCost <- adjAllocation + sum(demoOfContacted * strategyCostVector * strategyVector)
    
    return(totalCost)
  }
  ########################################
  
  ########################################
  AnnualStrategy <- function() {
    # Arguments: None.
    # Return:    A vector of the number of dogs in each demographic category that receive each treament.
    # Purpose:   At the start of each year, this function calculates the number of dogs that can be treated given costs, 
    #            demographics, and the specified strategy.
    
    # Calc the expected demographics associated with each strategy choice:
    currentAbundance   <- nrow(popMatrix)
    pupFemalesUnster   <- sum(popMatrix[, 'female'] == 1 & 
                                popMatrix[, 'puppy'] == 1 & 
                                popMatrix[, 'sterilized'] == 0)
    juvFemalesUnster   <- sum(popMatrix[, 'female'] == 1 & 
                                popMatrix[, 'adult'] == 0 &
                                popMatrix[, 'puppy'] == 0 &
                                popMatrix[, 'sterilized'] == 0)
    adultFemalesUnster <- sum(popMatrix[, 'female'] == 1 & 
                                popMatrix[, 'adult'] == 1 & 
                                popMatrix[, 'sterilized'] == 0)
    pupMalesUnster     <- sum(popMatrix[, 'female'] == 0 & 
                                popMatrix[, 'puppy'] == 1 & 
                                popMatrix[, 'sterilized'] == 0)
    juvMalesUnster     <- sum(popMatrix[, 'female'] == 0 & 
                                popMatrix[, 'adult'] == 0 &
                                popMatrix[, 'puppy'] == 0 &
                                popMatrix[, 'sterilized'] == 0)
    adultMalesUnster   <- sum(popMatrix[, 'female'] == 0 & 
                                popMatrix[, 'adult'] == 1 & 
                                popMatrix[, 'sterilized'] == 0)
    pupFemalesSter     <- sum(popMatrix[, 'female'] == 1 & 
                                popMatrix[, 'puppy'] == 1 &
                                popMatrix[, 'sterilized'] == 1)
    juvFemalesSter     <- sum(popMatrix[, 'female'] == 1 & 
                                popMatrix[, 'adult'] == 0 &
                                popMatrix[, 'puppy'] == 0 &
                                popMatrix[, 'sterilized'] == 1)
    adultFemalesSter   <- sum(popMatrix[, 'female'] == 1 & 
                                popMatrix[, 'adult'] == 1 &
                                popMatrix[, 'sterilized'] == 1)
    pupMalesSter       <- sum(popMatrix[, 'female'] == 0 & 
                                popMatrix[, 'puppy'] == 1 & 
                                popMatrix[, 'sterilized'] == 1)
    juvMalesSter       <- sum(popMatrix[, 'female'] == 0 & 
                                popMatrix[, 'adult'] == 0 &
                                popMatrix[, 'puppy'] == 0 &
                                popMatrix[, 'sterilized'] == 1)
    adultMalesSter     <- sum(popMatrix[, 'female'] == 0 & 
                                popMatrix[, 'adult'] == 1 & 
                                popMatrix[, 'sterilized'] == 1)
    dogDemoForStrategies <- c(pupMalesUnster + pupMalesSter,
                              pupFemalesUnster + pupFemalesSter,
                              adultMalesUnster + adultMalesSter,
                              adultFemalesUnster + adultFemalesSter,
                              juvMalesUnster + juvMalesSter,
                              juvFemalesUnster + juvFemalesSter,
                              pupMalesUnster, pupFemalesUnster,
                              adultMalesUnster, adultFemalesUnster, 
                              juvMalesUnster, juvFemalesUnster,
                              pupMalesUnster, pupFemalesUnster,
                              adultMalesUnster, adultFemalesUnster, 
                              juvMalesUnster, juvFemalesUnster,
                              pupMalesUnster + pupMalesSter,
                              pupFemalesUnster + pupFemalesSter,
                              adultMalesUnster + adultMalesSter,
                              adultFemalesUnster + adultMalesSter,
                              juvMalesUnster + juvMalesSter,
                              juvFemalesUnster + juvFemalesSter)
    names(dogDemoForStrategies) <- strategyNames
    
    # The following section of code determines how the annual budget is split between contacting/capturing dogs and 
    # treating them. It begins by allocating the entire budget to capture. Then it calculates the total cost of this 
    # given the treatment options that have been specified. If the total cost exceeds the available budget, it decreases 
    # the allocation to contact/capture and repeats the calculation. This process continues until the total cost is less 
    # than (or equal to) the total annual budget.
    
    if(annualBudget[j] > 0) {
      k = length(costSequence)
      totalCost = annualBudget[j] + 1
      while(totalCost > annualBudget[j]) {
        allocation = costSequence[k]
        totalCost = TotalCost(k, dogDemoForStrategies)
        k = k - 1
      }
      contactAllocation = costSequence[k + 1]
      # adjust contact cost for current population relative to k
      dogsContacted = contactMapping[k + 1] * 
        currentAbundance
      
      demoOfContacted <- dogsContacted * pmax((dogDemoForStrategies / nrow(popMatrix)), 0, na.rm=TRUE)
      strategyResults <- demoOfContacted * strategyVector
      
    } else {
      strategyResults = dogDemoForStrategies * 0
    }
    
    return(strategyResults) 
  }
  ########################################
  
  ########################################
  MortalityFunction <- function() {
    # Arguments: None.
    # Return:    An updated population matrix.
    # Purpose:   Induces out-migration, probabilistic mortality,
    #            old-age mortality, and censors population to carrying capacity.
    
    # Induce out-migration:
    emigDraw  <- runif(nrow(popMatrix))
    popMatrix <- popMatrix[emigDraw > emigrationProb, , drop=FALSE]
    
    # Induce standard and age-related mortality:
    n <- nrow(popMatrix)
    mortProbVector <- rep(adultMortalityProb, n)
    mortProbVector[popMatrix[, 'age'] <= maxJuvAge]  <- juvMortalityProb
    mortProbVector[popMatrix[, 'age'] <= maxPuppyAge] <- pupMortalityProb
    mortDraw <- runif(n)
    popMatrix <- popMatrix[mortDraw > mortProbVector, , drop=FALSE]
    popMatrix <- popMatrix[maxAge > popMatrix[, 'age'], , drop=FALSE]
    
    # Censor the population to carrying capacity:
    n <- nrow(popMatrix)
    if (n > 0) {
      mortProbVector <- rep(adultMortalityProb, n)
      mortProbVector[popMatrix[, 'age'] <= maxJuvAge]  <- juvMortalityProb
      mortProbVector[popMatrix[, 'age'] <= maxPuppyAge] <- pupMortalityProb
      survProbVector <- 1 - mortProbVector
      survivors <- sample(seq(1, n), min(carryingCap, n), prob=survProbVector, replace=FALSE)
      popMatrix <- popMatrix[survivors, , drop=FALSE]
    }
    
    return(popMatrix)
  }
  ########################################
  
  ########################################
  ReproductionFunction <- function(d) {
    # Arguments: Day of the year.
    # Return:    An updated population matrix.
    # Purpose:   Induces reproduction and adds the new puppies to
    #            the population.
    
    popMatrix[popMatrix[, 'age'] == (maxJuvAge + 1), 'adult'] <- 1
    popMatrix[popMatrix[, 'age'] == (maxPuppyAge +1), 'puppy'] <- 0
    fertFemales <- sum(popMatrix[, 'adult'] == 1 & 
                         popMatrix[, 'female'] == 1 &
                         popMatrix[, 'contracepted'] == 0 &
                         popMatrix[, 'sterilized'] == 0 &
                         (popMatrix[, 'exposed'] + popMatrix[, 'infective'] == 0 | popMatrix[, 'immune'] == 1))
    fertFemales <- max(0, fertFemales)
    litterDraw <- runif(fertFemales)
    puppies <- round(sum(litterDraw < litterProbability[d]) * meanLitterSize)
    newDogMatrix <- matrix(0, nrow=puppies, ncol=length(traitList))
    colnames(newDogMatrix) <- traitList
    newDogMatrix[, 'female'] <- sample(c(0, 1), puppies, 
                                       replace=TRUE, prob=c(1-femalePupProb, 
                                                            femalePupProb))
    newDogMatrix[, 'puppy'] <- 1
    newDogMatrix[, 'contactCost'] <- sample(marginalCost, nrow(newDogMatrix), replace=TRUE)
    newDogMatrix[, 'timeLimitExposed'] <- rgamma(nrow(newDogMatrix), shape=exposedTimeShape, rate=exposedTimeRate)
    newDogMatrix[, 'timeLimitInfective'] <- rgamma(nrow(newDogMatrix), shape=infectiveTimeShape, rate=infectiveTimeRate)
    popMatrix <- rbind(popMatrix, newDogMatrix)
    
    return(popMatrix)
  }
  ########################################
  
  ########################################
  ImmigrationFunction <- function() {
    # Arguments: None.
    # Return:    An updated population matrix.
    # Purpose:  Creates immigrant dogs and adds them to the population.
    
    # Calcuate number of new dogs to be added each day:
    newDogCount <- immigrantDogs%/%365
    probNewDog <- (immigrantDogs%%365) / 365
    newDogCount <- newDogCount + (runif(1) < probNewDog)
    
    # Create and add new dogs:
    newDogMatrix <- matrix(0, nrow=newDogCount, ncol=length(traitList))
    colnames(newDogMatrix) <- traitList
    newDogMatrix[, 'female'] <- sample(c(0, 1), newDogCount, 
                                       replace=TRUE, 
                                       prob=c(1-initialFracFemale, 
                                              initialFracFemale))
    newDogMatrix[, 'age'] <- sample(seq(0, maxAge), newDogCount, replace=TRUE)
    newDogMatrix[newDogMatrix[, 'age'] > maxJuvAge, 'adult'] <- 1
    newDogMatrix[newDogMatrix[, 'age'] <= maxPuppyAge, 'puppy'] <- 1
    newDogMatrix[, 'contactCost'] <- sample(marginalCost, nrow(newDogMatrix), replace=TRUE)
    newDogMatrix[, 'timeLimitExposed'] <- rgamma(nrow(newDogMatrix), shape=exposedTimeShape, rate=exposedTimeRate)
    newDogMatrix[, 'timeLimitInfective'] <- rgamma(nrow(newDogMatrix), shape=infectiveTimeShape, rate=infectiveTimeRate)
    popMatrix <- rbind(popMatrix, newDogMatrix)
    
    return(popMatrix)
  }
  ########################################
  
  ########################################
  DiseaseSpreadFunction <- function() {
    # Arguments: None.
    # Return:    An updated population matrix.
    # Purpose:   Induces disease transmission, exogenous disease introduction.
    
    # Exogenous transmission:
    if (d %in% pressureDays[[j]]) {
      temp <- rep(0, nrow(popMatrix))
      # Generate a vector of 0's and 1's that indicates which individuals are exposed:
      temp[sample(seq(1, nrow(popMatrix)), dogsPerMonthExposed)] <- 1
      # Change states if the individual can be moved to the exposed state:
      newExposed <- temp == 1 & popMatrix[, 'infective'] == 0 & popMatrix[, 'exposed'] == 0 &
        popMatrix[, 'immune'] == 0 & popMatrix[, 'vaccinated'] == 0
      popMatrix[newExposed, 'exposed']     <- 1
      popMatrix[newExposed, 'timeExposed'] <- 0
    }
    
    # Endogenous transmission:
    infectiveTimes <- popMatrix[popMatrix[, 'infective'] == 1, 'timeLimitInfective']
    if (length(infectiveTimes) > 0) { 
      size <- bitesPerRabidShape / infectiveTimes
      mu <- bitesPerRabidMean / infectiveTimes
      biteMatrix <- mapply(function(x, y){rnbinom(size=x, mu=y, n=1)}, x=size, y=mu)
      dailyRabidBites <- sum(biteMatrix)
    } else{
      dailyRabidBites <- 0
    }
    
    # Now we draw dogs randomly from population to be bitten:
    rowsBitten <- unique(sample(seq(1:nrow(popMatrix)), dailyRabidBites, replace=TRUE))
    bitten <- rep(0, nrow(popMatrix))
    bitten[rowsBitten] <- 1
    infectionDraw <- runif(nrow(popMatrix))
    # Treat dog as unbitten if did not actually aquire infection from bite:
    bitten[infectionDraw > probInfectionFromBite] <- 0
    # Take the dogs that received rabid bites and moved to exposed state if appropriate:
    newExposed <- bitten == 1 & popMatrix[, 'infective'] == 0 & popMatrix[, 'exposed'] == 0 &
      popMatrix[, 'immune'] == 0 & popMatrix[, 'vaccinated'] == 0
    popMatrix[newExposed, 'exposed']     <- 1
    popMatrix[newExposed, 'timeExposed'] <- 0
    
    return(popMatrix)
  } 
  ########################################
  
  ########################################
  DiseaseProgressionFunction <- function() {
    # Arguments: None.
    # Return:    An updated population matrix.
    # Purpose:   Induces transition from exposed and infective states.
    
    # Transition exposed to infective:
    newInfective <- popMatrix[, 'exposed'] == 1 & popMatrix[, 'timeExposed'] > popMatrix[, 'timeLimitExposed']
    recoverDraw <- runif(length(newInfective))
    # If a dog is both leaving exposed state and recovering, it gets a TRUE
    recover <- newInfective & recoverDraw < survivalProb
    # if a dog is both leaving the exposed state and moving to infective, it gets a TRUE
    infective <- newInfective & recoverDraw >= survivalProb
    popMatrix[infective, 'exposed']       <- 0
    popMatrix[infective, 'infective']     <- 1
    popMatrix[infective, 'timeInfective'] <- 0
    popMatrix[recover, 'exposed']         <- 0
    popMatrix[recover, 'immune']          <- 1
    
    # Transition infective to death:
    newDead <- popMatrix[, 'infective'] == 1 & popMatrix[, 'timeInfective'] > popMatrix[, 'timeLimitInfective']
    # Dog gets a TRUE if leaving infective, keep dogs with a FALSE
    popMatrix <- popMatrix[!newDead, , drop=FALSE]
    
    return(popMatrix)
  }
  ########################################
  
  ########################################
  CensusFunction <- function() {
    # Arguments: None.
    # Return: A vector of results.
    # Purpose: Calculate results that are recorded daily.
    
    censusVector['abundance'] <- nrow(popMatrix)
    censusVector['puppy'] <- sum(popMatrix[, 'age'] <= maxPuppyAge)
    censusVector['adult'] <- sum(popMatrix[, 'age'] > maxJuvAge)
    censusVector['females'] <- sum(popMatrix[, 'female'])
    censusVector['sterilized'] <- sum(popMatrix[, 'sterilized'])
    censusVector['femalesSterilized'] <- sum(popMatrix[, 'sterilized'] == 1 & popMatrix[, 'female'] == 1)
    censusVector['contracepted'] <- sum(popMatrix[, 'sterilized'])
    censusVector['femalesContracepted'] <- sum(popMatrix[, 'contracepted'] == 1 & 
                                                 popMatrix[, 'female'] == 1)
    censusVector['vaccinated'] <- sum(popMatrix[, 'vaccinated'])
    censusVector['immune'] <- sum(popMatrix[, 'immune'])
    censusVector['exposed'] <- sum(popMatrix[, 'exposed'])
    censusVector['infective'] <- sum(popMatrix[, 'infective'])
    bitesNonRabid <- bitesPerNonRabid * (censusVector['abundance'] - censusVector['infective'])
    bitesRabid <- bitesPerRabid * (censusVector['infective'])
    censusVector['PEPs'] <- PEPperNonRabidBite * bitesNonRabid + PEPperRabidBite * bitesRabid
    censusVector['lifeLoss'] <- lifeLossPerRabidBite * ((1 - PEPperRabidBite) * bitesRabid)
    
    return(censusVector)
  }
  ########################################
  
  
  ########################################
  ManagementFunction <- function(d, marginalCost, dailyBudget, totalSpending, totalContacted) {
    dailySpending <- 0
    if (dailyBudget[d] > 0) {
      count <- 0
      while (dailySpending < dailyBudget[d] & min(popMatrix[, 'contacted']) == 0) {
        
        # if there are uncontacted dogs left in the lowest marginal cost category, contact them first
        if (sum(popMatrix[, 'contacted'] == 0 & popMatrix[, 'contactCost'] == marginalCost[1]) > 0) {
          dogNumber <- sample(rep(which(popMatrix[, 'contacted'] == 0 & 
                                          popMatrix[, 'contactCost'] == marginalCost[1]), 2), 1)
          # now check for uncontacted in 2nd lowest marginal cost category
        } else if (sum(popMatrix[, 'contacted'] == 0 & popMatrix[, 'contactCost'] == marginalCost[2]) > 0) {
          dogNumber <- sample(rep(which(popMatrix[, 'contacted'] == 0 & 
                                          popMatrix[, 'contactCost'] == marginalCost[2]), 2), 1)
          # and for 2nd highest marginal cost category
        } else if (sum(popMatrix[, 'contacted'] == 0 & popMatrix[, 'contactCost'] == marginalCost[3]) > 0) {
          dogNumber <- sample(rep(which(popMatrix[, 'contacted'] == 0 & 
                                          popMatrix[, 'contactCost'] == marginalCost[3]), 2), 1)
          # and for highest marginal cost category
        } else if (sum(popMatrix[, 'contacted'] == 0 & popMatrix[, 'contactCost'] == marginalCost[4]) > 0) {
          dogNumber <- sample(rep(which(popMatrix[, 'contacted'] == 0 & 
                                          popMatrix[, 'contactCost'] == marginalCost[4]), 2), 1)
        } else {
          break
        }
        
        popMatrix[dogNumber, 'contacted'] <- 1
        totalContacted <- totalContacted + 1
        dailySpending <- dailySpending + as.numeric(popMatrix[dogNumber, 'contactCost'])
        
        if (popMatrix[dogNumber, 'female'] == 1) {
          # FEMALE management starts here
          
          if (popMatrix[dogNumber, 'puppy'] == 1) {
            # female PUPPY management here
            if (strategyVector['euthPuppyFemale'] == 1) {
              popMatrix <- popMatrix[!dogNumber, , drop=FALSE]
              dailySpending <- dailySpending + as.numeric(strategyCostVector['euthPuppyFemale'])
            } else {
              if (strategyVector['sterPuppyFemale'] == 1) {
                if (popMatrix[dogNumber, 'sterilized'] == 0) {
                  popMatrix[dogNumber, 'sterilized'] <- 1
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['sterPuppyFemale'])
                }
              } else if (strategyVector['contraPuppyFemale'] == 1) {
                if (popMatrix[dogNumber, 'sterilized'] == 0) {
                  # we won't contracept if dog has been sterilized, be we will even if already contracepted
                  popMatrix[dogNumber, 'contracepted'] <- 1
                  popMatrix[dogNumber, 'timeContra'] <- 0
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['contraPuppyFemale'])
                }
              }
              if (strategyVector['vaccPuppyFemale'] == 1) {
                if (popMatrix[dogNumber, 'vaccinated'] == 0) {
                  popMatrix[dogNumber, 'vaccinated'] <- 1
                  popMatrix[dogNumber, 'timeVacc'] <- 0
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['vaccPuppyFemale'])
                } else if (boosterGiven == TRUE) {
                  popMatrix[dogNumber, 'vaccinated'] <- 1
                  popMatrix[dogNumber, 'boosted'] <- 1
                  popMatrix[dogNumber, 'timeVacc'] <- 0
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['vaccPuppyFemale'])
                }
              }
            }
          } else if (popMatrix[dogNumber, 'adult'] == 1) {
            # female ADULT management here
            if (strategyVector['euthAdultFemale'] == 1) {
              popMatrix <- popMatrix[!dogNumber, , drop=FALSE]
              dailySpending <- dailySpending + as.numeric(strategyCostVector['euthAdultFemale'])
            } else {
              if (strategyVector['sterAdultFemale'] == 1) {
                if (popMatrix[dogNumber, 'sterilized'] == 0) {
                  popMatrix[dogNumber, 'sterilized'] <- 1
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['sterAdultFemale'])
                }
              } else if (strategyVector['contraAdultFemale'] == 1) {
                if (popMatrix[dogNumber, 'sterilized'] == 0) {
                  # we won't contracept if dog has been sterilized, be we will even if already contracepted
                  popMatrix[dogNumber, 'contracepted'] <- 1
                  popMatrix[dogNumber, 'timeContra'] <- 0
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['contraAdultFemale'])
                }
              }
              if (strategyVector['vaccAdultFemale'] == 1) {
                if (popMatrix[dogNumber, 'vaccinated'] == 0) {
                  popMatrix[dogNumber, 'vaccinated'] <- 1
                  popMatrix[dogNumber, 'timeVacc'] <- 0
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['vaccAdultFemale'])
                } else if (boosterGiven == TRUE) {
                  popMatrix[dogNumber, 'vaccinated'] <- 1
                  popMatrix[dogNumber, 'boosted'] <- 1
                  popMatrix[dogNumber, 'timeVacc'] <- 0
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['vaccAdultFemale'])
                }
              }
            }
          } else {
            # female JUVENILE management here
            if (strategyVector['euthJuvFemale'] == 1) {
              popMatrix <- popMatrix[!dogNumber, , drop=FALSE]
              dailySpending <- dailySpending + as.numeric(strategyCostVector['euthJuvFemale'])
            } else {
              if (strategyVector['sterJuvFemale'] == 1) {
                if (popMatrix[dogNumber, 'sterilized'] == 0) {
                  popMatrix[dogNumber, 'sterilized'] <- 1
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['sterJuvFemale'])
                }
              } else if (strategyVector['contraJuvFemale'] == 1) {
                if (popMatrix[dogNumber, 'sterilized'] == 0) {
                  # we won't contracept if dog has been sterilized, be we will even if already contracepted
                  popMatrix[dogNumber, 'contracepted'] <- 1
                  popMatrix[dogNumber, 'timeContra'] <- 0
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['contraJuvFemale'])
                }
              }
              if (strategyVector['vaccJuvFemale'] == 1) {
                if (popMatrix[dogNumber, 'vaccinated'] == 0) {
                  popMatrix[dogNumber, 'vaccinated'] <- 1
                  popMatrix[dogNumber, 'timeVacc'] <- 0
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['vaccJuvFemale'])
                } else if (boosterGiven == TRUE) {
                  popMatrix[dogNumber, 'vaccinated'] <- 1
                  popMatrix[dogNumber, 'boosted'] <- 1
                  popMatrix[dogNumber, 'timeVacc'] <- 0
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['vaccJuvFemale'])
                }
              }
            }
          } 
          
        } else {
          # MALE management starts here
          if (popMatrix[dogNumber, 'puppy'] == 1) {
            # male PUPPY management here
            if (strategyVector['euthPuppyMale'] == 1) {
              popMatrix <- popMatrix[!dogNumber, , drop=FALSE]
              dailySpending <- dailySpending + as.numeric(strategyCostVector['euthPuppyMale'])
            } else {
              if (strategyVector['sterPuppyMale'] == 1) {
                if (popMatrix[dogNumber, 'sterilized'] == 0) {
                  popMatrix[dogNumber, 'sterilized'] <- 1
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['sterPuppyMale'])
                }
              } else if (strategyVector['contraPuppyMale'] == 1) {
                if (popMatrix[dogNumber, 'sterilized'] == 0) {
                  # we won't contracept if dog has been sterilized, be we will even if already contracepted
                  popMatrix[dogNumber, 'contracepted'] <- 1
                  popMatrix[dogNumber, 'timeContra'] <- 0
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['contraPuppyMale'])
                }
              }
              if (strategyVector['vaccPuppyMale'] == 1) {
                if (popMatrix[dogNumber, 'vaccinated'] == 0) {
                  popMatrix[dogNumber, 'vaccinated'] <- 1
                  popMatrix[dogNumber, 'timeVacc'] <- 0
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['vaccPuppyMale'])
                } else if (boosterGiven == TRUE) {
                  popMatrix[dogNumber, 'vaccinated'] <- 1
                  popMatrix[dogNumber, 'boosted'] <- 1
                  popMatrix[dogNumber, 'timeVacc'] <- 0
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['vaccPuppyMale'])
                }
              }
            }
            
          } else if (popMatrix[dogNumber, 'adult'] == 1) {
            # male ADULT management here
            if (strategyVector['euthAdultMale'] == 1) {
              popMatrix <- popMatrix[!dogNumber, , drop=FALSE]
              dailySpending <- dailySpending + as.numeric(strategyCostVector['euthAdultMale'])
            } else {
              if (strategyVector['sterAdultMale'] == 1) {
                if (popMatrix[dogNumber, 'sterilized'] == 0) {
                  popMatrix[dogNumber, 'sterilized'] <- 1
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['sterAdultMale'])
                }
              } else if (strategyVector['contraAdultMale'] == 1) {
                if (popMatrix[dogNumber, 'sterilized'] == 0) {
                  # we won't contracept if dog has been sterilized, be we will even if already contracepted
                  popMatrix[dogNumber, 'contracepted'] <- 1
                  popMatrix[dogNumber, 'timeContra'] <- 0
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['contraAdultMale'])
                }
              }
              if (strategyVector['vaccAdultMale'] == 1) {
                if (popMatrix[dogNumber, 'vaccinated'] == 0) {
                  popMatrix[dogNumber, 'vaccinated'] <- 1
                  popMatrix[dogNumber, 'timeVacc'] <- 0
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['vaccAdultMale'])
                } else if (boosterGiven == TRUE) {
                  popMatrix[dogNumber, 'vaccinated'] <- 1
                  popMatrix[dogNumber, 'boosted'] <- 1
                  popMatrix[dogNumber, 'timeVacc'] <- 0
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['vaccAdultMale'])
                }
              }
            }
            
          } else {
            # male JUVENILE management here
            if (strategyVector['euthJuvMale'] == 1) {
              popMatrix <- popMatrix[!dogNumber, , drop=FALSE]
              dailySpending <- dailySpending + as.numeric(strategyCostVector['euthJuvMale'])
            } else {
              if (strategyVector['sterJuvMale'] == 1) {
                if (popMatrix[dogNumber, 'sterilized'] == 0) {
                  popMatrix[dogNumber, 'sterilized'] <- 1
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['sterJuvMale'])
                }
              } else if (strategyVector['contraJuvMale'] == 1) {
                if (popMatrix[dogNumber, 'sterilized'] == 0) {
                  # we won't contracept if dog has been sterilized, be we will even if already contracepted
                  popMatrix[dogNumber, 'contracepted'] <- 1
                  popMatrix[dogNumber, 'timeContra'] <- 0
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['contraJuvMale'])
                }
              }
              if (strategyVector['vaccJuvMale'] == 1) {
                if (popMatrix[dogNumber, 'vaccinated'] == 0) {
                  popMatrix[dogNumber, 'vaccinated'] <- 1
                  popMatrix[dogNumber, 'timeVacc'] <- 0
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['vaccJuvMale'])
                } else if (boosterGiven == TRUE) {
                  popMatrix[dogNumber, 'vaccinated'] <- 1
                  popMatrix[dogNumber, 'boosted'] <- 1
                  popMatrix[dogNumber, 'timeVacc'] <- 0
                  dailySpending <- dailySpending + as.numeric(strategyCostVector['vaccJuvMale'])
                }
              }
            }
          } 
        }
      }
    }
    
    return(list(popMatrix, totalContacted, dailySpending))
  }
  ########################################
  
  ########################################
  TimeFunction <- function() {
    # Arguments: None.
    # Return:    An updated population matrix.
    # Purpose:   Updates time-related columns in the population matrix.
    
    if (boosterGiven == FALSE) {
      # Turn off vaccinated and contracepted when time limit reached: 
      popMatrix[popMatrix[, 'timeVacc'] == timeVaccineEffective, 'vaccinated'] <- 0
      popMatrix[(popMatrix[, 'timeContra'] == timeContraEffectiveFemales & popMatrix[, 'female'] == 1), 
                'contracepted'] <- 0
      popMatrix[(popMatrix[, 'timeContra'] == timeContraEffectiveMales & popMatrix[, 'female'] == 0), 
                'contracepted'] <- 0
    } else {
      # Turn off vaccinated when time limit reached: 
      popMatrix[(popMatrix[, 'timeVacc'] == timeVaccineEffective & popMatrix[, 'boosted'] == 0), 'vaccinated'] <- 0
      
      # Turn off boosted and vaccinated when time limit reached:
      popMatrix[(popMatrix[, 'timeVacc'] == timeBoosterEffective & popMatrix[, 'boosted'] == 1), 'vaccinated'] <- 0
      popMatrix[(popMatrix[, 'timeVacc'] == timeBoosterEffective & popMatrix[, 'boosted'] == 1), 'boosted'] <- 0
      
      # Turn off contracepted when time limit reached:
      popMatrix[(popMatrix[, 'timeContra'] == timeContraEffectiveFemales & popMatrix[, 'female'] == 1), 
                'contracepted'] <- 0
      popMatrix[(popMatrix[, 'timeContra'] == timeContraEffectiveMales & popMatrix[, 'female'] == 0), 
                'contracepted'] <- 0
    }
    
    popMatrix[, 'age']           <- popMatrix[, 'age'] + 1
    popMatrix[, 'timeVacc']      <- popMatrix[, 'timeVacc'] + 1
    popMatrix[, 'timeContra']    <- popMatrix[, 'timeContra'] + 1
    popMatrix[, 'timeExposed']   <- popMatrix[, 'timeExposed'] + 1
    popMatrix[, 'timeInfective'] <- popMatrix[, 'timeInfective'] + 1
    
    return(popMatrix)
  }
  ########################################
  
  
  # inputs for initial population
  initialPopSize    <- paramCombos[c, 'initialPopSize']
  initialFracAdult  <- 0.61
  initialFracPup    <- 0.33
  initialFracFemale <- 0.38
  initialFracImmune <- paramCombos[c, 'initialFracImmune']
  initialFracContra <- 0.0   # UNUSED
  initialFracVacc   <- 0.0   # USER CHOICE
  initialFracSter   <- 0.0   # UNUSED
  
  # inputs for mortality
  maxJuvAge        <- 299   
  maxPuppyAge      <- 89
  maxAge           <- 4000
  carryingCap      <- initialPopSize * 1.2  # USER CHOICE
  pupAnnMortProb   <- 0.90  # Hampson et al. 2009
  juvAnnMortProb   <- 0.63  # Hampson et al. 2009
  adultAnnMortProb <- 0.32  # Hampson et al. 2009 
  emigrationProb   <- 0     # MISSING OR UNUSED
  
  # inputs for reproduction
  immigrantDogs        <- round(carryingCap * paramCombos[c, 'inMigrationFactor'])
  fecundityFactor      <- paramCombos[c, 'fecundityFactor']
  expLitterPer         <- 0.31 * fecundityFactor
  meanLitterSize       <- 4.4 * fecundityFactor            
  femalePupProb        <- 0.38
  fractionBirthPulse   <- 0.0              # MISSING or UNUSED         
  birthPulseVector     <- rep(0, 12)  
  birthPulseVector[1]  <- 0                # MISSING or UNUSED
  birthPulseVector[2]  <- 0                # MISSING or UNUSED
  birthPulseVector[3]  <- 0                # MISSING or UNUSED
  birthPulseVector[4]  <- 0                # MISSING or UNUSED  
  birthPulseVector[5]  <- 0                # MISSING or UNUSED
  birthPulseVector[6]  <- 0                # MISSING or UNUSED
  birthPulseVector[7]  <- 0                # MISSING or UNUSED
  birthPulseVector[8]  <- 0                # MISSING or UNUSED
  birthPulseVector[9]  <- 0                # MISSING or UNUSED
  birthPulseVector[10] <- 0                # MISSING or UNUSED  
  birthPulseVector[11] <- 0                # MISSING or UNUSED
  birthPulseVector[12] <- 0                # MISSING or UNUSED
  
  # inputs for disease
  monthInitIntroduction <- 1     # USER CHOICE
  monthsOfPressure      <- 2     # USER CHOICE
  dogsPerMonthExposed   <- ceiling(paramCombos[c,'epidemicSize']*initialPopSize/12)  # USER CHOICE
  timeLimitExposed      <- 22
  timeLimitInfective    <- 3
  survivalProb          <- paramCombos[c, 'survivalProb']
  meanFactor            <- paramCombos[c, 'meanFactor']
  bitesPerRabidMean     <- 2.15 * meanFactor
  varianceFactor        <- paramCombos[c, 'varianceFactor']
  bitesPerRabidShape     <- 1.33 * varianceFactor
  probInfectionFromBite <- 0.49
  
  exposedTimeShape      <- 1.08549138
  exposedTimeRate       <- 0.04919551
  infectiveTimeShape    <- 2.831788
  infectiveTimeRate     <- 0.9193612
  
  
  # inputs for benefits of management
  bitesPerNonRabid     <- 0  # UNUSED
  bitesPerRabid        <- 0  # UNUSED   
  PEPperNonRabidBite   <- 0  # UNUSED
  PEPperRabidBite      <- 0  # UNUSED
  costPerPEP           <- 0  # UNUSED
  lifeLossPerRabidBite <- 0  # UNUSED
  
  # inputs for treatment costs
  vaccineCost             <- 0  # UNUSED
  contraceptionCostFemale <- 0  # UNUSED
  contraceptionCostMale   <- 0  # UNUSED
  sterilizationCostFemale <- 0  # UNUSED
  sterilizationCostMale   <- 0  # UNUSED
  euthanasiaCost          <- 0  # UNUSED
  
  # inputs for effectiveness of contraception and vaccination
  timeVaccineEffective       <- 0  # UNUSED
  timeBoosterEffective       <- 0  # UNUSED
  timeContraEffectiveMales   <- 0  # UNUSED
  timeContraEffectiveFemales <- 0  # UNUSED
  
  # inputs for contact costs
  # note: 25, 50, 75, 100, mean 25% 50%, 75%, 100% of specified initial abundance
  contactCost25  <- 0  # UNUSED
  contactCost50  <- 0  # UNUSED
  contactCost75  <- 0  # UNUSED
  contactCost100 <- 0  # UNUSED
  
  # input for budget years 1-5    
  annualBudget     <- rep(0, 10)
  annualBudget[1]  <- 0  # UNUSED
  annualBudget[2]  <- 0  # UNUSED
  annualBudget[3]  <- 0  # UNUSED
  annualBudget[4]  <- 0  # UNUSED
  annualBudget[5]  <- 0  # UNUSED
  annualBudget[6]  <- 0  # UNUSED
  annualBudget[7]  <- 0  # UNUSED
  annualBudget[8]  <- 0  # UNUSED
  annualBudget[9]  <- 0  # UNUSED
  annualBudget[10] <- 0  # UNUSED
  
  # inputs for strategy
  # note: model assumes already sterilized dogs are not re-sterilized. Within the same year dogs will not be vaccinated 
  #       or contracepted twice. If dogs are re-contancted in a future year, they will be re-vaccinated or re-contracepted
  # note: contraception and sterilization cannot both equal 1 for same demographic
  # note: if euthanisia equal 1 for some demographic, all other treatments must equal zero  
  vaccPuppyMale     <- 0  # UNUSED
  vaccPuppyFemale   <- 0  # UNUSED
  vaccAdultMale     <- 0  # UNUSED
  vaccAdultFemale   <- 0  # UNUSED
  vaccJuvMale       <- 0  # UNUSED
  vaccJuvFemale     <- 0  # UNUSED
  contraPuppyMale   <- 0  # UNUSED
  contraPuppyFemale <- 0  # UNUSED
  contraAdultMale   <- 0  # UNUSED
  contraAdultFemale <- 0  # UNUSED
  contraJuvMale     <- 0  # UNUSED
  contraJuvFemale   <- 0  # UNUSED
  sterPuppyMale     <- 0  # UNUSED
  sterPuppyFemale   <- 0  # UNUSED
  sterAdultMale     <- 0  # UNUSED
  sterAdultFemale   <- 0  # UNUSED
  sterJuvMale       <- 0  # UNUSED
  sterJuvFemale     <- 0  # UNUSED
  euthPuppyMale     <- 0  # UNUSED
  euthPuppyFemale   <- 0  # UNUSED
  euthAdultMale     <- 0  # UNUSED
  euthAdultFemale   <- 0  # UNUSED
  euthJuvMale       <- 0  # UNUSED
  euthJuvFemale     <- 0  # UNUSED
  boosterGiven      <- TRUE  # UNUSED
  
  # inputs for management timing
  mgtMonthVector     <- rep(0, 12)
  mgtMonthVector[1]  <- 0  # UNUSED
  mgtMonthVector[2]  <- 1  # UNUSED
  mgtMonthVector[3]  <- 0  # UNUSED
  mgtMonthVector[4]  <- 0  # UNUSED
  mgtMonthVector[5]  <- 0  # UNUSED
  mgtMonthVector[6]  <- 0  # UNUSED
  mgtMonthVector[7]  <- 0  # UNUSED
  mgtMonthVector[8]  <- 0  # UNUSED
  mgtMonthVector[9]  <- 0  # UNUSED
  mgtMonthVector[10] <- 0  # UNUSED
  mgtMonthVector[11] <- 0  # UNUSED
  mgtMonthVector[12] <- 0  # UNUSED
  ########################################
  
  
  ########################################
  # Misc preliminary calculations and assignments:
  
  # Get total number of days in simulation:
  simulationEnd   <- 365 * simulationYears
  
  # A vector of month number for use in seasonal timing:
  monthSeries <- c(rep(1, 31), rep(2, 28), rep(3, 31), rep(4, 30), 
                   rep(5, 31), rep(6, 30), rep(7, 31), rep(8, 31),
                   rep(9, 30), rep(10, 31), rep(11, 30), rep(12, 31)) 
  
  monthFirstDays <- rep(c(match(1, monthSeries), match(2, monthSeries), match(3, monthSeries), match(4, monthSeries),
                          match(5, monthSeries), match(6, monthSeries), match(7, monthSeries), match(8, monthSeries),
                          match(9, monthSeries), match(10, monthSeries), match(11, monthSeries), match(12, monthSeries)),
                        simulationYears)
  
  # Get days of each year that disease will be introduced:
  pressureMonths <- seq(monthInitIntroduction, monthInitIntroduction + monthsOfPressure - 1)  
  pressureYears <- (pressureMonths - 1) %/% 12 + 1
  pressureDays <- list()
  for (i in 1:simulationYears) {
    if (sum(pressureYears == i) > 0) {
      pressureDays[[i]] <- monthFirstDays[pressureMonths[pressureYears == i]]
    } else {
      pressureDays[[i]] <- 0
    }
  }
  flush.console()
  
  # Calculate demographics of initial population:
  initialAdults     <- round(initialFracAdult * initialPopSize)
  initialSubAdults  <- initialPopSize - initialAdults
  initialPuppies    <- round(initialFracPup * initialSubAdults)
  initialJuveniles  <- initialSubAdults - initialPuppies
  
  # Calculate daily mortality probabilities:
  pupMortalityProb   <- 1 - (1 - pupAnnMortProb) ^ (1/365)
  juvMortalityProb   <- 1 - (1 - juvAnnMortProb) ^ (1/365)
  adultMortalityProb <- 1 - (1 - adultAnnMortProb) ^ (1/365)
  
  # Calculate daily litter probabilities:
  monthDayCount <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  if(sum(birthPulseVector) != 12 & sum(birthPulseVector != 0)) {
    peakDays <- sum(birthPulseVector * monthDayCount)
    peakProb <- (fractionBirthPulse * expLitterPer) / peakDays
    offPeakDays <- sum((!birthPulseVector) * monthDayCount)
    offPeakProb <- ((1-fractionBirthPulse) * expLitterPer) / offPeakDays
    litterProbability <- rep(offPeakProb, 365)
    for(m in 1:12) {
      if(birthPulseVector[m] == 1) {
        litterProbability[monthSeries == m] <- peakProb
      }
    }
  } else {
    litterProbability <- rep(expLitterPer / 365, 365)
  }
  
  # Calculate marginal costs of contact:
  marginalCost1 <- contactCost25 / (initialPopSize * 0.25)
  marginalCost2 <- contactCost50 / (initialPopSize * 0.25)
  marginalCost3 <- contactCost75 / (initialPopSize * 0.25)
  marginalCost4 <- contactCost100 / (initialPopSize * 0.25)
  marginalCost <- c(marginalCost1, marginalCost2, marginalCost3, marginalCost4)
  
  # List days that management will occur:
  mgtDayVector <- rep(0, 365)
  for(m in 1:12) {
    if(mgtMonthVector[m] == 1) {
      mgtDayVector[monthSeries == m] <- 1
    }
  }
  managementDays <- seq(1, 365)
  managementDays <- managementDays[mgtDayVector == 1]
  
  # A list of traits in the population matrix:
  traitList <- c('age', 'puppy', 'adult','female',
                 'sterilized', 'contracepted', 'timeContra',
                 'vaccinated', 'timeVacc',
                 'boosted', 'contacted', 'contactCost',
                 'exposed', 'timeExposed', 'timeLimitExposed',
                 'infective', 'timeInfective', 'timeLimitInfective',
                 'immune', 'month')
  
  # A list of results that will be tracked:
  censusSeries <- c('abundance', 'puppy', 'adult', 'females', 
                    'sterilized', 'femalesSterilized',
                    'contracepted', 'femalesContracepted', 
                    'vaccinated', 'immune', 'exposed', 'infective',
                    'PEPs', 'lifeLoss', 'newlyVaccinated') 
  censusVector <- rep(0, length(censusSeries))
  names(censusVector) <- censusSeries
  
  # Create a 3d array to store results:
  resultsMatrix <- array(data=NA, dim=c(simulationEnd, length(censusSeries), iterations))
  colnames(resultsMatrix) <- censusSeries
  
  # Create a vector of binary strategy indicators:
  strategyNames <- c('vaccPuppyMale', 'vaccPuppyFemale',
                     'vaccAdultMale', 'vaccAdultFemale', 
                     'vaccJuvMale', 'vaccJuvFemale',
                     'contraPuppyMale', 'contraPuppyFemale',
                     'contraAdultMale', 'contraAdultFemale',
                     'contraJuvMale', 'contraJuvFemale',
                     'sterPuppyMale', 'sterPuppyFemale',
                     'sterAdultMale', 'sterAdultFemale', 
                     'sterJuvMale', 'sterJuvFemale',
                     'euthPuppyMale', 'euthPuppyFemale',
                     'euthAdultMale', 'euthAdultFemale', 
                     'euthJuvMale', 'euthJuvFemale') 
  strategyVector <- c(vaccPuppyMale, vaccPuppyFemale,
                      vaccAdultMale, vaccAdultFemale, 
                      vaccJuvMale, vaccJuvFemale,
                      contraPuppyMale, contraPuppyFemale,
                      contraAdultMale, contraAdultFemale, 
                      contraJuvMale, contraJuvFemale,
                      sterPuppyMale, sterPuppyFemale,
                      sterAdultMale, sterAdultFemale, 
                      sterJuvMale, sterJuvFemale,
                      euthPuppyMale, euthPuppyFemale,
                      euthAdultMale, euthAdultFemale, 
                      euthJuvMale, euthJuvFemale)
  names(strategyVector) <- strategyNames
  
  # Create a cost vector to indicate unit cost of each strategy:
  strategyCostVector <- c(rep(vaccineCost, 6),
                          contraceptionCostMale, contraceptionCostFemale,
                          contraceptionCostMale, contraceptionCostFemale, 
                          contraceptionCostMale, contraceptionCostFemale,
                          sterilizationCostMale, sterilizationCostFemale,
                          sterilizationCostMale, sterilizationCostFemale, 
                          sterilizationCostMale, sterilizationCostFemale,
                          rep(euthanasiaCost, 6))
  names(strategyCostVector) <- strategyNames
  ########################################
  
  
  ########################################
  # Loop through iterations:
  for(i in 1:iterations) {
    print(paste('Running iteration', i))
    flush.console()
    popMatrix <- InitialPopulation()
    
    # Loop through years:
    for(j in 1:simulationYears) {
      # reset total spending, number of dogs contacted, and contacted indicator at start of year
      totalSpending <- 0
      totalContacted <- 0
      popMatrix[, 'contacted'] <- 0
      # get the daily budget for each day of year
      dailyBudget <- getDailyBudget(j)
      
      # Loop through days of the year
      for(d in 1:365) {
        popMatrix[, 'month'] <- monthSeries[d]
        resultsMatrix[(365 * (j-1) + d), ,i] <- CensusFunction()
        popMatrix <- MortalityFunction()
        popMatrix <- ReproductionFunction(d)
        popMatrix <- ImmigrationFunction()
        popMatrix <- DiseaseProgressionFunction()
        popMatrix <- DiseaseSpreadFunction()
        tempVacc <- sum(popMatrix[, 'vaccinated'])
        if (totalSpending < annualBudget[j]) {
          mgtReturnList <- ManagementFunction(d, marginalCost, dailyBudget, totalSpending, totalContacted)
          popMatrix <- mgtReturnList[[1]]
          totalContacted <- mgtReturnList[[2]]
          totalSpending <- totalSpending + mgtReturnList[[3]]
        }
        # Record new vaccinations:
        resultsMatrix[(365 * (j-1) + d), 'newlyVaccinated', i] <- sum(popMatrix[, 'vaccinated']) - tempVacc
        popMatrix <- TimeFunction()
      }  # close d for loop
    }  # close j for loop
  }  # close i for loop
  
  # Record % of iterations with no rabies at end of each year
  for (i in 1:simulationYears) {
    qElimination[i] <- sum(resultsMatrix[(i*365), 'infective', ] + resultsMatrix[(i*365), 'exposed', ] == 0)
    qNoDogsLeft[i] <- sum(resultsMatrix[(i*365), 'abundance', ] == 0)
    qElimWithDogsLeft[i] <- sum(resultsMatrix[(i*365), 'abundance', ] > 0 &
                                  (resultsMatrix[(i*365), 'infective', ] + 
                                     resultsMatrix[(i*365), 'exposed', ] == 0))
    qNoElimWithDogsLeft[i] <- sum(resultsMatrix[(i*365), 'abundance', ] > 0 &
                                    (resultsMatrix[(i*365), 'infective', ] + 
                                       resultsMatrix[(i*365), 'exposed', ] > 0))
  }
  resultRow <- c(qElimination, qNoDogsLeft, qElimWithDogsLeft, qNoElimWithDogsLeft)
  ########################################
  
}  # end for loop thru parameter combinations

# print time
print(proc.time() - startTime)

colnames(simulationData) <- c('qElim1',
                              'qExtinct1',
                              'qElimWithDogs1',
                              'qNoElimWithDogs1')

paramCombos <- paramCombos[, !names(paramCombos) %in% c("varianceFactor")] #remove column that is not relevant to further analysis

colnames(paramCombos) <- c('Population',
                          'Fecundity',
                          'Infectivity',
                          'Exposure_survival',
                          'Initial_immunity',
                          'In_migration',
                          'Incidence') #rename columns for regression for better understanding in regression analysis
write.csv(simulationData, 'simulationData_pop5400.csv')
write.csv(paramCombos, 'paramCombos_pop5400.csv')

