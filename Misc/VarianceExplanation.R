

# Bites across infectious period (per rabid dog) are a random variable following the binomial distribution.
# The following two parameters determine the shape and location of the distribution.
bitesPerRabidMean     <- 2.15  # parameter from Hampson paper (mean bites across infectious period)
bitesPerRabidShape    <- 1.33  # parameter from Hampson paper

# if we draw a million instances of this random value and get mean and variance:
mean(rnbinom(1000000, size=bitesPerRabidShape, mu=bitesPerRabidMean))  # equals ~2.15
var(rnbinom(1000000, size=bitesPerRabidShape, mu=bitesPerRabidMean))   # equals ~5.62


# In our model, we need daily bites, not bites across whole infectious period.
timeLimitInfective <- 3
# In a single day:
mean(rnbinom(1000000, size=bitesPerRabidShape/timeLimitInfective, mu=bitesPerRabidMean/timeLimitInfective))
var(rnbinom(1000000, size=bitesPerRabidShape/timeLimitInfective, mu=bitesPerRabidMean/timeLimitInfective))
# Now randomly draw daily bites for all days of the infectious period
df <- matrix(NA, nrow=1000000, ncol=timeLimitInfective) # each column gives the bites on a day, rows are the draws
for (i in 1:timeLimitInfective) {
  df[, i] <- rnbinom(1000000, size=bitesPerRabidShape/timeLimitInfective, mu=bitesPerRabidMean/timeLimitInfective)
}
# Now sum bites across all days of infectious period:
totalBites <- apply(df, 1, sum)
mean(totalBites) # equals ~2.15
var(totalBites) # equals ~5.62
# So our parameter conversion for daily values was correct.


# Ok, so what about variance?
# Variance depends on both 'mu' and 'size'. If you adjust the arguments below, you can see this.
var(rnbinom(1000000, size=1.33, mu=2.15))   # equals ~5.62

# Specifically, variance = mu + (mu^2)/size
# So in our case
2.15 + (2.15^2) / 1.33   # eqauls ~5.62

# If you want to adjust variance and keep the mean constant, use this formula:
# size = -(mu^2) / (mu-var)
# for Hampson's values:
-(2.15^2)/(2.15 - 5.62)  # size equals ~1.33
# double the variance:
-(2.15^2)/(2.15 - (5.62*2))  # size equals ~0.508
# half the variance:
-(2.15^2)/(2.15 - (5.62*0.5))  # size equals ~7.003
# So just plug in whatever variance you want in the abvove formulas and get the corresponding 'size' parameter.
# This new 'size' is what you use in rnbinom; keep mu the same.
var(rnbinom(1000000, size=-(2.15^2)/(2.15 - (5.62*0.5)), mu=2.15))  # var = ~2.81 (half of 5.62)




