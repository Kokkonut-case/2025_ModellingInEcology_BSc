##Basic Overview
#Installing Stuff | Run the code beginning from line 9 (if its going to be a rerun) otherwise it will produce an annoying error 
install.packages("ggplot2")
install.packages("popbio")

library(ggplot2)
library(popbio)

#General Values
age <- c(100, rep(0,3)) * 1.0 #vector so it might be more usable; guarantee its a numeric (integer destroyed my whole code)
max_ageclasses <- 4 #to limit and it being the value that defines x on the graph
years <- 50 #max years to also maybe limit the loop
growthRate <- 1.01537

#Survival Rates
SJuv <- 0.607 #survival of year 0-1
SYea1 <- 0.38211 #survival of year 1-2; upper Value
SYea2 <- 0.15989 #survival of year 1-2; lower value
SPRe1 <- 0.419244 #survival of years >= 2 (pre-reproductive); upper value
SPRe2 <- 0.231756 #survival of years >= 2 (reproductive); lower value
SRep <- 0.747 #survival of years >= 2 (reproductive)

#Fecundity
FYea <- 0.297395 #fecundity of year 1-2
FPre <- 0.498275 #fecundity of years >= 2 (pre-reproductive)
FRep <- 0.941145 #fecundity of years >= 2 (reproductive)

#Matrix
MatrixSF <- matrix(0.0, nrow = max_ageclasses, ncol = max_ageclasses) #4x4 matrix so I can just use max_ageclasses as value to limit 
MatrixSF[2,1] <- SJuv #add the different survival rates based on the position in the matrix (first value = row, second = column)
MatrixSF[3,2] <- SYea1
MatrixSF[4,2] <- SYea2
MatrixSF[3,3] <- SPRe1
MatrixSF[4,3] <- SPRe2
MatrixSF[4,4] <- SRep

MatrixSF[1,2] <- FYea #Fecundity rates
MatrixSF[1,3] <- FPre
MatrixSF[1,4] <- FRep


#Stable stage distribution
stable_dist <- popbio::stable.stage(MatrixSF)
startpop <- 100
age <- stable_dist * startpop

## Deterministic Model
AoT_det <- matrix(0.0, nrow = max_ageclasses, ncol = years + 1) #Age over time matrix; had to add 0.0 because otherwise it would all be 0
#rows = age classes & columns = time steps
AoT_det[,1] <- age #first value is empty because its not necessary

#Modelling
for (t in 1:years) { #t = time in every year (we will just let it loop)
  AoT_det[,t+1] <- MatrixSF %*% AoT_det[,t] #Age over time with time step +1 while using the survival matrix from above and matrix multiplication %*%
}

total_det <- colSums(AoT_det)

## Stochastic Model
mat_sims <- 1000
Pop_Matrix <- matrix(NA, nrow = years + 1, ncol = mat_sims) #new matrix to replicate the simulation 1k times and store its values
GrowthRateMatrix <- matrix(NA, nrow = years, ncol = mat_sims)

for (i in 1:mat_sims) { #outer loop with simulations
  AoT_stoch <- matrix(0.0, nrow = max_ageclasses, ncol = years + 1) #Age over time matrix; had to add 0.0 because otherwise it would all be 0
  AoT_stoch[,1] <- age #first value is empty because its not necessary
  
  for (t in 1:years) {   #inner loop for the time
    odds <- SJuv / (1 - SJuv)              # baseline odds of juvenile survival
    odds_var <- odds * rlnorm(1, 0, 0.3)   # stochasticity on odds
    P <- odds_var / (1 + odds_var)         # back to probability
    P <- P /SJuv
    
    Matrix_var <- MatrixSF
    Matrix_var[2,1] <- SJuv* P #[2,1] -> changing only the surv rate of the juveniles
    Matrix_var[2,1] <- min(Matrix_var[2,1], 1)
    
    AoT_stoch[, t + 1] <- Matrix_var %*% AoT_stoch[, t]
  }
  total_pop <- colSums(AoT_stoch)
  Pop_Matrix[, i] <- total_pop
  for (t in 1:years) { #adding growth rate
    GrowthRateMatrix[t, i] <- total_pop[t + 1] / total_pop[t]
  }
}

MeanPop <- rowMeans(Pop_Matrix, na.rm = TRUE) #getting the means for the rows | with years + 1
MeanGrowth <- rowMeans(GrowthRateMatrix, na.rm = TRUE) #based on years
total_stoch <- MeanPop

#Plotting
pop_dataframe <- data.frame(year = 0:years, Deterministic = total_det, Stochastic = total_stoch, growthRate = c(NA, MeanGrowth)) #create a new dataframe to use in the plot (with stochastic, deterministic and growth rate)

ggplot(pop_dataframe, aes(x = year)) +
  geom_line(aes(y = Deterministic, color = "Deterministic"), size = 1.2) +
  geom_line(aes(y = Stochastic, color = "Stochastic"), size = 1.2) +
  scale_color_manual(values = c("Deterministic" = "black", "Stochastic" = "purple")) + #custom colors
  labs(title = "Population over time with 1000 simulations", x = "Year", y = "Total Population Size", color = "Legend") +
  theme_minimal()

lambda_det <- rep(growthRate, years) #lambda for growthrate so I can use in the plot below
growth_dataframe <- data.frame(year = 1:years, Deterministic = lambda_det, Stochastic = MeanGrowth) #1:years because growth rate is between years

ggplot(growth_dataframe, aes(x = year)) +
  geom_line(aes(y = Deterministic, color = "Deterministic"), size = 1.2) +
  geom_line(aes(y = Stochastic, color = "Stochastic"), size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("Deterministic" = "black", "Stochastic" = "purple")) + #custom colors again
  labs(title = "Growth Rate Trajectory", x = "Year", y = expression(lambda), color = "Model") +
  theme_minimal()