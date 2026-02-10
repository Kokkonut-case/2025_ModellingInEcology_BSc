## Set your wd
setwd("/path/to/wd") # 

#############################################################
##################### LMC model #############################
#############################################################
###### Inspired by Hamilton´s unbeatable sex ratio (1967): (n-1)/2n ######
# n = number of egglaying females in the group
# calculate optimal sex ratio for different numbers of females
n_females <- 1:25
optimal_ratio <- (n_females-1)/(2*n_females)
# create the figure
plot_ratio <- plot(n_females, optimal_ratio, type = "l", ylim = c(0,0.5), xlim = c(0,25),
                   main = "LMC by Hamilton", xlab = "number of egg-laying females (n)", 
                   ylab = "optimal sex ratio", col = "black")

# MB: good commenting, chuck of code works, produces figure with a title and appropriate axis names 

#############################################################
##################### (1) FEMINIZATION ######################
#############################################################
T_FE <- 0.95 # transmission efficiency
beta_FE <- 0.9 # probability that infected males are feminized
time_steps <- 50 # number of generations
# dataframe to store frequencies
fe <- data.frame(
time = 1:time_steps,
f_FE = rep(NA, time_steps), # infected females
m_FE = rep(NA, time_steps), # infected males
f_U = rep(NA, time_steps), # uninfected females
m_U = rep(NA, time_steps) # uninfected males
)
# initial frequencies
eps <- 0.005 # small initial infection
fe[1, c("f_FE", "m_FE", "f_U", "m_U")] <- c( 
eps, # infected females (start very rare)
0, # infected males (almost none)
(1-eps)/2, # uninfected females (almost 0.5)
(1-eps)/2 # uninfected males (almost 0.5)
)
# MB: runs perfectly fine! As a side note, this more concise syntax also works: fe[1,2:5] <- c(eps,0,(1-eps)/2,(1-eps)/2)

#####
# Replicator equation of FE dynamics (R is not important) 
##### 
for(t in 1:(time_steps - 1)){
    # mean fitness
    omega_bar <- fe$f_FE[t] * T_FE * (1 + beta_FE) + fe$f_U[t]
    # replicator equations
    fe$f_FE[t+1] <- (fe$f_FE[t] * T_FE * (1 + beta_FE)) / omega_bar
    fe$m_FE[t+1] <- 0
    fe$f_U[t+1] <- fe$f_U[t] / omega_bar
    fe$m_U[t+1] <- fe$f_U[t] / omega_bar
}
#####
# Plot #
####
x11(); plot(fe$time, fe$f_FE, type = "l", ylim = c(0,1), main = "Feminization, Replicator Equation Dynamics", xlab = "Time (Generations)", ylab = "Frequency", lwd = 2)
lines(fe$time, fe$m_FE, lwd = 2, col = "red")
lines(fe$time, fe$f_U, lwd = 2, col = "darkgreen", lty=3)
lines(fe$time, fe$m_U, lwd = 2, col = "pink",lty=2)

# MB: runs without issue / produces a plot with axis labels and title BUT no color legend (I don't know which curve corresponds to what)
# missing green line --> because it overlaps with pink because recursions define f_U and m_U similarly --> mistake?

#############################################################
##################### (2) MALE KILLING ######################
#############################################################
T_MK <- 0.95 # transmission efficieny
R <- 0.3 # resource advantage for sisters due to male killing
time_steps <- 100 # number of generation
# dataframe to store frequencies
mk <- data.frame(
time = 1:time_steps,
f_MK = rep(NA, time_steps), # infected females
m_MK = rep(NA, time_steps), # infected males (killed)
f_U = rep(NA, time_steps), # uninfected females
m_U = rep(NA, time_steps) # uninfected males
)
# initial frequencies
eps <- 0.005 # small initial infection
mk[1, c("f_MK", "m_MK", "f_U", "m_U")] <- c(
eps, # rare infected females
0, # infected males are killed
(1-eps)/2, # uninfected females (almost 0.5)
(1-eps)/2 # uninfected males (almost 0.5)
)
#####
# Replicator equation of MK dynamics
#####
    for(t in 1:(time_steps - 1)){
    # mean fitness
    omega_bar <- mk$f_MK[t] * T_MK * (1 + R) + mk$f_U[t]
    # replicator equations
    mk$f_MK[t+1] <- (mk$f_MK[t] * T_MK * (1 + R)) / omega_bar
    mk$m_MK[t+1] <- 0
    mk$f_U[t+1] <- mk$f_U[t] / omega_bar
    mk$m_U[t+1] <- mk$f_U[t] / omega_bar
}

#####
# Plot #
####
x11(); plot(mk$time, mk$f_MK, type = "l", ylim = c(0,1), main = "Male killing, Replicator Equation Dynamics", 
xlab = "Time (Generations)", ylab = "Frequency", lwd = 2)
lines(mk$time, mk$f_U, lwd = 2, col = "darkgreen", lty=2)
lines(mk$time, mk$m_U, lwd = 2, col = "pink", lty=3)

# MB: same as above, runs well but is missing a color caption. Larger number of time steps = problem to compare FE vs. MK.
# Missing the green line --> overlapping with pink because recursions define f_u and m_u to be the same... bug?
# if you replace R = 0.3 by R = 0.9 and time_steps = 50, you obtain EXACTLY the same curve as with the feminization...