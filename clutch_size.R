getwd("path/to/wd") #our working directory was saved under our names, you can save it under your prefered folder and name.

library(ggplot2) #as we are ploting with the help of ggplot2, it needs to be loaded to the library

#we use R and ggplot to plot a function for the relationship between clutch size and adult mortality
adult_mortality <- function (c,k,m) 1/(1+exp(-k*(c-m))) #logistic function of adult mortality in dependance of clutch size
adult_mortality(2,3,5) #for c=2, k=3 and m=5, we have a mortality rate of 0.0001 for example
ggplot()+#graph the mortality rate function with the help of ggplot
  xlim(0,10)+
  geom_function(fun=adult_mortality, args=list(k=3, m=5))

#we produced a series of equations to calculate the fitness of both alleles,the relative fitness of the alleles and the proportion of the mutant allel at the next time step, following the steps of section 3.3.1 from Otto and Day (2007). each of this equations require the output of the previous equation for its calculation, which we can show in R calling a function within another function.
#using this knowledge we can code a function for proportion of the mutant allele at the next time step:p(t+1)=(p(t)*V_E)/(p(t)*V_E)+ (1-p(t))
#the relative fitness of allele E is calculated V_E =W_E/W_e
#keep in mind: clutchsize_adultmortality <- function (c,k,m) 1/(1+exp(-k*(c-m)))

#mj is the mortality rate for juveniles
#c is the resident and common clutch size (for allele e)
# b is the benefit parameter (how much bigger c is for allele E)
W_e <- function(c,k,m,mj)((1-adult_mortality(c,k,m))+c(1-mj)) #fitness W of allele e (resident strategy)
W_e(2,3,5,0.4) #for c=2, k=3, m=5 and mj=0.4 allele e has a fitness of 1.6, for example
W_E <- function(c,k,m,mj,b)((1-adult_mortality(c*b,k,m))+((c*b)*(1-mj)))#fitness W of allele E (invading strategy)
W_E(1,3,5,0.4,2) #for c=1, k=3, m=5, mj=0.4 and b=2 allele E has a fitness of 2.2, for example

V_E <- function(c,k,m,mj,b) W_E(c,k,m,mj,b)/W_e(c,k,m,mj) #relative fitness of allele E in relation to allele e.
V_E(1,3,5,0.4,2) #for c=1, k=3, m=5, mj=0.4 and b=2 allele E has a relative fitness of 1.4, for example

p_next_timestep <- function(c,k,m,mj,b,p)(p*V_E(c,k,m,mj,b))/((p*V_E(c,k,m,mj,b))+(1-p))
p_next_timestep (1,3,5,0.4,2,0.2) #for p=0.2, allele E takes in p+1 26% of the population
p_next_timestep (1,3,5,0.4,2,0.5) #for p=0.5, allele E takes in p+1 58% of the population
p_next_timestep (1,3,5,0.4,2,0.8) #for p=0.8, allele E takes in p+1 85% of the population

#with the help of ggplot we can now visualize how the future freq. of allele E depends on its current frequency. analyse this for with different values of mj.

straight_line <- function (l,p) l*p # a diagonal that shows no change between p on the last and next time step helps us analyse whether the frequency is increasing depending on whether the curve for set mj is below or above this diagonal.

ggplot() + #with the help of ggplot, we can visualize our function of p in the next time step. we chose as set parameters k=2 and m=1.
  xlim(0, 1)+
  ylim(0, 1)+
  geom_function(aes(colour="mj = 0.2"), fun = p_next_timestep,args = list(c = 1, b = 2, mj = 0.2, k = 2, m = 1),  linewidth=1.25)+
  geom_function(aes(colour="mj = 0.5"), fun = p_next_timestep,args = list(c = 1, b = 2, mj = 0.5, k = 2, m = 1),  linewidth=1.25)+
  geom_function(aes(colour="mj = 0.7"), fun = p_next_timestep,args = list(c = 1, b = 2, mj = 0.7, k = 2, m = 1),  linewidth=1.25)+
  geom_function(aes(colour="mj = 0.8"), fun = p_next_timestep,args = list(c = 1, b = 2, mj = 0.8, k = 2, m = 1),  linewidth=1.25)+
  geom_function(fun = straight_line, args=list(l=1), colour = "black",linetype=2, linewidth=1.25)+
  xlab('Current frequency of E')+
  ylab('Future frequency of E')+
  guides(colour=guide_legend(title= "Juvenile Mortality"))+
  theme_bw(base_size = 15)+
scale_colour_manual(values = c("#007191", "#62c8d3", "#f47a00", "#d31f11")) #Fig 6.

#we now create a loop that allows us to account not only for t+1, but for the next 
#10000 generations of the population

k<-2
m<-1
mj<-0.5
b<-2
p<-0.01
c<- 1

count <- 0 #from time step 0
while(count < 10){ #until the 10000th time step (meaning the 10000th generation)
  p<-p_next_timestep(c,k,m,mj,b,p) #p beggins at 0.01, and everytime the next step is accounted for, we start from p on the previous timestep --> p(t+1) = (p(t)*V_E(c,k,m,mj,b))/((p(t)*V_E(c,k,m,mj,b))+(1-p(t)))
  count <- count + 1 # keep accounting for the next p of allele E untill we get to 10000 generations
}

#we now create a data set with the help of expand.grid. this allows us to set different values to multiple rows that account for our variables and parameters

#our first data set (survivors) will have m=1 and k=3 as set parameters

survivors <- expand.grid(
  b = seq(1, 10, by = 0.5),
  mj = seq(0, 1, by = 0.05),
  c = 1,
  m = 1,
  k = 3,
  p = NA,
  ma = adult_mortality(c,k,m))

#we now run the loop so that, for each combination of the possible values, we obtain p at the next time step. the data set is automatically filled with the calculated p values. knowing the adult mortality values will help us analyse our heatmaps
#keep in mind: adult_mortality <- function (c,k,m) 1/(1+exp(-k*(c-m)))


for (i in 1:nrow(survivors)){ #for each variable in the number of rows of our data set
  c = survivors[i,3]
  m = survivors[i,4]
  k = survivors[i,5]
  b = survivors[i,1]
  mj = survivors[i,2]
  p = 0.01 #with p(0) being 0.01
  count <- 0 #beginning at time step 0
  while(count < 10000){ #until the 10000th time step (meaning the 10000th generation)
    p<-p_next_timestep(c,k,m,mj,b,p) 
    count <- count + 1 
  }
  survivors[i,6] <- p #so that p is automatically filled in
}

#we now plot a heat map that shows us if allele E invades, in relation to p on data set 1

ggplot(data = survivors, aes(mj,b,fill=p))+
  geom_tile(colour="black")+
  scale_fill_gradient(low="white", high="#feb326")+
  guides(fill=guide_colourbar(title= "p"))
         
#we repeat this process for our other data sets

#our second data set (HRS) will have m=-0.2 and k=1 as set parameters
HRS <- expand.grid(
  b = seq(1, 10, by = 0.5),
  mj = seq(0, 1, by = 0.05),
  c = 1,
  m = -0.2,
  k = 1,
  p = NA)

for (i in 1:nrow(HRS)){ 
  c = HRS[i,3]
  m = HRS[i,4]
  k = HRS[i,5]
  b = HRS[i,1]
  mj = HRS[i,2]
  p = 0.01 
  count <- 0 
  while(count < 10000){ 
    p<-p_next_timestep(c,k,m,mj,b,p) 
    count <- count + 1 
  }
  HRS[i,6] <- p 
}

ggplot(data = HRS, aes(mj,b,fill=p))+
  geom_tile(colour="black")+
  scale_fill_gradient(low="white", high="#1b7d4f")+
  guides(fill=guide_colourbar(title= "p"))

#our 3th data set (kestrels) will have m=5.5 and k=2 as set parameters
kestrels <- expand.grid(
  b = seq(1, 10, by = 0.5),
  mj = seq(0, 1, by = 0.05),
  c = 1,
  m = 5.5,
  k = 2,
  p = NA)

for (i in 1:nrow(kestrels)){ 
  c = kestrels[i,3]
  m = kestrels[i,4]
  k = kestrels[i,5]
  b = kestrels[i,1]
  mj = kestrels[i,2]
  p = 0.01 
  count <- 0 
  while(count < 10000){ 
    p<-p_next_timestep(c,k,m,mj,b,p) 
    count <- count + 1 
  }
  kestrels[i,6] <- p 
}

ggplot(data = kestrels, aes(mj,b,fill=p))+
  geom_tile(colour="black")+
  scale_fill_gradient(low="white", high="#e84d8a")+
  guides(fill=guide_colourbar(title= "p"))

#our 4th data set will have m=5 and k=10 as set parameters
b_mj_values4 <- expand.grid(
  b = seq(1, 10, by = 0.5),
  mj = seq(0, 1, by = 0.05),
  c = 1,
  m = 5,
  k = 10,
  p = NA)

for (i in 1:nrow(b_mj_values4)){ 
  c = b_mj_values4[i,3]
  m = b_mj_values4[i,4]
  k = b_mj_values4[i,5]
  b = b_mj_values4[i,1]
  mj = b_mj_values4[i,2]
  p = 0.01 
  count <- 0 
  while(count < 10000){ 
    p<-p_next_timestep(c,k,m,mj,b,p) 
    count <- count + 1 
  }
  b_mj_values4[i,6] <- p 
}

ggplot(data = b_mj_values4, aes(mj,b,fill=p))+
  geom_tile(colour="black")+
  scale_fill_gradient(low="white", high="#7f58af")+
  guides(fill=guide_colourbar(title= "p"))



#we now graph an adult mortality logic function for our different data set parameters, whicht allows us to analise the results of the heatmaps in relation to the form of the curves.


ggplot() + 
  xlim(0, 12)+
  ylim(0, 1)+
  geom_function(aes(colour="survivors"), fun = adult_mortality, args = list(k = 1, m = 3),  linewidth=1.25)+
  geom_function(aes(colour="HRS"), fun = adult_mortality,args = list(k = 1, m = -0.2),  linewidth=1.25)+
  geom_function(aes(colour="Kestrels"), fun = adult_mortality,args = list(k = 2, m = 5.5),  linewidth=1.25)+
  geom_function(aes(colour="big K"), fun = adult_mortality,args = list(k = 10, m = 5),  linewidth=1.25)+
  xlab("Clutch Size (C)")+
  ylab('Adult Mortality (ma)')+
  guides(colour=guide_legend(title= "Strategy"))+
  theme_bw(base_size = 15)+
scale_colour_manual(values = c("#7f58af", "#1b7d4f", "#e84d8a", "#feb326"))

#this concludes our code for the results seccion. deeper explanaition on the methods can be found on the method section of our paper, and analysis is found in the discusion seccion of our paper. thanks for reading! :)                                 