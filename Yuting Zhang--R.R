# Exercise 10-Q1
#Dynamic modeling used to explore cancer biology
#A population model to investigate evolution of drug resistance in tumors:
#Imagine a cancer cell in a tumor that spontaneously exhibited a mutation that confers drug resistance;
#The mutation does not have any positive or negative effects on growth rate of that sub-population when the cancer drug is absent;
#when the cancer drug is present the mutant sub-population grows at 50% of its growth rate in the absence of the drug and the non-mutant sub-population declines rapidly;

#in the absense of the cancer drug that the cells grow at a rate of 0.1 per day:
rN = 0.1
rM = 0.1

#the carrying capapcity (K) of the tumor is one million cells:
K=1e6

#Drug treatment of non-mutant cells results in a negative growth rate of -0.1;
dM = -0.1
#when the cancer drug is present the mutant sub-population grows at 50% of its growth rate in the absence of the drug:
dN = 0.5

#The mutation of a single cell occurred early in the tumor growth and when it occurred there were 100 total cells in the tumor:
N0 = 99
M0 = 1
T0 = 100

time = 500

#create vector to store N's and set initial N
Ns= numeric(length=time)
Ns[1]=N0
Ms= numeric(length=time)
Ms[1]=M0
Totals = numeric(length=time)
Totals[1]=T0

#simulate

for (t in 1:(time-1)){
  Ns[t+1] <- Ns[t]+rN*dN*Ns[t]*(1-(Ns[t]+Ms[t])/K)
  Ms[t+1] <- Ms[t]+rM*dM*Ms[t]*(1-(Ns[t]+Ms[t])/K)
  Totals[t+1] <- Ns[t+1] + Ms[t+1]
}

#plot simulation:
library(ggplot2)
simEvents <- data.frame(time=1:length(Totals),To = Totals)
ggplot(data=simEvents,aes(x=time,y=To))+geom_line()


# when time = 350 days, the growth of the two sub-populations in the tumor is equilibrium followed by drug treatment. 

