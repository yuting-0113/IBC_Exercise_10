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

#Drug treatment of non-mutant cells(Normal cells) results in a negative growth rate of -0.1;
dN = -0.1
#when the cancer drug is present the mutant sub-population grows at 50% of its growth rate in the absence of the drug:
dM = 0.05

#The mutation of a single cell occurred early in the tumor growth and when it occurred there were 100 total cells in the tumor:
N0 = 99
M0 = 1

#Set up the total timepoint:
time = 750

#create vector to store N's and set initial N
Ns= numeric(length=time)
Ns[1]=N0
Ms= numeric(length=time)
Ms[1]=M0
Drugs = character(length=time)
Drugs[1] <- "No drug"

#simulate

for (t in 1:(time-1)){
  if (t < 250){     #Cells are in the absent of drug treatment before Day250;
    Ns[t+1] <- Ns[t]+rN*Ns[t]*(1-(Ns[t]+Ms[t])/K)
    Ms[t+1] <- Ms[t]+rM*Ms[t]*(1-(Ns[t]+Ms[t])/K)
    Drugs[t+1] <- "No drug"
  }else{            #Cells are in the presence of drug treatment after Day250;
    Ns[t+1] <- Ns[t]+dN*Ns[t]*(1-(Ns[t]+Ms[t])/K)
    Ms[t+1] <- Ms[t]+dM*Ms[t]*(1-(Ns[t]+Ms[t])/K)
    Drugs[t+1] <- "drug treatment"
  }
}

#plot simulation:
library(ggplot2)

#Create dataframe for N and M:
NEvent <- data.frame(time=1:length(Ns),abundance=Ns, treatment=Drugs)
NEvent$treatment <- sub("^", "Normal,", NEvent$treatment )
MEvent <- data.frame(time=1:length(Ms),abundance=Ms, treatment=Drugs)
MEvent$treatment <- sub("^", "Mutant,", MEvent$treatment )
simEvents <- rbind(NEvent , MEvent)   #append the N with M data.

#plot:
ggplot(data=simEvents,aes(x=time,y=abundance,color=treatment))+
  geom_line() + 
  xlab("Time") + ylab("Cell number") 


# If treat the cells with drug at the Day 250, on Day650 equilibrium. 

