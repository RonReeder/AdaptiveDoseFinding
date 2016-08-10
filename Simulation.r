#  Note, you will need R 3.3.1 with rjags package, and you will also need to install a program called rjags4 outside of R: http://www.sourceforge.net/projects/mcmc-jags/files

N=15000
DurationOfTrialMonths = 50
Allocation = c(.25,.25,.25,.25)
ARMMeans=c(12,9,6,3)
ARMStandardDeviations=c(12,15,10,8)
Simulations = 1
set.seed(9)

MRSDistribution=matrix(,4,7)
MRSDistribution[1,] = c(.1,.1,.1,.1,.1,.1,.4)
MRSDistribution[2,] = c(.1,.1,.1,.1,.1,.4,.1)
MRSDistribution[3,] = c(.1,.1,.1,.4,.1,.1,.1)
MRSDistribution[4,] = c(.4,.1,.1,.1,.1,.1,.1)

UtilityWeighting = c(1,.91,.76,.65,.33, 0, 0) # for MRS

#Number of contols / actives in future phase three
m0=250
m1=250

library(rjags)


Allocate = function(ARMMeans, ARMStandardDeviations, Allocation, MRSDistribution){

	arm = sample(x=4,size=1,replace=TRUE, prob=Allocation) #treatment arm the subject is randomized to
	penumbra = rnorm(mean=ARMMeans[arm], sd=ARMStandardDeviations[arm], n=1) #penumbra outcome for each patient
	MRS = sample(x=7,size=1,replace=TRUE, prob=MRSDistribution[arm,])-1 #MRS outcome for each patient
	UtilityWeightedMRS = UtilityWeighting[MRS+1]

	c(
		arm
		, penumbra
		, UtilityWeightedMRS
	)

}

Fit = function(Subject, column=3){
	### Make dummy variables
	dose0 <- as.numeric(Patients[1:Subject,2]==1)
	dose1 <- as.numeric(Patients[1:Subject,2]==2)
	dose2 <- as.numeric(Patients[1:Subject,2]==3)
	dose3 <- as.numeric(Patients[1:Subject,2]==4)
	#cbind(Patients[1:Subject,2], dose0, dose1, dose2, dose3)

	df <- list(N=40, x=Patients[1:Subject,column], dose0=dose0, dose1=dose1, dose2=dose2, dose3=dose3)

	model1ForJAGS <- "model1.txt"
	cat("
	model {
	    for(i in 1:N){
	        x[i] ~ dnorm(mu[i],prec.x)  
	        mu[i] <- beta[1]*dose0[i] + beta[2]*dose1[i] + beta[3]*dose2[i] + beta[4]*dose3[i]
	    }
	    beta[1] ~ dnorm(10, prec.b)#prior for control
	    beta[2] ~ dnorm(beta[1], prec.b)
	    beta[3] ~ dnorm(beta[2], prec.b)
	    beta[4] ~ dnorm(beta[3], prec.b)
	    prec.b ~ dgamma(1/2, 1/2) #prior for the drift parameter - variation from dose to dose.  Note precision and not SD
	    prec.x ~ dgamma(1/2, 1/2) #prior for individual level pt variability.  Note precision and not SD
	    tau.b <- 1/sqrt(prec.b) 
	    tau.x <- 1/sqrt(prec.x)
	}
	",file=model1ForJAGS)


	inits.1 <- list(beta=c(1,1,1,1), prec.b=1, prec.x=1)
	jags1 <- jags.model("model1.txt",data=df,n.chains=1,n.adapt=1000, inits=inits.1)
	update(jags1, n.iter=10000)
	mcmc.samples <- coda.samples(jags1,c('beta','tau.b', 'tau.x'),10000)
  	
	Probabilities = apply(Samples[,2:4] == apply(Samples[,2:4], 1, min),2,mean)

	return(
		list(
			Samples = mcmc.samples[[1]][,1:4]
			, Probabilities = Probabilities
			, Allocation = c(.25,0.75*Probabilities[1]/sum(Probabilities[1:3]), 0.75*Probabilities[2]/sum(Probabilities[1:3]), 0.75*Probabilities[3]/sum(Probabilities[1:3]))
			, BetterThan1 = apply(Samples < Samples[,1],2,mean)
			, BestDose = 1+which.max(Probabilities)
		)
	)
}

#first dimension is for trial
#second is for variable
#third is for look
Results = array(,c(Simulations,22,4))
for (Trial in 1:Simulations)
{
	#Initialize Patient matrix
	Patients = matrix(,N,10)

	#Column 1 is the time from the start of the trial to the enrollment of the patient
	Patients[,1] = sort(runif(n=N, min = 0, max = DurationOfTrialMonths))

	LookIndex = 0
	for (Subject in 1:N){
		Patients[Subject, 2:4] = Allocate(ARMMeans=ARMMeans, ARMStandardDeviations=ARMStandardDeviations, Allocation = Allocation, MRSDistribution=MRSDistribution)
		
		if (Subject %in% c(40, 60,80, N)) {
			LookIndex = LookIndex + 1
			print(LookIndex)
			PenumbraModel=Fit(Subject, column = 3)
			Allocation = PenumbraModel$Allocation
			#Futile = (max(PenumbraModel$BetterThan1) < Cutoff)
			MRSModel=Fit(Subject,column = 4)

			#store results
				#	Number allocated to 1
				Results[Trial,1,LookIndex] = mean(Patients[1:Subject,2]==1)			
				#	Number allocated to 2
				Results[Trial, 2,LookIndex] = mean(Patients[1:Subject,2]==2)		
				#	Number allocated to 3
				Results[Trial, 3,LookIndex] = mean(Patients[1:Subject,2]==3)		
				#	Number allocated to 4
				Results[Trial, 4,LookIndex] = mean(Patients[1:Subject,2]==4)	
	
				#	Observed effect in dose 1 - 4
				Results[Trial, 5:8,LookIndex] =  tapply(Patients[1:Subject,3], Patients[1:Subject,2], mean)

				#	Fitted effect in dose 1 - 4
				Results[Trial, 9:12,LookIndex] = apply(PenumbraModel$Samples, 2, mean)
				
				#	Probability of being best
				Results[Trial, 13:15,LookIndex] = PenumbraModel$Probabilities

				#	Probability that best is better than control
				Results[Trial, 16,LookIndex] = max(PenumbraModel$BetterThan1)

				#	Allocation going forward
				Results[Trial, 17:20,LookIndex] = PenumbraModel$Allocation

				#	Probability that best (best based on penumbra) is better than control (better is based on MRS)
				p = MRSModel$BetterThan1[PenumbraModel$BestDose]
				Results[Trial, 21,LookIndex] = p

				#	Probability of success in a future phase 3 trial with utility weighted MRS outcome
				n0 = sum(Patients[1:Subject,2]==1) #number of controls				
				n1 = sum(Patients[1:Subject,2]==PenumbraModel$BestDose) #number in best dose (best by penumbra)
				Results[Trial, 22,LookIndex] = pnorm(
						(qnorm(p)*sqrt(1/n0+1/n1)-qnorm(.975)*sqrt(1/m0+1/m1))
							/sqrt(1/n0+1/n1+1/m0+1/m1)
					)
				
		}
	}
	
}

#Observed
Results[1,5:8,4]

#Fitted
Results[1,9:12,4]

#Observed allocation of subject
Results[1,1:4,4]