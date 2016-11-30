#  Note, you will need R 3.3.1 with rjags package, and you will also need to install a program called rjags4 outside of R: http://www.sourceforge.net/projects/mcmc-jags/files

#Speficy trial characteristics and sassumed distributions.
	#Maximum number of subjects in a trial.
		N=150

	#Duration of trial (not used in calculation because response is measured almost immediately)
		DurationOfTrialMonths = 50

	#Initial allocation (ranomization probabilities)	
		Allocation = c(.25,.25,.25,.25)

	#Number of simulations to run
		Simulations = 100

	#Distribution of penumbra outcomes in each arm
		ARMMeans=c(21,9,6,3)
		ARMStandardDeviations=c(12,15,10,8)

	#Distribution of MRS for each dosing group.
		MRSDistribution=matrix(,4,7)
		MRSDistribution[1,] = c(.09,.13,.09,.13,.16,.09,.31) #Placebo
		MRSDistribution[2,] = c(.12,.16,.11,.15,.15,.12,.19) #Dose 1 
		MRSDistribution[3,] = c(.15,.20,.10,.16,.15,.10,.14) #Dose 2
		MRSDistribution[4,] = c(.17,.24,.15,.09,.20,.04,.11) #Dose 3
		#MRSDistribution[1,] = c(0,1,0,0,0,0,0) #Placebo
		#MRSDistribution[2,] = c(0,1,0,0,0,0,0) #Dose 1 
		#MRSDistribution[3,] = c(0,1,0,0,0,0,0) #Dose 2
		#MRSDistribution[4,] = c(0,1,0,0,0,0,0) #Dose 3

	#Number of contols / actives in future phase three
		m0=25 #controls
		m1=25 #actives

	#Utility function for utility-weighted MRS.
		UtilityWeighting = c(1,.91,.76,.65,.33, 0, 0) 

	#The points at which interim evaluations are made (in terms of the number of subjects)
	#N is include so that a 'look' is also performed at the end of the trial.
		InterimLooks = c(40,80,120,N)

#Set seed for generating random data so that results are reproducable.
	set.seed(1)

#Load library for interface with RJags.
	library(rjags)

#This function is used to randomly determine the treatment a subject will receive as well as the subjects penumbra and MRS outcomes.
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

#	This function takes the Patient matrix (rows 1 through 'Subject'), fits a NDLR model, and estimates the probability that each active dose is the best
#	'Best' can mean in terms of the penumbra outcome (column = 3) or the modified MRS outcome (column = 4)
##### Output generated ######
	#Samples is a matrix of data generated based on the fitted distributions.
	#Probabilities is a vector of probabilities that each dose is the best.  This vector does not include placebeo, so it has only 3 elements.
	#BetterThan1 is a vector of length 4 (includes placebo) indicating the probability that each dose is better than placebo.  Note that this probability is zero for placebo.
	#BestDose is an integer representing which dose is estimated to be the best.  1 = placebo, etc.
		Fit = function(Patients,Subject, column, good){
			### Make dummy variables
			dose0 <- as.numeric(Patients[1:Subject,2]==1)
			dose1 <- as.numeric(Patients[1:Subject,2]==2)
			dose2 <- as.numeric(Patients[1:Subject,2]==3)
			dose3 <- as.numeric(Patients[1:Subject,2]==4)
			#cbind(Patients[1:Subject,2], dose0, dose1, dose2, dose3)

			df <- list(N=Subject, x=Patients[1:Subject,column], dose0=dose0, dose1=dose1, dose2=dose2, dose3=dose3)

			model1ForJAGS <- "model1.txt"
			cat("
			model {
			    for(i in 1:N){
			        x[i] ~ dnorm(mu[i],prec.x)  
			        mu[i] <- beta[1]*dose0[i] + beta[2]*dose1[i] + beta[3]*dose2[i] + beta[4]*dose3[i]
			    }
			    beta[1] ~ dnorm(20, prec.b)#prior for control
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
		  	Samples = mcmc.samples[[1]][,1:4]
		  	if (good == 'max'){Probabilities = apply(Samples[,2:4] == apply(Samples[,2:4], 1, max),2,mean)}
		  	if (good == 'min'){Probabilities = apply(Samples[,2:4] == apply(Samples[,2:4], 1, min),2,mean)}
		  	if (good == 'max'){BetterThan1 = apply(Samples > Samples[,1],2,mean)}
		  	if (good == 'min'){BetterThan1 = apply(Samples < Samples[,1],2,mean)}
			return(
				list(
					Samples = Samples
					, Probabilities = Probabilities
					, Allocation = c(.25,0.75*Probabilities[1]/sum(Probabilities[1:3]), 0.75*Probabilities[2]/sum(Probabilities[1:3]), 0.75*Probabilities[3]/sum(Probabilities[1:3]))
					, BetterThan1 = BetterThan1
					, BestDose = 1+which.max(Probabilities)
				)
			)
		}

#Now, we run several simulations of the trial in a loop.  The results will be store in the 'Results' array.
	#first dimension is for trial
	#second is for variable
	#third is for look
		Results = array(,c(Simulations,23,4))

	#Simulate trials
		pdf(file = "C:/Users/rreeder/Desktop/AdaptiveDoseFinding/TrialPlot.pdf", onefile = T)
		for (Trial in 1:Simulations)
		{
			#Initialize Patient matrix
				Patients = matrix(,N,4)

			#Column 1 is the time from the start of the trial to the enrollment of the patient
			#This variable isn't really used currently.
				Patients[,1] = sort(runif(n=N, min = 0, max = DurationOfTrialMonths))

			#The allocation probabilities are updated at interim looks.  LookIndex keeps track of which look we are on.
				LookIndex = 0

			#	Generate data for each subject.	
				for (Subject in 1:N){
					Patients[Subject, 2:4] = Allocate(ARMMeans=ARMMeans, ARMStandardDeviations=ARMStandardDeviations, Allocation = Allocation, MRSDistribution=MRSDistribution)
					
					#If we are at an interim look, update allocation probabilities, etc, and store them in the Results array.
						if (Subject %in% InterimLooks) {
							LookIndex = LookIndex + 1
							print(LookIndex)
							PenumbraModel=Fit(Patients = Patients, Subject=Subject, column = 3, good = 'min')
							Allocation = PenumbraModel$Allocation
							#Futile = (max(PenumbraModel$BetterThan1) < Cutoff)
							MRSModel=Fit(Patients = Patients, Subject=Subject,column = 4, good = 'max')

							#store results
								#	Proportion allocated to 1
								Results[Trial,1,LookIndex] = mean(Patients[1:Subject,2]==1)			
								#	Proportion allocated to 2
								Results[Trial, 2,LookIndex] = mean(Patients[1:Subject,2]==2)		
								#	Proportion allocated to 3
								Results[Trial, 3,LookIndex] = mean(Patients[1:Subject,2]==3)		
								#	Proportion allocated to 4
								Results[Trial, 4,LookIndex] = mean(Patients[1:Subject,2]==4)	
					
								#	Observed effect in dose 1 - 4
								Results[Trial, 5:8,LookIndex] =  tapply(Patients[1:Subject,3], factor(Patients[1:Subject,2], levels=1:4), mean)

								#	Fitted effect in dose 1 - 4
								Results[Trial, 9:12,LookIndex] = apply(PenumbraModel$Samples, 2, mean)
								
								#	Probability of being best
								Results[Trial, 13:15,LookIndex] = PenumbraModel$Probabilities

								#	Probability that best is better than control
								Results[Trial, 16,LookIndex] = max(PenumbraModel$BetterThan1)

								#	Allocation going forward
								Results[Trial, 17:20,LookIndex] = PenumbraModel$Allocation

								#	Probability that best (best based on penumbra) is better than control (better is also based on penumbra)
								ProbBestBetterThanControlPenumbra = PenumbraModel$BetterThan1[PenumbraModel$BestDose]
								Results[Trial, 21,LookIndex] = ProbBestBetterThanControlPenumbra

								#	Probability best (based on penumbra) is better than control (better is based on Weighted MRS)
								ProbBestBetterThanControlMRS = MRSModel$BetterThan1[PenumbraModel$BestDose]
								Results[Trial, 23,LookIndex] = ProbBestBetterThanControlMRS

								#	Probability of success in a future phase 3 trial with utility weighted MRS outcome
								n0 = sum(Patients[1:Subject,2]==1) #number of controls				
								n1 = sum(Patients[1:Subject,2]==PenumbraModel$BestDose) #number in best dose (best by penumbra)
								Results[Trial, 22,LookIndex] = pnorm(
										(qnorm(ProbBestBetterThanControlMRS)*sqrt(1/n0+1/n1)-qnorm(.975)*sqrt(1/m0+1/m1))
											/sqrt(1/n0+1/n1+1/m0+1/m1)
									)
								
							# plot
							if (Trial %in% c(1,2,3,4))
							{
								#MRS plot
								best <- c(NA, PenumbraModel$Probabilities)
								par(mar = c(4,4,4,4))
								x.spot <- barplot(Subject*Results[Trial,1:4,LookIndex], ylim = c(0,70), xlim = c(0, 4.75), names = c("Cntrl", "Dose1", "Dose2", "Dose3"))
								mtext(side = 2, "Sample Size", line = 2.5)	
								par(new = TRUE)
								plot(x.spot,  tapply(Patients[,3], factor(Patients[,2], levels = 1:4), mean), xlim = c(0, 4.75), ylim = c(0, 30), pch = 8, col = "red", 
										bty = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = 1.5)
								title(main = paste0("Simulated Trial: ", Trial, "; Number Enrolled: ", Subject))		
								axis(side = 4)	
								mtext(side = 4, "Mean Change in Penumbra Outcome", line = 2)	
								lines(x.spot, apply(PenumbraModel$Samples, 2, mean)[1:4], type = "b", lwd = 2)
								lines(x.spot, apply(PenumbraModel$Samples, 2, quantile, 0.025)[1:4], type = "b", lwd = 1, lty = 2)
								lines(x.spot, apply(PenumbraModel$Samples, 2, quantile, 0.975)[1:4], type = "b", lwd = 1, lty = 2)
								text(x = x.spot, y = rep(-5,4), best, xpd = TRUE)
								text(x = x.spot[1], y = -5, "Pr(Best)", xpd = TRUE)
								legend("topleft", legend = c(paste0("Pr(Best > Cntrl on Penumbra Change) = ", ProbBestBetterThanControlPenumbra),paste0("Pr(Phase III Success on Weighted MRS) = ",round(Results[Trial, 22,LookIndex],3))), bty = "n")

								#Penumbra Plot
								best <- c(NA, MRSModel$Probabilities)
								par(mar = c(4,4,4,4))
								x.spot <- barplot(Subject*Results[Trial,1:4,LookIndex], ylim = c(0,70), xlim = c(0, 4.75), names = c("Cntrl", "Dose1", "Dose2", "Dose3"))
								mtext(side = 2, "Sample Size", line = 2.5)	
								par(new = TRUE)
								plot(x.spot,  tapply(Patients[,4], factor(Patients[,2], levels = 1:4), mean), xlim = c(0, 4.75), ylim = c(0, 1), pch = 8, col = "red", 
										bty = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = 1.5)
								title(main = paste0("Simulated Trial: ", Trial, "; Number Enrolled: ", Subject))		
								axis(side = 4)	
								mtext(side = 4, "Utility-Weighted MRS", line = 2)	
								lines(x.spot, apply(MRSModel$Samples, 2, mean)[1:4], type = "b", lwd = 2)
								lines(x.spot, apply(MRSModel$Samples, 2, quantile, 0.025)[1:4], type = "b", lwd = 1, lty = 2)
								lines(x.spot, apply(MRSModel$Samples, 2, quantile, 0.975)[1:4], type = "b", lwd = 1, lty = 2)
								text(x = x.spot, y = rep(-0.15,4), best, xpd = TRUE)
								text(x = x.spot[1], y = -0.15, "Pr(Best)", xpd = TRUE)
								legend("topleft", legend = c(paste0("Pr(Best > Cntrl on MRS) = ", ProbBestBetterThanControlMRS),paste0("Pr(Phase III Success on Weighted MRS) = ",round(Results[Trial, 22,LookIndex],3))), bty = "n")

								
							}

						}
				}
		
		}
dev.off()

#Print a few things to verify results are sensible.
	#Observed
	Results[1,5:8,4]

	#Fitted
	Results[1,9:12,4]

	#Observed allocation of subject
	Results[1,1:4,4]

	#Prob best better than control on penumbra
	Results[, 21,4]

	#Prob success on Phase III based on MRS
	Results[, 22,4]

	#Prob best better than control on MRS
	Results[, 23,4]