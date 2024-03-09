# Name - Hrushikesh Vazurkar, University ID - S2550941

# Assignment Overview

# This project involves mapping a relationship between measured effect
# of human androgen receptors(AR) with the dosage of fungicide vinclozolin.
# It involves modelling the dosage of vinclozolin(dosage or x) as latent variable
# following second order random walk model, with priors tau and tau0 as gamma(4,0.04)
# and M being uniform(0,4000) respectively. 

# The various steps taken in this project are:
# 0. Read vin.txt file containing the observed effects for each dosage level in each experiment.
# 1. Create a second order random walk model as per the specifications in 
#    the JAGS file "M1s2550941.jags".
# 2. Check for convergence:
#       -> Effective Sample Sizes(ESS): ESS for each stochastic variable
#          must be at least greater than 25% of the total iterations.
#       -> Plotting min and max ESS(Validation): The graphs for minimum
#          and maximum effective sample values must be steady with 
#          no discernible patterns
# 3. Effect calculation from model output and effect vs dose_index graph
#       -> Observed effect: Refer to the observed effect values 
#          from the data file.
#       -> Calculated effect: effect[i,j] <- M[j] * mu[i]
#       -> 95% interval: Use effect[i,j] to calculate the 
#          credible intervals
#       -> Graph plotting: Use ggplot for effect vs dose_index, 
#          using observed, calculated effect and credible intervals
# 4. Modify the model such that there is a mu value for 
#    each random walk variable(dosage or x), as per specifications in
#    the M2s2550941.jags file.
# 5. Repeat step 3 for the modified model from step 4
# 6. Model comparison through DIC sampling

# INPUT: NA
# OUTPUT: vin_data - A list containing data from vin.txt and 
#         params - used as input to the jags file for modelling:
#           -> no of experiments, no of doses and observed effect.
# PURPOSE: Preliminary step - getting relevant data 
#          from file for downstream usage
readData <- function(){
  vin_data <- read.table("vin.txt", header = TRUE, sep = " ")
  no_of_experiments <- length(unique(vin_data$exper)); no_of_doses <- nrow(vin_data)/no_of_experiments
  params <- list(nexp = no_of_experiments, ndoses = no_of_doses, effect=matrix(vin_data$effect, nrow=no_of_doses, ncol=no_of_experiments))
  
  return (list(
    vin_data=vin_data,
    params=params
  ))
}

# Constants: burn-in and iteration values for both models
burn_in <- 5000; total_iterations <- 40000

# INPUT: jagsfile - name of the JAGS file containing the model specs and
#        params - the output from readData function, to be used as input in 
#        jags file
# OUTPUT: model - the updated model after burn-in and 
#         samples - used to calculate effect values for plotting.
# PURPOSE: Create models using coda library and generate model samples
#             -> Create model as per specified jagsfile
#             -> Get stochastic variables like M,mu,dosage(x),tau and tau0 
#                from model
createModel <- function(jagsfile, params){
  library(rjags)
  model <- jags.model(jagsfile, data=params)
  # Updating the model with burn-in iterations(which are discarded).
  # Allows the model to settle in the proper region(as per posterior distn.)
  update(model, burn_in)
  # Getting the samples from the model
  samples <- coda.samples(model, c("M", "mu", "x", "tau", "tau0"), n.iter = total_iterations)
  return (
    list(model=model, samples=samples)
  )
}

# NOTE: This function handles the effect calculation for both the models,
#       and has dedicated code to change mu valus depending on the model name.
# INPUT: modelName - the name of the model - either model1 or model2
# OUTPUT: No return value, just print plot of effect(observed and calculated)
#         versus Dose_index(1 to 9), with 95% credible intervals
# PURPOSE: Plotting the graph as per requirements
plotEffectGraph <- function(modelName, samples){
  
  # Calculated Effect calculations
  # OVERVIEW:
  #   What ?
  #   -> For each experiment, dosage combo -> allot effect = Mj * mui
  #   -> Hence, effect or output from the sample is stored as 3D matrix of 9 X 5 X total_iterations
  #   -> Created an effect_compressed -> with mean of each column stored at effect[i,j,] index
  #   Why ?
  #   -> We need mean values for each (experiment, dose), and we need 
  #   credible regions for each (experiment, dose). Using effect_compressed 
  #   for mean values and effect for credible region calculation.
  
  # effect_output: 3D matrix of 9X5Xtotal_iterations
  effect_output <- array(0, c(no_of_doses,no_of_experiments,total_iterations))
  for(j in 1:no_of_experiments){
    for(i in 1:no_of_doses){
      Mj <- model1_samples[[1]][,paste0('M[',j,']')] #M[j] from model sample
      # mu changes from model1 to model2.
      # if model1 -> mu[i], else if model2 -> mu[i,j]
      mu <- ifelse(modelName=="model1", samples[,paste0('mu[',i,']')],
                    samples[,paste0('mu[',i,',',j,']')])
      effect_output[i,j,] <- Mj * mu
    }
  }
  
  # effect_compressed, credible interval lower and upper: 
  # vector of 45(for ease of plotting)
  effect_output_compressed <- cred_lower <- cred_upper <- rep(0, no_of_doses*no_of_experiments)
  for(j in 1:no_of_experiments){
    for(i in 1:no_of_doses){
      # Converted from a 2D representation to vector for ease of plotting
      effect_output_compressed[(j-1)*9+i] <- mean(effect_output[i,j,])
      # lower and upper credible interval
      cred_lower[(j-1)*9+i] <- quantile(effect_output[i,j,],probs = 0.025)
      cred_upper[(j-1)*9+i] <- quantile(effect_output[i,j,],probs = 0.975)
    }
  }
  
  library(ggplot2)
  
  colors <- c("red","blue","green","purple","orange")
  # scaling colors across groups
  colScale <- scale_colour_manual(name = "colour_codes",values = colors)
  # scaling colors across fill -> for credible interval bands
  fillScale <- scale_fill_manual(name="colour_codes", values=colors)
  # dataframe for plotting
  data2 <- data.frame(
    x=rep(1:9), # dose_index: 1 to 9
    y1=vin_data$effect, # observed effect
    y2=effect_output_compressed, # calculated effect
    cred_lower=cred_lower, # lower and upper credible intervals
    cred_upper=cred_upper,
    colour_codes=rep(c("Experiment1","Experiment2","Experiment3","Experiment4","Experiment5"), each=9)
  )
  # line: observed effect, dot: calculated effect
  # ribbon: credible interval bands
  # observed effect
  observed <- ggplot(data = data2, aes(x=x, y=y1, color=colour_codes)) + geom_point()
  # observed effect + calculated
  observed_calculated <- observed + geom_line(aes(x=x, y=y2,color=colour_codes))
  # observed effect + calculated + credible regions
  observed_calculated_cred <- observed_calculated + geom_ribbon(aes(x=x,ymin=cred_lower,ymax=cred_upper,fill=colour_codes), alpha=0.2, linetype="dotted")
  # add formatting
  formatted <- observed_calculated_cred + colScale + fillScale + labs(title = paste("Effect vs Dose Index -",modelName), x = "Dose Index", y = "Effect")
  print(formatted)
}

# Step0: Read and extract data from vin.txt
data <- readData()
vin_data <- data$vin_data; no_of_experiments <- data$params$nexp;
no_of_doses <- data$params$ndoses; params <- data$params

# Step1: Create model1 from and obtain samples
sorw_model <- createModel("M1s2550941.jags", params)
model1 <- sorw_model$model; model1_samples <- sorw_model$samples; 

# Step2: Convergence Test
# Using coda library, get ESS for stochastic variable:
#   -> Print ESS: min_ESS_value > 25% of total_iterations and
#      max_ESS_value << total_iterations
#   -> Get max and min ESS value variables and plot them.
#      The graphs should have be steady without any discernible pattern
library(coda)
print("Effective Sample Size:")
print(effectiveSize(model1_samples))

# Non stochastic variables give constant values, not suited for ESS plots.
# Hence, discarding variables like mu from ESS values and plots
non_stochastic <- c(grep("mu",colnames(model1_samples[[1]]))) # find all mu values
ess <- effectiveSize(model1_samples)[-non_stochastic] # remove all mu values
max_ess <- names(ess)[which.max(ess)]; min_ess <- names(ess)[which.min(ess)]

# max and min ESS plots
print(plot(model1_samples[, max_ess], main = paste("Max: Trace Plot for -", max_ess)))
print(plot(model1_samples[, min_ess], main = paste("Min: Trace Plot for -", max_ess)))

print("Inference:")
cat(strwrap("The ESS values show a good indication of convergence, and the plots for max and min ESS are fairly steady with no discernible pattern, for the current burn-in and total iterations.", width=80), sep="\n")

# Step3: From model1 samples(output), calculate effect and plot
#        graph of Effect(observed and calculated) vs Dose_index(1 to 9) with
#        95% credible interval regions
plotEffectGraph("model1", model1_samples[[1]])

# Step4: Modified model2
mod_sorw_model <- createModel("M2s2550941.jags",params)
model2 <- mod_sorw_model$model; model2_samples <- mod_sorw_model$samples 

# Step5: From model2 samples(output), calculate effect and plot
#        graph of Effect(observed and calculated) vs Dose_index(1 to 9) with
#        95% credible interval regions
plotEffectGraph("model2", model2_samples[[1]])

# Step6: DIC sampling for model1 and model2 comparison
print("Step6: DIC sampling for model1 and model2 comparison")
dic_model <- jags.model("M1s2550941.jags",data=params,n.chains=7)
dic_model2 <- jags.model("M2s2550941.jags",data=params,n.chains=7)

print(dic.samples(dic_model,n.iter=total_iterations))
print(dic.samples(dic_model2,n.iter=total_iterations))

print("Inference from DIC Sampling Test:")
cat(strwrap("Since the mean and penalized deviance from model2 is lower than model1, model2 is more preferable than model1",width=80),sep="\n")