set.seed(123)#for reproducibility

# Set your working directory to the DMC folder  to load the model
wd <- here::here();

setwd("/Users/ale/ownCloud/3-R/DMC_190819/")
source("dmc/dmc.R")
load_model (dir_name="LBA",model_name="lba_B.R")

# Go back to the project's folder
setwd(wd)

pacman::p_load('tidyverse','magrittr')

# Set up the model -------------------------------------------------------------------------------

# Define factorial design 2 directions x 2 coherence levels x 2 angle difficulty levels
# s1/s2 left and right motion direction
# cohE/cohH easy/hard coherence
# angleE/angleH easy/hard decision boundary
factors <- list(S=c("s1","s2"),Coh=c("cohE","cohH"),Angle=c("angleE","angleH"))
# specify the responses
responses <- c("r1","r2") # keys o and p

# define a new factor M (Match) with correct and incorrect responses for each factor level
match.map <- list(M=list(s1="r1",s2="r2")) 

# define constants to allow parameter identification
const <- c(sd_v.false=1,st0=0)

# Specify how the factors map onto the parameters
# boundary  changes with manipulations, accum-rate faster for correct
p.map <- list(A="1",B=c("Coh","Angle"),t0="1",mean_v=c("M"),sd_v="M",st0="1")

model1 <- model.dmc(p.map,match.map = match.map,constants=const,
                   responses=responses,factors=factors)

p.vector  <- c(A=0.5, B.cohE.angleE = 3, B.cohE.angleH = 4, B.cohH.angleE = 5, B.cohH.angleH = 6,t0=.2,
               mean_v.true=2,mean_v.false=-1,sd_v.true=0.5)

check.p.vector(p.vector,model1)
print.cell.p(p.vector,model1)


# Load and preprocess the data ------------------

#data <- R.matlab::readMat(here::here('Data','lba_test_data.mat'))
data <- read.csv("AccrateRT.csv") %>% as.data.frame()


data <- simulate.dmc(p.vector,model1,n=1) 

# this needs to be changed 
# get indices of  correct  trials and remove implausible RTs
data <- data %>% mutate(Correct = ifelse((S=="s1" & R == "r1") | (S=="s2" & R == "r2"),"Correct","Wrong"))  %>%
        mutate(Plausibility = case_when(
                        Correct == "Correct" & RT < 0.25 ~ "Too quick",
                        Correct == "Correct" & RT > 4.5 ~ "Too slow",
                         TRUE ~ "Good")) %>%
        filter(Plausibility == "Good",.preserve = TRUE)

# remove outlier RTs (mean ± 2.5sd) on subject-by-subject basis
data <- data %>% group_by(s) %>% mutate(Plausibility = case_when(
  
  Correct == "Correct" & RT < (mean(RT)-2*sd(RT)) ~ "Outlier",
  Correct == "Correct" & RT > (mean(RT)+2*sd(RT)) ~ "Outlier",
  TRUE ~ "Good")) %>% 
  # filter out outliers and remove  the last two columns
  ungroup() %>% filter(Plausibility == "Good",.preserve = TRUE) %>%
  select(-c("Correct","Plausibility")) 


# Ensure that subject and factor columns are encoded as factors
data <- data %>% mutate( s = as.factor(s),
                         S = as.factor(S),
                         Coh   = as.factor(Coh),
                         Angle = as.factor(Angle),
                         R = as.factor(R))

# bind the data to the model object
data.model <- data.model.dmc(data, model)


# Fit the model non-hierarchically (i.e. treating subjects as fixed effects) ------
# number of cores available
mycores <- 4

# First we set  up our model fitting object
sampling_LBA <- h.samples.dmc(
  nmc = 400, # number of samples left AFTER thinning
  p.prior = p.prior, # list of priors
  thin = 10, #thinning to deal with autocorrelation
  data = data.model # note: Default number of chains is 3 * nparams
)

# Run the model using the automatic function h.RUN.dmc
sampling_LBA <- h.RUN.dmc(sampling_LBA, cores = mycores, verbose = TRUE)

# Save this sampling object
save(sampling_LBA, file =  "sampling_LBA.RData")





# Each subject can be thought of as a random effect. That is, rather than 
# thinking of each subject as being completely unrelated to all other subjects
# we treat them as coming from the same population. In practice this means that
# each subject's parameters are from a population distribution.

# This can be done with the "p.prior" (parameter prior) argument to 
# h.simulate.dmc where the prior.p.dmc function defines 
# distributions from which parameters can be sampled. We will use the default 
# distribution (normal), requiring a mean (p1) and standard deviation (p2) 
# parameter. Here we use p.vector to give the means and the same sd for 
# all parameters (a vector for p2 allows different sds for each parameter).

pop.mean <- c(A=0.5, B.cohE.angleE = 3, B.cohE.angleH = 4, B.cohH.angleE = 5, B.cohH.angleH = 6,
              mean_v.true=2,mean_v.false=-1,sd_v.true=0.5, t0=.2)
pop.scale <-c(A=.1, B.cohE.angleE = .1, B.cohE.angleH = .1, B.cohH.angleE = .1, B.cohH.angleH = .1,
              mean_v.true=.2,mean_v.false=.2,sd_v.true=0.1,t0=.05)

pop.prior <- prior.p.dmc(
  dists = rep("tnorm",9),
  p1=pop.mean,p2=pop.scale,
  lower=c(0,0,0,0,0,NA,NA,0,.1),upper=c(NA,NA,NA,NA,NA,NA,NA,NA,1)
)

##  Check population distributions
par(mfcol=c(2,5)); for (i in names(pop.prior)) plot.prior(i,pop.prior)


# Simulate some data: ns subjects, with n data point per design cell
raw.data1 <- h.simulate.dmc(model1, p.prior = pop.prior, n = 250, ns = 40)
data.model1 <- data.model.dmc(raw.data1, model1)

# Take a look at the first  subject's data
par(mfcol=c(4,2)) 
for (i in 1) { # First column=response left, Second column = response right. Rows = EasyCoh-EasyAngle, EasyCoh-HardAngle,HardCoh-EasyAngle, hard-hard
  plot.cell.density(data.cell=data.model1[[i]][data.model1[[i]]$Coh=="cohE" & data.model1[[i]]$Angle=="angleE" & data.model1[[i]]$S=="s1",],C="r1",xlim=c(1,6),show.mean = TRUE)
  plot.cell.density(data.cell=data.model1[[i]][data.model1[[i]]$Coh=="cohE" & data.model1[[i]]$Angle=="angleH" & data.model1[[i]]$S=="s1",],C="r1",xlim=c(1,6),show.mean = TRUE)
  plot.cell.density(data.cell=data.model1[[i]][data.model1[[i]]$Coh=="cohH" & data.model1[[i]]$Angle=="angleE" & data.model1[[i]]$S=="s1",],C="r1",xlim=c(1,6),show.mean = TRUE)
  plot.cell.density(data.cell=data.model1[[i]][data.model1[[i]]$Coh=="cohH" & data.model1[[i]]$Angle=="angleH" & data.model1[[i]]$S=="s1",],C="r1",xlim=c(1,6),show.mean = TRUE)

  plot.cell.density(data.cell=data.model1[[i]][data.model1[[i]]$Coh=="cohE" & data.model1[[i]]$Angle=="angleE" & data.model1[[i]]$S=="s2",],C="r2",xlim=c(1,6),show.mean = TRUE)
  plot.cell.density(data.cell=data.model1[[i]][data.model1[[i]]$Coh=="cohE" & data.model1[[i]]$Angle=="angleH" & data.model1[[i]]$S=="s2",],C="r2",xlim=c(1,6),show.mean = TRUE)
  plot.cell.density(data.cell=data.model1[[i]][data.model1[[i]]$Coh=="cohH" & data.model1[[i]]$Angle=="angleE" & data.model1[[i]]$S=="s2",],C="r2",xlim=c(1,6),show.mean = TRUE)
  plot.cell.density(data.cell=data.model1[[i]][data.model1[[i]]$Coh=="cohH" & data.model1[[i]]$Angle=="angleH" & data.model1[[i]]$S=="s2",],C="r2",xlim=c(1,6),show.mean = TRUE)
  
}  

# Take a look at parameters
ps <- round( attr(raw.data1, "parameters"), 2)
round(apply(ps,2,mean),2)


# accumulation rate  changes with manipulations and correct responses
p.map <- list(A="1",B= "1",t0="1",mean_v=c("Coh","Angle","M"),sd_v="M",st0="1")

model2 <- model.dmc(p.map,match.map = match.map,constants=const,
                    responses=responses,factors=factors)

p.vector  <- c(A=0.5, B = 3,  mean_v.cohE.angleE.true = 3, mean_v.cohE.angleH.true = 2.5, mean_v.cohH.angleE.true = 2.5, mean_v.cohH.angleH.true = 2,
               mean_v.cohE.angleE.false = 0, mean_v.cohE.angleH.false = 0, mean_v.cohH.angleE.false = 0, mean_v.cohH.angleH.false = 0,t0=.2,
               sd_v.true=0.5)

check.p.vector(p.vector,model2)
print.cell.p(p.vector,model1)


# Population distribution, Coherence and Angle effect on accumulation rate
pop.mean <- c(A=0.5, B = 3,  mean_v.cohE.angleE.true = 3, mean_v.cohE.angleH.true = 2.5, mean_v.cohH.angleE.true = 2.5, mean_v.cohH.angleH.true = 2,
              mean_v.cohE.angleE.false = 0, mean_v.cohE.angleH.false = 0, mean_v.cohH.angleE.false = 0, mean_v.cohH.angleH.false = 0,t0=.2,
              sd_v.true=0.5)
pop.scale <-c(A=.1, B = .1,  mean_v.cohE.angleE.true = .2, mean_v.cohE.angleH.true = .2, mean_v.cohH.angleE.true = .2, mean_v.cohH.angleH.true = .2,
  mean_v.cohE.angleE.false = .2, mean_v.cohE.angleH.false = .2, mean_v.cohH.angleE.false = .2, mean_v.cohH.angleH.false = .2, sd_v.true=0.1, t0=.05
  )

check.p.vector(pop.mean,model2)
check.p.vector(pop.scale,model2)



pop.prior <- prior.p.dmc(
  dists = rep("tnorm",12),
  p1=pop.mean,p2=pop.scale,
  lower=c(0,0,NA,NA,NA,NA,NA,NA,NA,NA,0,.1),upper=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1)
)

##  Check population distributions
par(mfcol=c(2,6)); for (i in names(pop.prior)) plot.prior(i,pop.prior)


# Simulate some data
raw.data2 <- h.simulate.dmc(model2, p.prior = pop.prior, n = 250, ns = 40)
data.model2 <- data.model.dmc(raw.data2, model2)

# Take a look at the first  subject's data
par(mfcol=c(4,2)) 
for (i in 1) { # First column=response left, Second column = response right. Rows = EasyCoh-EasyAngle, EasyCoh-HardAngle,HardCoh-EasyAngle, hard-hard
  plot.cell.density(data.cell=data.model2[[i]][data.model2[[i]]$Coh=="cohE" & data.model2[[i]]$Angle=="angleE" & data.model2[[i]]$S=="s1",],C="r1",xlim=c(1,6),show.mean = TRUE)
  plot.cell.density(data.cell=data.model2[[i]][data.model2[[i]]$Coh=="cohE" & data.model2[[i]]$Angle=="angleH" & data.model2[[i]]$S=="s1",],C="r1",xlim=c(1,6),show.mean = TRUE)
  plot.cell.density(data.cell=data.model2[[i]][data.model2[[i]]$Coh=="cohH" & data.model2[[i]]$Angle=="angleE" & data.model2[[i]]$S=="s1",],C="r1",xlim=c(1,6),show.mean = TRUE)
  plot.cell.density(data.cell=data.model2[[i]][data.model2[[i]]$Coh=="cohH" & data.model2[[i]]$Angle=="angleH" & data.model2[[i]]$S=="s1",],C="r1",xlim=c(1,6),show.mean = TRUE)
  plot.cell.density(data.cell=data.model2[[i]][data.model2[[i]]$Coh=="cohE" & data.model2[[i]]$Angle=="angleE" & data.model2[[i]]$S=="s2",],C="r2",xlim=c(1,6),show.mean = TRUE)
  plot.cell.density(data.cell=data.model2[[i]][data.model2[[i]]$Coh=="cohE" & data.model2[[i]]$Angle=="angleH" & data.model2[[i]]$S=="s2",],C="r2",xlim=c(1,6),show.mean = TRUE)
  plot.cell.density(data.cell=data.model2[[i]][data.model2[[i]]$Coh=="cohH" & data.model2[[i]]$Angle=="angleE" & data.model2[[i]]$S=="s2",],C="r2",xlim=c(1,6),show.mean = TRUE)
  plot.cell.density(data.cell=data.model2[[i]][data.model2[[i]]$Coh=="cohH" & data.model2[[i]]$Angle=="angleH" & data.model2[[i]]$S=="s2",],C="r2",xlim=c(1,6),show.mean = TRUE)
  
}  

# Take a look at parameters
ps <- round( attr(raw.data2, "parameters"), 2)
round(apply(ps,2,mean),2)



