set.seed(123)#for reproducibility
wd <- here::here();

setwd("/Users/ale/ownCloud/3-R/DMC_190819/")
source("dmc/dmc.R")
load_model (dir_name="LBA",model_name="lba_B.R")

setwd(wd)

pacman::p_load('tidyverse','magrittr')

# Define factorial design 2key x 2coh x 2angles
# s1/s2 left and right motion direction
# cohE/cohH easy/hard coherence
# angleE/angleH easy/hard decision boundary
factors <- list(S=c("s1","s2"),Coh=c("cohE","cohH"),Angle=c("angleE","angleH"))
responses <- c("r1","r2") # keys o and p

match.map <- list(M=list(s1="r1",s2="r2"))

const <- c(sd_v.false=1,st0=0)

# boundary  changes with manipulations, accum-rate faster for correct
p.map <- list(A="1",B=c("Coh","Angle"),t0="1",mean_v=c("M"),sd_v="M",st0="1")

model1 <- model.dmc(p.map,match.map = match.map,constants=const,
                   responses=responses,factors=factors)

p.vector  <- c(A=0.5, B.cohE.angleE = 3, B.cohE.angleH = 4, B.cohH.angleE = 5, B.cohH.angleH = 6,t0=.2,
               mean_v.true=2,mean_v.false=-1,sd_v.true=0.5)

check.p.vector(p.vector,model1)
print.cell.p(p.vector,model1)

data <- simulate.dmc(p.vector,model1,n=1)

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



