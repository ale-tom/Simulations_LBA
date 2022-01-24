set.seed(123)#for reproducibility

# Set your working directory to the DMC folder  to load the model
wd <- here::here();

#setwd("C:/Users/at07/ownCloud/3-R/DMC_190819/")
setwd("/imaging/rowe/users/at07/Matlab/Projects/CBU2020/LC/Model_fitting/")
source("dmc/dmc.R")
load_model (dir_name="LBA",model_name="lba_B.R")

# Go back to the project's folder
setwd(wd)

pacman::p_load('tidyverse','magrittr','slurmR')

# Setting required parameters
# remove old jobs
unlink(here::here('test', 'test-DorianSlurm'),recursive = TRUE)
# set up folder For slurm jobs
opts_slurmR$set_tmp_path(here::here("/test/"))
opts_slurmR$set_job_name("test-DorianSlurm")

# Load and preprocess the data -------------------------------------------------------------------
# this probably needs to be changed 

data <- read.csv(here::here('test',"Data/AccrateRT.csv"),na.strings=c("","NA")) %>% as.data.frame()

#data <- simulate.dmc(p.vector,model1,n=1) 

# get indices of  correct  trials and remove implausible RTs
data <- data %>% mutate(Correct = ifelse((S=="s1" & R == "r1") | (S=="s2" & R == "r2"),"Correct","Wrong"))  %>%
  filter(!is.na(R),.preserve = TRUE) %>%
  mutate(Plausibility = case_when(
    Correct == "Correct" & RT < 0.25 ~ "Too quick",
    Correct == "Correct" & RT > 4.5 ~ "Too slow",
    TRUE ~ "Good")) %>%
  filter(Plausibility == "Good",.preserve = TRUE)

# remove outlier RTs (mean ± 2.5sd) on  a subject-by-subject basis
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

#note i is the index for the model variant (total number of variants is 34)
LBAFit <- function(i,data) {
  # Set your working directory to the DMC folder  to load the model
  wd <- here::here();
  
  #setwd("C:/Users/at07/ownCloud/3-R/DMC_190819/")
  setwd("/imaging/rowe/users/at07/Matlab/Projects/CBU2020/LC/Model_fitting/")
  source("dmc/dmc.R")
  load_model (dir_name="LBA",model_name="lba_B.R")
  
  # Go back to the project's foldersavename =   sprintf('hsampling_LBA%i_LBA.RData',i)
  setwd(wd)
  
  pacman::p_load('tidyverse','magrittr')
  
  # prepare design matrix -------------------------------------------------------------------------
  #34 models
  levs <- list("1","Coh","Angle",c("Coh","Angle"))
  levs2 <- list(c(0,0),c(1,0),c(0,1),c(1,1))
  design_matrix <- expand.grid(c(1:3),c(1:3),c(1:3))
  design_matrix <- design_matrix[-1,] %>% rbind(expand.grid(c(1,4),c(1,4),c(1,4))) %>%  as.data.frame()
  colnames(design_matrix) <- c('B','mean_v','t0')
  
  
  
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
  # e.g. boundary  changes with manipulations, accum-rate faster for correct
  #p.map <- list(A="1",B=c("Coh","Angle"),t0="1",mean_v=c("M"),sd_v="M",st0="1")
  pB <-levs[design_matrix$B[i]] %>% unlist()
  pV <-levs[design_matrix$mean_v[i]] %>% unlist()
  pV <- if(pV=='1') 'M' else c(pV,'M')
  pT0 <- levs[design_matrix$t0[i]] %>% unlist()
  
  p.map <- list(A="1",B=pB,t0=pT0,mean_v=pV,sd_v="M",st0="1")
  
  model1 <- model.dmc(p.map,match.map = match.map,constants=const,
                      responses=responses,factors=factors)
  
  
  # Define the prior distributions -----------------------------------------------------------------
  
  vec_p <- function(param,factor1 = 1,factor2 = 1, value,plower,pupper)
  {
    if(!factor1 & !factor2){
      lbls <- param  
    }
    else if(factor1 & !factor2) {
      lbls <- do.call(paste, c(expand.grid(param, c("cohE","cohH")), list(sep='.')))}
    
    else if (!factor1 & factor2) {
      lbls <- do.call(paste, c(expand.grid(param,  c("angleE","angleH")), list(sep='.')))
    }
    else{
      lbls<- do.call(paste, c(expand.grid(param, c("cohE","cohH"), c("angleE","angleH")), list(sep='.')))
      
    }
    
    parv <- rep(value,length(lbls))
    names(parv) <- lbls
    plow <- rep(plower,length(parv))
    phigh <- rep(pupper,length(parv))
    
    vec <- list("parv"=parv,"lower"=plow,"upper"=phigh)
    return(vec)
    
  }
  
  vec_pV <- function(param,factor1 = 1,factor2 = 1, value_true,value_false,plower,pupper)
  {  
    true_false <- c("true","false")
    if(!factor1 & !factor2){
      lbls <- do.call(paste, c(expand.grid(param,true_false),list(sep ='.')))  
    }
    else if(factor1 & !factor2) {
      lbls <- do.call(paste, c(expand.grid(param, c("cohE","cohH"),true_false), list(sep='.')))
    }
    else if (!factor1 & factor2) {
      lbls <- do.call(paste, c(expand.grid(param,  c("angleE","angleH"),true_false), list(sep='.')))
    }
    else{
      lbls<- do.call(paste, c(expand.grid(param, c("cohE","cohH"), c("angleE","angleH"),true_false), list(sep='.')))
      
    }
    
    parv <- c(rep(value_true,length(lbls)/2),rep(value_false,length(lbls)/2))
    names(parv) <- lbls
    
    
    plow <- rep(plower,length(parv))
    phigh <- rep(pupper,length(parv))
    
    vec <- list("parv"=parv,"lower"=plow,"upper"=phigh)
    return(vec)
    
  }
  
  comp_vec <- function(vA,vB,vV,vsd,vt0) {
    pars<-c(vA$parv,vB$parv,vV$parv,vsd$parv,vt0$parv)
    upper<-c(vA$upper,vB$upper,vV$upper,vsd$upper,vt0$upper)
    lower<-c(vA$lower,vB$lower,vV$lower,vsd$lower,vt0$lower)
    comp <- list("pars" = pars,"upper"=upper,"lower"=lower)
    return(comp)
  }
  
  B <-levs2[design_matrix$B[i]] %>% unlist()
  V <-levs2[design_matrix$mean_v[i]] %>% unlist()
  T0 <- levs2[design_matrix$t0[i]] %>% unlist()
  
  pop_prior <- comp_vec(
    vec_p('A',0,0, .4,0,NA),
    vec_p('B',B[1],B[2], 4,0,NA),
    vec_pV('mean_v',V[1],V[2], 1,-1,NA,NA),
    vec_p('sd_v.true',0,0, .3,0,NA),
    vec_p('t0',T0[1],T0[2], .4,.1,1))
  pop_scale <- comp_vec(
    vec_p('A',0,0, .1,0,NA),
    vec_p('B',B[1],B[2], .2,0,NA),
    vec_pV('mean_v',V[1],V[2], .2,.2,NA,NA),
    vec_p('sd_v.true',0,0, .1,0,NA),
    vec_p('t0',T0[1],T0[2], .05,.1,1))
  
  
  
  
  #double check that parameters match with model
  pop.mean <- pop_prior$pars
  check.p.vector(pop.mean,model1)
  pop.scale <- pop_scale$pars
  check.p.vector(pop.scale,model1)
  
  nparams = length(pop.mean)
  paramnames = names(pop.mean)
  
  p.prior <- prior.p.dmc(dists = rep("tnorm",nparams),
                         p1=pop.mean,p2=pop.scale,
                         lower = pop_prior$lower,
                         upper = pop_prior$upper)
  
  
  
  
  # bind the data to the model object
  data.model <- data.model.dmc(data, model1)
  
 
  
  
  
  #Fit the model non-hierarchically (i.e.treating subjects as fixed effects ---------------------
  #number of cores available
  mycores <- 10

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
  savename =   sprintf('sampling_LBA%i_LBA.RData',i)
  save(sampling_LBA, file = here::here('test',savename))

  # 
  
  # # Fit the model hierarchically (treating subjects as a random effect) -----------------------------
  # 
  # # We will use the results of the non-hierarchical modelling to generate starting points for the
  # # hierarchical modelling
  # groupstart <- make.hstart(sampling_LBA) # group-level parameters
  # subjectstart <- make.theta1(sampling_LBA) # subject-level parameters
  # 
  # # If we perform hierarchical modelling, we assume that each parameter has a group-level
  # # distribution with some location (e.g. mean) and scale (e.g. SD). These group-level distributions
  # # then serve as priors for our parameters at the level of individual subjects.
  # # Thus, we need to define prior distributions for both the group-level mean and the group-level SD
  # # of each parameter. These are also known as hyperpriors.
  # 
  # # For the prior distributions for the location (i.e. mean) of our group-level distributions, we
  # # will simply recycle the priors we used for the fixed-effects modelling
  # mu.prior <- p.prior
  # 
  # # What's new with hierarchical modelling is we also have to define prior distributions for the
  # # scale (i.e. SD) of our group-level distributions
  # sigma.prior <- prior.p.dmc(
  #   # set the type of distribution: Gamma distribution for all parameters
  #   dists = rep(x = "gamma", times = nparams),
  #   # set the shape parameter (p1) of the distribution: Shape = 1, hence exponential distribution
  #   p1 = c(
  #     setNames(
  #       object = rep(x = 1, times = nparams),
  #       nm = paramnames
  #     )
  #   ),
  #   # set the scale parameter (p2) of the distribution
  #   p2 = c(
  #     setNames(
  #       object = rep(x = 1, times = nparams),
  #       nm = paramnames
  #     )
  #   )
  # )
  # 
  # # Lastly, we combine these into a hyper-prior list
  # pp.prior <- list(mu.prior, sigma.prior)
  # 
  # # Set up our model fitting object
  # hsampling_LBA <- h.samples.dmc(
  #   nmc = 400, # number of samples AFTER thinning
  #   pp.prior = pp.prior, # hyper-priors
  #   p.prior = mu.prior, # prior for each subject
  #   hstart.prior = groupstart,
  #   theta1 = subjectstart,
  #   thin = 10, # thinning to deal with autocorrelation
  #   data = data.model # note: Default number of chains is 3 * nparams
  # )
  # 
  # # Now we run the models using some automated functions from the DMC toolbox
  # hsampling_LBA <- h.run.unstuck.dmc(
  #   samples = hsampling_LBA, cores = mycores, report = 10,
  #   p.migrate = 0.05, h.p.migrate = 0.05)
  # hsampling_LAB <- h.run.converge.dmc(
  #   samples = h.samples.dmc(samples = hsampling_LBA, nmc = 120, thin = 25),
  #   nmc = 40, cores = mycores, report = 10, verbose = TRUE)
  # 
  # # after the model has converged, get some "clean" samples to be used for parameter estimation
  # hsampling_LBA <- h.run.dmc(
  #   h.samples.dmc(samples = hsampling_LBA, nmc = 500),
  #   report = 1, cores = mycores)
  # hsampling_CamCAN <- h.run.dmc(
  #   h.samples.dmc(samples = hsampling_LBA, nmc = 500),
  #   report = 1, cores = mycores)
  # 
  # savename =   sprintf('hsampling_LBA%i_LBA.RData',i)
  # save(hsampling_LBA, file = here::here('test',savename))
  # 
  
  
  
}

#submit jobs to cluster
# here I test only the first 5 variants 
job <- Slurm_lapply(c(1:5), LBAFit, data = data, njobs=5, mc.cores=10, plan = "collect")
#check jobs exit status 
status(job)
