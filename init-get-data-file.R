####### K. Garner, 2020
####### wrapper code for using wrangling and plotting functions for the variability/decision-making data
####### works for data output by: https://github.com/kel-github/variability-decision-making
####### this code will clean/sort the data, compute transition matrices by trial, measure the correlation between transition matrices/trials and plot
####### and compute the KL divergence for the final blocks of each session
rm(list=ls())

##### load packages, set wd, get custom functions
library(tidyverse)
library(cowplot)
library(wesanderson)
theme_set(theme_cowplot())
setwd("~/Dropbox/QBI/variability-learning/var-dec-data") # set the working directory to the *copied* data location
source("../analysis/variability_exp_behav_data_wrangling_functions.R")
source("../analysis/variability_exp_behav_data_plotting_functions.R")

#### load data
s <- c(1,2,3,4,5,6,9)
subs <- rep(s, each=2)
ses <- rep(c(1,2), times=length(s))
out.dat <- mapply(get_data, subs, ses, SIMPLIFY=FALSE)
out.dat <- do.call(rbind, out.dat)

### now data is loaded, set factors
out.dat$sub    = as.factor(out.dat$sub)
out.dat$sess[out.dat$sess < 1] = 2 # for a bug with sub 1's sess 2 output, need to check if this generalises across to other subs
out.dat$sess   = as.factor(out.dat$sess)
out.dat$cond   = as.factor(out.dat$cond)
levels(out.dat$cond) <- c("A", "B")
out.dat$block  = as.factor(out.dat$block)
out.dat$door_p = as.factor(out.dat$door_p)
out.dat$door   = as.factor(out.dat$door)

#######################################################################################################################################    
########################### PRODUCE INITIAL PLOT - PERSON vs WORLD) ################################################################### 
# Summarise long form
# get the probability of being on the door, given all the fixations of that participant
sum.dat.freq = out.dat %>% group_by(sub, sess, cond, block) %>%
                           count(door) %>%
                           mutate(look_freq = n/sum(n))
# get the probabililty for each door
sum.dat.door_p = out.dat %>% group_by(sub, sess, cond, door, block) %>%
                                      summarise(p = door_p[1])
all.sum.dat = inner_join(sum.dat.freq, sum.dat.door_p, by=c("sub", "sess", "cond", "block", "door"))
all.sum.dat$num_p <- NA
all.sum.dat$num_p[all.sum.dat$p == "0"] = 0
all.sum.dat$num_p[all.sum.dat$p == "0.25"] = 0.25
all.sum.dat$num_p[all.sum.dat$p == "0.083"] = 0.083

##### plot the frequency for any subject. to see their average looks over condition x session
##### use this line of code for data exploration
p = get.sub.level.plot(all.sum.dat[all.sum.dat$sub == "9",])

##########################################################################################################################################    
######################################################  COMPUTE TRANSITION MATRICES ######################################################  
########################################################################################################################################## 
se = 2
c = 2
subs = rep(s, each = se*c) 
sess = rep(c(1:2), each=se, times=length(s))
cond = rep(c("A", "B"), times=se*length(s))
dist.dat <- do.call(rbind, mapply(function(x,y,z) get_trans_sim_by_cond(data=out.dat[out.dat$sub == x & out.dat$sess == y & out.dat$cond == z, ], ndoors=16), subs, sess, cond, SIMPLIFY=FALSE))
plot.distances(dist.dat)
  
#######################################################################################################################################    
######################################################  COMPUTE KL VARIABLES ######################################################  
# KL things
# 1. functions to determine if a Dirichlet model approximates the experienced trials
#### get the door frequency data in desired format

##### THE BELOW NEEDS NEATENING UP!!!!
subjs = s
tsubs = length(subjs)
session = c(1, 2)
tconds = 2
tblocks = 8
conds = rep(c("A", "B"), each=tblocks, times = length(subs)*length(session))
blocks = rep(seq(1,tblocks,1), times = tconds*length(subs)*length(session))
ndoors = 16
subs = rep(subjs, each = length(session)*tconds*tblocks)
sessions = rep(session, times=tsubs*tconds*tblocks)
tmp = mapply(add.missing.doors, cond = conds, block = blocks, sub = subs, sess = sessions, MoreArgs=list(data=all.sum.dat, ndoors=ndoors), SIMPLIFY=FALSE)
tmp = do.call(rbind, tmp)
all.sum.dat = as.data.frame(all.sum.dat)
all.sum.dat = rbind(all.sum.dat, tmp)

subs = rep(s, each=length(session))
sessions = rep(session, times=tsubs)
sum.dat.door_p = mapply(get.summary.n.for.door, sub_num = subs, sess = sessions, SIMPLIFY=FALSE)
sum.dat.door_p = do.call(rbind, sum.dat.door_p)
sum.dat.for.kl = inner_join(all.sum.dat, sum.dat.door_p, by=c("sub", "sess", "cond", "block", "door"))
sum.dat.for.kl$door = as.factor(sum.dat.for.kl$door)
sum.dat.for.kl$sess = as.factor(sum.dat.for.kl$sess)
rm(tmp)

#### now add priors
ndoors = 16
doors = rep(c(1:ndoors), times = tconds*tsubs*length(session)) # doors * conditions
conds.to.calc = rep(c("A", "B"), each = ndoors, times=tsubs*length(session))
subs = rep(subjs, each=ndoors*tconds*length(session))
sess = rep(session, each=ndoors*tconds, times=length(session))
tblocks = 8

### ESTABLISH WORLD AND SUBJECT DIRICHLET PARAMS, GIVEN EXPERIENCE
priors.over.world.and.sub = mapply(get.priors, cond=conds.to.calc, door=doors, sub=subs, sess=sess, MoreArgs = list(data=sum.dat.for.kl, tblocks=8), SIMPLIFY = FALSE)
priors.over.world.and.sub = do.call(rbind, priors.over.world.and.sub)
priors.over.world.and.sub$block = as.factor(priors.over.world.and.sub$block)
priors.over.world.and.sub$door = as.factor(priors.over.world.and.sub$door)
priors.over.world.and.sub$sub = as.factor(priors.over.world.and.sub$sub)
priors.over.world.and.sub$sess = as.factor(priors.over.world.and.sub$sess)
sum.dat.for.kl = inner_join(sum.dat.for.kl, priors.over.world.and.sub, by=c("sub", "sess", "cond", "door", "block"))
sum.dat.for.kl$door = as.factor(sum.dat.for.kl$door)


# NOW COMPUTE KL DIVERGENCE BETWEEN BLOCKS FOR EACH PARTICIPANT AND SESSION
tblocks = 8
subs = rep(subjs, each = tconds*length(session)*tblocks)
conds = rep(c("A", "B"), times = tsubs*length(session), each=tblocks)
sess = rep(c("1", "2"), each=tblocks*tconds, times=length(s))
blocks = rep(c(1:8), times = length(subs))
kl.over.blocks = mapply(get.kl.div.blocks, cond = conds, block = blocks, sub = subs, sess = sess, MoreArgs=list(sum.dat.for.kl), SIMPLIFY=FALSE)
kl.over.blocks = do.call(rbind, kl.over.blocks)


tmp.dat = rbind(kl.over.blocks[kl.over.blocks$sess == 1 & kl.over.blocks$block == 1,],
                kl.over.blocks[kl.over.blocks$sess == 2 & kl.over.blocks$block == 8,])
tmp.dat %>% ggplot(aes(sess, D, fill=sess)) +
  geom_boxplot() +
  geom_line(aes(group=sub), position = position_dodge(0.2), alpha=0.5) +
  geom_point(aes(fill=sess,group=sub), position = position_dodge(0.2), alpha=0.5) +
  facet_wrap(~cond)


# this function computes p(x; a) 
get.p.given.a.param <- function(a, x){
  # COME BACK TO THIS
  # this function computes the dirichlet density function for observations
  # x, given alpha parameters a
  # implemented with help from: http://bariskurt.com/kullback-leibler-divergence-between-two-dirichlet-and-beta-distributions/
  # p(x; a) = gamma(sum(alphas(1:K))) / productk=1:K (gamma(a_k)) then multiply with the product of x_k^a_k-1  
  # 
  p_x_g_a = gamma(sum(a))/prod(gamma(a)) * prod(x^(a-1))
  p_x_g_a
}

#
























# now write function to compute KL Divergence between trial t and the model of the world
get.world.KL.over.trials <- function(data){
  
  # how to compute KL with zero entries
  # https://mathoverflow.net/questions/72668/how-to-compute-kl-divergence-when-pmf-contains-0s
  
  # get p of each door in each trial 
  tmp      = data %>% group_by(cond, block) %>% count(door) %>% mutate(p=n/sum(n))
  tmp.door = data %>% group_by(cond, block) %>% summarise(tgt = last(tgt_door))
  # for each trial in tmp, add missing doors, and add zeros, as no glances at them

  missing.dat = lapply(unique(tmp$t), add.missing.doors, tmp=tmp)   
  missing.dat = do.call(rbind, missing.dat)
  tmp = data.frame(tmp)
  tmp = rbind(tmp, missing.dat)
  tmp = inner_join(tmp, tmp.door, by=c("cond", "t"))
  
  #### now that data frame is set up, go through trials and update a 
  # dirichlet distribution for each context, for each block of trials (n=10)
  # save the parameters (alpha) for the next block
  
  # for each step, calculate the KL between the recursive prior for that context, and the 
  # prior for the other context - 
  # creates a new dataframe with the condition, and two KL values for each trial
  get.kl.over.trials <- function(tmp){
    
    cert_KLs   = c()
    uncert_KLs = c()
    cond_vec   = c()
    t_vec      = c()
    for (i in 1:length(unique(tmp$t))){
      

      if (i == 1){
        # use a to denote the prior counts. This will be updated in each
        # iteration after comparison to subject pattern below
        cert_a       = rep(1, times = 16)
        uncert_a     = cert_a
        sub_cert_a   = cert_a
        sub_uncert_a = cert_a
      }
      
      # now get the vector for the current trial
      get.trial = tmp[tmp$t == i, ]
      get.trial = get.trial[order(get.trial$door), ]
      
      # assign other variables
      t_vec[i] = get.trial$t[1]
      cond_vec[i] = get.trial$cond[1]
      # now get KLs for cert and uncert
      cert_KLs[i] = entropy::KL.Dirichlet(get.trial$n, cert_a, sub_cert_a, cert_a)
      uncert_KLs[i] = entropy::KL.Dirichlet(get.trial$n, uncert_a, sub_uncert_a, uncert_a)
      
      # now update the prior
      if (get.trial$cond[1] == "cert"){
        cert_a[get.trial$tgt[1]] = cert_a[get.trial$tgt[1]] + 1
        sub_cert_a = sub_cert_a + get.trial$n
      } else if (get.trial$cond[1] == "uncert") {
        uncert_a[get.trial$tgt[1]] = uncert_a[get.trial$tgt[1]] + 1
        sub_uncert_a = sub_uncert_a + get.trial$n
      }
     
    }
    
    kl.over.trials = data.frame(t = c(t_vec, t_vec), 
                                cond = c(cond_vec, cond_vec), 
                                KL_idx = rep(c(1,2), each = length(t_vec)),
                                KL = c(cert_KLs, uncert_KLs))
    
  }
  
  
  kl.dat = get.kl.over.trials(tmp)
  kl.dat$cond = as.factor(kl.dat$cond)
  levels(kl.dat$cond) = c("cert", "uncert")
  kl.dat$KL_idx = as.factor(kl.dat$KL_idx)
  levels(kl.dat$KL_idx) = c("cert", "uncert")
  kl.dat
}
  
kl.dat = get.world.KL.over.trials(out.dat)  
  

# now have data (for one sub), plot
plot.KL.dat <- function(kl.dat){
  
  
  p = kl.dat %>%
      ggplot(aes(t, KL, group=KL_idx, colour=KL_idx)) + 
      geom_line() +
      facet_wrap(~cond) +
      theme_classic() 
  
  p
}

plot.KL.dat(kl.dat)



##### long to wide
sum.dat.all.wide.freq = reshape2::dcast(all.sum.dat, sub ~ cond + block + door, value.var = c("look_freq")) 
write.csv(sum.dat.all.wide.freq, file="pilot1-subs21-22-loc-freq.csv")

sum.dat.all.wide.p = reshape2::dcast(all.sum.dat, sub ~ cond + block + door, value.var = c("p")) 
write.csv(sum.dat.all.wide.p, file="pilot1-subs21-22-loc-p.csv")                      

sum.dat.all.wide.rt = reshape2::dcast(all.sum.dat, sub ~ cond + block + door, value.var = c("rt")) 
write.csv(sum.dat.all.wide.rt, file="pilot1-subs21-22-loc-rt.csv") 

