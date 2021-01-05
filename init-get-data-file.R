rm(list=ls())
setwd("~/Dropbox/QBI/variability-learning/pilot-data")
source("variability_exp_behav_data_wrangling_functions.R")
source("variability_exp_behav_data_plotting_functions.R")
##### packages
library(dplyr)
library(ggplot2)
#### load data
subs <- c(1)
out.dat <- lapply(subs, data)
out.dat <- do.call(rbind, out.dat)

### now data is loaded, set factors
out.dat$sub    = as.factor(out.dat$sub)
out.dat$cond   = as.factor(out.dat$cond)
levels(out.dat$cond) <- c("cert", "uncert")
out.dat$block  = as.factor(out.dat$block)
out.dat$door_p = as.factor(out.dat$door_p)
out.dat$door   = as.factor(out.dat$door)

#######################################################################################################################################    
########################### PRODUCE INITIAL PLOT (LAST TWO BLOCKS - PERSON vs WORLD) ########################### 
# Summarise long form
sum.dat.freq = out.dat %>% group_by(sub, cond, block) %>%
                           count(door) %>%
                           mutate(look_freq = n/sum(n))
sum.dat.door_p = out.dat %>% group_by(sub, cond, door, block) %>%
                                      summarise(p = door_p[1])
all.sum.dat = inner_join(sum.dat.freq, sum.dat.door_p, by=c("sub", "cond", "block", "door"))
all.sum.dat$num_p[all.sum.dat$p == "0"] = 0
all.sum.dat$num_p[all.sum.dat$p == "0.25"] = 0.25
all.sum.dat$num_p[all.sum.dat$p == "0.083"] = 0.083

##### plot the frequency (ordered by common to not common)
p = get.sub.level.plot(all.sum.dat)

#######################################################################################################################################    
######################################################  COMPUTE KL VARIABLES ######################################################  
# KL things
# 1. functions to determine if a Dirichlet model approximates the experienced trials
#### get the door frequency data in desired format
conds = rep(c("cert", "uncert"), each=8)
blocks = rep(as.character(seq(1,8,1)), times = 2)
tmp = mapply(add.missing.doors, cond = conds, block = blocks, MoreArgs=list(data=all.sum.dat))
tmp = do.call(rbind, tmp)
all.sum.dat = as.data.frame(all.sum.dat)
all.sum.dat = rbind(all.sum.dat, tmp)
sum.dat.door_p = get.summary.n.for.door(1)
sum.dat.for.kl = inner_join(all.sum.dat, sum.dat.door_p, by=c("sub", "cond", "block", "door"))
sum.dat.for.kl$door = as.factor(sum.dat.for.kl$door)
rm(tmp)

#### now add priors
doors = rep(c(1:16), times = 2) # doors * conditions
conds.to.calc = rep(c("cert", "uncert"), each = 16)
tblocks = 8
priors.over.world.and.sub = mapply(get.priors, cond=conds.to.calc, door=doors, MoreArgs = list(sub = "1", data=sum.dat.for.kl, tblocks=8), SIMPLIFY = FALSE)
priors.over.world.and.sub = do.call(rbind, priors.over.world.and.sub)
priors.over.world.and.sub$block = as.factor(priors.over.world.and.sub$block)
priors.over.world.and.sub$door = as.factor(priors.over.world.and.sub$door)
sum.dat.for.kl = inner_join(sum.dat.for.kl, priors.over.world.and.sub, by=c("sub", "cond", "door", "block"))

conds = rep(c("cert", "uncert"), each = 8)
blocks = rep(as.character(seq(1,8,1)), times = 2)
kl.over.blocks = mapply(get.kl.div.blocks, cond = conds, block = blocks, MoreArgs=list(sum.dat.for.kl))

#### plot KL over blocks
kl.plot = plot.kls.from.list(kl.over.blocks)
kl.plot








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

