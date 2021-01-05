data <- function(sub_num){
  # this function reads in the behavioural and trial logs, and combines the two to provide a trial by trial dataframe
  # that contains the pre-defined door probabilities and participant response to them
  # input a numeric subject value (not zero padded)
  trials <- read.table(sprintf('sub-0%d/beh/sub-0%d_ses-1_task-iforage-v1_trls.tsv', sub_num, sub_num), header = TRUE)
  resps <- read.table(sprintf('sub-0%d/beh/sub-0%d_ses-1_task-iforage-v1_beh.tsv', sub_num, sub_num), header = TRUE)
  # remove all practice trials
  trials = trials %>% filter(t != 999) 
    
  resps  = resps  %>% filter(cond != 3 & open_d != 0) %>%
           filter(lead(door,1) != door)
  resps$t = resps$t - 5 # remove practice trials

  ### now assign condition blocks to trials
  trials$block = NA
  trials_per_block = 10 
  total_blocks = length(trials$cond[trials$cond==1])/trials_per_block
  trials$block[trials$cond == 1] = rep(c(1:total_blocks), each = length(trials$cond[trials$cond==1])/total_blocks)
  trials$block[trials$cond == 2] = rep(c(1:total_blocks), each = length(trials$cond[trials$cond==2])/total_blocks)
  ### and probability of location
  names(trials)[6] = "tgt_prob" # rename to make life easier
  
  ### now combine the two by sub and t
  sub.dat = inner_join(resps, trials, by=c("sub", "t", "sess", "cond"))
  sub.dat
}


add.missing.doors <- function(data, cond, block){
  # this function takes the summary data 'all.sum.dat', generated in init-get-data-file, and 
  # adds the missing doors from each trial, i.e. ones participants did not look at
  miss.doors = setdiff(levels(data$door[data$cond == cond]), unique(data$door[data$cond == cond & data$block == block]))
  if (sum(lengths(miss.doors)) > 0){
    sub = rep(as.character(data$sub[1]), times = length(miss.doors))
    this.cond = rep(cond, times = length(miss.doors))
    this.block = rep(block, times = length(miss.doors))
    door = miss.doors
    n = rep(0, times = length(miss.doors))
    look_freq = rep(0, times = length(miss.doors))
    num_p = rep(0, times = length(miss.doors))
    for (x in miss.doors){
      
      num_p[door == x] = data$num_p[data$cond == cond & data$door == x][1]
    }
    p = as.character(num_p)
    tmp = data.frame(sub = sub, cond = this.cond, block = this.block, door = door, n = n, look_freq = look_freq, num_p = num_p, p = p)
    tmp
    }
}

get.summary.n.for.door <- function(sub_num){
  # this function takes the sub number (sub_num), opens up the trial log for that participant, and
  # counts the frequencies with which each target appeared behind any given door, the output is a dataframe
  # that can be combined with the dataframe of subject saccade/door frequencies
  tmp.dat =  read.table(sprintf('pilot_b/sub-0%d_ses-1_task-iforage_trls.tsv', sub_num), header = TRUE)
  tmp.dat = tmp.dat[tmp.dat$t != 999, ]
  trials_per_block = 10
  total_blocks = length(tmp.dat$cond[tmp.dat$cond==1])/trials_per_block
  tmp.dat$block[tmp.dat$cond == 1] = rep(c(1:total_blocks), each = length(tmp.dat$cond[tmp.dat$cond==1])/total_blocks)
  tmp.dat$block[tmp.dat$cond == 2] = rep(c(1:total_blocks), each = length(tmp.dat$cond[tmp.dat$cond==2])/total_blocks)
  
  sum.dat.door_p = tmp.dat %>% group_by(sub, cond, block) %>%
    count(loc) %>%
    mutate(door_freq = n/sum(n))
  names(sum.dat.door_p)[4] = "door"
  sum.dat.door_p$sub = as.factor(sum.dat.door_p$sub)
  sum.dat.door_p$cond = as.factor(sum.dat.door_p$cond)
  levels(sum.dat.door_p$cond) = c("cert", "uncert")
  sum.dat.door_p$block = as.factor(sum.dat.door_p$block)
  sum.dat.door_p$door = as.factor(sum.dat.door_p$door)
  
  # add the other doors to each trial
  sum.dat.door_p = data.frame(sum.dat.door_p) # for rbinding with tmp
  door.test = as.factor(seq(1,16,1))
  for (i in levels(sum.dat.door_p$cond)){
    for (j in levels(sum.dat.door_p$block)){
      miss.doors = setdiff(door.test, sum.dat.door_p$door[sum.dat.door_p$cond == i & sum.dat.door_p$block == j])
      if (sum(lengths(miss.doors)) > 0){
        tmp  <-  data.frame(sub   = rep(as.character(sum.dat.door_p$sub[1]), length(miss.doors)),
                            cond  = rep(i,length(miss.doors)),
                            block = rep(j,length(miss.doors)),
                            door  = miss.doors,
                            n     = rep(0, length(miss.doors)),
                            door_freq = rep(0, length(miss.doors)))
        
        sum.dat.door_p = rbind(sum.dat.door_p, tmp)
      }
    }
  }
  sum.dat.door_p
}


get.priors <- function(data, sub, cond, door, tblocks){
  # taking the sum.dat.for.kl dataframe, this function adds the priors for the associated Dirichlet distribution
  world_prior=c()
  world_prior[1] = 1/16
  sub_prior=c()
  sub_prior[1] = 1/16
  block_count = c()
  block_count[1] = 1
  for (i in seq(2,tblocks,1)){
    world_prior[i] = sum(world_prior[i-1], data$door_freq[data$door == door & data$block== i-1 & data$cond == cond & data$sub == sub])
    sub_prior[i] = sum(sub_prior[i-1],  data$look_freq[data$door == door & data$block== i-1 & data$cond == cond & data$sub == sub])
    block_count[i] = i
  }
  out=data.frame(world=world_prior, 
                 subs_world=sub_prior, 
                 door = rep(door, times = length(world_prior)), 
                 block = block_count,
                 sub = rep(sub, length(world_prior)),
                 cond = rep(cond, length(world_prior)))
  out
}


get.kl.div.blocks <- function(data, block, cond){
  ### using the data frame 'sub.dat.for.kl', compute, for each block and condition, 
  ### the kl divergence between the data and the world probabilites
  
  current.dat = data[data$cond == cond & data$block == block, ]
  current.dat = current.dat[order(as.numeric(current.dat$door)),]
  # the below function (KL.Dirichlet) - uses these functions:
  # which appear to take the frequencies of the counts, and compute their probability
  # without computing the pdf as defined by the dirichlet distribution - so not using
  # https://github.com/cran/entropy/blob/master/R/KL.plugin.R
  # https://github.com/cran/entropy/blob/master/R/entropy.Dirichlet.R
#  kl.out      = with(current.dat, entropy::KL.Dirichlet(n.y, n.x, world, subs_world))
  
  # but see manual calculation for comparison, based on:
  # 
  
  world.alpha = with(current.dat, n.y + world) 
  sub.alpha = with(current.dat, n.x + subs_world)
  D = lgamma(sum(world.alpha)) - lgamma(sum(sub.alpha)) - sum(lgamma(world.alpha)) + sum(lgamma(sub.alpha)) + sum((world.alpha - sub.alpha)*(digamma(world.alpha)-digamma(sum(world.alpha))))
  # and this derivation on https://github.com/cran/Compositional/blob/master/R/kl.diri.R
  world.0 = sum(world.alpha) # ease for interpreting their eq
  sub.0 = sum(sub.alpha)
  f = sum( (world.alpha - sub.alpha) * ( digamma(world.alpha) - digamma(world.0) ) ) + sum( lgamma(sub.alpha) - lgamma(world.alpha) ) + lgamma(world.0) - lgamma(sub.0)
  out = list(f, D)
  out
}


