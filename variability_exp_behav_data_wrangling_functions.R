##########################################################################################################################################
############################ DATA WRANGLES ########################################################################################
##########################################################################################################################################

get_data <- function(sub_num, ses){
  # this function reads in the behavioural and trial logs, and combines the two to provide a trial by trial dataframe
  # that contains the pre-defined door probabilities and participant response to them
  # input a numeric subject value (not zero padded)
  trials <- read.table(sprintf('sub-0%d/beh/sub-0%d_ses-%d_task-iforage-v1_trls.tsv', sub_num, sub_num, ses), header = TRUE)
  resps <- read.table(sprintf('sub-0%d/beh/sub-0%d_ses-%d_task-iforage-v1_beh.tsv', sub_num, sub_num, ses), header = TRUE)
  # remove all practice trials
  trials = trials %>% filter(t != 999) 
  
  # remove the practice trials and the times where participants were not looking at a door,
  # then select the last of the multiple entries for each door, so that you get the total duration
  # that the door was open for (i.e. if the next entry for door ID does not match previous, and the next
  # depress duration is less than the current)
  resps  = resps  %>% filter(cond != 3 & open_d != 0) %>%
           filter(lead(door,1) != door | lead(depress_dur,1) < depress_dur )
  resps$t = resps$t - 5 # remove practice trials

  ### now assign condition blocks to trials
  trials$block = NA
  trials_per_block = 10 
  total_blocks = length(trials$cond[trials$cond==1])/trials_per_block
  trials$block[trials$cond == 1] = rep(c(1:total_blocks), each = length(trials$cond[trials$cond==1])/total_blocks)
  trials$block[trials$cond == 2] = rep(c(1:total_blocks), each = length(trials$cond[trials$cond==2])/total_blocks)
  ### and probability of location
  names(trials)[names(trials)=="prob"] = "tgt_prob" # rename to make life easier
  
  ### now combine the two by sub and t
  sub.dat = inner_join(resps, trials, by=c("sub", "t", "sess", "cond"))
  sub.dat
}


add.missing.doors <- function(data, cond, block, ndoors, sub, sess){
  # this function takes the summary data 'all.sum.dat', generated in init-get-data-file, and 
  # adds the missing doors from each trial, i.e. ones participants did not look at
  miss.doors = setdiff(as.factor(c(1:ndoors)), unique(data$door[data$sub == sub & data$sess == sess & data$cond == cond & data$block == block]))
  if (length(miss.doors) > 0){
    sub = rep(sub, times = length(miss.doors))
    sess = rep(sess, times = length(miss.doors))
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
    tmp = data.frame(sub = sub, sess = sess, cond = this.cond, block = this.block, door = door, n = n, look_freq = look_freq, num_p = num_p, p = p)
    tmp
    }
}

get.summary.n.for.door <- function(sub_num, sess){
  # this function takes the sub number (sub_num), opens up the trial log for that participant, and
  # counts the frequencies with which each target appeared behind any given door, the output is a dataframe
  # that can be combined with the dataframe of subject saccade/door frequencies
  tmp.dat =  read.table(sprintf('sub-0%d/beh/sub-0%d_ses-%d_task-iforage-v1_trls.tsv', sub_num, sub_num, sess), header = TRUE)
  tmp.dat = tmp.dat[tmp.dat$t != 999, ]
  trials_per_block = 10
  total_blocks = length(tmp.dat$cond[tmp.dat$cond==1])/trials_per_block
  tmp.dat$block[tmp.dat$cond == 1] = rep(c(1:total_blocks), each = length(tmp.dat$cond[tmp.dat$cond==1])/total_blocks)
  tmp.dat$block[tmp.dat$cond == 2] = rep(c(1:total_blocks), each = length(tmp.dat$cond[tmp.dat$cond==2])/total_blocks)
  
  sum.dat.door_p = tmp.dat %>% group_by(sub, sess, cond, block) %>%
    count(loc) %>%
    mutate(door_freq = n/sum(n))
  names(sum.dat.door_p)[names(sum.dat.door_p) == "loc"] = "door"
  sum.dat.door_p$sub = as.factor(sum.dat.door_p$sub)
  sum.dat.door_p$cond = as.factor(sum.dat.door_p$cond)
  levels(sum.dat.door_p$cond) = c("A", "B")
  sum.dat.door_p$block = as.factor(sum.dat.door_p$block)
  sum.dat.door_p$door = as.factor(sum.dat.door_p$door)
  
  # add the other doors to each trial
  sum.dat.door_p = data.frame(sum.dat.door_p) # for rbinding with tmp
  door.test = as.factor(seq(1,16,1))
  for (i in levels(sum.dat.door_p$cond)){
    for (j in levels(sum.dat.door_p$block)){
      miss.doors = setdiff(door.test, sum.dat.door_p$door[sum.dat.door_p$cond == i & sum.dat.door_p$block == j])
      if (sum(lengths(miss.doors)) > 0){
        tmp  <-  data.frame(sub   = rep(as.character(sub_num), length(miss.doors)),
                            sess  = rep(as.character(sess), length(miss.doors)),
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

##########################################################################################################################################
############################ KL BUSINESS ########################################################################################
##########################################################################################################################################

get.priors <- function(data, sub, sess, cond, door, tblocks){
  # taking the sum.dat.for.kl dataframe, this function adds the priors for the associated Dirichlet distribution
  world_prior=c()
  world_prior[1] = 1/16
  sub_prior=c()
  sub_prior[1] = 1/16
  block_count = c()
  block_count[1] = 1
  for (i in seq(2,tblocks,1)){
    world_prior[i] = sum(world_prior[i-1], data$door_freq[data$door == door & data$block== i-1 & data$cond == cond & data$sub == sub & data$sess == sess])
    sub_prior[i] = sum(sub_prior[i-1],  data$look_freq[data$door == door & data$block== i-1 & data$cond == cond & data$sub == sub & data$sess == sess])
    block_count[i] = i
  }
  out=data.frame(world=world_prior, 
                 subs_world=sub_prior, 
                 door = rep(door, times = length(world_prior)), 
                 block = block_count,
                 sub = rep(sub, length(world_prior)),
                 sess = rep(sess, length(world_prior)),
                 cond = rep(cond, length(world_prior)))
  out
}


get.kl.div.blocks <- function(data, block, cond, sub, sess){
  ### using the data frame 'sum.dat.for.kl', compute, for each block and condition, 
  ### the kl divergence between the subject data and the world probabilites
  current.dat = data[data$sub == sub & data$sess == sess & data$cond == cond & data$block == block, ]
  current.dat = current.dat %>% unique()
  current.dat$door = as.numeric(current.dat$door)
  current.dat = current.dat[order(current.dat$door),]

  world.alpha = with(current.dat, n.y + world) 
  sub.alpha = with(current.dat, n.x + subs_world)
  D = lgamma(sum(world.alpha)) - lgamma(sum(sub.alpha)) - sum(lgamma(world.alpha)) + sum(lgamma(sub.alpha)) + sum((world.alpha - sub.alpha)*(digamma(world.alpha)-digamma(sum(world.alpha))))
  # see this derivation on https://github.com/cran/Compositional/blob/master/R/kl.diri.R
  world.0 = sum(world.alpha) # ease for interpreting their eq
  sub.0 = sum(sub.alpha)
  f = sum( (world.alpha - sub.alpha) * ( digamma(world.alpha) - digamma(world.0) ) ) + sum( lgamma(sub.alpha) - lgamma(world.alpha) ) + lgamma(world.0) - lgamma(sub.0)
  out = data.frame(sub = sub, sess = sess, block = block, cond = cond, D = D, f = f)
  out
}


##########################################################################################################################################
############################ TRANSITION MATRICES ########################################################################################
##########################################################################################################################################

make_transition_matrix <- function(data, ndoors){
  # compute transition matrix for data (typically, 1 trials worth)
  out.mat <- matrix(data = 0, nrow=ndoors, ncol=ndoors)
  ntrials = length(data$sub)
  for (i in 2:ntrials){
    out.mat[ data$door[i], data$door[i-1] ] = out.mat[ data$door[i], data$door[i-1] ] + 1
  }
  out.mat = out.mat / sum(out.mat)
}

compute_mat_dist <- function(x, y){
  # compute similarity between two matrices, x & y
  # using d2 definition from: https://math.stackexchange.com/questions/507742/distance-similarity-between-two-matrices
  sqrt(sum((x-y)^2))
}

get_trans_sim_by_cond <- function(data, ndoors){
  # for one subject, session and condition, get the similarity between matrices across trials, and return as
  # a data frame
  trials = unique(data$t)
  out = data.frame(dist=rep(NA, length(trials-1)))
  for (i in 2:length(trials)){
    # a = rbind(tail(data[data$t == trials[i-1],],1), data[data$t == trials[i], ])
    # if (i == 2){
    #   b = data[data$t == trials[i - 1], ]
    # } else {
    #   b = rbind(tail(data[data$t == trials[i-2],],1), data[data$t == trials[i-1],])
    # }
    a = data[data$t == trials[i], ]
    b = data[data$t == trials[i-1], ]
    if (length(a$t) > 1 & length(b$t) > 1){
      a = make_transition_matrix(a, ndoors)
      b = make_transition_matrix(b, ndoors)
      out$dist[i-1] = compute_mat_dist(a, b)
    }
  }  
  # assign all other required variables to data frame
  out = out %>% drop_na()
  out$sub = data$sub[1]
  out$sess = data$sess[1]
  out$t = 1:length(out$dist)
  out$cond = data$cond[1]
  out
}

