get.sub.level.plot <- function(data){
  # NOTE - feed in one sub at a time, or amend to facet wrap over subs
  # first find and add missing values/variables
  for (i in levels(data$cond)){
    # add any doors missing from each block
    miss.doors = setdiff(levels(data$door[data$cond == i]), unique(data$door[data$cond == i]))
    if (sum(lengths(miss.doors)) > 0){
      tmp = c()
      tmp$sub[1:length(miss.doors)] = as.character(data$sub[1])
      tmp$cond[1:length(miss.doors)] = i
      tmp$block[1:length(miss.doors)] = as.character(data$block[1])
      tmp$door = miss.doors
      tmp$n[1:length(miss.doors)] = 0
      tmp$look_freq[1:length(miss.doors)] = 0
      tmp$num_p[1:length(miss.doors)] = 0
      data = rbind(data, tmp)
    }
  }
  
  plot.data = data %>% group_by(sess, cond, door) %>%
                   summarise(look_freq = mean(look_freq),
                   num_p = num_p[1])
  plot.data$door <- factor(plot.data$door, 
                    levels = as.character(c(1:16)))
  p <- 
    plot.data %>%
    ggplot(aes(door, look_freq, colour=cond, fill=cond)) + 
    geom_histogram(stat="identity", binwidth=16, alpha=0.5) +
    facet_wrap(~sess*cond) +
    geom_line(data=plot.data, aes(door, num_p, group=1), color=wesanderson::wes_palette("BottleRocket2")[c(4)]) +
    scale_color_manual(values=wesanderson::wes_palette("BottleRocket2")[c(1:2)]) +
    scale_fill_manual(values=wesanderson::wes_palette("BottleRocket2")[c(1:2)]) +
    theme_classic() 
  p
}

plot.kls.from.list <- function(data){
  # take the list output 'kl.over.blocks' and plot using same colour scheme as above
  plot.dat = data.frame(cond=rep(c("cert", "uncert"), each=8),
                           block=rep(seq(1,8,1), times=2),
                           kl=unlist(data[1,]))
  p <- plot.dat %>%
    ggplot(aes(block, kl, colour=cond, fill=cond)) +
    geom_line() + scale_colour_manual(values=wesanderson::wes_palette("BottleRocket2")[c(1:2)]) +
    theme_classic()
  
}
  
######################################################################################################
############################ plot transition distance densities ######################################
######################################################################################################
plot.distances <- function(data){
  data %>% ggplot(aes( y=dist, colour=sess, fill=sess)) +
    geom_density(alpha=0.4) +
    facet_wrap(~sub*cond) +
    theme(axis.title.x = element_text(face = "italic")) +
    scale_fill_manual(values=wes_palette("IsleofDogs1")[c(2:1)]) +
    scale_colour_manual(values=wes_palette("IsleofDogs1")[c(2:1)]) +
    theme_cowplot()
}

