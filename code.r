#The following is the code used for the  models and figure...
#...found in Hutchison et al, JAMA Open, 2023

#Raw Data is not available for privacy reasons, however for all STAN models, ...
#...the data list is created from the summary data found in the tables ...
#...or written results of the manuscript.

#Load libraries
library(tidyverse)
library(magrittr)
library(modelr)
library(ggplot2)
library(cowplot)
library(rstan)
library(rstanarm)
library(bayesplot)
library(rethinking)
library(tidybayes)
library(ggrepel)
library(RColorBrewer)
library(gganimate)
library(gtsummary)
library(gt)
library(ggdist)

#in-house functions
#function for deriving the shortest credible interval from the posterior
HDILow <- function(x, HDI = 0.9) {
  sortedPts <- sort(x)
  ciIdxInc <- ceiling(HDI * length(sortedPts) )
  nCIs <- length(sortedPts) - ciIdxInc 
  ciWidth <- rep(0, nCIs)
  for (i in 1:nCIs) {
    ciWidth[ i ] <- sortedPts[i + ciIdxInc] - sortedPts[i] }
  HDImin <- sortedPts[which.min(ciWidth)]
  HDImax <- sortedPts[which.min(ciWidth) + ciIdxInc] 
  return(HDImin)
}
HDIHigh <- function(x, HDI = 0.9) {
  sortedPts <- sort(x)
  ciIdxInc <- ceiling(HDI * length(sortedPts))
  nCIs <- length(sortedPts) - ciIdxInc
  ciWidth <- rep(0, nCIs)
  for (i in 1:nCIs) {
    ciWidth[i] <- sortedPts[i + ciIdxInc] - sortedPts[i]}
  HDImin <- sortedPts[which.min(ciWidth)]
  HDImax <- sortedPts[which.min(ciWidth) + ciIdxInc] 
  return(HDImax)
}

#more parsimoniuos z_score function than base R
z_score <- function(x) {
  out <- (x - mean(x)) / sd(x)
}

#function for creating data lists for each comparison to run STAN models
dat_create <- function(pre, post, pre_n, post_n, prior_mean, prior_sd) {

  #create group data based on manuscript table 2 or 3 (depending on comparison)
  pre <- data.frame(rep(c(0, 1), times = c((pre_n - pre), pre)))
  post <- data.frame(rep(c(0, 1), times = c((post_n - post), post)))
  
  #create group labels
  pre$group <- 2
  post$group <- 1
  colnames(post)[1] <- 'conc'
  colnames(pre)[1] <- 'conc'

  #combine
  conc <- rbind.data.frame(pre,post)
  
  #create a list for STAN
  dat_conc <- list(conc = conc$conc,
                   group = conc$group,
                   N = nrow(conc),
                   p_m = prior_mean, #prior mean (based on pre values: ...
                                     #...use inv_logit(x) to get pre value)
                   p_sd = prior_sd)  #prior sd (slight regularization)
  return(dat_conc)
}

#function for contrast creating with the posterior draws
#one input 'x' is the model sampled in STAN
contrast_create <- function(x) {
  #extract posterior draws
  post_draws <- as.data.frame(extract.samples(x))
  #extract the first two columns (a.1 and a.2)
  post_draws2 <- post_draws[c(1:2)]
  #apply the inverse logit function to the draws
  post_draws3 <- data.frame(apply(post_draws2, 2, inv_logit))
  #create a new column with the difference between the two draws
  post_draws3$diff <- post_draws3$a.2 - post_draws3$a.1
  return(post_draws3)
}

#for speed with STAN hmc modelling
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#compile the logistic regression STAN model used for all comparisons
logistic_mod <- stan_model("logistic_model.stan")

#####MODELLING#####

###Proportion Modelling (Table 2)###

##Open ice
dat_open <- dat_create(108, 167, 231, 455, -0.2, 0.3)
str(dat_open)
#fit the model
mod_open <- sampling(logistic_mod, data = dat_open,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_open <- contrast_create(mod_open)
#check the contrasts
precis(df_open, 2, prob = 0.9)

##Offensive zone
dat_offensive <- dat_create(79, 167, 230, 457, -0.65, 0.3)
str(dat_offensive)
#fit the model
mod_offensive <- sampling(logistic_mod, data = dat_offensive,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_offensive <- contrast_create(mod_offensive)
#check the contrasts
precis(df_offensive, 2, prob = 0.9)

##neutral zone
dat_neutral <- dat_create(45, 95, 230, 457, -1.4, 0.3)
str(dat_neutral)
#fit the model
mod_neutral <- sampling(logistic_mod, data = dat_neutral,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_neutral <- contrast_create(mod_neutral)
#check the contrasts
precis(df_neutral, 2, prob = 0.9)

##defensive zone
dat_defensive <- dat_create(106, 195, 230, 457, -0.15, 0.3)
str(dat_defensive)
#fit the model
mod_defensive <- sampling(logistic_mod, data = dat_defensive,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_defensive <- contrast_create(mod_defensive)
#check the contrasts
precis(df_defensive, 2, prob = 0.9)

##body checks (all)
dat_bca <- dat_create(173, 307, 231, 455, 1.1, 0.3)
str(dat_bca)
#fit the model
mod_bca <- sampling(logistic_mod, data = dat_bca,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_bca <- contrast_create(mod_bca)
#check the contrasts
precis(df_bca, 2, prob = 0.9)

##body checks to the head
dat_bch <- dat_create(113, 149, 231, 455, -0.08, 0.3)
str(dat_bch)
#fit the model
mod_bch <- sampling(logistic_mod, data = dat_bch,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_bch <- contrast_create(mod_bch)
#check the contrasts
precis(df_bch, 2, prob = 0.9)

##body checks to the head by the shoulder
dat_bcs <- dat_create(78, 91, 231, 455, -0.7, 0.3)
str(dat_bcs)
#fit the model
mod_bcs <- sampling(logistic_mod, data = dat_bcs,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_bcs <- contrast_create(mod_bcs)
#check the contrasts
precis(df_bcs, 2, prob = 0.9)

##body checks to the head by the glove
dat_bcg <- dat_create(8, 8, 231, 455, -3.3, 0.3)
str(dat_bcg)
#fit the model
mod_bcg <- sampling(logistic_mod, data = dat_bcg,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_bcg <- contrast_create(mod_bcg)
#check the contrasts
precis(df_bcg, 2, prob = 0.9)

##body checks to the head by the arm
dat_bcar <- dat_create(27, 50, 231, 455, -2, 0.3)
str(dat_bcar)
#fit the model
mod_bcar <- sampling(logistic_mod, data = dat_bcar,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_bcar <- contrast_create(mod_bcar)
#check the contrasts
precis(df_bcar, 2, prob = 0.9)

##lateral hits to the head
dat_lh2h <- dat_create(80, 61, 231, 455, -0.65, 0.3)
str(dat_lh2h)
#fit the model
mod_lh2h <- sampling(logistic_mod, data = dat_lh2h,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_lh2h <- contrast_create(mod_lh2h)
#check the contrasts
precis(df_lh2h, 2, prob = 0.9)

##First period
dat_first <- dat_create(117, 158, 231, 440, 0.03, 0.3)
str(dat_first)
#fit the model
mod_first <- sampling(logistic_mod, data = dat_first,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_first <- contrast_create(mod_first)
#check the contrasts
precis(df_first, 2, prob = 0.9) 

##puck possession
dat_puck <- dat_create(46, 77, 231, 456, -1.1, 0.3)
str(dat_puck)
#fit the model
mod_puck <- sampling(logistic_mod, data = dat_puck,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_puck <- contrast_create(mod_puck)
#check the contrasts
precis(df_puck, 2, prob = 0.9)

##penalty called
dat_penalty <- dat_create(60, 96, 220, 457, -1, 0.3)
str(dat_penalty)
#fit the model
mod_penalty <- sampling(logistic_mod, data = dat_penalty,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_penalty <- contrast_create(mod_penalty)
#check the contrasts
precis(df_penalty, 2, prob = 0.9)

##forwards
dat_forward <- dat_create(149, 281, 231, 457, 0.6, 0.3)
str(dat_forward)
#fit the model
mod_forward <- sampling(logistic_mod, data = dat_forward,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_forward <- contrast_create(mod_forward)
#check the contrasts
precis(df_forward, 2, prob = 0.9)

##defense
dat_defense <- dat_create(77, 146, 231, 457, -0.75, 0.3)
str(dat_defense)
#fit the model
mod_defense <- sampling(logistic_mod, data = dat_defense,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_defense <- contrast_create(mod_defense)
#check the contrasts
precis(df_defense, 2, prob = 0.9)

##goalie
dat_goalie <- dat_create(5, 30, 231, 457, -3.8, 0.3)
str(dat_goalie)
#fit the model
mod_goalie <- sampling(logistic_mod, data = dat_goalie,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_goalie <- contrast_create(mod_goalie)
#check the contrasts
precis(df_goalie, 2, prob = 0.9)

###Incidence Modelling (Table 3)###

##Overall concussion incidence
#create data list
dat_inc <- dat_create(301, 516, 4920, 6232, -2.75, 0.3)
str(dat_inc)
#fit the model
mod_inc <- sampling(logistic_mod, data = dat_inc, 
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_inc <- contrast_create(mod_inc)
#check the contrasts
precis(df_inc, 2, prob = 0.9)

##Open ice
dat_i_open <- dat_create(108, 167, 4920, 6232, -3.75, 0.3)
str(dat_i_open)
#fit the mod_iel
mod_i_open <- sampling(logistic_mod, data = dat_i_open,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_i_open <- contrast_create(mod_i_open)
#multiply by 100 to get the incidence per 100 games
df_i_open <- df_i_open * 100
#check the contrasts
precis(df_i_open, 2, prob = 0.9)

##Offensive zone
dat_i_offensive <- dat_create(79, 167, 4920, 6232, -4, 0.3)
str(dat_i_offensive)
#fit the mod_iel
mod_i_offensive <- sampling(logistic_mod, data = dat_i_offensive,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_i_offensive <- contrast_create(mod_i_offensive)
#multiply by 100 to get the incidence per 100 games
df_i_offensive <- df_i_offensive * 100
#check the contrasts
precis(df_i_offensive, 2, prob = 0.9)

##neutral zone
dat_i_neutral <- dat_create(45, 95, 4920, 6232, -4.75, 0.3)
str(dat_i_neutral)
#fit the mod_iel
mod_i_neutral <- sampling(logistic_mod, data = dat_i_neutral,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_i_neutral <- contrast_create(mod_i_neutral)
#multiply by 100 to get the incidence per 100 games
df_i_neutral <- df_i_neutral * 100
#check the contrasts
precis(df_i_neutral, 2, prob = 0.9)

##defensive zone
dat_i_defensive <- dat_create(106, 195, 4920, 6232, -3.75, 0.3)
str(dat_i_defensive)
#fit the mod_iel
mod_i_defensive <- sampling(logistic_mod, data = dat_i_defensive,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_i_defensive <- contrast_create(mod_i_defensive)
#multiply by 100 to get the incidence per 100 games
df_i_defensive <- df_i_defensive * 100
#check the contrasts
precis(df_i_defensive, 2, prob = 0.9)

##body checks (all)
dat_i_bca <- dat_create(173, 307, 4920, 6232, -3.25, 0.3)
str(dat_i_bca)
#fit the mod_iel
mod_i_bca <- sampling(logistic_mod, data = dat_i_bca,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_i_bca <- contrast_create(mod_i_bca)
#multiply by 100 to get the incidence per 100 games
df_i_bca <- df_i_bca * 100
#check the contrasts
precis(df_i_bca, 2, prob = 0.9)

##body checks to the head
dat_i_bch <- dat_create(113, 149, 4920, 6232, -3.75, 0.3)
str(dat_i_bch)
#fit the mod_iel
mod_i_bch <- sampling(logistic_mod, data = dat_i_bch,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_i_bch <- contrast_create(mod_i_bch)
#multiply by 100 to get the incidence per 100 games
df_i_bch <- df_i_bch * 100
#check the contrasts
precis(df_i_bch, 2, prob = 0.9)

##body checks to the head by the shoulder
dat_i_bcs <- dat_create(78, 91, 4920, 6232, -4.2, 0.3)
str(dat_i_bcs)
#fit the mod_iel
mod_i_bcs <- sampling(logistic_mod, data = dat_i_bcs,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_i_bcs <- contrast_create(mod_i_bcs)
#multiply by 100 to get the incidence per 100 games
df_i_bcs <- df_i_bcs * 100
#check the contrasts
precis(df_i_bcs, 2, prob = 0.9)

##body checks to the head by the glove
dat_i_bcg <- dat_create(8, 8, 4920, 6232, -6.5, 0.3)
str(dat_i_bcg)
#fit the mod_iel
mod_i_bcg <- sampling(logistic_mod, data = dat_i_bcg,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_i_bcg <- contrast_create(mod_i_bcg)
#multiply by 100 to get the incidence per 100 games
df_i_bcg <- df_i_bcg * 100
#check the contrasts
precis(df_i_bcg, 2, prob = 0.9)

##body checks to the head by the arm
dat_i_bcar <- dat_create(27, 50, 4920, 6232, -5.25, 0.3)
str(dat_i_bcar)
#fit the mod_iel
mod_i_bcar <- sampling(logistic_mod, data = dat_i_bcar,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_i_bcar <- contrast_create(mod_i_bcar)
#multiply by 100 to get the incidence per 100 games
df_i_bcar <- df_i_bcar * 100
#check the contrasts
precis(df_i_bcar, 2, prob = 0.9)

##lateral hits to the head
##body checks to the head by the arm
dat_i_lh2h <- dat_create(80, 61, 4920, 6232, -4, 0.3)
str(dat_i_lh2h)
#fit the mod_iel
mod_i_lh2h <- sampling(logistic_mod, data = dat_i_lh2h,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_i_lh2h <- contrast_create(mod_i_lh2h)
#multiply by 100 to get the incidence per 100 games
df_i_lh2h <- df_i_lh2h * 100
#check the contrasts
precis(df_i_lh2h, 2, prob = 0.9)

##First period
dat_i_first <- dat_create(117, 158, 4920, 6232, -3.75, 0.3)
str(dat_i_first)
#fit the mod_iel
mod_i_first <- sampling(logistic_mod, data = dat_i_first,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_i_first <- contrast_create(mod_i_first)
#multiply by 100 to get the incidence per 100 games
df_i_first <- df_i_first * 100
#check the contrasts
precis(df_i_first, 2, prob = 0.9) 

##puck possession
dat_i_puck <- dat_create(46, 77, 4920, 6232, -4.75, 0.3)
str(dat_i_puck)
#fit the mod_iel
mod_i_puck <- sampling(logistic_mod, data = dat_i_puck,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_i_puck <- contrast_create(mod_i_puck)
#multiply by 100 to get the incidence per 100 games
df_i_puck <- df_i_puck * 100
#check the contrasts
precis(df_i_puck, 2, prob = 0.9)

##penalty called
dat_i_penalty <- dat_create(60, 96, 4920, 6232, -4.25, 0.3)
str(dat_i_penalty)
#fit the mod_iel
mod_i_penalty <- sampling(logistic_mod, data = dat_i_penalty,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_i_penalty <- contrast_create(mod_i_penalty)
#multiply by 100 to get the incidence per 100 games
df_i_penalty <- df_i_penalty * 100
#check the contrasts
precis(df_i_penalty, 2, prob = 0.9)

##forwards
dat_i_forward <- dat_create(149, 281, 4920, 6232, -3.5, 0.3)
str(dat_i_forward)
#fit the mod_iel
mod_i_forward <- sampling(logistic_mod, data = dat_i_forward,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_i_forward <- contrast_create(mod_i_forward)
#multiply by 100 to get the incidence per 100 games
df_i_forward <- df_i_forward * 100
#check the contrasts
precis(df_i_forward, 2, prob = 0.9)

##defense
dat_i_defense <- dat_create(77, 146, 4920, 6232, -4.15, 0.3)
str(dat_i_defense)
#fit the mod_iel
mod_i_defense <- sampling(logistic_mod, data = dat_i_defense,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_i_defense <- contrast_create(mod_i_defense)
#multiply by 100 to get the incidence per 100 games
df_i_defense <- df_i_defense * 100
#check the contrasts
precis(df_i_defense, 2, prob = 0.9)

##goalie
dat_i_goalie <- dat_create(5, 30,  4920, 6232, -6.8, 0.3)
str(dat_i_goalie)
#fit the mod_iel
mod_i_goalie <- sampling(logistic_mod, data = dat_i_goalie,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_i_goalie <- contrast_create(mod_i_goalie)
#multiply by 100 to get the incidence per 100 games
df_i_goalie <- df_i_goalie * 100
#check the contrasts
precis(df_i_goalie, 2, prob = 0.9)

###Sensitivity Analyses###
##The lateral hit to the head comparison was repeated... 
##...with several different priors prepresenting different...
##...possible prior beliefs about the proportion of lateral hits to the head

#prior of 10%
dat_lh2h10 <- dat_create(80, 61, 231, 455, -2, 0.3)
str(dat_lh2h10)
#fit the model
mod_lh2h10 <- sampling(logistic_mod, data = dat_lh2h10,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_lh2h10 <- contrast_create(mod_lh2h10)
#check the contrasts
precis(df_lh2h10, 2, prob = 0.9)

#prior of 20%
dat_lh2h20 <- dat_create(80, 61, 231, 455, -1.5, 0.3)
str(dat_lh2h20)
#fit the model
mod_lh2h20 <- sampling(logistic_mod, data = dat_lh2h20,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_lh2h20 <- contrast_create(mod_lh2h20)
#check the contrasts
precis(df_lh2h20, 2, prob = 0.9)

#prior of 25%
dat_lh2h25 <- dat_create(80, 61, 231, 455, -1, 0.3)
str(dat_lh2h25)
#fit the model
mod_lh2h25 <- sampling(logistic_mod, data = dat_lh2h25,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_lh2h25 <- contrast_create(mod_lh2h25)
#check the contrasts
precis(df_lh2h25, 2, prob = 0.9)

#prior of 50%
dat_lh2h50 <- dat_create(80, 61, 231, 455, 0, 0.3)
str(dat_lh2h50)
#fit the model
mod_lh2h50 <- sampling(logistic_mod, data = dat_lh2h50,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_lh2h50 <- contrast_create(mod_lh2h50)
#check the contrasts
precis(df_lh2h50, 2, prob = 0.9)

#prior of 60%
dat_lh2h60 <- dat_create(80, 61, 231, 455, 0.5, 0.3)
str(dat_lh2h60)
#fit the model
mod_lh2h60 <- sampling(logistic_mod, data = dat_lh2h60,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_lh2h60 <- contrast_create(mod_lh2h60)
#check the contrasts
precis(df_lh2h60, 2, prob = 0.9)

#prior of 75%
dat_lh2h75 <- dat_create(80, 61, 231, 455, 1, 0.3)
str(dat_lh2h75)
#fit the model
mod_lh2h75 <- sampling(logistic_mod, data = dat_lh2h75,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_lh2h75 <- contrast_create(mod_lh2h75)
#check the contrasts
precis(df_lh2h75, 2, prob = 0.9)

###Poststratification Analysis###
#Post-stratification analysis for missing data incongruency...
#...between 2006-2010 (18.9% missing) and 2014-2019 (11.4% missing)...
#...if we matched 2006-2010 to 2014-2019, we would need to add ...
#...data to make up for the 7.5% more missing data...
#...So, there are 231 concussions with video data from 2006-2010...
#...out of a possible 301. To arrive at 11.4% missing data from 2006-2010...
#...(to match 2014-2019) we would have to add 22 concussion,...
#...reducing the 'no video' concussion number from 57 to 35...
#...since 35/301 is ~ 11.6% vs. 57/301 = 18.9%...
#...and to make up the correct proportion of events that were...
#...from lateral hits to the head, (34%) the constitution of these...
#...hits should be ~8 concussions from lateral hits to the head,...
#...and 14 concussions that were not lateral hits to the head.

#revised data list for poststratification analysis
dat_lh2hps <- dat_create(88, 61, 253, 455, 1, 0.3) #added as per above
str(dat_lh2hps)
#fit the model
mod_lh2hps <- sampling(logistic_mod, data = dat_lh2hps,
iter = 1000, chains = 4, cores = 4)
#create contrasts
df_lh2hps <- contrast_create(mod_lh2hps)
#check the contrasts
precis(df_lh2hps, 2, prob = 0.9)

###Figure in Manuscript###

##proportion plot prepping##
#create a data frame with the posterior contrasts from the models of proportion
diff_prop_df <- data.frame(df_open$diff, df_offensive$diff, df_neutral$diff, 
                      df_defensive$diff, df_bch$diff, df_bcs$diff, df_lh2h$diff, 
                      df_first$diff, df_puck$diff, df_penalty$diff)

#rename columns
colnames(diff_prop_df) <- c('Open Ice', 'Offensive','Neutral','Defensive',
                            'BC to Head (Total)','BC to Head (Shoulder)',
                            'Lateral Hit to Head', 'First Period',
                            'Puck Possession','Penalty Called')

#gather different groupings
#location
diff_prop_location <- gather(diff_prop_df[c(1:4)])
diff_prop_location$group <- "Location"

#mechanism
diff_prop_mechanism <- gather(diff_prop_df[c(5:7)])
diff_prop_mechanism$group <- 'Mechanism'

#Game Situation
diff_prop_game <- gather(diff_prop_df[c(8:10)])
diff_prop_game$group <- "Game"

#combine
contrast_prop_df <- rbind.data.frame(diff_prop_location, diff_prop_mechanism, 
                                      diff_prop_game)
colnames(contrast_prop_df) <-c('Characteristic',
        'Difference in Proportion: 2014 - 2019 vs. 2006 - 2010 (%)','group')
contrast_prop_df$group <- factor(contrast_prop_df$group,
levels = c('Game','Location','Mechanism'))

#reverse signs to match manuscript figure (for clarity of interpretation)
contrast_prop_df$`Difference in Proportion: 2014 - 2019 vs. 2006 - 2010 (%)` <-
- contrast_prop_df$`Difference in Proportion: 2014 - 2019 vs. 2006 - 2010 (%)`

#set theme for tidybayes
theme_set(theme_tidybayes() + panel_border())

#make plot
contrast_prop_fig <- contrast_prop_df %>%
ggplot(aes(y = `Characteristic`,
x = `Difference in Proportion: 2014 - 2019 vs. 2006 - 2010 (%)`,
fill = stat(x < 0 ))) +
geom_vline(xintercept = 0,alpha= 0.8,linetype =2 ) +
theme(strip.text.x = element_text(face = 'bold',size = 10),
axis.title.x = element_text(size = 12,face = 'bold'),
axis.title.y = element_blank(),
axis.text.y = element_text(face = 'bold',size = 10),
legend.position = 'none',
axis.text.x = element_text(face = 'bold',size = 10),
strip.text.y = element_blank())+
ylab("Characteristic")+ xlab('Difference in Proportion (percentage points)')+
stat_halfeye(.width = c(0.70,0.90)) +
facet_grid(group ~., scales="free",space = "free")+
scale_fill_manual(values = c("darkgray",'lightblue')) + 
labs(fill='Probability Mass < 0')
contrast_prop_fig

##incidence plot prepping##
#create a data frame with the posterior contrasts from the models for incidence
diff_inc_df <- data.frame(df_i_open$diff, df_i_offensive$diff, df_i_neutral$diff, 
                      df_i_defensive$diff, df_i_bch$diff, df_i_bcs$diff, 
                      df_i_lh2h$diff, df_i_first$diff, df_i_puck$diff, 
                      df_i_penalty$diff)
#rename
colnames(diff_inc_df) <- colnames(diff_prop_df)

#gather different groupings
#location
diff_inc_location <- gather(diff_inc_df[c(1:4)])
diff_inc_location$group <- "Location"

#mechanism
diff_inc_mechanism <- gather(diff_inc_df[c(5:7)])
diff_inc_mechanism$group <- 'Mechanism'

#Game Situation
diff_inc_game <- gather(diff_inc_df[c(8:10)])
diff_inc_game$group <- "Game"

#combine
contrast_inc_df <- rbind.data.frame(diff_inc_location, diff_inc_mechanism, 
                                    diff_inc_game)
colnames(contrast_inc_df) <-c('Characteristic',
       'Difference in Incidence: 2014 - 2019 vs. 2006 - 2010 (per 100 games)', 'group')
contrast_inc_df$group <- factor(contrast_inc_df$group,
levels = c('Game','Location','Mechanism'))

#reverse signs to match manuscript figure (for clarity of interpretation)
contrast_inc_df$`Difference in Incidence: 2014 - 2019 vs. 2006 - 2010 (per 100 games)` <-
- contrast_inc_df$`Difference in Incidence: 2014 - 2019 vs. 2006 - 2010 (per 100 games)`

#set theme for tidybayes
theme_set(theme_tidybayes() + panel_border())

#make plot
contrast_inc_fig <- contrast_inc_df %>%
ggplot(aes(y = `Characteristic`,
x = `Difference in Incidence: 2014 - 2019 vs. 2006 - 2010 (per 100 games)`,
fill = stat(x < 0 ))) +
geom_vline(xintercept = 0,alpha= 0.8,linetype =2 ) +
theme(strip.text.x = element_text(face = 'bold',size = 10),
axis.title.x = element_text(size = 12,face = 'bold'),
axis.title.y = element_blank(),
axis.text.y = element_text(face = 'bold',size = 10),
legend.position = 'none',
axis.text.x = element_text(face = 'bold',size = 10),
strip.text.y = element_blank())+
ylab("Characteristic")+ xlab('Difference in Incidence Rate (per 100 games)')+
stat_halfeye(.width = c(0.70,0.90)) +
facet_grid(group ~., scales="free",space = "free")+
scale_fill_manual(values = c("darkgray",'orange')) + labs(fill='Probability Mass < 0')
contrast_inc_fig