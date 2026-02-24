



####################################################################
##
## This script relies on data: "Data/Sonde_Clean/SV_clean-23072025.csv"
##                             "Data/Sonde_Clean/SV_clean-MLE23072025.csv"
##
## Script modified on 25 Jul 2025 by L. Genzoli (laurel.genzoli@gmail.com)
## This script models metabolism for 6 years for the full year at SV (Klamath)
##
## Script first calculates metab using the bayesisan solution for GPP, ER, and K600
## Second model run uses fixed K600 from first run, then estimates GPP and ER with MLE
##
##
####################################################################



################################################
## Code is for loading streamMetabolizer if its not loaded:
################################################
options(mc.cores = parallel::detectCores())
Reinstall unitted and streamMetabolizer if needed
remotes::install_github('appling/unitted', force = TRUE)
remotes::install_github("GLEON/LakeMetabolizer", force = TRUE)
remotes::install_github("USGS-R/streamMetabolizer", force = TRUE)


################################################
## Packages:
################################################

library(tidyverse)
library(lubridate)
library(rstan)
library(streamMetabolizer)

#setwd("R:/Blaszczak_Lab/Ongoing Projects/LG/Klamath_Metab")
#setwd('/Users/laurelgenzoli/Dropbox/2024_Projects/Mckinney_Fire_Metab/Mck-Fire-Metab-GH')
################################################



################################################
##
## BAYES RUN TO GET K600 ESTIMATES!!!!
##
################################################

## Data

SV.run.metab<- read_csv("Data/SV_clean-23072025.csv")%>% 
  select(DO.obs, solar.time, light, DO.sat, depth, discharge, temp.water)



################################################
## Function from Alice to define the nodes used for binning discharge:
## This is part of setting the specs.
## Set K600 priors in the function below!!!
## Used mean from Genzoli and Hall 2012 at SV
################################################

set_Q_nodes <- function(specs_poolBinned, discharge){
  Qrange = quantile(log(discharge),
                    probs = c(0.1, 0.9), na.rm = T)
  n = 4
  delta = (Qrange[2]-Qrange[1])/n
  while(delta > 1){
    n = n + 1
    delta <- (Qrange[2]-Qrange[1])/n}
  
  nodes <- seq(Qrange[1], Qrange[2], length.out = n)
  specs_poolBinned$K600_lnQ_nodes_centers <- nodes
  specs_poolBinned$K600_lnQ_nodes_meanlog <- rep(2.3, n) ## This is my K600 prior-mean
  specs_poolBinned$K600_lnQ_nodes_sdlog = 0.8 ## This is my K600 prior-sd (I made this wider)
  specs_poolBinned$K600_daily_sigma_sigma = 0.05 ## This is my K600 prior-sigma
  return(specs_poolBinned)}

## Check the priors:

values <- as.data.frame(rnorm(1000, mean = exp(2.3), sd = exp(0.8)))
colnames(values) <- "value"
ggplot(values, aes(value))+geom_histogram()


###############################################
#   Priors for K600: mean from Genzoli and Hall 2012 at SV (2.34)
#   mean.log.k = 2.34
#
#   read.table('/Users/laurelgenzoli/Dropbox/Klamath/Klamath_Revisions/Submit_EA/Klamath_Metab_2012_data.txt', header = T)%>%
#   filter(site == "SV")%>% 
#   summarize(mean.log.k = log(mean(K600pred, na.rm = T)))
###############################################



###############################################
## Set Parameters for sM
## Run model (this takes a loooooong time)
## Save the model output
################################################

name_poolBinned <- mm_name(type = "bayes", 
                           pool_K600 = "binned", 
                           err_obs_iid = T, err_proc_iid =T)
specs_poolBinned <- specs(name_poolBinned, 
                          burnin_steps=1000, saved_steps=1000,
                          n_cores=4, verbose=T) 

specs_poolBinned <- set_Q_nodes(specs_poolBinned, SV.run.metab$discharge)
plot(density(log(SV.run.metab$discharge), na.rm = T))
abline(v = specs_poolBinned$K600_lnQ_nodes_centers)

mm_SV <- metab(specs_poolBinned, data=SV.run.metab)
saveRDS(mm_SV, file = "Data/Metab_Output/Metab_Out/SV_BINNED.RDS") 

GPP.ER <- get_params(mm_SV, uncertainty = "ci")
mod.fits <- predict_DO(mm_SV)
fit.time <- as_tibble(as.list(get_fitting_time(mm_SV)))
metab.fit.rda <- get_fit(mm_SV)[-2]

write_csv(GPP.ER, "Data/Metab_Output/Metab_Out/SV_GPP_ER_BIN.csv")
write_csv(mod.fits, "Data/Metab_Output/Metab_Out/SV_fits_BIN.csv")
write_csv(fit.time, "Data/Metab_Output/Metab_Out/SV_fit-time_BIN.csv")
save(metab.fit.rda, file = "Data/Metab_Output/Metab_Out/SV_RDA_BIN.rda")






############################################
####                                    ####
####      MLE runs starts below!!       ####
####                                    ####
############################################


############################################
####
####    Run metab with MLE and hand-pooled K600
####    This section of the script was run after first half of metab_QA script
####    SV.run.dd comes from K600 estimates from Bayes run above
####
############################################


SV.run.metab<- read_csv("Data/SV_clean-MLE23072025.csv")%>% 
  select(DO.obs, solar.time, light, DO.sat, depth, temp.water)

SV.run.dd<- read_csv("Data/SV_clean-MLE23072025.csv")%>% 
  select(solar.time, K600.daily)%>%
  mutate(date = as.Date(solar.time))%>%
  group_by(date)%>%
  summarize(K600.daily = mean(K600.daily, na.rm = T))


mle.model <- metab_mle(specs(mm_name("mle")), 
                       data = SV.run.metab,
                       data_daily = SV.run.dd)

GPP.ER <- get_params(mle.model, uncertainty = "ci")
mod.fits <- predict_DO(mle.model)
metab.fit.rda <- get_fit(mle.model)[-2]

write_csv(GPP.ER, "Data/Metab_Output/MLE/SV_GPP_ER_BIN.csv")
write_csv(mod.fits, "Data/Metab_Output/MLE/SV_fits_BIN.csv")
save(metab.fit.rda, file = "Data/Metab_Output/MLE/SV_RDA_BIN.rda")



############################################
####
####    QA of MLE model output continues in metab_QA script
####
############################################


