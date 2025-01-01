###For KNB
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()

library(brms)
library(rstan)
library(tidybayes)

dat<-read.csv("flobuds_KNB.csv",header = TRUE)

dat$Light<-ifelse(dat$Light=="S",0,1) ## dummy variable treatments
dat$Force<-ifelse(dat$Force=="C",0,1)

dat<-filter(dat, !GEN.SPA %in% c("AME.SPP","BET.SPP", "BET.ALL","ACE.SAC")) ## remove species with no flowering


#####leaf budburst model #########
##################################
bb.int<-get_prior(budburst.9.~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = dat, family = gaussian()) # check default priors

mod.bb.int<-brm(budburst.9. ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+
                  (Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                data = dat, family = gaussian(),
                iter= 4000,
                warmup = 3000) 

if(FALSE){ ### diffus-ify priors to make sure they aren't influencing inference
  
  altpriorz<-c(prior(student_t(3,33,100),class="Intercept"),
               prior(student_t(3, 0, 100),class="sigma"),
               prior_string("normal(0,100)", class = "b")) ## set wider alternative priors to check they are weakly informative
  
  mod.bb.int.alprior<-brm(budburst.9. ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+
                            (Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                          data = dat, family = gaussian(),prior=altpriorz,chains=2,
                          iter= 4000,
                          warmup = 3000)  ##run model with alternative priors
  
  fixef(mod.bb.int.alprior) ## compare results
}
 fixef(mod.bb.int)


#####flowering model #########
##################################
flo.int<-get_prior(flo_day~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = dat, family = gaussian()) ## check default priors

mod.flo.int<-brm(flo_day~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                 data = dat, family = gaussian(),
                 iter= 4000,
                 warmup = 3000)  ##run model


if(FALSE){ # check flower model with diffuse priors too
mod.flo.int.priorz<-brm(flo_day~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                 data = dat, family = gaussian(),prior=altpriorz,
                 iter= 4000,
                 warmup = 3000,chains=2)  ##run model
fixef(mod.flo.int.priorz)}

fixef(mod.flo.int)
#####leafout model #########
##################################

lo.int<-get_prior(leaf_day.15.~Chill+Light+Force+Chill:Light+Chill:Force+Force:Light,data = dat, family = gaussian()) ## check default priors

mod.lo.int<-brm(leaf_day.15. ~ Chill+Light+Force+Chill:Light+Chill:Force+Force:Light+(Chill+Light+Force+Chill:Light+Chill:Force+Force:Light|GEN.SPA),
                data = dat, family = gaussian(),
                iter= 4000,
                warmup = 3000)   


#####################################
summary(mod.flo.int) ## check Rhats and N_Effs
summary(mod.bb.int)
summary(mod.lo.int)
########################

### check some mpre diagnostics
bayestestR::effective_sample(mod.flo.int,effects = c("all"))
bayestestR::effective_sample(mod.bb.int,effects = c("all"))
bayestestR::effective_sample(mod.lo.int,effects = c("all"))

pp_check(mod.bb.int,nsamples = 100)
pp_check(mod.flo.int, nsamples = 100)
pp_check(mod.bb.int, nsamples=100)


####################################

