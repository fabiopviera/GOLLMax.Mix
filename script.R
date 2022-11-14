
require(survival)
require(readxl)
require(gamlss)
require(gamlss.cens)


dat_ <- read_xlsx("Data/data_china.xlsx")

dat_$cens <- as.integer(dat_$cens)
dat_$age.c <- as.factor(dat_$age.c)
dat_$gender.c <- as.factor(dat_$gender.c)

head(dat_)

source("GAMLSS/GOLLMax.Mix.R")
source("GAMLSS/OLLMax.Mix.R")
source("GAMLSS/EMax.Mix.R")


fit.GOLLMax <- gamlss(Surv(time,cens)~gender.c+age.c,
                       tau.formula =~gender.c+age.c,
                       family=cens(GOLLMax.mixtura),
                       data=dat_,c.crit=0.0001,n.cyc=2000)

summary(fit.GOLLMax,type="qr")


fit.OLLMax <- gamlss(Surv(time,cens)~gender.c+age.c,
                   nu.formula =~gender.c+age.c,
                   family=cens(OLLMax.mixtura),
                       data=dat_,c.crit=0.0001,n.cyc=2000)

summary(fit.OLLMax,type="qr")


fit.EMax <- gamlss(Surv(time,cens)~gender.c+age.c,
                   nu.formula =~gender.c+age.c,
                   family=cens(EMax.mixtura),
                   data=dat_,c.crit=0.0001,n.cyc=2000)

summary(fit.EMax,type="qr")

