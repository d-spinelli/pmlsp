library(Matrix)
library(matrixcalc)
library(maxLik)
library(mvtnorm)
library(qrng)
library(spatialreg)
library(spdep)
library(numDeriv)
library(abind)

data(oldcol, package="spdep")

listw <- spdep::nb2listw(COL.nb, style="W")
COL.OLD$y<-as.numeric(COL.OLD$CRIME>35)
equation <- y~HOVAL+INC
xlag <-  y~PLUMB
set.seed(857489)

mod1 <- pmlsbp (equation, COL.OLD, W=listw, model = 'SAR',print.level = 2,grouping=7)
summary(mod1)
head(predict(mod1), digits=3)
mfx<-predict(mod1, 'margins')
head(mfx$Direct)
head(mfx$Indirect)
head(mfx$Total)
ape(mod1)

mod2 <- pmlsbp (equation, COL.OLD, W=listw, model = 'SAR',
                 mvtnorm_control=list(M=1e3, sim_type="qmc" , tol = .Machine$double.eps, fast = FALSE),
                 print.level = 2,grouping=7)
summary(mod2)


data("eurostat")
set.seed(12345)

eurostat$bin.emp <- ifelse(eurostat$employment >= 70, 1, 0)
eurostat$gdp.std <- as.numeric(scale(eurostat$gdp))

mod3 <- pmlsbp(bin.emp ~ 0 + gdp.std + isced_02+isced_34+isced_58 ,
               eurostat,
               W=eurostat.nb,
               zero.policy=T, M=eurostat.nb,
               grouping=2, model='SARAR')
summary(mod3)
ape(mod3)


eqn<- bin.emp ~ 0 + gdp.std + isced_02+isced_34+isced_58
mod4 <- pmlsbp(eqn  ,
               eurostat,
               W=eurostat.nb,
               mvtnorm_control=list(M=1e3, sim_type="qmc" , tol = .Machine$double.eps, fast = FALSE),
               zero.policy=T, W2=eurostat.nb,formula_xlag =bin.emp ~ 0 + gdp.std + isced_02+isced_34  ,
               grouping=6, model='SAR')
summary(mod4)
ape(mod4)
head(predict(mod4, 'prob'))

