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

mod1a <- pmlsbp (equation, COL.OLD, W=listw, model = 'SAR',
                 mvtnorm_control=list(M=1e3, sim_type="qmc" , tol = .Machine$double.eps, fast = FALSE),
                 print.level = 2,grouping=7)
summary(mod1a)

mod1b <- pmlsbp (equation, COL.OLD, W=listw, model = 'SAR',formula_xlag = xlag,
                 print.level = 2,grouping=7)
summary(mod1b)

mod2 <- pmlsbp (equation, COL.OLD, W=listw,M=listw, model = 'SARAR',print.level = 2,
                mvtnorm_control=list(M=1e3, sim_type="qmc" , tol = .Machine$double.eps, fast = FALSE),
                grouping=7,method='bhhh')
summary(mod2)
xlag <-  y~PLUMB
mod2a <- pmlsbp (equation, COL.OLD, W=listw,M=listw, model = 'SAR',print.level = 2,
                 formula_xlag = xlag, iterlim=500,
                 mvtnorm_control=list(M=1e3, sim_type="mc" , tol = .Machine$double.eps, fast = FALSE),
                 grouping=7,method='bhhh')
summary(mod2a)

