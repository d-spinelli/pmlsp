#' Fitting Spatial Probit Models by Partial Maximum Likelihood
#'
#' @description
#' Fits spatial autoregressive probit models using partial maximum likelihood
#' as proposed in Bille and Leorato (2020) <doi:10.1080/07474938.2019.1682314>.
#' Both Spatial Autoregressive (SAR) probit models and Spatial Autoregressive
#' Probit Models with Autoregressive Disturbances (SARAR) are supported.
#'
#' @details
#' The SARAR model is defined as
#' \deqn{
#' y^* = \rho W y + X \beta + W_2 Z \gamma + u
#' }{
#' y* = rho W y + X beta + W2 Z gamma + u
#' }
#' with
#' \deqn{
#' u = \lambda M u + \epsilon,
#' }{
#' u = lambda M u + epsilon
#' }
#' where \eqn{\epsilon \sim N(0,\sigma^2)} and the observed binary outcome is
#' \eqn{y = 1(y^* > 0)}. The SAR model is obtained by setting \eqn{\lambda = 0}.
#'
#' @usage
#' pmlsbp(
#'   formula, data, model = "SAR", grouping = 2,
#'   W = NULL,
#'   zero.policy = spatialreg::get.ZeroPolicyOption(),
#'   M = NULL, formula_xlag = NULL, W2 = NULL,
#'   method_inv = "solve",
#'   start = NULL, subset = NULL, na.action = na.fail,
#'   qu = Inf, iterlim = 1000,
#'   mvtnorm_control = list(M = 1e4, sim_type = "qmc",
#'                          tol = .Machine$double.eps, fast = FALSE),
#'   finalHessian = ifelse(method == "bhhh", method, TRUE),
#'   method = "bhhh", print.level = 2,
#'   vce.type = "asy",
#'   Conley = list(coords = NULL, LM = 2),
#'   nBoot = 1e3, spectral = FALSE,
#'   verbose = TRUE ,
#'   tol.solve = .Machine$double.eps,
#'   version = 0, ...
#' )
#'
#' @param formula An object of class \code{"formula"} specifying the model.
#' @param data An optional data frame containing the variables in the model.
#' @param model A character string specifying the model type; either
#'   \code{"SAR"} or \code{"SARAR"}.
#' @param grouping An integer defining the number of observations to include
#'   in the tuples used for estimation.
#' @param W A spatial weights object of class \code{"listw"}, for example
#'   created by \code{\link[spdep]{nb2listw}}.
#' @param zero.policy Logical; if \code{TRUE}, assigns zero to lagged values
#'   for observations without neighbors.
#' @param M A spatial weights object of class \code{"listw"} used for the
#'   disturbance process; relevant only if \code{model = "SARAR"}.
#' @param formula_xlag An optional \code{"formula"} specifying covariates to be
#'   spatially lagged.
#' @param W2 A spatial weights object used to construct spatially lagged
#'   covariates; relevant only if \code{formula_xlag} is not \code{NULL}.
#' @param method_inv Character string specifying the method used to invert
#'   \eqn{I - \rho W}; one of \code{"solve"}, \code{"chol"}, or \code{"fast"}.
#' @param start Numeric vector of starting values.
#' @param subset Optional vector specifying a subset of observations.
#' @param na.action A function specifying how missing values are handled.
#' @param qu Integer used only if \code{method_inv = "fast"}.
#' @param mvtnorm_control A list of control parameters for multivariate normal
#'   probability calculations; see \code{\link[mvtnorm]{lpmvnorm}}.
#' @param finalHessian Logical or character specifying how the final Hessian
#'   matrix is computed; see \code{\link[maxLik]{maxLik}}.
#' @param method Maximization method passed to \code{\link[maxLik]{maxLik}}.
#' @param print.level Integer controlling the amount of diagnostic output.
#' @param vce.type Character specifying the variance–covariance estimator;
#'   one of \code{"asy"}, \code{"bootstrap"}, or \code{"mConley"}.
#' @param Conley A list specifying options for the modified Conley estimator.
#' @param nBoot Integer specifying the number of bootstrap replications.
#' @param spectral Logical; if \code{TRUE}, uses spectral normalization of \eqn{W}.
#' @param verbose Logical; if \code{TRUE}, diagnostic messages are printed during
#'   model estimation. The default is \code{TRUE}. Setting \code{verbose = FALSE}
#'   suppresses intermediate output.
#' @param iterlim Integer specifying the maximum number of iterations.
#' @param tol.solve Numeric tolerance used when inverting matrices.
#' @param version Integer controlling the computation of bivariate normal
#'   integrals when \code{grouping = 2}.
#' @param ... Additional arguments passed to \code{\link[maxLik]{maxLik}}.
#'
#' @return
#' An object of class \code{"pmlsprobit"}, which is a list containing:
#' \describe{
#'   \item{\code{beta}}{Named vector of regression coefficients.}
#'   \item{\code{call}}{The matched function call.}
#'   \item{\code{code}}{Convergence code from \code{maxLik}.}
#'   \item{\code{estimate}}{Named vector of parameter estimates.}
#'   \item{\code{f}}{List of matrices used for prediction.}
#'   \item{\code{maximum}}{Value of the partial log-likelihood.}
#'   \item{\code{message}}{Convergence message from \code{maxLik}.}
#'   \item{\code{model}}{Model frame used in estimation.}
#'   \item{\code{model.type}}{Either \code{"SAR"} or \code{"SARAR"}.}
#'   \item{\code{rho}}{Estimated spatial autoregressive parameter.}
#'   \item{\code{slx}}{Logical indicating the presence of spatially lagged covariates.}
#'   \item{\code{start}}{Starting values used in estimation.}
#'   \item{\code{terms}}{Terms object.}
#'   \item{\code{vcov}}{Variance–covariance matrix of the estimates.}
#'   \item{\code{W}}{Spatial weights matrix \eqn{W}.}
#'   \item{\code{lambda}}{Estimated disturbance parameter (SARAR only).}
#'   \item{\code{M}}{Spatial weights matrix \eqn{M} (SARAR only).}
#'   \item{\code{W2}}{Spatial weights matrix \eqn{W_2} (if \code{slx = TRUE}).}
#' }
#'
#' @seealso
#' \code{\link{pmlsbp}}, \code{\link[maxLik]{maxLik}}
#'
#' @references
#' Bille', A. G., & Leorato, S. (2020).
#' Partial ML estimation for spatial autoregressive nonlinear probit models
#' with autoregressive disturbances.
#' \emph{Econometric Reviews}, 39(5), 437–475.
#' \doi{10.1080/07474938.2019.1682314}
#'
#' @examples
#' \donttest{
#' data(oldcol, package = "spdep")
#' dat <- COL.OLD
#' dat$y <- as.numeric(dat$CRIME > 35)
#'
#' listw <- spdep::nb2listw(COL.nb, style = "W")
#'
#' set.seed(857489)
#' mod <- pmlsbp(
#'   y ~ HOVAL + INC, data = dat, W = listw,
#'   model = "SAR", grouping = 7,
#'   na.action = na.omit, spectral = TRUE,
#'   iterlim = 1e5
#' )
#' }
#'
#' @export

pmlsbp <- function(formula,data, model="SAR", grouping=2, W=NULL,zero.policy =spatialreg::get.ZeroPolicyOption(),
                         M=NULL, formula_xlag=NULL, W2=NULL,  method_inv="solve",
                         start=NULL, subset=NULL, na.action=na.fail,qu=Inf, iterlim=1000,
                         mvtnorm_control=list(M=1e4, sim_type="qmc" , tol = .Machine$double.eps, fast = FALSE),
                         finalHessian=ifelse( method=="bhhh",method,TRUE),
                        method="bhhh",print.level=2, vce.type='asy' ,
                        Conley=list(coords=NULL,LM=2),nBoot=1e3 , spectral=FALSE, verbose=TRUE,
                        tol.solve= .Machine$double.eps, version=0, ...) {
  match.arg(model, c("SAR","SARAR"), several.ok = FALSE)
  match.arg(method_inv, c("chol","solve","fast"), several.ok = FALSE)
  stopifnot(is.numeric(qu))
  stopifnot(inherits(formula,"formula"))
  match.arg(vce.type, c("asy", "mConley", "bootstrap"), several.ok = F)



  if (is.null(mvtnorm_control$M))
    mvtnorm_control$M <- 1e4
  if (is.null(mvtnorm_control$tol))
    mvtnorm_control$tol <- .Machine$double.eps
  if (is.null(mvtnorm_control$fast))
    mvtnorm_control$fast <- FALSE
  if (is.null(mvtnorm_control$sim_type))
    mvtnorm_control$sim_type <- "qmc"
  stopifnot(mvtnorm_control$sim_type %in% c("mc","qmc"))

  #serve per chiamare model.frame con relative opzioni, senza generare errori
  call <- match.call()
  defaults  <- formals(sys.function())

  for (nm in names(defaults)) {
    if (is.null(call[[nm]]) & !is.null(defaults[[nm]])) {
      call[[nm]] <- defaults[[nm]]
    }
  }
  mmf <- match.call(expand.dots = FALSE)
  mmf$na.action <- na.action
  mmf$subset <- subset
  m <- match(c("formula", "data","na.action","subset"), names(mmf), 0L)
  mf <- mmf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)

  gav <- mf
  gav.formula <- all.vars(formula)

  if (!is.null(formula_xlag)) {
    mf.xlag <- mf
    mf.xlag$formula <- formula_xlag
    gav.formula <- c(gav.formula,all.vars(formula_xlag))
    mf.xlag <- eval(mf.xlag, parent.frame())
    mt.xlag <- attr(mf.xlag, "terms")
  }

  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  gav$formula <- reformulate(gav.formula)
  gav <- eval(gav,parent.frame())

  if (is.empty.model(mt))
    stop("'formula' is an empty model")
  y <- stats::model.response(mf, "any")
  if (is.character(y) | is.factor(y))
    stop ('the response variable must be numeric')
  stopifnot(length(unique(y))==2)
  y <- as.integer(as.logical(y))
  if (!all(y %in% c(0,1))) {
    stop('failed to convert the response variable into (0,1) values,
         please check the response variable')
  }
  X <- stats::model.matrix(mt, mf)
  if (nrow(X)!=length(y))
    stop("independent variables and response variable have different dimensions")
  if (is.null(subset))
    subset <- rep(T, nrow(X))
  if (is.null(W))
    stop("argument W is mandatory")
  if (!inherits(W,"matrix") & !inherits(W,"listw"))
    stop("argument W must be a matrix or a listw")
  if (inherits(W, "listw")) {
    if (!any(W$style==c("W","B","C")))
      stop("the style of W must be W or B")
    styleW<-W$style
    W <- subset(W, !(1:length(W$neighbours) %in% c(attr(mf, "na.action"), which(!subset))),zero.policy = zero.policy)
    W <- spdep::listw2mat (W)
  } else {
    W <- subset(W, !(1:dim(W)[2] %in% c(attr(mf, "na.action"), which(!subset))))
    }
  eig <-  list(W=range(eigen(W, symmetric = F)$values))
  if (styleW=='B' & spectral) {
    W <- W/max(abs(eig$W))
    eig <-  list(W=range(eigen(W, symmetric = F)$values))
  }
  if (nrow(W)!=length(y))
    stop("spatial matrix W and response variable have different dimensions")

  if (model=="SARAR" )   {
    if ( is.null(M) )
      stop("M is mandatory for SARAR model")
    if (!inherits(M,"matrix") & !inherits(M,"listw"))
      stop("M must be a matrix or a listw")
    if (inherits(M, "listw")) {
      if (!any(M$style==c("W","B","C")))
        stop("the style of M must be W or B")
      M <- subset(M, !(1:length(M$neighbours) %in% c(attr(mf, "na.action"), which(!subset))),zero.policy = zero.policy)
      M <- spdep::listw2mat(M)
    } else {
      M <- subset(M, !(1:dim(M)[2] %in% c(attr(mf, "na.action"), which(!subset))))
      }
    eig$M <- range(eigen(M, symmetric = F)$values)

    if (nrow(M)!=length(y))
      stop("spatial matrix M and response variable have different dimensions")
  }

  if (!is.null(formula_xlag)) {
    if ( is.null(W2) ) {
      W2 <- W
      styleW2 <- styleW
    }
    else {
      if (!inherits(W2,"matrix") & !inherits(W2,"listw"))
      stop("W2 must be a matrix or a listw")
    if (inherits(W2, "listw")) {
       if (!any(W2$style==c("W","B","C")))
        stop("the style of W2 must be W or B")
      styleW2 <-  W2$style
      W2 <- subset(W2, !(1:length(W2$neighbours) %in% c(attr(mf, "na.action"), which(!subset))),zero.policy = zero.policy)
      W2  <- spdep::listw2mat(W2)
    } else   {
      W2 <- subset(W2, !(1:dim(W2)[2] %in% c(attr(mf, "na.action"), which(!subset))))
      if (all(rowSums(W2)==1))
        styleW2 <-  "W"
        }
    }


    Xlag <- stats::model.matrix(mt.xlag, mf.xlag)
    attr.Xlag <- attributes(Xlag)
    attr.names <- names(attributes(Xlag))
    attr.names <- attr.names[attr.names != 'names']


    Xlag <- W2%*%stats::model.matrix(mt.xlag, mf.xlag)
    #cn <-  colnames(Xlag)
    if (styleW2=="W" && any(colnames(Xlag)=="(Intercept)")) {
       #cn <- cn[-which(colnames(Xlag)=="(Intercept)")]
       Xlag <- as.matrix(Xlag[,-which(colnames(Xlag)=="(Intercept)")])
       attr.Xlag$dimnames[[2]] <- (attr.Xlag$dimnames[[2]])[-which(attr.Xlag$dimnames[[2]]=="(Intercept)")]
       attr.Xlag$dim[2] <- attr.Xlag$dim[2]-1
       attr.Xlag$assign <- attr.Xlag$assign[-which(attr.Xlag$assign==0) ]
       #attr(attributes(mod$model$mf.xlag)[[2]], 'intercept')
    }

      attributes(Xlag)[attr.names] <- attr.Xlag[attr.names]
      attr(Xlag, 'dimnames')[[2]] <-   paste0("W2*",attr(Xlag, 'dimnames')[[2]])
    if (nrow(Xlag)!=length(y))
      stop("spatial lags and response variable have different dimensions")
    Xall <- cbind(X,Xlag)
    if (nrow(W2)!=length(y))
      stop("spatial matrix W2 and response variable have different dimensions")
  }
  else Xall <- X
  if (is.null(start)) {
    start <- suppressWarnings(glm(y~Xall+0, family = binomial(link = "probit"))$coefficients)
    names(start) <- colnames(Xall)

    start <- c(start,to_tilde(0, eig$W)[1])
    names(start)[length(start)] <- "rho"
    if (model=="SARAR") {
      start <- c(start,to_tilde(0, eig$M)[1])
      names(start)[length(start)] <- "lambda"
    }
  }
  if ("mConley" %in% vce.type) {
    if (is.null(Conley$coords))
      stop("Conley score requires coordinates data")
    stopifnot(is.data.frame(Conley$coords) | is.matrix(Conley$coords) )
    stopifnot(Conley$LM>1)
    Conley$coords <- subset(Conley$coords, !(1:nrow(Conley$coords) %in% c(attr(mf, "na.action"), which(!subset))))
  }


 groups <- group_matrices(y,Xall,grouping, M=mvtnorm_control$M, sim_type= mvtnorm_control$sim_type)

   if (method!="bobyqa") {
     m.maxlik <- match(c("method", "iterlim","tol","reltol","gradtol", "steptol","lambdatol","qac",
                         "qrtol", "marquardt_lambda0", "marquardt_lambdaStep" , "marquardt_maxLambda"  , "print.level" ,
                         "nm_alpha", "nm_beta", "nm_gamma", "sann_cand", "sann_temp","sann_tmax","sann_randomSeed","finalHessian","constraints" ), names(call), 0L)
     optimizer.call <- call[c(1L,m.maxlik)]
     optimizer.call[[1L]] <- quote(maxLik::maxLik)
     if (version==0 & grouping==2 & nrow(W)%%2==0 )
       optimizer.call$logLik <- ifelse(model=="SAR" ,quote(logLIK_SAR3), quote(logLIK_SARAR))
     else
       optimizer.call$logLik <- ifelse(model=="SAR" ,quote(logLIK_SAR), quote(logLIK_SARAR))
     optimizer.call$start <- quote(start)
     #optimizer.call$tol.solve <-quote(tol.solve)
   } else {
     optimizer.call <- call[c(1L)]
     optimizer.call[[1L]] <- quote(minqa::bobyqa)
     optimizer.call$par <- quote(start)
     #optimizer.call$feval <-quote(iterlim)
     optimizer.call$fn <- ifelse(model=="SAR" ,quote(logLIK_SAR), quote(logLIK_SARAR))
     optimizer.call$bobyqa <- quote(TRUE)
     optimizer.call$control <- quote(list(iprint=print.level , maxfun=iterlim) )
   }

   optimizer.call$X <-  quote(Xall)
   optimizer.call$y <-  quote(y)
   optimizer.call$eig <-  quote(eig)
   optimizer.call$W <-  quote(W)
   optimizer.call$qu <-  quote(qu)
   optimizer.call$method_inv <- quote(method_inv)
   optimizer.call$groups <- quote(groups)
   optimizer.call$mvtnorm_control <-  quote(mvtnorm_control)

   if (model!="SAR")
     optimizer.call$M <- quote(M)
   return <- eval(optimizer.call)


   return$rho <- to_natural(tail(return$estimate,1), eig$W )
   if (method=="bobyqa") {
     return$code<- return$ierr
     return$estimate <- return$par
   }

  return$groups <- groups
  return$start <- start

  obFun.call <- optimizer.call
  m.obFun <- match(c("X","y","eig","W","qu","method_inv","groups","mvtnorm_control","M"  ), names(optimizer.call), 0L)
  obFun.call <- optimizer.call[c(1L,m.obFun)]
  obFun.call[[1L]] <- ifelse(model=="SAR" ,quote(logLIK_SAR), quote(logLIK_SARAR))
  obFun.call$theta <- quote(return$estimate)


  obFun <- eval(obFun.call)

  return$beta <- return$estimate[1:ncol(Xall)]
  if (model=="SAR")  {
    return$rho <- to_natural(tail(return$estimate,1), eig$W )
    return$rho_tilde <- return$estimate[length(return$estimate)]
    return$estimate[length(return$estimate)] <- return$rho
    G <- diag(c(rep(1,length(return$estimate)-1), attr(return$rho,"jacobian")))
  }
  else {
    return$rho <- to_natural(return$estimate[length(return$estimate)-1], eig$W )
    return$lambda <- to_natural(tail(return$estimate,1), eig$M )
    return$rho_tilde <- return$estimate[length(return$estimate)-1]
    return$lambda_tilde <- return$estimate[length(return$estimate)]
    return$estimate[length(return$estimate)-1] <- return$rho
    return$estimate[length(return$estimate)] <- return$lambda
    G <- diag(c(rep(1,length(return$estimate)-2), attr(return$rho,"jacobian"), attr(return$lambda,"jacobian")))
  }

  ##maxlik uses rho_tilde for the maximization problem, to obtain the standard error the gradient must be adjusted to work with rho using the cain rule
if (method=='bobyqa') {
  return$gradientObs<-attr(obFun, 'gradient')
  return$hessian <- quote(-numericHessian(fn, t0=return$estimate, X=Xall, y = y, eig = eig, W = W, qu = qu, method_inv = method_inv,
                                    groups = groups, mvtnorm_control = mvtnorm_control, bobyqa=T))
  return$hessian[[2]][[2]] <- obFun.call[[1]]
  return$hessian<- eval(return$hessian)
}
  invH <- solve(as(return$hessian,"Matrix")) ##### devo tirare fuori hessiano
  V <-  list()
  if ("asy" %in% vce.type) {
    J <- crossprod(return$gradientObs)
    V <- G%*%invH%*%J%*%invH%*%t(G)
  }

  if ("mConley" %in% vce.type) {
      mean_coord <- aggregate(Conley$coords, list(return$groups$y), FUN=mean)
      knn <- spdep::knearneigh(mean_coord, k=Conley$LM-1)$nn     #matrix of knn for each site (pair of points)
      knn <-  cbind(matrix(1:nrow(knn), ncol=1), knn)
      KM <- 2*(1-abs(0:(Conley$LM-1)/Conley$LM))
      J <- matrix(0,ncol(return$gradientObs),ncol(return$gradientObs))
      for (j in 2:ncol(knn)) {
        J <- J+KM[j]*t(return$gradientObs)%*%return$gradientObs[knn[,j],]
      }
      J <-  (J+crossprod(return$gradientObs))
      V <- G%*%invH%*%J%*%invH%*%t(G)
    }     ## M, LM e tau?
  if ("bootstrap" %in% vce.type) { #bootstrap qui
    xb <- attributes(return$maximum)$f$Inv_rho%*%Xall%*%return$beta
    eta <- t(mvtnorm::rmvnorm(nBoot,sigma=attr(obFun, "f")$Sigma_rho))
    ystar <- t(tcrossprod(rep(1, ncol(eta)), xb) ) + eta
    ystar <- ifelse(ystar>0,1,0)
    bootstrap.call <- optimizer.call
    bootstrap.call$start <- quote(return$estimate)
    bootstrap.call$print.level <- quote(0)
    boot_estimate <- matrix(NA, nrow=nBoot, ncol=length(return$estimate))
    if (verbose) {
    cat("bootstrap iterations\n")
    pb <- txtProgressBar(min = 0, max = nBoot, initial = 0,style = 3)
           for (i in 1:nBoot) {
            setTxtProgressBar(pb,i)
            bootstrap.call$y <- quote(ystar[,i])
            bi <- try(eval(bootstrap.call), silent = T)
            if (!( "try-error" %in% class(bi)))
             boot_estimate[i,] <- bi$estimate
           }
    }
    boot_estimate <- na.omit(boot_estimate)
    V <- var(boot_estimate)
    #colnames(V$bootstrap) <- rownames(V$bootstrap) <- names(return$estimate)
    close(pb)
  }
  colnames(V) <- rownames(V) <- names(return$estimate)
  return$vcov <- V
  attr(return$vcov,"type") <- vce.type
  return$eig<-eig
  return$f <- attr(obFun, "f")
  return$model.type <- model
  return$slx <- !is.null(formula_xlag)
  return$W <- W
  return$terms <- mt
  return$model <- list(mf=mf)
  if (model=="SARAR")
    return$M <- M
  if (!is.null(formula_xlag)) {
    return$terms.xlag <- mt.xlag
    return$W2 <- W2
    return$model$mf.xlag <- mf.xlag
  }
  return$call <- call

  keep.results <- c('beta','estimate','call','rho','code','f','maximum','message','model','model.type','slx','start','terms','vcov','W', 'eig')
  if (!is.null(formula_xlag))  keep.results <- c(keep.results, 'W2','terms.xlag')
  if (model=="SARAR") keep.results <- c(keep.results,'M','lambda')
  ret<-return[keep.results]
  class(ret) <- "pmlsprobit"
  ret
  }

to_natural <- function(tilde, eig) {
  ret <-  eig[1]^(-1)+ ((eig[2]^(-1)-eig[1]^(-1))/(1+exp(-tilde)))
  attr(ret,"jacobian") <- (1/(ret-eig[1]^-1)-1/(ret-eig[2]^-1))^-1 ##drho_tilde/drho
  ret
}
to_tilde <- function(x,eig) {
  ret <-  -log((eig[2]^(-1)-eig[1]^(-1))/(x-eig[1]^(-1))-1)
  attr(ret,"jacobian") <-  1/(x-eig[1]^-1)-1/(x-eig[2]^-1) ##drho/drho-tilde (dret/dx)
  ret

}

#' Summarizing Partial Maximum Likelihood Spatial Probit Models
#'
#' @description
#' Provides a summary method for objects of class \code{"pmlsprobit"}.
#'
#' @param object An object of class \code{"pmlsprobit"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' An object of class \code{"summary.pmlsprobit"} containing:
#' \describe{
#'   \item{\code{estimate}}{Matrix of estimates, standard errors, z-values, and p-values.}
#'   \item{\code{rho}}{Estimated spatial autoregressive parameter.}
#'   \item{\code{lambda}}{Estimated spatial error parameter.}
#'   \item{\code{model.type}}{Type of spatial probit model fitted.}
#'   \item{\code{loglik}}{Log-likelihood at the optimum.}
#'   \item{\code{iteration}}{Number of iterations performed.}
#'   \item{\code{returnCode}}{Optimizer convergence code.}
#'   \item{\code{returnMessage}}{Message associated with the convergence code.}
#'   \item{\code{vcov}}{Variance–covariance matrix of the parameter estimates.}
#' }
#'
#' @seealso
#' Partial maximum likelihood estimation is implemented in
#' \code{\link{pmlsbp}}.
#'
#' @examples
#' \donttest{
#' data(oldcol, package = "spdep")
#' dat <- COL.OLD
#' dat$y <- as.numeric(dat$CRIME > 35)
#'
#' listw <- spdep::nb2listw(COL.nb, style = "W")
#'
#' set.seed(857489)
#' mod <- pmlsbp(
#'   y ~ HOVAL + INC, data = dat, W = listw,
#'   model = "SAR", grouping = 7,
#'   na.action = na.omit, spectral = TRUE,
#'   iterlim = 1e5
#' )
#' summary(mod)
#' }
#'
#'
#' @method summary pmlsprobit
#' @export
summary.pmlsprobit <- function(object, ... ) {
  if(!inherits(object, "pmlsprobit"))
    stop("'summary.pmlsprobit' called on a non-'pmlsprobit' object")
  vcov <- object$vcov

  se <- sqrt(diag(vcov))
  z <- object$estimate/se

  results <- cbind("Estimate"=object$estimate, "Std. error"=se, "z value"=z, "Pr(> z)" = 2*pnorm( -abs( z)) )
  summary <- list(  estimate=results,
                   rho=object$rho,
                   lambda=object$lambda,
                   model.type=object$model.type,
                   loglik=object$maximum,
                   iteration=object$iteration,
                   returnCode=object$code,
                   returnMessage=object$message ,
                   vcov=vcov)
  class(summary) <- "summary.pmlsprobit"
  summary
}

#' Print method for summary.pmlsprobit objects
#'
#' S3 method to print objects of class \code{summary.pmlsprobit}.
#'
#' @param x An object of class \code{summary.pmlsprobit}.
#' @param digits Number of significant digits to print for coefficients.
#' @param ... Additional arguments (currently ignored).
#'
#' @return No return value
#' @export
#' @method print summary.pmlsprobit
print.summary.pmlsprobit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  stopifnot(inherits(x, "summary.pmlsprobit"))

  cat("Partial Maximum Likelihood estimation\n")
  cat(x$model, " model\n")
  cat("Partial Log-Likelihood = ", x$loglik, "\n")

  if (!is.null(x$returnMessage) && nzchar(x$returnMessage)) {
    cat(x$returnMessage, "\n")
  }

  cat("---\n")

  if (!is.null(x$estimate)) {
    printCoefmat(x$estimate, digits = digits, ...)
  }

  invisible(x)
}

#' Variance–Covariance Matrix for Partial Maximum Likelihood Spatial Probit Models
#'
#' @description
#' Extracts the variance–covariance matrix from an object of class
#' \code{"pmlsprobit"}.
#'
#' @param object An object of class \code{"pmlsprobit"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' A numeric matrix containing the variance–covariance matrix of the parameter
#' estimates.
#'
#' @seealso
#' Partial maximum likelihood estimation is implemented in
#' \code{\link{pmlsbp}}.
#'
#' @method vcov pmlsprobit
#' @export
vcov.pmlsprobit <- function(object,...) {
  V<-object$vcov
  colnames(V) <- rownames(V) <- names(object$estimate)
  V
}

#' Coefficient Estimates for Partial Maximum Likelihood Spatial Probit Models
#'
#' @description
#' Extracts the estimated regression coefficients from an object of class
#' \code{"pmlsprobit"}.
#'
#' @param object An object of class \code{"pmlsprobit"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' A named numeric vector containing the coefficient estimates.
#'
#' @seealso
#' Partial maximum likelihood estimation is implemented in
#' \code{\link{pmlsbp}}.
#'
#' @method coef pmlsprobit
#' @export
coef.pmlsprobit <- function(object,...) {
  object$estimate
}

#' @describeIn pmlsbp prints the call and the coefficient estimates of a \code{"pmlsprobit"} object
#' @param x a \code{"pmlsprobit"} object
#' @param digits integer, used for number formatting the decimal digits
#' @return No return value
#' @export
print.pmlsprobit<-function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n",
                         collapse = "\n"), "\n\n", sep = "")
  if (length(coef.pmlsprobit(x))) {
    cat(paste0("Coefficients (",x$model.type," model) :\n"))
    print.default(format(coef.pmlsprobit(x), digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}


