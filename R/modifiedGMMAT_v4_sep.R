# library(Rcpp)
# source("modifiedDeseq2_glmm.R")
#sourceCpp("DESeq2_v1.cpp")
#sourceCpp("fitglmm.cpp")
glmmkin <- function(sizefactor,fit0,modelMatrix,weights,fixed,disp,ns, data = parent.frame(), kins = NULL, id, random.slope = NULL, groups = NULL, method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, taumin = 1e-5, taumax = 1e5, tauregion = 10, verbose = FALSE, ...) {
  call <- match.call()
  if(!is.null(kins) && !class(kins) %in% c("matrix", "list")) {
    if(is.null(attr(class(kins), "package"))) stop("Error: \"kins\" must be a matrix or a list.")
    else if(attr(class(kins), "package") != "Matrix") stop("Error: if \"kins\" is a sparse matrix, it must be created using the Matrix package.")
  }
  if(!method %in% c("REML", "ML"))
    stop("Error: \"method\" must be \"REML\" or \"ML\".")
  method.optim <- try(match.arg(method.optim, c("AI", "Brent", "Nelder-Mead")))
  if(class(method.optim) == "try-error")
    stop("Error: \"method.optim\" must be \"AI\", \"Brent\" or \"Nelder-Mead\".")
  if(method.optim == "AI" && method == "ML")
    stop("Error: method \"ML\" not available for method.optim \"AI\", use method \"REML\" instead.")
  if(method.optim == "Brent" && class(kins) == "list")
    stop("Error: method.optim \"Brent\" can only be applied in one-dimensional optimization, use a matrix for \"kins\".")
  if(method.optim != "AI" && ((!is.null(attr(class(kins), "package")) && attr(class(kins), "package") == "Matrix") || (class(kins) == "list" && any(sapply(kins, function(xx) !is.null(attr(class(xx), "package")) && attr(class(xx), "package") == "Matrix")))))
    stop("Error: sparse matrices can only be handled by method.optim \"AI\".")
  # if(class(family) != "family")
  #   stop("Error: \"family\" must be an object of class \"family\".")
  # if(!family$family %in% c("binomial", "gaussian", "Gamma", "inverse.gaussian", "poisson", "quasi", "quasibinomial", "quasipoisson"))
  #   stop("Error: \"family\" must be one of the following: binomial, gaussian, Gamma, inverse.gaussian, poisson, quasi, quasibinomial, quasipoisson.")
  # if(!is.null(groups)) {
  #   if(family$family != "gaussian") stop("Error: heteroscedastic linear mixed models are only applicable when \"family\" is gaussian.")
  #   if(method.optim != "AI") stop("Error: heteroscedastic linear mixed models are currently only implemented for method.optim \"AI\".")
  #   if(!groups %in% names(data)) stop("Error: \"groups\" must be one of the variables in the names of \"data\".")
  # }
  if(!id %in% names(data)) stop("Error: \"id\" must be one of the variables in the names of \"data\".")
  if(!is.null(random.slope)) {
    if(method.optim != "AI") stop("Error: random slope for longitudinal data is currently only implemented for method.optim \"AI\".")
    if(!random.slope %in% names(data)) stop("Error: \"random.slope\" must be one of the variables in the names of \"data\".")
  }
  if(!is.null(attr(class(kins), "package")) && attr(class(kins), "package") == "Matrix")  kins <- list(kins1 = kins)
  if(method.optim != "Brent" && class(kins) == "matrix") kins <- list(kins1 = kins)

  idx <- match(rownames(model.frame(formula = fixed, data = data, na.action = na.omit)), rownames(model.frame(formula = fixed, data = data, na.action = na.pass)))

  tmp.matrix<-model.frame(formula = fixed, data = data, na.action = na.omit)
  y<-t(unname(tmp.matrix[,1,drop=F]))
  #
  # if(ncol(tmp.matrix)==1){
  #   modelMatrix<-matrix(rep(1,nrow(tmp.matrix)),ncol=1)
  #   colnames(modelMatrix)<-"(Intercept)"
  # }else{
  #   modelMatrix<-tmp.matrix
  #   modelMatrix[,1]<-rep(1,nrow(tmp.matrix))
  #   colnames(modelMatrix)[1]<-"(Intercept)"
  # }
  # rownames(modelMatrix)<-seq(1:nrow(tmp.matrix))
  # <- glm(formula = fixed, data = data, family = family, ...)
  # fit0 <- fitNbinomGLMs(y=y,
  #                       sizefactor = sizefactor,
  #                       modelMatrix = modelMatrix,
  #                       alpha_hat=disp,
  #                       weight=weights[,,drop=FALSE])



  if(any(duplicated(data[idx, id]))) {
    cat("Duplicated id detected...\nAssuming longitudinal data with repeated measures...\n")

    if(method.optim == "Brent") {
      if(is.null(kins)) {
        kins <- diag(length(unique(data[idx, id])))
        rownames(kins) <- colnames(kins) <- unique(data[idx, id])
      } else stop("Error: method.optim \"Brent\" can only be applied to unrelated individuals in longitudinal data analysis.")
    } else {
      if(method.optim != "AI") kins[[length(kins) + 1]] <- diag(length(unique(data[idx,id])))
      else if(length(kins) > 0) kins[[length(kins) + 1]] <- Diagonal(n = length(unique(data[idx, id])))
      else kins <- list(kins1 = Diagonal(n = length(unique(data[idx, id]))))
      rownames(kins[[length(kins)]]) <- colnames(kins[[length(kins)]]) <- unique(data[idx, id])
    }
  } else if(!is.null(random.slope)) stop("Error: no duplicated \"id\" found, \"random.slope\" must be used for longitudinal data with duplicated \"id\".")
  if(class(kins) == "matrix") {
    match.idx1 <- match(data[idx, id], rownames(kins))
    match.idx2 <- match(data[idx, id], colnames(kins))
    if(any(is.na(c(match.idx1, match.idx2)))) stop("Error: kins matrix does not include all individuals in the data.")
    kins <- kins[match.idx1, match.idx2]
  } else if(class(kins) == "list") {
    for(i in 1:length(kins)) {
      match.idx1 <- match(data[idx, id], rownames(kins[[i]]))
      match.idx2 <- match(data[idx, id], colnames(kins[[i]]))
      if(any(is.na(c(match.idx1, match.idx2)))) stop("Error: kins matrix ", i, " does not include all individuals in the data.")
      kins[[i]] <- kins[[i]][match.idx1, match.idx2]
    }
  } else {
    if(!is.null(groups)) stop("Error: heteroscedastic linear models for unrelated observations have not been implemented.")
    y <- y[1,]
    n <- length(y)
    names(y)<-seq(1,length(y))
    offset <- log(sizefactor)
    if(is.null(offset)) offset <- rep(0, n)
    #family <- fit0$family
    beta<-fit0$betaMatrix
    eta<-modelMatrix%*%beta
    #eta <- fit0$linear.predictors
    mu<-as.vector(fit0$mu)
    mu.eta<-mu
    variance<-mu+disp*mu^2
    #mu <- fit0$fitted.values
    #mu.eta <- family$mu.eta(eta)
    Y <- eta - offset + (y - mu)/mu.eta
    sqrtW <- mu.eta/sqrt(1/as.vector(weights)*variance)
    X<-modelMatrix
    #X <- model.matrix(fit0)
    #alpha <- fit0$coef

    res <- y - mu
    tau <- 1
    Sigma_i <- diag(sqrtW^2/tau)
    rownames(Sigma_i) <- colnames(Sigma_i) <- rownames(X)
    fit <- list(theta=tau, n.groups=1, coefficients=beta, linear.predictors=eta, fitted.values=mu, Y=Y, X=X, P=NULL, residuals=res, scaled.residuals=res*as.vector(weights)/tau, cov=vcov(fit0), Sigma_i=Sigma_i,Sigma_iX=X*sqrtW^2/tau,converged=TRUE, call = call, id_include = data[idx, id])
    class(fit) <- "glmmkin"
    return(fit)
  }
  group.id <- if(is.null(groups)) rep(1, length(idx)) else data[idx, groups]
  time.var <- if(is.null(random.slope)) NULL else data[idx, random.slope]
  y=as.vector(y)
  fit <- glmmkin.fit(y=y,modelMatrix=modelMatrix,sizefactor=sizefactor,disp=disp,ns=ns,weights=weights,fit0, kins, time.var, group.id, method = method, method.optim = method.optim, maxiter = maxiter, tol = tol, taumin = taumin, taumax = taumax, tauregion = tauregion, verbose = verbose)
  fit$call <- call
  fit$id_include <- data[idx, id]
  class(fit) <- "glmmkin"
  return(fit)
}

glmmkin.fit <- function(y,modelMatrix,sizefactor,disp,ns,weights,fit0, kins, time.var, group.id, method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, taumin = 1e-5, taumax = 1e5, tauregion = 10, verbose = FALSE) {
  if(method.optim == "Brent") {
    fit <- glmmkin.brent(fit0, kins, method = method, maxiter = maxiter, tol = tol, taumin = taumin, taumax = taumax, tauregion = tauregion, verbose = verbose)
    if(fit$theta[2]/fit$theta[1] < 1.01 * tol) {
      warning("Variance estimate 0 observed, refitting model...", call. = FALSE)
      fit <- glmmkin.brent(fit0, kins, method = method, tau = 0, fixtau = 1, maxiter = maxiter, tol = tol, taumin = taumin, taumax = taumax, tauregion = tauregion, verbose = verbose)
    }
  } else {
    names(kins) <- paste("kins", 1:length(kins), sep="")
    if(method.optim == "AI") {
      group.unique <- unique(group.id)
      group.idx <- list()
      for(i in 1:length(group.unique)) group.idx[[i]] <- which(group.id == group.unique[i])
      covariance.idx <- fixrho.old <- fixrho.new <- NULL
      if(!is.null(time.var)) {
        q <- length(kins)
        covariance.idx <- matrix(0, q, 3)
        for(i in 1:q) kins[[q + i]] <- if(!is.null(attr(class(kins[[i]]), "package")) && attr(class(kins[[i]]), "package") == "Matrix") forceSymmetric(kins[[i]] * time.var + t(kins[[i]] * time.var)) else kins[[i]] * time.var + t(kins[[i]] * time.var)
        for(i in 1:q) {
          kins[[2*q + i]] <- if(!is.null(attr(class(kins[[i]]), "package")) && attr(class(kins[[i]]), "package") == "Matrix") forceSymmetric(t(kins[[i]] * time.var) * time.var) else t(kins[[i]] * time.var) * time.var
          covariance.idx[i, ] <- c(q + i, i, 2*q + i) + length(group.idx)
        }
        names(kins) <- paste("kins", 1:length(kins), sep="")
        fixrho.old <- rep(0, q)
      }
      fixtau.old <- rep(0, length(kins)+length(group.idx))
      fit <- glmmkin.ai(y=y,modelMatrix=modelMatrix,sizefactor=sizefactor,disp=disp,ns,weights=weights,fit0 = fit0, kins = kins, covariance.idx = covariance.idx, group.idx = group.idx, maxiter = maxiter, tol = tol, verbose = verbose)
      fixtau.new <- 1*(fit$theta < 1.01 * tol)
      if(!is.null(covariance.idx)) {
        fixtau.new[covariance.idx[, 1]] <- 0
        fixrho.new <- rep(0, nrow(covariance.idx))
        fixrho.idx <- apply(covariance.idx, 1, function(x) abs(fit$theta[x[1]]) > (1 - 1.01 * tol) * sqrt(fit$theta[x[2]] * fit$theta[x[3]]))
        fixrho.new[fixrho.idx] <- sign(fit$theta[covariance.idx[fixrho.idx, 1]])
      }
      while(any(fixtau.new != fixtau.old) || (!is.null(fixrho.new) && any(fixrho.new != fixrho.old))) {

        warning("Variance estimate on the boundary of the parameter space observed, refitting model...", call. = FALSE)
        fixtau.old <- fixtau.new
        if(!is.null(covariance.idx)) fixrho.old <- fixrho.new
        fit <- glmmkin.ai(y=y,modelMatrix=modelMatrix,sizefactor=sizefactor,disp=disp,ns,weights=weights,fit0 = fit0, kins = kins, covariance.idx = covariance.idx, group.idx = group.idx, fixtau = fixtau.old, fixrho = fixrho.old, maxiter = maxiter, tol = tol, verbose = verbose)
        fixtau.new <- 1*(fit$theta < 1.01 * tol)
        if(!is.null(covariance.idx)) {
          fixtau.new[covariance.idx[, 1]] <- 0
          fixrho.new <- rep(0, nrow(covariance.idx))
          fixrho.idx <- apply(covariance.idx, 1, function(x) abs(fit$theta[x[1]]) > (1 - 1.01 * tol) * sqrt(fit$theta[x[2]] * fit$theta[x[3]]))
          fixrho.new[fixrho.idx] <- sign(fit$theta[covariance.idx[fixrho.idx, 1]])
        }
      }
      if(!fit$converged) {warning("Average Information REML not converged", call. = FALSE)}
      # if(!fit$converged) {
      #   if(length(group.idx) != 1) stop("Error: Average Information REML not converged, cannot refit heteroscedastic linear mixed model using Brent or Nelder-Mead methods.")
      #   if(!is.null(time.var)) stop("Error: Average Information REML not converged, cannot refit random slope model for longitudinal data using Brent or Nelder-Mead methods.")
      #   if(length(kins) == 1) {
      #     warning("Average Information REML not converged, refitting model using Brent method...", call. = FALSE)
      #     fit <- glmmkin.brent(fit0, kins[[1]], method = method, maxiter = maxiter, tol = tol, taumin = taumin, taumax = taumax, tauregion = tauregion, verbose = verbose)
      #     if(fit$theta[2]/fit$theta[1] < 1.01 * tol) {
      #       warning("Variance estimate 0 observed, refitting model...", call. = FALSE)
      #       fit <- glmmkin.brent(fit0, kins[[1]], method = method, tau = 0, fixtau = 1, maxiter = maxiter, tol = tol, taumin = taumin, taumax = taumax, tauregion = tauregion, verbose = verbose)
      #     }
      #   } else {
      #     warning("Average Information REML not converged, refitting model using Nelder-Mead method...", call. = FALSE)
      #     fixtau.old <- rep(0, length(kins))
      #     fit <- glmmkin.nm(fit0, kins, method = method, maxiter = maxiter, tol = tol, verbose = verbose)
      #     fixtau.new <- 1*(fit$theta[-1]/fit$theta[1] < 1.01 * tol)
      #     while(any(fixtau.new != fixtau.old)) {
      #       warning("Variance estimate 0 observed, refitting model...", call. = FALSE)
      #       fixtau.old <- fixtau.new
      #       tau <- rep(1, length(kins))
      #       tau[which(fixtau.old == 1)] <- 0
      #       fit <- glmmkin.nm(fit0, kins, method = method, tau = tau, fixtau = fixtau.old, maxiter = maxiter, tol = tol, verbose = verbose)
      #       fixtau.new <- 1*(fit$theta[-1]/fit$theta[1] < 1.01 * tol)
      #     }
      #   }
      # }
    } else {
      fixtau.old <- rep(0, length(kins))
      fit <- glmmkin.nm(fit0, kins, method = method, maxiter = maxiter, tol = tol, verbose = verbose)
      fixtau.new <- 1*(fit$theta[-1]/fit$theta[1] < 1.01 * tol)
      while(any(fixtau.new != fixtau.old)) {
        warning("Variance estimate 0 observed, refitting model...", call. = FALSE)
        fixtau.old <- fixtau.new
        tau <- rep(1, length(kins))
        tau[which(fixtau.old == 1)] <- 0
        fit <- glmmkin.nm(fit0, kins, method = method, tau = tau, fixtau = fixtau.old, maxiter = maxiter, tol = tol, verbose = verbose)
        fixtau.new <- 1*(fit$theta[-1]/fit$theta[1] < 1.01 * tol)
      }
    }
  }
  return(fit)
}

glmmkin.ai <- function(y,modelMatrix,sizefactor,disp,ns,weights,fit0, kins, covariance.idx = NULL, group.idx, tau = rep(0, length(kins)+length(group.idx)), fixtau = rep(0, length(kins)+length(group.idx)), fixrho = NULL, maxiter = 500, tol = 1e-5, verbose = FALSE) {
  is.Matrix <- any(sapply(kins, function(xx) !is.null(attr(class(xx), "package")) && attr(class(xx), "package") == "Matrix"))
  y <- y
  n <- length(y)
  offset <- log(sizefactor)
  if(is.null(offset)) offset <- rep(0, n)
  #family <- fit0$family
  # if(!fix.intercept){
  #   eta <-as.vector( cbind(rep(1,nrow(modelMatrix)),modelMatrix)%*%(fit0$betaMatrix))
  #   beta <- as.vector(fit0$betaMatrix)[-1]
  # }else{
  eta <-as.vector( modelMatrix%*%(fit0$betaMatrix))
  beta <- as.vector(fit0$betaMatrix)
  #}
  mu <- as.vector(fit0$mu)
  mu.eta<-mu
  variance<-mu+disp*mu^2
  Y <- eta - offset + (y - mu)/mu.eta
  sqrtW <- mu.eta/sqrt(1/as.vector(weights)*variance)
  X <- modelMatrix
  p<-ncol(X)
  if(verbose) {
    cat("Fixed-effect coefficients:\n")
    print(beta)
  }
  #if(family$family %in% c("poisson", "binomial")) {
  tau[1] <- 1
  fixtau[1] <- 1
  #}
  q <- length(kins)
  ng <- length(group.idx)
  idxtau <- which(fixtau == 0)
  if(!is.null(covariance.idx)) idxtau2 <- intersect(covariance.idx[, 1], idxtau)
  q2 <- sum(fixtau == 0)
  if(q2 > 0) {
    tau[idxtau] <- rep(var(Y)/(q+ng), q2)
    if(!is.null(covariance.idx)) tau[idxtau2] <- 0
    diagSigma <- rep(0, n)
    for(i in 1:ng) diagSigma[group.idx[[i]]] <- tau[i]/sqrtW[group.idx[[i]]]^2
    for(i in 1:q) {tau[i+ng] <- tau[i+ng]/mean(diag(kins[[i]]))}
    # Sigma <- diag(diagSigma)
    # for(i in 1:q) {
    #   tau[i+ng] <- tau[i+ng]/mean(diag(kins[[i]]))
    #   Sigma <- Sigma + tau[i+ng]*kins[[i]]
    # }
    # #		if(is.Matrix) Sigma <- forceSymmetric(Sigma)
    #
    # Sigma_i <- chol2inv(chol(Sigma))
    debug1<<-"point 1_1"
    entry_Sigma_iXorY<<-  lapply(ns, function(sid){
      if(length(sid)==1){sub_Sigma<-matrix(diagSigma[sid])
      }else{sub_Sigma <- diag(diagSigma[sid])}
      #sub_Sigma <- diag(diagSigma[sid])
      #sub_kins<-matrix(1,nrow = length(sid),ncol=length(sid))
      sub_kins<<-lapply(kins,function(component){component[sid,sid,drop=F]})
      for(i in 1:q) sub_Sigma <<- sub_Sigma + tau[i+ng]*sub_kins[[i]]
      sub_Sigma_i <- chol2inv(chol(sub_Sigma))
      rm(sub_Sigma)
      sub_X<<-X[sid,,drop=F]
      sub_Sigma_iX <<- crossprod(sub_Sigma_i, sub_X)
      sub_Y<<-Y[sid]
      sub_Sigma_iY<<-crossprod(sub_Sigma_i,sub_Y)
      #return(list(sigma_i=sub_Sigma_i))
      return(cbind(sub_Sigma_iY,sub_Sigma_iX))
    })
    debug1<<-"point 1_2"
    Sigma_iX<-Reduce(rbind,entry_Sigma_iXorY)
    # Sigma_i_new<-bdiag(entry_Sigma_iXorY)
    # Sigma_iY<-as.matrix(crossprod(Sigma_i_new,Y))
    # Sigma_iX_new<-as.matrix(crossprod(Sigma_i_new, X))
    Sigma_iY<-Sigma_iX[,1,drop=F]
    Sigma_iX<-Sigma_iX[,-1,drop=F]
    entry_XSigma_iX<-lapply(ns, function(sid){
      sub_X<-X[sid,,drop=F]
      sub_XSigma_iX <- crossprod(sub_X, Sigma_iX[sid,,drop=F])
      return(sub_XSigma_iX)
    })
    cov <- chol2inv(chol(Reduce("+",entry_XSigma_iX)))

    Sigma_iXcov<- tcrossprod(Sigma_iX, cov)

    # entry_Sigma_i<-  sapply(ns, function(sid){
    #   if(length(sid)==1){sub_Sigma<-matrix(diagSigma[sid])
    #   }else{sub_Sigma <- diag(diagSigma[sid])}
    #   #sub_Sigma <- diag(diagSigma[sid])
    #   sub_kins<-matrix(1,nrow = length(sid),ncol=length(sid))
    #   for(i in 1:q) sub_Sigma <- sub_Sigma + tau[i+ng]*sub_kins
    #   sub_Sigma_i <- chol2inv(chol(sub_Sigma))
    #   rm(sub_Sigma)
    #   # sub_X<-X[sid,,drop=F]
    #   # sub_Sigma_iX <- crossprod(sub_Sigma_i, sub_X)
    #   # sub_Y<-Y[sid]
    #   # sub_Sigma_iY<-crossprod(sub_Sigma_i,sub_Y)
    #   return(list(sigma_i=sub_Sigma_i))
    #   #return(cbind(sub_Sigma_iY,sub_Sigma_iX))
    # })
    # Sigma_i<-bdiag(entry_Sigma_i)


    #rm(diagSigma)
    gc()
    #Sigma_iX <- crossprod(Sigma_i, X)

    #P <- Sigma_i - tcrossprod(tcrossprod(Sigma_iX, chol2inv(chol(crossprod(X, Sigma_iX)))), Sigma_iX)
    #rm(Sigma_i)
    #gc()
    #PY <- crossprod(P, Y)

    # XSigma_iX <- crossprod(X, Sigma_iX)
    # if(is.Matrix) XSigma_iX <- forceSymmetric(XSigma_iX)
    # cov <- chol2inv(chol(XSigma_iX))
    #Sigma_iXcov <- tcrossprod(Sigma_iX, cov)

    PY <- Sigma_iY - tcrossprod(Sigma_iX, t(crossprod(Sigma_iXcov, Y)))

    #PY <- crossprod(Sigma_i, Y) - tcrossprod(Sigma_iX, t(crossprod(Sigma_iXcov, Y)))
    tau0 <- tau
    score<-rep(NA,q2)
    for(i in 1:q2) {
      #if(idxtau[i] <= ng) tau[idxtau[i]] <- max(0, tau0[idxtau[i]] + tau0[idxtau[i]]^2 * (sum((PY/sqrtW)[group.idx[[idxtau[i]]]]^2) - sum((diag(P)/sqrtW^2)[group.idx[[idxtau[i]]]]))/n)
      if(idxtau[i] <= ng) tau[idxtau[i]] <- max(0, tau0[idxtau[i]] + tau0[idxtau[i]]^2 * (sum((PY/sqrtW)[group.idx[[idxtau[i]]]]^2) - sum(((diag(Sigma_i)-rowSums(Sigma_iX*Sigma_iXcov))/sqrtW^2)[group.idx[[idxtau[i]]]]))/n)
      else {
        #PAPY <- crossprod(P, crossprod(kins[[idxtau[i]-ng]], PY))
        #if(!is.null(covariance.idx)) tau[idxtau[i]] <- if(idxtau[i] %in% idxtau2) 0 else max(0, tau0[idxtau[i]] + tau0[idxtau[i]]^2 * (crossprod(Y, PAPY) - sum(P*kins[[idxtau[i]-ng]]))/n)
        #else tau[idxtau[i]] <- max(0, tau0[idxtau[i]] + tau0[idxtau[i]]^2 * (crossprod(Y, PAPY) - sum(P*kins[[idxtau[i]-ng]]))/n)
        APY <- crossprod(kins[[idxtau[i]-ng]], PY)
        #APY<-matrix(Reduce(c,sapply(ns,function(sid){rep(sum(PY[sid]),length(sid))})),ncol=1)

        entry_Sigma_iAPY<-  lapply(ns, function(sid){
          if(length(sid)==1){sub_Sigma<-matrix(diagSigma[sid])
          }else{sub_Sigma <- diag(diagSigma[sid])}
          #sub_Sigma <- diag(diagSigma[sid])
          sub_kins<-lapply(kins,function(component){component[sid,sid,drop=F]})
          for(comp in 1:q) sub_Sigma <- sub_Sigma + tau[comp+ng]*sub_kins[[comp]]
          sub_Sigma_i <- chol2inv(chol(sub_Sigma))
          rm(sub_Sigma)
          sum_Sigma_i<-sum(sub_Sigma_i*sub_kins[[i]])
          sub_Sigma_iAPY<-crossprod(sub_Sigma_i, APY[sid,,drop=F])
          return(cbind(sub_Sigma_iAPY,rep(sum_Sigma_i,length(sid))))
        })

        # entry_Sigma_iAPY<-  sapply(ns, function(sid){
        #   sub_Sigma_i <- Sigma_i_new[sid,sid]
        #   sub_Sigma_iAPY<-matrix(crossprod(sub_Sigma_i, APY[sid,]),ncol=1)
        #   return(sub_Sigma_iAPY)
        # })
        PAPY <- Reduce(rbind,entry_Sigma_iAPY)[,1,drop=F] - tcrossprod(Sigma_iX, t(crossprod(Sigma_iXcov, APY)))


        # APY <- crossprod(kins[[idxtau[i]-ng]], PY)
        # PAPY <- crossprod(Sigma_i, APY) - tcrossprod(Sigma_iX, t(crossprod(Sigma_iXcov, APY)))
        if(!is.null(covariance.idx)) tau[idxtau[i]] <- if(idxtau[i] %in% idxtau2) 0 else max(0, tau0[idxtau[i]] + tau0[idxtau[i]]^2 * (sum(Y* PAPY) - (sum(Sigma_i*kins[[idxtau[i]-ng]])-sum(Sigma_iX*crossprod(kins[[idxtau[i]-ng]],Sigma_iXcov))))/n)
        else {
          score[i] <- sum(Y * PAPY) - (sum(unique(Reduce(rbind,entry_Sigma_iAPY)[,2]))-sum(sapply(ns, function(sid){sum(Sigma_iX[sid,,drop=F]*crossprod(kins[[idxtau[i]-ng]][sid,sid],Sigma_iXcov[sid,,drop=F]))})))

          #score[i]<-sum(Y* PAPY) - (sum(Sigma_i*kins[[idxtau[i]-ng]])-sum(Sigma_iX*crossprod(kins[[idxtau[i]-ng]],Sigma_iXcov)))
          tau[idxtau[i]] <- max(0, tau0[idxtau[i]] + tau0[idxtau[i]]^2 * (score[i])/n)
        }

      }
    }
    #rm(P)
    #gc()
  }



  for (i in seq_len(maxiter)) {
    if(verbose) cat("\nIteration ", i, ":\n")
    beta0 <- beta
    tau0 <- tau
    #fit <- .Call(C_fitglmm_ai, Y, X, q, kins, ng, group.idx, sqrtW^2, tau, fixtau, tol)
    #if(is.Matrix) fit <- R_fitglmm_ai(Y, X, q, kins, ng, group.idx, sqrtW^2, tau, fixtau)
    # tic("R")
    # fit <- R_fitglmm_ai(Y, X, q, kins, ng, group.idx, sqrtW^2, tau, fixtau)
    # toc()
debug1<<-"point 1_3"
    fit <<- R_fitglmm_ai1(ns=ns,Y=Y, X=X, q=q, kins=kins, ng=ng, group.idx=group.idx, W=sqrtW^2, tau=tau, fixtau=fixtau)


    ##C++ version
    # tic("C")
    # fit <- fitglmm_ai(Y, X, q, kins, ng, group.idx, sqrtW^2, tau, fixtau)
    # toc()
    #else fit <- .Call(C_fitglmm_ai, Y, X, q, kins, ng, group.idx, sqrtW^2, tau, fixtau)

    if(q2 > 0) {
      debug1<<-"point 1_4"
      Dtau <- as.numeric(fit$Dtau)
      tau[idxtau] <- tau0[idxtau] + Dtau
      if(is.null(covariance.idx)) {
        tau[tau < tol & tau0 < tol] <- 0
        while(any(tau < 0)) {
          Dtau <- Dtau / 2
          tau[idxtau] <- tau0[idxtau] + Dtau
          tau[tau < tol & tau0 < tol] <- 0
        }
        tau[tau < tol] <- 0
      } else {
        fixrho.idx0 <- apply(covariance.idx, 1, function(x) abs(tau0[x[1]]) > (1 - 1.01 * tol) * sqrt(tau0[x[2]] * tau0[x[3]]))
        tau[-covariance.idx[, 1]][tau[-covariance.idx[, 1]] < tol & tau0[-covariance.idx[, 1]] < tol] <- 0
        if(any(fixrho != 0)) {
          idxrho <- which(fixrho != 0)
          tau[covariance.idx[idxrho, 1]] <- suppressWarnings(fixrho[idxrho] * sqrt(tau[covariance.idx[idxrho, 2]] * tau[covariance.idx[idxrho, 3]]))
        }
        fixrho.idx <- suppressWarnings(apply(covariance.idx, 1, function(x) abs(tau[x[1]]) > (1 - 1.01 * tol) * sqrt(tau[x[2]] * tau[x[3]])))
        tau[covariance.idx[fixrho.idx & fixrho.idx0, 1]] <- sign(tau[covariance.idx[fixrho.idx & fixrho.idx0, 1]]) * sqrt(tau[covariance.idx[fixrho.idx & fixrho.idx0, 2]] * tau[covariance.idx[fixrho.idx & fixrho.idx0, 3]])
        while(any(tau[-covariance.idx[, 1]] < 0) || any(apply(covariance.idx, 1, function(x) abs(tau[x[1]]) > sqrt(tau[x[2]] * tau[x[3]])))) {
          Dtau <- Dtau / 2
          tau[idxtau] <- tau0[idxtau] + Dtau
          tau[-covariance.idx[, 1]][tau[-covariance.idx[, 1]] < tol & tau0[-covariance.idx[, 1]] < tol] <- 0
          if(any(fixrho != 0)) {
            idxrho <- which(fixrho != 0)
            tau[covariance.idx[idxrho, 1]] <- suppressWarnings(fixrho[idxrho] * sqrt(tau[covariance.idx[idxrho, 2]] * tau[covariance.idx[idxrho, 3]]))
          }
          fixrho.idx <- suppressWarnings(apply(covariance.idx, 1, function(x) abs(tau[x[1]]) > (1 - 1.01 * tol) * sqrt(tau[x[2]] * tau[x[3]])))
          tau[covariance.idx[fixrho.idx & fixrho.idx0, 1]] <- sign(tau[covariance.idx[fixrho.idx & fixrho.idx0, 1]]) * sqrt(tau[covariance.idx[fixrho.idx & fixrho.idx0, 2]] * tau[covariance.idx[fixrho.idx & fixrho.idx0, 3]])
        }
        tau[-covariance.idx[, 1]][tau[-covariance.idx[, 1]] < tol] <- 0
        tau[covariance.idx[fixrho.idx, 1]] <- sign(tau[covariance.idx[fixrho.idx, 1]]) * sqrt(tau[covariance.idx[fixrho.idx, 2]] * tau[covariance.idx[fixrho.idx, 3]])
      }
    }

debug1<<-"point 1_5"
    cov <- as.matrix(fit$cov)

    beta <- as.numeric(fit$alpha)

    eta <- as.numeric(fit$eta) + offset
    if(verbose) {
      cat("Variance component estimates:\n")
      print(tau)
      cat("Fixed-effect coefficients:\n")
      print(beta)
    }
    mu <- exp(eta)
    mu.eta <- mu
    Y <- eta - offset + (y - mu)/mu.eta
    variance<-mu+disp*mu^2
    sqrtW <- mu.eta/sqrt(1/as.vector(weights)*variance)

    if(2*max(abs(beta - beta0)/(abs(beta) + abs(beta0) + tol), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) break
    if(max(abs(tau)) > tol^(-2)) {
      warning("Large variance estimate observed in the iterations, model not converged...", call. = FALSE)
      i <- maxiter
      break
    }
  }

  converged <- ifelse(i < maxiter, TRUE, FALSE)
  res <- y - mu
  res.var <- rep(1, n)
  for(i in 1:ng) res.var[group.idx[[i]]] <- tau[i]
  #if(!is.Matrix) fit$Sigma_i <- fit$Sigma_iX <- NULL
  #P<-Sigma_i-tcrossprod(fit$Sigma_iX,fit$Sigma_iXcov)
  if(!exists("Sigma_i")){
    n <- nrow(X)
    p <- ncol(X)
    q2 <- sum(fixtau == 0)

    diagSigma <- rep(0, n)
    for(i in 1:ng) diagSigma[group.idx[[i]]] <- tau[i]/sqrtW[group.idx[[i]]]^2

    entry_Sigma_i<-  lapply(ns, function(sid){
      if(length(sid)==1){sub_Sigma<-matrix(diagSigma[sid])
      }else{sub_Sigma <- diag(diagSigma[sid])}
      sub_kins<-lapply(kins,function(component){component[sid,sid,drop=F]})
      for(i in 1:q) sub_Sigma <- sub_Sigma + tau[i+ng]*sub_kins[[i]]
      sub_Sigma_i <- chol2inv(chol(sub_Sigma))
      rm(sub_Sigma)
      return(sub_Sigma_i)
    })
    names(entry_Sigma_i)<-names(ns)
     block_sum<-sapply(entry_Sigma_i,function(block){return(sum(as.matrix(block)))})
    # names(block_sum)<-names(entry_Sigma_i)
     Sigma_i<-bdiag(entry_Sigma_i)
     Sigma_iX <- crossprod(Sigma_i, X)
    # Sigma_iY <- as.vector(crossprod(Sigma_i, Y))
    # Sigma_iY2 <- as.vector(crossprod(Sigma_i, Y^2))
    # Sigma_icount <- as.vector(crossprod(Sigma_i, y))
    # Sigma_icount2 <- as.vector(crossprod(Sigma_i, y^2))
     XSigma_iX <- crossprod(X, Sigma_iX)
     XSigma_iX <- forceSymmetric(XSigma_iX)
     cov <- chol2inv(chol(XSigma_iX))
     Sigma_iXcov <- tcrossprod(Sigma_iX, cov)
    #
     rowSSigma_i<-as.vector(Matrix::rowSums(Sigma_i))

  }
  P<-as.matrix(Sigma_i-tcrossprod(Sigma_iX,Sigma_iXcov))
  #P1<<-as.matrix(Sigma_i-tcrossprod(fit$Sigma_iX,fit$Sigma_iXcov))
  #return(list(mu=mu,theta=tau, n.groups=ng, coefficients=beta, linear.predictors=eta, fitted.values=mu, Y=Y, X=X, P=P, residuals=res, scaled.residuals=res*as.vector(weights)/res.var, cov=cov, Sigma_iX=Sigma_iX, converged=converged,sqrtW=sqrtW,Sigma_iXcov=Sigma_iXcov,Sigma_iY=Sigma_iY,block_sum=block_sum,Sigma_icount=Sigma_icount,Sigma_iY2=Sigma_iY2,Sigma_icount2=Sigma_icount2,rowSSigma_i=rowSSigma_i,Sigma_i=Sigma_i))
  #return(list(mu=mu,theta=tau, n.groups=ng, coefficients=beta, linear.predictors=eta, fitted.values=mu, Y=Y, X=X, P=P, residuals=res, scaled.residuals=res*as.vector(weights)/res.var, cov=cov, Sigma_iX=Sigma_iX, converged=converged,sqrtW=sqrtW,Sigma_iXcov=Sigma_iXcov,Sigma_iY=Sigma_iY,block_sum=block_sum,Sigma_icount=Sigma_icount,rowSSigma_i=rowSSigma_i))
  return(list(mu=mu,theta=tau, n.groups=ng, coefficients=beta, linear.predictors=eta, fitted.values=mu, Y=Y, X=X, P=P, residuals=res, scaled.residuals=res*as.vector(weights)/res.var, cov=cov, Sigma_iX=Sigma_iX, converged=converged,sqrtW=sqrtW,rowSSigma_i=rowSSigma_i,block_sum=block_sum))
}







R_fitglmm_ai1 <- function(ns,Y, X, q, kins, ng, group.idx, W, tau, fixtau) {
  n <- nrow(X)
  p <- ncol(X)
  q2 <- sum(fixtau == 0)
  diagSigma <- rep(0, n)
  for(i in 1:ng) diagSigma[group.idx[[i]]] <- tau[i]/W[group.idx[[i]]]


  entry_Sigma_iXorY<-  lapply(ns, function(sid){
    if(length(sid)==1){sub_Sigma<-matrix(diagSigma[sid])
    }else{sub_Sigma <- diag(diagSigma[sid])}
    #sub_Sigma <- diag(diagSigma[sid])
    sub_kins<-lapply(kins,function(component){component[sid,sid,drop=F]})
    for(i in 1:q) sub_Sigma <- sub_Sigma + tau[i+ng]*sub_kins[[i]]
    sub_Sigma_i <- chol2inv(chol(sub_Sigma))
    rm(sub_Sigma)
    sub_X<-X[sid,,drop=F]
    sub_Sigma_iX <- crossprod(sub_Sigma_i, sub_X)
    sub_Y<-Y[sid]
    sub_Sigma_iY<-crossprod(sub_Sigma_i,sub_Y)
    #return(list(sigma_i=sub_Sigma_i))
    return(cbind(sub_Sigma_iY,sub_Sigma_iX))
  })

  Sigma_iX_new<-Reduce(rbind,entry_Sigma_iXorY)
  # Sigma_i_new<-bdiag(entry_Sigma_iXorY)
  # Sigma_iY<-as.matrix(crossprod(Sigma_i_new,Y))
  # Sigma_iX_new<-as.matrix(crossprod(Sigma_i_new, X))
  Sigma_iY<-Sigma_iX_new[,1,drop=F]
  Sigma_iX_new<-Sigma_iX_new[,-1,drop=F]

  entry_XSigma_iX<-lapply(ns, function(sid){
    sub_X<-X[sid,,drop=F]
    sub_XSigma_iX <- crossprod(sub_X, Sigma_iX_new[sid,,drop=F])
    return(sub_XSigma_iX)
  })
  cov_new <- chol2inv(chol(Reduce("+",entry_XSigma_iX)))

  Sigma_iXcov_new <- tcrossprod(Sigma_iX_new, cov_new)

  beta_new<-crossprod(cov_new,Reduce("+",lapply(ns, function(sid){
    sub_Y<-Y[sid]
    return(crossprod(Sigma_iX_new[sid,,drop=F], sub_Y))
  })))

  eta_new <- Y - diagSigma*(Sigma_iY-tcrossprod(Sigma_iX_new,t(beta_new)))


  if(q2 > 0) {
    idxtau <- which(fixtau == 0)
    PY <- Sigma_iY - tcrossprod(Sigma_iX_new, t(crossprod(Sigma_iXcov_new, Y)))
    wPY <- PY/W
    diagP <- ((1/diagSigma) - rowSums(Sigma_iX_new * Sigma_iXcov_new))/W
    AI <- matrix(NA, q2, q2)
    score <- rep(NA, q2)
    for(i in 1:q2) {
      if(idxtau[i] <= ng) {
        score[i] <- sum(wPY[group.idx[[idxtau[i]]]] * PY[group.idx[[idxtau[i]]]] - diagP[group.idx[[idxtau[i]]]])
        for(j in 1:i) {
          AI[i,j] <- sum(wPY[group.idx[[idxtau[i]]]] * crossprod(Sigma_i[group.idx[[idxtau[j]]], group.idx[[idxtau[i]]]], wPY[group.idx[[idxtau[j]]]])) - sum(crossprod(Sigma_iXcov[group.idx[[idxtau[i]]], ], wPY[group.idx[[idxtau[i]]]]) * crossprod(Sigma_iX[group.idx[[idxtau[j]]], ], wPY[group.idx[[idxtau[j]]]]))
          if(j != i) AI[j,i] <- AI[i,j]
        }
      } else {
        #APY<-matrix(Reduce(c,sapply(ns,function(sid){rep(sum(PY[sid]),length(sid))})),ncol=1)
        APY <- crossprod(kins[[idxtau[i]-ng]], PY)
        entry_Sigma_iAPY<-  lapply(ns, function(sid){
          if(length(sid)==1){sub_Sigma<-matrix(diagSigma[sid])
          }else{sub_Sigma <- diag(diagSigma[sid])}
          #sub_Sigma <- diag(diagSigma[sid])
          sub_kins<-lapply(kins,function(component){component[sid,sid,drop=F]})
          for(comp in 1:q) sub_Sigma <- sub_Sigma + tau[comp+ng]*sub_kins[[comp]]
          sub_Sigma_i <- chol2inv(chol(sub_Sigma))

          rm(sub_Sigma)
          #sum_Sigma_i<-sum(sub_Sigma_i)
          sum_Sigma_i<-sum(sub_Sigma_i*sub_kins[[i]])
          sub_Sigma_iAPY<-crossprod(sub_Sigma_i, APY[sid,,drop=F])
          return(cbind(sub_Sigma_iAPY,rep(sum_Sigma_i,length(sid))))
        })

        # entry_Sigma_iAPY<-  sapply(ns, function(sid){
        #   sub_Sigma_i <- Sigma_i_new[sid,sid]
        #   sub_Sigma_iAPY<-matrix(crossprod(sub_Sigma_i, APY[sid,]),ncol=1)
        #   return(sub_Sigma_iAPY)
        # })
        PAPY <- Reduce(rbind,entry_Sigma_iAPY)[,1,drop=F] - tcrossprod(Sigma_iX_new, t(crossprod(Sigma_iXcov_new, APY)))
        score[i] <- sum(Y * PAPY) - (sum(unique(Reduce(rbind,entry_Sigma_iAPY)[,2]))-sum(sapply(ns, function(sid){sum(Sigma_iX_new[sid,,drop=F]*crossprod(kins[[idxtau[i]-ng]][sid,sid],Sigma_iXcov_new[sid,,drop=F]))})))
        for(j in 1:i) {
          if(idxtau[j] <= ng) AI[i,j] <- sum(wPY[group.idx[[idxtau[j]]]] * PAPY[group.idx[[idxtau[j]]]])
          else AI[i,j] <- sum(PY * crossprod(kins[[idxtau[j]-ng]], PAPY))
          if(j != i) AI[j,i] <- AI[i,j]
        }


        # P=Sigma_i-Sigma_iX%*%t(Sigma_iXcov)
        # AI[i,i]=sum(diag(P%*%kins[[idxtau[i]-ng]]%*%P%*%kins[[idxtau[i]-ng]]))
      }
    }
    #print(AI)
    Dtau <- solve(AI, score)
    #print(Dtau)
    return(list(Dtau = Dtau, P = NULL, cov = cov_new, alpha = beta_new, eta = eta_new, Sigma_iX = Sigma_iX_new,Sigma_iXcov=Sigma_iXcov_new))
  }
  return(list(Dtau = NULL, P = NULL, cov = cov_new, alpha = beta_new, eta = eta_new, Sigma_iX = Sigma_iX_new,Sigma_iXcov=Sigma_iXcov_new))
}
