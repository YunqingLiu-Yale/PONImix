library(Matrix)
library(gamlss)
library(pscl)
library(MASS)
library(bbmle)
library(VGAM)
library(maxLik)
library(tictoc)

DEsingle_est_single <- function(counts, group, parallel = FALSE, BPPARAM = bpparam()){
  ngroup<-length(levels(group))
  # Preprocessing
  counts_norm <- counts

  # Memory management
  counts_norm <- Matrix(counts_norm, sparse = TRUE)
  gc()
  geneNum <- nrow(counts_norm)
  sampleNum <- ncol(counts_norm)

  # Function of testing homogeneity of two ZINB populations
  CallDE <- function(i){

    # Memory management
    if(i %% 100 == 0)
      gc()
    results_gene <-matrix(NA,nrow = 1,ncol=(5*ngroup+1))

    for(group_id in 1:ngroup){
      counts_group <- counts_norm[i, group == levels(group)[group_id]]
      if(sum(counts_group == 0) > 0){
        if(sum(counts_group == 0) == length(counts_group)){
          theta_1 <- 1
          mu_1 <- 0
          size_1 <- 1
          prob_1 <- size_1/(size_1 + mu_1)
        }else{
          options(show.error.messages = FALSE)
          zinb_try <- try(gamlssML(counts_group, family="ZINBI"), silent=TRUE)
          options(show.error.messages = TRUE)
          if('try-error' %in% class(zinb_try) | sum(is.na(vcov(zinb_try)))>0){
            zinb_try_twice <- try(zeroinfl(formula = counts_group ~ 1 | 1, dist = "negbin"), silent=TRUE)
            if('try-error' %in% class(zinb_try_twice)){
              if(!'try-error' %in% class(zinb_try)){
                pi_zi<-c(zinb_try$nu*1*(as.integer(counts_group)==0))
                pi_nb<-(1-zinb_try$nu)*dnbinom(counts_group,size =(1/zinb_try$sigma),mu=zinb_try$mu)
                weight_EM<-as.numeric(t(pi_nb)/t(pi_zi+pi_nb))
                zinb_try_GLM<-try(glm.nb(formula = counts_group ~ 1,weights=weight_EM), silent=TRUE)
              }
              if('try-error' %in% class(zinb_try_GLM) | 'try-error' %in% class(zinb_try)){
                print("MLE of ZINB failed!");
                results_gene[1,ncol(results_gene)] <- "ZINB failed!"
                return(results_gene)
              }else{
                zinb_1 <- zinb_try_GLM
                theta_1 <- zinb_try$nu
                mu_1 <- exp(zinb_1$coefficients);names(mu_1) <- NULL
                disp_1 <- 1/zinb_1$theta;names(disp_1) <- NULL
                beta_1 <- zinb_1$coefficients;names(beta_1) <- NULL
                se_beta_1<-sqrt(vcov(zinb_1)[1,1]);names(se_beta_1) <- NULL
              }
            }else{
              zinb_1 <- zinb_try_twice
              theta_1 <- plogis(zinb_1$coefficients$zero);names(theta_1) <- NULL
              mu_1 <- exp(zinb_1$coefficients$count);names(mu_1) <- NULL
              disp_1 <- 1/zinb_1$theta;names(disp_1) <- NULL
              beta_1 <- zinb_1$coefficients$count;names(beta_1) <- NULL
              se_beta_1 <- sqrt(vcov(zinb_1)[1,1]);names(se_beta_1) <- NULL
            }
          }else{
            zinb_1 <- zinb_try
            theta_1 <- zinb_1$nu;names(theta_1) <- NULL
            mu_1 <- zinb_1$mu;names(mu_1) <- NULL
            disp_1 <- zinb_1$sigma;names(disp_1) <- NULL
            beta_1 <- zinb_1$mu.coefficients;names(beta_1) <- NULL
            se_beta_1 <- sqrt(vcov(zinb_1)[1,1]);names(se_beta_1) <- NULL
          }
        }
      }else{
        op <- options(warn=2)
        nb_try <- try(glm.nb(formula = counts_group ~ 1), silent=TRUE)
        options(op)
        if('try-error' %in% class(nb_try)){
          nb_try_twice <- try(fitdistr(counts_group, "Negative Binomial"), silent=TRUE)
          if('try-error' %in% class(nb_try_twice)){
            nb_try_again <- try(mle2(counts_group~dnbinom(mu=exp(logmu),size=1/invk), data=data.frame(counts_group), start=list(logmu=0,invk=1), method="L-BFGS-B", lower=c(logmu=-Inf,invk=1e-8)), silent=TRUE)
            if('try-error' %in% class(nb_try_again)){
              nb_try_fourth <- try(glm.nb(formula = counts_group ~ 1), silent=TRUE)
              if('try-error' %in% class(nb_try_fourth)){
                print("MLE of NB failed!");
                results_gene[1,ncol(results_gene)] <- "NB failed!"
                return(results_gene)
              }else{
                nb_1 <- nb_try_fourth
                theta_1 <- 0
                mu_1 <- exp(nb_1$coefficients);names(mu_1) <- NULL
                disp_1 <- 1/nb_1$theta;names(disp_1) <- NULL
                beta_1 <- nb_1$coefficients;names(beta_1) <- NULL
                se_beta_1<-sqrt(vcov(nb_1)[1,1]);names(se_beta_1) <- NULL
              }
            }else{
              nb_1 <- nb_try_again
              theta_1 <- 0
              mu_1 <- exp(nb_1@coef["logmu"]);names(mu_1) <- NULL
              disp_1 <- nb_1@coef["invk"];names(disp_1) <- NULL
              beta_1 <- nb_1@coef["logmu"];names(beta_1) <- NULL
              se_beta_1<-sqrt(vcov(nb_1)[1,1]);names(se_beta_1) <- NULL
            }
          }else{
            nb_1 <- nb_try_twice
            theta_1 <- 0
            mu_1 <- nb_1$estimate["mu"];names(mu_1) <- NULL
            disp_1 <- 1/nb_1$estimate["size"];names(disp_1) <- NULL
            beta_1 <- log(mu_1);names(beta_1) <- NULL
            se_beta_1<-nb_1$sd["mu"];names(se_beta_1) <- NULL
          }
        }else{
          nb_1 <- nb_try
          theta_1 <- 0
          mu_1 <- exp(nb_1$coefficients);names(mu_1) <- NULL
          disp_1 <- 1/nb_1$theta;names(disp_1) <- NULL
          beta_1 <- nb_1$coefficients;names(beta_1) <- NULL
          se_beta_1<-sqrt(vcov(nb_1)[1,1]);names(se_beta_1) <- NULL
        }
      }
      ##########################################################
      results_gene[1,(1:5)+5*(group_id-1)] <- c(theta_1,mu_1,disp_1,beta_1,se_beta_1)
    }
    colnames(results_gene)<-c(paste(rep(c("theta","mu","disp","beta","se_beta"),ngroup),rep(1:ngroup,each=5),sep="_"),"Remark")
    return(results_gene)
  }
  # Call DEG gene by gene
  if(!parallel){
    results <- matrix(data=NA, nrow = geneNum, ncol = (5*ngroup+1), dimnames = list(row.names(counts_norm), c(paste(rep(c("theta","mu","disp","beta","se_beta"),ngroup),rep(1:ngroup,each=5),sep="_"),"Remark")))
    results <- as.data.frame(results)
    for(i in 1:geneNum){
      cat("\r",paste0("DEsingle is analyzing ", i," of ",geneNum," expressed genes"))
      results[i,] <- CallDE(i)
    }
  }else{
    message("DEsingle is analyzing ", geneNum, " expressed genes in parallel")
    results <- do.call(rbind, bplapply(1:geneNum, CallDE, BPPARAM = BPPARAM))
  }
  return(results)
}


BuildModelMatrix<-function(fixed,data){
  #fixed is formula; data is pheno.data
  if(all.vars(fixed)[1]!="count"){
    #print("The response name should be #count#...automatical change")
    if(is.na(all.vars(fixed)[2])){
      fixed<-count~1
      cov.formula<-~1
    }else{
      fixed<-as.formula(paste0("count~",paste0((all.vars(fixed))[-1],collapse = "+")))
      cov.formula<-reformulate(attr(delete.response(terms(fixed)), "term.labels"))
    }
  }else{
    if(is.na(all.vars(fixed)[2])){
      cov.formula<-~1
    }else{
      cov.formula<-reformulate(attr(delete.response(terms(fixed)), "term.labels"))
    }
  }
  tmp.matrix<-model.frame(formula = cov.formula, data = data, na.action = na.omit)
  if(ncol(tmp.matrix)==0){
    modelMatrix<-matrix(rep(1,nrow(tmp.matrix)),ncol=1)
    # modelMatrix_GLM<-modelMatrix
    colnames(modelMatrix)<-"(Intercept)"
    #fix.intercept<-T
  }else if(ncol(tmp.matrix)>0){
    #modelMatrix_GLM<-cbind(rep(1,nrow(tmp.matrix)),tmp.matrix)
    #modelMatrix<-tmp.matrix
    #fix.intercept<-F
    modelMatrix<-cbind(rep(1,nrow(tmp.matrix)),tmp.matrix)
    colnames(modelMatrix)[1]<-"(Intercept)"
  }else{print("Error in columns of Model Matrix");break();}
  rownames(modelMatrix)<-seq(1:nrow(tmp.matrix))
  modelMatrix<-as.matrix(modelMatrix)
  return(modelMatrix)
}



estimateDispMuGLM <- function(y,fixed,alpha_hat,weights,sizefactor,separateGLM,
                              data = parent.frame(),useWeights=TRUE
                              , minDisp=1e-8, kappa_0=1,
                              dispTol=1e-6, maxit=100, quiet=FALSE,
                              niter=20, linearMu=NULL,
                              minmu=0) {
  n_sample<-ncol(y)
  n_gene<-nrow(y)
  maxDisp <- max(10, ncol(y))
  alpha_hat <- alpha_hat_new <- alpha_init <- pmin(pmax(minDisp, alpha_hat), maxDisp)

  stopifnot(length(niter) == 1 & niter > 0)

  if (useWeights) {
    linearMu <- FALSE
  }

  fitidx <- rep(TRUE,nrow(y))
  mu <- matrix(0, nrow=nrow(y), ncol=ncol(y))

  beta.est <- matrix(0, nrow=nrow(y),ncol = length(all.vars(fixed)))
  beta.cov <- matrix(0, nrow=nrow(y),ncol = length(all.vars(fixed)))
  #tau.est <- rep(0, nrow=nrow(y))

  dispIter <- numeric(nrow(y))

  modelMatrix<-BuildModelMatrix(fixed,data)

  for (iter in seq_len(niter)) {



    dat_weight_disp<<-cbind(y[fitidx,,drop=FALSE],weights[fitidx,,drop=FALSE],alpha_hat[fitidx])
    fit_GLM<<-apply(dat_weight_disp, 1, function(yw){
      nsample<-(length(yw)-1)/2
      yw_y<<-yw[1:nsample]
      yw_w<<-yw[(nsample+1):(length(yw)-1)]
      yw_disp<<-yw[length(yw)]
      if(length(all.vars(fixed))==1){
        fit<-tryCatch(glm.nb(formula = yw_y~offset(log(sizefactor)),weights = yw_w,init.theta = (1/yw_disp)),error=function(e){NA})
      }else if(length(all.vars(fixed))==2){
        #fit<-tryCatch(glm.nb(formula = yw_y~data$disease+offset(log(sizefactor)),weights = yw_w,init.theta = (1/yw_disp)),error=function(e){NA})
        fit<-try(glm.nb(formula = yw_y~data$disease+offset(log(sizefactor)),weights = yw_w,init.theta = (1/yw_disp)),silent = T)
        if('try-error' %in% class(fit)){
          fit_try <- try(zeroinfl(formula = yw_y~data$disease+offset(log(sizefactor)) | 1, dist = "negbin"))

          if('try-error' %in% class(fit_try)){fit<-NA}else{
            fit_try$coefficients<- fit_try$coefficients$count
            fit_try$fitted.values<-exp(fit_try$coefficients[1]+data$disease*fit_try$coefficients[2])
            fit_try$iter=1
            fit=fit_try}
        }
      }else{
        stop("Cannot handle additional covariates")
      }
      if(!is.na(fit)){
        if(length(unname(sqrt(diag(vcov(fit)))))>2){betase=unname(sqrt(diag(vcov(fit))))[1:2]}else{betase=unname(sqrt(diag(vcov(fit))))}
        return(list(mu=fit$fitted.values,betaMatrix=unname(fit$coefficients),betaSE=betase,disp=(1/(fit$theta)),niter=fit$iter))
      }else{return(list(mu=rep(NA,length(yw_y)),betaMatrix=rep(NA,length(all.vars(fixed))),betaSE=rep(NA,length(all.vars(fixed))),disp=NA,niter=NA))}
    })
    Mu<-t(sapply(fit_GLM, function(l){return(unname(l$mu))}))
    Beta<-t(sapply(fit_GLM, function(l){unname(l$betaMatrix)}))
    Cov<-t(sapply(fit_GLM, function(l){unname(l$betaSE)}))
    Disp<<-as.vector(sapply(fit_GLM, function(l){unname(l$disp)}))
    fitIter<-as.vector(sapply(fit_GLM, function(l){unname(l$niter)}))

    #print(paste("fitidx=",sum(fitidx)))


    fitMu <- Mu

    fitMu[fitMu < minmu] <- minmu
    mu[fitidx,] <- fitMu

    beta.est[fitidx,]<-Beta
    beta.cov[fitidx,]<-Cov


    dispIter[fitidx] <- fitIter
    alpha_hat_new[fitidx] <- pmin(Disp, maxDisp)

    fitidx <- abs(log(alpha_hat_new) - log(alpha_hat)) > .005
    fitidx[is.na(fitidx)]<-F
    alpha_hat <- alpha_hat_new

    if (sum(fitidx) == 0) {break()}
  }

  dispGeneEst <- alpha_hat

  # if (niter == 1) {
  #   noIncrease <- dispRes$last_lp < dispRes$initial_lp + abs(dispRes$initial_lp)/1e6
  #   dispGeneEst[which(noIncrease)] <- alpha_init[which(noIncrease)]
  # }

  # dispGeneEstConv <- dispIter < maxit & !(dispIter == 1)
  #
  # refitDisp <- !dispGeneEstConv & dispGeneEst > minDisp*10
  #
  # tic("refit5")
  # if (sum(refitDisp) > 0) {
  #   print("refit")
  #   dispGrid <- fitDispGridWrapper(y = y[refitDisp,,drop=FALSE],
  #                                  x = modelMatrix,
  #                                  mu = mu[refitDisp,,drop=FALSE],
  #                                  logAlphaPriorMean = rep(0,sum(refitDisp)),
  #                                  logAlphaPriorSigmaSq = 1, usePrior = FALSE,
  #                                  weightsSEXP = weights[refitDisp,,drop=FALSE],
  #                                  useWeightsSEXP = TRUE)
  #   dispGeneEst[refitDisp] <- dispGrid
  # }
  # toc()

  dispGeneEst <- pmin(pmax(dispGeneEst, minDisp), maxDisp)

  # return(list(dispGeneEst=dispGeneEst,
  #             dispGeneIter=dispIter,
  #             mu=mu,Beta=beta.est,Tau=tau.est,Beta_cov=beta.cov))
  return(list(dispGeneEst=dispGeneEst,
              dispGeneIter=dispIter,
              mu=mu,Beta=beta.est,Beta_cov=beta.cov))
}

mse<-function(est,true){
  rowMeans(apply(est,1,function(x){(x-true)^2}))
}
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
fitNbinomGLMs <- function(y,sizefactor, modelMatrix, alpha_hat, weight, lambda_GLM,modelAsFormula=F,
                          renameCols=TRUE, betaTol=1e-8, maxit=100, useOptim=TRUE,
                          useQR=TRUE, forceOptim=FALSE, warnNonposVar=TRUE, minmu=0.5) {

  normalizationFactors <- matrix(rep(sizefactor,each=nrow(y)),
                                 ncol=ncol(y))
  if (length(alpha_hat) != nrow(y)) {
    stop("alpha_hat needs to be the same length as nrows(object)")
  }

  # set a wide prior for all coefficients
  if (missing(lambda_GLM)) {
    lambda_GLM <- rep(1e-6, ncol(modelMatrix))
  }
  weights<- weight
  useWeights <-TRUE
  # bypass the beta fitting if the model formula is only intercept and
  # the prior variance is large (1e6)
  # i.e., LRT with reduced ~ 1 and no beta prior
  justIntercept <- if (modelAsFormula) {
    modelFormula == formula(~ 1)
  } else {
    ncol(modelMatrix) == 1 & all(modelMatrix == 1)
  }
  if (justIntercept & all(lambda_GLM <= 1e-6)) {
    alpha <- alpha_hat
    betaConv <- rep(TRUE, nrow(y))
    betaIter <- rep(1,nrow(y))
    betaMatrix <-      matrix(log2(rowSums(weights*t(t(y)/sizefactor))
                                   /rowSums(weights)),ncol=1)
    mu <- normalizationFactors * as.numeric(2^betaMatrix)
    logLikeMat <- dnbinom(y, mu=mu, size=1/alpha, log=TRUE)
    logLike <-rowSums(weights*logLikeMat)

    #modelMatrix <- stats::model.matrix.default(~ 1, as.data.frame(colData(object)))
    modelMatrix <-modelMatrix
    colnames(modelMatrix) <- modelMatrixNames <- "Intercept"
    w <-weights * (mu^-1 + alpha)^-1

    xtwx <- rowSums(w)
    sigma <- xtwx^-1
    betaSE <- matrix(log2(exp(1)) * sqrt(sigma),ncol=1)
    hat_diagonals <- w * xtwx^-1;
    res <- list(logLike = logLike, betaConv = betaConv, betaMatrix = betaMatrix,
                betaSE = betaSE, mu = mu, betaIter = betaIter,
                modelMatrix=modelMatrix,
                nterms=1, hat_diagonals=hat_diagonals)
    return(res)
  }
  modelMatrix<-as.matrix(modelMatrix)
  modelMatrixNames<-colnames(modelMatrix)
  qrx <- qr(modelMatrix)
  # if full rank, estimate initial betas for IRLS below
  if (qrx$rank == ncol(modelMatrix)) {
    Q <- qr.Q(qrx)
    R <- qr.R(qrx)
    # y <- t(log(counts(object,normalized=TRUE) + .1))
    Y<-t(log(t(t(y)/sizefactor)+.1))
    beta_mat <- t(solve(R, t(Q) %*% Y))
  } else {
    if ("(Intercept)" %in% modelMatrixNames) {
      beta_mat <- matrix(0, ncol=ncol(modelMatrix), nrow=nrow(y))
      # use the natural log as fitBeta occurs in the natural log scale
      logBaseMean <- log(rowMeans(t(t(y)/sizefactor)))
      beta_mat[,which(modelMatrixNames == "(Intercept)")] <- logBaseMean
    } else {
      beta_mat <- matrix(1, ncol=ncol(modelMatrix), nrow=nrow(y))
    }
  }

  # here we convert from the log2 scale of the betas
  # and the beta prior variance to the log scale
  # used in fitBeta.
  # so we divide by the square of the
  # conversion factor, log(2)
  lambdaNatLogScale <- lambda_GLM / log(2)^2

  betaRes <- fitBetaWrapper(ySEXP = y, xSEXP = modelMatrix,
                            nfSEXP = normalizationFactors,
                            alpha_hatSEXP = alpha_hat,
                            beta_matSEXP = beta_mat,
                            lambdaSEXP = lambdaNatLogScale,
                            weightsSEXP = weights,
                            useWeightsSEXP = useWeights,
                            tolSEXP = betaTol, maxitSEXP = maxit,
                            useQRSEXP=useQR, minmuSEXP=minmu)

  # Note on deviance: the 'deviance' calculated in fitBeta() (C++)
  # is not returned in mcols(object)$deviance. instead, we calculate
  # the log likelihood below and use -2 * logLike.
  # (reason is that we have other ways of estimating beta:
  # above intercept code, and below optim code)

  mu <- normalizationFactors * t(exp(modelMatrix %*% t(betaRes$beta_mat)))
  dispersionVector <- rep(alpha_hat, times=ncol(y))
  logLike <- nbinomLogLike(y, mu, alpha_hat, weights, useWeights)

  # test for stability
  rowStable <- apply(betaRes$beta_mat,1,function(row) sum(is.na(row))) == 0

  # test for positive variances
  rowVarPositive <- apply(betaRes$beta_var_mat,1,function(row) sum(row <= 0)) == 0

  # test for convergence, stability and positive variances
  betaConv <- betaRes$iter < maxit

  # here we transform the betaMatrix and betaSE to a log2 scale
  betaMatrix <- log2(exp(1))*betaRes$beta_mat
  colnames(betaMatrix) <- modelMatrixNames
  colnames(modelMatrix) <- modelMatrixNames
  # warn below regarding these rows with negative variance
  betaSE <- log2(exp(1))*sqrt(pmax(betaRes$beta_var_mat,0))
  colnames(betaSE) <- paste0("SE_",modelMatrixNames)

  # switch based on whether we should also use optim
  # on rows which did not converge
  rowsForOptim <- if (useOptim) {
    which(!betaConv | !rowStable | !rowVarPositive)
  } else {
    which(!rowStable | !rowVarPositive)
  }

  if (forceOptim) {
    rowsForOptim <- seq_along(betaConv)
  }

  if (length(rowsForOptim) > 0) {
    # we use optim if didn't reach convergence with the IRLS code
    resOptim <- fitNbinomGLMsOptim(object,modelMatrix,lambda_GLM,
                                   rowsForOptim,rowStable,
                                   normalizationFactors,alpha_hat,
                                   weights,useWeights,
                                   betaMatrix,betaSE,betaConv,
                                   beta_mat,
                                   mu,logLike,minmu=minmu)
    betaMatrix <- resOptim$betaMatrix
    betaSE <- resOptim$betaSE
    betaConv <- resOptim$betaConv
    mu <- resOptim$mu
    logLike <- resOptim$logLike
  }

  stopifnot(!any(is.na(betaSE)))
  nNonposVar <- sum(rowSums(betaSE == 0) > 0)
  if (warnNonposVar & nNonposVar > 0) warning(nNonposVar,"rows had non-positive estimates of variance for coefficients")

  list(logLike = logLike, betaConv = betaConv, betaMatrix = betaMatrix,
       betaSE = betaSE, mu = mu, betaIter = betaRes$iter, modelMatrix=modelMatrix,
       nterms=ncol(modelMatrix), hat_diagonals=betaRes$hat_diagonals)
}

fitBetaWrapper <- function (ySEXP, xSEXP, nfSEXP, alpha_hatSEXP, contrastSEXP,
                            beta_matSEXP, lambdaSEXP, weightsSEXP, useWeightsSEXP,
                            tolSEXP, maxitSEXP, useQRSEXP, minmuSEXP) {
  if ( missing(contrastSEXP) ) {
    # contrast is not required, just give 1,0,0,...
    contrastSEXP <- c(1,rep(0,ncol(xSEXP)-1))
  }
  # test for any NAs in arguments
  arg.names <- names(formals(fitBetaWrapper))
  na.test <- sapply(mget(arg.names), function(x) any(is.na(x)))
  if (any(na.test)) stop(paste("in call to fitBeta, the following arguments contain NA:",
                               paste(arg.names[na.test],collapse=", ")))

  fitBeta(ySEXP=ySEXP, xSEXP=xSEXP, nfSEXP=nfSEXP, alpha_hatSEXP=alpha_hatSEXP,
          contrastSEXP=contrastSEXP, beta_matSEXP=beta_matSEXP,
          lambdaSEXP=lambdaSEXP, weightsSEXP=weightsSEXP, useWeightsSEXP=useWeightsSEXP,
          tolSEXP=tolSEXP, maxitSEXP=maxitSEXP, useQRSEXP=useQRSEXP, minmuSEXP=minmuSEXP)
}
nbinomLogLike <- function(counts, mu, disp, weights, useWeights) {
  if (is.null(disp)) return(NULL)
  if (useWeights) {
    rowSums(weights * matrix(dnbinom(counts,mu=mu,size=1/disp,
                                     log=TRUE),ncol=ncol(counts)))
  } else {
    rowSums(matrix(dnbinom(counts,mu=mu,size=1/disp,
                           log=TRUE),ncol=ncol(counts)))
  }
}



fitDispWrapper <- function (ySEXP, xSEXP, mu_hatSEXP, log_alphaSEXP, log_alpha_prior_meanSEXP,
                            log_alpha_prior_sigmasqSEXP, min_log_alphaSEXP, kappa_0SEXP,
                            tolSEXP, maxitSEXP, usePriorSEXP, weightsSEXP, useWeightsSEXP) {
  # test for any NAs in arguments
  arg.names <- names(formals(fitDispWrapper))
  na.test <- sapply(mget(arg.names), function(x) any(is.na(x)))
  if (any(na.test)) stop(paste("in call to fitDisp, the following arguments contain NA:",
                               paste(arg.names[na.test],collapse=", ")))
  fitDisp(ySEXP=ySEXP, xSEXP=xSEXP, mu_hatSEXP=mu_hatSEXP,
          log_alphaSEXP=log_alphaSEXP, log_alpha_prior_meanSEXP=log_alpha_prior_meanSEXP,
          log_alpha_prior_sigmasqSEXP=log_alpha_prior_sigmasqSEXP,
          min_log_alphaSEXP=min_log_alphaSEXP, kappa_0SEXP=kappa_0SEXP,
          tolSEXP=tolSEXP, maxitSEXP=maxitSEXP, usePriorSEXP=usePriorSEXP,
          weightsSEXP=weightsSEXP, useWeightsSEXP=useWeightsSEXP)
}

fitDispGridWrapper <- function(y, x, mu, logAlphaPriorMean, logAlphaPriorSigmaSq, usePrior,
                               weightsSEXP, useWeightsSEXP) {
  # test for any NAs in arguments
  arg.names <- names(formals(fitDispGridWrapper))
  na.test <- sapply(mget(arg.names), function(x) any(is.na(x)))
  if (any(na.test)) stop(paste("in call to fitDispGridWrapper, the following arguments contain NA:",
                               paste(arg.names[na.test],collapse=", ")))
  minLogAlpha <- log(1e-8)
  maxLogAlpha <- log(max(10, ncol(y)))
  dispGrid <- seq(from=minLogAlpha, to=maxLogAlpha, length=20)
  logAlpha <- fitDispGrid(ySEXP=y, xSEXP=x, mu_hatSEXP=mu, disp_gridSEXP=dispGrid,
                          log_alpha_prior_meanSEXP=logAlphaPriorMean,
                          log_alpha_prior_sigmasqSEXP=logAlphaPriorSigmaSq,
                          usePriorSEXP=usePrior,
                          weightsSEXP=weightsSEXP, useWeightsSEXP=useWeightsSEXP)$log_alpha
  exp(logAlpha)
}


EMcell_in_gene<-function(count.dat,depth,pheno.data,fixed=y~1,output_mu_matrix=F,separateGLM=F,initial=c(0.5,0.1,0.5),cutoff=FALSE,tolerance=rep(1e-3,4),max_iteration){
  if(class(count.dat) == "numeric"){
    count.dat<-matrix(count.dat,nrow = 1)
  }

  conv_iter=1
  n_sample<-ncol(count.dat)
  n_gene<-nrow(count.dat)


  if (missing(pheno.data)) {
    pheno.data<-data.frame(matrix(NA,nrow=n_sample,ncol=1))
  }

  if(!class(fixed)=="formula"){
    stop("Error: \"fixed\" must be a \"formula\".")
  }

  if(is.na(all.vars(fixed)[2])){
    if (missing(pheno.data)) {
      pheno.data<-data.frame(matrix(NA,nrow=n_sample,ncol=1))
    }
  }else{
    if (missing(pheno.data)) {
      stop("Error: Please input phenotype data \"pheno.data\".")
    }
    if(!(colnames(pheno.data) %in% all.vars(fixed))){
      stop("Error: Variables in formula must be in \"pheno.data\".")
    }
  }

  d<-depth

  #************************
  # depth<-colSums(count.dat)
  # d<-depth/median(depth)
  #************************

  #number of covariates
  p<-length(attr(delete.response(terms(fixed)), "term.labels"))



  if(class(initial) == "numeric"){
    initial<-matrix(initial,nrow = 1)
  }
  if(n_gene==nrow(initial)){
    pi0<-initial[,1]
    lambda0<-initial[,2]
    disp0<-initial[,3]
  }else{
    if(nrow(initial)!=1){print("Default: Using the first row of provided initial value")}
    pi0<-rep(initial[1,1],n_gene)
    lambda0<-rep(initial[1,2],n_gene)
    disp0<-rep(initial[1,3],n_gene)
  }


  #mu0<-rep(mean(dat[dat!=0]),n_gene)
  #mu0<-colMeans(dat[dat!=0])

  ###########################Let the mu become a matrix and update

  mu0<-matrix(apply(count.dat,1,function(x){mean(x[x!=0])}),ncol = n_sample,nrow = n_gene)

  #should be list of four vectors
  phi_new<-cbind(pi0,lambda0,disp0,mu0)
  phi_old<-matrix(0.1,nrow=n_gene,ncol=(3+n_sample))
  beta_est<-matrix(NA,nrow=n_gene,ncol=(p+1))
  #tau_est<-rep(NA,n_gene)
  beta_cov<-matrix(NA,nrow=n_gene,ncol=(p+1))
  #while (sum((phi_new-phi_old)^2)>1e-3)
  converge_idx <- rep(TRUE,nrow(count.dat))
  weights_output<-matrix(NA,nrow=n_gene,ncol=n_sample)
  while (sum(converge_idx) != 0) {

    dat<-count.dat[converge_idx,,drop=F]
    phi_old<-phi_new

    pi<-phi_old[converge_idx,1,drop=F]
    lambda<-phi_old[converge_idx,2,drop=F]
    disp<-phi_old[converge_idx,3,drop=F]
    mu<-phi_old[converge_idx,-(1:3),drop=F]

    r<-1/disp

    pi_poi<-apply(cbind(pi,lambda,dat), 1, function(x){x[1]*dpois(as.integer(x[-c(1,2)]),lambda = x[2])})
    pi_nb<-apply(cbind(pi,r,mu,dat), 1, function(x){(1-x[1])*dnbinom(x[(length(x)-n_sample+1):length(x)],size = x[2],mu=x[3:(n_sample+2)])})
    #E-step

    z1<-t(pi_poi)/t(pi_poi+pi_nb)
    weights_nb<-1-z1
    #M-step

    pi_new<-rowSums(z1)/n_sample

    lambda_tmp<-rowSums(z1*dat)/rowSums(z1)
    lambda_tmp[is.nan(lambda_tmp)]=0
    if(cutoff!=0){
      update_idx<-lambda_tmp<=cutoff
      lambda_new<-lambda
      lambda_new[update_idx]<-lambda_tmp[update_idx]
    }else{
      lambda_new<-lambda_tmp
    }

    result<-estimateDispMuGLM(y = dat,fixed=fixed,alpha_hat=disp,weights = weights_nb,sizefactor=d,data=pheno.data,separateGLM)

    disp_new<-result$dispGeneEst
    mu_new<-result$mu

    #if(zero_range((mu_norm[1,]/d))){mu_new<-(mu_norm[,1]/d[1])}else{print("Mu error");break();}
    phi_new[converge_idx,]<-cbind(pi_new,lambda_new,disp_new,mu_new)
    weights_output[converge_idx,]<-weights_nb
    beta_est[converge_idx,]<-result$Beta
    #tau_est[converge_idx]<-result$Tau
    beta_cov[converge_idx,]<-result$Beta_cov
    converge_idx <- (abs(phi_old[,1]-phi_new[,1])>tolerance[1])|(abs(phi_old[,2]-phi_new[,2])>tolerance[2])|(abs(phi_old[,3]-phi_new[,3])>tolerance[3])|(abs(phi_old[,4]-phi_new[,4])>tolerance[4])

    if(sum(rowSums(is.na(phi_new))>=(1+ncol(mu_new)))!=0){converge_idx[rowSums(is.na(phi_new))==(1+ncol(mu_new))]<-F}

    conv_iter=conv_iter+1
    #print(paste0("conv_iter=",conv_iter,",cong=",sum(converge_idx)))
    cat("\r",paste0("Iterations: ", conv_iter," Converged Genes: ",sum(!converge_idx)))
    converge_status<-!converge_idx
    show_conv_id=0

    if(conv_iter>=max_iteration){show_conv_id=which(converge_idx);print("Reach max iteration");break();}
  }
  #beta_cov from GLM is already square roots
  pvalue<-2*pnorm(abs(beta_est)/(beta_cov),lower.tail = F)
  #return(list(est=phi_new,niter=conv_iter,beta=beta_est,tau=tau_est,beta_cov=beta_cov,pvalue=pvalue,converge=converge_status))
  if(output_mu_matrix){
    return(list(est=phi_new,niter=conv_iter,beta=beta_est,beta_cov=beta_cov,pvalue=pvalue,weights=weights_output,converge=converge_status,conv_id=show_conv_id))
  }else{
    mu_new_single<-apply(phi_new[,-c(1:3)], 1, function(x){getmode(x/d)})
    phi_new_single<-cbind(phi_new[,c(1:3)],mu_new_single)
    return(list(est=phi_new_single,niter=conv_iter,beta=beta_est,beta_cov=beta_cov,pvalue=pvalue,weights=weights_output,converge=converge_status,conv_id=show_conv_id))
  }

}


library(Rcpp)
library(tictoc)
library(Matrix)

# source("modifiedDeseq2_v3.R")
# source("modifiedGMMAT.R")
# sourceCpp("fitglmm.cpp")
switch.sub<-function(G,subject,seed){
  set.seed(seed)
  tmp.sub<-subject
  perm.d<-G
  sub.disease<-paste(G,tmp.sub,sep = "_")
  num.d.sub<-sum(substring(names(table(sub.disease)),1,1)==1)
  num.c.sub<-length(table(sub.disease))-num.d.sub
  orig.label<-sample(names(table(sub.disease)),num.c.sub)

  switch.id<-(tmp.sub%in%(sub(".*_","",orig.label)))
  perm.d[switch.id]<-0
  perm.d[!(switch.id)]<-1
  pm.subinfo<-as.numeric(factor(tmp.sub))
  pm.pdata<-data.frame(disease=1*(perm.d==1))
  return(list(pdata=pm.pdata,subinfo=pm.subinfo))
}

scoretest<-function(x,mu,P,G,offset){
  Y <- log(mu) - offset + (x - mu)/mu
  score<-t(G)%*%P%*%Y
  fishervar=t(G)%*%P%*%G
  score.stat<-(score^2)/fishervar
  return(score.stat)
}

sep_GLMM_wald<-function(count_data,pheno.data,sizefactor,sub_info,separateGLM = T,internal_mat=T,
                        niter=20,cutoff=FALSE,nperm=0,tolerance=rep(1e-3,4),
                        output_mu_matrix=T,useWeights=TRUE
                        , minDisp=1e-8, kappa_0=1,
                        dispTol=1e-6, maxit=100,sep_sigma=F, quiet=FALSE,
                        linearMu=NULL,start.seed=NULL,
                        minmu=0,max_iter_per_gene=NULL,initial=NULL,kin.test=NULL,ref=NULL,score=F){

  tic("total")

  if(class(count_data) == "numeric"){
    count_data<-matrix(count_data,nrow = 1)
  }
  if (missing(sub_info)) {
    print("Default: All cells belong to one subject")
    sub_info=rep(1,n_sample)
    names(sub_info)<-rep(1,n_sample)
  }
  if(!is.numeric(sub_info)){tmp<-sub_info;sub_info<-as.numeric(factor(sub_info,levels=unique(sub_info)));names(sub_info)<-tmp;rm(tmp)}else{names(sub_info)<-sub_info}
  # else{
  #   count_data<-count_data[,order(sub_info),drop=F]
  #   sub_info<-sort(as.numeric(factor(sub_info)))
  # }

  reorder_subject<-order(sub_info)
  sub_info<-sort(sub_info)
  count_data<-count_data[,reorder_subject,drop=F]
  sizefactor<-sizefactor[reorder_subject]
  pheno.data<-pheno.data[reorder_subject,,drop=F]

  maxDisp <- max(10, ncol(count_data))

  n_sample<-ncol(count_data)
  n_gene<-nrow(count_data)
  score_stat<-rep(NA,n_gene)
  if(missing(max_iter_per_gene)){max_iteration=1000}else{max_iteration=max_iter_per_gene}
  if(missing(initial)){
    pi_tmp<-apply(count_data, 1, function(gene){sum(gene==0)/n_sample})

    lambda_tmp<-rep(0.1,n_gene)
    disp_tmp<-apply(count_data, 1, function(gene){abs((var(gene[gene>0])-mean(gene[gene>0]))/(mean(gene[gene>0])^2))})
    disp_tmp[is.na(disp_tmp)|disp_tmp<0]<-0
    initial<-cbind(pi_tmp,lambda_tmp,disp_tmp)
  }

  n_sub<-length(unique(sub_info))#number of subject
  sub_length=table(sub_info)
  if(sum(sub_length)!=n_sample){print("Sub_info error");break();}
  G_score<-as.numeric(pheno.data[,ncol(pheno.data)])
  if(length(table(G_score))>2)stop("Error: Test are only applicable to two groups comparison.")
  if(score){group<-factor(rep(1,length(pheno.data$disease)))}else{group<-factor(pheno.data$disease)  }
  ngroup<-length(levels(group))

  weights<- matrix(NA, nrow=n_gene, ncol=n_sample)
  mu <- matrix(0, nrow=n_gene, ncol=n_sample)
  beta.est <- matrix(0, nrow=n_gene,ncol = ngroup)
  beta.cov <- matrix(0, nrow=n_gene,ncol = ngroup)
  tau.est <- matrix(0, nrow=n_gene,ncol =ngroup)

  #tau.est <- rep(0, nrow=nrow(count_data))
  #score_stat<-matrix(NA,nrow=n_gene,ncol = ngroup)
  if(missing(ref)){
    ref=levels(group)[1]
  }else{stopifnot(ref%in%levels(group),print("Error: reference level is not in the dataset."))}
  internal_matrix<-list()







  for(group_idx in 1:ngroup){
    print(paste("The group ",group_idx," is ",levels(group)[group_idx]))
    fixed<-count~1
    cov.formula<-~1
    counts_group <- count_data[, group == levels(group)[group_idx]]
    group_id<-which(group == levels(group)[group_idx])

    #sub_info_group<-(as.numeric(factor(sub_info[group_id])))
    sub_info_group<-names(sub_info)[group_id]
    nsub_group<-length(unique(sub_info_group))
    nsample_group<-length(group_id)
    sub_length_group=table(sub_info_group)[unique(sub_info_group)]

    #    if(missing(kin.test)){
    print("Build up random effect matrix")


    randomMatrix<-matrix(0,nrow=nsample_group,ncol=nsub_group)
    for(i in 1:ncol(randomMatrix)){
      randomMatrix[,i]<-1*(sub_info_group==unique(sub_info_group)[i])
    }

    kin.test<-randomMatrix%*%t(randomMatrix)
    rownames(kin.test)=colnames(kin.test)=as.character(seq(1,nsample_group))
    #    }
    ns_group<-list()
    for(j in 1:nsub_group){
      sub<-unique(sub_info_group)[j]
      sub_id<-which(sub_info_group==sub)
      ns_group[[j]]<-sub_id
    }
    names(ns_group)<-unique(sub_info_group)
    tmp.matrix<-model.frame(formula = ~1, data = group[group_id], na.action = na.omit)
    modelMatrix<-matrix(rep(1,nrow(tmp.matrix)),ncol=1)
    colnames(modelMatrix)<-"(Intercept)"
    rownames(modelMatrix)<-seq(1:nrow(tmp.matrix))
    modelMatrix<-as.matrix(modelMatrix)
    fitidx <- rep(TRUE,nrow(counts_group))
    offset_score <- log(sizefactor[group_id])

    tic(paste0("Estimate GLM ",group_idx))
    #Need to inform which group is baseline
    GLM_res<-DEsingle_est_single(counts_group,factor(as.character(group[group_id])))
    toc()
    GLM_mu<-GLM_res$mu_1
    GLM_beta=GLM_res$beta_1
    GLM_disp<-GLM_res$disp_1
    GLM_est<-cbind(GLM_res$theta_1,GLM_res$disp_1,GLM_res$mu_1,GLM_res$beta_1,GLM_res$se_beta_1)
    colnames(GLM_est)<-paste(c("pi","disp","mu","beta","betaSE"),rep(group_idx,5),sep="_")
    rownames(GLM_est)<-rownames(counts_group)
    if(group_idx==1){summary_est<-GLM_est}else{summary_est<-cbind(summary_est,GLM_est)}


    pi_zi<-apply(cbind(GLM_res$theta_1,counts_group), 1, function(x){x[1]*1*(as.integer(x[-1])==0)})
    pi_nb<-apply(cbind(GLM_res$theta_1,1/GLM_res$disp_1,GLM_res$mu_1,counts_group), 1, function(x){(1-x[1])*dnbinom(x[-c(1:3)],size =x[2],mu=x[3])})
    weight_EM<-t(pi_nb)/t(pi_zi+pi_nb)

    # if(score){
    #   if(group_idx==1){
    #     tmp.GLM_mu<-GLM_mu
    #     tmp.GLM_beta=GLM_res$beta_1
    #     tmp.GLM_disp<-GLM_res$disp_1
    #   }else{
    #     tmp.GLM_mu<-cbind(tmp.GLM_mu,GLM_mu)
    #     tmp.GLM_beta=GLM_res$beta_1
    #     tmp.GLM_disp<-GLM_res$disp_1
    #   }
    #
    # }



      Mu<-matrix(NA,nrow=sum(fitidx),ncol = nsample_group)
      Beta<-matrix(NA,nrow=sum(fitidx),ncol = ncol(modelMatrix))
      if(class(kin.test)=="list"){
        Tau<-matrix(NA,nrow=sum(fitidx),ncol = (1+length(kin.test)))
      }else{
        Tau<-matrix(NA,nrow=sum(fitidx),ncol = 2)
      }

      Cov<-matrix(NA,nrow=sum(fitidx),ncol = ncol(modelMatrix))
      converge_glmm<-rep(NA,sum(fitidx))
      internal_matrix_group<-list()
      perm_score<-list()
      for(g in 1:sum(fitidx)){
        count=counts_group[fitidx,,drop=FALSE][g,]
        glmm_dat<-data.frame(count,modelMatrix,id=seq(1,nsample_group))
        colnames(glmm_dat)<-c("count",colnames(modelMatrix),"id")
        if((nsample_group>1) && is.vector(GLM_mu)){
          if(length(GLM_mu)!=n_gene){stop("Error: the number of gene error")}
          mu_start<-matrix(GLM_mu,nrow = n_gene,ncol=nsample_group)
        }else{mu_start<-GLM_mu}
        fit0<-list(betaMatrix=matrix(GLM_beta[fitidx],ncol=1)[g,],mu=mu_start[fitidx,,drop=FALSE][g,])
        disp=GLM_disp[fitidx][g]
        weight=weight_EM[fitidx,,drop=FALSE][g,]
        tic("GLMM model fitting")
        #tryCatch( {


        model.test<-glmmkin(sizefactor = sizefactor[group_id],fit0=fit0,modelMatrix=modelMatrix,weights = weight, fixed = fixed,disp = disp,ns=ns_group,data=glmm_dat,kins=kin.test,id="id")
        converge_glmm[g]<-model.test$converged
        #},error=function(e){converge_glmm[g]=FALSE;print(paste0("Cannot analyze gene",g))})
        if(is.na(converge_glmm[g])){converge_glmm[g]=FALSE}
        toc()

        if(converge_glmm[g]){


          #Tau[g]<-model.test$theta[2]
          Tau[g,]<-model.test$theta
          Mu[g,]<-model.test$mu
          Beta[g,]<-model.test$coefficients
          Cov[g,]<-diag(model.test$cov)
          P<-model.test$P
          scale.res<-model.test$scaled.residuals
          eta=model.test$linear.predictors





          #Tau[g]<-model.test$theta[2]


          if(score){
            score<-Matrix::crossprod(G_score,scale.res/(1+disp*model.test$mu))
            var<-diag(Matrix::crossprod(G_score,Matrix::crossprod(P,G_score)))
            score_stat[g]<-score^2/var

            if(nperm & score){
              tic("Permutation")
              if(missing(start.seed)){start.seed=1}
              seed.perm<-start.seed:(start.seed+nperm-1)
              perm_score[[g]]<-unlist(lapply(seed.perm,function(seed){
                G_perm<-(switch.sub(G_score,sub_info,seed)$pdata)[,1]
                score<-Matrix::crossprod(G_perm,scale.res/(1+disp*model.test$mu))
                var<-diag(Matrix::crossprod(G_perm,Matrix::crossprod(P,G_perm)))
                return(score^2/var)
              }))
              toc()
            }
          }

          # internal_matrix_group[[g]]<-list(model.test$Sigma_iX,model.test$Sigma_iXcov,model.test$sqrtW,scale.res,model.test$mu,model.test$linear.predictors,model.test$Y,model.test$Sigma_iY,model.test$block_sum,model.test$Sigma_icount,weight,fit0$mu,model.test$Sigma_icount2,model.test$Sigma_iY2,model.test$rowSSigma_i,model.test$Sigma_i)
          # names(internal_matrix_group[[g]])<-c("Sigma_iX","Sigma_iXcov","sqrtW","scale.res","mu_glmm","eta_glmm","Y","Sigma_iY","block_sum","Sigma_icount","weight_dropout","mu_glm","Sigma_icount2","Sigma_iY2","rowSSigma_i","Sigma_i")
          internal_matrix_group[[g]]<-list(model.test$Sigma_iX,model.test$Sigma_iXcov,model.test$sqrtW,scale.res,model.test$mu,model.test$linear.predictors,model.test$Y,model.test$Sigma_iY,model.test$block_sum,model.test$Sigma_icount,weight,fit0$mu,model.test$rowSSigma_i)
          names(internal_matrix_group[[g]])<-c("Sigma_iX","Sigma_iXcov","sqrtW","scale.res","mu_glmm","eta_glmm","Y","Sigma_iY","block_sum","Sigma_icount","weight_dropout","mu_glm","rowSSigma_i")

        }else{
          print(paste0("Cannot analyze stat",g))
        }
      }
      internal_matrix[[group_idx]]<-internal_matrix_group
      weights[,group_id]<-weight_EM
      alpha_hat <- alpha_hat_new <- alpha_init <- pmin(pmax(minDisp, GLM_disp), maxDisp)
      fitMu <- Mu

      fitMu[fitMu < minmu] <- minmu
      fitMu1<<-fitMu
      mu1<<-mu
      mu[fitidx,group_id] <- fitMu
      beta.est[fitidx,group_idx]<-Beta
      beta.cov[fitidx,group_idx]<-Cov
      tau.est[fitidx,group_idx]<-Tau[,2]



  }

  if(!score){
    ref_idx<-which(levels(group)==ref)
    comp_idx<-which(levels(group)!=ref)
    glmm_beta<-cbind(beta.est[,ref_idx],beta.est[,comp_idx]-beta.est[,ref_idx])
    glmm_betaSE<-cbind(sqrt(beta.cov[,ref_idx]),sqrt(rowSums(beta.cov)))
    glmm_wald<-glmm_beta^2/glmm_betaSE^2
    glm_beta=cbind(summary_est[,"beta_1"],summary_est[,"beta_2"]-summary_est[,"beta_1"])
    glm_betaSE<-cbind(summary_est[,"betaSE_1"],sqrt(rowSums(cbind(summary_est[,"betaSE_1"]^2,summary_est[,"betaSE_2"]^2))))
    glm_wald<-glm_beta^2/glm_betaSE^2
    glmm_score<-rep(NA,length(glmm_wald))
    perm_score<-list(NA)
  }else{
    glmm_beta<-beta.est
    glmm_betaSE<-sqrt(beta.cov)
    glmm_wald<-glmm_beta^2/glmm_betaSE^2
    glm_beta=summary_est[,"beta_1"]
    glm_betaSE<-summary_est[,"betaSE_1"]
    glm_wald<-glm_beta^2/glm_betaSE^2
    glmm_score<-score_stat
  }


  toc()
  return(list(Beta=glmm_beta,Tau=tau.est,Beta_SE=glmm_betaSE,GLMM_wald=glmm_wald,GLMM_score=glmm_score,GLM_est=summary_est,GLM_wald=glm_wald,covergence=converge_glmm,perm_score_test=perm_score,internal_matrix=internal_matrix))
}

