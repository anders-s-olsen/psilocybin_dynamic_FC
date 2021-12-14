permlme_function <- function(lme0.ML, lme1.ML, data = NULL, index.perm = NULL){ 
  
  
  ## ** 1- extract key quantities from input
  n.obs <- NROW(data) ## number of observations
  
  cluster <- getGroups(lme0.ML) ## to which cluster (patient) each observation belongs
  U.cluster <- levels(cluster)
  n.cluster <- length(U.cluster)
  index.cluster <- lapply(U.cluster, function(iCluster){which(cluster==iCluster)}) ## used to restaure proper ordering after tapply
  
  Y <- getResponse(lme1.ML) ## response
  name.Y <- all.vars(formula(lme1.ML))[[1]] ## name of the response variable
  X0 <- model.matrix(formula(lme0.ML),data) ## design matrix
  
  Omega0 <- getVarCov(lme0.ML, type = "marginal", individuals = levels(cluster)) ## residual variance-covariance matrix
  beta0 <- fixef(lme0.ML) ## regression coefficients
  
  ## ** 2- compute residuals
  Xbeta0 <- X0 %*% beta0
  residuals0 <- as.double(Y - X0 %*% beta0)
  
  ## ** 3- compute normalized residuals
  sqrtOmega0 <- lapply(Omega0,function(iOmega0){t(chol(iOmega0))})
  sqrtOmega0M1 <- lapply(sqrtOmega0,solve)
  
  residuals0N <- vector(length=n.obs, mode = "numeric")
  for(iCluster in 1:n.cluster){ ## iCluster <- 1
    residuals0N[index.cluster[[iCluster]]] <- sqrtOmega0M1[[iCluster]] %*% residuals0[index.cluster[[iCluster]]]
  }
  
  ## ** 4- estimate the distribution of the test statistics under the null
  # Anders edit: done elsewhere
  
  data.perm <- data
  
  ## permute residuals and fixed effects
  
  residuals0N.perm <- residuals0N[index.perm]
  
  ## rescale residuals
  for(iCluster in 1:n.cluster){ ## iCluster <- 1
    data.perm[[name.Y]][index.cluster[[iCluster]]] <- sqrtOmega0[[iCluster]] %*% residuals0N.perm[index.cluster[[iCluster]]]
  }
  ## add permuted fixed effects
  
  data.perm[[name.Y]] <- data.perm[[name.Y]] + Xbeta0[index.perm,,drop=FALSE]
  
  # compute new models using permuted observations
  lme0.permML <- try(update(lme0.ML, data = data.perm, method = "ML"), silent = TRUE)
  lme1.permML <- try(update(lme1.ML, data = data.perm, method = "ML"), silent = TRUE)
  if(inherits(lme0.permML,"try-error")||inherits(lme1.permML,"try-error")){
    LRT.stat <- NA
  }else{
    LRT.stat <- as.double(2*(logLik(lme1.permML)-logLik(lme0.permML)))
  }
  
  return(LRT.stat)
}

