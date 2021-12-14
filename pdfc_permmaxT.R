list.of.packages <- c("nlme","parallel","tictoc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library("nlme")
library("parallel")
library("tictoc")

source("pdfc_permlme_function.R")
dt <- read.csv("perm_in.csv")

nperm = 100000


# initialize output dataframe and corresponding counter
outdf <- data.frame(matrix(ncol = 8, nrow = sum(unique(dt$N_centroids))))
colnames(outdf) <- c("k","state","coefintercept","coefcov","coefCIcovlow","coefCIcovhigh","pval","pval_perm")
counter = 1


# Function that can be parallelized in x = permutations.
# The function permutes observations and runs the permlme_function
# the latter being an implementation of Lee2012 (Biometrics).
# dtk is the subset of data for a given k
LRTapply <- function(x,dtk,n.obs,e.lmeH0,e.lmeH1){
  
  index.perm <- sample(1:n.obs)
  LRT = 1:k
  for (state in 1:k){
    # make data subset specific to this state
    dtkc = dtk[dtk$Current_centroid==state,]
    LRT[state] <- permlme_function(e.lmeH0[[state]], e.lmeH1[[state]], data = dtkc, index.perm = index.perm)
  }
  return(LRT)
}


tic()
for (k in unique(dt$N_centroids)){
  
  # Make data subset specific to this k
  dtk <- dt[dt$N_centroids==k,]
  
  # Initialize variables because R
  e.lmeH0 <- as.list(1:k)
  e.lmeH1 <- as.list(1:k)
  LRT_init <- as.list(1:k)
  pval_init = as.list(1:k)
  perm_p = as.list(1:k)
  
  
  # Compute initial models
  for (state in 1:k){
    
    # make data subset specific to this state
    dtkc = dtk[dtk$Current_centroid==state,] 
    
    n.obs = nrow(dtkc)
    
    # Compute initial model
    e.lmeH0[[state]] <- nlme::lme(Frac_occ ~ 1, random =~1|Subject, data = dtkc,method="ML")
    e.lmeH1[[state]] <- nlme::lme(Frac_occ ~ 1 + cov1, random=~1|Subject, data = dtkc,method="ML")
    
    LRT_init[[state]] <- as.double(2*(logLik(e.lmeH1[[state]])-logLik(e.lmeH0[[state]])))
    pval_init[[state]] = pchisq(LRT_init[[state]],df=1,lower.tail = FALSE)
  }
  
  ######################### Run permutation testing for all states ################################
  LRT_list = mclapply(1:nperm,LRTapply,dtk=dtk,n.obs=n.obs,e.lmeH0=e.lmeH0,e.lmeH1=e.lmeH1,mc.cores=20)
  
  LRT_matrix = matrix(unlist(LRT_list),ncol=k,byrow=TRUE)
  
  # Find max statistic per permutation
  LRTmax <- apply(LRT_matrix,1,max)
  LRTmax = LRTmax[!is.na(LRTmax)]
  # find permutation p-value and paste information into output dataframe
  for (state in 1:k){
    
    perm_p[state] = sum(LRTmax>=LRT_init[[state]])/nperm
    
    outdf$k[counter] = k
    outdf$state[counter] = state
    outdf$coefintercept[counter]  = e.lmeH1[[state]]["coefficients"]$coefficients$fixed["(Intercept)"]
    outdf$coefcov[counter]  = e.lmeH1[[state]]["coefficients"]$coefficients$fixed["cov1"]
    
    # confidence intervals for lme fails when fit is singular
    coefci <- try(intervals(e.lmeH1[[state]]), silent = TRUE)
    if(inherits(coefci,"try-error")){
      outdf$coefCIcovlow[counter]  = NA
      outdf$coefCIcovhigh[counter]  = NA
    }else{
      outdf$coefCIcovlow[counter]  = coefci$fixed[2,1]
      outdf$coefCIcovhigh[counter]  = coefci$fixed[2,3]
    }
    
    outdf$pval[counter]  = pval_init[[state]]
    outdf$pval_perm[counter]  = perm_p[[state]]
    
    counter = counter + 1
  }
  print(paste0("Done with model k=",k," in"))
  toc()
}

write.csv(outdf,file="perm_out.csv",na="NaN")
toc()
