

### main classtpx function when a known theta matrix is provided by the user

class_topics_2 <- function(counts, 
                         K, 
                         theta_known, 
                         class_labs=NULL, 
                         method=c("omega.fix", "theta.fix", "theta.prior","no.fix"),
                         optional_theta=NULL,
                         shrink=TRUE,
                         shrink.method=c(1,2),
                         mash_user=NULL,
                         shape=NULL, 
                         initopics=NULL, 
                         tol=0.1, 
                         bf=FALSE, 
                         kill=2, 
                         ord=TRUE, verb=1, 
                         tmax=10000, wtol=10^(-4), 
                         qn=100, grp=NULL, admix=TRUE, 
                         prior_omega = NULL,
                         nonzero=FALSE, dcut=-10){
  
  
  K_classes <- dim(theta_known)[2]
  if(is.null(prior_omega)){
    prior_omega <- rep(1/K, K);
  }else{
    prior_omega <- (prior_omega/sum(prior_omega))
  }
  
  
}