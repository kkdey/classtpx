#######  Undocumented "class.tpx" utility functions #########

## ** Only referenced from topics.R

## check counts (can be an object from tm, slam, or a simple co-occurance matrix)
CheckCounts <- function(counts){
  if(class(counts)[1] == "TermDocumentMatrix"){ counts <- t(counts) }
  if(is.null(dimnames(counts)[[1]])){ dimnames(counts)[[1]] <- paste("doc",1:nrow(counts)) }
  if(is.null(dimnames(counts)[[2]])){ dimnames(counts)[[2]] <- paste("wrd",1:ncol(counts)) }
  empty <- row_sums(counts) == 0
  if(sum(which(empty==TRUE))){
    counts <- counts[!empty,]
    cat(paste("Removed", sum(empty), "blank documents.\n")) }
  return(as.simple_triplet_matrix(counts))
}
 
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


## Topic estimation and selection for a list of K values
class.tpxSelect <- function(X, K, 
                            known_samples, omega_known, theta_known,
                            initheta, alpha, tol, kill, verb,
                      admix=TRUE, prior_omega, grp=NULL, tmax=10000,
                      wtol=10^{-4}, qn=100, nonzero=FALSE, dcut=-10){

  ## check grp if simple mixture
  if(!admix){
    if(is.null(grp) || length(grp)!=nrow(X)){  grp <- rep(1,nrow(X)) }
    else{ grp <- factor(grp) }
  }

  ## return fit for single K
  if(length(K)==1 && bf==FALSE){
    if(verb){ cat(paste("Fitting the",K,"topic model.\n")) }
    fit <-  class.tpxfit(X=X, known_samples=known_samples, 
                         omega_known=omega_known, 
                         theta_known=theta_known, 
                         theta=initheta, alpha=alpha, tol=tol, 
                         verb=verb, admix=admix, 
                         prior_omega = prior_omega,
                         rp=grp, tmax=tmax, wtol=wtol, qn=qn)
    fit$D <- class.tpxResids(X=X, theta=fit$theta, omega=fit$omega, 
                             grp=grp, nonzero=nonzero)$D
    return(fit)
  }

  if(is.matrix(alpha)){ stop("Matrix alpha only works for fixed K") }
  
  if(verb){ cat(paste("Fit and Bayes Factor Estimation for K =",K[1]))
            if(length(K)>1){ cat(paste(" ...", max(K))) }
            cat("\n") }

  ## dimensions
  n <- nrow(X)
  p <- ncol(X)
  nK <- length(K)
    
  BF <- D <- NULL
  iter <- 0
  
  ## Null model log probability
  sx <- sum(X)
  qnull <- col_sums(X)/sx
  null <- sum( X$v*log(qnull[X$j]) ) - 0.5*(n+p)*(log(sx) - log(2*pi))
  
  ## allocate and initialize
  best <- -Inf
  bestfit <- NULL  
  
  ## loop over topic numbers
  for(i in 1:nK){
    
    ## Solve for map omega in NEF space
    fit <- class.tpxfit(X=X, known_samples=known_samples, omega_known=omega_known, 
                        theta_known=theta_known,theta=initheta, alpha=alpha, tol=tol, verb=verb,
                        admix=admix, prior_omega = prior_omega, grp=grp, tmax=tmax, wtol=wtol, qn=qn)
    
    BF <- c(BF, class.tpxML(X=X, theta=fit$theta, omega=fit$omega, alpha=fit$alpha, L=fit$L, dcut=dcut, admix=admix, grp=grp) - null)
    R <- class.tpxResids(X=X, theta=fit$theta, omega=fit$omega, grp=grp, nonzero=nonzero)
    D <- cbind(D, unlist(R$D))

    if(verb>0) cat(paste("log BF(", K[i], ") =", round(BF[i],2)))
    if(verb>1) cat(paste(" [ ", fit$iter,"steps, disp =",round(D[1,i],2)," ]\n")) else if(verb >0) cat("\n")
    
    if(is.nan(BF[i])){ 
      cat("NAN for Bayes factor.\n")
      return(bestfit)
      break} 
    
    if(BF[i] > best){ # check for a new "best" topic
      best <- BF[i]
      bestfit <- fit
    } else if(kill>0 && i>kill){ # break after kill consecutive drops
      if(prod(BF[i-0:(kill-1)] < BF[i-1:kill])==1) break }
    
    if(i<nK){
      if(!admix){ initheta <- class.tpxinit(X,2,K[i+1], alpha, 0) }
      else{ initheta <- class.tpxThetaStart(X, fit$theta, fit$omega, K[i+1]) }
    }
  }

  names(BF) <- dimnames(D)[[2]] <- paste(K[1:length(BF)]) 
 
  return(list(theta=bestfit$theta, omega=bestfit$omega, alpha=bestfit$alpha,
              BF=BF, D=D, K=K[which.max(BF)])) }

## theta initialization
class.tpxinit <- function(X, 
                          K1,
                          initheta, 
                          K_classes, 
                          method=c("omega.fix", "theta.fix", "no.fix"),
                          alpha, verb){
## initheta can be matrix, or c(nK, tmax, tol, verb)
  
  prior_omega = rep(1/(K1), (K1))
  
  if(is.matrix(initheta)){
    if(ncol(initheta)!=K1){ stop("mis-match between initheta and K.") }
    if(prod(initheta>0) != 1){ stop("use probs > 0 for initheta.") }
    return(class.normalizetpx(initheta, byrow=FALSE)) }

  if(is.matrix(alpha)){
    if(nrow(alpha)!=ncol(X) || ncol(alpha)!=K1){ stop("bad matrix alpha dimensions; check your K") }
    return(class.normalizetpx(alpha, byrow=FALSE)) }

  if(is.null(initheta)){ ilength <- K1-1 } else{ ilength <- initheta[1] }
  if(ilength < 1){ ilength <- 1 }

  ## set number of initial steps
  if(length(initheta)>1){ tmax <- initheta[2] }else{ tmax <- 3 }
  ## set the tolerance
  if(length(initheta)>2){ tol <- initheta[3] }else{ tol <- 0.5 }
  ## print option
  if(length(initheta)>3){ verb <- initheta[4] }else{ verb <- 0 }
  

  if(verb){ cat("Building initial topics") 
            if(verb > 1){ cat(" for K = ") } else{ cat("... ") } }
            
  nK <- length( Kseq <-  unique(ceiling(seq(2,K1,length=ilength))) )
  
  initheta <- tpxThetaStart(X, matrix(col_sums(X)/sum(X), ncol=1),
                            matrix(rep(1,nrow(X))), K1)
  
 # cat("no error till now \n")
  
  if(verb > 0)
    { cat("\n")
      print(list(Kseq=Kseq, tmax=tmax, tol=tol)) }
    
  ## loop over topic numbers
  for(i in 1:nK){

    ## Solve for map omega in NEF space
    fit <- class.tpxfit(X=X, known_samples=NULL,
                        omega_known=NULL, theta=initheta, K_classes=K_classes,
                        alpha=alpha, method=method, tol=tol, verb=verb, 
                        admix=TRUE, prior_omega = prior_omega,
                        grp=NULL, tmax=tmax, wtol=-1, qn=-1)
    if(verb>1){ cat(paste(Kseq[i],",", sep="")) }

  #    cat("kushal \n")
      if(i<nK){ initheta <- tpxThetaStart(X, fit$theta, fit$omega, Kseq[i+1]) }
      else{ initheta <- fit$theta }
    }
  if(verb){ cat("done.\n") }
  return(initheta)
}
               
## ** main workhorse function.  Only Called by the above wrappers.
## topic estimation for a given number of topics (taken as ncol(theta))
class.tpxfit <- function(X, known_samples, omega_known, theta, K_classes,
                         alpha, method=c("omega.fix", "theta.fix", "no.fix"),
                         tol, verb, admix, prior_omega, grp, tmax, wtol, qn)
{
  ## inputs and dimensions
  if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix") }
  K <- ncol(theta)
  n <- nrow(X)
  p <- ncol(X)
  m <- row_sums(X)
  if(is.null(alpha)){ alpha <- 1/(K*p) }
  if(is.matrix(alpha)){ if(nrow(alpha)!=p || ncol(alpha)!=K){ stop("bad matrix alpha dimensions") }}

  ## recycle these in tpcweights to save time
  xvo <- X$v[order(X$i)]
  wrd <- X$j[order(X$i)]-1
  doc <- c(0,cumsum(as.double(table(factor(X$i, levels=c(1:nrow(X)))))))
  
  omega <- matrix(0, n, K);
  
  ## Initialize
  if(method=="omega.fix") {
    X_unknown <- X[-(known_samples),];
    n_unknown <- nrow(X_unknown);
    xvo_unknown <- X_unknown$v[order(X_unknown$i)]
    wrd_unknown <- X_unknown$j[order(X_unknown$i)]-1
    doc_unknown <- c(0,cumsum(as.double(table(factor(X_unknown$i, levels=c(1:nrow(X_unknown)))))))
    
    omega_unknown <- class.tpxweights(n=n_unknown, p=p, xvo=xvo_unknown, wrd=wrd_unknown, doc=doc_unknown, start=class.tpxOmegaStart(X_unknown,theta), theta=theta)
   if(!admix){ omega_unknown <- matrix(apply(omega_unknown,2, function(w) tapply(w,grp,mean)), ncol=K) }
    omega[known_samples,] <- omega_known;
    omega[-(known_samples),] <- omega_unknown;
  }else{
    omega <- class.tpxweights(n=n, p=p, xvo=xvo, wrd=wrd, doc=doc, start=class.tpxOmegaStart(X,theta), theta=theta)
    if(!admix){ omega <- matrix(apply(omega,2, function(w) tapply(w,grp,mean)), ncol=K) }
  }

  ## tracking
  iter <- 0
  dif <- tol+1+qn
  update <- TRUE
  if(verb>0){
    cat("log posterior increase: " )
    digits <- max(1, -floor(log(tol, base=10))) }
  
  Y <- NULL # only used for qn > 0 
  Q0 <- col_sums(X)/sum(X)
  L <- class.tpxlpost(X=X, theta=theta, omega=omega, alpha=alpha, 
                      admix=admix, prior_omega=prior_omega, grp=grp) 
 # if(is.infinite(L)){ L <- sum( (log(Q0)*col_sums(X))[Q0>0] ) }
  
  ## Iterate towards MAP
  while( update  && iter < tmax ){ 

    ## sequential quadratic programming for conditional Y solution
    if(admix && wtol > 0){ 
      Wfit <- matrix(0, n, K);
      if(method=="omega.fix") {
        X_unknown <- X[-(known_samples),];
        Wfit_unknown <- class.tpxweights(n=nrow(X_unknown), p=ncol(X_unknown), xvo=xvo_unknown, wrd=wrd_unknown, doc=doc_unknown,
                                start=omega_unknown, theta=theta,  verb=0, nef=TRUE, wtol=wtol, tmax=20)
        Wfit[known_samples,] <- omega_known
        Wfit[-(known_samples),] <- Wfit_unknown}else{
          Wfit <- class.tpxweights(n=nrow(X), p=ncol(X), xvo=xvo, wrd=wrd, doc=doc,
                             start=omega, theta=theta,  verb=0, nef=TRUE, 
                             wtol=wtol, tmax=20);
        }}else{ Wfit <- omega }

    
    
    ## joint parameter EM update
    if(method=="theta.fix"){
      move1 <- class.tpxEM(X=X, m=m, theta=theta, omega=Wfit, 
                           alpha=alpha, admix=admix, 
                           prior_omega=prior_omega,
                           grp=grp)
      if(K_classes < K){
        move <- list(theta=cbind(theta[,1:K_classes],move1$theta[,((K_classes+1):K)]), omega=move1$omega)
      }else{
        move <- list(theta=theta, omega=move1$omega)
      }
    }else{
        move <- class.tpxEM(X=X, m=m, theta=theta, omega=Wfit, alpha=alpha, 
                           admix=admix, prior_omega=prior_omega, grp=grp)
    }
    
    
    ## quasinewton-newton acceleration
    if(method=="omega.fix"){
        move_unknown <- list(omega=move$omega[-(known_samples),], 
                             theta=move$theta)
        QNup <- class.tpxQN(move=move_unknown, Y=Y, X=X_unknown, 
                            alpha=alpha, verb=verb, admix=admix, 
                            prior_omega=prior_omega,
                            grp=grp, doqn=qn-dif)
        move_unknown <- QNup$move;
        omega_unknown <- move_unknown$omega;
        move$omega[-(known_samples),] <- omega_unknown;
        move$omega[known_samples,] <- omega_known;
        move$theta <- move_unknown$theta;
        QNup$L <-  class.tpxlpost(X=X, theta=move$theta, 
                                  omega=move$omega, alpha=alpha, 
                                  admix=admix, 
                                  prior_omega=prior_omega,
                                  grp=grp)
        } else if (method=="theta.fix"){
          QNup <- class.tpxQN(move=move, Y=Y, X=X, alpha=alpha, 
                              verb=verb, admix=admix, 
                              prior_omega=prior_omega,
                              grp=grp, 
                              doqn=qn-dif)
          move$omega <- QNup$move$omega;
        } else{
          QNup <- class.tpxQN(move=move, Y=Y, X=X, alpha=alpha, 
                              verb=verb, admix=admix, 
                              prior_omega=prior_omega,
                              grp=grp, 
                              doqn=qn-dif)
          move <- QNup$move;
      }
    Y <- QNup$Y
    
    if(QNup$L < L){  # happens on bad Wfit, so fully reverse
      if(verb > 10){ cat("_reversing a step_") }
      if(method=="theta.fix"){
        if(K_classes < K){
          move1 <- class.tpxEM(X=X, m=m, theta=theta, omega=omega, 
                               alpha=alpha, admix=admix, 
                               prior_omega = prior_omega, grp=grp)
          move <- list(theta=cbind(theta[,1:K_classes],move1$theta[,((K_classes+1):K)]), omega=move1$omega)
        }else{
          move <- list(theta=theta, omega=omega)
        }
      }else{
        move <- class.tpxEM(X=X, m=m, theta=theta, omega=omega, 
                            alpha=alpha, 
                            admix=admix, prior_omega=prior_omega,
                            grp=grp)
      }
      
      QNup$L <-  class.tpxlpost(X=X, theta=move$theta, 
                                omega=move$omega, alpha=alpha, 
                                admix=admix, 
                                prior_omega=prior_omega, 
                                grp=grp) }
      #L <-  class.tpxlpost(X=X, theta=theta, omega=omega, alpha=alpha, admix=admix, grp=grp) 
  
    ## calculate dif
    dif <- (QNup$L-L)
    #print(dif)
    L <- QNup$L
    
        
    ## check convergence
    if(abs(dif) < tol){
      if(sum(abs(theta-move$theta)) < tol){ update = FALSE } }

    ## print
    if(verb>0 && (iter-1)%%ceiling(10/verb)==0 && iter>0){
      cat( paste( round(abs(dif),digits), #" (", sum(abs(theta-move$theta)),")",
                 ", ", sep="") ) }
    
    ## heartbeat for long jobs
    if(((iter+1)%%1000)==0){ 
          cat(sprintf("p %d iter %d diff %g\n",
                nrow(theta), iter+1,round(dif))) }

    ## iterate
    iter <- iter+1
    theta <- move$theta
    omega <- move$omega
    
  }

  ## final log posterior
  L <- class.tpxlpost(X=X, theta=theta, omega=omega, 
                      alpha=alpha, admix=admix,
                      prior_omega=prior_omega,
                      grp=grp) 

  ## summary print
  if(verb>0){
    cat("done.")
    if(verb>1) { cat(paste(" (L = ", round(L,digits), ")", sep="")) }
    cat("\n")
  }
  
  out <- list(theta=theta, omega=omega, K=K, alpha=alpha, L=L, iter=iter)
  invisible(out) }

 
## ** called from topics.R (predict) and class.tpx.R
## Conditional solution for topic weights given theta
class.tpxweights <- function(n, p, xvo, wrd, doc, start, theta, verb=FALSE, nef=TRUE, wtol=10^{-5}, tmax=1000)
{
  K <- ncol(theta)
  start[start == 0] <- 0.1/K
  start <- start/rowSums(start) 
  omega <- .C("Romega",
              n = as.integer(n),
              p = as.integer(p),
              K = as.integer(K),
              doc = as.integer(doc),
              wrd = as.integer(wrd),
              X = as.double(xvo),
              theta = as.double(theta),
              W = as.double(t(start)),
              nef = as.integer(nef),
              tol = as.double(wtol),
              tmax = as.integer(tmax),
              verb = as.integer(verb),
              PACKAGE="classtpx")
  return(t(matrix(omega$W, nrow=ncol(theta), ncol=n))) }

## ** Called only in class.tpx.R

## single EM update. two versions: admix and mix
class.tpxEM <- function(X, m, theta, omega, alpha, admix, prior_omega, grp)
{
  n <- nrow(X)
  p <- ncol(X)
  K <- ncol(theta)

  if(admix){ Xhat <- (X$v/class.tpxQ(theta=theta, omega=omega, doc=X$i, wrd=X$j))*(omega[X$i,]*theta[X$j,])
             Zhat <- .C("Rzhat", n=as.integer(n), p=as.integer(p), K=as.integer(K), N=as.integer(nrow(Xhat)),
                         Xhat=as.double(Xhat), doc=as.integer(X$i-1), wrd=as.integer(X$j-1),
                        zj = as.double(rep(0,K*p)), zi = as.double(rep(0,K*n)), PACKAGE="classtpx")
             theta <- class.normalizetpx(matrix(Zhat$zj+alpha, ncol=K), byrow=FALSE)
             omega <- class.normalizetpx(matrix(Zhat$zi+prior_omega, ncol=K)) 
             } else{
    qhat <- class.tpxMixQ(X, omega, theta, grp, qhat=TRUE)$qhat
    ## EM update
    theta <- class.normalizetpx(tcrossprod_simple_triplet_matrix( t(X), t(qhat) ) + alpha, byrow=FALSE)
    omega <- class.normalizetpx(matrix(apply(qhat*m,2, function(x) tapply(x,grp,sum)), ncol=K)+prior_omega )  
    }
    
  return(list(theta=theta, omega=omega)) }

## Quasi Newton update for q>0 
class.tpxQN <- function(move, Y, X, alpha, verb, admix, 
                        prior_omega, grp, doqn)
{
  ## always check likelihood
  L <- class.tpxlpost(X=X, theta=move$theta, omega=move$omega,
                alpha=alpha, admix=admix, prior_omega=prior_omega,
                grp=grp) 
  move$omega[move$omega<=0] <- 1e-20;
  move$theta[move$theta<=0] <- 1e-20;
  move$omega <- class.normalizetpx(move$omega);
  move$theta <- class.normalizetpx(move$theta, byrow=FALSE);
  
  if(doqn < 0){ return(list(move=move, L=L, Y=Y)) }

  ## update Y accounting
  Y <- cbind(Y, class.tpxToNEF(theta=move$theta, omega=move$omega))
  if(ncol(Y) < 3){ return(list(Y=Y, move=move, L=L)) }
  if(ncol(Y) > 3){ warning("mis-specification in quasi-newton update; please report this bug.") }
  
  ## Check quasinewton secant conditions and solve F(x) - x = 0.
  U <- as.matrix(Y[,2]-Y[,1])
  V <- as.matrix(Y[,3]-Y[,2])
  sUU <- sum(U^2)
  sVU <- sum(V*U)
  Ynew <- Y[,3] + V*(sVU/(sUU-sVU)) 
  qnup <- class.tpxFromNEF(Ynew, n=nrow(move$omega),
                     p=nrow(move$theta), K=ncol(move$theta))

  ## check for a likelihood improvement
  Lqnup <- try(class.tpxlpost(X=X, theta=qnup$theta, omega=qnup$omega,
                        alpha=alpha, admix=admix, prior_omega=prior_omega, 
                        grp=grp), silent=TRUE)
  
  if(inherits(Lqnup, "try-error")){
    if(verb>10){ cat("(QN: try error) ") }
    return(list(Y=Y[,-1], move=move, L=L)) }
  
  if(verb>10){ cat(paste("(QN diff ", round(Lqnup-L,3), ")\n", sep="")) }
  
  if(Lqnup < L){
    return(list(Y=Y[,-1], move=move, L=L)) }
  else{
    L <- Lqnup
    Y <- cbind(Y[,2],Ynew)
    return( list(Y=Y, move=qnup, L=L) )
  }
}

  
## unnormalized log posterior (objective function)
class.tpxlpost <- function(X, theta, omega, alpha, admix=TRUE, 
                           prior_omega, grp=NULL)
{
  K <- ncol(theta);
  omega[omega<=0] <- 1e-20;
  theta[theta<=0] <- 1e-20;
  if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix.") }

  if(admix){ L <- sum( X$v*log(class.tpxQ(theta=theta, omega=omega, doc=X$i, wrd=X$j)) ) 
  }else{ L <- sum(class.tpxMixQ(X, omega, theta, grp)$lqlhd) }
  if(is.null(nrow(alpha))){ if(alpha != 0){ L <- L + sum(alpha*log(theta))  } } # unnormalized prior
  if(length(prior_omega)>1){
    L <- L + sum(log(omega)*(rep.row(prior_omega, dim(omega)[1])));
  }else{
    L <- L + sum(log(omega));
  }
  
  return(L) }

## log marginal likelihood
class.tpxML <- function(X, theta, omega, alpha, L, dcut, admix=TRUE, grp=NULL){
  ## get the indices
  K <- ncol(theta)
  p <- nrow(theta)
  n <- nrow(omega)

  ## return BIC for simple finite mixture model
  if(!admix){
    qhat <- class.tpxMixQ(X, omega, theta, grp, qhat=TRUE)$qhat
    ML <- sum(X$v*log(row_sums(qhat[X$i,]*theta[X$j,])))
    return( ML - 0.5*( K*p + (K-1)*n )*log(sum(X)) ) } 

  ML <- L  + lfactorial(K) # lhd multiplied by label switching modes

  ## block-diagonal approx to determinant of the negative log hessian matrix
  q <- class.tpxQ(theta=theta, omega=omega, doc=X$i, wrd=X$j)
  D <- class.tpxHnegDet(X=X, q=q, theta=theta, omega=omega, alpha=alpha)

  D[D < dcut] <- dcut
  ML <- ML - 0.5*sum( D  )   # -1/2 |-H|

  ML <- ML + (K*p + sum(omega>0.01))*log(2*pi)/2  # (d/2)log(2pi)
  if(is.null(nrow(alpha))){ # theta prior normalizing constant
    ML <- ML + K*( lgamma(p*(alpha+1)) - p*lgamma(alpha+1) )  }
  else{ ML <- ML + sum(lgamma(col_sums(alpha+1)) - col_sums(lgamma(alpha+1))) } # matrix version
  ## omega(phi) prior normalizing constant number of parameters
  ML <- ML +  sum(D[-(1:p)]>dcut)*( lfactorial(K) - K*lgamma( 1+1/K ) ) #
  
  return(ML) }

## find residuals for X$v
class.tpxResids <- function(X, theta, omega, grp=NULL, nonzero=TRUE)
{
  if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix.") }

  m <- row_sums(X)
  K <- ncol(theta)
  n <- nrow(X)
  phat <- sum(col_sums(X)>0)
  d <- n*(K-1) + K*( phat-1 )

  if(nrow(omega) == nrow(X)){
    qhat <- class.tpxQ(theta=theta, omega=omega, doc=X$i, wrd=X$j)
    xhat <- qhat*m[X$i]
  } else{
    q <- class.tpxMixQ(X=X, omega=omega, theta=theta, grp=grp, qhat=TRUE)$qhat
    qhat <- row_sums(q[X$i,]*theta[X$j,])
    xhat <- qhat*m[X$i] }

  if(nonzero || nrow(omega) < nrow(X)){
    ## Calculations based on nonzero counts
    ## conditional adjusted residuals
    e <- X$v^2 - 2*(X$v*xhat - xhat^2)
    s <- qhat*m[X$i]*(1-qhat)^{1-m[X$i]}
    r <- sqrt(e/s)
    df <- length(r)*(1-d/(n*phat))
    R <- sum(r^2) 
  }
  else{
    ## full table calculations
    e <- (X$v^2 - 2*X$v*m[X$i]*qhat)
    s <- m[X$i]*qhat*(1-qhat)
    fulltable <- .C("RcalcTau",
                 n = as.integer(nrow(omega)),
                 p = as.integer(nrow(theta)),
                 K = as.integer(ncol(theta)),
                 m = as.double(m),
                 omega = as.double(omega),
                 theta = as.double(theta),
                 tau = double(1), size=double(1),
                 PACKAGE="classtpx" )
    tau <- fulltable$tau 
    R <- sum(e/s) + tau
    df <-  fulltable$size - phat  - d
    r <- suppressWarnings(sqrt(e/s + tau))
    r[is.nan(r)] <- 0 ## should not happen, but can theoretically
  }
  
  ## collect and output
  sig2 <- R/df
  rho <- suppressWarnings(pchisq(R, df=df, lower.tail=FALSE))
  D <- list(dispersion=sig2, pvalue=rho, df=df)
  return( list(s=s, e=e, r=r, D=D) ) }

  
## fast initialization functions for theta (after increasing K) and omega (given theta)
class.tpxThetaStart <- function(X, theta, omega, K)
  {
    R <- class.tpxResids(X, theta=theta, omega=omega, nonzero=TRUE) 
    X$v <- R$e*(R$r>3) + 1/ncol(X)
    Kpast <- ncol(theta)
    Kdiff <- K-Kpast
    if(Kpast != ncol(omega) || Kpast >= K){ stop("bad K in class.tpxThetaStart") }
    initheta <- class.normalizetpx(Kpast*theta+rowMeans(theta), byrow=FALSE)
    n <- nrow(X)
    ki <- matrix(1:(n-n%%Kdiff), ncol=Kdiff)
    for(i in 1:Kdiff){ initheta <- cbind(initheta, (col_sums(X[ki[,i],])+1/ncol(X))/(sum(X[ki[,i],])+1)) }
    return( initheta )
  }

class.tpxOmegaStart <- function(X, theta)
  {
    if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix.") }
    omega <- try(tcrossprod_simple_triplet_matrix(X, solve(t(theta)%*%theta)%*%t(theta)), silent=TRUE )
    if(inherits(omega,"try-error")){ return( matrix( 1/ncol(theta), nrow=nrow(X), ncol=ncol(theta) ) ) }
    omega[omega <= 0] <- .5
    return(class.normalizetpx(omega, byrow=TRUE) )
  }


## fast computation of sparse P(X) for X>0
class.tpxQ <- function(theta, omega, doc, wrd){
  theta[theta==0] <- 1e-14;
  theta[theta==1] <- 1 - 1e-14
  omega[omega==0] <-  1e-14
  omega[omega==1] <- 1 - 1e-14
  
  omega <- class.normalizetpx(omega, byrow = TRUE)
  theta <- class.normalizetpx(theta, byrow = FALSE)
  
  if(length(wrd)!=length(doc)){stop("index mis-match in class.tpxQ") }
  if(ncol(omega)!=ncol(theta)){stop("theta/omega mis-match in class.tpxQ") }
  
  out <- .C("RcalcQ",
            n = as.integer(nrow(omega)),
            p = as.integer(nrow(theta)),
            K = as.integer(ncol(theta)),
            doc = as.integer(doc-1),
            wrd = as.integer(wrd-1),
            N = as.integer(length(wrd)),
            omega = as.double(omega),
            theta = as.double(theta),
            q = double(length(wrd)),
            PACKAGE="classtpx" )

  return( out$q ) }

## model and component likelihoods for mixture model
class.tpxMixQ <- function(X, omega, theta, grp=NULL, qhat=FALSE){
  if(is.null(grp)){ grp <- rep(1, nrow(X)) }
  K <- ncol(omega)
  n <- nrow(X)
  mixhat <- .C("RmixQ",
               n = as.integer(nrow(X)),
               p = as.integer(ncol(X)),
               K = as.integer(K),
               N = as.integer(length(X$v)),
               B = as.integer(nrow(omega)),
               cnt = as.double(X$v),
               doc = as.integer(X$i-1),
               wrd = as.integer(X$j-1),
               grp = as.integer(as.numeric(grp)-1),
               omega = as.double(omega),
               theta = as.double(theta),
               Q = double(K*n),
               PACKAGE="classtpx")
  ## model and component likelihoods
  lQ <- matrix(mixhat$Q, ncol=K)
  lqlhd <- log(row_sums(exp(lQ)))
  lqlhd[is.infinite(lqlhd)] <- -600 # remove infs
  if(qhat){
    qhat <- exp(lQ-lqlhd)
    ## deal with numerical overload
    infq <- row_sums(qhat) < .999
    if(sum(infq)>0){
      qhat[infq,] <- 0
      qhat[n*(apply(matrix(lQ[infq,],ncol=K),1,which.max)-1) + (1:n)[infq]] <- 1 }
  }
  return(list(lQ=lQ, lqlhd=lqlhd, qhat=qhat)) }

## negative log hessian block diagonal matrix for theta & omega
class.tpxHnegDet <- function(X, q, theta, omega, alpha){
  K <- ncol(theta)
  n <- nrow(omega)

  ## sparse Xij/Qij^2 
  Xq <- X
  Xq$v <- Xq$v/q^2
  
  ## negative 2nd derivitive matrices for theta
  HT <- tcrossprod_simple_triplet_matrix(t(Xq), apply(omega, 1, function(v) v%o%v ) )
  HT[,K*(0:(K-1))+1:K] <- HT[,K*(0:(K-1))+1:K] + alpha/theta^2 # will break for alpha<=1
  DT <- apply(HT, 1, class.tpxlogdet)

  ## ditto for omega
  HW <- matrix(.C("RnegHW",
                  n = as.integer(nrow(omega)),
                  p = as.integer(nrow(theta)),
                  K = as.integer(K-1),
                  omeg = as.double(omega[,-1]),
                  thet = as.double(theta[,-1]),
                  doc = as.integer(X$i-1),
                  wrd = as.integer(X$j-1),
                  cnt = as.double(X$v),
                  q = as.double(q),
                  N = as.integer(length(q)),
                  H = double(n*(K-1)^2),
                  PACKAGE="classtpx")$H,
               nrow=(K-1)^2, ncol=n)
  DW <- apply(HW, 2, class.tpxlogdet) 
  return( c(DT,DW) )  }

## functions to move theta/omega to and from NEF.  
class.tpxToNEF <- function(theta, omega){
  n <- nrow(omega)
  p <- nrow(theta)
  K <- ncol(omega)
  return(.C("RtoNEF",
            n=as.integer(n), p=as.integer(p), K=as.integer(K),
            Y=double((p-1)*K + n*(K-1)),
            theta=as.double(theta), tomega=as.double(t(omega)),
            PACKAGE="classtpx")$Y)
}

## 'From' NEF representation back to probabilities
class.tpxFromNEF <- function(Y, n, p, K){
  bck <- .C("RfromNEF",
            n=as.integer(n), p=as.integer(p), K=as.integer(K),
            Y=as.double(Y), theta=double(K*p), tomega=double(K*n),
            PACKAGE="classtpx")
  return(list(omega=t( matrix(bck$tomega, nrow=K) ), theta=matrix(bck$theta, ncol=K)))
}

## utility log determinant function for speed/stabilty
class.tpxlogdet <- function(v){
    v <- matrix(v, ncol=sqrt(length(v)))
    
    if( sum(zeros <- colSums(v)==0)!=0 ){
      cat("warning: boundary values in laplace approx\n")
      v <- v[-zeros,-zeros] }
   
    return(determinant(v, logarithm=TRUE)$modulus) }
