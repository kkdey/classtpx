library(ashr)
library(truncdist)

#' @title Main F-statistics Adaptive Shrinkage function
#' @description Takes observed F-statistics (fhat) and their degrees of freedom (df1, df2). Suppose fhat=alpha*F, where F is an F-distributed random variable.
#' We applies shrinkage to the F-stats, using Empirical Bayes methods, to compute shrunk estimates for alpha.
#' 
#' @param fhat  a p vector of F-statistics 
#' @param df1 the first degree of freedom for (F) distribution of fhat
#' @param df2 the second degree of freedom for (F) distribution of fhat
#' @param method specifies how ash is to be run. Can be "shrinkage" (if main aim is shrinkage) or "fdr" (if main aim is to assess fdr or fsr)
#' This is simply a convenient way to specify certain combinations of parameters: "shrinkage" sets pointmass=FALSE and prior="uniform";
#' "fdr" sets pointmass=TRUE and prior="nullbiased".
#' @param optmethod specifies optimization method used. Default is "mixIP", an interior point method, if REBayes is installed; otherwise an EM algorithm is used. The interior point method is faster for large problems (n>2000).
#' @param nullweight scalar, the weight put on the prior under "nullbiased" specification, see \code{prior}
#' @param randomstart logical, indicating whether to initialize EM randomly. If FALSE, then initializes to prior mean (for EM algorithm)
#' @param pointmass logical, indicating whether to use a point mass at zero as one of components for the mixture prior g
#' @param prior string, or numeric vector indicating Dirichlet prior on mixture proportions (defaults to "uniform", or (1,1...,1); also can be "nullbiased" (nullweight,1,...,1) to put more weight on first component), or "unit" (1/K,...,1/K) [for optmethod=mixVBEM version only]
#' @param mixsd vector of sds for underlying mixture components of prior g 
#' @param g the prior distribution for log(alpha) (usually estimated from the data; this is used primarily in simulated data to do computations with the "true" g)
#' @param gridmult the multiplier by which the default grid values for mixsd differ by one another. (Smaller values produce finer grids)
#' @param control A list of control parameters for the optmization algorithm. Default value is set to be   control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE). User may supply changes to this list of parameter, say, control=list(maxiter=10000,trace=TRUE)

#' @return fash returns a list with some or all of the following elements \cr
#' \item{fitted.g}{fitted prior g, an uniform mixture}
#' \item{fit}{a list, including the fitted prior g, the log-likelihood of model and mixture fitting convergence information}
#' \item{PosteriorMean.f}{A vector consisting the posterior mean of alpha}
#' \item{PosteriorMean.logf}{A vector consisting the posterior mean of log(alpha)}
#' \item{PositiveProb}{A vector of posterior probability that log(alpha) is positive}
#' \item{NegativeProb}{A vector of posterior probability that log(alpha) is negative}
#' \item{ZeroProb}{A vector of posterior probability that log(alpha) is zero}
#' \item{lfsr}{The local false sign rate}
#' \item{lfdr}{A vector of estimated local false discovery rate}
#' \item{qvalue}{A vector of q values}
#'
#' @export
#' @examples 
#' alpha = exp(rnorm(200,0,0.1))
#' fhat = alpha*rf(200,df1=10,df2=12)
#' fhat.fash = fash(fhat,df1=10,df2=12)
#' plot(fhat,beta.ash$PosteriorMean.f)
#' 
fash = function(fhat, df1, df2,
                method = c("fdr","shrink"),              
                oneside = FALSE,
                optmethod = c("mixIP","mixEM"),
                nullweight = 10,
                randomstart = FALSE,
                pointmass = TRUE,
                prior = c("nullbiased","uniform"),
                mixsd = NULL,
                g = NULL,
                gridmult = 2,
                control = list()){
  
  
  if(!missing(method)){
    method = match.arg(method) 
    if(method=="shrink"){
      if(missing(prior)){
        prior = "uniform"
      } else {
        warning("Specification of prior overrides default for method shrink")
      }
      if(missing(pointmass)){
        pointmass=FALSE
      } else {
        warning("Specification of pointmass overrides default for method shrink")
      }    
    }
    
    if(method=="fdr"){
      if(missing(prior)){
        prior = "nullbiased"
      } else {
        warning("Specification of prior overrides default for method fdr")
      }
      if(missing(pointmass)){
        pointmass=TRUE
      } else {
        warning("Specification of pointmass overrides default for method fdr")
      }
    }  
  }
  
  if(missing(optmethod)){
    if(require(REBayes,quietly=TRUE)){ #check whether REBayes package is present
      optmethod = "mixIP"
    } else{  #If REBayes package missing
      message("Due to absence of package REBayes, switching to EM algorithm")
      optmethod = "mixEM" #fallback if neither Rcpp or REBayes are installed
      message("Using vanilla EM; for faster performance install REBayes (preferred) or Rcpp")  
    }
  } else { #if optmethod specified
    optmethod = match.arg(optmethod)
  }
  
  if(!is.numeric(prior)){  prior = match.arg(prior)  } 
  if(gridmult<=1){  stop("gridmult must be > 1")  }
  
  logfhat = log(fhat)
  logfhat[df1==Inf | df2==Inf]=NA
  completeobs = (!is.na(fhat))
  n=sum(completeobs)
  
  #Handling control variables
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  if(n>50000){control.default$trace=TRUE}
  namc=names(control)
  if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput=modifyList(control.default, control)
  
  if(!is.null(g)){
    controlinput$maxiter = 0 # if g is specified, don't iterate the EM
    k = ncomp(g)
    null.comp=1 #null.comp not actually used unless randomstart true 
    prior = setprior(prior,k,nullweight,null.comp)
    if(randomstart){pi = initpi(k,n,null.comp,randomstart)
                    g$pi=pi} #if g specified, only initialize pi if randomstart is TRUE 
  }else {
    if(is.null(mixsd)){
      mixsd = autoselect.mixsd(logfhat[completeobs],df1,df2,gridmult)
    }
    if(pointmass){
      mixsd = c(0,mixsd)
    }    
    
    null.comp = which.min(mixsd) #which component is the "null"
    
    k = length(mixsd)
    prior = setprior(prior,k,nullweight,null.comp)
    pi = initpi(k,n,null.comp,randomstart)
    if(oneside==FALSE){
      g = unimix(pi,-mixsd,mixsd)
    }else{
      g = unimix(pi,rep(0,length(mixsd)),mixsd)
    }  
  }
  
  
  
  pi.fit = est_mixprop(logfhat,df1,df2,g,prior,optmethod, null.comp=null.comp,
                      control=controlinput)
  
  # ZeroProb/NegProb of log(f)!
  ####### TBD: lfsr when oneside==TRUE?
  ZeroProb = rep(0,length(logfhat))
  NegativeProb = rep(0,length(logfhat))
  ZeroProb[completeobs] = colSums(comppostprob_logf(pi.fit$g,logfhat[completeobs],df1,df2)[comp_sd(pi.fit$g)==0,,drop=FALSE])
  NegativeProb[completeobs] = cdf_post_logf(pi.fit$g, 0, logfhat[completeobs],df1,df2) - ZeroProb[completeobs]
  
  PosteriorMean.f = rep(0,length(logfhat))
  PosteriorMean.f[completeobs] = postmean_f(pi.fit$g,logfhat[completeobs],df1,df2)
  
  PosteriorMean.logf = rep(0,length(logfhat))
  PosteriorMean.logf[completeobs] = postmean_logf(pi.fit$g,logfhat[completeobs],df1,df2)
  #PosteriorSD[completeobs] = postsd(pi.fit$g,betahat[completeobs],sebetahat[completeobs],df)
  
  PositiveProb = 1- NegativeProb-ZeroProb
  lfsr = compute_lfsr(NegativeProb,ZeroProb)
  lfsra = compute_lfsra(PositiveProb,NegativeProb,ZeroProb) 
  lfdr = ZeroProb
  qvalue = qval.from.lfdr(lfdr)
  
  return(list(fitted.g=pi.fit$g, fit=pi.fit, 
              PosteriorMean.f=PosteriorMean.f, PosteriorMean.logf=PosteriorMean.logf,
              ZeroProb=ZeroProb, NegativeProb=NegativeProb, PositiveProb=PositiveProb,
              qvalue=qvalue, lfsr=lfsr, lfsra=lfsra, lfdr=lfdr))
  
}

# estimate g's mixture proportion
est_mixprop = function(logfhat,df1,df2,g,prior,optmethod,null.comp=1,control=list()){ 
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  namc=names(control)
  if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput=modifyList(control.default, control)
  
  pi.init = g$pi
  k = length(g$pi)
  n = length(logfhat)
  controlinput$tol = min(0.1/n,1.e-7) # set convergence criteria to be more stringent for larger samples
  
  if(controlinput$trace==TRUE){tic()}
  
  matrix_lik = t(compdens_conv_logf(g,logfhat,df1,df2))
  
  if (optmethod=="mixIP"){
    EMfit = mixIP(matrix_lik,prior,pi.init,control=controlinput)
  }else if (optmethod=="mixEM"){
    EMfit = mixEM(matrix_lik,prior,pi.init,control=controlinput)
  }
  
  if(!EMfit$converged & controlinput$maxiter>0){
    warning("EM algorithm in function mixEM failed to converge. Results may be unreliable. Try increasing maxiter and rerunning.")
  }
  
  pi = EMfit$pihat     
  penloglik = EMfit$B 
  converged = EMfit$converged
  niter = EMfit$niter
  
  loglik.final =  penloglik(pi,matrix_lik,1) #compute penloglik without penalty
  null.loglik = sum(log(matrix_lik[,null.comp]))  
  
  g$pi = pi
  if(controlinput$trace==TRUE){toc()}
  
  return(list(loglik=loglik.final,null.loglik=null.loglik,
              matrix_lik=matrix_lik,converged=converged,g=g))
}

normalize = function(x){return(x/sum(x))}

# penalized log-likelihood
penloglik = function(pi, matrix_lik, prior){
  pi = normalize(pmax(0,pi))
  m  = t(pi * t(matrix_lik)) # matrix_lik is n by k; so this is also n by k
  m.rowsum = rowSums(m)
  loglik = sum(log(m.rowsum))
  subset = (prior != 1.0)
  priordens = sum((prior-1)[subset]*log(pi[subset]))
  return(loglik+priordens)
}


compdens_conv_logf= function(m,x,v1,v2,FUN="+"){
  if(FUN!="+") stop("Error; compdens_conv not implemented for uniform with FUN!=+")
  compdens = t(pf(exp(outer(x,m$a,FUN="-")),df1=v1,df2=v2)-pf(exp(outer(x,m$b,FUN="-")),df1=v1,df2=v2))/(m$b-m$a)
  compdens[m$a==m$b,] = t(df(exp(outer(x,m$a,FUN="-")),df1=v1,df2=v2)*exp(outer(x,m$a,FUN="-")))[m$a==m$b,]
  return(compdens)
}

comppostprob_logf = function(m,x,v1,v2){
  tmp= (t(m$pi * compdens_conv_logf(m,x,v1,v2))/dens_conv_logf(m,x,v1,v2))
  ismissing = (is.na(x))
  tmp[ismissing,]=m$pi
  t(tmp)
}

dens_conv_logf = function(m,x,v1,v2,FUN="+"){
  colSums(m$pi * compdens_conv_logf(m,x,v1,v2,FUN))
}

cdf_post_logf=function(m,c,logfhat,v1,v2){
  colSums(comppostprob_logf(m,logfhat,v1,v2)*compcdf_post_logf(m,c,logfhat,v1,v2))
}

compcdf_post_logf=function(m,c,logfhat,v1,v2){
  k = length(m$pi)
  n = length(logfhat)
  tmp = matrix(1,nrow=k,ncol=n)
  tmp[m$a > c,] = 0
  subset = m$a<=c & m$b>c # subset of components (1..k) with nontrivial cdf
  if(sum(subset)>0){
    #pna = pf(exp(outer(logfhat,m$a[subset],FUN="-")), df1=v1, df2=v2)
    #pnc = pf(exp(outer(logfhat,rep(c,sum(subset)),FUN="-")), df1=v1, df2=v2)
    #pnb = pf(exp(outer(logfhat,m$b[subset],FUN="-")), df1=v1, df2=v2)
    pna = pf(exp(outer(m$a[subset],logfhat,FUN="-")), df1=v2, df2=v1)
    pnc = pf(exp(outer(rep(c,sum(subset)),logfhat,FUN="-")), df1=v2, df2=v1)
    pnb = pf(exp(outer(m$b[subset],logfhat,FUN="-")), df1=v2, df2=v1)
    tmp[subset,] = t((pnc-pna)/(pnb-pna))
  }
  subset = (m$a == m$b) #subset of components with trivial cdf
  tmp[subset,]= rep(m$a[subset] <= c,n)
  #Occasionally we would encounter issue such that in some entries pna[i,j]=pnb[i,j]=pnc[i,j]=0 or pna=pnb=pnc=1
  #Those are the observations with significant betahat(small sebetahat), resulting in pnorm() return 1 or 0
  #due to the thin tail property of normal distribution.(or t-distribution, although less likely to occur)
  #Then R would be dividing 0 by 0, resulting  in NA values
  #In practice, those observations would have 0 probability of belonging to those "problematic" components
  #Thus any sensible value in [0,1] would not matter much, as they are highly unlikely to come from those 
  #components in posterior distribution.
  #Here we simply assign the "naive" value as as (c-a)/(b-a)
  #As the component pdf is rather smaller over the region.
  tmpnaive=matrix(rep((c-m$a)/(m$b-m$a),length(logfhat)),nrow=k,ncol=n)
  tmp[is.nan(tmp)]= tmpnaive[is.nan(tmp)]
  tmp
}

postmean_f = function(m,logfhat,v1,v2){
  colSums(comppostprob_logf(m,logfhat,v1,v2) * comp_postmean_f(m,logfhat,v1,v2))
}

# Posterior Mean of alpha, not log(alpha)!
comp_postmean_f = function(m,logfhat,v1,v2){
  alpha = exp(outer(-logfhat, m$a,FUN="+"))
  beta = exp(outer(-logfhat, m$b, FUN="+"))   
  tmp = matrix(my_etruncf_vec(cbind(c(alpha),c(beta)),v2,v1), # here reverse v1 and v2!
               nrow=length(logfhat))
  tmp = exp(logfhat)*tmp
  
  ismissing = is.na(logfhat)
  tmp[ismissing,]= exp((m$a+m$b)/2)
  t(tmp)
}

# vectorized funtion to compute mean of truncated F-distribution on [a,b]
# a_b: first column is vector a, second column is vector b
my_etruncf_vec = function(a_b,v1,v2,log=FALSE){
  tmp = apply(a_b,1,my_etruncf,v1=v1,v2=v2,
              lowthres=qf(1e-10,v1,v2),
              highthres=qf(1-1e-10,v1,v2))
}

# compute mean of truncated F-distribution on [a,b]
my_etruncf= function(a_b,v1,v2,lowthres,highthres){
  a = a_b[1]
  b = a_b[2]
  if (a==b){
    tmp = a
  }else{
    if (b>highthres | a<lowthres){
      tmp = etruncf(df1=v1, df2=v2, a=a, b=b)
      if (is.na(tmp)){
        tmp = try(extrunc(spec="f",df1=v1, df2=v2, a=a, b=b), silent=TRUE)
      }
    }else{
      tmp = try(extrunc(spec="f", df1=v1, df2=v2, a=a, b=b), silent=TRUE)
    }
  }
  
  if (class(tmp)=="try-error"){
    tmp=NA
  }
  return(tmp) #deal with extreme case a=b
}

etruncf_num = function(x,df1,df2,a,b){
  c = 10^(-round((log10(df(a,df1,df2))+log10(df(b,df1,df2)))/2))
  c*x*df(x,df1=df1,df2=df2)
}

etruncf_denom = function(x,df1,df2,a,b){
  c = 10^(-round((log10(df(a,df1,df2))+log10(df(b,df1,df2)))/2))
  c*df(x,df1=df1,df2=df2)
}

etruncf = function(df1,df2,a,b){
  thresq = list(low=qf(1e-10,df1,df2), high=qf(1-1e-10,df1,df2))
  n = try(integrate(etruncf_num, lower=a, upper=b, df1=df1,df2=df2,a=a,b=b)$value,
          silent=TRUE)
  d = try(integrate(etruncf_denom, lower=a, upper=b, df1=df1,df2=df2,a=a,b=b)$value,
          silent=TRUE)
  
  if(class(n)=="try-error" | class(d)=="try-error"){
    if (a>thresq$high){
      return(a)
    }else if(b<thresq$low){
      return(b)
    }else if(pf(a,df1,df2)<1e-4 & pf(b,df1,df2)>(1-1e-4)){
      return(df2/(df2-2)) # return the mean of un-truncated F distribution
    }else{
      for (j in seq(10,4,by=-1)){
        if(b > qf(1-10^(-j),df1,df2)){b = qf(1-10^(-j),df1,df2)}
        if(a < qf(10^(-j),df1,df2)){a = qf(10^(-j),df1,df2)}
        n = try(integrate(etruncf_num,lower=a,upper=b,df1=df1,df2=df2,a=a,b=b)$value,silent=TRUE)
        d = try(integrate(etruncf_denom,lower=a,upper=b,df1=df1,df2=df2,a=a,b=b)$value,silent=TRUE)
        
        if(class(n)!="try-error" & class(d)!="try-error"){return(n/d)}
      }  
      ###### NOTE: what if still gets error?
      return(NA)
    }
  }else{
    return(n/d)
  } 
}

# Posterior Mean of log(alpha)
postmean_logf = function(m,logfhat,v1,v2){
  colSums(comppostprob_logf(m,logfhat,v1,v2) * comp_postmean_logf(m,logfhat,v1,v2))
}

comp_postmean_logf = function(m,logfhat,v1,v2){
  alpha = outer(-logfhat, m$a,FUN="+")
  beta = outer(-logfhat, m$b, FUN="+")   
  tmp = matrix(my_etrunclogf_vec(cbind(c(alpha),c(beta)),v2,v1),
               nrow=length(logfhat))
  tmp = logfhat+tmp
  
  ismissing = is.na(logfhat)
  tmp[ismissing,]= (m$a+m$b)/2
  t(tmp)
}

# vectorized funtion to compute mean of truncated log-F distribution on [a,b]
my_etrunclogf_vec = function(a_b,v1,v2,log=FALSE){
  tmp = apply(a_b,1,my_etrunclogf,v1=v1,v2=v2,
              lowthres=log(qf(1e-10,v1,v2)),
              highthres=log(qf(1-1e-10,v1,v2)))
}

my_etrunclogf= function(a_b,v1,v2,lowthres,highthres){
  a = a_b[1]
  b = a_b[2]
  if (a==b){
    tmp = a
  }else{
    if (a<=highthres & b>=lowthres){
      tmp = etrunclogf(df1=v1, df2=v2, a=a, b=b, adj=FALSE)
    }else{
      tmp = etrunclogf(df1=v1, df2=v2, a=a, b=b, adj=TRUE)
    }
    
  }
  return(tmp) #deal with extreme case a=b
}

etrunclogf_num = function(x,df1,df2,a,b){
  c = 10^(-round(min(log10(df(exp(a),df1,df2)*exp(a)),
                     log10(df(exp(b),df1,df2)*exp(b)))))
  c*x*df(exp(x),df1=df1,df2=df2)*exp(x)
}

etrunclogf_denom = function(x,df1,df2,a,b){
  c = 10^(-round(min(log10(df(exp(a),df1,df2)*exp(a)),
                     log10(df(exp(b),df1,df2)*exp(b)))))
  c*df(exp(x),df1=df1,df2=df2)*exp(x)
}

# x multiply by the density of truncated log-F distribution on (a,b) at x
xdtrunclogf = function(x,df1,df2,a,b){
  x*df(exp(x),df1=df1,df2=df2)*exp(x)/(pf(exp(b),df1,df2)-pf(exp(a),df1,df2))
}

etrunclogf = function(df1,df2,a,b,adj=FALSE){
  if (adj==TRUE){
    n = integrate(etrunclogf_num, lower=a, upper=b, df1=df1,df2=df2,a=a,b=b)$value
    d = integrate(etrunclogf_denom, lower=a, upper=b, df1=df1,df2=df2,a=a,b=b)$value
    return(n/d)
  }else{
    return(integrate(xdtrunclogf, lower=a, upper=b, df1=df1,df2=df2,a=a,b=b)$value)
  } 
}

comp_sd = function(m){
  (m$b-m$a)/sqrt(12)
}

compute_lfsr = function(NegativeProb,ZeroProb){
  ifelse(NegativeProb> 0.5*(1-ZeroProb),1-NegativeProb,NegativeProb+ZeroProb)
}

compute_lfsra = function(PositiveProb, NegativeProb,ZeroProb){
  ifelse(PositiveProb<NegativeProb,2*PositiveProb+ZeroProb,2*NegativeProb+ZeroProb)  
}  

# automatically choose a grid of sd's for g
autoselect.mixsd = function(logfhat,df1,df2,mult){
  sigmaamin = log(qf(0.85,df1,df2))/10 #so that the minimum is small compared with measurement precision
  if(all(logfhat^2<=log(qf(0.85,df1,df2))^2)){
    sigmaamax = 8*sigmaamin #to deal with the occassional odd case where this could happen; 8 is arbitrary
  }else{
    sigmaamax = 2*sqrt(max(logfhat^2-log(qf(0.85,df1,df2))^2)) #this computes a rough largest value you'd want to use, based on idea that sigmaamax^2 + sebetahat^2 should be at least betahat^2   
  }
  if(mult==0){
    return(c(0,sigmaamax/2))
  }else{
    npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult))
    return(mult^((-npoint):0) * sigmaamax)
  }
}

# set prior
setprior=function(prior,k,nullweight,null.comp){
  if(!is.numeric(prior)){
    if(prior=="nullbiased"){ # set up prior to favour "null"
      prior = rep(1,k)
      prior[null.comp] = nullweight #prior 10-1 in favour of null by default
    }else if(prior=="uniform"){
      prior = rep(1,k)
    }
  }
  if(length(prior)!=k | !is.numeric(prior)){
    stop("invalid prior specification")
  }
  return(prior)
}

# initialize pi
initpi = function(k,n,null.comp,randomstart){
  if(randomstart){
    pi = rgamma(k,1,1)
  } else {
    if(k<n){
      pi=rep(1,k)/n #default initialization strongly favours null; puts weight 1/n on everything except null
      pi[null.comp] = (n-k+1)/n #the motivation is data can quickly drive away from null, but tend to drive only slowly toward null.
    } else {
      pi=rep(1,k)/k
    }
  }    
  pi=normalize(pi)
  return(pi)
}

