

## A mechanism for feature/ variable selection using ash-shrinkage in classtpx 
## models


thetaSelect <- function(counts, known_samples, class_labs, shrink=TRUE)
{
  counts_class <- counts[known_samples,];
  mean_features <- apply(counts_class, 2, mean);
  FeatureSummary_class <- parallel::mclapply(1:dim(counts_class)[2], 
                                             function(l) {
                                               sd_element <- tapply(counts_class[,l], class_labs, sd);
                                               mean_element <- tapply(counts_class[,l], class_labs, mean);
                                               beta_element <- mean_element - mean_features[l];
                                               n.element <- as.numeric(table(class_labs));
                                               sebeta_element <- sd_element/sqrt(n.element);
                                               ll <- list("mean_element"=mean_element, "sd_element"=sd_element, "beta_element"=beta_element, "sebeta_element"=sebeta_element);
                                               return(ll)
                                             })
  
  if(!shrink){
    
    mean_class <- do.call(rbind, lapply(1:dim(counts_class)[2], function(l)
    {
      return(FeatureSummary_class[[l]]$mean_element)
    }))
    
    theta_class <- class.normalizetpx(mean_class+1e-20, byrow=FALSE)
    return(theta_class)
  }
  
  if(shrink){
    
    sebeta_class <- do.call(rbind, lapply(1:dim(counts_class)[2], function(l)
    {
      return(FeatureSummary_class[[l]]$sebeta_element)
    }))
    
    beta_class <- do.call(rbind, lapply(1:dim(counts_class)[2], function(l)
    {
      return(FeatureSummary_class[[l]]$beta_element)
    }))
    
    ash_beta_class <- do.call(cbind, lapply(unique(class_labs), 
                                            function(l) 
                                            {
                                              if(length(which(class_labs==l))==1){
                                                return(beta_class[,l])
                                              }else{
                                                return(suppressWarnings(ashr::ash(beta_class[,l], sebeta_class[,l], 
                                                                                  mixcompdist="normal")$PosteriorMean))
                                              }
                                              }));
    ash_mean_class <- ash_beta_class + mean_features;
    ash_theta_class <- class.normalizetpx(ash_mean_class+1e-20, byrow=FALSE)
    return(ash_theta_class)
  }
}