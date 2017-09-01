

## A mechanism for feature/ variable selection using ash-shrinkage in classtpx
## models


thetaSelect <- function(counts,
                        known_samples,
                        class_labs,
                        shrink=TRUE,
                        shrink.method=c(1,2),
                        nchunks=20,
                        trim= 0,
                        mash_user=NULL)
{
  if(!is.null(mash_user)){
    return(mash_user)
  }else{

  counts_class <- counts[known_samples,];
  if(!shrink){

    mean_features <- apply(counts_class, 2, function(x) return(mean(x, trim = trim)));
    
    FeatureSummary_class <- parallel::mclapply(1:dim(counts_class)[2],
                                               function(l) {
                                                 sd_element <- tapply(counts_class[,l], class_labs, sd);
                                                 central_element <- tapply(counts_class[,l], class_labs, function(x) return(mean(x, trim = trim)));
                                                 beta_element <- central_element - mean_features[l];
                                                 n.element <- as.numeric(table(class_labs));
                                                 sebeta_element <- sd_element/sqrt(n.element);
                                                 ll <- list("central_element"=central_element, "sd_element"=sd_element, "beta_element"=beta_element, "sebeta_element"=sebeta_element);
                                                 return(ll)
                                               })

    mean_class <- do.call(rbind, lapply(1:dim(counts_class)[2], function(l)
    {
      return(FeatureSummary_class[[l]]$central_element)
    }))

    theta_class <- class.normalizetpx(mean_class+1e-20, byrow=FALSE)
    return(theta_class)
  }

  if(shrink & shrink.method==2){

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
                                                                                  mixcompdist="normal")$result$PosteriorMean))
                                              }
                                              }));
    ash_mean_class <- ash_beta_class + mean_features;

    scale_clus <- ash_mean_class[,2]/ash_mean_class[,1];

#    chunks <- chunk(1:nrow(ash_mean_class), 20);

#    ash_theta_class_list <- parallel::mclapply(1:length(chunks),
#                                               function(l)
#                                               {
#                                                 out <- class.scallio(ash_mean_class[chunks[[l]],])$theta_class;
#                                                 return(out)
#                                               }, mc.cores=parallel::detectCores());

#    ash_theta_class <- matrix(0, nrow(ash_mean_class), ncol(ash_mean_class));



#    for(l in 1:length(chunks)){
#      ash_theta_class[chunks[[l]],] <- as.matrix(ash_theta_class_list[[l]])
#    }

#    scale_clus_scallio <- ash_theta_class[,2]/ash_theta_class[,1]
#      ash_theta_class <- class.normalizetpx(ash_theta_class, byrow=FALSE)
      ash_theta_class <- class.normalizetpx(ash_mean_class+1e-20, byrow=FALSE)
    return(ash_theta_class)
  }

  if(shrink & shrink.method==1){
    voom_class <- voom2(counts_class);
    mean_voom_features <- apply(voom_class, 2, mean);
    voom_class_adj <- voom_class - rep.row(mean_voom_features, dim(voom_class)[1])
    model_mat <- model.matrix(~as.factor(class_labs)-1)

    beta_class <- matrix(0, dim(voom_class)[2], dim(model_mat)[2]);
    sebeta_class <- matrix(0, dim(voom_class)[2], dim(model_mat)[2])
    for(k in 1:dim(model_mat)[2]){
      model_mat_temp <- cbind(model_mat[,k]);
      limma.obj <- limma::lmFit(t(voom_class_adj), design=matrix(model_mat_temp, nrow=dim(model_mat)[1]),
                                weights=limma::voom(t(counts_class))$weights)
  #    limma.obj <- limma::eBayes(limma.obj)
   #   mean_genes_limma <- apply(limma.obj$coefficients, 1, mean)
      beta_class[,k] <- as.matrix(limma.obj$coefficients[,1]);
      sebeta_class[,k] <- limma.obj$sigma*(as.matrix(limma.obj$stdev.unscaled[,1]));
    }

    ash_beta_class <- do.call(cbind, lapply(1:length(unique(class_labs)),
                                            function(l)
                                            {
                                              if(length(which(class_labs==l))==1){
                                                return(beta_class[,l])
                                              }else{
                                                return(suppressWarnings(ashr::ash(beta_class[,l], sebeta_class[,l],
                                                                                  mixcompdist="normal")$result$PosteriorMean))
                                              }
                                            }));

 #   voom_shrunk_class <- matrix(0, dim(counts_class)[1], dim(counts_class)[2])

#    for(i in 1:length(unique(class_labs))){
#      voom_shrunk_class[which(class_labs==unique(class_labs)[i]),] <-
#        voom_class[which(class_labs==unique(class_labs)[i]),] - rep.row(as.vector(beta_class[,i]), length(which(class_labs==unique(class_labs)[i])))
#      + rep.row(as.vector(ash_beta_class[,i]), length(which(class_labs==unique(class_labs)[i])));
#    }


 #   voom_shrunk_mean <- mean_voom_features + cbind.data.frame(beta_class);
#    voom_shrunk_class_2 <- matrix(0, dim(counts_class)[1], dim(counts_class)[2])

  voom_shrunk_class_2 <- matrix(0, dim(counts_class)[1], dim(counts_class)[2])

  for(i in 1:length(unique(class_labs))){
    voom_shrunk_class_2[which(class_labs==unique(class_labs)[i]),] <-
             voom_class[which(class_labs==unique(class_labs)[i]),] - rep.row(as.vector(beta_class[,i]), length(which(class_labs==unique(class_labs)[i]))) + rep.row(as.vector(ash_beta_class[,i]), length(which(class_labs==unique(class_labs)[i])));
  }


    lib_size_1 <- rowSums(counts_class);
    lib_size <- rep(mean(lib_size_1), length(lib_size_1));
    counts_shrunk_matrix <- (2^{voom_shrunk_class_2 - 6*log(10, base=2)})*(rep.col(lib_size_1+1, dim(voom_shrunk_class_2)[2])) - 0.5;
    counts_shrunk_matrix[counts_shrunk_matrix < 0]=1e-08;

    mean_counts_shrunk_class <- do.call(rbind, lapply(1:dim(counts_shrunk_matrix)[2], function(l)
    {
      mean_element <- tapply(counts_shrunk_matrix[,l], class_labs, mean);
      return(mean_element)
    }))

    ash_theta_class <- class.normalizetpx(mean_counts_shrunk_class+1e-20, byrow=FALSE);

    return(ash_theta_class)
  }
  }
}




#     voom_mean_class_1 <- apply(voom_class[which(class_labs==1),],2, mean)
#     voom_class_adj <- voom_class - rep.row(voom_mean_class_1, dim(voom_class)[1])
#
#     ash_beta_class <- matrix(0, dim(voom_class)[2], dim(model_mat)[2])
#     voom_shrunk_class <- matrix(0, dim(voom_class)[2], dim(model_mat)[2])
#
#     #beta_class <- matrix(0, dim(voom_class)[2], dim(model_mat)[2])
#
#  #   for(k in 1:dim(model_mat)[2]){
#       limma.obj <- limma::lmFit(t(voom_class_adj), model_mat)
#       limma.obj <- limma::eBayes(limma.obj)
#       mean_genes_limma <- apply(limma.obj$coefficients, 1, mean)
#       beta_class <- as.matrix(limma.obj$coefficients[,-1]);
#       sebeta_class <- limma.obj$sigma*(as.matrix(limma.obj$stdev.unscaled[,-1]));
#
#   #        if(length(which(class_labs==k))==1){
#   #         ash_beta_class <- beta_class;
#   #       }else{
#   #         ash_beta_class[,k] <- suppressWarnings(ashr::ash(beta_class, sebeta_class,
#   #                                    mixcompdist="normal")$PosteriorMean)
#   #       }
#   #       voom_shrunk_class[,k] <- limma.obj$coefficients[,1] + ash_beta_class[,k];
#   #    }
#
#
#
#     #    sebeta_class <- do.call(rbind, lapply(1:dim(voom_class)[2], function(l)
#     #    {
#     #      sd_element <- tapply(voom_class[,l], class_labs, sd);
#     #      n.element <- as.numeric(table(class_labs));
#     #      return(sd_element/sqrt(n.element))
#     #    }))
#
#
#      ash_beta_class <- do.call(cbind, lapply(1:length(unique(class_labs)[-1]),
#                                              function(l)
#                                              {
#                                                if(length(which(class_labs==(l+1)))==1){
#                                                  return(beta_class[,l])
#                                                }else{
#                                                  return(suppressWarnings(ashr::ash(beta_class[,l], sebeta_class[,l],
#                                                                                    mixcompdist="normal")$PosteriorMean))
#                                                }
#                                              }));
# #
#
#     voom_shrunk_mean <- limma.obj$coefficients[,1] + voom_mean_class_1 + cbind.data.frame(rep(0,dim(ash_beta_class)[1]), ash_beta_class);
#
#     voom_shrunk_class <- matrix(0, dim(counts_class)[1], dim(counts_class)[2])
#
#     for(i in 1:length(unique(class_labs))){
#       voom_shrunk_class[which(class_labs==unique(class_labs)[i]),] <- rep.row(as.vector(voom_shrunk_mean[,i]), length(which(class_labs==unique(class_labs)[i])));
#     }
#     lib_size <- rowSums(counts_class);
#
#     counts_shrunk_matrix <- (2^{voom_shrunk_class - 6*log(10, base=2)})*(rep.col(lib_size+1, dim(voom_shrunk_class)[2])) - 0.5
#     counts_shrunk_matrix[counts_shrunk_matrix < 0]=1e-08;
#
#     mean_counts_shrunk_class <- do.call(rbind, lapply(1:dim(counts_shrunk_matrix)[2], function(l)
#     {
#       mean_element <- tapply(counts_shrunk_matrix[,l], class_labs, mean);
#       return(mean_element)
#     }))
#
#     ash_theta_class <- class.normalizetpx(mean_counts_shrunk_class+1e-20, byrow=FALSE);
#     return(ash_theta_class)
#


#    chunks <- chunk(1:nrow(mean_counts_shrunk_class), 200);

#    ash_theta_class_list <- parallel::mclapply(1:length(chunks),
#                                function(l)
#                                {
#                                  out <- class.scallio(mean_counts_shrunk_class[chunks[[l]],])$theta_class;
#                                  return(out)
#                                }, mc.cores=parallel::detectCores());

#    ash_theta_class <- matrix(0, nrow(mean_counts_shrunk_class), ncol(mean_counts_shrunk_class));

#     for(l in 1:length(chunks)){
#      ash_theta_class[chunks[[l]],] <- as.matrix(ash_theta_class_list[[l]])
#     }

#
# }
#}


chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

class.scallio <- function(mean_counts_class){
  theta_init <- class.normalizetpx(mean_counts_class, byrow=FALSE);
  counter <- 0
  flag <- 0
  l <- 1
  while(l <= dim(mean_counts_class)[2]){
      scale_clus <- mean_counts_class[,-l]/mean_counts_class[,l];
    #  scale_clus[scale_clus <1e-09] = 1e-09
    #  scale_clus[scale_clus >1e+09] = 1e+09

      A <- rbind(as.matrix(t(scale_clus)), rep(1,dim(as.matrix(scale_clus))[1]));
      B <- rep(1,nrow(A));
      G <- diag(1, ncol(A));
      H <- rep(1e-06,ncol(A));

      solve.obj <- limSolve::lsei(A = A, B = B, G = G, H = H, verbose = FALSE)

      ash_theta_class <- matrix(0, length(scale_clus), nrow(A));
      ash_theta_class[,-(l)] <- A[-(nrow(A)),]*solve.obj$X;
      ash_theta_class[,l] <- solve.obj$X;

      scale_est <- ash_theta_class[,-l]/ash_theta_class[,l]
      na_indices <- which(is.na(rowSums(as.matrix(scale_est))));

      if(length(na_indices)==0){
        ll <- list("theta_class"=class.normalizetpx(ash_theta_class+1e-20, byrow=FALSE), "counter"=counter)
        return(ll)
      }

      if(l==1 & length(na_indices) > 0){
        temp <- ash_theta_class;
        counter <- counter + 1;
      }

      if(l !=1 & length(na_indices) > 0){
        scale_prev <- temp[,-l]/temp[,l];
        index1 <- which(is.na(rowSums(as.matrix(scale_prev))));
        scale_curr <- ash_theta_class[-l]/ash_theta_class[,l];
        index2 <- which(is.na(rowSums(as.matrix(scale_curr))));
        if(length(setdiff(index2, index1))==0){
        ash_theta_class[setdiff(index2, index1),] <- temp[setdiff(index2, index1),];
        }
        scale_est <- ash_theta_class[,-l]/ash_theta_class[,l]
        na_indices <- which(is.na(rowSums(as.matrix(scale_est))));
        if(length(na_indices)==0){
          ll <- list("theta_class"=class.normalizetpx(ash_theta_class+1e-20, byrow=FALSE), "counter"=counter)
          return(ll)
        }else{
        temp <- ash_theta_class;
        counter <- counter +1;
        }}
  }

  if(l==dim(mean_counts_class)[2]){
    scale_est <- ash_theta_class[,-l]/ash_theta_class[,l];
    na_indices <- which(is.na(rowSums(scale_est)));
    ash_theta_class[na_indices,] <- theta_init[na_indices,]
    ll <- list("theta_class"=class.normalizetpx(ash_theta_class+1e-20, byrow=FALSE), "counter"=counter)
    return(ll)
  }
}

voom2 <- function(counts){
  libsize.mat <- rep.col(rowSums(counts), dim(counts)[2]);
  voom.out <- log((counts+0.5), base=2) - log((libsize.mat+1), base=2)+ 6* log(10, base=2);
  return(voom.out)
}
