## Tools for manipulation of text count matrices ##

## converting count to frequency matrix
class.normalizetpx <- function(x, byrow=TRUE){
    if(byrow){ s <- slam::row_sums(x)
               s[s==0] <- 1
               return( x/s ) }
    else{
      s <- slam::col_sums(x)
      s[s==0] <- 1
      return(t(t(x)/s)) }
}

## converting a count/freq matrix to tfidf
stm_tfidf <- function(x){
  idf <- log( nrow(x) ) - log(slam::col_sums(x>0) + 1) 
  t( t(x) * idf )
}
    
## Dirichlet RNG
rdir <- function(n, alpha)
{
    x <- matrix(rgamma(length(alpha)*n,alpha),nrow=n,byrow=TRUE)
    return(t(x/rowSums(x))) }

class.voom_generator <- function(counts_class, 
                           class_labs,
                           doshrink=TRUE)
{
  voom_class <- voom2(counts_class);
  mean_voom_features <- apply(voom_class, 2, mean);
  voom_class_adj <- voom_class - rep.row(mean_voom_features, dim(voom_class)[1])
  model_mat <- model.matrix(~as.factor(class_labs)-1) 
  
  beta_class <- matrix(0, dim(voom_class)[2], dim(model_mat)[2]);
  sebeta_class <- matrix(0, dim(voom_class)[2], dim(model_mat)[2])
  for(k in 1:dim(model_mat)[2]){
    model_mat_temp <- cbind(model_mat[,k]);
    limma.obj <- limma::lmFit(t(voom_class_adj), model_mat_temp)
    limma.obj <- limma::eBayes(limma.obj)
    #   mean_genes_limma <- apply(limma.obj$coefficients, 1, mean)
    beta_class[,k] <- as.matrix(limma.obj$coefficients[,1]);
    sebeta_class[,k] <- limma.obj$sigma*(as.matrix(limma.obj$stdev.unscaled[,1]));
  }
  
  voom_mean_class <- t(mean_voom_features + cbind.data.frame(beta_class));
  
  if(doshrink==FALSE){
    ll <- list("voom_class"=voom_class,
               "voom_mean_class"=voom_mean_class)
    return(ll)
  }
  else{
        ash_beta_class <- do.call(cbind, lapply(1:length(unique(class_labs)), 
                                          function(l) 
                                          {
                                            if(length(which(class_labs==l))==1){
                                              return(beta_class[,l])
                                            }else{
                                              return(suppressWarnings(ashr::ash(beta_class[,l], sebeta_class[,l], 
                                                                                mixcompdist="normal")$PosteriorMean))
                                            }
                                          }));
  
      voom_shrunk_class_2 <- matrix(0, dim(counts_class)[1], dim(counts_class)[2])
  
      for(i in 1:length(unique(class_labs))){
          voom_shrunk_class_2[which(class_labs==unique(class_labs)[i]),] <- 
              voom_class[which(class_labs==unique(class_labs)[i]),] - rep.row(as.vector(beta_class[,i]), length(which(class_labs==unique(class_labs)[i]))) + rep.row(as.vector(ash_beta_class[,i]), length(which(class_labs==unique(class_labs)[i])));
        }
  
      voom_shrunk_mean_class <- t(mean_voom_features + cbind.data.frame(beta_class));
  
        ll <- list("voom_class"=voom_class,
                    "voom_shrunk_class"=voom_shrunk_class_2,
                    "voom_mean_class" = voom_mean_class,
                    "voom_shrunk_mean_class"=voom_shrunk_mean_class)
        return(ll)
  }
}

class.nonmodel_clust <- function(counts,
                                 known_samples,
                                 class_labs,
                                 method="svm",
                                 doshrink=FALSE)
{
  counts_class <- counts[known_samples,];
  out <- class.voom_generator(counts_class, class_labs, doshrink=doshrink)
  if(method=="svm"){
      if(doshrink==FALSE){
          data <- out$voom_class;
          mean_model <- out$voom_mean_class;
          sigma_model <- do.call(cbind, lapply(1:dim(data)[2], function(l) tapply(data[,l], class_labs, sd)))
    
          pooled_voom <- voom2(counts)
          test_voom <- pooled_voom[-(known_samples),];
    
          loglik_test <- matrix(0, dim(test_voom)[1], length(unique(class_labs)))
    
          for(k in 1:length(unique(class_labs))){
              loglik_test[,k] <- rowSums(t(apply(test_voom, 1, function(x) dnorm(x, mean_model[k,], sigma_model[k,], log=TRUE))))
          }
          return(loglik_test)
      }else{
            data <- out$voom_shrunk_class;
            mean_model <- out$voom_shrunk_mean_class;
            sigma_model <- do.call(cbind, lapply(1:dim(data)[2], function(l) tapply(data[,l], class_labs, sd)))
    
            pooled_voom <- voom2(counts)
            test_voom <- pooled_voom[-(known_samples),];
    
            loglik_test <- matrix(0, dim(test_voom)[1], length(unique(class_labs)))
    
            for(k in 1:length(unique(class_labs))){
                loglik_test[,k] <- rowSums(t(apply(test_voom, 1, function(x) dnorm(x, mean_model[k,], sigma_model[k,], log=TRUE))))
            }
          return(loglik_test)
      }
  }
}


class.model_clust <- function(counts, 
                              known_samples, 
                              class_labs,
                              dist= c("normal", "poisson", "negbinom"),
                              libadj=TRUE,
                              shrink=FALSE)
{
  counts_class <- counts[known_samples,]
  test_data <- counts[-(known_samples),];
  out <- class.voom_generator(counts_class, class_labs)
      if(dist=="normal"){
          if(shrink==FALSE){          
            data <- out$voom_class;
            mean_model <- out$voom_mean_class;
          }else{
            data <- out$voom_shrunk_class;
            mean_model <- out$voom_shrunk_mean_class;
          }
          sigma_model <- do.call(cbind, lapply(1:dim(data)[2], function(l) tapply(data[,l], class_labs, sd)))
    
          pooled_voom <- voom2(counts)
          test_voom <- pooled_voom[-(known_samples),];
    
          loglik_test <- matrix(0, dim(test_voom)[1], length(unique(class_labs)))
    
          for(k in 1:length(unique(class_labs))){
            loglik_test[,k] <- rowSums(t(apply(test_voom, 1, function(x) dnorm(x, mean_model[k,], sigma_model[k,], log=TRUE))))
            }
            return(loglik_test)
          }
  
      if(dist=="poisson"){
              lib_size_1 <- rowSums(counts[known_samples,]);
              lib_size <- rep(mean(lib_size_1), length(lib_size_1));
              if(shrink==FALSE){
                  if(libadj){
                       counts_shrunk_matrix <- (2^{out$voom_class - 6*log(10, base=2)})*(rep.col(lib_size+1, dim(out$voom_shrunk_mean_class)[2])) - 0.5;
                       counts_shrunk_matrix[counts_shrunk_matrix < 0]=1e-08;
                    }
                  if(!libadj){
                       counts_shrunk_matrix <- (2^{out$voom_class - 6*log(10, base=2)})*(rep.col(lib_size_1+1, dim(out$voom_shrunk_mean_class)[2])) - 0.5;
                       counts_shrunk_matrix[counts_shrunk_matrix < 0]=1e-08;
                  }
              }else{
                  if(libadj){
                       counts_shrunk_matrix <- (2^{out$voom_shrunk_class - 6*log(10, base=2)})*(rep.col(lib_size+1, dim(out$voom_shrunk_mean_class)[2])) - 0.5;
                       counts_shrunk_matrix[counts_shrunk_matrix < 0]=1e-08;
                  }
                  if(!libadj){
                       counts_shrunk_matrix <- (2^{out$voom_shrunk_class - 6*log(10, base=2)})*(rep.col(lib_size_1+1, dim(out$voom_shrunk_mean_class)[2])) - 0.5;
                       counts_shrunk_matrix[counts_shrunk_matrix < 0]=1e-08;
                  }
               }
              df <- data.frame(class_labs, counts_shrunk_matrix)
              model_lambda <- apply(counts_shrunk_matrix, 2, function(x) return (tapply(x, class_labs, mean)));
              cat("Loglik calc")
              loglik_test <- matrix(0, dim(test_data)[1], length(unique(class_labs)))
              for(k in 1:length(unique(class_labs))){
                loglik_test[,k] <- rowSums(t(apply(test_data, 1, function(x) dpois(x, as.numeric(model_lambda[k,]), log=TRUE))))
                }
                return(loglik_test)
            }
  
      if(dist=="negbinom"){
              lib_size_1 <- rowSums(counts[known_samples,]);
              lib_size <- rep(mean(lib_size_1), length(lib_size_1));
              if(shrink==FALSE){
                if(libadj){
                  counts_shrunk_matrix <- (2^{out$voom_class - 6*log(10, base=2)})*(rep.col(lib_size+1, dim(out$voom_shrunk_mean_class)[2])) - 0.5;
                  counts_shrunk_matrix[counts_shrunk_matrix < 0]=1e-08;
                }
                if(!libadj){
                  counts_shrunk_matrix <- (2^{out$voom_class - 6*log(10, base=2)})*(rep.col(lib_size_1+1, dim(out$voom_shrunk_mean_class)[2])) - 0.5;
                  counts_shrunk_matrix[counts_shrunk_matrix < 0]=1e-08;
                }
              }else{
                if(libadj){
                  counts_shrunk_matrix <- (2^{out$voom_shrunk_class - 6*log(10, base=2)})*(rep.col(lib_size+1, dim(out$voom_shrunk_mean_class)[2])) - 0.5;
                  counts_shrunk_matrix[counts_shrunk_matrix < 0]=1e-08;
                }
                if(!libadj){
                  counts_shrunk_matrix <- (2^{out$voom_shrunk_class - 6*log(10, base=2)})*(rep.col(lib_size_1+1, dim(out$voom_shrunk_mean_class)[2])) - 0.5;
                  counts_shrunk_matrix[counts_shrunk_matrix < 0]=1e-08;
                }
              }

              out_mean <- apply(counts_shrunk_matrix, 2, function(x) return (tapply(x, class_labs, mean)));
              out_var <- apply(counts_shrunk_matrix, 2, function(x) return (tapply(x, class_labs, var)));
    
              negben_p <- (out_mean/out_var);
              negben_p[negben_p > 1 ]=0.9999
              negben_p[negben_p < 0] = 0.0001
    
              negben_r <- ceiling(((negben_p)/(1-negben_p))*out_mean)
              
              cat("Loglik calc")
              loglik_test <- matrix(0, dim(test_data)[1], length(unique(class_labs)))
              for(k in 1:length(unique(class_labs))){
                    loglik_test[,k] <- rowSums(t(apply(test_data, 1, function(x) dnbinom(x, as.numeric(negben_r[k,]), as.numeric(negben_p[k,]), log=TRUE))))
              }
            return(loglik_test)
        }
}


chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

voom2 <- function(counts){
  libsize.mat <- rep.col(rowSums(counts), dim(counts)[2]);
  voom.out <- log((counts+0.5), base=2) - log((libsize.mat+1), base=2)+ 6* log(10, base=2);
  return(voom.out)
}
