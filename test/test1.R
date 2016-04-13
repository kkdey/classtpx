counts_class <- counts[known_samples,];
voom_class <- limma::voom(counts_class)$E;
model.matrix(~as.factor(class_labs)-1);
limma.obj <- limma::lmFit(t(voom_class), model.matrix(~as.factor(class_labs)-1))
limma::eBayes(limma.obj)
lm(voom_class[,100] ~ as.factor(class_labs)-1)

global_mean <- apply(limma.obj$coefficients, 1, mean)

beta_class <- limma.obj$coefficients - global_mean;

sebeta_class <- do.call(rbind, lapply(1:dim(voom_class)[2], function(l)
                          {
                              sd_element <- tapply(voom_class[,l], class_labs, sd);
                              n.element <- as.numeric(table(class_labs));
                              return(sd_element/sqrt(n.element))
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

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

voom_shrunk_mean <- global_mean + ash_beta_class;
voom_shrunk <- matrix(0, dim(counts_class)[1], dim(counts_class)[2])

for(i in 1:length(unique(class_labs))){
  voom_shrunk[which(class_labs==unique(class_labs)[i]),] <- rep.row(as.vector(voom_shrunk_mean[,i]), length(which(class_labs==unique(class_labs)[i])));
}
lib_size <- rowSums(counts_class);

counts_shrunk <- (2^{voom_shrunk - 6*log(10, base=2)})*(rep.col(lib_size+1, dim(voom_shrunk)[2])) -0.5
counts_shrunk[counts_shrunk < 0]=1e-08;

mean_counts_class <- do.call(rbind, lapply(1:dim(counts_shrunk)[2], function(l)
{
  mean_element <- tapply(counts_shrunk[,l], class_labs, mean);
  return(mean_element)
}))

theta_init <- class.normalizetpx(mean_counts_class, byrow=FALSE);

scale_clus_1 <- mean_counts_class[,1]/mean_counts_class[,2];
scale_clus_1[scale_clus_1 <1e-05] = 1e-05
scale_clus_1[scale_clus_1 >1e+05] = 1e+05

A <- rbind(scale_clus_1, rep(1,length(scale_clus_1)));
B <- rep(1,2);
G <- diag(1, length(scale_clus_1));
H <- rep(0, length(scale_clus_1));

solve.obj <- limSolve::lsei(A = A, B = B, G = G, H = H, verbose = FALSE)

ash_theta_class_1 <- matrix(0, length(scale_clus_1), length(unique(class_labs)));
ash_theta_class_1[,2] <- solve.obj$X;
ash_theta_class_1[,1] <- solve.obj$X * scale_clus_1;

scale_clus_12 <- ash_theta_class_1[,1]/ash_theta_class_1[,2];


scale_clus_2 <- mean_counts_class[,2]/mean_counts_class[,1];
scale_clus_2[scale_clus_2 <1e-05] = 1e-05
scale_clus_2[scale_clus_2 >1e+05] = 1e+05

A <- rbind(scale_clus_2, rep(1,length(scale_clus_2)));
B <- rep(1,2);
G <- diag(1, length(scale_clus_2));
H <- rep(0, length(scale_clus_2));

solve.obj <- limSolve::lsei(A = A, B = B, G = G, H = H, verbose = FALSE)

ash_theta_class_2 <- matrix(0, length(scale_clus), length(unique(class_labs)));
ash_theta_class_2[,1] <- solve.obj$X;
ash_theta_class_2[,2] <- solve.obj$X * scale_clus_2;

scale_clus_21 <- ash_theta_class_2[,2]/ash_theta_class_2[,1];

flag_scale <- 0;

if(length(intersect(which(is.na(scale_clus_12)), which(is.na(scale_clus_21))))==0){
  ash_theta_class_1[which(is.na(scale_clus_12)),] <- ash_theta_class_2[which(is.na(scale_clus_12)),];
}else{
  ash_theta_class_1[which(is.na(scale_clus_12)),] <- theta_init[which(is.na(scale_clus_12)),];
  flag_scale <- 1;
  warning("NAs occured in limsolve under both directional scaling")
}




ash_theta_class_1[,2]/ash_theta_class_1[,1]





