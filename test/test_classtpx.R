

############  test classtpx   ####################

library(maptpx)
library(classtpx)
counts <- get(load("../singlecell-clustering/output/GTEX_heart.rda"))
A_est <- get(load("../singlecell-clustering/output/A.total.14.rda"))
topic_clust = classtpx::class_topics(
  t(counts),
  K = 14,
  #known_samples = A.heart.12,
  #class_labs = rep(dim(data_heart)[2],"heart"),
  optional_theta = A_est[,1:14],
  method = "theta.fix",
  tol = 0.01,
  shrink = FALSE
)



method = "theta.fix"
optional_theta <- A_est
shrink=FALSE
mash_user=NULL
shape=NULL
initopics=NULL
tol=0.1
bf=FALSE
kill=2
ord=TRUE 
verb=1
tmax=10000
wtol=10^(-4)
qn=100
grp=NULL
admix=TRUE
prior_omega = NULL
trim = 0
nonzero=FALSE
dcut=-10

known_samples=NULL
class_labs=NULL
K <- dim(optional_theta)[2]
