
GhostKnockoff.prelim<-function(cor.G, M=5, method='sdp', max.size=500, corr_max=0.75, clusters=NULL, rep.index=NULL){
  n.G<-nrow(cor.G)
  Normal_50Studies<-matrix(rnorm(n.G*M*50),n.G*M,50)
  P.each<-matrix(0,n.G,n.G)

  eigen.fit<-eigen(cor.G)
  newEig <- ifelse(eigen.fit$values < 1e-5, 1e-5, eigen.fit$values)
  newMat <- eigen.fit$vectors %*% (newEig*t(eigen.fit$vectors))
  # normalize modified matrix eqn 6 from Brissette et al 2007
  newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))
  cor.G<-newMat
  cor.G<-as.matrix(forceSymmetric(cor.G))
  Sigma<-cor.G
  SigmaInv<-solve(Sigma)
  
  if(length(clusters)==0){
    #clustering to identify tightly linked variants, 
    #first apply hierarchical clustering to determine number of clusters
    corr_max<-corr_max
    Sigma.distance = as.dist(1 - abs(Sigma))
    if(ncol(Sigma)>1){
      fit = hclust(Sigma.distance, method="single")
      clusters = cutree(fit, h=1-corr_max)
    }else{clusters<-1}
  }

  if(method=='equi'){
    s<-create.solve_equi_M(Sigma,M=M)
    diag_s<-Matrix(diag(s,length(s)))
    clusters=NULL
  }
  if(method=='sdp'){
    s<-create.solve_sdp_M(Sigma,M=M)
    diag_s<-Matrix(diag(s,length(s)))
    clusters=NULL
  }
  if(method=='mvr'){ 
    s<-solve_MVR_s_with_julia(Sigma,M=M)
    diag_s<-Matrix(diag(s,length(s)))
    clusters=NULL
  }
  if(method=='me'){ 
    s<-solve_ME_s_with_julia(Sigma,M=M)
    diag_s<-Matrix(diag(s,length(s)))
    clusters=NULL
  }
  if(method=='asdp'){
    s<-create.solve_asdp_M(Sigma,M=M,max.size=max.size)
    diag_s<-Matrix(diag(s,length(s)))
    clusters=NULL
  }
  if(method=='group.equi'){
    temp.fit<-create.solve_group_equi_M(Sigma,M=M,clusters=clusters)
    diag_s<-temp.fit$S
    clusters<-temp.fit$clusters
  }
  if(method=='group.sdp'){
    temp.fit<-create.solve_group_sdp_M(Sigma,M=M,clusters=clusters)
    diag_s<-temp.fit$S
    clusters<-temp.fit$clusters
  }
  if(method=='group.sdp.pca'){
    temp.fit<-create.solve_group_sdp_pca_M(Sigma,M=M,clusters=clusters)
    diag_s<-temp.fit$S
    clusters<-temp.fit$clusters
  }
  if(method=='group.me.pca'){
    temp.fit<-create.solve_group_ME_pca_M(Sigma,M=M,clusters=clusters)
    diag_s<-temp.fit$S
    clusters<-temp.fit$clusters
  }
  if(method=='group.sdp.graph'){
    if(length(rep.index)==0){
      rep.index<-c()
      for(k in 1:max(clusters)){
        if(sum(clusters==k)==1){rep.index<-c(rep.index,which(clusters==k))}else{
          rep.index<-c(rep.index,which(clusters==k)[Get.rep(Sigma[clusters==k,clusters==k,drop=T],thres=0.9,method="R2.ratio")])
        }
      }
    }
    diag_s<-create.solve_group_graphSDP_M(Sigma,rep.index,M=M,clusters=clusters)$S
  }
  if(method=='group.me.graph'){
    if(length(rep.index)==0){
      rep.index<-c()
      for(k in 1:max(clusters)){
        if(sum(clusters==k)==1){rep.index<-c(rep.index,which(clusters==k))}else{
          rep.index<-c(rep.index,which(clusters==k)[Get.rep(Sigma[clusters==k,clusters==k,drop=T],thres=0.9,method="R2.ratio")])
        }
      }
    }
    diag_s<-create.solve_group_graphME_M(Sigma,rep.index,M=M,clusters=clusters)$S
  }

  if(sum(diag_s)==0){
    V.left<-matrix(0,n.G*M,n.G*M)
  }else{
    #Sigma_k<-2*diag_s - s*t(s*SigmaInv)
    Sigma_k<-2*diag_s - diag_s%*%SigmaInv%*%diag_s

    V.each<-Matrix(forceSymmetric(Sigma_k-diag_s))

    #random part of knockoff
    V<-matrix(1,M,M)%x%V.each+diag(1,M)%x%diag_s
    V<-Matrix(forceSymmetric(V))
    
    #diag(V)<-diag(V)+rep(s,M)
    V.left<-try(t(chol(V)),silent=T)
    if(class(V.left)=="try-error"){
      eigen.fit<-eigen(V)
      newEig <- ifelse(eigen.fit$values < 1e-5, 1e-5, eigen.fit$values)
      newMat <- eigen.fit$vectors %*% (newEig*t(eigen.fit$vectors))
      # normalize modified matrix eqn 6 from Brissette et al 2007
      newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))
      V<-newMat
      V.left<-t(chol(V))
    }
  }
  #P.each<-diag(1,length(s))-s*SigmaInv
  P.each<-diag(1,n.G)-diag_s%*%SigmaInv
  Normal_50Studies<-as.matrix(V.left%*%matrix(rnorm(ncol(V.left)*50),ncol(V.left),50))

  #is_posdef((M+1)/M*Sigma-diag_s)
  A.each<-Sigma-diag_s
  A<-matrix(1,M+1,M+1)%x%A.each+diag(1,M+1)%x%diag_s
  #is_posdef(A)
  A<-Matrix(forceSymmetric(A))
  
  
  A.left<-try(t(chol(A)),silent=T)
  if(class(A.left)[1]=="try-error"){
    #alpha<-0.05-min(eigen(A)$values)
    #A<-(1-alpha)*A+alpha*diag(nrow(A))
    eigen.fit<-eigen(A)
    newEig <- ifelse(eigen.fit$values < max(1e-5,2*min(abs(eigen.fit$values))), max(1e-5,2*min(abs(eigen.fit$values))), eigen.fit$values)
    newMat <- eigen.fit$vectors %*% (newEig*t(eigen.fit$vectors))
    # normalize modified matrix eqn 6 from Brissette et al 2007
    newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))
    A<-newMat
    #A.left<-chol(newMat)
    #svd.fit<-svd.A
    #u<-svd.fit$u
    #svd.fit$d[is.na(svd.fit$d)]<-0
    #cump<-cumsum(svd.fit$d)/sum(svd.fit$d)
    #n.svd<-which(cump>=0.999)[1]
    #if(is.na(n.svd)){n.svd<-nrow(A)}
    #svd.index<-intersect(1:n.svd,which(svd.fit$d!=0))
    #A.left<-t(sqrt(svd.fit$d[svd.index])*t(u[,svd.index,drop=F]))
  }
  A.left<-as.matrix(t(chol(A)))
  #svd.A<-svd(A)
  #svd.A<-svd(A)

  return(list(P.each=as.matrix(P.each), V.left=V.left, Normal_50Studies=as.matrix(Normal_50Studies), M=M, A=A, A.left=A.left, clusters=clusters, diag_s=diag_s, Sigma=Sigma))
}


GhostKnockoff.fit<-function(Zscore_0, N.effect, fit.prelim, method='susie',type='fdr',M.fwer=50){
  Zscore_0<-as.matrix(Zscore_0)
  N.effect<-as.vector(N.effect)
  Zscore_0[is.na(Zscore_0)]<-0
  N.effect[is.na(N.effect)]<-Inf
  
  M<-fit.prelim$M
  n.G<-nrow(Zscore_0)
  P.each<-fit.prelim$P.each
  Normal_50Studies<-fit.prelim$Normal_50Studies
  A<-as.matrix(fit.prelim$A)
  A.left<-fit.prelim$A.left
  V.left<-fit.prelim$V.left

  if(type=='fdr'){M.rep<-1}
  if(type=='fwer'){M.rep<-M.fwer}
  
  T_0<-list();T_k<-list()
  kappa<-list();tau<-list()
  for(m in 1:M.rep){
    #Normal_k<-matrix(Normal_50Studies[,m],nrow=n.G)
    Normal_k<-matrix(V.left%*%matrix(rnorm(ncol(V.left)),ncol(V.left),1),nrow=n.G)
    GK.Zscore_0<-Zscore_0
    GK.Zscore_k<-as.vector(P.each%*%GK.Zscore_0)+Normal_k
    
    if(method=='marginal'){
      T_0[[m]]<-(GK.Zscore_0)^2
      T_k[[m]]<-(GK.Zscore_k)^2
    }
    if(method=='lasso'){
      #calculate importance score
      r<-GK.Zscore_0/sqrt(N.effect)#sqrt(N.effect-1+GK.Zscore_0^2)
      r_k<-as.vector(GK.Zscore_k/sqrt(N.effect))#sqrt(N.effect-1+GK.Zscore_k^2))
      r_all<-as.matrix(c(r,r_k))
      
      nfold<-5
      nA<-N.effect*(nfold-1)/nfold;nB<-N.effect/nfold
      temp.left<-sqrt(nB/nA/N.effect)*as.matrix(A.left)
      r_all_A<-r_all+as.matrix(temp.left%*%matrix(rnorm(ncol(temp.left)),ncol(temp.left),1))
      r_all_B<-(r_all*N.effect-r_all_A*nA)/nB
      shrink=0.01#seq(0.05,1,0.05)
      beta.all<-c();parameter.set<-c()
      k<-1
      #temp.A<-(1-shrink[k])*A+diag(shrink[k],nrow(A))
      temp.A<-A+diag(shrink[k],nrow(A))
      fit.basil<-try(ghostbasil(temp.A, r_all_A, alpha=1, delta.strong.size = max(1,min(500,length(r_all_A)/20)), max.strong.size = nrow(temp.A),n.threads=1,use.strong.rule=F),silent=T)
      parameter.set<-rbind(parameter.set,cbind(fit.basil$lmdas,shrink[k]))
      beta.all<-cbind(beta.all,fit.basil$betas)
      
      Get.f<-function(x){x<-as.matrix(x);return(t(x)%*%r_all_B/sqrt(t(x)%*%temp.A%*%x))}
      f.lambda<-apply(beta.all,2,Get.f)
      f.lambda[is.na(f.lambda)]<--Inf
      #beta<-beta.all[,which.max(f.lambda)]
      parameter<-parameter.set[which.max(f.lambda),]
      temp.A<-(1-parameter[2])*A+diag(parameter[2],nrow(A))
      
      lambda.all<-fit.basil$lmdas
      lambda<-fit.basil$lmdas[which.max(f.lambda)]
      lambda.seq <- lambda.all[lambda.all > lambda]
      lambda.seq <- c(lambda.seq, lambda)

      fit.basil<-ghostbasil(temp.A, r_all,user.lambdas=lambda.seq, alpha=1, delta.strong.size = max(1,min(500,length(r_all_A)/20)), max.strong.size = nrow(temp.A),n.threads=1,use.strong.rule=F)
      beta<-fit.basil$betas[,ncol(fit.basil$betas)]
      
      T_0[[m]]<-abs(beta[1:n.G])
      T_k[[m]]<-abs(matrix(beta[-(1:n.G)],n.G,M))
    }
    if(method=='lasso.approx.lambda'){
      #calculate importance score
      r<-GK.Zscore_0/sqrt(N.effect)#sqrt(N.effect-1+GK.Zscore_0^2)
      r_k<-as.vector(GK.Zscore_k/sqrt(N.effect))#sqrt(N.effect-1+GK.Zscore_k^2))
      r_all<-as.matrix(c(r,r_k))
      
      lambda_max<-max(abs(rnorm(length(r_all))))/sqrt(N.effect)
      epsilon <- .0001
      K <- 100
      lambda.all <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon),
                                  length.out = K)), digits = 10)
      lambda<-lambda_max*0.6
      lambda.seq <- lambda.all[lambda.all > lambda]
      lambda.seq <- c(lambda.seq, lambda)
      
      temp.A<-A+0.01*diag(1,nrow(A))
      fit.basil<-ghostbasil(temp.A, r_all,user.lambdas=lambda.seq, alpha=1, delta.strong.size = max(1,min(500,length(r_all)/20)), max.strong.size = nrow(temp.A),n.threads=1,use.strong.rule=F)
      beta<-fit.basil$betas[,ncol(fit.basil$betas)]
      
      T_0[[m]]<-abs(beta[1:n.G])
      T_k[[m]]<-abs(matrix(beta[-(1:n.G)],n.G,M))
    }
    if(method=='susie'){
      fitted_rss <- suppressMessages(susie_rss(z=c(GK.Zscore_0,as.vector(GK.Zscore_k)), R=as.matrix(A), L = min(10,length(GK.Zscore_0))))
      fitted_vars<-summary(fitted_rss)$vars
      beta<-fitted_vars[order(fitted_vars[,1]),2]#*c(GK.Zscore_0,as.vector(GK.Zscore_k))^2
      T_0[[m]]<-abs(beta[1:n.G])
      T_k[[m]]<-abs(matrix(beta[-(1:n.G)],n.G,M))
    }
    MK.stat<-MK.statistic(T_0[[m]],T_k[[m]])
    kappa[[m]]<-MK.stat[,'kappa']
    tau[[m]]<-MK.stat[,'tau']
  }

  return(list(T_0=T_0,T_k=T_k,kappa=kappa,tau=tau))
}

GhostKnockoff.fwer<-function(kappa,tau,M,M.fwer=50,alpha=0.05,beta=1){
  v.factor<-M;alpha.factor<-M/v.factor
  rho<-0.25
  v0<-1*v.factor; eta<- compute_eta_given_M(M.fwer,alpha*alpha.factor/rho,beta)
  ## Initialization
  pi = rep(0,length(kappa[[1]]))
  ## Multiple knockoff runs
  for (m in 1:M.fwer){
    v <- floor(v0)+rbinom(1,1,v0-floor(v0))
    order_tau <- order(tau[[m]],decreasing = TRUE)
    sorted_kappa <- kappa[[m]][order_tau]
    negid <- which(sorted_kappa>0)
    if(v>0){
      if(length(negid)<v){
        S <- which(kappa==0)
        pi[S] <- pi[S]+1
      }else{
        TT <- negid[v]
        S <- which(sorted_kappa[1:TT]==0)
        S <- order_tau[S]
        pi[S] <- pi[S]+1
      }
    }
  }
  pi <- pi/M.fwer
  S <- which(pi>=eta)
  q<-rep(1,length(kappa[[1]]))
  q[S]<-0
  
  return(list(S=S,q=q))
}


Get.transfer.Tau<-function(kappa,tau,U,reveal_prop=0.5){
  revealed_id<-which(tau<=quantile(tau,reveal_prop))
  unrevealed_id<-(1:length(tau))[-revealed_id]
  U<-as.matrix(U)
  transfer.tau<-rep(0,length(tau))
  for(k in 1:length(unrevealed_id)){
    temp.y<-as.numeric(kappa[revealed_id]!=0)
    if(var(temp.y)==0){
      add_id<-unrevealed_id[which.min(tau[unrevealed_id])]#sample(unrevealed_id,1)
      revealed_id<-c(revealed_id,add_id)
      transfer.tau[add_id]<-k
      unrevealed_id<-(1:length(tau))[-revealed_id]
    }else{
      temp.x<-tau[revealed_id]
      temp.z<-U[revealed_id,,drop=F]
      temp.fit<-lm(temp.y~cbind(temp.x,temp.z))
      temp.x<-tau[unrevealed_id]
      temp.z<-U[unrevealed_id,,drop=F]
      fitted.pval = predict(temp.fit,data.frame(cbind(temp.x,temp.z)),type = "response")
      
      revealed_id<-c(revealed_id,unrevealed_id[which.max(fitted.pval)])
      transfer.tau[unrevealed_id[which.max(fitted.pval)]]<-k
      unrevealed_id<-(1:length(tau))[-revealed_id]
    }
  }
  return(transfer.tau)
}

## optimization function for FWER control
get_lp_bnd <- function(M,eta,beta){
  one <- rep(1,M+1)
  one_eta <- rep(0,M+1)
  supp <- (0:M)/M
  one_eta[supp>=eta] <- 1
  
  y <- Variable(M+1)
  obj <- Maximize(t(one_eta)%*% y)
  constraints <- list(
    y>=0,
    t(supp)%*%y==1,
    diff(y)<=(beta-1)*y[-(M+1)]
  )
  prob <- Problem(obj, constraints)
  res <- psolve(prob,solver = "ECOS")
  y <- res$getValue(y)
  opt_ratio <- t(one_eta) %*% y
  return(opt_ratio)
}

compute_eta_given_M <- function(M,alpha,beta){
  v <- 1 ## the parameter v 
  etalist <- seq(0.05, 1, by = 0.01) ## candidates for eta
  peta <- length(etalist)
  alpha_mat <- rep(0,peta) ## the matrix to store the FWER bound with each pair of (M,eta)
  for(i in 1:peta){
    res <- get_lp_bnd(M,etalist[i],beta) #change to 1/2
    alpha_mat[i] <- res * v
  }
  M_eta<-etalist[min(which(alpha_mat<=alpha))]
  return(M_eta)
}


#GhostKnockoff.filter<-function (T_0,T_k,clusters=1:length(T_0)){
#  T_0<-as.matrix(T_0);T_k<-as.matrix(T_k)
#  M<-ncol(T_k);Rej.Bound<-10000
  
#  T.temp<-cbind(T_0,T_k)
#  T.temp[is.na(T.temp)]<-0
  
#  which.max.alt<-function(x){
#    temp.index<-which(x==max(x))
#    if(length(temp.index)!=1){return(temp.index[2])}else{return(temp.index[1])}
#  }
#  kappa<-apply(T.temp,1,which.max.alt)-1
  
#  Get.OtherMedian<-function(x){median(x[-which.max(x)])}
#  tau<-apply(T.temp,1,max)-apply(T.temp,1,Get.OtherMedian)
  
#  b<-order(tau,decreasing=T)
#  c_0<-kappa[b]==0
  #calculate ratios for top Rej.Bound tau values
#  ratio<-c();temp_0<-0
#  for(i in 1:length(b)){
    #if(i==1){temp_0=c_0[i]}
#    if(i==1){exist.cluster<-0}else{exist.cluster<-unique(clusters[b[1:(i-1)]])}
#    temp_0<-temp_0+c_0[i]*(!clusters[b[i]]%in%exist.cluster)
#    temp_1<-(length(exist.cluster)+!clusters[b[i]]%in%exist.cluster)-temp_0
#    temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
#    ratio<-c(ratio,temp_ratio)
#    if(i>Rej.Bound){break}
#  }
  #calculate q values for top Rej.Bound values
#  q<-rep(1,length(tau));
#  if(length(which(tau[b]>0))!=0){
#    index_bound<-max(which(tau[b]>0))
#    for(i in 1:length(b)){
#      temp.index<-i:min(length(b),Rej.Bound,index_bound)
#      if(length(temp.index)==0){next}
#      q[b[i]]<-min(ratio[temp.index])*c_0[i]+1-c_0[i]
#      if(i>Rej.Bound){break}
#    }
#    q[q>1]<-1
#  }
#  
#  return(list(kappa=kappa,tau=tau,q=q))
#}


GhostKnockoff.filter<-function (T_0,T_k,clusters=1:length(T_0)){
  T_0<-as.matrix(T_0);T_k<-as.matrix(T_k)
  M<-ncol(T_k);Rej.Bound<-10000
  
  T.temp<-cbind(T_0,T_k)
  T.temp[is.na(T.temp)]<-0
  
  which.max.alt<-function(x){
    temp.index<-which(x==max(x))
    if(length(temp.index)!=1){return(temp.index[2])}else{return(temp.index[1])}
  }
  kappa<-apply(T.temp,1,which.max.alt)-1
  
  Get.OtherMedian<-function(x){median(x[-which.max(x)])}
  tau<-apply(T.temp,1,max)-apply(T.temp,1,Get.OtherMedian)
  
  b<-order(tau,decreasing=T)
  c_0<-kappa[b]==0
  #calculate ratios for top Rej.Bound tau values
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    #if(i==1){temp_0=c_0[i]}
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    G.factor<-max(table(clusters[b][1:i]))
    temp_ratio<-(1/M*G.factor+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  #calculate q values for top Rej.Bound values
  q<-rep(1,length(tau));
  if(length(which(tau[b]>0))!=0){
    index_bound<-max(which(tau[b]>0))
    for(i in 1:length(b)){
      temp.index<-i:min(length(b),Rej.Bound,index_bound)
      if(length(temp.index)==0){next}
      q[b[i]]<-min(ratio[temp.index])*c_0[i]+1-c_0[i]
      if(i>Rej.Bound){break}
    }
    q[q>1]<-1
  }
  
  return(list(kappa=kappa,tau=tau,q=q))
}



MK.q.byStat<-function (kappa,tau,M,clusters=1:length(kappa),Rej.Bound=10000){
  b<-order(tau,decreasing=T)
  c_0<-kappa[b]==0
  #calculate ratios for top Rej.Bound tau values
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    #if(i==1){temp_0=c_0[i]}
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    G.factor<-max(table(clusters[b][1:i]))
    temp_ratio<-(1/M*G.factor+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  #calculate q values for top Rej.Bound values
  q<-rep(1,length(tau));
  if(length(which(tau[b]>0))!=0){
    index_bound<-max(which(tau[b]>0))
    for(i in 1:length(b)){
      temp.index<-i:min(length(b),Rej.Bound,index_bound)
      if(length(temp.index)==0){next}
      q[b[i]]<-min(ratio[temp.index])*c_0[i]+1-c_0[i]
      if(i>Rej.Bound){break}
    }
    q[q>1]<-1
  }
  
  return(q)
}


create.solve_group_graphME_M <- function(Sigma, rep.index, M=1, clusters, gaptol=1e-6, maxit=1000, verbose=FALSE) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))
  # Convert the covariance matrix to a correlation matrix
  G = cov2cor(Sigma)
  p = dim(G)[1]
  
  # Check that the input matrix is positive-definite
  if (!is_posdef(G)) {
    warning('The covariance matrix is not positive-definite: knockoffs may not have power.', immediate.=T)
  }
  
  ini.clusters<-clusters[rep.index]
  ini.clusters<-match(ini.clusters,unique(ini.clusters))
  
  Sigma_rep<-Sigma[rep.index,rep.index,drop=F]
  inv.Sigma_rep<-solve(Sigma_rep)
  if(nrow(Sigma_rep)==1){S<-as.matrix(1)}else{
    S<-create.solve_group_ME_pca_M(Sigma_rep,M=M,clusters=ini.clusters)$S
  }
  
  diag_s<-matrix(0,nrow(Sigma),ncol(Sigma))
  diag_s[rep.index,rep.index]<-as.matrix(S)
  diag_s[rep.index,-rep.index]<-as.matrix(S%*%inv.Sigma_rep%*%Sigma[rep.index,-rep.index])
  diag_s[-rep.index,rep.index]<-as.matrix(t(diag_s[rep.index,-rep.index]))
  diag_s[-rep.index,-rep.index]<-as.matrix(Sigma[-rep.index,-rep.index]-
                                             Sigma[-rep.index,rep.index]%*%inv.Sigma_rep%*%Sigma[rep.index,-rep.index]+
                                             Sigma[-rep.index,rep.index]%*%inv.Sigma_rep%*%S%*%inv.Sigma_rep%*%Sigma[rep.index,-rep.index])
  #                                               Sigma[-rep.index,rep.index]%*%(inv.Sigma_rep-inv.Sigma_rep%*%S%*%inv.Sigma_rep)%*%Sigma[rep.index,-rep.index])
  diag_s<-Matrix(forceSymmetric(diag_s))
  
  S<-diag_s
  
  # Compensate for numerical errors (feasibility)
  # psd S
  psd = 0;
  s_eps = 1e-8;
  while (psd==0) {
    psd = is_posdef(as.matrix(S)*(1-s_eps)+diag(s_eps,nrow(S))) #change 2 to (M+1)/Ms
    if (!psd) {
      s_eps = s_eps*10
    }
  }
  S = S*(1-s_eps)+diag(s_eps,nrow(S))
  S<-Matrix(S)
  
  psd = 0;
  s_eps = 1e-8;
  while (psd==0) {
    psd = is_posdef((M+1)/M*G-S*(1-s_eps)) #change 2 to (M+1)/Ms
    if (!psd) {
      s_eps = s_eps*10
    }
  }
  S = S*(1-s_eps)
  S<-Matrix(S)
  
  # Scale back the results for a covariance matrix
  return(list(S=diag(Sigma)^(1/2)*t(diag(Sigma)^(1/2)*S),Sigma=diag(Sigma)^(1/2)*t(diag(Sigma)^(1/2)*G),clusters=clusters,obj.tr=sum(diag(S))))
}

create.solve_group_graphSDP_M <- function(Sigma, rep.index, M=1, clusters, gaptol=1e-6, maxit=1000, verbose=FALSE) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))
  # Convert the covariance matrix to a correlation matrix
  G = cov2cor(Sigma)
  p = dim(G)[1]
  
  # Check that the input matrix is positive-definite
  if (!is_posdef(G)) {
    warning('The covariance matrix is not positive-definite: knockoffs may not have power.', immediate.=T)
  }
  
  ini.clusters<-clusters[rep.index]
  ini.clusters<-match(ini.clusters,unique(ini.clusters))
  
  Sigma_rep<-Sigma[rep.index,rep.index]
  inv.Sigma_rep<-solve(Sigma_rep)
  S<-create.solve_group_sdp_pca_M(Sigma_rep,M=M,clusters=ini.clusters)$S
  
  diag_s<-matrix(0,nrow(Sigma),ncol(Sigma))
  diag_s[rep.index,rep.index]<-as.matrix(S)
  diag_s[rep.index,-rep.index]<-as.matrix(S%*%inv.Sigma_rep%*%Sigma[rep.index,-rep.index])
  diag_s[-rep.index,rep.index]<-as.matrix(t(diag_s[rep.index,-rep.index]))
  diag_s[-rep.index,-rep.index]<-as.matrix(Sigma[-rep.index,-rep.index]-
                                             Sigma[-rep.index,rep.index]%*%inv.Sigma_rep%*%Sigma[rep.index,-rep.index]+
                                             Sigma[-rep.index,rep.index]%*%inv.Sigma_rep%*%S%*%inv.Sigma_rep%*%Sigma[rep.index,-rep.index])
  #                                               Sigma[-rep.index,rep.index]%*%(inv.Sigma_rep-inv.Sigma_rep%*%S%*%inv.Sigma_rep)%*%Sigma[rep.index,-rep.index])
  diag_s<-Matrix(forceSymmetric(diag_s))
  
  S<-diag_s
  
  # Compensate for numerical errors (feasibility)
  # psd S
  psd = 0;
  s_eps = 1e-8;
  while (psd==0) {
    psd = is_posdef(as.matrix(S)*(1-s_eps)+diag(s_eps,nrow(S))) #change 2 to (M+1)/Ms
    if (!psd) {
      s_eps = s_eps*10
    }
  }
  S = S*(1-s_eps)+diag(s_eps,nrow(S))
  S<-Matrix(S)
  
  psd = 0;
  s_eps = 1e-8;
  while (psd==0) {
    psd = is_posdef((M+1)/M*G-S*(1-s_eps)) #change 2 to (M+1)/Ms
    if (!psd) {
      s_eps = s_eps*10
    }
  }
  S = S*(1-s_eps)
  S<-Matrix(S)
  
  # Scale back the results for a covariance matrix
  return(list(S=diag(Sigma)^(1/2)*t(diag(Sigma)^(1/2)*S),Sigma=diag(Sigma)^(1/2)*t(diag(Sigma)^(1/2)*G),clusters=clusters,obj.tr=sum(diag(S))))
}



solve_ME_s_with_julia <- function(Sigma, M=1) {
  # pass Sigma from R to Julia
  julia_assign("Sigma", Sigma)
  julia_assign("M", M)
  
  # solve for s vector in Julia (note: one must wrap `Symmetric()` keyword around Sigma)
  julia_command("s = solve_s(Symmetric(Sigma), :maxent, tol=0.000001, robust=false, m=Int(M));")
  
  # put s vector from Julia back to R
  s <- julia_eval("s")
  
  # return s vector
  return(s)
}

solve_MVR_s_with_julia <- function(Sigma, M=1) {
  # pass Sigma from R to Julia
  julia_assign("Sigma", Sigma)
  julia_assign("M", M)
  
  # solve for s vector in Julia (note: one must wrap `Symmetric()` keyword around Sigma)
  julia_command("s = solve_s(Symmetric(Sigma), :mvr, tol=0.000001, robust=false, m=Int(M));")
  
  # put s vector from Julia back to R
  s <- julia_eval("s")
  
  # return s vector
  return(s)
}

solve_group_ME_s_with_julia <- function(Sigma, M=1, groups) {
  # pass Sigma and group membership from R to Julia
  julia_assign("Sigma", Sigma)
  julia_assign("M", M)
  julia_assign("groups", groups)
  
  # solve for S matrix in Julia (note there are 2 outputs!!!)
  julia_command("S, gammas = solve_s_group(Sigma, groups, :maxent, niter=20, tol=0.005, robust=false, m=Int(M));")
  
  # put S matrix from Julia back to R
  S <- julia_eval("S")
  
  # return S vector
  return(S)
}

solve_group_MVR_s_with_julia <- function(Sigma, M=1, groups) {
  # pass Sigma and group membership from R to Julia
  julia_assign("Sigma", Sigma)
  julia_assign("M", M)
  julia_assign("groups", groups)
  
  # solve for S matrix in Julia (note there are 2 outputs!!!)
  julia_command("S, gammas = solve_s_group(Sigma, groups, :mvr, niter=20, tol=0.005, robust=false, m=Int(M));")
  
  # put S matrix from Julia back to R
  S <- julia_eval("S")
  
  # return S vector
  return(S)
}



create.solve_sdp_M <- function(Sigma, M=1, gaptol=1e-6, maxit=1000, verbose=FALSE) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))
  # Convert the covariance matrix to a correlation matrix
  G = cov2cor(Sigma)
  p = dim(G)[1]
  
  # Check that the input matrix is positive-definite
  if (!is_posdef(G)) {
    warning('The covariance matrix is not positive-definite: knockoffs may not have power.', immediate.=T)
  }
  
  # Convert problem for SCS
  
  # Linear constraints
  Cl1 = rep(0,p)
  Al1 = -Matrix::Diagonal(p)
  Cl2 = rep(1,p)
  Al2 = Matrix::Diagonal(p)
  
  # Positive-definite cone
  d_As = c(diag(p))
  As = Matrix::Diagonal(length(d_As), x=d_As)
  As = As[which(Matrix::rowSums(As) > 0),]
  Cs = c((M+1)/M*G) ##change from 2 to (M+1)/M
  
  # Assemble constraints and cones
  A = cbind(Al1,Al2,As)
  C = matrix(c(Cl1,Cl2,Cs),1)
  K=NULL
  K$s=p
  K$l=2*p #not sure if it should be changed - may be not as it is the dimention of the linear part.
  
  # Objective
  b = rep(1,p)
  
  # Solve SDP with Rdsdp
  OPTIONS=NULL
  OPTIONS$gaptol=gaptol
  OPTIONS$maxit=maxit
  OPTIONS$logsummary=0
  OPTIONS$outputstats=0
  OPTIONS$print=0
  if(verbose) cat("Solving SDP ... ")
  sol = Rdsdp::dsdp(A,b,C,K,OPTIONS)
  if(verbose) cat("done. \n")
  
  # Check whether the solution is feasible
  if( ! identical(sol$STATS$stype,"PDFeasible")) {
    warning('The SDP solver returned a non-feasible solution. Knockoffs may lose power.')
  }
  
  # Clip solution to correct numerical errors (domain)
  s = sol$y
  s[s<0]=0
  s[s>1]=1
  
  # Compensate for numerical errors (feasibility)
  if(verbose) cat("Verifying that the solution is correct ... ")
  psd = 0
  s_eps = 1e-8
  while ((psd==0) & (s_eps<=0.1)) {
    if (is_posdef((M+1)/M*G-diag(s*(1-s_eps),length(s)),tol=1e-9)) { ##change from 2 to (M+1)/M
      psd  = 1
    }
    else {
      s_eps = s_eps*10
    }
  }
  s = s*(1-s_eps)
  s[s<0]=0
  if(verbose) cat("done. \n")
  
  # Verify that the solution is correct
  if (all(s==0)) {
    warning('In creation of SDP knockoffs, procedure failed. Knockoffs will have no power.',immediate.=T)
  }
  
  # Scale back the results for a covariance matrix
  return(s*diag(Sigma))
}


create.solve_group_sdp_M <- function(Sigma, M=1, clusters, gaptol=1e-6, maxit=1000, verbose=FALSE) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))
  # Convert the covariance matrix to a correlation matrix
  G = cov2cor(Sigma)
  p = dim(G)[1]

  # Check that the input matrix is positive-definite
  if (!is_posdef(G)) {
    warning('The covariance matrix is not positive-definite: knockoffs may not have power.', immediate.=T)
  }

  #clustering to identify tightly linked variants
  #corr_max<-corr_max
  #Sigma.distance = as.dist(1 - abs(G))
  #if(ncol(G)>1){
  #  fit = hclust(Sigma.distance, method="single")
  #  clusters = cutree(fit, h=1-corr_max)
  #}else{clusters<-1}
  clusters.index<-match(unique(clusters),clusters)
  p.clusters<-max(clusters)
  cum.p.clusters<-cumsum(table(clusters))

  # Convert problem for SCS

  # Linear constraints
  Cl1 = rep(0,p.clusters)
  Al1 = -Matrix::Diagonal(p.clusters)
  Cl2 = rep(1,p.clusters)
  Al2 = Matrix::Diagonal(p.clusters)

  # Positive-definite cone
  #d_As = c(diag(p))
  #As = Matrix::Diagonal(length(d_As), x=d_As)
  #As = As[which(Matrix::rowSums(As) > 0),]
  As<-c()
  for(k in 1:max(clusters)){
    temp.G<-matrix(0,p,p)
    if(k==1){temp.index<-1:cum.p.clusters[k]}else{
      temp.index<-(cum.p.clusters[k-1]+1):cum.p.clusters[k]
    }
    temp.G[temp.index,temp.index]<-G[clusters==k,clusters==k,drop=F]
    As<-rbind(As,Matrix(c(temp.G),1,length(temp.G)))
  }
  G.clusters<-G[order(clusters),order(clusters)]
  Cs = c((M+1)/M*G.clusters) ##change from 2 to (M+1)/M

  # Assemble constraints and cones
  A = cbind(Al1,Al2,As)
  C = matrix(c(Cl1,Cl2,Cs),1)
  K=NULL
  K$s=p
  K$l=2*p.clusters #not sure if it should be changed - may be not as it is the dimention of the linear part.

  # Objective
  #b = rep(1,p.clusters)
  b <- table(clusters) #reweight clusters

  # Solve SDP with Rdsdp
  OPTIONS=NULL
  OPTIONS$gaptol=gaptol
  OPTIONS$maxit=maxit
  OPTIONS$logsummary=0
  OPTIONS$outputstats=0
  OPTIONS$print=0
  if(verbose) cat("Solving SDP ... ")
  sol = Rdsdp::dsdp(A,b,C,K,OPTIONS)
  if(verbose) cat("done. \n")

  # Check whether the solution is feasible
  if( ! identical(sol$STATS$stype,"PDFeasible")) {
    warning('The SDP solver returned a non-feasible solution. Knockoffs may lose power.')
  }

  # Clip solution to correct numerical errors (domain)
  s = sol$y
  s[s<0]=0
  s[s>1]=1

  S<-matrix(0,ncol(G),ncol(G))
  for(k in 1:max(clusters)){
    #print(s[k])
    temp.G<-G[clusters==k,clusters==k,drop=F]
    S[clusters==k,clusters==k]<-s[k]*temp.G
  }

  # Compensate for numerical errors (feasibility)
  psd = 0;
  s_eps = 1e-8;
  while (psd==0) {
    psd = is_posdef((M+1)/M*G-S*(1-s_eps)) #change 2 to (M+1)/Ms
    if (!psd) {
      s_eps = s_eps*10
    }
  }
  S = S*(1-s_eps)
  S<-Matrix(S)

  # Scale back the results for a covariance matrix
  return(list(S=diag(Sigma)^(1/2)*t(diag(Sigma)^(1/2)*S),clusters=clusters))
}


create.solve_group_sdp_pca_M <- function(Sigma, M=1, clusters, gaptol=1e-6, maxit=1000, verbose=FALSE) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))
  # Convert the covariance matrix to a correlation matrix
  G = cov2cor(Sigma)
  p = dim(G)[1]
  
  # Check that the input matrix is positive-definite
  if (!is_posdef(G)) {
    warning('The covariance matrix is not positive-definite: knockoffs may not have power.', immediate.=T)
  }
  
  #clustering to identify tightly linked variants
  clusters.index<-match(unique(clusters),clusters)
  p.clusters<-max(clusters)
  cum.p.clusters<-cumsum(table(clusters))
  
  n.clusters <- table(clusters) #reweight clusters
  cum.n.clusters<-c(0,cumsum(n.clusters))
  # Convert problem for SCS
  
  # Linear constraints
  Cl1 = rep(0,2*p)
  temp1<-Matrix::Diagonal(2*p)
  Al1 = -temp1
  
  Cl2 = rep(1,p)
  Al2<-matrix(0,p,p)
  As.B<-c() # Positive-definite cone
  #save block-wise pca
  m<-0;D<-list()
  for (k in 1:max(clusters)){
    temp.G<-matrix(0,p,p)
    if(k==1){temp.index<-1:cum.p.clusters[k]}else{
      temp.index<-(cum.p.clusters[k-1]+1):cum.p.clusters[k]
    }
    block.G<-G[clusters==k,clusters==k,drop=F]
    block.fit<-eigen(block.G)
    W_k<-c()
    for(l in 1:n.clusters[k]){
      m<-m+1
      D_kl<-block.fit$vectors[,l]%*%t(block.fit$vectors[,l])
      D[[m]]<-D_kl
      temp.G[temp.index,temp.index]<-D_kl
      As.B<-rbind(As.B,Matrix(c(temp.G),1,length(temp.G)))
      W_k<-rbind(W_k,diag(D_kl))
    }
    Al2[temp.index,temp.index]<-W_k
  }
  d_As.D = c(diag(p))
  As.D = Matrix::Diagonal(length(d_As.D), x=d_As.D)
  As.D = As.D[which(Matrix::rowSums(As.D) > 0),]
  As<-rbind(As.B,As.D)
  Al2<-rbind(Al2,diag(p))
  
  G.clusters<-G[order(clusters),order(clusters)]
  Cs = c((M+1)/M*G.clusters) ##change from 2 to (M+1)/M
  
  # Assemble constraints and cones
  A = cbind(Al1,Al2,As)
  C = matrix(c(Cl1,Cl2,Cs),1)
  K=NULL
  K$l=2*p +p #not sure if it should be changed - may be not as it is the dimention of the linear part.
  K$s=p
  
  # Objective
  b<-rep(1,2*p)
  #b<-c(rep(1/n.clusters,times=n.clusters),rep(1/n.clusters,times=n.clusters))
  
  # Solve SDP with Rdsdp
  OPTIONS=NULL
  OPTIONS$gaptol=gaptol
  OPTIONS$maxit=maxit
  OPTIONS$logsummary=0
  OPTIONS$outputstats=0
  OPTIONS$print=0
  if(verbose) cat("Solving SDP ... ")
  sol = Rdsdp::dsdp(A,b,C,K,OPTIONS)
  if(verbose) cat("done. \n")
  
  # Check whether the solution is feasible
  if( ! identical(sol$STATS$stype,"PDFeasible")) {
    warning('The SDP solver returned a non-feasible solution. Knockoffs may lose power.')
  }
  
  s = sol$y
  S<-matrix(0,ncol(G),ncol(G))
  for(k in 1:max(clusters)){
    block.G<-G[clusters==k,clusters==k,drop=F]
    block.fit<-eigen(block.G)
    for(l in 1:n.clusters[k]){
      m<-cum.n.clusters[k]+l
      D_kl<-D[[m]]
      if(l==1){D_k<-D_kl*s[m]}else{
        D_k<-D_k+D_kl*s[m]
      }
    }
    S[clusters==k,clusters==k]<-as.matrix(D_k+Matrix::Diagonal(n.clusters[k], x=s[(p+cum.n.clusters[k]+1):(p+cum.n.clusters[k+1])]))
  }
  
  # Clip solution to correct numerical errors (domain)
  diag(S)[diag(S)<0]=0
  diag(S)[diag(S)>1]=1
  
  # Compensate for numerical errors (feasibility)
  psd = 0;
  s_eps = 1e-8;
  while (psd==0) {
    psd = is_posdef((M+1)/M*G-S*(1-s_eps)) #change 2 to (M+1)/Ms
    if (!psd) {
      s_eps = s_eps*10
    }
  }
  S = S*(1-s_eps)
  S<-Matrix(S)
  
  # Scale back the results for a covariance matrix
  return(list(S=diag(Sigma)^(1/2)*t(diag(Sigma)^(1/2)*S),clusters=clusters))
}



create.solve_group_ME_pca_M <- function(Sigma, M=1, clusters, ini.S=NULL, pc.max=Inf, rep.max=Inf, rep.method="R2", gaptol=1e-6, maxit=1000, verbose=FALSE) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))
  # Convert the covariance matrix to a correlation matrix
  G = cov2cor(Sigma)
  p = dim(G)[1]
  
  # Check that the input matrix is positive-definite
  if (!is_posdef(G)) {
    warning('The covariance matrix is not positive-definite: knockoffs may not have power.', immediate.=T)
  }
  
  #clustering to identify tightly linked variants
  clusters.index<-match(unique(clusters),clusters)
  p.clusters<-max(clusters)
  n.clusters <- table(clusters) #reweight clusters
  
  # prepare predefined matrices
  V.pc<-c();V.rep<-c()
  for (j in 1:p.clusters){
    block.index<-which(clusters==j)
    block.G<-G[block.index,block.index,drop=F]
    #PCA
    if(pc.max>0){
      block.fit<-eigen(block.G)
      v<-matrix(0,p,min(n.clusters[j],pc.max))
      v[block.index,]<-block.fit$vectors[,1:min(n.clusters[j],pc.max)]
      V.pc<-cbind(V.pc,v)
    }
    if(rep.max>0){
      #ID
      if(rep.method=='ID'){
        A<-chol(block.G)
        block.fit<-rid(A,ncol(A),k=min(n.clusters[j],rep.max),rand=F,idx_only=T)
        v<-matrix(0,p,min(n.clusters[j],rep.max))
        v[block.index,]<-diag(1,nrow(block.G))[,block.fit$idx]
        V.rep<-cbind(V.rep,v)
      }
      #SubsetC
      if(rep.method=='R2'){
        block.fit<-subsetC(block.G,k=min(n.clusters[j],rep.max))
        v<-matrix(0,p,min(n.clusters[j],rep.max))
        v[block.index,]<-diag(1,nrow(block.G))[,block.fit$indices]
        V.rep<-cbind(V.rep,v)
      }
    }
  }
  V<-Matrix(cbind(V.pc,V.rep))
  temp<-apply(V,2,paste,collapse="")
  V<-V[,match(unique(temp),temp)]
  
  # initial matrix
  #alpha<-1e-8
  #Sigma<-(1-alpha)*Sigma+diag(alpha,p)
  eigen.fit<-eigen(Sigma)
  if(length(ini.S)!=0){S<-ini.S}else{
    #inv.Sigma<-solve(Sigma)
    #S<-diag(1/irlba(inv.Sigma,nv=1)$d,p)
    S<-diag(min(eigen.fit$values),p)
    #S<-diag(1e-8,p)
  }
  S<-Matrix(S)
  E<-Matrix((M+1)/M*Sigma-S)
  inv.E<-eigen.fit$vectors%*%as.matrix(1/((M+1)/M*eigen.fit$values-min(eigen.fit$values))*t(eigen.fit$vectors))#solve(E)
  inv.S<-Matrix(diag(1/S[1,1],p))#Matrix(diag(1/min(eigen.fit$values),p))
  obj<-sum(log((M+1)/M*eigen.fit$values-min(eigen.fit$values)))+M*p*log(min(eigen.fit$values))#log(det(E))+M*log(det(S))
  #obj<-0
  #print(obj)
  
  new.obj<-obj
  for(i in 1:maxit){
    change_obj<-0
    for (l in 1:ncol(V)){
      v<-V[,l,drop=F]
      alpha_S<-as.numeric(t(v)%*%inv.S%*%v)
      alpha_E<-as.numeric(t(v)%*%inv.E%*%v)
      delta_l<-max(min((M*alpha_S-alpha_E)/(M+1)/alpha_S/alpha_E,1/alpha_E),-1/alpha_S)
      S<-S+delta_l*v%*%t(v)
      
      change_obj_l<-log(1-delta_l*alpha_E+1e-8)+M*log(1+delta_l*alpha_S+1e-8)
      change_obj<-change_obj+change_obj_l
      
      temp<-t(v)%*%inv.E
      inv.E<-inv.E-t(temp)%*%((-delta_l^(-1)+alpha_E)^(-1)*temp)
      temp<-t(v)%*%inv.S
      inv.S<-inv.S-t(temp)%*%((delta_l^(-1)+alpha_S)^(-1)*temp)
    }
    new.obj<-new.obj+change_obj
    if(change_obj/abs(new.obj)<=0.01){break}
    print(new.obj)
  }
  
  # Compensate for numerical errors (feasibility)
  psd = 0;
  s_eps = 1e-8;
  while (psd==0) {
    psd = is_posdef((M+1)/M*G-S*(1-s_eps)) #change 2 to (M+1)/Ms
    if (!psd) {
      s_eps = s_eps*10
    }
  }
  S = S*(1-s_eps)
  S<-Matrix(S)
  
  # Scale back the results for a covariance matrix
  return(list(S=diag(Sigma)^(1/2)*t(diag(Sigma)^(1/2)*S),clusters=clusters,obj.tr=sum(diag(S)),obj.me=new.obj))
}


create.solve_equi_M <- function(Sigma,M=1) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))
  p = nrow(Sigma)
  tol = 1e-10
  # Convert the covariance matrix to a correlation matrix
  G = cov2cor(Sigma)

  # Check that the input matrix is positive-definite
  if (!is_posdef(G)) {
    stop('The covariance matrix is not positive-definite: cannot solve SDP',immediate.=T)
  }

  lambda_min = eigen(G, symmetric=T, only.values = T)$values[p]
  if (lambda_min<0) {
    stop('In creation of equi-correlated knockoffs, while computing the smallest eigenvalue of the
         covariance matrix. The covariance matrix is not positive-definite.')
  }

  s = rep(1, nrow(Sigma)) * min((M+1)/M*lambda_min, 1)

  # Compensate for numerical errors (feasibility)
  psd = 0;
  s_eps = 1e-8;
  while (psd==0) {
    psd = is_posdef((M+1)/M*G-diag(s*(1-s_eps),length(s)))
    if (!psd) {
      s_eps = s_eps*10
    }
  }
  s = s*(1-s_eps)

  # Scale back the results for a covariance matrix
  return(s*diag(Sigma))
}

create.solve_group_equi_M <- function(Sigma,M=1,clusters){
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))
  p = nrow(Sigma)
  tol = 1e-6
  # Convert the covariance matrix to a correlation matrix
  G = cov2cor(Sigma)
  # Check that the input matrix is positive-definite
  if (!is_posdef(G)) {
    stop('The covariance matrix is not positive-definite: cannot solve SDP',immediate.=T)
  }

  if (p>2) {
    #clustering to identify tightly linked variants
    converged=FALSE
    maxitr=1000
    while (!converged) {
      #temp.index<-cbind(1:length(clusters),clusters)
      #group.G<-G[temp.index[order(temp.index[,2]),1],temp.index[order(temp.index[,2]),1],drop=F]
      D<-matrix(0,ncol(G),ncol(G))
      for(k in 1:max(clusters)){
        temp.G<-G[clusters==k,clusters==k,drop=F]
        svd.temp.G<-svd(temp.G)
        temp.D<-svd.temp.G$u%*%((svd.temp.G$d^(-1/2))*t(svd.temp.G$v))
        D[clusters==k,clusters==k]<-temp.D
      }
      D<-Matrix(D)
      #group.D<-D[temp.index[order(temp.index[,2]),1],temp.index[order(temp.index[,2]),1],drop=F]
      #group.DGD<-group.D%*%group.G%*%group.D
      DGD<-D%*%G%*%D
      DGD<-as.matrix(forceSymmetric(DGD))
      
      lambda_min<-RSpectra::eigs(as.matrix(DGD), 1, which="SR", opts=list(retvec = FALSE, maxitr=maxitr, tol=1e-8))$values

      if (length(lambda_min)==1) {
        converged = TRUE
      } else {
        lambda_min = eigen(DGD, symmetric=T, only.values = T)$values[p]
        converged = TRUE
      }
      }
    } else {
      lambda_min = eigen(G, symmetric=T, only.values = T)$values[p]
    }

  if (lambda_min<0) {
    stop('In creation of equi-correlated knockoffs, while computing the smallest eigenvalue of the
         covariance matrix. The covariance matrix is not positive-definite.')
  }

  gamma<-min((M+1)/M*lambda_min, 1) #change 2 to M+1/M
  S<-matrix(0,ncol(G),ncol(G))
  for(k in 1:max(clusters)){
    temp.G<-G[clusters==k,clusters==k,drop=F]
    S[clusters==k,clusters==k]<-gamma*temp.G
  }

  # Compensate for numerical errors (feasibility)
  psd = 0;
  s_eps = 1e-8;
  while (psd==0) {
    psd = is_posdef((M+1)/M*G-S*(1-s_eps)) #change 2 to (M+1)/Ms
    if (!psd) {
      s_eps = s_eps*10
    }
  }
  S = S*(1-s_eps)

  S<-Matrix(S)
  # Scale back the results for a covariance matrix
  return(list(S=diag(Sigma)^(1/2)*t(diag(Sigma)^(1/2)*S),clusters=clusters))
}


create.solve_asdp_M <- function(Sigma, M=1, max.size=500, gaptol=1e-6, maxit=1000, verbose=FALSE) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))

  if(ncol(Sigma) <= max.size) return(create.solve_sdp_M(Sigma, M=M, gaptol=gaptol, maxit=maxit, verbose=verbose))

  # Approximate the covariance matrix as block diagonal
  if(verbose) cat(sprintf("Dividing the problem into subproblems of size <= %s ... ", max.size))
  cluster_sol = divide.sdp(Sigma, max.size=max.size)
  n.blocks = max(cluster_sol$clusters)
  if(verbose) cat("done. \n")

  # Solve the smaller SDPs corresponding to each block
  if(verbose) cat(sprintf("Solving %s smaller SDPs ... \n", n.blocks))
  s_asdp_list = list()
  if(verbose) pb <- utils::txtProgressBar(min = 0, max = n.blocks, style = 3)
  for(k in 1:n.blocks) {
    s_asdp_list[[k]] = create.solve_sdp_M(as.matrix(cluster_sol$subSigma[[k]]), M=M, gaptol=gaptol, maxit=maxit)
    if(verbose) utils::setTxtProgressBar(pb, k)
  }
  if(verbose) cat("\n")

  # Assemble the solutions into one vector of length p
  p = dim(Sigma)[1]
  idx_count = rep(1, n.blocks)
  s_asdp = rep(0,p)
  for( j in 1:p ){
    cluster_j = cluster_sol$clusters[j]
    s_asdp[j] = s_asdp_list[[cluster_j]][idx_count[cluster_j]]
    idx_count[cluster_j] = idx_count[cluster_j]+1
  }

  # Maximize the shrinkage factor
  if(verbose) cat(sprintf("Combinining the solutions of the %s smaller SDPs ... ", n.blocks))
  tol = 1e-9
  maxitr=1000
  gamma_range = c(seq(0,0.1,len=11)[-11],seq(0.1,1,len=10)) # change from 100 to 20 to make it accurate near 0 and scalable.
  #options(warn=-1)
  gamma_opt = gtools::binsearch( function(i) {
    G = (M+1)/M*Sigma - gamma_range[i]*diag(s_asdp)
    lambda_min = suppressWarnings(RSpectra::eigs(G, 1, which = "SR", opts = list(retvec = FALSE, maxitr=maxitr, tol=tol))$values)
    if (length(lambda_min)==0) {
      #lambda_min = 1  # Not converged
      # RSpectra::eigs did not converge. Using eigen instead."
      lambda_min = min(eigen(G)$values)
    }
    lambda_min
  }, range=c(1,length(gamma_range)) )
  s_asdp_scaled = gamma_range[min(gamma_opt$where)]*s_asdp
  options(warn=0)
  if(verbose) cat("done. \n")

  if(verbose) cat("Verifying that the solution is correct ... ")
  # Verify that the solution is correct
  if (!is_posdef((M+1)/M*Sigma-diag(s_asdp_scaled,length(s_asdp_scaled)))) {
    warning('In creation of approximate SDP knockoffs, procedure failed. Knockoffs will have no power.',immediate.=T)
    s_asdp_scaled = 0*s_asdp_scaled
  }
  if(verbose) cat("done. \n")

  # Return result
  s_asdp_scaled
}

#find representative variants per group
Get.group.rep<-function(Sigma,clusters,inv.Sigma=NULL,thres=0.75,search.method='subsetC',stop.method='R2.ratio'){
  if(length(inv.Sigma)==0 & stop.method=='R2.ratio'){inv.Sigma<-solve(Sigma)}
  rep.data<-c()
  for(j in 1:max(clusters)){
    #print(j)
    if(sum(clusters==j)==1){
      rep.data<-rbind(rep.data,cbind(j,which(clusters==j)))
    }else{
      cor.G<-Sigma[clusters==j,clusters==j]
      if(search.method=='ID'){
        #interpolative decomposition
        A<-chol(cor.G)
        temp.fit<-rid(A,ncol(A),rand=F,idx_only=T)
        index.all<-temp.fit$idx
      }
      if(search.method=='subsetC'){
        index.all<-subsetC(cor.G, k=nrow(cor.G), traceit=FALSE)$indices
      }
      index<-index.all[1]
      for(i in 1:(nrow(cor.G)-1)){
        #print(i)
        #for(i in 1:4){
        temp.A<-cor.G[index,index,drop=F]
        #pre-compute some matrices
        if(i==1){inv.A<-solve(temp.A)}
        B<-cor.G[index,(1:nrow(cor.G))[-index],drop=F]
        #representative residual R2
        R2.R<-colSums(B*inv.A%*%B)
        inv.AB<-inv.A%*%B
        
        if(stop.method=='R2.ratio'){
          #representative plus other groups R2
          index.O<-which(clusters!=j)
          index.OR<-c(which(clusters==j)[index],index.O)
          inv.A.OR<-inv.Sigma[index.OR,index.OR]-
            inv.Sigma[index.OR,-index.OR,drop=F]%*%solve(inv.Sigma[-index.OR,-index.OR])%*%t(inv.Sigma[index.OR,-index.OR,drop=F])
          B.OR<-Sigma[which(clusters==j)[-index],index.OR,drop=F]
          R2.OR<-rowSums(B.OR%*%inv.A.OR*B.OR)#diag(B%*%inv.A%*%t(B))
          
          if(mean(R2.R/R2.OR)>=thres){
            #print(mean(R2.R/R2.OR))
            break
          }
        }
        if(stop.method=='R2'){
          if(mean(R2.R)>=thres){
            #print(min(R2.R))
            break
          }
        }
        index.add<-index.all[i+1]
        b<-cor.G[index,index.add,drop=F]
        c<-cor.G[index.add,index.add,drop=F]
        #R<-as.numeric(solve(c-t(b)%*%inv.A%*%b))
        R<-solve(c-t(b)%*%inv.A%*%b)
        inv.Ab<-inv.A%*%b
        inv.A<-rbind(cbind(inv.A+inv.Ab%*%R%*%t(inv.Ab),-inv.Ab%*%R),cbind(-R%*%t(inv.Ab),R))
        #update results
        index<-c(index,index.add)
        
      }
      index<-which(clusters==j)[index]
      rep.data<-rbind(rep.data,cbind(j,as.numeric(index)))
    }
  }
  return(rep.data)
}



#find representative variants
Get.rep<-function(Sigma,all.R2,thres=0.1,search.method='subsetC',stop.method='R2.diff',clusters=1:nrow(Sigma)){
  cor.G<-Sigma
  if(search.method=='ID'){
    #interpolative decomposition
    A<-chol(cor.G)
    temp.fit<-rid(A,ncol(A),rand=F,idx_only=T)
    index.all<-temp.fit$idx
  }
  if(search.method=='subsetC'){
    index.all<-subsetC(cor.G, k=nrow(cor.G), traceit=FALSE)$indices
  }
  if(search.method=='greadyR2.diff'){
    k<-nrow(cor.G)
    k.group<-1:k;
    index.all<-which.max(colSums(abs(cor.G)))
    for(i in 1:(k-1)){
      temp.A<-cor.G[index.all,index.all,drop=F]
      if(i==1){inv.A<-solve(temp.A)}
      #pre-compute some matrices
      B<-G[index.all,(1:nrow(cor.G))[-index.all],drop=F]
      inv.R.all<-diag(cor.G)[(1:nrow(cor.G))[-index.all]]-colSums(B*inv.A%*%B)
      inv.AB<-inv.A%*%B
      #compute residual R2
      obj<-all.R2[-index.all]-(1-inv.R.all)
      index.add<-((1:nrow(cor.G))[-index.all])[which.min(obj)]
      
      b<-cor.G[index.all,index.add,drop=F]
      c<-cor.G[index.add,index.add,drop=F]
      R<-as.numeric(solve(c-t(b)%*%inv.A%*%b))
      inv.Ab<-inv.A%*%b
      inv.A<-rbind(cbind(inv.A+inv.Ab%*%t(inv.Ab)*R,-inv.Ab*R),c(-t(inv.Ab)*R,R))
      #update results
      index.all<-c(index.all,index.add)
    }
  }
  
  clusters.index.all<-unique(clusters[index.all])
  #index<-index.all[1]
  index<-which(clusters==clusters.index.all[1])
  #for(i in 1:(k-1)){
  for(i in 1:(length(clusters.index.all)-1)){
  #for(i in 1:4){
    temp.A<-cor.G[index,index,drop=F]
    #pre-compute some matrices
    if(i==1){inv.A<-solve(temp.A)}
    B<-cor.G[index,(1:nrow(cor.G))[-index],drop=F]
    #representative residual R2
    inv.R.all<-diag(cor.G)[(1:nrow(cor.G))[-index]]-colSums(B*inv.A%*%B)
    inv.AB<-inv.A%*%B
    if(stop.method=='R2.ratio'){
      if(min((1-inv.R.all)/all.R2[-index])>=1-thres){break}
    }
    if(stop.method=='R2.diff'){
      if(max(all.R2[-index]-(1-inv.R.all))<=thres){break}
    }
    index.add<-which(clusters==clusters.index.all[i+1])
    b<-cor.G[index,index.add,drop=F]
    c<-cor.G[index.add,index.add,drop=F]
    #R<-as.numeric(solve(c-t(b)%*%inv.A%*%b))
    R<-solve(c-t(b)%*%inv.A%*%b)
    inv.Ab<-inv.A%*%b
    inv.A<-rbind(cbind(inv.A+inv.Ab%*%R%*%t(inv.Ab),-inv.Ab%*%R),cbind(-R%*%t(inv.Ab),R))
    #update results
    index<-c(index,index.add)
  }
  rep.index<-as.numeric(index)
  
  return(rep.index)
}



# 2022 Dec. 28 by Zihuai He, to better construct blocks

divide.block<-function(Sigma,max.size,ini.clusters){
  if(max.size<max(table(ini.clusters))){
    warning('Maximum block size is less than the maximum cluster size; algorithm may not converge', immediate.=T)
  }
  p<-nrow(Sigma)
  blocks<-rep(0,p)
  index<-1:p
  for(j in 1:p){
    Sigma.distance = as.dist(1 - abs(Sigma[index,index]))
    fit = hclust(Sigma.distance, method="average")
    for(temp.corr_max in unique(c(seq(0,0.1,by=0.01),seq(0.1,0.9,by=0.05)))){
      temp.clusters = cutree(fit, h=1-temp.corr_max)
      include.index<-ini.clusters %in% unique(ini.clusters[index[which(temp.clusters == which.max(table(temp.clusters)))]])
      #include.index<-index[which(temp.clusters == which.max(table(temp.clusters)))]
      if(sum(include.index)<=max.size){break}
      #print(max(table(temp.clusters)))
    }
    if(j==1){inter_AvgCorr_max<-temp.corr_max}
    #print(temp.corr_max)
    blocks[include.index]<-j
    index<-which(blocks==0)
    #print(sum(blocks!=0))
    if(sum(index)<=max.size){
      blocks[index]<-max(blocks)+1
      break
    }
  }
  return(list(blocks=blocks,inter_AvgCorr_max=inter_AvgCorr_max))
}


divide.sdp <- function(Sigma, max.size) {
  # Convert the covariance matrix into a dissimilarity matrix
  # Add a small perturbation to stabilize the clustering in the case of highly symmetrical matrices
  p = ncol(Sigma)
  Eps = matrix(rnorm(p*p),p)*1e-6
  dissimilarity = 1 - abs(cov2cor(Sigma)+Eps)
  distance = as.dist(dissimilarity)

  # Hierarchical clustering
  fit = hclust(distance, method="complete")
  # Cut tree into clusters of size smaller than a threshold
  n.blocks.min = 1
  n.blocks.max = ncol(Sigma)
  for(it in 1:100) {
    n.blocks = ceiling((n.blocks.min+n.blocks.max)/2)
    clusters = cutree(fit, k=n.blocks)
    size = max(table(clusters))
    if(size <= max.size) {
      n.blocks.max = n.blocks
    }
    if(size >= max.size) {
      n.blocks.min = n.blocks
    }
    if(n.blocks.min == n.blocks.max) {
      break
    }
  }

  # Merge small clusters
  clusters.new = merge.clusters(clusters, max.size)
  while(sum(clusters.new != clusters)>0) {
    clusters = clusters.new
    clusters.new = merge.clusters(clusters, max.size)
  }
  clusters = clusters.new

  # Create covariance submatrices for each cluster
  subSigma = vector("list", max(clusters))
  for( k in 1:length(subSigma) ) {
    indices_k = clusters==k
    subSigma[[k]] = Sigma[indices_k,indices_k]
  }

  # Return the cluster assignments and the cluster covariance submatrices
  structure(list(clusters=clusters, subSigma=subSigma), class='knockoff.clusteredCovariance')
}


merge.clusters <- function(clusters, max.size) {
  cluster.sizes = table(clusters)
  clusters.new = rep(0, length(clusters))
  g = 1
  g.size = 0
  for(k in 1:max(clusters)) {
    if(g.size + cluster.sizes[k] > max.size) {
      g = g + 1
      g.size = 0
    }
    clusters.new[clusters==k] = g
    g.size = g.size + cluster.sizes[k]
  }
  return(clusters.new)
}


is_posdef = function(A, tol=1e-9) {
  p = nrow(matrix(A))

  if (p<500) {
    lambda_min = min(eigen(A,only.values=T)$values)
  }
  else {
    oldw <- getOption("warn")
    #options(warn = -1)
    lambda_min = suppressWarnings(RSpectra::eigs(A, 1, which="SM", opts=list(retvec = FALSE, maxitr=100, tol))$values)
    options(warn = oldw)
    if( length(lambda_min)==0) {
      # RSpectra::eigs did not converge. Using eigen instead."
      lambda_min = min(eigen(A,only.values=T)$values)
    }
  }
  return (lambda_min>tol*10)
}


#' @title Elastic net using summary statistics
#' @description Coordinate descent algorithm to solve:
#' 0.5 x'X'Xx - x'b + lambda1 ||x||_1 + 0.5 lambda2 ||x||_2^2
#' Function to get elastic net solutions given X, a reference panel, and
#' b, regression coefficients
elnetR <- function(lambda1, lambda2=0, X, b, thr=1e-4,
                   trace=0, maxiter=10000,
                   blocks=NULL,
                   x=NULL) {
  stopifnot(length(b) == ncol(X))
  diag <- colSums(X^2)

  if(length(lambda2) > 1) {
    nlambda2 <- length(lambda2)
    for(i in 1:nlambda2) {
      result <- elnetR(lambda1, lambda2[i], X, b, thr,
                       trace, maxiter, x)
      result <- list(fit=result, lambda2=lambda2[i])
      if(i == 1) Result <- rep(result, nlambda2) else
        Result[i] <- result

    }
    return(Result)
  }

  order <- order(lambda1, decreasing = T)
  lambda1a <- lambda1[order]
  conv <- lambda1a * NA
  len <- length(b)
  beta <- matrix(NA, len, length(lambda1))
  pred <- matrix(NA, nrow(X), length(lambda1))
  loss <- rep(NA, length(lambda1))
  fbeta <- loss

  if(is.null(x)) x <- b * 0.0 else {
    stopifnot(length(x) == len)
    x <- x + 0.0 # Making sure R creates a copy...
  }

  if(is.null(blocks)) {
    Blocks <- list(startvec=0, endvec=len - 1)
  } else {
    Blocks <- parseblocks(blocks)
    stopifnot(max(Blocks$endvec)==len - 1)
  }

  X <- as.matrix(X)
  yhat <- as.vector(X %*% x)

  for(i in 1:length(lambda1a)) {
    if(trace > 0) cat("lambda1: ", lambda1a[i], "\n")
    conv[i] <- repelnet(lambda1a[i], lambda2, diag, X, b,thr,x,yhat, trace-1,maxiter,
                        Blocks$startvec, Blocks$endvec)
    if(conv[i] != 1) stop("Not converging...")

    beta[,i] <- x
    pred[,i] <- yhat
    loss[i] <- sum(yhat^2) - 2* sum(b * x)
    fbeta[i] <- loss[i] + 2* sum(abs(x))*lambda1a[i] + sum(x^2)*lambda2
  }

  conv[order] <- conv
  beta[,order] <- beta
  pred[,order] <- pred
  loss[order] <- loss
  fbeta[order] <- fbeta

  return(list(lambda1=lambda1, lambda2=lambda2, beta=beta, conv=conv, pred=pred, loss=loss, fbeta=fbeta))
}

repelnet <- function(lambda1, lambda2, diag, X, r, thr, x, yhat, trace, maxiter, startvec, endvec) {
  .Call('_lassosum_repelnet', PACKAGE = 'lassosum', lambda1, lambda2, diag, X, r, thr, x, yhat, trace, maxiter, startvec, endvec)
}


step1 <- function(C,vlist=seq(ncol(C)),RSS0=sum(diag(C)),zero=1e-12){
  dC <- diag(C)
  rs <- colSums(C^2)/dC
  imax <- order(rs,decreasing=TRUE)[1]
  vmin <- sum(dC) - rs[imax]
  residC = C - outer(C[,imax],C[,imax],"*")/C[imax,imax]
  index = vlist[imax]
  izero = diag(residC) <= zero
  list(index = index, variance = vmin, R2 = 1-vmin/RSS0, C=residC[!izero,!izero],vlist=vlist[!izero])
}

check_early_stopping_rule <- function(rsq_l, rsq_m, rsq_u, cond_0_thresh=1e-2, cond_1_thresh=1e-2) 
{
  delta_u <- (rsq_u-rsq_m)
  delta_m <- (rsq_m-rsq_l)
  (delta_u < cond_0_thresh*rsq_u) && ((delta_m*rsq_u-delta_u*rsq_m) < cond_1_thresh*rsq_m*rsq_u)
}



subsetC.R2 <- function(C, k=NA, R2.thres=0.9, traceit=FALSE){
  ## C correlation matrix
  ## k subset size
  do.adaptive <- is.na(k)
  p <- ncol(C)
  if (do.adaptive) {
    k <- p-1
  }
  indices <- rep(0, k)
  RSS0 <- p
  R2 <- double(k)
  vlist = seq(p)
  for(i in 1:k){
    fit1 <- step1(C, RSS0=RSS0, vlist=vlist)
    indices[i] <- fit1$index
    C <- as.matrix(fit1$C)
    vlist <- fit1$vlist
    R2[i] <- fit1$R2
    if(traceit)cat(i, "index", fit1$index, "Variance Explained", fit1$variance,"R-squared",fit1$R2,"\n")
    
    if(fit1$R2>=R2.thres){break}
  }
  list(indices = indices[indices!=0], R2=R2[indices!=0])
}


subsetC <- function(C, k=NA, traceit=FALSE){
  ## C correlation matrix
  ## k subset size
  do.adaptive <- is.na(k)
  p <- ncol(C)
  if (do.adaptive) {
    k <- p-1
  }
  indices <- rep(0, k)
  RSS0 <- p
  R2 <- double(k)
  vlist = seq(p)
  for(i in 1:k){
    fit1 <- step1(C, RSS0=RSS0, vlist=vlist)
    indices[i] <- fit1$index
    C <- as.matrix(fit1$C)
    vlist <- fit1$vlist
    R2[i] <- fit1$R2
    if(traceit)cat(i, "index", fit1$index, "Variance Explained", fit1$variance,"R-squared",fit1$R2,"\n")
    
    # if there is at least 3 R2 values,
    # check early stopping rule
    if (do.adaptive && (i >= 3)) {
      rsq_u <- R2[i]
      rsq_m <- R2[i-1]
      rsq_l <- R2[i-2]
      if (check_early_stopping_rule(rsq_l, rsq_m, rsq_u)) {
        indices <- indices[1:i]
        R2 <- R2[1:i]
        break
      }
    }
  }
  list(indices = indices, R2=R2)
}


#step1 <- function(C,vlist=seq(ncol(C)),RSS0=sum(diag(C)),zero=1e-12){
#  dC <- diag(C)
#  rs <- colSums(C^2)/dC
#  imax <- order(rs,decreasing=TRUE)[1]
#  vmin <- sum(dC) - rs[imax]
#  residC = C - outer(C[,imax],C[,imax],"*")/C[imax,imax]
#  index = vlist[imax]
#  izero = diag(residC) <= zero
#  list(index = index, variance = vmin, R2 = 1-vmin/RSS0, C=residC[!izero,!izero],vlist=vlist[!izero])
#}


#subsetC <- function(C, k, traceit=FALSE){
  ## C correlation matrix
  ## k subset size
#  indices <- rep(0, k)
#  p <- ncol(C)
#  RSS0 <- p
#  R2 <- double(k)
#  vlist = seq(p)
#  for(i in 1:k){
#    fit1 <- step1(C, RSS0=RSS0, vlist=vlist)
#    indices[i] <- fit1$index
#    C <- as.matrix(fit1$C)
#    vlist <- fit1$vlist
#    R2[i] <- fit1$R2
#    if(traceit)cat(i, "index", fit1$index, "Variance Explained", fit1$variance,"R-squared",fit1$R2,"\n")
#  }
#  list(indices = indices, R2=R2)
#}
