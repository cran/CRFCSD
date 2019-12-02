.quadrature.rules <- function( recurrences, inner.products )
{
  
  np1 <- nrow( recurrences )
  n <- np1 - 1
  rules <- as.list( rep( NULL, n ) )
  monic.recurrences <- orthopolynom::monic.polynomial.recurrences( recurrences )
  matrices <- orthopolynom::jacobi.matrices( monic.recurrences )
  matrix.eigens <- lapply( matrices, eigen )
  roots <- orthopolynom::polynomial.roots( monic.recurrences )
  h.0 <- inner.products[1]
  for ( k in 1:n ) {
    values <- matrix.eigens[[k]]$values
    vectors <- matrix.eigens[[k]]$vectors
    x <- values
    w <- rep( 0, k )
    for ( j in 1:k ) {
      v.j <- vectors[1,j]
      w[j] <- h.0 * v.j * v.j
    }
    rule <- data.frame( cbind( x, w ) )
    names( rule ) <- c( "x", "w" )
    rules[[k]] <- rule
  }
  return( rules )
}


.ghrules=function(n,normalized=FALSE){
  r=orthopolynom::hermite.h.recurrences(n,normalized)
  ip=orthopolynom::hermite.h.inner.products(n)
  return(.quadrature.rules(r,ip))
}
.EMiter=function(lastpar,X,Z,Delta,n,ni,blC,r,myrules,lambda,R,Cauchy.pen,clustering){
  betadim=dim(X[[1]])[2]
  gammadim=dim(Z)[2]
  zeta=lastpar[1:(1+betadim+gammadim)]
  betagamma=lastpar[(2+betadim+gammadim):length(lastpar)]
  beta=betagamma[1:betadim]
  gamma=betagamma[(betadim+1):(betadim+gammadim)]
  theta=betagamma[betadim+gammadim+1]
  psi=betagamma[(2+betadim+gammadim):length(betagamma)]
  if(clustering){
  normalC=NormalConstant(myrules,beta,gamma,zeta,theta,Delta,betadim,gammadim,n,ni,X,Z,psi,blC,r)
  
  g0bvalue=g0bmat(myrules,beta,gamma,zeta,theta,Delta,betadim,gammadim,n,ni,X,Z,psi,blC,r)
  }
  else{
    myrule[,1]=rep(0,dim(myrules)[1])
    myrule[,2]=rep(1,dim(myrules)[1])
    normalC=rep(dim(myrules)[1],n)
    g0bvalue=matrix(1, nrow = dim(myrules)[1], ncol = n)
  }
  uijqmat=uijq(myrules,beta,gamma,zeta,theta,Delta,betadim,gammadim,n,ni,X,Z,psi,blC,r)
  zeta=optim(par = zeta,fn=Q1,rules=myrules,beta=beta,gamma=gamma,eta=zeta,
             theta=theta,Delta=Delta,betadim=betadim,gammadim=gammadim,n=n,ni=ni,X=X,Z=Z,psi=psi,blC=blC,r=r,g0bvalue=g0bvalue,
             normalconstant=normalC,uijqmat=uijqmat,Cauchyindex=Cauchy.pen,method = "BFGS")$par
  betagamma=optim(par = betagamma,fn=Q2,rules=myrules,beta=beta,gamma=gamma,eta=zeta,
                  theta=theta,Delta=Delta,betadim=betadim,gammadim=gammadim,n=n,ni=ni,X=X,Z=as.matrix(Z),psi=psi,blC=blC,r=r,g0bvalue=g0bvalue,
                  normalconstant=normalC,uijqmat=uijqmat,tunepar=lambda,R=R,Cauchyindex=Cauchy.pen,method = "BFGS")$par
  result=c(zeta,betagamma)
  return(result)
}
.EMest=function(X,Z,Delta,n,ni,C,r,quadnum,lambda,criterion,blC,myrules,Cauchy.pen,clustering){
  
  coefnum=dim(blC[[1]])[1]
  D=matrix(0, nrow = coefnum-2, ncol = coefnum)
  for(i in 1:(coefnum-2)){
    D[i,i:(i+2)]=c(1,-2,1)
  }
  R=(t(D)%*%D)
  betadim=dim(X[[1]])[2]
  gammadim=dim(Z)[2]
  lastpar=rep(0,2+2*(betadim+gammadim)+coefnum)
  difference=1
  while (difference>criterion) {
    outputpar=.EMiter(lastpar,X,Z,Delta,n,ni,blC,r,myrules,lambda,R,Cauchy.pen,clustering)
    difference=sum(abs((-outputpar+lastpar)/(lastpar)))
    lastpar=outputpar
  }
  
  return(outputpar)
}
.AICcomputepenpara=function(parest,X,Z,Delta,blC,r,n,ni,myrules,lambda,Cauchy.pen,clustering){
  coefnum=dim(blC[[1]])[1]
  D=matrix(0, nrow = coefnum-2, ncol = coefnum)
  for(i in 1:(coefnum-2)){
    D[i,i:(i+2)]=c(1,-2,1)
  }
  R=(t(D)%*%D)
  pardim=length(parest)
  result=matrix(0, nrow = pardim, ncol = pardim)
  betadim=dim(X[[1]])[2]
  gammadim=dim(Z)[2]
  zeta=parest[1:(1+betadim+gammadim)]
  betagamma=parest[(2+betadim+gammadim):length(parest)]
  beta=betagamma[1:betadim]
  gamma=betagamma[(betadim+1):(betadim+gammadim)]
  theta=betagamma[betadim+gammadim+1]
  psi=betagamma[(2+betadim+gammadim):length(betagamma)]
  if(clustering){
  g0bvalue=g0bmat(myrules,beta,gamma,zeta,theta,Delta,betadim,gammadim,n,ni,X,Z,psi,blC,r)
  normalC=NormalConstant(myrules,beta,gamma,zeta,theta,Delta,betadim,gammadim,n,ni,X,Z,psi,blC,r)
  }
  else{
    myrule[,1]=rep(0,dim(myrules)[1])
    myrule[,2]=rep(1,dim(myrules)[1])
    normalC=rep(dim(myrules)[1],n)
    g0bvalue=matrix(1, nrow = dim(myrules)[1], ncol = n)
  }
  quadnum=dim(myrules)[1]
  for (i in 1:n) {
    resulti=matrix(0, nrow = pardim, ncol = pardim)
    for (k in 1:quadnum) {
      
      firstdrivk=firstderiv_i(parest,myrules[k,1],Delta[[i]],X[[i]],Z[i,],ni[i],r,blC[[i]],betadim,gammadim)
      if(sum(is.nan( firstdrivk)>0)){
        firstdrivk=numDeriv::grad(loglikelihoodtest,x=parest,b=myrules[k,1],Delta=Delta[[i]],X=X[[i]],Z=Z[i,],
                                  ni=ni[i],r=r,blC=blC[[i]],betadim=betadim,gammadim=gammadim)
      }
      for (q in 1:quadnum) {
        
        firstdrivq=firstderiv_i(parest,myrules[q,1],Delta[[i]],X[[i]],Z[i,],ni[i],r,blC[[i]],betadim,gammadim)
        if(sum(is.nan( firstdrivq)>0)){
          firstdrivq=numDeriv::grad(loglikelihoodtest,x=parest,b=myrules[q,1],Delta=Delta[[i]],X=X[[i]],Z=Z[i,],
                                    ni=ni[i],r=r,blC=blC[[i]],betadim=betadim,gammadim=gammadim)
        }
        tempresult=as.matrix(firstdrivk)%*%t(as.matrix(firstdrivq))*myrules[k,2]*g0bvalue[k,i]*myrules[q,2]*g0bvalue[q,i]/(normalC[i])^2
        
        resulti=resulti+tempresult
        if(sum(is.nan(tempresult))>0){
          print(c(i,k,q))
        }
      }
    }
    result=result+resulti
  }
  if(Cauchy.pen){
    hessian=numDeriv::hessian(logcauchy,x=parest)
    result=result-hessian
    if(lambda!=0){
      Rbig=matrix(0, nrow = length(parest), ncol = length(parest))
      Rsecderiv=numDeriv::hessian(penaltyterm,x=psi,lambda=lambda,R=R)
      Rbig[(2*(betadim+gammadim)+3):(length(parest)),(2*(betadim+gammadim)+3):(length(parest))]=Rsecderiv
      df=sum(diag((result)%*%solve(result+Rbig)))
      log.likeli=Loglikelihood(parest,myrules,Delta,X,Z,n,ni,r,blC,betadim,gammadim,lambda,R,Cauchy.pen)
      AIC=df-log.likeli
      Informationmatrix=result+Rbig
      var=diag(solve(Informationmatrix,tol=1e-300))[1:(2*(betadim+gammadim+1))]}
    else{
      df=length(parest)
      log.likeli=Loglikelihood(parest,myrules,Delta,X,Z,n,ni,r,blC,betadim,gammadim,lambda,R,Cauchy.pen)
      AIC=df-log.likeli
      Informationmatrix=result
      var=diag(solve(Informationmatrix,tol=1e-300))[1:(2*(betadim+gammadim+1))]
    }
  }
  else{
    if(lambda!=0){
      Rbig=matrix(0, nrow = length(parest), ncol = length(parest))
      Rsecderiv=numDeriv::hessian(penaltyterm,x=psi,lambda=lambda,R=R)
      Rbig[(2*(betadim+gammadim)+3):(length(parest)),(2*(betadim+gammadim)+3):(length(parest))]=Rsecderiv
      df=sum(diag((result)%*%solve(result+Rbig)))
      log.likeli=Loglikelihood(parest,myrules,Delta,X,Z,n,ni,r,blC,betadim,gammadim,lambda,R)
      AIC=df-log.likeli
      Informationmatrix=result+Rbig
      var=diag(solve(Informationmatrix,tol=1e-300))[1:(2*(betadim+gammadim+1))]
    }
    else{
      df=length(parest)
      log.likeli=Loglikelihood(parest,myrules,Delta,X,Z,n,ni,r,blC,betadim,gammadim,lambda,R,Cauchy.pen)
      AIC=df-log.likeli
      Informationmatrix=result
      var=diag(solve(Informationmatrix,tol=1e-300))[1:(2*(betadim+gammadim+1))]
    }
  }
  finalresult=list()
  length(finalresult)=3
  finalresult[[1]]=log.likeli
  finalresult[[2]]=var
  finalresult[[3]]=AIC
  return(finalresult)
}
.EMitercure=function(lastpar,X,Z,Delta,n,ni,blC,r,myrules,lambda,R,Cauchy.pen,clustering){
  betadim=dim(X[[1]])[2]
  gammadim=dim(Z)[2]
  zeta=lastpar[1]
  betagamma=lastpar[2:length(lastpar)]
  beta=betagamma[1:betadim]
  gamma=betagamma[(betadim+1):(betadim+gammadim)]
  theta=betagamma[betadim+gammadim+1]
  psi=betagamma[(2+betadim+gammadim):length(betagamma)]
  if(clustering){
  normalC=NormalConstantcure(myrules,beta,gamma,zeta,theta,Delta,betadim,gammadim,n,ni,X,Z,psi,blC,r)
  g0bvalue=g0bmatcure(myrules,beta,gamma,zeta,theta,Delta,betadim,gammadim,n,ni,X,Z,psi,blC,r)
  }
  else{
    myrule[,1]=rep(0,dim(myrules)[1])
    myrule[,2]=rep(1,dim(myrules)[1])
    normalC=rep(dim(myrules)[1],n)
    g0bvalue=matrix(1, nrow = dim(myrules)[1], ncol = n)
  }
  uijqmat=uijqcure(myrules,beta,gamma,zeta,theta,Delta,betadim,gammadim,n,ni,X,Z,psi,blC,r)
  
  zeta=optim(par = zeta,fn=Q1cure,rules=myrules,beta=beta,gamma=gamma,eta=zeta,
             theta=theta,Delta=Delta,betadim=betadim,gammadim=gammadim,n=n,ni=ni,X=X,Z=Z,psi=psi,blC=blC,r=r,g0bvalue=g0bvalue,
             normalconstant=normalC,uijqmat=uijqmat,Cauchyindex=Cauchy.pen,method = "BFGS")$par
  betagamma=optim(par = betagamma,fn=Q2cure,rules=myrules,beta=beta,gamma=gamma,eta=zeta,
                  theta=theta,Delta=Delta,betadim=betadim,gammadim=gammadim,n=n,ni=ni,X=X,Z=as.matrix(Z),psi=psi,blC=blC,r=r,g0bvalue=g0bvalue,
                  normalconstant=normalC,uijqmat=uijqmat,tunepar=lambda,R=R,Cauchyindex=Cauchy.pen,method = "BFGS")$par
  result=c(zeta,betagamma)
  return(result)
}
.EMestcure=function(X,Z,Delta,n,ni,C,r,quadnum,lambda,criterion,blC,myrules,Cauchy.pen,clustering){
  
  coefnum=dim(blC[[1]])[1]
  D=matrix(0, nrow = coefnum-2, ncol = coefnum)
  for(i in 1:(coefnum-2)){
    D[i,i:(i+2)]=c(1,-2,1)
  }
  R=(t(D)%*%D)
  betadim=dim(X[[1]])[2]
  gammadim=dim(Z)[2]
  lastpar=rep(0.1,2+(betadim+gammadim)+coefnum)
  difference=1
  while (difference>criterion) {
    outputpar=.EMitercure(lastpar,X,Z,Delta,n,ni,blC,r,myrules,lambda,R,Cauchy.pen,clustering)
    difference=sum(abs((-outputpar+lastpar)/(lastpar)))
    lastpar=outputpar
  }
  
  return(outputpar)
}
.AICcomputepenparacure=function(parest,X,Z,Delta,blC,r,n,ni,myrules,lambda,Cauchy.pen,clustering){
  coefnum=dim(blC[[1]])[1]
  D=matrix(0, nrow = coefnum-2, ncol = coefnum)
  for(i in 1:(coefnum-2)){
    D[i,i:(i+2)]=c(1,-2,1)
  }
  R=(t(D)%*%D)
  pardim=length(parest)
  result=matrix(0, nrow = pardim, ncol = pardim)
  betadim=dim(X[[1]])[2]
  gammadim=dim(Z)[2]
  zeta=parest[1]
  betagamma=parest[2:length(parest)]
  beta=betagamma[1:betadim]
  gamma=betagamma[(betadim+1):(betadim+gammadim)]
  theta=betagamma[betadim+gammadim+1]
  psi=betagamma[(2+betadim+gammadim):length(betagamma)]
  if(clustering){
  g0bvalue=g0bmatcure(myrules,beta,gamma,zeta,theta,Delta,betadim,gammadim,n,ni,X,Z,psi,blC,r)
  normalC=NormalConstantcure(myrules,beta,gamma,zeta,theta,Delta,betadim,gammadim,n,ni,X,Z,psi,blC,r)
  }
  else{
    myrule[,1]=rep(0,dim(myrules)[1])
    myrule[,2]=rep(1,dim(myrules)[1])
    normalC=rep(dim(myrules)[1],n)
    g0bvalue=matrix(1, nrow = dim(myrules)[1], ncol = n)
  }
  quadnum=dim(myrules)[1]
  for (i in 1:n) {
    resulti=matrix(0, nrow = pardim, ncol = pardim)
    for (k in 1:quadnum) {
      
      firstdrivk=firstderiv_icure(parest,myrules[k,1],Delta[[i]],X[[i]],Z[i,],ni[i],r,blC[[i]],betadim,gammadim)
      if(sum(is.nan( firstdrivk)>0)){
        firstdrivk=numDeriv::grad(loglikelihoodtestcure,x=parest,b=myrules[k,1],Delta=Delta[[i]],X=X[[i]],Z=Z[i,],
                                  ni=ni[i],r=r,blC=blC[[i]],betadim=betadim,gammadim=gammadim)
      }
      for (q in 1:quadnum) {
        
        firstdrivq=firstderiv_icure(parest,myrules[q,1],Delta[[i]],X[[i]],Z[i,],ni[i],r,blC[[i]],betadim,gammadim)
        if(sum(is.nan( firstdrivq)>0)){
          firstdrivq=numDeriv::grad(loglikelihoodtestcure,x=parest,b=myrules[q,1],Delta=Delta[[i]],X=X[[i]],Z=Z[i,],
                                    ni=ni[i],r=r,blC=blC[[i]],betadim=betadim,gammadim=gammadim)
        }
        tempresult=as.matrix(firstdrivk)%*%t(as.matrix(firstdrivq))*myrules[k,2]*g0bvalue[k,i]*myrules[q,2]*g0bvalue[q,i]/(normalC[i])^2
        
        resulti=resulti+tempresult
        if(sum(is.nan(tempresult))>0){
          print(c(i,k,q))
        }
      }
    }
    result=result+resulti
  }
  if(Cauchy.pen){
    hessian=numDeriv::hessian(logcauchy,x=parest)
    result=result-hessian
    if(lambda!=0){
      Rbig=matrix(0, nrow = length(parest), ncol = length(parest))
      Rsecderiv=numDeriv::hessian(penaltyterm,x=psi,lambda=lambda,R=R)
      Rbig[(2*(betadim+gammadim)+3):(length(parest)),(2*(betadim+gammadim)+3):(length(parest))]=Rsecderiv
      df=sum(diag((result)%*%solve(result+Rbig)))
      log.likeli=Loglikelihoodcure(parest,myrules,Delta,X,Z,n,ni,r,blC,betadim,gammadim,lambda,R,Cauchy.pen)
      AIC=df-log.likeli
      Informationmatrix=result+Rbig
      var=diag(solve(Informationmatrix,tol=1e-300))[1:(betadim+gammadim+2)]}
    else{
      df=length(parest)
      log.likeli=Loglikelihoodcure(parest,myrules,Delta,X,Z,n,ni,r,blC,betadim,gammadim,lambda,R,Cauchy.pen)
      AIC=df-log.likeli
      Informationmatrix=result
      var=diag(solve(Informationmatrix,tol=1e-300))[1:(betadim+gammadim+2)]
    }
  }
  else{
    if(lambda!=0){
      Rbig=matrix(0, nrow = length(parest), ncol = length(parest))
      Rsecderiv=numDeriv::hessian(penaltyterm,x=psi,lambda=lambda,R=R)
      Rbig[(2*(betadim+gammadim)+3):(length(parest)),(2*(betadim+gammadim)+3):(length(parest))]=Rsecderiv
      df=sum(diag((result)%*%solve(result+Rbig)))
      log.likeli=Loglikelihoodcure(parest,myrules,Delta,X,Z,n,ni,r,blC,betadim,gammadim,lambda,R)
      AIC=df-log.likeli
      Informationmatrix=result+Rsecderiv
      var=diag(solve(Informationmatrix,tol=1e-300))[1:(betadim+gammadim+2)]
    }
    else{
      df=length(parest)
      log.likeli=Loglikelihoodcure(parest,myrules,Delta,X,Z,n,ni,r,blC,betadim,gammadim,lambda,R,Cauchy.pen)
      AIC=df-log.likeli
      Informationmatrix=result
      var=diag(solve(Informationmatrix,tol=1e-300))[1:(betadim+gammadim+2)]
    }
  }
  finalresult=list()
  length(finalresult)=3
  finalresult[[1]]=log.likeli
  finalresult[[2]]=var
  finalresult[[3]]=AIC
  return(finalresult)
}
.Data.trans=function(Rawdata,n_subject.cov,n_tooth.cov){
  n=length(unique(Rawdata[,1]))
  ni=sort(unique(Rawdata[,1]))
  mi=rep(0,n)
  for(i in 1:n){
    mi[i]=sum(Rawdata[,1]==ni[i])
  }
  CSTime=list()
  length(CSTime)=n
  X=list()
  length(X)=n
  Delta=list()
  length(Delta)=n
  for (i in 1:n) {
    CSTime[[i]]=Rawdata[,2][Rawdata[,1]==ni[i]]
    X[[i]]=as.matrix(Rawdata[Rawdata[,1]==ni[i],(3+n_subject.cov):(dim(Rawdata)[2]-1)])
    Delta[[i]]=Rawdata[,dim(Rawdata)[2]][Rawdata[,1]==ni[i]]
  }
  Z=matrix(0, nrow = n, ncol = n_subject.cov)
  
  for (i in 1:n) {
    
    Z[i,]=as.numeric(Rawdata[sum(mi[0:i]),3:(2+n_subject.cov)])
  }
  Z=as.matrix(Z)
  
  finalresult=list()
  length(finalresult)=6
  finalresult[[1]]=n
  finalresult[[2]]=mi
  finalresult[[3]]=CSTime
  finalresult[[4]]=Delta
  finalresult[[5]]=X
  finalresult[[6]]=Z
  return(finalresult)
}
.Datascale=function(Rawdata,n_subject.cov,n_tooth.cov){
  covdata=Rawdata[,3:(2+n_subject.cov+n_tooth.cov)]
  covdim=dim(covdata)[2]
  
  for(i in 1:covdim){
    if(is.numeric(covdata[,i])){
      covdata[,i]=scale(covdata[,i])[,1]
    }
    else{
      covdata[,i]=as.numeric(covdata[,i])-1
    }
  }
  finalresult=cbind(Rawdata[,1:2],covdata,Rawdata[,dim(Rawdata)[2]])
  return(finalresult)
}
.H.est=function(t,coefs,Rawdata,knots.num=2,degree=2){
  minCSTime=min(Rawdata[,2])
  maxCSTime=max(Rawdata[,2])
  knots <- seq(0,1  , length.out = (knots.num + 2))
  knots=knots[3:length(knots)-1]
  
  
  blC=t(splines2::ibs((t-(minCSTime-0.1))/(maxCSTime-minCSTime+0.2),knots = knots,degree=degree
                      ,Boundary.knots = c(0,1),intercept = TRUE))
  return(sum(blC*coefs))
  
}
.Surv.est=function(t,covariate,reg.est,coefs,Rawdata,n_subject.raw,n_within.raw,r,n_quad=30,knots.num=2,degree=2,cure.reg=TRUE){
  myrules=.ghrules(n_quad,normalized=FALSE)
  myrules=as.matrix(myrules[[n_quad]])
  Hhat=.H.est(t,coefs,Rawdata,knots.num,degree)
  Z=covariate[1:n_subject.raw]
  X=covariate[(n_subject.raw+1):length(covariate)]
  if(cure.reg){
    pi=1/(1+exp(-reg.est[1]-sum(reg.est[2:(1+n_within.raw)]*X)-sum(reg.est[(2+n_within.raw):(1+n_within.raw+n_subject.raw)])))
    beta=reg.est[(2+n_within.raw+n_subject.raw):(1+2*n_within.raw+n_subject.raw)]
    gamma=reg.est[(2+2*n_within.raw+n_subject.raw):(1+2*n_within.raw+2*n_subject.raw)]
    theta=exp(reg.est[length(reg.est)])
    if(r==0){
      S.sucp=sum(exp(-Hhat*exp(sum(X*beta)+sum(Z*gamma)+theta*myrules[,1]))*dnorm(myrules[,1],0,1)*exp(myrules[,1]^2)*myrules[,2])
    }
    else{
      S.sucp=sum((1+r*Hhat*exp(sum(X*beta)+sum(Z*gamma)+theta*myrules[,1]))^(-1/r)*dnorm(myrules[,1],0,1)*exp(myrules[,1]^2)*myrules[,2])
    }
    result=pi+(1-pi)*S.sucp
  }
  else{
    pi=reg.est[1]
    beta=reg.est[(2):(1+n_within.raw)]
    gamma=reg.est[(2+n_within.raw):(1+n_within.raw+n_subject.raw)]
    theta=exp(reg.est[length(reg.est)])
    if(r==0){
      S.sucp=sum(exp(-Hhat*exp(sum(X*beta)+sum(Z*gamma)+theta*myrules[,1]))*dnorm(myrules[,1],0,1)*exp(myrules[,1]^2)*myrules[,2])
    }
    else{
      S.sucp=sum((1+r*Hhat*exp(sum(X*beta)+sum(Z*gamma)+theta*myrules[,1]))^(-1/r)*dnorm(myrules[,1],0,1)*exp(myrules[,1]^2)*myrules[,2])
    }
    result=pi+(1-pi)*S.sucp
  }
  return(result)
}

.firstdriv_xi=function(omega,parest,rules,Delta,X,Z,n,ni,r,blC,betadim,gammadim){
  result=numDeriv::grad(testquadrature1,x=parest,rules=rules,Delta=Delta,X=X,Z=Z,n=n,
                        ni=ni,r=r,blC=blC,betadim=betadim,gammadim=gammadim,weight=omega)
  return(result)
}

.jacob_weight=function(omega,parest,rules,Delta,X,Z,n,ni,r,blC,betadim,gammadim){
  result=numDeriv::jacobian(.firstdriv_xi,x=omega,parest=parest,rules=rules,Delta=Delta,X=X,Z=Z,n=n,
                            ni=ni,r=r,blC=blC,betadim=betadim,gammadim=gammadim)
  return(result)
}






CSDfit=function(Rawdata,n_subject.raw,n_within.raw,r,n_quad=30,lambda=0,Cauchy.pen=TRUE,tolerance=1e-2,
                knots.num=2,degree=2,scale.numr=TRUE,cure.reg=TRUE,clustering=TRUE){
  myrules=.ghrules(n_quad,normalized=FALSE)
  myrules=as.matrix(myrules[[n_quad]])
  if(scale.numr==TRUE){
    Rawdata=.Datascale(Rawdata,n_subject.raw,n_within.raw)
    n_tooth.cov=n_within.raw
    n_subject.cov=n_subject.raw
    totaldata=.Data.trans(Rawdata,n_subject.cov,n_tooth.cov)
  }
  else{
    totaldata=.Data.trans(Rawdata,n_subject.raw,n_within.raw)
    n_tooth.cov=n_within.raw
    n_subject.cov=n_subject.raw
  }
  n=totaldata[[1]]
  ni=totaldata[[2]]
  C=totaldata[[3]]
  Delta=totaldata[[4]]
  X=totaldata[[5]]
  Z=totaldata[[6]]
  minCSTime=min(Rawdata[,2])
  maxCSTime=max(Rawdata[,2])
  blC <- list()
  length(blC) <- n
  knots <- seq(0,1  , length.out = (knots.num + 2))
  knots=knots[3:length(knots)-1]
  
  for (i in 1:n) {
    blC[[i]]=t(splines2::ibs((C[[i]]-(minCSTime-0.1))/(maxCSTime-minCSTime+0.2),knots = knots,degree=degree
                             ,Boundary.knots = c(0,1),intercept = TRUE))
  }
  
  if(cure.reg){
    par.est=.EMest(X,Z,Delta,n,ni,C,r,n_quad,lambda,tolerance,blC,myrules,Cauchy.pen,clustering)
    loglikeli.var.AIC=.AICcomputepenpara(par.est,X,Z,Delta,blC,r,n,ni,myrules,lambda,Cauchy.pen,clustering)
    reg.est=par.est[1:(2*n_subject.cov+2*n_tooth.cov+2)]
    reg.est[length(reg.est)]=exp(reg.est[length(reg.est)])
    reg.se=sqrt(loglikeli.var.AIC[[2]][1:(2*n_subject.cov+2*n_tooth.cov+2)])
    reg.se[length(reg.se)]=reg.est[length(reg.est)]*reg.se[length(reg.se)]
    z.stat=abs(reg.est/reg.se)
    p.est=2*(1-pnorm(z.stat))
    resultmat=matrix(0, nrow = 2*n_subject.cov+2*n_tooth.cov+2, ncol = 6)
    resultmat[,1]=reg.est
    resultmat[,2]=reg.se
    resultmat[,3]=z.stat
    resultmat[,4]=p.est
    resultmat[,5]=resultmat[,1]-1.96*resultmat[,2]
    resultmat[,6]=resultmat[,1]+1.96*resultmat[,2]
    resultmat[length(reg.est),5]=exp(par.est[length(reg.est)]-1.96*sqrt(loglikeli.var.AIC[[2]][length(reg.est)]))
    resultmat[length(reg.est),6]=exp(par.est[length(reg.est)]+1.96*sqrt(loglikeli.var.AIC[[2]][length(reg.est)]))
    resultmat=round(resultmat,2)
    resultmat=data.frame(resultmat)
    colnames(resultmat)=c("par.est","SE","Z","p-value","CI_lower","CI_upper")
    rownames(resultmat)[1]="Intercept"
    rownames(resultmat)[dim(resultmat)[1]]="Random_effect"
    for(k in 1:n_tooth.cov){
      rownames(resultmat)[1+k]=paste("Cure",colnames(Rawdata)[2+n_subject.cov+k],sep = "_")
    }
    for (k in 1:n_subject.cov) {
      rownames(resultmat)[1+k+n_tooth.cov]=paste("Cure",colnames(Rawdata)[2+k],sep = "_")
    }
    for(k in 1:n_tooth.cov){
      rownames(resultmat)[1+n_subject.cov+n_tooth.cov+k]=paste("Surv",colnames(Rawdata)[2+n_subject.cov+k],sep = "_")
    }
    for (k in 1:n_subject.cov) {
      rownames(resultmat)[1+2*n_tooth.cov+n_subject.cov+k]=paste("Surv",colnames(Rawdata)[2+k],sep = "_")
    }
  }
  else{
    par.est=.EMestcure(X,Z,Delta,n,ni,C,r,n_quad,lambda,tolerance,blC,myrules,Cauchy.pen,clustering)
    loglikeli.var.AIC=.AICcomputepenparacure(par.est,X,Z,Delta,blC,r,n,ni,myrules,lambda,Cauchy.pen,clustering)
    reg.est=par.est[1:(n_subject.cov+n_tooth.cov+2)]
    reg.est[length(reg.est)]=exp(reg.est[length(reg.est)])
    reg.se=sqrt(loglikeli.var.AIC[[2]])
    reg.se[1]=(1+exp(-reg.est[1]))^(-2)*exp(-reg.est[1])*reg.se[1]
    reg.est[1]=(1+exp(-reg.est[1]))^(-1)
    reg.se[length(reg.se)]=reg.est[length(reg.est)]*reg.se[length(reg.se)]
    z.stat=abs(reg.est/reg.se)
    p.est=2*(1-pnorm(z.stat))
    resultmat=matrix(0, nrow = n_subject.cov+n_tooth.cov+2, ncol = 6)
    resultmat[,1]=reg.est
    resultmat[,2]=reg.se
    resultmat[,3]=z.stat
    resultmat[,4]=p.est
    resultmat[,5]=resultmat[,1]-1.96*resultmat[,2]
    resultmat[,6]=resultmat[,1]+1.96*resultmat[,2]
    
    resultmat[length(reg.est),5]=exp(par.est[length(reg.est)]-1.96*sqrt(loglikeli.var.AIC[[2]][length(reg.est)]))
    resultmat[length(reg.est),6]=exp(par.est[length(reg.est)]+1.96*sqrt(loglikeli.var.AIC[[2]][length(reg.est)]))
    resultmat=round(resultmat,2)
    resultmat=data.frame(resultmat)
    colnames(resultmat)=c("par.est","SE","Z","p-value","CI_lower","CI_upper")
    rownames(resultmat)[1]="cure_rate"
    rownames(resultmat)[dim(resultmat)[1]]="Random_effect"
    for(k in 1:n_tooth.cov){
      rownames(resultmat)[1+k]=paste("Surv",colnames(Rawdata)[2+n_subject.cov+k],sep = "_")
    }
    for (k in 1:n_subject.cov) {
      rownames(resultmat)[1+k+n_tooth.cov]=paste("Surv",colnames(Rawdata)[2+k],sep = "_")
    }
    
  }
  

  output=list(parameter.est=resultmat,log_likelihood=loglikeli.var.AIC[[1]],AICvalue=loglikeli.var.AIC[[3]],
              coefs=round(exp(par.est[(n_subject.cov+n_tooth.cov+3):length(par.est)]),2))
  return(output)
}



boot.CSD=function(Rawdata,n_subject.raw,n_within.raw,r,boot.rep,seed.begin,
                  n_quad=30,lambda=0,Cauchy.pen=TRUE,tolerance=1e-2,
                  knots.num=2,degree=2,scale.numr=TRUE,cure.reg=TRUE,clustering=TRUE){
  myrules=.ghrules(n_quad,normalized=FALSE)
  myrules=as.matrix(myrules[[n_quad]])
  n_tooth.cov=n_within.raw
  n_subject.cov=n_subject.raw
  if(scale.numr==TRUE){
    Rawdata=.Datascale(Rawdata,n_subject.raw,n_within.raw)
    totaldata=.Data.trans(Rawdata,n_subject.cov,n_tooth.cov)
  }

  
  n=length(unique(Rawdata[,1]))
  bootnum=n
  bootresult=c()
  for(b in 1:boot.rep){
    
    set.seed(b+seed.begin)
    ind=sample(1:n,bootnum,replace = T)
    ni=sort(unique(Rawdata[,1]))
    mi=rep(0,n)
    for(i in 1:n){
      mi[i]=sum(Rawdata[,1]==ni[i])
    }
    Z=matrix(0, nrow = n, ncol = (n_subject.cov+1))
    
    for (i in 1:n) {
      Z[i,]=as.numeric(c(Rawdata[sum(mi[0:i]),1],Rawdata[sum(mi[0:i]),(3:(2+n_subject.cov))]))
    }
    ni=ni[ind]
    mi=rep(0,bootnum)
    for(i in 1:bootnum){
      mi[i]=sum(Rawdata[,1]==ni[i])
    }
    C=list()
    length(C)=bootnum
    X=list()
    length(X)=bootnum
    Delta=list()
    length(Delta)=bootnum
    for (i in 1:bootnum) {
      C[[i]]=Rawdata[Rawdata[,1]==ni[i],2]
      X[[i]]=as.matrix(as.numeric(Rawdata[Rawdata[,1]==ni[i],((3+n_subject.cov):(2+n_subject.cov+n_tooth.cov))]))
      Delta[[i]]=Rawdata[Rawdata[,1]==ni[i],dim(Rawdata)[2]]
    }
    newZ=matrix(0, nrow = bootnum, ncol = n_subject.cov)
    for(i in 1:bootnum){
      for(j in 1:n){
        if(Z[j,1]==ni[i]){
          newZ[i,]=Z[j,-1]
        }
      }
    }
    minCSTime=min(Rawdata[,2])
    maxCSTime=max(Rawdata[,2])
    blC <- list()
    length(blC) <- n
    knots <- seq(0,1  , length.out = (knots.num + 2))
    knots=knots[3:length(knots)-1]
    
    for (i in 1:bootnum) {
      blC[[i]]=t(splines2::ibs((C[[i]]-(minCSTime-0.1))/(maxCSTime-minCSTime+0.2),knots = knots,degree=degree
                               ,Boundary.knots = c(0,1),intercept = TRUE))
    }
    if(cure.reg){
    par.est=.EMest(X,newZ,Delta,bootnum,mi,C,r,n_quad,lambda,tolerance,blC,myrules,Cauchy.pen,clustering)
    
    }
    else{
      par.est=.EMestcure(X,newZ,Delta,bootnum,mi,C,r,n_quad,lambda,tolerance,blC,myrules,Cauchy.pen,clustering)
      par.est[1]=(1+exp(-par.est[1]))^(-1)
    }
    bootresult=rbind(bootresult,par.est)
  }
  bootresult=data.frame(bootresult)
  if(cure.reg){
    colnames(bootresult)[1]="Intercept"
    colnames(bootresult)[2+2*n_subject.cov+2*n_tooth.cov]="Random_effect"
    for(k in 1:n_tooth.cov){
      colnames(bootresult)[1+k]=paste("Cure",colnames(Rawdata)[2+n_subject.cov+k],sep = "_")
    }
    for (k in 1:n_subject.cov) {
      colnames(bootresult)[1+k+n_tooth.cov]=paste("Cure",colnames(Rawdata)[2+k],sep = "_")
    }
    for(k in 1:n_tooth.cov){
      colnames(bootresult)[1+n_subject.cov+n_tooth.cov+k]=paste("Surv",colnames(Rawdata)[2+n_subject.cov+k],sep = "_")
    }
    for (k in 1:n_subject.cov) {
      colnames(bootresult)[1+2*n_tooth.cov+n_subject.cov+k]=paste("Surv",colnames(Rawdata)[2+k],sep = "_")
    }
    for (k in 1:dim(blC[[1]])[1]) {
      colnames(bootresult)[2*n_tooth.cov+2*n_subject.cov+2+k]=paste("coef",k,sep="_")
    }
  }
  else{
    colnames(bootresult)[1]="Cure_rate"
    colnames(bootresult)[2+n_subject.cov+n_tooth.cov]="Random_effect"
    for(k in 1:n_tooth.cov){
      colnames(bootresult)[1+k]=paste("Surv",colnames(Rawdata)[2+n_subject.cov+k],sep = "_")
    }
    for (k in 1:n_subject.cov) {
      colnames(bootresult)[1+k+n_tooth.cov]=paste("Surv",colnames(Rawdata)[2+k],sep = "_")
    }
    for (k in 1:dim(blC[[1]])[1]) {
      colnames(bootresult)[n_tooth.cov+n_subject.cov+2+k]=paste("coef",k,sep="_")
    }
  }
  
 
  return(bootresult)
}




Surv.CI=function(t,x,boot.par,Rawdata,n_subject,n_within,r,CI.lev=0.95,n_quad=30,knots.num=2,degree=2,
                 lambda=0,Cauchy.pen=TRUE,tolerance=1e-2,scale.numr=TRUE,cure.reg=TRUE,clustering=TRUE){
  myrules=.ghrules(n_quad,normalized=FALSE)
  myrules=as.matrix(myrules[[n_quad]])
  if(scale.numr==TRUE){
    Rawdata=.Datascale(Rawdata,n_subject,n_within)
    n_tooth.cov=n_within
    n_subject.cov=n_subject
    totaldata=.Data.trans(Rawdata,n_subject.cov,n_tooth.cov)
  }
  else{
    totaldata=.Data.trans(Rawdata,n_subject,n_within)
    n_tooth.cov=n_within
    n_subject.cov=n_subject
  }
  n=totaldata[[1]]
  ni=totaldata[[2]]
  C=totaldata[[3]]
  Delta=totaldata[[4]]
  X=totaldata[[5]]
  Z=totaldata[[6]]
  minCSTime=min(Rawdata[,2])
  maxCSTime=max(Rawdata[,2])
  blC <- list()
  length(blC) <- n
  knots <- seq(0,1  , length.out = (knots.num + 2))
  knots=knots[3:length(knots)-1]
  t.num=length(t)
  result=c()
  for (i in 1:n) {
    blC[[i]]=t(splines2::ibs((C[[i]]-(minCSTime-0.1))/(maxCSTime-minCSTime+0.2),knots = knots,degree=degree
                             ,Boundary.knots = c(0,1),intercept = TRUE))
  }
  for(l in 1:t.num){
    if(cure.reg){
      par.est=.EMest(X,Z,Delta,n,ni,C,r,n_quad,lambda,tolerance,blC,myrules,Cauchy.pen,clustering)
      St=.Surv.est(t[l],x[l,],par.est[1:(2+2*n_subject.cov+2*n_tooth.cov)],
                   exp(par.est[(3+2*n_subject.cov+2*n_tooth.cov):length(par.est)]),Rawdata,n_subject,n_within,r,
                   n_quad,knots.num,degree,cure.reg = TRUE)
      boot.surv=c()
      for(b in 1:dim(boot.par)[1]){
        survb=.Surv.est(t[l],x[l,],as.numeric(boot.par[b,(1:(2+2*n_subject.cov+2*n_tooth.cov))]),
                        as.numeric(exp(boot.par[b,((3+2*n_subject.cov+2*n_tooth.cov):length(par.est))])),Rawdata,n_subject,n_within,r,
                        n_quad,knots.num,degree,cure.reg = TRUE)
        boot.surv=append(boot.surv,survb)
      }
    }
    else{
      par.est=.EMestcure(X,Z,Delta,n,ni,C,r,n_quad,lambda,tolerance,blC,myrules,Cauchy.pen,clustering)
      St=.Surv.est(t[l],x[l,],c((1+exp(-par.est[1]))^(-1),par.est[2:(2+n_subject.cov+n_tooth.cov)]),
                   exp(par.est[(3+n_subject.cov+n_tooth.cov):length(par.est)]),Rawdata,n_subject,n_within,r,
                   n_quad,knots.num,degree,cure.reg = FALSE)
      boot.surv=c()
      for(b in 1:dim(boot.par)[1]){
        survb=.Surv.est(t[l],x[l,],as.numeric(boot.par[b,(1:(2+n_subject.cov+n_tooth.cov))]),
                        as.numeric(exp(boot.par[b,((3+n_subject.cov+n_tooth.cov):length(par.est))])),Rawdata,n_subject,n_within,r,
                        n_quad,knots.num,degree,cure.reg = FALSE)
        boot.surv=append(boot.surv,survb)
      }
    }
    
    
    bound=quantile(boot.surv,probs = c((1-CI.lev)/2,(1+CI.lev)/2))
    result=rbind(result,c(St,bound))
  }
  result=data.frame(result)
  colnames(result)[1]="Pred_surv"
  colnames(result)[2]="CI_lower"
  colnames(result)[3]="CI_upper"
  return(result)
}


