#include <RcppArmadillo.h>
#include <iostream>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 

// [[Rcpp::export]]

arma::vec firstderiv_i(const arma::vec&parameters,const double&b,const arma::vec&Delta,
                       const arma::mat&X,const arma::vec&Z,const int&ni,const double&r,const arma::mat&blC,
                       const int&betadim,const int&gammadim){
  arma::vec result(parameters.n_elem);
  arma::vec covariate(betadim+gammadim);
  double eta0=parameters(0);
  arma::vec eta1=parameters.subvec(1,betadim);
  arma::vec eta2=parameters.subvec(1+betadim,betadim+gammadim);
  arma::vec beta=parameters.subvec(betadim+gammadim+1,betadim+gammadim+betadim);
  arma::vec gamma=parameters.subvec(betadim+gammadim+1+betadim,betadim+gammadim+betadim+gammadim);
  double theta=parameters(betadim+gammadim+betadim+gammadim+1);
  arma::vec psi=parameters.subvec(betadim+gammadim+betadim+gammadim+2,parameters.n_elem-1);
  arma::vec onesvector;
  arma::vec pi,si,expterm,Ht;
  onesvector.ones(ni);
  arma::vec exppilinear=exp(-onesvector*eta0-X*eta1-sum(Z%eta2)*onesvector);
  pi=pow(1+exppilinear,-1);
  expterm=exp(X*beta+sum(Z%gamma)*onesvector+std::exp(theta)*b*onesvector);
  Ht=trans(blC)*exp(psi);
  arma::vec pisi;
  
  if(r>0){
    si=pow(1+r*Ht%expterm,-1/r);
  }
  else{
    si=exp(-Ht%expterm);
  }
  pisi=pi+(onesvector-pi)%si;

  arma::vec piderivpart1,piderivpart2;
  piderivpart1=pow(1+exppilinear,-2)%exppilinear;
  arma::vec commonvec=Delta%pow(1-pi,-1)%(-piderivpart1)+(1-Delta)%pow(pi+(1-pi)%si,-1)%(1-si)%piderivpart1;
  result(0)=sum(commonvec);
  for(int j=0;j<ni;j++){
    result.subvec(1,betadim)=result.subvec(1,betadim)+trans(X.row(j))*commonvec(j);
    result.subvec(1+betadim,betadim+gammadim)=result.subvec(1+betadim,betadim+gammadim)+Z*commonvec(j);
  }
  arma::vec derivsi;
  if(r>0){
    derivsi=-pow(1+r*Ht%expterm,-1/r-1)%Ht%expterm;
  }
  else{
    derivsi=-exp(-Ht%expterm)%Ht%expterm;
  }
  arma::vec commonvec2=Delta%pow(1-si,-1)%(-derivsi)+(1-Delta)%pow(pi+(1-pi)%si,-1)%(1-pi)%derivsi;
  for(int j=0;j<ni;j++){
    result.subvec(betadim+gammadim+1,betadim+gammadim+betadim)=result.subvec(betadim+gammadim+1,betadim+gammadim+betadim)+
      trans(X.row(j))*commonvec2(j);
    result.subvec(betadim+gammadim+1+betadim,betadim+gammadim+betadim+gammadim)=result.subvec(betadim+gammadim+1+betadim,betadim+gammadim+betadim+gammadim)
      +Z*commonvec2(j);
    result(betadim+gammadim+betadim+gammadim+1)=sum(commonvec2*b*exp(theta));
  }
  arma::vec commonvec3,sipsideriv;
  if(r>0){
    sipsideriv=-pow(1+r*Ht%expterm,-1/r-1)%expterm;
  }
  else{
    sipsideriv=-exp(-Ht%expterm)%expterm;
  }
  commonvec3=Delta%pow(1-si,-1)%(-sipsideriv)+(1-Delta)%pow(pi+(1-pi)%si,-1)%(1-pi)%sipsideriv;
  for(int j=0;j<psi.n_elem;j++){
    result(betadim+gammadim+betadim+gammadim+2+j)=sum(exp(psi(j))*commonvec3%trans(blC.row(j)));
  }
  
  return result;
  
}


double g0b(const double&b,const arma::vec&beta,const arma::vec&gamma,const arma::vec&eta,const double&theta,
           const int&betadim,const int&gammadim,const arma::vec&Delta,
           const arma::mat&X,const arma::vec&Z,
           const arma::vec&psi,const arma::mat&blC,const double&r,const int&ni){
  double eta0=eta(0);
  arma::vec eta1=eta.subvec(1,betadim);
  arma::vec eta2=eta.subvec(betadim+1,betadim+gammadim);
  arma::vec pilinear=eta0+X*eta1+sum(trans(eta2)*Z);
  
  arma::vec pi=pow((1+exp(-pilinear)),-1);
  arma::vec S;
  arma::vec Slinear=X*beta+sum(trans(gamma)*Z)+std::exp(theta)*b;
  
  arma::vec H_t=trans(blC)*exp(psi);
  if(r==0){
    S=exp(-H_t%exp(Slinear));
  }
  else{
    S=pow(1+r*H_t%exp(Slinear),-1/r);
  }
  double result=1;
  for(int i=0;i<ni;i++){
    result=result*std::pow((1-pi(i))*(1-S(i)),Delta(i))*std::pow(pi(i)+(1-pi(i))*S(i),1-Delta(i));
  }
  result=result*R::dnorm4(b,0,1,false)*std::exp(b*b);
  return result;
}
// [[Rcpp::export]]
arma::mat g0bmat(const arma::mat&rules,const arma::vec&beta,const arma::vec&gamma,const arma::vec&eta,const double&theta,
                 const arma::field<arma::vec>&Delta,const int&betadim,const int&gammadim,
                 const int&n,const arma::vec&ni,
                 const arma::field<arma::mat>&X,
                 const arma::mat&Z,const arma::vec&psi,const arma::field<arma::mat>&blC,
                 const double&r){
  int order=rules.n_rows;
  arma::mat result(order,n);
  for(int q=0;q<order;q++){
    for(int i=0;i<n;i++){
      result(q,i)=g0b(rules(q,0),beta,gamma,eta,theta,betadim,gammadim,Delta(i),X(i),trans(Z.row(i)),psi,blC(i),r,ni(i));
    }
  }
  return result;
}




// [[Rcpp::export]]

arma::vec NormalConstant(const arma::mat&rules,const arma::vec&beta,const arma::vec&gamma,const arma::vec&eta,const double&theta,
                         const arma::field<arma::vec>&Delta,const int&betadim,const int&gammadim,
                         const int&n,const arma::vec&ni,
                         const arma::field<arma::mat>&X,
                         const arma::mat&Z,const arma::vec&psi,const arma::field<arma::mat>&blC,
                         const double&r){
  arma::vec result(n);
  int order=rules.n_rows;
  arma::vec functionvalue(order);
  for(int i=0;i<n;i++){
    for(int k=0;k<order;k++){
      functionvalue(k)=g0b(rules(k,0),beta,gamma,eta,theta,betadim,gammadim,Delta(i),X(i),trans(Z.row(i)),psi,blC(i),r,ni(i));
    }
    result(i)=sum(functionvalue%rules.col(1));
  }
  return result;
}




arma::vec ui(const double&b,const arma::vec&beta,const arma::vec&gamma,const arma::vec&eta,const double&theta,
             const int&betadim,const int&gammadim,const arma::vec&Delta,
             const arma::mat&X,const arma::vec&Z,
             const arma::vec&psi,const arma::mat&blC,const double&r,const int&ni){
  double eta0=eta(0);
  arma::vec eta1=eta.subvec(1,betadim);
  arma::vec eta2=eta.subvec(betadim+1,betadim+gammadim);
  arma::vec pilinear=eta0+X*eta1+sum(trans(eta2)*Z);
  arma::vec pi=pow((1+exp(-pilinear)),-1);
  arma::vec S;
  arma::vec Slinear=X*beta+sum(trans(gamma)*Z)+std::exp(theta)*b;
  arma::vec H_t=trans(blC)*exp(psi);
  if(r==0){
    S=exp(-H_t%exp(Slinear));
  }
  else{
    S=pow(1+r*H_t%exp(Slinear),-1/r);
  }
  arma::vec result;
  result=(1-Delta)%pi%pow(pi+(1-pi)%S,-1);
  return result;
}

// [[Rcpp::export]]
arma::field<arma::mat> uijq(const arma::mat&rules,const arma::vec&beta,const arma::vec&gamma,const arma::vec&eta,const double&theta,
                            const arma::field<arma::vec>&Delta,const int&betadim,const int&gammadim,
                            const int&n,const arma::vec&ni,
                            const arma::field<arma::mat>&X,
                            const arma::mat&Z,const arma::vec&psi,const arma::field<arma::mat>&blC,
                            const double&r){
  int order=rules.n_rows;
  arma::field<arma::mat> result(n);
  for(int i=0;i<n;i++){
    result(i).zeros(ni(i),order);
    for(int q=0;q<order;q++){
      result(i).col(q)=ui(rules(q,0),beta,gamma,eta,theta,betadim,gammadim,Delta(i),X(i),trans(Z.row(i)),psi,blC(i),r,ni(i));
    }
  }
  return result;
}

// [[Rcpp::export]]

double Q1(const arma::vec&parameters,const arma::mat&rules,const arma::vec&beta,const arma::vec&gamma,const arma::vec&eta,const double&theta,
          const arma::field<arma::vec>&Delta,const int&betadim,const int&gammadim,
          const int&n,const arma::vec&ni,
          const arma::field<arma::mat>&X,
          const arma::mat&Z,const arma::vec&psi,const arma::field<arma::mat>&blC,
          const double&r,const arma::mat&g0bvalue,const arma::vec&normalconstant,
          const arma::field<arma::mat>&uijqmat,const bool&Cauchyindex){
  arma::field<arma::vec> uij(n);
  double esteta0=parameters(0);
  arma::vec esteta1=parameters.subvec(1,betadim);
  arma::vec esteta2=parameters.subvec(betadim+1,betadim+gammadim);
  arma::vec pilinear;
  arma::vec pi;
  for(int i=0;i<n;i++){
    
    uij(i).zeros(ni(i));
    uij(i)=uijqmat(i)*(g0bvalue.col(i)%rules.col(1))/normalconstant(i);
  }
  double result=0;
  for(int i=0;i<n;i++){
    pilinear=esteta0+X(i)*esteta1+sum(trans(esteta2)*trans(Z.row(i)));
    pi=pow((1+exp(-pilinear)),-1);
    result=result+sum(uij(i)%log(pi)+(1-uij(i))%log(1-pi));
  }
  arma::vec parametersquare=pow(parameters,2);
  if(Cauchyindex){
    result=result-sum(log(1+parametersquare/6.25));
  }
  
  
  
  return -result;
}




// [[Rcpp::export]]
double Q2(const arma::vec&parameters,const arma::mat&rules,const arma::vec&beta,const arma::vec&gamma,const arma::vec&eta,const double&theta,
          const arma::field<arma::vec>&Delta,const int&betadim,const int&gammadim,
          const int&n,const arma::vec&ni,
          const arma::field<arma::mat>&X,
          const arma::mat&Z,const arma::vec&psi,const arma::field<arma::mat>&blC,
          const double&r,const arma::mat&g0bvalue,const arma::vec&normalconstant,
          const arma::field<arma::mat>&uijqmat,const double&tunepar,const arma::mat&R,const bool&Cauchyindex){
  arma::vec estbeta=parameters.subvec(0,betadim-1);
  arma::vec estgamma=parameters.subvec(betadim,betadim+gammadim-1);
  double esttheta=parameters(betadim+gammadim);
  arma::vec estpsi=parameters.subvec(betadim+gammadim+1,parameters.n_elem-1);
  arma::field<arma::mat> Smat(n);
  arma::vec Slinear;
  arma::vec H_t;
  int order=rules.n_rows;
  for(int i=0;i<n;i++){
    Smat(i).zeros(ni(i),order);
    H_t=trans(blC(i))*exp(estpsi);
    for(int q=0;q<order;q++){
      Slinear=X(i)*estbeta+sum(trans(estgamma)*trans(Z.row(i)))+std::exp(esttheta)*rules(q,0);
      if(r==0){
        Smat(i).col(q)=exp(-H_t%exp(Slinear));
      }
      else{
        Smat(i).col(q)=pow(1+r*H_t%exp(Slinear),-1/r);
      }
    }
  }
  for(int i=0;i<n;i++){
    for(int j=0;j<ni(i);j++){
      for(int q=0;q<order;q++){
        if(Smat(i)(j,q)>0.99999999999){
          Smat(i)(j,q)=0.99999999999;
        }
        else if(Smat(i)(j,q)<0.00000000001){
          Smat(i)(j,q)=0.00000000001;
        }
      }
    }
  }
  double part1=0;
  for(int i=0;i<n;i++){
    part1=part1+sum(Delta(i)%(log(1-Smat(i))*(g0bvalue.col(i)%rules.col(1))))/normalconstant(i);
  }
  double part2=0;
  for(int i=0;i<n;i++){
    part2=part2+sum((1-Delta(i))%((1-uijqmat(i))%log(Smat(i))*(g0bvalue.col(i)%rules.col(1))))/normalconstant(i);
  }
  double result;
  arma::vec parametersquare=pow(parameters.subvec(0,betadim+gammadim),2);
  arma::mat penalty;
  penalty=trans(estpsi)*R*estpsi;
  if(Cauchyindex){
    result=part1+part2-sum(log(1+parametersquare/6.25))-tunepar*penalty(0,0);
  }
  else{
    result=part1+part2-tunepar*penalty(0,0);
  }
  
  return -(result);
}







// [[Rcpp::export]]
double likelihoodfunc_i(const double&b,const arma::vec&parameters,
                       const arma::vec&Delta,const arma::mat&X,const arma::vec&Z,
                       const int&ni,const double&r,const arma::mat&blC,const int&betadim,const int&gammadim){
  int totaldim=parameters.n_elem;
  int zetadim=betadim+gammadim+1;
  double result=0;double S;arma::mat midresult;midresult.zeros(1,1);
  arma::vec covariate(betadim+gammadim);
  
  
  double uncurerate;
  for(int j=0;j<ni;j++){
    covariate.subvec(0,betadim-1)=trans(X.row(j));
    covariate.subvec(betadim,betadim+gammadim-1)=Z;
    
    uncurerate=1/(1+std::exp(-parameters(0)-sum(parameters.subvec(1,zetadim-1)%covariate)));
    
    if(r>0){
      S=pow(1+r*sum(trans(exp(parameters.subvec(zetadim+betadim+gammadim+1,totaldim-1)))*blC.col(j))
              *std::exp(sum(parameters.subvec(zetadim,zetadim+betadim+gammadim-1)%covariate)+std::exp(parameters(zetadim+betadim+gammadim))*b),-1/r);
    }
    else{
      S=std::exp(-sum(trans(exp(parameters.subvec(zetadim+betadim+gammadim+1,totaldim-1)))*blC.col(j))
                   *std::exp(sum(parameters.subvec(zetadim,zetadim+betadim+gammadim-1)%covariate)+std::exp(parameters(zetadim+betadim+gammadim))*b));
    }
    if(S>0.99999999999){
      S=0.99999999999;
    }
    if(uncurerate<std::pow(10,-30)){
      uncurerate=std::pow(10,-30);
    }
    double secondterm=uncurerate+(1-uncurerate)*S;
    if((uncurerate+(1-uncurerate)*S)<std::pow(10,-30)){
      
      secondterm=std::pow(10,-30);
    }
    result=result+Delta(j)*(log(1-S)+log(1-uncurerate))+(1-Delta(j))*std::log(secondterm);
  }
  
  // lambda=std::exp(parameters(zetadim+betadim)+sum(parameters.subvec(1+zetadim+betadim,betadim+etadim-2+zetadim)%Z)+parameters(zetadim+betadim+etadim-1)*(b));
  // for(int l=0;l<8;l++){
  //   
  //   tracprobability=tracprobability+std::pow(lambda,l+1)*inverfac(l);
  // }
  result=result+R::dnorm(b,0,1,true);
  result=std::exp(result);
  // result=result*std::pow(lambda,ni)/tracprobability*inverfac(ni-1);
  result=result*std::exp(b*b);
  return result;
}

// This function is used to caculate the likelihood value using Hermit quadrature.

// [[Rcpp::export]]
double Loglikelihood(const arma::vec&parameters,const arma::mat&rules,const arma::field<arma::vec>&Delta,
                       const arma::field<arma::vec>&X,const arma::mat&Z,const int&n,const arma::vec&ni,
                       const double&r,const arma::field<arma::mat>&blC,const int&betadim,const int&gammadim,
                       const double&tunepar,const arma::mat&R,const bool&Cauchyindex){
  int zetadim=betadim+gammadim+1;
  int order=rules.n_rows;double result=0;
  arma::vec weightvec;weightvec=rules.col(1);double term1;
  arma::vec functionvalue(order);
  arma::vec estpsi=parameters.subvec(zetadim+betadim+gammadim+1,parameters.n_elem-1);
  for(int i=0;i<n;i++){
    
    for(int k=0;k<order;k++){
      functionvalue(k)=likelihoodfunc_i(rules(k,0),parameters,
                    Delta(i),X(i),trans(Z.row(i)),ni(i),r,blC(i),betadim,gammadim);
    }
    
    term1=sum(functionvalue%weightvec);
    if(term1<std::pow(10,-30)){
      term1=std::pow(10,-30);}
    result=result+std::log(term1);
    
  }
  arma::mat penalty;
  penalty=trans(estpsi)*R*estpsi;
  arma::vec parametersquare=pow(parameters.subvec(0,zetadim+betadim+gammadim),2);
  if(Cauchyindex){
    result=result-sum(log(1+parametersquare/6.25))-tunepar*penalty(0,0);
  }
  else{
    result=result-tunepar*penalty(0,0);
  }
  return result;
}


// [[Rcpp::export]]
double loglikelihoodtest(const arma::vec&parameters,const double&b,
                         const arma::vec&Delta,const arma::mat&X,const arma::vec&Z,
                         const int&ni,const double&r,const arma::mat&blC,const int&betadim,const int&gammadim){
  int totaldim=parameters.n_elem;
  int zetadim=betadim+gammadim+1;
  double result=0;double S;arma::mat midresult;midresult.zeros(1,1);
  arma::vec covariate(betadim+gammadim);
  arma::vec facvec(8);facvec(0)=1;
  arma::vec inverfac(8);inverfac(0)=1;
  
  
  double uncurerate;
  for(int j=0;j<ni;j++){
    covariate.subvec(0,betadim-1)=trans(X.row(j));
    covariate.subvec(betadim,betadim+gammadim-1)=Z;
    
    uncurerate=1/(1+std::exp(-parameters(0)-sum(parameters.subvec(1,zetadim-1)%covariate)));
    
    if(r>0){
      S=pow(1+r*sum(trans(exp(parameters.subvec(zetadim+betadim+gammadim+1,totaldim-1)))*blC.col(j))
              *std::exp(sum(parameters.subvec(zetadim,zetadim+betadim+gammadim-1)%covariate)+std::exp(parameters(zetadim+betadim+gammadim))*b),-1/r);
    }
    else{
      S=std::exp(-sum(trans(exp(parameters.subvec(zetadim+betadim+gammadim+1,totaldim-1)))*blC.col(j))
                   *std::exp(sum(parameters.subvec(zetadim,zetadim+betadim+gammadim-1)%covariate)+std::exp(parameters(zetadim+betadim+gammadim))*b));
    }
    if(S>0.99999999999){
      S=0.99999999999;
    }
    // if(uncurerate<std::pow(10,-30)){
    //   uncurerate=std::pow(10,-30);
    // }
    double secondterm=uncurerate+(1-uncurerate)*S;
    if((uncurerate+(1-uncurerate)*S)<std::pow(10,-30)){
      
      secondterm=std::pow(10,-30);
    }
    result=result+Delta(j)*(log(1-S)+log(1-uncurerate))+(1-Delta(j))*std::log(secondterm);
  }
  
  return result;
}


// [[Rcpp::export]]
double logcauchy(const arma::vec&parameters){
  arma::vec parametersquare=pow(parameters.subvec(0,9),2);
  double result;
  result=-sum(log(1+parametersquare/6.25));
  return result;
}

// [[Rcpp::export]]
double penaltyterm(const arma::vec&psi,const double&lambda,const arma::mat&R){
  arma::mat penalty;
  penalty=exp(trans(psi))*R*exp(psi);
  double result;
  result=lambda*penalty(0,0);
  return(result);
}


// [[Rcpp::export]]

arma::vec firstderiv_icure(const arma::vec&parameters,const double&b,const arma::vec&Delta,
                           const arma::mat&X,const arma::vec&Z,const int&ni,const double&r,const arma::mat&blC,
                           const int&betadim,const int&gammadim){
  arma::vec result(parameters.n_elem);
  arma::vec covariate(betadim+gammadim);
  double eta0=parameters(0);
  arma::vec beta=parameters.subvec(1,betadim);
  arma::vec gamma=parameters.subvec(1+betadim,betadim+gammadim);
  double theta=parameters(betadim+gammadim+1);
  arma::vec psi=parameters.subvec(betadim+gammadim+2,parameters.n_elem-1);
  arma::vec onesvector;
  arma::vec pi,si,expterm,Ht;
  onesvector.ones(ni);
  arma::vec exppilinear=exp(-onesvector*eta0);
  pi=pow(1+exppilinear,-1);
  expterm=exp(X*beta+sum(Z%gamma)*onesvector+std::exp(theta)*b*onesvector);
  Ht=trans(blC)*exp(psi);
  arma::vec pisi;
  
  if(r>0){
    si=pow(1+r*Ht%expterm,-1/r);
  }
  else{
    si=exp(-Ht%expterm);
  }
  pisi=pi+(onesvector-pi)%si;
  arma::vec piderivpart1,piderivpart2;
  piderivpart1=pow(1+exppilinear,-2)%exppilinear;
  arma::vec commonvec=Delta%pow(1-pi,-1)%(-piderivpart1)+(1-Delta)%pow(pi+(1-pi)%si,-1)%(1-si)%piderivpart1;
  result(0)=sum(commonvec);
  arma::vec derivsi;
  if(r>0){
    derivsi=-pow(1+r*Ht%expterm,-1/r-1)%Ht%expterm;
  }
  else{
    derivsi=-exp(-Ht%expterm)%Ht%expterm;
  }
  arma::vec commonvec2=Delta%pow(1-si,-1)%(-derivsi)+(1-Delta)%pow(pi+(1-pi)%si,-1)%(1-pi)%derivsi;
  for(int j=0;j<ni;j++){
    result.subvec(1,betadim)=result.subvec(betadim+gammadim+1,betadim+gammadim+betadim)+
      trans(X.row(j))*commonvec2(j);
    result.subvec(1+betadim,betadim+gammadim)=result.subvec(betadim+gammadim+1+betadim,betadim+gammadim+betadim+gammadim)
      +Z*commonvec2(j);
    result(betadim+gammadim+1)=sum(commonvec2*b*exp(theta));
  }
  arma::vec commonvec3,sipsideriv;
  if(r>0){
    sipsideriv=-pow(1+r*Ht%expterm,-1/r-1)%expterm;
  }
  else{
    sipsideriv=-exp(-Ht%expterm)%expterm;
  }
  commonvec3=Delta%pow(1-si,-1)%(-sipsideriv)+(1-Delta)%pow(pi+(1-pi)%si,-1)%(1-pi)%sipsideriv;
  for(int j=0;j<psi.n_elem;j++){
    result(betadim+gammadim+2+j)=sum(exp(psi(j))*commonvec3%trans(blC.row(j)));
  }
  
  return result;
  
}


// [[Rcpp::export]]
double g0bcure(const double&b,const arma::vec&beta,const arma::vec&gamma,const double&eta0,const double&theta,
               const int&betadim,const int&gammadim,const arma::vec&Delta,
               const arma::mat&X,const arma::vec&Z,
               const arma::vec&psi,const arma::mat&blC,const double&r,const int&ni){
  
  arma::vec onevector;
  onevector.ones(ni);
  arma::vec pilinear=eta0*onevector;
  
  arma::vec pi=pow((1+exp(-pilinear)),-1);
  arma::vec S;
  arma::vec Slinear=X*beta+sum(trans(gamma)*Z)+std::exp(theta)*b;
  
  arma::vec H_t=trans(blC)*exp(psi);
  if(r==0){
    S=exp(-H_t%exp(Slinear));
  }
  else{
    S=pow(1+r*H_t%exp(Slinear),-1/r);
  }
  double result=1;
  for(int i=0;i<ni;i++){
    result=result*std::pow((1-pi(i))*(1-S(i)),Delta(i))*std::pow(pi(i)+(1-pi(i))*S(i),1-Delta(i));
  }
  result=result*R::dnorm4(b,0,1,false)*std::exp(b*b);
  return result;
}
// [[Rcpp::export]]
arma::mat g0bmatcure(const arma::mat&rules,const arma::vec&beta,const arma::vec&gamma,const double&eta0,const double&theta,
                     const arma::field<arma::vec>&Delta,const int&betadim,const int&gammadim,
                     const int&n,const arma::vec&ni,
                     const arma::field<arma::mat>&X,
                     const arma::mat&Z,const arma::vec&psi,const arma::field<arma::mat>&blC,
                     const double&r){
  int order=rules.n_rows;
  arma::mat result(order,n);
  for(int q=0;q<order;q++){
    for(int i=0;i<n;i++){
      result(q,i)=g0bcure(rules(q,0),beta,gamma,eta0,theta,betadim,gammadim,Delta(i),X(i),trans(Z.row(i)),psi,blC(i),r,ni(i));
    }
  }
  return result;
}




// [[Rcpp::export]]

arma::vec NormalConstantcure(const arma::mat&rules,const arma::vec&beta,const arma::vec&gamma,const double&eta0,const double&theta,
                             const arma::field<arma::vec>&Delta,const int&betadim,const int&gammadim,
                             const int&n,const arma::vec&ni,
                             const arma::field<arma::mat>&X,
                             const arma::mat&Z,const arma::vec&psi,const arma::field<arma::mat>&blC,
                             const double&r){
  arma::vec result(n);
  int order=rules.n_rows;
  arma::vec functionvalue(order);
  for(int i=0;i<n;i++){
    for(int k=0;k<order;k++){
      functionvalue(k)=g0bcure(rules(k,0),beta,gamma,eta0,theta,betadim,gammadim,Delta(i),X(i),trans(Z.row(i)),psi,blC(i),r,ni(i));
    }
    result(i)=sum(functionvalue%rules.col(1));
  }
  return result;
}



// [[Rcpp::export]]
arma::vec uicure(const double&b,const arma::vec&beta,const arma::vec&gamma,const double&eta0,const double&theta,
                 const int&betadim,const int&gammadim,const arma::vec&Delta,
                 const arma::mat&X,const arma::vec&Z,
                 const arma::vec&psi,const arma::mat&blC,const double&r,const int&ni){
  
  arma::vec onevector;
  onevector.ones(ni);
  arma::vec pilinear=eta0*onevector;
  arma::vec pi=pow((1+exp(-pilinear)),-1);
  arma::vec S;
  arma::vec Slinear=X*beta+sum(trans(gamma)*Z)+std::exp(theta)*b;
  arma::vec H_t=trans(blC)*exp(psi);
  if(r==0){
    S=exp(-H_t%exp(Slinear));
  }
  else{
    S=pow(1+r*H_t%exp(Slinear),-1/r);
  }
  arma::vec result;
  result=(1-Delta)%pi%pow(pi+(1-pi)%S,-1);
  return result;
}

// [[Rcpp::export]]
arma::field<arma::mat> uijqcure(const arma::mat&rules,const arma::vec&beta,const arma::vec&gamma,const double&eta0,const double&theta,
                                const arma::field<arma::vec>&Delta,const int&betadim,const int&gammadim,
                                const int&n,const arma::vec&ni,
                                const arma::field<arma::mat>&X,
                                const arma::mat&Z,const arma::vec&psi,const arma::field<arma::mat>&blC,
                                const double&r){
  int order=rules.n_rows;
  arma::field<arma::mat> result(n);
  for(int i=0;i<n;i++){
    result(i).zeros(ni(i),order);
    for(int q=0;q<order;q++){
      result(i).col(q)=uicure(rules(q,0),beta,gamma,eta0,theta,betadim,gammadim,Delta(i),X(i),trans(Z.row(i)),psi,blC(i),r,ni(i));
    }
  }
  return result;
}

// [[Rcpp::export]]

double Q1cure(const double&parameters,const arma::mat&rules,const arma::vec&beta,const arma::vec&gamma,const arma::vec&eta,const double&theta,
              const arma::field<arma::vec>&Delta,const int&betadim,const int&gammadim,
              const int&n,const arma::vec&ni,
              const arma::field<arma::mat>&X,
              const arma::mat&Z,const arma::vec&psi,const arma::field<arma::mat>&blC,
              const double&r,const arma::mat&g0bvalue,const arma::vec&normalconstant,
              const arma::field<arma::mat>&uijqmat,const bool&Cauchyindex){
  arma::field<arma::vec> uij(n);
  
  arma::vec pilinear;
  arma::vec pi;
  arma::vec onevector;
  for(int i=0;i<n;i++){
    
    uij(i).zeros(ni(i));
    uij(i)=uijqmat(i)*(g0bvalue.col(i)%rules.col(1))/normalconstant(i);
  }
  double result=0;
  for(int i=0;i<n;i++){
    onevector.ones(ni(i));
    pilinear=parameters*onevector;
    pi=pow((1+exp(-pilinear)),-1);
    result=result+sum(uij(i)%log(pi)+(1-uij(i))%log(1-pi));
  }
  if(Cauchyindex){
    result=result-(log(1+std::pow(parameters,2)/6.25));
  }
  return -result;
}

// [[Rcpp::export]]
double Q2cure(const arma::vec&parameters,const arma::mat&rules,const arma::vec&beta,const arma::vec&gamma,const arma::vec&eta,const double&theta,
              const arma::field<arma::vec>&Delta,const int&betadim,const int&gammadim,
              const int&n,const arma::vec&ni,
              const arma::field<arma::mat>&X,
              const arma::mat&Z,const arma::vec&psi,const arma::field<arma::mat>&blC,
              const double&r,const arma::mat&g0bvalue,const arma::vec&normalconstant,
              const arma::field<arma::mat>&uijqmat,const double&tunepar,const arma::mat&R,const bool&Cauchyindex){
  arma::vec estbeta=parameters.subvec(0,betadim-1);
  arma::vec estgamma=parameters.subvec(betadim,betadim+gammadim-1);
  double esttheta=parameters(betadim+gammadim);
  arma::vec estpsi=parameters.subvec(betadim+gammadim+1,parameters.n_elem-1);
  arma::field<arma::mat> Smat(n);
  arma::vec Slinear;
  arma::vec H_t;
  int order=rules.n_rows;
  for(int i=0;i<n;i++){
    Smat(i).zeros(ni(i),order);
    H_t=trans(blC(i))*exp(estpsi);
    for(int q=0;q<order;q++){
      Slinear=X(i)*estbeta+sum(trans(estgamma)*trans(Z.row(i)))+std::exp(esttheta)*rules(q,0);
      if(r==0){
        Smat(i).col(q)=exp(-H_t%exp(Slinear));
      }
      else{
        Smat(i).col(q)=pow(1+r*H_t%exp(Slinear),-1/r);
      }
    }
  }
  for(int i=0;i<n;i++){
    for(int j=0;j<ni(i);j++){
      for(int q=0;q<order;q++){
        if(Smat(i)(j,q)>0.99999999999){
          Smat(i)(j,q)=0.99999999999;
        }
        else if(Smat(i)(j,q)<0.00000000001){
          Smat(i)(j,q)=0.00000000001;
        }
      }
    }
  }
  double part1=0;
  for(int i=0;i<n;i++){
    part1=part1+sum(Delta(i)%(log(1-Smat(i))*(g0bvalue.col(i)%rules.col(1))))/normalconstant(i);
  }
  double part2=0;
  for(int i=0;i<n;i++){
    part2=part2+sum((1-Delta(i))%((1-uijqmat(i))%log(Smat(i))*(g0bvalue.col(i)%rules.col(1))))/normalconstant(i);
  }
  double result;
  arma::vec parametersquare=pow(parameters.subvec(0,betadim+gammadim),2);
  arma::mat penalty;
  penalty=trans(estpsi)*R*estpsi;
  if(Cauchyindex){
    result=part1+part2-sum(log(1+parametersquare/6.25))-tunepar*penalty(0,0);
  }
  else{
    result=part1+part2-tunepar*penalty(0,0);
  }
  
  return -(result);
}






// [[Rcpp::export]]
double likelihoodfunccure(const double&b,const arma::vec&parameters,
                          const arma::vec&Delta,const arma::mat&X,const arma::vec&Z,
                          const int&ni,const double&r,const arma::mat&blC,const int&betadim,const int&gammadim){
  int totaldim=parameters.n_elem;
  int zetadim=1;
  double result=0;double S;arma::mat midresult;midresult.zeros(1,1);
  arma::vec covariate(betadim+gammadim);
  
  
  double uncurerate;
  for(int j=0;j<ni;j++){
    covariate.subvec(0,betadim-1)=trans(X.row(j));
    covariate.subvec(betadim,betadim+gammadim-1)=Z;
    
    uncurerate=1/(1+std::exp(-parameters(0)));
    
    if(r>0){
      S=pow(1+r*sum(trans(exp(parameters.subvec(zetadim+betadim+gammadim+1,totaldim-1)))*blC.col(j))
              *std::exp(sum(parameters.subvec(zetadim,zetadim+betadim+gammadim-1)%covariate)+std::exp(parameters(zetadim+betadim+gammadim))*b),-1/r);
    }
    else{
      S=std::exp(-sum(trans(exp(parameters.subvec(zetadim+betadim+gammadim+1,totaldim-1)))*blC.col(j))
                   *std::exp(sum(parameters.subvec(zetadim,zetadim+betadim+gammadim-1)%covariate)+std::exp(parameters(zetadim+betadim+gammadim))*b));
    }
    if(S>0.99999999999){
      S=0.99999999999;
    }
    if(uncurerate<std::pow(10,-30)){
      uncurerate=std::pow(10,-30);
    }
    double secondterm=uncurerate+(1-uncurerate)*S;
    if((uncurerate+(1-uncurerate)*S)<std::pow(10,-30)){
      
      secondterm=std::pow(10,-30);
    }
    result=result+Delta(j)*(log(1-S)+log(1-uncurerate))+(1-Delta(j))*std::log(secondterm);
  }
  
  // lambda=std::exp(parameters(zetadim+betadim)+sum(parameters.subvec(1+zetadim+betadim,betadim+etadim-2+zetadim)%Z)+parameters(zetadim+betadim+etadim-1)*(b));
  // for(int l=0;l<8;l++){
  //   
  //   tracprobability=tracprobability+std::pow(lambda,l+1)*inverfac(l);
  // }
  result=result+R::dnorm(b,0,1,true);
  result=std::exp(result);
  // result=result*std::pow(lambda,ni)/tracprobability*inverfac(ni-1);
  result=result*std::exp(b*b);
  return result;
}

// This function is used to caculate the likelihood value using Hermit quadrature.

// [[Rcpp::export]]
double Loglikelihoodcure(const arma::vec&parameters,const arma::mat&rules,const arma::field<arma::vec>&Delta,
                         const arma::field<arma::vec>&X,const arma::mat&Z,const int&n,const arma::vec&ni,
                         const double&r,const arma::field<arma::mat>&blC,const int&betadim,const int&gammadim,
                         const double&tunepar,const arma::mat&R,const bool&Cauchyindex){
  int zetadim=1;
  int order=rules.n_rows;double result=0;
  arma::vec weightvec;weightvec=rules.col(1);double term1;
  arma::vec functionvalue(order);
  arma::vec estpsi=parameters.subvec(zetadim+betadim+gammadim+1,parameters.n_elem-1);
  for(int i=0;i<n;i++){
    
    for(int k=0;k<order;k++){
      functionvalue(k)=likelihoodfunccure(rules(k,0),parameters,Delta(i),X(i),trans(Z.row(i)),ni(i),r,blC(i),betadim,gammadim);
    }
    
    term1=sum(functionvalue%weightvec);
    if(term1<std::pow(10,-30)){
      term1=std::pow(10,-30);}
    result=result+std::log(term1);
    
  }
  arma::mat penalty;
  penalty=trans(estpsi)*R*estpsi;
  arma::vec parametersquare=pow(parameters.subvec(0,zetadim+betadim+gammadim),2);
  if(Cauchyindex){
    result=result-sum(log(1+parametersquare/6.25))-tunepar*penalty(0,0);
  }
  else{
    result=result-tunepar*penalty(0,0);
  }
  return result;
}








// [[Rcpp::export]]
double loglikelihoodtestcure(const arma::vec&parameters,const double&b,
                             const arma::vec&Delta,const arma::mat&X,const arma::vec&Z,
                             const int&ni,const double&r,const arma::mat&blC,const int&betadim,const int&gammadim){
  int totaldim=parameters.n_elem;
  int zetadim=1;
  double result=0;double S;arma::mat midresult;midresult.zeros(1,1);
  arma::vec covariate(betadim+gammadim);
  arma::vec facvec(8);facvec(0)=1;
  arma::vec inverfac(8);inverfac(0)=1;
  
  
  double uncurerate;
  for(int j=0;j<ni;j++){
    covariate.subvec(0,betadim-1)=trans(X.row(j));
    covariate.subvec(betadim,betadim+gammadim-1)=Z;
    
    uncurerate=1/(1+std::exp(-parameters(0)));    
    if(r>0){
      S=pow(1+r*sum(trans(exp(parameters.subvec(zetadim+betadim+gammadim+1,totaldim-1)))*blC.col(j))
              *std::exp(sum(parameters.subvec(zetadim,zetadim+betadim+gammadim-1)%covariate)+std::exp(parameters(zetadim+betadim+gammadim))*b),-1/r);
    }
    else{
      S=std::exp(-sum(trans(exp(parameters.subvec(zetadim+betadim+gammadim+1,totaldim-1)))*blC.col(j))
                   *std::exp(sum(parameters.subvec(zetadim,zetadim+betadim+gammadim-1)%covariate)+std::exp(parameters(zetadim+betadim+gammadim))*b));
    }
    if(S>0.99999999999){
      S=0.99999999999;
    }
    // if(uncurerate<std::pow(10,-30)){
    //   uncurerate=std::pow(10,-30);
    // }
    double secondterm=uncurerate+(1-uncurerate)*S;
    if((uncurerate+(1-uncurerate)*S)<std::pow(10,-30)){
      
      secondterm=std::pow(10,-30);
    }
    result=result+Delta(j)*(log(1-S)+log(1-uncurerate))+(1-Delta(j))*std::log(secondterm);
  }
  
  return result;
}



// [[Rcpp::export]]
double likelihoodfunc1(const double&b,const arma::vec&parameters,
                       const arma::vec&Delta,const arma::mat&X,const arma::vec&Z,
                       const int&ni,const double&r,const arma::mat&blC,const int&betadim,const int&gammadim){
  int totaldim=parameters.n_elem;
  int zetadim=betadim+gammadim+1;
  double result=0;double S,lambda;arma::mat midresult;midresult.zeros(1,1);
  double tracprobability=0;
  arma::vec covariate(betadim+gammadim);
  
  
  double uncurerate;
  for(int j=0;j<ni;j++){
    covariate.subvec(0,betadim-1)=trans(X.row(j));
    covariate.subvec(betadim,betadim+gammadim-1)=Z;
    
    uncurerate=1/(1+std::exp(-parameters(0)-sum(parameters.subvec(1,zetadim-1)%covariate)));
    
    if(r>0){
      S=pow(1+r*sum(trans(exp(parameters.subvec(zetadim+betadim+gammadim+1,totaldim-1)))*blC.col(j))
              *std::exp(sum(parameters.subvec(zetadim,zetadim+betadim+gammadim-1)%covariate)+std::exp(parameters(zetadim+betadim+gammadim))*b),-1/r);
    }
    else{
      S=std::exp(-sum(trans(exp(parameters.subvec(zetadim+betadim+gammadim+1,totaldim-1)))*blC.col(j))
                   *std::exp(sum(parameters.subvec(zetadim,zetadim+betadim+gammadim-1)%covariate)+std::exp(parameters(zetadim+betadim+gammadim))*b));
    }
    if(S>0.99999999999){
      S=0.99999999999;
    }
    if(uncurerate<std::pow(10,-30)){
      uncurerate=std::pow(10,-30);
    }
    double secondterm=uncurerate+(1-uncurerate)*S;
    if((uncurerate+(1-uncurerate)*S)<std::pow(10,-30)){
      
      secondterm=std::pow(10,-30);
    }
    result=result+Delta(j)*(log(1-S)+log(1-uncurerate))+(1-Delta(j))*std::log(secondterm);
  }
  
  
  result=result+R::dnorm(b,0,1,true);
  result=std::exp(result);
  result=result*std::exp(b*b);
  return result;
}


// [[Rcpp::export]]
double testquadrature1(const arma::vec&parameters,const arma::mat&rules,const arma::field<arma::vec>&Delta,
                       const arma::field<arma::vec>&X,const arma::mat&Z,const int&n,const arma::vec&ni,
                       const double&r,const arma::field<arma::mat>&blC,const int&betadim,const int&gammadim,
                       const arma::vec&weight){
  int zetadim=betadim+gammadim+1;
  int totaldim=parameters.n_elem;
  int order=rules.n_rows;double result=0;
  arma::vec weightvec;weightvec=rules.col(1);double term1;
  arma::vec functionvalue(order);
  for(int i=0;i<n;i++){
    
    for(int k=0;k<order;k++){
      functionvalue(k)=likelihoodfunc1(rules(k,0),parameters,
                    Delta(i),X(i),trans(Z.row(i)),ni(i),r,blC(i),betadim,gammadim);
    }
    
    term1=sum(functionvalue%weightvec);
    if(term1<std::pow(10,-30)){
      term1=std::pow(10,-30);}
    result=result+weight(i)*std::log(term1);
    
  }
  arma::vec parametersquare=pow(parameters.subvec(0,betadim+gammadim+zetadim),2);
  result=result-sum(log(1+parametersquare/6.25));
  return -result;
}



// [[Rcpp::export]]
arma::field<arma::mat> Maxeigen(const arma::mat&B){
  arma::vec eigval;
  arma::mat eigvec;
  arma::field<arma::mat> result(2);
  arma::eig_sym(eigval,eigvec,B);
  result(0)=eigval;
  result(1)=eigvec;
  return result;
}

