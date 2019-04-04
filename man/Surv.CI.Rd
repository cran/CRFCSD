\name{Surv.CI}
\alias{Surv.CI}
\title{
Predicating the survival probability and the corresponding confidence interval.
}
\description{
Surv.CI can predicate the survival probability and the corresponding confidence interval for a given set of time points and a given set of covariate values based on the bootstrap method. This uses the ouput of the function boot.CSD as an input. 
}
\usage{
Surv.CI(t,x,boot.par,Rawdata,n_subject,n_within,r,
CI.lev=0.95,n_quad=30,knots.num=2,degree=2,
lambda=0,Cauchy.pen=TRUE,tolerance=1e-2,scale.numr=TRUE,cure.reg=TRUE)
}
\arguments{
  \item{t}{The given set of time points. All time points should be between the minimum and maximum of the inspection time.}
  
  \item{x}{The given set of covariate vector. The first n_subject.raw covariates are subject (cluster) specific covariates. The rest of it are the within cluster covariates.}
  
  \item{boot.par}{The bootstrap estimation dataframe provided by the function boot.CSD}
  
  \item{Rawdata}{This is a dataframe of the current status data. The first column should be the index of the subject (cluster). The second column 
is the inspection time. The next n_subjec.raw columns are the subject (cluster-specifie) level covariates. Then the next n_within.raw columns are the within subject
covariates. The last column is the indicator of the event where  1 or 0 indicate if the event has or has not happened by the inspection time,  respectively. All the covariates are assumed to be either numerical or binary, and  our program automatically detects if a covariate is a binary or numerical variable.}

\item{n_subject}{The number of subject (cluster-specifie) level covariates.}

\item{n_within}{The number of within cluster covariates.}

\item{r}{The index of the Generalized odds ratio (GOR) model. This index is a non-negative number and it must be specified by the user. 
Here r=0 and 1 imply the proportional hazard and the proportional odds model, respectively. }

\item{CI.lev}{The confidence level. The default value is 0.95.}

\item{n_quad}{The number of Gauss-Hermite quadrature nodes used in numerical integration. The default value is 30.}

\item{lambda}{The tuning parameter of the roughness penalty used for estimating the non-parametric component of the GOR model. The default value is 0. One must use the roughness penalty when the number of basis functions in the non-parametric component of the GOR model is large.}

\item{Cauchy.pen}{logical. If TRUE, then we use Cauchy penalty on the regression parameters to reduce the samll sample bias. The default is TRUE. }

\item{tolerance}{This denotes the summation of the absolute values of the relative tolerance of all parameters in the model. It is used to define the convergence of the parameter estimates. The default value is 0.01.}

\item{knots.num}{The number of equidistant interior knots for the integrated B-spline approximation of the nonparametric component of the GOR model. The default value is 2.}

\item{degree}{The degree of integrated B-splines. The default value is 2.}

\item{scale.numr}{logical. If TRUE, then all numeric covariates (cluster specifie and within cluster) are scaled with mean zero and standard deviation one. The default value is TRUE.}

\item{cure.reg}{logical. TRUE and FALSE indicate modelling the cure rate part with covariates (including cluster and within cluster covariates) or not, respectively. The default value is TRUE.}

}

\value{
Surv.band returns the predicated survival probability and the corresponding confidence interval.
}



\examples{
data(PD)
t=c(45,50)
x=rbind(c(1,0,1,0),c(1,1,0,0))
data(boot.result)
CI_t=Surv.CI(t,x,boot.result,PD,3,1,0,n_quad=5)
}