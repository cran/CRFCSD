\name{boot.CSD}
\alias{boot.CSD}
\title{
Bootstrap estimation for the parameters.
}
\description{
For the purpose of the confidence band calculation, here we estimate the parametric and nonparametric components of the model- mixed cure GOR model based on the bootstrap resampling.
}

\usage{
  boot.CSD(Rawdata,n_subject.raw,n_within.raw,r,boot.rep,seed.begin,
                  n_quad=30,lambda=0,Cauchy.pen=TRUE,tolerance=1e-2,
                  knots.num=2,degree=2,scale.numr=TRUE,cure.reg=TRUE,clustering=TRUE)
}
\arguments{




\item{Rawdata}{This is a dataframe of the current status data. The first column should be the index of the subject (cluster). The second column 
is the inspection time. The next n_subjec.raw columns are the subject (cluster-specifie) level covariates. Then the next n_within.raw columns are the within subject
covariates. The last column is the indicator of the event where  1 or 0 indicate if the event has or has not happened by the inspection time,  respectively. All the covariates are assumed to be either numerical or binary, and  our program automatically detects if a covariate is a binary or numerical variable.}

\item{n_subject.raw}{The number of subject (cluster-specifie) level covariates.}

\item{n_within.raw}{The number of within cluster covariates.}

\item{r}{The index of the Generalized odds ratio (GOR) model. This index is a non-negative number and it must be specified by the user. 
Here r=0 and 1 imply the proportional hazard and the proportional odds model, respectively. }

\item{boot.rep}{The replications for the bootstrap (the number of bootstrap samples).}

\item{seed.begin}{The random seed to start the bootstrap sampling. If the user wants to use parallel computation, 
just set different seed.begin for each nodes, then combine the results.}

\item{n_quad}{The number of Gauss-Hermite quadrature nodes used in numerical integration. The default value is 30.}

\item{lambda}{The tuning parameter of the roughness penalty used for estimating the non-parametric component of the GOR model. The default value is 0. One must use the roughness penalty when the number of basis functions in the non-parametric component of the GOR model is large.}

\item{Cauchy.pen}{logical. If TRUE, then we use Cauchy penalty on the regression parameters to reduce the samll sample bias. The default is TRUE. }

\item{tolerance}{This denotes the summation of the absolute values of the relative tolerance of all parameters in the model. It is used to define the convergence of the parameter estimates. The default value is 0.01.}

\item{knots.num}{The number of equidistant interior knots for the integrated B-spline approximation of the nonparametric component of the GOR model. The default value is 2.}

\item{degree}{The degree of integrated B-splines. The default value is 2.}

\item{scale.numr}{logical. If TRUE, then all numeric covariates (cluster specifie and within cluster) are scaled with mean zero and standard deviation one. The default value is TRUE.}
\item{cure.reg}{logical. TRUE and FALSE indicate modelling the cure rate part with covariates (including cluster and within cluster covariates) or not, respectively. The default value is TRUE.}
   \item{clustering}{
   logical. TRUE and FALSE indicate assume there is clustering effect or not, respectively. The default value is TRUE.
   }
}
\value{
boot.CSD returns a dataframe, each row contains one parameter estimate based on one bootstrap sample.
}



\examples{
data(PD)
boot.result=boot.CSD(PD,3,1,0,3,10,n_quad=5)

}
