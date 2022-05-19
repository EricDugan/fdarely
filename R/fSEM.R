#' @title Functional Standard Error of Measurements
#'
#' @description Computation of the functional standard error of measurements (fSEM) for a multiple test session
#' scenario (with at least two sessions, test and retest), possibly in the presence of replicates. 
#'
#' @param data Data matrix collecting all individuals' curves evaluated on the same grid of \eqn{M} points. 
#' It is a matrix of dimension \eqn{J_{..}\times M}, where: 
#'     \eqn{J_{..} = \sum_{g=1}^G \sum_{i=1}^n J_{gi}}, 
#'     \eqn{G} is the number of test sessions (\eqn{G=2} in the test-retest scenario),
#'     \eqn{i} is the number of individuals,
#'     \eqn{J_{gi}} is the number of replicates for each individual and test session (\eqn{J_{gi} =1} in case of no replicates)
#'     \eqn{M} is the number of points of the grid.
#' 
#' @param individuals Factor with \eqn{J_{..}} elements coding the individuals 
#' (the number of levels correspond to the number of individuals)
#' 
#' @param test Factor with \eqn{J_{..}} elements coding the test sessions 
#' (the number of levels correspond to the number of test sessions, 2 in the test-retest setting)
#' 
#' @param alpha Significance level used for the pointwise confidence bands computation. The default is \code{alpha=0.05}.
#' 
#' @return A list containing the following components:
#' \item{SEM}{SEM value integrated over the data domain.}
#' \item{pointwise_SEM}{Pointwise SEM curve (evaluated on the same grid of \eqn{M} points provided by the user).}
#' \item{low_SEM}{Lower bound of the pointwise confidence interval for the SEM curve evaluated at level \code{alpha}.}
#' \item{up_SEM}{Upper bound of the pointwise confidence interval for the SEM curve evaluated at level \code{alpha}.}
#' 
#' @seealso See also \code{\link{fICC}} for the functional ICC. 
#'
#' @examples
#' # Examples with simulated data:
#' # EXAMPLE 1: reliable data
#' n = 10 # number of individuals
#' 
#' # parameters of the basis expansion
#' ncoef = 5 # number of basis coefficients
#' p = 100 # number of evaluation points
#' domain = seq(0,100,len=p) # domain
#' basis = create.bspline.basis(range(domain),nbasis=ncoef,norder=4) # b spline basis
#' # variability parameters:
#' sigmaW = c(0.5,1,1,1,0.5) # variability within individuals
#' sigmaB = c(2.5,5,5,5,2.5) # variability between individuals
#' 
#' # mean parameters :
#' # mean vector of basis coefficients
#' mucoef = c(0,5,0,-5,0) 
#' 
#' # generating data:
#' # coefficients of individual random effect
#' set.seed(1)
#' c_e_i = matrix(nrow=n,ncol=ncoef,data=rnorm(n*ncoef))
#' c_e_i = c_e_i*matrix(nrow=n,ncol=ncoef,data=sigmaB,byrow=TRUE) +
#'          matrix(nrow=n,ncol=ncoef,data=mucoef,byrow=TRUE)
#' 
#' # coefficients of test random effect
#' c_e_i1 = matrix(nrow=n,ncol=ncoef,data=rnorm(n*ncoef))
#' c_e_i1 = c_e_i1*matrix(nrow=n,ncol=ncoef,data=sigmaW,byrow=TRUE) 
#' 
#' # coefficients of retest random effect
#' c_e_i2 = matrix(nrow=n,ncol=ncoef,data=rnorm(n*ncoef))
#' c_e_i2 = c_e_i2*matrix(nrow=n,ncol=ncoef,data=sigmaW,byrow=TRUE) 
#' 
#' # generating functional data
#' f_e_i = fd(coef=t(c_e_i),basisobj=basis)
#' f_e_i1 = fd(coef=t(c_e_i1),basisobj=basis)
#' f_e_i2 = fd(coef=t(c_e_i2),basisobj=basis)
#' e_i = t(eval.fd(domain,f_e_i))
#' e_i1 = t(eval.fd(domain,f_e_i1))
#' e_i2 = t(eval.fd(domain,f_e_i2))
#' y1 = e_i + e_i1 # test
#' y2 = e_i + e_i2 # retest
#' 
#' data = rbind(y1,y2) # bind data
#' individual = c(1:n,1:n)
#' test = c(rep(1,n),rep(2,n))
#' 
#' # plot of simulated data
#' matplot(t(data),type='l',col=test,lty=individual)
#' 
#' # Computation of index
#' fSEM.result = fSEM(data,individual,test)
#' 
#' # Plot of pointwise index
#' plot(domain,fSEM.result$pointwise_SEM,type='l',
#'      ylim=range(c(fSEM.result$up_SEM,fSEM.result$low_SEM)))
#' lines(domain,fSEM.result$low_SEM,col=2)
#' lines(domain,fSEM.result$up_SEM,col=2)
#' 
#' 
#' 
#' # EXAMPLE 2: not reliable data (increased variability between test and retest)
#' n = 10 # number of individuals
#' 
#' # parameters of the basis expansion
#' ncoef = 5 # number of basis coefficients
#' p = 100 # number of evaluation points
#' domain = seq(0,100,len=p) # domain
#' basis = create.bspline.basis(range(domain),nbasis=ncoef,norder=4) # b spline basis
#' # variability parameters:
#' sigmaW = c(3,7,7,7,3) # variability within individuals
#' sigmaB = c(2.5,5,5,5,2.5) # variability between individuals
#' 
#' # mean parameters :
#' # mean vector of basis coefficients
#' mucoef = c(0,5,0,-5,0) 
#' 
#' # generating data:
#' set.seed(1)
#' # coefficients of individual random effect
#' c_e_i = matrix(nrow=n,ncol=ncoef,data=rnorm(n*ncoef))
#' c_e_i = c_e_i*matrix(nrow=n,ncol=ncoef,data=sigmaB,byrow=TRUE) + 
#'   matrix(nrow=n,ncol=ncoef,data=mucoef,byrow=TRUE)
#' 
#' # coefficients of test random effect
#' c_e_i1 = matrix(nrow=n,ncol=ncoef,data=rnorm(n*ncoef))
#' c_e_i1 = c_e_i1*matrix(nrow=n,ncol=ncoef,data=sigmaW,byrow=TRUE) 
#' 
#' # coefficients of retest random effect
#' c_e_i2 = matrix(nrow=n,ncol=ncoef,data=rnorm(n*ncoef))
#' c_e_i2 = c_e_i2*matrix(nrow=n,ncol=ncoef,data=sigmaW,byrow=TRUE) 
#' 
#' # generating functional data
#' f_e_i = fd(coef=t(c_e_i),basisobj=basis)
#' f_e_i1 = fd(coef=t(c_e_i1),basisobj=basis)
#' f_e_i2 = fd(coef=t(c_e_i2),basisobj=basis)
#' e_i = t(eval.fd(domain,f_e_i))
#' e_i1 = t(eval.fd(domain,f_e_i1))
#' e_i2 = t(eval.fd(domain,f_e_i2))
#' y1 = e_i + e_i1 # test
#' y2 = e_i + e_i2 # retest
#' 
#' data = rbind(y1,y2) # bind data
#' individual = c(1:n,1:n)
#' test = c(rep(1,n),rep(2,n))
#' 
#' # plot of simulated data
#' matplot(t(data),type='l',col=test,lty=individual)
#' 
#' # Computation of similarity
#' fSEM.result = fSEM(data,individual,test)
#' 
#' # Plot of pointwise index
#' plot(domain,fSEM.result$pointwise_SEM,type='l',
#'      ylim=range(c(fSEM.result$up_SEM,fSEM.result$low_SEM)))
#' lines(domain,fSEM.result$low_SEM,col=2)
#' lines(domain,fSEM.result$up_SEM,col=2)
#' 
#' @references
#' Pini, A., Markström, J., and Schelin, L. (2019): Test–retest reliability measures for curve data: 
#' an overview with recommendations and supplementary code, \emph{Sports Biomechanics} 21 (2): 179-200.
#' 
#' Schelin, L., Pini, A., Markström, J. L., Hager, C. K. (2021): Test-retest reliability of entire
#' time-series data from hip, knee and ankle kinematics and kinetics during one-leg hops
#' for distance: Analyses using integrated pointwise indices, \emph{Journal of Biomechanics}
#' Jul 19, 124:110546. 
#' 
#' Shrout, P. E., & Fleiss, J. L. (1979): Intraclass correlations: uses in assessing rater reliability. \emph{Psychological bulletin}, 86(2), 420.
#' 
#' McGraw, K. O., &  Wong S.P. (1996). Forming inferences about some intraclass correlation coefficients. \emph{Psychological methods} 1.1: 30.
#'
#' @export

fSEM = function(data,individuals,test,alpha=0.05){
  t = dim(data)[2]
  N = dim(data)[1]
  n = length(levels(factor(individuals)))
  nik = table(individuals, test)
  Ni = rowSums(nik)
  Nk = colSums(nik)
  K = length(Nk)
  
  ks = (N - sum(Ni^2)/N) / (n-1)
  kt = (N - sum(Nk^2)/N) / (K-1)
  
  pointwise_SEM = low_SEM = up_SEM = numeric(t)
  for(point in 1:t){
    data.all = data.frame(response = data[,point], individuals = factor(individuals), test=factor(test))
    anova = aov(response ~ individuals + test, data = data.all)
    
    MS = summary(anova)[[1]]$'Mean Sq'
    MSs = MS[1]
    MSt = MS[2]
    MSe = MS[3]
    
    pointwise_SEM[point] = sqrt((MSt+(kt-1)*MSe)/kt)
    c = 1/kt
    d = (kt-1)/kt
    g = (c*MSt+d*MSe)^2/((c*MSt)^2 / (K-1) + (d*MSe)^2 / (N-n-K-1) )
    
    low_SEM[point] = sqrt( g*pointwise_SEM[point]^2/qchisq(1-alpha/2,g) )
    up_SEM[point] = sqrt( g*pointwise_SEM[point]^2/qchisq(alpha/2,g) )
  }
  
  
  return(list(SEM=mean(pointwise_SEM), pointwise_SEM=pointwise_SEM,low_SEM=low_SEM,up_SEM=up_SEM))
}

