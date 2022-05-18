#' @title Functional Coefficient of Multiple Correlation
#'
#' @description Computation of the functional coefficient of multiple correlation (fCMC) for a test-retest
#' scenario. 
#'
#' @param data Data matrix collecting all individuals' curves evaluated on the same grid of \eqn{M} points. 
#' It is a matrix of dimension \eqn{2n\times M}, where: 
#'     \eqn{n} is the number of individuals, 
#'     \eqn{M} is the number of points of the grid.
#' 
#' @param individuals Factor with \eqn{2n} elements coding the individuals 
#' (the number of levels correspond to the number of individuals). Note that 
#' to compute CMC you must have \eqn{n} individuals, and two curves for each individual (corresponding to test and retest).
#' 
#' @param test Factor with \eqn{2n} elements coding the test sessions 
#' (the number of levels is 2 corresponding to  test and retest). Note that 
#' to compute CMC you must have \eqn{2} test sessions for each individual.
#' 
#' 
#' @return A list containing the following components:
#' \item{CMC}{CMC value computed by averaging the values of all individuals.}
#' \item{CMC_subject}{Vector of size \code{n} containing the CMC value for each individual. }
#' 
#' @seealso See also \code{\link{fICC}} for the functional ICC and \code{\link{fSEM}} for functional SEM. 
#'
#' @details 
#' Note that the CMC could be undefined for some individuals in the case of curves with a small ROM. 
#' In such cases, the common CMC value is computed as the average
#' among all individuals whose CMC is defined. 
#' 
#' @examples
#' # Example with simulated data:
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
#' fCMC(data,individual,test)
#' 
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
#' Ford, K. R., Myer, G. D., and Hewett, T. E. (2007). Reliability of landing 3D motion analysis: implications for longitudinal analyses. 
#' \emph{Medicine and Science in Sports and Exercise}, 39, 2021–2028.
#'
#' @export


fCMC = function(data,individuals,test){
  test = factor(test)
  individuals = factor(individuals)
  table_doe = table(individuals,test)
  if(sum(range(table_doe) == c(1,1)) != 2){
    stop('fCMC is only defined in a test-retest setting with no replicates. Number of individuals in test and retest must coincide')
  }
  
  data.all = data.frame(test,individuals,data)
  data.all = data.all[order(individuals),]
  testlev = levels(test)
  individualslev = levels(individuals)
  
  test = data.all[which(data.all$test==testlev[1]),-(1:2)]
  retest = data.all[which(data.all$test==testlev[2]),-(1:2)]
  
  t = dim(test)[2]
  n = dim(test)[1]
  G = 2 # number of trials per individual
  R = numeric(n)
  Y = abind::abind(test,retest,along=3) 
  mST = (test+retest)/2
  mS = rowMeans(mST) 
  for(subj in 1:n){
    den = sum((Y[subj,,] - mS[subj])^2 ) / (t*G-1)
    num = sum( (Y[subj,,] - cbind(mST[subj,],mST[subj,]))^2 ) / den / (t*(G-1)) 
    R[subj] = sqrt(1-num)
  }
  mR = mean(R,na.rm=TRUE) 
  return(list(CMC=mR,CMC_subject=R))
}

