#' @title Functional L2 Distance Between Test and Retest Curves
#'
#' @description Computation of the functional L2 distance for a test-retest
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
#' @param domain Vector of dimension \eqn{M} containing the grid points on the x axis where functional data are evaluated.
#' The grid must be equally spaced. 
#' 
#' @return A list containing the following components:
#' \item{dict}{Distance value computed by averaging the values of all individuals.}
#' \item{dist_subject}{Vector of size \code{n} containing the distance values for each individual. }
#' 
#' @seealso See also \code{\link{fICC}} for the functional ICC and \code{\link{fSEM}} for functional SEM. 
#'
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
#' # Computation of similarity
#' fDist(data,individual,test,domain)
#' 
#' 
#' @references
#' Pini, A., Markström, J., and Schelin, L. (2019): Test–retest reliability measures for curve data: 
#' an overview with recommendations and supplementary code, \emph{Sports Biomechanics}.
#' 
#'
#' @export
#' 

fDist = function(data,individuals,test,domain){
  test = factor(test)
  individuals = factor(individuals)
  table_doe = table(individuals,test)
  if(sum(range(table_doe) == c(1,1)) != 2){
    stop('fDist is only defined in a test-retest setting with no replicates. Number of individuals in test and retest must coincide')
  }
  
  data.all = data.frame(test,individuals,data)
  data.all = data.all[order(individuals),]
  testlev = levels(test)
  individualslev = levels(individuals)
  
  test = data.all[which(data.all$test==testlev[1]),-(1:2)]
  retest = data.all[which(data.all$test==testlev[2]),-(1:2)]
  
  t = dim(test)[2] 
  n = dim(test)[1]
  
  dt = domain[2] - domain[1]
  if(max(abs(diff(domain) - dt)) > (domain[length(domain)] - domain[1])/t*0.00001 ){
    stop("Evaluation points in domain must be equally spaced") }
  dist_subject = sqrt(dt * rowSums( (test - retest)^2)) 
  dist = mean(dist_subject) 
  return(list(dist=dist,dist_subject=dist_subject))
}
