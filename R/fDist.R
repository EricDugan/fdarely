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
