#' @title Functional Pearson Correlation Coefficient
#'
#' @description Computation of the functional Pearson correlation coefficient for a test-retest
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
#' \item{rho}{functional correlation coefficient value integrated over the data domain.}
#' \item{pointwise_rho}{Pointwise curve of the correlation coefficient (evaluated on the same grid of \eqn{M} points provided by the user). }
#' 
#' @seealso See also \code{\link{fICC}} for the functional ICC and \code{\link{fSEM}} for functional SEM. 
#'
#' @details 
#' Note that the CMC could be undefined for some individuals in the case of curves with a small ROM. 
#' In such cases, the common CMC value is computed as the average
#' among all individuals whose CMC is defined. 
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

fRho = function(data,individuals,test){
  test = factor(test)
  individuals = factor(individuals)
  table_doe = table(individuals,test)
  if(sum(range(table_doe) == c(1,1)) != 2){
    stop('fRho is only defined in a test-retest setting with no replicates. Number of individuals in test and retest must coincide')
  }
  
  data.all = data.frame(test,individuals,data)
  data.all = data.all[order(individuals),]
  testlev = levels(test)
  individualslev = levels(individuals)
  
  test = data.all[which(data.all$test==testlev[1]),-(1:2)]
  retest = data.all[which(data.all$test==testlev[2]),-(1:2)]
  
  t = dim(test)[2] 
  n = dim(test)[1]
  if(dim(retest)[2] != t){
    stop("Number of evaluation points in test and retest must coincide") 
  }
  rho_t = (diag(cor(test,retest))) 
  rho = mean(rho_t)
  return(list(rho=rho, pointwise_rho=rho_t)) 
}
