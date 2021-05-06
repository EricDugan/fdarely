#' @title Functional Similarity Between Test and Retest Curves
#'
#' @description Computation of the functional similarity for a test-retest
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
#' \item{sim}{Similarity value computed by averaging the values of all individuals.}
#' \item{sim_subject}{Vector of size \code{n} containing the similarity value for each individual. }
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

fSim = function(data,individuals,test){
  test = factor(test)
  individuals = factor(individuals)
  table_doe = table(individuals,test)
  if(sum(range(table_doe) == c(1,1)) != 2){
    stop('fSim is only defined in a test-retest setting with no replicates. Number of individuals in test and retest must coincide')
  }
  
  data.all = data.frame(test,individuals,data)
  data.all = data.all[order(individuals),]
  testlev = levels(test)
  individualslev = levels(individuals)
  
  test = data.all[which(data.all$test==testlev[1]),-(1:2)]
  retest = data.all[which(data.all$test==testlev[2]),-(1:2)]
  
  t = dim(test)[2] 
  n = dim(test)[1]
  
  sim_subject = rowSums(test * retest)/sqrt(rowSums(test^2)*rowSums(retest^2)) 
  sim = mean(sim_subject)
  return(list(sim=sim, sim_subject=sim_subject))
}
