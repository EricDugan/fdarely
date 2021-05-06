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
#' 
#' 
#' @references
#' Pini, A., Markström, J., and Schelin, L. (2019): Test–retest reliability measures for curve data: 
#' an overview with recommendations and supplementary code, \emph{Sports Biomechanics}.
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

