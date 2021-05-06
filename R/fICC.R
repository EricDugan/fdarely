#' @title Functional Intraclass Correlation Coefficient
#'
#' @description Computation of the functional intraclass correlation coefficient (fICC) for a multiple test session
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
#' 
#' @return A list containing the following components:
#' \item{ICC}{ICC value integrated over the data domain.}
#' \item{pointwise_ICC}{Pointwise ICC curve (evaluated on the same grid of \eqn{M} points provided by the user).}
#' \item{low_ICC}{Lower bound of the pointwise confidence interval for the ICC curve evaluated at level \code{alpha}.}
#' \item{up_ICC}{Upper bound of the pointwise confidence interval for the ICC curve evaluated at level \code{alpha}.}
#' 
#' @seealso See also \code{\link{fSEM}} for the functional SEM.
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


fICC = function(data,individuals,test,alpha=0.05){
  t = dim(data)[2]
  N = dim(data)[1]
  n = length(levels(factor(individuals)))
  nik = table(individuals, test)
  Ni = rowSums(nik)
  Nk = colSums(nik)
  K = length(Nk)
  
  pointwise_ICC = low_ICC = up_ICC = numeric(t)
  ks = (N - sum(Ni^2)/N) / (n-1)
  kt = (N - sum(Nk^2)/N) / (K-1)
  
  for(point in 1:t){
    data.all = data.frame(response = data[,point], individuals = factor(individuals), test=factor(test))
    anova = aov(response ~ individuals + test, data = data.all)
    MS = summary(anova)[[1]]$'Mean Sq'
    MSs = MS[1]
    MSt = MS[2]
    MSe = MS[3]
    pointwise_ICC[point] = (kt*(MSs-MSe)) / (ks*MSt + kt*MSs + (ks*kt-ks-kt)*MSe) #(MS[1]-MS[3])/(MS[1]+MS[3]+2/n*(MS[2]-MS[3]))
    a = ks*pointwise_ICC[point]/(kt*(1-pointwise_ICC[point]))
    b = 1 + ks*pointwise_ICC[point]*(kt-1)/(kt*(1-pointwise_ICC[point]))
    nu = (a*MSt+b*MSe)^2 / ((a*MSt)^2/(K-1) + (b*MSe)^2/(N-n-K+1 ))
    Fstar1 = qf(1-alpha/2,n-1,nu)
    Fstar2 = qf(1-alpha/2,nu,n-1)
    low_ICC[point] = kt*(MSs-Fstar1*MSe)/(kt*MSs + Fstar1*(K*MSt+(kt*ks-kt-ks)*MSe) )
    up_ICC[point] = kt*(Fstar2*MSs-MSe)/((ks*MSt+(kt*ks-kt-ks)*MSe) + Fstar2*kt*MSs)
  }
  
  return(list(ICC=mean(pointwise_ICC), pointwise_ICC=pointwise_ICC,low_ICC=low_ICC,up_ICC=up_ICC))
}
