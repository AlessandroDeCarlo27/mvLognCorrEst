#' @title Computation of Optimization functional
#' @description This function implements functional that must be minimized in order to
#' estimate indirect correlations. Given a starting correlation matrix with the indirect correlations to be estimated
#' encoded by NA, it automatically generates a matrix which depends by a vector, x, of unknown indirect correlations.
#' During the optimization task, the negative of the sum of the minimum eigenvaluein this matrix is minimized in order
#' to obtain indirect correlation estimates such that the correlation matrix is positive definite.
#' @param x vector of unknown variables (indirect correlations)
#' @param cbase Matrix object. It represents the input correlation matrix. Indirect correlations to be estimated must be
#' encoded by \code{NA} in this matrix.
#' @param var_optim Matrix object. The number of rows is equal to the number of unknown indirect correlations while the
#' columns are 4 and contain:
#' \itemize{
#'    \item{\code{var1}:} {numerical index of the first variable of the indirect correlation couple}
#'    \item{\code{var2}:} {numerical index of the second variable of the indirect correlation couple}
#'    \item{\code{lower}:} {lower bound for the range of indirect correlation between \code{var1} and \code{var2}}
#'    \item{\code{upper}:} {upper bound for the range of indirect correlation between \code{var1} and \code{var2}}
#'    }
#' @return The negative of the minimum eigenvalue of correlation matrix. This represents the quantity to minimize in
#' the optimization problem for estimating indirect correlations.
#' @author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}
#' @noRd


optim.fun <- function(x,cbase,var_optim){
    cbase[lower.tri(cbase,diag = F)] <- 0
    corr_eval <- cbase
    counter <- 1
    for (i in 1:nrow(var_optim)) {
        if(is.na(var_optim[i,3])||is.na(var_optim[i,4])){
            str2eval <- paste("corr_eval[",var_optim[i,1],",",var_optim[i,2],"]<- 0",sep="")
            eval(parse(text=str2eval))
        }else{
            str2eval <- paste("corr_eval[",var_optim[i,1],",",var_optim[i,2],"]<- x[",counter,"]",sep="")
            eval(parse(text=str2eval))
            counter <- counter+1
        }

    }
    corr_eval <- corr_eval+t(corr_eval)-diag(rep(1,ncol(corr_eval)))
    return(-min(eigen((corr_eval))$values))
}
