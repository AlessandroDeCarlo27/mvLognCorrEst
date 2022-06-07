#' @title Validation of input covariance matrix
#' @description function used to check if an input matrix represents a valid covariance
#' matrix. If these requisites are not respected, execution will be stopped.
#' @param covMatrix object that represents input covariance matrix to validate.
#' @note Eigenvalues are not considered because the package can handle not positive definite matrices too
#' by approximating them to their nearest positive definite.
#' @author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}
#' @noRd


validate_covMatrix <- function(covMatrix){
    if(!is.matrix(covMatrix)){
        stop("covMatrix must be a matrix object.")
    }

    if(!isSymmetric(covMatrix)){
        stop("covMatrix must be symmetric.")
    }

}
