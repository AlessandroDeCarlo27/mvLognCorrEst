#' @title Approximation to nearest positive definite for Covariance/Correlation Matrix
#' @description Function that test if a given input matrix object is positive definite or not.
#' If it is so, then input matrix is returned.
#' Otherwise, if \emph{flag_force} is \emph{TRUE}, input matrix is approximated to nearest positive definite, new matrix and
#' normalized Frobenius and Infinity norms are returned. If \emph{flag_force} is \emph{FALSE}, then a warning message is
#' displayed and input matrix is returned as output.
#' @importFrom Matrix nearPD norm
#' @param matr Matrix object representing a covariance/correlation matrix
#' @param name String with the name of matrix object representing a covariance/correlation matrix
#' @param flag_force Boolean flag. If \emph{TRUE}, input matrix is approximated
#' to nearest positive definite if necessary; if \emph{FALSE} approximation to nearest positive definite is not performed, a warning message is
#' displayed and input matrix is returned as output.
#' @param is_corr Boolean flag used to distinguish Correlation to Covariance matrices.
#' @return A list containing:\tabular{ll}{
#'    \code{matr} \tab Matrix object. It contains input correlation/covariance matrix, if the one passed in
#'                     input is positive definite or if \emph{flag_force} is set to \emph{FALSE}. If input matrix object is not
#'                     positive definite and  \emph{flag_force} is set to \emph{FALSE}, \code{matr} contains the nearest positive definite matrix
#'                     to input one.\cr
#'    \tab \cr
#'    \code{normF} \tab Scalar double. It represents the Frobenius Norm of the difference between input matrix and
#'                      its nearest PD, normalized for  the Frobenius Norm of input matrix. If input matrix is PD
#'                      or \emph{flag_force} is set to \emph{FALSE}, this element is not returned.\cr
#'    \tab \cr
#'    \code{normInf} \tab Scalar double. It represents the Infinity norm of the difference between input matrix and
#'                      its nearest PD, normalized for  the Infinity norm of input matrix. If input matrix is PD
#'                      or \emph{flag_force} is set to \emph{FALSE}, this element is not returned.\cr
#' }
#'@author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}
#' @seealso \code{\link[Matrix]{nearPD}}
#' @seealso \code{\link[Matrix]{norm}}
#'
#' @noRd



req_approx_to_PD<-function(matr,name,flag_force,is_corr){

    out_list <- list()

    if(min(eigen(matr)$values)<0){

        if(flag_force){
            warning(paste("Input",name,"is not positive definite. It will be approximated to its nearest positive definite. Set force_toPD=FALSE to avoid",sep=" "))
            new_matPD<- Matrix::nearPD(matr,corr = is_corr)
            out_list[["matr"]] <- as.matrix(new_matPD$mat)
            out_list[["normF"]] <- new_matPD$normF/Matrix::norm(matr,"f")
            out_list[["normInf"]] <- Matrix::norm(out_list[["matr"]]-matr,"i")/
                Matrix::norm(matr,"i")
        }else{
            warning(paste("Input",name,"is not positive definite. Set force_toPD=TRUE to approximate it to its nearest positive definite",sep=" "))
            out_list[["matr"]] <- matr
        }

    }else{
        out_list[["matr"]] <- matr
    }

    return(out_list)
}
