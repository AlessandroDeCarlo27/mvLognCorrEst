#' @title Convert Log-Normal parameters to Normal parameters
#' @description This function is used to convert mean vector and covariance matrix of a multivariate Log-Normal
#' distribution to mean vector and covariance matrix of the associated multivariate Normal distribution.
#' Two conditions are evaluated before computing the parameters of Normal distribution:
#' \enumerate{
#' \item For each couple of variables, \eqn{X1}, \eqn{X2}, it is tested if \eqn{\rho(X1,X2)} satisfies the condition
#' necessary to obtain a positive definite correlation/covariance matrix for the
#' Normal distribution underlying the Log-Normal one. This condition is also sufficient for 2x2 matrix [1].
#' If this condition is not satisfied, a warning message is displayed
#' but covariance matrix of the Normal distribution is still computed even if it is not positive definite.
#' \item For each couple of variables, \eqn{X1}, \eqn{X2}, it is tested if \eqn{\rho(X1,X2)*cv1*cv2 > -1}.
#' If this condition is not satisfied for all the couples the covariance matrix of the underlying Normal distribution cannot
#' be computed. In fact, the covariance between \eqn{X1}, \eqn{X2} of the Normal distribution, \eqn{\Sigma^2(X1,X2)},
#' is defined as \eqn{ln((\rho(X1,X2)*cv1*cv2 )+1)} [1]. In this case, a warning message will be displayed and a summary of
#' the performed tests is returned.
#' }
#'
#' @param mu Array object containing means of multivariate Log-Normal distribution.
#' @param covMatrix Matrix object containing covariance matrix of multivariate Log-Normal distribution.
#' @export
#' @return A list containing:\tabular{ll}{
#'    \code{muN} \tab Array object. It contains mean of the multivariate Normal distribution associated to
#'                    Log-Normal one.\cr
#'    \tab \cr
#'    \code{sigmaN} \tab Matrix object. It contains covariance matrix of the multivariate Normal distribution associated
#'                   to Log-Normal one.\cr
#'  \tab \cr
#'   \code{is_valid_logN_corrMat} \tab Boolean flag which indicates if the input correlation matrix satisfies
#'    condition 1.\cr
#'    \tab \cr
#'    \code{can_log_transf_covMat} \tab Boolean flag which indicates if the input correlation matrix satisfies
#'    condition 2. \cr
#'    \tab \cr
#'    \code{validation_res} \tab Matrix object which contains a brief summary of the tested conditions on the input
#'    correlation matrix. It is composed by N(=number of Log-Normal variable couples) rows and 8 columns:\cr
#'    \tab \itemize{
#'    \item{\code{var1}:} {numerical index of the first Log-Normal variable of the couple, \eqn{X1}}
#'    \item{\code{var2}:} {numerical index of the second Log-Normal variable of the couple, \eqn{X2}}
#'    \item{\code{lower}:} {lower bound for the range of Log-Normal correlation between \code{var1} and \code{var2}}
#'    \item{\code{upper}:} {upper bound for the range of Log-Normal correlation between \code{var1} and \code{var2}}
#'    \item{\code{tested_corr}:} {input correlation value between \code{var1} and \code{var2}}
#'    \item{\code{is_valid_bound}:} {numerical flag. If \emph{1} \code{tested_corr} is inside the range, \emph{0}
#'    otherwise}
#'    \item{\code{tested_for_logTransf}:} {value computed for testing condition 2 (\eqn{\rho(X1,X2)*cv1*cv2})}
#'    \item{\code{can_logTransf}:} {numerical flag. If \emph{1} covariance of normal distribution between
#'    \code{var1} and \code{var2} can be computed, \emph{0} otherwise.}
#'    } \cr
#' }
#' @note This function tries to compute the covariance matrix of the multivariate Normal distribution even if the
#' Log-Normal covariance matrix is not positive definite. User should check that input matrix is positive definite, however in this
#' situation \code{logn_to_normal} returns a warning message. Even if input matrix is positive definite
#' and conditions 1 and 2 are fulfilled, covariance matrix of multivariate Normal distribution could not be positive definite
#' because condition 1 is only necessary [1]. In order to fix
#' this problem, \code{nearPD} of \code{Matrix} package could be used.
#' @author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}
#' @references
#' [1] Henrique S. Xavier, Filipe B. Abdalla, Benjamin Joachimi, Improving lognormal models for
#' cosmological fields, Monthly Notices of the Royal Astronomical Society,
#' Volume 459, Issue 4, 11 July 2016, Pages 3693–3710,
#' \href{https://doi.org/10.1093/mnras/stw874}{https://doi.org/10.1093/mnras/stw874}
#' @seealso \code{\link[Matrix]{nearPD}}
#' @seealso \code{\link{get_logNcorr_bounds}}
#' @seealso \code{\link{validate_logN_corrMatrix}}
#'
#' @examples
#' #define correlations
#' corr<- diag(rep(1,4))
#' corr[1,4] <- 0.9
#' corr[4,1]<-corr[1,4]
#' corr[2,4] <- -0.3
#' corr[4,2] <- corr[2,4]
#' corr[3,2] <- -0.2
#' corr[2,3] <- corr[3,2]
#'
#' #define sd of variables
#' sd2 <- array(c(rep(1,4)))
#' #obtain covariance matrix
#' covMatrix2 <- sd2%*%t(sd2)*corr
#' #define mean vector
#' mu2 <- array(rep(2.5,4))
#' logn_to_normal(mu2,covMatrix2)
#'
#' #output:
#'
#' # $is_valid_logN_corrMat
#' # [1] TRUE
#'
#' # $can_log_transf_covMat
#' # [1] TRUE
#' #
#' # $validation_res
#' #       var1 var2     lower upper tested_corr is_valid_bound tested_for_logTransf can_logTransf
#' # [1,]    1    2 -0.862069     1         0.0              1                0.000             1
#' # [2,]    1    3 -0.862069     1         0.0              1                0.000             1
#' # [3,]    1    4 -0.862069     1         0.9              1                0.144             1
#' # [4,]    2    3 -0.862069     1        -0.2              1               -0.032             1
#' # [5,]    2    4 -0.862069     1        -0.3              1               -0.048             1
#' # [6,]    3    4 -0.862069     1         0.0              1                0.000             1
#'
#' # $muN
#' # [1] 0.8420807 0.8420807 0.8420807 0.8420807
#'
#' # $sigmaN
#' # [,1]        [,2]        [,3]        [,4]
#' # [1,] 0.1484200  0.00000000  0.00000000  0.13453089
#' # [2,] 0.0000000  0.14842001 -0.03252319 -0.04919024
#' # [3,] 0.0000000 -0.03252319  0.14842001  0.00000000
#' # [4,] 0.1345309 -0.04919024  0.00000000  0.14842001


logn_to_normal <- function(mu,covMatrix){
    #check input
    validate_covMatrix(covMatrix)
    validate_array("mu",mu,covMatrix,"covMatrix")

    if(min(eigen(covMatrix)$values)<0){
        warning("Log-Normal covariance matrix to be converted in Normal covariance matrix is not positive definite")
    }

    #check values of correlation matrix
    sdLn <- as.array(sqrt(diag(covMatrix)))
    corrMatrix <- covMatrix/(sdLn%*%t(sdLn))
    valid_res <- validate_logN_corrMatrix(mu,sdLn,corrMatrix)

    muN <- array(c(rep(0, length(mu))))
    names(muN) <- names(mu)
    sigmaN <- matrix(0L, ncol = length(mu), nrow = length(mu))
    colnames(sigmaN) <- colnames(covMatrix)
    row.names(sigmaN) <- row.names(sigmaN)

    if(!valid_res$can_log_transf_covMat){
        warning("Cannot compute covariance matrix of associated normal distribution. Exists at least a couple of variables such that constrained for logarithmic transformation is not valid. Check validation_res variable in output list for further details and the documentation. ")
        return(valid_res)

    }else{

        if(!valid_res$is_valid_logN_corrMat){
            warning("Correlations set for Log-Normal distributions do not respect theorical bounds. This condition it necessary (but not sufficient) for obtaining a positive definite covariance matrix for multivariate Normal distribution.")
        }


        #get covariace matrix of normal distribution associated to starting logN
        for (i in 1:length(mu)) {
            for (j in 1:length(mu)) {
                sigmaN[i, j] <- log((covMatrix[i, j] / (mu[i] * mu[j])) + 1)
            }
        }

        sd_logN <- diag(sigmaN)
        #get mean vector of normal distribution associated to starting logN
        muN <- log(mu) - 0.5 * sd_logN

        out_list <- list()
        out_list[["muN"]] <- muN
        out_list[["sigmaN"]] <- sigmaN
        out_list <- append(valid_res,out_list)
        return(out_list)
    }


}
