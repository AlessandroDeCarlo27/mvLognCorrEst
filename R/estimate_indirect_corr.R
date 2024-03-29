#' @title Estimation of Indirect Correlations
#' @description This function estimates indirect correlations starting from the incomplete
#' correlation matrix in input. Indirect correlations to be estimated must be indicated
#' by \code{NA} values in the input correlation matrix.
#' @export
#' @importFrom pracma fmincon
#' @importFrom Matrix nearPD
#' @param corrMatStart Matrix object which represents a correlation matrix. Indirect correlations to be
#' estimated must be indicated by \code{NA}. This matrix must be symmetric, thus it must contain at
#' least two \code{NA} values.
#' @param force_estimate Boolean flag. When this flag is \emph{TRUE}, if the obtained correlation matrix is not positive
#' definite, it is approximated to the nearest positive definite matrix based on the Frobenius norm. Matrix
#' approximation may alter fixed initial correlations. If this flag is set to \emph{FALSE} (default option), matrix approximation
#' is skipped and a warning message is returned.
#' @param widen_factor number between 0 and 1. If there is a unique path, the range for that indirect correlation is computed considering \eqn{cost +/- widen_factor*cost}, where \eqn{cost} is the cost of the unique existing path. Default value is 0.2.
#' @details Indirect correlations are estimated solving a constrained optimization problem. Starting
#' from the fixed correlations, a correlation graph is built. Then, for each couple of variables whose indirect
#' correlation is unknown (i.e. \code{NA} values), all the possible paths among them are considered (without
#' visiting a node more than once). The cost of each path is computed by multiplying the correlations along it.
#' The maximum and the minimum costs provide a reasonable range for the indirect correlation value.
#' If there is not any path between two nodes, that indirect correlation will not be estimated and it will be
#' automatically set to 0. If there is a unique path, the range for that indirect correlation is computed considering \eqn{cost +/- widen_factor*cost}, where \eqn{cost} is the cost of the unique existing path. The default value of \eqn{widen_factor} is 0.2.
#'
#' Given the bounds of indirect correlations, a constrained optimization problem is solved by minimizing the negative
#' of minimum eigenvalue of correlation matrix. The starting values for the indirect correlation values are set equal to
#' the middle-point of the computed bounds. If the estimated matrix is not positive
#' definite, user can force a second optimization step in which the previously obtained  matrix
#' is approximated to its nearest positive definite matrix. Frobenius norm is used to measure
#' distance between matrices. Note that this step may alter initial fixed correlations.
#'
#' An indirect correlation between two variables can be estimated only if they are linked by at least one path in
#' the correlation graph. If for all indirect correlations declared does not exist any path, this function
#' prints a warning message and plots the correlation graph to support the debug.
#' @return A list containing \tabular{ll}{
#'   \code{corrMatFinal} \tab Matrix object containing the final correlation matrix with indirect correlations estimated \cr
#'    \tab \cr
#'    \code{optim}\tab List of objects containing the outputs provided by the solver (\code{fmincon} of \code{pracma}
#'    package) used for the constrained optimization. It is returned if the optimization step converges
#'    to a positive definite matrix or if the optimization step fails and \code{force_estimate}
#'    is set to \code{FALSE}.
#'    \cr
#'    \tab \cr
#'    \code{optim1} \tab The same of \code{optim}. It is returned only when constrained optimization does not converge to
#'    a positive definite correlation matrix and \code{force_estimate} is set to \code{TRUE}. \cr
#'    \tab \cr
#'    \code{optim2} \tab List of objects containing the outputs provided by the function \code{nearPD} of \code{Matrix}
#'    package used to approximate the matrix obtained by solving the constrained optimization problem with the nearest positive
#'    definite correlation matrix. It is returned only when constrained optimization does not converge to a positive definite correlation matrix and \code{force_estimate} is set to \code{TRUE}.\cr
#'    \tab \cr
#'    \code{optimizationBounds} \tab{A matrix object with N(= number of indirect correlations) rows and 4 columns reporting:
#'          \itemize{
#'                \item{\code{var1}:} {numerical index of \eqn{X1}, the first variable of the indirect correlation couple}
#'                \item{\code{var2}:} {numerical index of \eqn{X2}, the second variable of the indirect correlation couple}
#'                \item{\code{lower}:} {lower bound for the range of indirect correlation between \code{var1} and \code{var2}}
#'                \item{\code{upper}:} {upper bound for the range of indirect correlation between \code{var1} and \code{var2}}.
#'                }
#'      }\cr
#'    }
#'
#' @author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}
#' @seealso \code{\link[Matrix]{nearPD}}
#' @seealso \code{\link[pracma]{fmincon}}
#' @seealso \code{\link{estimate_corr_bounds}}
#'
#' @examples
#' #define initial correlation structure
#' c_start <- diag(rep(1,10))
#' c_start[1,2] <- -0.6
#' c_start[1,3] <- -0.75
#' c_start[2,3] <-0.95
#' c_start[2,4] <- 0.75
#' c_start[2,6] <- -0.6
#' c_start[2,8] <- 0.75
#' c_start[3,4] <- 0.6
#' c_start[3,8] <-0.75
#' c_start[4,7] <- 0.6
#' c_start[4,8]<-0.75
#' c_start[5,7] <- -0.95
#' #symmetric correlation structure
#' c_start <- c_start+t(c_start)-diag(rep(1,ncol(c_start)))
#' #assign NA to indirect correlations to be estimated
#' c_start[c_start==0]<-NA
#' #names of variables
#' colnames(c_start)<- paste(rep("X",10),1:10,sep = "")
#' rownames(c_start) <- paste(rep("X",10),1:10,sep = "")
#' # plot initial correlation matrix
#' plot_graph_corr(c_start,"Graph of Initial Correlation Matrix")
#' r<-estimate_indirect_corr(c_start)
#' #see final output
#' plot_graph_corr(r$corrMatFinal,'Graph of Final Correlation Matrix')
#'
#'







estimate_indirect_corr <- function(corrMatStart,force_estimate=FALSE,widen_factor=0.2){

    #VALIDATE INPUT
    validate_corrMatrix(corrMatStart)

    #indirect correlations to estimate must be declared with NA inside the input matrix
    if(!any(is.na(corrMatStart))){
        stop("No indirect correlations to estimate are declared")
    }
    #bounds of correlation matrix
    bounds <- estimate_corr_bounds(corrMatStart,widen_factor)

    #test if all bounds are NA: this implies that all indirect correlations estimates will be set to 0
    if(all(is.na(bounds[,3]))&&all(is.na(bounds[,4]))){
        list_variables <- paste(paste(bounds[,1],bounds[,2],sep="--"),collapse = " ")
        warning(paste("Cannot estimate indirect correlations, paths between all indicated variables are missing:\n",list_variables,sep=""))
        plot_graph_corr(corrMatStart,'Independet Variables')
        matOpt <- get_corrMatrixOptim(NULL,corrMatStart,bounds)
        output_optim <- list()
        output_optim[["corrMatFinal"]]<-matOpt
        return(output_optim)
    }

    #set to 0 correlations to be estimated
    corrMatStart[is.na(corrMatStart)]<-0
    # get indices of couples for which an indirect effect exists according to
    #graph path analysis
    notNa_idx <- !(is.na(bounds[,3])&is.na(bounds[,4]))
    x0 <- bounds[notNa_idx,3]+(sign(bounds[notNa_idx,3])*
                                   ((abs(bounds[notNa_idx,4])-abs(bounds[notNa_idx,3]))/2))
    #solve constrained optimization problem
    r <- pracma::fmincon(x0,optim.fun,lb=bounds[notNa_idx,3],
                         ub=bounds[notNa_idx,4],cbase=corrMatStart,
                         var_optim=bounds)

    #compute correlation matrix with estimated indirect correlations
    corrMatfinal <- get_corrMatrixOptim(r$par,corrMatStart,bounds)
    output_optim <- list()
    output_optim[["optimizationBounds"]] <-bounds
    if(r$val<0){
        #optimization successfull
        output_optim[["optim"]]<-r
        output_optim[["corrMatFinal"]] <- corrMatfinal
        return(output_optim)
    }else{
        if(force_estimate){
            warning("Optimization step-1 failed. Matrix obtained is not positive-semidefinite.
                Matrix of step-1 will be approximated to its nearest positive definite correlation matrix.
                This may change fixed correlations.")
            output_optim[["optim1"]] <-r
            r2 <- Matrix::nearPD(corrMatfinal,corr = TRUE)
            output_optim[["optim2"]] <-r2
            output_optim[["corrMatFinal"]] <- as.matrix(r2$mat)
        }else{
            warning("Optimization step-1 failed. Matrix obtained is not positive-semidefinite.
                To approximate Matrix of step-1 to its nearest positive definite correlation matrix set force_estimate=TRUE")
            output_optim[["optim"]]<-r
            output_optim[["corrMatFinal"]] <- corrMatfinal
        }
    }

    return(output_optim)

}
