#' @title Computation of the maximum and the minimum costs of paths between two
#' nodes in a correlation graph
#' @description Given a correlation graph, a starting node, \eqn{X1} and a target node, \eqn{X2},
#' this function receives in input the correlation matrix associated to the considered graph and the list of all possible paths which link \eqn{X1} and \eqn{X2}
#' visiting each node at most once. Each path is represented by a vector which starts with \eqn{X1}, ends with \eqn{X2} and contains in
#' the middle the sequence of nodes to visit before reaching the target. For each path, the relative
#' cost is computed by multiplying the correlations met along the path. Then, their minimum and maximum values are returned in a 2-d
#' array. If there is not any path between two nodes, maximum and minimum values will be set to \code{NA}. If
#' there is a unique path, the bound is computed considering \eqn{cost +/- 0.5*cost}.
#' @param corrmat Matrix Object. It represents the correlation structure among variables
#' @param list_path List Object. It contains all the possible paths between two nodes of the graph.
#' @return Array object with minimum and maximum costs of the possible paths between a couple
#' of nodes
#' @note
#' In order to estimate indirect correlations, this function is used to compute thei upper and
#' lower bounds.
#' @author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}
#' @noRd

get_bounds <- function(corrmat,list_path){
    #if a path that links two variables does not exist the bounds are set to NA
    #in the optimization step it will be set to 0
    if (length(list_path)==0) {
        return(array(c(NA,NA)))
    }
    # if only a path is available, the range is built around the only value available
    if (length(list_path)==1) {
        single_value <- get_pathProduct(corrmat,as.vector(list_path[[1]]))
        if(single_value<0){
            max_r <- single_value+(abs(single_value)/5)
            min_r <- max(single_value-(abs(single_value)/5),-1) #prevent lower bounds < -1
        }else{
            max_r <- min(single_value+(single_value/5),1) #prevent upper bounds >1
            min_r <- single_value-(single_value/5)
        }
        return(array(c(min_r,max_r)))
    }

    path_values <- array(0L,dim=length(list_path))
    for (i in 1:length(list_path)) {
        path_values[i] <- get_pathProduct(corrmat,as.vector(list_path[[i]]))
    }

    return(array(c(min(path_values),max(path_values))))


}
