#' Naive computation for Energy Distance 
#'
#' Computes energy distance, and possibly a p-value. 
#' Suitable for multivariate data. Naive approach, quadratic in number
#' of observations.
#' 
#' @param X Matrix (or vector) of observations in first sample.
#' 
#' @param Y Matrix (or vector) of observations in second sample.
#' 
#' @param pval Boolean for whether to compute p-value or not. 
#' 
#' @param numperm Number of permutations. Default is \code{200}.
#'
#' @param seednum Seed number for generating permutations. Default is \code{0}, 
#'                which means seed is set randomly. For values larger than 
#'                \code{0}, results will be reproducible.
#' 
#' @details First checks number of columns (dimension) are equal. 
#'          Suppose matrix \eqn{X} has \eqn{n} rows and \eqn{d} columns, 
#'          and matrix \eqn{Y} has \eqn{m} rows; checks that \eqn{Y} 
#'          has \eqn{d} columns (if not, then throws error). 
#'          Then flattens matrices to vectors (or, if \eqn{d=1}, they are
#'          already vectors.
#'          Then calls C++ method. If the first sample has \eqn{n} 
#'          \eqn{d}-dimensional samples and the second sample has 
#'          \eqn{m} \eqn{d}-dimensional samples, then the algorithm
#'          computes the statistic in \eqn{O((n+m)^2)} time.
#'          
#'  Random seed is set for \code{std::mt19937} and \code{std::shuffle} in C++.
#'
#' @return A list with the following elements:
#'         \describe{
#'             \item{\code{pval}}{The p-value of the test, if it is  
#'                                computed (\code{pval=TRUE}). }
#'             \item{\code{stat}}{The statistic of the test, which
#'                                is always computed. }
#'          }
#'
#' @references
#'    Baringhaus L. and Franz C. (2004) "On a new multivariate two-sample test." 
#'    Journal of multivariate analysis 88(1):190-206
#' 
#'    Szekely G. J. and Rizzo M. L. (2004) "Testing for equal distributions in 
#'    high dimension." InterStat 5(16.10):1249-1272
#'
#' @examples
#'
#' X <- matrix(c(1:12), ncol=2, byrow=TRUE)
#' Y <- matrix(c(13:20), ncol=2, byrow=TRUE)
#' energydistList <- energydist(X=X, Y=Y, pval=FALSE)
#'
#' #computing p-value
#' energydistList <- energydist(X=X, Y=Y)
#' \donttest{
#' #computing p-value
#' #using 1000 permutations and seed 1 for reproducibility.
#' energydistList <- energydist(X=X, Y=Y, numperm=1000, seednum=1)
#' }
#'
#' @export
energydist <- function(X, Y, pval=TRUE, numperm=200, seednum=0){

    # check vectors/matrices are numeric
    if ( !(is.numeric(X)) ){
        stop("X needs to be numeric.")
    }
    if ( !(is.numeric(Y)) ){
        stop("Y needs to be numeric.")
    }

    # initialise, will update later
    nX <- 0
    dX <- 0
    nY <- 0
    dY <- 0
    Xvec <- c()
    Yvec <- c()

    # add checks for matrix/vectors
    if (is.matrix(X)){
        # get dimensions of matrices
        nX <- dim(X)[1]
        dX <- dim(X)[2]

        if (!(is.matrix(Y))){
            stop("If X is a matrix, Y must also be a matrix.")
        }

        nY <- dim(Y)[1]
        dY <- dim(Y)[2]

        # check dimensions compatible
        if (dX != dY){
            stop("Dimension (number of columns) of matrices need to be equal.")
        }

        # flatten to vectors; we use transpose
        # > X
        #      [,1] [,2] [,3]
        # [1,]    1    3    5
        # [2,]    2    4    6
        #
        # will be flattened to 
        # [1] 1 3 5 2 4 6
        # which is what we want

        Xvec <- as.vector(t(X))
        Yvec <- as.vector(t(Y))
    } else if (is.vector(X)){
        if (!(is.vector(Y))){
            stop("If X is a vector, Y must also be a vector.")
        }
        nX <- length(X)
        dX <- 1
        nY <- length(Y)
        dY <- 1
        Xvec <- X
        Yvec <- Y
    } else {
        stop("X must be a vector or a matrix.")
    }


    # if beta not positive, compute median heuristic (also from C++)
    # finally, compute MMD 
    energydistList <- list()
    if (pval){
        energyList <- energydist_pval_Rcpp(Xvec, Yvec, nX, dX, nY, dY, 
                                           numperm, seednum)
    } else {
        energyList <- energydist_Rcpp(Xvec, Yvec, nX, dX, nY, dY)
    }

    # check pval; if no pval, functions return -1, then changes to NA
    if (energyList$pval < 0){
        energyList$pval <- NA
    }
    return(energyList)
}

