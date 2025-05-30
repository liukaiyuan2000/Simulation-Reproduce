#' MEA-MMD: Multivariate Efficient Approximate Maximum Mean Discrepancy
#'
#' Computes maximum mean discrepancy statistics with Laplacian 
#' or Gaussian kernel. 
#' Suitable for multivariate data. Naive approach, quadratic in number
#' of observations.
#' 
#' @param X Matrix (or vector) of observations in first sample.
#' 
#' @param Y Matrix (or vector) of observations in second sample.
#' 
#' @param beta kernel parameter. Must be positive; if not, computes
#'             median heuristic in quadratic time for each projection. 
#'             Default value
#'             is \code{-0.1}, which will force median heuristic to be used.
#' 
#' @param pval Boolean for whether to compute p-value or not. 
#' 
#' @param type The type of projection used. Either \code{"proj"} for 
#'             random projections (default) or \code{"dist"} for interpoint
#'             distances.
#' 
#' @param numproj Number of projections (only used if \code{type="proj"}).
#'                Default is \code{20}.
#' 
#' @param nmethod Norm used for interpoint distances, if \code{type="dist"}.
#'                Needs to be either \code{2} (for two-norm, default) or 
#'                \code{1} (for one-norm).
#' 
#' @param distpval The p-value combination procedure if \code{type="dist"}.
#'                 Options are \code{"Hommel"} (default) or \code{"Fisher"}.
#'                 The Hommel method is preferred since the Type I error does 
#'                 not seem to be controlled if the Fisher method is used.
#' 
#' @param numperm Number of permutations. Default is \code{200}.
#'
#' @param seednum Seed number for generating permutations. Default is \code{0}, 
#'                which means seed is set randomly. For values larger than 
#'                \code{0}, results will be reproducible.
#' 
#' @return A list with the following elements:
#'         \describe{
#'             \item{\code{pval}}{The p-value of the test, if it is  
#'                                computed (\code{pval=TRUE}). Otherwise, 
#'                                it is set to \code{NA}.}
#'             \item{\code{stat}}{The statistic of the test, which
#'                                is only returned when \code{type="proj"},
#'                                otherwise it is set to \code{NA}.}
#'          }
#'
#' @references
#'    Bodenham, D. A., and Kawahara, Y. (2023)
#'    "euMMD: efficiently computing the MMD two-sample test statistic for 
#'    univariate data." Statistics and Computing 33.5 (2023): 110.
#' 
#' @examples
#' X <- matrix(c(1:12), ncol=2, byrow=TRUE)
#' Y <- matrix(c(13:20), ncol=2, byrow=TRUE)
#' # using the random projections method
#' mmdList <- meammd(X=X, Y=Y, pval=TRUE, type="proj", numproj=50)
#'
#' # using the method were distances are computed to the various points 
#' mmdList <- meammd(X=X, Y=Y, pval=TRUE, type="dist")
#'
#'
#' @export 
meammd <- function(X, Y, beta=-0.1, pval=TRUE, 
                   type=c("proj", "dist"), numproj=20, nmethod=c(2, 1),
                   distpval=c("Hommel", "Fisher"), numperm=200, seednum=0){

    # check vectors/matrices are numeric
    if ( !(is.numeric(X)) && !(is.matrix(X)) ){
        stop("X needs to be a numeric matrix.")
    }
    if ( !(is.numeric(Y)) && !(is.matrix(Y)) ){
        stop("Y needs to be a numeric matrix.")
    }

    # check kernel is correct
    type <- type[1]
    if ( (type != "proj") && (type != "dist") ){
        stop("type needs to be either 'proj' or 'dist'.")
    }

    # for distpval; default is Hommel
    pmethod <- 0
    if (type=="dist") {
        if (pval==FALSE) {
            stop("MEA-MMD-Dist does not return a statistic, only a p-value.")
        }
        distpval <- distpval[1]
        if ( (distpval != "Hommel") && (distpval != "Fisher") ) {
            stop("distpval needs to be either 'Hommel' or 'Fisher'.")
        }

        ## Temporarily block Fisher method
        if (distpval == "Fisher") {
            stop("'Fisher' not supported in this version.")
        }


        if (distpval != "Hommel"){
            pmethod <- 1
            # then will be Fisher
        }

        nmethod <- nmethod[1]
        if ( (nmethod != 1) && (nmethod != 2) ){
            stop("nmethod needs to be either 1 or 2")
        }
        if (nmethod==2){
            # C++ code expects 0 for 2 norm to be used
            nmethod <- 0
        }
    }

    nX <- dim(X)[1]
    dX <- dim(X)[2]
    nY <- dim(Y)[1]
    dY <- dim(Y)[2]

    # check dimensions compatible
    if (dX != dY){
        stop("Dimension (number of columns) of matrices need to be equal.")
    }


    # flatten to vectors;  we use transpose
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

    meammdList <- list()
    if (type=="proj"){
        if (pval){
            meammdList <- meammd_proj_pval_Rcpp(Xvec, 
                                                Yvec,
                                                nX, dX,
                                                nY, dY, 
                                                numperm, 
                                                numproj, 
                                                seednum, 
                                                beta)
        } else {
            stat <- meammd_proj_Rcpp(Xvec, 
                                     Yvec,
                                     nX, dX,
                                     nY, dY, 
                                     numproj, 
                                     seednum,
                                     beta)
            meammdList$stat <- stat
            meammdList$pval <- NA
        }
    } else {
        # type is dist
         pval <- meammd_dist_pval_Rcpp(Xvec, 
                                       Yvec,
                                       nX, dX,
                                       nY, dY, 
                                       numperm, 
                                       seednum, 
                                       beta, 
                                       pmethod, 
                                       nmethod)

            meammdList$stat <- NA
            meammdList$pval <- pval

    }
    # return list
    return(meammdList)
}
