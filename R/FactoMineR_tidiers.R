#' tidying methods for objects produced by \pkg{FactoMineR}
#'
#' These methods provide some summaries of objects produced by FactoMineR, such as MCA or PCA. 
#' 
#' @param x an object of class \code{MCA} or \code{PCA}
#' @param n number of dimensions to take into account
#' @param ... extra arguments (not used)
#' 
#' @name FactoMineR_tidiers
#' @examples
#' 
NULL

#' @method tidy MCA
#' @rdname FactoMineR_tidiers
#' @export

tidy.MCA <- function(x, n = 5, ...) {
    if (n > x$call$ncp) {
        warning("n larger than the numbers of dimensions computed. Reset to the  numbers of dimensions computed.")
        n <- x$call$ncp
    }
    unrowname(data.frame(
                term = rep(row.names(x$var$coord),  n),
                dimension = rep(1:n, length(row.names(x$var$coord))),
                coord = as.vector(x$var$coord[,1:n]),
                contrib = as.vector(x$var$contrib[,1:n]),
                cos2 = as.vector(x$var$cos2[,1:n]),
                v.test = as.vector(x$var$v.test[,1:n]),
                eta2 = as.vector(x$var$eta2[,1:n]),
                stringsAsFactors = FALSE
    ))

}

#' @method augment MCA
#' @rdname FactoMineR_tidiers
#' 
#' @param data original data this was fitted on
#' @param n number of dimensions to take into account
#' 
#' @template augment_NAs
#' 
#' @return \code{augment} returns one row for each original observation,
#' with columns (each prepended by a .) added. They are called .Dim followed by the numer of the dimension.
#' 
#' @export
augment.MCA <- function(x, data = x$call$X, n = 5, ...) {
    if (n > x$call$ncp) {
        warning("n larger than the numbers of dimensions computed. Reset to the  numbers of dimensions computed.")
        n <- x$call$ncp
    }
    
    if (!is.null(row.names(data))) {
        data$.rownames <- row.names(data)
    }

    data[, paste0(".", "Dim", 1:n)] <- x$ind$coord[, 1:n]
    
    return(data)
}

#' @method glance MCA
#' @rdname FactoMineR_tidiers
#' 
#' @param n number of dimensions to retain
#' @param ... extra arguments (not used)
#' 
#' @return \code{glance} returns one row with each column corresponding to the eigenvalue of the dimension
#' 
#' @export
glance.MCA <- function(x, n = 5) {
    if (n > x$call$ncp) {
        warning("n larger than the numbers of dimensions computed. Reset to the  numbers of dimensions computed.")
        n <- x$call$ncp
    }
    
    unrowname(tidyr::spread(data.frame(dim = row.names(x$eig)[1:n],
               value = x$eig$eigenvalue[1:n]), key = dim, value = value))
}

