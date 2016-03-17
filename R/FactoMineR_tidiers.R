#' tidying methods for objects produced by \pkg{FactoMineR}
#'
#' These methods provide some summaries of objects produced by FactoMineR, such as MCA or PCA. 
#' 
#' @param x an object of class \code{MCA} or \code{PCA}
#' @param var.sup supplementary variables. A data frame with the same number of rows as the original dataframes, and discrete variables
#' @param n number of dimensions to take into account
#' @param ... extra arguments (not used)
#' 
#' @name FactoMineR_tidiers
#' @examples
#' 
NULL

#' @method tidy MCA
#' @rdname FactoMineR_tidiers
#' @import dplyr tidyr
#' @export

tidy.MCA <- function(x, var.sup = NULL, n = 5, ...) {
    if (n > x$call$ncp) {
        warning("n larger than the numbers of dimensions computed. Reset to the  numbers of dimensions computed.")
        n <- x$call$ncp
    }
    tmp <- data.frame(
                term = rep(row.names(x$var$coord), n),
                dimension = rep(1:n, each = length(row.names(x$var$coord))),
                coord = as.vector(x$var$coord[,1:n]),
                contrib = as.vector(x$var$contrib[,1:n]),
                cos2 = as.vector(x$var$cos2[,1:n]),
                v.test = as.vector(x$var$v.test[,1:n]),
                eta2 = NA,
                stringsAsFactors = FALSE
    )
    tmp <- rbind(tmp, data.frame(
                term = rep(row.names(x$var$eta2), n),
                dimension = rep(1:n, each = length(row.names(x$var$eta2))),
                coord = NA,
                contrib = NA,
                cos2 = NA,
                v.test = NA,
                eta2 = as.vector(x$var$eta2[,1:n]),
                stringsAsFactors = FALSE)
    )
    
    if (!is.null(var.sup)) {
        if (nrow(as.data.frame(var.sup)) == nrow(x$call$X)) {
            varsups <- names(as.data.frame(var.sup))
            df <- cbind(as.data.frame(x$ind$coord), var.sup)
            tmp <- rbind(
                tmp,
                df %>% 
                    mutate(w = x$call$row.w) %>% 
                    gather_("variable", "term", varsups) %>% 
                    group_by(term) %>% 
                    summarise_each(funs(weighted.mean(., w = w)), 1:n) %>% 
                    gather(dimension, coord, -term) %>% 
                    mutate(dimension = as.integer(stringr::str_replace(dimension, "Dim ", ""))) %>% 
                    mutate(contrib = NA, cos2 = NA, v.test = NA, eta2 = NA)
            )

        }
    }
    
    unrowname(tmp)

}

#' @method augment MCA
#' @rdname FactoMineR_tidiers
#' 
#' @param data original data this was fitted on
#' @param ind.sup a dataframe with supplementary individuals. It must have the same (active) variables as the original dataframe
#' @param n number of dimensions to take into account
#' 
#' @template augment_NAs
#' 
#' @return \code{augment} returns one row for each original observation, and for the supplemnetary individuals if there are some,
#' with columns (each prepended by a .) added. They are called .Dim followed by the numer of the dimension.
#' 
#' @export
augment.MCA <- function(x, data = x$call$X, ind.sup = NULL, n = 5, ...) {
    if (n > x$call$ncp) {
        warning("n larger than the numbers of dimensions computed. Reset to the  numbers of dimensions computed.")
        n <- x$call$ncp
    }
    
    if (!is.null(row.names(data))) {
        data$.rownames <- row.names(data)
    }

    data[, paste0(".", "Dim", 1:n)] <- x$ind$coord[, 1:n]
    
    if (!is.null(ind.sup)) {
        if (identical(names(ind.sup), names(x$call$X[,x$call$quali]))) {
            data <- rbind(data, 
                          cbind(ind.sup, 
                                tab.disjonctif(ind.sup) %*% x$var$coord[,1:ncp]
                                )
                          )
        } else {
            warning("ind.sup has different variables names than the original data, not taken into account")
        }
    }
    
    
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

