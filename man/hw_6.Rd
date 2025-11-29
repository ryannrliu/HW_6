#' sparseNumeric: Sparse Numeric Vector Class and Methods
#'
#' Provides an S4 class \code{sparse_numeric} for storing numeric vectors in
#' sparse form and a collection of associated methods for arithmetic,
#' coercion, plotting, summary statistics, and standardization.
#'
#' The implementation performs operations using only the non-zero entries
#' whenever possible to avoid expanding to dense format.
#'
#' @docType package
#' @name sparseNumeric
#' @import methods
NULL

#' Sparse numeric vector S4 class
#'
#' An object of class \code{sparse_numeric} stores only the non-zero values of
#' a numeric vector, the 1-based positions of those values, and the full
#' length of the vector.
#'
#' @slot value numeric vector of non-zero values
#' @slot pos integer vector of 1-based positions for the corresponding values
#' @slot length integer length of the full (dense) vector
#' @exportClass sparse_numeric
#' @examples
#' x <- new("sparse_numeric", value = c(4, 3.2), pos = c(2L, 3L), length = 5L)
#' show(x)
setClass(
  Class = "sparse_numeric",
  slots = c(
    value = "numeric",
    pos = "integer",
    length = "integer"
  )
)

#' Validate a sparse_numeric object
#'
#' Internal validity method for \code{sparse_numeric} objects.
#'
#' @param object A \code{sparse_numeric} object.
setValidity("sparse_numeric", function(object) {
  if (length(object@length) != 1L || object@length < 0L)
    return("slot 'length' must be a single non-negative integer")
  if (!is.integer(object@pos))
    return("slot 'pos' must be integer")
  if (!is.numeric(object@value))
    return("slot 'value' must be numeric")
  if (length(object@value) != length(object@pos))
    return("slots 'value' and 'pos' must have the same length")
  if (length(object@pos) > 0L) {
    if (any(object@pos < 1L) || any(object@pos > object@length))
      return("entries of 'pos' must be within 1 and object@length")
    if (any(duplicated(object@pos)))
      return("entries of 'pos' must be unique")
  }
  TRUE
})

#' Coerce numeric to sparse_numeric
#'
#' Create a \code{sparse_numeric} object from a base numeric vector.
#' Zeros are omitted from the stored representation.
#'
#' @param from A numeric vector.
#' @return An object of class \code{sparse_numeric}.
#' @examples
#' as(c(0, 2, 0, 5), "sparse_numeric")
#' @export
setAs("numeric", "sparse_numeric", function(from) {
  nz <- which(from != 0)
  new("sparse_numeric",
      value = if (length(nz)) as.numeric(from[nz]) else numeric(0),
      pos = if (length(nz)) as.integer(nz) else integer(0),
      length = as.integer(length(from)))
})

#' Coerce sparse_numeric to numeric
#'
#' Expand a \code{sparse_numeric} back to a dense numeric vector.
#'
#' @param from A \code{sparse_numeric} object.
#' @return A numeric vector of length \code{from@length}.
#' @examples
#' s <- as(c(0,1,0,2), "sparse_numeric"); as(s, "numeric")
#' @export
setAs("sparse_numeric", "numeric", function(from) {
  out <- numeric(from@length)
  if (length(from@pos) > 0L) out[from@pos] <- from@value
  out
})

# internal helper: check lengths equal
.check_len <- function(x, y) {
  if (x@length != y@length) stop("vectors must have the same length")
}

#' Sparse addition generic
#'
#' Generic for addition on sparse vectors.
#' @param x,y operands
#' @export
setGeneric("sparse_add", function(x, y, ...) standardGeneric("sparse_add"))

#' Sparse subtraction generic
#' @param x,y operands
#' @export
setGeneric("sparse_sub", function(x, y, ...) standardGeneric("sparse_sub"))

#' Sparse multiplication (element-wise) generic
#' @param x,y operands
#' @export
setGeneric("sparse_mult", function(x, y, ...) standardGeneric("sparse_mult"))

#' Sparse crossproduct (dot product) generic
#' @param x,y operands
#' @export
setGeneric("sparse_crossprod", function(x, y, ...) standardGeneric("sparse_crossprod"))

#' Norm generic
#' @param x operand
#' @export
setGeneric("norm", function(x, ...) standardGeneric("norm"))

#' Standardize generic
#' @param x operand
#' @export
setGeneric("standardize", function(x, ...) standardGeneric("standardize"))

# Operator generics (allow S4 methods for primitive ops)
setGeneric("+", function(e1, e2) standardGeneric("+"))
setGeneric("-", function(e1, e2) standardGeneric("-"))
setGeneric("*", function(e1, e2) standardGeneric("*"))

#' Add two sparse_numeric objects
#'
#' Adds two sparse vectors and returns a \code{sparse_numeric} object.
#' This function operates on the non-zero entries only.
#'
#' @param x,y \code{sparse_numeric} objects
#' @return \code{sparse_numeric}
#' @examples
#' x <- as(c(0,1,0,2), "sparse_numeric")
#' y <- as(c(1,0,3,0), "sparse_numeric")
#' sparse_add(x, y)
#' @export
setMethod("sparse_add", signature(x = "sparse_numeric", y = "sparse_numeric"),
          function(x, y, ...) {
            .check_len(x, y)
            allpos <- sort(unique(c(x@pos, y@pos)))
            if (length(allpos) == 0L) return(new("sparse_numeric", value=numeric(0), pos=integer(0), length=x@length))

            xv <- yv <- numeric(length(allpos))

            mx <- match(allpos, x@pos)
            if (any(!is.na(mx))) xv[!is.na(mx)] <- x@value[mx[!is.na(mx)]]

            my <- match(allpos, y@pos)
            if (any(!is.na(my))) yv[!is.na(my)] <- y@value[my[!is.na(my)]]

            vals <- xv + yv
            keep <- which(vals != 0)
            if (length(keep) == 0L) return(new("sparse_numeric", value = numeric(0), pos = integer(0), length = x@length))

            new("sparse_numeric", value = vals[keep], pos = as.integer(allpos[keep]), length = x@length)
          })

#' Subtract two sparse_numeric objects
#'
#' @inheritParams sparse_add
#' @return \code{sparse_numeric}
#' @export
setMethod("sparse_sub", signature(x = "sparse_numeric", y = "sparse_numeric"),
          function(x, y, ...) {
            .check_len(x, y)
            allpos <- sort(unique(c(x@pos, y@pos)))
            if (length(allpos) == 0L) return(new("sparse_numeric", value=numeric(0), pos=integer(0), length=x@length))

            xv <- yv <- numeric(length(allpos))
            mx <- match(allpos, x@pos)
            if (any(!is.na(mx))) xv[!is.na(mx)] <- x@value[mx[!is.na(mx)]]

            my <- match(allpos, y@pos)
            if (any(!is.na(my))) yv[!is.na(my)] <- y@value[my[!is.na(my)]]

            vals <- xv - yv
            keep <- which(vals != 0)
            if (length(keep) == 0L) return(new("sparse_numeric", value = numeric(0), pos = integer(0), length = x@length))

            new("sparse_numeric", value = vals[keep], pos = as.integer(allpos[keep]), length = x@length)
          })

#' Multiply two sparse_numeric objects (element-wise)
#'
#' @inheritParams sparse_add
#' @return \code{sparse_numeric}
#' @export
setMethod("sparse_mult", signature(x = "sparse_numeric", y = "sparse_numeric"),
          function(x, y, ...) {
            .check_len(x, y)
            inter <- intersect(x@pos, y@pos)
            if (length(inter) == 0L) return(new("sparse_numeric", value = numeric(0), pos = integer(0), length = x@length))
            mx <- match(inter, x@pos)
            my <- match(inter, y@pos)
            vals <- x@value[mx] * y@value[my]
            keep <- which(vals != 0)
            if (length(keep) == 0L) return(new("sparse_numeric", value = numeric(0), pos = integer(0), length = x@length))
            new("sparse_numeric", value = vals[keep], pos = as.integer(inter[keep]), length = x@length)
          })

#' Sparse crossproduct (dot product)
#'
#' Compute the dot product between two sparse vectors (scalar numeric result).
#'
#' @inheritParams sparse_add
#' @return numeric scalar
#' @export
setMethod("sparse_crossprod", signature(x = "sparse_numeric", y = "sparse_numeric"),
          function(x, y, ...) {
            .check_len(x, y)
            inter <- intersect(x@pos, y@pos)
            if (length(inter) == 0L) return(0)
            mx <- match(inter, x@pos)
            my <- match(inter, y@pos)
            sum(x@value[mx] * y@value[my])
          })

# Methods to allow numeric operands (coerce numeric to sparse)
setMethod("sparse_add", c("sparse_numeric", "numeric"), function(x, y, ...) sparse_add(x, as(y, "sparse_numeric")))
setMethod("sparse_add", c("numeric", "sparse_numeric"), function(x, y, ...) sparse_add(as(x, "sparse_numeric"), y))
setMethod("sparse_sub", c("sparse_numeric", "numeric"), function(x, y, ...) sparse_sub(x, as(y, "sparse_numeric")))
setMethod("sparse_sub", c("numeric", "sparse_numeric"), function(x, y, ...) sparse_sub(as(x, "sparse_numeric"), y))
setMethod("sparse_mult", c("sparse_numeric", "numeric"), function(x, y, ...) sparse_mult(x, as(y, "sparse_numeric")))
setMethod("sparse_mult", c("numeric", "sparse_numeric"), function(x, y, ...) sparse_mult(as(x, "sparse_numeric"), y))
setMethod("sparse_crossprod", c("sparse_numeric", "numeric"), function(x, y, ...) sparse_crossprod(x, as(y, "sparse_numeric")))
setMethod("sparse_crossprod", c("numeric", "sparse_numeric"), function(x, y, ...) sparse_crossprod(as(x, "sparse_numeric"), y))

#' @rdname sparse_add
#' @exportMethod +
setMethod("+", c("sparse_numeric", "sparse_numeric"), function(e1, e2) sparse_add(e1, e2))
# support numeric + sparse and sparse + numeric
setMethod("+", c("sparse_numeric", "numeric"), function(e1, e2) sparse_add(e1, as(e2, "sparse_numeric")))
setMethod("+", c("numeric", "sparse_numeric"), function(e1, e2) sparse_add(as(e1, "sparse_numeric"), e2))

#' @rdname sparse_sub
#' @exportMethod -
setMethod("-", c("sparse_numeric", "sparse_numeric"), function(e1, e2) sparse_sub(e1, e2))
setMethod("-", c("sparse_numeric", "numeric"), function(e1, e2) sparse_sub(e1, as(e2, "sparse_numeric")))
setMethod("-", c("numeric", "sparse_numeric"), function(e1, e2) sparse_sub(as(e1, "sparse_numeric"), e2))

#' @rdname sparse_mult
#' @exportMethod *
setMethod("*", c("sparse_numeric", "sparse_numeric"), function(e1, e2) sparse_mult(e1, e2))
setMethod("*", c("sparse_numeric", "numeric"), function(e1, e2) sparse_mult(e1, as(e2, "sparse_numeric")))
setMethod("*", c("numeric", "sparse_numeric"), function(e1, e2) sparse_mult(as(e1, "sparse_numeric"), e2))

#' Show method for sparse_numeric
#'
#' Compact print for sparse vectors.
#' @param object sparse_numeric
#' @export
setMethod("show", "sparse_numeric", function(object) {
  cat("An object of class \"sparse_numeric\"\n")
  cat("  length :", object@length, "\n")
  if (length(object@pos) == 0L) {
    cat("  (all zeros)\n")
  } else {
    nshow <- min(10L, length(object@pos))
    idx <- seq_len(nshow)
    cat("  non-zero entries (first", nshow, "shown):\n")
    df <- data.frame(pos = object@pos[idx], value = object@value[idx])
    print(df, row.names = FALSE)
    if (length(object@pos) > nshow) cat("  ... (", length(object@pos) - nshow, "more)\n")
  }
})

#' Plot overlapping non-zero entries of two sparse vectors
#'
#' Creates a scatterplot of the overlapping non-zero positions for two sparse
#' vectors. If there are no overlapping positions, an empty plot with a title
#' is produced.
#'
#' @param x,y sparse_numeric objects
#' @param ... graphical parameters passed to \code{plot}
#' @export
setMethod("plot", signature(x = "sparse_numeric", y = "sparse_numeric"),
          function(x, y, ...) {
            .check_len(x, y)
            inter <- intersect(x@pos, y@pos)
            if (length(inter) == 0L) {
              plot.new()
              title(main = "No overlapping non-zero positions")
              return(invisible(NULL))
            }
            mx <- match(inter, x@pos); my <- match(inter, y@pos)
            xv <- x@value[mx]; yv <- y@value[my]
            plot(inter, xv, xlab = "position", ylab = "value",
                 main = "Overlapping non-zero elements", pch = 19, ...)
            points(inter, yv, pch = 17)
            legend("topright", legend = c("x", "y"), pch = c(19,17))
            invisible(NULL)
          })

#' Mean for sparse_numeric
#'
#' Calculate the mean of a sparse vector including zeros, without converting to
#' dense format.
#'
#' @param x sparse_numeric
#' @return numeric scalar
#' @export
setMethod("mean", "sparse_numeric", function(x, ...) {
  if (x@length == 0L) return(NaN)
  sum(x@value) / x@length
})

#' Norm (L2) for sparse_numeric
#'
#' Compute the Euclidean norm of the vector (sqrt of sum of squares).
#' Operates on the non-zero entries only.
#'
#' @param x sparse_numeric
#' @return numeric scalar
#' @export
setMethod("norm", "sparse_numeric", function(x, ...) {
  sqrt(sum(x@value^2))
})

#' Standardize a sparse_numeric vector
#'
#' Returns a \code{sparse_numeric} whose entries are (x - mean) / sd. This
#' function computes the mean and standard deviation without expanding to dense
#' format. Note that standardization typically produces many non-zero entries;
#' the returned object is still a \code{sparse_numeric} representation of the
#' standardized values (i.e. zeros omitted).
#'
#' @param x sparse_numeric
#' @return sparse_numeric
#' @examples
#' x <- as(c(1,0,3), "sparse_numeric"); standardize(x)
#' @export
setMethod("standardize", "sparse_numeric", function(x, ...) {
  n <- x@length
  if (n == 0L) {
    return(new("sparse_numeric", value = numeric(0), pos = integer(0), length = 0L))
  }
  m <- mean(x)
  nz <- length(x@pos)
  z <- n - nz
  # deviations for non-zero entries
  dev <- x@value - m
  ssq <- sum(dev^2) + z * (m^2)
  # sample standard deviation (n-1)
  if (n - 1L <= 0L) {
    sd <- 0
  } else {
    sd <- sqrt(ssq / (n - 1L))
  }
  if (sd == 0) {
    return(new("sparse_numeric", value = numeric(0), pos = integer(0), length = n))
  }
  # baseline for zero positions: (-m)/sd
  # compute standardized values at non-zero positions using dev/sd
  vals <- rep((-m) / sd, n)
  if (nz > 0L) vals[x@pos] <- dev / sd
  nzpos <- which(vals != 0)
  new("sparse_numeric", value = vals[nzpos], pos = as.integer(nzpos), length = n)
})
