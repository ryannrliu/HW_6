#' Sparse Numeric Vector Class and Methods
#'
#' This file defines the sparse_numeric S4 class and associated methods
#' including coercion, arithmetic, mean(), norm(), and standardize().
#'
#' @import methods
NULL

setClass(
  Class = "sparse_numeric",
  slots = c(
    value = "numeric",
    pos = "integer",
    length = "integer"
  )
)

#' @title Validate a sparse_numeric Object
#' @description Internal validity checking for sparse_numeric.
setValidity("sparse_numeric", function(object) {
  if (length(object@length) != 1L || object@length < 0L)
    return("slot 'length' must be a single non-negative integer")
  if (!is.integer(object@pos))
    return("slot 'pos' must be integer")
  if (!is.numeric(object@value))
    return("slot 'value' must be numeric")
  if (length(object@value) != length(object@pos))
    return("'value' and 'pos' must have same length")
  if (length(object@pos) > 0L) {
    if (any(object@pos < 1L | object@pos > object@length))
      return("positions out of bounds")
    if (any(duplicated(object@pos)))
      return("positions must be unique")
  }
  TRUE
})

#' @title Coerce numeric to sparse_numeric
#' @param from Numeric vector
#' @return A sparse_numeric object
#' @export
setAs("numeric", "sparse_numeric", function(from) {
  nz <- which(from != 0)
  new("sparse_numeric",
      value = if (length(nz)) from[nz] else numeric(),
      pos = if (length(nz)) as.integer(nz) else integer(),
      length = as.integer(length(from)))
})

#' @title Coerce sparse_numeric to numeric
#' @param from sparse_numeric object
#' @return Dense numeric vector
#' @export
setAs("sparse_numeric", "numeric", function(from) {
  out <- numeric(from@length)
  if (length(from@pos)) out[from@pos] <- from@value
  out
})

.check_len <- function(x, y) {
  if (x@length != y@length) stop("Lengths differ")
}

#' @title Sparse Addition
#' @param x,y sparse_numeric objects
#' @return sparse_numeric result
#' @export
setGeneric("sparse_add", function(x, y, ...) standardGeneric("sparse_add"))

setMethod("sparse_add", c("sparse_numeric", "sparse_numeric"), function(x, y) {
  .check_len(x, y)
  allpos <- sort(unique(c(x@pos, y@pos)))

  xv <- yv <- numeric(length(allpos))

  mx <- match(allpos, x@pos)
  xv[!is.na(mx)] <- x@value[mx[!is.na(mx)]]

  my <- match(allpos, y@pos)
  yv[!is.na(my)] <- y@value[my[!is.na(my)]]

  vals <- xv + yv
  keep <- vals != 0

  new("sparse_numeric", value=vals[keep], pos=as.integer(allpos[keep]), length=x@length)
})

#' @title Sparse Subtraction
#' @export
setGeneric("sparse_sub", function(x, y, ...) standardGeneric("sparse_sub"))

setMethod("sparse_sub", c("sparse_numeric", "sparse_numeric"), function(x, y) {
  .check_len(x, y)
  allpos <- sort(unique(c(x@pos, y@pos)))
  xv <- yv <- numeric(length(allpos))
  if (length(x@pos)) xv[match(allpos, x@pos, nomatch = 0)] <- x@value
  if (length(y@pos)) yv[match(allpos, y@pos, nomatch = 0)] <- y@value
  vals <- xv - yv
  keep <- which(vals != 0)
  new("sparse_numeric", value = vals[keep], pos = as.integer(allpos[keep]), length = x@length)
})

#' @title Sparse Multiplication
#' @export
setGeneric("sparse_mult", function(x, y, ...) standardGeneric("sparse_mult"))

setMethod("sparse_mult", c("sparse_numeric", "sparse_numeric"), function(x, y) {
  .check_len(x, y)
  inter <- intersect(x@pos, y@pos)
  if (!length(inter)) return(new("sparse_numeric", value = numeric(), pos = integer(), length = x@length))
  xv <- x@value[match(inter, x@pos)]
  yv <- y@value[match(inter, y@pos)]
  vals <- xv * yv
  keep <- which(vals != 0)
  new("sparse_numeric", value = vals[keep], pos = as.integer(inter[keep]), length = x@length)
})

#' @title Sparse Crossproduct
#' @export
setGeneric("sparse_crossprod", function(x, y, ...) standardGeneric("sparse_crossprod"))

setMethod("sparse_crossprod", c("sparse_numeric", "sparse_numeric"), function(x, y) {
  .check_len(x, y)
  inter <- intersect(x@pos, y@pos)
  if (!length(inter)) return(0)
  xv <- x@value[match(inter, x@pos)]
  yv <- y@value[match(inter, y@pos)]
  sum(xv * yv)
})

#' @title Show sparse_numeric
#' @export
setMethod("show", "sparse_numeric", function(object) {
  cat("sparse_numeric of length", object@length, "\n")
  if (!length(object@pos)) {
    cat("  (all zeros)\n")
  } else {
    print(data.frame(pos = object@pos, value = object@value))
  }
})

#' @title Plot Overlapping Nonzero Elements
#' @export
setMethod("plot", c("sparse_numeric", "sparse_numeric"), function(x, y, ...) {
  inter <- intersect(x@pos, y@pos)
  if (!length(inter)) {
    plot.new(); title("No overlap"); return()
  }
  xv <- x@value[match(inter, x@pos)]
  yv <- y@value[match(inter, y@pos)]
  plot(inter, xv, xlab = "position", ylab = "x value", pch = 19, ...)
  points(inter, yv, pch = 17, col = 2)
  legend("topright", legend = c("x", "y"), pch = c(19,17), col = c(1,2))
})

#' @title Mean for sparse_numeric
#' @description Computes mean without dense expansion.
#' @param x sparse_numeric
#' @return numeric
#' @export
setMethod("mean", "sparse_numeric", function(x, ...) sum(x@value) / x@length)

#' @title Norm of a Sparse Vector
#' @description Computes L2 norm.
#' @param x sparse_numeric
#' @return numeric
#' @export
setGeneric("norm", function(x, ...) standardGeneric("norm"))

setMethod("norm", "sparse_numeric", function(x, ...) sqrt(sum(x@value^2)))

#' @title Standardize a Sparse Vector
#' @description Returns (x - mean) / sd as sparse_numeric.
#' @export
setGeneric("standardize", function(x, ...) standardGeneric("standardize"))

setMethod("standardize", "sparse_numeric", function(x) {

  n <- x@length
  if (n == 0)
    return(new("sparse_numeric", value=numeric(), pos=integer(), length=0L))

  m <- mean(x)

  # compute SD
  nz <- length(x@pos)
  z <- n - nz
  dev <- x@value - m
  ssq <- sum(dev^2) + z * m^2

  sd <- sqrt(ssq / (n - 1))
  if (sd == 0)
    return(new("sparse_numeric", value=numeric(), pos=integer(), length=n))

  vals <- rep((-m) / sd, n)
  vals[x@pos] <- dev / sd

  nzpos <- which(vals != 0)
  new("sparse_numeric", value=vals[nzpos], pos=as.integer(nzpos), length=n)
})
setGeneric("+", function(e1,e2) standardGeneric("+"))
setGeneric("-", function(e1,e2) standardGeneric("-"))
setGeneric("*", function(e1,e2) standardGeneric("*"))
