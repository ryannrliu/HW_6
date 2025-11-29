library(testthat)
source("hw_6.R")



test_that("validObject works for correct object", {
  x <- new("sparse_numeric", value = c(1,2), pos = c(2L,4L), length = 5L)
  expect_true(validObject(x))
})


test_that("invalid length is detected", {
  x <- new("sparse_numeric", value = c(1), pos = c(1L), length = 1L)
  x@length <- 0L
  expect_error(validObject(x))
})


test_that("coercion numeric -> sparse and back works", {
  v <- c(0,3,0,4)
  s <- as(v, "sparse_numeric")
  expect_s4_class(s, "sparse_numeric")
  expect_equal(as(s, "numeric"), v)
})



test_that("sparse_add produces correct result", {
  x <- as(c(0,0,1,2), "sparse_numeric")
  y <- as(c(1,0,2,0), "sparse_numeric")
  res <- sparse_add(x, y)
  expect_s4_class(res, "sparse_numeric")
  expect_equal(as(res, "numeric"), c(1,0,3,2))
})


test_that("operator + works with numeric on either side", {
  x <- as(c(1,0,2), "sparse_numeric")
  y <- c(0,1,1)
  r1 <- x + as(y, "sparse_numeric")
  r2 <- as(x, "numeric") + y
  expect_equal(as(r1, "numeric"), r2)
  r3 <- as(x, "sparse_numeric") + y
  expect_equal(as(r3, "numeric"), r2)
})


test_that("sparse_sub and - operator work", {
  x <- as(c(1,2,0), "sparse_numeric")
  y <- as(c(1,1,1), "sparse_numeric")
  expect_equal(as(sparse_sub(x,y), "numeric"), c(0,1,-1))
  expect_equal(as(x - y, "numeric"), c(0,1,-1))
})


test_that("sparse_mult and * operator work", {
  x <- as(c(1,0,2,0), "sparse_numeric")
  y <- as(c(0,2,1,0), "sparse_numeric")
  expect_equal(as(sparse_mult(x,y), "numeric"), c(0,0,2,0))
  expect_equal(as(x * y, "numeric"), c(0,0,2,0))
})


test_that("mean computes correctly without densifying", {
  x <- as(c(0,0,3,0), "sparse_numeric")
  expect_equal(mean(x), 3/4)
})


test_that("norm computes L2 norm", {
  x <- as(c(3,4,0), "sparse_numeric")
  expect_equal(norm(x), 5)
})


test_that("sparse_crossprod equals dot product", {
  x <- as(c(1,0,2,0), "sparse_numeric")
  y <- as(c(0,2,1,0), "sparse_numeric")
  expect_equal(sparse_crossprod(x,y), sum(as(x, "numeric") * as(y, "numeric")))
})

test_that("standardize returns zeros for constant vector", {
  x <- as(rep(3, 5), "sparse_numeric")
  s <- standardize(x)
  expect_s4_class(s, "sparse_numeric")
  expect_equal(length(s@pos), 0)
  expect_equal(as(s, "numeric"), rep(0,5))
})


test_that("standardize produces expected values for a simple vector", {
  v <- c(1,0,3)
  x <- as(v, "sparse_numeric")
  s <- standardize(x)
  dense_s <- as(s, "numeric")
  target <- scale(v, center = TRUE, scale = TRUE)[,1]
  expect_equal(round(dense_s, 8), round(target, 8))
})

test_that("show prints without error", {
  x <- as(c(0,2,0), "sparse_numeric")
  expect_output(show(x), "sparse_numeric")
})


test_that("plot does not error on overlap and no-overlap", {
  x <- as(c(0,1,0), "sparse_numeric")
  y <- as(c(0,2,0), "sparse_numeric")
  expect_silent(plot(x, y))
  z <- as(c(1,0,0), "sparse_numeric")
  expect_silent(plot(x, z))
})

test_that("empty sparse vector is valid", {
  s <- new("sparse_numeric", value = numeric(), pos = integer(), length = 0L)
  expect_true(validObject(s))
})

test_that("duplicate positions trigger validity error", {
  expect_error(
    new("sparse_numeric", value = c(1,2), pos = c(2L,2L), length = 5L)
  )
})

test_that("out-of-range positions trigger error", {
  expect_error(
    new("sparse_numeric", value = 5, pos = 10L, length = 5L)
  )
})

test_that("negative length triggers error", {
  expect_error(
    new("sparse_numeric", value = 1, pos = 1L, length = -1L)
  )
})


### COERCION EDGE CASES -----------------------------------------------------

test_that("coercing empty numeric vector to sparse works", {
  v <- numeric(0)
  s <- as(v, "sparse_numeric")
  expect_equal(s@length, 0)
  expect_equal(length(s@pos), 0)
})

test_that("coercion preserves names (should drop silently)", {
  v <- c(a = 1, b = 0, c = 2)
  s <- as(v, "sparse_numeric")
  expect_equal(as(s, "numeric"), unname(v))
})


### ARITHMETIC EDGE CASES ---------------------------------------------------

test_that("adding two all-zero sparse vectors returns empty sparse", {
  x <- as(rep(0, 5), "sparse_numeric")
  y <- as(rep(0, 5), "sparse_numeric")
  res <- sparse_add(x, y)
  expect_equal(length(res@pos), 0)
  expect_equal(as(res, "numeric"), rep(0, 5))
})

test_that("subtraction when left side is all zeros", {
  x <- as(rep(0, 4), "sparse_numeric")
  y <- as(c(1,0,2,0), "sparse_numeric")
  res <- sparse_sub(x, y)
  expect_equal(as(res, "numeric"), c(-1,0,-2,0))
})

test_that("multiplication with no overlapping nonzero indices", {
  x <- as(c(1,0,0), "sparse_numeric")
  y <- as(c(0,2,0), "sparse_numeric")
  res <- sparse_mult(x, y)
  expect_equal(length(res@pos), 0)
})

test_that("crossprod of two zero vectors is zero", {
  x <- as(rep(0,10), "sparse_numeric")
  y <- as(rep(0,10), "sparse_numeric")
  expect_equal(sparse_crossprod(x, y), 0)
})


### SUMMARY STATISTICS -------------------------------------------------------

test_that("mean of empty sparse numeric is NaN", {
  x <- new("sparse_numeric", value = numeric(), pos = integer(), length = 0L)
  expect_true(is.nan(mean(x)))
})

test_that("norm of empty sparse vector is zero", {
  x <- new("sparse_numeric", value = numeric(), pos = integer(), length = 0L)
  expect_equal(norm(x), 0)
})


### STANDARDIZE EDGE CASES --------------------------------------------------

test_that("standardizing empty vector returns empty", {
  x <- new("sparse_numeric", value = numeric(), pos = integer(), length = 0L)
  s <- standardize(x)
  expect_equal(length(s@pos), 0)
  expect_equal(s@length, 0)
})

test_that("standardizing vector with one element returns zero", {
  x <- as(5, "sparse_numeric")
  s <- standardize(x)
  expect_equal(as(s, "numeric"), 0)
})


### SHOW & PLOT EDGE CASES --------------------------------------------------

test_that("show handles empty vector", {
  x <- new("sparse_numeric", value = numeric(), pos = integer(), length = 0L)
  expect_output(show(x))
})

test_that("plot handles completely zero vectors", {
  x <- as(rep(0,5), "sparse_numeric")
  y <- as(rep(0,5), "sparse_numeric")
  expect_silent(plot(x, y))
})
