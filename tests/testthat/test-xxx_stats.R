# unit tests for ALD_stats()

test_that("ALD_stats() returns the correct values", {
  # test 1
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 2, 3, 4, 5)
  expect_equal(ALD_stats(x, y), c(0, 0, 0, 0, 0))
  # test 2
  x <- c(1, 2, 3, 4, 5)
  y <- c(5, 4, 3, 2, 1)
  expect_equal(ALD_stats(x, y), c(4, 4, 4, 4, 4))
  # test 3
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 1, 1, 1, 1)
  expect_equal(ALD_stats(x, y), c(0, 0, 0, 0, 0))
  # test 4
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 2, 3, 4, 5)
  expect_equal(ALD_stats(x, y), c(0, 0, 0, 0, 0))
  # test 5
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 2, 3, 4, 5)
  expect_equal(ALD_stats(x, y), c(0, 0, 0, 0, 0))
  # test 6
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 2, 3, 4, 5)
  expect_equal(ALD_stats(x, y), c(0, 0, 0, 0, 0))
  # test 7
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 2, 3, 4, 5)
  expect_equal(ALD_stats(x, y), c(0, 0, 0, 0, 0))
})


# unit tests for IPD_stats()

test_that("IPD_stats() returns the correct values", {
  # test 1
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 2, 3, 4, 5)
  expect_equal(IPD_stats(x, y), c(0, 0, 0, 0, 0))
  # test 2
  x <- c(1, 2, 3, 4, 5)
  y <- c(5, 4, 3, 2, 1)
  expect_equal(IPD_stats(x, y), c(4, 4, 4, 4, 4))
  # test 3
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 1, 1, 1, 1)
  expect_equal(IPD_stats(x, y), c(0, 0, 0, 0, 0))
  # test 4
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 2, 3, 4, 5)
  expect_equal(IPD_stats(x, y), c(0, 0, 0, 0, 0))
  # test 5
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 2, 3, 4, 5)
  expect_equal(IPD_stats(x, y), c(0, 0, 0, 0, 0))
  # test 6
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 2, 3, 4, 5)
  expect_equal(IPD_stats(x, y), c(0, 0, 0, 0, 0))
  # test 7
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 2, 3, 4, 5)
  expect_equal(IPD_stats(x, y), c(0, 0, 0, 0, 0))
})


