#

# Mock Data
base_formula <- y ~ trt * age   # single covariate
formula1 <- y ~ X1 + trt * age  # EM and PF

mock_ipd <- data.frame(
  trt = c(0, 1),
  y = c(1, 0),
  age = c(30, 40),
  X1 = c(300, 400)
)

mock_ald <- data.frame(
  mean.age = 35,
  sd.age = 5,
  mean.X1 = 350,
  sd.X1 = 10,
  y.B.sum = 30,
  N.B = 100,
  y.C.sum = 20,
  N.C = 100
)

# Basic Functionality Test
test_that("simulate_ALD_pseudo_pop generates a data frame with correct dimensions", {
  
  result <- simulate_ALD_pseudo_pop(
    formula = base_formula,
    ipd = mock_ipd,
    ald = mock_ald,
    N = 1000)
  
  expect_type(result, type = "double")
  expect_equal(nrow(result), 1000)

  result <- simulate_ALD_pseudo_pop(
    formula = formula1,
    ipd = mock_ipd,
    ald = mock_ald,
    N = 1000)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1000)
})


# Output Validation
test_that("simulate_ALD_pseudo_pop generates expected columns", {
  
  result <- simulate_ALD_pseudo_pop(
    formula = base_formula,
    ipd = mock_ipd,
    ald = mock_ald,
    N = 1000)
  
  expect_true("age" %in% colnames(result))

  result <- simulate_ALD_pseudo_pop(
    formula = formula1,
    ipd = mock_ipd,
    ald = mock_ald,
    N = 1000)
  
  expect_true("age" %in% colnames(result))
  expect_true("X1" %in% colnames(result))
})
