# unit test functions in parse_formula.R and prep_data.R


test_that("test get_eff_mod_names", {
  expect_equal(get_eff_mod_names(y ~ x), character(0))
  expect_equal(get_eff_mod_names(y ~ x + z), character(0))
  expect_equal(get_eff_mod_names(y ~ x * z, trt_var = "x"), "z")
  expect_equal(get_eff_mod_names(y ~ x + x * z, trt_var = "x"), "z")
  expect_equal(get_eff_mod_names(y ~ x * z + c, trt_var = "x"), "z")
})

test_that("test get_treatment_name", {
  expect_equal(get_treatment_name(y ~ x), "x")
  expect_equal(get_treatment_name(y ~ x * z), "x")
  expect_equal(get_treatment_name(y ~ x:z), "x")
  expect_equal(get_treatment_name(y ~ x:s + x:z), "x")
  expect_equal(get_treatment_name(y ~ s:x + x:z), "x")
  expect_equal(get_treatment_name(y ~ s:x + z:x), "x")
  expect_equal(get_treatment_name(y ~ x + s:x), "x")
})

test_that("get_covariate_names", {
  expect_equal(get_covariate_names(y ~ x), "x")
  expect_equal(get_covariate_names(y ~ x + z), c("x", "z"))
  expect_equal(get_covariate_names(y ~ x * z), c("x", "z"))
  expect_equal(get_covariate_names(y ~ x * z + z), c("x", "z"))
})


test_that("prep_ald different formula formats", {
  
  # mock data
  data <- tibble::tribble(
    ~variable, ~trt, ~statistic, ~value,
    "X1",       NA,  "mean",     30,
    "X1",       NA,  "sd",       30,
    "X2",       NA,  "mean",     30,
    "X2",       NA,  "sd",       30,
    "X3",       NA,  "mean",     30,
    "X3",       NA,  "sd",       30,
    "X4",       NA,  "mean",     30,
    "X4",       NA,  "sd",       30,
    "y",       "B",  "sum",      30,
    "y",       "C",  "sum",      40,
    NA,        "B",  "N",        100,
    NA,        "C",  "N",        120
  )
  
  form <- as.formula("y ~ X3 + X4 + trt*X1 + trt*X2")
  
  expect_equal(prep_ald(form, data), data)
  
  form <- as.formula("y ~ X3 + X4 + trt + X1 + X2 + trt*X1 + trt*X2")
  
  expect_equal(prep_ald(form, data), data)
  
  form <- as.formula("y ~ X3 + X4 + trt*(X1 + X2)")
  
  expect_equal(prep_ald(form, data), data)
  
  form <- as.formula("y ~ X3 + X4 + trt:(X1 + X2)")
  
  expect_equal(prep_ald(form, data), data)
  
  form <- as.formula("y ~ X3 + X4 + trt:X1 + trt:X2")
  
  expect_equal(prep_ald(form, data), data)
  
  form <- as.formula("y ~ X3 + X4 + trt + trt:X1 + trt:X2")
  
  expect_equal(prep_ald(form, data), data)
})
