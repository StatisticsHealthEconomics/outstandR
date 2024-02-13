
# unit test functions in parse_formula.R


test_that("test get_effect_modifiers", {
  expect_equal(get_effect_modifiers(y ~ x), character(0))
  expect_equal(get_effect_modifiers(y ~ x + z), "z")
})

test_that("test get_treatment_name", {
  expect_equal(get_treatment_name(y ~ x), "x")
  expect_equal(get_treatment_name(y ~ x + z), "x")
  expec
})

test_that("get_mean_names", {
  expect_equal(get_mean_names(y ~ x), "y")
  expect_equal(get_mean_names(y ~ x + z), "y")
})

test_that("get_sd_names", {
  expect_equal(get_sd_names(y ~ x), "y")
  expect_equal(get_sd_names(y ~ x + z), "y")
})

test_that("get_covariate_names", {
  expect_equal(get_covariate_names(y ~ x), character(0))
  expect_equal(get_covariate_names(y ~ x + z), "z")
})

