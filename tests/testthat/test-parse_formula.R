
# unit test functions in parse_formula.R


test_that("test get_effect_modifiers", {
  expect_equal(get_effect_modifiers(y ~ x), character(0))
  expect_equal(get_effect_modifiers(y ~ x + z), character(0))
  expect_equal(get_effect_modifiers(y ~ x*z), "z") 
  expect_equal(get_effect_modifiers(y ~ x + x*z), "z") 
  expect_equal(get_effect_modifiers(y ~ x*z + c), "z") 
})

test_that("test get_treatment_name", {
  expect_warning(get_treatment_name(y ~ x), regexp = "Treatment name missing from formula.")
  suppressWarnings(expect_equal(get_treatment_name(y ~ x), NA_character_))
  expect_equal(get_treatment_name(y ~ x*z), "x")
})

test_that("get_mean_names", {
  expect_equal(get_mean_names(data.frame("mean.x"=NA, "x"=NA), "x"), "mean.x")
  expect_equal(get_mean_names(data.frame("mean.x"=NA), "x"), "mean.x")
  expect_warning(get_mean_names(data.frame("x"=NA), "x"), regexp = "No matching mean names found.")
  suppressWarnings(expect_equal(get_mean_names(data.frame("x"=NA), "x"), character(0)))
  suppressWarnings(expect_equal(get_mean_names(data.frame("mean."=NA, "x"=NA), "x"), character(0)))
})

test_that("get_sd_names", {
  expect_equal(get_sd_names(data.frame("sd.x"=NA, "x"=NA), "x"), "sd.x")
  expect_equal(get_sd_names(data.frame("sd.x"=NA), "x"), "sd.x")
  expect_warning(get_sd_names(data.frame("x"=NA), "x"), regexp = "No matching sd names found.")
  suppressWarnings(expect_equal(get_sd_names(data.frame("x"=NA), "x"), character(0)))
  suppressWarnings(expect_equal(get_sd_names(data.frame("sd."=NA, "x"=NA), "x"), character(0)))
})

test_that("get_covariate_names", {
  expect_equal(get_covariate_names(y ~ x), "x")
  expect_equal(get_covariate_names(y ~ x + z), c("x","z"))
  expect_equal(get_covariate_names(y ~ x*z), c("x","z"))
  expect_equal(get_covariate_names(y ~ x*z + z), c("x","z"))
})

