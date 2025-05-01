
# unit test functions in parse_formula.R and prep_data.R


test_that("test get_eff_mod_names", {
  expect_equal(get_eff_mod_names(y ~ x), character(0))
  expect_equal(get_eff_mod_names(y ~ x + z), character(0))
  expect_equal(get_eff_mod_names(y ~ x*z, trt_var = "x"), "z") 
  expect_equal(get_eff_mod_names(y ~ x + x*z, trt_var = "x"), "z") 
  expect_equal(get_eff_mod_names(y ~ x*z + c, trt_var = "x"), "z") 
})

test_that("test get_treatment_name", {
  expect_equal(get_treatment_name(y ~ x), "x")
  expect_equal(get_treatment_name(y ~ x*z), "x")
  expect_equal(get_treatment_name(y ~ x:z), "x")
  expect_equal(get_treatment_name(y ~ x:s + x:z), "x")
  expect_equal(get_treatment_name(y ~ s:x + x:z), "x")
  expect_equal(get_treatment_name(y ~ s:x + z:x), "x")
  expect_equal(get_treatment_name(y ~ x + s:x), "x")
})

test_that("get_mean_names", {
  expect_equal(unname(get_mean_names(data.frame("mean.x"=NA, "x"=NA), "x")), "mean.x")
  expect_equal(unname(get_mean_names(data.frame("mean.x"=NA), "x")), "mean.x")
  expect_warning(get_mean_names(data.frame("x"=NA), "x"), regexp = "No matching mean names found.")

  # suppressWarnings(expect_equal(get_mean_names(data.frame("x"=NA), "x"), character(0)))
  # suppressWarnings(expect_equal(get_mean_names(data.frame("mean."=NA, "x"=NA), "x"), character(0)))
})

test_that("get_sd_names", {
  expect_equal(unname(get_sd_names(data.frame("sd.x"=NA, "x"=NA), "x")), "sd.x")
  expect_equal(unname(get_sd_names(data.frame("sd.x"=NA), "x")), "sd.x")
  expect_warning(get_sd_names(data.frame("x"=NA), "x"), regexp = "No matching sd names found.")

  # suppressWarnings(expect_equal(get_sd_names(data.frame("x"=NA), "x"), character(0)))
  # suppressWarnings(expect_equal(get_sd_names(data.frame("sd."=NA, "x"=NA), "x"), character(0)))
})

test_that("get_covariate_names", {
  expect_equal(get_covariate_names(y ~ x), "x")
  expect_equal(get_covariate_names(y ~ x + z), c("x","z"))
  expect_equal(get_covariate_names(y ~ x*z), c("x","z"))
  expect_equal(get_covariate_names(y ~ x*z + z), c("x","z"))
})


test_that("prep_ald different formula formats", {
  
  # mock data
  data <- data.frame(
    mean.X1 = 0.6065438,
    sd.X1    = 0.3917644,
    mean.X2  = 0.6273452,
    sd.X2    = 0.4031580,
    mean.X3  = 0.5579042,
    sd.X3    = 0.4077155,
    mean.X4  = 0.6122686,
    sd.X4    = 0.4466750,
    y.C.sum  = 36,
    y.C.bar  = 0.5373134,
    y.C.sd   = 0.5023689,
    N.C      = 67,
    y.B.sum  = 35,
    y.B.bar  = 0.2631579,
    y.B.sd   = 0.4420122,
    N.B      = 133,
    stringsAsFactors = FALSE
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

