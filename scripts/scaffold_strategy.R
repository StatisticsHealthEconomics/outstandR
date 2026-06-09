
#' Scaffold a new outstandR strategy
#' 
#' Automates the file creation and editing steps required to add a new 
#' statistical modelling method to the package.
#' 
#' @param name The short, snake_case name of the new method (e.g., "custom_method")
scaffold_strategy <- function(name) {
  
  if (!requireNamespace("usethis", quietly = TRUE)) {
    stop("The 'usethis' package is required for scaffolding. Please install it.")
  }
  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("The 'cli' package is required for scaffolding. Please install it.")
  }
  
  cli::cli_h1("Scaffolding outstandR Strategy: {.val {name}}")
  
  # 1. Create the core calculation file
  calc_file <- paste0("calc_", name)
  cli::cli_alert_info("Creating core calculation file: {.file R/{calc_file}.R}")
  usethis::use_r(calc_file)
  
  # 2. Open the strategy constructor file
  cli::cli_alert_info("Opening {.file R/strategy_.R} for constructor definition...")
  usethis::edit_file("R/strategy_.R")
  
  # 3. Open the S3 registry file
  cli::cli_alert_info("Opening {.file R/calc_IPD_stats.R} to register S3 method...")
  usethis::edit_file("R/calc_IPD_stats.R")
  
  # 4. Open the print methods file
  cli::cli_alert_info("Opening {.file R/print-strategies.R} to add print method...")
  usethis::edit_file("R/print-strategies.R")
  
  cli::cli_alert_success("Scaffold complete! Follow the developer guide to fill in the functions.")
}

