## see the NICE Guidelines Technical Support Unit Meta-Analysis of Event Outcomes
## Guideline Methodology Document 3 Version 1 (January 2021)
## for how to convert from probability to different relative effect statistics

# # compute baseline probability in control group (P0)
# newdat <- ipd[ipd[[treat_nm]] == 0, ]
# P0 <- mean(predict(fit, newdata = newdat, type = "response"))
# 
# converted_effect <- convert_effect(
#   coef_est, link = strategy$family$link, to_scale = scale, P0)

# transformation functions
convert_lor_to_or <- function(lor, P0 = NULL) exp(lor)
convert_or_to_lor <- function(or, P0 = NULL) log(or)

convert_lrr_to_rr <- function(lrr, P0 = NULL) exp(lrr)
convert_rr_to_lrr <- function(rr, P0 = NULL) log(rr)

convert_or_to_rr <- function(or, P0) or / ((1 - P0) + (P0 * or))
convert_rr_to_or <- function(rr, P0) rr * ((1 - P0) / (1 - P0 * rr))

convert_rr_to_rd <- function(rr, P0) P0 * (rr - 1)
convert_rd_to_rr <- function(rd, P0) 1 + (rd / P0)

#
conversion_map <- function() {
  list(
    log_odds =
      list(odds_ratio = convert_lor_to_or),
    odds_ratio =
      list(log_odds = convert_or_to_lor,
           relative_risk = convert_or_to_rr),
    relative_risk =
      list(odds_ratio = convert_rr_to_or,
           log_relative_risk = convert_rr_to_lrr,
           risk_difference = convert_rr_to_rd),
    risk_difference =
      list(relative_risk = convert_rd_to_rr),
    delta_z = 
      list(odds_ratio = convert_delta_z_to_or),
  )
}

#' recursive function to find conversion path
#'
#' @examples
#' find_conversion_path("log_odds", "relative_risk")
#' 
find_conversion_path <- function(from, to, visited = c()) {
  
  conv_map <- conversion_map()
  
  # found target
  if (from == to) return(list(from))
  
  # Avoid revisiting nodes
  if (from %in% visited) return(NULL)
  
  # If there is no valid conversion from 'from', return NULL
  if (!from %in% names(conversion_map)) return(NULL)
  
  visited <- c(visited, from)  # mark as visited
  
  # try each possible conversion from 'from'
  for (next_step in names(conv_map[[from]])) {
    
    path <- find_conversion_path(next_step, to, visited)
    
    if (!is.null(path)) return(c(from, path))  # return first valid path found
  }
  
  return(NULL)  # no valid path found
}

# apply conversions along the found path
#
convert_effect <- function(value, link, to, P0) {
  
  from <- get_treatment_effect(link)  # scale
  conv_map <- conversion_map()
  
  path <- find_conversion_path(from, to)
  
  if (is.null(path)) stop("No valid conversion path found.")
  
  result <- value
  
  for (i in seq_along(path)[-1]) {
    prev <- path[i - 1]
    nxt <- path[i]
    result <- conv_map[[prev]][[nxt]](result, P0)
  }
  
  result
}

