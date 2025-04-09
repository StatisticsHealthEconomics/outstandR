# convert between effect sizes

library(metafor)

# Example aggregate summary statistics
n1 <- 100  # Sample size for group 1
n2 <- 120  # Sample size for group 2
mean1 <- 50  # Mean for group 1
mean2 <- 40  # Mean for group 2
sd1 <- 10  # Standard deviation for group 1
sd2 <- 8   # Standard deviation for group 2

# Compute log odds ratio (LOR) for continuous data
res <- escalc(measure = "OR", 
              m1i = mean1, sd1i = sd1, n1i = n1,
              m2i = mean2, sd2i = sd2, n2i = n2)

res

library(esc)

# Example data for summary measures
mean1 <- 50
mean2 <- 40
sd1 <- 10
sd2 <- 8
n1 <- 100
n2 <- 120

# Calculate log odds ratio (LOR) or risk ratio (RR)
result <- esc_mean_sd(grp1m = mean1, grp1sd = sd1, grp1n = n1,
                      grp2m = mean2, grp2sd = sd2, grp2n = n2,
                      es.type = "logit")  # 'logit' for log odds ratio

# Print result
result
