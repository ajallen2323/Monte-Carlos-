# Author:  Angelina Allen
# Date:    2024-05-07
# Purpose: Portfolio 3- Monte Carlos
#-------------------------------------------------------------------------------

# Load packages
suppressPackageStartupMessages(library("tidyverse")); theme_set(theme_bw())
library("knitr")
library(dbplyr)

n <- c(5, 10, 15, 20, 25, 30)

p <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

df_p <- data.frame(p = p)
df_n <- data.frame(n = n)
values <- merge(df_p, df_n)

# Function to calculate Monte Carlo expectation with standard error
mc_expectation <- function(x) {
  return(
    c(mean = mean(x),
      se   = sd(x)/sqrt(length(x)))
  )
}

#Frequentist

#Calculate MSE for each row
for (i in 1:nrow(values)){
  n <- values$n[i]
  p <- values$p[i]
  x <- rbinom(1e4, size = n, prob = p)
  theta_hat <- x / n
  values$freq_MSE[i] <- mc_expectation((theta_hat - p)^2)[[1]]
}


#Baysian

#Calculate MSE for each row
for (i in 1:nrow(values)){
  n <- values$n[i]
  p <- values$p[i]
  x <- rbinom(1e4, size = n, prob = p)
  theta_hat <- (0.5 + x) / (1 + n) 
  values$bays_MSE[i] <- mc_expectation((theta_hat - p)^2)[[1]]
}

MSEdata <- pivot_longer(values, cols = c(freq_MSE, bays_MSE), 
                          values_to = "MSE", names_to ="MSEtype")


ggplot(data = MSEdata,
            aes(x = p,
                y = MSE, 
                color = MSEtype)) + 
  geom_line() + 
  facet_grid(~n, scales = "free") +
  labs(x = "Probability (p)",
       y = "Mean Squared Error (MSE)",
       color = "Approach",
       title = "Comparison of Mean Squared Error (MSE) between Frequentist and
                Bayesian Methods") +
  scale_color_manual(values = c("red", "blue"), labels = c("Bayesian Approach", "Frequentist Approach")) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#MSE Monte Carlo uncertainty

mse_freq <- values$freq_MSE
mse_bays <- values$bays_MSE

mse_dif <- (mse_freq - mse_bays)^2

mc_expectation(mse_dif)[2]


#Frequentist

#Calculate Coverage for each row
alpha <- 0.05 
for (i in 1:nrow(values)){
  n <- values$n[i]
  p <- values$p[i]
  z <- -qnorm(alpha/2)
  d <- data.frame(x = rbinom(1e4,
                             size = n,
                             prob = p)) 
  ci <- d |>
    rowwise() |>
    mutate(
      # Wald interval
      p_hat = x/n,
      lowerbound = p_hat - (z * sqrt(p_hat*(1-p_hat)/n)),
      upperbound = p_hat + z * sqrt(p_hat*(1-p_hat)/n),
      Wald = lowerbound < p & p < upperbound
    ) |>
    select(Wald)
  values$freq_cov[i] <- sum(ci$Wald) / length(ci$Wald)
}

#Bayesian

#Calculate Coverage for each row

alpha <- 0.5

for (i in 1:nrow(values)){
  n <- values$n[i]
  p <- values$p[i]
  x <- rbinom(1e4, size = n, prob = p)
  b <- qbeta(c(alpha/2, 1-alpha/2), 0.5 + x, 0.5 + n - x)
  
  d <- data.frame(lowerbound = b[seq(1, length(b), by = 2)],
                  upperbound = b[seq(2, length(b), by = 2)]) |>
         mutate(
           within = lowerbound < p & p < upperbound)|>
         select(within)
  
  values$bays_cov[i] <- sum(d$within) / length(d$within)
}

Covdata <- pivot_longer(values, cols = c(freq_cov, bays_cov), 
                          values_to = "COV", names_to ="Covtype")

ggplot(data = Covdata,
            aes(x = p,
                y = COV, 
                color = Covtype)) + 
  geom_line() + 
  facet_grid(~n, scales = "free") +
  labs(x = "Probability (p)",
       y = "Coverage",
       color = "Approach",
       title = "Comparison of Coverage between Frequentist and
                Bayesian Methods") +
  scale_color_manual(values = c("red", "blue"), labels = c("Bayesian Approach", "Frequentist Approach")) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Coverage Monte Carlo uncertainty

cov_freq <- values$freq_cov
cov_bays <- values$bays_cov

cov_dif <- (cov_freq - cov_bays)^2

mc_expectation(cov_dif)[2]

