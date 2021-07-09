######### Home assignment 1 ######### 
# Philip Dahlqvist-Sjöberg
# 18-05-2020
# NOTE: Run all these libraries first.
library(ggplot2) # For nicer plots
library(extraDistr) # Beta-Binomial distributions
library(coda) # Optimize MCMC
library(mvtnorm) # Provide Multivariate Normal Density 
library(gridExtra) # Plot multiple ggplots in one frame
library(fGarch) # Skew Normal Distribution
library(dplyr)

# Function to write matrices to Latex code
write_matex <- function(x) {
  begin <- "\\begin{bmatrix}"
  end <- "\\end{bmatrix}"
  X <-
    apply(x, 1, function(x) {
      paste(
        paste(x, collapse = "&"),
        "\\\\"
      )
    })
  writeLines(c(begin, X, end))
}

###### 1 ######

# Gamma function used in true marginal distribution
gamma <- function(x){
  factorial(x-1)
}

# True marginal distribution
f_x <- function(x, a=2, b=4, n=16){
  choose(n,x)*(gamma(a+b)/(gamma(a)*gamma(b)))*((gamma(x+a)*gamma(n-x+b))/(gamma(a+b+n)))
}

# Gibbs sampler 1
# Conditional distribution for both x and y
f_x.y <- function(n=16, y){rbinom(1, n, y)}
f_y.x <- function(x, a=2, b=4, n=16){rbeta(1, x+a, n-x+b)}

# Number of iterations
m <- 500

# Lenght of Gibbs sequence
k <- 10

# Create empty output vectors
x <- rep(NA, k+1)
y <- rep(NA, k+1)
x_sample <- rep(NA, m)
y_sample <- rep(NA, m)

# Iterate the Gibbs sequences
for(j in 1:m){
  # Random starting value for y
  y[1] <- rbeta(1, 2, 4)
  
  for(i in 2:(k+1)){
    x[i] <- f_x.y(y = y[i-1])
    y[i] <- f_y.x(x = x[i])
  }
  x_sample[j] <- x[k+1]
  y_sample[j] <- y[k+1]
}


# Plot 1
# Create grid, type for ggplot to work and randomly generated data from Bet-Binom distrribution
grid <- c(seq(0,16), seq(0,16))
Function <- c(rep("True", 17), rep("Gibbs", 17))
f_x_freq <- rbbinom((1:m), 16, 2, 4)

# Create data frame of all input for ggplot
Freq <- c(as.data.frame(table(f_x_freq))$Freq, rep(0, 17- length(unique(f_x_freq))), 
          c(as.data.frame(table(x_sample))$Freq,rep(0, 17- length(unique(x_sample)) )))
data <- as.data.frame(cbind(Function, grid, Freq))
data[,2] <- as.numeric(as.character(data[,2]))
data[,3] <- as.numeric(as.character(data[,3]))

# Run first plot of two histgrams frequency
ggplot(data, aes(x= grid, y= Freq, fill = Function)) +
  theme_classic() +
  geom_bar(stat="identity", position = position_dodge(width = 0.9), width = 0.5, colour="black") +
  scale_x_continuous(limits = c(-0.4,16), breaks = seq(0,16)) +
  scale_y_continuous(limits = c(0,70), breaks = seq(0,70,10)) +
  labs(title = "", x = "", y = "")

# Plot 3
# Density function 2.11
f_x_density <- function(x){
  (1/m)*sum(dbinom(x, 16, y_sample))
}

# Create grid, type for ggplot to work and generated estimated marginal density
grid <- c(seq(0,16), seq(0,16))
Function <- c(rep("True", 17), rep("Gibbs", 17))
x_density <- sapply(grid[1:17], f_x_density)

# Calculate true marginal density and store as data frame
Freq <- c(sapply(seq(0,16), f_x), c(x_density))
data <- as.data.frame(cbind(Function, grid, Freq))
data[,2] <- as.numeric(as.character(data[,2]))
data[,3] <- as.numeric(as.character(data[,3]))

# Run second plot of two histograms density
ggplot(data, aes(x= grid, y= Freq, fill = Function)) +
  theme_classic() +
  geom_bar(stat="identity", position = position_dodge(width = 0.9), width = 0.5, colour="black") +
  scale_x_continuous(limits = c(-0.4,16), breaks = seq(0,16)) +
  scale_y_continuous(limits = c(0,0.12), breaks = seq(0,0.12, length = 7)) +
  labs(title = "", x = "", y = "")

# Calculate the error
error_non_optimized <- mean(abs(data[1:17,3] - data[18:34, 3]))

# Plot 5
# Write all three conditional distributions, x,y and n
f_x.y.n <- function(n, y){rbinom(1, n, y)}
f_y.x.n <- function(x, n, a = 2, b = 4){rbeta(1, x + a, n - x + b)}
f_n.x.y <- function(x, y, lambda = 16){x + rpois(1, ((1-y)*lambda))}

# Gibbs sampler 2: With three variables
# Number of iterations
m <- 500
# Length of Gibbs sequence
k <- 10

# Create empty vectors
x <- rep(NA, k+1)
y <- rep(NA, k+1)
n <- rep(NA, k+1)
x_sample <- rep(NA, m)
y_sample <- rep(NA, m)
n_sample <- rep(NA, m)

# Run iterations
for(j in 1:m){
  # Generate starting values
  y[1] <- rbeta(1, 2, 4)
  n[1] <- rpois(1, 16)
  
  for(i in 2:(k+1)){
    x[i] <- f_x.y.n(n = n[i-1], y = y[i-1])
    y[i] <- f_y.x.n(x = x[i], n = n[i-1])
    n[i] <- f_n.x.y(x = x[i], y = y[i])
  }
  x_sample[j] <- x[k+1]
  y_sample[j] <- y[k+1]
  n_sample[j] <- n[k+1]
}

# Density function 2.11
f_x_density <- function(x){
  (1/m)*sum(dbinom(x, n_sample, y_sample))
}

# Create grid, and neccessary calulations for ggplot
grid <- 0:max(unique(x_sample))
Freq <- sapply(grid, f_x_density)
Function <- rep("Gibbs", length(grid))
data <- as.data.frame(cbind(Function, grid, Freq))
data[,2] <- as.numeric(as.character(data[,2]))
data[,3] <- as.numeric(as.character(data[,3]))

# Create third plot, two histograms of marginal density 
ggplot(data, aes(x= grid, y= Freq, fill = Function)) +
  theme_classic() +
  geom_bar(stat="identity", position = position_dodge(width = 0.9), width = 0.5, colour="black") +
  scale_x_continuous(limits = c(-1,max(grid)), breaks = seq(2,max(grid), by = 4)) +
  scale_y_continuous(limits = c(0,0.12), breaks = seq(0,0.12, length = 7)) +
  labs(title = "", x = "", y = "")

# Optimize with Coda
# Set up values according to the paper. 
# k is length of Gibbs sequence, m is number of sequences
alpha <- 2; beta <- 4; k <- 4000; n <- 16; m <- 10

# Set up gibbs sampler.
y_gibbs <- NULL
x_gibbs <- NULL
for (i in 1:m) {
  y_k <- rep(NA, k)
  x_k <- rep(NA, k)
  y_cond <- rbeta(1, alpha, beta)
  for (j in 1:k) {
    x_k[j] <- rbinom(1, 16, y_cond)
    y_k[j] <- rbeta(1, x_k[j] + alpha, 16 - x_k[j] + beta)
    y_cond <- y_k[j]
  }
  y_gibbs <- cbind(y_gibbs, y_k)
  x_gibbs <- cbind(x_gibbs, x_k)
}

# Make inte mcmc objects for coda package.
x_gibbs_mcmc1 <- mcmc(cbind(x_gibbs[,1], y_gibbs[,1]))
x_gibbs_mcmc2 <- mcmc(cbind(x_gibbs[,2], y_gibbs[,2]))
x_gibbs_mcmc3 <- mcmc(cbind(x_gibbs[,3], y_gibbs[,3]))
x_gibbs_mcmc4 <- mcmc(cbind(x_gibbs[,4], y_gibbs[,4]))
x_gibbs_mcmc5 <- mcmc(cbind(x_gibbs[,5], y_gibbs[,5]))
x_gibbs_mcmc6 <- mcmc(cbind(x_gibbs[,6], y_gibbs[,6]))
x_gibbs_mcmc7 <- mcmc(cbind(x_gibbs[,7], y_gibbs[,7]))
x_gibbs_mcmc8 <- mcmc(cbind(x_gibbs[,8], y_gibbs[,8]))
x_gibbs_mcmc9 <- mcmc(cbind(x_gibbs[,9], y_gibbs[,9]))
x_gibbs_mcmc10 <- mcmc(cbind(x_gibbs[,10], y_gibbs[,10]))

# Make into mcmc list.
x_gibbs_mcmc <- mcmc.list(x_gibbs_mcmc1, x_gibbs_mcmc2,x_gibbs_mcmc3, x_gibbs_mcmc4, x_gibbs_mcmc5, x_gibbs_mcmc6,
                         x_gibbs_mcmc7, x_gibbs_mcmc8, x_gibbs_mcmc9, x_gibbs_mcmc10)

summary(x_gibbs_mcmc) # Order: x, y
plot(x_gibbs_mcmc)

# Autocorrelation.
acf(x_gibbs_mcmc1) # For X variable only
round(autocorr.diag(x_gibbs_mcmc),4) # For all variables

# Sample size adjusted for autocorrelation.
effectiveSize(x_gibbs_mcmc)

# Looks good after 1500 iterations. Shrink factor = 1 which is good! See p.254
gelman.diag(x_gibbs_mcmc)
gelman.plot(x_gibbs_mcmc) # Order: x, y

# Looks good.
geweke.diag(x_gibbs_mcmc)
geweke.plot(x_gibbs_mcmc)

heidel.diag(x_gibbs_mcmc)

# Optimum solution of Gibbs sampler for plot 3
# True marginal distribution
f_x <- function(x, a=2, b=4, n=16){
  choose(n,x)*(gamma(a+b)/(gamma(a)*gamma(b)))*((gamma(x+a)*gamma(n-x+b))/(gamma(a+b+n)))
}

# Conditional distributions
f_x.y <- function(n=16, y){rbinom(1, n, y)}
f_y.x <- function(x, a=2, b=4, n=16){rbeta(1, x+a, n-x+b)}

# Gibbs sampler 1: With two variables
burn <- 1500 # Convergence speed
m <- 500 # Length of sequence
k <- 10 # Lag for shinning sequence
e <- m * k + burn # Total number of iterations

# Set up values according to the paper. 
alpha = 2; beta = 4; n = 16

# Create empty vectors
x <- rep(NA, e+1)
y <- rep(NA, e+1)
x_sample <- rep(NA, m)
y_sample <- rep(NA, m)

# Starting value
y[1] <- rbeta(1, alpha, beta)

# Run iterations
for(i in 2:(e+1)){
  x[i] <- f_x.y(y = y[i-1])
  y[i] <- f_y.x(x = x[i])
}

# Remove the first 1000 (1001) iterations since convergence after 1000 iterations
x <- x[-(1:(burn+1))]
y <- y[-(1:(burn+1))]

# Systematic sampling for every k_th iterations (thinning)
x_sample <- x[seq(sample(1:k, 1), length(x), by = k)]
y_sample <- y[seq(sample(1:k, 1), length(y), by = k)]

# Density function 2.11
f_x_density <- function(x){
  (1/m)*sum(dbinom(x, 16, y_sample))
}

# Create grid and neccessary calculations for ggplot
grid <- c(seq(0,16), seq(0,16))
Function <- c(rep("True", 17), rep("Gibbs", 17))
x_density <- sapply(grid[1:17], f_x_density)

Freq <- c(sapply(seq(0,16), f_x), c(x_density))
data <- as.data.frame(cbind(Function, grid, Freq))
data[,2] <- as.numeric(as.character(data[,2]))
data[,3] <- as.numeric(as.character(data[,3]))

# Create plot of optimized estimated marginal density
ggplot(data, aes(x= grid, y= Freq, fill = Function)) +
  theme_classic() +
  geom_bar(stat="identity", position = position_dodge(width = 0.9), width = 0.5, colour="black") +
  scale_x_continuous(limits = c(-0.4,16), breaks = seq(0,16)) +
  scale_y_continuous(limits = c(0,0.12), breaks = seq(0,0.12, length = 7)) +
  labs(title = "", x = "", y = "")

# Calculate error
error_optimized <- mean(abs(data[1:17,3] - data[18:34, 3]))

# Compare error between optimized and un-optimized
round(c(error_non_optimized, error_optimized),4)

###### 2 ######

# Data
Hits <- seq(0,5)
Observations <- c(229, 211, 93, 35, 7, 1)
data <- as.data.frame(cbind(Hits, Observations))

# Alpha and Beta value, 0->0.5 due to gamma restrictions
alpha_beta = c(0.5,10,100)

# Posterior distribution
posterior = function (lambda, alpha, beta, s, n) {
  dgamma(lambda, alpha+s, beta+n)
}

# Prior distribution
prior = function (x, alpha, beta) {
  dgamma(x, alpha, beta)
}

# Create grid
grid <- seq(0,1.4,0.001)

# Change R graph settings
par(mfrow=c(3,3), xpd=NA, cex=0.4, mar=c(4.5, 3, 3, 3) + 0.1)

# Run all 9 plots iteratively
for (alpha in alpha_beta) {
  for (beta in alpha_beta) {
    # Calculate posterior values over grid
    gridpost <- sapply(grid, posterior, alpha=alpha, beta=beta, s=535, n=576)
    # Calculate prior values over grid
    gridprior <- sapply(grid, prior, alpha=alpha, beta=beta)
    # Plot posterior
    plot(grid, gridpost, type="l", col="blue", axes=FALSE,
         main = bquote(alpha == .(alpha) ~ "," ~ beta == .(beta)), xlab= bquote(theta), ylab="")
    # Add posterior axis
    axis(2,col="blue",las=1)
    axis(1,col="black",las=1)
    # Allow to plot over plot
    par(new=TRUE)
    # Plot prior 
    plot(grid, gridprior, type="l", col="red", axes=FALSE, xlab="", ylab="")
    # Add prior axis
    axis(4,col="red",las=1)
  }
}
# Add legend to last plot
legend(0.1, 3.5, legend=c("Posterior", "Prior"), col=c("blue", "red"), lty=1, cex=0.8 )
# Reset 
dev.off()
par(mfrow=c(1,1))

# Credible intevals
# Sum of y
s <- 535
# Total number of districts
n <- 576
# Expected value calulation function
exp_post <- function(alpha, beta){
  (alpha + s)/(beta + n)
}
# SD value calulation function
sd_post <- function(alpha, beta){
  (alpha + s)^(1/2)/(beta + n)
}
# Combine Expected and SD to output a interval
cred_interval <- function(alpha, beta){
  c(exp_post(alpha, beta)-1.96 * sd_post(alpha, beta), exp_post(alpha, beta)+ 1.96 * sd_post(alpha, beta))
}

alpha_beta <- c(0,10,100)
# Approximate equal tail credible interval
interval <- c()
for (alpha in alpha_beta) {
  for (beta in alpha_beta) {
    cred <- cred_interval(alpha, beta)
    interval <- rbind(interval, c(alpha, beta, cred))
  }
}

interval_app <- cbind(interval, c(interval[,4]-interval[,3]))
interval_app

alpha_beta <- c(0,10,100)
# Exact equal tail credible interval
interval <- c()
for (alpha in alpha_beta) {
  for (beta in alpha_beta) {
    cred <- c(qgamma(0.025, shape = (alpha + 535), rate = (beta + 576)),
              qgamma(0.975, shape = (alpha + 535), rate = (beta + 576)))
    interval <- rbind(interval, c(alpha, beta, cred))
  }
}

interval_exac <- cbind(interval, c(interval[,4]-interval[,3]))
interval_exac


# Highest Probability Density credible interval
library(HDInterval)
alpha_beta <- c(0,10,100)
interval <- c()
for (alpha in alpha_beta) {
  for (beta in alpha_beta) {
    cred <- hdi(qgamma, credMass = 0.95, shape = (alpha + 535), rate = (beta + 576))
    interval <- rbind(interval, c(alpha, beta, cred))
  }
}

interval_hpd <- cbind(interval, c(interval[,4]-interval[,3]))
interval_hpd

# Combine all interval types for alpha/beta = 0
interval_all <- round(rbind(interval_app[1,], interval_exac[1,], interval_hpd[1,]),4)

# Plot intervals
par(mfrow=c(1,1))
# Grid over 0.8, to 1.05
grid <- seq(0.8, 1.05, by = 0.001)
# Calculate interval from smallest to largest interval
x_grid <- seq(max(interval_all[,3]), max(interval_all[,4]), by = 0.001)
# Add extra value on both sides of the interval
x_coord <- c(max(interval_all[,3]), x_grid, max(interval_all[,4]))
# Add 0 on both sides of the y calculations, to make the polygon go to the x-axis
y_coord <- c(0, sapply(x_grid, dgamma, shape = 535, rate = 576), 0)
# PLot the posterior function
plot(grid, sapply(grid, dgamma, shape = 535, rate = 576), type = "l", col = "black",
     xlab = expression(theta), ylab = expression(paste(pi, '(', theta, '|y)')))
# Add coloring for the covered area
polygon(x_coord, y_coord, col = 'pink', lwd = 2,  border = 'violet')
# Add interval lines
abline(v=interval_all[,3:4], col = c(2,3,4,2,3,4), lty = c(1,2,3,1,2,3))
# Add a legend to display what line represent what interval
legend("topright", c("Approx.", "Exact", "HPD"), col = c(2,3,4), lty = c(1,2,3), cex = 0.8)

###### 3 ######

# Posterior distribution: log of equation (6) (pi)
post <- function(rho){
  (-3/2)*log(1 - rho^2) - n*log((1 - rho^2)^(1/2)) - sum(1/(2*(1 - rho^2))*(x^2-2*rho*x*y+y^2))
}

# Generate data
n <- 1000 # Number of generated observations
# Generate random data from the given distribution in equation (3), use seed to reproduce
set.seed(102);data <- as.data.frame(rmvnorm(n, mean = c(0,0), sigma = matrix(c(1, 0.4, 0.4, 1), ncol = 2))) 
names(data) <- c("x", "y")
# Take data as seperate vectors
x <- data[ , 1]
y <- data[ , 2]

# Metropolis algorithm
mh_symmetric <- function(f = 1, m = 10000, r = 0, a = 0.07, b = 0.07, t = 0.07){
  # Create empty output list 
  output <- list()
  
  # Gibbs sampler
  burn <- r # Possiblility to remove the first "r" observations in the chain
  e <- burn + m # Number of iterations, inncluding the burn in 
  
  # Initialize the chain. 
  rho <- 0 # as if there's no correlation at all.
  
  # Store the samples
  chain_rho <- rep(NA, e-burn)
  
  accepted_number <- 0 # Initial value
  for(i in 1:e){
    
    # Draw a value from the proposal distribution.
    if(f == 1){
      # Uniform(rho-0.07,rho+0.07); Equation 7
      rho_candidate <- runif(1, rho - a, rho + b)
    }else{
      # N(rho, t)
      rho_candidate <- rnorm(1, rho, t)
    }
    
    # Compute the acceptance probability, Equation 8 and Equation 6. 
    # We will do both equations in log domain here to avoid underflow. 
    accept <- post(rho_candidate) - post(rho)
    accept <- exp(accept)
    
    # Handle possible NaN 
    if(is.nan(accept)){
      accept <- 0
    }else{
      accept <- min(1, accept)
    }
    
    # Accept rho_candidate with probability accept.
    if (runif(1)<accept){
      rho <- rho_candidate
    accepted_number <- accepted_number+1
    }else{
      rho <- rho
    }
    # store without burn in 
    if(i >= burn){
      chain_rho[i-(burn)] <- rho
    }
  }
  output[[1]] <- chain_rho # Store rho chain
  output[[2]] <- accepted_number # Store number of accepted candidates
  output[[3]] <- cov(x, y) # Store true rho
  output[[4]] <- e - burn # Store number of iterations
  return(output)

}

# Run the algorith with selected input: m = number of iterations, r = number of burn ins, 
# a/b value in Uniform, t sd in Normal
# Set f = 1 for Uniform proposal, and 2 for normal proposal
chain_rho <- mh_symmetric(f = 1)
print(paste0("Acceptance ratio is ", (chain_rho[[2]]/chain_rho[[4]])))
print(paste0("Mean rho is ", round(mean(chain_rho[[1]]), 4)))
print(paste0("Sd for rho is ", round(sd(chain_rho[[1]]), 4)))
print(paste0("Compare with cov function: ", round(chain_rho[[3]], 4)))

# Create data frame for ggplot
iteration <- seq(1, length(chain_rho[[1]]))
chain_rho_df <- as.data.frame(cbind(iteration, chain_rho[[1]]))
names(chain_rho_df)[2] <- "chain_rho"

# Function of plot, to plot multiple models
plot_function <- function(type = 1){
  # Scatter plot
  p_scatter <- ggplot(data = data, aes(x,y)) + 
    theme_classic() +
    theme(panel.border = element_rect(linetype = "solid", fill = NA), 
          plot.margin=unit(c(0,1,0,0),"cm")) +
    geom_point(colour = "black", fill = "blue", shape = 21, size = 5) +
    scale_x_continuous(expand = c(0, 0), limits = c(-4, 4), breaks = seq(-4, 4, by = 1)) +
    scale_y_continuous(expand = c(0, 0), limits = c(-4, 4), breaks = seq(-4, 4, by = 1)) +
    labs(title = "", x = "", y = "")
  
  # Convergence plot
  p_converg <- ggplot(data = chain_rho_df, aes(iteration, chain_rho)) +
    theme_classic() +
    theme(panel.border = element_rect(linetype = "solid", fill = NA), 
          plot.margin=unit(c(0,1,0,0),"cm")) +
    geom_line(colour = "blue") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(iteration)), breaks = seq(0, max(iteration), by = 2000)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.6), breaks = seq(0, 0.6, by = 0.1)) +
    labs(title = "", x = "", y = expression(rho))
  
  # Have fixed x-axis or automatically adjustable
  limit_upper <- ifelse(type == 1, 0.6, max(chain_rho_df$chain_rho))
  limit_lower <- ifelse(type == 1, 0, min(chain_rho_df$chain_rho))
  # Histogram of rho form the chain 
  p_hist <- ggplot(data = chain_rho_df, aes(x = chain_rho)) +
    theme_classic() +
    theme(panel.border = element_rect(linetype = "solid", fill = NA), 
          plot.margin=unit(c(0,1,0,0),"cm")) +
    geom_histogram(colour = "black", fill = "blue", bins = 100) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1800), breaks = seq(0, 1800, by = 200)) +
    scale_x_continuous(expand = c(0, 0), limits = c(limit_lower, limit_upper), breaks = seq(limit_lower, limit_upper, by = 0.1)) +
    labs(title = "", x = expression(rho), y = "")
  
  # Three stacked plots vertically
  grid.arrange(p_scatter, p_converg, p_hist, nrow = 3)
}


# Type = 1 for fixed x-axis, else use 2
plot_function(1)

# Testing different uniform functions, equation (7)
vector <- seq(0, 0.5, by = 0.01) # Different values
# Function to output for the input value given 
loop_f <- function(i, func = 1){
  input <- vector[i]
  chain_rho <- mh_symmetric(f = func, a = input, b = input, t = input)
  output <- c((chain_rho[[2]]/chain_rho[[4]]),
    round(mean(chain_rho[[1]]), 4),
    round(sd(chain_rho[[1]]), 4),
    abs(round(chain_rho[[3]], 8)-round(mean(chain_rho[[1]]), 8)),
    round(i))
  names(output) <- c("Accaptance ratio", "Mean rho", "Sd rho", "Absolute error", "Iteration")
  
  return(output)
  
}

# Run algorithm for all values in vector. 
output <-  sapply(1:length(vector), loop_f, func = 1)
output_df <- as.data.frame(cbind(t(output), vector))
# Filter to see only those outputs with close target acceptance rate
output_df <- output_df %>%
  filter(`Accaptance ratio`<0.6, `Accaptance ratio`>0.4)
# From literature, out goal, for a 1-2 dimensional algorithm, is an acceptance rate of ~0.5. 

# Write data as latex code
write_matex(output_df)

# Asymmetrical proposal, Skew Normal Distribution
mh_a_symmetric <- function(m = 10000, r = 0, skew = 0.6, sigma_one = 0.04){
  output <- list() # Output multiple values / vectors
  
  # Gibbs sampler
  burn <- r # Burn in 
  e <- burn + m # Number of iterations including burn in 
  
  # Initialize the chain. 
  rho <- 0 # as if there's no correlation at all.
  
  # Store the samples
  chain_rho <- rep(NA, e-burn)
  
  accepted_number <- 0 # Initial value
  for(i in 1:e){
   
    # Create the proposal, here with non-symmetric distribution
    rho_candidate = rsnorm(1, rho, sigma_one, xi = skew)
    
    # Derive the posterior part. We will do both equations in log domain here to avoid underflow.
    accept <- post(rho_candidate)-post(rho)
    accept <- exp(accept)
    
    # Derive the full acceptance function
    accept <- accept/dsnorm(rho, rho_candidate, 
                           sigma_one, xi = skew)*dsnorm(rho_candidate, rho, sigma_one, xi = skew)
    # Handle NaN
    if(is.nan(accept)){
      accept <- 0
    }else{
    accept <- min(1, accept)
    }
    
    # Accept rho_candidate with probability accept.
    if (runif(1)<accept){
      rho <- rho_candidate
      accepted_number <- accepted_number+1
    }else{
      rho <- rho
    }
    # store whitout burn in 
    if(i >= burn){
      chain_rho[i-(burn)] <- rho
    }
  }
  output[[1]] <- chain_rho # Store rho chain
  output[[2]] <- accepted_number # Store number of accepted candidates
  output[[3]] <- cov(x, y) # Store true rho
  output[[4]] <- e - burn # Store number of iterations
  return(output)
  
}

# Run algorithm, slect skew parameter and sigma_one (sd)
chain_rho <- mh_a_symmetric(skew = 0.44, sigma_one = 0.05) # low skew converge to small values, large skew converge to large values
print(paste0("Acceptance ratio is ", (chain_rho[[2]]/chain_rho[[4]])))
print(paste0("Mean rho is ", round(mean(chain_rho[[1]]), 4)))
print(paste0("Sd for rho is ", round(sd(chain_rho[[1]]), 4)))
print(paste0("Compare with cov function: ", round(chain_rho[[3]], 4)))

# Create data frame to run plot function
iteration <- seq(1, length(chain_rho[[1]]))
chain_rho_df <- as.data.frame(cbind(iteration, chain_rho[[1]]))
names(chain_rho_df)[2] <- "chain_rho"

# Plot non-symmetric proposal
plot_function(1)

# Testing different skew values, sd is fixed at best value from the normal test (0.05)
vector <- seq(0, 2, by = 0.01) # Different values,
loop_f <- function(i){
  input <- vector[i]
  chain_rho <- mh_a_symmetric(skew = input, sigma_one = 0.05)
  output <- c((chain_rho[[2]]/chain_rho[[4]]),
              round(mean(chain_rho[[1]]), 4),
              round(sd(chain_rho[[1]]), 4),
              abs(round(chain_rho[[3]], 8)-round(mean(chain_rho[[1]]), 8)),
              round(i))
  names(output) <- c("Accaptance ratio", "Mean rho", "Sd rho", "Absolute error", "Iteration")
  
  return(output)
  
}

# Run algorithm for all values in vector. 
output <-  sapply(1:length(vector), loop_f)
output_df <- as.data.frame(cbind(t(output), vector))
# Filter those values with close target acceptance rate. 
output_df <- output_df %>%
  filter(`Accaptance ratio`<0.5005, `Accaptance ratio`>0.4995)

# Write to latex code
write_matex(output_df)

# Plot the optimal Skewed Normal function from running the test
plot(seq(0.3,0.6, by = 0.001), sapply(seq(0.3,0.6, by = 0.001), 
                 dsnorm, mean = 0.4, sd = 0.05, xi = 1.44), type = "l", 
     xlab = expression(rho), ylab = expression(paste('SN(',mu,",", sigma,",", alpha, ')')), col = "blue")

# Three streams of observations

# Posterior function given Jeffreys prior and 3-dim log-likelihood normal distribution
post <- function(rho1, rho2, rho3, data_matrix, n){
  sigma <- matrix(c(1, rho1, rho2, rho1, 1, rho3, rho2, rho3, 1), ncol = 3)
  sigma_solve <- solve(sigma)
  log_prior <- log(1/(det(sigma)^(4/2)))
  log_likelihood <- (-n/2)*log(det(sigma))-(1/2)*sum(apply(data_matrix, 1, function(row) t(row)%*%sigma_solve%*%(row)))
  output <- log_prior + log_likelihood
  return(output)
}

# Generate data
n <- 1000 # Number of generated observations
set.seed(102);data <- as.data.frame(rmvnorm(n, mean = c(0,0,0), sigma = matrix(c(1, 0.4, 0.4, 0.4, 1, 0.4, 0.4, 0.4, 1), ncol = 3))) # Generate random data from the given distribution in equation (3)
names(data) <- c("z1", "z2", "z3")
z1 <- data[ , 1]
z2 <- data[ , 2]
z3 <- data[ , 3]

# Change into matrix
data_matrix <- as.matrix(data)

# Metropolis algorithm
mh_3var_symmetric <- function(f = 1, m = 1000, r = 0, a = 0.07, b = 0.07, t = 0.07){
  # Empty output list
  output1 <- list()
  
  # Metropolis sampler
  burn <- r # Possiblility to remove the first "r" observations in the chain
  e <- burn + m # Number of iterations, inncluding the burn in 
  
  # Initialize the chain. 
  rho1 <- 0 # as if there's no correlation at all.
  rho2 <- rho1
  rho3 <- rho2
  
  # Store the samples
  chain_rho1 <- rep(NA, e-burn)
  chain_rho2 <- rep(NA, e-burn)
  chain_rho3 <- rep(NA, e-burn)
  
  accepted_number1 <- 0 # Initial value
  
  for(i in 1:e){
    
    # Draw a value from the proposal distribution.
    if(f == 1){
      # Uniform(rho-0.07,rho+0.07); Equation 7
      rho_candidate1 <- runif(1, rho1 - a, rho1 + b)
      rho_candidate2 <- runif(1, rho2 - a, rho2 + b)
      rho_candidate3 <- runif(1, rho3 - a, rho3 + b)
    }else{
      # Normal(rho, t)
      rho_candidate1 <- rnorm(1, rho1, t)
      rho_candidate2 <- rnorm(1, rho2, t)
      rho_candidate3 <- rnorm(1, rho3, t)
    }
    
    # Log-posterior
    accept1 <- exp(post(rho_candidate1, rho_candidate2,rho_candidate3,data_matrix,1000) - post(rho1, rho2, rho3,data_matrix,1000))
    
    # Handle possible NaN
    if(is.nan(accept1)){
      accept1 <- 0
    }else{
      accept1 <- min(1, accept1)
    }
    
    # Accept rho_candidate with probability accept.
    if(runif(1)<accept1){
      rho1 <- rho_candidate1
      rho2 <- rho_candidate2
      rho3 <- rho_candidate3
      accepted_number1 <- accepted_number1+1
    }else{
      rho1 <- rho1
      rho2 <- rho2
      rho3 <- rho3
    }
    
    # store without burn in 
    if(i >= burn){
      chain_rho1[i-(burn)] <- rho1
    }
    if(i >= burn){
      chain_rho2[i-(burn)] <- rho2
    }
    if(i >= burn){
      chain_rho3[i-(burn)] <- rho3
    }
  }
  
  output1[[1]] <- chain_rho1 # Store rho1 chain
  output1[[2]] <- accepted_number1 # Store number of accepted candidates
  output1[[3]] <- cov(z1, z2) # Store true rho1
  output1[[4]] <- e - burn # Store number of iterations
  
  output1[[5]] <- cov(z1, z3) # Store true rho2 
  output1[[6]] <- cov(z2, z3) # Store true rho3
  output1[[7]] <- chain_rho2 # Store rho2 chain
  output1[[8]] <- chain_rho3 # Store rho3 chain
  
  return(output1)
}

# f = 1 for Uniform, 2 for Normal. m = number of iterations, r = burn in
# a/b value in Uniform, t = sd in Normal
chain_rho <- mh_3var_symmetric(f = 1)

# Read out only the important stuff function
readout <- function(){
  print(paste0("Acceptance rate is ", chain_rho[[2]]/chain_rho[[4]]))
  print(paste0("rho1 is ", round(mean(chain_rho[[1]]),4), " vs true ", round(chain_rho[[3]], 4)))
  print(paste0("rho2 is ", round(mean(chain_rho[[7]]),4), " vs true ", round(chain_rho[[5]], 4)))
  print(paste0("rho3 is ", round(mean(chain_rho[[8]]),4), " vs true ", round(chain_rho[[6]], 4)))
}
# Run read out function
readout()

# Plot all three rho parameters trace plot
# Added blue line to represent the true rho for each parameter
par(mfrow=c(3,1))
plot(chain_rho[[1]], type = "l", xlab = "Iteration", ylab = expression(paste(rho,'1')))
abline(h=chain_rho[[3]], col = "blue")
plot(chain_rho[[7]], type = "l", xlab = "Iteration", ylab = expression(paste(rho,'2')))
abline(h=chain_rho[[5]], col = "blue")
plot(chain_rho[[8]], type = "l", xlab = "Iteration", ylab = expression(paste(rho,'3')))
abline(h=chain_rho[[6]], col = "blue")
par(mfrow=c(1,1))

# Try different input values for 3 var model. 
# NOTE: loop takes quite long time.
vector <- seq(0, 0.3, by = 0.01) # Different values, 
loop_f <- function(i, func = 1){
  input <- vector[i]
  chain_rho <- mh_3var_symmetric(f = func, a = input, b = input, t = input)
  output <- c(chain_rho[[2]]/chain_rho[[4]],
              
              round(mean(chain_rho[[1]]),4),
              round(chain_rho[[3]], 4),
              
              round(mean(chain_rho[[7]]),4),
              round(chain_rho[[5]], 4),
              
              round(mean(chain_rho[[8]]),4),
              round(chain_rho[[6]], 4),
              round(i))
  names(output) <- c("Accaptance ratio", "Mean rho1", "True rho1", "Mean rho2", "True rho2","Mean rho3", "True rho3","Iteration")
  
  return(output)
  
}

# NOTE: loop takes quite long time.
# Run algorithm for all values in vector. 
output <-  sapply(1:length(vector), loop_f, func = 1)
output_df <- as.data.frame(cbind(t(output), vector))
# Filter for values with target aceptance rate.
output_df <- output_df %>%
  filter(`Accaptance ratio` < 0.3, `Accaptance ratio`>0.15)

# Write to latex code
write_matex(output_df)

############
