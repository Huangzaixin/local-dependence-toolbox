###################################################################################################
# Example 3: Simulation data
###################################################################################################
# Author: Zaixin Huang
# Email: eric.huangzaixin@gmail.com
# This version 2024.11.05 
# The latest version can be downloaded from https://github.com/Huangzaixin/local-dependence-toolbox
###################################################################################################
install.packages("ggplot2") 
install.packages("fitdistrplus") 
install.packages("stats") 
install.packages("MASS") 
install.packages("ggExtra")

library(ggplot2)
library(fitdistrplus) 
library(stats)
library(MASS) 
require(copula)
library(ggExtra)

#######################  generate data from the mixture copula model ######################## 
randomNumber = 1000; 

# MixFG copula : Frank copula + Gumbel copula
mC <- mixCopula(list(frankCopula(1, dim=2), gumbelCopula(2, dim=2)), c(5,5)/10)
rndcopula <- rCopula(randomNumber, mC)

u <- rndcopula[,1]
v <- rndcopula[,2]

# save to csv file
write.csv(u, file = "data/simulated_u.csv", row.names = FALSE)
write.csv(v, file = "data/simulated_v.csv", row.names = FALSE)

u <- read.csv('data/simulated_u.csv', header = T)
v <- read.csv('data/simulated_v.csv', header = T)

# transform u, v into X, Y
x <- qgamma(u[,1], shape = 3, rate = 0.5)
y <- qnorm(v[,1], mean = 5, sd = 1.4)
xy <- data.frame(x = x, y = y)

#######################  scatter plot ######################
# region 1
reg1_x_min <- 0
reg1_x_max <- 2
reg1_y_min <- 0
reg1_y_max <- 3

# region 2
reg2_x_min <- 0
reg2_x_max <- 3.5
reg2_y_min <- 0
reg2_y_max <- 4.05

# region 3
reg3_x_min <- 2
reg3_x_max <- 9.8
reg3_y_min <- 3
reg3_y_max <- 7.2

# region 4
reg4_x_min <- 3.5
reg4_x_max <- 8.1
reg4_y_min <- 4.05
reg4_y_max <- 6.1

# region 5
reg5_x_min <- 8.1
reg5_x_max <- 22
reg5_y_min <- 6.1
reg5_y_max <- 10

line_5_1 <- data.frame( x = c(reg5_x_min, reg5_x_max), y = c(reg5_y_min, reg5_y_min))
line_5_2 <- data.frame( x = c(reg5_x_min, reg5_x_min), y = c(reg5_y_min, reg5_y_max))
line_5_3 <- data.frame( x = c(reg5_x_max, reg5_x_max), y = c(reg5_y_min, reg5_y_max))
line_5_4 <- data.frame( x = c(reg5_x_min, reg5_x_max), y = c(reg5_y_max, reg5_y_max))

# region 6
reg6_x_min <- 12.5
reg6_x_max <- 22
reg6_y_min <- 7.2
reg6_y_max <- 10

line_6_1 <- data.frame( x = c(reg6_x_min, reg6_x_max), y = c(reg6_y_min, reg6_y_min))
line_6_2 <- data.frame( x = c(reg6_x_min, reg6_x_min), y = c(reg6_y_min, reg6_y_max))
line_6_3 <- data.frame( x = c(reg6_x_max, reg6_x_max), y = c(reg6_y_min, reg6_y_max))
line_6_4 <- data.frame( x = c(reg6_x_min, reg6_x_max), y = c(reg6_y_max, reg6_y_max))

p <- ggplot(xy, aes(x = x, y = y)) + 
        geom_point(color = "black", fill = "lightblue", size = 2, shape = 21, stroke = 0.5) +
        # region 1
        geom_rect(aes(xmin = reg1_x_min, xmax = reg1_x_max, ymin = reg1_y_min, ymax = reg1_y_max),
                  color = "#EBB842", size = 0.5, fill = "NA", alpha = .02) +
        # region 2
        geom_rect(aes(xmin = reg2_x_min, xmax = reg2_x_max, ymin = reg2_y_min, ymax = reg2_y_max),
                  color = "#EBB842", size = 0.5, fill = "NA", alpha = .02) +
        # region 3
        geom_rect(aes(xmin = reg3_x_min, xmax = reg3_x_max, ymin = reg3_y_min, ymax = reg3_y_max),
                  color = "#ff5722", size = 0.5, fill = "NA", alpha = .02) +
        # region 4
        geom_rect(aes(xmin = reg4_x_min, xmax = reg4_x_max, ymin = reg4_y_min, ymax = reg4_y_max),
                  color = "#ff5722", size = 0.5, fill = "NA", alpha = .02) +
        # region 5 
        geom_segment(data = line_5_1, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), color = "#EBB842", size = 0.5) +
        geom_segment(data = line_5_2, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), color = "#EBB842", size = 0.5) +
        geom_segment(data = line_5_3, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), linetype = "longdash", color = "#EBB842", size = 0.5) +
        geom_segment(data = line_5_4, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), linetype = "longdash", color = "#EBB842", size = 0.5) +
        
        # region 6    
        geom_segment(data = line_6_1, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), color = "#EBB842", size = 0.5) +
        geom_segment(data = line_6_2, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), color = "#EBB842", size = 0.5) +
       
        labs(title = "", x = "X", y = "Y") +
        xlim(c(0, 24)) + ylim(c(0, 10)) + theme(aspect.ratio = 1) + 
        scale_x_continuous(breaks = seq(0, 24, by = 6)) +
        scale_y_continuous(breaks = seq(0, 10, by = 2.5)) + 
        theme(plot.title = element_text(hjust = 0.5)) + theme(plot.margin = unit(c(0, 0, 0, 0), "mm")) +
        theme(axis.text.x = element_text(size = 7.2), axis.text.y = element_text(size = 7.2)) +
        theme(axis.title=element_text(size = 8.5))

# densigram
p <- ggMarginal(p, type = "densigram", xparams = list(color = "#222", fill = "lightblue", size = 0.35, bins = 40), 
                yparams = list(color = "#222", fill = "lightblue", size = 0.35, bins = 40)) 

print(p)

ggsave("scatter_example_3.eps", plot=p, device="eps", width = 3.4, height = 3.4, units = "in")


###################### fit different distribution functions ###################### 
# X
# gamma
gamma_x <- fitdist(x, "gamma", "mle") 
gamma_x$loglik   # -2533.588, the largest log-likelihood value
# weibull
weibull_x <- fitdist(x, "weibull", "mle")  
weibull_x$loglik   # -2541.171
# lognormal
lognormal_x <- fitdist(x, "lnorm", "mle")  
lognormal_x$loglik  # -2578.149
# logisitic
logis_x <- fitdist(x, "logis", "mle")  
logis_x$loglik   # -2618.02

# Y
# normal
norm_y <- fitdistr(y, "normal")
norm_y$loglik   # -1735.25, the largest log-likelihood value
# t
t_y <- fitdistr(y, "t")
t_y$loglik   # -1735.57


#######################  histogram and probability density curve  ###################### 
# X
hist_x <- hist(x, breaks = 20, freq = FALSE, xlim = c(0,20), ylim = c(0,0.15),
                  main = "Histogram of X", xlab = "X", ylab = "density", col = "lightblue", border = "black")
index <- seq(0, 20, by = 0.01)
lines(index,dgamma(index, shape = gamma_x$estimate[1], rate = gamma_x$estimate[2]), type = "l", lty = 1, col = "red", lwd = 4)

# compute the 95% quantile of the gamma distribution
qgamma(0.95, shape = gamma_x$estimate[1], rate = gamma_x$estimate[2])

# compute the cumulative probability of the gamma distribution at a given value
quantile_reg1_x_min <- pgamma(reg1_x_min, shape = gamma_x$estimate[1], rate = gamma_x$estimate[2])
quantile_reg1_x_max <- pgamma(reg1_x_max, shape = gamma_x$estimate[1], rate = gamma_x$estimate[2])
quantile_reg2_x_min <- pgamma(reg2_x_min, shape = gamma_x$estimate[1], rate = gamma_x$estimate[2])
quantile_reg2_x_max <- pgamma(reg2_x_max, shape = gamma_x$estimate[1], rate = gamma_x$estimate[2])
quantile_reg3_x_min <- pgamma(reg3_x_min, shape = gamma_x$estimate[1], rate = gamma_x$estimate[2])
quantile_reg3_x_max <- pgamma(reg3_x_max, shape = gamma_x$estimate[1], rate = gamma_x$estimate[2])
quantile_reg4_x_min <- pgamma(reg4_x_min, shape = gamma_x$estimate[1], rate = gamma_x$estimate[2])
quantile_reg4_x_max <- pgamma(reg4_x_max, shape = gamma_x$estimate[1], rate = gamma_x$estimate[2])
quantile_reg5_x_min <- pgamma(reg5_x_min, shape = gamma_x$estimate[1], rate = gamma_x$estimate[2])
quantile_reg5_x_max <- 1
quantile_reg6_x_min <- pgamma(reg6_x_min, shape = gamma_x$estimate[1], rate = gamma_x$estimate[2])
quantile_reg6_x_max <- 1

# Y
hist_y <- hist(y, breaks = 10, freq = FALSE, xlim = c(0,15), ylim = c(0,0.3), 
                  main = "Histogram of Y", xlab = "Y", ylab = "density", col = "lightblue", border = "black")
index <- seq(0, 15, by = 0.01)
lines(index,dnorm(index, mean = norm_y$estimate[1], sd = norm_y$estimate[2]), type = "l", lty = 1, col = "red", lwd = 4)

# compute the 95% quantile of the normal distribution
qnorm(0.95, mean = norm_y$estimate[1], sd = norm_y$estimate[2])

# compute the cumulative probability of the normal distribution at a given value
quantile_reg1_y_min <- pnorm(reg1_y_min, mean = norm_y$estimate[1], sd = norm_y$estimate[2])
quantile_reg1_y_max <- pnorm(reg1_y_max, mean = norm_y$estimate[1], sd = norm_y$estimate[2])
quantile_reg2_y_min <- pnorm(reg2_y_min, mean = norm_y$estimate[1], sd = norm_y$estimate[2])
quantile_reg2_y_max <- pnorm(reg2_y_max, mean = norm_y$estimate[1], sd = norm_y$estimate[2])
quantile_reg3_y_min <- pnorm(reg3_y_min, mean = norm_y$estimate[1], sd = norm_y$estimate[2])
quantile_reg3_y_max <- pnorm(reg3_y_max, mean = norm_y$estimate[1], sd = norm_y$estimate[2])
quantile_reg4_y_min <- pnorm(reg4_y_min, mean = norm_y$estimate[1], sd = norm_y$estimate[2])
quantile_reg4_y_max <- pnorm(reg4_y_max, mean = norm_y$estimate[1], sd = norm_y$estimate[2])
quantile_reg5_y_min <- pnorm(reg5_y_min, mean = norm_y$estimate[1], sd = norm_y$estimate[2])
quantile_reg5_y_max <- 1
quantile_reg6_y_min <- pnorm(reg6_y_min, mean = norm_y$estimate[1], sd = norm_y$estimate[2])
quantile_reg6_y_max <- 1

##################### probability integral transformation #####################
cdf_x <- pgamma(x, shape = gamma_x$estimate[1], rate = gamma_x$estimate[2])
cdf_y <- pnorm(y, mean = norm_y$estimate[1], sd = norm_y$estimate[2])
uv_x_y <- data.frame(x = cdf_x, y = cdf_y)

##################### scatter plot of cdf_x and cdf_y ##################
p <- ggplot(uv_x_y, aes(x = x, y = y)) + 
        geom_point(color = "black", fill = "lightblue", size = 2, shape = 21, stroke = 0.5) + 
        # region 1
        geom_rect(aes(xmin = quantile_reg1_x_min, xmax = quantile_reg1_x_max, ymin = quantile_reg1_y_min, ymax = quantile_reg1_y_max),
                  color = "#EBB842", size = 0.5, fill = "NA", alpha = .02) +
        # region 2
        geom_rect(aes(xmin = quantile_reg2_x_min, xmax = quantile_reg2_x_max, ymin = quantile_reg2_y_min, ymax = quantile_reg2_y_max),
                  color = "#EBB842", size = 0.5, fill = "NA", alpha = .02) + 
        # region 3
        geom_rect(aes(xmin = quantile_reg3_x_min, xmax = quantile_reg3_x_max, ymin = quantile_reg3_y_min, ymax = quantile_reg3_y_max),
                  color = "#ff5722", size = 0.5, fill = "NA", alpha = .02) + 
        # region 4
        geom_rect(aes(xmin = quantile_reg4_x_min, xmax = quantile_reg4_x_max, ymin = quantile_reg4_y_min, ymax = quantile_reg4_y_max),
                  color = "#ff5722", size = 0.5, fill = "NA", alpha = .02) + 
        # region 5
        geom_rect(aes(xmin = quantile_reg5_x_min, xmax = 1, ymin = quantile_reg5_y_min, ymax = 1),
                  color = "#EBB842", size = 0.5, fill = "NA", alpha = .02) +
        # region 6
        geom_rect(aes(xmin = quantile_reg6_x_min, xmax = 1, ymin = quantile_reg6_y_min, ymax = 1),
                  color = "#EBB842", size = 0.5, fill = "NA", alpha = .02) +
        labs(title = "", x = expression(u), y = expression(v)) + 
        xlim(c(0, 1)) + ylim(c(0, 1)) + theme(aspect.ratio = 1) +
        theme(plot.title = element_text(hjust = 0.5)) + theme(plot.margin = unit(c(0, 0, 0, 0), "mm")) +
        theme(axis.text.x = element_text(size = 8.43), axis.text.y = element_text(size = 8.43)) +
        theme(axis.title = element_text(size = 9.2))
print(p)

ggsave("scatter_uv_example_3.eps", plot=p, device="eps", width = 3.1, height = 3.1, units = "in")

# save to csv file
write.csv(cdf_x, file = "data/u.csv", row.names = FALSE)
write.csv(cdf_y, file = "data/v.csv", row.names = FALSE)



 