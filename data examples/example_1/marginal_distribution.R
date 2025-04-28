###################################################################################################
# Example 1: Parkinson's disease data
###################################################################################################
# Author: Zaixin Huang
# Email: eric.huangzaixin@gmail.com
# This version 2025.03.07 
# The latest version can be downloaded from https://github.com/Huangzaixin/local-dependence-toolbox
###################################################################################################
install.packages("fitdistrplus") 
install.packages("univariateML")
install.packages("ggplot2")  
install.packages("ggExtra")
install.packages("VGAM")  
install.packages("MASS") 
install.packages('latex2exp')

library(fitdistrplus) 
library(ggplot2) 
library(ggExtra)
library(univariateML)   # mlinvgauss
library(VGAM) 
library(MASS)   # fitdistr, student's t distribution
library(latex2exp)

####################### read data ###################### 
pd_dataset <- read.csv("data/FAM171A2_alpha_syn.csv",header = T)

FAM171A2_data <- pd_dataset$FAM171A2   # CSF measured FAM171A2 level
alpha_syn_data <- pd_dataset$alpha.syn   # CSF measured total α-syn

####################### scatter plot ####################
reg1_FAM171A2_min <- 7.15
reg1_FAM171A2_max <- 7.6
reg1_alpha_syn_min <- 1630
reg1_alpha_syn_max <- 5500
line_1_1 <- data.frame( x = c(reg1_FAM171A2_min, reg1_FAM171A2_min), y = c(reg1_alpha_syn_min, reg1_alpha_syn_max))
line_1_2 <- data.frame( x = c(reg1_FAM171A2_min, reg1_FAM171A2_max), y = c(reg1_alpha_syn_max, reg1_alpha_syn_max))
line_1_3 <- data.frame( x = c(reg1_FAM171A2_max, reg1_FAM171A2_max), y = c(reg1_alpha_syn_min, reg1_alpha_syn_max))
line_1_4 <- data.frame( x = c(reg1_FAM171A2_min, reg1_FAM171A2_max), y = c(reg1_alpha_syn_min, reg1_alpha_syn_min))

reg2_FAM171A2_min <- 7.3
reg2_FAM171A2_max <- 7.65
reg2_alpha_syn_min <- 1380
reg2_alpha_syn_max <- 3100

reg3_FAM171A2_min <- 7.47
reg3_FAM171A2_max <- 7.84
reg3_alpha_syn_min <- 810
reg3_alpha_syn_max <- 2400

reg4_FAM171A2_min <- 7.56
reg4_FAM171A2_max <- 7.92
reg4_alpha_syn_min <- 700
reg4_alpha_syn_max <- 1830

reg5_FAM171A2_min <- 7.65
reg5_FAM171A2_max <- 8.04
reg5_alpha_syn_min <- 585
reg5_alpha_syn_max <- 1380
line_5_1 <- data.frame( x = c(reg5_FAM171A2_min, reg5_FAM171A2_min), y = c(reg5_alpha_syn_min, reg5_alpha_syn_max))
line_5_2 <- data.frame( x = c(reg5_FAM171A2_min, reg5_FAM171A2_max), y = c(reg5_alpha_syn_max, reg5_alpha_syn_max))
line_5_3 <- data.frame( x = c(reg5_FAM171A2_max, reg5_FAM171A2_max), y = c(reg5_alpha_syn_min, reg5_alpha_syn_max))
line_5_4 <- data.frame( x = c(reg5_FAM171A2_min, reg5_FAM171A2_max), y = c(reg5_alpha_syn_min, reg5_alpha_syn_min))

reg6_FAM171A2_min <- 7.70
reg6_FAM171A2_max <- 9.55
reg6_alpha_syn_min <- 200
reg6_alpha_syn_max <- 1175
line_6_1 <- data.frame( x = c(reg6_FAM171A2_min, reg6_FAM171A2_min), y = c(reg6_alpha_syn_min, reg6_alpha_syn_max))
line_6_2 <- data.frame( x = c(reg6_FAM171A2_min, reg6_FAM171A2_max), y = c(reg6_alpha_syn_max, reg6_alpha_syn_max))
line_6_3 <- data.frame( x = c(reg6_FAM171A2_max, reg6_FAM171A2_max), y = c(reg6_alpha_syn_min, reg6_alpha_syn_max))
line_6_4 <- data.frame( x = c(reg6_FAM171A2_min, reg6_FAM171A2_max), y = c(reg6_alpha_syn_min, reg6_alpha_syn_min))

p <- ggplot(pd_dataset, aes(x = FAM171A2_data, y = alpha_syn_data)) + 
    geom_point(color = "white", fill = "#1eb3b9", size = 2.1, shape = 21, stroke = 0.5) +
    # region 1
    geom_segment(data = line_1_1, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), linetype = "longdash", color = "#ea6458", size = 0.5) +
    geom_segment(data = line_1_2, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), linetype = "longdash", color = "#ea6458", size = 0.5) +
    geom_segment(data = line_1_3, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), color = "#ea6458", size = 0.5) +
    geom_segment(data = line_1_4, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), color = "#ea6458", size = 0.5) +
    # region 2
    geom_rect(aes(xmin = reg2_FAM171A2_min, xmax = reg2_FAM171A2_max, ymin = reg2_alpha_syn_min, ymax = reg2_alpha_syn_max),
              color = "#f4ad60", size = 0.5, fill = "NA", alpha = .02) +
    # region 3
    geom_rect(aes(xmin = reg3_FAM171A2_min, xmax = reg3_FAM171A2_max, ymin = reg3_alpha_syn_min, ymax = reg3_alpha_syn_max),
              color = "#ea6458", size = 0.5, fill = "NA", alpha = .02) +
    # region 4
    geom_rect(aes(xmin = reg4_FAM171A2_min, xmax = reg4_FAM171A2_max, ymin = reg4_alpha_syn_min, ymax = reg4_alpha_syn_max),
              color = "#f4ad60", size = 0.5, fill = "NA", alpha = .02) +
    # region 5
    geom_rect(aes(xmin = reg5_FAM171A2_min, xmax = reg5_FAM171A2_max, ymin = reg5_alpha_syn_min, ymax = reg5_alpha_syn_max),
              color = "#ea6458", size = 0.5, fill = "NA", alpha = .02) +
    # region 6
    geom_segment(data = line_6_1, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), color = "#f4ad60", size = 0.5) +
    geom_segment(data = line_6_2, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), color = "#f4ad60", size = 0.5) +
    geom_segment(data = line_6_3, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), linetype = "longdash", color = "#f4ad60", size = 0.5) +
    geom_segment(data = line_6_4, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), linetype = "longdash", color = "#f4ad60", size = 0.5) +
    labs(title = "", x = "CSF measured FAM171A2 level (RFU)", y = TeX("CSF measured total $\\alpha$-syn (pg/ml)")) +
    xlim(c(7.15, 9.55)) + ylim(c(200, 5500)) + theme(aspect.ratio = 1) + 
    scale_x_continuous(breaks = seq(7.15, 9.55, by = 0.6)) +
    scale_y_continuous(breaks = seq(200, 5500, by = 1325)) + 
    theme(plot.title = element_text(hjust = 0.5)) + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm")) +
    theme(axis.text.x = element_text(size = 7.2), axis.text.y = element_text(size = 7.2)) +
    theme(axis.title = element_text(size = 7.4))

# densigram
p <- ggMarginal(p, type = "densigram", xparams = list(color = "white", fill = "#1eb3b9", size = 0.35, bins = 40), 
                yparams = list(color = "white", fill = "#1eb3b9", size = 0.35, bins = 40)) 

print(p)

ggsave("scatter_example_1.eps", plot=p, device="eps", width = 3.2, height = 3.2, units = "in")


###################### different distribution functions ###################### 
# FAM171A2 level
# student's t
t_FAM171A2 <- fitdistr(FAM171A2_data, "t")
t_FAM171A2$loglik   # 174.93, the largest log-likelihood value
# gamma
gamma_FAM171A2 <- fitdist(FAM171A2_data, "gamma", "mle")  
gamma_FAM171A2$loglik   # 107.35
# lognormal
lognormal_FAM171A2 <- fitdist(FAM171A2_data, "lnorm", "mle")
lognormal_FAM171A2$loglik   # 111.59
# weibull
weibull_FAM171A2 <- fitdist(FAM171A2_data, "weibull", "mle")  
weibull_FAM171A2$loglik   # -122.70
# inverse gaussian
invgauss_FAM171A2 <- mlinvgauss(FAM171A2_data)
logLik(invgauss_FAM171A2)   # 111.42
# logisitic
logis_FAM171A2 <- fitdist(FAM171A2_data, "logis", "mle")  
logis_FAM171A2$loglik   # 157.93

# total α-syn
# student's t
t_alpha_syn <- fitdistr(alpha_syn_data, "t")
t_alpha_syn$loglik   # -3367.75
# gamma
gamma_alpha_syn <- fitdist(alpha_syn_data, "gamma", "mle")  
gamma_alpha_syn$loglik   # -3339.53
# weibull
weibull_alpha_syn <- fitdist(alpha_syn_data, "weibull", "mle")  
weibull_alpha_syn$loglik   # -3373.58
# lognormal
lognormal_alpha_syn <- fitdist(alpha_syn_data, "lnorm", "mle")  
lognormal_alpha_syn$loglik  # -3328.97, the largest log-likelihood value
# inverse gaussian
invgauss_alpha_syn <- mlinvgauss(alpha_syn_data)
logLik(invgauss_alpha_syn)  # -3329
# logisitic
logis_alpha_syn <- fitdist(alpha_syn_data, "logis", "mle")  
logis_alpha_syn$loglik  # -3374.47


####################### histogram and probability density curve ###################### 
# FAM171A2 level  
par(cex.axis = 0.7, cex.lab = 0.7, cex.main = 1)
hist_FAM171A2 <- hist(FAM171A2_data, breaks = 40, freq = FALSE, xlim = c(7,9.6), ylim = c(0,3.2),
                      main = "Histogram of CSF FAM171A2 level", xlab = "FAM171A2 level", 
                      ylab = "density", col = "lightblue", border = "black") 
x <- seq(7,9.6, by = 0.01)
lines(x,dgamma(x,shape = gamma_FAM171A2$estimate[1],rate = gamma_FAM171A2$estimate[2]), 
      type = "l", lty = 2, col = "black", lwd = 2)
lines(x,dlnorm(x,meanlog = lognormal_FAM171A2$estimate[1],sdlog = lognormal_FAM171A2$estimate[2]), 
      type = "l", lty = 3, col = "orange", lwd = 2)
lines(x,dweibull(x,shape = weibull_FAM171A2$estimate[1],scale = weibull_FAM171A2$estimate[2]), 
      type = "l", lty = 4, col = "#9ACD32", lwd = 2)
lines(x,dinv.gaussian(x, mu = coefficients(invgauss_FAM171A2)[1], lambda = coefficients(invgauss_FAM171A2)[2]), 
      type = "l", lty = 5, col = "blue", lwd = 2)
lines(x,dlogis(x, location = logis_FAM171A2$estimate[1], scale = logis_FAM171A2$estimate[2]), 
      type = "l", lty = 6, col = "purple", lwd = 2)
t_density <- dt((x - t_FAM171A2$estimate[1]) / t_FAM171A2$estimate[2], t_FAM171A2$estimate[3]) / t_FAM171A2$estimate[2]
lines(x, t_density, type = "l", lty = 1, col = "red", lwd = 3)
legend("topright", legend = c("t", "gamma", "lognormal", "weibull", "inverse gaussian", "logistic"),
       x.intersp = 0.2, y.intersp = 0.09, inset = c(-0.65, -0.35),
       col = c("red", "black", "orange", "#9ACD32","blue", "purple"),
       lty = c(1, 2, 3, 4, 5, 6), lwd = c(3, 2, 2, 2, 2, 2), bty = "n", cex = 0.57)

quantile_reg1_FAM171A2_min <- 0
quantile_reg1_FAM171A2_max <- pt((reg1_FAM171A2_max - t_FAM171A2$estimate[1])/t_FAM171A2$estimate[2], df = t_FAM171A2$estimate[3])
quantile_reg2_FAM171A2_min <- pt((reg2_FAM171A2_min - t_FAM171A2$estimate[1])/t_FAM171A2$estimate[2], df = t_FAM171A2$estimate[3])
quantile_reg2_FAM171A2_max <- pt((reg2_FAM171A2_max - t_FAM171A2$estimate[1])/t_FAM171A2$estimate[2], df = t_FAM171A2$estimate[3])
quantile_reg3_FAM171A2_min <- pt((reg3_FAM171A2_min - t_FAM171A2$estimate[1])/t_FAM171A2$estimate[2], df = t_FAM171A2$estimate[3])
quantile_reg3_FAM171A2_max <- pt((reg3_FAM171A2_max - t_FAM171A2$estimate[1])/t_FAM171A2$estimate[2], df = t_FAM171A2$estimate[3])
quantile_reg4_FAM171A2_min <- pt((reg4_FAM171A2_min - t_FAM171A2$estimate[1])/t_FAM171A2$estimate[2], df = t_FAM171A2$estimate[3])
quantile_reg4_FAM171A2_max <- pt((reg4_FAM171A2_max - t_FAM171A2$estimate[1])/t_FAM171A2$estimate[2], df = t_FAM171A2$estimate[3])
quantile_reg5_FAM171A2_min <- pt((reg5_FAM171A2_min - t_FAM171A2$estimate[1])/t_FAM171A2$estimate[2], df = t_FAM171A2$estimate[3])
quantile_reg5_FAM171A2_max <- pt((reg5_FAM171A2_max - t_FAM171A2$estimate[1])/t_FAM171A2$estimate[2], df = t_FAM171A2$estimate[3])
quantile_reg6_FAM171A2_min <- pt((reg6_FAM171A2_min - t_FAM171A2$estimate[1])/t_FAM171A2$estimate[2], df = t_FAM171A2$estimate[3])
quantile_reg6_FAM171A2_max <- 1


# total α-syn 
par(cex.axis = 0.7, cex.lab = 0.7, cex.main = 1)
hist_alpha_syn <- hist(alpha_syn_data, breaks = 40, freq = FALSE, xlim = c(0,5500), ylim = c(0,0.0009), 
                       main = TeX("Histogram of CSF total $\\alpha$-syn"), xlab = TeX("total $\\alpha$-syn level"), 
                       ylab = "density", col = "lightblue", border = "black")
x <- seq(0, 5500, by = 1)
lines(x,dgamma(x,shape = gamma_alpha_syn$estimate[1],rate = gamma_alpha_syn$estimate[2]), 
      type = "l", lty = 2, col = "black", lwd = 2)
lines(x,dweibull(x,shape = weibull_alpha_syn$estimate[1],scale = weibull_alpha_syn$estimate[2]), 
      type = "l", lty = 4, col = "#9ACD32", lwd = 2)
lines(x,dinv.gaussian(x, mu = coefficients(invgauss_alpha_syn)[1], lambda = coefficients(invgauss_alpha_syn)[2]), 
      type = "l", lty = 5, col = "blue", lwd = 2)
lines(x,dlogis(x, location = logis_alpha_syn$estimate[1], scale = logis_alpha_syn$estimate[2]), 
      type = "l", lty = 6, col = "purple", lwd = 2)
t_density <- dt((x - t_alpha_syn$estimate[1]) / t_alpha_syn$estimate[2], t_alpha_syn$estimate[3]) / t_alpha_syn$estimate[2]
lines(x, t_density, type = "l", lty = 3, col = "orange", lwd = 2)
lines(x,dlnorm(x,meanlog = lognormal_alpha_syn$estimate[1],sdlog = lognormal_alpha_syn$estimate[2]), 
      type = "l", lty = 1, col = "red", lwd = 3)
legend("topright", legend = c("lognormal", "gamma", "weibull", "inverse gaussian", "logistic", "t"),
       x.intersp = 0.2, y.intersp = 0.09, inset = c(-0.65, -0.35),
       col = c("red", "black", "#9ACD32", "blue", "purple", "orange"),
       lty = c(1, 2, 4, 5, 6, 3), lwd = c(3, 2, 2, 2, 2, 2), bty = "n", cex = 0.57)

quantile_reg1_alpha_syn_min <- plnorm(reg1_alpha_syn_min, meanlog = lognormal_alpha_syn$estimate[1], sdlog = lognormal_alpha_syn$estimate[2])
quantile_reg1_alpha_syn_max <- 1
quantile_reg2_alpha_syn_min <- plnorm(reg2_alpha_syn_min, meanlog = lognormal_alpha_syn$estimate[1], sdlog = lognormal_alpha_syn$estimate[2])
quantile_reg2_alpha_syn_max <- plnorm(reg2_alpha_syn_max, meanlog = lognormal_alpha_syn$estimate[1], sdlog = lognormal_alpha_syn$estimate[2])
quantile_reg3_alpha_syn_min <- plnorm(reg3_alpha_syn_min, meanlog = lognormal_alpha_syn$estimate[1], sdlog = lognormal_alpha_syn$estimate[2])
quantile_reg3_alpha_syn_max <- plnorm(reg3_alpha_syn_max, meanlog = lognormal_alpha_syn$estimate[1], sdlog = lognormal_alpha_syn$estimate[2])
quantile_reg4_alpha_syn_min <- plnorm(reg4_alpha_syn_min, meanlog = lognormal_alpha_syn$estimate[1], sdlog = lognormal_alpha_syn$estimate[2])
quantile_reg4_alpha_syn_max <- plnorm(reg4_alpha_syn_max, meanlog = lognormal_alpha_syn$estimate[1], sdlog = lognormal_alpha_syn$estimate[2])
quantile_reg5_alpha_syn_min <- plnorm(reg5_alpha_syn_min, meanlog = lognormal_alpha_syn$estimate[1], sdlog = lognormal_alpha_syn$estimate[2])
quantile_reg5_alpha_syn_max <- plnorm(reg5_alpha_syn_max, meanlog = lognormal_alpha_syn$estimate[1], sdlog = lognormal_alpha_syn$estimate[2])
quantile_reg6_alpha_syn_min <- 0
quantile_reg6_alpha_syn_max <- plnorm(reg6_alpha_syn_max, meanlog = lognormal_alpha_syn$estimate[1], sdlog = lognormal_alpha_syn$estimate[2])


##################### probability integral transformation #####################
cdf_FAM171A2 <- pt((FAM171A2_data - t_FAM171A2$estimate[1])/t_FAM171A2$estimate[2], df = t_FAM171A2$estimate[3])
cdf_alpha_syn <- plnorm(alpha_syn_data, meanlog = lognormal_alpha_syn$estimate[1], sdlog = lognormal_alpha_syn$estimate[2])
uv_FAM171A2_alpha_syn <- data.frame(x = cdf_FAM171A2, y = cdf_alpha_syn)


################## scatter plot of cdf_FAM171A2 and cdf_α-syn ##################
p <- ggplot(uv_FAM171A2_alpha_syn, aes(x = x, y = y)) + 
     geom_point(color = "white", fill = "#1eb3b9", size = 2.4, shape = 21, stroke = 0.5) + 
     # region 1
     geom_rect(aes(xmin = 0, xmax = quantile_reg1_FAM171A2_max, ymin = quantile_reg1_alpha_syn_min, ymax = 1),
               color = "#ea6458", size = 0.5, fill = "NA", alpha = .02) +
     # region 2
     geom_rect(aes(xmin = quantile_reg2_FAM171A2_min, xmax = quantile_reg2_FAM171A2_max, ymin = quantile_reg2_alpha_syn_min, ymax = quantile_reg2_alpha_syn_max),
               color = "#f4ad60", size = 0.5, fill = "NA", alpha = .02) + 
     # region 3
     geom_rect(aes(xmin = quantile_reg3_FAM171A2_min, xmax = quantile_reg3_FAM171A2_max, ymin = quantile_reg3_alpha_syn_min, ymax = quantile_reg3_alpha_syn_max),
               color = "#ea6458", size = 0.5, fill = "NA", alpha = .02) + 
     # region 4
     geom_rect(aes(xmin = quantile_reg4_FAM171A2_min, xmax = quantile_reg4_FAM171A2_max, ymin = quantile_reg4_alpha_syn_min, ymax = quantile_reg4_alpha_syn_max),
               color = "#f4ad60", size = 0.5, fill = "NA", alpha = .02) +
     # region 5
     geom_rect(aes(xmin = quantile_reg5_FAM171A2_min, xmax = quantile_reg5_FAM171A2_max, ymin = quantile_reg5_alpha_syn_min, ymax = quantile_reg5_alpha_syn_max),
               color = "#ea6458", size = 0.5, fill = "NA", alpha = .02) +
     # region 6
     geom_rect(aes(xmin = quantile_reg6_FAM171A2_min, xmax = 1, ymin = 0, ymax = quantile_reg6_alpha_syn_max),
               color = "#f4ad60", size = 0.5, fill = "NA", alpha = .02) +
     labs(title = "", x = expression(u^FAM171A2), y = expression(v^{alpha-syn})) + 
     xlim(c(0, 1)) + ylim(c(0, 1)) + theme(aspect.ratio = 1) +
     theme(plot.title = element_text(hjust = 0.5)) + theme(plot.margin = unit(c(0, 0, 0, 0), "mm")) +
     theme(axis.text.x = element_text(size = 8.49), axis.text.y = element_text(size = 8.49)) +
     theme(axis.title = element_text(size = 9.6))
print(p)

ggsave("scatter_uv_example_1.eps", plot=p, device="eps", width = 3.2, height = 3.2, units = "in")


####################### save to csv file #######################
write.csv(1-cdf_FAM171A2, file = "data/inv_u_FAM171A2.csv", row.names = FALSE)
write.csv(cdf_alpha_syn, file = "data/v_alpha_syn.csv", row.names = FALSE)


