###################################################################################################
# Example 2: COVID-19 data
###################################################################################################
# Author: Zaixin Huang
# Email: eric.huangzaixin@gmail.com
# This version 2024.11.05 
# The latest version can be downloaded from https://github.com/Huangzaixin/local-dependence-toolbox
###################################################################################################
install.packages("ggplot2")
install.packages("distributionsrd")
install.packages("ggExtra")

library(ggplot2)
library(distributionsrd)    # rightparetolognormal.mle
library(ggExtra)

#######################  read data  ######################## 
all_deaths_counts <- read.csv("data/us-counties_20210501.csv",header = T)
cases_data <- all_deaths_counts$cases
deaths_data <- all_deaths_counts$deaths
cases_deaths_data <- data.frame(x = cases_data, y = deaths_data)

#######################  scatter plot ######################
label_scientific <- function(l) { 
        l <- format(l, scientific = TRUE)
        l <- gsub("0e\\+00","0",l)
        l <- gsub("^(.*)e", "'\\1'e", l)
        l <- gsub("e\\+","e",l)
        l <- gsub("e", "%*%10^", l)
        l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
        parse(text=l)
}

# scatter plot 1: whole dataset
p <- ggplot(cases_deaths_data, aes(x = x, y = y)) + 
        geom_point(color = "black", fill = "lightblue", size = 1.7, shape = 21, stroke = 0.5) + 
        xlim(c(0, 600000)) + scale_x_continuous(limits = c(0, 600000),breaks = seq(from = 0, to = 600000, by = 200000),labels=label_scientific) +
        ylim(c(0, 12000)) + scale_y_continuous(limits = c(0, 12000),breaks = seq(from = 0, to = 12000, by = 3000),labels=label_scientific) +
        labs(title = "", x = "cases", y = "deaths") +
        theme(aspect.ratio = 1) + theme(plot.title = element_text(hjust = 0.5)) +
        theme(plot.title = element_text(hjust = 0.5)) + theme(plot.margin = unit(c(0, 3, 0, 0), "mm")) +
        theme(axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9))
print(p)


# scatter plot 2: points with particularly large values are not shown 
region1_xmax <- 1000
region1_ymax <- 18

region2_xmax <- 2000
region2_ymax <- 40

region3_xmax <- 4000
region3_ymax <- 80

region4_xmax <- 10000
region4_ymax <- 200

region5_xmin <- 1200
region5_xmax <- 40000
region5_ymin <- 22
region5_ymax <- 700

region6_xmin <- 20000
region7_xmin <- 32000
region8_xmin <- 45000

region6_ymin <- 350
region7_ymin <- 600
region8_ymin <- 800

line_1 <- data.frame( x = c(region6_xmin, region6_xmin), y = c(region6_ymin, 1400))
line_2 <- data.frame( x = c(region7_xmin, region7_xmin), y = c(region7_ymin, 1400))
line_3 <- data.frame( x = c(region8_xmin, region8_xmin), y = c(region8_ymin, 1400))

line_4 <- data.frame( x = c(region6_xmin, 80000), y = c(region6_ymin, region6_ymin))
line_5 <- data.frame( x = c(region7_xmin, 80000), y = c(region7_ymin, region7_ymin))
line_6 <- data.frame( x = c(region8_xmin, 80000), y = c(region8_ymin, region8_ymin))

line_7 <- data.frame( x = c(region6_xmin, 80000), y = c(1400, 1400))
line_8 <- data.frame( x = c(80000, 80000), y = c(region6_ymin, 1400))

p <- ggplot(cases_deaths_data, aes(x = x, y = y)) + 
        geom_point(color = "black", fill = "lightblue", size = 1.7, shape = 21, stroke = 0.5) + 
        # region 1
        geom_rect(aes(xmin = 0, xmax = region1_xmax, ymin = 0, ymax = region1_ymax), color = "#00adb5", size = 0.5, fill = "NA", alpha = .02) +
        # region 2
        geom_rect(aes(xmin = 0, xmax = region2_xmax, ymin = 0, ymax = region2_ymax), color = "#00adb5", size = 0.5, fill = "NA", alpha = .02) +
        # region 3
        geom_rect(aes(xmin = 0, xmax = region3_xmax, ymin = 0, ymax = region3_ymax), color = "#00adb5", size = 0.5, fill = "NA", alpha = .02) + 
        # region 4
        geom_rect(aes(xmin = 0, xmax = region4_xmax, ymin = 0, ymax = region4_ymax), color = "#00adb5", size = 0.5, fill = "NA", alpha = .02) + 
        # region 5
        geom_rect(aes(xmin = region5_xmin, xmax = region5_xmax, ymin = region5_ymin, ymax = region5_ymax), color = "#ff5722", size = 0.5, fill = "NA", alpha = .02) +
        # region 6、7、8
        geom_segment(data = line_1, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), color = "#00adb5", size = 0.5) +
        geom_segment(data = line_2, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), color = "#00adb5", size = 0.5) +
        geom_segment(data = line_3, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), color = "#00adb5", size = 0.5) +
        geom_segment(data = line_4, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), color = "#00adb5", size = 0.5) +
        geom_segment(data = line_5, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), color = "#00adb5", size = 0.5) +
        geom_segment(data = line_6, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), color = "#00adb5", size = 0.5) +
        geom_segment(data = line_7, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), color = "#00adb5", linetype = "longdash", size = 0.5) +
        geom_segment(data = line_8, aes(x = x[1], y = y[1], xend = x[2], yend = y[2]), color = "#00adb5", linetype = "longdash", size = 0.5) +
        
        xlim(c(0, 80000)) + scale_x_continuous(limits = c(0, 80000),breaks = seq(from = 0, to = 80000, by = 20000),labels=label_scientific) +
        ylim(c(0, 1400)) + scale_y_continuous(limits = c(0, 1400),breaks = seq(from = 0, to = 1400, by = 350),labels=label_scientific) +
        
        labs(title = "", x = "cases", y = "deaths") +
        theme(aspect.ratio = 1) + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.margin = unit(c(0, 3, 0, 0), "mm")) +
        theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) +
        theme(axis.title=element_text(size=8.5)) 

# densigram
p <- ggMarginal(p, type = "densigram", xparams = list(color = "#222", fill = "lightblue", size = 0.35, bins = 40), 
                yparams = list(color = "#222", fill = "lightblue", size = 0.35, bins = 40)) 

print(p)

ggsave("scatter_example_2.eps", plot=p, device="eps", width = 3.55, height = 3.55, units = "in")


#######################  fit lognormal-Pareto distribution ###################### 
est_lnor_pareto_cases <- rightparetolognormal.mle(cases_data)
est_lnor_pareto_deaths <- rightparetolognormal.mle(deaths_data)

#######################  histogram and probability density curve ###################### 
# cases
display_upper_bound <- 80000
hist_cases <- hist(cases_data, breaks = 4000, probability = TRUE, 
                   xlim = c(0,display_upper_bound), ylim = c(0,0.0003),
                   xlab = "cumulative cases", ylab = "density", 
                   col = "lightblue", border = "black", main = "Histogram of cases")

x <- seq(0, display_upper_bound, by = 1)
lines(x,drightparetolognormal(x, shape2 = est_lnor_pareto_cases$coefficients[3], 
                                 meanlog = est_lnor_pareto_cases$coefficients[1], 
                                 sdlog = est_lnor_pareto_cases$coefficients[2]), 
                                 type = "l", col = "red", lwd = 2) 

# deaths
display_upper_bound <- 1400
hist_deaths <- hist(deaths_data, breaks = 5000, probability = TRUE, 
                   xlim = c(0,display_upper_bound), ylim = c(0,0.02),
                   xlab = "cumulative deaths", ylab = "density", 
                   col = "lightblue", border = "black", main = "Histogram of deaths")
x <- seq(0, display_upper_bound, by = 1)
lines(x,drightparetolognormal(x, shape2 = est_lnor_pareto_deaths$coefficients[3], 
                              meanlog = est_lnor_pareto_deaths$coefficients[1], 
                              sdlog = est_lnor_pareto_deaths$coefficients[2]), 
                              type = "l", col = "red", lwd = 2) 

#######################  probability integral transformation  ###################### 
# cases
cdf_cases <- prightparetolognormal(cases_data,shape2 = est_lnor_pareto_cases$coefficients[3],
                                   meanlog = est_lnor_pareto_cases$coefficients[1],
                                   sdlog = est_lnor_pareto_cases$coefficients[2])

cases_95_quantile <- qrightparetolognormal(0.95,shape2 = est_lnor_pareto_cases$coefficients[3],
                                           meanlog = est_lnor_pareto_cases$coefficients[1],
                                           sdlog = est_lnor_pareto_cases$coefficients[2])

cdf_cases_region1_xmax <- prightparetolognormal(region1_xmax,shape2 = est_lnor_pareto_cases$coefficients[3],
                                       meanlog = est_lnor_pareto_cases$coefficients[1],
                                       sdlog = est_lnor_pareto_cases$coefficients[2])
cdf_cases_region2_xmax <- prightparetolognormal(region2_xmax,shape2 = est_lnor_pareto_cases$coefficients[3],
                                       meanlog = est_lnor_pareto_cases$coefficients[1],
                                       sdlog = est_lnor_pareto_cases$coefficients[2])
cdf_cases_region3_xmax <- prightparetolognormal(region3_xmax,shape2 = est_lnor_pareto_cases$coefficients[3],
                                       meanlog = est_lnor_pareto_cases$coefficients[1],
                                       sdlog = est_lnor_pareto_cases$coefficients[2])
cdf_cases_region4_xmax <- prightparetolognormal(region4_xmax,shape2 = est_lnor_pareto_cases$coefficients[3],
                                       meanlog = est_lnor_pareto_cases$coefficients[1],
                                       sdlog = est_lnor_pareto_cases$coefficients[2])
cdf_cases_region5_xmin <- prightparetolognormal(region5_xmin,shape2 = est_lnor_pareto_cases$coefficients[3],
                                       meanlog = est_lnor_pareto_cases$coefficients[1],
                                       sdlog = est_lnor_pareto_cases$coefficients[2])
cdf_cases_region5_xmax <- prightparetolognormal(region5_xmax,shape2 = est_lnor_pareto_cases$coefficients[3],
                                       meanlog = est_lnor_pareto_cases$coefficients[1],
                                       sdlog = est_lnor_pareto_cases$coefficients[2])
cdf_cases_region6_xmin <- prightparetolognormal(region6_xmin,shape2 = est_lnor_pareto_cases$coefficients[3],
                                       meanlog = est_lnor_pareto_cases$coefficients[1],
                                       sdlog = est_lnor_pareto_cases$coefficients[2])
cdf_cases_region7_xmin <- prightparetolognormal(region7_xmin,shape2 = est_lnor_pareto_cases$coefficients[3],
                                       meanlog = est_lnor_pareto_cases$coefficients[1],
                                       sdlog = est_lnor_pareto_cases$coefficients[2])
cdf_cases_region8_xmin <- prightparetolognormal(region8_xmin,shape2 = est_lnor_pareto_cases$coefficients[3],
                                       meanlog = est_lnor_pareto_cases$coefficients[1],
                                       sdlog = est_lnor_pareto_cases$coefficients[2])

# deaths
cdf_deaths <- prightparetolognormal(deaths_data,shape2 = est_lnor_pareto_deaths$coefficients[3],
                                    meanlog = est_lnor_pareto_deaths$coefficients[1],
                                    sdlog = est_lnor_pareto_deaths$coefficients[2])

deaths_95_quantile <- qrightparetolognormal(0.95,shape2 = est_lnor_pareto_deaths$coefficients[3],
                                            meanlog = est_lnor_pareto_deaths$coefficients[1],
                                            sdlog = est_lnor_pareto_deaths$coefficients[2])

cdf_deaths_region1_ymax <- prightparetolognormal(region1_ymax,shape2 = est_lnor_pareto_deaths$coefficients[3],
                                      meanlog = est_lnor_pareto_deaths$coefficients[1],
                                      sdlog = est_lnor_pareto_deaths$coefficients[2])
cdf_deaths_region2_ymax <- prightparetolognormal(region2_ymax,shape2 = est_lnor_pareto_deaths$coefficients[3],
                                      meanlog = est_lnor_pareto_deaths$coefficients[1],
                                      sdlog = est_lnor_pareto_deaths$coefficients[2])
cdf_deaths_region3_ymax <- prightparetolognormal(region3_ymax,shape2 = est_lnor_pareto_deaths$coefficients[3],
                                      meanlog = est_lnor_pareto_deaths$coefficients[1],
                                      sdlog = est_lnor_pareto_deaths$coefficients[2])
cdf_deaths_region4_ymax <- prightparetolognormal(region4_ymax,shape2 = est_lnor_pareto_deaths$coefficients[3],
                                       meanlog = est_lnor_pareto_deaths$coefficients[1],
                                       sdlog = est_lnor_pareto_deaths$coefficients[2])
cdf_deaths_region5_ymin <- prightparetolognormal(region5_ymin,shape2 = est_lnor_pareto_deaths$coefficients[3],
                                       meanlog = est_lnor_pareto_deaths$coefficients[1],
                                       sdlog = est_lnor_pareto_deaths$coefficients[2])
cdf_deaths_region5_ymax <- prightparetolognormal(region5_ymax,shape2 = est_lnor_pareto_deaths$coefficients[3],
                                       meanlog = est_lnor_pareto_deaths$coefficients[1],
                                       sdlog = est_lnor_pareto_deaths$coefficients[2])
cdf_deaths_region6_ymin <- prightparetolognormal(region6_ymin,shape2 = est_lnor_pareto_deaths$coefficients[3],
                                       meanlog = est_lnor_pareto_deaths$coefficients[1],
                                       sdlog = est_lnor_pareto_deaths$coefficients[2])
cdf_deaths_region7_ymin <- prightparetolognormal(region7_ymin,shape2 = est_lnor_pareto_deaths$coefficients[3],
                                       meanlog = est_lnor_pareto_deaths$coefficients[1],
                                       sdlog = est_lnor_pareto_deaths$coefficients[2])
cdf_deaths_region8_ymin <- prightparetolognormal(region8_ymin,shape2 = est_lnor_pareto_deaths$coefficients[3],
                                       meanlog = est_lnor_pareto_deaths$coefficients[1],
                                       sdlog = est_lnor_pareto_deaths$coefficients[2])


#######################  scatter plot for cdf_cases and cdf_deaths  ######################
uv_cases_deaths <- data.frame(x = cdf_cases, y = cdf_deaths)
p <- ggplot(uv_cases_deaths, aes(x = x, y = y)) + 
     geom_point(color = "black", fill = "lightblue", size = 1.7, shape = 21, stroke = 0.5) + 
     
     geom_rect(aes(xmin = 0, xmax = cdf_cases_region1_xmax, ymin = 0, ymax = cdf_deaths_region1_ymax), color = "#00adb5", size = 0.5, fill = "NA", alpha = .02) +
     geom_rect(aes(xmin = 0, xmax = cdf_cases_region2_xmax, ymin = 0, ymax = cdf_deaths_region2_ymax), color = "#00adb5", size = 0.5, fill = "NA", alpha = .02) +
     geom_rect(aes(xmin = 0, xmax = cdf_cases_region3_xmax, ymin = 0, ymax = cdf_deaths_region3_ymax), color = "#00adb5", size = 0.5, fill = "NA", alpha = .02) + 
     geom_rect(aes(xmin = 0, xmax = cdf_cases_region4_xmax, ymin = 0, ymax = cdf_deaths_region4_ymax), color = "#00adb5", size = 0.5, fill = "NA", alpha = .02) + 
        
     geom_rect(aes(xmin = cdf_cases_region5_xmin, xmax = cdf_cases_region5_xmax, ymin = cdf_deaths_region5_ymin, ymax = cdf_deaths_region5_ymax), color = "#ff5722", size = 0.5, fill = "NA", alpha = .02) +         
     
     geom_rect(aes(xmin = cdf_cases_region6_xmin, xmax = 1, ymin = cdf_deaths_region6_ymin, ymax = 1), color = "#00adb5", size = 0.5, fill = "NA", alpha = .02) + # #FED770
     geom_rect(aes(xmin = cdf_cases_region7_xmin, xmax = 1, ymin = cdf_deaths_region7_ymin, ymax = 1), color = "#00adb5", size = 0.5, fill = "NA", alpha = .02) +
     geom_rect(aes(xmin = cdf_cases_region8_xmin, xmax = 1, ymin = cdf_deaths_region8_ymin, ymax = 1), color = "#00adb5", size = 0.5, fill = "NA", alpha = .02) + 
     
     labs(title = "", x = expression(u^cases), y = expression(v^deaths)) +
     xlim(c(0, 1)) + ylim(c(0, 1)) + theme(aspect.ratio = 1) +
     theme(plot.title = element_text(hjust = 0.5)) + theme(plot.margin = unit(c(0, 3, 0, 0), "mm")) +
     theme(axis.text.x = element_text(size = 7.48), axis.text.y = element_text(size = 7.48)) +
     theme(axis.title=element_text(size=10))       
print(p)

ggsave("scatter_uv_example_2.eps", plot=p, device="eps", width = 3, height = 3, units = "in")


####################### save to csv file #######################
write.csv(cdf_cases, file = "data/u_cases.csv", row.names = FALSE)
write.csv(cdf_deaths, file = "data/v_deaths.csv", row.names = FALSE)


 