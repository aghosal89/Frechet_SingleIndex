
## This script creates the Global Fréchet regression model for the mortality distributions

# read libraries
library('latex2exp')
library(ggplot2)
library(tidyverse)

# Set working directory

# the following are the list of 40 countries for which mortality is considered in our model.
country_ <- read.csv("Countries_FSI.csv", header = T)[,-1]

# read the covariate data:
X_ctr<- read.csv("X_centerscale.csv", header = T)
X_ctr <- as.matrix(X_ctr[,-1])
rownames(X_ctr)<- country_

# Read the response data as quantiles
quant_all <- read.csv("quant_all.csv", header = T)
rownames(quant_all)<- country_
quant_all <- as.matrix(quant_all[,-1])

# length of quantiles 
m <- ncol(quant_all)

# support for quantiles
qSup <- seq(0,1, length.out=m)

# support for densities
dSup <- seq(20, 110, length.out=m)

# create vector for the regression coefficients
b<- matrix(NaN, nrow = ncol(X_ctr)+1, ncol = ncol(quant_all))
b_ucl95 <- b_lcl95 <- b

for (i in 1:m) {
  lm_temp<- lm(quant_all[,i]~X_ctr)
  b[,i]<- lm_temp$coefficients
  
  b_lcl95[,i]<- confint(lm_temp)[,1]
  b_ucl95[,i]<- confint(lm_temp)[,2]
}

# GDP YoY %age change
plot(qSup,b[2,],main= "GDP YoY %-age change" , ylim = c(-2.5, 1),
     xlab="Quantile,t", ylab=TeX(r"(${\beta}(t)\,\,{estimates}$)"), type="l")
lines(qSup, b_ucl95[2,], type ="l", col="green2")
lines(qSup, b_lcl95[2,], type ="l", col="blue2")
abline(h=0, pch=2, col="red")

# Current HealthCare Expenditure %of GDP 
plot(qSup,b[3,],main= "Current HCE %of GDP" , ylim = c(-2.0, 2),
     xlab="Quantile,t", ylab=TeX(r"(${\beta}(t)\,\,{estimates}$)"), type="l")
lines(qSup, b_ucl95[3,], type ="l", col="green2")
lines(qSup, b_lcl95[3,], type ="l", col="blue2")
abline(h=0, pch=2, col="red")

# CO2 Emissions metric tonnes per capita 
plot(qSup,b[4,],main= "CO2E metric tonnes PC" , ylim = c(-2.4, 1),
     xlab="Quantile,t", ylab=TeX(r"(${\beta}(t)\,\,{estimates}$)"), type="l")
lines(qSup, b_ucl95[4,], type ="l", col="green2")
lines(qSup, b_lcl95[4,], type ="l", col="blue2")
abline(h=0, pch=2, col="red")

# Infant Mortality per 1000 liver births 
plot(qSup,b[5,],main= "IM/1000 live births" , ylim = c(-2.8, 1),
     xlab="Quantile,t", ylab=TeX(r"(${\beta}(t)\,\,{estimates}$)"), type="l")
lines(qSup, b_ucl95[5,], type ="l", col="green2")
lines(qSup, b_lcl95[5,], type ="l", col="blue2")
abline(h=0, pch=2, col="red")

# Human Development Index 
plot(qSup,b[6,], main= "HDI", ylim = c(-.8, 6.5),
     xlab="Quantile,t", ylab=TeX(r"(${\beta}(t)\,\,{estimates}$)"), type="l")
lines(qSup, b_ucl95[6,], type ="l", col="green2")
lines(qSup, b_lcl95[6,], type ="l", col="blue2")
abline(h=0, pch=2, col="red")

# Integrate the coefficient functions over the quantiles

# GDP YoY %age change
th<- sapply(2:6, function(i) fdadensity:::trapzRcpp(X= qSup, Y=b[i,]))

th<- th/sqrt(sum(th^2))
th<- th[c(5, 2, 1, 4, 3)] # reorder according to Table 3

# Plot all coefficient function estimates together

betaDF <- data.frame(t = qSup, HDI = b[6, ], HCE = b[3, ], GDP = b[2, ], IM = b[5, ], CO2 = b[4, ]) %>% 
  gather("Covariate", "Beta", -t)
betaDF$Covariate = factor(betaDF$Covariate, 
                          levels = c("HDI", "HCE", "GDP", "IM", "CO2"))

png(filename = 'GlobalFrechetBetas.png', width = 3.5, height = 2.5, units = "in", res = 200)
p <- ggplot(data = betaDF) + 
  geom_line(mapping = aes(x = t, y = Beta, color = Covariate)) +
  scale_color_manual(values = c("purple3",
                                "green3",
                                "pink3",
                                "skyblue2",
                                "grey")) +
  labs(x = 't', y = TeX(r"(${\hat{\beta}}(t)$)")) + 
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6))
p
dev.off()




