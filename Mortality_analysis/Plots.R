## This script contains the codes to generate figures 7 - 10 in the paper. 

# set working directory~
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/FSI/Mortality_all")

## Codes for producing figure 7 in the document
###################################################
# plot all densities with top 6 and bottom 6 modes
###################################################

# read the mortality densities of all countries
density_all<- read.csv("density_all.csv", header= T)
country_ <- density_all[,1]
rownames(density_all) <- country_
density_all <- as.matrix(density_all[,-1])

# equidistant grid for densities 
dSup <- seq(20,110, length.out=101)

# equidistant grid for quantiles 
qSup <- seq(0,1, length.out = 101)

# Create plot
plot(dSup, density_all[2,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey', ylim = c(0, 0.045))
lines(dSup, density_all[3,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[8,], lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[9,],  xlab='Age',lwd=2, ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[10,], xlab='Age',lwd=2, ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[11,], xlab='Age',lwd=2, ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[12,], xlab='Age',lwd=2, ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[13,], xlab='Age',lwd=2, ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[14,], xlab='Age',lwd=2, ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[16,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[17,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[19,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[21,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[22,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[23,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[24,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[26,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[27,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[28,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[30,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[31,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[32,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[33,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[34,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[37,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[38,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[40,],lwd=2, xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[1,], 
      xlab='Age', ylab='Density', type='l', lwd=2, col=rgb(red = 0, green=0, blue = 1, alpha = 0.5))
lines(dSup, density_all[4,], lwd=2,
      xlab='Age', ylab='Density', type='l',col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5))
lines(dSup, density_all[5,], lwd=2,
      xlab='Age', ylab='Density', type='l',col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5))
lines(dSup, density_all[6,], lwd=2,
      xlab='Age', ylab='Density', type='l',col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
lines(dSup, density_all[7,], lwd=2,
      xlab='Age', ylab='Density', type='l',col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
lines(dSup, density_all[15,], lwd=2,
      xlab='Age', ylab='Density', type='l',col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
lines(dSup, density_all[18,], lwd=2,
      xlab='Age', ylab='Density', type='l',col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
lines(dSup, density_all[20,], lwd=2,
      xlab='Age', ylab='Density', type='l',col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5))
lines(dSup, density_all[25,], lwd=2,
      xlab='Age', ylab='Density', type='l',col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5))
lines(dSup, density_all[29,], lwd=2,
      xlab='Age', ylab='Density', type='l',col='lightgrey')
lines(dSup, density_all[35,], lwd=2,
      xlab='Age', ylab='Density', type='l',col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5))
lines(dSup, density_all[36,], lwd=2,
      xlab='Age', ylab='Density', type='l',col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5))
lines(dSup, density_all[39,], lwd=2,
      xlab='Age', ylab='Density', type='l',col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5))

# Add a legend to the plot
legend(20, 0.045, legend=c("Bottom 6", "Top 6"),
       col=c("red", "blue"), lty=c(1,1), cex=1.2,
       title="Modes of distribution", text.font=10, bg='lightblue')
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)




################################################
## Codes for producing figure 8 in the document
################################################

# This function computes the mode age of the mortality density over a grid
# Inputs:  1) dSup     : An equispaced grid of points of length m on the support of densities
#          2) densities: matrix of dimension rxm, each row represents a density on dSup. Total r densities.
# Outputs: a numeric vector of length r of modes of the r densities in the respective order.

mode_ages <- function(dSup, densities) {
  r<- nrow(densities)
  modes <- matrix(0, r,1)
  for (j in 1:r) {
    modes[j] <- dSup[which.max(densities[j,])]
  }
  return(modes)
}

# load library for frechet regression
library('frechet')

# read predicted densities from Global Frechet regression
gf_dpred<- read.csv("GF_Dpred.csv", header = T)[,-1]

# read predicted densities from Local Frechet regression with HDI covariate
lf_hdi_dpred <- read.csv("LF_HDI_Dpred.csv", header = T)[,-1]

# read predicted densities from Frechet Single index regression model 
fsi_dpred <- read.csv("FSI_Dpred.csv", header = T)[,-1]


# create the dataset for the plot
df_gf<- data.frame(Age=rep(dSup, 40), Method=rep("Global Frechet", 40*101), 
                   Density=as.vector(t(gf_dpred)), 
                   Country=rep(country_, each=101), 
                   Mode_Age= rep(mode_ages(dSup, gf_dpred), each=101))

df_lf_hdi<- data.frame(Age=rep(dSup, 40), Method=rep("LF:HDI", 40*101), 
                       Density=as.vector(t(lf_hdi_dpred)), 
                       Country=rep(country_, each=101), 
                       Mode_Age= rep(mode_ages(dSup, lf_hdi_dpred), each=101))

df_fsi <-data.frame(Age=rep(dSup, 40), Method=rep("FSI", 40*101), 
                    Density=as.vector(t(fsi_dpred)), 
                    Country=rep(country_, each=101), 
                    Mode_Age= rep(mode_ages(dSup,fsi_dpred), each=101))


df_all_actual <- data.frame(Age=rep(dSup, 40), Method= rep('Actual Mortality', 40*101), 
                            Density = as.vector(t(as.matrix(density_all))), 
                            Country=rep(country_,each=101), 
                            Mode_Age=rep(mode_ages(dSup,as.matrix(density_all)),each=101))

df_pred_plot<- rbind(df_all_actual, df_gf, df_lf_hdi, df_fsi)

library("ggplot2")
ggplot(data = df_pred_plot, aes(x=Age, y=Density, fill=Country )) +
  geom_path(aes(colour=(Mode_Age)), size=.6, alpha=1) +
  facet_wrap(~Method)+
  #ggtitle("Mortality distributions of 2013") +
  scale_color_gradientn(colours = terrain.colors(10)) +
  labs(colour = "Mode Age")
  
ggsave(file="Rplot_models_predictions.eps")









################################################
## Codes for producing figure 9 in the document
################################################

df_mspe_fold <- read.csv("MSPE_folds.csv", header = T)

df_mspe_fold$Model <- factor(df_mspe_fold$Model)
df_mspe_fold$log.MSPE <- log(df_mspe_fold$MSPE)
df_mspe_fold$sqrt.MSPE <- sqrt(df_mspe_fold$MSPE)

my.bp <- ggplot(data =df_mspe_fold, aes(y=log.MSPE, x=Model, fill=Model))
my.bp <- my.bp + geom_boxplot()
#my.bp <- my.bp + ggtitle('Distribution of MSPE for various models')
my.bp <- my.bp+ylab('log MSPE')+xlab("Models")
my.bp + theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) 
ggsave(file="Rplot_mspe_models_compare_log.eps")


# plots to compare MSPE of LF(HDI) FSI
lf_hdi_folds <- read.csv("LF_HDI_folds.csv", header= T)
lf_hce_folds <- read.csv("LF_HCE_folds.csv", header=T)
fsi_folds <- read.csv("FSI_MSPE_folds.csv", header = T)
gf_folds<- read.csv("GF_folds.csv", header = T)

hdi_fsi_ratio_log <- log(lf_hdi_folds[,2]/fsi_folds[,2])
hce_fsi_ratio_log <- log(lf_hce_folds[,2]/fsi_folds[,2])
gf_fsi_ratio_log <- log(gf_folds[,2]/fsi_folds[,2])

log_ratio_df <- rbind(data.frame(Comparison = "LF(HDI):FSI",MSPE.Ratio = hdi_fsi_ratio_log),
                      data.frame(Comparison = "LF(HCE):FSI", MSPE.Ratio = hce_fsi_ratio_log),
                      data.frame(Comparison = "GF:FSI", MSPE.Ratio = gf_fsi_ratio_log))

my.bp1 <- ggplot(data = log_ratio_df, aes(y=MSPE.Ratio, x=Comparison, fill= Comparison))
my.bp1 <- my.bp1 + geom_boxplot()
#my.bp <- my.bp + ggtitle('Distribution of MSPE for various models')
my.bp1 <- my.bp1 + ylab('log Ratio')+xlab("Models")
my.bp1 + theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) + geom_hline(yintercept=0, linetype="dashed")
ggsave(file="Rplot_mspe_compare_hdi_hce_fsi.eps")








#################################################
## Codes for producing figure 10 in the document
#################################################


# read the covariate data for all countries
X_ctr<- read.csv("X_centerscale.csv", header = T)
country_<- X_ctr[,1]
X_ctr <- X_ctr[,-1]
rownames(X_ctr)<- country_

# Read the response data as quantiles
quant_all <- read.csv("quant_all.csv", header = T)
rownames(quant_all)<- quant_all[,1]
quant_all <- quant_all[,-1]

# creating effects plot for HDI:
l= seq(min(X_ctr[,'HDI']), max(X_ctr[,'HDI']), length.out= 100)
covar_median_hdi_effect<- matrix(0, length(l), 5)
m=apply(X_ctr, 2, median)

thetahat <- read.csv('Theta_Hat.csv', header = T)
thetahat<- t(thetahat[,2])

h_fsi <- read.csv("FSI_bw.csv", header = T)[,2]

covar_median_hdi_effect<- cbind(GDP_yoy=rep(m[1],100), HC_exp=rep(m[2],100), 
                                CO2emission=rep(m[3],100), Infantm=rep(m[4],100))
ef_hdi<-cbind(covar_median_hdi_effect, HDI=l)
ef_hdi<- ef_hdi%*%t(thetahat)

X_ctr <- matrix(as.vector(as.matrix(X_ctr)),40, 5, byrow = F)
colnames(X_ctr) <- c('GDPC', 'HCE', 'CO2E', 'IM', 'HDI')

temp_hdi_effect<- LocDenReg(xin= as.matrix(X_ctr)%*%t(thetahat), qin=as.matrix(quant_all), 
                            xout=as.matrix(ef_hdi,40,1),
                            optns=list(bwReg= h_fsi, qSup=qSup, 
                                       dSup=dSup,lower=20, upper=110))

hdi_effect_data<- data.frame(Age=rep(dSup, 100), HDI.group=rep(factor(1:100),each=101),Density=as.vector(t(temp_hdi_effect$dout)),
                             HDI= rep(l, each=101))

ggplot(data = hdi_effect_data, aes(x=Age, y=Density, fill=HDI.group)) +
  geom_path(aes(colour=(HDI)), size=.9, alpha=1) +
  #ggtitle("Mortality Distributions Variation by partial HDI effect") +
  scale_color_gradientn(colours = terrain.colors(10)) +
  labs(colour = "Standard HDI") 

ggsave(file="Rplot_HDI_effect.eps")



# creating effects plot for HCE, healthcare expenditure:
l= seq(min(X_ctr[,'HCE']), max(X_ctr[,'HCE']), length.out= 100)
covar_median_hcexp_effect<- matrix(0, length(l), 5)
m=apply(X_ctr, 2, median)

covar_median_hcexp_effect<- cbind(GDP_yoy=rep(m[1],100), 
                                CO2emission=rep(m[3],100), Infantm=rep(m[4],100), HDI=rep(m[5],100))
ef_hcexp<-cbind(covar_median_hcexp_effect, HCE=l)
ef_hcexp<- ef_hcexp%*%t(thetahat)

X_ctr <- matrix(as.vector(as.matrix(X_ctr)),40, 5, byrow = F)
colnames(X_ctr) <- c('GDPC', 'CO2E', 'IM', 'HDI',"HCE")

temp_hce_effect<- LocDenReg(xin= X_ctr%*%t(thetahat), qin=as.matrix(quant_all), 
                            xout=as.matrix(ef_hcexp,40,1),
                            optns=list(bwReg= h_fsi, qSup=qSup, 
                                       dSup=dSup,lower=20, upper=110))

hce_effect_data<- data.frame(Age=rep(dSup, 100), HCE.group=rep(factor(1:100),each=101),Density=as.vector(t(temp_hce_effect$dout)),
                             HCE= rep(l, each=101))

ggplot(data = hce_effect_data, aes(x=Age, y=Density, fill=HCE.group)) +
  geom_path(aes(colour=(HCE)), size=.9, alpha=1) +
  #ggtitle("Mortality Distributions Variation by partial HCE effect") +
  scale_color_gradientn(colours = terrain.colors(10)) +
  labs(colour = "Standard HCE") 

ggsave(file="Rplot_HCE_effect.eps")




