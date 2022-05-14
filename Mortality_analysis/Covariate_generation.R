############
# script #1
############

## Setting up the covariate data 

# following is the list of countries on which mortality data is available from mortality.org for
# which the covariate data is available from various sources as mentioned in the paper, the
# variables that had most impact on human mortality in the countries considered are: 

# 1) GDP year-on-year percentage change (GDPC)
# 2) Healthcare expenditure as percentage of GDP (HCE)
# 3) CO2 emissions in tonnes per capita (CO2E)
# 4) Infant mortality rate per 1000 live births (IM)
# 5) Human Development Index (HDI)

# Aritra's working directory
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/FSI/Mortality_all")

# the following are the list of 40 countries for which mortality is considered in our model.
country_ <- read.csv("Countries_FSI.csv", header = T)[,-1]

# Define n, the number of observations (countries)
n<- length(country_)

# the year for which we analyze mortality
Yr <-2013

# getting covariate data~
# GDP year-on-year percentage change~
gdpc <- read.csv("GDP_YoY_perc.csv", stringsAsFactors =T)  
# selecting data for the countries considered above
gdpc = subset(gdpc, gdpc[,c("Country.Name")] %in% country_)  
# selecting for given year for the year 2013
colnames(gdpc) <-c('Country',colnames(gdpc[,-1]))
gdpc_yr <- gdpc[,paste('X',Yr,sep ="")]       

# Country government expenditure on healthcare as percentage of GDP
hce <-read.csv("HC_exp_percGDP.csv", header=T)  
# selecting the data for given year for above countries
hce = subset(hce, hce[,c("Country.Name")] %in% country_)      
colnames(hce) <- c('Country', colnames(hce[,-1]))
# selecting data for the year 2013 
hce_yr<-hce[,paste('X',Yr,sep ="")]  
# Missing data imputed for 'Hong Kong SAR, China' for year 2013.
hce_yr[18]<- mean(hce_yr[-18])  

# The CO2 emissions tonnes per capita for each country~
co2e <- read.csv("CO2_emissions_pc.csv", header=T)       
# selecting the data for given year for above countries
co2e = subset(co2e, co2e[,c("Country.Name")] %in% country_) 
colnames(co2e) <- c('Country', colnames(co2e[,-1]))
# selecting data for the year 2013 
co2e_yr<- co2e[,paste('X',Yr,sep ="")]      

# Infant mortality per 1000 births
im <- read.csv('infmort.csv', header = T)
# selecting the data for given year for above countries
im = subset(im, im[,c("Country.Name")] %in% country_)       
colnames(im) <- c('Country', colnames(im[,-1]))
# selecting data for the year 2013 
im_yr<- im[,paste('X',Yr,sep ="")]
# Missing data imputed for 'Hong Kong SAR, China' for year 2013.
im_yr[18]<- mean(im_yr[-18]) 


# Human Development Index data 
hdi <- read.csv("HumanDevelopmentIndex(HDI).csv", header = T)  
# uploading the dataset for HDI manually
hdi <- hdi[c(9,10,17,26,16,32,165,35,45,65,47,160,56,60,61,180,67,75,42,76,82,77,83,
             84,86,91,101,102,95,122,128,123,138,139,142,155,156,164,178,181),]

colnames(hdi) <- c('Country', colnames(hdi[,-1]))
# selecting the for given year for above countries
hdi_yr <- as.numeric(hdi[,paste('X',Yr, sep ="")])      

# creating the design matrix for given countries and year
X <- data.frame(GDPC= gdpc_yr, HCE =hce_yr, CO2E=co2e_yr, IM= im_yr, HDI= hdi_yr)

# rearranging the countries in order they appear in world bank dataset
X_ctr <- scale(X) # centering and scaling the data 

# define p, the number of covariates
p<- 5

# creating the design matrix for storing and further analysis
X_ctr<- matrix(as.vector(X_ctr), n, p, byrow = F)
row.names(X_ctr)<- country_  
colnames(X_ctr) <- c("GDPC", "HCE", "CO2E", "IM", "HDI")

write.csv(X_ctr, "X_centerscale.csv")

