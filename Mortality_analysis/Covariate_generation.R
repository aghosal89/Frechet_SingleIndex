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
setwd("~/MATLAB/DistData/Mortality_all")

#setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/FSI/Mortality_all")

# the following are the list of 40 countries for which mortality is considered in our model.
country_ <- read.csv("Countries_FSI.csv", header = T)[,-1]

# Define n, the number of observations (countries)
n<- length(country_)

# the year for which we analyze mortality
Yr <-2013

# getting covariate data~
# GDP year-on-year percentage change~
gdpc <- read.csv("API_NY.GDP.MKTP.KD.ZG_DS2_en_csv_v2_4498512.csv", stringsAsFactors =T)  
# selecting data for the countries considered above
gdpc = subset(gdpc, gdpc[,1] %in% country_)  
# selecting data for the year 2013
colnames(gdpc) <-c('Country',colnames(gdpc[,-1]))
gdpc_yr <- gdpc[,paste('X',Yr,sep ="")]       

# Country government expenditure on healthcare as percentage of GDP
hce <-read.csv("API_SH.XPD.CHEX.GD.ZS_DS2_en_csv_v2_4499032.csv", header=T)  
# selecting the data for given year for above countries
hce = subset(hce, hce[,1] %in% country_)      
colnames(hce) <- c('Country', colnames(hce[,-1]))
# selecting data for the year 2013 
hce_yr<-hce[,paste('X',Yr,sep ="")]  

# The CO2 emissions tonnes per capita for each country~
co2e <- read.csv("API_EN.ATM.CO2E.PC_DS2_en_csv_v2_4498382.csv", header=T)       
# selecting the data for given year for above countries
co2e = subset(co2e, co2e[,1] %in% country_) 
colnames(co2e) <- c('Country', colnames(co2e[,-1]))
# selecting data for the year 2013 
co2e_yr<- co2e[,paste('X',Yr,sep ="")]      


# Infant mortality per 1000 births
im<- read.csv("UNICEF-CME_DF_2021_WQ-1.0-download.csv", header = TRUE)
# selecting data for the year 2013 
im_yr<- subset(im, im$REF_DATE.Reference.Date==2013.5 & 
                 im$SERIES_NAME.Series.Name=="UN_IGME: UN IGME estimate")
im_yr<- im_yr$OBS_VALUE.Observation.Value

# Human Development Index data 
hdi <- read.csv("HDR21-22_Composite_indices_complete_time_series.csv", header = T)  

country_2<- country_
country_2[9]<- "Czechia"
country_2[25]<- "Korea (Republic of)"
country_2[35]<- "Slovakia"

# uploading the dataset for HDI manually
hdi <- subset(hdi, hdi$country %in% country_2)
# selecting the for given year for above countries
hdi_yr <- as.numeric(hdi$hdi_2013)      


# creating the design matrix for given countries and year
X <- data.frame(GDPC= gdpc_yr, HCE =hce_yr, CO2E=co2e_yr, IM= im_yr, HDI= hdi_yr)
rownames(X)<- country_

# rearranging the countries in order they appear in world bank dataset
X_ctr <- scale(X) # centering and scaling the data 

# define p, the number of covariates
p<- 5

# creating the design matrix for storing and further analysis
X_ctr<- matrix(as.vector(X_ctr), n, p, byrow = F)
row.names(X_ctr)<- country_  
colnames(X_ctr) <- c("GDPC", "HCE", "CO2E", "IM", "HDI")

write.csv(X_ctr, "X_centerscale.csv")

