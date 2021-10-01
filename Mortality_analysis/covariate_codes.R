
################################
# Setting up the covariate data 
################################

# creating the array of countries on which mortality data is available~
country <- c("Australia","Austria","Belarus","Belgium","Bulgaria","Canada","Chile","Croatia","Czech Republic","Denmark",
            "Estonia","Finland",'France',"Germany","Greece",'Hong Kong SAR, China',"Hungary","Iceland","Ireland","Israel",
            'Italy',"Japan","Latvia","Lithuania","Luxembourg","Netherlands",'New Zealand','Norway',"Poland","Portugal",
            'Korea, Rep.','Russian Federation','Slovak Republic',"Slovenia","Spain","Sweden","Switzerland",'Taiwan',
            'United Kingdom', 'United States','Ukraine')

length(country)

# the year for which we analyze mortality
Yr <-2013

#Setting working directory
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/LocalFrechetSpherical/Mortality_all")

#getting covariate data~
# GDP year-on-year percentage change~
gdp <- read.csv("GDP_YoY_perc.csv", stringsAsFactors =T)                 # gdp year-on-year percentage change 
gdp = subset(gdp, gdp[,c("Country.Name")] %in% country)     # selecting data for the countries considered above
country_<- gdp[,c("Country.Name")]                         # storing the counties in order they appear in covariate dataset
colnames(gdp) <- c('Country',colnames(gdp[,-1]))
gdp_yr <- gdp[,paste('X',Yr,sep ="")]                            # selecting the for given year for above countries
#View(gdp)

write.csv(country_, "Countries_FSI.csv")

# Country government expenditure on healthcare as percentage of GDP
hcexp <- read.csv("HC_exp_percGDP.csv", header=T)                    # Healthcare expenditure percentage of gdp 
hcexp = subset(hcexp, hcexp[,c("Country.Name")] %in% country_)      # selecting the for given year for above countries
colnames(hcexp) <- c('Country', colnames(hcexp[,-1]))
hcexp_yr<- hcexp[,paste('X',Yr,sep ="")]                                 # selecting data for the countries considered above
hcexp_yr[18]<- mean(hcexp_yr[-18])  # since the data correcponding to 'Hong Kong SAR, China' is missing for year 2013.
#View(hcexp)

# The CO2 emissions tonnes per capita for each country~
co2em <- read.csv("CO2_emissions_pc.csv", header=T)       # damage in USD caused by atmospheric CO2
co2em = subset(co2em, co2em[,c("Country.Name")] %in% country_) # selecting data for the countries considered above
colnames(co2em) <- c('Country', colnames(co2em[,-1]))
co2em_yr<- co2em[,paste('X',Yr,sep ="")]      # selecting the for given year for above countries
#View(co2em)

# Infant mortality per 1000 births
infmort <- read.csv('infmort.csv', header = T)
infmort = subset(infmort, infmort[,c("Country.Name")] %in% country_)       # selecting data for the countries considered above
colnames(infmort) <- c('Country', colnames(infmort[,-1]))
infmort_yr<- infmort[,paste('X',Yr,sep ="")]                               # selecting the for given year for above countries
infmort_yr[18]<- mean(infmort_yr[-18]) # the data for 'Hong Kong SAR, China' was missing hence replaced with mean


# Human Development Index data
hdi <- read.csv("HumanDevelopmentIndex(HDI).csv", header = T)  # reading the dataset for HDI                
hdi <- hdi[c(9,10,17,26,16,32,165,35,45,65,47,160,56,60,61,180,67,75,42,76,82,77,83,
             84,86,91,101,102,95,122,128,123,138,139,142,155,156,164,178,181),]  # manually selecting data for the countries considered above

colnames(hdi) <- c('Country', colnames(hdi[,-1]))
hdi_yr <- as.numeric(hdi[,paste('X',Yr, sep ="")])      # selecting the for given year for above countries

# creating the design matrix for given countries and year
X <- data.frame(GDP= gdp_yr, HCexp =hcexp_yr, CO2em=co2em_yr, Infmt= infmort_yr, HDI= hdi_yr)

# rearranging the countries in order they appear in world bank dataset
X_ctr <- scale(X) # centering and scaling the data 
X_ctr<- matrix(as.vector(X_ctr), 40 ,5, byrow = F)
row.names(X_ctr)<- country_  
colnames(X_ctr) <- c("GDP", "HCexp", "CO2em", "Infmt", "HDI")

# writing the scaled data for further use:
write.csv(X_ctr, "X_centerscale.csv")

