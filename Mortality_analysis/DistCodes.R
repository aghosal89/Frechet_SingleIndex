
# loading necessary library~
library('frechet') 
library('ggplot2')
library('neldermead')
library('foreach')
library('doParallel')
library('fdadensity')
library('rlang')

# creating the array of countries on which mortality data is available~
country <- c("Australia","Austria","Belarus","Belgium","Bulgaria","Canada","Chile","Croatia","Czech Republic","Denmark",
            "Estonia","Finland",'France',"Germany","Greece",'Hong Kong SAR, China',"Hungary","Iceland","Ireland","Israel",
            'Italy',"Japan","Latvia","Lithuania","Luxembourg","Netherlands",'New Zealand','Norway',"Poland","Portugal",
            'Korea, Rep.','Russian Federation','Slovak Republic',"Slovenia","Spain","Sweden","Switzerland",'Taiwan',
            'United Kingdom', 'United States','Ukraine')

length(country)

# the year for which we analyze mortality
Yr <-2013

#Setting working directory (Modified according to the user)
setwd("~/MATLAB/DistData/AllCountries_lt")

#getting covariate data~

gdp <- read.csv("GDP_YOY_perc.csv", stringsAsFactors =T)                 # gdp year-on-year percentage change 
#gdp = subset(gdp, gdp[,c("Country.Name")] %in% country)     # selecting data for the countries considered above
gdp = subset(gdp, gdp[,c("ï..Country.Name")] %in% country)     # selecting data for the countries considered above
country_<- gdp[,c("ï..Country.Name")]                         # storing the counties in order they appear in covariate dataset
gdp <- gdp[,paste('X',Yr,sep ="")]                            # selecting the for given year for above countries

unemp <- read.csv("Unemp_perc.csv", header=T)                     # unemplyment percentage data 
unemp = subset(unemp, unemp[,c("ï..Country.Name")] %in% country)  # selecting data for the countries considered above
unemp<- unemp[,paste('X',Yr,sep ="")]                             # selecting the for given year for above countries

hcexp <- read.csv("HC_exp_percGDP.csv", header=T)                    # Healthcare expenditure percentage of gdp 
hcexp = subset(hcexp, hcexp[,c("ï..Country.Name")] %in% country)      # selecting the for given year for above countries
hcexp<- hcexp[,paste('X',Yr,sep ="")]                                 # selecting data for the countries considered above
hcexp[18]<- mean(hcexp[-18])

co2d <- read.csv("CO2_damage.csv", header=T)                      # damage in USD caused by atmospheric CO2
co2d = subset(co2d, co2d[,c("ï..Country.Name")] %in% country)       # selecting data for the countries considered above
co2d<- co2d[,paste('X',Yr,sep ="")]                               # selecting the for given year for above countries

infmort <- read.csv("infmort.csv", header = T)                          # death of infants under age 1 year for every 1000 live births
infmort <- subset(infmort, infmort[,c("Country.Name")] %in% country)    # selecting data for the countries considered above
infmort <- infmort[,paste('X',Yr, sep ="")]                              # selecting the for given year for above countries
infmort[18]<- mean(infmort[-18])

# creating the design matrix for given countries and year
cvr <- data.frame(cbind(gdp, unemp, hcexp, co2d, infmort))

# rearranging the countries in order they appear in world bank dataset
row.names(cvr)<- country_  
colnames(cvr)<- c("GDP.YoY.%","Unemp.%","HC.Expen.%of.GDP","CO2damage.USD","Infant.Mort.1000.births")

# The function to create the data vectors for computation of quantiles/density.
# each age is repeated the number of times as the deaths observed at that age.
data_vec<- function(x,f) {
  l<- 0    # initiate the vector by 0                  
  for (i in 1:length(x)) {        
    l<- c(l, rep(x[i], f[i]))   # repeat the age for number of deaths that year and create vector in the increasing order
  }
  l=l[-1] # remove the 0 as the first element 
  return(as.numeric(l))
}

#The data obtained from https://mortality.org/ is such that age of individuals alive beyond 110 years is recorded as 110+. 
# the following function removes the + sign at age 110 for further numerical analysis.
Remove_plus <- function(data) {
  n = nrow(data)  
  for (i in 1:n) {
    if(data[i, "Age"]=="110+") {
      data[i, "Age"]=110
    }
  }
  data$Age<- as.numeric(data$Age)
  return(data)
}

#################
# Data analysis #
#################

# in the following segment we gather the response data by the following steps:
#  1) read the lifetable corresponding to a specific country
#  2) remove '+' sign from age '110+' category using 'Remove_plus' function above.
#  3) gather the 'age-at-death' for age range [20,110].
#  4) subset data for the chosen year. 
#  5) create dataset for gathering quantile function using 'data_vec' function above


# Mortality of Australia:
lt_australia <- read.csv("lt_australia.csv", header = T) #Reading the mortality data
lt_australia<- Remove_plus(lt_australia)
lt_australia <- subset(lt_australia, lt_australia$Age %in% seq(20,110)) #subsetting data by age-group
lt_australia <- lt_australia[,c("Year", "Age", "dx")] #subset the death data ~
lt_australia<- subset(lt_australia, lt_australia[,1] %in% Yr)
unique(lt_australia$Year)
#View(lt_australia)

myls_australia = data_vec(lt_australia[,"Age"], lt_australia[,"dx"])

# Mortality of Austria:
lt_austria <- read.csv("lt_austria.csv", header = T) #Reading the mortality data
lt_austria<- Remove_plus(lt_austria)
lt_austria <- subset(lt_austria, lt_austria$Age %in% seq(20,110)) #subsetting data by age-group
lt_austria <- lt_austria[,c("Year", "Age", "dx")] #subset the death data ~
lt_austria<- subset(lt_austria, lt_austria[,1] %in% Yr)
unique(lt_austria$Year)
#View(lt_austria)

myls_austria <- data_vec(lt_austria[,"Age"], lt_austria[,"dx"])  

# Mortality of Belarus:
lt_belarus <- read.csv("lt_belarus.csv", header = T) #Reading the mortality data
lt_belarus<- Remove_plus(lt_belarus)
lt_belarus <- subset(lt_belarus, lt_belarus$Age %in% seq(20,110)) #subsetting data by age-group
lt_belarus <- lt_belarus[,c("Year", "Age", "dx")] #subset the death data ~
lt_belarus<- subset(lt_belarus, lt_belarus[,1] %in% Yr)
#View(lt_belarus)
unique(lt_belarus$Year)

myls_belarus <- data_vec(lt_belarus[,"Age"], lt_belarus[,"dx"])  

# Mortality of Belgium:
lt_belgium <- read.csv("lt_belgium.csv", header = T) #Reading the mortality data
lt_belgium<- Remove_plus(lt_belgium)
lt_belgium <- subset(lt_belgium, lt_belgium$Age %in% seq(20,110)) #subsetting data by age-group
lt_belgium <- lt_belgium[,c("Year", "Age", "dx")] #subset the death data ~
lt_belgium<- subset(lt_belgium, lt_belgium[,1] %in% Yr)
#View(lt_belgium)
unique(lt_belgium$Year)

myls_belgium <- data_vec(lt_belgium[,"Age"], lt_belgium[,"dx"])  

# Mortality of bulgaria:
lt_bulgaria <- read.csv("lt_bulgaria.csv", header = T) #Reading the mortality data
lt_bulgaria<- Remove_plus(lt_bulgaria)
lt_bulgaria <- subset(lt_bulgaria, lt_bulgaria$Age %in% seq(20,110)) #subsetting data by age-group
lt_bulgaria <- lt_bulgaria[,c("Year", "Age", "dx")] #subset the death data ~
lt_bulgaria<- subset(lt_bulgaria, lt_bulgaria[,1] %in% Yr)
unique(lt_bulgaria$Year)
#View(lt_bulgaria)

myls_bulgaria <- data_vec(lt_bulgaria[,"Age"], lt_bulgaria[,"dx"])  

# Mortality of canada:
lt_canada <- read.csv("lt_canada.csv", header = T) #Reading the mortality data
lt_canada<- Remove_plus(lt_canada)
lt_canada <- subset(lt_canada, lt_canada$Age %in% seq(20,110)) #subsetting data by age-group
lt_canada <- lt_canada[,c("Year", "Age", "dx")] #subset the death data ~
lt_canada<- subset(lt_canada, lt_canada[,1] %in% Yr)
unique(lt_canada$Year)
#View(lt_canada)

myls_canada <- data_vec(lt_canada[,"Age"], lt_canada[,"dx"])  


# Mortality of Chile:
lt_chile <- read.csv("lt_chile.csv", header = T) #Reading the mortality data
lt_chile<- Remove_plus(lt_chile)
lt_chile <- subset(lt_chile, lt_chile$Age %in% seq(20,110)) #subsetting data by age-group
lt_chile <- lt_chile[,c("Year", "Age", "dx")] #subset the death data ~
lt_chile<- subset(lt_chile, lt_chile[,1] %in% Yr)
unique(lt_chile$Year)
#View(lt_chile)

myls_chile <- data_vec(lt_chile[,"Age"], lt_chile[,"dx"])  


# Mortality of Croatia:
lt_croatia <- read.csv("lt_croatia.csv", header = T) #Reading the mortality data
lt_croatia<- Remove_plus(lt_croatia)
lt_croatia <- subset(lt_croatia, lt_croatia$Age %in% seq(20,110)) #subsetting data by age-group
lt_croatia <- lt_croatia[,c("Year", "Age", "dx")] #subset the death data ~
lt_croatia<- subset(lt_croatia, lt_croatia[,1] %in% Yr)
unique(lt_croatia$Year)
#View(lt_croatia)

myls_croatia <- data_vec(lt_croatia[,"Age"], lt_croatia[,"dx"])  



# Mortality of Czech Republic:
lt_czech <- read.csv("lt_czech.csv", header = T) #Reading the mortality data
lt_czech<- Remove_plus(lt_czech)
lt_czech <- subset(lt_czech, lt_czech$Age %in% seq(20,110)) #subsetting data by age-group
lt_czech <- lt_czech[,c("Year", "Age", "dx")] #subset the death data ~
lt_czech<- subset(lt_czech, lt_czech[,1] %in% Yr)
unique(lt_czech$Year)
#View(lt_czech)

myls_czech <- data_vec(lt_czech[,"Age"], lt_czech[,"dx"])  


# Mortality of denmark:
lt_denmark <- read.csv("lt_denmark.csv", header = T) #Reading the mortality data
lt_denmark<- Remove_plus(lt_denmark)
lt_denmark <- subset(lt_denmark, lt_denmark$Age %in% seq(20,110)) #subsetting data by age-group
lt_denmark <- lt_denmark[,c("Year", "Age", "dx")] #subset the death data ~
lt_denmark<- subset(lt_denmark, lt_denmark[,1] %in% Yr)
unique(lt_denmark$Year)
#View(lt_denmark)

myls_denmark <- data_vec(lt_denmark[,"Age"], lt_denmark[,"dx"])  


# Mortality of estonia:
lt_estonia <- read.csv("lt_estonia.csv", header = T) #Reading the mortality data
lt_estonia<- Remove_plus(lt_estonia)
lt_estonia <- subset(lt_estonia, lt_estonia$Age %in% seq(20,110)) #subsetting data by age-group
lt_estonia <- lt_estonia[,c("Year", "Age", "dx")] #subset the death data ~
lt_estonia<- subset(lt_estonia, lt_estonia[,1] %in% Yr)
unique(lt_estonia$Year)
#View(lt_estonia)

myls_estonia <- data_vec(lt_estonia[,"Age"], lt_estonia[,"dx"])  

# Mortality of Finland:
lt_finland <- read.csv("lt_finland.csv", header = T) #Reading the mortality data
lt_finland<- Remove_plus(lt_finland)
lt_finland <- subset(lt_finland, lt_finland$Age %in% seq(20,110)) #subsetting data by age-group
lt_finland <- lt_finland[,c("Year", "Age", "dx")] #subset the death data ~
lt_finland<- subset(lt_finland, lt_finland[,1] %in% Yr)
unique(lt_finland$Year)
#View(lt_finland)

myls_finland <- data_vec(lt_finland[,"Age"], lt_finland[,"dx"])  

# Mortality of france:
lt_france <- read.csv("lt_france.csv", header = T) #Reading the mortality data
lt_france<- Remove_plus(lt_france)
lt_france <- subset(lt_france, lt_france$Age %in% seq(20,110)) #subsetting data by age-group
lt_france <- lt_france[,c("Year", "Age", "dx")] #subset the death data ~
lt_france<- subset(lt_france, lt_france[,1] %in% Yr)
unique(lt_france$Year)
#View(lt_france)

myls_france <- data_vec(lt_france[,"Age"], lt_france[,"dx"])  

# Mortality of denmark:
lt_germany <- read.csv("lt_germany.csv", header = T) #Reading the mortality data
lt_germany<- Remove_plus(lt_germany)
lt_germany <- subset(lt_germany, lt_germany$Age %in% seq(20,110)) #subsetting data by age-group
lt_germany <- lt_germany[,c("Year", "Age", "dx")] #subset the death data ~
lt_germany<- subset(lt_germany, lt_germany[,1] %in% Yr)
unique(lt_germany$Year)
#View(lt_germany)

myls_germany <- data_vec(lt_germany[,"Age"], lt_germany[,"dx"])  


# Mortality of greece:
lt_greece <- read.csv("lt_greece.csv", header = T) #Reading the mortality data
lt_greece<- Remove_plus(lt_greece)
lt_greece <- subset(lt_greece, lt_greece$Age %in% seq(20,110)) #subsetting data by age-group
lt_greece <- lt_greece[,c("Year", "Age", "dx")] #subset the death data ~
lt_greece<- subset(lt_greece, lt_greece[,1] %in% Yr)
unique(lt_greece$Year)
#View(lt_greece)

myls_greece <- data_vec(lt_greece[,"Age"], lt_greece[,"dx"])  

# Mortality of hongkong:
lt_hongkong <- read.csv("lt_hongkong.csv", header = T) #Reading the mortality data
lt_hongkong<- Remove_plus(lt_hongkong)
lt_hongkong <- subset(lt_hongkong, lt_hongkong$Age %in% seq(20,110)) #subsetting data by age-group
lt_hongkong <- lt_hongkong[,c("Year", "Age", "dx")] #subset the death data ~
lt_hongkong<- subset(lt_hongkong, lt_hongkong[,1] %in% Yr)
unique(lt_hongkong$Year)
#View(lt_hongkong)

myls_hongkong <- data_vec(lt_hongkong[,"Age"], lt_hongkong[,"dx"])  


# Mortality of hungary:
lt_hungary<- read.csv("lt_hungary.csv", header = T) #Reading the mortality data
lt_hungary<- Remove_plus(lt_hungary)
lt_hungary <- subset(lt_hungary, lt_hungary$Age %in% seq(20,110)) #subsetting data by age-group
lt_hungary <- lt_hungary[,c("Year", "Age", "dx")] #subset the death data ~
lt_hungary<- subset(lt_hungary, lt_hungary[,1] %in% Yr)
unique(lt_hungary$Year)
#View(lt_hungary)

myls_hungary <- data_vec(lt_hungary[,"Age"], lt_hungary[,"dx"])  


# Mortality of iceland:
lt_iceland<- read.csv("lt_iceland.csv", header = T) #Reading the mortality data
lt_iceland<- Remove_plus(lt_iceland)
lt_iceland <- subset(lt_iceland, lt_iceland$Age %in% seq(20,110)) #subsetting data by age-group
lt_iceland <- lt_iceland[,c("Year", "Age", "dx")] #subset the death data ~
lt_iceland<- subset(lt_iceland, lt_iceland[,1] %in% Yr)
unique(lt_iceland$Year)
#View(lt_iceland)

myls_iceland <- data_vec(lt_iceland[,"Age"], lt_iceland[,"dx"])


# Mortality of ireland:
lt_ireland<- read.csv("lt_ireland.csv", header = T) #Reading the mortality data
lt_ireland<- Remove_plus(lt_ireland)
lt_ireland <- subset(lt_ireland, lt_ireland$Age %in% seq(20,110)) #subsetting data by age-group
lt_ireland <- lt_ireland[,c("Year", "Age", "dx")] #subset the death data ~
lt_ireland<- subset(lt_ireland, lt_ireland[,1] %in% Yr)
unique(lt_ireland$Year)
#View(lt_ireland)

myls_ireland<- data_vec(lt_ireland[,"Age"], lt_ireland[,"dx"])  


# Mortality of israel:
lt_israel<- read.csv("lt_israel.csv", header = T) #Reading the mortality data
lt_israel<- Remove_plus(lt_israel)
lt_israel <- subset(lt_israel, lt_israel$Age %in% seq(20,110)) #subsetting data by age-group
lt_israel <- lt_israel[,c("Year", "Age", "dx")] #subset the death data ~
lt_israel<- subset(lt_israel, lt_israel[,1] %in% Yr)
unique(lt_israel$Year)
#View(lt_israel)

myls_israel<- data_vec(lt_israel[,"Age"], lt_israel[,"dx"])  


# Mortality of italy:
lt_italy<- read.csv("lt_italy.csv", header = T) #Reading the mortality data
lt_italy<- Remove_plus(lt_italy)
lt_italy <- subset(lt_italy, lt_italy$Age %in% seq(20,110)) #subsetting data by age-group
lt_italy <- lt_italy[,c("Year", "Age", "dx")] #subset the death data ~
lt_italy<- subset(lt_italy, lt_italy[,1] %in% Yr)
unique(lt_italy$Year)
#View(lt_italy)

myls_italy<- data_vec(lt_italy[,"Age"], lt_italy[,"dx"])  


# Mortality of japan:
lt_japan <- read.csv("lt_japan.csv", header = T) #Reading the mortality data
lt_japan<- Remove_plus(lt_japan)
lt_japan <- subset(lt_japan, lt_japan$Age %in% seq(20,110)) #subsetting data by age-group
lt_japan <- lt_japan[,c("Year", "Age", "dx")] #subset the death data ~
lt_japan<- subset(lt_japan, lt_japan[,1] %in% Yr)
unique(lt_japan$Year)
#View(lt_japan)

myls_japan <- data_vec(lt_japan[,"Age"], lt_japan[,"dx"])  


# Mortality of latvia:
lt_latvia <- read.csv("lt_latvia.csv", header = T) #Reading the mortality data
lt_latvia<- Remove_plus(lt_latvia)
lt_latvia <- subset(lt_latvia, lt_latvia$Age %in% seq(20,110)) #subsetting data by age-group
lt_latvia <- lt_latvia[,c("Year", "Age", "dx")] #subset the death data ~
lt_latvia<- subset(lt_latvia, lt_latvia[,1] %in% Yr)
unique(lt_latvia$Year)
#View(lt_latvia)

myls_latvia <- data_vec(lt_latvia[,"Age"], lt_latvia[,"dx"])  


# Mortality of lithuania:
lt_lithuania<- read.csv("lt_lithuania.csv", header = T) #Reading the mortality data
lt_lithuania<- Remove_plus(lt_lithuania)
lt_lithuania <- subset(lt_lithuania, lt_lithuania$Age %in% seq(20,110)) #subsetting data by age-group
lt_lithuania <- lt_lithuania[,c("Year", "Age", "dx")] #subset the death data ~
lt_lithuania<- subset(lt_lithuania, lt_lithuania[,1] %in% Yr)
unique(lt_lithuania$Year)
#View(lt_lithuania)

myls_lithuania<- data_vec(lt_lithuania[,"Age"], lt_lithuania[,"dx"])  


# Mortality of luxembourg:
lt_luxembourg<- read.csv("lt_luxembourg.csv", header = T) #Reading the mortality data
lt_luxembourg<- Remove_plus(lt_luxembourg)
lt_luxembourg <- subset(lt_luxembourg, lt_luxembourg$Age %in% seq(20,110)) #subsetting data by age-group
lt_luxembourg <- lt_luxembourg[,c("Year", "Age", "dx")] #subset the death data ~
lt_luxembourg<- subset(lt_luxembourg, lt_luxembourg[,1] %in% Yr)
unique(lt_luxembourg$Year)
#View(lt_luxembourg)

myls_luxembourg<- data_vec(lt_luxembourg[,"Age"], lt_luxembourg[,"dx"])  


# Mortality of netherlands:
lt_netherlands<- read.csv("lt_netherlands.csv", header = T) #Reading the mortality data
lt_netherlands<- Remove_plus(lt_netherlands)
lt_netherlands <- subset(lt_netherlands, lt_netherlands$Age %in% seq(20,110)) #subsetting data by age-group
lt_netherlands <- lt_netherlands[,c("Year", "Age", "dx")] #subset the death data ~
lt_netherlands<- subset(lt_netherlands, lt_netherlands[,1] %in% Yr)
unique(lt_netherlands$Year)
#View(lt_netherlands)

myls_netherlands<- data_vec(lt_netherlands[,"Age"], lt_netherlands[,"dx"])  

# Mortality of newz:
lt_newz<- read.csv("lt_newz.csv", header = T) #Reading the mortality data
lt_newz<- Remove_plus(lt_newz)
lt_newz <- subset(lt_newz, lt_newz$Age %in% seq(20,110)) #subsetting data by age-group
lt_newz <- lt_newz[,c("Year", "Age", "dx")] #subset the death data ~
lt_newz<- subset(lt_newz, lt_newz[,1] %in% Yr)
unique(lt_newz$Year)
#View(lt_newz)

myls_newz<- data_vec(lt_newz[,"Age"], lt_newz[,"dx"])  

# Mortality of Norway:
lt_norway<- read.csv("lt_norway.csv", header = T) #Reading the mortality data
lt_norway<- Remove_plus(lt_norway)
lt_norway <- subset(lt_norway, lt_norway$Age %in% seq(20,110)) #subsetting data by age-group
lt_norway <- lt_norway[,c("Year", "Age", "dx")] #subset the death data ~
lt_norway<- subset(lt_norway, lt_norway[,1] %in% Yr)
unique(lt_norway$Year)
#View(lt_norway)

myls_norway<- data_vec(lt_norway[,"Age"], lt_norway[,"dx"])  


# Mortality of Poland:
lt_poland<- read.csv("lt_poland.csv", header = T) #Reading the mortality data
lt_poland<- Remove_plus(lt_poland)
lt_poland <- subset(lt_poland, lt_poland$Age %in% seq(20,110)) #subsetting data by age-group
lt_poland <- lt_poland[,c("Year", "Age", "dx")] #subset the death data ~
lt_poland<- subset(lt_poland, lt_poland[,1] %in% Yr)
unique(lt_poland$Year)
#View(lt_poland)

myls_poland<- data_vec(lt_poland[,"Age"], lt_poland[,"dx"])  


# Mortality of Portugal:
lt_portugal<- read.csv("lt_portugal.csv", header = T) #Reading the mortality data
lt_portugal<- Remove_plus(lt_portugal)
lt_portugal <- subset(lt_portugal, lt_portugal$Age %in% seq(20,110)) #subsetting data by age-group
lt_portugal <- lt_portugal[,c("Year", "Age", "dx")] #subset the death data ~
lt_portugal<- subset(lt_portugal, lt_portugal[,1] %in% Yr)
unique(lt_portugal$Year)
#View(lt_portugal)

myls_portugal<- data_vec(lt_portugal[,"Age"], lt_portugal[,"dx"])  


# Mortality of korea:
lt_korea<- read.csv("lt_korea.csv", header = T) #Reading the mortality data
lt_korea<- Remove_plus(lt_korea)
lt_korea <- subset(lt_korea, lt_korea$Age %in% seq(20,110)) #subsetting data by age-group
lt_korea <- lt_korea[,c("Year", "Age", "dx")] #subset the death data ~
lt_korea<- subset(lt_korea, lt_korea[,1] %in% Yr)
unique(lt_korea$Year)
#View(lt_korea)

myls_korea<- data_vec(lt_korea[,"Age"], lt_korea[,"dx"])  


# Mortality of russia:
lt_russia<- read.csv("lt_russia.csv", header = T) #Reading the mortality data
lt_russia<- Remove_plus(lt_russia)
lt_russia <- subset(lt_russia, lt_russia$Age %in% seq(20,110)) #subsetting data by age-group
lt_russia <- lt_russia[,c("Year", "Age", "dx")] #subset the death data ~
lt_russia<- subset(lt_russia, lt_russia[,1] %in% Yr)
unique(lt_russia$Year)
#View(lt_russia)

myls_russia<- data_vec(lt_russia[,"Age"], lt_russia[,"dx"])  


# Mortality of slovakia:
lt_slovakia<- read.csv("lt_slovakia.csv", header = T) #Reading the mortality data
lt_slovakia<- Remove_plus(lt_slovakia)
lt_slovakia <- subset(lt_slovakia, lt_slovakia$Age %in% seq(20,110)) #subsetting data by age-group
lt_slovakia <- lt_slovakia[,c("Year", "Age", "dx")] #subset the death data ~
lt_slovakia<- subset(lt_slovakia, lt_slovakia[,1] %in% Yr)
unique(lt_slovakia$Year)
#View(lt_slovakia)

myls_slovakia<- data_vec(lt_slovakia[,"Age"], lt_slovakia[,"dx"])  


# Mortality of slovenia:
lt_slovenia<- read.csv("lt_slovenia.csv", header = T) #Reading the mortality data
lt_slovenia<- Remove_plus(lt_slovenia)
lt_slovenia <- subset(lt_slovenia, lt_slovenia$Age %in% seq(20,110)) #subsetting data by age-group
lt_slovenia <- lt_slovenia[,c("Year", "Age", "dx")] #subset the death data ~
lt_slovenia<- subset(lt_slovenia, lt_slovenia[,1] %in% Yr)
unique(lt_slovenia$Year)
#View(lt_slovenia)

myls_slovenia<- data_vec(lt_slovenia[,"Age"], lt_slovenia[,"dx"])  


# Mortality of spain:
lt_spain<- read.csv("lt_spain.csv", header = T) #Reading the mortality data
lt_spain<- Remove_plus(lt_spain)
lt_spain <- subset(lt_spain, lt_spain$Age %in% seq(20,110)) #subsetting data by age-group
lt_spain <- lt_spain[,c("Year", "Age", "dx")] #subset the death data ~
lt_spain<- subset(lt_spain, lt_spain[,1] %in% Yr)
unique(lt_spain$Year)
#View(lt_spain)

myls_spain<- data_vec(lt_spain[,"Age"], lt_spain[,"dx"])  


# Mortality of sweden:
lt_sweden<- read.csv("lt_sweden.csv", header = T) #Reading the mortality data
lt_sweden<- Remove_plus(lt_sweden)
lt_sweden <- subset(lt_sweden, lt_sweden$Age %in% seq(20,110)) #subsetting data by age-group
lt_sweden <- lt_sweden[,c("Year", "Age", "dx")] #subset the death data ~
lt_sweden<- subset(lt_sweden, lt_sweden[,1] %in% Yr)
unique(lt_sweden$Year)
#View(lt_sweden)

myls_sweden<- data_vec(lt_sweden[,"Age"], lt_sweden[,"dx"])  


# Mortality of switzerland:
lt_switzerland<- read.csv("lt_switzerland.csv", header = T) #Reading the mortality data
lt_switzerland<- Remove_plus(lt_switzerland)
lt_switzerland <- subset(lt_switzerland, lt_switzerland$Age %in% seq(20,110)) #subsetting data by age-group
lt_switzerland <- lt_switzerland[,c("Year", "Age", "dx")] #subset the death data ~
lt_switzerland<- subset(lt_switzerland, lt_switzerland[,1] %in% Yr)
unique(lt_switzerland$Year)
#View(lt_switzerland)

myls_switzerland<- data_vec(lt_switzerland[,"Age"], lt_switzerland[,"dx"])  


# Mortality of taiwan:
lt_taiwan<- read.csv("lt_Taiwan.csv", header = T) #Reading the mortality data
lt_taiwan<- Remove_plus(lt_taiwan)
lt_taiwan <- subset(lt_taiwan, lt_taiwan$Age %in% seq(20,110)) #subsetting data by age-group
lt_taiwan <- lt_taiwan[,c("Year", "Age", "dx")] #subset the death data ~
lt_taiwan<- subset(lt_taiwan, lt_taiwan[,1] %in% Yr)
unique(lt_taiwan$Year)
#View(lt_taiwan)

myls_taiwan<- data_vec(lt_taiwan[,"Age"], lt_taiwan[,"dx"])  


# Mortality of UK:
lt_UK<- read.csv("lt_UK.csv", header = T) #Reading the mortality data
lt_UK<- Remove_plus(lt_UK)
lt_UK <- subset(lt_UK, lt_UK$Age %in% seq(20,110)) #subsetting data by age-group
lt_UK <- lt_UK[,c("Year", "Age", "dx")] #subset the death data ~
lt_UK<- subset(lt_UK, lt_UK[,1] %in% Yr)
unique(lt_UK$Year)
#View(lt_UK)

myls_UK<- data_vec(lt_UK[,"Age"], lt_UK[,"dx"])  


# Mortality of USA:
lt_USA<- read.csv("lt_USA.csv", header = T) #Reading the mortality data
lt_USA<- Remove_plus(lt_USA)
lt_USA <- subset(lt_USA, lt_USA$Age %in% seq(20,110)) #subsetting data by age-group
lt_USA <- lt_USA[,c("Year", "Age", "dx")] #subset the death data ~
lt_USA<- subset(lt_USA, lt_USA[,1] %in% Yr)
unique(lt_USA$Year)
#View(lt_USA)

myls_USA<- data_vec(lt_USA[,"Age"], lt_USA[,"dx"])  


# Mortality of ukraine:
lt_ukraine<- read.csv("lt_Ukraine.csv", header = T) #Reading the mortality data
lt_ukraine<- Remove_plus(lt_ukraine)
lt_ukraine <- subset(lt_ukraine, lt_ukraine$Age %in% seq(20,110)) #subsetting data by age-group
lt_ukraine <- lt_ukraine[,c("Year", "Age", "dx")] #subset the death data ~
lt_ukraine<- subset(lt_ukraine, lt_ukraine[,1] %in% Yr)
unique(lt_ukraine$Year)
#View(lt_ukraine)

myls_ukraine<- data_vec(lt_ukraine[,"Age"], lt_ukraine[,"dx"])  

qSup = seq(0,1, length.out = 101)  # the grid of quantiles for mortality distribution 
qin = qSup 
lq<- length(qSup)
dSup<- seq(20,110, length.out = 101)


#################################
# Main Algorithm for computation
#################################

# 2D vectors on unit circle using polar coordinates

## Euclidean norm ##
e_norm <- function(x) {
  return(sqrt(sum(x^2)))
}

#> country_

myls <- list()
#  "Australia"                            "Austria"                    "Belgium"              
myls$australia <- myls_australia; myls$austria <- myls_austria; myls$belgium <- myls_belgium
#  "Bulgaria"                         "Belarus"                     "Canada"               
myls$bulgaria <- myls_bulgaria; myls$belarus<- myls_belarus; myls$canada <- myls_canada
#  "Switzerland"                            "Chile"             "Czech Republic"       
myls$switzerland <- myls_switzerland;myls$chile <- myls_chile; myls$czech <- myls_czech
#  "Germany"                      "Denmark"                        "Spain"            
myls$germany <- myls_germany; myls$denmark <- myls_denmark; myls$spain <- myls_spain
#  "Estonia"                            "Finland"               "France"              
myls$estonia <- myls_estonia; myls$finland <- myls_finland; myls$france <- myls_france
#  "United Kingdom"         "Greece"                "Hong Kong SAR, China"               
myls$uk <- myls_UK; myls$greece <- myls_greece; myls$hongkong <- myls_hongkong
#  "Croatia"                     "Hungary"          "Ireland"                      
myls$croatia <- myls_croatia; myls$hungary <- myls_hungary; myls$ireland <- myls_ireland
#   "Iceland"                   "Israel"               "Italy"                              
myls$iceland <- myls_iceland; myls$israel <- myls_israel; myls$italy <- myls_italy
#  "Japan"                       "Korea, Rep."          "Lithuania"                   
myls$japan <- myls_japan; myls$korea <- myls_korea; myls$lithuania <- myls_lithuania
#  "Luxembourg"                         "Latvia"               "Netherlands"                      
myls$luxembourg <- myls_luxembourg; myls$latvia <- myls_latvia; myls$netherlands <- myls_netherlands
#  "Norway"                "New Zealand"                "Poland"                           
myls$norway <- myls_norway; myls$newz <- myls_newz; myls$poland <- myls_poland;
# "Portugal"                   "Russian Federation"           "Slovak Republic"              
myls$portugal <- myls_portugal; myls$russia <- myls_russia; myls$slovakia <- myls_slovakia;
#  "Slovenia"                        "Sweden"                   "Ukraine"              
myls$slovenia <- myls_slovenia; myls$sweden <- myls_sweden; myls$ukraine <- myls_ukraine; 
#  "United States"
myls$usa <- myls_USA;


#density_all <- vector(mode = 'list', length = 40)
density_all <- matrix(NA, 40, 101)
for (i in 1:40) {
  density_all[i,] <- CreateDensity(y= myls[[i]], optns= list(nRegGrid = 101))$y
}

quant_all <- matrix(NA, nrow = 40, ncol = length(qSup))
for (i in 1:40) {
  quant_all[i,] <- dens2quantile(dens= density_all[i,], dSup=dSup, qSup = qSup, useSplines = F)
}


### The following function computes cartesian coordinates from polar coordinates
# Inputs:  eta: angles in radian
#            r: radius

polar2cart <- function(eta, r) {
  s <- length(eta)
  theta <- vector(length = s+1)
  
  if(s==4) {
    theta[1]<- r*cos(eta[1])*cos(eta[2])*cos(eta[3])*cos(eta[4])
    theta[2]<- r*cos(eta[1])*cos(eta[2])*cos(eta[3])*sin(eta[4])
    theta[3]<- r*cos(eta[1])*cos(eta[2])*sin(eta[3])
    theta[4]<- r*cos(eta[1])*sin(eta[2])
    theta[5]<- r*sin(eta[1])
    return(theta)
  }
  
  if(s==3) {
  theta[1]<- r*cos(eta[1])*cos(eta[2])*cos(eta[3])
  theta[2]<- r*cos(eta[1])*cos(eta[2])*sin(eta[3])
  theta[3]<- r*cos(eta[1])*sin(eta[2])
  theta[4]<- r*sin(eta[1])
  return(theta)
  
  }
  else if(s==2) {
    theta[1]<- r*cos(eta[1])*cos(eta[2])
    theta[2]<- r*cos(eta[1])*sin(eta[2])
    theta[3]<- r*sin(eta[1])
    return(theta)
  }
  
  else if(s==1) {
    theta[1]<- r*cos(eta[1])
    theta[2]<- r*sin(eta[1])
    return(theta)
  }
  
}


# The following function computes polar coordinates from the cartesian 
# INPUTS:  x is the cartesian coordinate

cart2polar <- function(x) {
  p= length(x)                  # number of covariates
  eta =  vector(length = p-1)   # vector for storing the radian angles
  r <- sqrt(sum(x^2))         # find the norm
  if(p==2) {
    eta[1]= atan(x[2]/x[1])
  }
  
  if(p==4) {
    eta[1] = atan(x[4]/sqrt(sum((c(x[1],x[2],x[3]))^2)))
    eta[2] = atan(x[3]/sqrt(sum((c(x[1],x[2]))^2)))
    eta[3] = atan(x[2]/x[1])
  }
  
  if(p==5) {
    eta[1] = atan(x[5]/sqrt(sum((c(x[1],x[2],x[3],x[4]))^2)))
    eta[2] = atan(x[4]/sqrt(sum((c(x[1],x[2],x[3]))^2)))
    eta[3] = atan(x[3]/sqrt(sum((c(x[1],x[2]))^2)))
    eta[4] = atan(x[2]/x[1])
  }
  return(c(r,eta))
}


######### Local Frechet on each covariate ###########
X_ctr <- scale(as.matrix(cvr)) # centering and scaling the data 

# for gdp YoY % change
lf_gdp<- LocDenReg(xin=as.vector(X_ctr[,"GDP.YoY.%"]), yin=myls, optns =list(qSup = qSup,
                                                               lower=20, upper=110))
metric_v <- matrix(0,nrow(lf_gdp$qin), 1)
for(i in 1:nrow(lf_gdp$qin)) {
  metric_v[i,] <- W_metric(lf_gdp$qout[i,], lf_gdp$qin[i,], qSup)
}

mean(metric_v)
# 2.455191
# 2.437389

# for unemployment %
lf_unemp<- LocDenReg(xin=as.vector(X_ctr[,"Unemp.%"]), yin=myls, optns= list(qSup = qSup,
                                                                           lower=20, upper=110, infSupport=F))
metric_v <- matrix(0,nrow(lf_unemp$qin),1 )
for(i in 1:nrow(lf_unemp$qin)) {
  metric_v[i,] <- W_metric(lf_unemp$qout[i,], lf_unemp$qin[i,], qSup)
}
mean(metric_v)
# 2.305016
# 2.310186

# for health care expenditure %
lf_hcexp<- LocDenReg(xin=as.vector(X_ctr[,"HC.Expen.%of.GDP"]), yin=myls, optns= list(qSup = qSup,
                                                                           lower=20, upper=110))
metric_v <- matrix(0,nrow(lf_hcexp$qin),1 )
for(i in 1:nrow(lf_hcexp$qin)) {
  metric_v[i,] <- W_metric(lf_hcexp$qout[i,], lf_hcexp$qin[i,], qSup)
}

mean(metric_v)
# 1.766682

# for CO2 damage %
lf_co2d<- LocDenReg(xin=as.vector(X_ctr[,"CO2damage.USD"]), yin=myls, optns= list(qSup = qSup, lower=20, upper=110))

metric_v <- matrix(0, nrow(lf_co2d$qin),1)
for(i in 1:nrow(lf_co2d$qin)) {
  metric_v[i,] <- W_metric(lf_co2d$qout[i,], lf_co2d$qin[i,], qSup)
}
mean(metric_v)
# 2.630886

# for Infant mortality per 1000 births
lf_infmort<-LocDenReg(xin=as.vector(X_ctr[,"Infant.Mort.1000.births"]), yin=myls, optns=list(qSup=qSup, lower=20, upper=110))

metric_v <- matrix(0, nrow(lf_infmort$qin), 1)
for(i in 1:nrow(lf_infmort$qin)) {
  metric_v[i,] <- W_metric(lf_infmort$qout[i,], lf_infmort$qin[i,], qSup)
}

mean(metric_v)
# 1.831657


####### Global Frechet on all covariates ########
gf_all <-GloDenReg(xin = as.matrix(X_ctr), yin= myls,optns= list(qSup = qSup, lower=20, upper=110))

metric_v <- matrix(0, nrow(gf_all$qin),1)
for(i in 1:nrow(gf_all$qin)) {
  metric_v[i,] <- W_metric(gf_all$qout[i,], gf_all$qin[i,], qSup)
}

mean(metric_v)
# 1.516844
# 1.513268

# creating cost function for nelder-mead algorithm 
# Inputs:
# 1) Y_list: a list of (Response corresponding to single year, given countries)
# 2) x: nxp design matrix (covariate values corresponding to a single year)
# 3) eta: coordinates of the point in radian
# 4) qsup: equispaced grid of points on support of distribution of cdf

f<- seq(-pi+(pi/15), pi-(pi/15), length.out = 2)
start <- as.matrix(expand.grid(f,f,f,f))

# function to select the first coordinate positive~
start_pos<- matrix(NA, 1, 4)

for (i in 1:nrow(start)) {
  temp = polar2cart(start[i,], 1)
  if(temp[1]>0.02) {
    start_pos<- rbind(start_pos, temp)
  }
}

start_pos <- start_pos[-1,]
rownames(start_pos)<- c(1:nrow(start_pos))

# cost function for Nelder-Mead algorithm
# Inputs:  1) din  : densities of 'age-at-death' distribution on the support dSup
#          2) x   : covariate matrix
#          3) eta : the polar coordinates to be optimized.
#          4) dSup: support for the age-of-death distribution
#          5) h   : bandwidth for frechet regression
qin <- quant_all
x <- X_ctr
eta <- start_pos[1,]
h <- h[1]

nm_cost <- function(qin, x, eta, qSup, h) { 
  
  theta = polar2cart(eta,1)     # cartesian coordinates from polar coordinates
  x_c <- as.matrix(x) %*% theta   # creating the single index
  
  # creating the local frechet object by optimizing over the Wasserstein space
  Y_pred <- LocDenReg(xin=as.vector(x_c), qin=qin, optns= list(bwReg= h,
                        qSup = qSup, dSup=dSup, lower=20, upper=110))    
  
  # creating the matrix for storing the MSPE 
  metric_v <- matrix(0, dim(Y_pred$qin)[1], 1) 
  for(k in 1:nrow(metric_v)) {
  #get Wasserstein distance between the response and its estimate over each observation
    metric_v[k,1]<- dist4den(list(x= qSup, y= Y_pred$qout[k,]), list(x= qSup, y=Y_pred$qin[k,]), 
                             fctn_type = 'quantile') 
  }
  zx <- apply(metric_v, 2, mean) # computing MSPE over all observations
  return(zx)
}














#######################################################
# Choosing bandwidth by leave-one-out Cross-Validation
#######################################################

# modifying the cost function for leave-one-out CV

# cost function for Nelder-Mead algorithm
# Inputs:  1) i   : index for the leave-one-out.
#          2) qin : quantiles of 'age-at-death' distribution on the support qSup
#          2) x   : covariate matrix
#          3) eta : the polar coordinates to be optimized.
#          4) dSup: support for the age-of-death distribution.
#          5) h   : bandwidth for frechet regression.

h=0.8
nm_cost_L1OCV <- function(i, qin, x, eta, qSup, h) { 

  theta = polar2cart(eta,1)      # cartesian coordinates from polar coordinates
  x_in<-as.matrix(x[-i,])%*%theta   # creating the single index excluding i-th observation
  
  x_out <- x[i,]%*%theta  # creating the single index for i-th observation
  
  q_in <- qin[-i,] # the quantiles excluding i-th observation
  q_out <- qin[i,] # the quantiles for the i-th observation
  # creating the local frechet object by optimizing over the Wasserstein space
  Y_pred <- LocDenReg(xin=as.vector(x_in), xout = x_out, qin=q_in, 
                      optns= list(bwReg= h, qSup = qSup, dSup=dSup, 
                                  lower=20, upper=110))    
  
  #get Wasserstein distance between the response and its estimate for i-th observation
  pe<- dist4den(list(x= qSup, y= Y_pred$qout), list(x= qSup, y=q_out), fctn_type = 'quantile') 
  return(pe)
}


# to find the bandwidth range for analysis~
h_max = max(apply(X_ctr, 1, e_norm))

## to find the lowest possible value for bandwidth h
metric_v_temp <- matrix(NA, 40 ,1)
mv <-matrix(NA, 40 ,1)

for (j in 1:40) {
  for (i in 1:40) {
    if(i!=j)
      mv[i]<- e_norm(X_ctr[j,] - X_ctr[i,])
  }
  metric_v_temp[j] = min(mv[-j])
}

h_min <- min(metric_v_temp)

# the sequence of bandwidths to optimize over ~
h = exp(seq(log(h_min),log(h_max), length.out = 30))


# function for leave-one-out cross validation for i-th observation
L1O_CV <- function(i, quant_all, X_ctr, start_pos, qSup, h) {
  # matrix to store the MSPE for a given bandwidth and a given eta. 
  MSPE_m = matrix(NA, nrow = nrow(start_pos), length(h))
  
  for (j in 1:length(h)) {
    print(paste('bw',j,sep='_'))
    ftemp <- function(eta) nm_cost_L1OCV(i, quant_all, X_ctr, eta, qSup, h[j])
    
    for (k in 1:nrow(start_pos)) {
      opt <- optimset(TolFun= 1e-10) # decide on tolerance for the algorithm
      run <- fminbnd(ftemp, as.vector(c(start_pos[k,1],start_pos[k,2],
                start_pos[k,3],start_pos[k,4])), c(-pi,-pi,-pi,-pi), 
                c(pi,pi,pi,pi), opt)
      MSPE_m[k, j] <- as.vector(neldermead.get(this=run, key ='fopt'))
    }
  }
  return(MSPE_m)
}

# eta_opt <- eta_mincost[which.min(min_cr), ]
# theta_opt <- polar2cart(eta_opt, 1)

## Leave-one-out cv

# loading necessary libraries again~
library('frechet') 
library('ggplot2')
library('neldermead')
library('foreach')
library('doParallel')
library('fdadensity')
library('rlang')
library('vctrs')

# setting up parallel backend to use many processors
cores <- detectCores()
cl <- makeCluster(cores[1]-1) # not to overload the server
registerDoParallel(cl)


finalMatrix <- foreach(i=1:40, .combine = sum) %dopar% {
    
  library('frechet')
  library('rlang')
  library('neldermead')
  tempMatrix <- L1O_CV(i, quant_all, X_ctr, start_pos, qSup, h)
    
  #nm_cost_L1OCV(i, quant_all, X_ctr, eta, qSup, h[1])
    
  tempMatrix
}

# stop cluster
stopCluster(cl)




###### creating 30 folds for 5-fold cross-validation 
fold <- matrix(NA, 30, 5)

set.seed(1121)
for (i in 1:nrow(fold)) {
  fold[i,] <- sample(1:40, 5, replace = F)
}

q_fold <- density_all[-fold[1,],]






















# Unused codes warehouse
########## For 2008 predicted~

xout<- seq(min(x_c), max(x_c), length.out =50)
  
#t1<-proc.time()
Y_pred <- LocDenReg(xin=x_c, yin=myls[[10]], xout =xout, optns=list(dSup=dSup, bwReg= 2, qSup= qSup))
Y_pred
  
x <- Y_pred$din[[1]]$x
x
  
  #For the in-sample
  z = matrix(0, 28, 101)
  for (i in 1:nrow(z)) {
    z[i,]<- Y_pred$dout[i,]
  }
  
  #For the out-sample
  z=Y_pred$dout
  z
  
  d<- cbind(data.frame(xout), as.data.frame(z))
  d<- d[order(d[,1]),]
  d
  
  #y = x_c
  
  fig1 <- plot_ly(x=seq(20,110,length.out = 101), y= d$xout, z=as.matrix(d[,-c(1,2)]), type = "heatmap")
  fig1 <- fig1 %>% layout(
    xaxis = list(title="Age"),
    yaxis = list(title= "Covariate Score")
  )
  fig1
  
########## For 2008 actual~
  
x<- matrix(0, 20, 101)
for (i in 1:nrow(x)) {
    x[i,]<- Y_pred$din[[i]]$x
  }
  
  y = x_c
  z = matrix(0, 20, 101)
  for (i in 1:nrow(z)) {
    z[i,]<- Y_pred$din[[i]]$y
  }
  
  fig1 <- plot_ly(
    type = 'surface',
    contours = list(
      x = list(show = TRUE, start = 1.5, end = 2, size = 0.04, color = 'white'),
      z = list(show = TRUE, start = 0.5, end = 0.8, size = 0.05)),
    x = ~x,
    y = ~y,
    z = ~z)
  fig1 <- fig1 %>% layout(
    title = "Mortality distribution of 2008",
    scene = list(
      xaxis = list(nticks = 20, title="Age"),
      zaxis = list(nticks = 4),
      yaxis = list(title= "Covariate Score"),
      camera = list(eye = list(x = 0, y = -1, z = 0.6)),
      aspectratio = list(x = .9, y = .8, z = 0.2)))
  fig1

#############################################
## Code for finding h via cross-validation ##
#############################################
  
  # s<- matrix(sample(seq(1,28), replace = F), 7, 4) # sampled matrix 
  # mse_s<- matrix(0, nrow(s), 1)
  
  #for (i in 1:nrow(beta)) {
  #  for (j in 1:nrow(s)) {
  #    r<- as.numeric(s[j,]) #select the counties to be removed at first run
  #    myls_temp<- vector("list", length =21) #temporary list object to store the subset data for all dates 
  #    for (l in 1:length(myls_temp)) { 
  #      myls_temp[[l]]<- myls[[l]][-r] #Removing the data for countries from each year
  #    }
  #    # removing country data from all the covariates
  #    gdp_temp<- gdp[,-r] 
  #    unemp_temp<- unemp[,-r]
  #    
  #    # For the grid of parameters~
  #    for (v in 1:nrow(beta)) {
  #      #Get the score variable
  #      x_c<- beta[v,1]*matrix(as.numeric(c(as.matrix(unemp_temp))),nrow(unemp_temp),ncol(unemp_temp),byrow=F) + 
  #        beta[v,2]*matrix(as.numeric(c(as.matrix(gdp_temp))),nrow(gdp_temp),ncol(gdp_temp),byrow=F) 
  #      
  #      metric_v<- matrix(0, nrow(x_c), 1) 
  #      
  #      for (t in 1:nrow(x_c)) {
  #        xctemp<- x_c[t,]
  #        yls<- myls_temp[[t]]
  #        metric_v[t,]<-criterion_row(xctemp, yls, h[i])
  #      }
  #      MSE_vec[v,]<- mean(metric_v)      
  #    }
  #    
  #    beta_opt<- beta[which.min(MSE_vec), ] # get the optimum beta using the bigger rest of the data of countries
  #    # get the covariates for the removed countries, then from the optimal beta above, get the score variable, then get the response of
  # the removed variables, then compute the model predictions, then the Wasserstein metric for predicted from the true, 
  #    gdp_n<- gdp[,r] 
  #    unemp_n<- unemp[,r]
  #    myls_n<- vector("list",length= 25) #temporary list object to store the subset data for all dates 
  #    for (l in 1:length(myls_n)) { 
  #      myls_n[[l]]<- myls[[l]][r] #Removing the data for countries from each year
  #    }
  #    
  #    x_c<- beta_opt[1]*matrix(as.numeric(c(as.matrix(unemp_n))),nrow(unemp_n),ncol(unemp_n),byrow=F) + 
  #      beta[2]*matrix(as.numeric(c(as.matrix(gdp_n))),nrow(gdp_n),ncol(gdp_n),byrow=F) #Get the score variable
  #    
  #    metric_v<- matrix(0, nrow(x_c), 1) #
  
  #    for (t in 1:nrow(x_c)) {
  #      xctemp<- x_c[t,]
  #      yls<- myls_n[[t]]
  #      metric_v[t,]<-criterion_row(xctemp, yls, h[i])
  #    }
  #    mse_s[j,] <- mean(metric_v)
  #  }
  #  #mse_h[i,]<- mean(mse_s)
  
  
  #################################
  ## Computing the distributions ##
  #################################
  #  choosing a sequence of bandwidths
  #h <- seq(0.9, 1.5, by = 0.02)
  
  # choice of lag
  # lag = 0
  
  #  for(i in 1:length(h)) {
  #    for(j in 1:nrow(unemp)) {
  #x_c<- beta[i,1]*as.numeric(unemp[j,]) + beta[i,2]*as.numeric(gdp[j,])
  #      Y_pred<- LocDenReg(xin = x, yin= myls[[j+1]], optns = list(bwReg=2, 
  #                qSup = qSup,dSup=dSup))
  #      metric_v<- matrix(0, dim(Y_pred$qin)[1], 1)
  
  #      for(k in 1:nrow(metric_v)) {
  #        metric_v[k,]<- W_metric(Y_pred$qout[k,], Y_pred$qin[k,], Y_pred$qSup)
  #      }
  #    }
  #    MSE_vec[i,] <- mean(metric_v)
  #  }
  
  #  betaopt0 <- beta[which.min(MSE_vec),]
  #  betaopt0
  
