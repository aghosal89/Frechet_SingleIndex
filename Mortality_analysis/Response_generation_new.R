############
# Script #2
############

## In this script we set up the response data

# In the following segment we create the response data as objects in L^2-Wasserstein space to 
# run our FSI model.

# The data is provided in mortality.org website in the form of age-at-death for a given year 
# and a given country for the age [0,110+], we read the lifetable corresponding to a specific 
# country for year 2013, remove '+' sign from age '110+' category, gather the 'age-at-death' 
# for age range [20,110], subset data for the given year, create a single vector
# where the age is repeated as many times as the number of deaths in the dataset, and finally
# create quantile an array of quatiles of the vector of a given length, say m.

# Set working directory

# define m, number of support points for quantiles/densities.
m<- 101

# the year for which we analyze mortality.
Yr <- 2013

# this function computes the euclidean norm of the vector x:
e_norm <- function(x) sqrt(sum(x^2))

# This function creates the data vectors for computation of quantiles/densities
# Inputs: 1) x: vector of ages.
#         2) f: number of deaths for age in x.
# Output: a vector with ages repeated as many times as the number of deaths.

data_vec<- function(x,f) {
  l<- 0                     
  # creating vector where the age is repeated the number of deaths that year in the increasing order
  for (i in 1:length(x)) {        
    l<- c(l, rep(x[i], f[i]))   
  }
  # remove the initial 0 as the first element
  l=l[-1]  
  return(as.numeric(l))
}

# In 'mortality.org', the age of individuals above 110 are recorded as '110+'. 
# This function removes the '+' sign for further computation.
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


# Mortality of Australia:
#Reading the mortality data
lt_australia <- read.csv("lt_australia.csv", header = T) 
# remove '+' sign from 110 age
lt_australia<- Remove_plus(lt_australia)
# subsetting data by age-group
lt_australia <- subset(lt_australia, lt_australia$Age %in% seq(20,110)) 
# subset the death data by year 2013
lt_australia <- lt_australia[,c("Year", "Age", "dx")] 
lt_australia<- subset(lt_australia, lt_australia[,1] %in% Yr)
# creating data vector to create the quantile functions
myls_australia = data_vec(lt_australia[,"Age"], lt_australia[,"dx"])

# for the countries below we repeat the same procedute we did for 'Australia'
# Mortality of Austria:
lt_austria <- read.csv("lt_austria.csv", header = T) 
lt_austria<- Remove_plus(lt_austria)
lt_austria <- subset(lt_austria, lt_austria$Age %in% seq(20,110)) 
lt_austria <- lt_austria[,c("Year", "Age", "dx")] 
lt_austria<- subset(lt_austria, lt_austria[,1] %in% Yr)
myls_austria <- data_vec(lt_austria[,"Age"], lt_austria[,"dx"])  

# Mortality of Belarus:
lt_belarus <- read.csv("lt_belarus.csv", header = T) 
lt_belarus<- Remove_plus(lt_belarus)
lt_belarus <- subset(lt_belarus, lt_belarus$Age %in% seq(20,110)) 
lt_belarus <- lt_belarus[,c("Year", "Age", "dx")] 
lt_belarus<- subset(lt_belarus, lt_belarus[,1] %in% Yr)
myls_belarus <- data_vec(lt_belarus[,"Age"], lt_belarus[,"dx"])  

# Mortality of Belgium:
lt_belgium <- read.csv("lt_belgium.csv", header = T) 
lt_belgium<- Remove_plus(lt_belgium)
lt_belgium <- subset(lt_belgium, lt_belgium$Age %in% seq(20,110)) 
lt_belgium <- lt_belgium[,c("Year", "Age", "dx")] 
lt_belgium<- subset(lt_belgium, lt_belgium[,1] %in% Yr)
myls_belgium <- data_vec(lt_belgium[,"Age"], lt_belgium[,"dx"])  

# Mortality of bulgaria:
lt_bulgaria <- read.csv("lt_bulgaria.csv", header = T) 
lt_bulgaria<- Remove_plus(lt_bulgaria)
lt_bulgaria <- subset(lt_bulgaria, lt_bulgaria$Age %in% seq(20,110)) 
lt_bulgaria <- lt_bulgaria[,c("Year", "Age", "dx")] 
lt_bulgaria<- subset(lt_bulgaria, lt_bulgaria[,1] %in% Yr)
myls_bulgaria <- data_vec(lt_bulgaria[,"Age"], lt_bulgaria[,"dx"])  

# Mortality of canada:
lt_canada <- read.csv("lt_canada.csv", header = T) 
lt_canada<- Remove_plus(lt_canada)
lt_canada <- subset(lt_canada, lt_canada$Age %in% seq(20,110)) 
lt_canada <- lt_canada[,c("Year", "Age", "dx")] 
lt_canada<- subset(lt_canada, lt_canada[,1] %in% Yr)
myls_canada <- data_vec(lt_canada[,"Age"], lt_canada[,"dx"])  

# Mortality of Chile:
lt_chile <- read.csv("lt_chile.csv", header = T) 
lt_chile<- Remove_plus(lt_chile)
lt_chile <- subset(lt_chile, lt_chile$Age %in% seq(20,110)) 
lt_chile <- lt_chile[,c("Year", "Age", "dx")] 
lt_chile<- subset(lt_chile, lt_chile[,1] %in% Yr)
myls_chile <- data_vec(lt_chile[,"Age"], lt_chile[,"dx"])  

# Mortality of Croatia:
lt_croatia <- read.csv("lt_croatia.csv", header = T) 
lt_croatia<- Remove_plus(lt_croatia)
lt_croatia <- subset(lt_croatia, lt_croatia$Age %in% seq(20,110)) 
lt_croatia <- lt_croatia[,c("Year", "Age", "dx")] 
lt_croatia<- subset(lt_croatia, lt_croatia[,1] %in% Yr)
myls_croatia <- data_vec(lt_croatia[,"Age"], lt_croatia[,"dx"])  

# Mortality of Czech Republic:
lt_czech <- read.csv("lt_czech.csv", header = T) 
lt_czech<- Remove_plus(lt_czech)
lt_czech <- subset(lt_czech, lt_czech$Age %in% seq(20,110)) 
lt_czech <- lt_czech[,c("Year", "Age", "dx")] 
lt_czech<- subset(lt_czech, lt_czech[,1] %in% Yr)
myls_czech <- data_vec(lt_czech[,"Age"], lt_czech[,"dx"])  

# Mortality of denmark:
lt_denmark <- read.csv("lt_denmark.csv", header = T) 
lt_denmark<- Remove_plus(lt_denmark)
lt_denmark <- subset(lt_denmark, lt_denmark$Age %in% seq(20,110)) 
lt_denmark <- lt_denmark[,c("Year", "Age", "dx")] 
lt_denmark<- subset(lt_denmark, lt_denmark[,1] %in% Yr)
myls_denmark <- data_vec(lt_denmark[,"Age"], lt_denmark[,"dx"])  

# Mortality of estonia:
lt_estonia <- read.csv("lt_estonia.csv", header = T) 
lt_estonia<- Remove_plus(lt_estonia)
lt_estonia <- subset(lt_estonia, lt_estonia$Age %in% seq(20,110)) 
lt_estonia <- lt_estonia[,c("Year", "Age", "dx")] 
lt_estonia<- subset(lt_estonia, lt_estonia[,1] %in% Yr)
myls_estonia <- data_vec(lt_estonia[,"Age"], lt_estonia[,"dx"])  

# Mortality of Finland:
lt_finland <- read.csv("lt_finland.csv", header = T) 
lt_finland<- Remove_plus(lt_finland)
lt_finland <- subset(lt_finland, lt_finland$Age %in% seq(20,110)) 
lt_finland <- lt_finland[,c("Year", "Age", "dx")] 
lt_finland<- subset(lt_finland, lt_finland[,1] %in% Yr)
myls_finland <- data_vec(lt_finland[,"Age"], lt_finland[,"dx"])  

# Mortality of france:
lt_france <- read.csv("lt_france.csv", header = T) 
lt_france<- Remove_plus(lt_france)
lt_france <- subset(lt_france, lt_france$Age %in% seq(20,110)) 
lt_france <- lt_france[,c("Year", "Age", "dx")] 
lt_france<- subset(lt_france, lt_france[,1] %in% Yr)
myls_france <- data_vec(lt_france[,"Age"], lt_france[,"dx"])  

# Mortality of denmark:
lt_germany <- read.csv("lt_germany.csv", header = T) 
lt_germany<- Remove_plus(lt_germany)
lt_germany <- subset(lt_germany, lt_germany$Age %in% seq(20,110)) 
lt_germany <- lt_germany[,c("Year", "Age", "dx")] 
lt_germany<- subset(lt_germany, lt_germany[,1] %in% Yr)
myls_germany <- data_vec(lt_germany[,"Age"], lt_germany[,"dx"])  


# Mortality of greece:
lt_greece <- read.csv("lt_greece.csv", header = T) 
lt_greece<- Remove_plus(lt_greece)
lt_greece <- subset(lt_greece, lt_greece$Age %in% seq(20,110)) 
lt_greece <- lt_greece[,c("Year", "Age", "dx")] 
lt_greece<- subset(lt_greece, lt_greece[,1] %in% Yr)
myls_greece <- data_vec(lt_greece[,"Age"], lt_greece[,"dx"])  

# Mortality of hungary:
lt_hungary<- read.csv("lt_hungary.csv", header = T) 
lt_hungary<- Remove_plus(lt_hungary)
lt_hungary <- subset(lt_hungary, lt_hungary$Age %in% seq(20,110)) 
lt_hungary <- lt_hungary[,c("Year", "Age", "dx")] 
lt_hungary<- subset(lt_hungary, lt_hungary[,1] %in% Yr)
myls_hungary <- data_vec(lt_hungary[,"Age"], lt_hungary[,"dx"])  

# Mortality of iceland:
lt_iceland<- read.csv("lt_iceland.csv", header = T) 
lt_iceland<- Remove_plus(lt_iceland)
lt_iceland <- subset(lt_iceland, lt_iceland$Age %in% seq(20,110)) 
lt_iceland <- lt_iceland[,c("Year", "Age", "dx")] 
lt_iceland<- subset(lt_iceland, lt_iceland[,1] %in% Yr)
myls_iceland <- data_vec(lt_iceland[,"Age"], lt_iceland[,"dx"])


# Mortality of ireland:
lt_ireland<- read.csv("lt_ireland.csv", header = T) 
lt_ireland<- Remove_plus(lt_ireland)
lt_ireland <- subset(lt_ireland, lt_ireland$Age %in% seq(20,110)) 
lt_ireland <- lt_ireland[,c("Year", "Age", "dx")] 
lt_ireland<- subset(lt_ireland, lt_ireland[,1] %in% Yr)
myls_ireland<- data_vec(lt_ireland[,"Age"], lt_ireland[,"dx"])  

# Mortality of israel:
lt_israel<- read.csv("lt_israel.csv", header = T) 
lt_israel<- Remove_plus(lt_israel)
lt_israel <- subset(lt_israel, lt_israel$Age %in% seq(20,110)) 
lt_israel <- lt_israel[,c("Year", "Age", "dx")] 
lt_israel<- subset(lt_israel, lt_israel[,1] %in% Yr)

myls_israel<- data_vec(lt_israel[,"Age"], lt_israel[,"dx"])  

# Mortality of italy:
lt_italy<- read.csv("lt_italy.csv", header = T) 
lt_italy<- Remove_plus(lt_italy)
lt_italy <- subset(lt_italy, lt_italy$Age %in% seq(20,110)) 
lt_italy <- lt_italy[,c("Year", "Age", "dx")] 
lt_italy<- subset(lt_italy, lt_italy[,1] %in% Yr)

myls_italy<- data_vec(lt_italy[,"Age"], lt_italy[,"dx"])  

# Mortality of japan:
lt_japan <- read.csv("lt_japan.csv", header = T) 
lt_japan<- Remove_plus(lt_japan)
lt_japan <- subset(lt_japan, lt_japan$Age %in% seq(20,110)) 
lt_japan <- lt_japan[,c("Year", "Age", "dx")] 
lt_japan<- subset(lt_japan, lt_japan[,1] %in% Yr)

myls_japan <- data_vec(lt_japan[,"Age"], lt_japan[,"dx"])  

# Mortality of latvia:
lt_latvia <- read.csv("lt_latvia.csv", header = T) 
lt_latvia<- Remove_plus(lt_latvia)
lt_latvia <- subset(lt_latvia, lt_latvia$Age %in% seq(20,110)) 
lt_latvia <- lt_latvia[,c("Year", "Age", "dx")] 
lt_latvia<- subset(lt_latvia, lt_latvia[,1] %in% Yr)

myls_latvia <- data_vec(lt_latvia[,"Age"], lt_latvia[,"dx"])  


# Mortality of lithuania:
lt_lithuania<- read.csv("lt_lithuania.csv", header = T) 
lt_lithuania<- Remove_plus(lt_lithuania)
lt_lithuania <- subset(lt_lithuania, lt_lithuania$Age %in% seq(20,110)) 
lt_lithuania <- lt_lithuania[,c("Year", "Age", "dx")] 
lt_lithuania<- subset(lt_lithuania, lt_lithuania[,1] %in% Yr)

myls_lithuania<- data_vec(lt_lithuania[,"Age"], lt_lithuania[,"dx"])  


# Mortality of luxembourg:
lt_luxembourg<- read.csv("lt_luxembourg.csv", header = T) 
lt_luxembourg<- Remove_plus(lt_luxembourg)
lt_luxembourg <- subset(lt_luxembourg, lt_luxembourg$Age %in% seq(20,110)) 
lt_luxembourg <- lt_luxembourg[,c("Year", "Age", "dx")] 
lt_luxembourg<- subset(lt_luxembourg, lt_luxembourg[,1] %in% Yr)

myls_luxembourg<- data_vec(lt_luxembourg[,"Age"], lt_luxembourg[,"dx"])  


# Mortality of netherlands:
lt_netherlands<- read.csv("lt_netherlands.csv", header = T) 
lt_netherlands<- Remove_plus(lt_netherlands)
lt_netherlands <- subset(lt_netherlands, lt_netherlands$Age %in% seq(20,110)) 
lt_netherlands <- lt_netherlands[,c("Year", "Age", "dx")] 
lt_netherlands<- subset(lt_netherlands, lt_netherlands[,1] %in% Yr)

myls_netherlands<- data_vec(lt_netherlands[,"Age"], lt_netherlands[,"dx"])  

# Mortality of New Zealand:
lt_newz<- read.csv("lt_newz.csv", header = T) 
lt_newz<- Remove_plus(lt_newz)
lt_newz <- subset(lt_newz, lt_newz$Age %in% seq(20,110)) 
lt_newz <- lt_newz[,c("Year", "Age", "dx")] 
lt_newz<- subset(lt_newz, lt_newz[,1] %in% Yr)
myls_newz<- data_vec(lt_newz[,"Age"], lt_newz[,"dx"])  

# Mortality of Norway:
lt_norway<- read.csv("lt_norway.csv", header = T) 
lt_norway<- Remove_plus(lt_norway)
lt_norway <- subset(lt_norway, lt_norway$Age %in% seq(20,110)) 
lt_norway <- lt_norway[,c("Year", "Age", "dx")] 
lt_norway<- subset(lt_norway, lt_norway[,1] %in% Yr)
myls_norway<- data_vec(lt_norway[,"Age"], lt_norway[,"dx"])  


# Mortality of Poland:
lt_poland<- read.csv("lt_poland.csv", header = T) 
lt_poland<- Remove_plus(lt_poland)
lt_poland <- subset(lt_poland, lt_poland$Age %in% seq(20,110)) 
lt_poland <- lt_poland[,c("Year", "Age", "dx")] 
lt_poland<- subset(lt_poland, lt_poland[,1] %in% Yr)

myls_poland<- data_vec(lt_poland[,"Age"], lt_poland[,"dx"])  

# Mortality of Portugal:
lt_portugal<- read.csv("lt_portugal.csv", header = T) 
lt_portugal<- Remove_plus(lt_portugal)
lt_portugal <- subset(lt_portugal, lt_portugal$Age %in% seq(20,110)) 
lt_portugal <- lt_portugal[,c("Year", "Age", "dx")] 
lt_portugal<- subset(lt_portugal, lt_portugal[,1] %in% Yr)

myls_portugal<- data_vec(lt_portugal[,"Age"], lt_portugal[,"dx"])  


# Mortality of korea:
lt_korea<- read.csv("lt_korea.csv", header = T) 
lt_korea<- Remove_plus(lt_korea)
lt_korea <- subset(lt_korea, lt_korea$Age %in% seq(20,110)) 
lt_korea <- lt_korea[,c("Year", "Age", "dx")] 
lt_korea<- subset(lt_korea, lt_korea[,1] %in% Yr)

myls_korea<- data_vec(lt_korea[,"Age"], lt_korea[,"dx"])  


# Mortality of russia:
lt_russia<- read.csv("lt_russia.csv", header = T) 
lt_russia<- Remove_plus(lt_russia)
lt_russia <- subset(lt_russia, lt_russia$Age %in% seq(20,110)) 
lt_russia <- lt_russia[,c("Year", "Age", "dx")] 
lt_russia<- subset(lt_russia, lt_russia[,1] %in% Yr)

myls_russia<- data_vec(lt_russia[,"Age"], lt_russia[,"dx"])  


# Mortality of slovakia:
lt_slovakia<- read.csv("lt_slovakia.csv", header = T) 
lt_slovakia<- Remove_plus(lt_slovakia)
lt_slovakia <- subset(lt_slovakia, lt_slovakia$Age %in% seq(20,110)) 
lt_slovakia <- lt_slovakia[,c("Year", "Age", "dx")] 
lt_slovakia<- subset(lt_slovakia, lt_slovakia[,1] %in% Yr)

myls_slovakia<- data_vec(lt_slovakia[,"Age"], lt_slovakia[,"dx"])  


# Mortality of slovenia:
lt_slovenia<- read.csv("lt_slovenia.csv", header = T) 
lt_slovenia<- Remove_plus(lt_slovenia)
lt_slovenia <- subset(lt_slovenia, lt_slovenia$Age %in% seq(20,110)) 
lt_slovenia <- lt_slovenia[,c("Year", "Age", "dx")] 
lt_slovenia<- subset(lt_slovenia, lt_slovenia[,1] %in% Yr)

myls_slovenia<- data_vec(lt_slovenia[,"Age"], lt_slovenia[,"dx"])  


# Mortality of spain:
lt_spain<- read.csv("lt_spain.csv", header = T) 
lt_spain<- Remove_plus(lt_spain)
lt_spain <- subset(lt_spain, lt_spain$Age %in% seq(20,110)) 
lt_spain <- lt_spain[,c("Year", "Age", "dx")] 
lt_spain<- subset(lt_spain, lt_spain[,1] %in% Yr)

myls_spain<- data_vec(lt_spain[,"Age"], lt_spain[,"dx"])  


# Mortality of sweden:
lt_sweden<- read.csv("lt_sweden.csv", header = T) 
lt_sweden<- Remove_plus(lt_sweden)
lt_sweden <- subset(lt_sweden, lt_sweden$Age %in% seq(20,110)) 
lt_sweden <- lt_sweden[,c("Year", "Age", "dx")] 
lt_sweden<- subset(lt_sweden, lt_sweden[,1] %in% Yr)

myls_sweden<- data_vec(lt_sweden[,"Age"], lt_sweden[,"dx"])  


# Mortality of switzerland:
lt_switzerland<- read.csv("lt_switzerland.csv", header = T) 
lt_switzerland<- Remove_plus(lt_switzerland)
lt_switzerland <- subset(lt_switzerland, lt_switzerland$Age %in% seq(20,110)) 
lt_switzerland <- lt_switzerland[,c("Year", "Age", "dx")] 
lt_switzerland<- subset(lt_switzerland, lt_switzerland[,1] %in% Yr)

myls_switzerland<- data_vec(lt_switzerland[,"Age"], lt_switzerland[,"dx"])  


# Mortality of UK:
lt_UK<- read.csv("lt_UK.csv", header = T) 
lt_UK<- Remove_plus(lt_UK)
lt_UK <- subset(lt_UK, lt_UK$Age %in% seq(20,110)) 
lt_UK <- lt_UK[,c("Year", "Age", "dx")] 
lt_UK<- subset(lt_UK, lt_UK[,1] %in% Yr)
myls_UK<- data_vec(lt_UK[,"Age"], lt_UK[,"dx"])  


# Mortality of USA:
lt_USA<- read.csv("lt_USA.csv", header = T) 
lt_USA<- Remove_plus(lt_USA)
lt_USA <- subset(lt_USA, lt_USA$Age %in% seq(20,110)) 
lt_USA <- lt_USA[,c("Year", "Age", "dx")] 
lt_USA<- subset(lt_USA, lt_USA[,1] %in% Yr)

myls_USA<- data_vec(lt_USA[,"Age"], lt_USA[,"dx"])  


# Mortality of ukraine:
lt_ukraine<- read.csv("lt_Ukraine.csv", header = T) 
lt_ukraine<- Remove_plus(lt_ukraine)
lt_ukraine <- subset(lt_ukraine, lt_ukraine$Age %in% seq(20,110)) 
lt_ukraine <- lt_ukraine[,c("Year", "Age", "dx")] 
lt_ukraine<- subset(lt_ukraine, lt_ukraine[,1] %in% Yr)

myls_ukraine<- data_vec(lt_ukraine[,"Age"], lt_ukraine[,"dx"])  

# the grid of quantiles 
qSup <- seq(0,1, length.out=m)   

# the grid on the densities
dSup <- seq(20,110, length.out=m)  


#################################
# Main Algorithm for computation
#################################

# creating the list of age-at-death vectors for countries  

myls <- list()
myls$australia <-myls_australia; myls$austria <- myls_austria; myls$belgium <- myls_belgium
myls$bulgaria <- myls_bulgaria; myls$belarus<- myls_belarus; myls$canada <- myls_canada
myls$switzerland <- myls_switzerland; myls$chile <- myls_chile; myls$czech <- myls_czech
myls$germany <- myls_germany; myls$denmark <- myls_denmark; myls$spain <- myls_spain
myls$estonia <- myls_estonia; myls$finland <- myls_finland; myls$france <- myls_france
myls$uk <- myls_UK; myls$greece <- myls_greece;
myls$croatia <- myls_croatia; myls$hungary <- myls_hungary; myls$ireland <- myls_ireland
myls$iceland <- myls_iceland; myls$israel <- myls_israel; myls$italy <- myls_italy
myls$japan <- myls_japan; myls$korea <- myls_korea; myls$lithuania <- myls_lithuania
myls$luxembourg<-myls_luxembourg; myls$latvia<-myls_latvia; myls$netherlands<-myls_netherlands
myls$norway <- myls_norway; myls$newz <- myls_newz; myls$poland <- myls_poland;
myls$portugal <- myls_portugal; myls$russia <- myls_russia; myls$slovakia <- myls_slovakia;
myls$slovenia <- myls_slovenia; myls$sweden <- myls_sweden; myls$ukraine <- myls_ukraine; 
myls$usa <- myls_USA;


# define n as the number of observations (countries)
n<- length(myls)

# read necessary libraries for estimating densities and quantiles 
library('frechet') 
library('fdadensity')

# list of countries were saved separately in a .csv file in the exact order appearing above
# in the script 'covariate_generation.R'

# here we read the list of those countries 
country_ <- read.csv("Countries_FSI.csv", header = T)
country_ <- country_[,2]

# creating mortality densities for each country 
density_all <- matrix(NA, n, m)
for (i in 1:n) {
  density_all[i,] <- CreateDensity(y= myls[[i]], optns= list(nRegGrid = m))$y
}

density_all <- data.frame(density_all)
rownames(density_all) <- country_

# save densities for further use 
write.csv(density_all, "density_all.csv")

# creating the quantiles for the mortality data 
quant_all <- matrix(NA, nrow = n, ncol = length(qSup))
for (i in 1:n) {
  quant_all[i,]<-dens2quantile(dens = as.numeric(density_all[i,]), 
                                 dSup=dSup, qSup = qSup, useSplines = F)
}

quant_all <- data.frame(quant_all)
rownames(quant_all) <- country_

# save the quantiles for further use 
write.csv(quant_all, "quant_all.csv")

