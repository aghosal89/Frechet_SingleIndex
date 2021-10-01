
# function for euclidean norm~
e_norm <- function(x) sqrt(sum(x^2))

# The function to create the data vectors for computation of quantiles/density
data_vec<- function(x,f) {
  l<- 0    # initiate the vector by 0                  
  for (i in 1:length(x)) {        
    l<- c(l, rep(x[i], f[i]))   # repeat the age for number of deaths that year and create vector in the increasing order
  }
  l=l[-1] # remove the 0 as the first element 
  return(as.numeric(l))
}

#Remove the + sign at age 110
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


# function to compute mode age of each density
mode_ages <- function(dSup, densities) {
  r<- nrow(densities)
  modes <- matrix(0, r,1)
  for (j in 1:r) {
    modes[j] <- dSup[which.max(densities[j,])]
  }
  return(modes)
}



