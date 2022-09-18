
############################
## The following codes can 
## be used to generate the 
## folds for analysis:
############################

# Output:

# 30x10 matrix whose rows are indices for observations to be included 
# in the testing partition of the sample. The rest of observations
# are to consttitute the training partition.

Folds_new<- matrix(0, 30, 10)
set.seed(1989)

for(i in 1:nrow(Folds_new)) {
  Folds_new[i,]<- sample(39,10, replace = FALSE)
}

# set working directory

# save the data on folds
write.csv(Folds_new, "Folds_new.csv")

