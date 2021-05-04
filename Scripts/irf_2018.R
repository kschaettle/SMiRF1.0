args = commandArgs()
'%ni%' <- Negate('%in%')

library(magrittr)
#library(httr)
library(dplyr)
library(data.table)
#httr::set_config(config(ssl_verifypeer = 0L))  # only need this if you're getting ssl errors from R     #for consistency

all_data <- read.csv(args[3], header=TRUE)
all_data <- all_data[sample(nrow(all_data)),]  #shuffle
all_data <- all_data[1:(nrow(all_data)*as.numeric(args[5])),]   #sub fraction of observations  

ncols <- length(all_data[1, ])

yield_col_id <- c('Yield')    #this is the target column for regression
library(ggplot2)
library(iRF)
#source('surfacePlot.R')        #from iRF package; must be in same directory

all_data_backup <- all_data

holdout_limits <- c(621880, 621920)  #longitudinal limits of test data
buffer_size <- 20    #buffer around strip to exclude from training data
holdout_limits2 <- c(holdout_limits[1]-buffer_size, holdout_limits[2]+buffer_size)

infile <- all_data
tr1 <- rownames(infile)[(infile[,1] > holdout_limits[1])]
all_temp <- data.frame(infile[tr1,])
test_rows <- rownames(all_temp)[(all_temp[,1] < holdout_limits[2])]

tr2 <- rownames(infile)[(infile[,1] < holdout_limits2[1])]
tr3 <- rownames(infile)[(infile[,1] > holdout_limits2[2])]
train_rows <- c(tr2, tr3)

#set up training matrix and testing matrix
train_mat  <- all_data[train_rows,]
test_mat   <- all_data[test_rows,]
all_data   <- rbind(train_mat, test_mat)
colnames(all_data) <- colnames(all_data_backup)
all_data[1,]

keep_columns <- c('Conductivity', 'Phosphorus','Potassium','Zinc','Sulfur','Boron','Magnesium','Manganese','Copper','OrganicMatter')
keep_columns2 <- c('Conductivity2018', 'Phosphorus2018','Potassium2018','Zinc2018','Sulfur2018','Boron2018','Magnesium2018','Manganese2018','Copper2018','OrganicMatter2018') 		   
X_train <- all_data[train_rows, keep_columns] 
X_train <- data.matrix(X_train)

X_test  <- all_data[test_rows, keep_columns]
X_test  <- data.matrix(X_test)
Y_train <- data.matrix(all_data[train_rows, yield_col_id])
Y_test  <- data.matrix(all_data[test_rows , yield_col_id])

num_train <- length(X_train[,1])      #number of training points
length(Y_test)

X2 <- data.matrix(all_data[, keep_columns2])

x <- data.matrix(rbind(X_train, X_test))
y <- data.matrix(rbind(Y_train, Y_test))
y <- as.numeric(y)

varnames <- colnames(x)

n <- nrow(x)
p <- ncol(x)
train.id <- 1:num_train   #since we re-organized all_data
total_iterations <- 10

t1 <- Sys.time()
f <- iRF(x=x[train.id,], y=y[train.id], n.iter=total_iterations, n.bootstrap=1, n.core=4, get.prevalence=FALSE, int.sign=TRUE)
options(tibble.print_max=Inf)

colnames(X2) <- varnames
rf <- f$rf.list[[total_iterations]]  #new random forest model
prediction <- predict(rf, X2)
prediction
outmat = cbind(all_data, prediction)
outmat

write.csv(data.frame(outmat), file='zzziRF_2018_residuals.csv')
