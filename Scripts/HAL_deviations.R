args = commandArgs()
#library(session)
#restore.session(file=".RSession_mars")
library(hal9001)
library(pdp)


n <- 0

interact_names <- c("Conductivity_pH","Slope_pH","pH_OrganicMatter","pH_Manganese","pH_CEC","Phosphorus_pH","pH_Boron","pH_Magnesium","pH_Copper","pH_Zinc", "Phosphorus_Potassium")

interact_names <- interact_names[as.numeric(args[4]):as.numeric(args[5])]
infile <- read.csv(args[3], header=TRUE)
infile <- infile[sample(nrow(infile)), ]

keep_fraction <- 0.05    #fraction of data points to use for making surfaces
copies <- 100    #number of times to re-sample for making the surfaces


p3 <- c(0)
for (interact_name in interact_names){
	n <- n+1
	interact_vector = c(interact_name)
	split_string <- strsplit(interact_vector, '_')
	print(n)
	print(split_string[[1]])

    for (copy_iter in 1:copies){
	    print(copy_iter)
	    infile_new <- infile[sample(nrow(infile)), ]
	    infile_new <- infile_new[1:(keep_fraction*nrow(infile)), ]

	    holdout_limits <- c(621880, 621920)  #longitudinal limits of test data
	    buffer_size <- 20    #buffer around strip to exclude from training data
	    holdout_limits2 <- c(holdout_limits[1]-buffer_size, holdout_limits[2]+buffer_size)

	    tr1 <- rownames(infile_new)[(infile_new[,1] > holdout_limits[1])]
	    all_temp <- data.frame(infile_new[tr1,])
	    test_rows <- rownames(all_temp)[(all_temp[,1] < holdout_limits[2])]
	    infile_test <- infile_new[test_rows, ]

	    tr2 <- rownames(infile_new)[(infile_new[,1] < holdout_limits2[1])]
	    tr3 <- rownames(infile_new)[(infile_new[,1] > holdout_limits2[2])]
	    train_rows <- c(tr2, tr3)
	    infile_train <- infile_new[train_rows, ]

    features  <- split_string[[1]]

    keep_cols_reduced <- c(features, 'Yield')

    infile_reduced <- infile_new[, keep_cols_reduced]
    infile_train_reduced <- infile_train[, keep_cols_reduced]
    infile_test_reduced  <- infile_test[ , keep_cols_reduced]

    degree = length(features)

    X=data.matrix(infile_train_reduced[,c(features)])
    Y=data.matrix(infile_train_reduced[,'Yield'])

    HAL_model = fit_hal(X=X, Y=Y)

    grid_resolution=30
    if (degree==4){grid_resolution=20}   #to save computation

    seq1_limits = c(min(infile[,features[1]]), max(infile[, features[1]]))
    seq2_limits = c(min(infile[,features[2]]), max(infile[, features[2]]))
    seq1 <- seq(seq1_limits[1], seq1_limits[2], by=(seq1_limits[2]-seq1_limits[1])/(grid_resolution-1) )
    seq2 <- seq(seq2_limits[1], seq2_limits[2], by=(seq2_limits[2]-seq2_limits[1])/(grid_resolution-1) )

    X2 <- data.frame(matrix(0, ncol=2, nrow=(length(seq1)*length(seq2))))

    for (x1 in 1:length(seq1)){
	    for (x2 in 1:length(seq2)){
	        my_position <- (x1-1)*length(seq2) + x2
                X2[my_position,1] = seq1[x1]
	        X2[my_position,2] = seq2[x2]
		         }
		    }

    predicted_yield <- predict(HAL_model, new_data=X2)

    if (copy_iter ==1){
    p3 <- X2
    }
    p3 <- cbind(p3, predicted_yield)      #add the column

    } #end copy iter

    #now we take the point-wise standard deviation
    p4 <- p3[, 1:(degree+1)]  #create copy of right size
    
    for (row_iter in 1:nrow(p4)){   #now we assign the pointwise standard deviations
           p4[row_iter, degree+1] <- sd(p3[row_iter, (degree+1):ncol(p3)])
	   }

    write.csv(p4, file=paste(interact_name,"_HAL_deviations", ".csv", sep=""))
 
    } #end for interact names
