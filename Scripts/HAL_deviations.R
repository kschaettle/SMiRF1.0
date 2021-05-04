args = commandArgs()
#library(session)
#restore.session(file=".RSession_mars")
library(hal9001)
library(pdp)


n <- 0
#c("Slope_Conductivity_pH_OrganicMatter","Slope_Conductivity_pH_Manganese","Conductivity_pH_Manganese_OrganicMatter","Slope_pH_Manganese_OrganicMatter","Slope_Conductivity_pH_CEC","Slope_Conductivity_Phosphorus_pH","Slope_Conductivity_pH_Boron","Conductivity_pH_CEC_OrganicMatter","Slope_Conductivity_pH_Magnesium","Slope_pH_CEC_OrganicMatter","Slope_Conductivity_pH_Copper","Conductivity_Phosphorus_pH_OrganicMatter","Conductivity_pH_Manganese_CEC","Slope_Phosphorus_pH_OrganicMatter","Slope_pH_Manganese_CEC","Conductivity_pH_Boron_OrganicMatter","Conductivity_Phosphorus_pH_Manganese","Conductivity_pH_Boron_Manganese","Slope_Phosphorus_pH_Manganese","Slope_pH_Boron_OrganicMatter","Slope_pH_Boron_Manganese","Conductivity_pH_Magnesium_OrganicMatter","Slope_pH_Magnesium_OrganicMatter","Slope_pH_Copper_OrganicMatter","Conductivity_pH_Magnesium_Manganese","Slope_pH_Magnesium_Manganese","pH_Manganese_CEC_OrganicMatter","Conductivity_pH_Manganese_Copper","Slope_pH_Manganese_Copper","Phosphorus_pH_Manganese_OrganicMatter","Conductivity_pH_Boron_CEC","Conductivity_Phosphorus_pH_CEC","Slope_Phosphorus_pH_CEC","Slope_pH_Boron_CEC","Slope_Conductivity_pH_Zinc","pH_Boron_Manganese_OrganicMatter","Conductivity_Phosphorus_pH_Boron","Slope_Phosphorus_pH_Boron","pH_Magnesium_Manganese_OrganicMatter","Conductivity_pH_Magnesium_CEC","Conductivity_Phosphorus_pH_Magnesium","Slope_pH_Copper_CEC","Conductivity_pH_Boron_Magnesium","Phosphorus_pH_CEC_OrganicMatter","Slope_Phosphorus_pH_Copper","Slope_pH_Zinc_OrganicMatter","Slope_pH_Zinc_Manganese","pH_Copper_CEC_OrganicMatter","Slope_Conductivity_Potassium_pH","Conductivity_Phosphorus_Potassium_pH")
	
#c("Slope_Conductivity_pH","Conductivity_pH_OrganicMatter","Slope_pH_OrganicMatter","Conductivity_pH_Manganese","Slope_pH_Manganese","pH_Manganese_OrganicMatter","Conductivity_pH_CEC","Slope_pH_CEC","Conductivity_Phosphorus_pH","Conductivity_pH_Boron","Slope_Phosphorus_pH","Slope_pH_Boron","Conductivity_pH_Magnesium","pH_CEC_OrganicMatter","Slope_pH_Magnesium","Conductivity_pH_Copper","Slope_pH_Copper","Phosphorus_pH_OrganicMatter","pH_Manganese_CEC","Phosphorus_pH_Manganese","pH_Boron_OrganicMatter","pH_Boron_Manganese","pH_Magnesium_OrganicMatter","pH_Magnesium_Manganese","pH_Boron_CEC","Conductivity_pH_Zinc","Slope_pH_Zinc","Phosphorus_pH_Boron","pH_Zinc_OrganicMatter","pH_Zinc_Manganese") 
	
interact_names <- c("Conductivity_pH","Slope_pH","pH_OrganicMatter","pH_Manganese","pH_CEC","Phosphorus_pH","pH_Boron","pH_Magnesium","pH_Copper","pH_Zinc", "Phosphorus_Potassium")

#c("Conductivity_Lime","Conductivity_Lime_Boron","Conductivity_Lime_CEC","Conductivity_Lime_Copper","Conductivity_Lime_Manganese","Conductivity_Lime_OrganicMatter","Conductivity_Lime_pH","Conductivity_Lime_Phosphorus","Conductivity_Lime_Potash","Conductivity_Lime_Potash_Copper","Conductivity_Lime_Potash_Manganese","Conductivity_Lime_Potash_OrganicMatter","Conductivity_Lime_Potash_pH","Conductivity_Lime_Potash_TSP","Conductivity_Lime_TSP","Conductivity_Lime_TSP_Boron","Conductivity_Lime_TSP_Manganese","Conductivity_Lime_TSP_OrganicMatter","Conductivity_Lime_TSP_pH","Conductivity_Lime_TSP_Phosphorus","Lime_Boron","Lime_Copper","Lime_Manganese","Lime_OrganicMatter","Lime_pH","Lime_Phosphorus","Lime_Potash","Lime_Potash_CEC","Lime_Potash_OrganicMatter","Lime_Potash_pH","Lime_Potash_TSP","Lime_Potash_TSP_Manganese","Lime_Potash_TSP_OrganicMatter","Lime_Potash_TSP_pH","Lime_TSP","Lime_TSP_Copper","Lime_TSP_Manganese","Lime_TSP_OrganicMatter","Lime_TSP_pH","Lime_TSP_Phosphorus","Slope_Conductivity_Lime","Slope_Conductivity_Lime_Boron","Slope_Conductivity_Lime_CEC","Slope_Conductivity_Lime_Copper","Slope_Conductivity_Lime_Manganese","Slope_Conductivity_Lime_OrganicMatter","Slope_Conductivity_Lime_pH","Slope_Conductivity_Lime_Phosphorus","Slope_Conductivity_Lime_Potash","Slope_Conductivity_Lime_TSP","Slope_Lime","Slope_Lime_Boron","Slope_Lime_Boron_OrganicMatter","Slope_Lime_CEC","Slope_Lime_Copper","Slope_Lime_Manganese","Slope_Lime_OrganicMatter","Slope_Lime_pH","Slope_Lime_Phosphorus","Slope_Lime_Potash","Slope_Lime_Potash_Boron","Slope_Lime_Potash_Manganese","Slope_Lime_Potash_OrganicMatter","Slope_Lime_Potash_TSP","Slope_Lime_TSP","Slope_Lime_TSP_CEC","Slope_Lime_TSP_Manganese","Slope_Lime_TSP_OrganicMatter","Slope_Lime_TSP_pH","Slope_Lime_TSP_Phosphorus")
interact_names <- interact_names[as.numeric(args[4]):as.numeric(args[5])]
infile <- read.csv(args[3], header=TRUE)
infile <- infile[sample(nrow(infile)), ]

keep_fraction <- 0.05
copies <- 25


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
