args = commandArgs()
#library(session)
#restore.session(file=".RSession_mars")
library(hal9001)
library(pdp)

n <- 0
interact_names <- c("NIR_Blue_NDVI_pH","Slope_Conductivity_NIR_Manganese","Slope_Conductivity_NIR_RedEdge","Slope_Conductivity_NIR_Red","Slope_Conductivity_NIR_Green","Slope_Conductivity_NIR_Blue","Slope_Conductivity_NIR_NDVI","Slope_NIR_RedEdge_Red","Slope_NIR_RedEdge_Green","Slope_NIR_RedEdge_Blue","Slope_NIR_RedEdge_NDVI","Slope_NIR_Red_Manganese","Slope_NIR_Red_Green","Slope_NIR_Red_Blue","Slope_NIR_Red_NDVI","Slope_NIR_Green_Blue","Slope_NIR_Green_NDVI","Slope_NIR_Blue_pH","Slope_NIR_Blue_Manganese","Slope_NIR_Blue_NDVI","Slope_NIR_NDVI_OrganicMatter","Conductivity_NIR_RedEdge_Red","Conductivity_NIR_RedEdge_Blue","Conductivity_NIR_Red_pH","Conductivity_NIR_Red_Manganese","Conductivity_NIR_Red_Blue","Conductivity_NIR_Blue_pH","Conductivity_NIR_Blue_Manganese","Conductivity_NIR_Blue_NDVI","NIR_RedEdge_pH_CEC","NIR_RedEdge_Red_Blue","NIR_RedEdge_Red_NDVI","NIR_RedEdge_Blue_pH","NIR_RedEdge_Blue_Manganese","NIR_RedEdge_Blue_NDVI","NIR_Red_pH_OrganicMatter","NIR_Red_Blue_NDVI")
	
interact_names <- interact_names[as.numeric(args[4]):as.numeric(args[5])]
infile <- read.csv(args[3], header=TRUE)
infile_old <- infile[,]
infile <- infile[sample(nrow(infile)), ]
keep_fraction <- 0.001  #0.0025
infile <- infile[1:(keep_fraction*nrow(infile)), ]


holdout_limits <- c(621880, 621920)  #longitudinal limits of test data
buffer_size <- 20    #buffer around strip to exclude from training data
holdout_limits2 <- c(holdout_limits[1]-buffer_size, holdout_limits[2]+buffer_size)

tr1 <- rownames(infile)[(infile[,1] > holdout_limits[1])]
all_temp <- data.frame(infile[tr1,])
test_rows <- rownames(all_temp)[(all_temp[,1] < holdout_limits[2])]
infile_test <- infile[test_rows, ]

tr2 <- rownames(infile)[(infile[,1] < holdout_limits2[1])]
tr3 <- rownames(infile)[(infile[,1] > holdout_limits2[2])]
train_rows <- c(tr2, tr3)
infile_train <- infile[train_rows, ]


for (interact_name in interact_names){
    interact_vector = c(interact_name)
    split_string <- strsplit(interact_vector, '_')
    n <- n+1
    print(n)
    print(split_string)

    features  <- split_string[[1]]
    feature_1 <- split_string[[1]][1]
    feature_2 <- split_string[[1]][2]

    keep_cols_reduced <- c(features, 'Yield')

    infile_reduced <- infile[, keep_cols_reduced]
    infile_train_reduced <- infile_train[, keep_cols_reduced]
    infile_test_reduced  <- infile_test[ , keep_cols_reduced]

    degree = length(features)

    X=data.matrix(infile_train_reduced[,c(features)])
    Y=data.matrix(infile_train_reduced[,'Yield'])

    print(nrow(X))

    HAL_model = fit_hal(X=X, Y=Y, yolo=FALSE, max_degree=2, n_folds=4)
    	#take out fit_type argument for default
        #fit_type="lassi"


    grid_resolution=30
    if (degree==4){grid_resolution=20}   #to save computation

    seq1_limits = c(min(infile_old[,features[1]]), max(infile_old[, features[1]]))
    seq2_limits = c(min(infile_old[,features[2]]), max(infile_old[, features[2]]))
    seq3_limits = c(min(infile_old[,features[3]]), max(infile_old[, features[3]]))
    seq4_limits = c(min(infile_old[,features[4]]), max(infile_old[, features[4]]))
    seq1 <- seq(seq1_limits[1], seq1_limits[2], by=(seq1_limits[2]-seq1_limits[1])/(grid_resolution-1) )
    seq2 <- seq(seq2_limits[1], seq2_limits[2], by=(seq2_limits[2]-seq2_limits[1])/(grid_resolution-1) )
    seq3 <- seq(seq3_limits[1], seq3_limits[2], by=(seq3_limits[2]-seq3_limits[1])/(grid_resolution-1) ) 
    seq4 <- seq(seq4_limits[1], seq4_limits[2], by=(seq4_limits[2]-seq4_limits[1])/(grid_resolution-1) )

    X2 <- data.frame(matrix(0, ncol=4, nrow=(length(seq1)*length(seq2)*length(seq3)*length(seq4))))

    for (x1 in 1:length(seq1)){
	    for (x2 in 1:length(seq2)){
		    for (x3 in 1:length(seq3)){
			    for (x4 in 1:length(seq4)){
	        my_position <- (x1-1)*length(seq2)*length(seq3)*length(seq4) + (x2-1)*length(seq3)*length(seq4) + (x3-1)*length(seq4) + x4
                X2[my_position,1] = seq1[x1]
	        X2[my_position,2] = seq2[x2]
		X2[my_position,3] = seq3[x3]
		X2[my_position,4] = seq4[x4]
		         }
		    }
	    }
    }

    predicted_yield <- predict(HAL_model, new_data=X2)

    p3 <- cbind(X2, predicted_yield)

    write.csv(p3, file=paste(interact_name,"_HAL", ".csv", sep=""))

    X3 <- infile_old[, features]
    predicted_yield2 <- predict(HAL_model, new_data=X3)
    p4 <- cbind(infile_old, predicted_yield2)
    write.csv(p4, file=paste(interact_name, "_pointshal.csv", sep=""))


    }
