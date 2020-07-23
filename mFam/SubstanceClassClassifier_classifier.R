
######################################################################################
## classifiers

#calc_auprc <- function(model, data){
#  index_class2 <- data$Class == "m"
#  index_class1 <- data$Class == "p"
#  predictions <- predict(model, data, type = "prob")
#  pr.curve(predictions$Class2[index_class2],
#           predictions$Class2[index_class1],
#           curve = TRUE)
#}
#if(FALSE){
#auprcSummary <- function(data, lev = NULL, model = NULL){
#  index_class1 <- data$obs == "p"
#  index_class2 <- data$obs == "m"
#  
#  the_curve <- pr.curve(scores.class0 = data$p[index_class2], scores.class1 = data$p[index_class1], curve = FALSE)
#  
#  out <- the_curve$auc.integral
#  names(out) <- "AUPRC"
#  
#  out
#}
#twoClassSummary2 <- function(data, lev = NULL, model = NULL) 
#{
#  lvls <- levels(data$obs)
#  if (length(lvls) > 2) 
#    stop(paste("Your outcome has", length(lvls), "levels. The twoClassSummary() function isn't appropriate."))
#  requireNamespaceQuietStop("ModelMetrics")
#  if (!all(levels(data[, "pred"]) == lvls)) 
#    stop("levels of observed and predicted data do not match")
#  rocAUC <- ModelMetrics::auc(ifelse(data$obs == lev[2], 0, 
#                                     1), data[, lvls[1]])
#  out <- c(rocAUC, sensitivity(data[, "pred"], data[, "obs"], 
#                               lev[1]), specificity(data[, "pred"], data[, "obs"], lev[2]))
#  names(out) <- c("ROC", "Sens", "Spec")
#  out
#}
#}

caret_classifier <- list(
  train = function(matrix_train, classes_pm_train, modelName, classWeights){
    if(FALSE){
      ######################################################
      ## add class to matrix and use (inefficient->overhead) formula + data parameters of train
      
      ## add class factor column
      matrix_train_caret <- matrix_train
      rownames(matrix_train_caret) <- seq_len(nrow(matrix_train_caret))
      matrix_train_caret <- cbind(as.data.frame(as.matrix(matrix_train_caret)), rep(x = 0, times = nrow(matrix_train_caret)))
      colnames(matrix_train_caret)[[ncol(matrix_train_caret)]] <- "class"
      posRows <- classes_pm_train=="+"
      matrix_train_caret[ posRows, "class"] <- "p"
      matrix_train_caret[!posRows, "class"] <- "m"
      matrix_train_caret[, "class"] <- as.factor(matrix_train_caret[, "class"])
      
      ## rename mass column names
      colnamesHere <- colnames(matrix_train_caret)
      neutralLosses <- grepl(x = colnamesHere, pattern = "^\\-\\d+(\\.\\d+)?$")
      fragments     <- grepl(x = colnamesHere, pattern = "^\\d+(\\.\\d+)?$")
      colnamesHere[neutralLosses] <- gsub(x = colnamesHere[neutralLosses], pattern = "^-", replacement = "m")
      colnamesHere[fragments]     <- paste("p", colnamesHere[fragments], sep = "")
      colnames(matrix_train_caret) <- colnamesHere
    }
    
    ######################################################
    ## treat data and outcome independently and use x + y parameters of train for the sake of efficiency
    
    ## to matrix
    matrix_train_caret2 <- as.matrix(matrix_train)
    rownames(matrix_train_caret2) <- seq_len(nrow(matrix_train_caret2))
    classes_pm_train_caret2 <- classes_pm_train
    ## add class column
    posRows <- classes_pm_train_caret2=="+"
    classes_pm_train_caret2[ posRows] <- "p"
    classes_pm_train_caret2[!posRows] <- "m"
    class <- as.factor(classes_pm_train_caret2)
    ## make valid column names
    colnamesHere <- colnames(matrix_train_caret2)
    neutralLosses <- grepl(x = colnamesHere, pattern = "^\\-\\d+(\\.\\d+)?$")
    fragments     <- grepl(x = colnamesHere, pattern = "^\\d+(\\.\\d+)?$")
    colnamesHere[neutralLosses] <- gsub(x = colnamesHere[neutralLosses], pattern = "^-", replacement = "m")
    colnamesHere[fragments]     <- paste("p", colnamesHere[fragments], sep = "")
    colnames(matrix_train_caret2) <- colnamesHere
    
    ################################################################
    ## train control and class weights
    useRoc <- FALSE
    
    ctrl <- caret::trainControl(
      #method = "none",
      ## cross validation
      method = "repeatedcv", 
      number = 5,
      repeats = 1, 
      ## tuning
      search = "random",
      ## misc
      classProbs = TRUE, 
      summaryFunction = ifelse(test = useRoc, yes = caret::twoClassSummary, no = caret::prSummary), 
      #summaryFunction = auprcSummary, 
      #summaryFunction = caret::prSummary, 
      #summaryFunction = caret::twoClassSummary, 
      verboseIter = FALSE,
      allowParallel = TRUE
    )
    
    #modelProps <- caret::getModelInfo()[[modelName]]
    
    if(classWeights){
      #weight_p <- (1 / table(class)["p"]) * 0.5
      #weight_m <- (1 / table(class)["m"]) * 0.5
      weight_p <- table(class)["m"]
      weight_m <- table(class)["p"]
      class_weights <- ifelse(test = class == "p", yes = weight_p, no = weight_m)
    } else {
      class_weights <- NULL
    }
    
    ################################################################
    ## do
    #set.seed(1)
    #startTime <<- Sys.time()
    log <- capture.output({
    classifier <- caret::train(
      #class ~ .,
      #data = matrix_train_caret,
      x = matrix_train_caret2,
      y = class,
      weights = class_weights,
      method = modelName,
      tuneLength = 10,
      trControl = ctrl,
      metric = ifelse(test = useRoc, yes = "ROC", no = "AUC"), 
      #preProc = c("center", "scale")
      #tuneGrid = expand.grid(k = c(5, 11, 21, 25))
      preProc = NULL
    )
    })
    #endTime <- Sys.time()
    #rocs <- classifier$results$ROC
    
    ## Variable Importance
    #varImp(object=classifier)
    #importanceDf <- varImp(object=classifier)[["importance"]]
    
    return(classifier)
  },
  classify = function(classifier, matrix_test){
    matrix_test_caret <- matrix_test
    matrix_test_caret <- as.data.frame(as.matrix(matrix_test_caret))
    
    ## rename mass column names
    colnamesHere <- colnames(matrix_test_caret)
    neutralLosses <- grepl(x = colnamesHere, pattern = "^\\-\\d+(\\.\\d+)?$")
    fragments     <- grepl(x = colnamesHere, pattern = "^\\d+(\\.\\d+)?$")
    colnamesHere[neutralLosses] <- gsub(x = colnamesHere[neutralLosses], pattern = "^-", replacement = "m")
    colnamesHere[fragments]     <- paste("p", colnamesHere[fragments], sep = "")
    colnames(matrix_test_caret) <- colnamesHere
    
    ## https://www.analyticsvidhya.com/blog/2016/12/practical-guide-to-implement-machine-learning-with-caret-package-in-r-with-practice-problem/
    #predictions<-predict.train(object=model_gbm,testSet[,predictors],type="raw")
    
    #plsClasses <- predict(classifier, newdata = matrix_test)
    log <- capture.output({
      plsProbs <- predict(classifier, newdata = matrix_test_caret, type = "prob")
    })
    scores <- plsProbs$p #pred$predict[, "+", 1]
    
    return(scores)
  }
)
caret_classifier_stacked <- list(
  train = function(matrix_train, classes_pm_train, modelNames){
    ## TODO class weights
    
    ######################################################
    ## treat data and outcome independently and use x + y parameters of train for the sake of efficiency
    matrix_train_caret2 <- as.matrix(matrix_train)
    rownames(matrix_train_caret2) <- seq_len(nrow(matrix_train_caret2))
    classes_pm_train_caret2 <- classes_pm_train
    posRows <- classes_pm_train_caret2=="+"
    classes_pm_train_caret2[ posRows] <- "p"
    classes_pm_train_caret2[!posRows] <- "m"
    class <- as.factor(classes_pm_train_caret2)
    
    colnamesHere <- colnames(matrix_train_caret2)
    neutralLosses <- grepl(x = colnamesHere, pattern = "^\\-\\d+(\\.\\d+)?$")
    fragments     <- grepl(x = colnamesHere, pattern = "^\\d+(\\.\\d+)?$")
    colnamesHere[neutralLosses] <- gsub(x = colnamesHere[neutralLosses], pattern = "^-", replacement = "m")
    colnamesHere[fragments]     <- paste("p", colnamesHere[fragments], sep = "")
    colnames(matrix_train_caret2) <- colnamesHere
    
    ################################################################
    ## do
    #ctrl <- trainControl(
    #  method = "repeatedcv", 
    #  number = 5,
    #  repeats = 1, 
    #  classProbs = TRUE, 
    #  summaryFunction = twoClassSummary, 
    #  verboseIter = FALSE,
    #  allowParallel = TRUE
    #)
    useRoc <- FALSE
    
    ctrl <- caret::trainControl(
      #method = "none",
      method = "repeatedcv", 
      number = 5,
      repeats = 1, 
      classProbs = TRUE, 
      summaryFunction = ifelse(test = useRoc, yes = caret::twoClassSummary, no = caret::prSummary), 
      #summaryFunction = auprcSummary, 
      #summaryFunction = caret::prSummary, 
      #summaryFunction = caret::twoClassSummary, 
      verboseIter = FALSE,
      allowParallel = TRUE,
      savePredictions="final",
      #index=createResample(matrix_train_caret$class, 5)
      index=createResample(class, 5)
    )
    
    #caret::getModelInfo()[[modelName]]
    
#    if(FALSE){
#    fastModels2 <- c(
#      "slda",
#      #"lda2",
#      "lda",
#      "pam",
#      #"simpls",
#      "pls",
#      #"kernelpls",
#      #"widekernelpls",
#      "binda",
#      #"gamLoess",
#      "OneR",
#      "rpart2",
#      "svmLinear",
#      "rpart1SE",
#      "glm",
#      "gam",
#      "C5.0Tree",
#      "C5.0Rules"
#    )
#    }
    #modelNames <- paste(fastModels2, collapse = "|")
    modelNameVector <- strsplit(x = modelNames, split = "\\|")[[1]]
    
    #set.seed(1)
    #startTime <<- Sys.time()
    
    #sink("/dev/null")
    log <- capture.output({
      model_list <- caretList(
        x = matrix_train_caret2,
        y = class,
        #tuneLength = 10,
        trControl = ctrl,
        metric = ifelse(test = useRoc, yes = "ROC", no = "AUC"), 
        #preProc = c("center", "scale")
        preProc = NULL,
        methodList=modelNameVector
      )
    })
    #sink()
    
    ## ROC per model
    #unlist(lapply(X = model_list, FUN = function(x){max(x$results$ROC)}))
    
#    if(FALSE){
#    for(modelName in modelNameVector){
#      print(modelName)
#      tryCatch(expr = {
#        ctrl <- caret::trainControl(
#          #method = "none",
#          method = "repeatedcv", 
#          number = 5,
#          repeats = 1, 
#          classProbs = TRUE, 
#          summaryFunction = caret::twoClassSummary, 
#          verboseIter = FALSE,
#          allowParallel = FALSE
#        )
#        
#        #caret::getModelInfo()[[modelName]]
#        
#        #set.seed(1)
#        #startTime <<- Sys.time()
#        classifier <- caret::train(
#          #class ~ .,
#          #data = matrix_train_caret,
#          x = matrix_train_caret2,
#          y = class,
#          method = modelName,
#          tuneLength = 10,
#          trControl = ctrl,
#          metric = "ROC",
#          #preProc = c("center", "scale")
#          preProc = NULL
#        )
#        print(classifier$results$ROC)
#      },error = function(e){
#        print(paste("error", e))
#      },warning = function(w){
#        print(paste("warning", w))
#      }
#      )
#    }
#    }
    
    
    #endTime <- Sys.time()
    #rocs <- classifier$results$ROC
    
    ## check prediction correlations
    #tmp <- modelCor(resamples(model_list))
    #tmp[tmp < 0.97 & tmp > -0.97] <- 0
    #tmp
    
#    if(FALSE){
#    ## ensembl
#    greedy_ensemble <- caretEnsemble(
#      all.models = model_list, 
#      metric = ifelse(test = useRoc, yes = "ROC", no = "AUC"), 
#      trControl = trainControl(
#        number=5,
#        summaryFunction=twoClassSummary,
#        classProbs=TRUE
#      ))
#    #summary(greedy_ensemble)
#    }
    
    ## stacking
    #sink("/dev/null")
    log <- capture.output({
    stacked_ensemble <- caretStack(
      model_list,
      #method="glm",
      method="gbm",
      #method="rf",
      metric = ifelse(test = useRoc, yes = "ROC", no = "AUC"), 
      trControl=trainControl(
        method="boot",
        number=5,
        savePredictions="final",
        classProbs=TRUE,
        summaryFunction=twoClassSummary
      )
    )
    })
    #sink()
    
    ## Variable Importance
    #varImp(object=model_gbm)
    
    return(stacked_ensemble)
  },
  classify = function(classifier, matrix_test){
    matrix_test_caret <- matrix_test
    matrix_test_caret <- as.data.frame(as.matrix(matrix_test_caret))
    
    ## rename mass column names
    colnamesHere <- colnames(matrix_test_caret)
    neutralLosses <- grepl(x = colnamesHere, pattern = "^\\-\\d+(\\.\\d+)?$")
    fragments     <- grepl(x = colnamesHere, pattern = "^\\d+(\\.\\d+)?$")
    colnamesHere[neutralLosses] <- gsub(x = colnamesHere[neutralLosses], pattern = "^-", replacement = "m")
    colnamesHere[fragments]     <- paste("p", colnamesHere[fragments], sep = "")
    colnames(matrix_test_caret) <- colnamesHere
    
    ## https://www.analyticsvidhya.com/blog/2016/12/practical-guide-to-implement-machine-learning-with-caret-package-in-r-with-practice-problem/
    #predictions<-predict.train(object=model_gbm,testSet[,predictors],type="raw")
    
    #p <- as.data.frame(predict(model_list, newdata=matrix_test_caret))
    #print(p)
    
    
    #plsClasses <- predict(classifier, newdata = matrix_test)
    #sink("/dev/null")
    log <- capture.output({
    plsProbs <- predict(classifier, newdata = matrix_test_caret, type = "prob")
    })
    #sink()
    scores <- plsProbs #pred$predict[, "+", 1]
    
    return(scores)
  }
)

splsda_classifier <- list(
  train = function(matrix_train, classes_pm_train){
#    if(FALSE){
#      posRows <- classes_pm_train=="+"
#      colSumPos <- Matrix::colSums(x = matrix_train[ posRows, ])
#      #colSumNeg <- Matrix::colSums(x = matrix_train[!posRows, ])
#      colFrequencyPos <- colSumPos / sum( posRows)
#      selectedColumns <- colFrequencyPos >= 0.05
#      
#      matrix_train2 <- matrix_train
#      matrix_train2[, !selectedColumns] <- 0
#    }
    
    rownames(matrix_train) <- seq_len(nrow(matrix_train))
    
    #caret_splsda = mixOmics::splsda(X = dataFrame2, Y = groupLabels, ncomp = numberOfComponents, scale = FALSE)
    classifier     <- mixOmics::splsda(X = matrix_train, Y = classes_pm_train, scale = F, ncomp = 1)
    
    return(classifier)
  },
  classify = function(classifier, matrix_test){
    #pred <- predict(classifier, matrix_test, probability=TRUE)
    #scores <- attr(pred, "probabilities")
    #scores <- scores[, "+"]
    
    colnames(matrix_test) <- as.character(colnames(matrix_test))
    rownames(matrix_test) <- seq_len(nrow(matrix_test))
    pred <- predict(classifier, matrix_test, dist = "max.dist")
    scores <- pred$predict[, "+", 1]
    
    return(scores)
  }
)

svm_classifier <- list(
  train = function(matrix_train, classes_pm_train){
    posRows <- classes_pm_train=="+"
    colSumPos <- Matrix::colSums(x = matrix_train[ posRows, ])
    #colSumNeg <- Matrix::colSums(x = matrix_train[!posRows, ])
    colFrequencyPos <- colSumPos / sum( posRows)
    selectedColumns <- colFrequencyPos >= 0.05
    
    matrix_train2 <- matrix_train
    matrix_train2[, !selectedColumns] <- 0
    
    #matrix_train3 <- cbind(matrix_train2, rep(x = 0, times = nrow(matrix_train2)))
    #colnames(matrix_train3)[[ncol(matrix_train3)]] <- "class"
    #matrix_train3[posRows, "class"] <- 1
    
    #classifier <- svm(matrix_train2, factor(classes_pm_train), kernel='linear')
    #parametersHere <- tune(svm, train.x=matrix_train2, train.y=factor(classes_pm_train), kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(.05, .1, .25,.5,1,2)))
    classifier     <- svm(
      x = matrix_train2, y = factor(classes_pm_train), probability=TRUE, 
      kernel="radial"#, cost=parametersHere$best.parameters$cost, gamma=parametersHere$best.parameters$gamma
    )
    
    #trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
    #svm_Linear <- caret::train(class ~., data = matrix_train3, method = "svmLinear", trControl=trctrl, preProcess = c("center", "scale"), tuneLength = 10)
    #svm_model <- svm(class ~ ., data=matrix_train3)
    #svm.model <- svm(matrix_train3[,"class"] ~., data=matrix_train3, kernel="linear")
    
    return(classifier)
  },
  classify = function(classifier, matrix_test){
    pred <- predict(classifier, matrix_test, probability=TRUE)
    scores <- attr(pred, "probabilities")
    scores <- scores[, "+"]
    return(scores)
  }
)
lda_classifier <- list(
  train = function(matrix_train, classes_pm_train){
#    if(FALSE){
#      posRows <- classes_pm_train=="+"
#      colSumPos <- Matrix::colSums(x = matrix_train[ posRows, ])
#      #colSumNeg <- Matrix::colSums(x = matrix_train[!posRows, ])
#      colFrequencyPos <- colSumPos / sum( posRows)
#      selectedColumns <- colFrequencyPos >= 0.05
#      
#      matrix_train2 <- matrix_train[, selectedColumns]
#      #matrix_train2[, !selectedColumns] <- 0
#    }
    #matrix_train3 <- cbind(matrix_train2, rep(x = 0, times = nrow(matrix_train2)))
    #colnames(matrix_train3)[[ncol(matrix_train3)]] <- "class"
    #matrix_train3[posRows, "class"] <- 1
    
    ## train and predict
    ## Error variable 105 appears to be constant within groups
    ## Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'
    #classifier  <- lda(x = matrix_train, grouping=classes_pm_train, method = "t")
    #classifier  <- lda(x = matrix_train, grouping=classes_pm_train)
    
    scalingFactor <- 1
    repeatComputation <- TRUE
    while(repeatComputation){
      classifier  <- tryCatch(
        {
          result <- lda(x = matrix_train * scalingFactor, grouping=classes_pm_train)
          repeatComputation <- FALSE
          return(result)
        }, error = function(e) {
          if(e == "Error in La.svd(x, nu, nv): error code 1 from Lapack routine 'dgesdd'\n"){
            print(paste(e, scalingFactor))
            scalingFactor <<- scalingFactor + 1
          }
        }
      )
    }
    
    #model  <- lda(x = scale(dist(matrix_train)), grouping=classes_pm_train)
    #model2  <- lda(formula=classes_pm_train~., data=as.data.frame(as.matrix(matrix_train)))
    
    return(classifier)
  },
  classify = function(classifier, matrix_test){
    #matrix_test2 <- matrix_test[, selectedColumns]
    prediction   <- predict(object=classifier, newdata=matrix_test)
    #predicted_classes_pm <- as.character(prediction$class)
    
    #scores_pos <- prediction$posterior[, "+"]
    #scores_neg <- prediction$posterior[, "-"]
    
    #if(!ratio){
    #  scores <- scores_pos
    #} else {
    #  scores <- scores_pos / scores_neg
    #}
    
    scores <- prediction$posterior[, "+"]
    
    return(scores)
  }
)

colSumsRatio_classifier <- list(
  train = function(matrix_train, classes_pm_train){
    posRows <- classes_pm_train=="+"
    colSumPos <- Matrix::colSums(x = matrix_train[ posRows, ])
    colSumNeg <- Matrix::colSums(x = matrix_train[!posRows, ])
    
    pseudoCount <- 1
    colSums_PosNeg <- log2( ((colSumPos+pseudoCount) / sum( posRows)) / ((colSumNeg+pseudoCount) / sum(!posRows)) )
    #colSums_PosNeg <- ((colSumPos+pseudoCount) / sum( posRows)) / ((colSumNeg+pseudoCount) / sum(!posRows)) -1
    colSums_PosNeg[colSums_PosNeg < 0] <- 0 ## only pos fragments
    #colSums_PosNeg[abs(colnames(matrix_train)) < 1] <- 0 ## no precursor-NLs
    #names(colSums_PosNeg) <- colnames(matrix_train)
    
    ## storage
    colSums_PosNeg <- unname(colSums_PosNeg)
    
    return(colSums_PosNeg)
  },
  classify = function(classifier, matrix_test){
    ## convert matrix format
    dgTMatrix <- as(matrix_test, "dgTMatrix")
    matrixRows <- dgTMatrix@i
    matrixCols <- dgTMatrix@j
    matrixVals <- dgTMatrix@x
    stMatrix_test <- simple_triplet_matrix(
      i    = matrixRows + 1, 
      j    = matrixCols + 1, 
      v    = matrixVals, 
      nrow = nrow(matrix_test), 
      ncol = ncol(matrix_test)
    )
    
    ## score test matrix row by row
    scores <- rowapply_simple_triplet_matrix(x = stMatrix_test, FUN = function(x){
      sum(x * classifier)
    })
    
    return(scores)
  }
)
colSumsPos_classifier <- list(
  train = function(matrix_train, classes_pm_train){
    posRows <- classes_pm_train=="+"
    colSumPos <- Matrix::colSums(x = matrix_train[ posRows, ])
    colSumNeg <- Matrix::colSums(x = matrix_train[!posRows, ])
    
    colSums_PosNeg <- (colSumPos / sum( posRows)) - (colSumNeg / sum(!posRows))
    colSums_PosNeg[colSums_PosNeg < 0] <- 0 ## only pos fragments
    #colSums_PosNeg[abs(colnames(matrix_train)) < 1] <- 0 ## no precursor-NLs
    #names(colSums_PosNeg) <- colnames(matrix_train)
    
    ## storage
    colSums_PosNeg <- unname(colSums_PosNeg)
    
    return(colSums_PosNeg)
  },
  classify = function(classifier, matrix_test){
    ## convert matrix format
    dgTMatrix <- as(matrix_test, "dgTMatrix")
    matrixRows <- dgTMatrix@i
    matrixCols <- dgTMatrix@j
    matrixVals <- dgTMatrix@x
    stMatrix_test <- simple_triplet_matrix(
      i    = matrixRows + 1, 
      j    = matrixCols + 1, 
      v    = matrixVals, 
      nrow = nrow(matrix_test), 
      ncol = ncol(matrix_test)
    )
    
    ## score test matrix row by row
    scores <- rowapply_simple_triplet_matrix(x = stMatrix_test, FUN = function(x){
      sum(x * classifier)
    })
    
    return(scores)
  }
)
colSumsPosOnly_classifier <- list(
  train = function(matrix_train, classes_pm_train){
    posRows <- classes_pm_train=="+"
    colSumPos <- Matrix::colSums(x = matrix_train[ posRows, ])
    #colSumNeg <- Matrix::colSums(x = matrix_train[!posRows, ])
    
    colSums_PosNeg <- (colSumPos / sum( posRows))
    #colSums_PosNeg <- (colSumPos / sum( posRows)) - (colSumNeg / sum(!posRows))
    #colSums_PosNeg[colSums_PosNeg < 0] <- 0 ## only pos fragments
    #colSums_PosNeg[abs(colnames(matrix_train)) < 1] <- 0 ## no precursor-NLs
    #names(colSums_PosNeg) <- colnames(matrix_train)
    
    ## storage
    colSums_PosNeg <- unname(colSums_PosNeg)
    
    return(colSums_PosNeg)
  },
  classify = function(classifier, matrix_test){
    ## convert matrix format
    dgTMatrix <- as(matrix_test, "dgTMatrix")
    matrixRows <- dgTMatrix@i
    matrixCols <- dgTMatrix@j
    matrixVals <- dgTMatrix@x
    stMatrix_test <- simple_triplet_matrix(
      i    = matrixRows + 1, 
      j    = matrixCols + 1, 
      v    = matrixVals, 
      nrow = nrow(matrix_test), 
      ncol = ncol(matrix_test)
    )
    
    ## score test matrix row by row
    scores <- rowapply_simple_triplet_matrix(x = stMatrix_test, FUN = function(x){
      sum(x * classifier)
    })
    
    return(scores)
  }
)

## \text{score}(s) = \sum_{f \in s} P(f|s\text{ is in class}) - P(f|s\text{ is in not class})
## P(f|s\text{ is in class}) = \displaystyle \frac{\displaystyle\sum_{s\in S_+} \delta(f \in s)}{|S_+|}
## P(f|s\text{ is not in class}) = \displaystyle \frac{\displaystyle\sum_{s\in S_-} \delta(f \in s)}{|S_-|}
colSums_classifier <- list(
  train = function(matrix_train, classes_pm_train){
    posRows <- classes_pm_train=="+"
    colSumPos <- Matrix::colSums(x = matrix_train[ posRows, ])
    colSumNeg <- Matrix::colSums(x = matrix_train[!posRows, ])
    
    colSums_PosNeg <- (colSumPos / sum( posRows)) - (colSumNeg / sum(!posRows))
    #colSums_PosNeg[abs(colnames(matrix_train)) < 1] <- 0 ## no precursor-NLs
    #names(colSums_PosNeg) <- colnames(matrix_train)
    
    ## storage
    colSums_PosNeg <- unname(colSums_PosNeg)
    
    return(colSums_PosNeg)
  },
  classify = function(classifier, matrix_test){
    ## convert matrix format
    dgTMatrix <- as(matrix_test, "dgTMatrix")
    matrixRows <- dgTMatrix@i
    matrixCols <- dgTMatrix@j
    matrixVals <- dgTMatrix@x
    stMatrix_test <- simple_triplet_matrix(
      i    = matrixRows + 1, 
      j    = matrixCols + 1, 
      v    = matrixVals, 
      nrow = nrow(matrix_test), 
      ncol = ncol(matrix_test)
    )
    
    ## score test matrix row by row
    scores <- rowapply_simple_triplet_matrix(x = stMatrix_test, FUN = function(x){
      sum(x * classifier)
    })
    
    return(scores)
  }
)
colSumsLog_classifier <- list(
  train = function(matrix_train, classes_pm_train){
    posRows <- classes_pm_train=="+"
    colSumPos <- Matrix::colSums(x = matrix_train[ posRows, ]) + max(.Machine$double.neg.eps, .Machine$double.eps)#1/sum( posRows)#1E-6#1
    colSumNeg <- Matrix::colSums(x = matrix_train[!posRows, ]) + max(.Machine$double.neg.eps, .Machine$double.eps)#1/sum(!posRows)#1E-6#1
    
    colSums_PosNeg <- log(colSumPos / sum( posRows)) - log(colSumNeg / sum(!posRows))
    #colSums_PosNeg[abs(colnames(matrix_train)) < 1] <- 0 ## no precursor-NLs
    #names(colSums_PosNeg) <- colnames(matrix_train)
    
    return(colSums_PosNeg)
  },
  classify = function(classifier, matrix_test){
    ## convert matrix format
    dgTMatrix <- as(matrix_test, "dgTMatrix")
    matrixRows <- dgTMatrix@i
    matrixCols <- dgTMatrix@j
    matrixVals <- dgTMatrix@x
    stMatrix_test <- simple_triplet_matrix(
      i    = matrixRows + 1, 
      j    = matrixCols + 1, 
      v    = matrixVals, 
      nrow = nrow(matrix_test), 
      ncol = ncol(matrix_test)
    )
    
    ## score test matrix row by row
    scores <- rowapply_simple_triplet_matrix(x = stMatrix_test, FUN = function(x){
      sum(x * classifier)
    })
    
    return(scores)
  }
)
tensorflow_classifier <- list(
  train = function(matrix_train, classes_pm_train, topology){
    ####################################################
    ## parameters
    printProgress <- FALSE
    
    ## model parameters
    numberOfMS2PeakGroups <- ncol(matrix_train)
    
    #numLayers = 0
    #numNeurons_l = NULL
    #dropout_l = NULL
    numLayers        = topology$numLayers
    numNeurons_l = topology$numNeurons_l
    dropout_l         = topology$dropout_l
    numberOfInputNeurons = numberOfMS2PeakGroups
    numberOfOutputNeurons = 2
    
    ## train parameters
    meanCenter        <- FALSE
    maximumNumberOfTrainingIterationsWithoutImprovement <- 10
    maximumNumberOfTrainingIterations <- 1000
    learningRate <- 0.2
    
    numberOfSpectra_train     <- length(classes_pm_train)
    numberOfMS2PeakGroupsHere <- ncol(matrix_train)
    matrix_train     <- as.matrix(matrix_train)
    
    if(meanCenter)
      matrix_train <- as.matrix(apply(X = matrix_train, MARGIN = 2, FUN = function(x){
        x - mean(x = x)
      }))
    
    ####################################################
    ## class weights
    negToPosRatio <- sum(classes_pm_train=="-") / sum(classes_pm_train=="+")
    
    resultMultiplier <- vector(mode = "numeric", length = numberOfSpectra_train)
    resultMultiplier[classes_pm_train=="+"] <- negToPosRatio
    resultMultiplier[classes_pm_train=="-"] <- 1.
    
    ## classes with weigths each of 1
    classes_pm_train2 <- integer(length = numberOfSpectra_train)
    classes_pm_train2[classes_pm_train == "+"] <- 1L
    classes_pm_train2[classes_pm_train == "-"] <- 0L
    classes_pm_train2 <- matrix(data = c(classes_pm_train2, 1L - classes_pm_train2), nrow = numberOfSpectra_train)
    
    ####################################################
    ## Create the model
    #variableScope <- tf$VariableScope(reuse = TRUE)
    variableScope <- tf$VariableScope(reuse = FALSE)
    tf$reset_default_graph()
    
    modelList <- createTFModel(numLayers = numLayers, numNeurons_l = numNeurons_l, dropout_l = dropout_l, numberOfInputNeurons = numberOfInputNeurons, numberOfOutputNeurons = numberOfOutputNeurons, variableScope=variableScope)
    dropout <- modelList$dropout
    x  <- modelList$x
    y  <- modelList$y
    y_ <- modelList$y_
    dropoutTensorNames <- modelList$dropoutTensorNames
    multiLayerPerceptronStructure <- modelList$multiLayerPerceptronStructure
    
    ##########################
    ## define accuracy and auc
    correct_prediction <- tf$equal(tf$argmax(y, 1L), tf$argmax(y_, 1L), name = "correct_prediction")
    accuracy <- tf$reduce_mean(tf$cast(correct_prediction, tf$float32), name = "accuracy")
    auc <- tf$metrics$auc(labels = y_, predictions = y, name = "AUC")[[2]]
    
    ##########################
    ## Define loss and optimizer
    resultMult  <- tf$constant(value = c(negToPosRatio, 1), dtype = tf$float32, name = "many_to_one_multiplier")
    cross_entropy <- tf$reduce_mean(-tf$reduce_sum(tf$multiply(tf$multiply(y_, log(y)), resultMult), reduction_indices=1L), name = "weighted_cross_entropy")
    multiLayerPerceptronStructure[length(multiLayerPerceptronStructure)+1] <- paste("loss_function", "weighted_cross_entropy", sep = "=")
    
    ##########################
    # regularization: This is a good beta value to start with
#    if(FALSE){
#      regularizer_beta = 0.01
#      regularizer = tf$nn$l2_loss(tf$concat(values = list(inputToHidden_weights,tf$transpose(hiddenToOutput_weights)), axis = 0L))
#      cross_entropy = tf$reduce_mean(cross_entropy + regularizer_beta * regularizer)
#    }
    
    #train_step <- tf$train$AdamOptimizer(learning_rate=learningRate, name = "AdamOptimizer")$minimize(loss = cross_entropy)
    train_step <- tf$train$GradientDescentOptimizer(learning_rate = learningRate, name = "GradientDescentOptimizer")$minimize(loss = cross_entropy)
    multiLayerPerceptronStructure[length(multiLayerPerceptronStructure)+1] <- paste("optimizer", "GradientDescentOptimizer", sep = "=")
    multiLayerPerceptronStructure[length(multiLayerPerceptronStructure)+1] <- paste("learning_rate", learningRate, sep = "=")
    
    ####################################################
    ## Create session and initialize variables
    #if(!exists(x = "sess"))   sess <- tf$InteractiveSession()
    #if(!exists(x = "sess"))   sess <- tf$Session()
    sess <- tf$Session()
    
    sess$run(tf$global_variables_initializer())
    sess$run(tf$local_variables_initializer())
    
    bestAccuracy <- -1
    iterationOfBestAccuracy <- -1
    startTime <- Sys.time()
    accuracies <- numeric()
    
    dropoutTensorNames <- paste(dropoutTensorNames, 0, sep = ":")
    #list1 <- list(x = matrix_train)
    #list2 <- list(x = matrix_train, y_ = classes_pm_train2)
    list1 <- list()
    list2 <- list()
    list1[[x$name ]] = matrix_train
    list2[[x$name ]] = matrix_train
    list2[[y_$name]] = classes_pm_train2
    if(dropout){
      list3 <- list2
      list1[dropoutTensorNames] <- 1.0
      list2[dropoutTensorNames] <- 1.0
      list3[dropoutTensorNames] <- 0.5
      
      ##feed_dict_train  <- dict(x = matrix_train, y_ = classes_pm_train2, "hidden_layer0_dropout_tensor:0" = tf$constant(value = 1.0, dtype = tf$float32, shape = shape(1), name = "dropout_tmp"))
      ##feed_dict_train  <- dict(x = matrix_train, y_ = classes_pm_train2, "hidden_layer0_dropout_tensor:0" = 1.0)
      #feed_dict_train  <- dict(x = matrix_train, y_ = classes_pm_train2, keep_prob = 1.0)
      #feed_dict_train2 <- dict(x = matrix_train, keep_prob = 1.0)
      #feed_dict_train3 <- dict(x = matrix_train, y_ = classes_pm_train2, keep_prob = 0.5)
      feed_dict_train  <- dict(list2)
      feed_dict_train2 <- dict(list1)
      feed_dict_train3 <- dict(list3)
    } else {
      feed_dict_train  <- dict(list2)
      feed_dict_train2 <- dict(list1)
      feed_dict_train3 <- feed_dict_train
    }
    
    ####################################################
    ## train
    for (i in 1:maximumNumberOfTrainingIterations) {
      train_accuracy       <- accuracy$eval(feed_dict = feed_dict_train, session = sess)
      train_cross_entropy <- cross_entropy$eval(feed_dict = feed_dict_train, session = sess)
      train_auc           <- auc$eval(feed_dict = feed_dict_train, session = sess)
      accuracies[[length(accuracies) + 1]] <- train_accuracy
      if(train_accuracy > bestAccuracy){
        bestAccuracy <- train_accuracy
        iterationOfBestAccuracy <- i
      } else {
        accuracyVariance <- var(tail(accuracies, 10))
        #if(is.na(accuracyVariance)) accuracyVariance <- 1.
        if(i - iterationOfBestAccuracy > maximumNumberOfTrainingIterationsWithoutImprovement & accuracyVariance < 1e-6)
          break;
      }
      if(printProgress)  cat(paste(
        classIdx, "/", numberOfClasses,
        "class", class, 
        "; step", i, 
        "; accuracy tr",  round(train_accuracy, digits = 4), 
        "; cross entropy", round(train_cross_entropy, digits = 4),
        "; auc", round(train_auc, digits = 4),
        ifelse(iterationOfBestAccuracy == i, " Increase", ""), 
        "\n"
      ))
      
      train_step$run(feed_dict = feed_dict_train, session = sess)
    }
    endTime <- Sys.time()
    time <- difftime(endTime, startTime, units = "secs")[[1]]
    
    multiLayerPerceptronStructure[length(multiLayerPerceptronStructure)+1] <- paste("number_of_training_steps", i, sep = "=")
    multiLayerPerceptronStructure[length(multiLayerPerceptronStructure)+1] <- paste("training_accuracy", bestAccuracy, sep = "=")
    multiLayerPerceptronStructure[length(multiLayerPerceptronStructure)+1] <- paste("training_time(s)", time, sep = "=")
    
    ## close session
    #sess$close()
    
    modelList$multiLayerPerceptronStructure <- multiLayerPerceptronStructure
    
    resultList <- list(
      session = sess,
      modelList = modelList
    )
    
    return(resultList)
  },
  classify = function(classifier, matrix_test){
    list <- classifier
    
    sess <- list$session
    dropout <- list$modelList$dropout
    y <- list$modelList$y
    x <- list$modelList$x
    dropoutTensorNames <- list$modelList$dropoutTensorNames
    
    numberOfSpectra_test  <- nrow(matrix_test)
    numberOfMS2PeakGroupsHere <- ncol(matrix_test)
    
    meanCenter        <- FALSE
    matrix_test      <- as.matrix(matrix_test)
    
    if(meanCenter)
      matrix_test  <- as.matrix(apply(X = matrix_test,  MARGIN = 2, FUN = function(x){
        x - mean(x = x)
      }))
    
    
    dropoutTensorNames <- paste(dropoutTensorNames, 0, sep = ":")
    
    listTmp <- list()
    listTmp[[x$name]] = matrix_test
    if(dropout){
      listTmp[dropoutTensorNames] <- 1.0
      feed_dict_test2  <- dict(listTmp)
    } else {
      feed_dict_test2  <- dict(listTmp)
    }
    
    ####################################################
    ## metrics
    #correct_predictions <- correct_prediction$eval(feed_dict = feed_dict_test, session = sess)
    probabilities <- y$eval(feed_dict = feed_dict_test2, session = sess)
    #predictions_test <- tf$argmax(input = y, axis = 1L, name = "predictions_test")$eval(feed_dict = feed_dict_test2, session = sess)#, keep_prob = 1.0), session = sess)
    #predictions_test <- as.integer(predictions_test)
    
    predicted_scores <- probabilities[, 1]
    #predicted_scores <- probabilities[, 1] - probabilities[, 2]
    #predicted_scores <- probabilities[, 1] / probabilities[, 2]
    
    sess$close()
    
    return(predicted_scores)
  }
)

######################################################################################
## classifiers_old
predict_LDA <- function(matrix_train, classes_pm_train, matrix_test, classes_pm_test, ratio){
  ## train and predict
  model  <- lda(x = matrix_train, grouping=classes_pm_train)
  #model  <- lda(x = scale(dist(matrix_train)), grouping=classes_pm_train)
  #model2  <- lda(formula=classes_pm_train~., data=as.data.frame(as.matrix(matrix_train)))
  prediction   <- predict(object=model, newdata=matrix_test)
  predicted_classes_pm <- as.character(prediction$class)
  
  scores_pos <- prediction$posterior[, "+"]
  scores_neg <- prediction$posterior[, "-"]
  
  if(!ratio){
    scores <- scores_pos
  } else {
    scores <- scores_pos / scores_neg
  }
  
  return(scores)
  #return(predicted_classes_pm)
}
predict_Correlation <- function(matrix_train, classes_pm_train, matrix_test, classes_pm_test, ratio, corMethod = c("pearson", "kendall", "spearman"), linkage = c("single", "average", "centroid")){
  if(corMethod == "kendall")
    stop("Operation not supported")
  
  correlations <- cor(x = t(as.matrix(matrix_test)), y = t(as.matrix(matrix_train)), method = corMethod)
  
#  if(FALSE){
#    ## classes
#    switch(linkage,
#           "single"={
#             ## single linkage
#             predicted_classes_pm <- unlist(apply(X = correlations, MARGIN = 1, FUN = function(x){
#               if(all(is.na(x)))
#                 return("-")
#               else{
#                 maxIdx <- which.max(x)
#                 class <- classes_pm_train[maxIdx]
#                 return(class)
#               }
#             }))
#           },
#           "average"={
#             ## average linkage
#             predicted_classes_pm <- unlist(apply(X = correlations, MARGIN = 1, FUN = function(x){
#               if(all(is.na(x)))
#                 return("-")
#               else{
#                 plus  <- x[classes_pm_train=="+"]
#                 minus <- x[classes_pm_train=="-"]
#                 if(mean(plus, na.rm = TRUE) > mean(minus, na.rm = TRUE))
#                   #if(max(plus) > max(minus))
#                   #if(mean(plus) > mean(minus))
#                   return("+")
#                 else
#                   return("-")
#               }
#             }))
#           },
#           "centroid"={
#             ## centroid linkage
#             distPlus  <- as.matrix(dist(matrix_train[classes_pm_train=="+", ]))
#             distMinus <- as.matrix(dist(matrix_train[classes_pm_train=="-", ]))
#             centroidPlus  <- which.min(apply(X = distPlus , MARGIN = 1, FUN = sum))
#             centroidMinus <- which.min(apply(X = distMinus, MARGIN = 1, FUN = sum))
#             centroidPlus  <- which(classes_pm_train=="+")[[centroidPlus]]
#             centroidMinus <- which(classes_pm_train=="-")[[centroidMinus]]
#             
#             predicted_classes_pm <- unlist(apply(X = correlations, MARGIN = 1, FUN = function(x){
#               if(all(is.na(x)))
#                 return("-")
#               else{
#                 #if(x[[centroidPlus]] > x[[centroidMinus]])
#                 if(x[[centroidPlus]] > mean(x[classes_pm_train=="-"], na.rm = TRUE))
#                   #if(max(plus) > max(minus))
#                   #if(mean(plus) > mean(minus))
#                   return("+")
#                 else
#                   return("-")
#               }
#             }))
#           },
#           stop(paste("Unknown linkage (", linkage, ")!", sep = ""))
#    )
#  }
  
  posItems <- classes_pm_train=="+"
  negItems <- classes_pm_train=="-"
  switch(linkage,
         "single"={
           ## single linkage
           scoresPosNeg <- apply(X = correlations, MARGIN = 1, FUN = function(x){
             if(all(is.na(x)))
               return(c("+" = 0, "-" = 1))
             else{
               maxPos <- max(x[posItems], na.rm = TRUE)
               maxNeg <- max(x[negItems], na.rm = TRUE)
               return(c("+" = maxPos, "-" = maxNeg))
             }
           })
         },
         "average"={
           ## average linkage
           scoresPosNeg <- apply(X = correlations, MARGIN = 1, FUN = function(x){
             if(all(is.na(x)))
               return(c("+" = 0, "-" = 1))
             else{
               meanPos <- mean(posItems, na.rm = TRUE)
               meanNeg <- mean(negItems, na.rm = TRUE)
               return(c("+" = meanPos, "-" = meanNeg))
             }
           })
         },
         "centroid"={
           ## centroid linkage
           distPlus  <- as.matrix(dist(matrix_train[posItems, ]))
           #distMinus <- as.matrix(dist(matrix_train[negItems, ]))
           centroidPlus  <- which.min(apply(X = distPlus , MARGIN = 1, FUN = sum))
           #centroidMinus <- which.min(apply(X = distMinus, MARGIN = 1, FUN = sum))
           centroidPlus  <- which(posItems)[[centroidPlus]]
           #centroidMinus <- which(negItems)[[centroidMinus]]
           
           scoresPosNeg <- apply(X = correlations, MARGIN = 1, FUN = function(x){
             if(all(is.na(x)))
               return(c("+" = 0, "-" = 1))
             else{
               centroidPos <- x[[centroidPlus]]
               meanNeg <- mean(x[negItems], na.rm = TRUE)
               return(c("+" = centroidPos, "-" = meanNeg))
             }
           })
         },
         stop(paste("Unknown linkage (", linkage, ")!", sep = ""))
  )
  
  scores_pos <- scoresPosNeg["+", ]
  scores_neg <- scoresPosNeg["-", ]
  
  if(!ratio){
    scores <- scores_pos
  } else {
    scores <- scores_pos / scores_neg
  }
  
  return(scores)
  #return(predicted_classes_pm)
}
predict_RDA <- function(matrix_train, classes_pm_train, matrix_test, classes_pm_test){
  stop("Operation not supported")
  
  ## train and predict
  model  <- rda(x = matrix_train, grouping=classes_pm_train)
  prediction   <- predict(object=model, newdata=matrix_test)
  predicted_classes_pm <- as.character(prediction$class)
  
  return(predicted_classes_pm)
}
predict_SVM <- function(matrix_train, classes_pm_train, matrix_test, classes_pm_test, ratio){
  model <- svm(x = matrix_train, y = as.factor(classes_pm_train), scale = F, probability = T)
  prediction <- predict(model, matrix_test, decision.values=TRUE, probability=TRUE)
  probabilities <- attr(prediction, "probabilities")
  #predicted_classes_pm <- unname(apply(probabilities, 1, function(x) names(which.max(x))))
  
  scores_pos <- probabilities[, "+"]
  scores_neg <- probabilities[, "-"]
  
  if(!ratio){
    scores <- scores_pos
  } else {
    scores <- scores_pos / scores_neg
  }
  
  return(scores)
  #return(predicted_classes_pm)
}
predict_SOM <- function(matrix_train, classes_pm_train, matrix_test, classes_pm_test){
  stop("Operation not supported")
  model <- som(X = scale(matrix_train), grid=somgrid(length(unique(classes_pm_train)),2,"hexagonal"), rlen=100)
  prediction <- predict(model, newdata=scale(matrix_test), trainX=scale(matrix_train), trainY=as.factor(classes_pm_train))
  predicted_classes_pm <- as.character(prediction$prediction)
  
  return(predicted_classes_pm)
}
predict_XYF <- function(matrix_train, classes_pm_train, matrix_test, classes_pm_test){
  stop("Operation not supported")
  model <- xyf(X = scale(matrix_train), Y = as.factor(classes_pm_train), grid=somgrid(length(unique(classes_pm_train)),2,"hexagonal"), rlen=100)
  prediction <- predict(model, newdata=scale(matrix_test))
  predicted_classes_pm <- as.character(prediction$prediction)
  
  return(predicted_classes_pm)
}
predict_NeuralNet <- function(matrix_train, classes_pm_train, matrix_test, classes_pm_test){
  stop("Operation not supported")
  #nn_df <- cbind(matrix_train, classes_pm_train)
  model <- nnet(classes_pm_train ~ . , data=as.matrix(matrix_train), size=10, rang=0.6, decay=5e-4, maxit=10, MaxNWts=80000)
  prediction <- predict(model, newdata=matrix_train, type="class")
  predicted_classes_pm <- as.character(prediction)
  
  return(predicted_classes_pm)
}

predict_ColSums <- function(matrix_train, classes_pm_train, matrix_test, classes_pm_test){
  posRows <- classes_pm_train=="+"
  #colSumPos <- apply(X = matrix_train[ posRows, ], MARGIN = 2, FUN = sum)
  #colSumNeg <- apply(X = matrix_train[!posRows, ], MARGIN = 2, FUN = sum)
  colSumPos <- Matrix::colSums(x = matrix_train[ posRows, ])
  colSumNeg <- Matrix::colSums(x = matrix_train[!posRows, ])
  
  colSums <- (colSumPos / sum( posRows)) - (colSumNeg / sum(!posRows))
  
#  if(FALSE){
#    ## top ten
#    names(colSums) <- colnames(matrix_train)
#    tail(x = sort(colSums), n = 10)
#  }
  
  #scores <- apply(X = matrix_test, MARGIN = 1, FUN = function(x){
  #  sum(x * colSums)
  #})
  
  dgTMatrix <- as(matrix_test, "dgTMatrix")
  matrixRows <- dgTMatrix@i
  matrixCols <- dgTMatrix@j
  matrixVals <- dgTMatrix@x
  stMatrix <- simple_triplet_matrix(i = matrixRows + 1, j = matrixCols + 1, v = matrixVals, nrow=nrow(matrix_test), ncol=ncol(matrix_test))
  
  scores <- rowapply_simple_triplet_matrix(x = stMatrix, FUN = function(x){
    sum(x * colSums)
  })
  
  return(scores)
}
predict_CosinusDistance <- function(matrix_train, classes_pm_train, matrix_test, classes_pm_test, ratio){
  ## sum(a*b) / (sqrt(sum(a*a)) * sqrt(sum(b*b)))
  posRows <- classes_pm_train=="+"
  colSumPos <- apply(X = matrix_train[ posRows, ], MARGIN = 2, FUN = sum)
  colSumNeg <- apply(X = matrix_train[!posRows, ], MARGIN = 2, FUN = sum)
  
  scores_pos <- apply(X = matrix_test, MARGIN = 1, FUN = function(x){
    sum(x * colSumPos) / (sqrt(sum(x * x)) * sqrt(sum(colSumPos * colSumPos)))
  })
  scores_neg <- apply(X = matrix_test, MARGIN = 1, FUN = function(x){
    sum(x * colSumNeg) / (sqrt(sum(x * x)) * sqrt(sum(colSumNeg * colSumNeg)))
  })
  
  if(!ratio){
    scores <- scores_pos
  } else {
    scores <- scores_pos / scores_neg
  }
  
  return(scores)
}
predict_Prod <- function(matrix_train, classes_pm_train, matrix_test, classes_pm_test, ratio){
  ## sum(a*norm(b))
  
  posRows <- classes_pm_train=="+"
  colSumPos <- apply(X = matrix_train[ posRows, ], MARGIN = 2, FUN = sum)
  colSumNeg <- apply(X = matrix_train[!posRows, ], MARGIN = 2, FUN = sum)
  
  colSumPos <- colSumPos / sum(colSumPos)
  colSumNeg <- colSumNeg / sum(colSumNeg)
  
  scores_pos <- apply(X = matrix_test, MARGIN = 1, FUN = function(x){
    sum(x * colSumPos)
  })
  scores_neg <- apply(X = matrix_test, MARGIN = 1, FUN = function(x){
    sum(x * colSumNeg)
  })
  
  if(!ratio){
    scores <- scores_pos
  } else {
    scores <- scores_pos / scores_neg
  }
  
  return(scores)
}
predict_Jaccard <- function(matrix_train, classes_pm_train, matrix_test, classes_pm_test, ratio){
  ## sum(intersection(a,b))/sum(union(a,b))
  
  posRows <- classes_pm_train=="+"
  colSumPos <- apply(X = matrix_train[ posRows, ], MARGIN = 2, FUN = sum)
  colSumNeg <- apply(X = matrix_train[!posRows, ], MARGIN = 2, FUN = sum)
  
  scores_pos <- apply(X = matrix_test, MARGIN = 1, FUN = function(x){
    intersection <- colSumPos > 0 & x > 0
    union        <- colSumPos > 0 | x > 0
    sum(intersection) / sum(union)
  })
  scores_neg <- apply(X = matrix_test, MARGIN = 1, FUN = function(x){
    intersection <- colSumNeg > 0 & x > 0
    union        <- colSumNeg > 0 | x > 0
    sum(intersection) / sum(union)
  })
  
  if(!ratio){
    scores <- scores_pos
  } else {
    scores <- scores_pos / scores_neg
  }
  
  return(scores)
}
predict_JaccardWeighted <- function(matrix_train, classes_pm_train, matrix_test, classes_pm_test, ratio){
  ## sum(intersection(aW,bW))/sum(union(aW,bW))
  
  posRows <- classes_pm_train=="+"
  colSumPos <- apply(X = matrix_train[ posRows, ], MARGIN = 2, FUN = sum)
  colSumNeg <- apply(X = matrix_train[!posRows, ], MARGIN = 2, FUN = sum)
  
  scores_pos <- apply(X = matrix_test, MARGIN = 1, FUN = function(x){
    intersection <- colSumPos > 0 & x > 0
    union        <- colSumPos > 0 | x > 0
    sum(colSumPos[intersection]) / sum(colSumPos[union])
  })
  scores_neg <- apply(X = matrix_test, MARGIN = 1, FUN = function(x){
    intersection <- colSumNeg > 0 & x > 0
    union        <- colSumNeg > 0 | x > 0
    sum(colSumNeg[intersection]) / sum(colSumNeg[union])
  })
  
  if(!ratio){
    scores <- scores_pos
  } else {
    scores <- scores_pos / scores_neg
  }
  
  return(scores)
}

classify_lda <- function(){
  # Select compounds that separates the groups most using Wilk's lambda criterion
  model_wilks <- greedy.wilks(X=feat_list, grouping=species, niveau=0.0001, na.action=na.omit)
  model_wilks$formula
  
  # PLDA using Wilks
  model_plda <- lda(formula=model_wilks$formula, data=feat_list)
  model_plda <- predict(object=model_plda, newdata=feat_list)
  
  pdf(paste(f.pol_char(),"/","species/plda_wilks.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
  plot(model_plda$x[,c(1:2)], type="n")
  mtext(side=3, line=2, "Species: PLDA using Wilks lambda criterion", cex=1, font=2)
  points(model_plda$x[,c(1:2)], pch=16, col=species_samples_colors, cex=1.4)
  #text(model_plda$x[,c(1:2)], labels=rownames(model_plda$x), col=species_samples_colors, cex=0.4)
  legend("topleft", bty="n", pch=rep(16,length(levels(species))), col=species_colors,
         pt.cex=1, cex=0.8, y.intersp=0.8, text.width=0.5, legend=levels(species)) 
  dev.off()
}
