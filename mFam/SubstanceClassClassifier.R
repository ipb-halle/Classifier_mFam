
##############################################################################################################################################
## main program
runTest <- function(parameterSetAll){
  ##################################################################################################
  ## parameters
  
  ##############################################################################
  ## general parameters
  
  ## I/O
  progress               = as.logical(parameterSetAll$progress)
  outDetail              = as.logical(parameterSetAll$outDetail)
  resultFolderForClassifiers           = parameterSetAll$resultFolderForClassifiers
  resultFolderForMetfamilyProjects     = parameterSetAll$resultFolderForMetfamilyProjects
  writeClassifiers       = as.logical(parameterSetAll$writeClassifiers)
  writeResults           = as.logical(parameterSetAll$writeResults)
  outSuffix              = parameterSetAll$outSuffix
  maximumNumberOfScores  = as.integer(parameterSetAll$maximumNumberOfScores)
  removeRareFragments    = as.logical(parameterSetAll$removeRareFragments)
  builtClassifiers       = as.logical(parameterSetAll$builtClassifiers)
  writeMetFamilyProjects = as.logical(parameterSetAll$writeMetFamilyProjects)
  calculateRMDstatistics = as.logical(parameterSetAll$calculateRMDstatistics)
  ## classes
  mergeDuplicatedSpectra              = as.logical(parameterSetAll$mergeDuplicatedSpectra)
  takeSpectraWithMaximumNumberOfPeaks = as.logical(parameterSetAll$takeSpectraWithMaximumNumberOfPeaks)
  #considerAlternativeSubstanceClasses = parameterSetAll$considerAlternativeSubstanceClasses
  #processSubstanceClasses             = parameterSetAll$processSubstanceClasses
  #thisSubstanceClass                  = parameterSetAll$thisSubstanceClass
  #thisSubstituent                     = parameterSetAll$thisSubstituent
  classOfClass = parameterSetAll$classOfClass
  thisClass = parameterSetAll$thisClass
  ## repeated random subsampling validation
  minimumNumberOfPosSpectraPerClass = as.integer(parameterSetAll$minimumNumberOfPosSpectraPerClass)
  minimumNumberOfNegSpectraPerClass = as.integer(parameterSetAll$minimumNumberOfNegSpectraPerClass)
  numberOfDataSetDecompositions     = as.integer(parameterSetAll$numberOfDataSetDecompositions)
  proportionTraining                = as.numeric(parameterSetAll$proportionTraining)
  fragmentColumnSelectionMethod     = parameterSetAll$fragmentColumnSelectionMethod
  minimumProportionOfPositiveFragments = as.numeric(parameterSetAll$minimumProportionOfPositiveFragments)
  ## postprocessing
  unitResolution                    = as.logical(parameterSetAll$unitResolution)
  #######################################################################
  ## cache
  reParseMsp              = as.logical(parameterSetAll$reParseMsp)
  reProcessAnno           = as.logical(parameterSetAll$reProcessAnno)
  reProcessLibrary        = as.logical(parameterSetAll$reProcessLibrary)
  reProcessFragmentMatrix = as.logical(parameterSetAll$reProcessFragmentMatrix)
  
  #######################################################################
  ## constraints
  
  ## library one
  annoFile    = parameterSetAll$annoFile
  thisLibrary = parameterSetAll$thisLibrary
  thisMethod  = parameterSetAll$thisMethod
  
  thisParameterSet <- NULL
  if("thisParameterSet" %in% names(parameterSetAll)) thisParameterSet <- parameterSetAll$thisParameterSet
  splitClasses = FALSE
  if("splitClasses" %in% names(parameterSetAll)) splitClasses <- parameterSetAll$splitClasses
  splitClassesParameterSet <- parameterSetAll$splitClassesParameterSet
  if("computeFragmentFisherTest" %in% names(parameterSetAll)) computeFragmentFisherTest <- parameterSetAll$computeFragmentFisherTest
  computeFragmentFisherTest <- FALSE
  
  
  
  
  #######################################################################
  ## input files and MS parameters
  #returnObj <- getLibraryData()
  #libraries                 <- returnObj$libraries
  #annoFiles                 <- returnObj$annoFiles
  #parameterSets             <- returnObj$parameterSets
  #allowedInstrumentTypeSets <- returnObj$allowedInstrumentTypeSets
  libraries                 <- thisLibrary
  annoFiles                 <- annoFile
  parameterSets             <- highResolutionParameterSet
  allowedInstrumentTypeSets <- "all"
  
  #librariesOriginal <- libraries
  #annoFilesOriginal <- annoFiles
  #parameterSetsOriginal <- parameterSets
  #allowedInstrumentTypeSetsOriginal <- allowedInstrumentTypeSets
  
  if(unitResolution)
    allowedInstrumentTypeSets <- as.list(rep(x = "all", times = length(allowedInstrumentTypeSets)))
  
  
  
  #######################################################################
  ## sanity checks
  if(all(mergeDuplicatedSpectra, takeSpectraWithMaximumNumberOfPeaks))
    stop("error")
  if(!("thisParameterSet" %in% names(parameterSetAll))) stop("No parameterSet provided")
  if(classOfClass %in% c("ChEBI_Identified", "ChEBI_Predicted"))
    stop(classOfClass)
  if(is.null(thisLibrary)) stop("No library provided")
  
  print("============================================================ BREAKPOINT")
  print(length(parameterSetAll))
  print(length(libraries))
  print(length(annoFiles))
  print(length(parameterSets))
  print(length(allowedInstrumentTypeSets))
  
  ## sanity checks
##  if(
##    length(libraries) != length(annoFiles) | 
##    length(libraries) != length(parameterSets) | 
##    length(libraries) != length(allowedInstrumentTypeSets)
##  )
##    stop("Lengths not equal")
  
  #######################################################################
  ## print params
  print("### folders")
  print(paste("resultFolderForClassifiers: ", resultFolderForClassifiers, sep = ""))
  print(paste("resultFolderForMetfamilyProjects: ", resultFolderForMetfamilyProjects, sep = ""))
  ## I/O
  print("### IO")
  print(paste("progress: ", progress, sep = ""))
  print(paste("outDetail: ", outDetail, sep = ""))
  print(paste("writeClassifiers: ", writeClassifiers, sep = ""))
  print(paste("writeResults: ", writeResults, sep = ""))
  print(paste("outSuffix: ", outSuffix, sep = ""))
  print(paste("maximumNumberOfScores: ", maximumNumberOfScores, sep = ""))
  print(paste("removeRareFragments: ", removeRareFragments, sep = ""))
  print(paste("builtClassifiers: ", builtClassifiers, sep = ""))
  print(paste("writeMetFamilyProjects: ", writeMetFamilyProjects, sep = ""))
  print(paste("calculateRMDstatistics: ", calculateRMDstatistics, sep = ""))
  ## classes
  print("### classes")
  print(paste("mergeDuplicatedSpectra: ", mergeDuplicatedSpectra, sep = ""))
  print(paste("takeSpectraWithMaximumNumberOfPeaks: ", takeSpectraWithMaximumNumberOfPeaks, sep = ""))
  print(paste("classOfClass: ", classOfClass, sep = ""))
  print(paste("thisClass: ", thisClass, sep = ""))
  ## repeated random subsampling validation
  print("### data")
  print(paste("minimumNumberOfPosSpectraPerClass: ", minimumNumberOfPosSpectraPerClass, sep = ""))
  print(paste("minimumNumberOfNegSpectraPerClass: ", minimumNumberOfNegSpectraPerClass, sep = ""))
  print(paste("numberOfDataSetDecompositions: ", numberOfDataSetDecompositions, sep = ""))
  print(paste("proportionTraining: ", proportionTraining, sep = ""))
  print(paste("fragmentColumnSelectionMethod: ", fragmentColumnSelectionMethod, sep = ""))
  print(paste("minimumProportionOfPositiveFragments: ", minimumProportionOfPositiveFragments, sep = ""))
  ## postprocessing
  print("### postprocessing")
  print(paste("unitResolution: ", unitResolution, sep = ""))
  ## cache
  print("### cache")
  print(paste("reParseMsp: ", reParseMsp, sep = ""))
  print(paste("reProcessAnno: ", reProcessAnno, sep = ""))
  print(paste("reProcessLibrary: ", reProcessLibrary, sep = ""))
  print(paste("reProcessFragmentMatrix: ", reProcessFragmentMatrix, sep = ""))
  ## constraints
  print("### constraints")
  print(paste("thisLibrary: ", thisLibrary, sep = ""))
  print(paste("thisMethod: ", thisMethod, sep = ""))
  print(paste("thisParameterSet: ", paste(names(thisParameterSet), thisParameterSet, sep = "=", collapse = "; "), sep = ""))
  print(paste("splitClasses: ", splitClasses, sep = ""))
  print(paste("splitClassesParameterSet: ", paste(names(splitClassesParameterSet), splitClassesParameterSet, sep = "=", collapse = "; "), sep = ""))
  print(paste("computeFragmentFisherTest: ", computeFragmentFisherTest, sep = ""))
  
  #######################################################################
  ## input and output
  if(!dir.exists(resultFolderForClassifiers))
    dir.create(resultFolderForClassifiers)
  if(!dir.exists(resultFolderForMetfamilyProjects))
    dir.create(resultFolderForMetfamilyProjects)
  
  #######################################################################
  ## constraints
  
  structureFormats <- list(InChI="InChI", InChIKey="InChIKey", SMILES="SMILES")
  #classOfClass <- paste("ChemOnt", ifelse(processSubstanceClasses, "SubstanceClass", "Substituent"), sep = "_")
  
  #######################################################################
  ## result dataFrame
  colNamesResults <- c(
    "Library", "Number of spectra", "Substance class", "Number of positive spectra", "Method", "Parameters", "AlgorithmName",
    "AUC", "AUC-PR", "TPR for FPR = 5%", "Sn for Max sum Sn+Sp", "Sp for Max sum Sn+Sp", "TNR for FNR = 5%", "time_t (s)", "time_p (s)", "AUCs", "AUC-PRs"
  )
  
  #resultDf <- as.data.frame(matrix(ncol = length(colNames)))
  #colnames(resultDf) <- colNamesResults
  
  #######################################################################
  ## analysis methods
  applicableCaretModelNames <- getApplicableCaretModels()
  
  if(all(!is.null(thisMethod), grepl(x = thisMethod, pattern = "^method=caretStacked"))){
    modelNames <- strsplit(x = thisMethod, split = "modelNames=")[[1]][[2]]
  }else{
    modelNames <- NA
  }
  
  methods <- list(
    ## scores
    "CosinusDistance"  = predict_CosinusDistance,
    "ColSums"          = colSums_classifier,
    "ColSumsPos"       = colSumsPos_classifier,
    "ColSumsPosOnly"   = colSumsPosOnly_classifier,
    "ColSumsRatio"     = colSumsRatio_classifier,
    "ColSumsLog"       = colSumsLog_classifier,
    "Prod"             = predict_Prod,
    "Jaccard"          = predict_Jaccard,
    "JaccardWeighted"  = predict_JaccardWeighted,
    ## classes
    "LDA"              = lda_classifier,
    "SPLSDA"           = splsda_classifier,
    "Correlation"      = predict_Correlation,
    #"RDA"              = predict_RDA,
    "SVM"              = svm_classifier,
    #"SOM"              = predict_SOM,
    #"XYF"              = predict_XYF,
    #"NeuralNet"        = predict_NeuralNet,
    "TensorFlow"        = tensorflow_classifier,
    "caret"             = caret_classifier,
    "caretStacked"      = caret_classifier_stacked
  )
  
  params <- list(
    ## scores
    "CosinusDistance"  = list(
      list(ratio=FALSE, smoothIntensities=FALSE),
      list(ratio=FALSE, smoothIntensities=TRUE),
      list(ratio=TRUE, smoothIntensities=FALSE),
      list(ratio=TRUE, smoothIntensities=TRUE)
    ),
    "ColSums"          = list(
      list(smoothIntensities=FALSE),
      list(smoothIntensities=TRUE)
    ),
    "ColSumsPos"          = list(
      list(smoothIntensities=FALSE),
      list(smoothIntensities=TRUE)
    ),
    "ColSumsPosOnly"          = list(
      list(smoothIntensities=FALSE),
      list(smoothIntensities=TRUE)
    ),
    "ColSumsRatio"          = list(
      list(smoothIntensities=FALSE),
      list(smoothIntensities=TRUE)
    ),
    "ColSumsLog"        = list(
      list(smoothIntensities=FALSE),
      list(smoothIntensities=TRUE)
    ),
    "Prod"             = list(
      list(ratio=FALSE, smoothIntensities=FALSE),
      list(ratio=FALSE, smoothIntensities=TRUE),
      list(ratio=TRUE, smoothIntensities=FALSE),
      list(ratio=TRUE, smoothIntensities=TRUE)
    ),
    "Jaccard"          = list(
      list(ratio=FALSE, smoothIntensities=FALSE),
      list(ratio=FALSE, smoothIntensities=TRUE),
      list(ratio=TRUE, smoothIntensities=FALSE),
      list(ratio=TRUE, smoothIntensities=TRUE)
    ),
    "JaccardWeighted"  = list(
      list(ratio=FALSE, smoothIntensities=FALSE),
      list(ratio=FALSE, smoothIntensities=TRUE),
      list(ratio=TRUE, smoothIntensities=FALSE),
      list(ratio=TRUE, smoothIntensities=TRUE)
    ),
    ## classes
    "LDA"              = list(
      list(smoothIntensities=FALSE),
      list(smoothIntensities=TRUE)
    ),
    "SPLSDA"              = list(
      list(smoothIntensities=FALSE),
      list(smoothIntensities=TRUE)
    ),
    "Correlation"      = list(
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[2]], ratio=FALSE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[2]], ratio=FALSE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[2]], ratio=TRUE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[2]], ratio=TRUE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[3]], ratio=FALSE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[3]], ratio=FALSE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[3]], ratio=TRUE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[1]], linkage = c("single", "average", "centroid")[[3]], ratio=TRUE, smoothIntensities=TRUE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=FALSE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=TRUE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=FALSE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=TRUE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=FALSE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=TRUE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=FALSE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=TRUE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=FALSE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=TRUE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=FALSE),
      #list(corMethod = c("pearson", "kendall", "spearman")[[2]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[1]], ratio=FALSE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[1]], ratio=TRUE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[2]], ratio=FALSE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[2]], ratio=FALSE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[2]], ratio=TRUE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[2]], ratio=TRUE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[3]], ratio=FALSE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[3]], ratio=FALSE, smoothIntensities=TRUE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[3]], ratio=TRUE, smoothIntensities=FALSE),
      list(corMethod = c("pearson", "kendall", "spearman")[[3]], linkage = c("single", "average", "centroid")[[3]], ratio=TRUE, smoothIntensities=TRUE)
    ),
    #"RDA"              = list(
    #  list()
    #),
    "SVM"              = list(
      #list(ratio=FALSE, smoothIntensities=FALSE),
      #list(ratio=TRUE,  smoothIntensities=FALSE),
      #list(ratio=FALSE, smoothIntensities=TRUE),
      #list(ratio=TRUE,  smoothIntensities=TRUE)
      list(smoothIntensities=FALSE),
      list(smoothIntensities=TRUE)
    ),
    #"SOM"              = list(
    #  list()
    #),
    #"XYF"              = list(
    #  list()
    #),
    #"NeuralNet"        = list(
    #  list()
    #),
    "TensorFlow"       = list(
      list(topology=list(numLayers = 0, numNeurons_l = NULL, dropout_l = NULL), smoothIntensities=FALSE),
      list(topology=list(numLayers = 1, numNeurons_l =   4, dropout_l = FALSE), smoothIntensities=FALSE),
      list(topology=list(numLayers = 1, numNeurons_l =  20, dropout_l = FALSE), smoothIntensities=FALSE),
      list(topology=list(numLayers = 1, numNeurons_l = 100, dropout_l = FALSE), smoothIntensities=FALSE),
      list(topology=list(numLayers = 1, numNeurons_l =   4, dropout_l = TRUE ), smoothIntensities=FALSE),
      list(topology=list(numLayers = 1, numNeurons_l =  20, dropout_l = TRUE ), smoothIntensities=FALSE),
      list(topology=list(numLayers = 1, numNeurons_l = 100, dropout_l = TRUE ), smoothIntensities=FALSE),
      list(topology=list(numLayers = 2, numNeurons_l = c(100, 20), dropout_l = FALSE), smoothIntensities=FALSE),
      list(topology=list(numLayers = 2, numNeurons_l = c( 20,  4), dropout_l = FALSE), smoothIntensities=FALSE),
      list(topology=list(numLayers = 2, numNeurons_l = c( 5,   4), dropout_l = FALSE), smoothIntensities=FALSE),
      list(topology=list(numLayers = 2, numNeurons_l = c(100, 20), dropout_l = TRUE ), smoothIntensities=FALSE),
      list(topology=list(numLayers = 2, numNeurons_l = c( 20,  4), dropout_l = TRUE ), smoothIntensities=FALSE),
      list(topology=list(numLayers = 2, numNeurons_l = c( 5,   4), dropout_l = TRUE ), smoothIntensities=FALSE)
    ),
    ## classWeights
    "caret"              = unlist(lapply(X = applicableCaretModelNames, FUN = function(model){
      apply(X = as.data.frame(expand.grid("smoothIntensities" = c(T,F), "classWeights" = c(T,F))), MARGIN = 1, FUN = function(params){
        list(smoothIntensities=params[["smoothIntensities"]], classWeights=params[["classWeights"]], modelName=model)
      })
      #lapply(X = c(TRUE, FALSE), FUN = function(si){
      #  list(smoothIntensities=si, modelName=model)
      #})
    }), recursive = FALSE),
    "caretStacked" = list(
      list(smoothIntensities=FALSE, modelNames=modelNames),
      list(smoothIntensities=TRUE,  modelNames=modelNames)
    )
  )
  
  algorithms <- list()
  for(methodIdx in seq_along(methods)){
  	print(paste0("Running classifier method ", methodIdx, "..."))
    methodName <- names(methods)[[methodIdx]]
    paramsLists <- params[[methodName]]
    for(paramsListIdx in seq_along(paramsLists)){
      paramsList <- paramsLists[[paramsListIdx]]
      paramsString <- paste(names(paramsList), paramsList[names(paramsList)], sep = "=", collapse = ", ")
      algoName <- paste("method=", methodName, "; ", paramsString, "", sep = "")
      
      algorithms[[length(algorithms) + 1]] <- list(
        "method" = methods[[methodIdx]],
        "methodName" = methodName,
        "params" = paramsList,
        "paramsString" = paramsString,
        "algoName" = algoName
      )
    }
  }
  #print(unlist(lapply(X = algorithms, FUN = function(x){x$algoName})))
  
  #######################################################################
  ## selections
  isMetFamilyProject <- FALSE
  if(!is.null(thisLibrary)){
    isMetFamilyProject <- FALSE #all(length(thisLibrary) == 1, file.exists(thisLibrary))
    
 #   if(isMetFamilyProject){
 #     libraryIdx   <- 1
 #     
 #     libraries     <- thisLibrary
 #     #annoFiles     <- NULL
 #     parameterSets <- list(highResolutionParameterSet)
 #     allowedInstrumentTypeSets <- NULL
 #     
 #     fileSpectra  <- thisLibrary
 #     annoFile     <- NULL
 #     parameterSet <- highResolutionParameterSet
 #     allowedInstrumentTypes <- NULL
 #   } else {
      ##########################################################################
      ## select library
      idx <- which(grepl(pattern = paste(thisLibrary, "$", sep = ""), x = libraries))
      libraries     <- libraries[idx]
      #annoFiles     <- annoFiles[idx]
      annoFiles     <- annoFile
      parameterSets <- parameterSets[idx]
      allowedInstrumentTypeSets <- allowedInstrumentTypeSets[idx]
      
      libraryIdx   <- 1
      fileSpectra  <- libraries[[1]]
      #annoFile     <- annoFiles[[1]]
      parameterSet <- parameterSets[[1]]
      allowedInstrumentTypes <- allowedInstrumentTypeSets[[1]]
      
      if(!is.null(thisParameterSet)){
        parameterSet <- thisParameterSet
        parameterSets[[1]] <- thisParameterSet
      }
#    }
  }
  #numberOfLibraries <- length(libraries)
  
  #if(isMetFamilyProject & writeMetFamilyProjects){
  #  writeMetFamilyProjects <- FALSE
  #  print("### Warning ### writeMetFamilyProjects not possible for input MetFamily project?")
  #}
  
  if(!is.null(thisMethod)){
    ##################################
    ## select method
    algoNames <- unlist(lapply(X = algorithms, FUN = function(x){x$algoName}))
    idx <- which(algoNames == thisMethod)
    algorithms <- algorithms[idx]
    methodIdx  <- 1
    algorithm <- algorithms[[1]]
    
    if(algorithm$methodName == "caret"){
      caretModelProperties <- caret::getModelInfo()[[algorithm$params$modelName]]
      if(algorithm$params$classWeights & !("Accepts Case Weights" %in% caretModelProperties$tags))
        stop("Attempt to use weights without being supported by the model")
    }
  }
  numberOfAlgorithms <- length(algorithms)
  
  #######################################################################
  ## run
  for(libraryIdx in seq_along(libraries)){
    if(isMetFamilyProject){
      fileSpectra  <- libraries[[libraryIdx]]
      annoFile     <- NULL
      parameterSet <- parameterSets[[libraryIdx]]
      allowedInstrumentTypes <- NULL
    } else {
      fileSpectra  <- libraries[[libraryIdx]]
      annoFile     <- annoFiles[[libraryIdx]]
      parameterSet <- parameterSets[[libraryIdx]]
      allowedInstrumentTypes <- allowedInstrumentTypeSets[[libraryIdx]]
    }
    ## annoFile <- "/home/htreutle/Downloads/MetSWATH/MONA/190523_MSMS_HR_someScaffolds.tsv"
    
    ## get library data
    libraryName <- basename(fileSpectra)
    
#    if(FALSE){
#      ## reparse msp
#      processSpectraAndAnnotation(fileSpectra, parameterSet, allowedInstrumentTypes, annoFile, structureFormats, TRUE, progress)
#      next
#    }
    
    cat(paste("Using library", libraryName, "\n"))
    
#=BREAKPOINT========================================================================================================  
#=BREAKPOINT========================================================================================================  
#=BREAKPOINT========================================================================================================  
    if(isMetFamilyProject){
      dataList <- readClusterDataFromProjectFile(file = parameterSetAll$thisLibrary, progress = NA)
      #dataList$annoArrayOfLists
      
      numberOfSpectra                                      <- dataList$numberOfPrecursors
      #spectraList                                          <- dataList$spectraList
      #annoTable                                            <- dataList$annoTable
      #if(exists("spectraList")) rm(spectraList)
      if(exists("annoTable"))   rm(annoTable)
      spectraList <- lapply(X = seq_len(dataList$numberOfPrecursors), FUN = function(idx){
        indeces <- which(dataList$featureMatrix[idx, ] != 0)
        indeces <- indeces[dataList$fragmentMasses[indeces] > 0]
        ms2Peaks_mzs  <- dataList$fragmentMasses[indeces]
        ms2Peaks_ints <- dataList$featureMatrix[idx, indeces]
        return(list(
          "name"           = dataList$dataFrameInfos$"Metabolite name"[[idx]],
          "ms1Int"         = -1,
          "rt"             = as.numeric(dataList$dataFrameInfos$"RT"[[idx]]),
          "mz"             = as.numeric(dataList$dataFrameInfos$"m/z"[[idx]]),
          "metName"        = dataList$dataFrameInfos$"Metabolite name"[[idx]],
          "adduct"         = dataList$dataFrameInfos$"Adduct ion name"[[idx]],
          "quantMass"      = -1,
          "compoundClass"  = dataList$dataFrameInfos$Annotation[[idx]],
          "instrumentType" = "Unknown",
          "inchi"          = NA,
          "inchiKey"       = NA,
          "smiles"         = NA,
          "peakNumber"     = length(ms2Peaks_mzs),
          "ms2Peaks_mz"    = ms2Peaks_mzs,
          "ms2Peaks_int"   = ms2Peaks_ints,
          "spectrumString" = paste(ms2Peaks_mzs, ms2Peaks_ints, sep = " ", collapse = ";")#,
          #"entryInterval"  = NULL
        ))
      })
      
      ## duplicated structures
      duplicatedStructures                                 <- rep(x = FALSE, times = numberOfSpectra)
      duplicatedSpectrumIndecesToRemove_maxPeakCount       <- vector(mode = "integer")
      duplicatedSpectrumIndecesToRemove_merge              <- vector(mode = "integer")
      numberOfSpectraUnique                                <- numberOfSpectra
      ## matrix
      matrixOriginal                                       <- dataList$featureMatrix
      fragmentMasses                                       <- dataList$fragmentMasses
      numberOfMS2PeakGroups                                <- ncol(matrixOriginal)
      spectraCount_fragment                                <- apply(X = matrixOriginal, MARGIN = 2, FUN = function(column){sum(column>0)})
      
      ## substance classes
      substanceclasses                                     <- vector(mode = "character", length = numberOfSpectra)
      allSubstanceClasses                                  <- dataList$annoArrayOfLists
      substanceclassesWithSuperClassesWithEnoughSpectra_sc <- vector(mode = "character", length = numberOfSpectra)
      substanceclassesWithSuperClassesWithEnoughSpectra_asc <- vector(mode = "character", length = numberOfSpectra)
      ## substituents
      substituents                                         <- vector(mode = "character", length = numberOfSpectra)
      substituentsWithEnoughSpectra_subst                  <- vector(mode = "character", length = numberOfSpectra)
      scaffolds                        = vector(mode = "character", length = numberOfSpectra)
      scaffoldsWithEnoughSpectra_subst = vector(mode = "character", length = numberOfSpectra)
      
      presentClasses <- sort(unique(unlist(allSubstanceClasses)))
      
      substanceclassesWithSuperClassesWithEnoughSpectra_sc <- NULL
      substanceclassesWithSuperClassesWithEnoughSpectra_asc <- presentClasses
      scaffoldsWithEnoughSpectra_subst <- NULL
      
      if(all(!is.null(thisClass), !(thisClass %in% presentClasses))) stop(paste("Class '", thisClass, "' not in the set of available classes:", paste(presentClasses, collapse = "; ")))
    } else {
      paramsHash <- digest::sha1(algo = "crc32", x = unlist(c(
        parameterSet=sort(unlist(parameterSet)), 
        allowedInstrumentTypes=sort(unlist(allowedInstrumentTypes)), 
        annoFile=basename(annoFile), 
        #annoFile=gsub(x = annoFile, pattern = "^/ifs/data/", replacement = "/vol/"), 
        structureFormats=sort(unlist(structureFormats)), 
        progress=progress, 
        minimumNumberOfPosSpectraPerClass=minimumNumberOfPosSpectraPerClass, 
        minimumNumberOfNegSpectraPerClass=minimumNumberOfNegSpectraPerClass#,
        #mergeDuplicatedSpectra=mergeDuplicatedSpectra,
        ##considerAlternativeSubstanceClasses=considerAlternativeSubstanceClasses,
        #unitResolution = unitResolution
      )))
      #paramsHash_old <- digest::sha1(algo = "crc32", x = unlist(c(
      #  parameterSet=sort(unlist(parameterSet)), 
      #  allowedInstrumentTypes=sort(unlist(allowedInstrumentTypes)), 
      #  annoFile=annoFile, 
      #  #annoFile=gsub(x = annoFile, pattern = "^/ifs/data/", replacement = "/vol/"), 
      #  structureFormats=sort(unlist(structureFormats)), 
      #  progress=progress, 
      #  minimumNumberOfPosSpectraPerClass=minimumNumberOfPosSpectraPerClass, 
      #  minimumNumberOfNegSpectraPerClass=minimumNumberOfNegSpectraPerClass#,
      #  #mergeDuplicatedSpectra=mergeDuplicatedSpectra,
      #  ##considerAlternativeSubstanceClasses=considerAlternativeSubstanceClasses,
      #  #unitResolution = unitResolution
      #)))
      #print(paste(paramsHash_old, ">", paramsHash))
      print(paste("Hash ", paramsHash, ": ", paste(sort(unlist(parameterSet)), collapse = "; "), "__", paste(sort(unlist(allowedInstrumentTypes)), collapse = "; "), "__", basename(annoFile), "__", paste(sort(unlist(structureFormats)), collapse = "; "), "__", progress, "__", minimumNumberOfPosSpectraPerClass, "__", minimumNumberOfNegSpectraPerClass, sep = ""))
      if(paramsHash == "44069439") {
        paramsHash <- "" ## backward compatibility
      } else {
        paramsHash <- paste(paramsHash, "_", sep = "")
      }
      
      fileProcessedLibraryRDataOld <- paste( substr(x = fileSpectra, start = 1, stop = regexpr(pattern = "\\.[a-zA-Z]{3,4}$", text = fileSpectra)[[1]] - 1), 
                                             "_",
                                             #"AltSC=", considerAlternativeSubstanceClasses, "_",
                                             "Class=", classOfClass, "_",
                                             "MergeSpectra=", mergeDuplicatedSpectra, "_",
                                             "unitResolution=", unitResolution, "_",
                                             "processed.RData", 
                                             sep = "")
      fileProcessedLibraryRData <- paste( substr(x = fileSpectra, start = 1, stop = regexpr(pattern = "\\.[a-zA-Z]{3,4}$", text = fileSpectra)[[1]] - 1), 
                                          "_",
                                          #"AltSC=", considerAlternativeSubstanceClasses, "_",
                                          #"Class=", classOfClass, "_",
                                          "MergeSpectra=", mergeDuplicatedSpectra, "_",
                                          "unitResolution=", unitResolution, "_",
                                          paramsHash,
                                          "processed.RData", 
                                          sep = "")
      #if(file.exists(fileProcessedLibraryRData2)){
      #  file.rename(from = fileProcessedLibraryRData2, to = fileProcessedLibraryRData)
      #}
      #next
      
      if(file.exists(fileProcessedLibraryRData) & reProcessLibrary){
        file.remove(fileProcessedLibraryRData)
      }
      
      if(file.exists(fileProcessedLibraryRData)){
        print(paste("load", fileProcessedLibraryRData))
        load(file = fileProcessedLibraryRData)
      } else {
        #if(file.exists())
        print(paste("create", fileProcessedLibraryRData))
        returnObj <- processLibrary(
          fileSpectra=fileSpectra, 
          parameterSet=parameterSet, 
          allowedInstrumentTypes=allowedInstrumentTypes, 
          annoFile=annoFile, 
          structureFormats=structureFormats, 
          progress=progress, 
          minimumNumberOfPosSpectraPerClass=minimumNumberOfPosSpectraPerClass, 
          minimumNumberOfNegSpectraPerClass=minimumNumberOfNegSpectraPerClass,
          mergeDuplicatedSpectra=mergeDuplicatedSpectra,
          #considerAlternativeSubstanceClasses=considerAlternativeSubstanceClasses,
          unitResolution = unitResolution,
          reParseMsp = reParseMsp,
          reProcessAnno = reProcessAnno,
          reProcessFragmentMatrix = reProcessFragmentMatrix
        )
        
        print(fileProcessedLibraryRData)
        save(file = fileProcessedLibraryRData, returnObj)
      }
      
      numberOfSpectra                                      <- returnObj$numberOfSpectra
      spectraList                                          <- returnObj$spectraList
      annoTable                                            <- returnObj$annoTable
      ## duplicated structures
      duplicatedStructures                                 <- returnObj$duplicatedStructures
      duplicatedSpectrumIndecesToRemove_maxPeakCount       <- returnObj$duplicatedSpectrumIndecesToRemove_maxPeakCount
      duplicatedSpectrumIndecesToRemove_merge              <- returnObj$duplicatedSpectrumIndecesToRemove_merge
      numberOfSpectraUnique                                <- returnObj$numberOfSpectraUnique
      ## matrix
      matrixOriginal                                       <- returnObj$fragmentMatrix
      fragmentMasses                                       <- returnObj$fragmentMasses
      numberOfMS2PeakGroups                                <- returnObj$numberOfMS2PeakGroups
      spectraCount_fragment                                <- returnObj$spectraCount_fragment
      ## substance classes
      substanceclasses                                     <- returnObj$substanceclasses
      allSubstanceClasses                                  <- returnObj$allSubstanceClasses
      substanceclassesWithSuperClassesWithEnoughSpectra_sc <- returnObj$substanceclassesWithSuperClassesWithEnoughSpectra_sc
      substanceclassesWithSuperClassesWithEnoughSpectra_asc <- returnObj$substanceclassesWithSuperClassesWithEnoughSpectra_asc
      ## substituents
      substituents                                         <- returnObj$substituents
      substituentsWithEnoughSpectra_subst                  <- returnObj$substituentsWithEnoughSpectra_subst
      scaffolds                        = returnObj$scaffolds
      scaffoldsWithEnoughSpectra_subst = returnObj$scaffoldsWithEnoughSpectra_subst
      rm(returnObj)
    }
    #=BREAKPOINT========================================================================================================  
    #=BREAKPOINT========================================================================================================  
    #=BREAKPOINT========================================================================================================  
    
    if(length(duplicatedSpectrumIndecesToRemove_maxPeakCount) != length(duplicatedSpectrumIndecesToRemove_merge)) stop("1:50 bug there")
    
    set.seed(1)
    spectraList <- lapply(X = spectraList, FUN = function(spectrum){
      if(is.na(spectrum$rt)){
        spectrum$rt <- runif(n = 1, min = 1000, max = 1000000)
      }
      return(spectrum)
    })
    precursorMzHere <- unlist(lapply(X = spectraList, FUN = function(x){x$mz}))
    precursorRtHere <- unlist(lapply(X = spectraList, FUN = function(x){x$rt}))
    #precursorRtHere[is.na(precursorRtHere)] <- -1
    #precursorRtHere <- seq_along(spectraListHere) ## no duplicates
    precursorLabels <- paste(precursorMzHere, precursorRtHere, sep = " / ")
    while(any(duplicated(precursorLabels))){
      precursorRtHere <- precursorRtHere + runif(n = length(precursorRtHere), min = 0.01, max = 0.1)
      precursorLabels <- paste(precursorMzHere, precursorRtHere, sep = " / ")
      for(idx in seq_along(spectraList))
        spectraList[[idx]]$rt <- precursorRtHere[[idx]]
    }
    
    
    if(mergeDuplicatedSpectra){
      duplicatedSpectrumIndecesToRemove <- duplicatedSpectrumIndecesToRemove_merge
    } else {
      duplicatedSpectrumIndecesToRemove <- duplicatedSpectrumIndecesToRemove_maxPeakCount
    }
    
    spectraToRetain <- rep(x = TRUE, times = numberOfSpectra)
    spectraToRetain[duplicatedSpectrumIndecesToRemove] <- FALSE
    
#    if(FALSE){
#      matrixMergedSpectra <- matrixOriginal[spectraToRetain, ]
#      rownames(matrixMergedSpectra) <- inchiKeysToInchiKeysBlock1(annoTable$InChIKey[spectraToRetain])
#      file <- gsub(x = fileSpectra, pattern = "\\.msp$", replacement = "_fragmentMatrix.tsv")
#      write.table(x = as.matrix(matrixMergedSpectra), file = file, sep = "\t")
#      #plot()
#    }
    
    ## load splitted classes
    paramsHash2 <- digest::sha1(algo = "crc32", x = unlist(splitClassesParameterSet))
    fileSpectraShort <- ifelse(test = isMetFamilyProject, yes = fileSpectra, no = substr(x = fileSpectra, start = 1, stop = regexpr(pattern = "\\.[a-zA-Z]{3,4}$", text = fileSpectra)[[1]] - 1))
    
    fileSplitClassRData <- paste(fileSpectraShort , "_", paramsHash2, "_splitClasses.RData", sep = "")
    if(splitClasses & file.exists(fileSplitClassRData)){
      print(paste("load", fileSplitClassRData))
      load(file = fileSplitClassRData)
    } else
      results__subSets    = list()
    
    ## classifier field
    #results__spectrum_class <- array(data = list(), dim = c(numberOfSpectraUnique))
    results__class      = list()
    #results__class2     = list()
    results__algo_class = array(data = list(), dim = c(numberOfAlgorithms))
    classifiers__algo_class = array(data = list(), dim = c(numberOfAlgorithms))
    for(methodIdx in seq_along(algorithms)){
      results__algo_class    [[methodIdx]] <- list()
      classifiers__algo_class[[methodIdx]] <- list()
    }
    #results__class3     = vector(mode = "numeric")
    
    classes <- NULL # vector of classes with enough spectra
    trueClasses <- NULL # list of vectors of true classes of the spectra
    
    switch(classOfClass,
           "ChemOnt|MainClass"={
             classes <- substanceclassesWithSuperClassesWithEnoughSpectra_sc
             trueClasses <- as.list(substanceclasses)
           },
           "ChemOnt|AllClasses"={
             classes <- substanceclassesWithSuperClassesWithEnoughSpectra_asc
             trueClasses <- allSubstanceClasses
           },
           "ChemOnt|Substituent"={
             classes <- substituentsWithEnoughSpectra_subst
             trueClasses <- substituents
           },
           "Scaffold"={
             classes <- scaffoldsWithEnoughSpectra_subst
             trueClasses <- scaffolds
           },
           stop(paste("Unknown class (", classOfClass, ")!", sep = ""))
    )
    ## selected class
    if(!is.null(thisClass)){
      idx <- which(grepl(pattern = paste(thisClass, "$", sep = ""), x = classes))
      if(length(idx) != 1){
        #stop(paste("Class not there:", thisClass))
        print("")
        print("##################################################################")
        print(paste("Class not there:", thisClass, "in", paste(classes, collapse = "; ")))
        print("##################################################################")
        print("")
        return()
      }
      classes <- classes[idx]
      #selectedClass <- substanceclassesWithSuperClassesWithEnoughSpectra_sc[[1]]
    }
    numberOfClasses <- length(classes)
    
    #######################################################################
    ## iterate classes
    rmdStatistics <- list()
    
    ## calc smiles stuff
#    if(calculateRMDstatistics){
#      # "obabel -i inchi -o smi AllInchis_new.inchi > AllInchis_new.smi"
#      tmpFileIn  <- "/home/htreutle/Downloads/tmp/inchisTmp.txt"
#      tmpFileOut <- "/home/htreutle/Downloads/tmp/smilesTmp.txt"
#      inchis <- sort(unique(unlist(lapply(X = spectraList, FUN = function(x){x$inchi}))))
#      writeLines(text = inchis, con = tmpFileIn)
#      #cmd <- paste("obabel -i inchi -o smi ", tmpFileIn, " > ", tmpFileOut, sep = "")
#      #output <- system(command = cmd, intern = TRUE)
#      smiles <- readLines(con = tmpFileOut)
#      
#      for(idx in seq_along(spectraList)){
#        #if(spectraList[[idx]]$smiles == "NA")
#          spectraList[[idx]]$smiles <- smiles[[which(inchis == spectraList[[idx]]$inchi)]]
#      }
#    }
    
    selectedSpectra_class <- lapply(X = classes, FUN = function(classHere){
      classRegExHere <- gsub(x = classHere,      pattern = "\\\\", replacement = "\\\\\\\\")
      classRegExHere <- gsub(x = classRegExHere, pattern = "\\[",  replacement = "\\\\[")
      classRegExHere <- gsub(x = classRegExHere, pattern = "\\]",  replacement = "\\\\]")
      
      unlist(lapply(X = trueClasses, FUN = function(x){any(stri_startswith_fixed(pattern = classRegExHere, str = x))}))
    })
    names(selectedSpectra_class) <- classes
    
    print(paste(length(classes), " classes with ", sum(unlist(selectedSpectra_class)), " members, e.g. ", paste(classes[1:5], collapse = "; "), sep = ""))
    
    #classIdx = which(classes == "InChI=1S/C7H14N2/c1-6-2-8-4-7(1)5-9-3-6/h6-9H,1-5H2")
    
    #classIdx <- 1
    #classIdx <- 865 Flav neg all
    #classIdx <- 821 Flav neg all
    
    #classIdx <- 243
    #classIdx <- 353 ## free neg Cinnamic acids and derivatives
    #classIdx <- 335 ## MONA neg Cinnamic acids and derivatives
    #classIdx <- 765 ## free MS/MS neg
    #classIdx <- 359 ## free MS/MS
    #classIdx <- 421 ## MoNA pos
    #classIdx <- 340 ## MoNA neg
    #classIdx <- 1022 ## NIST pos
    #classIdx <- 1004 ## NIST pos
    #classIdx <- 636 ## free pos
    
    #for(classIdx in 1:10){
    #for(classIdx in 1056:numberOfClasses){
    for(classIdx in seq_along(classes)){
      class <- classes[[classIdx]]
      
      #selectedSpectra <- unlist(lapply(X = trueClasses, FUN = function(x){any(stri_startswith_fixed(pattern = classRegEx, str = x))}))
      selectedSpectra <- selectedSpectra_class[[classIdx]]
      
      numberOfPositiveSpectra <- sum(selectedSpectra)
      numberOfNegativeSpectra <- sum(!selectedSpectra)
      #numberOfNegativeSpectra <- numberOfSpectra - numberOfPositiveSpectra
      
      numberOfPositiveSpectraUnique <- numberOfPositiveSpectra - sum(duplicatedStructures[ selectedSpectra])
      numberOfNegativeSpectraUnique <- numberOfNegativeSpectra - sum(duplicatedStructures[!selectedSpectra])
      
      cat(paste(libraryName, classIdx, " / ", numberOfClasses, class, " in ", numberOfPositiveSpectraUnique, " vs ", numberOfNegativeSpectraUnique, "\n"))
      
      if(numberOfNegativeSpectra != (numberOfSpectra - numberOfPositiveSpectra))  stop("oh oh")
      
      ## select spectra of class
      classes_pm <- vector(mode = "character", length = numberOfSpectra)
      classes_pm[ selectedSpectra] <- "+"
      classes_pm[!selectedSpectra] <- "-"
      
      #######################################################################
      ## sort out fragments in the negatives with too less spectra for efficiency reasons
      minimumNumberOfSpectraPerFragment <- 3
      if(removeRareFragments){
        if(as.numeric(numberOfMS2PeakGroups) * as.numeric(sum(selectedSpectra)) < 1E6){
          positiveFragments <- apply(X = matrixOriginal[selectedSpectra, ], MARGIN = 2, FUN = function(x){sum(x!=0)}) > 0
        } else {
          positiveFragments <- Matrix::colSums(x = matrixOriginal[selectedSpectra, ] != 0) > 0
        }
        fragmentsToRemove <- (spectraCount_fragment <= minimumNumberOfSpectraPerFragment) & (!positiveFragments)
      } else {
        fragmentsToRemove <- logical(length = numberOfMS2PeakGroups)
      }
      
      #print(paste(numberOfMS2PeakGroups-sum(fragmentsToRemove), "/", numberOfMS2PeakGroups, "fragment groups retained:"))
      matrix <- matrixOriginal[, !fragmentsToRemove]
      rownames(matrix) <- paste(seq_len(nrow(matrix)), rownames(matrix), sep = "__")
      fragmentMasses <- fragmentMasses[!fragmentsToRemove]
      numberOfMS2PeakGroups <- length(fragmentMasses)
      spectraCount_fragment <- spectraCount_fragment[!fragmentsToRemove]
      
      ## XXX fix: remove duplicated fragment columns
      duplicatedColNames <- duplicated(colnames(matrix))
      for(colName in unique(colnames(matrix)[duplicatedColNames])){
        colIndeces <- which(colnames(matrix)==colName)
        rowSums <- apply(X = as.matrix(matrix[, colIndeces]), MARGIN = 1, FUN = sum)
        matrix[, colIndeces[[1]]] <- rowSums
        matrix <- matrix[, -colIndeces[-1]]
      }
      fragmentMasses <- fragmentMasses[!duplicatedColNames]
      numberOfMS2PeakGroups <- length(fragmentMasses)
      spectraCount_fragment <- spectraCount_fragment[!duplicatedColNames]
      
      ## complete matrix and classes for complete classifier
      matrix_train_all     <- matrix
      classes_pm_train_all <- classes_pm
      if(length(duplicatedSpectrumIndecesToRemove) > 0){
        matrix_train_all     <- matrix_train_all    [-duplicatedSpectrumIndecesToRemove, ]
        classes_pm_train_all <- classes_pm_train_all[-duplicatedSpectrumIndecesToRemove]
      }
      spectraToRetainPos <- spectraToRetain &  selectedSpectra
      spectraToRetainNeg <- spectraToRetain & !selectedSpectra
      
      #######################################################################
      ## writeMetFamilyProjects
      if(writeMetFamilyProjects | (splitClasses & is.null(results__subSets[[class]]))){
        ## idxHere <- which(unlist(lapply(spectraListHere, function(x){x$name})) == "Cyanidin-3-O-rhamnoside cation")
        spectraListHere <- spectraList[spectraToRetain]
        spectraListHere <- spectraListHere[classes_pm_train_all=="+"]
        
        libraryName2 <- gsub(x = libraryName, pattern = "^\\d\\d\\d\\d\\-\\d\\d\\-\\d\\d_\\d\\d:\\d\\d:\\d\\d_", replacement = "")
        resultFolderForMetfamilyProjectsHere <- paste(resultFolderForMetfamilyProjects, "/MetFamily_", libraryName2, sep = "")
        
        file <- createMetFamilyClassFile(class=class, numberOfPrecursors=length(spectraListHere), resultFolder=resultFolderForMetfamilyProjectsHere)
        print(paste(file.exists(file), file))
        if(!file.exists(file)){
          if(!dir.exists(resultFolderForMetfamilyProjectsHere))
            dir.create(path = resultFolderForMetfamilyProjectsHere)
          writeMetFamilyProject(spectraListHere = spectraListHere, classHere = tail(x = strsplit(x = class, split = "; ")[[1]], n=1), parameterSet = parameterSet, file = file)
          print(paste("Written MetFamily project", file))
        }
      }## writeMetFamilyProjects
      
      #######################################################################
      ## split classes
      if(splitClasses & is.null(results__subSets[[class]]) & length(spectraList[spectraToRetainPos]) <= splitClassesParameterSet$maximumNumberOfPrecursors){# && writeMetFamilyProjects){
        spectraListHereMF <- spectraList[spectraToRetain]
        spectraListHereMF <- spectraListHereMF[classes_pm_train_all=="+"]
        
        spectraListHere <- spectraList[spectraToRetainPos]
        file <- createMetFamilyClassFile(class=class, numberOfPrecursors=length(spectraListHere), resultFolder=resultFolderForMetfamilyProjectsHere)
        
        #source("/home/htreutle/Code/Java/MetFamily/R/R_packages.R")
        
        #library("mzR")
        #library("xcms")
#        library("matrixStats")
#        library("Matrix")
#        library("tools")
#        library("stringr")
#        library("slam")
        
#        source("/home/htreutle/Code/Java/MetFamily/R/DataProcessing.R")
#        source("/home/htreutle/Code/Java/MetFamily/R/Analysis.R")
#        source("/home/htreutle/Code/Java/MetFamily/R/Plots.R")
#        source("/home/htreutle/Code/Java/MetFamily/R/TreeAlgorithms.R")
        
        dataList <- readClusterDataFromProjectFile(file = file, progress = NA)
        distanceMatrix <- calculateDistanceMatrix(dataList = dataList, filter = seq_len(dataList$numberOfPrecursors), distanceMeasure = splitClassesParameterSet$distanceMeasure, progress = NA)
        clusterDataList <- calculateCluster(dataList = dataList, filterObj = list(filter=distanceMatrix$filter), distanceMatrix = distanceMatrix$distanceMatrix, method = splitClassesParameterSet$clusterMethod, distanceMeasure = splitClassesParameterSet$distanceMeasure, progress = NA)
        
        rootIndex <- length(clusterDataList$cluster$height)
        df <- analyzeTreeForFrequentFragments(clusterDataList = clusterDataList, nodeIdx = rootIndex, minimumNumberOfChildren = splitClassesParameterSet$minimumNumberOfChildren)
        df <- df[df$NodeIdx > 0,]
        
        ## XXX remove
        proportionOfCompounds <- sum(df$numberOfChildren) / dataList$numberOfPrecursors
        print(paste(dataList$numberOfPrecursors, nrow(df), proportionOfCompounds))
        
        #inchikeysHere <- unlist(lapply(X = spectraListHere, FUN = function(x){x$inchiKey}))
        labelsHere <- unlist(lapply(X = spectraListHere, FUN = function(x){paste(x$mz, x$rt, sep = " / ")}))
        #labelsPresent <- gsub(x = gsub(x = dataList$precursorLabels, pattern = " ", replacement = ""), pattern = "/", replacement = " / ")
        labelsPresent <- paste(dataList$dataFrameInfos$`m/z`, dataList$dataFrameInfos$RT, sep = " / ")
        precursorIndecesHere <- lapply(X = df$NodeIdx, FUN = function(nodeIdx){getPrecursorSetFromTreeSelection(clusterDataList = clusterDataList, clusterLabel = nodeIdx)})
        spectraToRetainPosSubSets <- lapply(X = precursorIndecesHere, FUN = function(precursorIndeces){
          spectraToRetainPosSubSetHere <- vector(mode = "logical", length = numberOfPositiveSpectraUnique)
          spectraToRetainPosSubSetHere[precursorIndeces] <- TRUE
          #spectraToRetainPosSubSetHere <- spectraToRetainPosSubSetHere[match(x = inchikeysHere, table = dataList$dataFrameInfos$InChIKey)]
          spectraToRetainPosSubSetHere <- spectraToRetainPosSubSetHere[match(x = labelsHere, table = labelsPresent)]
          
          spectraToRetainPosSubSet <- spectraToRetainPos
          spectraToRetainPosSubSet[spectraToRetainPosSubSet] <- spectraToRetainPosSubSetHere
          return(spectraToRetainPosSubSet)
        })
        
        if(length(spectraToRetainPosSubSets) > 0){
          ids <- sapply(X = as.character(seq_along(spectraToRetainPosSubSets)), FUN = function(id){paste(paste(rep(x = "0", times = (floor(log10(length(spectraToRetainPosSubSets))) + 1 - nchar(id))), collapse = ""), id, sep = "")})
          names(spectraToRetainPosSubSets) <- paste(class, "; ", "SubSet_", ids, sep = "")
        }
        results__subSets[[class]] <- spectraToRetainPosSubSets
      }
      
      ## TODO 999
      #next
      
      #plot(sort(unlist(lapply(results__subSets, length))))
      
      if(!builtClassifiers){
        selectedSpectraWithSubSets <- list()
      } else {
        if(length(results__subSets[[class]]) > 0){
          selectedSpectraWithSubSets <- c(setNames(object = list(selectedSpectra), nm = class), results__subSets[[class]])
        } else {
          selectedSpectraWithSubSets <- list(class = selectedSpectra)
          names(selectedSpectraWithSubSets) <- class
        }
      }
      
      #######################################################################
      ## relative mass defect statistic
      if(calculateRMDstatistics){
        print(paste(classIdx, "/", length(classes), class))
        spectraListHere <- spectraList[-duplicatedSpectrumIndecesToRemove]
        spectraListHere <- spectraListHere[classes_pm_train_all=="+"]
        
        inchisHere <- unlist(lapply(X = spectraListHere, FUN = function(x){x$inchi}))
        smilesHere <- unlist(lapply(X = spectraListHere, FUN = function(x){x$smiles}))
        
        #sumFormulasHere <- unlist(lapply(X = strsplit(x = inchisHere, split = "/"), FUN = function(x){x[[2]]}))
        #smilesHere <- sapply(X = inchisHere, FUN = webchem::cs_inchi_smiles)
        
        sp <- get.smiles.parser()
        molecule <- parse.smiles(smilesHere)
        
        nothing <- sapply(X = molecule, FUN = convert.implicit.to.explicit)
        nothing <- sapply(X = molecule, FUN = do.aromaticity)
        nothing <- sapply(X = molecule, FUN = do.typing)
        nothing <- sapply(X = molecule, FUN = do.isotopes)
        
        formula <- sapply(X = molecule, FUN = get.mol2formula, charge=0)
        #formula <- sapply(X = formula, FUN = function(x){x@string})
        these <- sapply(X = formula, FUN = function(x){
          all(x@isotopes[, "isoto"] %in% c("C","H","N","O","P","S"))
        })
        
        molecule <- molecule[these]
        
        if(length(molecule) > 0){
          #exactMass <- sapply(X = formula, FUN = function(x){x@mass})
          exactMass <- sapply(X = molecule, FUN = get.exact.mass)
          
          nominalMass <- sapply(X = molecule, FUN = function(y){sum(sapply(X = sapply(X = get.atoms(object = y), FUN = get.atomic.number), FUN = function(x){
            atomNominalMass <- NA
            switch(as.character(x),
                   "6"={# C
                     atomNominalMass <- 12
                   },
                   "1"={# H
                     atomNominalMass <- 1
                   },
                   "7"={# N
                     atomNominalMass <- 14
                   },
                   "8"={# O
                     atomNominalMass <- 16
                   },
                   "15"={# P
                     atomNominalMass <- 31
                   },
                   "16"={# S
                     atomNominalMass <- 32
                   },
                   "35"={# Br
                     atomNominalMass <- 79
                   },
                   "17"={# Cl
                     atomNominalMass <- 35
                   },
                   "9"={# F
                     atomNominalMass <- 19
                   },
                   "14"={# Si
                     atomNominalMass <- 28
                   },
                   "12"={# Mg
                     atomNominalMass <- 24
                   },
                   "53"={# I
                     atomNominalMass <- 191
                   },
                   "11"={# Na
                     atomNominalMass <- 23
                   },
                   "5"={# B
                     atomNominalMass <- 10
                   },
                   {## unknown state
                     stop(paste("Error: unknown atom number", x))
                   }
            )## end switch
            return(atomNominalMass)
          }))})
          
          #precursorMzHere <- unlist(lapply(X = spectraListHere, FUN = function(x){x$mz}))
          #relativeMassDefect <- (precursorMzHere - round(precursorMzHere)) / precursorMzHere * 1E6
          relativeMassDefect <- (exactMass - nominalMass) / exactMass * 1E6
          
          rmdStatistics[[tail(n=1, x = strsplit(x = class, split = "; ")[[1]])]] <- relativeMassDefect
        }
      }
      
      for(subSetIdx in seq_along(selectedSpectraWithSubSets)){
        class <- names(selectedSpectraWithSubSets)[[subSetIdx]]
        selectedSpectra <- selectedSpectraWithSubSets[[subSetIdx]]
        
        numberOfPositiveSpectra <- sum(selectedSpectra)
        numberOfNegativeSpectra <- sum(!selectedSpectra)
        #numberOfNegativeSpectra <- numberOfSpectra - numberOfPositiveSpectra
        
        numberOfPositiveSpectraUnique <- numberOfPositiveSpectra - sum(duplicatedStructures[ selectedSpectra])
        numberOfNegativeSpectraUnique <- numberOfNegativeSpectra - sum(duplicatedStructures[!selectedSpectra])
        
        cat(paste(libraryName, classIdx, " / ", numberOfClasses, class, " in ", numberOfPositiveSpectraUnique, " vs ", numberOfNegativeSpectraUnique, "\n"))
        
        if(numberOfNegativeSpectra != (numberOfSpectra - numberOfPositiveSpectra))  stop("oh oh")
        
        ## select spectra of class
        classes_pm <- vector(mode = "character", length = numberOfSpectra)
        classes_pm[ selectedSpectra] <- "+"
        classes_pm[!selectedSpectra] <- "-"
        
        #######################################################################
        ## sort out fragments in the negatives with too less spectra for efficiency reasons
        minimumNumberOfSpectraPerFragment <- 3
        if(removeRareFragments){
          if(as.numeric(numberOfMS2PeakGroups) * as.numeric(sum(selectedSpectra)) < 1E6){
            positiveFragments <- apply(X = matrixOriginal[selectedSpectra, ], MARGIN = 2, FUN = function(x){sum(x!=0)}) > 0
          } else {
            positiveFragments <- Matrix::colSums(x = matrixOriginal[selectedSpectra, ] != 0) > 0
          }
          fragmentsToRemove <- (spectraCount_fragment <= minimumNumberOfSpectraPerFragment) & (!positiveFragments)
        } else {
          fragmentsToRemove <- logical(length = numberOfMS2PeakGroups)
        }
        
        #print(paste(numberOfMS2PeakGroups-sum(fragmentsToRemove), "/", numberOfMS2PeakGroups, "fragment groups retained:"))
        matrix <- matrixOriginal[, !fragmentsToRemove]
        rownames(matrix) <- paste(seq_len(nrow(matrix)), rownames(matrix), sep = "__")
        fragmentMasses <- fragmentMasses[!fragmentsToRemove]
        numberOfMS2PeakGroups <- length(fragmentMasses)
        spectraCount_fragment <- spectraCount_fragment[!fragmentsToRemove]
        
        ## XXX fix: remove duplicated fragment columns
        duplicatedColNames <- duplicated(colnames(matrix))
        for(colName in unique(colnames(matrix)[duplicatedColNames])){
          colIndeces <- which(colnames(matrix)==colName)
          rowSums <- apply(X = as.matrix(matrix[, colIndeces]), MARGIN = 1, FUN = sum)
          matrix[, colIndeces[[1]]] <- rowSums
          matrix <- matrix[, -colIndeces[-1]]
        }
        fragmentMasses <- fragmentMasses[!duplicatedColNames]
        numberOfMS2PeakGroups <- length(fragmentMasses)
        spectraCount_fragment <- spectraCount_fragment[!duplicatedColNames]
        
        ## complete matrix and classes for complete classifier
        matrix_train_all     <- matrix
        classes_pm_train_all <- classes_pm
        if(length(duplicatedSpectrumIndecesToRemove) > 0){
          matrix_train_all     <- matrix_train_all    [-duplicatedSpectrumIndecesToRemove, ]
          classes_pm_train_all <- classes_pm_train_all[-duplicatedSpectrumIndecesToRemove]
        }
        spectraToRetainPos <- spectraToRetain &  selectedSpectra
        spectraToRetainNeg <- spectraToRetain & !selectedSpectra
        
        
        
      #if(builtClassifiers){
        #######################################################################
        ## analyze
        #positiveSpectraIndeces <- seq_len(numberOfSpectra)[ selectedSpectra]
        #negativeSpectraIndeces <- seq_len(numberOfSpectra)[!selectedSpectra]
        
        ###########################################
        ## built classifier(s)
        
        ## repeated random subsampling validation
        set.seed(1)
        spectrumId_train_i <- list()
        spectrumId_test_i  <- list()
        for(dataSetIdx in seq_len(numberOfDataSetDecompositions)){
          #print(paste("Data set", dataSetIdx, "/", n))
          
          ## remove duplicated compounds
          if(!takeSpectraWithMaximumNumberOfPeaks & !mergeDuplicatedSpectra){
            inchiKeysBlock1    <- unlist(lapply(X = strsplit(x = annoTable[, "InChIKey"], split = "\\-"), FUN = function(x){x[[1]]}))
            duplicatedStructures <- duplicated(x = inchiKeysBlock1)
            duplicatedInchiKeysBlock1 <- unique(inchiKeysBlock1[duplicatedStructures])
            
            duplicatedSpectrumIndecesToRemove <- unlist(lapply(X = duplicatedInchiKeysBlock1, FUN = function(x){
              indeces <- which(inchiKeysBlock1 %in% x)
              indexToRetain <- sample(x = length(indeces), size = 1)
              indecesOut <- indeces[-indexToRetain]
              return(indecesOut)
            }))
          } else {
            ## duplicatedSpectrumIndecesToRemove is set fix
          }
          
          spectrumId_here <- seq_len(numberOfSpectra)
          selectedSpectraHere <- selectedSpectra
          if(length(duplicatedSpectrumIndecesToRemove) > 0){
            spectrumId_here     <- spectrumId_here    [-duplicatedSpectrumIndecesToRemove]
            selectedSpectraHere <- selectedSpectraHere[-duplicatedSpectrumIndecesToRemove]
          }
          
          numberOfSpectraHere <- length(selectedSpectraHere)
          numberOfPositiveSpectraHere <- sum(selectedSpectraHere)
          numberOfNegativeSpectraHere <- numberOfSpectraHere - numberOfPositiveSpectraHere
          positiveSpectraIndecesHere <- spectrumId_here[ selectedSpectraHere]
          negativeSpectraIndecesHere <- spectrumId_here[!selectedSpectraHere]
          
          ## sets
          numberOfPositives_train <- round(numberOfPositiveSpectraHere * proportionTraining)
          #numberOfPositives_test  <- numberOfPositiveSpectraHere - numberOfPositives_train
          positives_train <- sample(seq_len(numberOfPositiveSpectraHere), size = numberOfPositives_train, replace = FALSE)
          positives_test  <- setdiff(seq_len(numberOfPositiveSpectraHere), positives_train)
          
          numberOfNegatives_train <- round(numberOfNegativeSpectraHere * proportionTraining)
          #numberOfNegatives_test  <- numberOfNegativeSpectraHere - numberOfNegatives_train
          negatives_train <- sample(seq_len(numberOfNegativeSpectraHere), size = numberOfNegatives_train, replace = FALSE)
          negatives_test  <- setdiff(seq_len(numberOfNegativeSpectraHere), negatives_train)
          
          #numberOfSpectra_train <- numberOfPositives_train + numberOfNegatives_train
          #numberOfSpectra_test  <- numberOfPositives_test  + numberOfNegatives_test
          
          spectrumId_train  <- c(positiveSpectraIndecesHere[positives_train], negativeSpectraIndecesHere[negatives_train])
          spectrumId_test   <- c(positiveSpectraIndecesHere[positives_test ], negativeSpectraIndecesHere[negatives_test ])
          
          ## box
          spectrumId_train_i[[dataSetIdx]] <- spectrumId_train
          spectrumId_test_i [[dataSetIdx]] <- spectrumId_test
          
          if(length(spectrumId_train) == 0 | length(spectrumId_test) == 0)
            stop("oh")
        }## dataset
        
        #################################################################################################################
        #################################################################################################################
        for(methodIdx in seq_along(algorithms)){
          algorithm <- algorithms[[methodIdx]]
          methodFunctions   <- algorithm$method
          methodName        <- algorithm$methodName
          parameters_train  <- algorithm$params
          paramsString      <- algorithm$paramsString
          algoName          <- algorithm$algoName
          if(outDetail)
            cat(paste("Using algo", algoName, "\n"))
          
          smoothIntensities <- parameters_train$smoothIntensities
          parameters_train$smoothIntensities <- NULL
          
          #classifierName <- gsub(x = paste("library=", libraryName, "; Class=", classOfClass, "; AltSC=", considerAlternativeSubstanceClasses, "; ", algoName, sep = ""), pattern = "; ", replacement = "_")
          classifierName <- gsub(x = paste("library=", libraryName, "; Class=", classOfClass, "; ", algoName, sep = ""), pattern = "; ", replacement = "_")
          classifier <- list()
          classifier$classifierName          <- classifierName
          classifier$numberOfSpectra         <- numberOfSpectraUnique
          classifier$numberOfPositiveSpectra <- numberOfPositiveSpectraUnique
          classifier$numberOfNegativeSpectra <- numberOfNegativeSpectraUnique
          classifier$algorithm               <- algorithm
          classifier$class                   <- class
          classifier$fragmentMasses          <- fragmentMasses
          classifier$classOfClass            <- classOfClass
          
          ## I/O
          classifier$maximumNumberOfScores   <- maximumNumberOfScores
          classifier$removeRareFragments     <- removeRareFragments
          ## replicates
          classifier$mergeDuplicatedSpectra              <- mergeDuplicatedSpectra
          classifier$takeSpectraWithMaximumNumberOfPeaks <- takeSpectraWithMaximumNumberOfPeaks
          ## classes
          #classifier$considerAlternativeSubstanceClasses <- considerAlternativeSubstanceClasses
          #classifier$processSubstanceClasses <- processSubstanceClasses
          #classifier$classOfClass <- classOfClass
          ## repeated random subsampling validation
          classifier$minimumNumberOfPosSpectraPerClass <- minimumNumberOfPosSpectraPerClass
          classifier$minimumNumberOfNegSpectraPerClass <- minimumNumberOfNegSpectraPerClass
          classifier$numberOfDataSetDecompositions     <- numberOfDataSetDecompositions
          classifier$proportionTraining      <- proportionTraining
          ## postprocessing
          classifier$unitResolution <- unitResolution
          
          ## "AUC", "PPV for FDR = 5%", "Sn for Max sum Sn+Sp", "Sp for Max sum Sn+Sp", "Sn for FNR = 5%", "time (s)"
          auc_roc_j <- list()
          auc_pr_j <- list()
          tpr_j <- list()
          max_Sn_j <- list()
          max_Sp_j <- list()
          tnr_j <- list()
          time_t_j <- list()
          time_p_j <- list()
          predictedScores_pos_j <- list()
          predictedScores_neg_j <- list()
          
          for(dataSetIdx in seq_len(numberOfDataSetDecompositions)){
            spectrumId_train <- spectrumId_train_i[[dataSetIdx]]
            spectrumId_test  <- spectrumId_test_i [[dataSetIdx]]
            
            parameters_train$matrix_train     <- matrix    [spectrumId_train, ]
            parameters_train$classes_pm_train <- classes_pm[spectrumId_train]
            
            parameters_test <- list()
            parameters_test$matrix_test      <- matrix    [spectrumId_test , ]
            classes_pm_test                  <- classes_pm[spectrumId_test ]
            
            if(!smoothIntensities){
              #parameters$matrix_train [parameters$matrix_train  != 0] <- 1
              #parameters$matrix_test  [parameters$matrix_test   != 0] <- 1
              parameters_train$matrix_train@x[parameters_train$matrix_train@x != 0] <- 1
              parameters_test$matrix_test@x  [parameters_test$matrix_test@x   != 0] <- 1
            }
            if(grepl(x = algoName, pattern = "^method=((LDA)|(SPLSDA)|(caret)).*$")){
              ## shrink the columns for the classifier to the relevant one in the positive set
              selectedColumns <- selectColumnsForClassification(parameters_train$matrix_train, parameters_train$classes_pm_train, fragmentColumnSelectionMethod, minimumProportionOfPositiveFragments)
              numberOfSelectedColumns <- sum(selectedColumns)
              if(numberOfSelectedColumns == 0) stop("Number of selected columns is zero in the positive set")
              
              parameters_train$matrix_train <- parameters_train$matrix_train[, selectedColumns, drop=FALSE]
              parameters_test$matrix_test   <- parameters_test$matrix_test  [, selectedColumns, drop=FALSE]
            }
            
            #matrix_train     <- parameters_train$matrix_train
            #classes_pm_train <- parameters_train$classes_pm_train
            #matrix_test      <- parameters_test $matrix_test
            #classes_pm_test  <- parameters_test $classes_pm_test
            
            if(outDetail)
              cat(paste(libraryName, classIdx, " / ", numberOfClasses, class, " in ", numberOfPositiveSpectraUnique, " vs ", numberOfNegativeSpectraUnique, methodName, paste("(", methodIdx, " / ", numberOfAlgorithms, ")", sep = ""), "smooth", smoothIntensities))
            
            ########################## do
            startTime_t <- Sys.time()
            classifierHere <- tryCatch(
              {
                do.call(what = methodFunctions$train,    args = parameters_train)
              }, error = function(e) {
                print("")
                print(paste("#################### train"))
                print(paste("####################", e))
                print(paste("####################"))
                stop(e)
              }
            )
            parameters_test$classifier <- classifierHere
            endTime_t <- Sys.time()
            startTime_p <- Sys.time()
            predicted_scores <- tryCatch(
              {
                do.call(what = methodFunctions$classify, args = parameters_test)
              }, error = function(e) {
                print("")
                print(paste("#################### prediction"))
                print(paste("####################", e))
                print(paste("####################"))
                stop(e)
              }
            )
            endTime_p <- Sys.time()
            time_t_here <- difftime(endTime_t, startTime_t, units = "secs")[[1]]
            time_p_here <- difftime(endTime_p, startTime_p, units = "secs")[[1]]
            
            ########################## performance
            if(length(predicted_scores) != length(classes_pm_test) | length(predicted_scores) != length(spectrumId_test)) stop("Number of predicted scores inconsistent")
            returnObj   <- evaluatePerformance(predicted_scores = predicted_scores, classes_pm_test = classes_pm_test, spectrumId_test=spectrumId_test)
            
            auc_roc_here <- returnObj$auc_roc
            auc_pr_here  <- returnObj$auc_pr
            tpr_here     <- returnObj$tpr
            max_Sn_here  <- returnObj$max_Sn
            max_Sp_here  <- returnObj$max_Sp
            tnr_here     <- returnObj$tnr
            predictedScores_pos_here <- returnObj$predictedScores_pos
            predictedScores_neg_here <- returnObj$predictedScores_neg
            spectrumId_test_neg_here <- returnObj$spectrumId_test_neg
            spectrumId_test_pos_here <- returnObj$spectrumId_test_pos
            names(predictedScores_neg_here) <- spectrumId_test_neg_here
            names(predictedScores_pos_here) <- spectrumId_test_pos_here
            
            ########################## out
            lineOut <- paste("\tAUC-ROC", auc_roc_here, "\tAUC-PR", auc_pr_here, "\tTPR", tpr_here, "\tTNR", tnr_here, "\tmaxSn", max_Sn_here, "\tmaxSp", max_Sp_here, "\ttime_t", time_t_here, "\ttime_p", time_p_here, "\n")
            if(outDetail)
              cat(lineOut)
            
            auc_roc_j[[dataSetIdx]] <- auc_roc_here
            auc_pr_j [[dataSetIdx]] <- auc_pr_here
            tpr_j    [[dataSetIdx]] <- tpr_here
            max_Sn_j [[dataSetIdx]] <- max_Sn_here
            max_Sp_j [[dataSetIdx]] <- max_Sp_here
            tnr_j    [[dataSetIdx]] <- tnr_here
            time_t_j [[dataSetIdx]] <- time_t_here
            time_p_j [[dataSetIdx]] <- time_p_here
            predictedScores_pos_j[[dataSetIdx]] <- predictedScores_pos_here
            predictedScores_neg_j[[dataSetIdx]] <- predictedScores_neg_here
          }## dataset decompositions
          
          auc_roc <- median(unlist(as.numeric(auc_roc_j)))
          auc_pr  <- median(unlist(as.numeric(auc_pr_j )))
          tpr     <- median(unlist(as.numeric(tpr_j    )))
          max_Sn  <- median(unlist(as.numeric(max_Sn_j )))
          max_Sp  <- median(unlist(as.numeric(max_Sp_j )))
          tnr     <- median(unlist(as.numeric(tnr_j    )))
          time_t  <- median(unlist(as.numeric(time_t_j )))
          time_p  <- median(unlist(as.numeric(time_p_j )))
          
          ## aggregate validation runs
          resultRow <- c(
            "Library"=libraryName, "Number of spectra"=numberOfSpectraHere, "Substance class"=class, "Number of positive spectra"=numberOfPositiveSpectraUnique, 
            "Method"=methodName, "Parameters"=paramsString, "AlgorithmName"=algoName, 
            "AUC"=auc_roc, "AUC-PR"=auc_pr, "TPR for FPR = 5%"=tpr, "Sn for Max sum Sn+Sp"=max_Sn, "Sp for Max sum Sn+Sp"=max_Sp, "TNR for FNR = 5%"=tnr, 
            "time_t (s)"=time_t, "time_p (s)"=time_p, 
            "AUCs"=paste(as.numeric(unlist(auc_roc_j)), collapse = ", "), "AUC-PRs"=paste(as.numeric(unlist(auc_pr_j)), collapse = ", ")
          )
          
          if(writeClassifiers){
            ####################################################
            ## performance
            classifier$AUC    <- auc_roc
            classifier$AUC_PR <- auc_pr
            classifier$TPR_for_FPR_of_5Percent <- tpr
            classifier$TNR_for_FNR_of_5Percent <- tnr
            
            ####################################################
            ## complete classifier
            parameters_train_all <- algorithm$params
            smoothIntensities <- parameters_train_all$smoothIntensities
            parameters_train_all$smoothIntensities <- NULL
            
            parameters_train_all$matrix_train     <- matrix_train_all
            parameters_train_all$classes_pm_train <- classes_pm_train_all
            
            if(!smoothIntensities)
              parameters_train_all$matrix_train@x[parameters_train_all$matrix_train@x != 0] <- 1
            
            if(grepl(x = algoName, pattern = "^method=((LDA)|(SPLSDA)|(caret)).*$")){
              selectedColumns <- selectColumnsForClassification(parameters_train$matrix_train, parameters_train$classes_pm_train, fragmentColumnSelectionMethod, minimumProportionOfPositiveFragments)
              numberOfSelectedColumns <- sum(selectedColumns)
              if(numberOfSelectedColumns == 0) stop("Number of selected columns is zero")
              
              parameters_train_all$matrix_train <- parameters_train_all$matrix_train[, selectedColumns]
              classifier$fragmentMassSelection          <- which(selectedColumns)
            } else {
              classifier$fragmentMassSelection          <- NULL
            }
            
            classifierAll <- do.call(what = methodFunctions$train,    args = parameters_train_all)
            classifier$classifier <- classifierAll
            
            ####################################################
            ## fragment statistics
            
            ## TODO Chi-squared test
            ## https://www.ncbi.nlm.nih.gov/pubmed/16518876
            ## ... counting the number of occurrences for each m/z value, for both the good and the bad spectra. 
            ## If an m/z value is significantly over-represented in either the good or bad spectra, the chi-square test will reveal this.
            
            returnObj <- fragmentStatistics(matrix = matrix_train_all, classes = classes_pm_train_all)
            classifier$frequentFragments       <- returnObj$frequentFragments
            classifier$characteristicFragments <- returnObj$characteristicFragments
            classifier$minimumFrequency        <- returnObj$minimumFrequency
            classifier$frequentFragmentsMeanIntensity       <- returnObj$frequentFragmentsMeanIntensity
            classifier$characteristicFragmentsMeanIntensity <- returnObj$characteristicFragmentsMeanIntensity
            
            if(grepl(x = algoName, pattern = "^method=caret; .*$")){
              ## TODO caret vs caretStacked?
              #thisMethod <- "method=caret; smoothIntensities=TRUE, modelName=C5.0"
              #thisMethod <- "method=caretStacked; smoothIntensities=FALSE, modelNames=slda|lda|pam|pls|binda|OneR|rpart2|svmLinear|rpart1SE|glm|gam|C5.0Tree|C5.0Rules"
              importanceDf <- varImp(object=classifierAll)[["importance"]]
              importance <- importanceDf$"p"
              if(!is.null(importance)){
                features <- rownames(importanceDf)
                features <- gsub(x = features, pattern = "^[mp]", replacement = "")
                names(importance) <- features
                importance <- importance[importance >= max(importance) * 0.2]
                importance <- sort(x = importance, decreasing = TRUE)
              } else {
                importance <- NA
              }
              classifier$importantFragments <- importance
            } else {
              classifier$importantFragments <- NA
            }
            ## TODO Chi-squared Test of Independence?
            
            #results__class_cfi <- c(results__class_cfi, classifier$characteristicFragmentsMeanIntensity[[1]])
            
            ####################################################
            ## compute quantiles and downsample to max 10000 values
            positiveScores <- sort(unlist(predictedScores_pos_j))
            negativeScores <- sort(unlist(predictedScores_neg_j))
            
            quantiles <- seq(from = 0, to = 1, by = 0.001)
            quantilesValuesPositive <- sapply(X = quantiles, FUN = function(x){
              idx <- length(positiveScores) * x
              if(idx == 0)
                idx = 1
              positiveScores[[idx]]
            })
            quantilesValuesNegative <- sapply(X = quantiles, FUN = function(x){
              idx <- length(negativeScores) * x
              if(idx == 0)
                idx = 1
              negativeScores[[idx]]
            })
            
            if(length(positiveScores) > maximumNumberOfScores)
              positiveScores <- positiveScores[round(seq(from = 1, to = length(positiveScores), length.out = maximumNumberOfScores))]
            if(length(negativeScores) > maximumNumberOfScores)
              negativeScores <- negativeScores[round(seq(from = 1, to = length(negativeScores), length.out = maximumNumberOfScores))]
            
            classifier$quantiles <- quantiles
            classifier$quantilesValuesPositive <- quantilesValuesPositive
            classifier$quantilesValuesNegative <- quantilesValuesNegative
            
            classifier$positiveScores <- positiveScores
            classifier$negativeScores <- negativeScores
            
            ####################################################
            ## tests for alternative substance classes
            
            ## overrepresented substance classes in the negatives by BiNChe?
            ## Fisher's_exact_test for alternative (and different) substance classes
            
            ##############################
            ## collect pos and neg spectra
            selectedSpectraPos <- which(selectedSpectra)
            selectedSpectraNeg <- which(!selectedSpectra)
            if(length(duplicatedSpectrumIndecesToRemove) > 0){
              selectedSpectraPos <- which( selectedSpectra & spectraToRetain)
              selectedSpectraNeg <- which(!selectedSpectra & spectraToRetain)
            }
            
            ##############################
            ## low vs high scores
            
            ## score threshold
            fdrThreshold <- 0.05
            negativeQuantile <- quantilesValuesNegative[[which.min(abs(quantiles  - (1 - fdrThreshold)))]]
            
            ## pos scores
            predictedScores_pos_j_tmp <- unlist(predictedScores_pos_j)
            table_pos_lowScores  <- table(names(sort(predictedScores_pos_j_tmp[predictedScores_pos_j_tmp <  negativeQuantile])))
            table_pos_highScores <- table(names(sort(predictedScores_pos_j_tmp[predictedScores_pos_j_tmp >= negativeQuantile])))
            
            features_pos_all    <- sort(unique(as.numeric(c(names(table_pos_highScores), names(table_pos_lowScores)))))
            features_pos_isHigh <- sapply(X = features_pos_all, FUN = function(id){
              high <- ifelse(test = id %in% names(table_pos_highScores), yes = table_pos_highScores[names(table_pos_highScores)==id], no = 0)
              low  <- ifelse(test = id %in% names(table_pos_lowScores),  yes = table_pos_lowScores [names(table_pos_lowScores )==id], no = 0)
              return(high - low >= 2)
            })
            features_pos_high <- features_pos_all[ features_pos_isHigh]
            features_pos_low  <- features_pos_all[!features_pos_isHigh]
            
            ## neg scores
            predictedScores_neg_j_tmp <- unlist(predictedScores_neg_j)
            table_neg_lowScores  <- table(names(sort(predictedScores_neg_j_tmp[predictedScores_neg_j_tmp <  negativeQuantile])))
            table_neg_highScores <- table(names(sort(predictedScores_neg_j_tmp[predictedScores_neg_j_tmp >= negativeQuantile])))
            
            features_neg_all    <- sort(unique(as.numeric(c(names(table_neg_highScores), names(table_neg_lowScores)))))
            features_neg_isHigh <- sapply(X = features_neg_all, FUN = function(id){
              high <- ifelse(test = id %in% names(table_neg_highScores), yes = table_neg_highScores[names(table_neg_highScores)==id], no = 0)
              low  <- ifelse(test = id %in% names(table_neg_lowScores),  yes = table_neg_lowScores [names(table_neg_lowScores )==id], no = 0)
              return(high - low >= 2)
            })
            features_neg_high <- features_neg_all[ features_neg_isHigh]
            features_neg_low  <- features_neg_all[!features_neg_isHigh]
            
            ##############################
            ## confusion matrix -> p-value
            spectra_InClass_pos_highScore  <- length(features_pos_high)
            spectra_InClass_pos_lowScore   <- length(features_pos_low)
            spectra_OutClass_pos_highScore <- length(features_neg_high)
            spectra_OutClass_pos_lowScore  <- length(features_neg_low)
            challenge.df_pos = matrix(c(spectra_InClass_pos_highScore,spectra_InClass_pos_lowScore,spectra_OutClass_pos_highScore,spectra_OutClass_pos_lowScore), nrow = 2)
            pValue <- fisher.test(challenge.df_pos)[["p.value"]]
            
            ##############################
            ## test all classes for under-/overrepresentedness
            if(computeFragmentFisherTest){
              alternativeSubstanceClasses_class  <- NULL
              alternativeSubstanceClasses_pValue <- NULL
              differentSubstanceClasses_class    <- NULL
              differentSubstanceClasses_pValue   <- NULL
              
              ## TODO classIdx2 is broken for subSets
              for(classIdx2 in seq_along(classes)){
                ##############################
                ## selected spectra
                #class3 <- classes[[classIdx2]]
                #class3 <- gsub(x = class3, pattern = "\\\\", replacement = "\\\\\\\\")
                #class3 <- gsub(x = class3, pattern = "\\[", replacement = "\\\\[")
                #class3 <- gsub(x = class3, pattern = "\\]", replacement = "\\\\]")
                #selectedSpectraHere2   <- unlist(lapply(X = trueClasses, FUN = function(x){any(stri_startswith_fixed(pattern = class3, str = x))}))
                selectedSpectraHere2 <- selectedSpectra_class[[classIdx2]]
                selectedSpectraInClass <- which(selectedSpectraHere2)
                if(length(duplicatedSpectrumIndecesToRemove) > 0)
                  selectedSpectraInClass <- which(selectedSpectraHere2 & spectraToRetain)
                
                ##############################
                ## confusion matrix -> p-value
                spectra_InClass_highScore  <- sum(  features_neg_high %in% selectedSpectraInClass )
                spectra_InClass_lowScore   <- sum(  features_neg_low  %in% selectedSpectraInClass )
                spectra_OutClass_highScore <- sum(!(features_neg_high %in% selectedSpectraInClass))
                spectra_OutClass_lowScore  <- sum(!(features_neg_low  %in% selectedSpectraInClass))
                challenge.df = matrix(c(spectra_InClass_highScore,spectra_InClass_lowScore,spectra_OutClass_highScore,spectra_OutClass_lowScore), nrow = 2)
                pValue <- fisher.test(challenge.df)[["p.value"]]
                if(
                  pValue < 0.05 &
                  spectra_InClass_highScore / spectra_OutClass_highScore > spectra_InClass_lowScore / spectra_OutClass_lowScore
                ){
                  ## overrepresented in high-scoring spectra
                  #print(paste(pValue, classIdx2, classes[[classIdx2]]))
                  alternativeSubstanceClasses_class  <- c(alternativeSubstanceClasses_class,  classes[[classIdx2]])
                  alternativeSubstanceClasses_pValue <- c(alternativeSubstanceClasses_pValue, pValue)
                } else {
                  if(pValue < 0.05){
                    ## underrepresented in high-scoring spectra
                    #print(paste(spectra_InClass_highScore, spectra_OutClass_highScore, spectra_InClass_lowScore, spectra_OutClass_lowScore))
                    differentSubstanceClasses_class  <- c(differentSubstanceClasses_class,  classes[[classIdx2]])
                    differentSubstanceClasses_pValue <- c(differentSubstanceClasses_pValue, pValue)
                  }
                }
              }
              
              ## box
              alternativeSubstanceClasses  <- alternativeSubstanceClasses_pValue
              differentSubstanceClasses    <- differentSubstanceClasses_pValue
              names(alternativeSubstanceClasses) <- alternativeSubstanceClasses_class
              names(differentSubstanceClasses)   <- differentSubstanceClasses_class
              
              # install.packages("fitdistrplus")
              ## plot a Cullen and Frey graph
              #fitdistrplus::descdist(variance_vector)
              ## Kolmogorov-Smirnoff-Test
              #ks.test
              #}
              classifier$alternativeSubstanceClasses <- alternativeSubstanceClasses
              classifier$differentSubstanceClasses   <- differentSubstanceClasses
            }
            
            ## plot scores distribution
#            if(FALSE){
#              predictedScores_neg_here2 <- unlist(predictedScores_neg_j)
#              predictedScores_pos_here2 <- unlist(predictedScores_pos_j)
#              predictedScores_neg_here2_sorted <- sort(predictedScores_neg_here2)
#              predictedScores_pos_here2_sorted <- sort(predictedScores_pos_here2)
#              nNeg <- length(predictedScores_neg_here2)
#              nPos <- length(predictedScores_pos_here2)
#              
#              quantileIdx <- which(mapply(function(x, y) {isTRUE(all.equal(x, y))}, classifier$quantiles, 1 - fdrThreshold))
#              if(length(quantileIdx) == 0)  stop(paste("Quantile for fdrThreshold", fdrThreshold, "not there"))
#              scoreThreshold <- classifier$quantilesValuesNegative[[quantileIdx]]
#              
#              
#              classTmp <- tail(x = strsplit(x = class, split = "; ")[[1]], n=1)
#              file <- paste("/home/htreutle/Events/181018 FoSem/ScoreDistribution_", classTmp, ".png", sep = "")
#              widthInInch     <- 5
#              heigthInInch    <- 5
#              resolutionInDPI <- 600
#              widthInPixel    <- widthInInch  * resolutionInDPI
#              heightInPixel   <- heigthInInch * resolutionInDPI
#              png(filename = file, width = widthInPixel, height = heightInPixel, res = resolutionInDPI, bg = "white")
#              
#              main = paste(classTmp, "; AUC=", format(x = auc_roc, digits=4), sep = "")
#              plot(x = c(0,1), y = c(0,1), type="n", xlab="Score", ylab="Cumulative distribution", main=main)
#              par(new=TRUE)
#              plot(x = predictedScores_pos_here2_sorted, y = (1:nPos)/nPos, type="s", col = "blue", xlab="", ylab="", axes=F)
#              par(new=TRUE)
#              plot(x = predictedScores_neg_here2_sorted, y = (1:nNeg)/nNeg, type="s", col = "red", xlab="", ylab="", axes=F)
#              segments(x0 = scoreThreshold, x1 = scoreThreshold, y0 = 0, y1 = 1)
#              text(x = scoreThreshold+0.015, y = 0.5, labels = paste("FPR=", format(x = fdrThreshold, digits=4), "; TPR=", format(x = tpr, digits=4), sep = ""), adj = c(0,0.5))
#              #segments(x0 = 0, x1 = 1, y0 = 0.95, y1 = 0.95)
#              legend(x = 0.55, y = 0.20, legend = c("Foreground", "Background"), col = c("blue", "red"), lty=1)
#              dev.off()
#            }
          }
          ####################################################
          ## store classifier and results
          results__algo_class    [[methodIdx]][[class]] <- resultRow
          classifiers__algo_class[[methodIdx]][[class]] <- classifier
        }## method
      }## builtClassifiers
    }## class
    
#    ## TODO remove redundant pos sets
#    if(FALSE){
#      #classifierFile      <- "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2019-05-09_08:40:15_180912_MoNA-export-LC-MS-MS_Positive_Mode_processed.msp__Classifier.RData"
#      #fileSplitClassRData <-                               "/vol/Workspace/SEB/Bioinformatics_RG/htreutle_data/Spectral_libraries/MONA/180912/180912_MoNA-export-LC-MS-MS_Positive_Mode_processed_5ece85cd_splitClasses.RData"
#      #load(classifierFile)
#      #load(fileSplitClassRData)
#      classesWithSubSets <- names(classifiers_class)
#      isSubSet <- grepl(x = classesWithSubSets, pattern = "; SubSet_\\d+$")
#      subSetClasses <- classesWithSubSets[isSubSet]
#      
#      selectedSpectra_class2 <- lapply(X = selectedSpectra_class, FUN = function(selectedSpectraHere){return(selectedSpectraHere & spectraToRetain)})
#      redundantSubSets_toClass   <- unlist(sapply(X = classes, FUN = function(class){lapply(X = results__subSets[[class]], FUN = function(results__subSet){
#        is.logical(all.equal(current = selectedSpectra_class2[[class]], target = results__subSet))
#      })}))
#      
#      results__subSets_unlist <- unlist(x = results__subSets, recursive = F, use.names = F)
#      redundantSubSets_toSubSetsTmp <- duplicated(results__subSets_unlist, recursive = F)
#      subSetClassIndeces <- unique(unlist(sapply(X = results__subSets_unlist[redundantSubSets_toSubSetsTmp], FUN = function(results__subSet){
#        #print(sum(results__subSet))
#        #results__subSet <- results__subSets_unlist[redundantSubSets_toSubSetsTmp][[5]]
#        #indeces <- which(unlist(lapply(X = results__subSets_unlist, FUN = function(results__subSets_unlistHere){ is.logical(all.equal(target = results__subSets_unlistHere, current = results__subSet))})))
#        indeces <- which(unlist(lapply(X = results__subSets_unlist, FUN = function(results__subSets_unlistHere){ all(results__subSets_unlistHere == results__subSet)})))
#        subSetClassesHereWithSubSet <- subSetClasses[indeces]
#        subSetClassesHere <- gsub(x = subSetClassesHereWithSubSet, pattern = "; SubSet_\\d+$", replacement = "")
#        #print(length(indeces))
#        subSetClassesHereRedundant <- sapply(X = seq_along(subSetClassesHere), FUN = function(subSetClassIdx){any(grepl(x = subSetClassesHere[-subSetClassIdx], pattern = subSetClassesHere[[subSetClassIdx]]))})
#        subSetClassIndeces <- match(x = subSetClassesHereWithSubSet[subSetClassesHereRedundant], table = subSetClasses)
#        return(subSetClassIndeces)
#      })))
#      redundantSubSets_toSubSets <- logical(length = length(subSetClasses))
#      redundantSubSets_toSubSets[subSetClassIndeces] <- TRUE
#      
#      redundantSubSets <- redundantSubSets_toClass | redundantSubSets_toSubSets
#      print(paste(sum(redundantSubSets_toClass), " / ", length(subSetClasses) , " equal to the class; ", sum(redundantSubSets_toSubSets), " present in more general class(es), ", (length(subSetClasses) - sum(redundantSubSets)), " / ", length(subSetClasses), " SubSets remain for ", length(classes), " classes", sep = ""))
#      
#      redundantClasses <- logical(length = length(classesWithSubSets))
#      redundantClasses[isSubSet] <- redundantSubSets
#      
#      ## remove
#      classifiers_class <- classifiers_class[!redundantClasses]
#      # results__subSets ## leave as is
#      #for(idx in seq_along(classifiers__algo_class)){
#      #  if(length(classifiers__algo_class[[idx]] != length(redundantClasses))) stop("ui")
#      #  #classifiers__algo_class[[idx]]
#      #}
#    }
    
    ## TODO 999
    if(all(splitClasses, !file.exists(fileSplitClassRData))){#, !is.null(thisClass)))
      print(fileSplitClassRData)
      save(file = fileSplitClassRData, results__subSets)
    }
#    if(FALSE){
#      ## all pos: 153/1283 no subSets
#      tmp <- lapply(X = classes, FUN = function(class){
#        inClass <- sum(selectedSpectra_class[[class]] & spectraToRetain)
#        numberOfSubSets <- length(results__subSets[[class]])
#        numberOfSpectraInSubSets <- sum(unlist(results__subSets[[class]]))
#        c(inClass, numberOfSubSets)
#      })
#      xe <- unlist(lapply(X = classes, FUN = function(class){  sum(selectedSpectra_class[[class]] & spectraToRetain)  }))
#      ye <- unlist(lapply(X = classes, FUN = function(class){  length(results__subSets[[class]])    }))
#      ze <- unlist(lapply(X = classes, FUN = function(class){  sum(unlist(lapply(selectedSpectra_class[[class]], FUN = function(vec){vec & spectraToRetain})))  }))
#      plot(xe, ye, log="x")
#      plot(xe, ze, log="xy")
#      lines(x = c(1,max(xe)), y = c(1,max(xe)))
#      plot(seq(from=0, to = max(ye)), sapply(X = seq(from=0, to = max(ye)), FUN = function(c){sum(ye==c)}), log="y", xlab="Number of subsets", ylab="Number of classes", main=gsub(x = gsub(x = libraryName, pattern = "^\\d\\d\\d\\d\\-\\d\\d\\-\\d\\d_\\d\\d:\\d\\d:\\d\\d_", replacement = ""), pattern = "\\.msp$", replacement = ""))
#    }
    
    #######################################################################
    #######################################################################
    ## out
    timeStamp <- gsub(pattern = "[ ]", x = Sys.time(), replacement = "_")
    if(writeClassifiers | writeResults){
      #setupName <- paste(
      #  classifiers_class[[1]]$classifierName, "_",
      #  "NC=", numberOfClasses, "_",
      #  ifelse(all(allowedInstrumentTypes=="all"), "IT=ALL",""), "_",
      #  "MS=", mergeDuplicatedSpectra, "_",
      #  ifelse(unitResolution, "UR=TRUE",""),
      #  sep = ""
      #)
      
      ## path       80
      ## timestamp  20
      ## filetype   11
      ## library    8+x+5 x=80
      ## classifier 80
      ## properties 30
      ## extension   5
      ##            239+x / 255
      libraryName2 <- libraryName
      suffixHere   <- ifelse(test = !is.null(outSuffix), yes = paste(outSuffix, "_", sep = ""), no = "_")
      classifierFile <- paste(resultFolderForClassifiers, "/", timeStamp, "_", libraryName2, "_", suffixHere, "Classifier", ".RData", sep = "")
      resultFile     <- paste(resultFolderForClassifiers, "/", timeStamp, "_", libraryName2, "_", suffixHere, "Results",    ".tsv",   sep = "")
      propertiesFile <- paste(resultFolderForClassifiers, "/", timeStamp, "_", libraryName2, "_", suffixHere, "Properties", ".txt",   sep = "")
      
      ## write results and classifier(s)
      for(methodIdx in seq_along(algorithms)){
        ## get data
        classifiers_class <- classifiers__algo_class[[methodIdx]]
        results__class    <- results__algo_class    [[methodIdx]]
        
        if(writeResults){
          resultDf <- as.data.frame(t(matrix(data = unlist(results__class), ncol = length(results__class))), stringsAsFactors = FALSE)
          colnames(resultDf) <- colNamesResults
          
          ## write results
          write.table(file = resultFile,     x = resultDf, row.names = FALSE, sep = "\t")
          print(resultFile)
        }
        if(writeClassifiers){
          ###################################################################################################
          ## properties file
          properties <- list()
          
          ## library
          properties$library <- libraryName
          ## I/O
          properties$maximumNumberOfScores   <- maximumNumberOfScores
          properties$removeRareFragments     <- removeRareFragments
          ## replicates
          properties$mergeDuplicatedSpectra              <- mergeDuplicatedSpectra
          properties$takeSpectraWithMaximumNumberOfPeaks <- takeSpectraWithMaximumNumberOfPeaks
          ## classes
          #properties$considerAlternativeSubstanceClasses <- considerAlternativeSubstanceClasses
          #properties$processSubstanceClasses <- processSubstanceClasses
          properties$classOfClass <- classOfClass
          ## repeated random subsampling validation
          properties$minimumNumberOfPosSpectraPerClass <- minimumNumberOfPosSpectraPerClass
          properties$minimumNumberOfNegSpectraPerClass <- minimumNumberOfNegSpectraPerClass
          properties$numberOfDataSetDecompositions     <- numberOfDataSetDecompositions
          properties$proportionTraining      <- proportionTraining
          ## postprocessing
          properties$unitResolution <- unitResolution
          
          ## algo
          algorithm <- algorithms[[methodIdx]]
          #methodFunctions   <- algorithm$method
          properties$methodName        <- algorithm$methodName
          #properties$parameters_train  <- algorithm$params
          properties$paramsString      <- algorithm$paramsString
          properties$algoName          <- algorithm$algoName
          
          ## write results
          save(       file = classifierFile, classifiers_class)
          #classifiers_class <- classifiers_class[!grepl(x = names(classifiers_class), pattern = "; SubSet_\\d+$")]
          #save(classifiers_class, file = gsub(x = classifierFile, pattern = "\\.msp__Classifier\\.RData$", replacement = paste(".msp__", "NoSubSets_", "Classifier.RData", sep = "")))
          
          writeLines( con  = propertiesFile, text = paste(names(properties), properties, sep = " = "))
          print(classifierFile)
          print(propertiesFile)
        }
      }## methods
    }## write something
  }## libraries
  
  #######################################################################
  #######################################################################
  ## misc
  
  ## AUC-ROC vs AUC-PR
#  if(FALSE){
#    resultFile <- "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-10-19_11:29:59_2018-02-13_09_14_10_neg_11328_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv"
#    resultTable <- read.table(file = resultFile, header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
#    plot(x = resultTable$AUC, resultTable$AUC.PR, xlim=c(0.5,1), ylim=c(0.,1))
#    lines(x = c(0.5, 1), y = c(0., 1))
#  }
  ## AUC vs AUC plot
#  if(FALSE){
#    ## considerAlternativeSubstanceClasses ?
#    resultDf1 <- read.table(file = "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=MoNA-export-GC-MS_-_EI_-_Positive.msp_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv",
#                            header = T, sep = "\t", comment.char = "", stringsAsFactors = F, check.names = F)
#    resultDf2 <- read.table(file = "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=mainlibNIST2014_withINCHI.MSP_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv",
#                            header = T, sep = "\t", comment.char = "", stringsAsFactors = F, check.names = F)
#    
#    resultDf2 <- resultDf2[resultDf2$`Substance class` %in% resultDf$`Substance class`, ]
#    
#    plot(x = resultDf$AUC, y = resultDf3$AUC, xlim=c(0.5,1), ylim=c(0.5,1), xlab="AUC Without alternative SC", ylab="AUC Without alternative SC")
#    segments(0.5,0.5,1,1)
#    
#    plot(sort(resultDf$AUC), col="red", ylim=c(0.5,1))
#    points(sort(resultDf1$AUC), col="blue")
#    points(sort(resultDf2$AUC), col="green")
#    legend(x = 1, y = 1, legend = c("1 vs 2", "1", "2"), col = c("red", "blue", "green"), lty = c(1,1,1))
#  }
  ## method vs method plot
#  if(FALSE){
#    ## fetch mean AUC per method
#    for(methodName in methodNames){
#      for(smoothIntensities in c(T,F)){
#        these <- resultDf$Method==methodName & resultDf$SmoothIntensities==as.character(smoothIntensities)
#        print(paste(
#          sum(as.numeric(resultDf$AUC[these]), na.rm = TRUE) / sum(these, na.rm = TRUE), 
#          sum(as.numeric(resultDf$`time (s)`[these]), na.rm = TRUE) / sum(these, na.rm = TRUE), 
#          methodName, smoothIntensities
#        ))
#      }
#    }
#  }
  ## combine results
#  if(FALSE){
    ## combine library results
    #libraryNames <- c("GC_EI_pos", "LC_ESI_neg", "LC_MSMS_neg", "LC_MSMS_pos")
#    libraryNames <- c(
#      "170713_Weizmass_neg",
#      "170713_Weizmass_pos",
#      "Wiley_and_Wiley_10",
#      "mainlibNIST2014_GC",
#      "MoNA_GC_EI_Pos",
#      "MoNA_LC_ESI_Neg",
#      "MoNA_LC_MSMS_Neg",
#      "MoNA_LC_MSMS_Pos",
#      "UNPD_ISDB_R_LC_Neg"
#    )
#    files <- c(
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=170713_Shahaf_lib_spectra_neg.msp_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=170713_Shahaf_lib_spectra_neg.msp_Class=ChemOnt_Substituent_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=170713_Shahaf_lib_spectra_pos.msp_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=170713_Shahaf_lib_spectra_pos.msp_Class=ChemOnt_Substituent_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=GCMS_DB-Wiley_and_Wiley_10.msp_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=GCMS_DB-Wiley_and_Wiley_10.msp_Class=ChemOnt_Substituent_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=mainlibNIST2014_withINCHI.MSP_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=mainlibNIST2014_withINCHI.MSP_Class=ChemOnt_Substituent_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=MoNA-export-GC-MS_-_EI_-_Positive.msp_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=MoNA-export-GC-MS_-_EI_-_Positive.msp_Class=ChemOnt_Substituent_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=MoNA-export-LC-MS_-_ESI_-_Negative.msp_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=MoNA-export-LC-MS_-_ESI_-_Negative.msp_Class=ChemOnt_Substituent_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=MoNA-export-LC-MS_-_MSMS_-_Negative.msp_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=MoNA-export-LC-MS_-_MSMS_-_Negative.msp_Class=ChemOnt_Substituent_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=MoNA-export-LC-MS_-_MSMS_-_Positive.msp_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=MoNA-export-LC-MS_-_MSMS_-_Positive.msp_Class=ChemOnt_Substituent_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=UNPD_ISDB_R.mgf_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_library=UNPD_ISDB_R.mgf_Class=ChemOnt_Substituent_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE.tsv"
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_GC_EI_pos_SVM_one.tsv",
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_LC_ESI_Neg_SVM_one.tsv",
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_LC_MSMS_neg_SVM_one.tsv",
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_LC_MSMS_pos_SVM_one.tsv"
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_GC_EI_Pos_Substituents_Detailed.tsv",
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_LC_ESI_Neg_Substituents_Detailed.tsv",
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_LC_MSMS_Neg_Substituents_Detailed.tsv",
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_LC_MSMS_Pos_Substituents_Detailed.tsv"
#    )
#    
#    if(TRUE){
#      files <- files[seq(from=1, to=length(files), by=2)]
#    } else {
#      files <- files[seq(from=2, to=length(files), by=2)]
#    }
#    
#    librarySelection <- c(T,T,F,F,F,F,T,T,F)
#    files        <- files       [librarySelection]
#    libraryNames <- libraryNames[librarySelection]
#    
#    
#    libraryNames <- c(
#      "MoNA_GC_EI_Pos",
#      "MoNA_LC_MSMS_Pos",
#      "MoNA_LC_MSMS_Neg"
#    )
#    
#    files <- c(
#      ## GC-EI
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-08_08:55:24_Results_library=MoNA-export-GC-MS_-_EI_-_Positive.msp_Class=ChemOnt_Substituent_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE_NumClasses=512_mergedSpectra.tsv",
#      ## LC-pos
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-08_08:51:16_Results_library=MoNA-export-LC-MS_-_MSMS_-_Positive.msp_Class=ChemOnt_Substituent_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE_NumClasses=351_mergedSpectra.tsv",
#      ## LC-neg
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-08_08:45:32_Results_library=MoNA-export-LC-MS_-_MSMS_-_Negative.msp_Class=ChemOnt_Substituent_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE_NumClasses=286_mergedSpectra.tsv"
#    )
#    
#    
#    libraryNames <- c(
#      #"MoNA",
#      "Merge_neg",
#      "NIST2017_neg",
#      
#      #"Respect_pos",
#      "Respect_neg",
#      #"GNPS_pos",
#      "GNPS_neg",
#      #"IPB_pos",
#      "IPB_neg",
#      #"WeizMass_pos",
#      "WeizMass_neg",
#      #"MoNA_pos",
#      "MoNA_neg"
#    )
#    files <- c(
#      ## GC-EI
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-06_15:30:43_Results_library=MoNA-export-GC-MS_-_EI_-_Positive.msp_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE_NumClasses=502_mergedSpectra.tsv",
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-06_22:02:17_Results_library=mainlibNIST2014_withINCHI.MSP_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE_NumClasses=1766_mergedSpectra.tsv",
#      ## LC-pos
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-06_15:26:28_Results_library=170713_Shahaf_lib_spectra_pos.msp_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE_NumClasses=282_mergedSpectra.tsv",
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-06_15:23:24_Results_library=MoNA-export-LC-MS_-_MSMS_-_Positive.msp_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE_NumClasses=374_mergedSpectra.tsv",
#      ## LC-neg
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-06_15:18:15_Results_library=170713_Shahaf_lib_spectra_neg.msp_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE_NumClasses=269_mergedSpectra.tsv",
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-06_15:16:54_Results_library=MoNA-export-LC-MS_-_MSMS_-_Negative.msp_Class=ChemOnt_SubstanceClass_AltSC=TRUE_method=ColSums_smoothIntensities=FALSE_NumClasses=310_mergedSpectra.tsv"
#      
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-15_08:34:16_MoNA-export-LC-MS_-_MSMS_-_Negative.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-14_15_52_58_2018-02-09_15_00_25_MSMS_neg_106843_MoNA_WeizMASS_NIST2017_IPB_GNPS_Respect.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-15_15_22_21_2018-02-13_09_14_10_neg_91215_nist_msms.MSP_Results.tsv",
#      
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-20_09:23:12_2018-02-13_09_14_10_pos_2768_MSMS-Respect-Curated-Pos.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-20_09_21_21_2018-02-13_09_14_10_neg_1628_MSMS-Respect-Curated-Neg.msp_Results.tsv",
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-20_09:20:10_2018-02-13_09_14_10_pos_3771_MSMS-GNPS-Curated-Pos.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-20_09_02_41_2018-02-13_09_14_10_neg_153_MSMS-GNPS-Curated-Neg.msp_Results.tsv",
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-20_09:02:03_2018-02-13_09_14_10_pos_111_Leibniz positive mode_20170801.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-20_09:01:49_2018-02-13_09_14_10_neg_292_Leibniz negative mode_20170801.msp_Results.tsv",
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-20_08:39:25_2018-02-13_09_14_10_pos_2214_lib_spectra_pos_new.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-20_08:35:32_2018-02-13_09_14_10_neg_2227_lib_spectra_neg_new.msp_Results.tsv",
#      #"/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-20_08:32:47_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-02-20_07:34:28_2018-02-13_09_14_10_neg_11328_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv"
#    )
#    
#    libraryNames <- c(
#      "Merge_pos",
#      "WeizMass_pos",
#      "MoNA_pos",
#      "Merge_neg",
#      "WeizMass_neg",
#      "MoNA_neg"
#    )
#    files <- c(
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-03-12_15:58:45_2018-02-09_15_00_25_MSMS_pos_418622_MoNA_WeizMASS_NIST2017_IPB_GNPS_Respect.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-03-12_14:01:16_2018-02-13_09_14_10_pos_2214_lib_spectra_pos_new.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-03-12_14:00:05_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-03-12_13:59:06_2018-02-09_15_00_25_MSMS_neg_106843_MoNA_WeizMASS_NIST2017_IPB_GNPS_Respect.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-03-12_13:42:15_2018-02-13_09_14_10_neg_2227_lib_spectra_neg_new.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/2018-03-12_13:41:33_2018-02-13_09_14_10_neg_11328_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv"
#    )
#    
#    libraryNames <- c(
#      "merge_ColSums; sI=FALSE",
#      "merge_ColSums; sI=TRUE",
#      "merge_ColSumsPos; sI=FALSE",
#      "merge_ColSumsPos; sI=TRUE",
#      "merge_ColSumsPosOnly; sI=FALSE",
#      "merge_ColSumsPosOnly; sI=TRUE",
#      "merge_ColSumsRatio; sI=FALSE",
#      "merge_ColSumsRatio; sI=TRUE",
#      "merge_ColSumsLog; sI=FALSE",
#      "merge_ColSumsLog; sI=TRUE",
#      "max_ColSums; sI=FALSE",
#      "max_ColSums; sI=TRUE",
#      "max_ColSumsPos; sI=FALSE",
#      "max_ColSumsPos; sI=TRUE",
#      "max_ColSumsPosOnly; sI=FALSE",
#      "max_ColSumsPosOnly; sI=TRUE",
#      "max_ColSumsRatio; sI=FALSE",
#      "max_ColSumsRatio; sI=TRUE",
#      "max_ColSumsLog; sI=FALSE",
#      "max_ColSumsLog; sI=TRUE"
#    )
#    files <- c(
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_13:48:33_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_13:49:25_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_13:50:15_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_13:51:06_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_13:51:47_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_13:52:29_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_13:53:21_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_13:54:12_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_13:55:04_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_13:55:57_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_13:56:49_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_13:57:39_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_13:58:30_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_13:59:20_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_14:00:01_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_14:00:44_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_14:01:36_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_14:02:28_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_14:03:21_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_scaffolds/2018-06-08_14:04:15_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv"
#    )
#    files <- c(
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_08:10:36_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_08:16:59_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_08:23:43_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_08:30:21_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_08:35:46_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_08:41:10_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_08:47:54_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_08:54:21_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_09:00:56_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_09:07:18_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_09:13:17_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_09:19:15_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_09:25:16_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_09:31:19_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_09:36:20_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_09:41:23_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_09:47:50_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_09:54:29_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_10:00:51_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv",
#      "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/colSumsTest_ChemOntAll/2018-06-13_10:07:26_2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp_Results.tsv"
#    )
#    
#    if(length(files) != length(libraryNames)) stop("lengths not equal")
#    
#    takeSubstanceClasses <- TRUE
#    
#    fileOut <- paste("/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/tmp/Results_combined_", ifelse(takeSubstanceClasses, "SC", "Subst"), "_", length(files), "Libraries.tsv", sep = "")
#    #fileOut <- paste("/home/htreutle/Events/180208 SEB seminar/pics/Results_combined_", ifelse(takeSubstanceClasses, "SC", "Subst"), "_", length(files), "Libraries.tsv", sep = "")
#    #fileOut <- paste("/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_combined_", ifelse(takeSubstanceClasses, "SC", "Subst"), "_", length(files), "Libraries.tsv", sep = "")
#    #fileOut <- "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/Results_combined_Subst_4Libraries.tsv"
#    
#    tables <- list()
#    for(idx in seq_along(files))
#      tables[[idx]] <- read.table(file = files[[idx]], header = T, sep = "\t", comment.char = "", check.names = F, stringsAsFactors = F, quote = "\"")
#    scs <- NULL
#    for(idx in seq_along(files))
#      scs <- c(scs, tables[[idx]][, "Substance class"])
#    scs <- unique(scs)
#    scs <- sort(scs)
#    
#    resultDfAll <- data.frame(
#      "Substance class" = scs, 
#      stringsAsFactors = F, check.names = F
#    )
#    ## empty columns for spectra count and AUC
#    for(idx in seq_along(files))
#      resultDfAll = cbind(
#        resultDfAll, 
#        rep(x = "", times = nrow(resultDfAll)), 
#        rep(x = "", times = nrow(resultDfAll)), 
#        stringsAsFactors = F
#    )
#    
#    colnames(resultDfAll) <- c(
#      ifelse(test = takeSubstanceClasses, yes = "Substance class", no = "Substituent"), 
#      paste("Spectra", libraryNames), 
#      paste("AUC", libraryNames)
#    )
#    
#    for(sc in scs){
#      for(idx in seq_along(files)){
#        from = which(tables[[idx]][, "Substance class"] == sc)
#        to   = which(resultDfAll  [, "Substance class"] == sc)
#        if(length(from) == 0){
#          num <- "-"
#          auc <- "-"
#        } else {
#          num <- paste(tables[[idx]][,"Number of positive spectra"][from], tables[[idx]][, "Number of spectra"][from], sep = " / ")
#          auc <- tables[[idx]][, "AUC"][from]
#        }
#        
#        resultDfAll[to, 1 + idx] <- num
#        resultDfAll[to, 1 + length(tables) + idx] <- auc
#      }
#    }
#    
#    write.table(file = fileOut, x = resultDfAll, row.names = FALSE, sep = "\t")
#    print(fileOut)
#    
#    if(FALSE)
#      for(idx in seq_along(files)){
#        fileOut2 <- gsub(pattern = ".tsv", replacement = paste("_", libraryNames[[idx]], "hist.png", sep = ""), x = paste(fileOut))
#        
#        widthToHeightRatio = 16/10;
#        size = 6
#        resolution <- 300
#        width <- size * widthToHeightRatio
#        height <- size
#        
#        fileName <- gsub(pattern = ".tsv$", replacement = "_hist_2.png", x = files[[idx]])
#        png(filename = fileName, res = resolution, width = width * resolution, height = height * resolution)
#        
#        hist(as.numeric(resultDfAll[, paste("AUC", libraryNames[[idx]])]), main = paste("SC", libraryNames[[idx]]), xlab = "AUC", breaks = seq(from=0, to=1, by=0.05))
#        
#        dev.off()
#      }
#    #for(idx in seq_along(files))
#    #  write.table(file = files[[idx]], x = tables[[idx]], row.names = FALSE, sep = "\t")
#  }
}

##############################################################################################################################################
## MetFamily
createMetFamilyClassFile <- function(class, numberOfPrecursors, resultFolder){
  classFile <- gsub(x = class, pattern = "; ", replacement = "__")
  classFile <- gsub(x = classFile, pattern = ", ", replacement = "_")
  classFile <- gsub(x = classFile, pattern = " ", replacement = "_")
  classFile <- gsub(x = classFile, pattern = "/", replacement = "_")
  dataSetPrefix <- paste(classFile, "___", numberOfPrecursors, sep = "")
  #file   <- paste(resultFolder, "/", dataSetPrefix, ".csv",    sep = "")
  file   <- paste(resultFolder, "/", dataSetPrefix, ".csv.gz",    sep = "")
  
  ## class string length: 10 * (2+4) + 58
  replaceIdx <- 1
  while(nchar(file) >= 260){
    classFileTokens <- strsplit(x = classFile, split = "__")[[1]]
    
    if(replaceIdx > length(classFileTokens))
      break
    
    classFileTokens[[replaceIdx]] <- paste(substr(x = classFileTokens[[replaceIdx]], start = 1, stop = 3), ".", sep = "")
    classFile <- paste(classFileTokens, collapse = "__")
    
    dataSetPrefix <- paste(classFile, "___", numberOfPrecursors, sep = "")
    file   <- paste(resultFolder, "/", dataSetPrefix, ".csv",    sep = "")
    
    replaceIdx <- replaceIdx+1
  }
  return(file)
}
writeMetFamilyProject <- function(spectraListHere, classHere, parameterSet, file){
  precursorMzHere <- unlist(lapply(X = spectraListHere, FUN = function(x){x$mz}))
  precursorRtHere <- unlist(lapply(X = spectraListHere, FUN = function(x){x$rt}))
  #precursorRtHere[is.na(precursorRtHere)] <- -1
  #precursorRtHere <- seq_along(spectraListHere) ## no duplicates
  precursorLabels <- paste(precursorMzHere, precursorRtHere, sep = " / ")
  while(any(duplicated(precursorLabels))){
    print("### WARN renaming precursors for MetFamily project")
    precursorRtHere <- precursorRtHere + runif(n = length(precursorRtHere), min = 0.01, max = 0.1)
    precursorLabels <- paste(precursorMzHere, precursorRtHere, sep = " / ")
    for(idx in seq_along(spectraListHere))
      spectraListHere[[idx]]$rt <- precursorRtHere[[idx]]
  }
  
  inchisHere    <- unlist(lapply(X = spectraListHere, FUN = function(x){x$inchi}))
  smilesHere    <- unlist(lapply(X = spectraListHere, FUN = function(x){x$smiles}))
  inchikeysHere <- unlist(lapply(X = spectraListHere, FUN = function(x){x$inchiKey}))
  
  ## return(!grepl(x = annotationValue, pattern = "(, )|(=)|(^$)|(^[ \t\n\r\v\f]+$)"))
  #classHere <- class
  classHere <- gsub(x = classHere, pattern = ", ", replacement = ",")
  classHere <- gsub(x = classHere, pattern = "=",  replacement = "_")
  
  metaboliteFamiliesHere <- rep(x = classHere, times = length(spectraListHere))
  uniqueMetaboliteFamiliesHere <- sort(unique(metaboliteFamiliesHere))
  metaboliteFamilyColorsHere <- rep(x = "black", times = length(uniqueMetaboliteFamiliesHere))
  
  parameterSet$mzDeviationAbsolute_mapping <- 0.01
  parameterSet$maximumRtDifference <- 0.01
  
  resultObj <- convertToProjectFile2(
    filePeakMatrix = NULL, 
    spectraList = spectraListHere, precursorMz = precursorMzHere, precursorRt = precursorRtHere, 
    metaboliteFamilies = metaboliteFamiliesHere, uniqueMetaboliteFamilies = uniqueMetaboliteFamiliesHere, metaboliteFamilyColors = metaboliteFamilyColorsHere, 
    furtherProperties = list("InChI" = inchisHere, "SMILES" = smilesHere, "InChIKey" = inchikeysHere),
    parameterSet = parameterSet, 
    progress = NA
  )
  
  ## to lines
  matrixRows <- resultObj$matrixRows
  matrixCols <- resultObj$matrixCols
  matrixVals <- resultObj$matrixVals
  
  matrixRows <- c(matrixRows, 1)
  matrixCols <- c(matrixCols, 1)
  matrixVals <- c(matrixVals, serializeParameterSet(parameterSet))
  
  ## convert matrix to dataframe
  numberOfRows    <- max(matrixRows)
  numberOfColumns <- max(matrixCols)
  lines <- vector(mode = "character", length = numberOfRows)
  for(rowIdx in seq_len(numberOfRows)){
    indeces <- matrixRows == rowIdx
    tokens  <- vector(mode = "character", length = numberOfColumns)
    tokens[matrixCols[indeces]] <- matrixVals[indeces]
    lines[[rowIdx]] <- paste(tokens, collapse = "\t")
  }
  
  #writeLines(text = lines, con = file)
  gz1 <- gzfile(description = file, open = "w")
  writeLines(text = lines, con = gz1)
  close(gz1)
}

##############################################################################################################################################
## data and models
highResolutionParameterSet <- list(
  ## parse msp
  minimumIntensityOfMaximalMS2peak                  = 000,
  minimumProportionOfMS2peaks                       = 0.05,
  neutralLossesPrecursorToFragments                 = TRUE,
  neutralLossesFragmentsToFragments                 = FALSE,
  doPrecursorDeisotoping                            = FALSE,
  ## built fragment matrix
  mzDeviationAbsolute_grouping                      = 0.01,
  mzDeviationInPPM_grouping                         = 20,
  doMs2PeakGroupDeisotoping                         = TRUE,
  mzDeviationAbsolute_ms2PeakGroupDeisotoping       = 0.01,
  mzDeviationInPPM_ms2PeakGroupDeisotoping          = 10,
  proportionOfMatchingPeaks_ms2PeakGroupDeisotoping = 0.9
)
gcResolutionParameterSet <- list(
  ## parse msp
  minimumIntensityOfMaximalMS2peak                  = 000,
  minimumProportionOfMS2peaks                       = 0.05,
  neutralLossesPrecursorToFragments                 = FALSE,
  neutralLossesFragmentsToFragments                 = FALSE,
  doPrecursorDeisotoping                            = FALSE,
  ## built fragment matrix
  mzDeviationAbsolute_grouping                      = 0.3,
  mzDeviationInPPM_grouping                         = 0,
  doMs2PeakGroupDeisotoping                         = TRUE,
  mzDeviationAbsolute_ms2PeakGroupDeisotoping       = 0.3,
  mzDeviationInPPM_ms2PeakGroupDeisotoping          = 0,
  proportionOfMatchingPeaks_ms2PeakGroupDeisotoping = 0.9
)

getApplicableCaretModels <- function(verbose.output = FALSE){
  models <- caret::getModelInfo()
  modelNames <- names(models)
  applicableModelNames <- vector(mode = "character", length = 0)
  for(modelName in modelNames){
    list <- models[[modelName]]
    
    if(verbose.output)  print(paste(modelName, "class", ("Classification" %in% list$type), "prob", !is.null(list$prob)))
    if(!("Classification" %in% list$type))
      next;
    if(is.null(list$prob))
      next;
    applicableModelNames <- c(applicableModelNames, modelName)
  }
  
  packagesList <- lapply(X = models[applicableModelNames], FUN = "[", "library")
  #packagesNotThere <- c("adaptDA", "mxnet")
  packagesNotThere <- c("FCNN4R", "adaptDA", "mxnet", "sparsediscrim")
  excludedModels <- unlist(lapply(X = packagesList, function(x){any(packagesNotThere %in% x)}))
  applicableModelNames <- applicableModelNames[!excludedModels]
  
  return(applicableModelNames)
}
#classifierOverview <- function(){
#  folder <- "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis/"
#  resultObj <- getAvailableClassifiers(folder)
#  availableClassifiersDf <- resultObj$availableClassifiersDf
#  availableClassifiersDfProperties <- resultObj$availableClassifiersDfProperties
#  
#  rownames(availableClassifiersDf) <- basename(rownames(availableClassifiersDf))
#  
#  timeStamp <- gsub(pattern = "[ ]", x = Sys.time(), replacement = "_")
#  folder <- "/vol/Workspace/SEB/Bioinformatics_RG/htreutle_data/Spectral_libraries/"
#  file <- paste(folder, timeStamp, "_LibraryOverview.tsv", sep = "")
#  
#  write.table(x = availableClassifiersDf, file = file, sep = "\t")
#  
#  availableClassifiersDf
#}
processLibrary <- function(
  fileSpectra, parameterSet, allowedInstrumentTypes, annoFile, structureFormats, progress, 
  minimumNumberOfPosSpectraPerClass, minimumNumberOfNegSpectraPerClass, mergeDuplicatedSpectra, #considerAlternativeSubstanceClasses, 
  unitResolution,
  reParseMsp, reProcessAnno, reProcessFragmentMatrix
){
  ## rm(list=setdiff(ls(), c("fileSpectra", "parameterSet", "allowedInstrumentTypes", "annoFile", "structureFormats", "progress", "minimumNumberOfPosSpectraPerClass", "minimumNumberOfNegSpectraPerClass", "processSpectraAndAnnotation")))
  #######################################################################
  ## parse spectra file and annotate
  returnObj <- processSpectraAndAnnotation(fileSpectra, parameterSet, allowedInstrumentTypes, annoFile, structureFormats, reParseMsp, reProcessAnno, progress)
  spectraList     = returnObj$spectraList
  numberOfSpectra = returnObj$numberOfSpectra
  precursorMz     = returnObj$precursorMz
  precursorRt     = returnObj$precursorRt
  inchis          = returnObj$inchis
  annoTable       = returnObj$annoTable
  
#  if(FALSE){
#    subSet <- 5000
#    if(!is.na(subSet)){
#      spectraList     = spectraList[1:subSet]
#      numberOfSpectra = subSet
#      precursorMz     = precursorMz[1:subSet]
#      precursorRt     = precursorRt[1:subSet]
#      inchis          = inchis[1:subSet]
#      annoTable       = annoTable[1:subSet, ]
#    }
#  }
  
  #spectrumId <- paste(precursorMz, precursorRt)
  #spectrumId <- seq_len(numberOfSpectra)
  
  ######################################################################
  ## annotate spectra
  inchiKeys          <- annoTable[, "InChIKey"]
  inchiKeysBlock1    <- unlist(lapply(X = strsplit(x = inchiKeys, split = "\\-"), FUN = function(x){x[[1]]}))
  
  #######################################################################
  ## duplicated
  print("Duplicated spectra")
  duplicatedStructures <- duplicated(x = inchiKeysBlock1)
  #duplicatedStructures <- duplicated(x = inchis)
  uniqueNumberOfSpectra <- numberOfSpectra - sum(duplicatedStructures)
  duplicatedInchiKeysBlock1 <- unique(inchiKeysBlock1[duplicatedStructures])
  
  #######################################################################
  ## built fragment matrix
  print("Fragment matrix")
  paramsHash <- digest::sha1(algo = "crc32", x = unlist(c(
    ## parse
    minimumIntensityOfMaximalMS2peak = parameterSet$minimumIntensityOfMaximalMS2peak, 
    minimumProportionOfMS2peaks = parameterSet$minimumProportionOfMS2peaks, 
    neutralLossesPrecursorToFragments = parameterSet$neutralLossesPrecursorToFragments,
    neutralLossesFragmentsToFragments = parameterSet$neutralLossesFragmentsToFragments,
    ## matrix
    mzDeviationAbsolute_grouping = parameterSet$mzDeviationAbsolute_grouping, 
    mzDeviationInPPM_grouping = parameterSet$mzDeviationInPPM_grouping, 
    doMs2PeakGroupDeisotoping = parameterSet$doMs2PeakGroupDeisotoping, 
    mzDeviationAbsolute_ms2PeakGroupDeisotoping = parameterSet$mzDeviationAbsolute_ms2PeakGroupDeisotoping,
    mzDeviationInPPM_ms2PeakGroupDeisotoping = parameterSet$mzDeviationInPPM_ms2PeakGroupDeisotoping,
    proportionOfMatchingPeaks_ms2PeakGroupDeisotoping = parameterSet$proportionOfMatchingPeaks_ms2PeakGroupDeisotoping
  )))
  #print(paste("Hash ", paramsHash, ": ", paste(unlist(c(minimumIntensityOfMaximalMS2peak = parameterSet$minimumIntensityOfMaximalMS2peak, minimumProportionOfMS2peaks = parameterSet$minimumProportionOfMS2peaks, neutralLossesPrecursorToFragments = parameterSet$neutralLossesPrecursorToFragments, neutralLossesFragmentsToFragments = parameterSet$neutralLossesFragmentsToFragments, mzDeviationAbsolute_grouping = parameterSet$mzDeviationAbsolute_grouping, mzDeviationInPPM_grouping = parameterSet$mzDeviationInPPM_grouping, doMs2PeakGroupDeisotoping = parameterSet$doMs2PeakGroupDeisotoping, mzDeviationAbsolute_ms2PeakGroupDeisotoping = parameterSet$mzDeviationAbsolute_ms2PeakGroupDeisotoping, mzDeviationInPPM_ms2PeakGroupDeisotoping = parameterSet$mzDeviationInPPM_ms2PeakGroupDeisotoping, proportionOfMatchingPeaks_ms2PeakGroupDeisotoping = parameterSet$proportionOfMatchingPeaks_ms2PeakGroupDeisotoping)), collapse = "__"), sep = ""))
  if(paramsHash == "af46fbb2"){paramsHash <- ""}
  else                        {paramsHash <- paste("_", paramsHash, sep = "")}
  
  
  fileFragmentMatrixRData <- paste(
    substr(x = fileSpectra, start = 1, stop = regexpr(pattern = "\\.[a-zA-Z]{3,4}$", text = fileSpectra)[[1]] - 1), 
    ifelse(all(allowedInstrumentTypes=="all"), "_allowedInstrumentTypes=ALL",""), "_", 
    "fragmentMatrix", 
    paramsHash,
    ".RData", 
  sep = "")
  
  if(file.exists(fileFragmentMatrixRData) & reProcessFragmentMatrix){
    file.remove(fileFragmentMatrixRData)
  }
  
  if(file.exists(fileFragmentMatrixRData)){
    print(paste("load", fileFragmentMatrixRData))
    load(file = fileFragmentMatrixRData)
  } else {
    #returnObj <- builtMatrixOld(
    tmpParameterSet <- parameterSet[c("mzDeviationAbsolute_grouping", 
                                      "mzDeviationInPPM_grouping", 
                                      "doMs2PeakGroupDeisotoping", 
                                      "mzDeviationAbsolute_ms2PeakGroupDeisotoping",
                                      "mzDeviationInPPM_ms2PeakGroupDeisotoping",
                                      "proportionOfMatchingPeaks_ms2PeakGroupDeisotoping")]
    print(paste("Building FragmentMatrix with parameters:", paste(names(tmpParameterSet), tmpParameterSet, sep = "=", collapse = "; ")))
    returnObj <- builtMatrix(
    #returnObj <- builtMatrix_big(
      spectraList = spectraList, 
      mzDeviationAbsolute_grouping = parameterSet$mzDeviationAbsolute_grouping, 
      mzDeviationInPPM_grouping = parameterSet$mzDeviationInPPM_grouping, 
      doMs2PeakGroupDeisotoping = parameterSet$doMs2PeakGroupDeisotoping, 
      mzDeviationAbsolute_ms2PeakGroupDeisotoping = parameterSet$mzDeviationAbsolute_ms2PeakGroupDeisotoping,
      mzDeviationInPPM_ms2PeakGroupDeisotoping = parameterSet$mzDeviationInPPM_ms2PeakGroupDeisotoping,
      proportionOfMatchingPeaks_ms2PeakGroupDeisotoping = parameterSet$proportionOfMatchingPeaks_ms2PeakGroupDeisotoping, 
      progress = progress
    )
    
    print(fileFragmentMatrixRData)
    save(file = fileFragmentMatrixRData, returnObj)
  }
  #rm(spectraList)
  
  ######################################################################
  ## fragment matrix
  fragmentMatrix <- returnObj$matrix
  fragmentMasses <- returnObj$fragmentMasses
  numberOfMS2PeakGroups <- returnObj$numberOfMS2PeakGroups
  if(is.null(fragmentMatrix@Dimnames[[1]])){
    fragmentMatrix@Dimnames[[1]] <- paste(precursorMz, precursorRt, sep = " / ") #TODO: mz / rt / x for uniqueness...?
    fragmentMatrix@Dimnames[[2]] <- fragmentMasses
  }
  
  rm(returnObj)
  
  print("Spectra count for fragments")
  spectraCount_fragment <- Matrix::colSums(fragmentMatrix != 0)
  
  ## take duplicated spectrum with maximum number of peaks
  print("Select duplicated spectra")
  #if(takeSpectraWithMaximumNumberOfPeaks){
  fragmentCounts <- Matrix::rowSums(fragmentMatrix != 0)
  duplicatedSpectrumIndecesToRemove_maxPeakCount <- unlist(lapply(X = duplicatedInchiKeysBlock1, FUN = function(x){
    indeces <- which(inchiKeysBlock1 %in% x)
    peakCount <- fragmentCounts[indeces]
    indexToRetain <- which.max(peakCount)
    indecesOut <- indeces[-indexToRetain]
    return(indecesOut)
  }))
  #}
  rm(fragmentCounts)
  
#  if(FALSE){
#    cpdNames <- unlist(lapply(X = spectraList, FUN = function(x){x$name}))
#    cpdNames[grepl(pattern = "rhamnoside", x = cpdNames)]
#    inchiKeyBlock1Here <- unique(inchiKeysBlock1[grepl(pattern = "Cyanidin-3-O-rhamnoside cation", x = cpdNames)])
#    # USWXMMRFOWNEOR
#    colSumsHere <- colSums(fragmentMatrix[inchiKeysBlock1 == inchiKeyBlock1Here, ] > 0)
#    paste(colSumsHere[colSumsHere > 0], fragmentMasses[colSumsHere > 0])
#    idx <- which(duplicatedInchiKeysBlock1 == inchiKeyBlock1Here)
#  }
  
  duplicatedSpectrumIndecesToRemove_merge <- integer()
  if(mergeDuplicatedSpectra){
    ## save.image("/newhome/htreutle/Downloads/temp/MergeBug.RData")
    ## load("/newhome/htreutle/Downloads/temp/MergeBug.RData")
    
    print("Merge and remove duplicated spectra")
    start <- Sys.time()
    
    #######################################################################
    ## merge duplicated
    precursorIndecesToRemove <- vector(mode = "integer", length = 0)
    #for(idx in 1:50) {
    for(idx in seq_along(duplicatedInchiKeysBlock1)) {
    #unlist(sapply(X = duplicatedInchiKeysBlock1, FUN = function(x){
      #print(idx)
      
      if(all(length(duplicatedInchiKeysBlock1)>=10, (idx %% (as.integer(length(duplicatedInchiKeysBlock1)/10))) == 0))
        print(paste("Merge and remove duplicated spectra ", idx, " / ", length(duplicatedInchiKeysBlock1), sep = ""))
      
      indeces   <- which(inchiKeysBlock1 == duplicatedInchiKeysBlock1[[idx]])
      if(length(indeces) <= 1)  stop(indeces)
      
      ###############################
      ## perf
      #mergedSpectrum  <- apply(X = fragmentMatrix[indeces, ], MARGIN = 2, FUN = max)
      fragmentMatrixPart <- as.matrix(fragmentMatrix[indeces, ])
      
      ## take maximum fragment intensity
      mergedSpectrum  <- apply(X = fragmentMatrixPart, MARGIN = 2, FUN = max)
      ## perf
      ###############################
      
      #peakCounts      <- apply(X = fragmentMatrixPart != 0, MARGIN = 1, FUN = sum)
      #mergedSpectrum  <- fragmentMatrixPart[which.max(peakCounts), ]
      
      ms2Peaks_mz     <- fragmentMasses[mergedSpectrum != 0]
      ms2Peaks_int    <- mergedSpectrum[mergedSpectrum != 0]
      isFragment      <- ms2Peaks_mz > 0
      
      #mergedEntry <- spectraList[[indeces[[1]]]]
      #if(FALSE){
      mergedEntry <- suppressWarnings(list(
        "name"           = names(which.max(table( sapply(X = indeces, FUN = function(x){spectraList[[x]]$name          }) ))),
        "ms1Int"         = median(                sapply(X = indeces, FUN = function(x){spectraList[[x]]$ms1Int        }) ),
        "rt"             = median(                sapply(X = indeces, FUN = function(x){spectraList[[x]]$rt            }) ),
        "mz"             = median(                sapply(X = indeces, FUN = function(x){spectraList[[x]]$mz            }) ),
        "metName"        = names(table(           sapply(X = indeces, FUN = function(x){spectraList[[x]]$metName       }) ))[[1]],
        "adduct"         = names(table(           sapply(X = indeces, FUN = function(x){spectraList[[x]]$adduct        }) ))[[1]],
        "quantMass"      = median(                sapply(X = indeces, FUN = function(x){spectraList[[x]]$quantMass     }) ),
        "compoundClass"  = paste(sort(unique(unlist(strsplit(x = sapply(X = indeces, FUN = function(x){spectraList[[x]]$compoundClass}), split = "; ")))), collapse = "; "),
        "instrumentType" = names(which.max(table( sapply(X = indeces, FUN = function(x){spectraList[[x]]$instrumentType}) ))),
        "inchi"          = names(which.max(table( sapply(X = indeces, FUN = function(x){spectraList[[x]]$inchi         }) ))),
        "inchiKey"       = names(which.max(table( sapply(X = indeces, FUN = function(x){spectraList[[x]]$inchiKey      }) ))),
        "smiles"         = names(which.max(table( sapply(X = indeces, FUN = function(x){spectraList[[x]]$smiles        }) ))),
        "peakNumber"     = length(ms2Peaks_mz),
        "ms2Peaks_mz"    = ms2Peaks_mz,
        "ms2Peaks_int"   = ms2Peaks_int,
        "spectrumString" = paste(ms2Peaks_mz[isFragment], ms2Peaks_int[isFragment], sep = " ", collapse = ";")
      ))
      #}
      
      ## update
#      if(TRUE){
        representativeIdx <- indeces[[1]]
        indecesToRemove   <- indeces[-1]
        spectraList  [[representativeIdx]]  <- mergedEntry
        fragmentMatrix[representativeIdx, ] <- mergedSpectrum
#      } else {
#        indecesToRemove   <- vector(mode = "integer", length = 0)
#        for(index in indeces){
#          spectraList   [[index]] <- mergedEntry
#          fragmentMatrix[index, ] <- mergedSpectrum
#        }
#      }
      
      precursorIndecesToRemove <- c(precursorIndecesToRemove, indecesToRemove)
    }
    
#    if(FALSE){
#    fragmentMatrix <- fragmentMatrix[-precursorIndecesToRemove, ]
#    spectraList <- spectraList[-precursorIndecesToRemove]
#    
#    numberOfSpectra = numberOfSpectra - length(precursorIndecesToRemove)
#    #precursorMz     = precursorMz[-precursorIndecesToRemove]
#    #precursorRt     = precursorRt[-precursorIndecesToRemove]
#    precursorMz     = unlist(lapply(X = spectraList, FUN = function(x){x$mz}))
#    precursorRt     = unlist(lapply(X = spectraList, FUN = function(x){x$rt}))
#    inchis          = inchis[-precursorIndecesToRemove]
#    annoTable       = annoTable[-precursorIndecesToRemove, ]
#    
#    fragmentMatrix@Dimnames[[1]] <- paste(precursorMz, precursorRt, sep = " / ") #TODO: mz / rt / x for uniqueness...?
#    #spectrumId <- paste(precursorMz, precursorRt)
#    #spectrumId <- seq_len(numberOfSpectra)
#    
#    #######################################################################
#    ## duplicated
#    duplicatedStructures <- duplicatedStructures[-precursorIndecesToRemove]
#    uniqueNumberOfSpectra <- numberOfSpectra - sum(duplicatedStructures)
#    duplicatedInchiKeysBlock1 <- unique(inchiKeysBlock1[duplicatedStructures])
#    
#    duplicatedSpectrumIndecesToRemove_merge <- integer()
#    }
    duplicatedSpectrumIndecesToRemove_merge <- precursorIndecesToRemove
    
    end <- Sys.time()
    print(paste("Merge and remove duplicated spectra: ", difftime(time1 = end, time2 = start, units = "secs"), "s", sep = ""))
  }## mergeDuplicatedSpectra
  rm(inchiKeys)
  rm(inchiKeysBlock1)
  
  #######################################################################
  ## annotations
  substanceclasses   <- annoTable[, "CHEMONT_name"]
  alternativeParents <- annoTable[, "AlternativeParents_IDs"]
  substituents       <- annoTable[, "Substituents"]
  scaffolds          <- annoTable[, "Scaffolds"]
  
  substituentsList <- strsplit(x = substituents, split = "; ")
  scaffoldsList    <- strsplit(x = scaffolds,    split = "; ")
  #rm(annoTable)
  
  #######################################################################
  ## substituents
  print("Substituents")
  uniqueSubstituents_subst <- sort(unique(unlist(str_split(string = substituents, pattern = "; "))))
  if(sum(uniqueSubstituents_subst == "") > 0)
    uniqueSubstituents_subst <- uniqueSubstituents_subst[-which(uniqueSubstituents_subst == "")]
  
  ## remove C-atom numbers
  uniqueSubstituents_substTmp <- uniqueSubstituents_subst
  uniqueSubstituents_substTmp <- gsub(x = uniqueSubstituents_substTmp, pattern = "\\\\", replacement = "\\\\\\\\")
  uniqueSubstituents_substTmp <- gsub(x = uniqueSubstituents_substTmp, pattern = "\\[", replacement = "\\\\[")
  uniqueSubstituents_substTmp <- gsub(x = uniqueSubstituents_substTmp, pattern = "\\]", replacement = "\\\\]")
  uniqueSubstituents_substTmp <- gsub(x = uniqueSubstituents_substTmp, pattern = "\\[", replacement = "\\\\(")
  uniqueSubstituents_substTmp <- gsub(x = uniqueSubstituents_substTmp, pattern = "\\]", replacement = "\\\\)")
  substituentsTmp <- substituents[!duplicatedStructures]
  
  start <- Sys.time()
  numberOfSpectra_subst <- unlist(lapply(X = uniqueSubstituents_substTmp, FUN = function(x){
    sum(
      stri_startswith_fixed(pattern = paste(x, "; ", sep = ""), str = substituentsTmp) |
        stri_endswith_fixed( pattern = paste("; ", x, sep = ""), str = substituentsTmp) |
        (stri_startswith_fixed(pattern = x, str = substituentsTmp) & stri_endswith_fixed( pattern = x, str = substituentsTmp)) |
        stri_detect_fixed(pattern = paste("; ", x, "; ", sep = ""), str = substituentsTmp)
    )
  }))
  end <- Sys.time()
  print(paste("Subst: ", difftime(time1 = end, time2 = start, units = "secs"), "s", sep = ""))
  
  rm(substituentsTmp)
  
  theseSpectra_subst <- 
    (numberOfSpectra_subst >= minimumNumberOfPosSpectraPerClass) & 
    ((uniqueNumberOfSpectra - numberOfSpectra_subst) >= minimumNumberOfNegSpectraPerClass)
  substituentsWithEnoughSpectra_subst <- uniqueSubstituents_subst[theseSpectra_subst]
  numberOfSpectra_subst <- numberOfSpectra_subst[theseSpectra_subst]
  ## TODO consider for substituents the classes-path to root class analog to substance classes?
  
  ##################################
  ## alternative parents...
  allSubstanceClasses <- processChemOnt(substanceclasses = substanceclasses, alternativeParents = alternativeParents)
  
  print("Substance classes")
  start <- Sys.time()
  
  returnObj <- countSubstanceClasses(substanceclasses = substanceclasses[!duplicatedStructures], allSubstanceClasses = allSubstanceClasses[!duplicatedStructures])
  numberOfSpectra_sc  <- returnObj$"numberOfSpectra_sc"
  numberOfSpectra_asc <- returnObj$"numberOfSpectra_asc"
  substanceclassesWithSuperClasses_sc <- returnObj$"substanceclassesWithSuperClasses_sc"
  
  theseSpectra_sc <- 
    (numberOfSpectra_sc >= minimumNumberOfPosSpectraPerClass) & 
    ((uniqueNumberOfSpectra - numberOfSpectra_sc) >= minimumNumberOfNegSpectraPerClass) &
    !grepl(pattern = "^Organic compounds$", x = substanceclassesWithSuperClasses_sc)
  theseSpectra_asc <- 
    (numberOfSpectra_asc >= minimumNumberOfPosSpectraPerClass) & 
    ((uniqueNumberOfSpectra - numberOfSpectra_asc) >= minimumNumberOfNegSpectraPerClass) &
    !grepl(pattern = "^Organic compounds$", x = substanceclassesWithSuperClasses_sc)
  substanceclassesWithSuperClassesWithEnoughSpectra_sc  <- substanceclassesWithSuperClasses_sc[theseSpectra_sc ]
  substanceclassesWithSuperClassesWithEnoughSpectra_asc <- substanceclassesWithSuperClasses_sc[theseSpectra_asc]
  rm(numberOfSpectra_sc)
  rm(numberOfSpectra_asc)
  #numberOfSpectra_sc <- numberOfSpectra_sc[theseSpectra_sc]
  
  end <- Sys.time()
  print(paste("Substance classes: ", difftime(time1 = end, time2 = start, units = "secs"), "s", sep = ""))
  
  ##################################
  ## scaffolds
  print("Scaffolds")
  start <- Sys.time()
  
  scaffoldsListTmp <- scaffoldsList[!duplicatedStructures]
  uniqueScaffolds <- sort(unique(unlist(scaffoldsListTmp)))
  numberOfSpectra_scf <- sapply(X = uniqueScaffolds, FUN = function(x){
    sum(unlist(lapply(X = scaffoldsListTmp, FUN = function(y){
      any(y == x)
    })))
  })
  theseSpectra_scf <- 
    (numberOfSpectra_scf >= minimumNumberOfPosSpectraPerClass) & 
    ((uniqueNumberOfSpectra - numberOfSpectra_scf) >= minimumNumberOfNegSpectraPerClass)
  scaffoldsWithEnoughSpectra_subst  <- uniqueScaffolds[theseSpectra_scf]
  
  end <- Sys.time()
  print(paste("Scaffolds: ", difftime(time1 = end, time2 = start, units = "secs"), "s", sep = ""))
  
  ######################################################################
  ## unitResolution
  if(unitResolution){
    print("Unit resolution")
    start <- Sys.time()
    
    fragmentMassesRounded        <- round( fragmentMasses)
    duplicatedFragmentMasses_bool <- duplicated(fragmentMassesRounded)
    duplicatedFragmentMasses <- unique(fragmentMassesRounded[duplicatedFragmentMasses_bool])
    #fragmentMassesNew        <- unique(fragmentMassesRounded)
    fragmentMassesNew        <- fragmentMassesRounded[!duplicatedFragmentMasses_bool]
    
    fragmentIndecesToRemove <- vector(mode = "integer", length = 0)
    for(idx in seq_along(duplicatedFragmentMasses)){
      if(all(length(duplicatedFragmentMasses)>=10, (idx %% (as.integer(length(duplicatedFragmentMasses)/10))) == 0))
        print(paste("Unit resolution ", idx, " / ", length(duplicatedFragmentMasses), sep = ""))
      
      indeces <- which(fragmentMassesRounded == duplicatedFragmentMasses[[idx]])
      indexRepresentant <- indeces[[1]]
      indecesToRemove   <- indeces[-1]
      
      fragmentMatrix[, indexRepresentant] <- Matrix::rowSums(x = fragmentMatrix[, indeces])
      
      fragmentIndecesToRemove <- c(fragmentIndecesToRemove, indecesToRemove)
    }
    
    if(length(fragmentIndecesToRemove) > 0){
      fragmentMatrix <- fragmentMatrix[, -fragmentIndecesToRemove]
      fragmentMasses <- fragmentMassesNew
      numberOfMS2PeakGroups <- length(fragmentMasses)
      fragmentMatrix@Dimnames[[2]] <- fragmentMasses
      spectraCount_fragment <- Matrix::colSums(fragmentMatrix != 0)
    }
    
    end <- Sys.time()
    print(paste("Unit resolution: ", difftime(time1 = end, time2 = start, units = "secs"), "s", sep = ""))
  }## unitResolution
  
  ######################################################################
  ## wrap
  returnObj <- list(
    numberOfSpectra = numberOfSpectra,
    spectraList = spectraList,
    annoTable = annoTable,
    ## duplicated structures
    duplicatedStructures                           = duplicatedStructures,
    duplicatedSpectrumIndecesToRemove_maxPeakCount = duplicatedSpectrumIndecesToRemove_maxPeakCount,
    duplicatedSpectrumIndecesToRemove_merge        = duplicatedSpectrumIndecesToRemove_merge,
    numberOfSpectraUnique                          = uniqueNumberOfSpectra,
    ## matrix
    fragmentMatrix        = fragmentMatrix,
    fragmentMasses        = fragmentMasses,
    numberOfMS2PeakGroups = numberOfMS2PeakGroups,
    spectraCount_fragment = spectraCount_fragment,
    ## substance classes
    substanceclasses                                      = substanceclasses,
    substanceclassesWithSuperClassesWithEnoughSpectra_sc  = substanceclassesWithSuperClassesWithEnoughSpectra_sc,
    ## all substance classes
    allSubstanceClasses                                   = allSubstanceClasses,
    substanceclassesWithSuperClassesWithEnoughSpectra_asc = substanceclassesWithSuperClassesWithEnoughSpectra_asc,
    ## substituents
    substituents                        = substituentsList,
    substituentsWithEnoughSpectra_subst = substituentsWithEnoughSpectra_subst,
    ## scaffolds
    scaffolds                        = scaffoldsList,
    scaffoldsWithEnoughSpectra_subst = scaffoldsWithEnoughSpectra_subst
  )
  return(returnObj)
}## processLibrary
processChemOnt <- function(substanceclasses, alternativeParents){
  ##################################
  ## alternative parents... see OBO.R
  ## iterate substance classes with path to taxonomy root
  library("ontologyIndex")
  file <- "./data/ChemOnt_2_1.obo"
  oboIndex <- get_OBO(file = file)
  
  #alternativeParents <- "CHEMONTID:0003458; CHEMONTID:0004557; CHEMONTID:0000323; CHEMONTID:0003940; CHEMONTID:0000469; CHEMONTID:0004150"
  
  alternativeParentIDsLists  <- strsplit(x = alternativeParents, split = "; ")
  alternativeParentIDsUnique <- unique(unlist(alternativeParentIDsLists))
  alternativeParentPathNamesUnique <- unlist(lapply(X = alternativeParentIDsUnique, FUN = function(y){
    ancestorPathIDs   <- get_ancestors(ontology = oboIndex, terms = y)
    ancestorPathNames <- unlist(lapply(X = ancestorPathIDs, FUN = function(x){
      get_term_property(ontology = oboIndex, property_name = "name", term = x)
    }))
    ancestorPathNames <- paste(ancestorPathNames[2:length(ancestorPathNames)], collapse = "; ")
    return(ancestorPathNames)
  }))
  alternativeParentPathNamesLists <- lapply(X = alternativeParentIDsLists, function(x){
    alternativeParentPathNames <- alternativeParentPathNamesUnique[alternativeParentIDsUnique %in% x]
  })
  
  allSubstanceClasses <- mapply(c, as.list(substanceclasses), alternativeParentPathNamesLists, SIMPLIFY=FALSE)
  rm(oboIndex)
  
  return(allSubstanceClasses)
}
countSubstanceClasses <- function(substanceclasses, allSubstanceClasses){
  uniqueSubstanceClasses <- unique(unlist(allSubstanceClasses))
  
  substanceclassesWithSuperClasses_sc <- sort(unique(unlist(
    lapply(X = strsplit(x = uniqueSubstanceClasses, split = "; "), FUN = function(x){
      x <- x[x != "NA"]
      sapply(X = seq_len(length(x)), FUN = function(y){
        paste(x[1:y], collapse = "; ")
      })
    })
  )))
  
  substanceclassesToIgnore <- c(
    "Organic compounds"
  )
  substanceclassesWithSuperClasses_sc <- substanceclassesWithSuperClasses_sc[-(substanceclassesWithSuperClasses_sc %in% substanceclassesToIgnore)]
  
  substanceclassesWithSuperClasses_sc2 <- substanceclassesWithSuperClasses_sc
  substanceclassesWithSuperClasses_sc2 <- gsub(x = substanceclassesWithSuperClasses_sc2, pattern = "\\\\", replacement = "\\\\\\\\")
  substanceclassesWithSuperClasses_sc2 <- gsub(x = substanceclassesWithSuperClasses_sc2, pattern = "\\[", replacement = "\\\\[")
  substanceclassesWithSuperClasses_sc2 <- gsub(x = substanceclassesWithSuperClasses_sc2, pattern = "\\]", replacement = "\\\\]")
  substanceclassesWithSuperClasses_sc2 <- gsub(x = substanceclassesWithSuperClasses_sc2, pattern = "\\(", replacement = "\\\\(")
  substanceclassesWithSuperClasses_sc2 <- gsub(x = substanceclassesWithSuperClasses_sc2, pattern = "\\)", replacement = "\\\\)")
  
  ## sort out substance classes with too few spectra
  #substanceclassesTmp <- substanceclasses#[!duplicatedStructures]
  #alternativeParentPathNamesListsTmp <- allSubstanceClasses#[!duplicatedStructures]
  numberOfSpectra_sc <- unlist(lapply(X = substanceclassesWithSuperClasses_sc2, FUN = function(x){
    #sum(grepl(pattern = paste("^", x, sep = ""), x = substanceclassesTmp))
    directParent      <- stri_startswith_fixed(pattern = x, str = substanceclasses)
    return(sum(directParent))
  }))
  numberOfSpectra_asc <- unlist(lapply(X = substanceclassesWithSuperClasses_sc2, FUN = function(x){
    #directParent      <- stri_startswith_fixed(pattern = x, str = substanceclasses)
    alternativeParent <- unlist(lapply(X = allSubstanceClasses, FUN = function(y){any(stri_startswith_fixed(pattern = x, str = y))}))
    #return(sum(directParent | alternativeParent))
    return(sum(alternativeParent))
  }))
  return(list("numberOfSpectra_sc" = numberOfSpectra_sc, "numberOfSpectra_asc" = numberOfSpectra_asc, "substanceclassesWithSuperClasses_sc" = substanceclassesWithSuperClasses_sc))
}
processSpectraAndAnnotation <- function(fileSpectra, parameterSet, allowedInstrumentTypes, annoFile, structureFormats, reParseMsp, reProcessAnno, progress){
  #source("/home/htreutle/Code/Java/MetFam_util/SubstanceClassClassifier.R")
  
  #######################################################################
  ## parse msp
  paramsHash <- digest::sha1(algo = "crc32", x = unlist(c(
    minimumIntensityOfMaximalMS2peak = parameterSet$minimumIntensityOfMaximalMS2peak, 
    minimumProportionOfMS2peaks = parameterSet$minimumProportionOfMS2peaks, 
    neutralLossesPrecursorToFragments = parameterSet$neutralLossesPrecursorToFragments,
    neutralLossesFragmentsToFragments = parameterSet$neutralLossesFragmentsToFragments
  )))
  #if(paramsHash == "da3c656b") paramsHash <- ""
  #else                         
  paramsHash <- paste("_", paramsHash, sep = "")
  
  
  fileSpectraRData <- paste(substr(x = fileSpectra, start = 1, stop = regexpr(pattern = "\\.[a-zA-Z]{3,4}$", text = fileSpectra)[[1]] - 1), paramsHash, ".RData", sep = "")
  
  if(file.exists(fileSpectraRData) & reParseMsp){
    file.remove(fileSpectraRData)
  }
  
  if(file.exists(fileSpectraRData)){
    print(paste("load", fileSpectraRData))
    load(file = fileSpectraRData)
  } else {
    #returnObj <- parseMSP(
    returnObj <- parseMSP_big(
      fileSpectra = fileSpectra, 
      minimumIntensityOfMaximalMS2peak = parameterSet$minimumIntensityOfMaximalMS2peak, 
      minimumProportionOfMS2peaks = parameterSet$minimumProportionOfMS2peaks, 
      neutralLossesPrecursorToFragments = parameterSet$neutralLossesPrecursorToFragments,
      neutralLossesFragmentsToFragments = parameterSet$neutralLossesFragmentsToFragments,
      progress = progress
    )
    print(fileSpectraRData)
    save(file = fileSpectraRData, returnObj)
  }
  
  spectraList     <- returnObj$spectraList
  numberOfSpectra <- returnObj$numberOfSpectra
  precursorMz     <- returnObj$precursorMz
  precursorRt     <- returnObj$precursorRt
  rm(returnObj)
  
  #######################################################################
  ## select instruments (e.g. sort out low resolution)
#  if (FALSE) {
#  if(all(length(allowedInstrumentTypes) == 1, all(allowedInstrumentTypes == "all"))) {
#    isAllowedInstrumentType <- rep(x = TRUE, times = length(spectraList))
#  }else{
#    isAllowedInstrumentType <- unlist(lapply(X = spectraList, FUN = function(x){x$instrument %in% allowedInstrumentTypes}))
#  }
#  print(paste(sum(!isAllowedInstrumentType), "/", numberOfSpectra, "spectra arise from other instruments,", sum(isAllowedInstrumentType), "remain"))
#  
#  spectraList     <- spectraList[isAllowedInstrumentType]
#  numberOfSpectra <- length(spectraList)
#  precursorMz     <- precursorMz[isAllowedInstrumentType]
#  precursorRt     <- precursorRt[isAllowedInstrumentType]
#  }
  #######################################################################
  ## sort out unknown inchis
  noInchi    <- unlist(lapply(X = spectraList, FUN = function(x){x$inchi == ""}))
  noInchiKey <- unlist(lapply(X = spectraList, FUN = function(x){x$inchiKey == ""}))
  noSmiles   <- unlist(lapply(X = spectraList, FUN = function(x){x$smiles == ""}))
  noStructure <- noInchi & noInchiKey & noSmiles
  #print(paste(sum(noInchi), "/", numberOfSpectra, "spectra have no InChI,", sum(!noInchi), "remain"))
  print(paste(sum(noStructure), "/", numberOfSpectra, "spectra have no structure information,", sum(!noStructure), "remain"))
  
  #spectraList     <- spectraList[!noInchi]
  #numberOfSpectra <- length(spectraList)
  #precursorMz     <- precursorMz[!noInchi]
  #precursorRt     <- precursorRt[!noInchi]
  spectraList     <- spectraList[!noStructure]
  numberOfSpectra <- length(spectraList)
  precursorMz     <- precursorMz[!noStructure]
  precursorRt     <- precursorRt[!noStructure]
  
  inchis    <- unlist(lapply(X = spectraList, FUN = function(x){x$inchi}))
  inchiKey  <- unlist(lapply(X = spectraList, FUN = function(x){x$inchiKey}))
  smiles    <- unlist(lapply(X = spectraList, FUN = function(x){x$smiles}))
  
  numberOfInchis    <- sum(inchis   != "" & inchis   != "NA")
  numberOfInchiKeys <- sum(inchiKey != "" & inchiKey != "NA")
  numberOfSmiles    <- sum(smiles   != "" & smiles   != "NA")
  
  usedStructureFormat <- NULL
  structures <- NULL
  if(numberOfSmiles == max(c(numberOfInchis, numberOfInchiKeys, numberOfSmiles))){
    usedStructureFormat <- structureFormats$SMILES
    structures <- smiles
  }
  if(numberOfInchiKeys == max(c(numberOfInchis, numberOfInchiKeys, numberOfSmiles))){
    usedStructureFormat <- structureFormats$InChIKey
    structures <- inchiKey
  }
  if(numberOfInchis == max(c(numberOfInchis, numberOfInchiKeys, numberOfSmiles))){
    usedStructureFormat <- structureFormats$InChI
    structures <- inchis
  }
  
  #######################################################################
  ## annotate substance class
  annoFileName <- gsub(x = basename(annoFile), pattern = "\\.tsv$", replacement = "")
  #fileAnnoTable <- paste(substr(x = fileSpectra, start = 1, stop = regexpr(pattern = "\\.[a-zA-Z]{3,4}$", text = fileSpectra)[[1]] - 1), "_annotationTable.tsv", sep = "")
  fileAnnoTable <- paste(substr(x = fileSpectra, start = 1, stop = regexpr(pattern = "\\.[a-zA-Z]{3,4}$", text = fileSpectra)[[1]] - 1), "_annotationTable_", annoFileName, paramsHash, ".tsv", sep = "")
  
  if(file.exists(fileAnnoTable) & reProcessAnno){
    file.remove(fileAnnoTable)
  }
  
  if(file.exists(fileAnnoTable)){
    annoTable <- read.table(file = fileAnnoTable, sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
  } else {
    annoTable <- read.table(file = annoFile,      sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
  }
  
  if(!(usedStructureFormat %in% colnames(annoTable)))
    stop(paste(usedStructureFormat, "not supported"))
  
  ## sort out spectra without annotation entry
#  if(FALSE){
#    structureWithEntry <- structures %in% annoTable[, usedStructureFormat]
#    print(paste(sum(!structureWithEntry), "/", numberOfSpectra, "spectra have no annotation entry,", sum(structureWithEntry), "remain"))
#    
#    spectraList     <- spectraList     [structureWithEntry]
#    numberOfSpectra <- length(spectraList)
#    precursorMz     <- precursorMz     [structureWithEntry]
#    precursorRt     <- precursorRt     [structureWithEntry]
#    inchis          <- inchis          [structureWithEntry]
#    structures      <- structures      [structureWithEntry]
#    
#    ## map annotations to spectra
#    mapping <- match(x = structures, table = annoTable[, usedStructureFormat])
#    annoTable <- annoTable[mapping, ]
#    annoTable[is.na(annoTable)] <- ""
#  } else {
    switch(usedStructureFormat,
           "InChI"={
             structuresPrefixes      <- inchisToTruncatedInchis(structures)
             #structuresPrefixesThere <- inchisToTruncatedInchis(annoTable[, usedStructureFormat])
             structuresPrefixesThere <- structuresPrefixes
             mapping <- match(x = structuresPrefixes, table = structuresPrefixesThere)
             structureWithEntry <- !is.na(mapping)
             
             print(paste(sum(!structureWithEntry), "/", numberOfSpectra, "spectra have no annotation entry,", sum(structureWithEntry), "remain"))
             
             spectraList     <- spectraList     [structureWithEntry]
             numberOfSpectra <- length(spectraList)
             precursorMz     <- precursorMz     [structureWithEntry]
             precursorRt     <- precursorRt     [structureWithEntry]
             inchis          <- inchis          [structureWithEntry]
             structures      <- structures      [structureWithEntry]
             
             ## map annotations to spectra
             #mapping <- match(x = structures, table = annoTable[, usedStructureFormat])
             #structuresPrefixes <- unlist(lapply(X = strsplit(x = structures, split = "-"), FUN = function(x){x[[1]]}))
             #mapping <- match(x = structuresPrefixes, table = structuresPrefixesThere)
             mapping <- mapping      [structureWithEntry]
             
             annoTable <- annoTable[mapping, ]
             annoTable[is.na(annoTable)] <- ""
           },
           "InChIKey"={
             structuresPrefixes <- inchiKeysToInchiKeysBlock1(structures)
             structuresPrefixesThere <- inchiKeysToInchiKeysBlock1(annoTable[, usedStructureFormat])
             mapping <- match(x = structuresPrefixes, table = structuresPrefixesThere)
             structureWithEntry <- !is.na(mapping)
             
             print(paste(sum(!structureWithEntry), "/", numberOfSpectra, "spectra have no annotation entry,", sum(structureWithEntry), "remain"))
             
             spectraList     <- spectraList     [structureWithEntry]
             numberOfSpectra <- length(spectraList)
             precursorMz     <- precursorMz     [structureWithEntry]
             precursorRt     <- precursorRt     [structureWithEntry]
             inchis          <- inchis          [structureWithEntry]
             structures      <- structures      [structureWithEntry]
             
             ## map annotations to spectra
             #mapping <- match(x = structures, table = annoTable[, usedStructureFormat])
             #structuresPrefixes <- unlist(lapply(X = strsplit(x = structures, split = "-"), FUN = function(x){x[[1]]}))
             #mapping <- match(x = structuresPrefixes, table = structuresPrefixesThere)
             mapping <- mapping      [structureWithEntry]
             
             annoTable <- annoTable[mapping, ]
             annoTable[is.na(annoTable)] <- ""
           },
           "SMILES"={
             ## business as usual
             structureWithEntry <- structures %in% annoTable[, usedStructureFormat]
             print(paste(sum(!structureWithEntry), "/", numberOfSpectra, "spectra have no annotation entry,", sum(structureWithEntry), "remain"))
             
             spectraList     <- spectraList     [structureWithEntry]
             numberOfSpectra <- length(spectraList)
             precursorMz     <- precursorMz     [structureWithEntry]
             precursorRt     <- precursorRt     [structureWithEntry]
             inchis          <- inchis          [structureWithEntry]
             structures      <- structures      [structureWithEntry]
             
             ## map annotations to spectra
             mapping <- match(x = structures, table = annoTable[, usedStructureFormat])
             annoTable <- annoTable[mapping, ]
             annoTable[is.na(annoTable)] <- ""
           },
           stop(paste("Unknown structure format (", usedStructureFormat, ")!", sep = ""))
    )
#  }
  
  ## inchis with anno
  #annoTable[is.na(annoTable[, "CHEMONT_name"]), "CHEMONT_name"] <- "NA"
  structureWithAnno <- !(annoTable[, "CHEMONT_name"] == "NA" | annoTable[, "CHEMONT_name"] == "")
  print(paste(sum(!structureWithAnno, na.rm = TRUE) + sum(is.na(structureWithAnno)), "/", numberOfSpectra, "spectra have no annotation,", sum(structureWithAnno), "remain"))
  
  spectraList     <- spectraList[structureWithAnno]
  numberOfSpectra <- length(spectraList)
  precursorMz     <- precursorMz[structureWithAnno]
  precursorRt     <- precursorRt[structureWithAnno]
  inchis          <- inchis     [structureWithAnno]
  annoTable       <- annoTable  [structureWithAnno, ]
  
  ## inchis with inchiKey
  inchiKeys          <- annoTable[, "InChIKey"]
  structureWithInchiKey <- inchiKeys != ""
  print(paste(sum(!structureWithInchiKey, na.rm = TRUE) + sum(is.na(structureWithInchiKey)), "/", numberOfSpectra, "spectra have no inchiKey,", sum(structureWithInchiKey), "remain"))
  
  spectraList     <- spectraList[structureWithInchiKey]
  numberOfSpectra <- length(spectraList)
  precursorMz     <- precursorMz[structureWithInchiKey]
  precursorRt     <- precursorRt[structureWithInchiKey]
  inchis          <- inchis     [structureWithInchiKey]
  annoTable       <- annoTable  [structureWithInchiKey, ]
  
  if(!file.exists(fileAnnoTable)){
    write.table(x = annoTable, file = fileAnnoTable, sep = "\t", row.names = FALSE)
  }
  
  returnObj <- list(
    spectraList     = spectraList,
    numberOfSpectra = numberOfSpectra,
    precursorMz     = precursorMz,
    precursorRt     = precursorRt,
    inchis          = inchis,
    annoTable       = annoTable
  )
  return(returnObj)
}

##############################################################################################################################################
## util
inchiKeysToInchiKeysBlock1 <- function(inchiKeys){
  inchiKeysSplitted <- strsplit(x = inchiKeys, split = "-")
  inchiKeysBlock1 <- unlist(lapply(X = inchiKeysSplitted, FUN = function(x){
    if(length(x) == 0){
      return("")
    } else {
      return(x[[1]])
    }
  }))
  return(inchiKeysBlock1)
}
inchisToTruncatedInchis <- function(inchis){
  ## remove inchi parts: version / formula / c / h / p / q / b / t/m / s / i/h / f / r
	#print(inchis)
	#print(typeof(inchis))
	#if (is.null(inchis)) return("C")
  inchisSplitted <- strsplit(x = inchis, split = "/")
  structurePrefixes <- unlist(lapply(X = inchisSplitted, FUN = function(x){
    toTruncate <- grepl(pattern = "^[btmsifr].+", x = x)
    if(any(toTruncate)){
      x <- x[1:(min(which(toTruncate)) - 1)]
    }
    return(paste(x, collapse = "/"))
  }))
  return(structurePrefixes)
}
selectColumnsForClassification <- function(matrix_train, classes_pm_train, fragmentColumnSelectionMethod, minimumProportionOfPositiveFragments){
  ## selecting columns
  posRows <- classes_pm_train=="+"
  colSumPos <- Matrix::colSums(x = matrix_train[ posRows, ])
  #colSumNeg <- Matrix::colSums(x = matrix_train[!posRows, ])
  colFrequencyPos <- colSumPos / sum( posRows)
  
  switch(fragmentColumnSelectionMethod,
         "AbsoluteProportion"={
           selectedColumns <- colFrequencyPos >= minimumProportionOfPositiveFragments
         },
         "ProportionOfSumOfFragmentPresence"={
           ## proportionOfSumOfFragmentPresence
           proportionOfSumOfFragmentPresence <- minimumProportionOfPositiveFragments#0.5 ## high means less fragments
           tmp1 <- sort(colFrequencyPos)
           tmp2 <- cumsum(tmp1)
           idx1 <- min(which(tmp2>sum(colFrequencyPos)*(1-proportionOfSumOfFragmentPresence)))
           idx2 <- which(names(colFrequencyPos)==names(tmp2)[[idx1]])
           threshValue <- colFrequencyPos[[idx2]]
           selectedColumns <- colFrequencyPos >= threshValue
         },
         "ProportionOfHighestFragmentPresence"={
           ## proportionOfHighestFragmentPresence
           proportionOfHighestFragmentPresence <- minimumProportionOfPositiveFragments#0.05
           threshValue <- max(colFrequencyPos) * proportionOfHighestFragmentPresence
           selectedColumns <- colFrequencyPos >= threshValue
         },
         stop(paste("Unknown selection method (", columnSelectionMethod, ")!", sep = ""))
  )
  
  return(selectedColumns)
}
fragmentStatistics <- function(matrix, classes){
  ## make binary
  #matrix@x[matrix@x != 0] <- 1
  
  ## remove MS/MS-precursor-m/z - MS-precursor-m/z
  #matrix <- matrix[, abs(colnames(matrix)) >= 1]
  ## pos (and neg) statistics for frequent fragments
  posRows <- classes=="+"
  colSums_pos <- Matrix::colSums(x = matrix[ posRows, ] > 0) / sum( posRows)
  colSums_neg <- Matrix::colSums(x = matrix[!posRows, ] > 0) / sum(!posRows)
  ## pos/neg statistics for characteristic fragments
  colSums_posNeg <- colSums_pos - colSums_neg
  ## annotate fragment masses
  names(colSums_pos)    <- colnames(matrix)
  names(colSums_posNeg) <- colnames(matrix)
  ## select
  minimumFrequency <- 0.01
  colSums_pos    <- rev(sort(colSums_pos   ))
  colSums_posNeg <- rev(sort(colSums_posNeg))
  colSums_pos    <- colSums_pos   [colSums_pos    >= minimumFrequency]
  colSums_posNeg <- colSums_posNeg[colSums_posNeg >= minimumFrequency]
  
  ## mean fragment intensity
  meanFragmentIntensity    <- Matrix::colSums(x = matrix            ) / Matrix::colSums(x = matrix             > 0)
  meanFragmentIntensityPos <- Matrix::colSums(x = matrix[ posRows, ]) / Matrix::colSums(x = matrix[ posRows, ] > 0)
  meanFragmentIntensityNeg <- Matrix::colSums(x = matrix[!posRows, ]) / Matrix::colSums(x = matrix[!posRows, ] > 0)
  
  columnIndeces_pos    <- match(x = names(colSums_pos), table = colnames(matrix))
  meanIntensity_pos    <- meanFragmentIntensity[columnIndeces_pos]
  columnIndeces_posNeg <- match(x = names(colSums_posNeg), table = colnames(matrix))
  meanIntensity_posNeg <- meanFragmentIntensity[columnIndeces_posNeg]
  
  ## wrap
  returnObj <- list()
  returnObj$minimumFrequency        <- minimumFrequency
  returnObj$frequentFragments       <- colSums_pos
  returnObj$characteristicFragments <- colSums_posNeg
  returnObj$frequentFragmentsMeanIntensity       <- meanIntensity_pos
  returnObj$characteristicFragmentsMeanIntensity <- meanIntensity_posNeg
  
  return(returnObj)
}

##############################################################################################
evaluatePerformance <- function(predicted_scores, classes_pm_test, spectrumId_test){
  ## remove NA and replace Inf stuff
  isOut <- which(is.na(predicted_scores))
  isInf <- which(is.infinite(predicted_scores))
  
  if(length(isInf) > 0){
    ## max plus 1
    predicted_scores[isInf] <- max(predicted_scores[-isInf], na.rm = TRUE) + 1
  }
  if(length(isOut) > 0){
    predicted_scores <- predicted_scores[-isOut]
    classes_pm_test  <- classes_pm_test [-isOut]
  }
  
  ## auc
  pred <- prediction(predictions = predicted_scores, labels = classes_pm_test)
  perf <- performance(prediction.obj = pred, measure = "tpr", x.measure = "fpr")
  auc_roc  <- performance(prediction.obj = pred, measure = "auc")@y.values[[1]]
  
  ## area under precision recall curve
  auc_pr   <- PRROC::pr.curve(scores.class0 = predicted_scores[classes_pm_test=="+"], scores.class1 = predicted_scores[classes_pm_test=="-"])#, curve = T)
  auc_pr   <- auc_pr$auc.integral
  #auc_roc2 <- PRROC::roc.curve(scores.class0 = predicted_scores[classes_pm_test=="+"], scores.class1 = predicted_scores[classes_pm_test=="-"])#, curve = T)
  
#  if(FALSE){
#    
#    ## auc
#    plot(perf, main=auc_roc)
#    segments(0,1,1,1, lty = 2)
#    segments(1,0,1,1, lty = 2)
#    
#    ## scores
#    ## plot cumulative scores distribution
#    posScores <- sort(predicted_scores[classes_pm_test=="+"])
#    negScores <- sort(predicted_scores[classes_pm_test=="-"])
#    plot(NA, xlim=range(c(posScores, negScores)), ylim=c(0,1))
#    lines(x = negScores, y = 1:length(negScores) / length(negScores), col="red")
#    lines(x = posScores, y = 1:length(posScores) / length(posScores), col="blue")
#    legend(x = 0.7, y = 0.15, legend = c("neg", "pos"), lwd = 1, col = c("red", "blue"))
#    #legend(x = min(predicted_scores), y = 1, legend = c("neg", "pos"), lwd = 1, col = c("red", "blue"))
#    
#    ## Histograms absolute
#    hist(negScores, col=rgb(1,0,0,0.5), breaks = 20, xlim=range(c(posScores, negScores)), ylim=c(0,max(hist(negScores, breaks = 20, plot = F)$counts)), xlab="score", ylab="absolute frequency")
#    hist(posScores, col=rgb(0,0,1,0.5), breaks = 20, add=T)
#    box()
#    
#    ## Histograms relative
#    h1 <- hist(negScores, breaks = 20, plot = F)
#    h2 <- hist(posScores, breaks = 20, plot = F)
#    h1$counts <- h1$counts / sum(h1$counts)
#    h2$counts <- h2$counts / sum(h2$counts)
#    
#    plot(NA, xlim=range(c(posScores, negScores)), ylim=c(0,max(c(h1$counts, h2$counts))), xlab="score", ylab="relative frequency")
#    plot(h1, col=rgb(1,0,0,0.5), add = T)
#    plot(h2, col=rgb(0,0,1,0.5), add = T)
#  }
  
  ## TPR for FPR = 5%
  perf <- performance(prediction.obj = pred, measure = "tpr", x.measure = "fpr")
  tpr <- perf@y.values[[1]]
  fpr <- perf@x.values[[1]]
  
  order <- order(fpr)
  tpr <- tpr[order]
  fpr <- fpr[order]
  idx <- which(fpr <= 0.05)
  tpr <- tpr[[max(idx)]]
  
  ## TNR for FNR = 5%
  perf <- performance(prediction.obj = pred, measure = "tnr", x.measure = "fnr")
  tnr <- perf@y.values[[1]]
  fnr <- perf@x.values[[1]]
  
  order <- order(fnr)
  tnr <- tnr[order]
  fnr <- fnr[order]
  idx <- which(fnr <= 0.05)
  tnr <- tnr[[max(idx)]]
  
  ## Sn and Sp for Max sum Sn+Sp
  perf <- performance(prediction.obj = pred, measure = "tpr", x.measure = "tnr")
  sn <- perf@y.values[[1]]
  sp <- perf@x.values[[1]]
  sumSnSp <- sn + sp
  idx <- which.max(sumSnSp)
  
  max_Sn <- sn[[idx]]
  max_Sp <- sp[[idx]]
  
  predictedScores_pos <- predicted_scores[classes_pm_test == "+"]
  predictedScores_neg <- predicted_scores[classes_pm_test == "-"]
  spectrumId_test_pos <- spectrumId_test [classes_pm_test == "+"]
  spectrumId_test_neg <- spectrumId_test [classes_pm_test == "-"]
  
  returnObj <- list()
  returnObj$auc_roc <- auc_roc
  returnObj$auc_pr  <- auc_pr
  returnObj$tpr     <- tpr
  returnObj$max_Sn  <- max_Sn
  returnObj$max_Sp  <- max_Sp
  returnObj$tnr     <- tnr
  returnObj$predictedScores_pos <- predictedScores_pos
  returnObj$predictedScores_neg <- predictedScores_neg
  returnObj$spectrumId_test_pos <- spectrumId_test_pos
  returnObj$spectrumId_test_neg <- spectrumId_test_neg
  
  return(returnObj)
}
