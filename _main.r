# ---------- Preparations ----------
# Load libraries
library(SparseM)
library(slam)
library(matrixStats)
library(Matrix)
library(plotrix)
library(MASS)
library(RMassBank)
library(rinchi)
library(stringr)
library(stringi)
library(ROCR)
library(PRROC)
library(klaR)
library(e1071)
library(kohonen)
library(nnet)
library(rda)
library(caret)
library(caretEnsemble)
library(squash)
library(RCurl)
library(ontologyIndex)
library(tools)


# Set encoding of files and set locale
options(encoding="UTF-8")
Sys.setlocale(category="LC_ALL", locale="C")



# ---------- Global variables ----------
# Working directory
#setwd("/usr/local/share/mFam/")
setwd("~/Desktop/Projekte/Habilitation/de.NBI/Classifier_mFam/")



# ---------- Load MetFamily ----------
# Load MetFamily
source("MetFamily/FragmentMatrixFunctions.R")
source("MetFamily/Annotation.R")
source("MetFamily/DataProcessing.R")

# Load MFam Classifier
source("MFam/SubstanceClassClassifier.R")
source("MFam/SubstanceClassClassifier_classifier.R")



# ---------- MetFamily and MFam Parameter ----------
# Cache
reParseMsp              <- FALSE
reProcessAnno           <- FALSE
reProcessFragmentMatrix <- FALSE
reProcessLibrary        <- FALSE

# I/O
resultFolderForClassifiers       <- "./data/Classifier_ROC_Analysis"
resultFolderForMetfamilyProjects <- "./data/MetFamily_class_projects"

progress <- FALSE
outDetail <- FALSE
maximumNumberOfScores <- 10000
removeRareFragments <- FALSE

builtClassifiers       <- TRUE
writeMetFamilyProjects <- FALSE ####TRUE
calculateRMDstatistics <- FALSE
computeFragmentFisherTest <- FALSE

writeClassifiers <- TRUE
writeResults     <- TRUE

outSuffix <- ""

# Replicates
mergeDuplicatedSpectra              <- TRUE
takeSpectraWithMaximumNumberOfPeaks <- FALSE

# Classes
classOfClass <- c(
	"ChemOnt|MainClass",   # 1
	"ChemOnt|AllClasses",  # 2
	"ChemOnt|Substituent", # 3
	"ChEBI|Identified",    # 4
	"ChEBI|Predicted",     # 5
	"Scaffold"             # 6
)[[2]]
#thisClass <- "Organic compounds; Lipids and lipid-like molecules; Prenol lipids; Terpene lactones; Sesquiterpene lactones"
thisClass <- NULL

# Repeated random subsampling validation
minimumNumberOfPosSpectraPerClass <- 10
minimumNumberOfNegSpectraPerClass <- 10
numberOfDataSetDecompositions <- 10
proportionTraining <- 0.7
#fragmentColumnSelectionMethod <- "AbsoluteProportion"
#fragmentColumnSelectionMethod <- "ProportionOfSumOfFragmentPresence"
fragmentColumnSelectionMethod <- "ProportionOfHighestFragmentPresence"
minimumProportionOfPositiveFragments <- 0.05

# Postprocessing
unitResolution <- FALSE

# Constraints
# Library one
#thisLibrary    <- NULL
#thisLibrary    <- "./data/180912_MoNA-export-LC-MS-MS_Negative_Mode_processed.msp"
#thisLibrary <- "./data/2018-02-13_09_14_10_neg_11328_MoNA-export-LC-MS-MS_Spectra.msp"
thisLibrary <- "./data/2018-02-13_09_14_10_pos_21908_MoNA-export-LC-MS-MS_Spectra.msp"
annoFile <- "./data/190523_MSMS_HR_someScaffolds_completed.tsv"

splitClasses <- FALSE
splitClassesParameterSet <- list(
	minimumNumberOfChildren = 5,
	maximumNumberOfPrecursors = 1000,
	distanceMeasure = "Jaccard (intensity-weighted)",
	#distanceMeasure = "Jaccard",
	clusterMethod = "ward.D"
)

thisMethod <- "method=ColSumsPos; smoothIntensities=FALSE"
#thisMethod <- "method=caret; smoothIntensities=FALSE, modelName=binda"

thisParameterSet <- list(
	## parse msp
	minimumIntensityOfMaximalMS2peak                  = 000,
	minimumProportionOfMS2peaks                       = 0.05,
	neutralLossesPrecursorToFragments                 = TRUE,
	neutralLossesFragmentsToFragments                 = FALSE,
	## built fragment matrix
	mzDeviationAbsolute_grouping                      = 0.01,
	mzDeviationInPPM_grouping                         = 10,
	#doMs2PeakGroupDeisotoping                         = TRUE,
	doMs2PeakGroupDeisotoping                         = FALSE,
	mzDeviationAbsolute_ms2PeakGroupDeisotoping       = 0.01,
	mzDeviationInPPM_ms2PeakGroupDeisotoping          = 10,
	proportionOfMatchingPeaks_ms2PeakGroupDeisotoping = 0.9
)

# Box parameters
parameterSetAll <- list(
	# I/O
	progress              = progress,
	outDetail             = outDetail,
	resultFolderForClassifiers          = resultFolderForClassifiers,
	resultFolderForMetfamilyProjects = resultFolderForMetfamilyProjects,
	writeClassifiers      = writeClassifiers,
	writeResults          = writeResults,
	outSuffix             = outSuffix,
	maximumNumberOfScores = maximumNumberOfScores,
	removeRareFragments   = removeRareFragments,
	builtClassifiers      = builtClassifiers,
	writeMetFamilyProjects = writeMetFamilyProjects,
	calculateRMDstatistics = calculateRMDstatistics,
	
	# Classes
	mergeDuplicatedSpectra              = mergeDuplicatedSpectra,
	takeSpectraWithMaximumNumberOfPeaks = takeSpectraWithMaximumNumberOfPeaks,
	classOfClass = classOfClass,
	thisClass = thisClass,
	splitClasses = splitClasses,
	splitClassesParameterSet = splitClassesParameterSet,
	##computeFragmentFisherTest         = computeFragmentFisherTest,
	
	# Repeated random subsampling validation
	minimumNumberOfPosSpectraPerClass = minimumNumberOfPosSpectraPerClass,
	minimumNumberOfNegSpectraPerClass = minimumNumberOfNegSpectraPerClass,
	numberOfDataSetDecompositions     = numberOfDataSetDecompositions,
	proportionTraining                = proportionTraining,
	fragmentColumnSelectionMethod     = fragmentColumnSelectionMethod,
	minimumProportionOfPositiveFragments = minimumProportionOfPositiveFragments,
	
	# Postprocessing
	unitResolution                    = unitResolution,
	
	# Cache
	reParseMsp = reParseMsp,
	reProcessAnno = reProcessAnno,
	reProcessLibrary = reProcessLibrary,
	reProcessFragmentMatrix = reProcessFragmentMatrix,
	
	# Constraints
	# Library one
	annoFile = annoFile,
	thisLibrary = thisLibrary,
	thisParameterSet = thisParameterSet,
	thisMethod  = thisMethod
)



# ---------- MFam Classifier ----------
runTest(parameterSetAll = parameterSetAll)

