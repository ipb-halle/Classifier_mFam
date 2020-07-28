#!/usr/bin/env Rscript

# Setup R error handling to go to stderr
options(show.error.messages=F, error=function() { cat(geterrmessage(), file=stderr()); q("no",1,F) } )

# Set locales and encoding
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
loc <- Sys.setlocale(category="LC_ALL", locale="C")
options(encoding="UTF-8")

# Set options
options(stringAsfactors=FALSE, useFancyQuotes=FALSE)



# ---------- Preparations ----------
# Parallelization
library(parallel)
nSlaves <- parallel::detectCores(all.tests=FALSE, logical=FALSE)
library(doMC)
registerDoMC(nSlaves)

# Load libraries
library(GA)



# ---------- Arguments and user variables ----------
#args <- list()
#args[1] <- "~/Desktop/Projekte/Habilitation/de.NBI/mFam-Classifier"
#args[2] <- "data/2018-02-13_pos_21908_MoNA_Spectra.msp"
#args[3] <- "data/2019-05-23_Scaffolds.tsv"
#args[4] <- "data/ga_output"

# Take in trailing command line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
	print("Error! No or not enough arguments given.")
	print("Usage: $0 working_dir classifier_library annotation_file output_dir")
	quit(save="no", status=1, runLast=FALSE)
}

# Working directory
setwd(as.character(args[1]))

# mFam compound library
mFam_file <- as.character(unlist(args[2]))

# Scaffolds
scaff <- as.character(unlist(args[3]))

# Output dir
out_dir <- as.character(unlist(args[4]))



# ---------- Read MSP ----------
read_msp <- function(msp_file) {
	# Default variables
	minimumIntensityOfMaximalMS2peak <- 0
	minimumProportionOfMS2peaks <- 0.05
	neutralLossesPrecursorToFragments <- TRUE
	neutralLossesFragmentsToFragments <- FALSE
	offset <- 0
	
	# Read file
	msp_lines <- readLines(con=msp_file)
	
	# Numbers
	msp_lines_length <- length(msp_lines)
	numberOfMS2PeaksOriginal <- 0
	numberOfMS2PeaksWithNeutralLosses <- 0
	numberOfTooHeavyFragments <- 0
	numberOfMS2PeaksAboveThreshold <- 0
	numberOfMS2PeaksBelowThreshold <- 0
	numberOfSpectraDiscardedDueToNoPeaks <- 0
	numberOfSpectraDiscardedDueToMaxIntensity <- 0
	numberOfSpectraDiscardedDueToTooHeavy <- 0
	
	# Does file start with empty line?
	endOfRecord <- TRUE
	if (msp_lines_length > 0)
		if (nchar(trimws(msp_lines[[1]])) > 0)
			endOfRecord <- FALSE
	
	# Check for patterns
	isBI <- grepl(pattern="^BEGIN IONS$", x=msp_lines)
	isName <- grepl(pattern="(^Name:)|(^NAME:)", x=msp_lines)
	isNAme <- grepl(pattern="^NAME=", x=msp_lines)
	isNAME <- grepl(pattern="^TITLE=", x=msp_lines)
	isRT <- grepl(pattern="^RETENTIONTIME:", x=msp_lines)
	isRt <- grepl(pattern="^retention time:", x=msp_lines)
	isRts <- grepl(pattern="^RTINSECONDS=", x=msp_lines)
	isMZ <- grepl(pattern="^PRECURSORMZ:", x=msp_lines)
	isMz <- grepl(pattern="^precursor m/z:", x=msp_lines)
	isTotEm <- grepl(pattern="^total exact mass:", x=msp_lines)
	isEMass <- grepl(pattern="(^exact mass:)|(^EXACT_MASS:)", x=msp_lines)
	isPMass <- grepl(pattern="^PEPMASS=", x=msp_lines)
	isE2Mass <- grepl(pattern="^EXACTMASS=", x=msp_lines)
	isMetN <- grepl(pattern="^METABOLITENAME:", x=msp_lines)
	isAddN <- grepl(pattern="(^ADDUCTIONNAME:)|(^Adductionname:)", x=msp_lines)
	isScanN <- grepl(pattern="^SCANNUMBER:", x=msp_lines)
	isMIon <- grepl(pattern="^MODELION:", x=msp_lines)
	isPrety <- grepl(pattern="^precursor type:", x=msp_lines)
	isPretyp <- grepl(pattern="^PRECURSORTYPE:", x=msp_lines)
	isNumP <- grepl(pattern="^Num Peaks:", x=msp_lines)
	isPeak <- grepl(pattern="^\\d+(((\\.)|(,))\\d+)?[ \t]\\d+(((\\.)|(,))\\d+)?$", x=msp_lines)
	isCoCl <- grepl(pattern="^compound class:", x=msp_lines)
	isInty <- grepl(pattern="^instrument type:", x=msp_lines)
	isIntyp <- grepl(pattern="^INSTRUMENTTYPE:", x=msp_lines)
	isIntype <- grepl(pattern="^SOURCE_INSTRUMENT=", x=msp_lines)
	isInt <- grepl(pattern="^INSTRUMENT:", x=msp_lines)
	isInchi <- grepl(pattern="(^InChI:)|(^InChI=)|(^INCHI=)|(^INCHI:)", x=msp_lines)
	isInchiKey <- grepl(pattern="(^InChIKey:)|(^INCHIKEY:)|(^InChIKey=)|(^INCHIKEY=)|(^INCHIAUX=)", x=msp_lines)
	isSmiles <- grepl(pattern="(^SMILES:)|(^SMILES=)", x=msp_lines)
	
	someStrings <- trimws(c(
		substring(text=msp_lines[isMZ], first=nchar("RETENTIONTIME:") + 1), 
		substring(text=msp_lines[isMZ], first=nchar("PRECURSORMZ:") + 1), 
		substring(text=msp_lines[isMz], first=nchar("precursor m/z:") + 1)
	))
	decimalDelimiterIsComma <- ifelse(test=length(someStrings) > 0, yes=all(grepl(x=someStrings, pattern="^(\\d+,\\d+$)|(^\\d+$)")), no=FALSE)
	if (decimalDelimiterIsComma){
		msp_lines_02 <- gsub(pattern=",", replacement=".", x=gsub(pattern="\\.", replacement="", x=msp_lines))
	} else {
		msp_lines_02 <- msp_lines
	}
	
	# Extract
	parsedName <- trimws(substring(text=msp_lines, first=nchar("Name:") + 1))
	parsedNAme <- trimws(substring(text=msp_lines, first=nchar("NAME=") + 1))
	parsedNAME <- trimws(substring(text=msp_lines, first=nchar("TITLE=") + 1))
	parsedRT <- as.numeric(trimws(substring(text=msp_lines_02, first=nchar("RETENTIONTIME:") + 1)))
	parsedRt <- as.numeric(unlist(lapply(X=strsplit(x=  trimws(substring(text=msp_lines_02, first=nchar("retention time:") + 1)), split=" "), FUN=function(x){if (length(x)==0) return(NA) else return(x[[1]])})))
	parsedRts <- as.numeric(trimws(substring(text=msp_lines_02, first=nchar("RTINSECONDS=") + 1)))
	parsedMZ <- as.numeric( trimws(substring(text=msp_lines_02, first=nchar("PRECURSORMZ:") + 1)))
	parsedMz <- as.numeric(trimws(substring(text=msp_lines_02, first=nchar("precursor m/z:") + 1)))
	parsedTotEm <- as.numeric(trimws(substring(text=msp_lines_02, first=nchar("total exact mass:") + 1)))
	parsedEMass <- as.numeric(trimws(substring(text=msp_lines_02, first=nchar("exact mass:") + 1)))
	parsedE2Mass <- as.numeric(trimws(substring(text=msp_lines_02, first=nchar("EXACTMASS=") + 1)))
	parsedPMass <- as.numeric(trimws(substring(text=msp_lines_02, first=nchar("PEPMASS=") + 1)))
	parsedMetN <- trimws(substring(text=msp_lines, first=nchar("METABOLITENAME:") + 1))
	parsedAddN <- trimws(substring(text=msp_lines, first=nchar("Adductionname:") + 1))
	parsedScanN <- as.numeric(trimws(substring(text=msp_lines_02, first=nchar("SCANNUMBER:") + 1)))
	parsedMIon <- as.numeric(trimws(substring(text=msp_lines_02, first=nchar("MODELION:") + 1)))
	parsedPrety <- trimws(substring(text=msp_lines, first=nchar("precursor type:") + 1))
	parsedPretyp <- trimws(substring(text=msp_lines, first=nchar("PRECURSORTYPE:") + 1))
	parsedNumP <- as.numeric(trimws(substring(text=msp_lines_02, first=nchar("Num Peaks:") + 1)))
	
	parsedTokensTmp <- strsplit(x=trimws(msp_lines_02), split="[ \t]")
	parsedms2Peaks_mz <- as.numeric(unlist(lapply(X=parsedTokensTmp, FUN=function(x){if (length(x)<1) return(NA) else return(x[[1]])})))
	parsedms2Peaks_int <- as.numeric(unlist(lapply(X=parsedTokensTmp, FUN=function(x){if (length(x)<2) return(NA) else return(x[[2]])})))
	
	parsedCoCl <- trimws(substring(text=msp_lines, first=nchar("compound class:") + 1))
	parsedInty <- trimws(substring(text=msp_lines, first=nchar("instrument type:") + 1))
	parsedIntype <- trimws(substring(text=msp_lines, first=nchar("SOURCE_INSTRUMENT=") + 1))
	parsedIntyp <- trimws(substring(text=msp_lines, first=nchar("INSTRUMENTTYPE:") + 1))
	parsedInt <- trimws(substring(text=msp_lines, first=nchar("INSTRUMENT:") + 1))
	parsedInchi <- trimws(substring(text=msp_lines, first=nchar("InChI:") + 1))
	parsedInchiKey <- trimws(substring(text=msp_lines, first=nchar("InChIKey:") + 1))
	parsedSmiles <- trimws(substring(text=msp_lines, first=nchar("SMILES:") + 1))
	
	# Entry line intervals
	entryBorders   <- c(which(isName | isNAME | isBI), length(msp_lines)+1)
	entryIntervals <- matrix(data=unlist(lapply(X=seq_len(length(entryBorders) - 1), FUN=function(x){c(entryBorders[[x]], entryBorders[[x+1]] - 1)})), nrow=2)
	
	numberOfSpectraOriginal <- length(entryBorders)
	
	# Parse
	spectraList <- apply(X=entryIntervals, MARGIN=2, FUN=function(x) {
		msp_lines2 <- msp_lines[x[[1]]:x[[2]]]
		msp_lines_022 <- msp_lines_02[x[[1]]:x[[2]]]
		parsedName2 <- parsedName[x[[1]]:x[[2]]]
		parsedNAme2 <- parsedNAme[x[[1]]:x[[2]]]
		parsedNAME2 <- parsedNAME[x[[1]]:x[[2]]]
		parsedRT2 <- parsedRT[x[[1]]:x[[2]]]
		parsedRt2 <- parsedRt[x[[1]]:x[[2]]]
		parsedRts2 <- parsedRts[x[[1]]:x[[2]]]
		parsedMZ2 <- parsedMZ[x[[1]]:x[[2]]]
		parsedMz2 <- parsedMz[x[[1]]:x[[2]]]
		parsedTotEm2 <- parsedTotEm[x[[1]]:x[[2]]]
		parsedEMass2 <- parsedEMass[x[[1]]:x[[2]]]
		parsedE2Mass2 <- parsedE2Mass[x[[1]]:x[[2]]]
		parsedPMass2 <- parsedPMass[x[[1]]:x[[2]]]
		parsedMetN2 <- parsedMetN[x[[1]]:x[[2]]]
		parsedAddN2 <- parsedAddN[x[[1]]:x[[2]]]
		parsedScanN2 <- parsedScanN[x[[1]]:x[[2]]]
		parsedMIon2 <- parsedMIon[x[[1]]:x[[2]]]
		parsedPrety2 <- parsedPrety[x[[1]]:x[[2]]]
		parsedPretyp2 <- parsedPretyp[x[[1]]:x[[2]]]
		parsedNumP2 <- parsedNumP[x[[1]]:x[[2]]]
		parsedms2Peaks_mz2 <- parsedms2Peaks_mz[x[[1]]:x[[2]]]
		parsedms2Peaks_int2 <- parsedms2Peaks_int[x[[1]]:x[[2]]]
		parsedCoCl2 <- parsedCoCl[x[[1]]:x[[2]]]
		parsedInty2 <- parsedInty[x[[1]]:x[[2]]]
		parsedIntype2 <- parsedIntype[x[[1]]:x[[2]]]
		parsedIntyp2 <- parsedIntyp[x[[1]]:x[[2]]]
		parsedInt2 <- parsedInt[x[[1]]:x[[2]]]
		parsedInchi2 <- parsedInchi[x[[1]]:x[[2]]]
		parsedInchiKey2 <- parsedInchiKey[x[[1]]:x[[2]]]
		parsedSmiles2 <- parsedSmiles[x[[1]]:x[[2]]]
		
		isName2 <- isName[x[[1]]:x[[2]]]
		isNAme2 <- isName[x[[1]]:x[[2]]]
		isNAME2 <- isNAME[x[[1]]:x[[2]]]
		isRT2 <- isRT[x[[1]]:x[[2]]]
		isRt2 <- isRt[x[[1]]:x[[2]]]
		isRts2 <- isRts[x[[1]]:x[[2]]]
		isMZ2 <- isMZ[x[[1]]:x[[2]]]
		isMz2 <- isMz[x[[1]]:x[[2]]]
		isTotEm2 <- isTotEm[x[[1]]:x[[2]]]
		isEMass2 <- isEMass[x[[1]]:x[[2]]]
		isE2Mass2 <- isE2Mass[x[[1]]:x[[2]]]
		isPMass2 <- isPMass[x[[1]]:x[[2]]]
		isMetN2 <- isMetN[x[[1]]:x[[2]]]
		isAddN2 <- isAddN[x[[1]]:x[[2]]]
		isScanN2 <- isScanN[x[[1]]:x[[2]]]
		isMIon2 <- isMIon[x[[1]]:x[[2]]]
		isPrety2 <- isPrety[x[[1]]:x[[2]]]
		isPretyp2 <- isPretyp[x[[1]]:x[[2]]]
		isNumP2 <- isNumP[x[[1]]:x[[2]]]
		isPeak2 <- isPeak[x[[1]]:x[[2]]]
		isCoCl2 <- isCoCl[x[[1]]:x[[2]]]
		isInty2 <- isInty[x[[1]]:x[[2]]]
		isIntype2 <- isIntype[x[[1]]:x[[2]]]
		isIntyp2 <- isIntyp[x[[1]]:x[[2]]]
		isInt2 <- isInt[x[[1]]:x[[2]]]
		isInchi2 <- isInchi[x[[1]]:x[[2]]]
		isInchiKey2 <- isInchiKey[x[[1]]:x[[2]]]
		isSmiles2 <- isSmiles[x[[1]]:x[[2]]]
		
		name <- NULL
		ms1Int <- NA
		rt <- NULL
		mz <- NULL
		metName <- "Unknown"
		adduct <- "Unknown"
		scanNumber <- NA
		quantMass <- NA
		peakNumber <- NA
		ms2Peaks_mz  <- vector(mode="numeric")
		ms2Peaks_int <- vector(mode="numeric")
		compoundClass <- "Unknown"
		instrumentType <- "Unknown"
		inchi <- ""
		inchiKey <- ""
		smiles <- ""
		
		if (any(isName2)) name <- parsedName2[[which(isName2)[[1]]]]
		if (any(isNAme2)) name <- parsedNAme2[[which(isNAme2)[[1]]]]
		if (any(isNAME2)) name <- parsedNAME2[[which(isNAME2)[[1]]]]
		
		if (any(isRT2)) rt <- parsedRT2[[which(isRT2)[[1]]]]
		if (any(isRt2) & any(is.null(rt), is.na(rt))) rt <- parsedRt2[[which(isRt2)[[1]]]]
		if (any(isRts2) & any(is.null(rt), is.na(rt))) rt <- parsedRts2[[which(isRts2)[[1]]]]
		
		if (any(isMZ2)) mz <- parsedMZ2[[which(isMZ2)[[1]]]]
		if (any(isMz2) & any(is.null(mz), is.na(mz), ifelse(test=is.null(mz), yes=FALSE, no=mz%%1==0))) mz <- parsedMz2[[which(isMz2)[[1]]]]
		
		if (any(isTotEm2) & any(is.null(mz), is.na(mz), ifelse(test=is.null(mz), yes=FALSE, no=mz%%1==0)) & any(isPrety2)){
			mzTmp <- parsedTotEm2[[which(isTotEm2)[[1]]]]
			pit <- parsedPrety2[[which(isPrety2)[[1]]]]
			switch(pit,
				   "[M-H]-"={ mz <- mzTmp - 1.008 },
				   "[M-H]" ={ mz <- mzTmp - 1.008 },
				   "[M+H]+"={ mz <- mzTmp + 1.008 },
				   "[M+H]" ={ mz <- mzTmp + 1.008 },
				   "[M+Na]+"= { mz <- mzTmp + 22.9898 },
				   "[M+Na]"={ mz <- mzTmp + 22.9898 }
			)
		}
		if (any(isEMass2) & any(is.null(mz), is.na(mz), ifelse(test=is.null(mz), yes=FALSE, no=mz%%1==0))){
			mz <- parsedEMass2[[which(isEMass2)[[1]]]]
		}
		if (any(isPMass2) & any(is.null(mz), is.na(mz), ifelse(test=is.null(mz), yes=FALSE, no=mz%%1==0))){
			if (!is.na(parsedPMass2[[which(isPMass2)[[1]]]])){
				mz <- parsedPMass2[[which(isPMass2)[[1]]]]
			} else {
				tmp <- as.numeric(strsplit(x=trimws(substring(text=msp_lines_022[[which(isPMass2)[[1]]]], first=nchar("PEPMASS=") + 1)), split="[\t ]")[[1]])
				mz <- tmp[[1]]
				if (length(tmp) > 1)
					ms1Int <- tmp[[2]]
			}
		}
		if (any(isE2Mass2) & any(is.null(mz), is.na(mz), ifelse(test=is.null(mz), yes=FALSE, no=mz%%1==0))){
			mz <- parsedE2Mass2[[which(isE2Mass2)[[1]]]]
		}
		
		if (any(isMetN2)) metName <- parsedMetN2[[which(isMetN2)[[1]]]]
		
		if (any(isAddN2) & adduct == "Unknown") adduct <- parsedAddN2[[which(isAddN2)[[1]]]]
		if (any(isPrety2) & adduct == "Unknown") adduct <- parsedPrety2[[which(isPrety2)[[1]]]]
		if (any(isPretyp2) & adduct == "Unknown") adduct <- parsedPretyp2[[which(isPretyp2)[[1]]]]
		
		if (any(isScanN2)) scanNumber <- parsedScanN2[[which(isScanN2)[[1]]]]
		if (any(isMIon2)) quantMass <- parsedMIon2[[which(isMIon2)[[1]]]]
		
		if (any(isNumP2)) peakNumber <- parsedNumP2[[which(isNumP2)[[1]]]]
		
		if (any(isPeak2)) {
			ms2Peaks_mz  <- parsedms2Peaks_mz2 [which(isPeak2)]
			ms2Peaks_int <- parsedms2Peaks_int2[which(isPeak2)]
		}
		
		if (any(isCoCl2)) compoundClass <- parsedCoCl2[[which(isCoCl2)[[1]]]]
		
		if (any(isInty2)) instrumentType <- parsedInty2[[which(isInty2)[[1]]]]
		if (any(isIntype2)) instrumentType <- parsedIntype2[[which(isIntype2)[[1]]]]
		if (any(isIntyp2)) instrumentType <- parsedIntyp2[[which(isIntyp2)[[1]]]]
		if (all(any(isInt2), instrumentType %in% c("Unknown", "NA"))) instrumentType <- parsedInt2[[which(isInt2)[[1]]]]
		
		if (any(isInchi2)) inchi <- parsedInchi2[[which(isInchi2)[[1]]]]
		if (any(isInchiKey2)) inchiKey <- parsedInchiKey2[[which(isInchiKey2)[[1]]]]
		if (any(isSmiles2)) smiles <- parsedSmiles2[[which(isSmiles2)[[1]]]]
		
		if (is.null(rt)) rt <- 0
		
		if (is.null(mz)) mz <- max(ms2Peaks_mz)
		ms2Peaks_mz_original  <- ms2Peaks_mz
		ms2Peaks_int_original <- ms2Peaks_int
		
		numberOfMS2PeaksOriginal <- numberOfMS2PeaksOriginal + length(ms2Peaks_mz)
		if (length(ms2Peaks_mz) == 0) numberOfSpectraDiscardedDueToNoPeaks <- numberOfSpectraDiscardedDueToNoPeaks + 1
		
		# Remove fragments with mass greater than precursor
		numberOfTooHeavyFragmentsHere <- 0
		if (all(!is.null(mz), !is.na(mz))) {
			tooHeavy <- ms2Peaks_mz > mz
			ms2Peaks_mz <- ms2Peaks_mz[!tooHeavy]
			ms2Peaks_int <- ms2Peaks_int[!tooHeavy]
			numberOfTooHeavyFragmentsHere <- sum(tooHeavy)
			
			if (length(ms2Peaks_mz) == 0 & numberOfTooHeavyFragmentsHere > 0) numberOfSpectraDiscardedDueToTooHeavy <- numberOfSpectraDiscardedDueToTooHeavy + 1
		}
		numberOfTooHeavyFragments <- numberOfTooHeavyFragments + numberOfTooHeavyFragmentsHere
		
		# Filter for MS2 peak intensity
		peakNumber <- length(ms2Peaks_mz)
		if (peakNumber > 0) {
			# Filter for ms2 peak intensity relative to maximum peak intensity in spectrum
			maximumIntensity <- max(ms2Peaks_int)
			
			# Spectrum
			if (maximumIntensity >= minimumIntensityOfMaximalMS2peak) {
				intensityThreshold <- maximumIntensity * minimumProportionOfMS2peaks
				fragmentsAboveThreshold <- ms2Peaks_int >= intensityThreshold
				
				ms2Peaks_mz  <- ms2Peaks_mz[fragmentsAboveThreshold]
				ms2Peaks_int <- ms2Peaks_int[fragmentsAboveThreshold]
				numberOfMS2PeaksAboveThreshold <- numberOfMS2PeaksAboveThreshold + sum(fragmentsAboveThreshold)
				numberOfMS2PeaksBelowThreshold <- numberOfMS2PeaksBelowThreshold + sum(!fragmentsAboveThreshold)
				# No spectrum
			} else {
				numberOfMS2PeaksBelowThreshold <- numberOfMS2PeaksBelowThreshold + length(ms2Peaks_mz)
				numberOfSpectraDiscardedDueToMaxIntensity <- numberOfSpectraDiscardedDueToMaxIntensity + 1
				ms2Peaks_mz <- vector(mode="numeric")
				ms2Peaks_int <- vector(mode="numeric")
			}
		}
		
		# Normalize MS2 peaks to maximum=1
		peakNumber <- length(ms2Peaks_mz)
		if (peakNumber > 0) {
			max <- max(ms2Peaks_int)
			ms2Peaks_int <- ms2Peaks_int / max
		}
		
		# Add neutral losses
		if (peakNumber > 0) {
			# Neutral losses with regard to precursor
			if (all(!is.null(mz), !is.na(mz), neutralLossesPrecursorToFragments)) {
				ms2PeaksNLPF_mz <- ms2Peaks_mz - as.numeric(mz)
				ms2PeaksNLPF_int <- ms2Peaks_int
			} else {
				ms2PeaksNLPF_mz <- vector(mode="numeric")
				ms2PeaksNLPF_int <- vector(mode="numeric")
			}
			
			# Neutral losses with regard to fragments
			if (neutralLossesFragmentsToFragments) {
				m_mz  <- outer(X=ms2Peaks_mz, Y=ms2Peaks_mz, FUN=function(x,y) {x-y})
				m_int <- outer(X=ms2Peaks_int, Y=ms2Peaks_int, FUN=function(x,y) {(x+y) / 2})
				upper <- upper.tri(x=m_mz)
				ms2PeaksNLFF_mz <- m_mz[upper]
				ms2PeaksNLFF_int <- m_int[upper]
			} else {
				ms2PeaksNLFF_mz <- vector(mode="numeric")
				ms2PeaksNLFF_int <- vector(mode="numeric")
			}
			
			ms2Peaks_mz <- c(ms2Peaks_mz, ms2PeaksNLPF_mz, ms2PeaksNLFF_mz)
			ms2Peaks_int <- c(ms2Peaks_int, ms2PeaksNLPF_int, ms2PeaksNLFF_int)
		}
		
		# Precursor mz
		if (all(!is.null(mz), !is.na(mz))) {
			mz <- round(as.numeric(mz), digits=4)
		} else {
			if (!is.na(quantMass)){
				mz <- quantMass
			} else {
				if (!is.na(scanNumber)){
					mz <- scanNumber
				} else {
					mz <- max(ms2Peaks_mz)
				}
			}
		}
		
		# Spectrum string
		spectrumString <- paste(ms2Peaks_mz_original, ms2Peaks_int_original, sep=" ", collapse=";")
		
		# Spectrum set
		spectrumItem <- list(
			name=name,
			ms1Int=ms1Int,
			rt=rt,
			mz=mz,
			metName=metName,
			adduct=adduct,
			quantMass=quantMass,
			compoundClass=compoundClass,
			instrumentType=instrumentType,
			inchi=inchi,
			inchiKey=inchiKey,
			smiles=smiles,
			peakNumber=length(ms2Peaks_mz),
			ms2Peaks_mz =ms2Peaks_mz,
			ms2Peaks_int=ms2Peaks_int,
			spectrumString=spectrumString,
			entryInterval=x + offset
		)
		if (spectrumItem$peakNumber > 0) {
			numberOfMS2PeaksWithNeutralLosses <- numberOfMS2PeaksWithNeutralLosses + spectrumItem$peakNumber
			return(spectrumItem)
		} else {
			return(NULL)
		}
	})
	
	# Remove NULL entries
	spectraList[unlist(lapply(X=spectraList, FUN=is.null))] <- NULL
	
	numberOfSpectra <- length(spectraList)
	
	# Postprocess
	precursorMz <- unlist(lapply(X=spectraList, FUN=function(x) {
		if (is.null(x)) return(NULL)
		else return(x$mz)
	}))
	precursorRt <- unlist(lapply(X=spectraList, FUN=function(x) {
		if (is.null(x)) return(NULL)
		else return(x$rt)
	}))
	
	returnObj <- list()
	returnObj$fileSpectra <- NA
	returnObj$spectraList <- spectraList
	returnObj$numberOfSpectra <- numberOfSpectra
	returnObj$numberOfSpectraOriginal <- numberOfSpectraOriginal
	returnObj$numberOfMS2PeaksOriginal <- numberOfMS2PeaksOriginal
	returnObj$numberOfMS2PeaksWithNeutralLosses <- numberOfMS2PeaksWithNeutralLosses
	returnObj$numberOfMS2PeaksAboveThreshold <- numberOfMS2PeaksAboveThreshold
	returnObj$numberOfMS2PeaksBelowThreshold <- numberOfMS2PeaksBelowThreshold
	returnObj$numberOfTooHeavyFragments <- numberOfTooHeavyFragments
	returnObj$numberOfSpectraDiscardedDueToNoPeaks <- numberOfSpectraDiscardedDueToNoPeaks
	returnObj$numberOfSpectraDiscardedDueToMaxIntensity <- numberOfSpectraDiscardedDueToMaxIntensity
	returnObj$numberOfSpectraDiscardedDueToTooHeavy <- numberOfSpectraDiscardedDueToTooHeavy
	returnObj$precursorMz <- precursorMz
	returnObj$precursorRt <- precursorRt
	
	return(returnObj)
}



# ---------- Write MSP ----------
write_msp <- function(msp_filename, msp_lib, indices=NULL) {
	# Length / number of spectra
	msp_lib_len <- length(msp_lib[[2]])
	
	# Write all spectra
	if (is.null(indices)) {
		indices <- rep(1, msp_lib_len)
	}
	
	# Loop through spectra
	msp_file <- NULL
	for (msp_entry in 1:msp_lib_len) {
		# Do not write spectra with index of 0
		if (indices[msp_entry] == 1) {
			msp_file[length(msp_file)+1] <- paste0("NAME: ", msp_lib[[2]][[msp_entry]]$name)
			msp_file[length(msp_file)+1] <- paste0("RETENTIONTIME: ", msp_lib[[2]][[msp_entry]]$rt)
			msp_file[length(msp_file)+1] <- paste0("PRECURSORMZ: ", msp_lib[[2]][[msp_entry]]$mz)
			msp_file[length(msp_file)+1] <- paste0("PRECURSORTYPE: ", msp_lib[[2]][[msp_entry]]$adduct)
			msp_file[length(msp_file)+1] <- paste0("IONMODE: ", msp_lib[[2]][[msp_entry]]$adduct)
			msp_file[length(msp_file)+1] <- paste0("INCHIKEY: ", msp_lib[[2]][[msp_entry]]$inchikey)
			msp_file[length(msp_file)+1] <- paste0("INCHI: ", msp_lib[[2]][[msp_entry]]$inchi)
			msp_file[length(msp_file)+1] <- paste0("SMILES: ", msp_lib[[2]][[msp_entry]]$smiles)
			msp_file[length(msp_file)+1] <- paste0("INTENSITY: ", msp_lib[[2]][[msp_entry]]$ms1Int)
			msp_file[length(msp_file)+1] <- paste0("INSTRUMENTTYPE: ", msp_lib[[2]][[msp_entry]]$instrumentType)
			msp_file[length(msp_file)+1] <- paste0("Num Peaks: ", msp_lib[[2]][[msp_entry]]$peakNumber)
			for (peak in unlist(strsplit(msp_lib[[2]][[msp_entry]]$spectrumString, ";"))) {
				msp_file[length(msp_file)+1] <- gsub(x=peak, pattern=" ", replacement="\t")
			}
			msp_file[length(msp_file)+1] <- ""	
		}
	}
	
	# Save file
	writeLines(text=msp_file, con=msp_filename)
}



# ---------- Evaluaton function of Genetic Algorithm fitness ----------
mFam_ga_fitness <- function(x) {
	# Options
	exec_docker <- FALSE
	
	# Define fitness score
	fitness_score <- 0
	
	# Create output dir
	ga_run_id <- ga_run_id + 1
	#ga_out_dir <- paste0(out_dir,"/ga_",sprintf("%08s",as.character(ga_run_id)))
	ga_out_dir <- paste0(out_dir,"/ga_", formatC(as.numeric(ga_run_id), width=8, format="d", flag="0"))
	dir.create(path=ga_out_dir, recursive=TRUE, mode="0755")
	
	print(ga_out_dir)
	x <- round(runif(mFam_num_spectra,0,1), 0)
	print(x)
	
	# Apply binary vector to mFam library
	write_msp(msp_filename=paste0(ga_out_dir,"/","lib.msp"), msp_lib=mFam_lib, indices=x)
	gen_mFam_lib <- read_msp(paste0(ga_out_dir,"/","lib.msp"))
	
	# Execute mFam classifier
	if (exec_docker == FALSE) {
		gen_mFam_exec <- system2(command="./galaxy/mFam_train_classifier.r",
								 args=c(getwd(),
								 	   paste0(ga_out_dir,"/lib.msp"),
								 	   paste0(getwd(),"/data/2019-05-23_Scaffolds.tsv"),
								 	   paste0(ga_out_dir)),
								 wait=TRUE, timeout=0)
	} else {
		gen_mFam_exec <- system2(command="docker",
								 args=c("run",
								 	   "-ti",
								 	   "-v",
								 	   paste0(getwd(),":",getwd()),
								 	   "ipb-halle/mfam-classifier",
								 	   "/usr/local/bin/mFam_train_classifier.sh",
								 	   paste0(ga_out_dir,"/lib.msp"),
								 	   paste0(getwd(),"/data/2019-05-23_Scaffolds.tsv"),
								 	   paste0(ga_out_dir,"/out_obj.rdata"),
								 	   paste0(ga_out_dir,"/out_obj.tsv"),
								 	   paste0(ga_out_dir,"/out_obj.txt")),
								 wait=TRUE, timeout=0)
	}
	if ((gen_mFam_exec != 0) | (gen_mFam_exec == 127)) {
		print(paste0("Warning. GA run #", ga_out_dir, " exited with errors."))
		fitness_score <- 0
		return(fitness_score)
	}
	
	# Evaluate mFam classifier results
	gen_mFam_results_file <- list.files(ga_out_dir, pattern="*Results.tsv", recursive=FALSE, full.names=TRUE)
	if (length(gen_mFam_results_file) == 0) {
		fitness_score <- 0
		return(fitness_score)
	} else {
		gen_mFam_results <- read.table(file=gen_mFam_results_file, sep='\t', header=TRUE, stringsAsFactors=TRUE, fill=TRUE)
		score_auc_pr <- sum(gen_mFam_results$AUC.PR)
		score_num_classes <- nrow(gen_mFam_results)
		score_tpr <- sum(gen_mFam_results$TPR.for.FPR...5.)
	}
	
	# Calculate fitness score
	fitness_score <- score_auc_pr / score_num_classes
	
	print(paste0("GA #", ga_out_dir, ", Fitness Score: ", fitness_score))
	
	return(fitness_score)
}



# ---------- Read mFam library ----------
mFam_lib <- read_msp(mFam_file)
mFam_num_spectra <- mFam_lib$numberOfSpectra
ga_run_id <- 0

# ---------- Perform Genetic Algorithm here ----------
model_ga <- ga(type="binary",          # Optimization data type
			   fitness=mFam_ga_fitness,# Fitness function
			   elitism=3,              # Number of best individuals (compounds) to pass to next iteration
			   pmutation=1/100,        # Mutation rate probability
			   popSize=10,             # Number of individuals / solutions
			   nBits=mFam_num_spectra, # Total number of variables in compound matrix
			   run=10,                 # Max iterations without improvement (stopping criteria)
			   maxiter=10,             # Max iterations
			   monitor=gaMonitor,      # Monitoring of intermediate results: plot | gaMonitor | FALSE
			   keepBest=TRUE,          # Keep the best solution at the end
			   parallel=TRUE           # Allow parallel procesing
)

# Print summary on Genetic Algorithm
plot(model_ga, main="Genetic Algorithm Performance", cex.points=0.9, col=c("dodgerblue4", "dodgerblue3",  adjustcolor("dodgerblue2", alpha.f=0.1)), pch=c(19,16), lty=c(1,2), legend=TRUE, grid=graphics:::grid)
summary(model_ga)
best_ga_solution <- as.numeric(model_ga@solution[1,])
#knapsack[best_ga_solution == 1, ]
#cat(best_ga_solution %*% knapsack$survivalpoints)
#cat(best_ga_solution %*% knapsack$weight)


