###################################################
##  Main Control Script to process model output  ##
###################################################

###########################################################################################
## RUN (do not change anything below this line)

library(rwrfhydro)
library(data.table)

load(maskFile)
source("util_FUNC.R")


# Model Reads 
if (readMod | readForc) {
        source("read_MODELOUT.R")
}

# Obs
if (calcStats | createPlots) {
	if (!is.null(AMFfile) & amfProc & exists("ptgeo.amf")) {
		if (file.exists(AMFfile)) {
			load(AMFfile)
			obsAmfData <- subset(obsAmfData, obsAmfData$site_id %in% ptgeo.amf$id)
		} else {
			stop(paste("Ameriflux obs file specified but does not exist:", AMFfile))
		}
	}
	if (!is.null(SNOfile) & (snoProc | swePlot) & exists("ptgeo.sno")) {
        	if (file.exists(SNOfile)) {
                	load(SNOfile)
			obsSnoData <- subset(obsSnoData, obsSnoData$site_id %in% ptgeo.sno$id)
		} else {
                	stop(paste("SNOTEL obs file specified but does not exist:", SNOfile))
       		}
	}
	if (!is.null(METfile) & metProc & exists("ptgeo.met")) {
        	if (file.exists(METfile)) {
                	load(METfile)
			obsMetData <- subset(obsMetData, obsMetData$site_id %in% ptgeo.met$id)
		} else {
                	stop(paste("MET obs file specified but does not exist:", METfile))
        	}
	}
	if (!is.null(STRfile) & ( (calcStats & strProc) | (createPlots & (accflowPlot | hydroPlot | flowswePlot | flowlsmPlot)) ) ) {
		obsStrData_FINAL <- data.frame()
		obsStrMeta_FINAL <- data.frame()
		# Gridded routing case w/ subset
		if ( !reachRting & exists("stid2gageList") ) {
			if (is.list(stid2gageList)) {
				gageList <- data.frame(st_id=names(stid2gageList), site_no=unlist(stid2gageList), stringsAsFactors=FALSE)
			} else {
				gageList <- stid2gageList
			}
		# Reach routing case w/ subset
		} else if (reachRting) {
			if ( calcStats & strProc & exists("statsLink2gage") & !is.null(statsLink2gage) ) {
				gageList <- statsLink2gage[,c("link","site_no")]
			} else if ( createPlots & (accflowPlot | hydroPlot | flowswePlot | flowlsmPlot) & exists("plotLink2gage") & !is.null(plotLink2gage) ) {
				gageList <- plotLink2gage[,c("link","site_no")]
			} else {
				gageList <- subset(rtLinks[,c("link","site_no")], !(rtLinks$site_no == ''))
			}
		# No subset
		} else {
			gageList <- NULL
		} 
		for (i in STRfile) {
			if (file.exists(i)) {
				load(i)
				if (exists("obsStrData.map")) obsStrData <- remapData(obsStrData, obsStrData.map)
				if (exists("obsStrMeta.map")) obsStrMeta <- remapData(obsStrMeta, obsStrMeta.map)
				if (is.data.table(obsStrData)) {
					message("Obs are in data.table format")
					obsStrData_FINAL <- data.table()
                                        if ( !is.null(gageList) ) {
                                                obsStrData_TMP <- obsStrData[site_no %in% unique(gageList$site_no),]
                                                obsStrMeta_TMP <- subset(obsStrMeta, obsStrMeta$site_no %in% unique(gageList$site_no))
                                        } else {
                                                obsStrData_TMP <- obsStrData
                                                obsStrMeta_TMP <- obsStrMeta
                                        }
                                        obsStrData_FINAL <- rbindlist(list(obsStrData_FINAL, obsStrData_TMP), fill=TRUE)
                                        obsStrMeta_FINAL <- plyr::rbind.fill(obsStrMeta_FINAL, obsStrMeta_TMP)
				} else {
					if ( !is.null(gageList) ) { 
						obsStrData_TMP <- subset(obsStrData, obsStrData$site_no %in% unique(gageList$site_no))
						obsStrMeta_TMP <- subset(obsStrMeta, obsStrMeta$site_no %in% unique(gageList$site_no))
					} else {
						obsStrData_TMP <- obsStrData
						obsStrMeta_TMP <- obsStrMeta
					}
					obsStrData_FINAL <- plyr::rbind.fill(obsStrData_FINAL, obsStrData_TMP)
                        		obsStrMeta_FINAL <- plyr::rbind.fill(obsStrMeta_FINAL, obsStrMeta_TMP)
				}
			} else {
				stop(paste("Streamflow obs file specified but does not exist:", STRfile))
			}
		}
		obsStrData <- obsStrData_FINAL
		obsStrMeta <- obsStrMeta_FINAL
		if ( reachRting & !is.null(gageList) ) {
			obsStrData <- plyr::join(obsStrData, gageList, by="site_no")
		}
		rm(obsStrData_FINAL, obsStrMeta_FINAL, obsStrData_TMP, obsStrMeta_TMP)
	}
}

# Stats Calculations
if (calcStats & (strProc | snoProc | amfProc | metProc)) {
	message("Calculating stats")
	if (is.null(modReadFileOut)) {
		if (file.exists(modReadFileIn)) {
			load(modReadFileIn)
		}
	} else {
		if (is.null(modReadFileIn)) {
			if (file.exists(modReadFileOut)) {
				load(modReadFileOut)
			}
		} else {
			if (file.exists(modReadFileIn)) {
				load(modReadFileIn)
			}
		}
	}
	if (metProc) {
        	if (is.null(forcReadFileOut)) {
                	if (file.exists(forcReadFileIn)) {
                        	load(forcReadFileIn)
                	}
        	} else {
                	if (is.null(forcReadFileIn)) {
                        	if (file.exists(forcReadFileOut)) {
                                	load(forcReadFileOut)
                        	}
                	} else {
                        	if (file.exists(forcReadFileIn)) {
                                	load(forcReadFileIn)
                        	}
                	}
        	}
	}
	source("calc_PERFSTATS.R")
}

# Plots
if (createPlots) {
	if (accflowPlot | hydroPlot | accprecipPlot | 
			flowswePlot | flowlsmPlot | swePlot | 
			strBiasMap | strCorrMap | 
			snosweErrMap | snoprecipErrMap |
			amfetErrMap | amfetCorrMap) {
        	message("Generating plots")
		if (strBiasMap | strCorrMap | 
			snosweErrMap | snoprecipErrMap |
			amfetErrMap | amfetCorrMap) {
			load(statsFileOut)
		}
        	if (is.null(modReadFileOut)) {
                	if (file.exists(modReadFileIn)) {
                        	load(modReadFileIn)
                	}
        	} else {
                	if (is.null(modReadFileIn)) {
                        	if (file.exists(modReadFileOut)) {
                                	load(modReadFileOut)
                        	}
                	} else {
                        	if (file.exists(modReadFileIn)) {
                                	load(modReadFileIn)
                        	}
                	}
        	}
	}
	if (metPlot) {
                if (is.null(forcReadFileOut)) {
                        if (file.exists(forcReadFileIn)) {
                                load(forcReadFileIn)
                        }
                } else {
                        if (is.null(forcReadFileIn)) {
                                if (file.exists(forcReadFileOut)) {
                                        load(forcReadFileOut)
                                }
                        } else {
                                if (file.exists(forcReadFileIn)) {
                                        load(forcReadFileIn)
                                }
                        }
                }
        }
        source("calc_PLOTS.R")
}

# EXIT
#quit("no")

