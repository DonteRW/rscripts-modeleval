###################################################
##                GENERATE PLOTS                 ##
###################################################


################## General Setup ##################

source("PlotBasin.R")
source("PlotSnotel.R")
source("PlotMaps.R")
library(ggplot2)

# Multi-core
parallelFlag <- FALSE
if (ncores>1) {
        library(doParallel)
        cl <- makeForkCluster(ncores)
        registerDoParallel(cl)
        parallelFlag <- TRUE
}

## -----------------------------------------------------------------------
# Initial configuration setup

# Setup palettes

#lineColors <- c(scales::alpha("dodgerblue", 0.8), scales::alpha("darkorange1", 0.7), "olivedrab", "chocolate", "darkmagenta")
lineColors <- c(scales::alpha("darkorange1", 0.7), scales::alpha("dodgerblue", 0.8), scales::alpha("olivedrab", 0.8), 
       scales::alpha("darkmagenta", 0.8), scales::alpha("chocolate", 0.8))
lineTyp <- 1
#lineWd <- c(1,3,1)
lineWd <- c(3,3,1,1,1)

# Sequential palettes
seqColPurp5 <- c('#edf8fb', '#b3cde3', '#8c96c6', '#8856a7', '#810f7c')
seqColGrn5 <- c("#f7f7f7", "#ffffcc", "#c2e699", "#78c679", "#238443")

# Divergent palettes
divColBluWhtRed6 <- c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020", "#800000")

divColBluYelRed6 <- c('#2c7bb6', '#abd9e9', '#ffffbf', '#fdae61', '#d7191c', '#800000')
divColBluYelRed7 <- c('#26466D', '#2c7bb6', '#abd9e9', '#ffffbf', '#fdae61', '#d7191c', '#800000')

divColRedYelBlu6 <- c('#d7191c', '#fdae61', '#ffffbf', '#abd9e9', '#2c7bb6', '#162252')
divColRedYelBlu7 <- c('#800000', '#d7191c', '#fdae61', '#ffffbf', '#abd9e9', '#2c7bb6', '#162252')

divColBrnWhtGrn5 <- c('#a6611a', '#dfc27d', '#f5f5f5', '#80cdc1', '#018571')
divColBrnWhtGrn6 <- c('#5E2605', '#a6611a', '#dfc27d', '#f5f5f5', '#80cdc1', '#018571')

# Setup plot directories

dir.create(writePlotDir, showWarnings = FALSE)
dir.create(paste0(writePlotDir,"/HYDRO"), showWarnings = FALSE)
dir.create(paste0(writePlotDir,"/CLIM"), showWarnings = FALSE)
dir.create(paste0(writePlotDir,"/SNOW"), showWarnings = FALSE)
dir.create(paste0(writePlotDir,"/MAPS"), showWarnings = FALSE)
dir.create(paste0(writePlotDir,"/BOXPLOTS"), showWarnings = FALSE)

if (is.null(colParseHtml) | !exists("plotLink2gage")) {
	parseList <- "global"
} else {
	parseList <- unique(plotLink2gage[,colParseHtml])
}

print(parseList)
foreach(parseGrp=parseList, .packages = c("data.table")) %dopar% {
#for (parseGrp in parseList) {
	if (parseGrp == "global") { 
		if (!is.null(plotLink2gage)) {
        		gageList <- plotLink2gage
		} else {
        		if (!exists("gageList") & exists("rtLinks")) {
                		gageList<-subset(rtLinks[c("link","site_no")], rtLinks$gages!="")
                		#names(link2gage)[names(link2gage)=="gages"]<-"site_no"
        		}
		}
	} else {
		gageList <- subset(plotLink2gage, plotLink2gage[,colParseHtml] == parseGrp)
	}
print("In dopar")
# Initialize HTML

if (writeHtml) {
	library(knitr)
	library(pander)
	library(xtable)
	if (accflowPlot | hydroPlot | flowlsmPlot) {
		writeLines('```{r set-options, echo=FALSE, cache=FALSE}\noptions(width=1600)\nopts_chunk$set(comment = "", warning = FALSE, message = FALSE, echo = TRUE, tidy = FALSE, size="small")\n```', 
			con=paste0(writePlotDir,"/HYDRO/plots_hydro_", parseGrp, ".Rmd"))
		cat('# MODEL OUTPUT: HYDROLOGY\n', file=paste0(writePlotDir,"/HYDRO/plots_hydro_", parseGrp, ".Rmd"), append=TRUE)
		cat(paste0('##', parseGrp, '\n'), file=paste0(writePlotDir,"/HYDRO/plots_hydro_", parseGrp, ".Rmd"), append=TRUE)
	}
	if (accprecipPlot) {
                writeLines('```{r set-options, echo=FALSE, cache=FALSE}\noptions(width=1600)\nopts_chunk$set(comment = "", warning = FALSE, message = FALSE, echo = TRUE, tidy = FALSE, size="small")\n```', 
			con=paste0(writePlotDir,"/CLIM/plots_climate_", parseGrp, ".Rmd"))
                cat('# MODEL OUTPUT: CLIMATE\n', file=paste0(writePlotDir,"/CLIM/plots_climate_", parseGrp, ".Rmd"), append=TRUE)
		cat(paste0('##', parseGrp, '\n'), file=paste0(writePlotDir,"/CLIM/plots_climate_", parseGrp, ".Rmd"), append=TRUE)
	}
	if (flowswePlot | swePlot) {
                writeLines('```{r set-options, echo=FALSE, cache=FALSE}\noptions(width=1600)\nopts_chunk$set(comment = "", warning = FALSE, message = FALSE, echo = TRUE, tidy = FALSE, size="small")\n```', 
			con=paste0(writePlotDir,"/SNOW/plots_snow_", parseGrp, ".Rmd"))
                cat('# MODEL OUTPUT: SNOW\n', file=paste0(writePlotDir,"/SNOW/plots_snow_", parseGrp, ".Rmd"), append=TRUE)
		cat(paste0('##', parseGrp, '\n'), file=paste0(writePlotDir,"/SNOW/plots_snow_", parseGrp, ".Rmd"), append=TRUE)
	}
        if (strBiasMap | strCorrMap | snosweErrMap | snoprecipErrMap | amfetErrMap | amfetCorrMap) {
                writeLines('```{r set-options, echo=FALSE, cache=FALSE}\noptions(width=1600)\nopts_chunk$set(comment = "", warning = FALSE, message = FALSE, echo = TRUE, tidy = FALSE, size="small")\n```', 
			con=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"))
                cat('# MODEL OUTPUT: MAPS\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
		cat(paste0('##', parseGrp, '\n'), file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
        }
	if (strBiasEcoBoxplot | strCorrEcoBoxplot) {
                writeLines('```{r set-options, echo=FALSE, cache=FALSE}\noptions(width=1600)\nopts_chunk$set(comment = "", warning = FALSE, message = FALSE, echo = TRUE, tidy = FALSE, size="small")\n```', 
                        con=paste0(writePlotDir,"/BOXPLOTS/plots_boxplots_", parseGrp, ".Rmd"))
                cat('# MODEL OUTPUT: BOXPLOTS\n', file=paste0(writePlotDir,"/BOXPLOTS/plots_boxplots_", parseGrp, ".Rmd"), append=TRUE)
		cat(paste0('##', parseGrp, '\n'), file=paste0(writePlotDir,"/BOXPLOTS/plots_boxplots_", parseGrp, ".Rmd"), append=TRUE)
        }
}

## -----------------------------------------------------------------------
# Data prep

# Get needed geo info
ncid <- ncdf4::nc_open(geoFile)
geoDX <- ncdf4::ncatt_get(ncid,varid=0,'DX')$value
ncdf4::nc_close(ncid)
hydDX <- geoDX/aggfact


# Setup data lists
if (accprecipPlot | flowswePlot | flowlsmPlot) {
	modLdasout_BAS <- list(native=subset(modLdasout[["native"]], modLdasout[["native"]]$fileGroup=="ldasout.basgeo"),
                       snoday=subset(modLdasout[["snoday"]], modLdasout[["snoday"]]$fileGroup=="ldasout.basgeo"),
                       utcday=subset(modLdasout[["utcday"]], modLdasout[["utcday"]]$fileGroup=="ldasout.basgeo"))
}

if (swePlot) {
	modLdasout_SNO <- list(native=subset(modLdasout[["native"]], modLdasout[["native"]]$fileGroup=="ldasout.sno"),
                       snoday=subset(modLdasout[["snoday"]], modLdasout[["snoday"]]$fileGroup=="ldasout.sno"),
                       utcday=subset(modLdasout[["utcday"]], modLdasout[["utcday"]]$fileGroup=="ldasout.sno"))
}


## -----------------------------------------------------------------------
# Generate Plots

# Accumulated Flow
if (accflowPlot) {
message("Generating accumulated flow plots...")
# Setup
accflowList <- list()
if (is.null(accflowTags)) accflowTags <- unique(modFrxstout$tag)
for (i in 1:length(accflowTags)) {
        accflowList[[i]] <- subset(modFrxstout, modFrxstout$tag==accflowTags[i])
}
accflowColors <- lineColors[1:length(accflowList)]
accflowTypes <- rep(lineTyp, length(accflowList))
accflowWidths <- rep(lineWd, length(accflowList))
# Loop plots
for (n in names(gage2basinList)) {
        png(paste0(writePlotDir, "/HYDRO/accstrflow_", n, ".png"), width=2100, height=1350, res=225)
        PlotAccFlow(n, modDfs=accflowList,
                        obs=obsStrData.dy,
                        stdate=accflowStartDate,
                        enddate=accflowEndDate,
                        labMods=accflowTags,
                        labObs="Observed",
                        lnCols=accflowColors,
                        lnTyps=accflowTypes, lnWds=accflowWidths,
                        labTitle=paste0("Accumulated Flow: ", n, " (", obsStrMeta$site_name[obsStrMeta$site_no==n], ")"), obsCol="cumqvol_mm")
        dev.off()
}
if (writeHtml) {
	cat('## Accumulated Flow Plots\n', file=paste0(writePlotDir,"/HYDRO/plots_hydro_", parseGrp, ".Rmd"), append=TRUE)
	for (n in names(gage2basinList)) {
		cat(paste0("```{r accflow_", n, ", fig.width = 12, fig.height = 6, out.width='700', out.height='350', echo=FALSE}\n"), 
			file=paste0(writePlotDir,"/HYDRO/plots_hydro_", parseGrp, ".Rmd"), append=TRUE)
		plottxt <- knitr::knit_expand(text='PlotAccFlow("{{n}}", modDfs=accflowList,
                        obs=obsStrData.dy,
                        stdate=accflowStartDate,
                        enddate=accflowEndDate,
                        labMods=accflowTags,
                        labObs="Observed",
                        lnCols=accflowColors,
                        lnTyps=accflowTypes, lnWds=accflowWidths,
                        labTitle=paste0("Accumulated Flow: ", "{{n}}", " (", obsStrMeta$site_name[obsStrMeta$site_no=="{{n}}"], ")"), obsCol="cumqvol_mm")\n')
		cat(plottxt, file=paste0(writePlotDir,"/HYDRO/plots_hydro_", parseGrp, ".Rmd"), append=TRUE)
		cat('```\n', file=paste0(writePlotDir,"/HYDRO/plots_hydro_", parseGrp, ".Rmd"), append=TRUE)
	}
}
}

# Hydrographs
if (hydroPlot) {
message("Generating hydrograph plots...")
# Setup
hydroList <- list()
if (reachRting) {
	if (is.null(hydroTags)) hydroTags <- unique(modChrtout$tag)
        if (exists("gageList")) {
        	gageNames <- unique(gageList$link)
        } else {
                gageNames <- unique(obsStrData$link)
        }
        idCol <- "link"
} else {
	if (is.null(hydroTags)) hydroTags <- unique(modFrxstout$tag)
        gageNames <- names(gage2basinList)
        idCol <- "site_no"
}
for (i in 1:length(hydroTags)) {
	if (!reachRting) {
		if (hydroPlotDaily) {
        		hydroList[[i]] <- subset(modFrxstout.d, modFrxstout.d$tag==hydroTags[i])
		} else {
			hydroList[[i]] <- subset(modFrxstout, modFrxstout$tag==hydroTags[i])
		}
	} else {
		if (hydroPlotDaily) {
			hydroList[[i]] <- subset(modChrtout.d, modChrtout.d$tag==hydroTags[i])
		} else {
			hydroList[[i]] <- subset(modChrtout, modChrtout$tag==hydroTags[i])
		}
	}
}
hydroColors <- lineColors[1:length(hydroList)]
hydroTypes <- rep(lineTyp, length(hydroList))
hydroWidths <- rep(lineWd, length(hydroList))
# Loop plots
for (n in gageNames) {
	if (idCol == "site_no") {
		siteId <- n
		plotTitle <- paste0("Streamflow: ", n, " (", obsStrMeta$site_name[obsStrMeta$site_no==n], ")")
	} else if (idCol =="link") {
		siteId <- subset(rtLinks$site_no, rtLinks$link==n)
		plotTitle <- paste0("Streamflow: ", subset(rtLinks$site_no, rtLinks$link==n), 
			" (", obsStrMeta$site_name[obsStrMeta$site_no==subset(rtLinks$site_no, rtLinks$link==n)], ")")
	}
        png(paste0(writePlotDir, "/HYDRO/hydrogr_", siteId, ".png"), width=2100, height=1350, res=225)
        PlotFlow(n, modDfs=hydroList,
                        obs=obsStrData,
                        labMods=hydroTags,
                        labObs="Observed",
                        lnCols=hydroColors,
                        lnWds=hydroWidths,
                        labTitle=plotTitle,
                        stdate=hydroStartDate, enddate=hydroEndDate, 
			obsCol="q_cms", idCol=idCol, ymaxPerc=0.99)
        dev.off()
}
if (writeHtml) {
        cat('## Hydrographs\n', file=paste0(writePlotDir,"/HYDRO/plots_hydro_", parseGrp, ".Rmd"), append=TRUE)
        for (n in gageNames) {
                cat(paste0("```{r hydro_", n, ", fig.width = 12, fig.height = 6, out.width='700', out.height='350', echo=FALSE}\n"), 
			file=paste0(writePlotDir,"/HYDRO/plots_hydro_", parseGrp, ".Rmd"), append=TRUE)
                plottxt <- knitr::knit_expand(text='plotTitle <- paste0("Streamflow: ", subset(rtLinks$site_no, rtLinks$link=={{n}}),   
                        " (", obsStrMeta$site_name[obsStrMeta$site_no==subset(rtLinks$site_no, rtLinks$link=={{n}})], ")");
			PlotFlow("{{n}}", modDfs=hydroList,
                        obs=obsStrData,
                        labMods=hydroTags,
                        labObs="Observed",
                        lnCols=hydroColors,
                        lnWds=hydroWidths,
                        labTitle=plotTitle,
                        stdate=hydroStartDate, enddate=hydroEndDate, 
			obsCol="q_cms", idCol=idCol, ymaxPerc=0.99)\n')
                cat(plottxt, file=paste0(writePlotDir,"/HYDRO/plots_hydro_", parseGrp, ".Rmd"), append=TRUE)
                cat('```\n', file=paste0(writePlotDir,"/HYDRO/plots_hydro_", parseGrp, ".Rmd"), append=TRUE)
        }
}
}

# Accumulated Precip
if (accprecipPlot) {
message("Generating accumulated precip plots...")
# Setup
accprecipList <- list()
if (is.null(accprecipTags)) accprecipTags <- unique(modLdasout_BAS[["native"]]$tag)
for (i in 1:length(accprecipTags)) {
        accprecipList[[i]] <- subset(modLdasout_BAS[["native"]], modLdasout_BAS[["native"]]$tag==accprecipTags[i])
}
accprecipColors <- lineColors[1:length(accprecipList)]
accprecipTypes <- rep(lineTyp, length(accprecipList))
accprecipWidths <- rep(lineWd, length(accprecipList))
# Loop plots
for (n in names(gage2basinList)) {
        png(paste0(writePlotDir, "/CLIM/accprecip_", n, ".png"), width=2100, height=1350, res=225)
        PlotAccPrecip(n, modDfs=accprecipList,
                        stdate=accprecipStartDate,
                        enddate=accprecipEndDate,
                        labMods=accprecipTags,
                        lnCols=accprecipColors,
                        lnTyps=accprecipTypes, lnWds=accprecipWidths,
                        labTitle=paste0("Accumulated Precip: ", n, " (", obsStrMeta$site_name[obsStrMeta$site_no==n], ")"))
        dev.off()
}
if (writeHtml) {
        cat('## Accumulated Basin-Mean Precip Plots\n', file=paste0(writePlotDir,"/CLIM/plots_climate_", parseGrp, ".Rmd"), append=TRUE)
        for (n in names(gage2basinList)) {
                cat(paste0("```{r accprecip_", n, ", fig.width = 12, fig.height = 6, out.width='700', out.height='350', echo=FALSE}\n"), 
				file=paste0(writePlotDir,"/CLIM/plots_climate_", parseGrp, ".Rmd"), append=TRUE)
                plottxt <- knitr::knit_expand(text='PlotAccPrecip("{{n}}", modDfs=accprecipList,
                        stdate=accprecipStartDate,
                        enddate=accprecipEndDate,
                        labMods=accprecipTags,
                        lnCols=accprecipColors,
                        lnTyps=accprecipTypes, lnWds=accprecipWidths,
                        labTitle=paste0("Accumulated Precip: ", "{{n}}", " (", obsStrMeta$site_name[obsStrMeta$site_no=="{{n}}"], ")"))\n')
                cat(plottxt, file=paste0(writePlotDir,"/CLIM/plots_climate_", parseGrp, ".Rmd"), append=TRUE)
                cat('```\n', file=paste0(writePlotDir,"/CLIM/plots_climate_", parseGrp, ".Rmd"), append=TRUE)
        }
}
}

# Flow and basin mean SWE
if (flowswePlot) {
message("Generating flow + basin SWE plots...")
# Setup
flowsweStrList <- list()
flowsweLsmList <- list()
if (is.null(flowsweTags)) flowsweTags <- unique(modLdasout_BAS[["native"]]$tag)
for (i in 1:length(flowsweTags)) {
        flowsweStrList[[i]] <- subset(modFrxstout, modFrxstout$tag==flowsweTags[i])
        flowsweLsmList[[i]] <- subset(modLdasout_BAS[["native"]], modLdasout_BAS[["native"]]$tag==flowsweTags[i])
}
flowsweColors <- lineColors[1:length(flowsweStrList)]
flowsweTypes <- rep(lineTyp, length(flowsweStrList))
flowsweWidths <- rep(lineWd, length(flowsweStrList))
# Loop plots
for (n in names(gage2basinList)) {
        png(paste0(writePlotDir, "/SNOW/flowswe_", n, ".png"), width=2100, height=1350, res=225)
        PlotFlowSwe(n, modDfs=flowsweStrList,
                        lsmDfs=flowsweLsmList,
                        obs=obsStrData.dy,
                        labMods=flowsweTags,
                        lnCols=flowsweColors,
                        lnWds=flowsweWidths,
                        labTitle=paste0("Streamflow with Basin-Mean SWE: ", n),
                        stdate=flowsweStartDate, enddate=flowsweEndDate)
        dev.off()
}
if (writeHtml) {
        cat('## Streamflow & Basin-mean SWE Plots\n', file=paste0(writePlotDir,"/SNOW/plots_snow_", parseGrp, ".Rmd"), append=TRUE)
        for (n in names(gage2basinList)) {
                cat(paste0("```{r flowswe_", n, ", fig.width = 12, fig.height = 6, out.width='700', out.height='350', echo=FALSE}\n"), 
			file=paste0(writePlotDir,"/SNOW/plots_snow_", parseGrp, ".Rmd"), append=TRUE)
                plottxt <- knitr::knit_expand(text='PlotFlowSwe("{{n}}", modDfs=flowsweStrList,
                        lsmDfs=flowsweLsmList,
                        obs=obsStrData.dy,
                        labMods=flowsweTags,
                        lnCols=flowsweColors,
                        lnWds=flowsweWidths,
                        labTitle=paste0("Streamflow with Basin-Mean SWE: ", "{{n}}"),
                        stdate=flowsweStartDate, enddate=flowsweEndDate)\n')
                cat(plottxt, file=paste0(writePlotDir,"/SNOW/plots_snow_", parseGrp, ".Rmd"), append=TRUE)
                cat('```\n', file=paste0(writePlotDir,"/SNOW/plots_snow_", parseGrp, ".Rmd"), append=TRUE)
        }
}
}

# Flow and basin mean LSM runoff
if (flowlsmPlot) {
message("Generating flow + basin runoff plots...")
# Setup
if (is.null(flowlsmTags)) {
        if (exists("modFrxstout")) {
                flowlsmTags <- unique(modFrxstout$tag)
                gageNames <- names(gage2basinList)
                idCol <- "site_no"
        } else if (exists("modChrtout")) {
                flowlsmTags <- unique(modChrtout$tag)
                if (!is.null(gageList)) {
                        gageNames <- unique(gageList$link)
                } else {
                        gageNames <- unique(obsStrData$link)
                }
                idCol <- "link"
        }
}
flowlsmColors <- lineColors[1:2]
flowlsmTypes <- rep(lineTyp, 2)
flowlsmWidths <- rep(lineWd, 2)
# Loop plots
for (i in 1:length(flowlsmTags)) {
        if (exists("modFrxstout")) {
                strDf <- subset(modFrxstout, modFrxstout$tag==flowlsmTags[i])
        } else if (exists("modChrtout")) {
                strDf <- subset(modChrtout, modChrtout$tag==flowlsmTags[i])
        } else {
                stop()
        }
        lsmDf <- subset(modLdasout_BAS[["native"]], modLdasout_BAS[["native"]]$tag==flowlsmTags[i])
	ts <- as.integer(difftime(lsmDf$POSIXct[2],lsmDf$POSIXct[1], units="secs"))
	for (n in gageNames) {
		ngage <- ifelse(idCol=="link", as.integer(subset(gageList$site_no, gageList$link==n)), n)
		ngageChar <- ifelse(nchar(as.character(ngage))==7, paste0("0", as.character(ngage)), as.character(ngage))
        	png(paste0(writePlotDir, "/HYDRO/flowlsm_", n, ".png"), width=2100, height=1350, res=225)
        	PlotFlowLsm(n, modDf=strDf, lsmDf=lsmDf, 
                        obs=obsStrData,
                        labMods=flowlsmTags,
                        labObs="Observed",
                        lnCols=flowlsmColors,
                        lnWds=flowlsmWidths,
                        labTitle=paste0("Streamflow: ", ngageChar, " (", obsStrMeta$site_name[obsStrMeta$site_no==ngageChar], ")"),
                        stdate=flowlsmStartDate, enddate=flowlsmEndDate, obsCol="q_cms", idCol=idCol,
			tsSecs=ts, areaSqKm=mskgeo.areaList[[as.character(ngage)]]*geoDX/1000, ngage=ngage)
        	dev.off()
	}
}
if (writeHtml) {
for (i in 1:length(flowlsmTags)) {
        if (exists("modFrxstout")) {
                strDf <- subset(modFrxstout, modFrxstout$tag==flowlsmTags[i])
        } else if (exists("modChrtout")) {
                strDf <- subset(modChrtout, modChrtout$tag==flowlsmTags[i])
        } else {
                stop()
        }
        lsmDf <- subset(modLdasout_BAS[["native"]], modLdasout_BAS[["native"]]$tag==flowlsmTags[i])
        ts <- as.integer(difftime(lsmDf$POSIXct[2],lsmDf$POSIXct[1], units="secs"))
        cat(paste0('## Streamflow & Basin-mean LSM Runoff Plots:', i, '\n'), file=paste0(writePlotDir,"/HYDRO/plots_hydro_", parseGrp, ".Rmd"), append=TRUE)
        for (n in gageNames) {
                cat(paste0("```{r flowlsm_", n, ", fig.width = 12, fig.height = 6, out.width='700', out.height='350', echo=FALSE}\n"),
                        file=paste0(writePlotDir,"/HYDRO/plots_hydro_", parseGrp, ".Rmd"), append=TRUE)
                plottxt <- knitr::knit_expand(text='ngage <- ifelse(idCol=="link", 
					as.integer(subset(gageList$site_no, gageList$link=="{{n}}")), "{{n}}");
			ngageChar <- ifelse(nchar(as.character(ngage))==7, paste0("0", as.character(ngage)), as.character(ngage));
			PlotFlowLsm("{{n}}", modDf=strDf,
			lsmDf=lsmDf,
                        obs=obsStrData,
                        labMods=flowlsmTags,
                        labObs="Observed",
                        lnCols=flowlsmColors,
                        lnWds=flowlsmWidths,
                        labTitle=paste0("Streamflow: ", ngageChar, " (", obsStrMeta$site_name[obsStrMeta$site_no==ngageChar], ")"),
                        stdate=flowlsmStartDate, enddate=flowlsmEndDate, obsCol="q_cms", idCol=idCol,
			tsSecs=ts, areaSqKm=mskgeo.areaList[[as.character(ngage)]]*geoDX/1000, ngage=ngage)\n')
                cat(plottxt, file=paste0(writePlotDir,"/HYDRO/plots_hydro_", parseGrp, ".Rmd"), append=TRUE)
                cat('```\n', file=paste0(writePlotDir,"/HYDRO/plots_hydro_", parseGrp, ".Rmd"), append=TRUE)
        }
}
}
}
   
# SWE
if (swePlot) {
message("Generating SWE plots...")
# Setup
sweList <- list()
if (is.null(sweTags)) sweTags <- unique(modLdasout_SNO[["native"]]$tag)
for (i in 1:length(sweTags)) {
        sweList[[i]] <- subset(modLdasout_SNO[["native"]], modLdasout_SNO[["native"]]$tag==sweTags[i])
}
sweColors <- lineColors[1:length(sweList)]
sweTypes <- rep(lineTyp, length(sweList))
sweWidths <- rep(lineWd, length(sweList))
# Loop plots
sites <- unique(modLdasout_SNO[["native"]]$statArg)
for (n in sites) {
  png(paste0(writePlotDir, "/SNOW/swe_", n, ".png"), width=2100, height=1350, res=225)
  PlotSwe(n, modDfs=sweList,
                obs=obsSnoData, obsmeta=obsSnoMeta,
                labMods=sweTags,
                lnCols=sweColors,
                lnWds=sweWidths,
                precCol.obs="CumPrec_mm", precCol.mod="ACCPRCP",
                sweCol.obs="SWE_mm", sweCol.mod="SNEQV", fact=1, snowh=FALSE,
                labTitle="Accumulated Precipitation and SWE",
                stdate=sweStartDate, enddate=sweEndDate)
  dev.off()
}
if (writeHtml) {
        cat('## Station SWE Plots\n', file=paste0(writePlotDir,"/SNOW/plots_snow_", parseGrp, ".Rmd"), append=TRUE)
        for (n in sites) {
                cat(paste0("```{r swe_", n, ", fig.width = 12, fig.height = 6, out.width='700', out.height='350', echo=FALSE}\n"), 
			file=paste0(writePlotDir,"/SNOW/plots_snow_", parseGrp, ".Rmd"), append=TRUE)
                plottxt <- knitr::knit_expand(text='PlotSwe("{{n}}", modDfs=sweList,
                	obs=obsSnoData, obsmeta=obsSnoMeta,
                	labMods=sweTags,
                	lnCols=sweColors,
                	lnWds=sweWidths,
                	precCol.obs="CumPrec_mm", precCol.mod="ACCPRCP",
                	sweCol.obs="SWE_mm", sweCol.mod="SNEQV", fact=1, snowh=FALSE,
                	labTitle="Accumulated Precipitation and SWE",
                	stdate=sweStartDate, enddate=sweEndDate)\n')
                	cat(plottxt, file=paste0(writePlotDir,"/SNOW/plots_snow_", parseGrp, ".Rmd"), append=TRUE)
                	cat('```\n', file=paste0(writePlotDir,"/SNOW/plots_snow_", parseGrp, ".Rmd"), append=TRUE)
        }
}
}

# MET
if (metPlot) {
	message("Generating MET plots...")
	# Setup
	modLdasin_MET <- subset(modLdasin[["utcday"]], modLdasin[["utcday"]]$fileGroup=="ldasin.met")
	if (is.null(metTags)) metTags <- unique(modLdasin_MET$tag)
	metSites <- unique(modLdasin_MET$statArg)
	for (i in metTags) {
		modLdasin_MET_TAG <- subset(modLdasin_MET, modLdasin_MET$tag==i)
		# Loop Sites
		for (n in metSites) {
  			print(n)
  			# Temperature
  			png(paste0(writePlotDir, "/CLIM/met_temp_", i, "_", n, ".png"), width=1350, height=2100, res=225)
  			PlotMet(obs=obsMetData.dy,
                        	mod=modLdasin_MET_TAG,
                        	site=n,
                        	obsVars=c("Tmean_K", "Tmax_K", "Tmin_K"),
                        	modVars=c("T2D_mean", "T2D_max", "T2D_min"),
                        	lnLabs=c("Mean Temp (C)", "Max Temp (C)", "Min Temp (C)"),
                        	title=paste0(obsMetMeta$site_name[obsMetMeta$site_id==n], ":\nDaily Temperature"),
                        	xLab="", adj=(-273.15), 
				stdate=metStartDate, enddate=metEndDate)
  			dev.off()
			if (writeHtml) {
			        cat('## MET Station Plots\n', file=paste0(writePlotDir,"/CLIM/plots_climate_", parseGrp, ".Rmd"), append=TRUE)
                		cat(paste0("```{r met_", i, "_", n, ", fig.width = 13.5, fig.height = 21, out.width='1350', out.height='2100', echo=FALSE}\n"),
                        		file=paste0(writePlotDir,"/CLIM/plots_climate_", parseGrp, ".Rmd"), append=TRUE)
                		plottxt <- knitr::knit_expand(text='PlotMet(obs=obsMetData.dy,
                                	mod=modLdasin_MET_TAG,
                                	site="{{n}}",
                                	obsVars=c("Tmean_K", "Tmax_K", "Tmin_K"),
                                	modVars=c("T2D_mean", "T2D_max", "T2D_min"),
                                	lnLabs=c("Mean Temp (C)", "Max Temp (C)", "Min Temp (C)"),
                                	title=paste0(obsMetMeta$site_name[obsMetMeta$site_id=="{{n}}"], ":\nDaily Temperature"),
                                	xLab="", adj=(-273.15), 
                                	stdate=metStartDate, enddate=metEndDate)\n')
                        	cat(plottxt, file=paste0(writePlotDir,"/CLIM/plots_climate_", parseGrp, ".Rmd"), append=TRUE)
                        	cat('```\n', file=paste0(writePlotDir,"/CLIM/plots_climate_", parseGrp, ".Rmd"), append=TRUE)
        		}
  			# SW Radiation
  			png(paste0(writePlotDir, "/CLIM/met_swrad_", i, "_", n, ".png"), width=1350, height=2100, res=225)
  			PlotMet(obs=obsMetData.dy,
                        	mod=modLdasin_MET_TAG,
                        	site=n,
                        	obsVars=c("SWRad_mean", "SWRad_max", "SWRad_min"),
                        	modVars=c("SWDOWN_mean", "SWDOWN_max", "SWDOWN_min"),
                        	lnLabs=c("Mean Rad (W/m2)", "Max Rad (W/m2)", "Min Rad (W/m2)"),
                        	title=paste0(obsMetMeta$site_name[obsMetMeta$site_id==n], ":\nDaily Shortwave Radiation"),
                        	xLab="",
				stdate=metStartDate, enddate=metEndDate)
  			dev.off()
  			# Wind
  			png(paste0(writePlotDir, "/CLIM/met_wind_", i, "_", n, ".png"), width=1350, height=2100, res=225)
  			PlotMet(obs=obsMetData.dy,
                        	mod=modLdasin_MET_TAG,
                        	site=n,
                        	obsVars=c("Wind_mean", "Wind_max", "Wind_min"),
                        	modVars=c("Wind_mean", "Wind_max", "Wind_min"),
                        	lnLabs=c("Mean Speed (m/s)", "Max Speed (m/s)", "Min Speed (m/s)"),
                        	title=paste0(obsMetMeta$site_name[obsMetMeta$site_id==n], ":\nDaily Wind Speed"),
                        	xLab="",
				stdate=metStartDate, enddate=metEndDate)
  			dev.off()
  			# Humidity
  			png(paste0(writePlotDir, "/CLIM/met_relhum_", i, "_", n, ".png"), width=1350, height=2100, res=225)
  			PlotMet(obs=obsMetData.dy,
                        	mod=modLdasin_MET_TAG,
                        	site=n,
                        	obsVars=c("RH_mean", "RH_max", "RH_min"),
                        	modVars=c("RelHum_mean", "RelHum_max", "RelHum_min"),
                        	lnLabs=c("Mean RH (0-1)", "Max RH (0-1)", "Min RH (0-1)"),
                        	title=paste0(obsMetMeta$site_name[obsMetMeta$site_id==n], ":\nRelative Humidity"),
                        	xLab="",
				stdate=metStartDate, enddate=metEndDate)
  			dev.off()
  			# Pressure
  			png(paste0(writePlotDir, "/CLIM/met_press_", i, "_", n, ".png"), width=1350, height=2100, res=225)
  			PlotMet(obs=obsMetData.dy,
                        	mod=modLdasin_MET_TAG,
                        	site=n,
                        	obsVars=c("SurfPressmean_Pa", "SurfPressmax_Pa", "SurfPressmin_Pa"),
                        	modVars=c("PSFC_mean", "PSFC_max", "PSFC_min"),
                        	lnLabs=c("Mean Press (kPa)", "Max Press (kPa)", "Min Press (kPa)"),
                        	title=paste0(obsMetMeta$site_name[obsMetMeta$site_id==n], ":\nSurface Pressure"),
                        	xLab="", mult=0.001,
				stdate=metStartDate, enddate=metEndDate)
  			dev.off()
		}
	}
}

# BOXPLOTS

if (strBiasEcoBoxplot | strCorrEcoBoxplot) {

	message("Generating ECOREG boxplots...")

	# Setup ecoregion data
	stats_str.eco <- plyr::join(stats_str, site2region, by="site_no")
	stats_str.eco$L1_INT <- as.integer(stats_str.eco$NA_L1CODE)
	stats_str.eco$L3FACT <- factor(stats_str.eco$US_L3NAME)
	stats_str.eco$L3FACT <- reorder(stats_str.eco$L3FACT, stats_str.eco$L1_INT)
	stats_str.eco$L1FACT <- factor(stats_str.eco$NA_L1NAME)
	stats_str.eco$L1FACT <- reorder(stats_str.eco$L1FACT, stats_str.eco$L1_INT)

	# Setup subsets
	if (exists("trustThresh") & !is.null(trustThresh)) {
		stats_str.eco <- plyr::join(stats_str.eco, trustGages, by="site_no")
		stats_str.eco <- subset(stats_str.eco, stats_str.eco$fractPerfect > trustThresh | is.na(stats_str.eco$fractPerfect))
	}

        # Setup tags
        if (is.null(strEcoTags)) strEcoTags <- unique(stats_str$tag)
        if (is.null(strEcoSeas)) strEcoSeas <- unique(stats_str$seas)

	# Setup html
        if (writeHtml) {
                cat('## Streamflow Ecoregion Boxplots\n', file=paste0(writePlotDir,"/BOXPLOTS/plots_boxplots_", parseGrp, ".Rmd"), append=TRUE)
                strEco.ggList <- list()
                strEco.tblList <- list()
        }

	stats_str.sd <- data.frame()
	stats_str.medians <- data.frame()
	# Loop
        for (i in strEcoTags) {
                for (j in strEcoSeas) {
			stats_str.eco.sub <- subset(stats_str.eco, stats_str.eco$tag==i & stats_str.eco$seas==j)
			# Subset for completeness threshold
			if (exists("nThresh") & !is.null(nThresh)) {
				stats_str.eco.sub <- subset(stats_str.eco.sub, stats_str.eco.sub$dy_n>=nThresh*max(stats_str.eco.sub$dy_n, na.rm=TRUE))
			}
			# Subset for mimimum flow threshold
			if (exists("qThresh") & !is.null(qThresh)) {
				tmp <- subset(stats_qmean, stats_qmean$typ=="Obs" & stats_qmean$tag==i & stats_qmean$seas==j)
				stats_str.eco.sub <- plyr::join(stats_str.eco.sub, tmp[,c("site_no", "qmean")], by="site_no")
				stats_str.eco.sub$qmean_cfs <- stats_str.eco.sub$qmean/(0.3048^3)
				stats_str.eco.sub <- subset(stats_str.eco.sub, stats_str.eco.sub$qmean>qThresh)
			}
			# Subset for maximum bias threshold
			if (exists("biasThresh") & !is.null(biasThresh)) {
				stats_str.eco.sub <- subset(stats_str.eco.sub, abs(stats_str.eco.sub$dy_bias)<biasThresh)
			}

			# Bias
			if (strBiasEcoBoxplot) {

                                give.n <- function(x){return(c(y=-120, label = length(x))) }
                                gg_dybias.l1 <- ggplot(stats_str.eco.sub, aes(x=L1FACT, y=dy_bias, color=L1FACT)) +
                                        geom_boxplot(outlier.size=1.5) +
                                        theme(axis.text.x=element_text(angle=60, hjust=1, vjust=1, size=8)) +
                                        geom_hline(aes(yintercept=0), col='black', linetype=2) +
                                        coord_cartesian(ylim = c(-130,300)) +
                                        scale_y_continuous(minor_breaks = seq(-100, 300, 50)) +
                                        stat_summary(fun.data = give.n, geom = "text", size=2, color='black') +
                                        theme(axis.title.x = element_blank()) +
                                        ylab("Bias (%)") +
                                        scale_color_brewer(palette = "Set1")
                                ggsave(filename=paste0(writePlotDir, "/BOXPLOTS/boxplot_dybias_ecol1_", i, "_", j, ".png"),
                                        plot=gg_dybias.l1, units="in", width=15, height=6, dpi=300)
                                if (writeHtml) {
                                        strEco.ggList <- c(strEco.ggList, list(gg_dybias.l1))
                                }

				give.n <- function(x){return(c(y=-120, label = length(x))) }
				gg_dybias.l3 <- ggplot(stats_str.eco.sub, aes(x=L3FACT, y=dy_bias, color=L1FACT)) + 
					geom_boxplot(outlier.size=1.5) + 
					theme(axis.text.x=element_text(angle=60, hjust=1, vjust=1, size=8)) + 
					geom_hline(aes(yintercept=0), col='black', linetype=2) + 
					coord_cartesian(ylim = c(-130,500)) + 
					scale_y_continuous(minor_breaks = seq(-100, 500, 50)) + 
					stat_summary(fun.data = give.n, geom = "text", size=2, color='black') + 
					theme(axis.title.x = element_blank()) + 
					ylab("Bias (%)") + 
					scale_color_brewer(palette = "Set1")
				ggsave(filename=paste0(writePlotDir, "/BOXPLOTS/boxplot_dybias_ecol3_", i, "_", j, ".png"), 
					plot=gg_dybias.l3, units="in", width=15, height=6, dpi=300)
                        	if (writeHtml) {
                                	strEco.ggList <- c(strEco.ggList, list(gg_dybias.l3))
				}
			}

			# Correlation
			if (strCorrEcoBoxplot) {

                                give.n <- function(x){return(c(y=0.03, label = length(x))) }
                                gg_dycor.l1 <- ggplot(stats_str.eco.sub, aes(x=L1FACT, y=dy_cor, color=L1FACT)) +
                                        geom_boxplot(outlier.size=1) +
                                        theme(axis.text.x=element_text(angle=60, hjust=1, vjust=1, size=8)) +
                                        geom_hline(aes(yintercept=0.6), col='black', linetype=2) +
                                        coord_cartesian(ylim = c(0,1)) +
                                        scale_y_continuous(minor_breaks = seq(0, 1, 0.1), breaks = seq(0, 1, 0.2)) +
                                        stat_summary(fun.data = give.n, geom = "text", size=2, color='black') +
                                        theme(axis.title.x = element_blank()) +
                                        ylab("Daily Correlation") +
                                        scale_color_brewer(palette = "Set1")
                                ggsave(filename=paste0(writePlotDir, "/BOXPLOTS/boxplot_dycor_ecol1_", i, "_", j, ".png"),
                                        plot=gg_dycor.l1, units="in", width=15, height=6, dpi=300)
                                if (writeHtml) {
                                        strEco.ggList <- c(strEco.ggList, list(gg_dycor.l1))
                                }

				give.n <- function(x){return(c(y=0.03, label = length(x))) }
				gg_dycor.l3 <- ggplot(stats_str.eco.sub, aes(x=L3FACT, y=dy_cor, color=L1FACT)) + 
					geom_boxplot(outlier.size=1) + 
					theme(axis.text.x=element_text(angle=60, hjust=1, vjust=1, size=8)) + 
					geom_hline(aes(yintercept=0.6), col='black', linetype=2) + 
					coord_cartesian(ylim = c(0,1)) + 
					scale_y_continuous(minor_breaks = seq(0, 1, 0.1), breaks = seq(0, 1, 0.2)) + 
					stat_summary(fun.data = give.n, geom = "text", size=2, color='black') + 
					theme(axis.title.x = element_blank()) + 
					ylab("Daily Correlation") + 
					scale_color_brewer(palette = "Set1")
				ggsave(filename=paste0(writePlotDir, "/BOXPLOTS/boxplot_dycor_ecol3_", i, "_", j, ".png"), 
					plot=gg_dycor.l3, units="in", width=15, height=6, dpi=300)
                                if (writeHtml) {
                                        strEco.ggList <- c(strEco.ggList, list(gg_dycor.l3))
                                }
			}

			# Min 10 Bias
			if (strBiasMin10EcoBoxplot) {

                                give.n <- function(x){return(c(y=-200, label = length(x))) }
                                gg_dybias_min10.l1 <- ggplot(stats_str.eco.sub, aes(x=L1FACT, y=min10_bias, color=L1FACT)) +
                                        geom_boxplot(outlier.size=1.5) +
                                        theme(axis.text.x=element_text(angle=60, hjust=1, vjust=1, size=8)) +
                                        geom_hline(aes(yintercept=0), col='black', linetype=2) +
                                        coord_cartesian(ylim = c(-230,2000)) +
                                        scale_y_continuous(minor_breaks = seq(-200, 2000, 100)) +
                                        stat_summary(fun.data = give.n, geom = "text", size=2, color='black') +
                                        theme(axis.title.x = element_blank()) +
                                        ylab("Bias in Low (10%) Flows (%)") +
                                        scale_color_brewer(palette = "Set1")
                                ggsave(filename=paste0(writePlotDir, "/BOXPLOTS/boxplot_dybias_min10_ecol1_", i, "_", j, ".png"),
                                        plot=gg_dybias_min10.l1, units="in", width=15, height=6, dpi=300)
                                if (writeHtml) {
                                        strEco.ggList <- c(strEco.ggList, list(gg_dybias_min10.l1))
                                }

				give.n <- function(x){return(c(y=-200, label = length(x))) }
				gg_dybias_min10.l3 <- ggplot(stats_str.eco.sub, aes(x=L3FACT, y=min10_bias, color=L1FACT)) + 
					geom_boxplot(outlier.size=1.5) + 
					theme(axis.text.x=element_text(angle=60, hjust=1, vjust=1, size=8)) + 
					geom_hline(aes(yintercept=0), col='black', linetype=2) + 
					coord_cartesian(ylim = c(-230,2000)) + 
					scale_y_continuous(minor_breaks = seq(-200, 2000, 100)) + 
					stat_summary(fun.data = give.n, geom = "text", size=2, color='black') + 
					theme(axis.title.x = element_blank()) + 
					ylab("Bias in Low (10%) Flows (%)") + 
					scale_color_brewer(palette = "Set1")
				ggsave(filename=paste0(writePlotDir, "/BOXPLOTS/boxplot_dybias_min10_ecol3_", i, "_", j, ".png"), 
					plot=gg_dybias_min10.l3, units="in", width=15, height=6, dpi=300)
                                if (writeHtml) {
                                        strEco.ggList <- c(strEco.ggList, list(gg_dybias_min10.l3))
                                }
			}

                        # Max 10 Bias
                        if (strBiasMax10EcoBoxplot) {

                                give.n <- function(x){return(c(y=-170, label = length(x))) }
                                gg_dybias_max10.l1 <- ggplot(stats_str.eco.sub, aes(x=L1FACT, y=max10_bias, color=L1FACT)) +
                                        geom_boxplot(outlier.size=1.5) +
                                        theme(axis.text.x=element_text(angle=60, hjust=1, vjust=1, size=8)) +
                                        geom_hline(aes(yintercept=0), col='black', linetype=2) +
                                        coord_cartesian(ylim = c(-180,400)) +
                                        scale_y_continuous(minor_breaks = seq(-150, 400, 50)) +
                                        stat_summary(fun.data = give.n, geom = "text", size=2, color='black') +
                                        theme(axis.title.x = element_blank()) +
                                        ylab("Bias in High (10%) Flows (%)") +
                                        scale_color_brewer(palette = "Set1")
                                ggsave(filename=paste0(writePlotDir, "/BOXPLOTS/boxplot_dybias_max10_ecol1_", i, "_", j, ".png"),
                                        plot=gg_dybias_max10.l1, units="in", width=15, height=6, dpi=300)
                                if (writeHtml) {
                                        strEco.ggList <- c(strEco.ggList, list(gg_dybias_max10.l1))
                                }

                                give.n <- function(x){return(c(y=-170, label = length(x))) }
                                gg_dybias_max10.l3 <- ggplot(stats_str.eco.sub, aes(x=L3FACT, y=max10_bias, color=L1FACT)) +
                                        geom_boxplot(outlier.size=1.5) +
                                        theme(axis.text.x=element_text(angle=60, hjust=1, vjust=1, size=8)) +
                                        geom_hline(aes(yintercept=0), col='black', linetype=2) +
                                        coord_cartesian(ylim = c(-180,400)) +
                                        scale_y_continuous(minor_breaks = seq(-150, 400, 50)) +
                                        stat_summary(fun.data = give.n, geom = "text", size=2, color='black') +
                                        theme(axis.title.x = element_blank()) +
                                        ylab("Bias in High (10%) Flows (%)") +
                                        scale_color_brewer(palette = "Set1")
                                ggsave(filename=paste0(writePlotDir, "/BOXPLOTS/boxplot_dybias_max10_ecol3_", i, "_", j, ".png"),
                                        plot=gg_dybias_max10.l3, units="in", width=15, height=6, dpi=300)
                                if (writeHtml) {
                                        strEco.ggList <- c(strEco.ggList, list(gg_dybias_max10.l3))
                                }
                        }

			# Cumulative stats
		        stats_str.sd.sub <- aggregate(stats_str.eco.sub[,1:57], by=list(L3_INT=stats_str.eco.sub$L3_INT, 
                                tag=stats_str.eco.sub$tag, seas=stats_str.eco.sub$seas), sd, na.rm=TRUE)
        		stats_str.medians.sub <- aggregate(stats_str.eco.sub[,1:57], by=list(L3_INT=stats_str.eco.sub$L3_INT,
                                tag=stats_str.eco.sub$tag, seas=stats_str.eco.sub$seas), median, na.rm=TRUE)
			stats_str.sd <- rbind(stats_str.sd, stats_str.sd.sub)
			stats_str.medians <- rbind(stats_str.medians, stats_str.medians.sub)
                        if (writeHtml & statsEcoTables) {
                                        strEco.tblList <- c(strEco.tblList, list(stats_str.sd.sub), list(stats_str.medians.sub))
                        }

		}
	}

        if (writeHtml) {
                for (i in 1:length(strEco.ggList)) {
                        # Plot
                        cat(paste0("```{r strecoboxplot_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE}\n"),
                                file=paste0(writePlotDir,"/BOXPLOTS/plots_boxplots_", parseGrp, ".Rmd"), append=TRUE)
                        plottxt <- knitr::knit_expand(text='strEco.ggList[[{{i}}]]\n')
                        cat(plottxt, file=paste0(writePlotDir,"/BOXPLOTS/plots_boxplots_", parseGrp, ".Rmd"), append=TRUE)
                        cat('```\n', file=paste0(writePlotDir,"/BOXPLOTS/plots_boxplots_", parseGrp, ".Rmd"), append=TRUE)
                        # Table
                        if (statsEcoTables) {
                                cat(paste0("```{r strecoboxplottbl_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE, results='asis'}\n"),
                                file=paste0(writePlotDir,"/BOXPLOTS/plots_boxplots_", parseGrp, ".Rmd"), append=TRUE)
                                #tbltxt <- knitr::knit_expand(text='pandoc.table(strBias.tblList[[{{i}}]], style = "simple", split.table=160)\n')
                                tbltxt <- knitr::knit_expand(text='print(xtable(strEco.tblList[[{{i}}]]), type="html", comment=FALSE)\n')
                                cat(tbltxt, file=paste0(writePlotDir,"/BOXPLOTS/plots_boxplots_", parseGrp, ".Rmd"), append=TRUE)
                                cat('```\n', file=paste0(writePlotDir,"/BOXPLOTS/plots_boxplots_", parseGrp, ".Rmd"), append=TRUE)
                        }
                }
           }

	# Export aggregate stats
	stats_str.sd_tmp <- stats_str.sd
	stats_str.sd_tmp[is.na(stats_str.sd_tmp)]<-(-1e+30)
	write.table(stats_str.sd_tmp, file=paste0(writePlotDir, "/BOXPLOTS/ecol3_statsstr_sd.txt"), sep="\t", row.names=FALSE)

	stats_str.medians_tmp <- stats_str.medians
	stats_str.medians_tmp[is.na(stats_str.medians_tmp)]<-(-1e+30)
	write.table(stats_str.medians_tmp, file=paste0(writePlotDir, "/BOXPLOTS/ecol3_statsstr_medians.txt"), sep="\t", row.names=FALSE)
}


# MAPS

# Initialize for maps
if (strBiasMap | strCorrMap | snosweErrMap | snoprecipErrMap | amfetErrMap | amfetCorrMap) {
	library(ggplot2)
	library(ggmap)
	library(gridExtra)
	if (strBiasMap | strCorrMap) {
		if (reachRting & exists("modChrtout")) {
			modStrout <- modChrtout
		} else if (exists("modFrxstout")) {
			modStrout <- modFrxstout
		}
		# Setup date ranges
		stdate_stats_PRINT <- ifelse(is.null(stdate_stats), ifelse(exists("modStrout"), min(modStrout$POSIXct), ""), stdate_stats)
		enddate_stats_PRINT <- ifelse(is.null(enddate_stats), ifelse(exists("modStrout"), max(modStrout$POSIXct), ""), enddate_stats)
		stdate_stats_sub_PRINT <- ifelse(is.null(stdate_stats_sub), ifelse(exists("modStrout"), min(modStrout$POSIXct), ""), stdate_stats_sub)
		enddate_stats_sub_PRINT <- ifelse(is.null(enddate_stats_sub), ifelse(exists("modStrout"), max(modStrout$POSIXct), ""), enddate_stats_sub)
		statsDateList_STR <- list("Full" = paste0(format(as.POSIXct(stdate_stats_PRINT, origin="1970-01-01 00:00.00 UTC", tz="UTC"), "%Y-%m-%d %H:%M"), 
			" to ", format(as.POSIXct(enddate_stats_PRINT, origin="1970-01-01 00:00.00 UTC", tz="UTC"), "%Y-%m-%d %H:%M")), 
			"Sub" = paste0(format(as.POSIXct(stdate_stats_sub_PRINT, origin="1970-01-01 00:00.00 UTC", tz="UTC"), "%Y-%m-%d %H:%M"), 
                        " to ", format(as.POSIXct(enddate_stats_sub_PRINT, origin="1970-01-01 00:00.00 UTC", tz="UTC"), "%Y-%m-%d %H:%M")))
	}
	if (snosweErrMap | snoprecipErrMap |amfetErrMap | amfetCorrMap) {
                stdate_stats_PRINT <- ifelse(is.null(stdate_stats), ifelse(exists("modLdasout"), min(modLdasout$POSIXct), ""), stdate_stats)
                enddate_stats_PRINT <- ifelse(is.null(enddate_stats), ifelse(exists("modLdasout"), max(modLdasout$POSIXct), ""), enddate_stats)
                stdate_stats_sub_PRINT <- ifelse(is.null(stdate_stats_sub), ifelse(exists("modLdasout"), min(modLdasout$POSIXct), ""), stdate_stats_sub)
                enddate_stats_sub_PRINT <- ifelse(is.null(enddate_stats_sub), ifelse(exists("modLdasout"), max(modLdasout$POSIXct), ""), enddate_stats_sub)
                statsDateList_LDAS <- list("Full" = paste0(format(as.POSIXct(stdate_stats_PRINT, origin="1970-01-01 00:00.00 UTC", tz="UTC"), "%Y-%m-%d %H:%M"),
                        " to ", format(as.POSIXct(enddate_stats_PRINT, origin="1970-01-01 00:00.00 UTC", tz="UTC"), "%Y-%m-%d %H:%M")),
                        "Sub" = paste0(format(as.POSIXct(stdate_stats_sub_PRINT, origin="1970-01-01 00:00.00 UTC", tz="UTC"), "%Y-%m-%d %H:%M"),
                        " to ", format(as.POSIXct(enddate_stats_sub_PRINT, origin="1970-01-01 00:00.00 UTC", tz="UTC"), "%Y-%m-%d %H:%M")))
	}
	# Setup map
	geoMap <- SetupMap(geoFile)
}

# STRFLOW Bias Maps
if (strBiasMap) {
	message("Generating STRFLOW Bias error map...")
	# Setup
	if (is.null(strBiasTags)) strBiasTags <- unique(stats_str$tag)
	if (is.null(strBiasSeas)) strBiasSeas <- unique(stats_str$seas)
        if (plotMapDaily) {colorVar <- "dy_bias"} else {colorVar <- "t_bias"}
        if (plotMapDaily) {exclVar <- "dy_n"} else {exclVar <- "t_n"}
        if (exists("trustThresh") & !is.null(trustThresh)) {
		stats_str <- plyr::join(stats_str, trustGages, by="site_no")
		stats_str <- subset(stats_str, stats_str$fractPerfect > trustThresh | is.na(stats_str$fractPerfect))
	}
	if (writeHtml) {
        	cat('## Streamflow Bias Maps\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
		strBias.ggList <- list()
		strBias.freqList <- list()
		strBias.histList <- list()
		strBias.tblList <- list()
	}
	for (i in strBiasTags) {
        	for (j in strBiasSeas) {
                	tbltmp <- subset(stats_str, stats_str$tag==i & stats_str$seas==j & stats_str$site_no %in% unique(gageList$site_no))
			tbltmp <- plyr::join(tbltmp, subset(stats_qmean, stats_qmean$tag==i & stats_qmean$seas==j & stats_qmean$typ=="Obs"), by="site_no")
			if (!is.null(qThresh)) tbltmp <- subset(tbltmp, tbltmp$qmean >= qThresh)
                	gg <- PlotMapErrors(geoMap, tbltmp,
                        	plotTitle="Modeled Streamflow Bias at USGS Gages",
				plotSubTitle=paste0(i, ", ", statsDateList_STR[[j]]),
                        	sizeVar="qmean", colorVar=colorVar,
                        	sizeLab="Mean\nFlowrate\n(cms)", colorLab="Bias (%)",
				minThreshSize=0, maxThreshSize=200,
				minThreshCol=(-100), maxThreshCol=100,
				minPtsize=1.5, maxPtsize=10,
				exclVar=exclVar, exclThresh=nThresh*max(tbltmp[,exclVar]),
				colBreaks=divColRedYelBlu7, 
                        	valBreaks=c(-Inf, -100, -60, -20, 20, 60, 100, Inf))
                	ggplot2::ggsave(filename=paste0(writePlotDir, "/MAPS/str_bias_map_", i, "_", j, ".png"),
                        	plot=gg[[1]], units="in", width=8, height=6, dpi=300)
                        ggplot2::ggsave(filename=paste0(writePlotDir, "/MAPS/str_bias_hist_", i, "_", j, ".png"),
                                plot=gg[[3]], units="in", width=6, height=4, dpi=300)
			if (writeHtml) {
				strBias.ggList <- c(strBias.ggList, list(gg[[1]]))
                                # Convert freq count to percentage
				freqtbl <- data.frame(gg[[2]]/sum(gg[[2]])*100)
				names(freqtbl) <- c("Bias", "Frequency")
				freqtbl$Frequency <- sprintf("%.1f %%", freqtbl$Frequency)
				strBias.freqList <- c(strBias.freqList, list(freqtbl))
                                strBias.histList <- c(strBias.histList, list(gg[[3]]))
				if (statsMapTables) {
					tbltmp <- tbltmp[order(tbltmp$site_no),]
                                	tbltmp <- data.frame(site_no=tbltmp$site_no, lat=tbltmp$lat, lon=tbltmp$lon, n=tbltmp$t_n, 
                                                bias=tbltmp$t_bias, mae=tbltmp$t_mae,
                                                corr=tbltmp$t_cor, daily_corr=tbltmp$dy_cor, mo_corr=tbltmp$mo_cor,
                                                nse=tbltmp$t_nse, daily_nse=tbltmp$dy_nse, mo_nse=tbltmp$mo_nse)
					strBias.tblList <- c(strBias.tblList, list(tbltmp))
				}
			}
		}
	}
        if (writeHtml) {
                for (i in 1:length(strBias.ggList)) {
                        # Map
                        cat(paste0("```{r strbiasmap_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE}\n"),
                                file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        plottxt <- knitr::knit_expand(text='strBias.ggList[[{{i}}]]\n')
                        cat(plottxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        # Freq Histogram
                        cat(paste0("```{r strbiashist_", i, ", fig.width = 6, fig.height = 4, out.width='600', out.height='400', echo=FALSE}\n"),
                                file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
			plottxt <- knitr::knit_expand(text='strBias.histList[[{{i}}]]\n')
                        cat(plottxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
			cat(paste0("```{r strbiasfreq_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE, results='asis'}\n"),
                                file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        freqtxt <- knitr::knit_expand(text='print(xtable(strBias.freqList[[{{i}}]]), type="html", comment=FALSE)\n')
                        cat(freqtxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
			# Table
			if (statsMapTables) {
				cat(paste0("```{r strbiastbl_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE, results='asis'}\n"),
                         	       file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
				#tbltxt <- knitr::knit_expand(text='pandoc.table(strBias.tblList[[{{i}}]], style = "simple", split.table=160)\n')
				tbltxt <- knitr::knit_expand(text='print(xtable(strBias.tblList[[{{i}}]]), type="html", comment=FALSE)\n')
				cat(tbltxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
				cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
			}
		}
           }
}

# STRFLOW Correlation Maps
if (strCorrMap) {
	message("Generating STRFLOW Corr error map...")
	# Setup
	if (is.null(strCorrTags)) strCorrTags <- unique(stats_str$tag)
	if (is.null(strCorrSeas)) strCorrSeas <- unique(stats_str$seas)
        if (plotMapDaily) {colorVar <- "dy_cor"} else {colorVar <- "t_cor"}
        if (plotMapDaily) {exclVar <- "dy_n"} else {exclVar <- "t_n"}
	if (writeHtml) {
        	cat('## Streamflow Correlation Maps\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
        	strCorr.ggList <- list()
		strCorr.freqList <- list()
		strCorr.histList <- list()
		strCorr.tblList <- list()
	}
	for (i in strCorrTags) {
        	for (j in strCorrSeas) {
			tbltmp <- subset(stats_str, stats_str$tag==i & stats_str$seas==j & stats_str$site_no %in% unique(gageList$site_no))
			tbltmp <- plyr::join(tbltmp, subset(stats_qmean, stats_qmean$tag==i & stats_qmean$seas==j & stats_qmean$typ=="Obs"), by="site_no")
			if (!is.null(qThresh)) tbltmp <- subset(tbltmp, tbltmp$qmean >= qThresh)
                	gg <- PlotMapErrors(geoMap, tbltmp,
                        	plotTitle="Modeled Streamflow Correlation at USGS Gages",
				plotSubTitle=paste0(i, ", ", statsDateList_STR[[j]]),
                        	sizeVar="qmean", colorVar=colorVar,
                        	sizeLab="Mean\nFlowrate\n(cms)", colorLab="Daily\nCorrelation",
				colorLow="orange", colorMid="yellow", colorHigh="cyan4",
				minThreshSize=0, maxThreshSize=200,
                                minThreshCol=0, maxThreshCol=1,
				minPtsize=1.5, maxPtsize=10,
				exclVar=exclVar, exclThresh=nThresh*max(tbltmp[,exclVar]),
                                colBreaks=seqColPurp5,
                                valBreaks=c(-1, 0.2, 0.4, 0.6, 0.8, 1.0))
                	ggplot2::ggsave(filename=paste0(writePlotDir, "/MAPS/str_corr_map_", i, "_", j, ".png"),
                        	plot=gg[[1]], units="in", width=8, height=6, dpi=300)
                        ggplot2::ggsave(filename=paste0(writePlotDir, "/MAPS/str_corr_hist_", i, "_", j, ".png"),
                                plot=gg[[3]], units="in", width=6, height=4, dpi=300)
                	if (writeHtml) {
                        	strCorr.ggList <- c(strCorr.ggList, list(gg[[1]]))
                                # Convert freq count to percentage
                                freqtbl <- data.frame(gg[[2]]/sum(gg[[2]])*100)
                                names(freqtbl) <- c("Correlation", "Frequency")
                                freqtbl$Frequency <- sprintf("%.1f %%", freqtbl$Frequency)
				strCorr.freqList <- c(strCorr.freqList, list(freqtbl))
				strCorr.histList <- c(strCorr.histList, list(gg[[3]]))
				if (statsMapTables) {
					tbltmp <- tbltmp[order(tbltmp$site_no),]
                                	tbltmp <- data.frame(site_no=tbltmp$site_no, lat=tbltmp$lat, lon=tbltmp$lon, n=tbltmp$t_n, 
						bias=tbltmp$t_bias, mae=tbltmp$t_mae,
						corr=tbltmp$t_cor, daily_corr=tbltmp$dy_cor, mo_corr=tbltmp$mo_cor, 
						nse=tbltmp$t_nse, daily_nse=tbltmp$dy_nse, mo_nse=tbltmp$mo_nse)
                                	strCorr.tblList <- c(strCorr.tblList, list(tbltmp))
                		}
			}
          	}
   	}
	if (writeHtml) {
        	for (i in 1:length(strCorr.ggList)) {
			# Map
                	cat(paste0("```{r strcorrmap_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE}\n"),
                                file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        plottxt <- knitr::knit_expand(text='strCorr.ggList[[{{i}}]]\n')
                	cat(plottxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        # Histogram
                        cat(paste0("```{r strcorrhist_", i, ", fig.width = 6, fig.height = 4, out.width='600', out.height='400', echo=FALSE}\n"),
                                file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        plottxt <- knitr::knit_expand(text='strCorr.histList[[{{i}}]]\n')
                        cat(plottxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        cat(paste0("```{r strcorrfreq_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE, results='asis'}\n"),
                                file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        freqtxt <- knitr::knit_expand(text='print(xtable(strCorr.freqList[[{{i}}]]), type="html", comment=FALSE)\n')
                        cat(freqtxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
			# Table
			if (statsMapTables) {
                        	cat(paste0("```{r strcorrtbl_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE, results='asis'}\n"),
                                	file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        	#tbltxt <- knitr::knit_expand(text='pandoc.table(strCorr.tblList[[{{i}}]], style = "simple", split.table=160)\n')
				tbltxt <- knitr::knit_expand(text='print(xtable(strCorr.tblList[[{{i}}]]), type="html", comment=FALSE)\n')
                        	cat(tbltxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
				cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
			}
                }
        }
}

# SNOTEL SWE Maps
if (snosweErrMap) {
	message("Generating SNOTEL SWE error map...")
	# Setup
	if (is.null(snosweErrTags)) snosweErrTags <- unique(stats_ldasout_sno$tag)
	if (is.null(snosweErrSeas)) snosweErrSeas <- unique(stats_ldasout_sno$seas)
        if (writeHtml) {
                cat('## SNOTEL SWE Error Maps\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                snosweErr.ggList <- list()
		snosweErr.freqList <- list()
		snosweErr.histList <- list()
		snosweErr.tblList <- list()
        }
	for (i in snosweErrTags) {
		for (j in snosweErrSeas) {
			tbltmp <- subset(stats_ldasout_sno, stats_ldasout_sno$tag==i & stats_ldasout_sno$var=="SWE" & stats_ldasout_sno$seas==j)
			gg <- PlotMapErrors(geoMap, tbltmp,
                        	plotTitle="Modeled SWE Errors at SNOTEL Stations",
				plotSubTitle=paste0(i, ", ", statsDateList_LDAS[[j]]),
                        	sizeVar="t_mae", colorVar="t_bias",
                        	sizeLab="Mean Abs\nError (mm)", colorLab="Bias (%)",
                                minThreshSize=0, maxThreshSize=100,
                                minThreshCol=(-100), maxThreshCol=100,
                                minPtsize=2, maxPtsize=8,
				exclVar="t_n", exclThresh=nThresh*max(stats_ldasout_sno$t_n),
                                colBreaks=divColBluYelRed6,
                                valBreaks=c(-Inf, -25, -10, 10, 25, 100, Inf))
			ggplot2::ggsave(filename=paste0(writePlotDir, "/MAPS/sno_sweerr_map_", i, "_", j, ".png"), 
				plot=gg[[1]], units="in", width=8, height=6, dpi=300)
                        ggplot2::ggsave(filename=paste0(writePlotDir, "/MAPS/sno_sweerr_hist_", i, "_", j, ".png"),
                                plot=gg[[3]], units="in", width=6, height=4, dpi=300)
                        if (writeHtml) {
                                snosweErr.ggList <- c(snosweErr.ggList, list(gg[[1]]))
				snosweErr.freqList <- c(snosweErr.freqList, list(gg[[2]]))
				snosweErr.histList <- c(snosweErr.histList, list(gg[[3]]))
				if (statsMapTables) {
					tbltmp <- tbltmp[order(tbltmp$site_id),]
                                	tbltmp <- data.frame(site_no=tbltmp$site_id, site_name=tbltmp$site_name, 
						lat=tbltmp$lat, lon=tbltmp$lon, n=tbltmp$t_n, 
                                                bias=tbltmp$t_bias, mae=tbltmp$t_mae,
                                                corr=tbltmp$t_cor, daily_corr=tbltmp$dy_cor, mo_corr=tbltmp$mo_cor,
                                                nse=tbltmp$t_nse, daily_nse=tbltmp$dy_nse, mo_nse=tbltmp$mo_nse)
                                	snosweErr.tblList <- c(snosweErr.tblList, list(tbltmp))
				}
                        }
		}
	}
        if (writeHtml) {
                for (i in 1:length(snosweErr.ggList)) {
			# Map
                        cat(paste0("```{r snosweerrmap_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE}\n"),
                                file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        plottxt <- knitr::knit_expand(text='snosweErr.ggList[[{{i}}]]\n')
                        cat(plottxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        # Histogram
                        cat(paste0("```{r snosweerrhist_", i, ", fig.width = 6, fig.height = 4, out.width='600', out.height='400', echo=FALSE}\n"),
                                file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        plottxt <- knitr::knit_expand(text='snosweErr.histList[[{{i}}]]\n')
                        cat(plottxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        # Table
			if (statsMapTables) {
                        	cat(paste0("```{r snosweerrtbl_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE, results='asis'}\n"),
                                file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        	#tbltxt <- knitr::knit_expand(text='pandoc.table(snosweErr.tblList[[{{i}}]], style = "simple", split.table=160)\n')
				tbltxt <- knitr::knit_expand(text='print(xtable(snosweErr.tblList[[{{i}}]]), type="html", comment=FALSE)\n')
                        	cat(tbltxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
				cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
			}
                }
        }
}

# SNOTEL Precip Maps
if (snoprecipErrMap) {
	message("Generating SNOTEL precip error map...")
	# Setup
	if (is.null(snoprecipErrTags)) snoprecipErrTags <- unique(stats_ldasout_sno$tag)
	if (is.null(snoprecipErrSeas)) snoprecipErrSeas <- unique(stats_ldasout_sno$seas)
        if (writeHtml) {
                cat('## SNOTEL Precipitation Error Maps\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                snoprecipErr.ggList <- list()
		snoprecipErr.freqList <- list()
		snoprecipErr.histList <- list()
		snoprecipErr.tblList <- list()
        }
	for (i in snoprecipErrTags) {
        	for (j in snoprecipErrSeas) {
			tbltmp <- subset(stats_ldasout_sno, stats_ldasout_sno$tag==i & stats_ldasout_sno$seas==j & stats_ldasout_sno$var=="Precip")
                	gg <- PlotMapErrors(geoMap, tbltmp,
                        	plotTitle="Modeled Precipitation Errors at SNOTEL Stations",
				plotSubTitle=paste0(i, ", ", statsDateList_LDAS[[j]]),
                        	sizeVar="t_mae", colorVar="t_bias",
                        	sizeLab="Mean Abs\nError (mm)", colorLab="Bias (%)",
                                minThreshSize=0, maxThreshSize=3,
                                minThreshCol=(-100), maxThreshCol=100,
                                minPtsize=2, maxPtsize=8,
				exclVar="t_n", exclThresh=nThresh*max(stats_ldasout_sno$t_n),
                                colBreaks=divColBluYelRed6,
                                valBreaks=c(-Inf, -25, -10, 10, 25, 100, Inf))
                	ggplot2::ggsave(filename=paste0(writePlotDir, "/MAPS/sno_preciperr_map_", i, "_", j, ".png"),
                        	plot=gg[[1]], units="in", width=8, height=6, dpi=300)
                        ggplot2::ggsave(filename=paste0(writePlotDir, "/MAPS/sno_preciperr_hist_", i, "_", j, ".png"),
                                plot=gg[[3]], units="in", width=6, height=4, dpi=300)
                        if (writeHtml) {
                                snoprecipErr.ggList <- c(snoprecipErr.ggList, list(gg[[1]]))
				snoprecipErr.freqList <- c(snoprecipErr.freqList, list(gg[[2]]))
				snoprecipErr.histList <- c(snoprecipErr.histList, list(gg[[3]]))
				if (statsMapTables) {
					tbltmp <- tbltmp[order(tbltmp$site_id),]
                                	tbltmp <- data.frame(site_no=tbltmp$site_id, site_name=tbltmp$site_name, 
                                                lat=tbltmp$lat, lon=tbltmp$lon, n=tbltmp$t_n,
                                                bias=tbltmp$t_bias, mae=tbltmp$t_mae,
                                                corr=tbltmp$t_cor, daily_corr=tbltmp$dy_cor, mo_corr=tbltmp$mo_cor,
                                                nse=tbltmp$t_nse, daily_nse=tbltmp$dy_nse, mo_nse=tbltmp$mo_nse)
                                	snoprecipErr.tblList <- c(snoprecipErr.tblList, list(tbltmp))
				}
                        }
                }
        }
        if (writeHtml) {
                for (i in 1:length(snoprecipErr.ggList)) {
			# Map
                        cat(paste0("```{r snopreciperrmap_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE}\n"),
                                file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        plottxt <- knitr::knit_expand(text='snoprecipErr.ggList[[{{i}}]]\n')
			cat(plottxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        # Histogram
                        cat(paste0("```{r snopreciperrhist_", i, ", fig.width = 6, fig.height = 4, out.width='600', out.height='400', echo=FALSE}\n"),
                                file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        plottxt <- knitr::knit_expand(text='snoprecipErr.histList[[{{i}}]]\n')
                        cat(plottxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        # Table
			if (statsMapTables) {
                        	cat(paste0("```{r snopreciperrtbl_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE, results='asis'}\n"),
                                	file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        	#tbltxt <- knitr::knit_expand(text='pandoc.table(snoprecipErr.tblList[[{{i}}]], style = "simple", split.table=160)\n')
				tbltxt <- knitr::knit_expand(text='print(xtable(snoprecipErr.tblList[[{{i}}]]), type="html", comment=FALSE)\n')
                        	cat(tbltxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
				cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
			}
                }
        }
}

# AMERIFLUX ET Maps
if (amfetErrMap) {
        message("Generating AMERIFLUX ET error map...")
        # Setup
        if (is.null(amfetErrTags)) amfetErrTags <- unique(stats_ldasout_amf$tag)
        if (is.null(amfetErrSeas)) amfetErrSeas <- unique(stats_ldasout_amf$seas)
        if (writeHtml) {
                cat('## AMERIFLUX ET Error Maps\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                amfetErr.ggList <- list()
		amfetErr.freqList <- list()
		amfetErr.histList <- list()
                amfetErr.tblList <- list()
        }
        for (i in amfetErrTags) {
                for (j in amfetErrSeas) {
			tbltmp <- subset(stats_ldasout_amf, stats_ldasout_amf$tag==i & stats_ldasout_amf$seas==j & stats_ldasout_amf$var=="ET")
                        gg <- PlotMapErrors(geoMap, tbltmp,
                                plotTitle="Modeled ET Errors at Ameriflux Stations",
                                plotSubTitle=paste0(i, ", ", statsDateList_LDAS[[j]]),
                                sizeVar="dy_mae", colorVar="dy_bias",
                                sizeLab="Mean Daily\nAbsolute\nError\n(mm)", colorLab="Bias (%)",
                                minThreshSize=0, maxThreshSize=3,
                                minThreshCol=(-100), maxThreshCol=100,
                                minPtsize=2, maxPtsize=8,
                                exclVar="dy_n", exclThresh=100,
                                colBreaks=divColBluYelRed6,
                                valBreaks=c(-Inf, -25, -10, 10, 25, 100, Inf))
                        ggplot2::ggsave(filename=paste0(writePlotDir, "/MAPS/amf_eterr_map_", i, "_", j, ".png"),
                                plot=gg[[1]], units="in", width=8, height=6, dpi=300)
                        ggplot2::ggsave(filename=paste0(writePlotDir, "/MAPS/amf_eterr_hist_", i, "_", j, ".png"),
                                plot=gg[[3]], units="in", width=6, height=4, dpi=300)
                        if (writeHtml) {
                                amfetErr.ggList <- c(amfetErr.ggList, list(gg[[1]]))
				amfetErr.freqList <- c(amfetErr.freqList, list(gg[[2]]))
				amfetErr.histList <- c(amfetErr.histList, list(gg[[3]]))
				if (statsMapTables) {
                                	tbltmp <- tbltmp[order(tbltmp$site_id),]
                                	tbltmp <- data.frame(site_no=tbltmp$site_id, site_name=tbltmp$name,
                                                lat=tbltmp$lat, lon=tbltmp$lon, n=tbltmp$t_n,
                                                bias=tbltmp$t_bias, mae=tbltmp$t_mae,
                                                corr=tbltmp$t_cor, daily_corr=tbltmp$dy_cor, mo_corr=tbltmp$mo_cor,
                                                nse=tbltmp$t_nse, daily_nse=tbltmp$dy_nse, mo_nse=tbltmp$mo_nse)
                                	amfetErr.tblList <- c(amfetErr.tblList, list(tbltmp))
				}
                        }
                }
        }
        if (writeHtml) {
                for (i in 1:length(amfetErr.ggList)) {
                        # Map
                        cat(paste0("```{r amfeterrmap_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE}\n"),
                                file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        plottxt <- knitr::knit_expand(text='amfetErr.ggList[[{{i}}]]\n')
                        cat(plottxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        # Histogram
                        cat(paste0("```{r amfeterrhist_", i, ", fig.width = 6, fig.height = 4, out.width='600', out.height='400', echo=FALSE}\n"),
                                file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        plottxt <- knitr::knit_expand(text='amfetErr.histList[[{{i}}]]\n')
                        cat(plottxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        # Table
			if (statsMapTables) {
                        	cat(paste0("```{r amfeterrtbl_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE, results='asis'}\n"),
                                	file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        	#tbltxt <- knitr::knit_expand(text='pandoc.table(amfetErr.tblList[[{{i}}]], style = "simple", split.table=160)\n')
                        	tbltxt <- knitr::knit_expand(text='print(xtable(amfetErr.tblList[[{{i}}]]), type="html", comment=FALSE)\n')
                        	cat(tbltxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        	cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
			}
                }
        }
}

if (amfetCorrMap) {
        message("Generating AMERIFLUX ET correlation map...")
        # Setup
        if (is.null(amfetCorrTags)) amfetCorrTags <- unique(stats_ldasout_amf$tag)
        if (is.null(amfetCorrSeas)) amfetCorrSeas <- unique(stats_ldasout_amf$seas)
        if (writeHtml) {
                cat('## AMERIFLUX ET Correlation Maps\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                amfetCorr.ggList <- list()
		amfetCorr.freqList <- list()
		amfetCorr.histList <- list()
                amfetCorr.tblList <- list()
        }
        for (i in amfetCorrTags) {
                for (j in amfetCorrSeas) {
			tbltmp <- subset(stats_ldasout_amf, stats_ldasout_amf$tag==i & stats_ldasout_amf$seas==j & stats_ldasout_amf$var=="ET")
                        gg <- PlotMapErrors(geoMap, tbltmp,
                                plotTitle="Modeled ET Correlation at Ameriflux Stations",
                                plotSubTitle=paste0(i, ", ", statsDateList_LDAS[[j]]),
                                sizeVar="dy_cor", colorVar="dy_cor",
                                sizeLab="Daily\nCorrelation", colorLab="Daily\nCorrelation",
                                colorLow="orange", colorMid="yellow", colorHigh="cyan4",
                                minThreshSize=0, maxThreshSize=1,
                                minThreshCol=0, maxThreshCol=1,
                                minPtsize=0.5, maxPtsize=6,
                                exclVar="dy_n", exclThresh=100,
                                colBreaks=seqColPurp5,
                                valBreaks=c(-1, 0.2, 0.4, 0.6, 0.8, 1.0))
			ggplot2::ggsave(filename=paste0(writePlotDir, "/MAPS/amf_etcorr_map_", i, "_", j, ".png"),
                                plot=gg[[1]], units="in", width=8, height=6, dpi=300)
                        ggplot2::ggsave(filename=paste0(writePlotDir, "/MAPS/amf_etcorr_hist_", i, "_", j, ".png"),
                                plot=gg[[3]], units="in", width=6, height=4, dpi=300)
                        if (writeHtml) {
                                amfetCorr.ggList <- c(amfetCorr.ggList, list(gg[[1]]))
				amfetCorr.freqList <- c(amfetCorr.freqList, list(gg[[2]]))
				amfetCorr.histList <- c(amfetCorr.histList, list(gg[[3]]))
				if (statsMapTables) {
                                	tbltmp <- tbltmp[order(tbltmp$site_id),]
                                	tbltmp <- data.frame(site_no=tbltmp$site_id, site_name=tbltmp$name,
                                                lat=tbltmp$lat, lon=tbltmp$lon, n=tbltmp$t_n,
                                                bias=tbltmp$t_bias, mae=tbltmp$t_mae,
                                                corr=tbltmp$t_cor, daily_corr=tbltmp$dy_cor, mo_corr=tbltmp$mo_cor,
                                                nse=tbltmp$t_nse, daily_nse=tbltmp$dy_nse, mo_nse=tbltmp$mo_nse)
                                	amfetCorr.tblList <- c(amfetCorr.tblList, list(tbltmp))
				}
                        }
                }
        }
        if (writeHtml) {
                for (i in 1:length(amfetCorr.ggList)) {
                        # Map
                        cat(paste0("```{r amfetcorrmap_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE}\n"),
                                file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        plottxt <- knitr::knit_expand(text='amfetCorr.ggList[[{{i}}]]\n')
                        cat(plottxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        # Histogram
                        cat(paste0("```{r amfetcorrhist_", i, ", fig.width = 6, fig.height = 4, out.width='600', out.height='400', echo=FALSE}\n"),
                                file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        plottxt <- knitr::knit_expand(text='amfetCorr.histList[[{{i}}]]\n')
                        cat(plottxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        # Table
			if (statsMapTables) {
                        	cat(paste0("```{r amfetcorrtbl_", i, ", fig.width = 8, fig.height = 6, out.width='800', out.height='600', echo=FALSE, results='asis'}\n"),
                                	file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        	#tbltxt <- knitr::knit_expand(text='pandoc.table(amfetCorr.tblList[[{{i}}]], style = "simple", split.table=160)\n')
                        	tbltxt <- knitr::knit_expand(text='print(xtable(amfetCorr.tblList[[{{i}}]]), type="html", comment=FALSE)\n')
                        	cat(tbltxt, file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
                        	cat('```\n', file=paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), append=TRUE)
			}
                }
        }
}


# Output HTML
if (writeHtml) {
	if (accflowPlot | hydroPlot | flowlsmPlot) {
		knit2html(paste0(writePlotDir,"/HYDRO/plots_hydro_", parseGrp, ".Rmd"), paste0(writePlotDir,"/HYDRO/plots_hydro_", parseGrp, ".html"))
		file.remove("plots_hydro.md")
	}
	if (accprecipPlot) {	
		knit2html(paste0(writePlotDir,"/CLIM/plots_climate_", parseGrp, ".Rmd"), paste0(writePlotDir,"/CLIM/plots_climate_", parseGrp, ".html"))
		file.remove("plots_climate.md")
	}
	if (flowswePlot | swePlot) {	
		knit2html(paste0(writePlotDir,"/SNOW/plots_snow_", parseGrp, ".Rmd"), paste0(writePlotDir,"/SNOW/plots_snow_", parseGrp, ".html"))
		file.remove("plots_snow.md")
	}
	if (strBiasMap | strCorrMap | snosweErrMap | snoprecipErrMap | amfetErrMap | amfetCorrMap) {	
		knit2html(paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".Rmd"), paste0(writePlotDir,"/MAPS/plots_stats_", parseGrp, ".html"))
		file.remove("plots_stats.md")
	}
        if ( strBiasEcoBoxplot | strCorrEcoBoxplot | strBiasEcoBoxplot ) {
                knit2html(paste0(writePlotDir,"/BOXPLOTS/plots_boxplots_", parseGrp, ".Rmd"), paste0(writePlotDir,"/BOXPLOTS/plots_boxplots_", parseGrp, ".html"))
                file.remove("plots_boxplots.md")
        }
	unlink("figure", recursive=TRUE)
}

}

## ------------------------------------------------------------------------
# Cleanup

if (ncores > 1) {stopCluster(cl)}
#save.image("rscripts_working.Rdata")
proc.time()    

