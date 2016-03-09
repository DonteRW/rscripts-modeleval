###################################################
##  Main Control Script to process model output  ##
###################################################

################## General Setup ##################

## Number of cores to use? Must have the doParallel package installed
ncores <- 8

## Specify the high-resolution routing domain file
hydFile <- '/glade/p/ral/RHAP/gochis/WRF_Hydro_code/WRF-Hydro_NCEP_test_version_Oct_27_2015/DOMAIN/Fulldom_hires_netcdf_file_1km.nc'

## Specify the low-resolution geogrid file
geoFile <- '/d6/adugger/WRF_Hydro/CONUS_3km/DOMAIN/geo_em.d01.nc.conus_3km_nlcd11'

## Specify the aggregation factor between hydrogrid and geogrid
aggfact <- 1

## Specify location of .Rdata file containing pre-processed mask objects
maskFile <- '/d6/adugger/WRF_Hydro/CONUS_IOC/DOMAIN/gagesII_MASKS.Rdata'

## Specify whether the model run used NHD reach-based routing (otherwise gridded routing assumed)
# If TRUE, mask file should also contain rtLinks dataframe.
reachRting <- TRUE

## Temp directory to write intermediate files
tmpDir <- '/d6/adugger/WRF_Hydro/NFIE/ANALYSIS/'

tmplist <- c("05087500", "10329500", "09107000", "09497980", "08340500", "08068780", "08095300", "06218500", "06352000", "06446500", "06441500", "06696980",
             "06888000", "06930000", "05411850", "03456500", "03237500", "07331300", "07291000", "07292500", "02372250", "02315000", "02297155", "02314500",
             "02092500", "02011400", "01532000", "01144000", "04057510", "12013500", "14325000", "13011900", "13339500", "14231000", "11264500", "11138500",
             "07148400", "02147500", "02322000", "02321500", "02322800")

################## Observations ###################

## Path to Ameriflux data .Rdata file
AMFfile <- "/glade/p/ral/RHAP/adugger/CONUS_IOC/OBS/AMF/obs_AMF_1998_current.Rdata" 

## Path to SNOTEl data .Rdata file
SNOfile <- "/glade/p/ral/RHAP/adugger/CONUS_IOC/OBS/SNOTEL/obs_SNOTEL_1998_current.Rdata"

## Path to meteorological station data .Rdata file
METfile <- NULL

## Path to streamflow data .Rdata file
#STRfile <- "/d6/adugger/WRF_Hydro/NFIE/OBS/obs_USGS_gagesII_2015.Rdata"
STRfile <- "/d6/adugger/WRF_Hydro/NFIE/OBS/obsStrData_gagesII_nfie2015_DV.Rdata"
#STRfile <- "/d6/adugger/WRF_Hydro/NFIE/OBS/obsStrData_gagesII_nfie2015_IV.Rdata"


################ Model Output Reads ###############

## Read model output?
readMod <- FALSE

## If TRUE, specify the following to read in model output:

        # Specify the model run output directory or directories
	modPathList <- c('/glade/p/ral/RHAP/adugger/CONUS_IOC/IOC_calib_runs/terr_rtg/v1.1_w_calib_no_oCONUS_no_res_new_route_link_Dec_26_2015')

        # Specify tags to identify the model run or runs (should be 1:1 with number of model output directories)
	modTagList <- c('SPINUP5YR_Full_Routing_v1.1')

        # Specify the output .Rdata file to create
        modReadFileOut <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/ANALYSIS/151228_conus_gagesII_su2010fullrtv11_modelout_STR.Rdata'
        # Append to existing file? FALSE = create new file (or overwrite existing!)
        modAppend <- FALSE

	# Select what aggregations/imports to run:

		# Basin means and imports
		readBasinLdasout <- FALSE
		readBasinRtout <- FALSE
		readGwout <- FALSE
		readFrxstout <- FALSE

		# Channel routing
		readChrtout <- TRUE
			# Read only links with gages?
			readChrtout_GAGES <- FALSE
			# Read specified subset? Provide object with link and site_no columns
			#readLink2gage <- read.table("/glade/p/ral/RHAP/adugger/CONUS_IOC/DOMAIN/link2gage_gagesII.txt", sep="\t", header=TRUE, colClasses=c("integer", "character"))

		# Snotel sites
		readSnoLdasout <- FALSE

		# Ameriflux sites
		readAmfLdasout <- FALSE

		# MET sites
		readMetLdasout <- FALSE

	# Subset LDASOUT variables?
	varsLdasoutSUB <- FALSE
	varsLdasoutNFIE <- FALSE
	varsLdasoutIOC0 <- TRUE

	# Specify start and end dates if you do NOT want to read all files
	readModStart <- NULL
	readModEnd <- NULL


################## Forcing Reads ##################

## Read forcing data?
readForc <- FALSE

## If TRUE, specify the following:

	# Specify the path to the forcing data
	forcPathList <- c('/glade/p/ral/RHAP/gochis/WRF_Hydro_code/WRF-Hydro_NCEP_test_version_Oct_27_2015/forcing/') 

        # Specify tags to identify the forcings (should be 1:1 with number of model forcing directories)
	forcTagList <- c('NLDAS2-Downscaled')

	# Specify the forcing output .Rdata file to create
	forcReadFileOut <- '/glade/p/ral/RHAP/adugger/CONUS_IOC/ANALYSIS/conus_nldas_forcings_2010.Rdata'
        # Append to existing file? FALSE = create new file (or overwrite existing!)
        forcAppend <- FALSE

	# Select what aggregations/imports to run:

		# Basin means
		readBasinLdasin <- FALSE

		# SNOTEL sites
		readSnoLdasin <- TRUE

		# Ameriflux sites
		readAmfLdasin <- TRUE

		# MET sites
		readMetLdasin <- FALSE

        # Specify start and end dates if you do NOT want to read all files
        readForcStart <- as.POSIXct("2010-01-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC")
        readForcEnd <- as.POSIXct("2010-12-31 23:59", format="%Y-%m-%d %H:%M", tz="UTC")


############# Model Performance Stats #############

## Calculate stats?
calcStats <- FALSE

	## Calculate streamflow performance stats?
	strProc <- TRUE
		# Read specified subset? Provide object with link and site_no columns
                statsLink2gage <- read.table("/d6/adugger/WRF_Hydro/NFIE/OBS/link2gage_gagesII_CORRECTED.txt",
                                        sep="\t", header=TRUE, colClasses=c("integer","character"))
                # Calculate daily stats?
                strProcDaily <- TRUE

	## Calculate SNOTEL performance stats?
	snoProc <- FALSE

	## Calculate Ameriflux performance stats?
	amfProc <- FALSE

	## Calculate MET performance stats?
	metProc <- FALSE

## If any are TRUE, specify the following:

	# If the raw data read .Rdata file exists (vs. created above), specify the file
	modReadFileIn <- '/d6/adugger/WRF_Hydro/NFIE/ANALYSIS/160126_nfie_modout_STR_w_ioc_mrms.Rdata'

        # Specify the stats output .Rdata file to create
        statsFileOut <- '/d6/adugger/WRF_Hydro/NFIE/ANALYSIS/160126_nfie_stats_STR_w_ioc_mrms_DVs.Rdata'

	# Range dates for main stats
	stdate_stats <- as.POSIXct("2015-05-23 00:00", format="%Y-%m-%d %H:%M", tz="UTC")
	enddate_stats <- as.POSIXct("2015-09-30 23:59", format="%Y-%m-%d %H:%M", tz="UTC")

	# Range dates for seasonal stats (e.g., spring)
	stdate_stats_sub <- as.POSIXct("2015-06-01 00:00", format="%Y-%m-%d %H:%M", tz="UTC")
	enddate_stats_sub <- as.POSIXct("2015-09-30 23:59", format="%Y-%m-%d %H:%M", tz="UTC")

	# Write stats tables?
	writeStatsFile <- TRUE
	# If TRUE, specify output directory
	writeDir <- '/d6/adugger/WRF_Hydro/NFIE/ANALYSIS/160126_nfie_gagesII_w_ioc_mrms_DV_PLOTS_ALL'



################### Plotting ######################

## Create plots and/or maps?
createPlots <- TRUE

	## Specify .Rdata file containing model time series data for plotting (if needed)
	plotModFile <- '/d6/adugger/WRF_Hydro/NFIE/ANALYSIS/160126_nfie_modout_STR_w_ioc_mrms.Rdata'

	## Specify .Rdata file containing model statistics for plotting (if needed)
	plotStatsFile <- '/d6/adugger/WRF_Hydro/NFIE/ANALYSIS/160126_nfie_stats_STR_w_ioc_mrms_DVs.Rdata'

	## Specify .Rdata file containing model forcings for plotting (if needed)
        plotForcFile <- NULL

	## Create HTML files?
	writeHtml <- TRUE

	## Specify plot output directory
	writePlotDir <- '/d6/adugger/WRF_Hydro/NFIE/ANALYSIS/160126_nfie_gagesII_w_ioc_mrms_DV_PLOTS_ALL'

	## Plot specified subset? Provide object with link and site_no columns. HYDRO PLOTS ONLY!
        plotLink2gage <- read.table("/d6/adugger/WRF_Hydro/NFIE/OBS/link2gage_gagesII_CORRECTED.txt",
                                        sep="\t", header=TRUE, colClasses=c("integer","character"))
	#plotLink2gage <- subset(plotLink2gage, plotLink2gage$site_no %in% tmplist)
	#plotLink2gage <- plotLink2gage[order(plotLink2gage$site_no),]

	######### TIME SERIES PLOTS ###########

	## Generate accumulated flow plots?
	accflowPlot <- FALSE

		# Specify which run tags to plot
		accflowTags <- NULL

		# Specify start date
		accflowStartDate <- as.POSIXct("2014-04-01", format="%Y-%m-%d", tz="UTC")

		# Specify end date
		accflowEndDate <- NULL

	## Generate hydrographs?
	hydroPlot <- TRUE

        	# Specify which run tags to plot
        	hydroTags <- c('NFIE_No_Routing_MRMS', 'IOCv1.1_Full_Routing_MRMS', 'IOCv1.1_No_Routing_MRMS', 'IOCv1.1_No_Routing_NLDAS')

        	# Specify start date
        	hydroStartDate <- as.POSIXct("2015-05-23", format="%Y-%m-%d", tz="UTC")
        
        	# Specify end date
        	hydroEndDate <- as.POSIXct("2015-09-30", format="%Y-%m-%d", tz="UTC")

                # Plot daily values?
                hydroPlotDaily <- TRUE

	## Generate accumulated precip plots?
	accprecipPlot <- FALSE

        	# Specify which run tags to plot
        	accprecipTags <- NULL
        
        	# Specify start date
        	accprecipStartDate <- as.POSIXct("2013-10-01", format="%Y-%m-%d", tz="UTC")
        
        	# Specify end date
        	accprecipEndDate <- NULL

	## Generate Streamflow and Basin-mean SWE plots?
	flowswePlot <- FALSE

        	# Specify which run tags to plot
        	flowsweTags <- NULL
        
        	# Specify start date
        	flowsweStartDate <- as.POSIXct("2013-10-01", format="%Y-%m-%d", tz="UTC")
        
        	# Specify end date
        	flowsweEndDate <- NULL

        ## Generate Streamflow and Basin-mean LSM Runoff plots?
        flowlsmPlot <- FALSE

                # Specify which run tags to plot
                flowlsmTags <- NULL

                # Specify start date
                flowlsmStartDate <- NULL

                # Specify end date
                flowlsmEndDate <- NULL

	## Generate SWE station plots?
	swePlot <- FALSE

        	# Specify which run tags to plot
        	sweTags <- NULL

        	# Specify start date
        	sweStartDate <- as.POSIXct("2013-10-01", format="%Y-%m-%d", tz="UTC")

        	# Specify end date
        	sweEndDate <- NULL

        ## Generate MET station plots?
        metPlot <- FALSE

                # Specify which run tags to plot
                metTags <- NULL

                # Specify start date
                metStartDate <- as.POSIXct("2014-10-01", format="%Y-%m-%d", tz="UTC")

                # Specify end date
                metEndDate <- NULL


	########### MAPS #############

	## Threshold for "completeness" of sites to include. Multiplier on the max n in the set.
	#  Set to 0 to plot all.
	nThresh <- 0.75

	## Trusted gages table (if used, otherwise NULL). 
	#  The table should have a "site_no" column and a "fractPerfect" column with the "trust" metric (0-1).
	trustGages <- read.table("/d6/adugger/WRF_Hydro/NFIE/OBS/trustGages.txt",
                                        sep="\t", header=TRUE, colClasses=c("character", "numeric"))
	# Threshold for "trusted" gages. Set to 0 to plot all.
	trustThresh <- 0.75

	# Plot daily stats?
	plotMapDaily <- TRUE

	## Generate STRFLOW bias maps?
	strBiasMap <- TRUE

        	# Specify which run tags to plot
        	strBiasTags <- NULL

        	# Specify which run seasons to plot
        	strBiasSeas <- NULL

	## Generate STRFLOW correlation maps?
	strCorrMap <- TRUE

        	# Specify which run tags to plot
        	strCorrTags <- NULL

        	# Specify which run seasons to plot
        	strCorrSeas <- NULL

	## Generate SNOTEL SWE error maps?
	snosweErrMap <- FALSE

        	# Specify which run tags to plot
        	snosweErrTags <- NULL

        	# Specify which run seasons to plot
        	snosweErrSeas <- NULL

	## Generate SNOTEL Precip error maps?
	snoprecipErrMap <- FALSE

        	# Specify which run tags to plot
        	snoprecipErrTags <- NULL

        	# Specify which run seasons to plot
        	snoprecipErrSeas <- NULL

        ## Generate SNOTEL SWE error maps?
        amfetErrMap <- FALSE

                # Specify which run tags to plot
                amfetErrTags <- NULL

                # Specify which run seasons to plot
                amfetErrSeas <- NULL

        ## Generate SNOTEL SWE error maps?
        amfetCorrMap <- FALSE
        
                # Specify which run tags to plot
                amfetCorrTags <- NULL

                # Specify which run seasons to plot
                amfetCorrSeas <- NULL

	## Include summary stats tables?
	statsMapTables <- FALSE


###########################################################################################
## RUN (do not change anything below this line)

source("run_CONFIG.R")

