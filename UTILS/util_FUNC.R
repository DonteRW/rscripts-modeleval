# Basin mean function
basin_avg <- function(myvar, mskvar, minValid=-1e+29) {
   myvar[which(myvar<minValid)]<-NA
   sum(mskvar*myvar, na.rm=TRUE)/sum(mskvar, na.rm=TRUE)
 }

# Filename to date conversion functions
rt2dt <- function(x) {as.POSIXct(unlist(strsplit(x, "[.]"))[1], format="%Y%m%d%H%M", tz="UTC")}
ldas2dt <- function(x) {as.POSIXct(unlist(strsplit(x, "[.]"))[1], format="%Y%m%d%H", tz="UTC")}

# Subset file list by dates
subDates <- function(filesList, startDate, endDate, func) {
	tmp <- basename(filesList)
	tmp <- as.POSIXct(apply(as.data.frame(tmp), 1, func), origin='1970-01-01 00:00.00 UTC', tz="UTC")              
	if (!is.null(startDate) & !is.null(endDate)) {
		filesList <- filesList[tmp >= startDate & tmp <= endDate]
        } else if (!is.null(startDate) & is.null(endDate)) {
		filesList <- filesList[tmp >= startDate]
	} else if (is.null(startDate) & !is.null(endDate)) {
                filesList <- filesList[tmp <= endDate]
	}
	filesList
}

# Subset object by dates
subDf <- function(df, stdate=NULL, enddate=NULL) {
  # Subset
  if (!is.null(stdate) & !is.null(enddate)) {
    df <- subset(df, df$POSIXct >= stdate & df$POSIXct <= enddate)
  }
  if (!is.null(stdate) & is.null(enddate)) {
    df <- subset(df, df$POSIXct >= stdate)
  }
  if (is.null(stdate) & !is.null(enddate)) {
    df <- subset(df, df$POSIXct <= enddate)
  }
  df
}

# Remap data
remapData <- function(inObj, mapObj) {
first <- TRUE
for (i in names(mapObj)) {
	if ( mapObj[[i]] %in% names(inObj) ) {
		if (is.data.table(inObj)) {
			out <- data.table(x=inObj[,mapObj[[i]], with=FALSE])
                        names(out) <- i
                        if(first) {outDf <- out} else {outDf <- cbind(outDf, out)}
                        first <- FALSE
		} else {
			out <- data.frame(x=inObj[,mapObj[[i]]], stringsAsFactors=FALSE)
			names(out) <- i
			if(first) {outDf <- out} else {outDf <- cbind(outDf, out)}
			first <- FALSE
		}
	}
}
outDf
}
