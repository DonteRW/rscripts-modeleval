PlotAccPrecipFlow <- function(n, str1=modFrxstout_wy2015_NLDAS2dwnsc_fullrtng,
                        lsm1=modLdasout_wy2015_NLDAS2dwnsc_fullrtng_BAS,
                        obsstr=obsStr.dy,
                        str2=modFrxstout_wy2015_NLDAS2dwnsc_NSSL_fullrtng,
                        lsm2=modLdasout_wy2015_NLDAS2dwnsc_NSSL_fullrtng_BAS) {
  str1 <- subset(str1, str1$site_no==n)
  lsm1 <- subset(lsm1, lsm1$site_no==n)
  obsstr <- subset(obsstr, obsstr$site_no==n)
  str2 <- subset(str2, str2$site_no==n)
  lsm2 <- subset(lsm2, lsm2$site_no==n)
  with(lsm1, plot(POSIXct, ACCPRCP, typ='l', col='darkmagenta'))
  with(str1, lines(POSIXct, ACCFLOW, col='deepskyblue'))
  with(lsm1, lines(POSIXct, ACCECAN+ACCEDIR+ACCETRAN, col='darkolivegreen3'))
  with(lsm2, lines(POSIXct, ACCPRCP, col='darkmagenta', lty=2))
  with(str2, lines(POSIXct, ACCFLOW, col='deepskyblue', lty=2))
  with(lsm2, lines(POSIXct, ACCECAN+ACCEDIR+ACCETRAN, col='darkolivegreen3', lty=2))
  with(obsstr, lines(POSIXct, cumqvol_mm, col='blue', lwd=2, lty=1))
}


PlotAccPrecip <- function(n, modDfs,
                        labMods=NULL,
                        labTitle="Accumulated Precipitation",
                        lnCols=NULL, lnWds=NULL, lnTyps=NULL,
                        stdate=NULL,
                        enddate=NULL,
                        modCol="ACCPRCP") {
  # Parse type of input for model data (dataframe or list of multiple dataframes)
  if (is.data.frame(modDfs)) {
        str1 <- modDfs
        strcnt <- 1
  } else if (is.list(modDfs)) {
        str1 <- modDfs[[1]]
        strcnt <- length(modDfs)
  } else {
        stop("modDfs must be a dataframe or a list of dataframes")
  }
  if (is.data.table(str1)) str1<-data.frame(str1)
  # Subset by dates
  if (is.null(stdate)) stdate <- min(str1$POSIXct)
  if (is.null(enddate)) enddate <- max(str1$POSIXct)
  str1 <- subset(str1, str1$statArg==n & str1$POSIXct>=stdate & str1$POSIXct<=enddate)
  # Calculate maximum y val for plot limits
  ymax <- max(str1[nrow(str1),modCol]-str1[1,modCol], na.rm=TRUE)
  if (!is.data.frame(modDfs) & is.list(modDfs) & length(modDfs)>1) {
        for (stri in modDfs) {
		if (is.data.table(stri)) stri<-data.frame(stri)
                stri <- subset(stri, stri$statArg==n & stri$POSIXct>=stdate & stri$POSIXct<=enddate)
                ymax <- max(ymax, stri[nrow(stri),modCol]-stri[1,modCol], na.rm=TRUE)
                }
        }
  # Set colors, widths, types
  if (is.null(lnCols)) lnCols <- sample(colours(), strcnt)
  if (is.null(lnWds)) lnWds <- rep(1, strcnt)
  if (is.null(lnTyps)) lnTyps <- rep(1, strcnt)
  # Set labels
  if (is.null(labMods)) labMods <- paste0("Model", 1:strcnt)
  # Create plot
  plot(str1$POSIXct, str1[,modCol]-str1[1,modCol], typ='l', ylim=c(0, ymax),
        xlim=c(stdate, enddate), col=lnCols[1], lty=lnTyps[1], lwd=lnWds[1],
        xlab="", ylab="Accumulated precipitation (mm)", cex.axis=1.2, cex.lab=1.2)
  if (!is.data.frame(modDfs) & is.list(modDfs) & length(modDfs)>1) {
        for (j in 2:length(modDfs)) {
                stri <- modDfs[[j]]
		if (is.data.table(stri)) stri<-data.frame(stri)
                stri <- subset(stri, stri$statArg==n & stri$POSIXct>=stdate & stri$POSIXct<=enddate)
                lines(stri$POSIXct, stri[,modCol]-stri[1,modCol], col=lnCols[j], lty=lnTyps[j], lwd=lnWds[j])
                }
        }
  title(labTitle, cex.main=1.6)
  legend("topleft", c(labMods[1:strcnt]),
                lty=c(lnTyps[1:strcnt]), lwd=c(lnWds[1:strcnt]),
                col=c(lnCols[1:strcnt]), cex=1.2)
}


PlotAccFlow <- function(n, modDfs, obs,
                        labMods=NULL,
                        labObs="Observed",
                        labTitle="Accumulated Flow",
                        lnCols=NULL, lnWds=NULL, lnTyps=NULL,
                        stdate=NULL,
                        enddate=NULL,
                        modCol="ACCFLOW", obsCol="cumqvol_mm") {
  # Parse type of input for model data (dataframe or list of multiple dataframes)
  if (is.data.frame(modDfs)) {
        str1 <- modDfs
        strcnt <- 1
  } else if (is.list(modDfs)) {
        str1 <- modDfs[[1]]
        strcnt <- length(modDfs)
  } else {
        stop("modDfs must be a dataframe or a list of dataframes")
  }
  if (is.data.table(str1)) str1<-data.frame(str1)
  # Subset by dates
  if (is.null(stdate)) stdate <- min(str1$POSIXct)
  if (is.null(enddate)) enddate <- max(str1$POSIXct)
  str1 <- subset(str1, str1$site_no==n & str1$POSIXct>=stdate & str1$POSIXct<=enddate)
  obs <- subset(obs, obs$site_no==n & obs$POSIXct>=stdate & obs$POSIXct<=enddate)
  # Calculate maximum y val for plot limits
  ymax <- max(str1[nrow(str1),modCol]-str1[1,modCol], obs[nrow(obs),obsCol]-obs[1,obsCol], na.rm=TRUE)
  if (!is.data.frame(modDfs) & is.list(modDfs) & length(modDfs)>1) {
        for (stri in modDfs) {
		if (is.data.table(stri)) stri<-data.frame(stri)
                stri <- subset(stri, stri$site_no==n & stri$POSIXct>=stdate & stri$POSIXct<=enddate)
                ymax <- max(ymax, stri[nrow(stri),modCol]-stri[1,modCol], na.rm=TRUE)
                }
        }
  # Set colors, widths, types
  if (is.null(lnCols)) lnCols <- sample(colours(), strcnt)
  if (is.null(lnWds)) lnWds <- rep(1, strcnt)
  if (is.null(lnTyps)) lnTyps <- rep(1, strcnt)
  # Set labels
  if (is.null(labMods)) labMods <- paste0("Model", 1:strcnt)
  # Create plot
  plot(str1$POSIXct, str1[,modCol]-str1[1,modCol], typ='l', ylim=c(0, ymax),
        xlim=c(stdate, enddate), col=lnCols[1], lty=lnTyps[1], lwd=lnWds[1],
        xlab="", ylab="Accumulated flow (mm)", cex.axis=1.2, cex.lab=1.2)
  if (!is.data.frame(modDfs) & is.list(modDfs) & length(modDfs)>1) {
        for (j in 2:length(modDfs)) {
                stri <- modDfs[[j]]
		if (is.data.table(stri)) stri<-data.frame(stri)
                stri <- subset(stri, stri$site_no==n & stri$POSIXct>=stdate & stri$POSIXct<=enddate)
                lines(stri$POSIXct, stri[,modCol]-stri[1,modCol], col=lnCols[j], lty=lnTyps[j], lwd=lnWds[j])
                }
        }
  lines(obs$POSIXct, obs[,obsCol]-obs[1,obsCol], col='black', lwd=2, lty=1)
  title(labTitle, cex.main=1.6)
  legend("topleft", c(labMods[1:strcnt], labObs),
                lty=c(lnTyps[1:strcnt],1), lwd=c(lnWds[1:strcnt],2),
                col=c(lnCols[1:strcnt], 'black'), cex=1.2)
}


PlotFlow <- function(n, modDfs, obs,
                        labMods=NULL,
                        labObs="Observed",
                        labTitle="Streamflow",
                        lnCols=NULL, lnWds=NULL, lnTyps=NULL,
                        stdate=NULL,
                        enddate=NULL,
                        modCol="q_cms", obsCol="mean_qcms",
			idCol="site_no", ymaxPerc=1.0,
			legAlign="topleft") {
  # Parse type of input for model data (dataframe or list of multiple dataframes)
  if (is.data.frame(modDfs)) {
        str1 <- modDfs
        strcnt <- 1
  } else if (is.list(modDfs)) {
        str1 <- modDfs[[1]]
        strcnt <- length(modDfs)
  } else {
        stop("modDfs must be a dataframe or a list of dataframes")
  }
  if (is.data.table(str1) & is.data.table(obs)) {
        # Subset by dates
        if (is.null(stdate)) stdate <- min(str1$POSIXct)
        if (is.null(enddate)) enddate <- max(str1$POSIXct)
        str1 <- str1[get(idCol)==n & POSIXct>=stdate & POSIXct<=enddate,]
        obs <- obs[get(idCol)==as.integer(n) & POSIXct>=stdate & POSIXct<=enddate,]
        # Calculate maximum y val for plot limits
        ymax <- max(quantile(str1[,get(modCol)], ymaxPerc, na.rm=TRUE), quantile(obs[,get(obsCol)], ymaxPerc, na.rm=TRUE), na.rm=TRUE)
        if (!is.data.frame(modDfs) & is.list(modDfs) & length(modDfs)>1) {
                for (stri in modDfs) {
                        stri <- stri[get(idCol)==n & POSIXct>=stdate & POSIXct<=enddate,]
                        ymax <- max(ymax, quantile(stri[,get(modCol)], ymaxPerc, na.rm=TRUE), na.rm=TRUE)
                        }
                }
        # Set colors, widths, types
        if (is.null(lnCols)) lnCols <- sample(colours(), strcnt)
        if (is.null(lnWds)) lnWds <- rep(1, strcnt)
        if (is.null(lnTyps)) lnTyps <- rep(1, strcnt)
        # Set labels
        if (is.null(labMods)) labMods <- paste0("Model", 1:strcnt)
        # Create plot
        if (!is.data.frame(modDfs) & is.list(modDfs) & length(modDfs)>1) {
        	plot(str1$POSIXct, str1[,get(modCol)], typ='l', ylim=c(0, ymax),
                	xlim=c(stdate, enddate), col=lnCols[1], lty=lnTyps[1], lwd=0,
                	xlab="", ylab="Streamflow (m3/s)", cex.axis=1.2, cex.lab=1.2)
                for (j in 1:length(modDfs)) {
                        stri <- modDfs[[j]]
                        stri <- stri[get(idCol)==n & POSIXct>=stdate & POSIXct<=enddate,]
                        lines(stri$POSIXct, stri[,get(modCol)], col=lnCols[j], lty=lnTyps[j], lwd=lnWds[j])
                        }
        } else {
        	plot(str1$POSIXct, str1[,get(modCol)], typ='l', ylim=c(0, ymax),
                	xlim=c(stdate, enddate), col=lnCols[1], lty=lnTyps[1], lwd=lnWds[1],
               		xlab="", ylab="Streamflow (m3/s)", cex.axis=1.2, cex.lab=1.2)
	}
        lines(obs$POSIXct, obs[,get(obsCol)], col='black', lwd=2, lty=1)
        title(labTitle, cex.main=1.6)
        legend(legAlign, c(labMods[1:strcnt], labObs),
                lty=c(lnTyps[1:strcnt],1), lwd=c(lnWds[1:strcnt],2),
                col=c(lnCols[1:strcnt], 'black'), cex=1.2,
                bg="white")
  } else {
  	if (is.data.table(str1)) str1<-data.frame(str1)
  	# Subset by dates
  	if (is.null(stdate)) stdate <- min(str1$POSIXct)
  	if (is.null(enddate)) enddate <- max(str1$POSIXct)
  	str1 <- subset(str1, str1[,idCol]==n & str1$POSIXct>=stdate & str1$POSIXct<=enddate)
  	if (is.data.table(obs)) {
        	obs <- obs[get(idCol)==as.integer(n) & POSIXct>=stdate & POSIXct<=enddate,]
  		obs <- data.frame(obs)
  	} else {
  		obs <- subset(obs, obs[,idCol]==n & obs$POSIXct>=stdate & obs$POSIXct<=enddate)
  	}
  	# Calculate maximum y val for plot limits
  	ymax <- max(quantile(str1[,modCol], ymaxPerc, na.rm=TRUE), quantile(obs[,obsCol], ymaxPerc, na.rm=TRUE), na.rm=TRUE)
  	if (!is.data.frame(modDfs) & is.list(modDfs) & length(modDfs)>1) {
        	for (stri in modDfs) {
			if (is.data.table(stri)) stri<-data.frame(stri)
        		stri <- subset(stri, stri[,idCol]==n & stri$POSIXct>=stdate & stri$POSIXct<=enddate)
			ymax <- max(ymax, quantile(stri[,modCol], ymaxPerc, na.rm=TRUE), na.rm=TRUE)
                	}
        	}
  	# Set colors, widths, types
  	if (is.null(lnCols)) lnCols <- sample(colours(), strcnt)
  	if (is.null(lnWds)) lnWds <- rep(1, strcnt)
  	if (is.null(lnTyps)) lnTyps <- rep(1, strcnt)
  	# Set labels
  	if (is.null(labMods)) labMods <- paste0("Model", 1:strcnt)
  	# Create plot
  	plot(str1$POSIXct, str1[,modCol], typ='l', ylim=c(0, ymax),
        	xlim=c(stdate, enddate), col=lnCols[1], lty=lnTyps[1], lwd=0,
        	xlab="", ylab="Streamflow (m3/s)", cex.axis=1.2, cex.lab=1.2)
  	if (!is.data.frame(modDfs) & is.list(modDfs) & length(modDfs)>1) {
        	for (j in 1:length(modDfs)) {
                	stri <- modDfs[[j]]
			if (is.data.table(stri)) stri<-data.frame(stri)
                	stri <- subset(stri, stri[,idCol]==n & stri$POSIXct>=stdate & stri$POSIXct<=enddate)
                	lines(stri$POSIXct, stri[,modCol], col=lnCols[j], lty=lnTyps[j], lwd=lnWds[j])
                	}
        	}
  	lines(obs$POSIXct, obs[,obsCol], col='black', lwd=2, lty=1)
  	title(labTitle, cex.main=1.6)
  	legend("topleft", c(labMods[1:strcnt], labObs),
                lty=c(lnTyps[1:strcnt],1), lwd=c(lnWds[1:strcnt],2),
                col=c(lnCols[1:strcnt], 'black'), cex=1.2,
                bg="white")
	}
}



PlotFlowSwe <- function(n, modDfs, lsmDfs, obs,
                        labMods=NULL,
                        labObs="Observed",
                        labTitle="Streamflow with Basin-Mean SWE",
                        lnCols=NULL, lnWds=NULL, lnTyps=NULL,
                        stdate=NULL,
                        enddate=NULL,
                        modCol="q_cms", lsmCol="SNEQV", obsCol="mean_qcms") {
  # Parse type of input for model data (dataframe or list of multiple dataframes)
  if (is.data.frame(modDfs)) {
        str1 <- modDfs
        strcnt <- 1
  } else if (is.list(modDfs)) {
        str1 <- modDfs[[1]]
        strcnt <- length(modDfs)
  } else {
        stop("modDfs must be a dataframe or a list of dataframes")
  }
  if (is.data.table(str1)) str1<-data.frame(str1)
  # Parse type of input for model data (dataframe or list of multiple dataframes)
  if (is.data.frame(lsmDfs)) {
        lsm1 <- lsmDfs
        lsmcnt <- 1
  } else if (is.list(lsmDfs)) {
        lsm1 <- lsmDfs[[1]]
        lsmcnt <- length(lsmDfs)
  } else {
        stop("lsmDfs must be a dataframe or a list of dataframes")
  }
  if (is.data.table(lsm1)) lsm1<-data.frame(lsm1)
  # Subset by dates
  if (is.null(stdate)) stdate <- min(str1$POSIXct)
  if (is.null(enddate)) enddate <- max(str1$POSIXct)
  str1 <- subset(str1, str1[idCol]==n & str1$POSIXct>=stdate & str1$POSIXct<=enddate)
  lsm1 <- subset(lsm1, lsm1$statArg==n & lsm1$POSIXct>=stdate & lsm1$POSIXct<=enddate)
  obs <- subset(obs, obs[idCol]==n & obs$POSIXct>=stdate & obs$POSIXct<=enddate)
  # Calculate maximum y val for plot limits
  ymax <- max(str1[,modCol], obs[,obsCol], na.rm=TRUE)
  ymax_lsm <- max(lsm1[,lsmCol], na.rm=TRUE)
  if (!is.data.frame(modDfs) & is.list(modDfs) & length(modDfs)>1) {
        for (stri in modDfs) {
		if (is.data.table(stri)) stri<-data.frame(stri)
                stri <- subset(stri, stri[idCol]==n & stri$POSIXct>=stdate & stri$POSIXct<=enddate)
                ymax <- max(ymax, stri[,modCol], na.rm=TRUE)
                }   
        }   
  if (!is.data.frame(lsmDfs) & is.list(lsmDfs) & length(lsmDfs)>1) {
        for (lsmi in lsmDfs) {
		if (is.data.table(lsmi)) lsmi<-data.frame(lsmi)
                lsmi <- subset(lsmi, lsmi$statArg==n & lsmi$POSIXct>=stdate & lsmi$POSIXct<=enddate)
                ymax_lsm <- max(ymax_lsm, lsmi[,lsmCol], na.rm=TRUE)
                }
        }
  # Set colors, widths, types
  if (is.null(lnCols)) lnCols <- sample(colours(), strcnt)
  if (is.null(lnWds)) lnWds <- rep(1, strcnt)
  if (is.null(lnTyps)) lnTyps <- rep(1, strcnt)
  # Set labels
  if (is.null(labMods)) labMods <- paste0("Model", 1:strcnt)
  # Create plot
  par(mar=c(5,4,4,5)+.1)
  plot(str1$POSIXct, str1[,modCol], typ='l', ylim=c(0, ymax),
        xlim=c(stdate, enddate), col=lnCols[1], lty=lnTyps[1], lwd=lnWds[1],
        xlab="", ylab="Streamflow (m3/s)", cex.axis=1.2, cex.lab=1.2)
  if (!is.data.frame(modDfs) & is.list(modDfs) & length(modDfs)>1) {
        for (j in 2:length(modDfs)) {
                stri <- modDfs[[j]]
		if (is.data.table(stri)) stri<-data.frame(stri)
                stri <- subset(stri, stri[idCol]==n & stri$POSIXct>=stdate & stri$POSIXct<=enddate)
                lines(stri$POSIXct, stri[,modCol], col=lnCols[j], lty=lnTyps[j], lwd=lnWds[j])
                }   
        }   
  lines(obs$POSIXct, obs[,obsCol], col='black', lwd=2, lty=1)
  par(new=TRUE)
  plot(lsm1$POSIXct, lsm1[,lsmCol], col=lnCols[1], lty=2, lwd=lnWds[1], typ='l',
        ylim=c(0, ymax_lsm),
        xaxt="n", yaxt="n", xlab="", ylab="")
  axis(4)
  mtext("SWE (mm)", side=4, line=3)
  if (!is.data.frame(lsmDfs) & is.list(lsmDfs) & length(lsmDfs)>1) {
        for (j in 2:length(lsmDfs)) {
                lsmi <- lsmDfs[[j]]
		if (is.data.table(lsmi)) lsmi<-data.frame(lsmi)
                lsmi <- subset(lsmi, lsmi$statArg==n)
                lines(lsmi$POSIXct, lsmi[,lsmCol], col=lnCols[j], lty=2, lwd=lnWds[j])
                }
        }
  title(labTitle, cex.main=1.6)
  legend("topleft", c(labMods[1:strcnt], labObs, "", "Streamflow (m3/s)", "SWE (mm)"),
                lty=c(lnTyps[1:strcnt],1,1,1,2), lwd=c(lnWds[1:strcnt],2,1,2,2),
                col=c(lnCols[1:strcnt], 'black','white','grey40','grey40'), cex=1.2,
                bg="white")
}


PlotFlowLsm <- function(n, modDf, lsmDf, obs, 
                        labMods=NULL,
                        labObs="Observed",
                        labTitle="Streamflow",
                        lnCols=NULL, lnWds=NULL, lnTyps=NULL,
                        stdate=NULL,
                        enddate=NULL,
                        modCol="q_cms", obsCol="mean_qcms",
                        idCol="site_no", tsSecs=86400, areaSqKm, ngage=NULL) {
  # Parse type of input for model data (dataframe or list of multiple dataframes)
  str1 <- modDf
  lsm1 <- lsmDf
  if (is.data.table(str1)) str1<-data.frame(str1)
  if (is.data.table(lsm1)) lsm1<-data.frame(lsm1)
  # Subset by dates
  if (is.null(stdate)) stdate <- min(str1$POSIXct)
  if (is.null(enddate)) enddate <- max(str1$POSIXct)
  ngage <- ifelse(!is.null(ngage), ngage, n)
  str1 <- subset(str1, str1[idCol]==n & str1$POSIXct>=stdate & str1$POSIXct<=enddate)
  lsm1 <- subset(lsm1, lsm1["statArg"]==ngage & lsm1$POSIXct>=stdate & lsm1$POSIXct<=enddate)
  if (is.data.table(obs)) {
        obs <- obs[get(idCol)==n & POSIXct>=stdate & POSIXct<=enddate,]
        obs <- data.frame(obs)
  } else {
        obs <- subset(obs, obs[idCol]==n & obs$POSIXct>=stdate & obs$POSIXct<=enddate)
  }
  # Calculate maximum y val for plot limits
  ymax <- max(str1[,modCol], obs[,obsCol], na.rm=TRUE)
  # Set colors, widths, types
  if (is.null(lnCols)) lnCols <- sample(colours(), 2)
  if (is.null(lnWds)) lnWds <- rep(1, 2)
  if (is.null(lnTyps)) lnTyps <- rep(1, 2)
  # Set labels
  if (is.null(labMods)) labMods <- paste0("Model", 1)
  # Create plot
  plot(str1$POSIXct, str1[,modCol], typ='l', ylim=c(0, ymax),
        xlim=c(stdate, enddate), col=lnCols[1], lty=lnTyps[1], lwd=lnWds[1],
        xlab="", ylab="Streamflow (m3/s)", cex.axis=1.2, cex.lab=1.2)
  with(lsm1, lines(POSIXct, (DEL_UGDRNOFF+DEL_SFCRNOFF)/tsSecs/1000*areaSqKm*1000*1000, 
	col=lnCols[2], lty=lnTyps[2], lwd=lnWds[2]))
  lines(obs$POSIXct, obs[,obsCol], col='black', lwd=2, lty=1)
  title(labTitle, cex.main=1.6)
  legend("topleft", c(labMods, "LSM Runoff", labObs),
                lty=c(lnTyps, 1), lwd=c(lnWds, 2),
                col=c(lnCols, 'black'), cex=1.2,
                bg="white")
}

