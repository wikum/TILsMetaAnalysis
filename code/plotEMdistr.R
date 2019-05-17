###########################################################################
#### A function to plot a normal distribution from EM-algorithm results
#### Arguments: the expression vector, the mclust output, and the sample name
### It returns split points between distributions as side effect

plotEMdistr <- function(exprDat, emOut, name,
			upper, 
			quant=c(0.25, 0.75),
                        forXlab="Log2 expression",
                        forTitle="Log2 expression", plot=TRUE){
	## Require mclust
	require(mclust)

	## Get rid of NA
	exprDat <- na.omit(exprDat)

	## Compute density for all data
	densExpDat <- density(exprDat)

	## Compute density for normals
	mySeq <- seq(min(exprDat, na.rm=T), max(exprDat, na.rm=T), length=1000)
	means <- emOut$parameters$mean
	sds <- sqrt(emOut$parameters$variance$sigmasq)
	props <- emOut$parameters$pro
	ppp <- mapply(mean=means, sd=sds, props=props, MoreArgs=list(x=mySeq),
		      FUN= function(x, mean, sd, props) { dnorm(x, mean, sd) * props })

	## Define y-axis limits
	if (missing(upper)) {
		upper <- max(c(densExpDat$y, ppp))
	}

    ## To define split line:
	if (ncol(ppp) > 1) {
		sel <- which( diff(apply(ppp, 1, which.max)) > 0 )
		mySplitLines <- mySeq[sel]
		## If only one distribution
	} else {
		mySplitLines <- quantile(exprDat, quant)
	}

	
	if(plot){
	## Plot
	plot(densExpDat, lwd=2, ylim=c(0, upper), las=0,
	     xlab=paste(name, forXlab),
	     main=paste(name, forTitle))

	## Add background and grid
	rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey90")
	grid(col="white")

	## Replot
	lines(densExpDat, lwd=2, ylim=c(0, upper))

	## Add lines
	require(RColorBrewer)

	## cols <- c("blue", "red")
	cols <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Set2"))[1:ncol(ppp)]

	## print(dim(ppp))
	sapply(1:ncol(ppp), function(i, x, y, cols, ...) {
		lines(x , y[,i], col=cols[i], lwd=2, lty=4)
		abline(v=mySplitLines, lty=2,
		       ## col=rgb(0.75, 0.4, 0.5, 0.75))
		       col=rgb(0.5, 0.5, 0.5, 0.9))
	}, y=ppp, x=mySeq, cols=cols, splitLine)
	
	}
	## Return split lines
	
	mySplitLines
}


