##############################################
# R function for aggregating multiple groups'
# sienaGOF-results.
# v. March 25, 2015
# written by Christian Steglich
##############################################

sienaGroupGOF <- function(goflist,agg.type='tests') {
# goflist should be a list of sienaGOF objects that
#    are obtained from independent groups and that are
#    compatible in the widest sense of the word;
# agg.type can be 'tests' or 'data', see below.
	n.gr <- length(goflist)
	n.sims <- unlist(lapply(1:n.gr,function(i){
		length(goflist[[i]][[1]]$SimulatedTestStat)
	}))
	if (var(rep(n.sims,2))>0) stop(
		'combined groups have different numbers of GOF statistics'
	)
	n.sims <- min(n.sims)
	if (agg.type=='tests') {
	# aggregation of MH distances over independent groups
		obs.mhd <- unlist(lapply(1:n.gr,function(i){
			goflist[[i]][[1]]$ObservedTestStat
		}))
		mhd <- sum(obs.mhd)
		# NOTE that obs.mhd is actually the chi-squared-
		# distributed square of Mahalnobis' distance!
		sim.mhd <- apply(
			do.call(cbind,lapply(1:n.gr,function(i){
				goflist[[i]][[1]]$SimulatedTestStat})),
			MARGIN=1,
			FUN=function(x){sqrt(sum(x^2))}
		)
		p.value <- sum(mhd<=sim.mhd)/n.sims
	} else if (agg.type=='data') {
	# calculation of MH distance after pooling data of groups
		obs.stats <- colSums(
			do.call(rbind,lapply(1:n.gr,function(i){
				goflist[[i]][[1]]$Observations}))
		)
		sim.stats <- apply(
			array(unlist(lapply(1:n.gr,function(i){
					goflist[[i]][[1]]$Simulations
				})),dim=c(n.sims,length(obs.stats),n.gr)),
			MARGIN=c(1,2),
			FUN=sum
		)
		e.stats <- colMeans(sim.stats)
		cen.sim.stats <- scale(sim.stats,scale=FALSE)
		cov.inv <- MASS::ginv(cov(sim.stats))
		MH <- function(x) {x %*% cov.inv %*% x}
		mhd <- c(MH(obs.stats-e.stats))
		sim.mhd <- apply(cen.sim.stats,1,MH)
		p.value <- sum(mhd<=sim.mhd)/n.sims
	} else stop(	paste(
		'agg.type of unknown type',as.character(agg.type)
	))
	return(data.frame(
		p.value=p.value,
		MHD=mhd,
		agg.type=agg.type,
		nGroups=n.gr,
		nSims=n.sims
	))
}

plot.sienaGroupGOF <- function (goflist,center=FALSE,scale=FALSE,
violin=TRUE,key=NULL,perc=.05,period=1,main=NULL,ylab=NULL) {
# plot function for data aggregation version of groupGOF
	n.gr <- length(goflist)
	n.sims <- unlist(lapply(1:n.gr,function(i){
		length(goflist[[i]][[1]]$SimulatedTestStat)

	}))
	if (var(rep(n.sims,2))>0) stop(
		'combined groups have different numbers of GOF statistics'
	)
	n.sims <- min(n.sims)
	obs.stats <- t(colSums(
		do.call(rbind,lapply(1:n.gr,function(i){
			goflist[[i]][[1]]$Observations}))
	))
	sim.stats <- apply(
		array(unlist(lapply(1:n.gr,function(i){
				goflist[[i]][[1]]$Simulations
			})),dim=c(n.sims,length(obs.stats),n.gr)),
		MARGIN=c(1,2),
		FUN=sum
	)
	gofobject <- goflist[[1]] # exemplarily
	gofobject[[1]]$Observations <- obs.stats
	gofobject[[1]]$Simulations <- sim.stats
	gofobject[[1]]$p <- sienaGroupGOF(goflist,agg.type='data')$p.value
	return(plot(gofobject,center = center,scale=scale,
		violin=violin,key=key,perc=perc,
		period=period,main=main,ylab=ylab))
}
