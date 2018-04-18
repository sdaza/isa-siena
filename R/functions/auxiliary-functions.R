###################################
# auxiliary functions
###################################

igraphNetworkExtraction <- function(i, data, sims, period, groupName, varName){
   dimsOfDepVar<- attr(data[[groupName]]$depvars[[varName]], "netdims")
   missings <- is.na(data[[groupName]]$depvars[[varName]][,,period]) |
               is.na(data[[groupName]]$depvars[[varName]][,,period+1])
   if (is.null(i)) {
 # sienaGOF wants the observation:
     original <- data[[groupName]]$depvars[[varName]][,,period+1]
     original[missings] <- 0
     returnValue <- graph.adjacency(original)
   }
   else
   {
     missings <- graph.adjacency(missings)
 #sienaGOF wants the i-th simulation:
     returnValue <- graph.difference(
     graph.empty(dimsOfDepVar) +
         edges(t(sims[[i]][[groupName]][[varName]][[period]][,1:2])),
              missings)
   }
   returnValue
 }

CliqueCensus <- function (i, obsData, sims, period, groupName, varName, levls = 1:5){
     x <- networkExtraction(i, obsData, sims, period, groupName, varName)
     cc0 <- sna::clique.census(x, mode='graph', tabulate.by.vertex = FALSE,
                               enumerate=FALSE)[[1]]
     cc <- 0*levls
     names(cc) <- as.character(levls)
     levels.used <- as.numeric(intersect(names(cc0), names(cc)))
     cc[levels.used] <- cc0[levels.used]
     cc
 }

EigenvalueDistribution <- function (i, data, sims, period, groupName, varName,
                         levls=c(seq(0,1,by=0.125)), cumulative=TRUE){
   require(igraph)
   x <- igraphNetworkExtraction(i, data, sims, period, groupName, varName)
   a <- igraph::evcent(x)$vector
   a[is.na(a)] <- Inf
   lel <- length(levls)
   if (cumulative)
   {
     cdi <- sapply(2:lel, function(i){sum(a<=levls[i])})
   }
   else
   {
     cdi <- sapply(2:lel, function(i){
                   sum(a<=levls[i]) - sum(a <= levls[i-1])})
   }
   names(cdi) <- as.character(levls[2:lel])
   cdi
  }

TriadCensus <- function(i, data, sims, wave, groupName, varName, levls=1:16){
     x <- networkExtraction(i, data, sims, wave, groupName, varName)
   if (network::network.edgecount(x) <= 0){x <- symmetrize(x)}
     # because else triad.census(x) will lead to an error
     tc <- sna::triad.census(x)[1,levls]
     # names are transferred automatically
     tc
 }

 CliqueCensus<-function (i, obsData, sims, period, groupName, varName, levls = 1:5){
       require(sna)
       x <- networkExtraction(i, obsData, sims, period, groupName, varName)
       cc0 <- sna::clique.census(x, mode='graph', tabulate.by.vertex = FALSE,
                                 enumerate=FALSE)[[1]]
       cc <- 0*levls
       names(cc) <- as.character(levls)
       levels.used <- as.numeric(intersect(names(cc0), names(cc)))
       cc[levels.used] <- cc0[levels.used]
       cc
   }

GeodesicDistribution <- function (i, data, sims, period, groupName,
                         varName, levls=c(1:5,Inf), cumulative=TRUE, ...) {
   x <- networkExtraction(i, data, sims, period, groupName, varName)
   a <- sna::geodist(symmetrize(x))$gdist
   if (cumulative)
   {
     gdi <- sapply(levls, function(i){ sum(a<=i) })
   }
 else
   {
     gdi <- sapply(levls, function(i){ sum(a==i) })
   }
   names(gdi) <- as.character(levls)
   gdi
 }

MoranGeary <- function(i, data, sims, wave, groupName, varName, levls=1:2){
x <- network::as.sociomatrix(networkExtraction(i, data, sims, wave, groupName, varName[1]))
z <- behaviorExtraction(i,data,sims,wave,groupName,varName[2])
n <- length(z)
z.ave <- mean(z,na.rm=TRUE)
numerator <- n*sum(x*outer(z-z.ave,z-z.ave),na.rm=TRUE)
denominator <- sum(x,na.rm=TRUE)*sum((z-z.ave)^2,na.rm=TRUE)
res <- numerator/denominator
numerator <- (n-1)*sum(x*(outer(z,z,FUN='-')^2),na.rm=TRUE)
denominator <- 2*sum(x,na.rm=TRUE)*sum((z-z.ave)^2,na.rm=TRUE)
res[2] <- numerator/denominator
names(res) <- c("Moran","Geary")
return(res)
}

MoranGOF <- function(i, data, sims, wave, groupName, varName, levls=1:2){
x <- network::as.sociomatrix(networkExtraction(i, data, sims, wave, groupName, varName[1]))
z <- behaviorExtraction(i,data,sims,wave,groupName,varName[2])
n <- length(z)
z.ave <- mean(z,na.rm=TRUE)
numerator <- n*sum(x*outer(z-z.ave,z-z.ave),na.rm=TRUE)
denominator <- sum(x,na.rm=TRUE)*sum((z-z.ave)^2,na.rm=TRUE)
res <- numerator/denominator
names(res) <- c("Moran")
return(res)
}

Moran <- function(x, mynetwork) {
net <- network::as.sociomatrix(mynetwork)
x.ave <- mean(x, na.rm = TRUE)
n <- length(x)
numerator <- n * sum ( net * outer(x - x.ave, x - x.ave), na.rm = TRUE )
denominator <- sum(net, na.rm = TRUE) * sum ((x - x.ave)^2, na.rm = TRUE)
res <- numerator / denominator
return(res)
}

Geary <- function(x, mynetwork) {
net <- network::as.sociomatrix(mynetwork)
x.ave <- mean(x, na.rm = TRUE)
n <- length(x)
numerator <- (n - 1) * sum(net * (outer(x, x, FUN = '-')^2 ), na.rm = TRUE)
denominator <- 2 * sum(net, na.rm=TRUE) * sum((x - x.ave)^2 ,na.rm=TRUE)
res <- numerator / denominator
return(res)
}

# get moran and geary measures
getAutoCorr <- function(model, period, varnames, groupName = "Data1",
                      observed = FALSE, numsim = NULL) {
output <- list()
nsim <- 1:length(model$sims)
if (!is.null(numsim)) {
  if (numsim > 0) nsim <- sample(nsim, numsim)
}

for (j in seq_along(period)) {
  print(paste0("Computing Period = ", j,"..."))
  moran <- NA; geary <- NA
    if (observed) {
      net <- networkExtraction(NULL, model$f, model$sims,
                               period = period[j], groupName= groupName,
                               varName= varnames[1])
      beh <- behaviorExtraction(NULL, model$f, model$sims,
                                period = period[j], groupName= groupName,
                                varName= varnames[2])
      moran <-  Moran(beh, net)
      geary <-  Geary(beh, net)
    }
    else if (observed == FALSE) {
      for (i in seq_along(nsim)) {
        net <- networkExtraction(nsim[i], model$f, model$sims,
                                 period = period[j], groupName= groupName,
                                 varName= varnames[1])
        beh <- behaviorExtraction(nsim[i], model$f, model$sims,
                                 period = period[j], groupName= groupName,
                                 varName= varnames[2])
         moran[i] <-  Moran(beh, net)
         geary[i] <-  Geary(beh, net)

      }
    }
  output[[j]] <- data.frame(moran, geary)
  result <- data.table::rbindlist(output, idcol = "time")
}
return(result)
}

# convergence function
siena07ToConvergence <- function(alg, dat, eff, ans0 = NULL, nodes = 1, ...) {

numr <- 0
ans <- siena07(alg, data = dat, effects = eff, prevAns = ans0, batch = TRUE,
               verbose = FALSE, useCluster = TRUE, nbrNodes = nodes,
               returnDeps = TRUE, silent = FALSE) # the first run

repeat {
numr <- numr + 1
tm <- ans$tconv.max
cat(numr, tm,"\n")
if (tm < 0.25 & numr > 0) {break}
if (tm > 8) {break}
# count number of repeated runs
# convergence indicator
# report how far we are
# success
# divergence without much hope
# of returning to good parameter values
if (numr > 10) {break}  # now it has lasted too long
ans <- siena07(alg, data = dat, effects = eff, prevAns = ans, batch = TRUE,
               verbose = FALSE, useCluster = TRUE, nbrNodes = nodes, returnDeps = TRUE,
               silent = FALSE)
}
 if (tm > 0.25)
 {
    cat("Warning: convergence inadequate.\n")
 }

ans

}