########################################################
# rsiena with measurement update variation using 30 replicates
# author: sebastian daza
# version: 1.00
#######################################################

#+ clean workspace and load libraries
# rm(list=ls(all=TRUE))

# libraries
library(data.table)
library(RSienaTest)
library(texreg)
library(igraph)
library(sna)
library(sdazar)
# library(colorout)

# linux servers
setwd("/home/s/sdaza/00projects/siena/")
source("/home/s/sdaza/00projects/siena/R/functions/auxiliary-functions.R")
results_path <- "R/results/timing/"

# anylogic experiment
# vr_selection = new boolean[] {false, true};
# vr_influence = new boolean[] {true, false};
# vr_behavior = new int[] {0,0}; // uniform
# vr_radius = new boolean[] {true, true};

scenarios <- c("I", "S") # iterations
specifications <- c("5", "10", "20")

iterations <- 1:length(scenarios)
specs <- 1:length(specifications)

iterations
#+ set up values ISA
agents <- 200 # define number of agents to create network
max_rep <- 100 # max number of replications
wi <- 500; le <- 500 # dimensions of space

# define number of cores
# (nodes <- floor(parallel::detectCores() / 3))
nodes <- 5

# create R object
loop <- data.table(expand.grid(iter = iterations, sp = specs))

loop
# auxiliar vector to update waves
v <- c(0,1,2,3)

#+ siena setup
sim <- 1000 # number of simulations

# read data (b, n, d, p)
load("R/data/measurementUpdate.Rdata")

# p # paramaters

# track path
track_path <- paste0(results_path, "track.Rdata")

# track function
keeptrack <- function(track_path) {

  if ( file.exists(track_path) ) {
    load(track_path) } else {
       track <- data.table() }

  # look for files
  mpattern <- paste0("rr_iter_", scen, "_spec_", spec)
  files <- list.files(path = results_path, pattern = mpattern)
  areplicates <- NULL

  if (length(files) > 0) {
  pat <- "(.+)(rep_)([0-9]+)(.*)"
  areplicates <- sort(as.numeric(sapply(files, function(x) sub(pat, "\\3", x))))
  }

  if (length(areplicates) > 0) {
    track <- rbind(track,
    data.table(iter = niter, replicate = areplicates,
      running = 0, finished = 1))
    track <- unique(track)
    save(track, file = track_path)
  }

  # define replicates to be done
  rreplicates <- NULL
  if (nrow(track) > 0 ) {
  rreplicates <- track[running == 1 & finished == 0, replicate]
  }
  sreplicates <- unique(c(areplicates, rreplicates))
  treplicates <- 1:max_rep
  if (is.null(sreplicates)) {
    rep <- treplicates
  } else {
    rep <- treplicates[!treplicates %in% sreplicates] # pending replicates
  }
  save(track, file = track_path)
  return(rep)
}

###########################################################
#+ loop through iterations
###########################################################

for (j in 1:nrow(loop)) {

# get scenario (iteration) and siena specification
scen <- loop[j, iter]
spec <- loop[j, sp]

# define waves based on iteration, this is only thing changing
if (spec == 1) { nwaves <- v * 1 + 36  } # 5
else if  (spec == 2) { nwaves <- v * 2 + 36  } # 10
else if (spec == 3) { nwaves <- v * 5 + 36  } # 20

# adjust accordingly
print(paste0(":::::::::::::::: ", paste0(scen, "-",  spec)))

# name iteration
niter <- paste0(scen, "-", spec)

rep <- keeptrack(track_path)
if (length(rep) == 0) { next }

# loop through replications
for (h in rep) {

# update replicates
rep <- keeptrack(track_path)

# rewrite h
if (length(rep) > 0 ) { h <- rep[1] } else { next }

# start replicate
print(paste0(". . . starting iteration ", scen, " spec ", spec, " rep ", h))

# update track
load(track_path)
track <- rbind(track, data.table(iter = niter,
  replicate = h, running = 1, finished = 0))
setkey(track, iter, replicate)
save(track, file = track_path)

# extract data
tn <- n[iteration == scen & replication == h] # network
tb <- b[iteration == scen & replication == h] # behavior
# td <- d[iteration == scen & replication == h] # position

#+ loop to create networks and get behavior
gnet <- list()
net <- list()
attr <- list()

# loop to get data in the right format
for (i in 1:length(nwaves)) {
    t <- as.data.frame(tn[measurement == nwaves[i], 1:2])
    g <- graph.data.frame(t, vertices = 1:agents, directed = TRUE)
    m <- get.adjacency(g, sparse = FALSE)
    gnet[[i]] <- g
    net[[i]] <- m
    attr[[i]] <- tb[measurement == nwaves[i], list(id, mybehavior, radius, tendency)]
}

remove(t, g, m)

#  dependent behavior variable
at <- rbindlist(attr, idcol = "time")
at[, mybehavior := as.numeric(mybehavior)]
at <- dcast(at, id + radius + tendency ~ time, value.var = c("mybehavior"), )
beh <- as.matrix(at[, as.character(seq_along(nwaves)), with = FALSE])

# use the corresponding specification and save model
#+ set data for network data
subnet <- net[seq_along(nwaves)]
net_array <- array(unlist(subnet), dim = c(nrow(net[[1]]), ncol(net[[1]]),
                                           length(subnet)))
network <- sienaDependent(net_array)

# behavior
beh <- sienaDependent(beh, type = "behavior")

##################
# covariates
#################

# radius covariate
radius <- coCovar(log(at$radius))

# tendency
tendency <- coCovar(ifelse(at$tendency == TRUE, 1, 0)) # FIXME: change other syntaxes

##########################################
# create data object by specification
##########################################

myData <- sienaDataCreate(network, beh, tendency)

# create effects object for model specification
myEffects <- getEffects(myData)

# varying density by wave
myEffects <- includeTimeDummy(myEffects, density, timeDummy = "all")

# structural network effects
myEffects <- includeEffects(myEffects, recip, outAct, outPop, inPop,
                       fix = FALSE, test = FALSE, include = TRUE)
myEffects <- includeEffects(myEffects, cycle3,
                            fix = FALSE, test = FALSE, include = TRUE)
myEffects <- includeEffects(myEffects, transTrip, transRecTrip,
                            fix = FALSE, test = FALSE, include = TRUE)
myEffects <- setEffect(myEffects, gwespFF, parameter = 69)

# selection (similarity)
myEffects <- includeEffects(myEffects, egoX, altX, simX,
                            interaction1 = "beh", fix = FALSE,
                            test = FALSE, include = TRUE)

# influence (assimilation)
myEffects <- includeEffects(myEffects, name = "beh", avAlt,
                            interaction1 = "network",
                            fix = FALSE, test = FALSE, include = TRUE)

myEffects <- includeEffects(myEffects, effFrom,
                            interaction1 = 'tendency', name = "beh",
                            fix = FALSE, test = FALSE, include = TRUE)

#+ define algorithm
# n3 suggested to be 3000 for publication
myalgorithm <- sienaAlgorithmCreate(nsub = 4, n3 = sim)

# save results
model <- list()

model[[as.character(h)]] <- tryCatch(
  siena07ToConvergence(myalgorithm, dat = myData,
    eff = myEffects, nodes = nodes), error = function (e) e)

rr <- paste0(results_path, "rr_iter_", scen, "_spec_", spec, "_rep_", h, ".rds")
tryCatch(saveRDS(model, file = rr), error = function (e) print("Error saving models!"))

print(screenreg(model[[as.character(h)]]))
print(paste0(". . . ending iteration ", scen, " spec ", spec, " rep ", h))

# update track
load(track_path)
track[iter == niter & replicate == h, c('running', 'finished') := list(0, 1)]
track <- unique(track)
save(track, file = track_path)

} # end loop through replications

} # end loop through iteration/specs

###########################################################
