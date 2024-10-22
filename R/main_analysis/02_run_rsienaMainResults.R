########################################################
# rsiena model using replicates
# uniform behavior distribution
# author: sebastian daza
# version: 1.00
#######################################################

# run using multiple linux servers

# libraries
library(data.table)
library(RSienaTest)
library(texreg)
library(igraph)
library(sna)
library(sdazar)
library(colorout)
library(rlist)
library(SGCS) # for compute toroidal distances

# adjust paths
setwd("/home/s/sdaza/00projects/siena/")
source("/home/s/sdaza/00projects/siena/R/functions/auxiliary-functions.R")
results_path <- "R/results/main/"

#+ define scenarios and specifications for loop

# Reference from Anylogic experiment setup (Java)
# selection = new boolean[] {false, false, true, false, false, true};
# influence = new boolean[] {false, true, false, false, true, false};
# behavior = new int[] {0,0,0,0,0,0};
# radius = new boolean[] {false, false, false, true, true, true};

# experiment information (important)
scenarios <- c("B", "I", "S","B+CR", "I+CR", "S+CR") # iterations
specifications <- c("N+S", "N+I", "N+S+I", "N+S+I+D", "N+S+I+R",
                    "S+I+R+D", "N+S+I+R+D")

iterations <- 1:length(scenarios)
specs <- 1:length(specifications)

# create R object
loop <- data.table(expand.grid(iter = iterations, sp = specs))

g <- list()
g[[1]] <- loop[iter %in% c(1,4) & sp == 3]
g[[2]]<- loop[iter == 2 & sp %in% c(1,3,4,5,6,7)]
g[[3]] <- loop[iter == 3 & sp %in% c(2,3,4,5,6,7)]
g[[4]] <- loop[iter == 5 & sp %in% c(3,4)]
g[[5]] <- loop[iter == 6 & sp %in% c(3,4)]

grid <- rbindlist(g) # grid to create loop
grid
rm(loop, g)

#+ set up values ISA
agents <- 200 # define number of agents to create network

# define number of cores
# (nodes <- floor(parallel::detectCores() / 3))
nodes <- 5

wi <- 500; le <- 500 # dimensions of space
max_rep <- 100 # max number of replicates FIXME: adjust later
nwaves <- c(6,7,8,9) # four waves

#+ siena setup
sim <- 1000 # number of simulations

# read data generated by Anyloic (b, n, d, p)
load("R/data/mainResults.Rdata")

p # explore parameters parameters

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
#+ loop trough scenario and specification and replication
##########################################################

for (j in 1:nrow(grid)) {
# for (j in c(1,3,4,5,8)) { # testing

  # get scenario (iteration) and siena specification
  scen <- grid[j, iter]
  spec <- grid[j, sp]

 # adjust accordingly
  print(paste0(":::::::::::::::: ", paste0(scen, "-",  spec)))

  # name iteration
  niter <- paste0(scen, "-", spec)

  rep <- keeptrack(track_path)
  if (length(rep) == 0) { next }

  ###############################
  # loop through replications
  ###############################

  for (h in rep) {

  # update replicates
  rep <- keeptrack(track_path)

  # rewrite h
  if (length(rep) > 0 ) { h <- rep[1] } else { next }

  # start replicate
  print(paste0(". . . starting iteration ", scen, " spec ", spec, " rep ", h))

  # update track
  load(track_path)
  track <- rbind(track,
    data.table(iter = niter, replicate = h, running = 1, finished = 0))
  setkey(track, iter, replicate)
  save(track, file = track_path )

  # select  data
  tn <- n[iteration == scen & replication == h] # network
  tb <- b[iteration == scen & replication == h] # behavior
  td <- d[iteration == scen & replication == h] # position

  #+ loop to create networks and get behavior
  gnet <- list()
  net <- list()
  attr <- list()
  dist <- list()

  # loop to get data in the right format
  for (i in 1:length(nwaves)) {
      t <- as.data.frame(tn[measurement == nwaves[i], 1:2])
      g <- graph.data.frame(t, vertices = 1:agents, directed = TRUE)
      m <- get.adjacency(g, sparse = FALSE)
      gnet[[i]] <- g
      net[[i]] <- m
      attr[[i]] <- tb[measurement == nwaves[i], list(id, mybehavior, radius, tendency)]
      xy <- list(x  = td[measurement == nwaves[i], x], y = td[measurement == nwaves[i], y])
      pp <- as.ppp(xy, c(0, wi, 0, le)) # convert to class ppp used by the spatstat package
      dpp <- log(pairwise_distances(pp, toroidal = TRUE)) # log of toroidal distances
      dm <-  matrix(0, agents, agents)
      dm[lower.tri(dm, diag = FALSE)] <- dpp
      dm[upper.tri(dm)] <- t(dm)[upper.tri(dm)] # create the full matrix
      dist[[i]] <- as.matrix(dm)
      isSymmetric(dist[[i]]) # check if matrix is symmetric
  }

  remove(t, g, m, xy, pp) # remove temp objects

  #  dependent behavior variable
  at <- rbindlist(attr, idcol = "time")
  prop.table(table(at[time == 1, .(tendency)], useNA = 'ifany'))
  prop.table(table(at[time == 3, .(tendency)], useNA = 'ifany'))
  at[, mybehavior := as.numeric(mybehavior)]
  at <- dcast(at, id + radius + tendency ~ time, value.var = c("mybehavior"), )

  beh <- as.matrix(at[, as.character(seq_along(nwaves)), with = FALSE])

  # use the corresponding specification and save results

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

  # distance
  subdist <- dist[1:(length(nwaves)-1)]
  if (length(nwaves) > 2) {
        dist_array <- array(unlist(subdist),
                        dim = c(nrow(dist[[1]]), ncol(dist[[1]]), length(subdist)))
        distance <- varDyadCovar(dist_array)
        }
        if (length (nwaves) == 2) {
        distance <- coDyadCovar(subdist[[1]])
        }

  ##########################################
  # create data object by specification
  ##########################################

  if (spec %in% 1:3) {
        myData <- sienaDataCreate(network, beh, tendency)
  }
  if (spec == 4) {
        myData <- sienaDataCreate(network, beh, distance, tendency)
        }
  if (spec == 5) {
        myData <- sienaDataCreate(network, beh, radius, tendency)
  }
  if (spec %in% 6:7) {
        myData <- sienaDataCreate(network, beh, radius, distance, tendency)
  }

  # create effects object for model specification
  myEffects <- getEffects(myData)

  # density coefficient varies by wave
  myEffects <- includeTimeDummy(myEffects, density, timeDummy = "all")

  # structural network effects
  if (spec != 6) {
  myEffects <- includeEffects(myEffects, recip, outAct, outPop, inPop,
                         fix = FALSE, test = FALSE, include = TRUE)
  myEffects <- includeEffects(myEffects, cycle3,
                              fix = FALSE, test = FALSE, include = TRUE)
  myEffects <- includeEffects(myEffects, transTrip, transRecTrip,
                              fix = FALSE, test = FALSE, include = TRUE)
  myEffects <- setEffect(myEffects, gwespFF, parameter = 120) # versus 69
  }

  # distance
  if (spec %in% c(4,6,7)) {
  myEffects <- includeEffects(myEffects, X, interaction1 = "distance",
                            fix = FALSE, test = FALSE,  include = TRUE)
  }

  # radius
  if (spec %in% c(5,6,7)) {
  myEffects <- includeEffects(myEffects, egoX, interaction1 = "radius",
                            fix = FALSE, test = FALSE,  include = TRUE)
  }

  # interaction distance - radius
  # if (spec %in% c(6,7)) {
  # myEffects  <- includeInteraction(myEffects, X, egoX,
  #                                interaction1 = c("distance", "radius"),
  #                                fix = FALSE, test = FALSE, include = TRUE)
  # }

  # selection (similarity)
  if (spec %in% c(1,3,4,5,6,7)) {
  myEffects <- includeEffects(myEffects, egoX, altX, simX,
    # similarity
    # myEffects <- includeEffects(myEffects, egoX, altX, altSqX, egoXaltX, # interaction
                              interaction1 = "beh", fix = FALSE,
                              test = FALSE, include = TRUE)
  }

  # influence (assimilation)
  if (spec %in% c(2,3,4,5,6,7)) {
  myEffects <- includeEffects(myEffects, name = "beh", avAlt,
    # average alter
    # myEffects <- includeEffects(myEffects, name = "beh", avSim, # average similiary
                              interaction1 = "network",
                              fix = FALSE, test = FALSE, include = TRUE)

  }

  # for all specifications: adjust for behavior tendency
  myEffects <- includeEffects(myEffects, effFrom,
    interaction1 = 'tendency', name = "beh", fix = FALSE, test = FALSE, include = TRUE) # TODO: add this line to the rest of sintaxis

  #+ define algorithm
  myalgorithm <- sienaAlgorithmCreate(nsub = 4, n3 = sim) # NOTE:  n3=1000 for siena simulations and GOF

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
