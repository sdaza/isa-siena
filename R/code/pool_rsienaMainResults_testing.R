###################################################################
# pool rsiena main results, compute coverage rate and standardized bias
# author: sebastian daza
# version: 1.00
###################################################################

#+ clean workspace and load  libraries
# rm(list=ls(all=TRUE))

# libraries
library(sdazar)
library(RSienaTest)
library(texreg)
library(igraph)
library(sna)
# library(colorout)
library(ggplot2)
library(metafor)
library(xtable)
library(rlist) # to manage model lists

#+ my laptop paths
setwd("/Users/sdaza/Desktop/model/")

#+ linstat paths
# setwd("/home/s/sdaza/00projects/siena/")
# source("/home/s/sdaza/00projects/siena/R/functions/auxiliary-functions.R")
# source("/home/s/sdaza/00projects/siena/R/functions/sienaGroupGOF.R")

#+ experiment information (important)
scenarios <- c("B", "I", "S","B+CR", "I+CR", "S+CR") # iterations
specifications <- c("N+S", "N+I", "N+S+I", "N+S+I+D", "N+S+I+R",
                    "S+I+R+D", "N+S+I+R+D")

# sienaeffects <- c("SS", "A", "SA")

iterations <- 1:length(scenarios)
specs <- 1:length(specifications)
# seffects <- 1:length(sienaeffects)

#+ define scenarios and specifications for loop
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

# grid <- grid[c(1,2,3,4,5)] # for testing
# grid

# set max number of replicates
max_rep <- 30

#+ list to save latex rows
latex_rows <- list()

#+ start loop by scenario and Rsiena specification

j <- 1

scen <- grid[j, iter]
spec <- grid[j, sp]

# load(paste0("R/results/rr_iter_", scen, "_spec_", spec, ".Rdata"))
load(paste0("R/results/main/rr_iter_", scen, "_spec_", spec, ".Rdata"))

# check convergence
rep <- length(models)
conv <- NULL
for (i in 1:rep) {
  conv <- c(conv, ifelse(is.null(models[[i]]$tconv.max), NA, models[[i]]$tconv.max))
}

hist(conv)

# remove models didn't converge
models[which(conv > .25 | is.na(conv))] <- NULL

# get one replication and explore effects
length(models)
print(screenreg(models[[1]]))


# extract coefficients and standard errors
estimates <- matrix(nr = rep, nc = 2) # two parameters (selection and influence)
sterrors <- matrix(nr = rep, nc = 2)

# order is important
colnames(estimates) <- c("selection", "influence")
rownames(estimates) <- 1:rep

colnames(sterrors) <- c("selection", "influence")
rownames(sterrors) <- 1:rep

# estimates
# sterrors

for (i in 1:rep) {
  myEffects <- which(models[[i]]$effects$effectName %in% c("beh similarity",
  # myEffects <- which(models[[i]]$effects$effectName %in% c("beh ego x beh alter",
                                                           # "beh average similarity",
                                                           "beh average alter"))
  estimates[i, ] <- models[[i]]$theta[myEffects]
  sterrors[i, ] <- sqrt(diag(models[[i]]$covtheta))[myEffects]

}

estimates
sterrors

# known parameter
parameter <- 0

# getting z-values for 95% of confidence (or nominal coverage rate)
z <- 1.96

# get vectors of estimates and standard errors (to reduce cluttering)
selection  <- estimates[, "selection"]
selection_se  <- sterrors[, "selection"]
influence  <- estimates[, "influence"]
influence_se  <- sterrors[, "influence"]

selection
influence
#+ create different estimates depending on scenario

# no selection influence or selection
covrate_s <- mean(
  ifelse(
  ((selection - z * selection_se) < parameter) &
  ((selection + z * selection_se) > parameter),
  1, 0)
)

test <- data.table(cbind(selection, selection_se))
test[, lower := selection - z * selection_se]
test[, upper := selection + z * selection_se]
test[, check := ifelse(lower < parameter & upper > parameter, 1, 0)]

metaresults <- list()
for (i in 1:length(myEffects)) {
  metaresults[[i]] <- rma(yi = estimates[,i], sei = sterrors[,i], method = "REML")
}

covrate_s
(covrate_i <- mean(ifelse(influence - z * influence_se < parameter
                        & influence + z * influence_se > parameter, 1, 0)))

(sbias_s <- ( mean(selection) - parameter ) / sd( selection))
(sbias_i <- ( mean(influence) - parameter ) / sd( influence))

meta_s <- metaresults[[1]]$b[1]
meta_s_se <- metaresults[[1]]$se[1]

meta_i <- metaresults[[2]]$b[1]
meta_i_se <- metaresults[[2]]$se[1]



#####################################
# end
#####################################
