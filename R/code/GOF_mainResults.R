##############################################################
# creating GOF plots specific iterations from main results
# author: sebastian daza
# version: 1.00
############################################################

#+ clean workspace and load  libraries
rm(list=ls(all=TRUE))

# libraries
library(sdazar)
library(RSienaTest)
library(igraph)
library(sna)
library(colorout)
library(rlist) # to manage lists

#+ linux server paths
setwd("/home/s/sdaza/00projects/siena/")
source("/home/s/sdaza/00projects/siena/R/functions/auxiliary-functions.R")
source("/home/s/sdaza/00projects/siena/R/functions/sienaGroupGOF.R")
source("/home/s/sdaza/00projects/siena/R/functions/extract_texreg.R")

#+ experiment information (important)
scenarios <- c("B", "I", "S","B+CR", "I+CR", "S+CR") # iterations
specifications <- c("N+S", "N+I", "N+S+I", "N+S+I+D", "N+S+I+R",
                    "S+I+S+R+D", "N+S+I+R+D")

# sienaeffects <- c("SS", "A", "SA")

iterations <- 1:length(scenarios)
specs <- 1:length(specifications)
# seffects <- 1:length(sienaeffects)

# max replication for group GOF
max_rep <- 100

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

# grid
grid <- grid[c(4,8)] # only two cases

#+ start loop by scenario and siena specification
for (j in 1:nrow(grid)) {

scen <- grid[j, iter]
spec <- grid[j, sp]

load(paste0("R/results/main/rr_iter_", scen, "_spec_", spec, ".Rdata"))

# check convergence
rep <- length(models)
conv <- NULL
for (i in 1:rep) {
conv <- c(conv, ifelse(is.null(models[[i]]$tconv.max), NA, models[[i]]$tconv.max))
}

# remove models didn't converge
models[which(conv > .25)] <- NULL

#+ create list to save gofs
lid <- list()
lod <- list()
lgd <- list()
ltc <- list()
lbh <- list()
lau <- list()

#+ loop

for (i in 1:max_rep) {

      print(paste0("::::::: iteration ", j, " replicate ", i, " :::::::::::::"))

      lbh[[i]] <- sienaGOF(models[[i]], BehaviorDistribution,
                  verbose = TRUE, join = TRUE, varName = "beh")

      lid[[i]] <- sienaGOF(models[[i]], verbose = TRUE, varName= "network",
                           IndegreeDistribution, join= TRUE, cumulative = TRUE)

      lod[[i]] <- sienaGOF(models[[i]], verbose=TRUE, varName = "network",
                           OutdegreeDistribution, join = TRUE, cumulative = TRUE)

      lgd[[i]] <- sienaGOF(models[[i]], verbose=TRUE, varName = "network",
                           GeodesicDistribution, join = TRUE, cumulative = TRUE)


      ltc[[i]] <- sienaGOF(models[[i]], TriadCensus, verbose = TRUE, join = TRUE,
                          varName = "network")

      lau[[i]] <- sienaGOF(models[[i]], MoranGeary, verbose =TRUE, join = FALSE, varName = c("network","beh"))

}

# plots
savepdf(paste0("R/results/plots/bh_", j), 13, 16)
print(plot.sienaGroupGOF(lbh))
dev.off()

savepdf(paste0("R/results/plots/od_", j), 13, 16)
print(plot.sienaGroupGOF(lod))
dev.off()

savepdf(paste0("R/results/plots/id_", j), 13, 16)
print(plot.sienaGroupGOF(lid))
dev.off()

savepdf(paste0("R/results/plots/gd_", j), 13, 16)
print(plot.sienaGroupGOF(lgd))
dev.off()

savepdf(paste0("R/results/plots/tc_", j), 13, 16)
print(plot.sienaGroupGOF(ltc))
dev.off()

savepdf(paste0("R/results/plots/au_", j), 13, 16)
print(plot.sienaGroupGOF(lau))
dev.off()

}

# end script

