########################################################
# creating GOF plots
# author: sebastian daza
# version: 1.00
#######################################################

#+ clean workspace and load  libraries
rm(list=ls(all=TRUE))

# libraries
library(sdazar)
library(RSienaTest)
library(texreg)
library(igraph)
library(sna)
library(colorout)
library(ggplot2)
library(metafor)
library(xtable)
library(rlist) # to manage model lists
library(foreach)
library(doParallel)
registerDoParallel(cores = 5) # define number of cores

#+ linstat paths
setwd("/home/s/sdaza/00projects/siena/")
source("/home/s/sdaza/00projects/siena/R/functions/auxiliary-functions.R")
source("/home/s/sdaza/00projects/siena/R/functions/sienaGroupGOF.R")
source("/home/s/sdaza/00projects/siena/R/functions/extract_texreg.R")
results_path <- "R/results/main/"


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
grid <- grid[c(4,8)] # only select two examples

#+ start loop by scenario and siena specification
for (j in 1:nrow(grid)) {

scen <- grid[j, iter]
spec <- grid[j, sp]

# combine files
mpattern <- paste0("rr_iter_", scen, "_spec_", spec)
files <- list.files(path = results_path, pattern = mpattern)

if (length(files) == 0 ) { next }

files <- paste0(results_path, files)
models <- do.call('list.merge', lapply(files, readRDS))

# clean working space and save files
# file.remove(files)
# rr <- paste0(results_path, "rr_iter_", scen, "_spec_", spec, ".Rdata")
# save(models, file = rr)

# check convergence
rep <- sort(as.numeric(names(models)))

conv <- NULL
for (i in rep) {
conv <- c(conv, ifelse(is.null(models[[as.character(i)]]$tconv.max), NA, models[[as.character(i)]]$tconv.max))
}

# remove models didn't converge
models[as.character(which(conv > .25 | is.na(conv)))] <- NULL
rep <- sort(as.numeric(names(models))) # update replication list

#+ to do it separately by gof to optimize memory use

print("########## behavior")
goflist <- foreach(i = 1:max_rep) %dopar% {
   library(RSienaTest)
   print(paste0("replicate ", i))
   sienaGOF(models[[i]], BehaviorDistribution, levls=1:10,
            verbose = FALSE, join = TRUE, varName = "beh")
}

savepdf(paste0("R/results/plots/bh_", j))
print(plot.sienaGroupGOF(goflist))
dev.off()
remove(goflist)

# indegree
print("########## indegree")
goflist <- foreach(i = 1:max_rep) %dopar% {
   library(RSienaTest)
   print(paste0("replicate ", i))

   sienaGOF(models[[i]], verbose = FALSE, varName= "network",
                     IndegreeDistribution, join= TRUE, cumulative = TRUE)
}

savepdf(paste0("R/results/plots/id_", j))
print(plot.sienaGroupGOF(goflist))
dev.off()
remove(goflist)

# outdegree
print("########## outdegree")
goflist <- foreach(i = 1:max_rep) %dopar% {
   library(RSienaTest)
   print(paste0("replicate ", i))
   sienaGOF(models[[i]], verbose = FALSE, varName = "network",
                        OutdegreeDistribution, join = TRUE, cumulative = TRUE)
}

savepdf(paste0("R/results/plots/od_", j))
print(plot.sienaGroupGOF(goflist))
dev.off()
remove(goflist)

# geodesic distance
print("########## geodesic dist")
goflist <- foreach(i = 1:max_rep) %dopar% {
   library(RSienaTest)
   print(paste0("replicate ", i))
   sienaGOF(models[[i]], verbose = FALSE, varName = "network",
                     GeodesicDistribution, join = TRUE, cumulative = TRUE)
}

savepdf(paste0("R/results/plots/gd_", j))
print(plot.sienaGroupGOF(goflist))
dev.off()
remove(goflist)

# triad census
print("########## triad census")
goflist <- foreach(i = 1:max_rep) %dopar% {
   library(RSienaTest)
   print(paste0("replicate ", i))
   return(sienaGOF(models[[i]], TriadCensus, verbose = FALSE, join = TRUE,
                    varName = "network"))
}

savepdf(paste0("R/results/plots/tc_", j))
print(plot.sienaGroupGOF(goflist, center = TRUE, scale = TRUE))
dev.off()
remove(goflist)

# autocorrelation
print("########## autocorrelation")
goflist <- foreach(i = 1:max_rep) %dopar% {
   library(RSienaTest)
   print(paste0("replicate ", i))
   sienaGOF(models[[i]], MoranGeary, verbose =TRUE, join = TRUE,
                     varName = c("network","beh"))
}

savepdf(paste0("R/results/plots/au_", j))
print(plot.sienaGroupGOF(goflist))
dev.off()
remove(goflist)

}

# end of script
