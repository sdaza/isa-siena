########################################################
# assess and define the number of replicates
# author: sebastian daza
# version: 1.00
#######################################################

#+ libraries
library(data.table)
library(igraph)
library(sna)
library(sdazar)
library(ggplot2)
library(ggthemes)

gpath <- "/Users/sdaza/Desktop/model/"
folder_output <- "replicates2000"

source(paste0(gpath,"R/functions/auxiliary-functions.R"))

# compute network statistics

#+ set values
agents <- 200 # define number of agents to create network
iteration <- 1:2 # only two, influence and selection
lab_iteration <- c("I: Only Influence", "S: Only Selection")
nwaves <- c(6:9) # waves to consider
max_rep <- 2000

#+ define paths
fpath <- paste0(gpath, "output/", folder_output) # files

# get files
behavior.files <- list.files(path = fpath, pattern = "behavior")
b <- lapply(paste0(fpath, "/", behavior.files), fread, sep=",")
b <- rbindlist(b)

network.files <- list.files(path = fpath, pattern = "network")
n <- lapply(paste0(fpath, "/", network.files), fread, sep=",")
n <- rbindlist(n)

# categorize behavior variable
b$mybehavior <- cut(b$behavior, breaks = 10, labels = 1:10)
table(b$mybehavior, useNA = "ifany")

#+ loop through iterations

for (j in iteration) {

  #+ create variables
  lmoranb <- list()
  rep <- 1:max_rep

  #+ loop through replications
  for (h in rep) {

  print(paste0(". . . iteration ", j, " replication ", h))

  # extract data
  tn <- n[iteration == j & replication == h]
  tb <- b[iteration == j & replication == h]
  # td <- d[iteration == j & replication == h]

  #+ loop to create networks and get behavior
  gnet <- list()
  net <- list()
  attr <- list()

  # loop to get objects
  for (i in 1:length(nwaves)) {
      t <- as.data.frame(tn[measurement == nwaves[i], 1:2])
      g <- graph.data.frame(t, vertices = 1:agents, directed = TRUE)
      m <- get.adjacency(g, sparse = FALSE)
      gnet[[i]] <- g
      net[[i]] <- m
      attr[[i]] <- tb[measurement == nwaves[i], list(id, mybehavior, radius, tendency)]
  }

  moranb <- NA

  for (i in 1:length(nwaves)) {
    bh <- as.numeric(attr[[i]][, mybehavior])
    rd <- as.numeric(attr[[i]][, radius])
    moranb[i] <- Moran(bh, net[[i]])
  }

  lmoranb[[h]] <- mean(moranb)

  } # end loop through replications

  # coefficient of variation per replication
  mb <- unlist(lmoranb)
  mb
  CV <- function(x) (sd(x)/mean(x))

  print(paste0("Length Moran's I: ", length(mb)))


  # get samples and plot
  cvs <- list()
  for (i in 2:500) { # size of samples
    s <- sample(rep, i, replace = FALSE) # with 2000 replicates, max sample of 1/4
    cvs[[i - 1]] <-  CV(mb[s])
  }

  cvs <- data.frame(cvs = unlist(cvs), n = 2:500)
  assign(paste0("test", j), cvs)

  savepdf(paste0(gpath, "R/results/plots/cv_",j))
    print(ggplot(cvs, aes(y = cvs, x = n)) +  ylim(0, 0.5) +
    geom_point()  + theme_light() +
    labs(title = lab_iteration[j], x = "\nSample sizes", y =  "CV Moran's I\n") +
    geom_vline(xintercept = 100, linetype = 3))
  dev.off()

} # end loop through iterations
