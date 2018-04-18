########################################################
# ISA descriptives
# author: sebastian daza
# 2017-05-05
#######################################################

#+ libraries
library(data.table)
library(RSienaTest)
library(texreg)
library(igraph)
library(sna)
library(sdazar)
# library(colorout)

# library(intergraph)
# library(GGally)

# linstat
# setwd("/home/s/sdaza/00projects/saom/")
# source("/home/s/sdaza/00projects/saom/R/auxiliary-functions.R")

setwd("/Users/sdaza/Desktop/model")
source("R/functions/auxiliary-functions.R")

#+ compute statistics

# experiments from anylogic
# uniform distribution

# selection = new boolean[] {false, false, true, false, false, true};
# influence = new boolean[] {false, true, false, false, true, false};
# behavior = new int[] {0,0,0,0,0,0};
# radius = new boolean[] {false, false, false, true, true, true};

# define iterations and some ISA parameters

#+ set values
agents <- 200 # define number of agents to create network
iteration <- 1:6 # ISA iterations
nwaves <- c(6:9) # waves to include
max_rep <- 100 # number of replicates

#+ load data
load("R/data/mainResults.Rdata")

# checks on behavior variable
print(setkey(b[, .(min = min(behavior), max = max(behavior)),
             by = mybehavior], mybehavior)) # right

# hist(b[mybehavior == 1, behavior]) # distribution of 1 and 10 is skew
# hist(b[mybehavior == 10, behavior])
# hist(b[mybehavior == 2, behavior])
# hist(b[mybehavior == 8, behavior])
# hist(b[mybehavior == 4, behavior])
# hist(b[mybehavior == 6, behavior])

table(b$iteration) # 6
# table(b$measurement)

# matrix to save results
cnames <- c("Density", "Degree", "Reciprocity", "Transitivity",
            "Moran's Behavior", "Moran's Radius", "Jaccard", "Stability Behavior")


est <- matrix(NA, ncol = length(cnames), nrow = max(iteration))
sds <- matrix(NA, ncol = length(cnames), nrow = max(iteration))

#+ loop through iterations and replicates

for (j in iteration) {

  #+ create variables
  ldens <- list()
  lavdegree <- list()
  lrecip <- list()
  ltrans <- list()
  lmoranb <- list()
  lmoranr <- list()
  ljaccard <- list()
  lstbeh <- list()

  rep <- 1:max_rep

  #+ loop through replications
  for (h in rep) {

  print(paste0(". . . iteration ", j, " replication ", h))

  # extract data (network and behavior)
  tn <- n[iteration == j & replication == h]
  tb <- b[iteration == j & replication == h]

  #+ loop to create networks and get behavior
  gnet <- list()
  net <- list()
  attr <- list()

  # loop to get network and attribute data
  for (i in 1:length(nwaves)) {
      t <- as.data.frame(tn[measurement == nwaves[i], 1:2])
      g <- graph.data.frame(t, vertices = 1:agents, directed = TRUE)
      m <- get.adjacency(g, sparse = FALSE)
      gnet[[i]] <- g
      net[[i]] <- m
      attr[[i]] <- tb[measurement == nwaves[i], list(id, mybehavior, radius, tendency)]
  }

  dens <- NA; avdegree <- NA; trans <- NA
  recip <- NA; moranb <- NA; moranr <- NA

  for (i in 1:length(nwaves)) {
    bh <- as.numeric(attr[[i]][, mybehavior])
    rd <- as.numeric(attr[[i]][, radius])
    dens[i] <- igraph::edge_density(gnet[[i]])
    avdegree[i] <- mean(igraph::degree(gnet[[i]]))
    trans[i] <- igraph::transitivity(gnet[[i]])
    recip[i] <- igraph::reciprocity(gnet[[i]])
    moranb[i] <- Moran(bh, net[[i]])
    moranr[i] <- Moran(rd, net[[i]])
  }

  # stability behavior
  stability <- NA; bdistance <- NA; jaccard <- NA; ndistance <- NA

  for (i in 1:(length(nwaves)-1)) {
    temp1 <- copy(attr[[i]])
    temp2 <- copy(attr[[i + 1]])
    change <- table(temp1$mybehavior, temp2$mybehavior)
    stability[i] <- sum(diag(change)) / sum(change)
    bdistance[i] <- sum(change[row(change) != col(change)])
    n1 <- net[[i]]
    n2 <- net[[i + 1]]
    t <- table(n1, n2)
    jaccard[i] <-  t[2, 2] / (t[2, 1] + t[1, 2] + t[2, 2])
    ndistance[i] <- t[2, 1] + t[1, 2]
  }

  ldens[[h]] <- mean(dens)
  lavdegree[[h]] <- mean(avdegree)
  lrecip[[h]] <- mean(recip)
  ltrans[[h]] <- mean(trans)
  lmoranb[[h]] <- mean(moranb)
  lmoranr[[h]] <- mean(moranr)
  ljaccard[[h]] <- mean(jaccard)
  lstbeh[[h]] <- mean(stability)


} # end loop through replications

listnames <- c("ldens", "lavdegree", "lrecip", "ltrans", "lmoranb", "lmoranr",
         "ljaccard", "lstbeh")

# mean through replications
for (i in 1:length(listnames)) {
  est[j,i] <- round(mean(unlist(get(listnames[i]))), digits = 2)
  sds[j,i] <- round(sd(unlist(get(listnames[i]))), digits = 2)
}

} # end loop through iterations

#+ create table
fest <- as.vector(format(est, nsmall = 2))
fsds <- paste0("(",as.vector(format(sds, nsmall = 2)),")")
M <- matrix(rbind(fest, fsds), nrow = length(iteration) * 2)

# define scenario labels
scenarios <- c("B", " ", "I", " ", "S", " ",
               "B+CR", " ",  "I+CR", " ", "S+CR", " ")

# list to save table rows
latexrows <- list()

for (i in 1:(max(iteration)*2)) {
  latexrows[[i]] <- paste0(scenarios[i], " & ",  paste0(M[i,], collapse = " & "), " \\\\ \n ")
}

latexrows <- lapply(latexrows, function(y) gsub("\\(  NA\\)", "--", y))
latexrows <- lapply(latexrows, function(y) gsub("  NaN", "--", y))
sep <- list("\\addlinespace[10pt] \n")

latexrows <- c(latexrows[1:6], sep, latexrows[7:12])

top <- paste0("
  % ISA descriptive table
  \\renewcommand{\\arraystretch}{1.2} % distance between rows
  \\begin{table}[htp]
  \\begin{threeparttable}
  \\setlength{\\tabcolsep}{2.5pt}
  \\centering
  \\scriptsize
  \\caption{ISA Scenario Descriptives (",
  max_rep, " replicates)} \\label{desc}
  \\begin{tabular}{lcccccccc}
  \\hline
  \\addlinespace
   & \\multicolumn{8}{c}{Statistics}  \\\\
  \\cmidrule{2-9}
  ISA Scenario  & Density & Degree (avg) & Reciprocity & Transitivity &
  Moran's I Behavior & Moran's I Radius  & Jaccard & Stability Behavior \\\\
  \\addlinespace
  \\hline
  \\addlinespace \n")

bottom <- "\\addlinespace
  \\hline
  \\end{tabular}
     \\begin{tablenotes}[flushleft]
        \\scriptsize
          \\item Standard errors in parentheses.
          \\item B = Baseline; S = Selection; I = Influence; CR = Constant radius.
      \\end{tablenotes}
    \\end{threeparttable}
  \\end{table}"

cat(top, do.call(paste, latexrows), bottom,
    file = "latex/tables/isa_descriptive.tex")

#+ finish table
