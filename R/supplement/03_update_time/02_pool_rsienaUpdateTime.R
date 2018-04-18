########################################################
# pool rsiena results varying measurement update
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
# library(colorout)
library(ggplot2)
library(metafor)
library(xtable)
library(rlist) # to manage model lists

#+ linstat paths
setwd("/home/s/sdaza/00projects/siena/")
source("/home/s/sdaza/00projects/siena/R/functions/sienaGroupGOF.R")
results_path <- "R/results/timing/"

#+ my laptop paths
# setwd("/Users/sdaza/Desktop/model/")

#+ experiment information (important)
scenarios <- c("I", "S") # iterations
specifications <- c("5", "10", "20")

iterations <- 1:length(scenarios)
update <- 1:length(specifications)

max_rep <- 100

# create R object
loop <- data.table(expand.grid(iter = iterations, up = update))
loop[, sp := 3] # rsiena specification 3

#+ siena setup
sim <- 1000 # number of simulations

###########################################################
#+ loop through iterations
###########################################################

#+ list to save latex rows
latex_rows <- list()

#+ start loop by scenario and siena specification

for (j in 1:nrow(loop)) {

scen <- loop[j, iter]
up <- loop[j, up]
spec <- loop[j, sp]

# combine files
mpattern <- paste0("rr_iter_", scen, "_spec_", up)
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

trep <- 1:max_rep

print(paste0("Missing replications (", max_rep - length(rep), ") iter ", scen,  " spec ", up, ": ",
  paste0(trep[!trep %in% rep], collapse = ", ")))

# extract coefficients and standard errors
estimates <- matrix(nr = length(rep), nc = 2) # two parameters (selection and influence)
sterrors <- matrix(nr = length(rep), nc = 2)

# order is important
colnames(estimates) <- c("selection", "influence")
rownames(estimates) <- rep

colnames(sterrors) <- c("selection", "influence")
rownames(sterrors) <- rep

# estimates
# sterrors

for (i in as.character(rep)) {
  myEffects <- which(models[[i]]$effects$effectName %in%
    c("beh similarity", "beh average alter"))
  estimates[i, ] <- models[[i]]$theta[myEffects]
  sterrors[i, ] <- sqrt(diag(models[[i]]$covtheta))[myEffects]

}

# meta analysis (using package metafor)
metaresults <- list()
for (i in 1:length(myEffects)) {
  metaresults[[i]] <- rma(yi = estimates[,i], sei = sterrors[,i], method = "REML")
}


#+ assessment measures

# known parameter
parameter <- 0

# getting z-values for 95% of confidence (or nominal coverage rate)
z <- 1.96

# get vectors of estimates and standard errors (to reduce cluttering)
selection  <- estimates[, "selection"]
selection_se  <- sterrors[, "selection"]
influence  <- estimates[, "influence"]
influence_se  <- sterrors[, "influence"]

#+ create different estimates depending on scenario

# only influence
if (scen == 1) {
(covrate_s <- mean(ifelse(selection - z * selection_se < parameter
                        & selection + z * selection_se > parameter, 1, 0)))

(sbias_s <- ( mean(selection) - parameter ) / sd( selection))

(pos_i <- mean(ifelse(influence - z * influence_se > 0 &
                    influence + z * influence_se > 0, 1, 0)))

meta_s <- metaresults[[1]]$b[1]
meta_s_se <- metaresults[[1]]$se[1]

meta_i <- metaresults[[2]]$b[1]
meta_i_se <- metaresults[[2]]$se[1]
}


# only selection
if (scen == 2) {

(covrate_i <- mean(ifelse(influence - z * influence_se < parameter
                        & influence + z * influence_se > parameter, 1, 0)))

(sbias_i <- ( mean(influence) - parameter ) / sd( influence))

(pos_s <- mean(ifelse(selection - z * selection_se > 0 &
                    selection + z * selection_se > 0, 1, 0)))

meta_s <- metaresults[[1]]$b[1]
meta_s_se <- metaresults[[1]]$se[1]

meta_i <- metaresults[[2]]$b[1]
meta_i_se <- metaresults[[2]]$se[1]
}

# rows for latex table

# only influence
if (scen == 1) {
latex_rows[[j]] <- paste0(
  scenarios[scen], " & N+S+I & & ",
  format(round(covrate_s, 2), nsmall = 2), " & ",
  format(round(sbias_s, 2), nsmall = 2), " & - & ",
  format(round(meta_s, 2), nsmall = 2) , " (",
  format(round(meta_s_se, 2), nsmall = 2), ") & & - & - & ",
  format(round(pos_i, 2), nsmall = 2), " & ",
  format(round(meta_i, 2), nsmall = 2) , " (",
  format(round(meta_i_se, 2), nsmall = 2), ") \\\\ \n"
  )
}

# only selection
if (scen == 2 ) {
latex_rows[[j]] <- paste0(
   scenarios[scen], " & N+S+I & & - & - &",
   format(round(pos_s, 2), nsmall = 2), " & ",
   format(round(meta_s, 2), nsmall = 2) , " (",
   format(round(meta_s_se, 2), nsmall = 2), ") & & ",
   format(round(covrate_i, 2), nsmall = 2), " & ",
   format(round(sbias_i, 2), nsmall = 2), " & - & ",
   format(round(meta_i, 2), nsmall = 2) , " (", format(round(meta_i_se, 2), nsmall = 2), ") \\\\ \n"
  )
}

} # end of the experiment analysis


############################
#+ create latex table
#############################

top <- paste0("
\\renewcommand{\\arraystretch}{1.2}
% results table
\\begin{table}[htp]
\\begin{threeparttable}
\\setlength{\\tabcolsep}{5pt}
\\centering
\\scriptsize
\\caption{Coverage Rate, Bias and Siena Estimates by \\newline Scenario, Specification and Update Time (", max_rep, " replicates)} \\label{suppl:update_time}
\\begin{tabular}{lllcccclcccc}
\\hline
\\addlinespace
\\multicolumn{3}{c}{}  & \\multicolumn{9}{c}{Estimates (Siena)}  \\\\
\\addlinespace
\\cmidrule{4-12}
\\addlinespace
 & & & \\multicolumn{4}{c}{Selection} &  & \\multicolumn{4}{c}{Influence} \\\\
\\addlinespace
\\cmidrule{1-2} \\cmidrule{4-7} \\cmidrule{8-12}
ISA Scenario & Siena Specification  &  & Coverage  & Bias & $Pr(\\beta_s>0)$ & RE Estimate &  & Coverage  & Bias & $Pr(\\beta_i>0)$ & RE Estimate \\\\
\\addlinespace
\\hline\n")

bottom <- "\\addlinespace
\\addlinespace
\\hline
\\end{tabular}
   \\begin{tablenotes}[flushleft]
      \\scriptsize
      \\item ISA scenario: B = Baseline, I = Influence, S = Selection, CR = Constant radius.
      \\item Siena specification: N = Network structural effects, S = Selection, I = Influence,
      D = Distance, R = Radius.
      % \\item Performance statistics: $\\beta_s$ = Behavior Similarity, $\\beta_i$ = Behavior Average Similarity.
       \\item Performance statistics: $\\beta_s$ = Behavior Similarity, $\\beta_i$ = Behavior Average Alter.
      % \\item Performance statistics: $\\beta_s$ = Behavior Ego-Alter Interaction, $\\beta_i$ = Behavior Average Alter.
      \\item $Pr(\\beta>0)$ = Probability confidence interval includes only positive values.
      \\item RE estimate = average estimate of the true effect using a random-effects model.
    \\end{tablenotes}
  \\end{threeparttable}
\\end{table}"


sep <- list("\\addlinespace
             \\addlinespace
             \\multicolumn{12}{l}{\\textbf{Update Time = 5 days}}\\\\
             \\addlinespace
             \\addlinespace\n",
              "\\addlinespace
             \\addlinespace
              \\multicolumn{12}{l}{\\textbf{Update Time = 10 days}}\\\\
              \\addlinespace
              \\addlinespace\n",
              "\\addlinespace
              \\addlinespace
              \\multicolumn{12}{l}{\\textbf{Update Time = 20 days}}\\\\
              \\addlinespace
              \\addlinespace\n",
              " \\addlinespace
                \\addlinespace
                \\hdashline\n")

clean_latex_rows <- list.clean(c(
                           sep[1], latex_rows[1:2], sep[4],
                           sep[2], latex_rows[3:4], sep[4],
                           sep[3], latex_rows[5:6]
                           ))

# create and save table
cat(top, do.call(paste, clean_latex_rows) , bottom,
    file = paste0("latex/tables/results_timing", ".tex"))

# end script
