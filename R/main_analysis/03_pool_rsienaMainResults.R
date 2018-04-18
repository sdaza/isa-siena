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
library(colorout)
library(ggplot2)
library(metafor)
library(xtable)
library(rlist) # to manage model lists

#+ my laptop paths
# setwd("/Users/sdaza/Desktop/model/")

#+ linstat paths
setwd("/home/s/sdaza/00projects/siena/")
source("/home/s/sdaza/00projects/siena/R/functions/auxiliary-functions.R")
source("/home/s/sdaza/00projects/siena/R/functions/sienaGroupGOF.R")
results_path <- "R/results/main/"

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
print(grid)
rm(loop, g)

# grid <- grid[c(1,2,3,4,5)] # for testing

# set max number of replicates
max_rep <- 100

#+ list to save latex rows
latex_rows <- list()

#+ start loop by scenario and Rsiena specification
for (j in 1:nrow(grid)) {

# for (j in 1) {

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

  trep <- 1:max_rep
  print(paste0("Missing replications (", max_rep - length(rep), ") iter ", scen,  " spec ", spec, ": ",
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

  # no selection influence or selection
  if (scen %in% c(1,4)) {
    covrate_s <- mean(
      ifelse(
      ((selection - z * selection_se) < parameter) &
      ((selection + z * selection_se) > parameter),
      1, 0)
    )

  (covrate_i <- mean(ifelse(influence - z * influence_se < parameter
                          & influence + z * influence_se > parameter, 1, 0)))

  (sbias_s <- ( mean(selection) - parameter ) / sd( selection))
  (sbias_i <- ( mean(influence) - parameter ) / sd( influence))

  meta_s <- metaresults[[1]]$b[1]
  meta_s_se <- metaresults[[1]]$se[1]

  meta_i <- metaresults[[2]]$b[1]
  meta_i_se <- metaresults[[2]]$se[1]

  }

  # only influence
  if (scen %in% c(2,5) & spec %in% c(3:7)) {

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

  if (scen %in% c(2,5) & spec == 1) {
  (covrate_s <- mean(ifelse(selection - z * selection_se < parameter
                          & selection + z * selection_se > parameter, 1, 0)))

  (sbias_s <- ( mean(selection) - parameter ) / sd( selection))

  meta_s <- metaresults[[1]]$b[1]
  meta_s_se <- metaresults[[1]]$se[1]

  }

  # only selection
  if (scen %in% c(3,6) & spec %in% c(3:7)) {

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

  if (scen %in% c(3,6) & spec == 2) {

  (covrate_i <- mean(ifelse(influence - z * influence_se < parameter
                          & influence + z * influence_se > parameter, 1, 0)))

  (sbias_i <- ( mean(influence) - parameter ) / sd( influence))

  meta_i <- metaresults[[1]]$b[1]
  meta_i_se <- metaresults[[1]]$se[1]
  }


  # rows for latex table

  # none (influence or selection)
  if (scen %in% c(1,4)) {
  latex_rows[[j]] <- paste0(
    scenarios[scen], " & ", specifications[spec], " & & ",
    format(round(covrate_s, 2), nsmall = 2), " & ",
    format(round(sbias_s, 2), nsmall = 2), " & - & ",
    format(round(meta_s, 2), nsmall = 2) , " (",
    format(round(meta_s_se, 2), nsmall = 2), ") & & ",
    format(round(covrate_i, 2), nsmall = 2), " & ",
    format(round(sbias_i, 2), nsmall = 2), " & - & ",
    format(round(meta_i, 2), nsmall = 2) , " (",
    format(round(meta_i_se, 2), nsmall = 2), ") \\\\ \n"
    )
  }

  # only influence
  if (scen %in% c(2,5) & spec %in% c(3:7)) {
  latex_rows[[j]] <- paste0(
    scenarios[scen], " & ", specifications[spec], " & & ",
    format(round(covrate_s, 2), nsmall = 2), " & ",
    format(round(sbias_s, 2), nsmall = 2), " & - & ",
    format(round(meta_s, 2), nsmall = 2) , " (",
    format(round(meta_s_se, 2), nsmall = 2), ") & & - & - & ",
    format(round(pos_i, 2), nsmall = 2), " & ",
    format(round(meta_i, 2), nsmall = 2) , " (",
    format(round(meta_i_se, 2), nsmall = 2), ") \\\\ \n"
    )
  }

  if (scen %in% c(2,5) & spec == 1) {
  latex_rows[[j]] <- paste0(
    scenarios[scen], " & ", specifications[spec], " & & ",
    format(round(covrate_s, 2), nsmall = 2), " & ",
    format(round(sbias_s, 2), nsmall = 2), " & - & ",
    format(round(meta_s, 2), nsmall = 2) , " (",
    format(round(meta_s_se, 2), nsmall = 2), ") & & - & - & ",
    " - & ",
    " - \\\\ \n"
    )
  }


  # only selection
  if (scen %in% c(3,6) & spec %in% c(3:7)) {
  latex_rows[[j]] <- paste0(
     scenarios[scen], " & ", specifications[spec], " & & - & - &",
     format(round(pos_s, 2), nsmall = 2), " & ",
     format(round(meta_s, 2), nsmall = 2) , " (",
     format(round(meta_s_se, 2), nsmall = 2), ") & & ",
     format(round(covrate_i, 2), nsmall = 2), " & ",
     format(round(sbias_i, 2), nsmall = 2), " & - & ",
     format(round(meta_i, 2), nsmall = 2) , " (", format(round(meta_i_se, 2), nsmall = 2), ") \\\\ \n"
    )
  }

  if (scen %in% c(3,6) & spec == 2) {
  latex_rows[[j]] <- paste0(
     scenarios[scen], " & ", specifications[spec], " & & - & - & ",
     "- & ",
     " -  & & ",
     format(round(covrate_i, 2), nsmall = 2), " & ",
     format(round(sbias_i, 2), nsmall = 2), " & - & ",
     format(round(meta_i, 2), nsmall = 2) , " (", format(round(meta_i_se, 2), nsmall = 2), ") \\\\ \n"
    )
  }

  # for testing!
  # testing1 <- estimates
  # testing2 <- sterrors

} # end of the experiment analysis


# more testing!
# screenreg(models[[5]])
# estimates[5,]
# sterrors[5,]

# maxrep <- length(models)
# savepdf("testing1")
# metafor::forest(metaresults[[1]], slab = 1:maxrep, cex = 0.80, alim = c(-1,1),
#     xlab =  "Selection Coefficients (log scale)", xlim = c(-1.5,2), addcred = TRUE)
#   text(-1.5, maxrep + 2, "Replication", pos = 4, font = 2, cex = 0.80)
#   text(1.3, maxrep + 2, "Estimates [95% CI]", pos = 4, font = 2, cex = 0.80)
# dev.off()

# savepdf("testing2")
# metafor::forest(metaresults[[2]], slab = 1:maxrep, cex = 0.80, alim = c(-1,1),
#     xlab =  "Influence Coefficients (log scale)", xlim = c(-1.5,2), addcred = TRUE)
#   text(-1.5, maxrep + 2, "Replication", pos = 4, font = 2, cex = 0.80)
#   text(1.3, maxrep + 2, "Estimates [95% CI]", pos = 4, font = 2, cex = 0.80)
# dev.off()


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
  \\caption{Coverage Rate, Bias and Siena Estimates by Scenario and Specification (", max_rep, " replicates)} \\label{results_main}
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
  \\hline
  \\addlinespace
  \\addlinespace \n")



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
        % \\item Performance statistics: $\\beta_s$ = Behavior Similarity, $\\beta_i$ = Behavior Average Alter.
        \\item Performance statistics: $\\beta_s$ = Behavior Ego-Alter Interaction, $\\beta_i$ = Behavior Average Alter.
        \\item $Pr(\\beta>0)$ = Probability confidence interval includes only positive values.
        \\item RE estimate = average estimate of the true effect using a random-effects model.
      \\end{tablenotes}
    \\end{threeparttable}
  \\end{table}"

sep <- list("\\addlinespace
            \\addlinespace
             \\hdashline
             \\addlinespace
             \\addlinespace",
             "\\addlinespace
             \\addlinespace
             \\addlinespace
             \\hline
             \\hline
             \\addlinespace
             \\addlinespace
             \\addlinespace")

clean_latex_rows <- list.clean(c(
                           latex_rows[1], sep[1],
                           latex_rows[3:8], sep[1],
                           latex_rows[9:14], sep[2],
                           latex_rows[2], sep[1],
                           latex_rows[15:16], sep[1],
                           latex_rows[17:18]
                           ))

# to explore individual iterations
# clean_latex_rows <- list.clean(latex_rows) # FIXME: comment out this line for final table

# create and save table
cat(top, do.call(paste, clean_latex_rows), bottom,
  file = "latex/tables/results_main_iter_08.tex") # NOTE: adjust

#####################################
# end
#####################################
