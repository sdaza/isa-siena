########################################################
# read anylogic output
# author: sebastian daza
# version: 1.00
#######################################################

# libraries
library(data.table)

# path
path <- "/Users/sdaza/Desktop/model/"
fpath <- "/Users/sdaza/Desktop/model/output/socialNorms/"

#+ read files
behavior.files <- list.files(path = fpath, pattern = "behavior")
b <- lapply(paste0(fpath, "/", behavior.files), fread, sep=",")
b <- rbindlist(b)

network.files <- list.files(path = fpath, pattern = "network")
n <- lapply(paste0(fpath, "/", network.files), fread, sep=",")
n <- rbindlist(n)

# no need to load these files
parameter.files <- list.files(path = fpath, pattern = "parameter")
p <- lapply(paste0(fpath, "/", parameter.files), fread, sep=",")
p <- rbindlist(p)

print(replications <- max(p$replication))
p <- p[, replication := NULL][!duplicated(p)]
setkey(p, iteration)
print(p)

# position (distance)
position.files <- list.files(path = fpath, pattern = "position")
d <- lapply(paste0(fpath, "/", position.files), fread, sep=",")
d <- rbindlist(d)

# categorize behavior variable
b$mybehavior <- cut(b$behavior, breaks = 10, labels = 1:10,
  ordered_result = TRUE)
table(b$mybehavior, useNA = "ifany")

# save files
save(b, n, d, p, file = paste0(path, "R/data/socialNorms.Rdata"))
