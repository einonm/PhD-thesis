# set a default repository for packages
local({r <- getOption("repos")
       r["CRAN"] <- "http://www.stats.bris.ac.uk/R/"
       options(repos=r)
})

if (!require("BiocManager", quietly = TRUE))
    install.packages('BiocManager', dependencies=T)

BiocManager::install("qvalue")
