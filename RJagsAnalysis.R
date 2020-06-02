####################################################
#
# Author: Stefano Piacentini
# email:  stefano.piacentini@uniroma1.it
# Last Update: 02/06/2020
#
####################################################

####################################################

# Function for generation of pseudodata from the
# background template.
#
# expratio: E / E0
# stringer: string for identification of the exposure 

genPseudoData <- function(expratio = 1, stringer = '') {
    
    exposure <- 6786 * expratio

    bkg <- read.table(paste("./histograms/bkg1.txt", sep = ""))
    
    n   <- length(bkg[,1])
    
    x       <- bkg[,1]
    bkg     <- bkg[,2]
    norm    <- sum(bkg)
    rB_true <- norm
    
    
    datan <- 0
    for (i in 1:n) {
        lambda   <- exposure * rB_true * bkg[i] / norm
        datan[i] <- rpois(1, lambda)
    }
    
    df <- data.frame(x, datan)
    
    filename <- paste("./histograms/pseudodata1_", stringer, ".txt", sep = '')
    
    write.table(df, filename, sep = "\t", row.names = FALSE, col.names = FALSE)
    
    cat(sprintf("%s saved.\n",filename))
    
}

####################################################

# Function for jags analysis
#
# mass: mass of the DM candidate in MeV [string]
# process:  0 -> Migdal
#           1 -> Bremsstrahlung
#           2 -> Nuclear Recoil
# levelM:   12  -> n = 1,2
#           123 -> n = 1,2,3
# mod: name of the bug model (w/o file extension)
# save_plot: TRUE if you want to save the MCMC chain
#            plot into a png file.
# expratio: E / E0
# stringer: string for identification of the exposure 


RJagsAnalysis<- function(mass, process, levelM, nIter = 1e+4, mod = "signal",
                        expratio = 1,
                        stringer = '',
                        save_plot = TRUE
                        ) {
    
    knownModels = c("signal") # Add model here if you need to use other models
    
    if(!(mod %in% knownModels)) {
        cat(sprintf("ERROR. Uknown model.\n"))
        return(-1)
    }
    modello    <-  paste("Models/", mod, ".bug", sep = "")
    
    if (process == 0) {
        proc <- "Migdal"
    } else if (process ==1) {
        proc <- "BR"
    } else {
        proc <- "NuclearRecoil"
    }
    
    if (proc != "Migdal" && levelM == 123) {
        cat(sprintf("ERROR. The convention is that levelM must be 12 for BR and NR.\n"))
        return(-1)
    }
    
    if (levelM == 123) {
        levelM <- "n"
    } else {
        levelM <- ""
    }
    
    exposure <- 6786 * expratio

    signal <- read.table(paste("./histograms/", proc, "/sig_", mass, levelM, ".txt", sep = ""))
    
    bkg <- read.table("./histograms/bkg1.txt")
        
    x    <- bkg[,1]
    norm <- sum(bkg[,2])
    rB_true <- norm
    n <- length(bkg[,1])
    
    t <- exposure
    
    sig <- signal[,2]
    
    data <- read.table(paste("./histograms/pseudodata1_", stringer, ".txt", sep = ""))
    data <- data[,2]
    
    # Factor scaling with the exposure for bkg uncertainties 
    f <- 1/sqrt(expratio) 
    
    fitInput <- NULL
    
    fitInput$n         <- n
    fitInput$x         <- round(data)
    fitInput$t         <- t
    fitInput$bkg   <- bkg[,2]/norm

    # Error reproducing a background MC number of events ~ 1e+5
    err_rel <- sqrt(fitInput$bkg[length(fitInput$bkg)]/fitInput$bkg) * 0.03
    for(i in 1:n) {
        if(x[i] < 7) {
            err_rel[i] <- 0.3
        }
    }
    fitInput$bkg_e     <- f * fitInput$bkg * err_rel
    fitInput$rB_true   <- rB_true
    fitInput$srb       <- rB_true
    fitInput$sig       <- sig
        
    inits <- list(rB=rB_true, rS=0.0)

    jm <- jags.model(modello, fitInput, inits, quiet = TRUE)
    
    update(jm, 10000)
    
    to.sample = c('rB', 'rS')
    
    catena <- coda.samples(jm, to.sample, n.iter=nIter)  # sampling 

    # print(summary(catena))
    if(save_plot) {
        
        dir  <-  paste("plots/pseudo", stringer, sep = "")
        dir.create(dir, showWarnings = FALSE)
        dir  <-  paste(dir, "/", mod, sep = "")
        dir.create(dir, showWarnings = FALSE)
        dir  <-  paste(dir, "/", proc, sep = "")
        dir.create(dir, showWarnings = FALSE)

        file <-  paste(dir, "/", "ris_", mass, "_", levelM,".png", sep = "")


        png(file,
            width= 1600, height = 300*length(to.sample))

        plot(catena)

        dev.off()
    }
    
    return(catena)
    
}