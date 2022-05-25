#Load External
##Libraries
library(AnaCoDa)
library(Biostrings)
library(deming)
library(dplyr)
library(forcats)
library(ggplot2)
library(Hmisc)
library(knitr)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)

###Library help information
anacodaHelp <- ??AnaCoDa

##Define key parameters

roundInitial <- 1
roundMax <- 6
roundMy <- roundInitial

runModel <- TRUE #true if you want to fit model, false if you want to process objects already produced by model
makeFitPlots <- TRUE
comparisonCategories <- c("Mutation", "Selection")
compareToPrevious <- TRUE
makeComparisonPlots <- TRUE

saveParameters <- TRUE
saveSmallObjects <- TRUE
maxSmallObjectSize <- 5E5

adaptive.ratio <- rep(c(rep(c(1,0), 2), 0), 5) ##Alternate between adaptive and non-adaptive then end with a second non-adaptive round
samples <- rep(5, length(adaptive.ratio))*100 #steps determined by samples*thin
thinning <- 10
steps <- samples * thinning
adaptive.width <- 200/thinning
adaptive.steps <- steps * adaptive.ratio


whichModel <- "ROC"
sampleGenome <- TRUE
sampleSeed <- 6272001
geneSampleSize <- 100 #number of genes to use in fitting

with.phi <- FALSE
initialSphi <- 1
est.hyper = FALSE
fasta.file <- "orf_genomic_R64-3-1_20210421.fasta"
ncores = 3



##Set up objects

inputRestartFile <- paste0("Restart/rstart.round", roundMy - 1, ".rst")
outputRestartFile <- paste0("Restart/rstart.round", roundMy, ".rst")


### Create Genome Object
genome <- initializeGenomeObject(file = fasta.file)

geneIndex <- 1:length(getNames(genome))
                      

### If Sampling Genome
if(sampleGenome) {
  set.seed(sampleSeed)
  geneIndex <- sort(sample(geneIndex, geneSampleSize, replace = FALSE))
  #create reduced genome
  genome <- genome$getGenomeForGeneIndices(geneIndex, simulated = FALSE)
  rm(geneIndex)
  set.seed(NULL) # reset seed
}


genomeIDs <- getNames(genome) #really get IDs
genomeLength <- length(genomeIDs)

## Do Multiple Rounds of Model Fitting
startTime <- format(Sys.time(), "%Y-%m-%d_%H.%M")

roundMy <- roundInitial

while(roundMy <= roundMax) {
  
  inputRestartFile <- paste0("Restart/rstart.round", roundMy - 1, ".rst")
  outputRestartFile <- paste0("Restart/rstart.round", roundMy, ".rst")

if(runModel) {
  if(roundMy ==1) {
    ##initial set up
    ##initialize parameter object
    ## set initial phi value
    parameter <- 
      initializeParameterObject(
        genome = genome,
        model = whichModel,
        sphi = initialSphi,
        num.mixtures = 1,
        gene.assignment = rep(1, genomeLength)
      )
  
} else { ##anything other than first round
  parameter <- initializeParameterObject(
               genome = genome,
               init.with.restart.file = input.restart.file,
               model = whichModel
  )
}

## Initialize Model Object
model <- initializeModelObject(parameter = parameter, model = whichModel, with.phi = with.phi)

## Initialize MCMC and Model Object
mcmc <- initializeMCMCObject(samples = samples[roundMy],
                             thinning = thinning,
                             adaptive.width = adaptive.width, 
                             est.expression = TRUE,
                             est.csp = TRUE,
                             est.hyper = est.hyper)
mcmc$setStepsToAdapt(adaptive.steps[roundMy])

##Set up restart files
setRestartSettings(mcmc = mcmc, filename = output.restart.file, samples = adaptive.width*10, write.multiple = FALSE)
##run mcmc on genome with parameter using model
print(paste0("Starting MCMC round: ", roundMy))

sys.runtime <- system.time(
              runMCMC(mcmc = mcmc, genome = genome, model = model, ncores = ncores))

sys.runtime <- data.frame(Value=names(sys.runtime), Time = as.vector(sys.runtime))


print(paste0("Finished MCMC round: ", roundMy))

#Process Output
print("Begin processing output")

print("Save parameter objects")
outputFile <- paste0("R.Objects/mcmc.round-", roundMy, ".Rdata")
writeParameterObject(parameter = parameter, 
                     file = outputFile)

print("Save MCMC object")
outputFile <- paste0("R.Objects/mcmc.round-", roundMy, ".Rdata")
writeMCMCObject(mcmc = mcmc,
                file = outputFile)

##END OF IF(RUNMODEL)
}else{
print("Not running model, loading AnaCoDa objects")
print(paste0("Round: ", roundMy))

outputFile <- paste0("R.Objects/parameter.round-", roundMy, ".Rdata")

parameter <- loadParameterObject(outputFile)

outputFile <- paste0("R.Objects/mcmc.round-", roundMy, ".Rdata")
mcmc <- loadMCMCObject(outputFile)

model <- initializeModelObject(parameter = parameter, model = whichModel, with.phi = with.phi)
}

##Export parameters to csv
if(saveParameters) {
  print("Save and write summaries of phi to memory")
  
  print("Write summaries of CSP")
  ## save output to csv file
  getCSPEstimates(parameter = parameter, filename = paste0("Parameters/csp.round-", roundMy), mixture = 1,
   samples = samples)
  
  print("Save estimates of phi")
  phiVals <- getExpressionEstimates(
              parameter = parameter,
              gene.index = c(1:genomeLength),
              samples = samples,
              quantiles = c(0.025, 0.5, 0.975))
  
  phiTibble <- bind_cols(isoform.ID = genomeIDs, phiVals)
  write.table(phiTibble,
              file = paste0("Parameters/phi.round-", roundMy, ".csv"),
              sep = ",",
              col.names = TRUE,
              quote = FALSE,
              row.names = FALSE)
}
  
##Make plots
if(makeFitPlots) {
  outputFile <- paste0("Graphs/csp.traces.round-", roundMy, ".pdf")
  print("Printing traces")
  trace <- parameter$getTraceObject()
  
  
  pdf("Trace Plots: Rhoba")
  plot(x = trace, what = "Mutation", samples = 1000)
  plot(x = trace, what = "Selection", samples = 1000)
  plot(x = trace, what = "AcceptanceRatio", samples = 1000)
  dev.off()
  
}
roundMy <- roundMy + 1 
}

##Save small objects
if(saveSmallObjects) {
  objectSize <- object.size(get(x));
 # if(objectSize < maxSmallObjectSize) {
   # objectSize
 # }
    
#} else {
   # NA}
  