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
keep.fractions <- rep(1,length(samples))
samples.to.keep <- keep.fractions * samples #saves relative to last sample
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
initial.divergence <- 4

createGenomeObject <- function(fasta.files, gene.index = NULL, observed.expression.files = NULL, ...){

#Set up objects

inputDirectory <- "~/Documents/Research/gilchrist-lab-eb/input"
outputDirectory <- "~/Documents/Research/gilchrist-lab-eb/output/"

inputRestartFile <- paste0(outputDirectory, "Restart/rstart.round.", roundMy - 1, ".rst")
outputRestartFile <- paste0(outputDirectory, "Restart/rstart.round.", roundMy, ".rst")


#Genome object 
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

##Do multiple rounds of model fitting

roundMy <- 6

while(roundMy <= roundMax) {
  
  if(runModel) {
    if(roundMy == 1) {
      ##initial set up
      ##initialize parameter object
      ## set initial phi value
      parameter <- 
        initializeParameterObject(
          genome = genome,
          model = whichModel,
          sphi = initialSphi,
          num.mixtures = 1,
          gene.assignment = rep(1, genomeLength))
          
          divergence.iteration <- initial.divergence
        
    }
    
    else { ##anything other than first round
      parameter <- initializeParameterObject(
        genome = genome,
        init.with.restart.file = inputRestartFile,
        model = whichModel
        divergence.iteration <- 0
      )
    }
      

  
  
##Initialize Model Object
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
  setRestartSettings(mcmc = mcmc, filename = outputRestartFile, samples = adaptive.width*10, write.multiple = FALSE)
  
#Run MCMC on genome with parameter using model
  print(paste0 ("Starting MCMC round: ", roundMy))
  
runMCMC(mcmc = mcmc, genome = genome, model = model, ncores = ncores, divergence.iteration = divergence.iteration)

print(paste0("Finished MCMC round: ", roundMy))

##Process Output
print("Being processing output")
print ("Write Restart File? Done automatically")

print("Save parameter object")
outputFile <- paste0(outputDirectory, "R.Objects/parameter.round-", roundMy, ".Rdata")
writeParameterObject(parameter = parameter,
                     file = outputFile)

print("Save MCMC Object")
outputFile <- paste0(outputDirectory, "R.Objects/mcmc.round-", roundMy, ".Rdata")
writeMCMCObject(mcmc = mcmc,
                file = outputFile)

#END IF(RUN MODEL)
} else {
  
  print("Not running model, loading AnaCoDa objects instead")
  print(paste0("Round: ", roundMy))
  
  outputFile <- paste0(outputDirectory, "R.Objects/parameter.round-", roundMy, ".Rdata")
  parameter <- loadParameterObject(outputFile)
  
  outputFile <- paste0(outputDirectory, "R.Objects/mcmc.round-", roundMy, ".Rdata")
  mcmc <- loadMCMCObject(outputFile)
  
  model <- initializeModelObject(parameter = parameter, model = whichModel, with.phi = with.phi)
  
}

#Save parameters

if(saveParameters) {
  print("Save and Write Summaries of phi to memory")
  if(samples.to.keep[roundMy] > 100) {
    
    print("Write summaries of CSP")
    ##Write output to csv file
    getCSPEstimates(parameter = parameter, filename = paste0(outputDirectory, "Parameters/csp.round-", roundMy), mixture = 1, samples = samples.to.keep[roundMy], relative.to.optimal.codon = FALSE)
    
    print("Save Estimates of phi")
    phi.vals <- getExpressionEstimates(
      parameter = parameter, 
      gene.index = c(1:genomeLength),
      samples = samples.to.keep[roundMy],
      quantiles = c(0.025, 0.5, 0.975))
    
    phi.tibble <- bind_cols(isoform.ID = genomeIDs, phi.vals)
    write.table(phi.tibble,
                file = paste0(outputDirectory, "Parameters/phi.round-", roundMy, ".csv"),
                sep = ",",
                col.names = TRUE,
                quote = FALSE,
                row.names = FALSE)
  }
}

#Make plots
  if(makeFitPlots) {
    
    print("Printing traces")
    trace <- parameter$getTraceObject()
    
    outputFile <- paste0(outputDirectory, "Graphs/csp.traces.round-", roundMy, ".pdf")
    
    pdf(outputFile, title = "Codon Specific Traces: Reference Yeast")
    plot(x = trace, what = "Mutation", samples = 1000)
    plot(x = trace, what = "Selection", samples = 1000)
    plot(x = trace, what = "AcceptanceRatio", samples = 1000)
    dev.off()
    
   
  }
  
  roundMy <- roundMy + 1
}
