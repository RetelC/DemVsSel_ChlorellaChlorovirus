#/bin/R
#############################
## Calculate Euclidean distances
## Cas Retel, cas. retel@eawag.ch
## 2019.10.14
#############################
## Standalone Rscript to calculate Euclidean distance between two
## population sequencing (Poolseq) samples. Writes a line of text 
## that identifies the two populations compared, and the distance between
## them, to a user-specified output file. 
## Inputs: 
##   infile = path to text file in the .sync format
##     (see https://sourceforge.net/p/popoolation2/wiki/Tutorial/), 
##      containing an arbitrary number of populations (length-1 character)
##   nsamp1, nsamp2 = two identifyers of which populations
##     from the input file to use (both length-1 numeric)
##   mac = minimum allele count: Variants with an alternative allele count 
##     of less than this number are removed from analysis (length-1 numeric)
##   mincov = minimum coverage: Positions that are covered by less than 
##     this number of reads are removed from analysis (length-1 numeric)
##   n.report = print progress every time this number of lines are completed, 
##     using base::message() functionality (length-1 numeric)
##   outfile = path to write output file to (length-1 character)

## Meant to be called from the command line using Rscript. Example usage: 
## Rscript calc_eucdist_from_sync_paired.R \
## ppl_dir/Pools2018-1.repls.q0.sync 1 2 2 10 100 ${outfile}

rm(list=ls())
## read in arguments
arguments <- list(
  infile = "/Users/reteladmin/Documents/HVInt/Project2_dem_vs_sel/genom/Virus2018-1/ppl_dir/Virus2018.repls.ds1k.sync",
  nsamp1 = 1, nsamp2 = 2,
  mac = 3, mincov = 200,
  n.report = 20,
  outfile = "/Users/reteladmin/Documents/HVInt/Project2_dem_vs_sel/genom/Virus2018-1/eucdist_dir/foursamps.eucdist"
)

## check input arguments formats
arguments <- as.list(commandArgs(trailingOnly=TRUE))
if(length(arguments) != 7){
  stop(paste0(
    "calc_eucdist_from_sync: ", length(arguments), 
    " arguments were provided where 7 were expected"
  ))
}
names(arguments) <- c(
  "infile", "nsamp1", "nsamp2", "mac", "mincov", "n.report", "outfile"
)
for(i in 2:6){
  arguments[[i]] <- as.numeric(arguments[[i]])
}

if((!is.numeric(arguments[["nsamp1"]])) & (length(arguments[["nsamp1"]]) == 1)){
  stop("calc_eucdist_from_sync: argument {nsamp1} is not a length-one numeric")
}
if((!is.numeric(arguments[["nsamp2"]])) & (length(arguments[["nsamp2"]]) == 1)){
  stop("calc_eucdist_from_sync: argument {nsamp2} is not a length-one numeric")
}
if((!is.numeric(arguments[["mac"]])) & (length(arguments[["mac"]]) == 1)){
  stop("calc_eucdist_from_sync: argument {mac} is not a length-one numeric")
}
if((!is.numeric(arguments[["mincov"]])) & (length(arguments[["mincov"]]) == 1)){
  stop("calc_eucdist_from_sync: argument {mincov} is not a length-one numeric")
}
if((!is.numeric(arguments[["n.report"]])) & (length(arguments[["n.report"]]) == 1)){
  stop("calc_eucdist_from_sync: argument {n.report} is not a length-one numeric")
}

## load required packages and custom functions
pkgs <- list("magrittr")
sapply(pkgs, require, character.only=T)
## ! two scripts below are available at ! 
## ! https://github.com/RetelC/TempDynamics_HostVirCoevol/ ! 
## ! dir_scripts should be changed to the path where they are locally stored !
dir_scripts <- "~/Documents/HVInt/scripts/"
source(paste0(dir_scripts, 'fs_syncCalculations.R'))
source(paste0(dir_scripts, 'writeReadSync.R'))

## report on progress
message(c(
  "###### calc_eucdist_from_sync initiated ######\n", 
  paste0(names(arguments), " : ", unlist(arguments), "\n")
))


## read in input file
snc.inp <- readSync(arguments[["infile"]], header_v = F)
n.inp <- nrow(snc.inp)
# head(snc.inp)

snc.samps <- snc.inp[, c(1:3, 3 + c(arguments[["nsamp1"]], arguments[["nsamp2"]]))]
rm(snc.inp)
## assess how many positions can be reliably evaluated
snc.samps.cov <- apply(snc.samps[, 4:5], 2, sumAF)
log.eval <- apply(snc.samps.cov, 1, (function(x) all(x > arguments[["mincov"]])))
n.eval <- sum(log.eval, na.rm=TRUE)

snc.cov <- snc.samps[log.eval, ]
rm(snc.samps)
## report on progress
message(c(
  "###### .sync file read in ######\n", 
  "number of input rows : ", n.inp, "\n", 
  "number of rows with coverage above {mincov}: ", n.eval, "\n"
))

## biallelic loci are loci that have an allele count of at least 
## {mac} for exactly two alleles. 
# log.biallelic <- apply(
#   snc.cov, 1, (function(row) checkBiallelic(snc_row = row, mac=arguments[["mac"]]))
# )
# n.biallelic <- sum(log.biallelic)
## multiallelic loci are loci that have an allele count of at least {mac} 
## for at least two alleles: 
log.multiallelic <- apply(
  snc.cov, 1, (function(row) checkMultiallelic(snc_row = row, mac=arguments[["mac"]]))
)
n.multiallelic <- sum(log.multiallelic)

## save up some space, report on progress
snc.variable <- snc.cov[log.multiallelic, ]
rm(snc.cov)
message(c(
  "###### filtered monoallelic sites ######\n", 
  "number of multiallelic sites : ", n.multiallelic, "\n"
))


## calculate major allele frequencies
## fill in allele frequencies
maf.variable <- matrix(0, nrow=n.multiallelic, ncol=ncol(snc.variable) - 3)
message(c(
  "--- calculating allele frequency ---"
))
for(i in 1:n.multiallelic){
  if((i %% arguments[["n.report"]]) == 0){
    message(c("--- ", i, " lines completed ---"))
  }
  maf.variable[i, ] <- as.numeric(calcMajorAF(snc.variable[i, ])[, 4:5])
}

## calculate euclidean distance
edist <- dist(t(maf.variable))

message(c(
  "Euclidean distance between ", colnames(snc.variable[4]), 
         " and ", colnames(snc.variable[5]), " :\t", edist, "\n"
))
writeLines(
  paste0("Euclidean distance between ", colnames(snc.variable[4]), 
         " and ", colnames(snc.variable[5]), " :\t", edist), 
  con=arguments[["outfile"]]
)
message(c(
  "###### DONE ######"
))
