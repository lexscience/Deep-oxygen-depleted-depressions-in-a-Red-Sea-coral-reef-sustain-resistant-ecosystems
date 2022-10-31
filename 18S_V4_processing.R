### Blue Hole Project Sequence Data Analsysis 18S SSU rRNA gene - V4 ###
### Larissa Fruehe, October 2022###
run<-commandArgs()[7]
ncpus<-as.numeric(commandArgs()[8])

### starting with demultiplexed MiSeq data PE 2x300bp for ~420bp long V4 region of 18S SSU rRNA gene ### 


# 1) set file paths
path <- "~/Documents/work/BlueHoles/dataanalysis/v4_newtry/" 
path.demultiplex <- file.path(path, "demultiplex")
  if(!dir.exists(path.demultiplex)) dir.create(path.demultiplex)
path.filt <- file.path(path, "filtered")
  if(!dir.exists(path.filt)) dir.create(path.filt)
path.cut <- file.path(path, "cutadapt")
  if(!dir.exists(path.cut)) dir.create(path.cut)

list.files(path.demultiplex)


# 2) load libraries and dependancies
library(devtools)
library(dada2)
library(ShortRead)
library(Biostrings)
library(beepr)
library(doParallel)


# 3) check if all samples are named correctly
fnFs <- list.files(path.demultiplex, pattern = "_R1_001.fastq.gz", full.names = TRUE)
fnRs <- list.files(path.demultiplex, pattern = "_R2_001.fastq.gz", full.names = TRUE)
all.equal(sapply(strsplit(basename(fnFs), "."), `[`, 1),
          sapply(strsplit(basename(fnRs), "."), `[`, 1)) # TRUE


filtFs <- file.path(path.filt, gsub("L001","filt",basename(fnFs)))
filtRs <- file.path(path.filt, gsub("L001","filt",basename(fnRs)))

FWD <- "CCAGCASCYGCGGTAATTCC" 
REV <- "ACTTTCGTTCTTGATYRATGA"  

#to ensure we have the right primers, and the correct orientation of the primers on the reads, we will verify the presence and orientation of these primers in the data.
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

#The presence of ambiguous bases (Ns) in the sequencing reads makes accurate mapping of short primer sequences difficult. Next we are going to “pre-filter” the sequences just to remove those with Ns, but perform no other filtering.
fnFs.filtN <- file.path(path.filt, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path.filt, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 2, multithread = TRUE)
beep()

##### remove primers with cutadapt#######

cutadapt <- "/Users/fruehel/.local/bin/cutadapt" #cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R


if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, #-n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], #output files
                             fnFs.filtN[i], fnRs.filtN[i], #input files
                             "-m", 20)) #remove too short or empty sequences
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
########################################################################
# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

readFastq(cutFs[1:2])


# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)
plotQualityProfile(cutFs[1:24])
plotQualityProfile(cutRs[1:24])

### filter and trim
filtFs <- file.path(path.filt, "filtered", basename(cutFs))
filtRs <- file.path(path.filt, "filtered", basename(cutRs))
