#### Improvements & Pipelining #### input ND ####
#### 09.03.2023 LF ####

library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(phyloseq)
library(plyr)
library(beepr)

# Set the working directory to folder where sample files (.txt) are located.
setwd("demultiplex")
# Extract the headers of the sample table in the file.
table <- read.csv(file="nameslist.csv", header = T)
# Extract the "Descrription" column which contains the actual names of the sample
new.names <- table$newname
# Extract target .txt files and rename to Description name.
old.names <- table$oldname
file.rename(from=old.names, to=new.names)

# set file paths
### change main path accordingly 
path <- "~/Documents/work/BlueHoles/dataanalysis/sequencing_data/CO1/" 
path.demultiplex <- file.path(path, "demultiplex")
if(!dir.exists(path.demultiplex)) dir.create(path.demultiplex)
path.filt <- file.path(path, "filtered")
if(!dir.exists(path.filt)) dir.create(path.filt)
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)

list.files(path.demultiplex)

# generate matched lists of the forward and reverse read files, as well as parsing out the sample name
fnFs <- sort(list.files(path.demultiplex, pattern = "R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path.demultiplex, pattern = "R2.fastq.gz", full.names = TRUE))

# set primers for each dataset
#V3V4
FWD <- "CCTACGGGNGGCWGCAG"
REV <- "GACTACHVGGGTATCTAATCC" 

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

# Calculate number of reads containing forward and reverse primer sequences (considering all possible primer orientations. Only exact matches are found.).
# Only one set of paired end fastq.gz files will be checked (first sample in this case).
# This is is sufficient, assuming all the files were created using the same library preparation.

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))


##### remove primers with cutadapt#######

cutadapt <- "/Users/fruehel/.local/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R
  
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
    system2(cutadapt, args = c("-e 0.07 --discard-untrimmed", R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                               fnFs[i], fnRs[i], "--info-file=info.tsv", # input files
                               "-m", 20)) #remove too short or empty sequences
  }


rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2.fastq.gz", full.names = TRUE))

if(length(cutFs) == length(cutRs)) print("Forward and reverse files match. Go forth and explore")
if(length(cutFs) != length(cutRs)) stop("Forward and reverse files do not match. Better go back and have a check")

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "-")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

if(length(cutFs) <= 20) {
  fwd_qual_plots<-plotQualityProfile(cutFs) + 
    scale_x_continuous(breaks=seq(0,300,20)) + 
    scale_y_continuous(breaks=seq(0,40,5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_hline(yintercept = 30)
  rev_qual_plots<-plotQualityProfile(cutRs) + 
    scale_x_continuous(breaks=seq(0,300,20)) + 
    scale_y_continuous(breaks=seq(0,40,5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_hline(yintercept = 30)
} else {
  rand_samples <- sample(size = 20, 1:length(cutFs)) # grab 20 random samples to plot
  fwd_qual_plots <- plotQualityProfile(cutFs[rand_samples]) + 
    scale_x_continuous(breaks=seq(0,300,20)) + 
    scale_y_continuous(breaks=seq(0,40,5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_hline(yintercept = 30)
  rev_qual_plots <- plotQualityProfile(cutRs[rand_samples]) + 
    scale_x_continuous(breaks=seq(0,300,20)) + 
    scale_y_continuous(breaks=seq(0,40,5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_hline(yintercept = 30)
}
fwd_qual_plots
rev_qual_plots

jpeg(file="Qual_forward.jpg",res=300, width=15, height=8, units="in")
fwd_qual_plots
dev.off()
jpeg(file="Qual_reverse.jpg",res=300, width=15, height=8, units="in")
rev_qual_plots
dev.off()

# Assign filenames to the fastq.gz files of filtered and trimmed reads.

filtFs <- file.path(path, "filtered", basename(cutFs))
filtRs <- file.path(path, "filtered", basename(cutRs))

###Parameters vary for all sets
#V3V4 270/210 and 3/3
#V4 200/200 and 2/2 
#CO1 220/215

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, trimLeft = 5, 
                     truncLen = c(210,160),  maxN = 0, maxEE = 1,
                     truncQ = 2, minLen = 50, rm.phix = TRUE, 
                     compress = TRUE)
out <- readRDS("Track.rds")
saveRDS(out, file.path(path,"Track.rds"))
# Check if file names match

sample.namesF <- sapply(strsplit(basename(filtFs), "-"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "-"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(identical(sample.namesF, sample.namesR)) {print("Files are still matching.....congratulations")
} else {stop("Forward and reverse files do not match.")}
names(filtFs) <- sample.namesF
names(filtRs) <- sample.namesR

set.seed(100) # set seed to ensure that randomized steps are replicatable
errF <- learnErrors(filtFs, multithread=T)
errR <- learnErrors(filtRs, multithread=T)

pdf("Errorrates_fwd.pdf",width=5,height=5)
plotErrors(errF, nominalQ=TRUE)
dev.off()
pdf("Errorrates_rvs.pdf",width=5,height=5)
plotErrors(errR, nominalQ=TRUE)
dev.off()
saveRDS(errF, file.path(path,"errF.rds"))
saveRDS(errR, file.path(path,"errR.rds"))
beep()

# Apply the dada2's core sequence-variant inference algorithm:
# Set pool = pseudo", see https://urldefense.com/v3/__https://benjjneb.github.io/dada2/pool.html__;!!Nmw4Hv0!26LJAPXOS0UsplQMzcuviNShtzF3TlBAxf0UpgL6R8nseGHufTvEwEE-tMLRxCrPegutropwZaAR5GnCfIwLsNO-xX7sGy0wsednpQ0$ 
dadaFs <- dada(filtFs, err=errF, multithread=T, pool="pseudo")
names(dadaFs) <- sample.namesF
saveRDS(dadaFs, "dadaF.rds")
dadaRs <- dada(filtRs, err=errR, multithread=T, pool="pseudo")
names(dadaFs) <- sample.namesR
saveRDS(dadaRs, "dadaR.rds")
beep()

mergers <- readRDS("mergers.rds")
mergers <- mergePairs(dadaFs, filtFs, dadaRs,filtRs, 
                      minOverlap=20, maxMismatch=2, verbose=TRUE)
saveRDS(mergers, "mergers.rds")

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
saveRDS(seqtab, "seqtab.rds")
seqtab <- readRDS("seqtab.rds")

#get seq length for filtering
sl <- hist(nchar(getSequences(seqtab)), main="Distribution of sequence lengths", breaks=50)
axis(1,at=seq(300,500,by=10))
jpeg(file="Seqlength.jpg",res=450, width=15, height=8, units="in")
sl
dev.off()

#V4 375-380
#V3V4 400-430

seqtab_filtered <- seqtab[,nchar(colnames(seqtab)) %in% seq(375,380)]
saveRDS(seqtab_filtered, "seqtab_filtered.rds")
seqtab_filtered <- readRDS("/seqtab_filtered.rds")
seqtab_nochim <- removeBimeraDenovo(seqtab_filtered, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab_nochim, "seqtab_nochim.rds")

getN <- function(x)sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab_filtered), rowSums(seqtab_nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "seqtab", "seqtab_filt", "seqtab_nochim")
rownames(track) <- sample.namesF
track
write.csv(track,"SeqOV.csv")

# ### different tax assignments 

# taxa with Midori download LONGEST 
taxa_long <- assignTaxonomy(seqtab_nochim, "MIDORI2_LONGEST_NUC_GB255_CO1_DADA2.fasta.gz", multithread=TRUE)
saveRDS(taxa_long, "Taxa_long.rds")

# taxa with Midori download UNIQ
taxa_uniq <- assignTaxonomy(seqtab_nochim, "MIDORI2_UNIQ_NUC_GB255_CO1_DADA2.fasta.gz", multithread=TRUE)
saveRDS(taxa_uniq, "Taxa_uniq.rds")

# taxa with Midori download LONGEST SP
taxa_long_sp <- assignTaxonomy(seqtab_nochim, "MIDORI2_LONGEST_NUC_SP_GB255_CO1_DADA2.fasta.gz", multithread=TRUE)
saveRDS(taxa_long_sp, "Taxa_long_sp.rds")

# taxa with Midori download UNIQ SP
taxa_uniq_sp <- assignTaxonomy(seqtab_nochim, "MIDORI2_UNIQ_NUC_SP_GB255_CO1_DADA2.fasta.gz", multithread=TRUE)
saveRDS(taxa_uniq_sp, "Taxa_uniq_sp.rds")

# taxa with Midori download UNIQ SP
taxa_total_raw <- assignTaxonomy(seqtab_nochim, "MIDORI2_TOTAL_NUC_GB255_CO1_RAW.fasta.gz", multithread=TRUE)
saveRDS(taxa_total_raw, "Taxa_total_raw.rds")

##Chimera removal, taxonomic assignment and phyloseq object will be created on IBEX cluster
