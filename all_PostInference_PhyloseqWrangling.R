####phyloseq post inference wrangling

library(phyloseq)
library(ggplot2)

library("devtools")
devtools::install_github("benjjneb/dada2")

path <- "/CO1"
setwd(path)

seqtab_nochim <- readRDS("seqtab_nochim.rds")
names <- read.csv2("BOLD.csv", head=T)
colnames(seqtab_nochim) <- names$ID
taxa <- as.matrix(read.csv2("SOFT.csv", head=T, row.names = 1))
row.names(seqtab_nochim)
samdf <- read.csv2("samdf.csv", row.names = 1, head=T)



ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(taxa))


dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
saveRDS(ps, "ps_CO1_SOFT.rds")

unname(taxa)


plot_richness(ps, x="Category", measures=c("Shannon", "Simpson"))

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="Category", title="Bray NMDS")

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:15]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Category", fill="Order") + facet_wrap(~Category, scales="free_x")

df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

a.ord <- ordinate(ps, "PCoA", "bray")
plot_ordination(ps, a.ord, type="samples", color="Category")


#####export files to further work on 

OTU1 = as(otu_table(ps), "matrix")
# transpose if necessary
if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)

taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

Taxa1 = as(tax_table(ps), "matrix")
# transpose if necessary
if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
# Coerce to data.frame
Taxadf = as.data.frame(Taxa1)
OTUdf <- t(OTUdf)
write.table(OTUdf,file="PS_ASVs_SOFT.csv",sep=",")
write.table(Taxadf,file="PS_Taxa_SOFT.csv",sep=",")

#### export sequence fasta file from seqtab_nochim
seqtab <- readRDS("V4/seqtab.rds")
seqtab_nochim <- readRDS("V4/seqtab_nochim.rds")

asv_seqs <- colnames(seqtab_nochim)
asv_headers <- vector(dim(seqtab_nochim)[2], mode="character")
for (i in 1:dim(seqtab_nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "FinSeq_COI.fasta")




