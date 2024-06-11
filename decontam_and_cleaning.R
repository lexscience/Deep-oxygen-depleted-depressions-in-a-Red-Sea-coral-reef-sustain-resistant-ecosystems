##Decontam final version##


#### PREVALENCE METHOD ####
BiocManager::install("decontam")
library(phyloseq)
library(ggplot2)
library(decontam)
library(dplyr)

ps <- readRDS("16s/")
(sample_data(ps))

sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev05$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")


contamdf.prev05_list <- contamdf.prev05 %>% filter(contaminant == "FALSE")
noncontamASVs <- as.list(rownames(contamdf.prev05_list))

asv_table <- read.csv2("unclean_COI_SOFT.csv", head=T)
asv_table_clean <- asv_table %>% filter(asv %in% noncontamASVs)
#dismiss E and C sample
excluded_vars <- c("E1N1_0", "C1N1_0")
asv_table_clean <- asv_table_clean %>% select(-one_of(excluded_vars))

# dismiss non target taxonomy hit sequences
asv_table_clean <- asv_table_clean %>% 
  # filter(Kingdom == "Bacteria") %>% 
  # filter(Order != "Chloroplast") %>% 
  # filter(Family != "Mitochrondia")
  
  #filter(Kingdom == "Eukaryota") %>%
  #filter(Supergroup != "NA") %>%
  #filter(Genus != "Zooxanthella") %>% 
  #filter(Class != "Metazoa") %>% 
  #filter(Class != "Fungi")
  
  # filter(Class != "Insecta") %>%
  # filter(Class != "Arachnida") %>%
  # filter(Class != "Aves") %>%
  # filter(Phylum != "Actinobacteria") %>%
  # filter(Phylum != "Proteobacteria") %>%
  # filter(Phylum != "No Match")

  filter(superkingdom != "Bacteria") %>%
  filter(class != "Insecta") %>%
  filter(class != "Arachnida") %>%
  filter(phylum != "Streptophyta")

write.csv(asv_table_clean, "clean_COI_SOFT.csv")

## remove singletons and doubletons in Excel
## check unclean table for all ASVs die in blanks vorkommen -> in excel tabelle und dann nochmal manuell cleanen

contams <- read.csv("V9/contams.csv", head=F)
contams <- as.list(contams$V1)

asv_table_clean <- read.csv("V9/clean_V9.csv", head=T)

asv_table_clean <- asv_table_clean %>% filter(!asv %in% contams)

write.csv(asv_table_clean, "clean_V9.csv")
