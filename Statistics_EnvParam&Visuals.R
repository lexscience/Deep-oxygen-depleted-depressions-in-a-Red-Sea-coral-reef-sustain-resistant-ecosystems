#Klein et al 2024#
#Non-sequence data analyses#

#Pearson correlations analysing the relationship between environmental parameters#

#Load data file#
df <- read.csv("/Users/.../Supplementary_Data_9.csv", header = TRUE)
head(df[, c("Site", "Sample_depth", "O2_umol_kg_seabird", "DIC_umol_kg", "pH_total", "pH_total_log")], 50)
str(df)
#Pearson correlation analysing dissolved inorganic carbon and dissolved oxgen#
res <- cor.test(df$DIC_umol_kg, df$O2_umol_kg_seabird, 
                method = "pearson")
#view results
res

#Pearson correlation analysing pH total (log) and dissolved oxgen#
res2 <- cor.test(df$pH_total_log, df$O2_umol_kg_seabird, 
                method = "pearson")
#view results
res2

#Mann-Whitney test analysing fish swimming speeds between site 2 and the open-water reference sites#

#Load data file#
df3 <-read.csv("/Users/.../Supplementary_Data_10.csv", header = TRUE)
head(df3[, c("Swim_sp_reference_sites", "Swim_sp_depression_sites")], 30)
str(df3)

#Mann-Whitney test, where alternate hypothesis is reference > depression
res3 <- wilcox.test(df3$Swim_sp_reference_sites, df3$Swim_sp_depression_sites, alternative = "greater", paired = FALSE)

#view results
res3
