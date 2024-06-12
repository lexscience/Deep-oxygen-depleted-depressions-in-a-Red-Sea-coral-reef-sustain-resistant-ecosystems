####################################
######      Data Analysis     ######
####################################
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(tibble)

table <- read.table("assembledTables_Final/clean_V3V4.csv", head = T, sep=",")
env <- read.csv("assembledTables_Final/samdf.csv", head=T, sep=",") %>% 
  select(Sample, Site)

matrix_all<-t(table[,2:22]) 
colnames(matrix_all) <- table$asv

taxonomy <- table[,c(1,23:29)]

relab_ab<-decostand(matrix_all,"total")
relab_df<-data.frame(relab_ab) %>% rownames_to_column("Sample") %>%
  gather("asv","ra",-Sample)

relab_df_env <- left_join(relab_df, env) %>%
  select(Site, Sample, asv, ra) %>%
  group_by(Site, Sample, asv) %>%
  summarise_each(funs(sum))

#BH1
D1_relab_df_env <- relab_df_env %>% filter(Site == "D1")
#BH2
D2_relab_df_env <- relab_df_env %>% filter(Site == "D2")
#BH3
S1_relab_df_env <- relab_df_env %>% filter(Site == "S1")
#BH4
S2_relab_df_env <- relab_df_env %>% filter(Site == "S2")


################## PHYLUM 

relab_df_taxonomy <-left_join(D1_relab_df_env,taxonomy) %>%
  select(ra, Sample, Phylum) %>%
  group_by(Sample, Phylum, Site) %>%
  summarize_each(funs(sum))

# relab_df_taxonomy_unass <- relab_df_taxonomy %>% 
#   mutate(Phylum2=ifelse(Phylum=="","unassigned bacteria",as.character(Phylum)))%>% 
#   group_by(Station, Category, Site, Phylum2) %>%
#   summarize(ra2=sum(ra)) %>%
#   as.data.frame()

filtered_df <- relab_df_taxonomy %>%
  group_by(Phylum) %>%
  summarize(ra2=mean(ra)) %>% 
  as.data.frame()

attach(filtered_df)
filtered_df_s <- filtered_df[order(-ra2),]
detach(filtered_df)
filtered_df_s_top15 <- filtered_df_s[1:20,1]


relab_df_taxonomy_plot_andere <- relab_df_taxonomy %>%
  mutate(Phylum2=ifelse(Phylum %in% filtered_df_s_top15, Phylum, "Other")) %>% 
  group_by(Sample, Phylum2) %>%
  summarize(ra2=sum(ra)) %>%
  as.data.frame() %>%
  left_join(env) %>%
  group_by(Sample, Phylum2) %>%
  summarize(ra3=mean(ra2))

colourCount = length(unique(relab_df_taxonomy_plot_andere$Phylum2))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

#BH1
relab_df_taxonomy_plot_andere$Sample <- factor(relab_df_taxonomy_plot_andere$Sample,levels=c("D1N4_45", "D1N3_90", "D1N2_125", "D1N1_419"), 
                                               labels = c("45 m", "90 m", "125 m", "419 m"))
#BH2
relab_df_taxonomy_plot_andere$Sample <- factor(relab_df_taxonomy_plot_andere$Sample, levels=c("D2N5_5", "D2N1_40", "D2N2_125", "D2N3_200", "D2N6_325", "D2N4_620"), 
                                               labels = c("5 m", "40 m", "125 m", "200 m", "325 m", "620 m"))
#BH3
relab_df_taxonomy_plot_andere$Sample <- factor(relab_df_taxonomy_plot_andere$Sample, levels=c("S1N5_3", "S1N4_10", "S1N3_20", "S1N2_30", "S1N1_45"), 
                                               labels = c("3 m", "10 m", "20 m", "30 m", "45 m"))
#BH4
relab_df_taxonomy_plot_andere$Sample <- factor(relab_df_taxonomy_plot_andere$Sample,  levels=c("S2N6_5", "S2N1_5", "S2N5_15", "S2N4_25", "S2N3_40", "S2N2_50"), 
                                               labels = c("3 m", "5 m", "15 m", "25 m", "40 m", "50 m"))
##PLOT Supergroup LEVEL
plot_Phylum2 <-ggplot(relab_df_taxonomy_plot_andere, aes(x=Sample, y=ra3, fill=Phylum2)) +
  geom_bar(stat="identity", position="fill") +
  #facet_grid(~Site,space="free", scales = "free") +
  theme(axis.text.x=element_text(size=18),
        axis.ticks.x=element_blank(), 
        strip.text.x = element_text(size = 18, colour = "grey"))+
  guides(fill=guide_legend(title="Phylum", ncol=1)) +
  theme_classic() +
  labs(x="Sample", y="Relative Read Abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = getPalette(colourCount))

plot_Phylum2
ggsave("Taxonomy_BH4.pdf", plot_Phylum2, height = 10, width = 15)

################## Class 

relab_df_taxonomy<-left_join(D1_relab_df_env,taxonomy) %>%
  select(ra, Site, Class) %>%
  group_by(Sample, Site, Class) %>%
  summarize_each(funs(sum))

# relab_df_taxonomy_unass <- relab_df_taxonomy %>% 
#   mutate(Phylum2=ifelse(Phylum=="","unassigned bacteria",as.character(Phylum)))%>% 
#   group_by(Station, Category, Site, Phylum2) %>%
#   summarize(ra2=sum(ra)) %>%
#   as.data.frame()

filtered_df <- relab_df_taxonomy %>%
  group_by(Class) %>%
  summarize(ra2=mean(ra)) %>% 
  as.data.frame()

attach(filtered_df)
filtered_df_s <- filtered_df[order(-ra2),]
detach(filtered_df)
filtered_df_s_top15 <- filtered_df_s[1:15,1]

relab_df_taxonomy_plot_andere <- relab_df_taxonomy %>%
  mutate(Class2=ifelse(Class %in% filtered_df_s_top15, Class, "Other")) %>% 
  group_by(Sample, Class2) %>%
  summarize(ra2=sum(ra)) %>%
  as.data.frame() %>%
  left_join(env) %>%
  group_by(Sample, Class2) %>%
  summarize(ra3=mean(ra2))

colourCount = length(unique(relab_df_taxonomy_plot_andere$Class2))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

##PLOT Supergroup LEVEL
plot_Class2 <-ggplot(relab_df_taxonomy_plot_andere, aes(x=Sample, y=ra3, fill=Class2)) +
  geom_bar(stat="identity", position="fill") +
  #facet_grid(.~ Sample,space="free_x", switch="y") +
  theme(axis.text.x=element_text(size=18),
        axis.ticks.x=element_blank(), 
        strip.text.x = element_text(size = 18, colour = "grey"))+
  guides(fill=guide_legend(title="Class", ncol=1)) +
  theme_classic() +
  labs(x="Site", y="Relative Read Abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = getPalette(colourCount))

plot_Class2

################## Order 

relab_df_taxonomy<-left_join(D1_relab_df_env,taxonomy) %>%
  select(ra, Site, Order) %>%
  group_by(Sample, Site, Order) %>%
  summarize_each(funs(sum))


# relab_df_taxonomy_unass <- relab_df_taxonomy %>% 
#   mutate(Phylum2=ifelse(Phylum=="","unassigned bacteria",as.character(Phylum)))%>% 
#   group_by(Station, Category, Site, Phylum2) %>%
#   summarize(ra2=sum(ra)) %>%
#   as.data.frame()

filtered_df <- relab_df_taxonomy %>%
  group_by(Order) %>%
  summarize(ra2=mean(ra)) %>% 
  as.data.frame()

attach(filtered_df)
filtered_df_s <- filtered_df[order(-ra2),]
detach(filtered_df)
filtered_df_s_top15 <- filtered_df_s[1:15,1]

relab_df_taxonomy_plot_andere <- relab_df_taxonomy %>%
  mutate(Order2=ifelse(Order %in% filtered_df_s_top15, Order, "Other")) %>% 
  group_by(Sample, Order2) %>%
  summarize(ra2=sum(ra)) %>%
  as.data.frame() %>%
  left_join(env) %>%
  group_by(Sample, Order2) %>%
  summarize(ra3=mean(ra2))

colourCount = length(unique(relab_df_taxonomy_plot_andere$Order2))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

##PLOT Supergroup LEVEL
plot_Order2 <-ggplot(relab_df_taxonomy_plot_andere, aes(x=Sample, y=ra3, fill=Order2)) +
  geom_bar(stat="identity", position="fill") +
  #facet_grid(.~ Site,space="free_x", switch="y") +
  theme(axis.text.x=element_text(size=18),
        axis.ticks.x=element_blank(), 
        strip.text.x = element_text(size = 18, colour = "grey"))+
  guides(fill=guide_legend(title="Order", ncol=1)) +
  theme_classic() +
  labs(x="Site", y="Relative Read Abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = getPalette(colourCount))

plot_Order2

library(ggpubr)
V3V4_all <- ggarrange(plot_Phylum2, plot_Class2, plot_Order2, ncol = 1, nrow = 3, labels = c("A", "B", "C"))
ggsave("Bacteria_Taxonomy.pdf",V3V4_all,width=13,height=18)


################### filter for Dinoflagellates

table <- read.table("assembledTables_Final/clean_V4.csv", head = T, sep=",")
env <- read.csv("assembledTables_Final/samdf.csv", head=T, sep=",") %>% 
  select(Sample, Site)

table_dinos <- table %>% filter(Division=="Dinoflagellata")

matrix_all<-t(table_dinos[,2:22])
colnames(matrix_all) <- table_dinos$asv

taxonomy <- table_dinos[,c(1,23:29)]

relab_ab<-decostand(matrix_all,"total")
relab_df<-data.frame(relab_ab) %>% rownames_to_column("Sample") %>%
  gather("asv","ra",-Sample)

relab_df_env <- left_join(relab_df, env) %>%
  select(Site, Sample, asv, ra) %>%
  group_by(Site, Sample, asv) %>%
  summarise_each(funs(sum))

#BH1
D1_relab_df_env <- relab_df_env %>% filter(Site == "D1")
#BH2
D2_relab_df_env <- relab_df_env %>% filter(Site == "D2")
#BH3
S1_relab_df_env <- relab_df_env %>% filter(Site == "S1")
#BH4
S2_relab_df_env <- relab_df_env %>% filter(Site == "S2")

relab_df_taxonomy <-left_join(S2_relab_df_env,taxonomy) %>%
  select(ra, Sample, Order) %>%
  group_by(Sample, Order, Site) %>%
  summarize_each(funs(sum))

filtered_df <- relab_df_taxonomy %>%
  group_by(Order) %>%
  summarize(ra2=mean(ra)) %>% 
  as.data.frame()

attach(filtered_df)
filtered_df_s <- filtered_df[order(-ra2),]
detach(filtered_df)
filtered_df_s_top15 <- filtered_df_s[1:10,1]


relab_df_taxonomy_plot_andere <- relab_df_taxonomy %>%
  mutate(Order2=ifelse(Order %in% filtered_df_s_top15, Order, "Other")) %>% 
  group_by(Sample, Order2) %>%
  summarize(ra2=sum(ra)) %>%
  as.data.frame() %>%
  left_join(env) %>%
  group_by(Sample, Order2) %>%
  summarize(ra3=mean(ra2))

colourCount = length(unique(relab_df_taxonomy_plot_andere$Order2))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

#BH1
relab_df_taxonomy_plot_andere$Sample <- factor(relab_df_taxonomy_plot_andere$Sample,levels=c("D1N4_45", "D1N3_90", "D1N2_125", "D1N1_419"), 
                                               labels = c("45 m", "90 m", "125 m", "419 m"))
#BH2
relab_df_taxonomy_plot_andere$Sample <- factor(relab_df_taxonomy_plot_andere$Sample, levels=c("D2N5_5", "D2N1_40", "D2N2_125", "D2N3_200", "D2N6_325", "D2N4_620"), 
                                               labels = c("5 m", "40 m", "125 m", "200 m", "325 m", "620 m"))
#BH3
relab_df_taxonomy_plot_andere$Sample <- factor(relab_df_taxonomy_plot_andere$Sample, levels=c("S1N5_3", "S1N4_10", "S1N3_20", "S1N2_30", "S1N1_45"), 
                                               labels = c("3 m", "10 m", "20 m", "30 m", "45 m"))
#BH4
relab_df_taxonomy_plot_andere$Sample <- factor(relab_df_taxonomy_plot_andere$Sample,  levels=c("S2N6_5", "S2N1_5", "S2N5_15", "S2N4_25", "S2N3_40", "S2N2_50"), 
                                               labels = c("3 m", "5 m", "15 m", "25 m", "40 m", "50 m"))
##PLOT Supergroup LEVEL
plot_Phylum2 <-ggplot(relab_df_taxonomy_plot_andere, aes(x=Sample, y=ra3, fill=Order2)) +
  geom_bar(stat="identity", position="fill") +
  #facet_grid(~Site,space="free", scales = "free") +
  theme(axis.text.x=element_text(size=18),
        axis.ticks.x=element_blank(), 
        strip.text.x = element_text(size = 18, colour = "grey"))+
  guides(fill=guide_legend(title="Dinoflagellata Division", ncol=1)) +
  theme_classic() +
  labs(x="Sample", y="Relative Read Abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = getPalette(colourCount))

plot_Phylum2
ggsave("Taxonomy_BH4_Dino.pdf", plot_Phylum2, height = 10, width = 15)
