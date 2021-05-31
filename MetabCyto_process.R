library(xlsx)
library(Rtsne)
library(clusterProfiler)
library(ggplot2)
library(ggdendro)
options(stringsAsFactors = F)

# Prepare Metabolomica Data of COVID-19 (, non-COVID-19) Patients and Healthy Controls ---------------------------

# Reference files
m_c.tar <- read.xlsx(paste0(wd1, "metabolite.HMDB.KEGG.HULAB.xlsx"),  # not provided in Supp files
                     sheetIndex = 2, startRow = 1, header = T, as.data.frame = T)
m_c.untar <- read.xlsx(paste0(wd1, "metabolite.HMDB.KEGG.HULAB.xlsx"), # not provided
                       sheetIndex = 3, startRow = 1, header = T, as.data.frame = T)
m_c.all <- rbind(data.frame("Metabolite" = m_c.tar$Metabolite, "KEGG" = m_c.tar$KEGG),
                 data.frame("Metabolite" = m_c.untar$Metabolite, "KEGG" = m_c.untar$KEGG))
# metabolits in each pathways
path_metab <- read.table('path.cpd.txt', header = F, sep = '\t') # not provided in Supp files
path_metab$V2 <- gsub("cpd:", "", path_metab$V2)
# metabolic pathways
path_intro <- read.xlsx("has_MetabolicPathway.xlsx",  # not provided in Supp files
                        sheetIndex = 2, startRow = 1, as.data.frame = T, header = F, check.names = T)

# Targeted metabolomics data of COVID-19 patients and healthy controls
df.tar.metab <- read.csv(paste0(wd2, "COVID-19.targeted.QCmad-TIC.csv"), # SuppData 3
                         sep = ",", header = T, row.names = 1, check.names = F)
df.tar.metab <- df.tar.metab[, -grep("arnitine", colnames(df.tar.metab))]
cpd.tar <- m_c.tar$KEGG[match(colnames(df.tar.metab), m_c.tar$Metabolite)]

# Untargeted metabolomics data of COVID-19 patients and healthy controls
df.untar.metab <- read.csv(paste0(wd2, "COVID-19.untargeted.QCmad-TIC.csv"), # SuppData 3
                           header = T, sep = ",", row.names = 1, check.names = F)
cpd.untar <- m_c.untar$KEGG[match(colnames(df.untar.metab), m_c.untar$Metabolite)]
df.untar.metab <- df.untar.metab[, -na.omit(match(cpd.tar[which(cpd.tar != "-")], cpd.untar))]
cpd.untar <- cpd.untar[-na.omit(match(cpd.tar[which(cpd.tar != "-")], cpd.untar))]

# Targeted metabolomics data of non-COVID-19 patients
df.noncovid <- read.xlsx("COVID-19.targeted.nonCOVID.QCmad-TIC.xlsx", # SuppData 3
                         sheetIndex = 1, startRow = 1, row.names = 1, header = T, as.data.frame = T, check.names = F)
df.noncovid <- df.noncovid[, 1:20]
df.noncovid <- data.frame(t(df.noncovid), check.names = F)

# for COVID-19 patients and healthy controls
df.metab <- cbind(df.tar.metab, df.untar.metab)
# COVID-19 patients and  healthy controls
df.metab <- df.metab[c(1:37, 96:112), ]

# for COVID-19, non-COVID-19 patients and healthy controls
df.metab <- rbind(df.tar.metab[c(1:37, 96:112),],
                  df.noncovid[, match(colnames(df.tar.metab), colnames(df.noncovid))])

# t-SNE Plot of COVID-19 (, non-COVID-19) patients and Healthy Controls ---------------------------

df.PCA <- log10(data.frame(df.metab))
# t-SNE analysis
tsne_out <- Rtsne(
  df.PCA,
  normalize = T,
  dims = 2,
  pca = F,
  perplexity = 5,
  theta = 0.01,
  max_iter = 1000
  #random_state = 9999
)
tsne.result <- as.data.frame(tsne_out$Y, row.names = rownames(df.metab))

# for COVID-19 patients and healthy controls
tsne.result$group <- rep(c("Severe", "Mild", "Healthy"), c(23, 14, 17))
# for COVID-19, non-COVID-19 patients and healthy controls
# tsne.result$group <- rep(c("Severe", "Mild", "Healthy", "non-COVID-19"), c(23, 14, 17, 20))
tsne.result$name <- rownames(tsne.result)

# t-SNE plot
ggplot(tsne.result, aes(x = V1, y = V2, fill = group))+
  geom_point(size = 2.3, shape = 21, color = "black") +
  
  # for COVID-19 patients and healthy controls
  scale_fill_manual(values = alpha(c("Healthy" = "#6fe7dd",
                                     "Mild" = "#3490de", "Severe" = "#6639a6"), 1)) +
  # for COVID-19, non-COVID-19 patients and healthy controls
  # scale_fill_manual(values = alpha(c("Healthy" = "#6fe7dd", "Mild" = "#3490de", 
  #                                    "Severe" = "#6639a6", "non-COVID-19" = "grey60"), 1)) +
  
  stat_ellipse(level = 0.95, show.legend = T, linetype = 2, color = "grey20")+
  labs(x = "tSNE-1", y = "tSNE-2", fill = "Group") + 
  theme_bw() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 1) +
  theme(axis.title.x = element_text(size = 8, hjust = 0.5, face = "plain"),
        axis.title.y = element_text(size = 8, face = "plain"), 
        axis.text.y = element_text(size = 8, face = "plain", colour = "black"),
        axis.text.x = element_text(size = 8, face = "plain", colour = "black")) +
  theme(legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8), 
        legend.spacing.y = unit(0.2, 'cm'), 
        legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'))

# Significantly Altered Metabolites between H/M/S ---------------------------

# for COVID-19 patients and healthy controls
df.metab$Degree <- rep(c("Severe", "Mild", "Healthy"), c(23, 14, 17))
# for COVID-19, non-COVID-19 patients and healthy controls
# df.metab$Degree <- rep(c("Severe", "Mild", "Healthy", "non-COVID-19"), c(23, 14, 17, 20))

# compare two groups
type1 <- "non-COVID-19"
type2 <- "Healthy"
df.test <- rbind(df.metab[which(df.metab$Degree == type1),], 
                 df.metab[which(df.metab$Degree == type2),])
matab_num <- ncol(df.test) - 1
test_results <- data.frame(Metabolites = colnames(df.test)[1:matab_num])
# Mann-Whitney U test
test_results$wilcox <- apply(df.test[, 1:matab_num], 2,
                             function(x) unlist(wilcox.test(as.numeric(x) ~ df.test$Degree, data = df.test)[3]))
# adjust p value using BH method
test_results$wilcox_BH <- p.adjust(test_results$wilcox, method = "BH")
# Log2 (FC)
cn <- paste("FC_", type1, "/", type2, sep = "")
test_results[cn] <- apply(df.test[,1:matab_num], 2, 
                          function(x) 
                            mean(as.numeric(x[which(df.test$Degree == type1)]))/
                            mean(as.numeric(x[which(df.test$Degree == type2)])))
test_results$LOG2FC <- log2(test_results[,4])

# KEGG enrichment analysis ---------------------------

group.name <- c("MvsH", "SvsH", "SvsM")
SMH_kegg <- data.frame()
for (i in c(1:2)) {
  metab.wilcox <- read.xlsx(paste0(wd3, "Serum.metabolite.MannWhitneyU.xlsx"),  # Source Data
                            sheetIndex = i, startRow = 1, as.data.frame = T, header = T, check.names = F)
  metab.wilcox$cpd <- m_c.all$KEGG[match(metab.wilcox$Metabolite, m_c.all$Metabolite)]
  # Significantly altered metabolites
  metab.wilcox <- metab.wilcox[metab.wilcox$Utest_BH < 0.05 & 
                                 abs(metab.wilcox$LOG2FC) > 0.3219281, ]
  # KEGG enrichment analyses
  x <- enricher(metab.wilcox$cpd, 
                TERM2GENE = path_metab, TERM2NAME = path_intro,
                minGSSize = 2, pvalueCutoff = 0.1, pAdjustMethod = "BH")
  kegg_table <- na.omit(as.data.frame(x))
  kegg_table <- kegg_table[kegg_table$Count > 2,]
  
  # calculate the absolute mean log2(FC) of metabolites in each pathway
  path_avelog2FC <- function(x){
    metabs <- unlist(strsplit(x[8], "[/]"))
    mean(abs(metab.wilcox[match(metabs, metab.wilcox$cpd),5]))
  }
  kegg_table$metabs_mean <- apply(kegg_table, 1, path_avelog2FC)
  # calculate fold enrichment
  kegg_table$FoldEnrich <- apply(kegg_table, 1, 
                                 function(x) as.numeric(unlist(strsplit(x[3], '[/]'))[1])/
                                   as.numeric(unlist(strsplit(x[3], '[/]'))[2])*
                                   as.numeric(unlist(strsplit(x[4], '[/]'))[2])/
                                   as.numeric(unlist(strsplit(x[4], '[/]'))[1]))
  # add group
  kegg_table$Group <- group.name[i]
  # combine
  SMH_kegg <- rbind(SMH_kegg, kegg_table)
}

path_all <- data.frame(row.names = unique(SMH_kegg$Description), 
                       MH_abslog2FC = SMH_kegg$metabs_mean[match(unique(SMH_kegg$Description), 
                                                                 SMH_kegg$Description[SMH_kegg$Group == "MvsH"])], 
                       SH_abslog2FC = SMH_kegg$metabs_mean[match(unique(SMH_kegg$Description), 
                                                                 SMH_kegg$Description[SMH_kegg$Group == "SvsH"])])
# order the pathways
path_all$MH_abslog2FC[which(is.na(path_all$MH_abslog2FC))] <- 0
path.ord <- hclust(dist(path_all, method = "euclidean"), method = "average")
path.ord <- rownames(path_all)[path.ord$order]
SMH_kegg$Description <- factor(SMH_kegg$Description, levels = rev(path.ord))

# point plot of KEGG enrichment results
ggplot(SMH_kegg, aes(Group, Description)) +
  geom_point(aes(fill = metabs_mean, size = -log10(p.adjust)), color = "black", shape = 21) +
  scale_size(range = c(1, 3.2), breaks = c(2,5,8)) +
  scale_fill_viridis_c(option = "A", direction = -1, begin = 0.4, breaks = c(0.5, 1.5, 2.5)) +
  theme_dendro() + 
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  coord_fixed(ratio = 0.7) +
  labs(x = "", y = "", title = "", 
       fill = "abs. Log2 (Fold change)", size = "-Log10 (P value)") +
  theme(axis.text.x = element_text(size = 6, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 6, face = "plain", colour = "black")) +
  theme(legend.text = element_text(size = 6, face = "plain", colour = "black"), 
        legend.title = element_text(size = 6, face = "plain", colour = "black"), 
        legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm"))

# Metabolites Altered along Disease Severity ---------------------------

metab.mhu.bind <- data.frame(row.names = colnames(df.metab))
for (i in c(1:3)) {
  metab.mhu <- read.xlsx(paste0(wd3, "Serum.metabolite.MannWhitneyU.xlsx"), # Source dara
                         sheetIndex = i, startRow = 1, as.data.frame = T, header = T, check.names = F)
  metab.mhu$Group <- group.name[i]
  metab.mhu.bind <- data.frame(metab.mhu.bind, 
                               metab.mhu[match(rownames(metab.mhu.bind), metab.mhu$Metabolite), -1])
}

# increased along disease severity
up.up <- metab.mhu.bind[which(metab.mhu.bind$LOG2FC > 0 & 
                                metab.mhu.bind$LOG2FC.2 > 0), ]
# decreased along disease severity
down.down <- metab.mhu.bind[which(metab.mhu.bind$LOG2FC < 0 & 
                                    metab.mhu.bind$LOG2FC.2 < 0), ]
metab.consis <- rbind(up.up, down.down)
# significant in SH or SM
metab.consis <- metab.consis[metab.consis$BH.1 < 0.05 | metab.consis$BH.2 < 0.05, ]
# add KEGG number
metab.consis$cpd <- m_c.all$KEGG[match(rownames(metab.consis), m_c.all$Metabolite)]

# KEGG enrichment analysis
x <- enricher(metab.consis$cpd, TERM2GENE = path_metab, TERM2NAME = path_intro,
              minGSSize = 2, pvalueCutoff = 0.1, pAdjustMethod = "BH")
kegg_table <- na.omit(as.data.frame(x))
kegg_table <- kegg_table[kegg_table$Count > 2, ]
# calculate fold enrichment
kegg_table$FoldEnrich <- apply(kegg_table, 1, 
                               function(x) as.numeric(unlist(strsplit(x[3], '[/]'))[1])/
                                 as.numeric(unlist(strsplit(x[3], '[/]'))[2])*
                                 as.numeric(unlist(strsplit(x[4], '[/]'))[2])/
                                 as.numeric(unlist(strsplit(x[4], '[/]'))[1]))

# set pathway order
kegg_table$Description <- factor(kegg_table$Description, levels = rev(kegg_table$Description))

# plot enrichment reslut
ggplot(kegg_table, aes(Description, -log10(p.adjust), fill = FoldEnrich)) + 
  geom_bar(stat = "identity", width = 0.7, color = "grey10") + 
  coord_flip() +
  scale_fill_viridis_c(option = "A", begin = 0.3, end = 0.9, breaks = c(10,15,20))+
  theme_classic() +
  labs(fill = "Fold enrichment", x="", y="-Log10 (P value)") + 
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(plot.title = element_text(vjust = 6, hjust = 0.5,  size = 8, face = "plain"))  +
  theme(axis.title.x =element_text(vjust = -2, hjust = 0.55, size = 8, face = "plain"),
        axis.title.y = element_text(vjust = 5, size = 8, face="plain")) +
  theme(axis.text.x = element_text(size = 8, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 8, face = "plain", colour = "black")) +
  theme(legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8),
        legend.key.height = unit(0.3, 'cm'), 
        legend.key.width = unit(0.3, 'cm'))

# Prepare Data of COVID-19 Patients and Healthy Controls ---------------------------

# Cytokine data of COVID-19 patients and healthy controls
df.cyto <- read.xlsx(paste0(wd2, "细胞因子数据.R.xlsx"), # SuppData 4
                     sheetIndex = 1, startRow = 1, header = T, row.names = 1, 
                     as.data.frame = T, check.names = F)
df.cyto <- df.cyto[1:54, c(25:72)]

# Significantly Altered Cytokines between H/M/S ---------------------------

df.cyto$Degree <- rep(c("Severe", "Mild", "Healthy"), c(23, 14, 17))
# compare two groups
type1 <- "Severe"
type2 <- "Healthy"
df.test <- rbind(df.cyto[which(df.cyto$Degree == type1),], 
                 df.cyto[which(df.cyto$Degree == type2),])
matab_num <- ncol(df.test) - 1
test_results <- data.frame(Metabolites = colnames(df.test)[1:matab_num])
# Mann-Whitney U test
test_results$wilcox <- apply(df.test[, 1:matab_num], 2,
                             function(x) unlist(wilcox.test(as.numeric(x) ~ df.test$Degree, data = df.test)[3]))
# adjust p value using BH method
test_results$wilcox_BH <- p.adjust(test_results$wilcox, method = "BH")
# Log2 (FC)
cn <- paste("FC_", type1, "/", type2, sep = "")
test_results[cn] <- apply(df.test[,1:matab_num], 2, 
                          function(x) 
                            mean(as.numeric(x[which(df.test$Degree == type1)]))/
                            mean(as.numeric(x[which(df.test$Degree == type2)])))
test_results$LOG2FC <- log2(test_results[,4])

# Point Plot of Cytokines ---------------------------

cyto.mhu.bind <- data.frame(row.names = colnames(df.cyto))
for (i in c(1:2)) {
  cyto.mhu <- read.xlsx("Serum.cytokine.MannWhitneyU.xlsx", # Source data
                        sheetIndex = i, startRow = 1, as.data.frame = T, header = T, check.names = F)
  cyto.mhu$Group <- group.name[i]
  cyto.mhu.bind <- rbind(cyto.mhu.bind, 
                         cyto.mhu)
}
cyto.mhu.bind$type <- ifelse(cyto.mhu.bind$BH < 0.05, "Significant", "Not Significant")
cyto.mhu.bind$type <- factor(cyto.mhu.bind$type, levels = c("Significant", "Not Significant"))

cyto_all <- data.frame(row.names = unique(cyto.mhu.bind$cytokine), 
                       MH_log2FC = cyto.mhu.bind$LOG2FC[match(unique(cyto.mhu.bind$cytokine), 
                                                              cyto.mhu.bind$cytokine[cyto.mhu.bind$Group == "MvsH"])], 
                       SH_log2FC = cyto.mhu.bind$LOG2FC[match(unique(cyto.mhu.bind$cytokine), 
                                                              cyto.mhu.bind$cytokine[cyto.mhu.bind$Group == "MvsH"])])
# sort the cytokines
cyto.ord <- hclust(dist(cyto_all, method = "euclidean"), method = "average")
cyto.ord <- rownames(cyto_all)[cyto.ord$order]
cyto.mhu.bind$cytokine <- factor(cyto.mhu.bind$cytokine, levels = rev(cyto.ord))

# point plot
ggplot(cyto.mhu.bind, aes(Group, cytokine)) +
  geom_point(aes(fill = LOG2FC, size = -log10(BH), color = type, shape = type)) + 
  scale_fill_gradientn(colours = colorRampPalette(c("#3d67a3", "white", "#ce1020"))(99)[-c(47:49, 50, 51:53)], 
                       limits = c(-6,6), breaks = c(-6, -3, 0, 3, 6)) +
  scale_shape_manual(values = c(21, 23)) +
  scale_color_manual(values = c("black", "grey50")) +
  scale_size_continuous(range = c(1, 3.5), breaks = c(1, 3, 5, 7)) +
  theme_dendro() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  coord_fixed(ratio = 0.7) +
  labs(fill = "Log2 (Fold change)", size = "-Log10 (P value)", shape = "", color = "") +
  theme(axis.text.x = element_text(size = 8, face = "plain", colour = "black", angle = 45), 
        axis.text.y = element_text(size = 8, face = "plain", colour = "black")) +
  theme(legend.text = element_text(size = 8, face = "plain", colour = "black"), 
        legend.title = element_text(size = 8, face = "plain", colour = "black"), 
        legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm"))   
        
