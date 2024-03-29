library(xlsx)
library(clusterProfiler)
library(ggplot2)
library(ggpubr)
library(ggdendro)
library(ggthemes)
library(igraph)
library(ggraph)
library(circlize)
options(stringsAsFactors = F)

# Prepare Metabolomica and Cytokine Data of COVID-19 Patients and Healthy Controls ---------------------------

# Reference files
m_c.tar <- read.xlsx("metabolite.HMDB.KEGG.HULAB.xlsx", # not provided in Supp files
                     sheetIndex = 2, startRow = 1, header = T, as.data.frame = T)
m_c.untar <- read.xlsx("metabolite.HMDB.KEGG.HULAB.xlsx", # not provided in Supp files
                       sheetIndex = 3, startRow = 1, header = T, as.data.frame = T)
m_c.all <- rbind(data.frame("Metabolite" = m_c.tar$Metabolite, "KEGG" = m_c.tar$KEGG),
                 data.frame("Metabolite" = m_c.untar$Metabolite, "KEGG" = m_c.untar$KEGG))
# metabolits in each pathways
path_metab <- read.table('path.cpd.txt', header = F, sep = '\t') # not provided in Supp files
path_metab$V2 <- gsub("cpd:", "", path_metab$V2)
# metabolic pathways
path_intro <- read.xlsx("has_MetabolicPathway.xlsx", # not provided in Supp files
                        sheetIndex = 2, startRow = 1, as.data.frame = T, header = F, check.names = T)

# Targeted metabolomics data of COVID-19 patients and healthy controls
df.tar.metab <- read.csv("COVID-19.targeted.QCmad-TIC.csv",  # SuppData 3
                         sep = ",", header = T, row.names = 1, check.names = F)
df.tar.metab <- df.tar.metab[, -grep("arnitine", colnames(df.tar.metab))]
cpd.tar <- m_c.tar$KEGG[match(colnames(df.tar.metab), m_c.tar$Metabolite)]

# Untargeted metabolomics data of COVID-19 patients and healthy controls
df.untar.metab <- read.csv("COVID-19.untargeted.QCmad-TIC.csv", # SuppData 3
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

# Cytokine data of COVID-19 patients and healthy controls
df.cyto <- read.xlsx("Cytokine.R.xlsx", # SuppData 4
                     sheetIndex = 1, startRow = 1, header = T, row.names = 1, 
                     as.data.frame = T, check.names = F)
df.cyto <- df.cyto[1:54, c(25:72)]

# clinical information
df.cli <- read.xlsx("Cytokine.R.xlsx", # SuppData 1
                    sheetIndex = 1, startRow = 1, header = T, row.names = 1, 
                    as.data.frame = T, check.names = F)
df.cli <- df.cli[1:54, 4:5]
df.cli$Age <- as.numeric(df.cli$Age)

# conmine clinical, cytokine, metabolite data
df.regress <- data.frame(df.cli, df.cyto[match(rownames(df.cli), rownames(df.cyto)), ], 
                         df.metab[match(rownames(df.cli), rownames(df.metab)), ])

# Linear Regression Analysis of Cytokines and Metabolites ---------------------------

# function for regression analysis
CytoPath <- function(x){
  cyto.metab.corr <- data.frame() 
  lm.metab.list <- data.frame() 
  
  mc <- make.names(colnames(df.metab))
  for (m in c(1:length(mc))) {
    fmla <- as.formula(paste(x," ~ ", paste(c("Gender", "Age", mc[m]), collapse = " + ")))
    rg.result <- summary(lm(formula = fmla, df.regress.part))
    rg.coef <- as.data.frame(rg.result$coefficients)
    rg.coef$BH <- p.adjust(rg.coef$`Pr(>|t|)`, method = "BH")
    if (rg.coef$BH[4] < 0.1) {
      lm.metab.list <- rbind(lm.metab.list, rg.coef[4,])
    }
    # lm.metab.list <- rbind(lm.metab.list, rg.coef[4,])
  }
  lm.metab.list$Cyto <- rep(x, nrow(lm.metab.list))
  lm.metab.list$Metab <- rownames(lm.metab.list)
  lm.metab.list$cpd <- m_c.all$KEGG[match(lm.metab.list$Metab, make.names(m_c.all$Metabolite))]
  cyto.metab.corr <- rbind(cyto.metab.corr, lm.metab.list)
  
  # KEGG enrichment analysis
  enrichRes <- enricher(lm.metab.list$cpd[lm.metab.list$cpd != "-" & lm.metab.list$BH < 0.1],
                        TERM2GENE = path_metab, TERM2NAME = path_intro,
                        minGSSize = 2, pvalueCutoff = 0.1, pAdjustMethod = "BH")
  kegg_table <- na.omit(as.data.frame(enrichRes))

  if(nrow(kegg_table) > 0){
    # calculate absolute mean log2FC
    path_aveT <- function(y){
      metabs <- unlist(strsplit(y[8], "[/]"))
      mean(abs(lm.metab.list[match(metabs, lm.metab.list$cpd), 3]))
    }
    kegg_table$t_mean <- apply(kegg_table, 1, path_aveT)
    # calculate fold enrichment
    kegg_table$FoldEnrich <- apply(kegg_table, 1,
                                   function(x) as.numeric(unlist(strsplit(x[3], '[/]'))[1])/
                                     as.numeric(unlist(strsplit(x[3], '[/]'))[2])*
                                     as.numeric(unlist(strsplit(x[4], '[/]'))[2])/
                                     as.numeric(unlist(strsplit(x[4], '[/]'))[1]))
    kegg_table$cyto <- x
  }
  # save results
  lab <- list(kegg = kegg_table, corr = cyto.metab.corr)
  return(lab)
}

group.name <- c("S", "M", "H")
samp.ind <- list(S = c(1:23), M = c(24:37), H = c(38:54))
# cytokine-metabolite
cyto_path.all <- data.frame()
cyto_metab.c.all <- data.frame()
for (s in c(1:3)) {
  df.regress.part <- df.regress[samp.ind[[group.name[s]]], ]
  # log10 transform
  df.regress.part[3:303] <- log10(df.regress.part[3:303])
  # linear regression analysis
  cyto_path.d <- data.frame()
  cyto_metab.corr <- data.frame()
  for (i in c(1:ncol(df.cyto))) {
    x <- colnames(data.frame(df.cyto)[i])
    cyto_path.d <- rbind(cyto_path.d, CytoPath(x)$kegg)
    cyto_metab.corr <- rbind(cyto_metab.corr, CytoPath(x)$corr)
  }
  cyto_path.d$Group <- group.name[s]
  cyto_metab.corr$Group <- group.name[s]
  
  # combine all sample types
  cyto_path.all <- rbind(cyto_path.all, cyto_path.d)
  cyto_metab.c.all <- rbind(cyto_metab.c.all, cyto_metab.corr)
}
cyto_metab.c.all$Cyto <- colnames(df.cyto)[match(cyto_metab.c.all$Cyto, make.names(colnames(df.cyto)))]
# cyto_path.all$cyto <- colnames(df.cyto)[match(cyto_path.all$cyto, make.names(colnames(df.cyto)))]
cyto_metab.c.all$Metab <- colnames(df.metab)[match(cyto_metab.c.all$Metab, make.names(colnames(df.metab)))]

# Show the Correlations between Cytokines and Metabolic Pathways ---------------------------

# significantly altered cytokines in severe patients compared with healthy controls
diff.cyto <- read.xlsx('Serum.cytokine.MannWhitneyU.xlsx', # Source data
                       sheetIndex = 2, startRow = 1, header = T, as.data.frame = T, check.names = F)
diff.cyto <- diff.cyto[diff.cyto$BH < 0.05, ]

# range(cyto_path.all$t_mean)
cyto_path.df <- cyto_path.all[cyto_path.all$cyto %in% make.names(diff.cyto$cytokine) & 
                               cyto_path.all$Group == "S", ]
cyto_path.df <- cyto_path.df[cyto_path.df$Count > 1, ]

if (T) {
  path_corr <- data.frame(row.names = unique(cyto_path.df$Description))
  for (i in unique(cyto_path.df$cyto)) {
    name <- cyto_path.df[which(cyto_path.df$cyto == i),]
    path_corr[i] <- name$t_mean[match(rownames(path_corr), name$Description)]
  }
  # Na replace
  NaReplace <- function(x){
    sub <- which(is.na(x))
    x[sub] <- 0
    x
  }
  path_corr.na <- data.frame(apply(path_corr, 2, NaReplace))
  # sort pathways
  path.ord <- hclust(dist(path_corr.na, method = "euclidean"), 
                     method = "average")
  path.ord <- rownames(path_corr.na)[path.ord$order]
  cyto_path.df$Description <- factor(cyto_path.df$Description, levels = rev(path.ord))
  # sort cytokines
  cyto.ord <- hclust(dist(t(path_corr.na), method = "euclidean"), 
                     method = "complete")
  cyto.ord <- colnames(path_corr.na)[cyto.ord$order]
  cyto_path.df$cyto <- factor(cyto_path.df$cyto, levels = cyto.ord)
}

# point plot of correlations between
range(cyto_path.df$t_mean)
ggplot(cyto_path.df, aes(cyto, Description)) +
  geom_point(aes(fill = t_mean, size = -log10(p.adjust)), color = "black", shape = 21) +
  scale_size(range = c(1, 3), breaks = c(2, 4, 6, 8)) +  # c(2, 4, 6, 8) c(2, 4, 6)
  scale_fill_viridis_c(direction = -1, begin = 0.2, option = "A", 
                       breaks = c(2, 3, 4, 5, 6)) +  # c(2.5, 3.0, 3.5, 4.0) c(2, 3, 4, 5)
  theme_dendro() + 
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  labs(fill = "abs. T statistics", size = "-Log10 (P value)") +
  theme(axis.text.x = element_text(size = 6, face = "plain", colour = "black", 
                                   angle = 45, vjust = 1, hjust = 1), 
        axis.text.y = element_text(size = 6, face = "plain", colour = "black")) +
  theme(legend.text = element_text(size = 6, face = "plain", colour = "black"), 
        legend.title = element_text(size = 6, face = "plain", colour = "black"),
        legend.spacing.y = unit(0.3, 'cm'), 
        legend.key.size = unit(0.3, 'cm'))

# Correlation Network by igraph ---------------------------

aim.cyto <- c("IL-6", "IL-1a", "IL-1b", "IL-18", "IP-10", 
              "IFN-g", "GM-CSF", "M-CSF", "G-CSF", "IL-12(p40)", "IL-10", 
              "MIP-1b", "MIP-1a", "MCP-3")

cyto_metab.c.p <- cyto_metab.c.all[cyto_metab.c.all$Group == "S" & cyto_metab.c.all$Cyto %in% aim.cyto, ]
cyto_metab.c.p <- cyto_metab.c.p[-c(grep("PC", cyto_metab.c.p$Metab), grep("LPE", cyto_metab.c.p$Metab), 
                                    grep("FFA", cyto_metab.c.p$Metab), grep("arnitine", cyto_metab.c.p$Metab)), ]

if (T) {
  # nodes
  num.cyto <- unique(cyto_metab.c.p$Cyto)
  num.metab <- unique(cyto_metab.c.p$Metab)
  nodes <- data.frame(name = c(num.cyto, num.metab), 
                      cpd = c(rep("NA", length(num.cyto)), 
                              cyto_metab.c.p$cpd[match(num.metab, cyto_metab.c.p$Metab)]),
                      type = c(rep("cyto", length(num.cyto)), 
                               rep("metab", length(num.metab))))
  # links
  links <- data.frame(from = cyto_metab.c.p$Cyto, 
                      to = cyto_metab.c.p$Metab, 
                      cpd = cyto_metab.c.p$cpd,
                      weight = abs(cyto_metab.c.p$`t value`), 
                      direct = ifelse(cyto_metab.c.p$`t value` > 0, "pos", "neg"))
  # label specific pathways
  path <- c("Arginine biosynthesis", "Arginine and proline metabolism",  
            "Nicotinate and nicotinamide metabolism", "Tryptophan metabolism", 
            "Purine metabolism", "Pyrimidine metabolism",
            "Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)", 
            "Primary bile acid biosynthesis", "Cysteine and methionine metabolism", 
            "Valine, leucine and isoleucine biosynthesis", "Alanine, aspartate and glutamate metabolism", 
            "Glycine, serine and threonine metabolism", "Phenylalanine, tyrosine and tryptophan biosynthesis")
  
  path.cpd <- path_intro[match(path, path_intro$X2),]
  path.cpd <- path_metab[path_metab$V1 %in% path.cpd$X1,]
  for (i in c(1:nrow(path.cpd))) {
    path.cpd$path[i] <- path_intro$X2[which(path_intro$X1 == path.cpd$V1[i])]
  }
  nodes$path <- path.cpd$path[match(nodes$cpd, path.cpd$V2)]
  nodes$path[is.na(nodes$path)] <- nodes$type[is.na(nodes$path)]
  # combine the pathways
  nodes$path[which(nodes$path == "Arginine biosynthesis" | 
                     nodes$path == "Arginine and proline metabolism")] <- "Arginine metabolism"
  nodes$path[which(nodes$path == "Alanine, aspartate and glutamate metabolism" | 
                     nodes$path == "Glycine, serine and threonine metabolism" | 
                     nodes$path == "Phenylalanine, tyrosine and tryptophan biosynthesis")] <- "Amino acid metabolism"
}

# construct the correlation network
net <- graph_from_data_frame(links, vertices = nodes, directed = F)
# set color
V(net)$vertex.color <- c("#dddddd", "#c8d5b9")[1 + (V(net)$type == "metab")]
E(net)$color <- c("#25667b", "#eae222")[1 + (E(net)$direct == "pos")]

# cluster the network
clg <- cluster_fast_greedy(net, modularity = T, membership = T, weights = E(net)$weight)
# plot each subgraph
table(clg$membership)
sub <- V(net)[clg$membership == 4]
if (T) {
  g <- induced_subgraph(net, sub)
  V(g)$path
  
  color <- list("cyto" = "#dddddd", "metab" = "#c8d5b9", 
                "Arginine metabolism" = "#a6cee3", 
                "Citrate cycle (TCA cycle)" = "#1f78b4", "Glycolysis / Gluconeogenesis" = "#c05555"
                "Nicotinate and nicotinamide metabolism" = "#15983f", "Tryptophan metabolism" = "#ff7f00", 
                "Purine metabolism" = "#e31a1c", "Pyrimidine metabolism" = "#fdbf6f", 
                "Primary bile acid biosynthesis" = "#fb9a99",
                "Cysteine and methionine metabolism" = "#b2df8a", 
                "Valine, leucine and isoleucine biosynthesis" = "#cab2d6",
                "Amino acid metabolism" = "#557571", )
  dir_color <- list("pos" = "#eae222", "neg" = "#25667b")
  
  ggraph(g, layout = 'dh') + 
    geom_edge_fan(aes(color = direct), show.legend = FALSE) +  # edge_width = weight
    geom_node_point(aes(fill = path), shape = 21, size = 3, color = "grey40") + 
    geom_node_text(aes(label = name), size = 1.8) +
    scale_edge_color_manual(values = alpha(dir_color, 0.7)) + 
    scale_fill_manual(values = color) + 
    scale_edge_width(range = c(0.4, 1.1)) + 
    guides(size = F, fill = F) +  
    theme_graph() 
}

# Correlations between Cytokines and Metabolites in Specific Pathways ---------------------------

path <- c("Arginine biosynthesis", "Arginine and proline metabolism")
path <- c("Purine metabolism")
path <- c("Tryptophan metabolism","Nicotinate and nicotinamide metabolism")

path.cpd <- path_intro[match(path, path_intro$X2),]
path.cpd <- path_metab[path_metab$V1 %in% path.cpd$X1,]
for (i in c(1:nrow(path.cpd))) {
  path.cpd$path[i] <- path_intro$X2[which(path_intro$X1 == path.cpd$V1[i])]
}

corrRes.df <- cyto_metab.c.all[cyto_metab.c.all$Group == "S" & cyto_metab.c.all$cpd %in% path.cpd$V2,]
corrRes.df$path <- path.cpd$path[match(corrRes.df$cpd, path.cpd$V2)]
ex <- c("Glyoxylic acid", "Fumaric acid", "Guanidineacetic acid", "Pyruvate", "S-Adenosylmethionine")
ex <- c("Allantoic acid", "Glutamine", "Glycine", "N-Formiminoglycine")
ex <- c("1-Methyl-4-pyridone-3-carboxamide", "6-Hydroxypseudooxynicotine",  
        "5-Hydroxyindoleacetic acid", "1-Methyl-5-carboxylamide-2-pyridone", "Aspartic acid", "Fumaric acid", "Melatonin", "Pyruvate")
corrRes.df <- corrRes.df[!corrRes.df$Metab %in% ex,]

if (T) {
  # nodes
  xlim <- data.frame()
  for (i in unique(corrRes.df$Cyto)) {
    xlim <- rbind(xlim, data.frame(name = i, start = 0, end = length(which(corrRes.df$Cyto == i)), color = "#dddddd"))
  }
  for (i in unique(corrRes.df$Metab)) {
    xlim <- rbind(xlim, data.frame(name = i, start = 0, end = length(which(corrRes.df$Metab == i)), color = "#94c9e4"))
  }
  
  # links
  links <- corrRes.df[, c(6:8,3)]
  links$cytoS <- rep(NA, nrow(links))
  links$cytoE <- rep(NA, nrow(links))
  links$metabS <- rep(NA, nrow(links))
  links$metabE <- rep(NA, nrow(links))
  for (i in unique(links$Cyto)) {
    links$cytoE[which(links$Cyto == i)] <- c(1:nrow(links[which(links$Cyto == i),]))
  }
  for (i in unique(links$Metab)) {
    links$metabE[which(links$Metab == i)] <- c(1:nrow(links[which(links$Metab == i),]))
  }
  links$cytoS <- links$cytoE - 1
  links$metabS <- links$metabE - 1
  links$color <- ifelse(links$`t value` > 0, alpha("#eae222",.8), alpha("#25667b", .8)) 
}

# circos plot
circos.clear()
circos.par(start.degree = 0, gap.degree = 1, 
           track.margin = c(0.01, 0.01), cell.padding = c(0.02, 0, 0.02, 0))
circos.genomicInitialize(xlim, plotType = "NULL")
circos.genomicLabels(xlim, labels = xlim$name, niceFacing, side = "outside",
                     padding = 0.1, connection_height = convert_height(3, "mm"), 
                     cex = 0.5, line_lwd = 0.6)
circos.genomicTrackPlotRegion(xlim, ylim = c(0,1), track.height = 0.05, bg.col = xlim$col, 
                              bg.border = "black", bg.lwd = 0.5)
for (n in c(1:nrow(links))) {
  circos.link(links[n, 1], links[n, 5:6], links[n, 2], links[n, 7:8], 
              col = alpha(links[n, 9], .6), border = NA, h = 0.7) #intreval to interval
}

# Specific Cytokines ---------------------------
cyto_metab.c.s <- cyto_metab.c.all[cyto_metab.c.all$Cyto == "IP-10", ]
cyto_metab.c.s <- cyto_metab.c.s[-c(grep("PC", cyto_metab.c.s$Metab), grep("LPE", cyto_metab.c.s$Metab), 
                                    grep("FFA", cyto_metab.c.s$Metab), grep("arnitine", cyto_metab.c.s$Metab)), ]

ggplot(cyto_metab.c.s, aes(Group, Metab)) + 
  geom_point(aes(fill = `t value`, size = -log10(BH)), shape = 21) +
  theme_dendro() +
  scale_size(range = c(1, 2.4), breaks = c(1, 2, 3)) + 
  scale_fill_gradientn(colors = colorRampPalette(c("#3d67a3", "white", "#ce1020"))(100), 
                       limits = c(-7, 7), breaks = c(-6, -3, 0, 3, 6)) + 
  theme_dendro() + 
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  labs(fill = "T statistics", size = "-Log10 (P value)") +
  theme(axis.text.x = element_text(size = 6, face = "plain", colour = "black", 
                                   angle = 45, vjust = 1, hjust = 1), 
        axis.text.y = element_text(size = 6, face = "plain", colour = "black")) +
  theme(legend.text = element_text(size = 6, face = "plain", colour = "black"), 
        legend.title = element_text(size = 6, face = "plain", colour = "black"),
        legend.spacing.y = unit(0.3, 'cm'), 
        legend.key.size = unit(0.3, 'cm'))

# Specific Metabolites ---------------------------

metab <- c("Arginine", "Aspartic acid", "Ornithine", "Citrulline", "Proline", "Urea", 
           "Tryptophan", "Kynurenine", "Kynurenic acid", "NMN", "Nicotinate", "Nicotinamide", "NAD", 
           "Hypoxanthine", "Xanthine", "Guanine", "Adenine", "Adenosine", 
           "Xanthosine", "XMP", "AMP", "GMP", "AICAR")
cyto <- c("IL-6", "IL-1a", "IL-1b", "TNF-a", "GM-CSF", "G-CSF", "M-CSF", "IFN-g", "IL-12(p40)", "IL-10", 
          "IP-10", "IFN-a2", "MIP-1b", "MCP-1", "MIP-1a", "MCP-3", "IL-18", "IL-1ra")

cyto_metab.c.s <- cyto_metab.c.all[cyto_metab.c.all$Metab %in% metab & cyto_metab.c.all$Group == "M", ]
cyto_metab.c.s <- cyto_metab.c.s[cyto_metab.c.s$Cyto %in% cyto, ]
cyto_metab.c.s$sig <- ifelse(cyto_metab.c.s$BH < 0.1, "sig", "not sig")
cyto_metab.c.s$sig <- factor(cyto_metab.c.s$sig, levels = c("sig", "not sig"))
cyto_metab.c.s$Metab <- factor(cyto_metab.c.s$Metab, levels = metab)
cyto_metab.c.s$Cyto <- factor(cyto_metab.c.s$Cyto, levels = rev(cyto))

ggplot(cyto_metab.c.s, aes(Metab, Cyto)) + 
  geom_point(aes(fill = `t value`, size = -log10(BH), color = sig), shape = 21) +
  scale_size(range = c(1, 2.8), breaks = c(1, 2, 3)) + 
  scale_color_manual(values = c("black", "grey80")) +
  scale_fill_gradientn(colors = colorRampPalette(c("#3d67a3", "white", "#ce1020"))(100), 
                       limits = c(-6, 6), breaks = c(-5, -2.5, 0, 2.5, 5)) + 
  theme_bw() +
  # theme_dendro() + 
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  labs(fill = "T statistics", size = "-Log10 (P value)") +
  theme(axis.text.x = element_text(size = 6, face = "plain", colour = "black", 
                                   angle = 45, vjust = 1, hjust = 1), 
        axis.text.y = element_text(size = 6, face = "plain", colour = "black")) +
  theme(legend.text = element_text(size = 6, face = "plain", colour = "black"), 
        legend.title = element_text(size = 6, face = "plain", colour = "black"),
        legend.spacing.y = unit(0.1, 'cm'), 
        legend.key.size = unit(0.3, 'cm'))
        
