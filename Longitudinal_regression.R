library(xlsx)
library(Mfuzz)
library(ggplot2)
library(ggdendro)
library(clusterProfiler)
library(mgcv)
library(circlize)
library(RColorBrewer)
options(stringsAsFactors = F)

# Prepare Data ---------------------------

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

# Targeted metabolomics data
df.tar.metab <- read.csv("COVID-19.targeted.QCmad-TIC.csv", # SuppData 3
                         sep = ",", header = T, row.names = 1, check.names = F)
df.tar.metab <- df.tar.metab[, -grep("arnitine", colnames(df.tar.metab))]
cpd.tar <- m_c.tar$KEGG[match(colnames(df.tar.metab), m_c.tar$Metabolite)]

# Untargeted metabolomics data
df.untar.metab <- read.csv("COVID-19.untargeted.QCmad-TIC.csv",  # SuppData 3
                           header = T, sep = ",", row.names = 1, check.names = F)
cpd.untar <- m_c.untar$KEGG[match(colnames(df.untar.metab), m_c.untar$Metabolite)]
df.untar.metab <- df.untar.metab[, -na.omit(match(cpd.tar[which(cpd.tar != "-")], cpd.untar))]
cpd.untar <- cpd.untar[-na.omit(match(cpd.tar[which(cpd.tar != "-")], cpd.untar))]

# combine targeted and untargeted metabolomics data
df.metab <- cbind(df.tar.metab, df.untar.metab)
# write.xlsx(df.metab, paste0(wd2, "COVID-19.metab_all.QCmad-TIC.xlsx"),
#            sheetName = "all", row.names = T, col.names = T, append = T)

# Cytokine data from healthy controls (hc)
df.cyto.hc <- read.xlsx(paste0(wd2, "Cytokine.R.xlsx"), # SuppData 4
                        sheetIndex = 1, startRow = 1, header = T, row.names = 1, 
                        as.data.frame = T, check.names = F)
df.cyto.hc <- df.cyto.hc[38:54, 23:72]

# Cytokine data from follow-up (fu) patients 
df.cyto.fu <- read.xlsx(paste0(wd2, "Cytokine.R.xlsx"), # SuppData 3
                          sheetIndex = 2, startRow = 1, header = T, row.names = 1, 
                        as.data.frame = T, check.names = F)
df.cyto.fu <- df.cyto.fu[rownames(df.cyto.fu) != "S053",]

# days after symptoms onset
days.so <- data.frame("day" = df.cyto.fu$SampleDay, row.names = rownames(df.cyto.fu)) # SuppData 1
days.so$point <- (days.so$day - 1) %/% 3
days.so$point[which(days.so$point == 0 | days.so$point > 11)] <- NA
df.cyto.fu <- df.cyto.fu[!is.na(days.so$point), 27:76]
# df.cyto.fu <- df.cyto.fu[, 29:76]

# combine cytokine data from healthy controls and follow-up patients
df.cyto <- rbind(df.cyto.hc, df.cyto.fu)

# Mean Levels of Cytokines and Metabolites from HC and FU-patients ---------------------------

# Calculate mean levels of cytokines in each time point
df.ave.cyto <- data.frame(row.names = colnames(df.cyto.fu))
# use cytokine levels from hc as baseline
df.ave.cyto$H <- apply(df.cyto.hc, 2, function(x) mean(x))
for (i in c(1:11)) {
  tp <- paste0("p_", sprintf("%02d", i))
  sam <- rownames(days.so)[which(days.so$point == i)]
  df.ave.cyto[tp] <- apply(df.cyto.fu[sam,], 2, function(x) mean(x))
}

# Calculate mean levels of metabolites in each time point
df.ave.metab <- data.frame(row.names = colnames(df.metab))
# use metabolite levels from hc as baseline
df.ave.metab$H <- apply(df.metab[96:112,], 2, function(x) mean(x))
for (i in c(1:11)) {
  tp <- paste0("p_", sprintf("%02d", i))
  sam <- rownames(days.so)[which(days.so$point == i)]
  df.ave.metab[tp] <- apply(df.metab[sam,], 2, function(x) mean(x))
}

# combine mean cytokine and metabolite data
df.ave <- rbind(df.ave.cyto, df.ave.metab)

# Levels of Cytokines and Metabolites from HC and FU-patients ---------------------------

# Levels of all analytes in healthy controls and follow-up patients
df.analyte <- cbind(df.cyto, df.metab[match(rownames(df.cyto), rownames(df.metab)),])
df.analyte$point <- paste0("p_", 
                           sprintf("%02d", 
                                   days.so$point[match(rownames(df.analyte), rownames(days.so))]))
df.analyte$point[df.analyte$point == "p_NA"] <- "H"

# Significantly altered analytes in follow-up patients
# H point_01 point_03 point_05 point_07 point_09 point_11 
analyst.df <- data.frame()
for (i in c(c(1:11))) {
  tp <- paste0("p_", sprintf("%02d", i))
  if(T){
    df.test <- rbind(df.analyte[which(df.analyte$point == tp),], 
                     df.analyte[which(df.analyte$point == "H"),])
    matab_num <- ncol(df.test)-1
    test_results <- data.frame(Metabolites = colnames(df.test)[3:matab_num])
    # Mann-Whitney U test
    test_results$wilcox <- apply(df.test[, 3:matab_num], 2, 
                                 function(x) unlist(wilcox.test(x ~ df.test$point, data = df.test)[3]))
    # adjust p value using BH method
    test_results$wil_BH <- p.adjust(test_results$wilcox, method = "BH") 
    # Log2 (FC)
    test_results["FC"] <- apply(df.test[, 3:matab_num], 2, 
                                function(x) 
                                  mean(as.numeric(x[which(df.test$point == tp)]))/
                                  mean(as.numeric(x[which(df.test$point == "H")])))
    test_results$LOG2FC <- log2(test_results[,4])
    # add KEGG number to metbaolite
    test_results$cpd <- m_c.all$KEGG[match(test_results$Metabolites, m_c.all$Metabolite)]
    # time point information
    test_results$point <- tp
  }
  
  analyst.df <- rbind(analyst.df, test_results)
}

# analyst.union <- analyst.df[analyst.df$wil_BH < 0.1,]
analyst <- unique(analyst.df$Metabolites[analyst.df$wil_BH < 0.1])

# Mfuzz cluster ---------------------------

df.ave <- df.ave[rownames(df.ave) %in% analyst,]

# transform data frame to ExpressionSet object
df.ave.expr <- as.matrix(df.ave)
df.ave.expr <- ExpressionSet(df.ave.expr)
df.ave.expr <- standardise(df.ave.expr)

# Dmin(df.ave.expr, m = 1.43, crange = seq(3,10), repeats = 10, visu = TRUE)
# cselection(df.ave.expr, m=1.25, crange=seq(2,10), repeats=5, visu=TRUE)

# mufzz cluster
cl <- mfuzz(df.ave.expr, c = 4, m = 1.43) 
mfuzz.plot2(df.ave.expr, cl = cl, mfrow = c(2,2), x11 = F, ax.col = "black", cex = 3,
            xlab = "Days after symptoms onset", ylab = "Relative abundance", 
            centre = T, centre.col = "grey20", centre.lwd = 1.6,
            colo = colorRampPalette(c("yellow", "#ce1020"))(100), min.mem = 0)
# cl$size
# cl$cluster[cl$cluster == 2]

ms <- data.frame(cl$membership)
# analytes in each cluster and save the results
clust.analyte <- data.frame()
for (i in c(1:4)) {
  clust.analyte <- rbind(clust.analyte, 
                         data.frame(Cluster = cl$cluster[cl$cluster == i]))
}
ms$Cluster <- clust.analyte$Cluster[match(rownames(ms), rownames(clust.analyte))]

# Heatmap of Cytokines and Metabolites in Each Cluster  ---------------------------

# choose plot cytokines or metabolites
cm <- df.ave.metab
# choose cluster number
clust_num <- 3

if (T) {
  clust.cm <- cm[rownames(cm) %in% rownames(clust.analyte)[clust.analyte$Cluster == clust_num],]
  clust.cm <- apply(clust.cm, 2, function(x) log2(x/clust.cm$H))
  clust.cm.plot <- data.frame()
  for (i in c(1:12)) {
    clust.cm.plot <- rbind(clust.cm.plot, 
                           data.frame("Analyte" = rownames(clust.cm), 
                                      "Level" = clust.cm[,i], 
                                      "Point" = colnames(clust.cm)[i]))
  }
  # hclust and set the order of cytokines
  analyte.ord <- hclust(dist(clust.cm, method = "euclidean"), method = "average")  # plot(analyte.ord)
  analyte.ord <- analyte.ord$order
  clust.cm.plot$Analyte <- factor(clust.cm.plot$Analyte, 
                                  level = rev(rownames(clust.cm)[analyte.ord]))
  
  # plot heatmap
  ggplot(clust.cm.plot, aes(Point, Analyte)) +
    geom_tile(aes(fill = Level), color = "black", size = 0.4, height = 1, width = 1) + 
    scale_fill_gradientn(colours = c(colorRampPalette(c("#3d67a3", "white","#ce1020"))(100)), 
                         limits = c(-7, 7)) +
    theme_dendro() + 
    theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
    coord_fixed(ratio = 1) +
    labs(x = "", y = "", title = "", fill = "Log2 (Fold change)") +
    theme(axis.title.x = element_text(size = 6, face="plain", colour = "black"),
          axis.title.y = element_text(size = 6, face="plain", colour = "black"), 
          axis.text.x = element_text(size = 6, face = "plain", colour = "black", angle = 45), 
          axis.text.y = element_text(size = 6, face = "plain", colour = "black")) +
    theme(legend.title = element_text(size = 6, face = "plain", colour = "black"), 
          legend.text = element_text(size = 6, face = "plain", colour = "black"), 
          legend.key.height = unit(0.3, 'cm'), 
          legend.key.width = unit(0.3, 'cm'))
  
}

# KEGG Enrichment Analyses for Metabolites in Each Cluster ---------------------------

clust.kegg <- data.frame()
for (clust_num in c(1:4)) {
  clust.metab <- df.ave.metab[rownames(df.ave.metab) %in% 
                                rownames(clust.analyte)[clust.analyte$Cluster == clust_num],] 
  clust.metab$cpd <- m_c.all$KEGG[match(rownames(clust.metab), m_c.all$Metabolite)]
  
  # pathways with FDR < 0.05 & Counts > 2
  x <- enricher(clust.metab$cpd, TERM2GENE = path_metab, TERM2NAME = path_intro,
                minGSSize = 2, pvalueCutoff = 0.05, pAdjustMethod = "fdr")
  kegg_table <- na.omit(as.data.frame(x))
  kegg_table <- kegg_table[which(kegg_table$Count > 2),]
  # calculate fold enrichment
  kegg_table$FoldEnrich <- apply(kegg_table, 1, 
                                 function(x) 
                                   as.numeric(unlist(strsplit(x[3], '[/]'))[1])/
                                   as.numeric(unlist(strsplit(x[3], '[/]'))[2])*
                                   as.numeric(unlist(strsplit(x[4], '[/]'))[2])/
                                   as.numeric(unlist(strsplit(x[4], '[/]'))[1]))
  
  kegg_table$Cluster <- clust_num
  clust.kegg <- rbind(clust.kegg, kegg_table)
}

# point plot of enriched pathways
ggplot(clust.kegg, aes(Cluster, Description)) +
  geom_point(aes(size = FoldEnrich, fill = -log10(p.adjust)), color = "black", shape = 21) +
  scale_size(range = c(1, 3), breaks = c(10, 20, 30)) +
  scale_fill_viridis_c(begin = 0.2, option = "D", breaks = c(4, 6, 8, 10)) +
  theme_bw() + 
  theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 3, 
        panel.grid = element_blank()) +
  labs(x = "", y = "", title = "", size = "Fold enrichment", fill = "-Log10 (P value)") +
  theme(axis.title.x = element_text(size = 6, face="plain", colour = "black"),
        axis.title.y = element_text(size = 6, face="plain", colour = "black"), 
        axis.text.x = element_text(size = 6, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 6, face = "plain", colour = "black")) +
  theme(legend.title = element_text(size = 6, face = "plain", colour = "black"), 
        legend.text = element_text(size = 6, face = "plain", colour = "black"), 
        legend.key.height = unit(0.3, 'cm'), 
        legend.key.width = unit(0.3, 'cm'))

# Generalized Addictive Model (GAM) ---------------------------

df.gam <- data.frame(df.analyte, check.names = T)
metab_name <- data.frame(Check = colnames(df.gam), 
                         Name = colnames(df.analyte))
# patients's information
info <- read.xlsx("细胞因子数据.R.xlsx", # SuppData 4
                  sheetIndex = 1, startRow = 1, header = T, row.names = 1, as.data.frame = T, check.names = F)
df.gam$Age <- as.numeric(info$Age[match(rownames(df.gam), rownames(info))])
df.gam$Gender <- info$Gender[match(rownames(df.gam), rownames(info))]
df.gam$point <- days.so$point[match(rownames(df.gam), rownames(days.so))]
df.gam$point[is.na(df.gam$point)] <- 0

# log10-transformed cytokine levels
df.gam[, 3:50] <- log10(df.gam[, 3:50])

gamRes.df <- data.frame()
# cytokines 3:50; metabolites 51:303
for (m in colnames(df.gam)[3:303]) {
  gamRes <- gam(df.gam[[m]] ~ Gender + Age + s(point), data = df.gam)
  gamResum <- summary(gamRes)
  # save statistical results
  gamRes.df <- rbind(gamRes.df, 
                     data.frame(t(as.matrix(c(gamResum$p.table[1,], gamResum$p.table[2,], 
                                  gamResum$p.table[3,], gamResum$s.table[1,], 
                                  gamResum$r.sq))), check.names = F))
  
  # create PDF file
  pdf(paste0("cyto.time/", m, ".pdf"), onefile = F, width = 2.41, height = 2.41)
  # plot GAM result
  plot <- ggplot(data = gamRes$model, aes_string(x = names(gamRes$model)[4], y = names(gamRes$model)[1])) +
    stat_smooth(method = "gam", formula = y ~ s(x), se = T, level = 0.95, show.legend=T,
                color = "grey20", fill = "grey85", lwd = .3) +
    geom_point(color = "black", size = 0.8) +
    stat_summary(fun = "mean", geom = "line", colour = "#3490de", lwd = .3) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = .3, color = "#3490de", lwd = .3) +
    geom_hline(yintercept = mean(gamRes$model[which(gamRes$model[,4] == 0), 1]),
               lty = 2, col = "#3490de", lwd = 0.3) +
    theme_classic() +
    theme(aspect.ratio = 0.7) +
    labs(x = "Days after symptoms onset", y = paste(m, "pg/mL", sep = " "),
         title = paste("R2 = ", signif(summary(gamRes)$r.sq, 2),
                       " P =", signif(summary(gamRes)$s.table[[4]], 2))) +
    theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 0.67) +
    theme(title = element_text(size = 7, vjust = 1, hjust = 1)) +
    theme(axis.title.x = element_text(size = 8, hjust = 0.5, face = "plain", colour = "black"),
          axis.title.y = element_text(size = 8, hjust = 0.5, face = "plain", colour = "black"),
          axis.text.y = element_text(size = 8, face = "plain", colour = "black"),
          axis.text.x = element_text(size = 8, face = "plain", colour = "black"))
  # save and exit
  print(plot)
  dev.off()
}
rownames(gamRes.df) <- metab_name$Name[match(colnames(df.gam)[3:303], metab_name$Check)]

# Cytokine-Metabolite Regression Analyses ---------------------------

df.gam <- data.frame(df.analyte[18:69, 3:303], check.names = T)
colnames(df.gam) <- make.names(colnames(df.gam))
df.gam <- log10(df.gam)

# patients's information
info <- read.xlsx("Cytokine.R.xlsx", sheetIndex = 1, startRow = 1, header = T, 
                  row.names = 1, as.data.frame = T, check.names = F)
df.gam$Age <- as.numeric(info$Age[match(rownames(df.gam), rownames(info))])
df.gam$Gender <- info$Gender[match(rownames(df.gam), rownames(info))]
df.gam$point <- days.so$point[match(rownames(df.gam), rownames(days.so))]

# use samples with time point from 1 to 7
df.gam <- df.gam[which(df.gam$point <= 7 & !is.na(df.gam$point)), ]

# linear regression
lmRes.df <- data.frame()
for (m in colnames(df.gam)[1:48]) {
  for (i in colnames(df.gam)[49:301]) {
    fmla <- as.formula(paste(m," ~ ", paste(c("Gender", "Age", i), collapse= " + ")))
    lmRes <- lm(formula = fmla, data = df.gam)
    lmRes <- data.frame(coef(summary(lmRes)))
    lmRes.df <- rbind(lmRes.df, data.frame(Estimate = lmRes$Estimate[4], 
                                           Std.err = lmRes$Std..Error[4], 
                                           t.value =  lmRes$t.value[4], 
                                           p = lmRes$Pr...t..[4], 
                                           BH = p.adjust(lmRes$Pr...t.., method = "BH")[4],
                                           Cyto = m, 
                                           Metab = i))
  }
}
lmRes.df$Metab <- colnames(df.metab)[match(lmRes.df$Metab, make.names(colnames(df.metab)))]
lmRes.df$Cyto <- colnames(df.cyto)[match(lmRes.df$Cyto, make.names(colnames(df.cyto)))]

# Use Circlize to Show Correlations between Cytokines and Metabolites ---------------------------

lmRes.df.n <- lmRes.df[lmRes.df$BH < 0.05,]
# add KEGG number
lmRes.df.n$cpd <- m_c.all$KEGG[match(lmRes.df.n$Metab, m_c.all$Metabolite)]

# aim pathways
path <- c("Arginine biosynthesis", "Arginine and proline metabolism",  
          "Nicotinate and nicotinamide metabolism", "Tryptophan metabolism", 
          "Purine metabolism", "Pyrimidine metabolism", 
          "Citrate cycle (TCA cycle)", "Primary bile acid biosynthesis", 
          "Valine, leucine and isoleucine biosynthesis", "Cysteine and methionine metabolism")
path.cpd <- path_intro[match(path, path_intro$X2),]
path.cpd <- path_metab[path_metab$V1 %in% path.cpd$X1,]
for (i in c(1:nrow(path.cpd))) {
  path.cpd$path[i] <- path_intro$X2[which(path_intro$X1 == path.cpd$V1[i])]
}  

# aim cytokines
cyto <- c("IL-6", "IP-10",  "M-CSF")
# cyto <- c("MCP-3", "MIP-1b",  "MIP-1a", "TNF-a", "G-CSF")
# cyto <- c("IL-18", "IL-1a",  "IL-1b", "IFN-g", "IL-10", "MCP-1", "IFN-a2")

# simplify correlation data frame
corrRes.df <- data.frame()
for (c in cyto) {
  corrRes <- lmRes.df.n[which(lmRes.df.n$Cyto == c & lmRes.df.n$cpd != "-"), ]
  corrRes.df <- rbind(corrRes.df, corrRes)
}
corrRes.df <- corrRes.df[corrRes.df$cpd %in% path.cpd$V2,]
corrRes.df$path <- path.cpd$path[match(corrRes.df$cpd, path.cpd$V2)]

ex <- c("1-Methyl-4-pyridone-3-carboxamide", "1-Methyl-5-carboxylamide-2-pyridone", 
        "5-Hydroxyindoleacetaldehyde", "3-Methyl pyruvic acid",
         "2-Isopropylmaleate", "Iminoaspartic acid", 
        "Glycinamide ribonucleotide", "1-Methylpyrrolinium")
corrRes.df <- corrRes.df[!corrRes.df$Metab %in% ex,]
corrRes.df$path[corrRes.df$path %in% 
                  c("Arginine biosynthesis ","Arginine and proline metabolism ")] <- "Arginine metabolism "
# nodes
xlim <- data.frame()
for (i in unique(corrRes.df$Cyto)) {
  xlim <- rbind(xlim, 
                data.frame(name = i, start = 0, 
                           end = length(which(corrRes.df$Cyto == i))))
}
for (i in unique(corrRes.df$Metab)) {
  xlim <- rbind(xlim, 
                data.frame(name = i, start = 0, 
                           end = length(which(corrRes.df$Metab == i))))
}

# links
links <- corrRes.df[, c(6:7, 9, 3)]
links$cytoE <- rep(NA, nrow(links))
links$metabE <- rep(NA, nrow(links))
for (i in unique(links$Cyto)) {
  links$cytoE[which(links$Cyto == i)] <- c(1:nrow(links[which(links$Cyto == i), ]))
}
links$cytoS <- links$cytoE - 1
for (i in unique(links$Metab)) {
  links$metabE[which(links$Metab == i)] <- c(1:nrow(links[which(links$Metab == i), ]))
}
links$metabS <- links$metabE - 1
# direction color
links$col <- ifelse(links$t.value > 0, alpha("#eae222", .6), alpha("#25667b", .6)) 
# pathway color
links$pcol <- factor(links$path, labels = brewer.pal(9, "Paired"))
# analyte color
n <- length(cyto)
r <- nrow(xlim)
xlim$col[(n + 1):r] <- as.character(links$pcol[match(xlim$name[(n + 1):r], links$Metab)])
xlim[(n + 1):r,] <- xlim[order(xlim$col[(n + 1):r])+(n),]

# circos plot
# par(mfrow = c(1,1))
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
  circos.link(links[n, 1], links[n, c(7, 5)], links[n, 2], links[n, c(8, 6)], 
              col = alpha(links[n, 9], .6), border = NA, h = 0.7) #intreval to interval
}

