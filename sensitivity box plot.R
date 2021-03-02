rm(list = ls())
library(tidyverse)
LFC <- read_csv("D:/Winnie/Depmap/raw data/primary-screen-replicate-collapsed-logfold-change.csv") %>% as.data.frame()
treatment.info <- read_csv("D:/Winnie/Depmap/raw data/primary-screen-replicate-collapsed-treatment-info.csv")%>% as.data.frame()

LFC <- LFC %>% column_to_rownames("X1")
colnames(LFC) <- treatment.info$broad_id[match(colnames(LFC), treatment.info$column_name)]

cell_line_info <- read_csv("D:/Winnie/Depmap/raw data/sample_info.csv") %>% as.data.frame()
cell_line_info <- select(cell_line_info, DepMap_ID, stripped_cell_line_name, primary_disease, primary_or_metastasis)

LFC$cancer_type <- cell_line_info$primary_disease[match(rownames(LFC), cell_line_info$DepMap_ID)]
LFC <- LFC %>% select(cancer_type, everything())

Drug_LFC <- LFC %>% select(cancer_type, 
                           treatment.info$broad_id[match("duloxetine", treatment.info$name)])
Drug_LFC <-  rename(Drug_LFC,"LFC.value" =`BRD-K71103788-003-12-0` )
Drug_LFC <- Drug_LFC %>% rownames_to_column("gene") %>% 
  filter(cancer_type != "Engineered" |!is.na(cancer_type)) %>% 
  column_to_rownames("gene")
Drug_LFC <- Drug_LFC %>% na.omit()
#Drug_AUC$PorM <- cell_line_info$primary_or_metastasis[match(Drug_AUC$depmap_id, cell_line_info$DepMap_ID)]
#Drug_AUC <- Drug_AUC %>% filter(cancer_type != "UNABLE TO CLASSIFY")
ggplot(iris, aes(x = reorder(Species, Sepal.Width, FUN = median), y = Sepal.Width)) + geom_boxplot()
ggplot(iris, aes(x = Species, y = Sepal.Width)) + geom_boxplot()

Drug_LFC %>% 
  ggplot(aes(x = reorder(cancer_type,LFC.value, median),
             y= LFC.value, fill = cancer_type)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(size=8, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_blank(),legend.position = "none",
        axis.title.y = element_text(size = 9)) +
  ylab("PRISM Primary Screening LFC")

cancer <- unique(Drug_LFC$cancer_type)
lineage.enrich <- data.frame(Cancer = cancer)
rownames(lineage.enrich) <- cancer
for (i in 1:length(cancer)){
TP <- 0
FP <- 0
FN <- 0
TN <- 0
cancer_name <- cancer[i]
cancer_LFC <- Drug_LFC %>% filter(cancer_type == cancer_name)
for (a in 1:length(cancer_LFC$cancer_type)) {
if (cancer_LFC$LFC.value[a] < 0.3) {TP <- TP+1} 
  else {FN <- FN +1}
}
non_cancer_LFC <- Drug_LFC %>% filter(cancer_type != cancer_name)
for (b in 1:length(non_cancer_LFC$cancer_type)) {
  if (non_cancer_LFC$LFC.value[b] < 0.3) {FP <- FP+1} 
  else {TN <- TN +1}
}
lineage.enrich[cancer_name,"adj.p"] <- p.adjust(fisher.test(matrix(c(TP, FP, FN, TN), nrow =2), alternative = "less")$p.value, "BH") 
}
#Biomarker
cancer_LFC <- Drug_LFC %>% rownames_to_column("gene") %>% 
  filter(cancer_type == "Colon/Colorectal Cancer") %>% 
  column_to_rownames("gene")
gene_signature <- read_csv("D:/Winnie/Depmap/raw data/CCLE_mutation_Winnie_full.csv")
gene_signature <- gene_signature %>% column_to_rownames("X1")
cancer_signature <- gene_signature %>% select(rownames(cancer_LFC))

biomarker <- data.frame(gene = rownames(cancer_signature), wilcox= NA)
rownames(biomarker) <- biomarker$gene
for (i in 1: length(biomarker$gene)){
  gene <- cancer_signature[i,] %>% t() %>% as.data.frame()
  colnames(gene) <- "status"
  WT <- rownames(gene %>% rownames_to_column("gene") %>% 
                   filter(status == "0") %>% 
                   column_to_rownames("gene"))
  Mut <- rownames(gene %>% rownames_to_column("gene") %>% 
                    filter(status == "1") %>% 
                    column_to_rownames("gene"))
  WT.lfc <- cancer_LFC$LFC.value[match(WT,rownames(cancer_LFC))] %>% na.omit()
  Mut.lfc <- cancer_LFC$LFC.value[match(Mut, rownames(cancer_LFC))] %>% na.omit()
  if(length(Mut.lfc) >= 1) {
    wilcox <- wilcox.test(WT.lfc, Mut.lfc, alternative = "greater")
    biomarker[i, "wilcox"] <- wilcox[["p.value"]]
  } else {
    biomarker[i, "wilcox"] <- NA
  }
  if (length(Mut.lfc) >= 2){
    ttest <- t.test(WT.lfc, Mut.lfc, alternative = "greater")
    biomarker[i, "Studentt"] <- ttest[["p.value"]]
  } else {
    biomarker[i, "Studentt"] <- NA
  }
  biomarker[i, "Median"] <- median(Mut.lfc)
}
biomarker <- biomarker %>% na.omit()
significance <- filter(biomarker, wilcox <= 0.05, Median < 0.3)
#write.csv(significance,"D:/Winnie/Depmap/Duloxetine_CRC_biomarker.csv")
#Function to draw box plot
boxplot <- function(a){
  biomarker_signature <- cancer_signature[a,] %>% t() %>% as.data.frame()
  cancer_LFC[,'Type'] <- biomarker_signature[,a][match(rownames(cancer_LFC), rownames(biomarker_signature))] %>% as.character()
  cancer_LFC %>%
    ggplot( aes(x=Type, y=LFC.value, fill = Type )) +
    geom_boxplot(width = 0.5) +
    geom_dotplot(binaxis = 'y',
                 dotsize = 1,
                 stackdir = 'center') +
    scale_x_discrete(labels=c("WT", "Mutation")) + 
    labs(y="PRISM Primary screening LFC", x=a) +
    geom_text(aes(label= ifelse(cancer_LFC$Type == 1,as.character(cell),'')),
              hjust=-0.2,vjust=0.5, size = 3)
}

#Breast cancer
TNBC <- cell_line_info %>% filter(lineage_sub_subtype == "ERneg_HER2neg")

cancer_LFC <- Drug_LFC %>% rownames_to_column("gene") %>% 
  filter(rownames(Drug_LFC) %in% TNBC$DepMap_ID) %>% 
  column_to_rownames("gene")
#gene_signature <- read_csv("D:/Winnie/Depmap/raw data/CCLE_mutation_Winnie_full.csv")
#gene_signature <- gene_signature %>% column_to_rownames("X1")
cancer_signature <- gene_signature %>% select(rownames(cancer_LFC))

biomarker <- data.frame(gene = rownames(cancer_signature), wilcox= NA)
rownames(biomarker) <- biomarker$gene
for (i in 1: length(biomarker$gene)){
  gene <- cancer_signature[i,] %>% t() %>% as.data.frame()
  colnames(gene) <- "status"
  WT <- rownames(gene %>% rownames_to_column("gene") %>% 
                   filter(status == "0") %>% 
                   column_to_rownames("gene"))
  Mut <- rownames(gene %>% rownames_to_column("gene") %>% 
                    filter(status == "1") %>% 
                    column_to_rownames("gene"))
  WT.lfc <- cancer_LFC$LFC.value[match(WT,rownames(cancer_LFC))] %>% na.omit()
  Mut.lfc <- cancer_LFC$LFC.value[match(Mut, rownames(cancer_LFC))] %>% na.omit()
  if(length(Mut.lfc) >= 1) {
    wilcox <- wilcox.test(WT.lfc, Mut.lfc, alternative = "greater")
    biomarker[i, "wilcox"] <- wilcox[["p.value"]]
  } else {
    biomarker[i, "wilcox"] <- NA
  }
  if (length(Mut.lfc) >= 2){
    ttest <- t.test(WT.lfc, Mut.lfc, alternative = "greater")
    biomarker[i, "Studentt"] <- ttest[["p.value"]]
  } else {
    biomarker[i, "Studentt"] <- NA
  }
  biomarker[i, "Median"] <- median(Mut.lfc)
}
biomarker <- biomarker %>% na.omit()
significance <- filter(biomarker, wilcox <= 0.05, Median < 0.3)
write.csv(significance,"D:/Winnie/Depmap/Duloxetine_BC_biomarker.csv")
cancer_LFC$cell <- cell_line_info$stripped_cell_line_name[match(rownames(cancer_LFC), cell_line_info$DepMap_ID)]


boxplot("CELSR2")


#Primary or Metastasis

Drug_LFC$PorM <- cell_line_info$primary_or_metastasis[match(rownames(Drug_LFC), cell_line_info$DepMap_ID)]
Drug_LFC$PorM <- factor(Drug_LFC$PorM,levels=c("Primary", "Metastasis",NA))
Drug_LFC %>% ggplot(aes(x = as.factor(reorder(cancer_type, LFC.value, median)), LFC.value, fill = PorM)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(size=8, angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size=9),
        axis.title.x = element_blank(), legend.title = element_blank(),
        legend.key.size = unit(0.5,"line"), legend.position = "top") +
  ylab("Afatinib PRISM Primary Screening LFC")

Drug_LFC %>% filter(cancer_type == "Colon/Colorectal Cancer") %>% 
  ggplot(aes(x = PorM, LFC.value, fill = PorM)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(size=8),
        axis.title.y = element_text(size=9),
        axis.title.x = element_text(size=9), legend.title = element_blank(),
        legend.key.size = unit(0.5,"line"), legend.position = "top") +
  ylab("Afatinib PRISM Primary Screening LFC") +
  xlab("Colon/Colorectal Cancer") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
wilcox.test(Drug_LFC %>% filter(cancer_type == "Colon/Colorectal Cancer" & PorM == "Primary") %>% pull(`LFC.value`),
            Drug_LFC %>% filter(cancer_type == "Colon/Colorectal Cancer" & PorM == "Metastasis") %>% pull(`LFC.value`))
