
# Import required libraries
library(ggplot2)
library(dplyr)
library(readxl)
library(dplyr)
library(DESeq2)
library(edgeR)
library(enrichR)

# Count Matrix
# folder where all xlsx are extracted
files <- list.files("GSE309456_RAW", pattern="\\.xlsx$", full.names=TRUE)
length(files)
#Extracting gene ID and total exon reads
read_one <- function(f) {
  df <- read_xlsx(f)
  df2 <- df %>%
    select(gene = `Name`,
           count = `Total exon reads`) %>%
    filter(!is.na(gene))

  v <- df2$count
  names(v) <- df2$gene
  return(v)
}

counts_list <- lapply(files, read_one)

# keep common genes across all samples
common_genes <- Reduce(intersect, lapply(counts_list, names))
counts_list <- lapply(counts_list, function(v) v[common_genes])

count_mat <- do.call(cbind, counts_list)
rownames(count_mat) <- common_genes

sample_names <- gsub("^GSM[0-9]+_|_GE_\\.xlsx$", "", basename(files))
colnames(count_mat) <- sample_names
count_mat

write.csv(count_mat, "counts_matrix.csv")
colnames(count_mat)

# Reading Count Data
count_mat <- as.matrix(read.csv("counts_matrix.csv",
                                row.names = 1,
                                check.names = FALSE))

mode(count_mat) <- "numeric"

# CREATE METADATA
condition <- ifelse(grepl("p62KO", colnames(count_mat), ignore.case = TRUE),
                    "p62KO",
                    "MOC2")

coldata <- data.frame(
  row.names = colnames(count_mat),
  condition = factor(condition)
)

# Set reference
coldata$condition <- relevel(coldata$condition, ref = "MOC2")

# DESEQ2 ANALYSIS
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = coldata,
                              design = ~ condition)

dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)

# PCA plot 
vsd <- vst(dds, blind = FALSE)

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot") +
  theme_minimal()

ggsave("PCA_plot.png", width = 6, height = 5)

res <- results(dds, contrast = c("condition", "p62KO", "MOC2"))

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

res_df$padj[is.na(res_df$padj)] <- 1

lfc_threshold <- log2(1.5)

sig_deseq <- subset(res_df,
                    padj < 0.1 &
                      abs(log2FoldChange) > lfc_threshold)

up_deseq <- sig_deseq[sig_deseq$log2FoldChange > 0, ]
down_deseq <- sig_deseq[sig_deseq$log2FoldChange < 0, ]

write.csv(res_df, "DESeq2_all_genes.csv", row.names = FALSE)
write.csv(sig_deseq, "DESeq2_significant.csv", row.names = FALSE)
write.csv(up_deseq, "DESeq2_UP.csv", row.names = FALSE)
write.csv(down_deseq, "DESeq2_DOWN.csv", row.names = FALSE)

# Volcano Plot for DESeq2

res_df$significance <- "Not Significant"
res_df$significance[res_df$padj < 0.1 &
                      abs(res_df$log2FoldChange) > lfc_threshold] <- "Significant"

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.6) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  ggtitle("DESeq2 Volcano Plot")

ggsave("DESeq2_Volcano.png", width = 6, height = 5)

# edgeR Analysis

library(edgeR)
group <- factor(condition)
y <- DGEList(counts = count_mat, group = group)
keep <- filterByExpr(y)
y <- y[keep,, keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~ group)
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)

edgeR_res <- topTags(qlf, n = Inf)$table
edgeR_res$gene <- rownames(edgeR_res)

edgeR_res$FDR[is.na(edgeR_res$FDR)] <- 1

edgeR_sig <- subset(edgeR_res,
                    FDR < 0.1 &
                      abs(logFC) > lfc_threshold)

up_edgeR <- edgeR_sig[edgeR_sig$logFC > 0, ]
down_edgeR <- edgeR_sig[edgeR_sig$logFC < 0, ]

write.csv(edgeR_res, "edgeR_all_genes.csv", row.names = FALSE)
write.csv(edgeR_sig, "edgeR_significant.csv", row.names = FALSE)
write.csv(up_edgeR, "edgeR_UP.csv", row.names = FALSE)
write.csv(down_edgeR, "edgeR_DOWN.csv", row.names = FALSE)

# Volcano Plot for edgeR
edgeR_res$significance <- "Not Significant"
edgeR_res$significance[edgeR_res$FDR < 0.1 &
                         abs(edgeR_res$logFC) > lfc_threshold] <- "Significant"

ggplot(edgeR_res, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = significance), alpha = 0.6) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  ggtitle("edgeR Volcano Plot")

ggsave("edgeR_Volcano.png", width = 6, height = 5)

# COMMON genes between edgeR and DESeq2

common_up <- intersect(up_deseq$gene, up_edgeR$gene)
common_down <- intersect(down_deseq$gene, down_edgeR$gene)
common_all <- intersect(sig_deseq$gene, edgeR_sig$gene)

length(common_all)

# GO Enrichment Analysis
# 1.EnrichR
# 2.DAVID

library(enrichR)

dbs <- c("GO_Biological_Process_2023",
         "GO_Molecular_Function_2023",
         "GO_Cellular_Component_2023")

enrich_up <- enrichr(common_up, dbs)
enrich_down <- enrichr(common_down, dbs)

go_up_bp <- enrich_up[["GO_Biological_Process_2023"]]
go_down_bp <- enrich_down[["GO_Biological_Process_2023"]]

# Filtering Top GO Terms
go_up_top <- go_up_bp %>%
  filter(Adjusted.P.value < 0.05) %>%
  arrange(Adjusted.P.value) %>%
  head(10)

go_down_top <- go_down_bp %>%
  filter(Adjusted.P.value < 0.05) %>%
  arrange(Adjusted.P.value) %>%
  head(10)

# Barplot GO UP genes

ggplot(go_up_top,
       aes(x = reorder(Term, -log10(Adjusted.P.value)),
           y = -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  coord_flip() +
  theme_minimal() +
  ggtitle("GO Biological Process (UP Genes)")

# Barplot GO down genes
ggplot(go_down_top,
       aes(x = reorder(Term, -log10(Adjusted.P.value)),
           y = -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity", fill = "darkred") +
  coord_flip() +
  theme_minimal() +
  ggtitle("GO Biological Process (DOWN Genes)")

write.csv(go_up_top, "GO_UP_Top10.csv", row.names = FALSE)
write.csv(go_down_top, "GO_DOWN_Top10.csv", row.names = FALSE)

# DAVID Analysis
write.table(common_up, "DAVID_common_UP.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(common_down, "DAVID_common_DOWN.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Files Downloaded from the DAVID official Database

up_data <- read.csv("DAVID_UP.csv")
down_data <- read.csv("DAVID_DOWN.csv")

# Filetering the significant terms FDR < 0.05

up_sig <- up_data %>% filter(FDR < 0.05)
down_sig <- down_data %>% filter(FDR < 0.05)

plot_go <- function(data, category_name, title_name, color_fill) {
  
  df <- data %>%
    filter(Category == category_name) %>%
    arrange(FDR) %>%
    head(10)
  
  df$logFDR <- -log10(df$FDR)
  
  ggplot(df,
         aes(x = reorder(Term, logFDR),
             y = logFDR)) +
    geom_bar(stat = "identity", fill = color_fill) +
    coord_flip() +
    theme_minimal() +
    labs(title = title_name,
         x = "GO Term",
         y = "-log10(FDR)")
}


# DAVID Up gene list and plots

plot_go(up_sig, "GOTERM_BP_DIRECT",
        "UP Genes - GO Biological Process", "steelblue")

plot_go(up_sig, "GOTERM_MF_DIRECT",
        "UP Genes - GO Molecular Function", "steelblue")

plot_go(up_sig, "GOTERM_CC_DIRECT",
        "UP Genes - GO Cellular Component", "steelblue")


# DAVID Down gene list and plots

plot_go(down_sig, "GOTERM_BP_DIRECT",
        "DOWN Genes - GO Biological Process", "firebrick")

plot_go(down_sig, "GOTERM_MF_DIRECT",
        "DOWN Genes - GO Molecular Function", "firebrick")

plot_go(down_sig, "GOTERM_CC_DIRECT",
        "DOWN Genes - GO Cellular Component", "firebrick")
