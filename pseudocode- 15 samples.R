##INSTALL + LOAD PACKAGES
BiocManager::install(version = "3.21")
BiocManager::install("phyloseq")
library("devtools")
install_github("biobakery/maaslin3")
BiocManager::install("microbiome")
BiocManager::install("ALDEx2")
install.packages("tidyverse")
install.packages('vegan',
                 repos = c('https://vegandevs.r-universe.dev','https://cloud.r-project.org'))
install_github("statdivlab/radEmu")
install.packages("randomcoloR")

library(phyloseq)
library(maaslin3)
library(ALDEx2)
library(ggplot2)
library(tidyr)
library(stringr)
library(vegan)
library(stats)
library(radEmu)
library(randomcoloR)

##SET INPUT DIRECTORY
##ENA metadata TSV file
# define input directory
input_dir <- "/Users/brindha/Documents/R- microbiome analysis/input"
setwd(input_dir)

##BUILD METADATA FROM ENA TSV
#TSV file (ENA run table) 
tsv_file <- file.path(input_dir, "ena_sra-run_20260127-1506.tsv")

meta0 <- read.delim(tsv_file, stringsAsFactors = FALSE)

#extract D##S or D##P from description
code <- sub(".*(D[0-9]+[SP]).*", "\\1", meta0$description)

disease_status <- ifelse(grepl("S$", code), "healthy",
                         ifelse(grepl("P$", code), "periodontitis", NA))

meta <- data.frame(
  sample_id = meta0$accession,
  subject_id = code,
  forward_read = paste0(meta0$accession, "_1.fastq.gz"),
  reverse_read = paste0(meta0$accession, "_2.fastq.gz"),
  disease_status = disease_status,
  stringsAsFactors = FALSE
)

write.csv(meta, file.path(input_dir, "metadata_wholesome.csv"), row.names = FALSE)
View(meta)

##RUN KNEADDATA (QC + HUMAN READ REMOVAL)
## using human Bowtie2 database
#Samples
samples <- meta$sample_id

#Output folder
qc_dir <- "/Users/brindha/Documents/R- microbiome analysis/kneaddata_out"
dir.create(qc_dir, showWarnings = FALSE)

#Human DB
human_db <- "/Users/brindha/kneaddata_db"
human_db <- path.expand(human_db)

#kneadData
kneaddata_bin <- "/opt/homebrew/Caskroom/miniforge/base/envs/kneaddata_env/bin/kneaddata"
trimmomatic_bin <- "/opt/homebrew/Caskroom/miniforge/base/envs/kneaddata_env/share/trimmomatic-0.40-0/trimmomatic.jar"
trimmomatic_dir <- "/opt/homebrew/Caskroom/miniforge/base/envs/kneaddata_env/share/trimmomatic-0.40-0"
bowtie2_bin <- "/opt/homebrew/Caskroom/miniforge/base/envs/kneaddata_env/bin/bowtie2"


samples_ok <- samples[
  file.exists(file.path(input_dir, paste0(samples, "_1.fastq.gz"))) &
    file.exists(file.path(input_dir, paste0(samples, "_2.fastq.gz")))
]

for (s in samples_ok) {
  
  r1 <- file.path(input_dir, paste0(s, "_1.fastq.gz"))
  r2 <- file.path(input_dir, paste0(s, "_2.fastq.gz"))
  
  cmd <- paste(
    kneaddata_bin,
    "-i1", shQuote(r1),
    "-i2", shQuote(r2),
    "-o",  shQuote(qc_dir),
    "-db", shQuote(human_db),
    "--trimmomatic", shQuote(trimmomatic_dir),
    "--bowtie2", shQuote(bowtie2_bin),
    "--trimmomatic-options", shQuote("SLIDINGWINDOW:4:20 MINLEN:50"),
    "--bowtie2-options", shQuote("--very-sensitive --dovetail"),
    "--remove-intermediate-output",
    collapse = " "
  )
  
  message("Running: ", s)
  message(cmd)
  
  system(cmd)
}


##PARSE KRAKEN REPORTS INTO GENUS TABLE
## genus-by-sample count matrix
setwd("/Users/brindha/Documents/R- microbiome analysis/kneaddata_out")

#installed packages(c("tidyverse", "phyloseq"))

input_dir <- "/Users/brindha/Documents/R- microbiome analysis/input"

reports <- list.files(input_dir, pattern = "\\.kraken\\.report$", full.names = TRUE)

read_report <- function(f) {
  x <- read.table(f, sep = "\t", stringsAsFactors = FALSE, quote = "", fill = TRUE)
  colnames(x) <- c("pct", "reads_clade", "reads_taxon", "rank", "taxid", "name")
  x$name <- trimws(x$name)
  x <- x[x$rank == "G", c("name", "reads_clade")]
  x
}

tables <- lapply(reports, read_report)
names(tables) <- sub("\\.sub\\.kraken\\.report$", "", basename(reports))

all_taxa <- sort(unique(unlist(lapply(tables, \(x) x$name))))

mat <- sapply(tables, function(x) {
  v <- setNames(x$reads_clade, x$name)
  v[all_taxa] <- v[all_taxa]
  v[is.na(v)] <- 0
  v
})

genus_table <- data.frame(taxon = all_taxa, mat, check.names = FALSE)
write.csv(genus_table,
          file.path(input_dir, "kraken_genus_table.csv"),
          row.names = FALSE)

genus_table <- read.csv(file.path(input_dir, "kraken_genus_table.csv"))
View(genus_table)
write.csv(meta, file.path(input_dir, "kraken_genus_table.csv"), row.names = FALSE)

##BUILD RELATIVE ABUNDANCE TABLE
## so samples with different read depths can be compared
input_dir <- "/Users/brindha/Documents/R- microbiome analysis/input"

tax <- read.csv(file.path(input_dir, "kraken_genus_table.csv"), check.names = FALSE)
meta <- read.csv(file.path(input_dir, "metadata_wholesome.csv"))

# rownames for taxa table
rownames(tax) <- tax$taxon
tax$taxon <- NULL

# keep only shared samples
shared <- intersect(colnames(tax), meta$sample_id)
tax <- tax[, shared, drop = FALSE]
meta <- meta[match(shared, meta$sample_id), ]

tax_rel <- sweep(tax, 2, colSums(tax), "/")

##ALPHA DIVERSITY + STAT TEST
## if diversity differs between healthy vs periodontitis
shannon <- vegan::diversity(t(tax_rel), index = "shannon")
meta$shannon <- shannon

wilcox.test(shannon ~ disease_status, data = meta)

boxplot(shannon ~ disease_status, data = meta,
        xlab = "Group", ylab = "Shannon diversity")

library(ggplot2)

ggplot(meta, aes(x = disease_status, y = shannon, fill = disease_status)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  labs(x = "Group", y = "Shannon diversity") +
  theme_classic() +
  theme(legend.position = "none")

## BETA DIVERSITY (BRAY-CURTIS) + PCoA
## and visualize clustering by disease group
library(vegan)
library(ggplot2)

# Bray-Curtis distance (samples must be rows)
bray <- vegdist(t(tax_rel), method = "bray")

# PCoA
pcoa <- cmdscale(bray, k = 2, eig = TRUE)

pcoa_df <- data.frame(
  sample_id = rownames(pcoa$points),
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2]
)

# merge with metadata
pcoa_df <- merge(pcoa_df, meta, by = "sample_id")

# plot
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = disease_status)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = "PCoA1", y = "PCoA2")

##PER-SAMPLE METRICS FILE
## Purpose: Save alpha diversity metrics for each sample
library(vegan)

# relative abundance table already made: tax_rel
# meta already aligned

meta$shannon <- diversity(t(tax_rel), index = "shannon")
meta$simpson <- diversity(t(tax_rel), index = "simpson")
meta$richness <- specnumber(t(tax_rel))

write.csv(meta, file.path(input_dir, "per_sample_metrics.csv"), row.names = FALSE)
View(meta)

##HEATMAP OF ALL GENERA (ALL SAMPLES)
## Purpose: Visual overview of genus abundance patterns
library(tidyr)
library(ggplot2)

tax_rel_long <- as.data.frame(tax_rel)
tax_rel_long$genus <- rownames(tax_rel_long)

tax_long <- pivot_longer(
  tax_rel_long,
  cols = -genus,
  names_to = "sample_id",
  values_to = "rel_abundance"
)

p_heat <- ggplot(tax_long, aes(x = sample_id, y = genus, fill = rel_abundance)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Sample", y = "Genus", fill = "Rel. Abundance")

ggsave(file.path(input_dir, "all_genera_heatmap.png"),
       plot = p_heat, width = 12, height = 8)

##INDIVIDUAL SAMPLE GENUS BARPLOTS
## Purpose: Create one plot per sample showing all genera
## present in that sample (individual analysis output)

library(ggplot2)

outdir <- file.path(input_dir, "all_genera_per_sample_plots")
dir.create(outdir, showWarnings = FALSE)

for (s in colnames(tax_rel)) {
  
  df <- data.frame(
    genus = rownames(tax_rel),
    rel_abundance = tax_rel[, s],
    stringsAsFactors = FALSE
  )
  
  df <- df[df$rel_abundance > 0, ]
  df <- df[order(df$rel_abundance, decreasing = TRUE), ]
  
  p <- ggplot(df, aes(x = reorder(genus, rel_abundance), y = rel_abundance)) +
    geom_col() +
    coord_flip() +
    theme_classic() +
    labs(title = s, x = "Genus", y = "Relative abundance")
  
  ggsave(file.path(outdir, paste0(s, "_all_genera.png")),
         plot = p, width = 8, height = 10)
}

##EXPORT FINAL TABLES (COUNTS + RELATIVE ABUNDANCE)
write.csv(
  data.frame(genus = rownames(tax), tax, check.names = FALSE),
  file.path(input_dir, "genus_counts_table.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(genus = rownames(tax_rel), tax_rel, check.names = FALSE),
  file.path(input_dir, "genus_relative_abundance_table.csv"),
  row.names = FALSE
)

##EXPORT LONG-FORMAT RELATIVE ABUNDANCE TABLE
library(tidyr)

tax_rel_long <- as.data.frame(tax_rel)
tax_rel_long$genus <- rownames(tax_rel_long)

tax_long <- pivot_longer(
  tax_rel_long,
  cols = -genus,
  names_to = "sample_id",
  values_to = "rel_abundance"
)

write.csv(tax_long,
          file.path(input_dir, "genus_relative_abundance_long.csv"),
          row.names = FALSE)

View(tax_long)


##merge with metadata
##merge metadata + per-sample metrics
meta <- read.csv(file.path(input_dir, "metadata_wholesome.csv"))
metrics <- read.csv(file.path(input_dir, "per_sample_metrics.csv"))

meta2 <- merge(meta, metrics, by = "sample_id", all.x = TRUE)

write.csv(meta2,
          file.path(input_dir, "metadata_plus_metrics.csv"),
          row.names = FALSE)


## MERGE METADATA INTO LONG ABUNDANCE TABLE
tax_long <- read.csv(file.path(input_dir, "genus_relative_abundance_long.csv"))
meta <- read.csv(file.path(input_dir, "metadata_wholesome.csv"))

tax_long2 <- merge(tax_long, meta, by.x = "sample_id", by.y = "sample_id")

write.csv(tax_long2,
          file.path(input_dir, "genus_long_plus_metadata.csv"),
          row.names = FALSE)

head(meta2)
head(tax_long2)

##CLEAN METADATA FOR MAASLIN3
meta2 <- read.csv(file.path(input_dir, "metadata_plus_metrics.csv"))
View(meta2)
tax_long2 <- read.csv(file.path(input_dir, "genus_long_plus_metadata.csv"))
View(tax_long2)

meta_clean <- meta2[, c(
  "sample_id",
  "subject_id.x",
  "disease_status.x",
  "forward_read.x",
  "reverse_read.x",
  "shannon",
  "simpson",
  "richness"
)]

colnames(meta_clean) <- c(
  "sample_id",
  "subject_id",
  "disease_status",
  "forward_read",
  "reverse_read",
  "shannon",
  "simpson",
  "richness"
)

# drop samples without diversity values
meta_clean <- meta_clean[!is.na(meta_clean$shannon), ]

write.csv(meta_clean,
          file.path(input_dir, "per_sample_metrics_clean.csv"),
          row.names = FALSE)

View(meta_clean)

tax_rel <- read.csv(file.path(input_dir, "genus_relative_abundance_table.csv"),
                    check.names = FALSE)

rownames(tax_rel) <- tax_rel$genus
tax_rel$genus <- NULL

meta_clean <- read.csv(file.path(input_dir, "per_sample_metrics_clean.csv"))
rownames(meta_clean) <- meta_clean$sample_id

shared <- intersect(colnames(tax_rel), rownames(meta_clean))
tax_rel <- tax_rel[, shared, drop = FALSE]
meta_clean <- meta_clean[shared, , drop = FALSE]

##MAASLIN3 (ALL GENERA)
outdir <- file.path(input_dir, "maaslin3_genus_results")
dir.create(outdir, showWarnings = FALSE)

fit <- maaslin3(
  input_data = tax_rel,
  input_metadata = meta_clean,
  output = outdir,
  fixed_effects = c("disease_status"),
  normalization = "NONE",
  transform = "LOG"
)

outdir <- file.path(input_dir, "maaslin3_genus_results")

res_all <- read.delim(file.path(outdir, "all_results.tsv"))
res_sig <- read.delim(file.path(outdir, "significant_results.tsv"))

View(res_all)
res_all <- res_all[order(res_all$qval_individual), ]
head(res_all, 20)

res_all <- res_all[order(res_all$pval_individual), ]
head(res_all, 20)

##SALIVA MICROBES FILTERING
bad <- c("Homo", "Gammaretrovirus")

tax_rel_saliva <- tax_rel[!rownames(tax_rel) %in% bad, , drop = FALSE]

present <- rowSums(tax_rel_saliva > 0) >= 2
tax_rel_saliva <- tax_rel_saliva[present, , drop = FALSE]


write.csv(
  data.frame(genus = rownames(tax_rel_saliva), tax_rel_saliva, check.names = FALSE),
  file.path(input_dir, "genus_relative_abundance_saliva_filtered.csv"),
  row.names = FALSE
)

View(tax_rel_saliva)

##RELOAD SALIVA TABLE + COMPUTE SALIVA METRICS
tax_rel_saliva <- read.csv(
  file.path(input_dir, "genus_relative_abundance_saliva_filtered.csv"),
  check.names = FALSE
)

rownames(tax_rel_saliva) <- tax_rel_saliva$genus
tax_rel_saliva$genus <- NULL
tax_rel_saliva <- as.matrix(tax_rel_saliva)

library(vegan)

meta <- read.csv(file.path(input_dir, "metadata_wholesome.csv"))

shared <- intersect(colnames(tax_rel_saliva), meta$sample_id)
tax_rel_saliva <- tax_rel_saliva[, shared, drop = FALSE]
meta <- meta[match(shared, meta$sample_id), ]

meta$shannon  <- diversity(t(tax_rel_saliva), index = "shannon")
meta$simpson  <- diversity(t(tax_rel_saliva), index = "simpson")
meta$richness <- specnumber(t(tax_rel_saliva))

write.csv(meta,
          file.path(input_dir, "per_sample_metrics_saliva.csv"),
          row.names = FALSE)

##SALIVA PCoA + PERMANOVA
library(vegan)
library(ggplot2)

bray <- vegdist(t(tax_rel_saliva), method = "bray")
pcoa <- cmdscale(bray, k = 2, eig = TRUE)

pcoa_df <- data.frame(
  sample_id = rownames(pcoa$points),
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2]
)

pcoa_df <- merge(pcoa_df, meta, by = "sample_id")

ggplot(pcoa_df, aes(x = PC1, y = PC2, color = disease_status)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = "PCoA1 (Bray-Curtis)", y = "PCoA2 (Bray-Curtis)")

adonis2(bray ~ disease_status, data = meta, permutations = 999)

ggsave(file.path(input_dir, "saliva_pcoa_braycurtis.png"),
       width = 7, height = 5)


## MAASLIN3 (SALIVA ONLY)
library(maaslin3)

outdir2 <- file.path(input_dir, "maaslin3_saliva_results")
dir.create(outdir2, showWarnings = FALSE)

# tax_rel_saliva = genera x samples (relative abundance)
# meta = metadata (must match samples)

shared <- intersect(colnames(tax_rel_saliva), meta$sample_id)
tax_rel_saliva2 <- tax_rel_saliva[, shared, drop = FALSE]
meta2 <- meta[match(shared, meta$sample_id), ]
rownames(meta2) <- meta2$sample_id

fit2 <- maaslin3(
  input_data = tax_rel_saliva2,
  input_metadata = meta2,
  output = outdir2,
  fixed_effects = c("disease_status"),
  normalization = "NONE",
  transform = "LOG"
)

res_all2 <- read.delim(file.path(outdir2, "all_results.tsv"))
res_all2 <- res_all2[order(res_all2$qval_individual), ]
View(res_all2)

outdir2 <- file.path(input_dir, "maaslin3_saliva_results")

res_all2 <- read.delim(file.path(outdir2, "all_results.tsv"))
res_all2 <- res_all2[order(res_all2$pval_individual), ]

##DIAGNOSTIC MODEL (LOGISTIC REGRESSION + 5-FOLD CV)
## Purpose: Prototype saliva diagnostic classifier for
## periodontitis vs healthy using genus features
library(pROC)

# tax_rel_saliva: genera x samples
X <- t(tax_rel_saliva)   # samples x genera

# outcome
y <- meta$disease_status
y <- factor(y, levels = c("healthy", "periodontitis"))

# keep only complete samples
keep <- !is.na(y)
X <- X[keep, , drop = FALSE]
y <- y[keep]


X_log <- log10(X + 1e-6)

set.seed(1)

n <- nrow(X_log)
folds <- sample(rep(1:5, length.out = n))

pred_prob <- rep(NA, n)

for (k in 1:5) {
  
  train <- folds != k
  test  <- folds == k
  
  df_train <- data.frame(y = y[train], X_log[train, , drop = FALSE])
  df_test  <- data.frame(X_log[test, , drop = FALSE])
  
  fit <- glm(y ~ ., data = df_train, family = binomial())
  
  pred_prob[test] <- predict(fit, newdata = df_test, type = "response")
}

roc_obj <- roc(y, pred_prob, levels = c("healthy", "periodontitis"))
auc(roc_obj)
plot(roc_obj, main = "Saliva diagnostic prototype (5-fold CV)")

##SAVE PREDICTIONS TABLE
## Purpose: Save per-sample predicted probability output
pred_table <- data.frame(
  sample_id = rownames(X_log),
  true_status = y,
  predicted_prob_periodontitis = pred_prob,
  stringsAsFactors = FALSE
)

write.csv(pred_table,
          file.path(input_dir, "saliva_classifier_predictions.csv"),
          row.names = FALSE)

View(pred_table)


##CONFUSION MATRIX + BASIC PERFORMANCE (THR=0.5)
pred_class <- ifelse(pred_prob >= 0.5, "periodontitis", "healthy")
pred_class <- factor(pred_class, levels = c("healthy", "periodontitis"))

table(True = y, Predicted = pred_class)

conf_mat <- table(True = y, Predicted = pred_class)
conf_mat

accuracy <- sum(diag(conf_mat)) / sum(conf_mat)

sensitivity <- conf_mat["periodontitis", "periodontitis"] / sum(conf_mat["periodontitis", ])
specificity <- conf_mat["healthy", "healthy"] / sum(conf_mat["healthy", ])

accuracy
sensitivity
specificity

##BEST THRESHOLD FROM ROC
##Find optimal probability cutoff instead of 0.5

best <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
best

thr <- as.numeric(best["threshold"])
thr
length(thr)
## NOTE: Your code below uses pred_prob2/y2/ok but does not
## define them before using. That is why you previously got errors.

##SAVE FINAL FIGURES
pred_class_best <- ifelse(pred_prob2 >= thr, "periodontitis", "healthy")
pred_class_best <- factor(pred_class_best, levels = c("healthy", "periodontitis"))

length(pred_class_best)

conf_mat2 <- table(True = y2, Predicted = pred_class_best)
conf_mat2

accuracy2 <- sum(diag(conf_mat2)) / sum(conf_mat2)
sensitivity2 <- conf_mat2["periodontitis", "periodontitis"] / sum(conf_mat2["periodontitis", ])
specificity2 <- conf_mat2["healthy", "healthy"] / sum(conf_mat2["healthy", ])

perf <- data.frame(
  threshold = thr,
  accuracy = accuracy2,
  sensitivity = sensitivity2,
  specificity = specificity2
)

write.csv(perf,
          file.path(input_dir, "saliva_classifier_performance.csv"),
          row.names = FALSE)

diagnostic_table2 <- data.frame(
  sample_id = rownames(X_log)[ok],
  true_status = y2,
  predicted_prob_periodontitis = pred_prob2,
  predicted_class = pred_class_best,
  threshold_used = thr,
  stringsAsFactors = FALSE
)

write.csv(diagnostic_table2,
          file.path(input_dir, "saliva_diagnostic_table_best_threshold.csv"),
          row.names = FALSE)

View(diagnostic_table2)


library(ggplot2)

p1 <- ggplot(meta, aes(x = disease_status, y = shannon, fill = disease_status)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Group", y = "Shannon diversity")

ggsave(file.path(input_dir, "FIG1_shannon_boxplot.png"),
       plot = p1, width = 6, height = 4)

p2 <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = disease_status)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = "PCoA1 (Bray-Curtis)", y = "PCoA2 (Bray-Curtis)")

ggsave(file.path(input_dir, "FIG2_pcoa_braycurtis.png"),
       plot = p2, width = 6, height = 4)

png(file.path(input_dir, "FIG3_ROC_curve.png"), width = 800, height = 600)
plot(roc_obj, main = "Saliva diagnostic prototype (5-fold CV)")
dev.off()

write.csv(as.data.frame.matrix(conf_mat2),
          file.path(input_dir, "FIG4_confusion_matrix.csv"))

browseURL(input_dir)

png(file.path(input_dir, "FIG3_ROC_curve_clean.png"), width = 800, height = 600)

plot(
  roc_obj,
  legacy.axes = TRUE,
  print.auc = TRUE,
  main = "Saliva diagnostic prototype (5-fold CV)"
)

# Correct baseline diagonal
abline(a = 0, b = 1, lty = 2, col = "gray40")

dev.off()

browseURL(file.path(input_dir, "FIG3_ROC_curve_clean.png"))


------------------------------------------------------------------------
#oral health#
------------------------------------------------------------------------
library(readr)
library(dplyr)
library(tidyr)

# INPUT DIRECTORY
input_dir <- "/Users/brindha/Documents/R- microbiome analysis/input"
setwd(input_dir)

# INPUT FILES
results_path   <- "all_results.tsv"
abundance_path <- "genus_relative_abundance_saliva_filtered.csv"
metadata_path  <- "per_sample_metrics_saliva.csv"

# CHECK FILES EXIST
stopifnot(file.exists(results_path))
stopifnot(file.exists(abundance_path))
stopifnot(file.exists(metadata_path))

# LOAD DATA
res  <- read_tsv(results_path, show_col_types = FALSE)
abun <- read_csv(abundance_path, show_col_types = FALSE)
meta <- read_csv(metadata_path, show_col_types = FALSE)

# SELECT PERIODONTITIS-ASSOCIATED GENERA (p <= 0.05)
sig <- res %>%
  filter(
    model == "abundance",
    metadata == "disease_status",
    value == "periodontitis",
    !is.na(pval_individual),
    pval_individual <= 0.05
  ) %>%
  select(feature, coef) %>%
  distinct()

if (nrow(sig) == 0) stop("No genera passed p <= 0.05. Cannot compute health score.")

# LONG FORMAT ABUNDANCE TABLE
abun_long <- abun %>%
  pivot_longer(
    cols = -genus,
    names_to = "sample_id",
    values_to = "rel_abundance"
  )

# SCORE
scores <- abun_long %>%
  inner_join(sig, by = c("genus" = "feature")) %>%
  mutate(weighted = rel_abundance * coef) %>%
  group_by(sample_id) %>%
  summarise(periodontitis_linear_score = sum(weighted), .groups = "drop")

# HEALTH SCORE (0–100)
s_min <- min(scores$periodontitis_linear_score)
s_max <- max(scores$periodontitis_linear_score)

scores <- scores %>%
  mutate(
    scaled = (periodontitis_linear_score - s_min) / (s_max - s_min + 1e-12),
    health_score_0_100 = 100 * (1 - scaled)
  ) %>%
  select(sample_id, periodontitis_linear_score, health_score_0_100)

# MERGE WITH METADATA + EXPORT
final <- meta %>%
  select(sample_id, disease_status) %>%
  left_join(scores, by = "sample_id") %>%
  arrange(desc(health_score_0_100))

write_csv(final, "saliva_health_score_results.csv")

print(final)
cat("\nSaved: saliva_health_score_results.csv\n")



library(dplyr)

tax_long <- read.csv("genus_long_plus_metadata.csv")

# Example: paper_taxa should be replaced by the taxa from the paper
paper_taxa <- c("Porphyromonas", "Treponema", "Fusobacterium")

tax_long %>%
  filter(genus %in% paper_taxa) %>%
  group_by(genus, disease_status) %>%
  summarise(mean_abundance = mean(rel_abundance), .groups="drop") %>%
  arrange(genus, disease_status)



-----------------------------------------------------------------------
#orah health#
-----------------------------------------------------------------------
  
## ============================================================
## ORAL HEALTH SCORE PIPELINE (MaAsLin3 → Weighted Score → ROC)
## Works with: all_results.tsv + genus_relative_abundance.csv + metadata.csv
## ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

## -----------------------------
## 1) Paths (EDIT)
## -----------------------------
input_dir <- "/Users/brindha/Documents/R- microbiome analysis/input"
out_dir   <- "/Users/brindha/Documents/R- microbiome analysis/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(input_dir)

results_path   <- "all_results.tsv"
abundance_path <- "genus_relative_abundance_saliva_filtered.csv"
metadata_path  <- "per_sample_metrics_saliva.csv"

## -----------------------------
## 2) Load files
## -----------------------------
res  <- read_tsv(results_path, show_col_types = FALSE)
abun <- read_csv(abundance_path, show_col_types = FALSE)
meta <- read_csv(metadata_path, show_col_types = FALSE)

## -----------------------------
## 3) Basic cleaning (not strict)
## -----------------------------
meta <- meta %>%
  rename_with(tolower) %>%
  mutate(
    sample_id = as.character(sample_id),
    disease_status = factor(disease_status)
  )

## abundance file must have first column = genus
abun <- abun %>%
  rename_with(tolower)

if (!("genus" %in% colnames(abun))) {
  stop("Abundance file must contain a column named 'genus'.")
}

## -----------------------------
## 4) Extract significant genera (periodontitis-associated)
## -----------------------------
sig <- res %>%
  filter(
    model == "abundance",
    metadata == "disease_status",
    value == "periodontitis",
    !is.na(pval_individual),
    pval_individual <= 0.05
  ) %>%
  select(feature, coef, pval_individual) %>%
  distinct() %>%
  arrange(pval_individual)

write_csv(sig, file.path(out_dir, "significant_genera_maaslin3.csv"))

if (nrow(sig) == 0) stop("No genera passed p <= 0.05. Cannot compute score.")


## 5) Convert abundance table to long format

abun_long <- abun %>%
  pivot_longer(
    cols = -genus,
    names_to = "sample_id",
    values_to = "rel_abundance"
  ) %>%
  mutate(
    sample_id = as.character(sample_id),
    rel_abundance = as.numeric(rel_abundance)
  )


## 6) Compute weighted Periodontitis Linear Score

scores <- abun_long %>%
  inner_join(sig, by = c("genus" = "feature")) %>%
  mutate(weighted = rel_abundance * coef) %>%
  group_by(sample_id) %>%
  summarise(periodontitis_linear_score = sum(weighted, na.rm = TRUE), .groups = "drop")


## 7) Convert to 0–100 Oral Health Score

s_min <- min(scores$periodontitis_linear_score, na.rm = TRUE)
s_max <- max(scores$periodontitis_linear_score, na.rm = TRUE)

scores <- scores %>%
  mutate(
    scaled = (periodontitis_linear_score - s_min) / (s_max - s_min + 1e-12),
    oral_health_score_0_100 = 100 * (1 - scaled)
  ) %>%
  select(sample_id, periodontitis_linear_score, oral_health_score_0_100)

write_csv(scores, file.path(out_dir, "oral_health_scores.csv"))


## 8) Merge with metadata

final <- meta %>%
  select(sample_id, disease_status) %>%
  left_join(scores, by = "sample_id") %>%
  arrange(desc(oral_health_score_0_100))

write_csv(final, file.path(out_dir, "saliva_health_score_results.csv"))

print(final)


## 9) Plot: Health score by disease status

p1 <- ggplot(final, aes(x = disease_status, y = oral_health_score_0_100, fill = disease_status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  labs(
    title = "Oral Health Score (0–100) by Disease Status",
    x = "",
    y = "Oral Health Score (higher = healthier-like)"
  )

ggsave(file.path(out_dir, "oral_health_score_boxplot.png"), p1, width = 6, height = 4)


## 10) Plot: Ranked samples

p2 <- ggplot(final, aes(x = reorder(sample_id, oral_health_score_0_100),
                        y = oral_health_score_0_100,
                        fill = disease_status)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Ranked Oral Health Score per Sample",
    x = "Sample",
    y = "Oral Health Score (0–100)"
  )

ggsave(file.path(out_dir, "oral_health_score_ranked.png"), p2, width = 7, height = 5)


## 11) Statistical test (Wilcoxon)

if (length(unique(final$disease_status)) == 2) {
  wtest <- wilcox.test(oral_health_score_0_100 ~ disease_status, data = final)
  capture.output(wtest, file = file.path(out_dir, "oral_health_score_wilcox.txt"))
}


## 12) OPTIONAL: ROC + AUC (no extra packages)
##     Healthy vs Periodontitis

roc_df <- final %>%
  filter(!is.na(oral_health_score_0_100)) %>%
  mutate(
    y = ifelse(disease_status == "periodontitis", 1, 0),
    score = 100 - oral_health_score_0_100
  )

if (length(unique(roc_df$y)) == 2) {
  
  thresholds <- sort(unique(roc_df$score))
  roc_curve <- lapply(thresholds, function(t) {
    pred <- ifelse(roc_df$score >= t, 1, 0)
    tp <- sum(pred == 1 & roc_df$y == 1)
    fp <- sum(pred == 1 & roc_df$y == 0)
    tn <- sum(pred == 0 & roc_df$y == 0)
    fn <- sum(pred == 0 & roc_df$y == 1)
    
    tpr <- tp / (tp + fn + 1e-12)
    fpr <- fp / (fp + tn + 1e-12)
    
    data.frame(threshold = t, TPR = tpr, FPR = fpr)
  }) %>% bind_rows()
  
  roc_curve <- roc_curve %>% arrange(FPR, TPR)
  
  # AUC by trapezoid rule
  auc <- sum(diff(roc_curve$FPR) * (head(roc_curve$TPR, -1) + tail(roc_curve$TPR, -1)) / 2)
  
  write_csv(roc_curve, file.path(out_dir, "roc_curve.csv"))
  
  p3 <- ggplot(roc_curve, aes(x = FPR, y = TPR)) +
    geom_line(linewidth = 1) +
    geom_abline(linetype = "dashed") +
    theme_bw() +
    labs(
      title = paste0("ROC Curve (AUC = ", round(auc, 3), ")"),
      x = "False Positive Rate",
      y = "True Positive Rate"
    )
  
  ggsave(file.path(out_dir, "roc_curve.png"), p3, width = 5, height = 5)
  
  cat("\nROC AUC =", round(auc, 3), "\n")
}

cat("\nSaved outputs to:\n", out_dir, "\n")
