library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)

all <- args[1]
pcr <- args[2]
plot_dir <- args[3]

all_basename <- basename(all)
pcr_basename <- basename(pcr)

sample_name <- sub("^all_", "", all_basename)
sample_name <- sub(".rpkm$", "", sample_name)


# all  <- "./../rpkm_calculation/nookaew_cm/RPKM/all_nookaew_cm.rpkm"
# pcr <- "./../rpkm_calculation/nookaew_cm/RPKM/pcr_nookaew_cm.rpkm"
# plot_dir <- "./../rpkm_calculation/nookaew_cm/Plots"

# all  <- "./../rpkm_calculation/hes_star/RPKM/all_hes_star.rpkm"
# pcr <- "./../rpkm_calculation/hes_star/RPKM/pcr_hes_star.rpkm"
# plot_dir <- "./../rpkm_calculation/hes_star/Plots"

# all  <- "./../rpkm_calculation/ebna_hisat/RPKM/all_ebna_hisat.rpkm"
# pcr <- "./../rpkm_calculation/ebna_hisat/RPKM/pcr_ebna_hisat.rpkm"

data <- read.csv(all, header = FALSE, col.names = c("gene", "rpkm", "chr"), sep = "\t", stringsAsFactors = FALSE)
data2 <- read.table(pcr, header = FALSE, col.names = c("gene", "rpkm", "chr"), sep = "\t", stringsAsFactors = FALSE)

Q1_all <- quantile(data$rpkm, 0.25)
Q3_all <- quantile(data$rpkm, 0.75)
IQR_all <- Q3_all - Q1_all
outlier_threshold_all <- Q3_all + 1.5 * IQR_all

Q1_pcr <- quantile(data2$rpkm, 0.25)
Q3_pcr <- quantile(data2$rpkm, 0.75)
IQR_pcr <- Q3_pcr - Q1_pcr
outlier_threshold_pcr <- Q3_pcr + 1.5 * IQR_pcr

top_outliers_all <- data %>%
  filter(rpkm > outlier_threshold_all) %>%
  arrange(desc(rpkm)) %>%
  head(10)
print("top_outliers_all")
print(top_outliers_all)

top_outliers_pcr <- data2 %>%
  filter(rpkm > outlier_threshold_pcr) %>%
  arrange(desc(rpkm)) %>%
  head(10)
print("top_outliers_pcr")
print(top_outliers_pcr)

if ("I" %in% data2$chr) {
  chromosome_levels <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "Mito")
} else {
  chromosome_levels <- c(as.character(1:22), "X", "Y", "MT")
}

data_sorted <- data %>%
  mutate(
    chr = factor(chr, levels = chromosome_levels),
    chr_order = as.numeric(chr),
    chr_pos = 1:n(),
    jittered_chr = as.numeric(chr) + runif(n(), min = -0.3, max = 0.3)
  )

data_sorted_pcr <- data2 %>%
  mutate(
    chr = factor(chr, levels = chromosome_levels),
    chr_order = as.numeric(chr),
    chr_pos = 1:n(),
    jittered_chr = as.numeric(chr) + runif(n(), min = -0.3, max = 0.3)
  )

top_outliers_all <- data_sorted %>%
  filter(gene %in% top_outliers_all$gene)

top_outliers_pcr <- data_sorted_pcr %>%
  filter(gene %in% top_outliers_pcr$gene)

common_genes <- intersect(top_outliers_pcr$gene, top_outliers_all$gene)

highlight_genes_uniq <- data_sorted_pcr %>%
  filter(gene %in% common_genes)

highlight_genes_all <- data_sorted %>%
  filter(gene %in% common_genes)

plot_title <- paste(sample_name, "\nRPKM Values of PCR Reads")
g_pcr <- ggplot(data_sorted_pcr, aes(x = jittered_chr, y = rpkm, color = factor(chr_order %% 2))) +
  geom_point(alpha = 0.5) +
  geom_text_repel(
    data = top_outliers_pcr,
    aes(label = gene),
    color = "red",
    box.padding = 0.6,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c("0" = "blue", "1" = "darkblue")) +
  scale_x_continuous(
    breaks = 1:length(chromosome_levels),
    labels = chromosome_levels
  ) +
  # scale_y_log10(
  #   breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000),
  #   labels = c("0.00001", "0.0001", "0.001", "0.01", "0.1", "1", "10", "100", "1000", "10000", "100000")
  # ) +
  labs(title = plot_title,
       x = "Chromosome",
       y = "RPKM") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    plot.title = element_text(size = 30, hjust = 0.5)
  )



plot_title <- paste(sample_name, "\nRPKM Values of All Reads")
g_all <- ggplot(data_sorted, aes(x = jittered_chr, y = rpkm, color = factor(chr_order %% 2))) +
  geom_point(alpha = 0.5) +
  # geom_text_repel(
  #   data = highlight_genes_all,
  #   aes(label = gene),
  #   color = "green",
  #   box.padding = 0.8,
  #   max.overlaps = Inf
  # ) +
  geom_text_repel(
    data = top_outliers_all,
    aes(label = gene),
    color = "red",
    box.padding = 0.6,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c("0" = "blue", "1" = "darkblue")) +
  scale_x_continuous(
    breaks = 1:length(chromosome_levels),
    labels = chromosome_levels
  ) +
  # scale_y_log10(
  #   breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000),
  #   labels = c("0.00001", "0.0001", "0.001", "0.01", "0.1", "1", "10", "100", "1000", "10000", "100000")
  # ) +
  labs(title = plot_title,
       x = "Chromosome",
       y = "RPKM") +
  theme_minimal() +
  scale_fill_manual(values = c("All Reads" = "blue", "PCR Unique Reads" = "red")) +  # Set fill color
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    plot.title = element_text(size = 30, hjust = 0.5)
  )

# Combine plots
plo <- g_all + g_pcr + plot_layout(ncol = 2)

hist_file <- file.path(plot_dir, "rpkm_mplot.png")
ggsave(hist_file, plo, width = 16, height = 10, dpi = 300)


plot_title <- paste(sample_name, "\nAll Reads vs PCR Unique Reads")
data$Dataset <- "All Reads"
data2$Dataset <- "PCR Unique Reads"
combined_data <- rbind(data, data2)
top_outliers_combined <- rbind(
  top_outliers_all %>% mutate(Dataset = "All Reads"),
  top_outliers_pcr %>% mutate(Dataset = "PCR Unique Reads")
)
# q <- ggplot(combined_data, aes(x = Dataset, y = rpkm, fill = Dataset)) +
#   geom_boxplot(alpha = 0.7, outlier.color = "red", outlier.size = 2) +
#   scale_y_log10(
#     breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000),
#     labels = c("0.00001", "0.0001", "0.001", "0.01", "0.1", "1", "10", "100", "1000", "10000", "100000")
#   ) +
#   scale_fill_manual(values = c("All Reads" = "blue", "PCR Unique Reads" = "red")) +  # Correct label for "PCR Unique Reads"
#   labs(
#     title = plot_title,
#     x = "Dataset",
#     y = "RPKM (log scale)"
#   ) +
#   theme_minimal() +
#   theme(
#     legend.position = "none",
#     axis.text.x = element_text(size = 15, color = "black"),
#     axis.text.y = element_text(size = 15, color = "black"),
#     axis.title.x = element_text(size = 17, face = "bold"),
#     axis.title.y = element_text(size = 17, face = "bold"),
#     plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
#     panel.grid.major = element_line(color = "gray", linetype = "dashed"),
#     panel.grid.minor = element_blank()
#   ) + 
#   coord_flip()
#
# # Save the plot
# hist_file <- file.path(plot_dir, "rpkm_hist.png")
# ggsave(hist_file, q, width = 12, height = 8, dpi = 300)

plot_title <- paste(sample_name, "\nAll Reads vs PCR Unique Reads")
p <- ggplot(combined_data, aes(x = rpkm, fill = Dataset)) +
  geom_density(alpha = 0.7, position = "identity") +  # Add histogram layer
  scale_x_log10(
    breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000),
    labels = c("0.00001", "0.0001", "0.001", "0.01", "0.1", "1", "10", "100", "1000", "10000", "100000")
  ) +
  scale_fill_manual(values = c("All Reads" = "blue", "PCR Unique Reads" = "red")) +
  labs(
    title = plot_title,
    x = "RPKM (log scale)",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 25),  # Increase legend text size
    legend.title = element_text(size = 30, face = "bold"),  # Increase legend title size
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black"),
    axis.title.x = element_text(size = 28, face = "bold"),
    axis.title.y = element_text(size = 28, face = "bold"),
    plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "gray", linetype = "dashed"),
    panel.grid.minor = element_blank()
  )

hist_file <- file.path(plot_dir, "rpkm_density.png")
ggsave(hist_file, p, width = 12, height = 8, dpi = 300)
