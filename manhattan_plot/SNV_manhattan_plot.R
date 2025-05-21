library(topr)
library(ggplot2)
# Read GWAS data from a file
known <- read.table("./initial_632656.txt", header = TRUE, stringsAsFactors = FALSE)
novel <- read.table("./overlaying_2838.txt", header = TRUE, stringsAsFactors = FALSE)
manhattan_plot <- manhattan(list(known, novel), color=c("#648FFF","#FFB000"), legend_labels=c("All", "Overlapping with sORFs"), label_fontface = "italic", scale=0.9, nudge_y = 20, sign_thresh = 5e-8, protein_coding_only = TRUE)
# Save the plot as PNG
ggsave("./Manhattan.png", plot = manhattan_plot, width = 40, height = 15, units = "cm", dpi = 300)
