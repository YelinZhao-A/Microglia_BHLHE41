# Identify genes that has the same expression pattern in both bulk data sets
library(readxl)
library(dplyr)
dir.create("output", showWarnings = FALSE)

# read and clean data####
bulk1 <- read_excel("Merged_data_reorganized.xlsx", sheet = "DS1_DKO data") %>% 
  select(Gene_Symbol,`logFC`,`Adjust P value`, Dataset ) %>% 
  filter(`Adjust P value` < 0.05) %>% filter(abs(`logFC`)>0.5)
bulk2 <- read_excel("Merged_data_reorganized.xlsx", sheet = "DS2_own") %>% 
  select(Gene_Symbol,`Fold change`,`Adjust P value`, Dataset ) %>% 
  filter(`Adjust P value` < 0.05) 
bulk2$logFC <- sign(bulk2$`Fold change`) * log2(abs(bulk2$`Fold change`))
bulk2 <- bulk2 %>% filter(abs(logFC)>0.5)

merged <- merge(bulk1, bulk2, by = "Gene_Symbol", all = TRUE)

# make a scatter plot####
library(ggplot2)
library(ggrepel)
genes_to_label <- c("Cd22","Cacng1","Upk1b",'Ptprg','Plaur','Ccr5','Fcgr2b','Ciita')
p=ggplot(merged, aes(x = logFC.x, y = logFC.y)) +
  geom_point() +  # Plot all the points
  geom_point(data = subset(merged, Gene_Symbol %in% genes_to_label), 
             aes(x = logFC.x, y = logFC.y), 
             color = "red", size = 2) +  # Mark points for genes in the list as red
  geom_label_repel(data = subset(merged, Gene_Symbol %in% genes_to_label), 
                   aes(label = Gene_Symbol), 
                   color = "red", 
                   size = 3, 
                   box.padding = 0.5,  # Increased padding to move label farther
                   point.padding = 0.1, # Increased padding between the point and label
                   segment.color = "red",  # Line connecting label to point
                   segment.size = 0.5) +  # Set the thickness of the line
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add dashed line at y = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Add dashed line at x = 0
  labs(title = "Scatter Plot of logFC for Genes in Two Bulk Data Sets",
       x = "logFC (Bhlhe40/Bhlhe41 double knockout (DKO) microglia (Podleśny-Drabiniok A et al. Nat Commun 2024))",
       y = "logFC (Bhlhe41 knockout microglia)") +
  theme(plot.title = element_text(hjust = 0.5))
pdf("output/Scatter_plot_bulk_data.pdf", width = 8, height = 6)
print(p)
dev.off()

