#
# For the NRPB 2014 Workshop: "Hands-On Training In Methylation Sequencing Analysis"
# I-Hsuan Lin, 2014/12/19
#
options(width=130)
library(data.table)
library(ggplot2)
library(scales)
library(gridExtra)

path = "../Output/"

# Load data files
hmr = data.frame(fread(paste(path,"h1_imr90_hmr_coverage.bed", sep = "")))
pmd = data.frame(fread(paste(path,"h1_imr90_pmd_coverage.bed", sep = "")))

# Add column names
names(hmr) = c("chr","left","right","H1_length","H1_coverage","IMR90_length","IMR90_coverage")
names(pmd) = c("chr","left","right","H1_length","H1_coverage","IMR90_length","IMR90_coverage")

# Define the 23 chromosome to compute
chr = c(paste("chr",seq(1:22),sep = ""),"chrX")

# Remove regions from non-standard chromosomes
hmr = hmr[hmr$chr %in% chr,]
pmd = pmd[pmd$chr %in% chr,]

# Re-arrange chromosome order
pmd$chr = as.factor(pmd$chr)
pmd$chr = factor(pmd$chr, levels(pmd$chr)[c(1,12,16:22,2:11,13:15,23)])

# Calculate the merged HMR and PMD lengths
hmr$merged_length = hmr$right - hmr$left
pmd$merged_length = pmd$right - pmd$left

# Calculate the overlapped differences and define HMR conservation levels
hmr$diff = hmr$H1_coverage - hmr$IMR90_coverage
hmr$cut_diff = cut(hmr$diff, breaks=c(-1.1,-0.9999,-0.7,-0.2,0.2,0.7,0.9999,1))
levels(hmr$cut_diff) = c("IMR90-specific","IMR90 >> H1","IMR90 > H1","IMR90 ~~ H1","IMR90 < H1","IMR90 << H1","H1-specific")

# Calculate the overlapped differences and define PMD conservation levels
pmd$diff = pmd$H1_coverage - pmd$IMR90_coverage
pmd$cut_diff = cut(pmd$diff, breaks=c(-1.1,-0.9999,0,0.9999,1))
levels(pmd$cut_diff) = c("IMR90-specific","IMR90 > H1","IMR90 < H1","H1-specific")

hmr_length = data.frame(
	name = c( rep("H1", length(hmr[hmr$H1_length != 0,]$H1_length)), rep("IMR90", length(hmr[hmr$IMR90_length != 0,]$IMR90_length))),
	length = c( hmr[hmr$H1_length != 0,]$H1_length, hmr[hmr$IMR90_length != 0,]$IMR90_length)
)

pmd_length = data.frame(
	name = c( rep("H1", length(pmd[pmd$H1_length != 0,]$H1_length)), rep("IMR90", length(pmd[pmd$IMR90_length != 0,]$IMR90_length))),
	length = c( pmd[pmd$H1_length != 0,]$H1_length, pmd[pmd$IMR90_length != 0,]$IMR90_length)
)

# Display HMR and PMD size distribution using box plots
plot_hmr = ggplot(hmr_length, aes(x = name, y = length, fill = name)) +
geom_boxplot(outlier.shape = 4, size = 0.5, color = "black") + theme_bw() +
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels=trans_format("log10",math_format(10^.x))) +
guides(fill = guide_legend(keywidth = 2, keyheight = 2)) +
theme(legend.position="none",
axis.title.x = element_blank(), 
axis.title.y = element_text(face = "bold", size = 16),
axis.text.x = element_text(size = 14), 
axis.text.y = element_text(size = 14),
plot.title = element_text(face = "bold", size = 18), panel.margin = unit(0.2, "in")) +
labs(title = "HMR Lengths", y = "Length (bp)") +
annotate("text", x = 1, y = 10, size = 4, label = paste("n =", dim(hmr[hmr$H1_length != 0,])[1])) +
annotate("text", x = 2, y = 10, size = 4, label = paste("n =", dim(hmr[hmr$IMR90_length != 0,])[1]))

plot_pmd = ggplot(pmd_length, aes(x = name, y = length, fill = name)) +
geom_boxplot(outlier.shape = 4, size = 0.5, color = "black") + theme_bw() +
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels=trans_format("log10",math_format(10^.x))) +
guides(fill = guide_legend(keywidth = 2, keyheight = 2)) +
theme(legend.position="none",
axis.title.x = element_blank(), 
axis.title.y = element_text(face = "bold", size = 16),
axis.text.x = element_text(size = 14), 
axis.text.y = element_text(size = 14),
plot.title = element_text(face = "bold", size = 18), panel.margin = unit(0.2, "in")) +
labs(title = "PMD Lengths", y = "Length (bp)") +
annotate("text", x = 1, y = 1000, size = 4, label = paste("n =", dim(pmd[pmd$H1_length != 0,])[1])) +
annotate("text", x = 2, y = 1000, size = 4, label = paste("n =", dim(pmd[pmd$IMR90_length != 0,])[1]))

# Display conservation levels between H1 and IMR90 using box plots
plot_hmr_lengths = ggplot(hmr, aes(x = cut_diff, y = merged_length)) + 
geom_boxplot(outlier.shape = 4, size = 0.5, color = "black") + theme_bw() +
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels=trans_format("log10",math_format(10^.x))) +
theme(legend.position="none",
axis.title.x = element_blank(),
axis.title.y = element_text(face = "bold", size = 16),
axis.text.x = element_text(size = 14, angle = 35, hjust=1, vjust = 1), 
axis.text.y = element_text(size = 14),
plot.title = element_text(face = "bold", size = 16), panel.margin = unit(0.2, "in")) +
labs(title = "Size of Merged HMRs of\nVarious Conservation Levels", y = "Length (bp)")

plot_pmd_lengths = ggplot(pmd, aes(x = chr, y = merged_length/1000000)) + 
geom_boxplot(outlier.shape = 4, size = 0.1, color = "black") + theme_bw() +
scale_y_continuous(labels = comma) +
facet_wrap(~ cut_diff, scales = "free_y", ncol = 1) + theme(legend.position="none",
axis.title.x = element_blank(),
axis.title.y = element_text(face = "bold", size = 14),
axis.text.x = element_text(size = 8, angle = 35, hjust=1, vjust = 1), 
axis.text.y = element_text(size = 8),
strip.text.x = element_text(face = "bold", size = 10),
strip.background = element_blank(),
plot.title = element_text(face = "bold", size = 16), panel.margin = unit(0.2, "in")) +
labs(title = "Size of Merged PMDs of\nVarious Conservation Levels", y = "Length (Mbp)")

pdf("Figure4.pdf", width=12, height=6, pointsize=12)
grid.arrange(plot_hmr, plot_pmd, nrow = 1)
grid.arrange(plot_hmr_lengths, plot_pmd_lengths, nrow = 1)
dev.off()
