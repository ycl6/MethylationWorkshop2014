#
# For the NRPB 2014 Workshop: "Hands-On Training In Methylation Sequencing Analysis"
# I-Hsuan Lin, 2014/12/19
#
options(width=130)
library(data.table)
library(reshape)
library(ggplot2)
library(scales)
library(gridExtra)

path = "../Output/"

# Load data files
genomic.h1 = data.frame(fread(paste(path,"hg18.500bins.h1.meth", sep = "")))
genomic.imr90 = data.frame(fread(paste(path,"hg18.500bins.imr90.meth", sep = "")))

# Add column names
data.genomic = cbind(genomic.h1, genomic.imr90[,4:5])
names(data.genomic) = c("chr","left","right","cg_h1","H1","cg_imr90","IMR90")

# Define the 23 chromosome to compute
chr = c(paste("chr",seq(1:22),sep = ""),"chrX")

# Remove regions that have too few CpG data and regions from non-standard chromosomes
data.genomic = data.genomic[data.genomic$cg_h1 > 9 & data.genomic$cg_imr90 > 9 & data.genomic$chr %in% chr,]

# Re-arrange chromosome order
data.genomic$chr = as.factor(data.genomic$chr)
data.genomic$chr = factor(data.genomic$chr, levels(data.genomic$chr)[c(1,12,16:22,2:11,13:15,23)])

# Re-arrange data format with melt
data.genomic = (melt(data.genomic, id = c("chr","left","right","cg_h1","cg_imr90")))

# Line plots
plot_chr = ggplot(data.genomic, aes(x = left/1000000, y = value, color = variable)) + 
facet_wrap(~ chr, ncol = 4, scales = "free_x") +
geom_line(size = 0.2) + scale_x_continuous(labels = comma) + 
guides(fill = guide_legend(title="CpG methylation", keywidth = 1, keyheight = 1, override.aes = list(size = 1))) +
theme_bw() + theme(legend.position="top",
axis.title.x = element_text(face = "bold", size = 16),
axis.title.y = element_text(face = "bold", size = 16),
axis.text.x = element_text(size = 10), 
axis.text.y = element_text(size = 10),
legend.title = element_text(face = "bold", size = 16), legend.text = element_text(size = 16),
legend.key = element_blank(), plot.title = element_text(face = "bold", size = 18),
strip.text.x = element_text(face = "bold", size = 16),
strip.background = element_blank(),
plot.title = element_text(face = "bold", size = 18), panel.margin = unit(0.2, "in")) +
labs(title = "Chromosomal Methylation Levels", x = "Positions (Mbp)", y = "Methylation")

pdf("Figure1.pdf", width=15, height=10, pointsize=12)
plot_chr
dev.off()
