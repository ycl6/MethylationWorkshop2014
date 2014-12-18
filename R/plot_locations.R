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
exon.h1 = data.frame(fread(paste(path,"gencode.v3c.exon_merged.h1.meth", sep = "")))
exon.imr90 = data.frame(fread(paste(path,"gencode.v3c.exon_merged.imr90.meth", sep = "")))

intergenic.h1 = data.frame(fread(paste(path,"gencode.v3c.intergenic.h1.meth", sep = "")))
intergenic.imr90 = data.frame(fread(paste(path,"gencode.v3c.intergenic.imr90.meth", sep = "")))

intron.h1 = data.frame(fread(paste(path,"gencode.v3c.intron_merged.h1.meth", sep = "")))
intron.imr90 = data.frame(fread(paste(path,"gencode.v3c.intron_merged.imr90.meth", sep = "")))

promoter.h1 = data.frame(fread(paste(path,"gencode.v3c.promoter_merged.h1.meth", sep = "")))
promoter.imr90 = data.frame(fread(paste(path,"gencode.v3c.promoter_merged.imr90.meth", sep = "")))

utr.h1 = data.frame(fread(paste(path,"gencode.v3c.UTR.h1.meth", sep = "")))
utr.imr90 = data.frame(fread(paste(path,"gencode.v3c.UTR.imr90.meth", sep = "")))

# Add column names
data.exon = cbind(exon.h1, exon.imr90[,4:5])
names(data.exon) = c("chr","left","right","cg_h1","meth_h1","cg_imr90","meth_imr90")

data.intergenic = cbind(intergenic.h1, intergenic.imr90[,4:5])
names(data.intergenic) = c("chr","left","right","cg_h1","meth_h1","cg_imr90","meth_imr90")

data.intron = cbind(intron.h1, intron.imr90[,4:5])
names(data.intron) = c("chr","left","right","cg_h1","meth_h1","cg_imr90","meth_imr90")

data.promoter = cbind(promoter.h1, promoter.imr90[,4:5])
names(data.promoter) = c("chr","left","right","cg_h1","meth_h1","cg_imr90","meth_imr90")

data.utr = cbind(utr.h1, utr.imr90[,7:8])
names(data.utr) = c("chr","left","right","utr","id","strand","cg_h1","meth_h1","cg_imr90","meth_imr90")

# Define the 23 chromosome to compute
chr = c(paste("chr",seq(1:22),sep = ""),"chrX")

# Remove regions that have too few CpG data and regions from non-standard chromosomes
data.exon = data.exon[data.exon$cg_h1 > 9 & data.exon$cg_imr90 > 9 & data.exon$chr %in% chr,]
data.intergenic = data.intergenic[data.intergenic$cg_h1 > 9 & data.intergenic$cg_imr90 > 9 & data.intergenic$chr %in% chr,]
data.intron = data.intron[data.intron$cg_h1 > 9 & data.intron$cg_imr90 > 9 & data.intron$chr %in% chr,]
data.promoter = data.promoter[data.promoter$cg_h1 > 9 & data.promoter$cg_imr90 > 9 & data.promoter$chr %in% chr,]
data.utr = data.utr[data.utr$cg_h1 > 9 & data.utr$cg_imr90 > 9 & data.utr$chr %in% chr,]

# Define chromosome type, chr1 to chr22 are autosomes
data.exon$chrType = "autosome"
data.exon[data.exon$chr == "chrX",]$chrType = "chrX"
data.intergenic$chrType = "autosome"
data.intergenic[data.intergenic$chr == "chrX",]$chrType = "chrX"
data.intron$chrType = "autosome"
data.intron[data.intron$chr == "chrX",]$chrType = "chrX"
data.promoter$chrType = "autosome"
data.promoter[data.promoter$chr == "chrX",]$chrType = "chrX"
data.utr$chrType = "autosome"
data.utr[data.utr$chr == "chrX",]$chrType = "chrX"

# Classify methylation levels into 4 categories
breaks = c(-1,0.2,0.5,0.8,1)
breaksNames = c("<20%","20%-50%","50%-80%",">80%")

data.exon$cut_h1 = cut(data.exon$meth_h1, breaks=breaks)
levels(data.exon$cut_h1) = breaksNames
data.exon$cut_imr90 = cut(data.exon$meth_imr90, breaks=breaks)
levels(data.exon$cut_imr90) = breaksNames
data.intergenic$cut_h1 = cut(data.intergenic$meth_h1, breaks=breaks)
levels(data.intergenic$cut_h1) = breaksNames
data.intergenic$cut_imr90 = cut(data.intergenic$meth_imr90, breaks=breaks)
levels(data.intergenic$cut_imr90) = breaksNames
data.intron$cut_h1 = cut(data.intron$meth_h1, breaks=breaks)
levels(data.intron$cut_h1) = breaksNames
data.intron$cut_imr90 = cut(data.intron$meth_imr90, breaks=breaks)
levels(data.intron$cut_imr90) = breaksNames
data.promoter$cut_h1 = cut(data.promoter$meth_h1, breaks=breaks)
levels(data.promoter$cut_h1) = breaksNames
data.promoter$cut_imr90 = cut(data.promoter$meth_imr90, breaks=breaks)
levels(data.promoter$cut_imr90) = breaksNames
data.utr$cut_h1 = cut(data.utr$meth_h1, breaks=breaks)
levels(data.utr$cut_h1) = breaksNames
data.utr$cut_imr90 = cut(data.utr$meth_imr90, breaks=breaks)

# Join various regions
data = rbind(data.exon[,8:10], data.intergenic[,8:10], data.intron[,8:10], data.promoter[,8:10], data.utr[,11:13])
data$location = c(rep("exon",dim(data.exon)[1]), rep("intergenic",dim(data.intergenic)[1]), rep("intron",dim(data.intron)[1]), rep("promoter",dim(data.promoter)[1]), rep("utr",dim(data.utr)[1]))
data = cbind(as.data.frame(prop.table(table(data$location, data$cut_h1),1)), as.data.frame(prop.table(table(data$location, data$cut_imr90),1))[,3])
names(data) = c("location","methylation","H1","IMR90")

# Re-arrange data format with melt
data = melt(data, id = c("location","methylation"))

# Bar plots
plot = ggplot(data, aes(x = location, y = value, fill = variable)) + facet_wrap(~ methylation, nrow = 1) +
geom_bar(stat = "identity", width=.5, position = "dodge") + theme_bw() +
guides(fill = guide_legend(title="CpG methylation", keywidth = 1, keyheight = 1, override.aes = list(size = 1))) +
theme(legend.position="top",
axis.title.x = element_blank(),
axis.title.y = element_text(face = "bold", size = 16),
axis.text.x = element_text(size = 14, angle = 35, hjust=1, vjust = 1),
axis.text.y = element_text(size = 14),
legend.title = element_text(face = "bold", size = 16), legend.text = element_text(size = 16),
legend.key = element_blank(), plot.title = element_text(face = "bold", size = 18),
strip.text.x = element_text(face = "bold", size = 16),
strip.background = element_blank(),
plot.title = element_text(face = "bold", size = 18), panel.margin = unit(0.2, "in")) +
labs(title = "Overview of Methylation Levels", y = "Proportion")

# Scattered plots
plot.exon = ggplot(data.exon, aes(x = meth_h1, y = meth_imr90, color = chrType)) +
geom_point(data = data.exon[data.exon$chrType == "autosome",], size = 0.6, alpha = 0.5) +
geom_point(data = data.exon[data.exon$chrType == "chrX",], size = 1.2) +
scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,1)) +
facet_wrap(~ chrType) + theme_bw() + geom_abline(intercept = 0, slope = 1, size = 1) +
theme(legend.position = "none",
axis.title.x = element_text(face = "bold", size = 16),
axis.title.y = element_text(face = "bold", size = 16),
axis.text.x = element_text(size = 14),
axis.text.y = element_text(size = 14),
strip.text.x = element_text(face = "bold", size = 16),
strip.background = element_blank(),
plot.title = element_text(face = "bold", size = 18), panel.margin = unit(0.2, "in")) +
labs(title = "Methylation at Exons", x = "H1",y = "IMR90")

plot.intron = ggplot(data.intron, aes(x = meth_h1, y = meth_imr90, color = chrType)) +
geom_point(data = data.intron[data.intron$chrType == "autosome",], size = 0.6, alpha = 0.5) +
geom_point(data = data.intron[data.intron$chrType == "chrX",], size = 1.2) +
scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,1)) +
facet_wrap(~ chrType) + theme_bw() + geom_abline(intercept = 0, slope = 1, size = 1) +
theme(legend.position = "none",
axis.title.x = element_text(face = "bold", size = 16),
axis.title.y = element_text(face = "bold", size = 16),
axis.text.x = element_text(size = 14),
axis.text.y = element_text(size = 14),
strip.text.x = element_text(face = "bold", size = 16),
strip.background = element_blank(),
plot.title = element_text(face = "bold", size = 18), panel.margin = unit(0.2, "in")) +
labs(title = "Methylation at Introns", x = "H1",y = "IMR90")

plot.promoter = ggplot(data.promoter, aes(x = meth_h1, y = meth_imr90, color = chrType)) +
geom_point(data = data.promoter[data.promoter$chrType == "autosome",], size = 0.6, alpha = 0.5) +
geom_point(data = data.promoter[data.promoter$chrType == "chrX",], size = 1.2) +
scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,1)) +
facet_wrap(~ chrType) + theme_bw() + geom_abline(intercept = 0, slope = 1, size = 1) +
theme(legend.position = "none",
axis.title.x = element_text(face = "bold", size = 16),
axis.title.y = element_text(face = "bold", size = 16),
axis.text.x = element_text(size = 14),
axis.text.y = element_text(size = 14),
strip.text.x = element_text(face = "bold", size = 16),
strip.background = element_blank(),
plot.title = element_text(face = "bold", size = 18), panel.margin = unit(0.2, "in")) +
labs(title = "Methylation at Promoters", x = "H1",y = "IMR90")

plot.utr = ggplot(data.utr, aes(x = meth_h1, y = meth_imr90, color = chrType)) +
geom_point(data = data.utr[data.utr$chrType == "autosome",], size = 0.6, alpha = 0.5) +
geom_point(data = data.utr[data.utr$chrType == "chrX",], size = 1.2) +
scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,1)) +
facet_grid(utr ~ chrType) + theme_bw() + geom_abline(intercept = 0, slope = 1, size = 1) +
theme(legend.position = "none",
axis.title.x = element_text(face = "bold", size = 16),
axis.title.y = element_text(face = "bold", size = 16),
axis.text.x = element_text(size = 14),
axis.text.y = element_text(size = 14),
strip.text.x = element_text(face = "bold", size = 16),
strip.text.y = element_text(face = "bold", size = 16),
strip.background = element_blank(),
plot.title = element_text(face = "bold", size = 18), panel.margin = unit(0.2, "in")) +
labs(title = "Methylation at UTRs", x = "H1",y = "IMR90")

plot.intergenic = ggplot(data.intergenic, aes(x = meth_h1, y = meth_imr90, color = chrType)) +
geom_point(data = data.intergenic[data.intergenic$chrType == "autosome",], size = 0.6, alpha = 0.5) +
geom_point(data = data.intergenic[data.intergenic$chrType == "chrX",], size = 1.2) +
scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,1)) +
facet_wrap(~ chrType) + theme_bw() + geom_abline(intercept = 0, slope = 1, size = 1) +
theme(legend.position = "none",
axis.title.x = element_text(face = "bold", size = 16),
axis.title.y = element_text(face = "bold", size = 16),
axis.text.x = element_text(size = 14),
axis.text.y = element_text(size = 14),
strip.text.x = element_text(face = "bold", size = 16),
strip.background = element_blank(),
plot.title = element_text(face = "bold", size = 18), panel.margin = unit(0.2, "in")) +
labs(title = "Methylation at Intergenic Regions", x = "H1",y = "IMR90")

pdf("Figure2a.pdf", width=10, height=5, pointsize=12)
plot
dev.off()

pdf("Figure2b.pdf", width=20, height=10, pointsize=12)
grid.arrange(plot.exon, plot.intron, plot.promoter, plot.intergenic, ncol = 2)
dev.off()

pdf("Figure2c.pdf", width=10, height=10, pointsize=12)
plot.utr
dev.off()
