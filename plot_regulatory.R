#
# For the NRPB 2014 Workshop: "Hands-On Training In Methylation Sequencing Analysis"
# I-Hsuan Lin, 2014/12/19
#
options(width=130, digits=4)
library(data.table)
library(reshape)
library(ggplot2)
library(scales)
library(gridExtra)

path = "../Output/"

# Load data files
cgi.h1 = data.frame(fread(paste(path,"cpgIslandExt.h1.meth", sep = "")))
cgi.imr90 = data.frame(fread(paste(path,"cpgIslandExt.imr90.meth", sep = "")))

tfbs.h1 = data.frame(fread(paste(path,"wgEncodeRegTfbsClustered.h1.meth", sep = "")))
tfbs.imr90 = data.frame(fread(paste(path,"wgEncodeRegTfbsClustered.imr90.meth", sep = "")))

rmsk.h1 = data.frame(fread(paste(path,"hg18.rmskRM327.h1.meth", sep = "")))
rmsk.imr90 = data.frame(fread(paste(path,"hg18.rmskRM327.imr90.meth", sep = "")))

# Add column names
data.cgi = cbind(cgi.h1, cgi.imr90[,4:5])
names(data.cgi) = c("chr","left","right","cg_h1","meth_h1","cg_imr90","meth_imr90")

data.tfbs = cbind(tfbs.h1, tfbs.imr90[,5:6])
names(data.tfbs) = c("chr","left","right","tf","cg_h1","meth_h1","cg_imr90","meth_imr90")

data.rmsk = cbind(rmsk.h1, rmsk.imr90[,5:6])
names(data.rmsk) = c("chr","left","right","repClass","cg_h1","meth_h1","cg_imr90","meth_imr90")

# Define the 23 chromosome to compute
chr = c(paste("chr",seq(1:22),sep = ""),"chrX")

# Remove regions that have too few CpG data and regions from non-standard chromosomes
data.cgi = data.cgi[data.cgi$cg_h1 > 9 & data.cgi$cg_imr90 > 9 & data.cgi$chr %in% chr,]
data.tfbs = data.tfbs[data.tfbs$cg_h1 > 9 & data.tfbs$cg_imr90 > 9 & data.tfbs$chr %in% chr,]
data.rmsk = data.rmsk[data.rmsk$cg_h1 > 9 & data.rmsk$cg_imr90 > 9 & data.rmsk$chr %in% chr,]

# Classify methylation levels into 4 categories
breaks = c(-1,0.2,0.5,0.8,1)
breaksNames = c("<20%","20%-50%","50%-80%",">80%")

data.cgi$cut_h1 = cut(data.cgi$meth_h1, breaks=breaks)
levels(data.cgi$cut_h1) = breaksNames
data.cgi$cut_imr90 = cut(data.cgi$meth_imr90, breaks=breaks)
levels(data.cgi$cut_imr90) = breaksNames

data.tfbs$cut_h1 = cut(data.tfbs$meth_h1, breaks=breaks)
levels(data.tfbs$cut_h1) = breaksNames
data.tfbs$cut_imr90 = cut(data.tfbs$meth_imr90, breaks=breaks)
levels(data.tfbs$cut_imr90) = breaksNames

data.rmsk$cut_h1 = cut(data.rmsk$meth_h1, breaks=breaks)
levels(data.rmsk$cut_h1) = breaksNames
data.rmsk$cut_imr90 = cut(data.rmsk$meth_imr90, breaks=breaks)
levels(data.rmsk$cut_imr90) = breaksNames

pt =  prop.table(table(data.cgi$cut_h1, data.cgi$cut_imr90))*100
print(pt)

#table(data.tfbs$tf)
#
#   BAF155    BAF170      BATF    BCL11A      BCL3   BHLHE40      Brg1     CEBPB     c-Fos     c-Jun     c-Myc      CTCF       EBF
#    13814      8316      1476      1783      2992       374      3221      1102      2259      1689      6786      5215      6105
#    Egr-1      ERRA     FOSL2     FOXP2      GABP        GR     GRp20      HEY1     HNF4A      HSF1      Ini1     Input      IRF4
#     1559       905      3612       972      2603      2574       315     10768      1577       594      6376       286      3580
#     JunD       Max     NF-E2      NFKB      Nrf1      NRSF      p300  PAX5-C20  PAX5-N19      Pbx3     PGC1A      Pol3    POU2F2
#     3592      7743       218      2523      1671      1931      2970      5565      3003      2217       150       106      7231
#     PU.1     Rad21      RXRA Sin3Ak-20      SIX5       SP1    SREBP1    SREBP2       SRF     STAT1     STAT2      TAF1     TCF12
#     3989      1639      1632      5625       874      4668       475         1       245       762       282      9101      4529
#    USF-1     XRCC4    ZBTB33
#     3214         1       585

# Remove SREBP2 and XRCC4 that has too few sites, and re-arrange data format with melt
data1 = melt(data.tfbs[!data.tfbs$tf %in% c("SREBP2","XRCC4"),c(4,6,8)], id = "tf")
levels(data1$variable) = c("H1","IMR90")

data2 = melt(data.rmsk[,c(4,6,8)], id = "repClass")
levels(data2$variable) = c("H1","IMR90")

# Box plots
plot1 = ggplot(data1, aes(x = tf, y = value, fill = variable)) +
geom_boxplot(outlier.shape = NA, size = 0.5, color = "black") + theme_bw() +
guides(fill = guide_legend(keywidth = 2, keyheight = 2)) +
theme(legend.position="top",
axis.title.x = element_blank(), 
axis.title.y = element_text(face = "bold", size = 16),
axis.text.x = element_text(size = 10, angle = 35, hjust=1, vjust = 1), 
axis.text.y = element_text(size = 16),
legend.title = element_blank(), legend.text = element_text(size = 16),
legend.key = element_blank(), plot.title = element_text(face = "bold", size = 18),
plot.title = element_text(face = "bold", size = 18), panel.margin = unit(0.2, "in")) +
labs(title = "Methylation Levels at TFBS", y = "Methylation")

plot2 = ggplot(data2, aes(x = repClass, y = value, fill = variable)) +
geom_boxplot(outlier.shape = NA, size = 0.5, color = "black") + theme_bw() +
guides(fill = guide_legend(keywidth = 2, keyheight = 2)) +
theme(legend.position="top",
axis.title.x = element_blank(), 
axis.title.y = element_text(face = "bold", size = 16),
axis.text.x = element_text(size = 16), 
axis.text.y = element_text(size = 16),
legend.title = element_blank(), legend.text = element_text(size = 16),
legend.key = element_blank(), plot.title = element_text(face = "bold", size = 18),
plot.title = element_text(face = "bold", size = 18), panel.margin = unit(0.2, "in")) +
labs(title = "Methylation Levels at Repeating Elements", y = "Methylation") +
annotate("text", x = 1, y = 0.05, size = 5, label = paste("n =", table(data2$repClass)[1]/2)) +
annotate("text", x = 2, y = 0.05, size = 5, label = paste("n =", table(data2$repClass)[2]/2)) +
annotate("text", x = 3, y = 0.05, size = 5, label = paste("n =", table(data2$repClass)[3]/2)) +
annotate("text", x = 4, y = 0.05, size = 5, label = paste("n =", table(data2$repClass)[4]/2)) +
annotate("text", x = 5, y = 0.05, size = 5, label = paste("n =", table(data2$repClass)[5]/2))

pdf("Figure3.pdf", width=18, height=6, pointsize=12)
plot1
plot2
dev.off()
