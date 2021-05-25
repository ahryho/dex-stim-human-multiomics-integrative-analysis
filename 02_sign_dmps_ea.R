library(data.table)
library(ggplot2)
library(RColorBrewer)

methyl.src.data.pre <- "~/bio/datasets/methylation/20_DMA/02_dmp/dmps_significant_annotated.txt"

dmps.anno.fn      <- paste0(methyl.src.data.pre, "02_dmp/dmp_bcc_pcs_anno.csv")
dmps.anno.df      <- fread(dmps.anno.fn, sep = "\t")

dmps.sign.anno.fn <- paste0(methyl.src.data.pre, "02_dmp/dmps_significant_annotated.txt")
dmps.anno.sign.df <- fread(dmps.sign.anno.fn)

str(dmps.anno.sign.df)

dmp.anno.subdf <- dmps.anno.sign.df[,c("Name","chr","pos", "UCSC_RefGene_Group", "UCSC_RefGene_Name", "Regulatory_Feature_Group", "Relation_to_Island")]
methyl_anno <- dmp.anno.subdf

ggplot(methyl_anno, aes(x=chr)) + geom_bar()
# replace all empty cells with NA
ggplot(methyl_anno, aes(x=Relation_to_Island)) + geom_bar()

island_association <- ggplot(methyl_anno, aes(Relation_to_Island)) +
  geom_bar(aes(fill=Regulatory_Feature_Group)) +
  labs(title="Island and Association") +
  scale_fill_brewer(palette="Paired", na.value="grey") + ylab("Count") +
  theme_classic()

#############################
# DENSITY
meth <- methyl_anno
# put the chromosomes in the good order: chr1, chr2, chr22, chrX
goodChrOrder <- paste("chr",c(1:22,"X","Y"),sep="")
meth$chr <- factor(meth$chr,levels=goodChrOrder)

# Plot the densities of snps in the bed file for each chr seperately
meth_Density<-ggplot(meth) + 
  geom_histogram(aes(x=pos),binwidth=1e6) + # pick a binwidth that is not too small 
  facet_wrap(~ chr,ncol=2) + # seperate plots for each chr, x-scales can differ from chr to chr
  ggtitle("Density of methylated CpG across hg19 - autosomes") +
  xlab("Position in the genome") + 
  ylab("SNP density") + 
  theme_bw() # I prefer the black and white theme


#############################

# # plot all annotations
# all_anno <- unlist(strsplit(as.character(methyl_anno$UCSC_RefGene_Group), split = "\\;" ))
# dt <- as.data.table(all_anno, 1:length(all_anno))
# ggplot(dt, aes(all_anno)) + geom_bar() +
#   labs(title="Annotations of cg sites - multiple annotations possible",
#        x="Annotations", y="Count") + theme_classic() + Emy_Style
# length(all_anno)
source("./exploratory_analysis/0_exploratory_analysis_functions.R")

# only the first annotation
one_anno_df <- methyl_anno
one_anno_df$UCSC_RefGene_Group <- gsub(";.*", "", one_anno_df$UCSC_RefGene_Group)
colourCount = length(unique(one_anno_df$Regulatory_Feature_Group))
getPalette = colorRampPalette(brewer.pal(5, "Set2"))

p1 <- ggplot(one_anno_df, aes(UCSC_RefGene_Group)) +
  geom_bar(aes(fill=Regulatory_Feature_Group)) +
  ylim(0,3e+05) +
  labs(title="Annotation and Association\none annotation") +
  scale_fill_manual(values=getPalette(colourCount), na.value="grey") + ylab("Count") +
  theme_classic() + Emy_Style


# Regulatory Feature Group vs Relation to Island
p2 <- ggplot(methyl_anno, aes(Regulatory_Feature_Group)) +
  geom_bar(aes(fill=Relation_to_Island)) +
  scale_fill_brewer(palette="RdYlBu") +
  labs(y="Count", title="Association and Island") +
  theme_classic() + Emy_Style


###############################################
mylist <- strsplit(as.character(methyl_anno$UCSC_RefGene_Group), split = "\\;" )
names(mylist) <- methyl_anno$cgSite
g <- plyr::ldply(mylist, rbind)

m_g[19727210:19727227,]
m_g<- melt(g, id = ".id") %>% as.data.table
first_pos <- m_g[variable ==1 ,]
other_pos <- m_g[variable !=1 ,]
other_pos <- other_pos[complete.cases(other_pos),]
m_g <- rbind(first_pos,other_pos)
dim(m_g)
table(m_g$variable)

m_g_m <- merge(as.data.table(m_g), methyl_anno[,.(cgSite, Regulatory_Feature_Group)],
               by.x=".id", by.y="cgSite")


p3 <- ggplot(m_g_m, aes(x=value, fill=Regulatory_Feature_Group)) + 
  geom_bar(aes(fill=Regulatory_Feature_Group)) +
  ylim(0,3e+05) +
  labs(title="Annotation and Association\nmultiple annotations",
       y="Count",
       x="UCSC_RefGene_Group") +
  scale_fill_manual(values=getPalette(colourCount), na.value="grey") +
  theme_classic() + Emy_Style

###############################################

sub <- m_g_m[value == "TSS1500" | value =="TSS200",]
p4 <- ggplot(sub, aes(x=value, fill=Regulatory_Feature_Group)) +
  geom_bar(aes(fill=Regulatory_Feature_Group)) +
  scale_fill_manual(values=getPalette(colourCount), na.value="grey") +
  labs(title="Transcription Start Site and Association",
       y="Count",
       x="Transcription Start Site") +
  theme_classic() + Emy_Style

library(cowplot)
mylegend <- get_legend(p1)

pdf("./exploratory_analysis/cgSites_Annotation.pdf", width=10, height=10)
meth_Density
plotarrange <- plot_grid(p1 + theme(legend.position = "none"), 
                         p3 + theme(legend.position = "none"))
plot_grid(plotarrange, mylegend, rel_widths = c(2, .5))
p4;p2
island_association
dev.off()

###############################################
# number of genes and cg Sites associated with promoter and TS200/TS1500

dt_TSS200_Promoter <- methyl_anno[grepl('TSS200',methyl_anno$UCSC_RefGene_Group)][Regulatory_Feature_Group == "Promoter_Associated"]
length(unique(dt_TSS200_Promoter$cgSite))
length(unique(dt_TSS200_Promoter$gene))

dt_TSS1500_Promoter <- methyl_anno[grepl('TSS1500',methyl_anno$UCSC_RefGene_Group)][Regulatory_Feature_Group == "Promoter_Associated"]
length(unique(dt_TSS1500_Promoter$cgSite))
length(unique(dt_TSS1500_Promoter$gene))

dt_TSS_Promoter <- methyl_anno[grepl('TSS1500',methyl_anno$UCSC_RefGene_Group) | grepl('TSS200',methyl_anno$UCSC_RefGene_Group),] %>% 
  .[Regulatory_Feature_Group == "Promoter_Associated"]
length(unique(dt_TSS_Promoter$cgSite))
length(unique(dt_TSS_Promoter$gene))


##############################################
# plot the distances from gene
ggplot(db_gene_cg[distance < 100,], aes(x=distance)) + geom_histogram(bins=10)
ggplot(db_gene_cg[distance < 100000,], aes(x=distance)) + geom_histogram(bins=10)
ggplot(db_gene_cg[distance < 100,], aes(x=distance)) + geom_histogram(bins=10)

ggplot(db_gene_cg, aes(y=distance)) + geom_boxplot()
summary(db_gene_cg[, distance])



####################################################
####################################################
####################################################
numgenes <- db_gene_cg[, .N, by=gene]
ggplot(numgenes, aes(x=N)) + geom_bar() +
  labs(x="Number of methylation sites", title="Number of Methylation Sites per Gene") 

db_gene_cg[, .N, by=distance][order(distance)]
distance <- db_gene_cg[distance < 100,]
ggplot(distance, aes(x=distance)) + geom_bar() + labs(x="distance to gene")

####################################################
db_ilmn_gene <- na.omit(db_ilmn_gene)
db_ilmn_gene[, .N, by=Gene][order(N)]

numgenes <- db_ilmn_gene[, .N, by=Gene]
ggplot(numgenes, aes(x=N)) + geom_bar() +
  labs(x="Number of ILMN", title="Number of ILMN identifier per Gene") 
####################################################
####################################################
####################################################

# compare the methylanno results to Janine's reannotation
db_gene_cg <- read.table("../Dex-eQTM/ReAnnotation/450K_Janine_annotation_withoutBadProbes_Final_and_rs_ch.txt",
                         header=T) %>% setDT
db_ilmn_gene <- read.table("../Dex-eQTM/ReAnnotation/v3_v4_sharedContent_QC.txt",
                           header=T) %>% setDT
str(db_gene_cg)

plot(db_gene_cg$Chr)


#############################
# DENSITY
meth_janine <- db_gene_cg
meth_janine$Chr <- as.factor(meth_janine$Chr)
# put the chromosomes in the good order: chr1, chr2, chr22, chrX
goodChrOrder <- paste(c(1:22),sep="")
meth_janine$Chr <- factor(meth_janine$Chr,levels=goodChrOrder)


# compare to the other one
load("Methyl_annno_illumina.rda")
meth_illumina <- methyl_anno
meth_illumina$chr <- as.factor(gsub("chr", "", meth_illumina$chr))
m_ill <- meth_illumina[,.(cgSite, chr, pos)]
m_ill <- m_ill[, db := "Illumina"]

m_j <- meth_janine[, .(Probe_Id, Chr, P_start)]
m_j <- m_j[, db := "Janine"]
colnames(m_j) <- c("cgSite", "chr", "pos", "db")

coll_m <- rbind(m_j, m_ill)


# Plot the densities of snps in the bed file for each chr seperately
meth_Density<-ggplot(coll_m) + 
  geom_histogram(aes(x=pos, fill=db),binwidth=1e6) + # pick a binwidth that is not too small 
  facet_wrap(~ chr,ncol=2) + # seperate plots for each chr, x-scales can differ from chr to chr
  ggtitle("Density of methylated CpG across hg19 - autosomes") +
  xlab("Position in the genome") + 
  ylab("SNP density") + 
  theme_bw() # I prefer the black and white theme

pdf("Methylation_density_comparison.pdf", height=10, width=10)
meth_Density
dev.off()
