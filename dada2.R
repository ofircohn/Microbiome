#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2")
  rm(list = ls())
  cat("\f")  
library(dada2)
library(Biostrings)
library(ggplot2)
library(stringr)
library(DECIPHER)
library(readr)
library(ShortRead)

raw="C:/Users/Ofir Cohn/Projects/Microbiome/baseline_CC/Raw_data/" #path to files cleaned of primers, reverse primers, and LSU/SSU sequences
dir="C:/Users/Ofir Cohn/Projects/Microbiome/baseline_CC/" #folder am working in
#dir.create(dir,"dada2/")
#fwd_primer <- "GTGCCAGCMGCCGCGGTAA"
#rev_primer <- "GGACTACHVGGGTWTCTAAT"
#fwd_primer_rev <- as.character(reverseComplement(DNAStringSet(fwd_primer)))
#rev_primer_rev <- as.character(reverseComplement(DNAStringSet(rev_primer)))
# count_primers <- function(primer, filename) {
#   num_hits <- vcountPattern(primer, sread(readFastq(filename)), fixed = FALSE)
#   return(sum(num_hits > 0))
# }
# fwd_files <- sort(list.files(raw, pattern = "R1", full.names = TRUE)) 
# rev_files <- sort(list.files(raw, pattern = "R2", full.names = TRUE))
# count_primers(fwd_primer, fwd_files[[5]])
# count_primers(rev_primer, rev_files[[5]])
names_F=list.files(raw,"_R1")
names_R=list.files(raw,"_R2")
fnFs=paste0(raw,names_F)
fnRs=paste0(raw,names_R)
sample.names = sub("_", "-", names_F) # replace _ with -
filtFs <- paste0(dir,"dada2/",sample.names,"_R1.fastq.gz") #make vector of names for filrtered files. add ".gz" if need to compress. 
filtRs <- paste0(dir,"dada2/",sample.names,"_R2.fastq.gz") #make vector of names for filrtered files. add ".gz" if need to compress.

# quality
plotQualityProfile(fnFs[5])
plotQualityProfile(fnRs[5])
# filter
out <- filterAndTrim(fnFs, filtFs,fnRs,filtRs, maxN = 0, maxEE=c(2,2),
                     truncQ = 2,truncLen=c(240,200), trimLeft=c(19, 20) ,rm.phix = TRUE, 
                     compress =TRUE, verbose=TRUE)
out=as.data.frame(out) #format out for easy viewing
out$per.survived<-out$reads.out/out$reads.in
out$SampID=gsub(".fastq.gz","",rownames(out))
write.table(out,paste0(dir,"summary_dadafilt.txt"), row.names=F,sep="\t",quote=F)

# LearnErrors
errF=learnErrors(filtFs,multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
dev.off()
errR=learnErrors(filtRs,multithread=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()
print(dada2:::checkConvergence(errF))
print(dada2:::checkConvergence(errR))

# Dereplicate
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- names_F
names(derepRs) <- names_R

# Infer sample composition
dadaFs=dada(derepFs,err=errF, multithread = TRUE, pool="pseudo")
dadaRs=dada(derepRs,err=errR, multithread = TRUE, pool="pseudo")
print(dadaFs)
print(dadaRs)

# Merge forward/reverse reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])
# Remove chimeras
seqtab <- makeSequenceTable(mergers) #OTU matrix
seqtab.nochim <- removeBimeraDenovo(seqtab, multithread=TRUE, verbose=TRUE)
print(sort(rowSums(seqtab.nochim)/rowSums(seqtab)))
sum(seqtab.nochim)/sum(seqtab) #in this case we barely lost any in terms of abundance
# Overview of counts throughout pipline
getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(row.names=sample.names, dada2_input=out[,1],
                          filtered=out[,2], dada_f=sapply(dadaFs, getN),
                          dada_r=sapply(dadaRs, getN), merged=sapply(mergers, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/out[,1]*100, 1))
write.table(summary_tab,paste0(dir, "read-count-tracking.tsv"), quote=FALSE, sep="\t", col.names=NA)

# extraction
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta,paste0(dir, "ASVs.fa"))
# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, paste0(dir,"ASVs_counts.tsv"), sep="\t", quote=F, col.names=NA)
# Assign taxa
taxa.ref<-"C:/Users/Ofir Cohn/Projects/Microbiome/baseline_CC/silva_nr99_v138.1_train_set.fa.gz"
species.ref<-"C:/Users/Ofir Cohn/Projects/Microbiome/baseline_CC/silva_species_assignment_v138.1.fa.gz"
taxa <- assignTaxonomy(seqtab.nochim, taxa.ref, multithread=TRUE)
taxa <- addSpecies(taxa, species.ref, allowMultiple=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
#taxa.80 <- assignTaxonomy(seqtab.nochim, taxa.ref, minBoot=80)
taxa_data<-as.data.frame(taxa)
taxa_data$asvID<-paste0("ASV_",seq(1:nrow(taxa)))
taxa_data$seq<-rownames(taxa)
rownames(taxa_data)<-NULL
write.csv(taxa_data,paste0(dir,"taxa.csv"))
#
write.csv(seqtab.nochim,paste0(dir,"seqtab.csv"), col.names=NA)

# metadata = read_excel(paste0(dir,"samples.xlsx"))
# samdf <- data.frame(FASTQ = rownames(seqtab.nochim))
# samdf$Sample = sub("_", "-", samdf$FASTQ)
# samdf$Sample <- gsub("\\_.*", "", samdf$Sample)
# metadata = merge(samdf,metadata, by = "Sample", all = T)
# rownames(metadata) <- metadata$FASTQ.x
# write.csv(metadata,paste0(dir,"samples.csv"), col.names=NA)

# pyloseq visualization
rm(list = ls())
cat("\f")  
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
# read files
dir="C:/Users/Ofir Cohn/Projects/Microbiome/baseline_CC/" #folder am working in
taxmat = as.matrix(read.csv(paste0(dir,"taxa.csv"),row.names = 9))[,c(-1,-9)]
otumat = as.matrix(read.table(paste0(dir,"ASVs_counts.tsv")))
samples = read.csv(paste0(dir,"samples.csv"),row.names = 1)
# build phlyo 
ps <- phyloseq(otu_table(otumat, taxa_are_rows=T), 
               sample_data(samples), 
               tax_table(taxmat))
# remove midDNA, NA's and Controls
ps <- subset_taxa(ps, !(Family %in% c("Mitochondria"))) # remove mitDNA
ps <- subset_taxa(ps, (Kingdom %in% c("Bacteria")))
ps = subset_taxa(ps, Phylum != "NA")
ps = subset_taxa(ps, Family != "NA")
#ps <- filter_taxa(ps, function(x){sum(x > 0) > 1}, prune = TRUE)

#a = as.data.frame(tax_table(ps))
# plot with NC
plot_richness(ps, x="Sample", measures=c("Shannon", "Simpson"), color="Group")
#
ps_clean <- subset_samples(ps, Group!="CT")
plot_richness(ps_clean, x="Sample", measures=c("Shannon", "Simpson"), color="Group")
plot_richness(ps_clean, x="Sample", measures=c("Shannon", "Simpson"), color="Baseline_T") +
  scale_colour_gradientn(colours=c('blue','white','red'))

#### prevelance table
prevelancedf = apply(X = otu_table(ps_clean),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})
prevelancedf = data.frame(Prevalence = prevelancedf,
                          TotalAbundance = taxa_sums(ps_clean),
                          tax_table(ps_clean))
prevelancedf[1:10,]
# low prevelance/abundance phylum and subset them out.
summary_prevalence  = plyr::ddply(prevelancedf, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})
# Using the table above, determine the phyla to filter
sum(summary_prevalence$total_abundance)*0.0001
table(summary_prevalence$total_abundance/sum(summary_prevalence$total_abundance) <= 0.0001)
phyla2Filter <- summary_prevalence$Phylum[summary_prevalence$total_abundance/sum(summary_prevalence$total_abundance) <= 0.0001]
# Filter entries with unidentified Phylum.
ps.1 = subset_taxa(ps_clean, !Phylum %in% phyla2Filter)
summary_prevalence <- summary_prevalence[!summary_prevalence$Phylum %in% phyla2Filter,]
summary_prevalence
ps.1

# Subset to the remaining phyla by prevelance.
prevelancedf1 = subset(prevelancedf, Phylum %in% get_taxa_unique(ps.1, taxonomic.rank = "Phylum"))

ggplot(prevelancedf1, aes(TotalAbundance,Prevalence / nsamples(ps.1),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.10, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")


ggplot(prevelancedf1, aes(TotalAbundance, Prevalence / nsamples(ps.1),color=Family)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Genus) + theme(legend.position="none")

#  Define prevalence threshold as 10% of total samples ~ set of replicates
prevalenceThreshold = 0.05 * nsamples(ps.1)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevelancedf1)[(prevelancedf1$Prevalence >= prevalenceThreshold)]
length(keepTaxa)
ps.2 = prune_taxa(keepTaxa, ps.1)
ps.2

# Agglomerate taxa at the Genus level (combine all with the same name)
# keeping all taxa without genus level assignment

length(get_taxa_unique(ps.2, taxonomic.rank = "Genus"))
length(get_taxa_unique(ps.2, taxonomic.rank = "Family"))
length(get_taxa_unique(ps.2, taxonomic.rank = "Order"))

ps.3 = tax_glom(ps.2, "Genus", NArm = TRUE)
ps.3 = tax_glom(ps.2, "Family", NArm = TRUE)
ps.3
## out of curiosity how many "reads" does this leave us at???
sum(colSums(otu_table(ps.3)))
ntaxa(ps.2)
ntaxa(ps.3)

# Now lets filter out samples (outliers and low performing samples)
# Do some simple ordination looking for outlier samples, first we variance stabilize 
# the data with a log transform, the perform PCoA using bray's distances
library(ggfortify)
df <- iris[1:4]
pca_res <- prcomp(df, scale. = TRUE)

autoplot(pca_res)

summary(dd)

ee = as.data.frame(dd$x[,1:2])

comp = dd$rotation
dd %>% biplot(cex=0.5)  

ggplot() + 
  geom_point(data = ee, aes(x=PC1, y=PC2))


logt  = transform_sample_counts(ps.3, function(x) log(1 + x) )



emptycols <- sapply(aaaaa, function (k) all(is.na(k)))
df <- aaaaa[!emptycols]

out.pcoa.logt <- ordinate(logt, method = "PCoA", distance = "bray")
out.pcoa.logt <- ordinate(logt, method = "MDS", distance = "bray")
evals <- out.pcoa.logt$values$Eigenvalues
plot_ordination(logt, out.pcoa.logt, type = "Samples", 
                color = "Baseline_T", shape = "Group") + labs(col = "T") +
  coord_fixed(sqrt(evals[2] / evals[1]))+
  geom_point(size=3)+
  scale_colour_gradientn(colours=c('blue','white','red'))

# Show taxa proportions per sample (quickplot)
library(gridExtra)
grid.arrange(nrow = 3,
             qplot(as(otu_table(logt),"matrix")[, "T15_9_S8_L001_R1_001.fastq.gz"], geom = "histogram", bins=30) +
               xlab("Relative abundance"),
             
             qplot(as(otu_table(logt),"matrix")[, "T17_3_S6_L001_R1_001.fastq.gz"], geom = "histogram", bins=30) +
               xlab("Relative abundance"),
             
             qplot(as(otu_table(logt),"matrix")[, "T15_1_S5_L001_R1_001.fastq.gz"], geom = "histogram", bins=30) +
               xlab("Relative abundance")
)
# if you needed to remove candidate outliers, can use the below to remove sample 
#s16sV3V5.pruned <- prune_samples(sample_names(s16sV3V5.3) != c("sample1","sample2"), s16sV3V5.3)

# Look for low perfroming samples
qplot(colSums(otu_table(ps.3)),bins=50) + xlab("Logged counts-per-sample")
ps.4 <- prune_samples(sample_sums(ps.3)>=10000, ps.3)
ps.4


ps.4 <- prune_samples(sample_sums(ps.2)>=1000, ps.2)

# We transform microbiome count data to account for differences in library size
# variance, scale, etc.
# RLE - is the scaling factor method proposed by Anders and Huber (2010). 
# We call it "relative log expression", as median library is calculated 
# from the geometric mean of all columns and the median ratio of each sample 
# to the median library is taken as the scale factor.
## for Firmictures
plot_abundance = function(physeq, meta, title = "",
                          Facet = "Family", Color = "Family"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Proteobacteria"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = meta,y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}

# transform counts into "abundances"
library(vegan)
library(edgeR)
ps.4.ra = transform_sample_counts(ps.4, function(x){x / sum(x)})
# transform counts into "hellinger standardized counts"
ps.4hell <- ps.4
otu_table(ps.4hell) <-otu_table(decostand(otu_table(ps.4hell), method = "hellinger"), taxa_are_rows=TRUE)

# RLE counts
ps.4RLE <- ps.4.ra
RLE_normalization <- function(phyloseq){
  prior.count = 1
  count_scale = median(sample_sums(phyloseq))
  m = as(otu_table(phyloseq), "matrix")
  d = DGEList(counts=m, remove.zeros = FALSE)
  z = calcNormFactors(d, method="RLE")
  y <- as.matrix(z)
  lib.size <- z$samples$lib.size * z$samples$norm.factors
  ## rescale to median sample count
  out <- round(count_scale * sweep(y,MARGIN=2, STATS=lib.size,FUN = "/"))
  dimnames(out) <- dimnames(y)
  out
}
otu_table(ps.4RLE) <- otu_table(RLE_normalization(ps.4), taxa_are_rows=TRUE)
ps.4RLElog = transform_sample_counts(ps.4RLE, function(x){ log2(x +1)})
ps.4log <- transform_sample_counts(ps.4, function(x) log(1 + x))

plotOriginal = plot_abundance(ps.4, "Sample", title="original")
plotRelative = plot_abundance(ps.4.ra, "Sample", title="relative")
plotHellinger = plot_abundance(ps.4hell, "Sample", title="Hellinger")
plotRLE = plot_abundance(ps.4RLE, "Sample", title="LogRLE")
plotLogRLE = plot_abundance(ps.4RLElog, "Sample", title="LogRLE")

plotLog = plot_abundance(ps.4RLE, "Sample", title="Log")

# Combine each plot into one graphic.
grid.arrange(nrow = 5, plotOriginal, plotRelative, plotHellinger, plotLogRLE,plotLog)

# Graphical Summaries
plot_richness(ps.4,x="Sample", measures=c("Observed","Chao1"),color = "Baseline_T")+
  scale_colour_gradientn(colours=c('blue','white','red'))
plot_richness(ps.4,x="Sample", measures=c("Observed","Chao1"),color = "Baseline_R")+
  scale_colour_gradientn(colours=c('blue','white','red'))
plot_richness(ps.4RLE, x="Sample", measures=c("Observed", "Chao1"), color="Baseline_T")+
  scale_colour_gradientn(colours=c('blue','white','red'))

plot_bar(ps.4.ra, x="Sample", y="Abundance", fill="Phylum")
plot_bar(ps.4RLElog, x="Sample", y="Abundance", fill="Class")
plot_bar(ps.2, x="Sample", y="Abundance", fill="order")

# Other Richness measures, "Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher" try some of these others.
er <- estimate_richness(ps.4RLE, measures=c("Chao1", "Shannon"))
res.aov <- aov(er$Shannon ~ Baseline_R, data = as(sample_data(ps.4RLE),"data.frame"))
summary(res.aov)

# Subset dataset by phylum
Proteobacteria = subset_taxa(ps.4RLE, Phylum=="Proteobacteria")
title = "plot_bar; Proteobacteria-only"
plot_bar(Proteobacteria, "Sample", "Abundance", "Class", title=title)
prop  = transform_sample_counts(ps.4, function(x) x / sum(x) )
keepTaxa <- ((apply(otu_table(prop) >= 0.005,1,sum,na.rm=TRUE) > 2) | (apply(otu_table(prop) >= 0.05, 1, sum,na.rm=TRUE) > 0))
table(keepTaxa)
ps.4RLE_trim <- prune_taxa(keepTaxa,ps.4RLE)
plot_heatmap(ps.4RLE_trim, "PCoA", distance="bray", sample.label="Sample", taxa.label="Genus", low="#FFFFCC", high="#000033", na.value="white")
plot_net(ps.4RLE_trim, maxdist=0.4, color="Sample", shape="Group")
hell.tip.labels <- as(get_variable(ps.4RLE, "Sample"), "character")
# This is the actual hierarchical clustering call, specifying average-linkage clustering
d <- phyloseq::distance(ps.4RLE_trim, method="bray", type="samples")
d <- phyloseq::distance(ps.4RLE, method="bray", type="samples")

RLE.hclust     <- hclust(d, method="average")
plot(RLE.hclust)
#Lets write out a plot
pdf("My_dendro.pdf", width=7, height=7, pointsize=8)
plot(RLE.hclust)
dev.off()


# Ordination
v4.RLE.ord <- ordinate(ps.4RLE, "NMDS", "bray")
v4.RLE.ord <- ordinate(ps.4, "NMDS", "bray")

p1 = plot_ordination(ps.4, v4.RLE.ord, type="taxa", color="Phylum", title="taxa")
print(p1)
p1 + facet_wrap(~Phylum, 5)

p1 = plot_ordination(ps.4, v4.RLE.ord, type="taxa", color="Family", title="taxa")
print(p1)
p1 + facet_wrap(~Family, 5)



# Differential Abundances
m = as(otu_table(ps.4), "matrix")
m = as(otu_table(logt), "matrix")

# Define gene annotations (`genes`) as tax_table
taxonomy = tax_table(logt, errorIfNULL=FALSE)
if( !is.null(taxonomy) ){
  taxonomy = data.frame(as(taxonomy, "matrix"))
}
# Now turn into a DGEList
d = DGEList(counts=m, genes=taxonomy, remove.zeros = TRUE)

## reapply filter
prop  = transform_sample_counts(ps.4, function(x) x / sum(x) )
keepTaxa <- ((apply(otu_table(prop) >= 0.005,1,sum,na.rm=TRUE) > 2) | (apply(otu_table(prop) >= 0.05, 1, sum,na.rm=TRUE) > 0))
table(keepTaxa)

d <- d[keepTaxa,]

# Calculate the normalization factors
z = calcNormFactors(d, method="RLE")
# Check for division by zero inside `calcNormFactors`
if( !all(is.finite(z$samples$norm.factors)) ){
  stop("Something wrong with edgeR::calcNormFactors on this data,
       non-finite $norm.factors, consider changing `method` argument")
}

plotMDS(z, col = as.numeric(factor(sample_data(ps.4)$Baseline_T)), labels = sample_names(ps.4), cex=0.5)

# Creat a model based on Treatment and depth
mm <- model.matrix( ~Group, data=data.frame(as(sample_data(ps.4),"matrix"))) # specify model with no intercept for easier contrasts
mm
y <- voom(d, mm, plot = T)

fit <- lmFit(y, mm)
head(coef(fit))

# single contrast comparing Timepoint 5 - 20
contr <- makeContrasts(TimpointT2vT1 = "TimepointT2",
                       levels = colnames(coef(fit)))



theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=3) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))





################# microbiome biomarker
rm(list = ls())
cat("\f")  
library(microbiomeMarker)
library(phyloseq)
dir="C:/Users/Ofir Cohn/Projects/Microbiome/baseline_CC/" #folder am working in
taxmat = as.matrix(read.csv(paste0(dir,"taxa.csv"),row.names = 9))[,c(-1,-9)]
otumat = as.matrix(read.table(paste0(dir,"ASVs_counts.tsv")))
samples = read.csv(paste0(dir,"samples.csv"),row.names = 1)
# build phlyo 
ps <- phyloseq(otu_table(otumat, taxa_are_rows=T), 
               sample_data(samples), 
               tax_table(taxmat))
ps <- subset_samples(ps, Group!="CT")
ps <- subset_taxa(ps, (Kingdom %in% c("Bacteria")))
ps = subset_taxa(ps, Phylum != "NA")
ps = subset_taxa(ps, Family != "NA")

tg_welch <- run_test_two_groups(
  ps,
  transform = "log10p",
  norm = "RLE",
  taxa_rank = "all",
  group = "R_cut",
  method = "t.test"
)

head(marker_table(tg_welch))

var.test(aaa[,1:5], aaa[6:10], alternative = "two.sided")


# small example phyloseq object for test
ps_small <- phyloseq::subset_taxa(
  ps,
  Phylum %in% c("Firmicutes", "Proteobacteria")
)
mm_lr <- run_sl(
  ps_small,
  group = "R_cut",
  nfolds = 2,
  nrepeats = 1,
  taxa_rank = "Order",
  top_n = 15,
  norm = "TSS",
  method = "LR",
)
marker_table(mm_lr)

pht <- run_posthoc_test(ps, group = "R_cut")
pht
extract_posthoc_res(pht, "Staphylococcales")[[1]]
pht[[1]]


p_abd <- plot_abundance(tg_welch, group = "R_cut")
p_abd
plot_heatmap(tg_welch, transform = "log10p", group = "Order")
plot_ef_bar(tg_welch)
plot_sl_roc(mm_lr, group = "Order")

p_pht <- plot_postHocTest(pht, feature = "p__Bacteroidetes|g__Bacteroides")
p_pht