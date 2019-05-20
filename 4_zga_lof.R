# TRANSCRIPTIONAL EFFECT OF VARIOUS LOSS-OF-FUNCTIONS (LOF) ON ZYGOTIC GENOME ACTIVATION (ZGA) IN XENOPUS TROPICALIS
# Author: G.E. Gentsch
#
#
# Differential expression was performed using the likelihood ratio test of DESeq2
# Reference: Love et al. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2
#
# EXON and INTRON read counts from RNA-Seq experiments (exp1-5) in biological triplicates at indicated developmental stages
# Paired-end reads were aligned to the X. tropicalis genome assembly v7.1 using STAR v2.5.3a with default settings (Supplementary Table 1) and a revised version of gene models v7.2 to improve mapping accuracy across splice junctions.
# The alignments were sorted by read name using the sort function of samtools v1.3.1.
# EXON and INTRON counts (-t 'exon;intron') were extracted from unstranded (-s 0) alignment files using VERSE v0.1.569 in featureCounts (default) mode (-z 0). Intron coordinates were adjusted to exclude any overlap with exon annotation.
#
# EXPERIMENTS:
#
# Experiment 1 (exp1), n=3: Alfa-amanitin at stage 8+, 9+ and 10+
# Experiment 1 controls: uninjected embryos (uni)
#
# Experiment 2 and 3 (exp23) performed simultaneously, n=3: Wnt (exp2), Nodal (exp3) or BMP (exp3) LOF at stage 8+, 9+ and 10+
# - controls: control MO (cMO) and uninjected embryox (Experiment 2 only) and DMSO-treated embryos(Experiment 3 only)
#
# Experiment 4 (exp4): mPouV and mPouV/Sox3 LOF at stage 8+, 9+ and 10+
# - controls: uninjected embryos (uni)
#
# Experiment 5 (exp5): mVegT LOF at stage 8+, 9+ and 10+
# - controls: standard control MO (cMO) and uninjected (uni) embryos
#
# Spatial information (spatio): Regional expression across the animal-vegetal and dorso-ventral axes
# Reference for experiment 6: Blitz et al. (2017) A catalog of Xenopus tropicalis transcription factors and their regional expression in the early gastrula stage embryo
#
#
# Set working directory accordingly to execute zga_lof.R.
# Leave folder structure containing read counts unchanged.

#
# Installing (if required) and uploading several R packages
required.pkg <- c( "DESeq2", "scales", "ggplot2", "GenomicFeatures", "GenomicRanges", "rtracklayer", "extrafont", "GOstats", "GSEABase", "GO.db", "limma" )

for ( pkg in required.pkg ) {
    if ( pkg %in% rownames(installed.packages()) == FALSE) {
        install.packages( pkg ) }
    if ( pkg %in% rownames(.packages()) == FALSE ) {
        library(pkg, character.only = TRUE) }
}


# Auxiliary functions
source("./R_utils/name_dissector.R")
# functions: f1(), f2(), f3(), f4() and f5()
source("./R_utils/colour_schemes.R")
source("./R_utils/plot_2xLOF.R")
# function: plotDifSp()



# Extract gene ID and gene name from gff3 file
# ---------------------------------------------------

gff3                <- import.gff3( "./xenTro71/xenTro71.gff3.gz" )
gene.gff3           <- gff3[ which( elementMetadata( gff3 )[ ,"type" ] == "gene" ), ]
id.name             <- data.frame( elementMetadata( gene.gff3 )$ID, elementMetadata( gene.gff3 )$Name )
colnames(id.name)   <- c("ID","Name")
id.name$label       <- ifelse( id.name$ID %in% id.name$Name, paste( id.name$ID,"|", sep=""), paste( id.name$ID, id.name$Name, sep="|") )

# Extract gene ID from first EXON/INTRON count tables to add gene name (if available)
exon.1              <- read.table( "./RNA/exp1/uni_st8p_frogA_GG201.exon.txt", header=F, col.names=c( "ID", "counts" ) )
intron.1            <- read.table("./RNA/exp1/uni_st8p_frogA_GG201.intron.txt", header=F, col.names=c( "ID","counts" ) )

# Keep order of read count list by using plyr::join
id.exon             <- plyr::join(exon.1,id.name,by="ID")[,c(1,4)]
id.intron           <- plyr::join(intron.1,id.name,by="ID")[,c(1,4)]


# Differential expression analysis using DESeq2
# ----------------------------------------------------
#
# Enter read counts (separate for EXON [e] and INTRON [i] counts) per experiment and set factors
#
# Experiment 1
directory                   <- "./RNA/exp1"
sampleExon                  <- grep("exon",list.files(directory),value=T)
sampleIntron                <- grep("intron",list.files(directory),value=T)
cond                        <- sapply(sampleExon, f1)
st                          <- sapply(sampleExon, f2)
frog                        <- sapply(sampleExon, f3)
exonTable                   <- data.frame(sampleName=sampleExon,filename=sampleExon,condition=cond,stage=st,replicate=frog)
intronTable                 <- data.frame(sampleName=sampleIntron,filename=sampleIntron,condition=cond,stage=st,replicate=frog)
dds.e.exp1                  <- DESeqDataSetFromHTSeqCount(sampleTable=exonTable, directory=directory, design = ~ stage + condition)
dds.i.exp1                  <- DESeqDataSetFromHTSeqCount(sampleTable=intronTable, directory=directory, design = ~ stage + condition)
rownames(dds.e.exp1)        <- id.exon$label
rownames(dds.i.exp1)        <- id.intron$label
#
colData(dds.e.exp1)$stage   <- factor(colData(dds.e.exp1)$stage, levels=c("st8p","st9p","st10p"))
colData(dds.e.exp1)$condition <- factor(colData(dds.e.exp1)$condition, levels=c("uni","amanitin"))
colData(dds.e.exp1)$replicate <- factor(colData(dds.e.exp1)$replicate, levels=c("frogA","frogB","frogC"))
dds.e.exp1$stage            <- relevel(dds.e.exp1$stage,"st8p","st9p")
dds.e.exp1$condition        <- relevel(dds.e.exp1$condition,"uni")
colData(dds.i.exp1)$stage   <- factor(colData(dds.i.exp1)$stage, levels=c("st8p","st9p","st10p"))
colData(dds.i.exp1)$condition <- factor(colData(dds.i.exp1)$condition, levels=c("uni","amanitin"))
colData(dds.i.exp1)$replicate <- factor(colData(dds.i.exp1)$replicate, levels=c("frogA","frogB","frogC"))
dds.i.exp1$stage            <- relevel(dds.i.exp1$stage,"st8p","st9p")
dds.i.exp1$condition        <- relevel(dds.i.exp1$condition,"uni")

#
# Experiment 2 and 3
directory                   <- "./RNA/exp2_3"
sampleExon                  <- grep( "exon", list.files( directory ), value=TRUE )
sampleIntron                <- grep( "intron", list.files( directory ), value=TRUE )
cond                        <- sapply( sampleExon, f1 )
st                          <- sapply( sampleExon, f2 )
frog                        <- sapply( sampleExon, f3 )
exonTable                   <- data.frame( sampleName=sampleExon, filename=sampleExon, condition=cond, stage=st, replicate=frog )
intronTable                 <- data.frame( sampleName=sampleIntron, filename=sampleIntron, condition=cond, stage=st, replicate=frog )
dds.e.exp23                 <- DESeqDataSetFromHTSeqCount( sampleTable=exonTable, directory=directory, design = ~ stage + condition )
dds.i.exp23 <- DESeqDataSetFromHTSeqCount( sampleTable=intronTable, directory=directory, design = ~ stage + condition )
rownames( dds.e.exp23 ) <- id.exon$label
rownames( dds.i.exp23 ) <- id.intron$gene

colData( dds.e.exp23 )$stage <- factor( colData( dds.e.exp23 )$stage, levels=c( "st8p","st9p","st10p" ))
colData( dds.e.exp23 )$condition <- factor( colData( dds.e.exp23 )$condition, levels=c( "uni", "cMO", "dmso", "BmpLOF", "NodalLOF", "WntLOF" ) )
colData( dds.e.exp23 )$replicate <- factor( colData( dds.e.exp23 )$replicate, levels=c( "frogA", "frogB", "frogC" ) )
dds.e.exp23$stage           <- relevel( dds.e.exp23$stage, "st8p", "st9p" )
dds.e.exp23$condition       <- relevel( dds.e.exp23$condition, "uni", "cMO", "dmso" )
colData( dds.i.exp23 )$stage <- factor( colData( dds.i.exp23 )$stage, levels=c( "st8p", "st9p", "st10p" ) )
colData( dds.i.exp23 )$condition <- factor( colData( dds.i.exp23 )$condition, levels=c( "uni", "cMO", "dmso", "BmpLOF", "NodalLOF", "WntLOF" ) )
colData( dds.i.exp23 )$replicate <- factor( colData( dds.i.exp23 )$replicate, levels=c( "frogA", "frogB", "frogC" ) )
dds.i.exp23$stage           <- relevel( dds.i.exp23$stage, "st8p", "st9p" )
dds.i.exp23$condition       <- relevel( dds.i.exp23$condition, "uni", "cMO", "dmso" )

#
# Experiment 4
directory                   <- "./RNA/exp4"
sampleExon                  <- grep( "exon", list.files(directory), value=TRUE )
sampleIntron                <- grep("intron",list.files(directory), value=TRUE )
cond                        <- sapply( sampleExon, f1 )
st                          <- sapply( sampleExon, f2 )
frog                        <- sapply( sampleExon, f3 )
exonTable                   <- data.frame( sampleName=sampleExon, filename=sampleExon, condition=cond, stage=st, replicate=frog)
intronTable                 <- data.frame( sampleName=sampleIntron, filename=sampleIntron, condition=cond, stage=st, replicate=frog)
dds.e.exp4                  <- DESeqDataSetFromHTSeqCount( sampleTable=exonTable, directory=directory, design = ~ stage + condition)
dds.i.exp4                  <- DESeqDataSetFromHTSeqCount( sampleTable=intronTable, directory=directory, design = ~ stage + condition)
rownames(dds.e.exp4)        <- id.exon$label
rownames(dds.i.exp4)        <- id.intron$label

colData(dds.e.exp4)$stage   <- factor( colData(dds.e.exp4)$stage, levels=c("st8p","st9p","st10p") )
colData(dds.e.exp4)$condition <- factor( colData(dds.e.exp4)$condition, levels=c("uni","mPouVkd","mPSkd") )
colData(dds.e.exp4)$replicate <- factor( colData(dds.e.exp4)$replicate, levels=c("frogC","frogD","frogE") )
dds.e.exp4$stage            <- relevel( dds.e.exp4$stage,"st8p","st9p")
dds.e.exp4$condition        <- relevel( dds.e.exp4$condition,"uni")
colData(dds.i.exp4)$stage   <- factor( colData(dds.i.exp4)$stage, levels=c("st8p","st9p","st10p"))
colData(dds.i.exp4)$condition <- factor( colData(dds.i.exp4)$condition, levels=c("uni","mPouVkd","mPSkd") )
colData(dds.i.exp4)$replicate <- factor( colData(dds.i.exp4)$replicate, levels=c("frogC","frogD","frogE") )
dds.i.exp4$stage            <- relevel( dds.i.exp4$stage, "st8p","st9p" )
dds.i.exp4$condition        <- relevel( dds.i.exp4$condition, "uni" )


# Experiment 5
directory                   <- "./RNA/exp5"
sampleExon                  <- grep( "exon", list.files(directory), value=TRUE )
sampleIntron                <- grep( "intron", list.files(directory), value=TRUE )
cond                        <- sapply( sampleExon, f1 )
st                          <- sapply( sampleExon, f2 )
frog                        <- sapply( sampleExon, f3 )
exonTable                   <- data.frame( sampleName=sampleExon, filename=sampleExon, condition=cond, stage=st, replicate=frog )
intronTable                 <- data.frame( sampleName=sampleIntron, filename=sampleIntron, condition=cond, stage=st, replicate=frog )
dds.e.exp5                  <- DESeqDataSetFromHTSeqCount( sampleTable=exonTable, directory=directory, design = ~ stage + condition )
dds.i.exp5                  <- DESeqDataSetFromHTSeqCount( sampleTable=intronTable, directory=directory, design = ~ stage + condition )
rownames( dds.e.exp5 )      <- id.exon$label
rownames( dds.i.exp5 )      <- id.intron$label

colData(dds.e.exp5)$stage   <- factor(colData(dds.e.exp5)$stage, levels=c("st8p","st9p","st10p"))
colData(dds.e.exp5)$condition <- factor(colData(dds.e.exp5)$condition, levels=c("uni","cMO","mVegTkd"))
colData(dds.e.exp5)$replicate <- factor(colData(dds.e.exp5)$replicate, levels=c("frogA","frogE","frogF"))
dds.e.exp5$stage            <- relevel(dds.e.exp5$stage,"st8p","st9p")
dds.e.exp5$condition        <- relevel(dds.e.exp5$condition,"uni","cMO")
colData(dds.i.exp5)$stage   <- factor(colData(dds.i.exp5)$stage, levels=c("st8p","st9p","st10p"))
colData(dds.i.exp5)$condition <- factor(colData(dds.i.exp5)$condition, levels=c("uni","cMO","mVegTkd"))
colData(dds.i.exp5)$replicate <- factor(colData(dds.i.exp5)$replicate, levels=c("frogA","frogE","frogF"))
dds.i.exp5$stage            <- relevel(dds.i.exp5$stage,"st8p","st9p")
dds.i.exp5$condition        <- relevel(dds.i.exp5$condition,"uni","cMO")

# Spatial information
directory                   <- "./RNA/spatial_info"
sampleExon                  <- grep( "exon",list.files(directory), value=TRUE )
sampleIntron                <- grep( "intron",list.files(directory), value=TRUE )
cond                        <- sapply( sampleExon, f1 )
st                          <- sapply( sampleExon, f2 )
frog                        <- sapply( sampleExon, f3 )
exonTable                   <- data.frame( sampleName=sampleExon, filename=sampleExon, condition=cond, stage=st, replicate=frog )
intronTable                 <- data.frame( sampleName=sampleIntron, filename=sampleIntron, condition=cond, stage=st, replicate=frog )
dds.e.spatio                <- DESeqDataSetFromHTSeqCount( sampleTable=exonTable, directory=directory, design = ~ condition )
dds.i.spatio                <- DESeqDataSetFromHTSeqCount( sampleTable=intronTable, directory=directory, design = ~ condition )
rownames( dds.e.spatio )    <- id.exon$label
rownames( dds.i.spatio )    <- id.intron$label

colData( dds.e.spatio )$stage <- factor( colData( dds.e.spatio )$stage, levels=c("st10"))
colData( dds.e.spatio )$condition <- factor( colData( dds.e.spatio )$condition, levels=c("we","ac","vg","dmz","vmz"))
colData( dds.e.spatio )$replicate <- factor( colData( dds.e.spatio )$replicate, levels=c("frogA","frogB"))
dds.e.spatio$stage          <- relevel( dds.e.spatio$stage, "st10" )
dds.e.spatio$condition      <- relevel( dds.e.spatio$condition, "we" )
colData( dds.i.spatio )$stage <- factor( colData( dds.i.spatio )$stage, levels=c( "st10" ) )
colData( dds.i.spatio )$condition <- factor( colData(dds.i.spatio )$condition, levels=c( "we", "ac", "vg", "dmz", "vmz" ) )
colData( dds.i.spatio )$replicate <- factor( colData(dds.i.spatio )$replicate, levels=c( "frogA", "frogB" ) )
dds.i.spatio$stage          <- relevel( dds.i.spatio$stage, "st10" )
dds.i.spatio$condition      <- relevel( dds.i.spatio$condition, "we" )


# Estimates the size factors using the "median ratio method" described by Equation 5 in Anders and Huber (2010)
# Reference: Anders and Huber (2010) Differential expression analysis for sequence count data
dds.e.exp1                  <- estimateSizeFactors( dds.e.exp1 )
dds.e.exp23                 <- estimateSizeFactors( dds.e.exp23 )
dds.e.exp4                  <- estimateSizeFactors( dds.e.exp4 )
dds.e.exp5                  <- estimateSizeFactors( dds.e.exp5 )
dds.e.spatio                <- estimateSizeFactors( dds.e.spatio )

# INTRON count size estimation: use size factor from EXON counts
sizeFactors( dds.i.exp1 )   <- sizeFactors( dds.e.exp1 )
sizeFactors( dds.i.exp23 )  <- sizeFactors( dds.e.exp23 )
sizeFactors( dds.i.exp4 )   <- sizeFactors( dds.e.exp4 )
sizeFactors( dds.i.exp5 )   <- sizeFactors( dds.e.exp5 )
sizeFactors( dds.i.spatio ) <- sizeFactors( dds.e.spatio )


# obtain dispersion estimates for Negative Binomial (NB) distributed data
# :EXON
dds.e.exp1                  <- estimateDispersions( dds.e.exp1 )
dds.e.exp23                 <- estimateDispersions( dds.e.exp23 )
dds.e.exp4                  <- estimateDispersions( dds.e.exp4 )
dds.e.exp5                  <- estimateDispersions( dds.e.exp5 )
dds.e.spatio                <- estimateDispersions( dds.e.spatio )
# :INTRON
dds.i.exp1                  <- estimateDispersions( dds.i.exp1 )
dds.i.exp23                 <- estimateDispersions( dds.i.exp23 )
dds.i.exp4                  <- estimateDispersions( dds.i.exp4 )
dds.i.exp5                  <- estimateDispersions( dds.i.exp5 )
dds.i.spatio                <- estimateDispersions( dds.i.spatio )

#
# Subset data for LIKELIHOOD RATIO TEST
# -------------------------------------------------

# Experiment 1
# :EXON
dds.e.ama.exp1              <- dds.e.exp1[ , which( dds.e.exp1$condition == "uni" | dds.e.exp1$condition == "amanitin" ) ]
dds.e.ama.exp1$condition    <- relevel( dds.e.ama.exp1$condition,"uni" )
dds.e.ama.exp1$condition    <- droplevels( dds.e.ama.exp1$condition )
# :INTRON
dds.i.ama.exp1              <- dds.i.exp1[ , which( dds.i.exp1$condition == "uni" | dds.i.exp1$condition == "amanitin" ) ]
dds.i.ama.exp1$condition    <- relevel( dds.i.ama.exp1$condition,"uni" )
dds.i.ama.exp1$condition    <- droplevels( dds.i.ama.exp1$condition )

# Experiment 2 and 3
# :EXON
dds.e.bmp.exp23             <- dds.e.exp23[ , which( dds.e.exp23$condition == "dmso" | dds.e.exp23$condition == "BmpLOF" ) ]
dds.e.bmp.exp23$condition   <- relevel(dds.e.bmp.exp23$condition,"dmso")
dds.e.bmp.exp23$condition   <- droplevels( dds.e.bmp.exp23$condition )
dds.e.nodal.exp23           <- dds.e.exp23[ , which( dds.e.exp23$condition == "dmso" | dds.e.exp23$condition == "NodalLOF" ) ]
dds.e.nodal.exp23$condition <- relevel(dds.e.nodal.exp23$condition,"dmso")
dds.e.nodal.exp23$condition <- droplevels( dds.e.nodal.exp23$condition )
dds.e.wnt.exp23             <- dds.e.exp23[ , which( dds.e.exp23$condition == "uni" | dds.e.exp23$condition == "WntLOF" ) ]
dds.e.wnt.exp23$condition   <- relevel( dds.e.wnt.exp23$condition,"uni" )
dds.e.wnt.exp23$condition   <- droplevels( dds.e.wnt.exp23$condition )
dds.e.cMO.exp23             <- dds.e.exp23[ , which( dds.e.exp23$condition == "uni" | dds.e.exp23$condition == "cMO" ) ]
dds.e.cMO.exp23$condition   <- relevel( dds.e.cMO.exp23$condition,"uni" )
dds.e.cMO.exp23$condition   <- droplevels( dds.e.cMO.exp23$condition )
# :INTRON
dds.i.bmp.exp23             <- dds.i.exp23[ , which( dds.i.exp23$condition == "dmso" | dds.i.exp23$condition == "BmpLOF" ) ]
dds.i.bmp.exp23$condition   <- relevel( dds.i.bmp.exp23$condition,"dmso" )
dds.i.bmp.exp23$condition   <- droplevels( dds.i.bmp.exp23$condition )
dds.i.nodal.exp23           <- dds.i.exp23[ , which( dds.i.exp23$condition == "dmso" | dds.i.exp23$condition == "NodalLOF" ) ]
dds.i.nodal.exp23$condition <- relevel( dds.i.nodal.exp23$condition,"dmso" )
dds.i.nodal.exp23$condition <- droplevels( dds.i.nodal.exp23$condition )
dds.i.wnt.exp23             <- dds.i.exp23[ , which( dds.i.exp23$condition == "uni" | dds.i.exp23$condition == "WntLOF" ) ]
dds.i.wnt.exp23$condition   <- relevel(dds.i.wnt.exp23$condition,"uni")
dds.i.wnt.exp23$condition   <- droplevels( dds.i.wnt.exp23$condition )
dds.i.cMO.exp23             <- dds.i.exp23[ , which( dds.i.exp23$condition == "uni" | dds.i.exp23$condition == "cMO" ) ]
dds.i.cMO.exp23$condition   <- relevel( dds.i.cMO.exp23$condition,"uni" )
dds.i.cMO.exp23$condition   <- droplevels( dds.i.cMO.exp23$condition )

# Experiment 4
# :EXON
dds.e.mPouV.exp4            <- dds.e.exp4[ , which( dds.e.exp4$condition == "uni" | dds.e.exp4$condition == "mPouVkd" ) ]
dds.e.mPouV.exp4$condition  <- relevel( dds.e.mPouV.exp4$condition,"uni" )
dds.e.mPouV.exp4$condition  <- droplevels( dds.e.mPouV.exp4$condition )
dds.e.mPS.exp4              <- dds.e.exp4[ , which( dds.e.exp4$condition == "uni" | dds.e.exp4$condition == "mPSkd" ) ]
dds.e.mPS.exp4$condition    <- relevel( dds.e.mPS.exp4$condition,"uni" )
dds.e.mPS.exp4$condition    <- droplevels( dds.e.mPS.exp4$condition )
# :INTRON
dds.i.mPouV.exp4            <- dds.i.exp4[ , which( dds.i.exp4$condition == "uni" | dds.i.exp4$condition == "mPouVkd" ) ]
dds.i.mPouV.exp4$condition  <- relevel( dds.i.mPouV.exp4$condition,"uni" )
dds.i.mPouV.exp4$condition  <- droplevels( dds.i.mPouV.exp4$condition )
dds.i.mPS.exp4              <- dds.i.exp4[ , which( dds.i.exp4$condition == "uni" | dds.i.exp4$condition == "mPSkd" ) ]
dds.i.mPS.exp4$condition    <- relevel( dds.i.mPS.exp4$condition, "uni" )
dds.i.mPS.exp4$condition    <- droplevels( dds.i.mPS.exp4$condition )

# Experiment 5
# :EXON
dds.e.vegt.exp5             <- dds.e.exp5[ , which( dds.e.exp5$condition == "uni" | dds.e.exp5$condition == "mVegTkd" )]
dds.e.vegt.exp5$condition   <- relevel( dds.e.vegt.exp5$condition, "uni" )
dds.e.vegt.exp5$condition   <- droplevels( dds.e.vegt.exp5$condition )
dds.e.cMO.exp5              <- dds.e.exp5[ , which( dds.e.exp5$condition == "uni" | dds.e.exp5$condition == "cMO" )]
dds.e.cMO.exp5$condition    <- relevel( dds.e.cMO.exp5$condition,"uni" )
dds.e.cMO.exp5$condition    <- droplevels( dds.e.cMO.exp5$condition )
# :INTRON
dds.i.vegt.exp5             <- dds.i.exp5[ , which( dds.i.exp5$condition == "uni" | dds.i.exp5$condition == "mVegTkd" )]
dds.i.vegt.exp5$condition   <- relevel( dds.i.vegt.exp5$condition,"uni" )
dds.i.vegt.exp5$condition   <- droplevels( dds.i.vegt.exp5$condition )
dds.i.cMO.exp5              <- dds.i.exp5[ , which( dds.i.exp5$condition == "uni" | dds.i.exp5$condition == "cMO" )]
dds.i.cMO.exp5$condition    <- relevel( dds.i.cMO.exp5$condition,"uni" )
dds.i.cMO.exp5$condition    <- droplevels( dds.i.cMO.exp5$condition )

# Spatial Expression
# :EXON
dds.e.av.spatio             <- dds.e.spatio[ , which( dds.e.spatio$condition == "ac" | dds.e.spatio$condition == "vg" ) ]
dds.e.av.spatio$condition   <- relevel( dds.e.av.spatio$condition,"vg" )
dds.e.av.spatio$condition   <- droplevels( dds.e.av.spatio$condition )
dds.e.dv.spatio             <- dds.e.spatio[ , which( dds.e.spatio$condition == "dmz" | dds.e.spatio$condition == "vmz" ) ]
dds.e.dv.spatio$condition   <- relevel( dds.e.dv.spatio$condition, "vmz" )
dds.e.dv.spatio$condition   <- droplevels( dds.e.dv.spatio$condition )
# :INTRON
dds.i.av.spatio             <- dds.i.spatio[ , which( dds.i.spatio$condition == "ac" | dds.i.spatio$condition == "vg" ) ]
dds.i.av.spatio$condition   <- relevel( dds.i.av.spatio$condition, "vg" )
dds.i.av.spatio$condition   <- droplevels( dds.i.av.spatio$condition )
dds.i.dv.spatio             <- dds.i.spatio[ , which( dds.i.spatio$condition == "dmz" | dds.i.spatio$condition == "vmz" ) ]
dds.i.dv.spatio$condition   <- relevel( dds.i.dv.spatio$condition,"vmz")
dds.i.dv.spatio$condition   <- droplevels( dds.i.dv.spatio$condition )


#
# LIKELIHOOD RATIO TEST (LRT)
# --------------------------------------------

# Experiment 1
# :EXON
dds.e.ama.exp1              <- nbinomLRT( dds.e.ama.exp1, reduced = ~ stage )
# :INTRON
dds.i.ama.exp1              <- nbinomLRT( dds.i.ama.exp1, reduced = ~ stage )


# Experiment 2 and 3
# :EXON
dds.e.bmp.exp23             <- nbinomLRT( dds.e.bmp.exp23, reduced = ~ stage )
dds.e.nodal.exp23           <- nbinomLRT( dds.e.nodal.exp23, reduced = ~ stage )
dds.e.wnt.exp23             <- nbinomLRT( dds.e.wnt.exp23, reduced = ~ stage )
dds.e.cMO.exp23             <- nbinomLRT( dds.e.cMO.exp23, reduced = ~ stage )
# :INTRON
dds.i.bmp.exp23             <- nbinomLRT( dds.i.bmp.exp23, reduced = ~ stage )
dds.i.nodal.exp23           <- nbinomLRT( dds.i.nodal.exp23, reduced = ~ stage )
dds.i.wnt.exp23             <- nbinomLRT( dds.i.wnt.exp23, reduced = ~ stage )
dds.i.cMO.exp23             <- nbinomLRT( dds.i.cMO.exp23, reduced = ~ stage )

# Experiment 4
# :EXON
dds.e.mPouV.exp4            <- nbinomLRT( dds.e.mPouV.exp4, reduced = ~ stage )
dds.e.mPS.exp4              <- nbinomLRT( dds.e.mPS.exp4, reduced = ~ stage )
# :INTRON
dds.i.mPouV.exp4            <- nbinomLRT( dds.i.mPouV.exp4, reduced = ~ stage )
dds.i.mPS.exp4              <- nbinomLRT( dds.i.mPS.exp4, reduced = ~ stage )

# Experiment 5
# :EXON
dds.e.vegt.exp5             <- nbinomLRT( dds.e.vegt.exp5, reduced = ~ stage )
dds.e.cMO.exp5              <- nbinomLRT( dds.e.cMO.exp5, reduced = ~ stage )
# :INTRON
dds.i.vegt.exp5             <- nbinomLRT( dds.i.vegt.exp5, reduced = ~ stage )
dds.i.cMO.exp5              <- nbinomLRT( dds.i.cMO.exp5, reduced = ~ stage )

# Spatial information
# :EXON
dds.e.av.spatio             <- nbinomLRT( dds.e.av.spatio, reduced = ~ 1 )
dds.e.dv.spatio             <- nbinomLRT( dds.e.dv.spatio, reduced = ~ 1 )
# :INTRON
dds.i.av.spatio             <- nbinomLRT( dds.i.av.spatio, reduced = ~ 1 )
dds.i.dv.spatio             <- nbinomLRT( dds.i.dv.spatio, reduced = ~ 1 )


# LIKELIHOOD RATIO TEST results (switch off Cook's cutoff and independent filtering)
# Experiment 1
# :EXON
res.dds.e.ama.exp1          <- results( dds.e.ama.exp1, cooksCutoff=FALSE, independentFiltering=FALSE )
# :INTRON
res.dds.i.ama.exp1          <- results( dds.i.ama.exp1, cooksCutoff=FALSE, independentFiltering=FALSE )

# Experiment 2 and 3
# :EXON
res.dds.e.bmp.exp23         <- results( dds.e.bmp.exp23, cooksCutoff=FALSE, independentFiltering=FALSE )
res.dds.e.nodal.exp23       <- results( dds.e.nodal.exp23, cooksCutoff=FALSE, independentFiltering=FALSE )
res.dds.e.wnt.exp23         <- results( dds.e.wnt.exp23, cooksCutoff=FALSE, independentFiltering=FALSE )
res.dds.e.cMO.exp23         <- results( dds.e.cMO.exp23, cooksCutoff=FALSE, independentFiltering=FALSE )
# :INTRON
res.dds.i.bmp.exp23         <- results( dds.i.bmp.exp23, cooksCutoff=FALSE, independentFiltering=FALSE )
res.dds.i.nodal.exp23       <- results( dds.i.nodal.exp23, cooksCutoff=FALSE, independentFiltering=FALSE )
res.dds.i.wnt.exp23         <- results( dds.i.wnt.exp23, cooksCutoff=FALSE, independentFiltering=FALSE )
res.dds.i.cMO.exp23         <- results( dds.i.cMO.exp23, cooksCutoff=FALSE, independentFiltering=FALSE )

# Experiment 4
# :EXON
res.dds.e.mPouV.exp4        <- results( dds.e.mPouV.exp4, cooksCutoff=FALSE, independentFiltering=FALSE )
res.dds.e.mPS.exp4          <- results( dds.e.mPS.exp4, cooksCutoff=FALSE, independentFiltering=FALSE )
# :INTRON
res.dds.i.mPouV.exp4        <- results( dds.i.mPouV.exp4, cooksCutoff=FALSE, independentFiltering=FALSE )
res.dds.i.mPS.exp4          <- results( dds.i.mPS.exp4, cooksCutoff=FALSE, independentFiltering=FALSE )

# Experiment 5
# :EXON
res.dds.e.vegt.exp5         <- results( dds.e.vegt.exp5, cooksCutoff=FALSE, independentFiltering=FALSE )
res.dds.e.cMO.exp5          <- results( dds.e.cMO.exp5, cooksCutoff=FALSE, independentFiltering=FALSE )
# :INTRON
res.dds.i.vegt.exp5         <- results( dds.i.vegt.exp5, cooksCutoff=FALSE, independentFiltering=FALSE )
res.dds.i.cMO.exp5          <- results( dds.i.cMO.exp5, cooksCutoff=FALSE, independentFiltering=FALSE )

# Spatial information
# :EXON
res.dds.e.av.spatio         <- results( dds.e.av.spatio, cooksCutoff=FALSE, independentFiltering=FALSE )
res.dds.e.dv.spatio         <- results( dds.e.dv.spatio, cooksCutoff=FALSE, independentFiltering=FALSE )
# :INTRON
res.dds.i.av.spatio         <- results( dds.i.av.spatio, cooksCutoff=FALSE, independentFiltering=FALSE )
res.dds.i.dv.spatio         <- results( dds.i.dv.spatio, cooksCutoff=FALSE, independentFiltering=FALSE )

#
#
#
#
# Differential ZGA: calculate fold changes from stages when gene is sensitive to alfa-amanitin (>50% down-regulation) and adjusted LRT p-value <= 10%
# --------------------------------------------------------------------------------------------------------------------------------------------------------

# Extract normalised EXON and INTRON read counts from "DESeqDataSet" using accessor "counts" and convert them to separate data frames
norm.df <- function(DESeqDS) {
    DF <- as.data.frame( counts( DESeqDS, normalized=TRUE ) )
    colnames( DF ) <- paste( colnames( DESeqDS ), ".ddsNorm", sep="" )
    return( DF )
}

# :EXON
norm.dds.e.exp1             <- norm.df( dds.e.exp1 )
norm.dds.e.exp23            <- norm.df( dds.e.exp23 )
norm.dds.e.exp4             <- norm.df( dds.e.exp4 )
norm.dds.e.exp5             <- norm.df( dds.e.exp5 )

# :INTRON
norm.dds.i.exp1             <- norm.df( dds.i.exp1 )
norm.dds.i.exp23            <- norm.df( dds.i.exp23 )
norm.dds.i.exp4             <- norm.df( dds.i.exp4 )
norm.dds.i.exp5             <- norm.df( dds.i.exp5 )


# Calculate mean from biological replicates
means.exp <- function( normDF ) {
    condition <- unique( sapply( colnames( normDF ), f1.2) )
    means <- sapply( condition, function(x) cbind( apply( normDF[ , grep( paste("^",x,sep=""), colnames(normDF)) ], 1, mean ) ) )
    ei <- strsplit( deparse( substitute( normDF ) ), ".", fixed=T)[[1]][3]
    exp <- strsplit( deparse( substitute( normDF ) ), ".", fixed=T)[[1]][4]
    colnames(means) <- paste0( gsub("_",".",colnames(means)), ".", ei, ".", exp )
    means <- as.data.frame(means)
    row.names(means) <- if(ei == "e") id.exon$label else id.intron$label
    return(means)
}

# :EXON
means.e.exp1                <- means.exp( norm.dds.e.exp1 )
means.e.exp23               <- means.exp( norm.dds.e.exp23 )
means.e.exp4                <- means.exp( norm.dds.e.exp4 )
means.e.exp5                <- means.exp( norm.dds.e.exp5 )

# :INTRON
means.i.exp1                <- means.exp( norm.dds.i.exp1 )
means.i.exp23               <- means.exp( norm.dds.i.exp23 )
means.i.exp4                <- means.exp( norm.dds.i.exp4 )
means.i.exp5                <- means.exp( norm.dds.i.exp5 )

means.e.allExp              <- cbind( means.e.exp1, means.e.exp23, means.e.exp4, means.e.exp5 )
means.i.allExp              <- cbind( means.i.exp1, means.i.exp23, means.i.exp4, means.i.exp5 )


# Calculate fold change in percentages
attach(means.e.allExp)
attach(means.i.allExp)

# Experiment 1
# :EXON
amanitin.st8p.uni.e.exp1    <- 100 / uni.st8p.e.exp1 * amanitin.st8p.e.exp1
amanitin.st9p.uni.e.exp1    <- 100 / uni.st9p.e.exp1 * amanitin.st9p.e.exp1
amanitin.st10p.uni.e.exp1   <- 100 / uni.st10p.e.exp1 * amanitin.st10p.e.exp1
# :INTRON
amanitin.st8p.uni.i.exp1    <- 100 / uni.st8p.i.exp1 * amanitin.st8p.i.exp1
amanitin.st9p.uni.i.exp1    <- 100 / uni.st9p.i.exp1 * amanitin.st9p.i.exp1
amanitin.st10p.uni.i.exp1   <- 100 / uni.st10p.i.exp1 * amanitin.st10p.i.exp1

# Experiment 2 and 3
# :EXON
BmpLOF.st8p.dmso.e.exp23    <- 100 / dmso.st8p.e.exp23 * BmpLOF.st8p.e.exp23
BmpLOF.st9p.dmso.e.exp23    <- 100 / dmso.st9p.e.exp23 * BmpLOF.st9p.e.exp23
BmpLOF.st10p.dmso.e.exp23   <- 100 / dmso.st10p.e.exp23 * BmpLOF.st10p.e.exp23
NodalLOF.st8p.dmso.e.exp23  <- 100 / dmso.st8p.e.exp23 * NodalLOF.st8p.e.exp23
NodalLOF.st9p.dmso.e.exp23  <- 100 / dmso.st9p.e.exp23 * NodalLOF.st9p.e.exp23
NodalLOF.st10p.dmso.e.exp23 <- 100 / dmso.st10p.e.exp23 * NodalLOF.st10p.e.exp23
WntLOF.st8p.uni.e.exp23     <- 100 / uni.st8p.e.exp23 * WntLOF.st8p.e.exp23
WntLOF.st9p.uni.e.exp23     <- 100 / uni.st9p.e.exp23 * WntLOF.st9p.e.exp23
WntLOF.st10p.uni.e.exp23    <- 100 / uni.st10p.e.exp23 * WntLOF.st10p.e.exp23
cMO.st8p.uni.e.exp23        <- 100 / uni.st8p.e.exp23 * cMO.st8p.e.exp23
cMO.st9p.uni.e.exp23        <- 100 / uni.st9p.e.exp23 * cMO.st9p.e.exp23
cMO.st10p.uni.e.exp23       <- 100 / uni.st10p.e.exp23 * cMO.st10p.e.exp23
# :INTRON
BmpLOF.st8p.dmso.i.exp23    <- 100 / dmso.st8p.i.exp23 * BmpLOF.st8p.i.exp23
BmpLOF.st9p.dmso.i.exp23    <- 100 / dmso.st9p.i.exp23 * BmpLOF.st9p.i.exp23
BmpLOF.st10p.dmso.i.exp23   <- 100 / dmso.st10p.i.exp23 * BmpLOF.st10p.i.exp23
NodalLOF.st8p.dmso.i.exp23  <- 100 / dmso.st8p.i.exp23 * NodalLOF.st8p.i.exp23
NodalLOF.st9p.dmso.i.exp23  <- 100 / dmso.st9p.i.exp23 * NodalLOF.st9p.i.exp23
NodalLOF.st10p.dmso.i.exp23 <- 100 / dmso.st10p.i.exp23 * NodalLOF.st10p.i.exp23
WntLOF.st8p.uni.i.exp23     <- 100 / uni.st8p.i.exp23 * WntLOF.st8p.i.exp23
WntLOF.st9p.uni.i.exp23     <- 100 / uni.st9p.i.exp23 * WntLOF.st9p.i.exp23
WntLOF.st10p.uni.i.exp23    <- 100 / uni.st10p.i.exp23 * WntLOF.st10p.i.exp23
cMO.st8p.uni.i.exp23        <- 100 / uni.st8p.i.exp23 * cMO.st8p.i.exp23
cMO.st9p.uni.i.exp23        <- 100 / uni.st9p.i.exp23 * cMO.st9p.i.exp23
cMO.st10p.uni.i.exp23       <- 100 / uni.st10p.i.exp23 * cMO.st10p.i.exp23

# Experiment 4
# :EXON
mPouVLOF.st8p.uni.e.exp4    <- 100 / uni.st8p.e.exp4 * mPouVkd.st8p.e.exp4
mPouVLOF.st9p.uni.e.exp4    <- 100 / uni.st9p.e.exp4 * mPouVkd.st9p.e.exp4
mPouVLOF.st10p.uni.e.exp4   <- 100 / uni.st10p.e.exp4 * mPouVkd.st10p.e.exp4
mPSLOF.st8p.uni.e.exp4      <- 100 / uni.st8p.e.exp4 * mPSkd.st8p.e.exp4
mPSLOF.st9p.uni.e.exp4      <- 100 / uni.st9p.e.exp4 * mPSkd.st9p.e.exp4
mPSLOF.st10p.uni.e.exp4     <- 100 / uni.st10p.e.exp4 * mPSkd.st10p.e.exp4
# :INTRON
mPouVLOF.st8p.uni.i.exp4    <- 100 / uni.st8p.i.exp4 * mPouVkd.st8p.i.exp4
mPouVLOF.st9p.uni.i.exp4    <- 100 / uni.st9p.i.exp4 * mPouVkd.st9p.i.exp4
mPouVLOF.st10p.uni.i.exp4   <- 100 / uni.st10p.i.exp4 * mPouVkd.st10p.i.exp4
mPSLOF.st8p.uni.i.exp4      <- 100 / uni.st8p.i.exp4 * mPSkd.st8p.i.exp4
mPSLOF.st9p.uni.i.exp4      <- 100 / uni.st9p.i.exp4 * mPSkd.st9p.i.exp4
mPSLOF.st10p.uni.i.exp4     <- 100 / uni.st10p.i.exp4 * mPSkd.st10p.i.exp4

# Experiment 5
# :EXON
mVegTLOF.st8p.uni.e.exp5    <- 100 / uni.st8p.e.exp5 * mVegTkd.st8p.e.exp5
mVegTLOF.st9p.uni.e.exp5    <- 100 / uni.st9p.e.exp5 * mVegTkd.st9p.e.exp5
mVegTLOF.st10p.uni.e.exp5   <- 100 / uni.st10p.e.exp5 * mVegTkd.st10p.e.exp5
cMO.st8p.uni.e.exp5         <- 100 / uni.st8p.e.exp5 * cMO.st8p.e.exp5
cMO.st9p.uni.e.exp5         <- 100 / uni.st9p.e.exp5 * cMO.st9p.e.exp5
cMO.st10p.uni.e.exp5        <- 100 / uni.st10p.e.exp5 * cMO.st10p.e.exp5
# :INTRON
mVegTLOF.st8p.uni.i.exp5    <- 100 / uni.st8p.i.exp5 * mVegTkd.st8p.i.exp5
mVegTLOF.st9p.uni.i.exp5    <- 100 / uni.st9p.i.exp5 * mVegTkd.st9p.i.exp5
mVegTLOF.st10p.uni.i.exp5   <- 100 / uni.st10p.i.exp5 * mVegTkd.st10p.i.exp5
cMO.st8p.uni.i.exp5         <- 100 / uni.st8p.i.exp5 * cMO.st8p.i.exp5
cMO.st9p.uni.i.exp5         <- 100 / uni.st9p.i.exp5 * cMO.st9p.i.exp5
cMO.st10p.uni.i.exp5        <- 100 / uni.st10p.i.exp5 * cMO.st10p.i.exp5

detach( means.e.allExp )
detach( means.i.allExp )


# Assemble expression percentages in single data frame
# :EXON
perc.e.coln <- c( ls( pattern="[089]p.uni.e.exp[1245]" ),
            ls( pattern="[089]p.cMO.e.exp[1245]" ),
            ls( pattern="[089]p.dmso.e.exp[1245]" ) )

perc.e.allExp <- sapply( perc.e.coln, function(x) cbind(get(x)) )
row.names( perc.e.allExp ) <- id.exon$label

# :INTRON
perc.i.coln <- c( ls( pattern="[089]p.uni.i.exp[1245]" ),
            ls( pattern="[089]p.cMO.i.exp[1245]" ),
            ls( pattern="[089]p.dmso.i.exp[1245]" ) )

perc.i.allExp <- sapply( perc.i.coln, function(x) cbind(get(x)) )
row.names(perc.i.allExp) <- id.intron$label

pm.e.allExp <- cbind(perc.e.allExp,means.e.allExp)
pm.i.allExp <- cbind(perc.i.allExp,means.i.allExp)


# Filter out 'low-read-count' genes using EXON counts as follows:
# min. 10 reads at least at one stage over all controls from each experiment.
# Reads were DESeq2-normalised (but no regularised log shrinkage, which moderates fold change, see DESeq2 vignette, p. 17)

index.1                 <- apply( means.e.exp1, 1, function(x) any( x[ grep("uni",colnames(means.e.exp1)) ] >= 10 ) )
index.23a               <- apply( means.e.exp23, 1, function(x) any( x[ grep("uni",colnames(means.e.exp23)) ] >= 10 ) )
index.23b               <- apply( means.e.exp23, 1, function(x) any( x[ grep("dmso",colnames(means.e.exp23)) ] >= 10 ) )
index.4                 <- apply( means.e.exp4, 1, function(x) any( x[ grep("uni",colnames(means.e.exp4)) ] >= 10 ) )
index.5                 <- apply( means.e.exp5, 1, function(x) any( x[ grep("uni",colnames(means.e.exp5)) ] >= 10 ) )

index.u                 <- index.1 + index.23a + index.23b + index.4 + index.5

pm.e.allExp.m10         <- pm.e.allExp[ index.u==5, ]
pm.e.allExp.m10.m       <- as.matrix( pm.e.allExp.m10 )


# Introduce NAs into stage-specific differential expression analysis for each stage if any control is below 10 norm. reads and convert to NA at stage 8 where alfa-amanitin >66.666%
# This is to remove unreliable gene expression data and expression levels at stage 8+ independent of alfa-amanitin.

st8.ctrls               <- colnames( means.e.allExp )[ c( grep("uni.st8p", colnames(means.e.allExp)), grep("dmso.st8p", colnames( means.e.allExp ))) ]

index.st8               <- apply( pm.e.allExp.m10, 1, function(x) any( x[st8.ctrls] < 10 ) )
pm.e.allExp.m10.m[ index.st8, grep( "st8p", colnames(perc.e.allExp) ) ] <- NA
index.st8.ama66         <- which(pm.e.allExp.m10[,"amanitin.st8p.uni.e.exp1"] > 2/3*100 )
pm.e.allExp.m10.m[ index.st8.ama66, grep( "st8p", colnames(perc.e.allExp) ) ] <- NA

st9.ctrls               <- colnames( means.e.allExp )[ c(grep("uni.st9p", colnames( means.e.allExp )), grep("dmso.st9p", colnames( means.e.allExp )))]
index.st9               <- apply( pm.e.allExp.m10, 1, function(x) any( x[st9.ctrls] < 10 ) )
pm.e.allExp.m10.m[ index.st9 , grep( "st9p", colnames(perc.e.allExp) ) ] <- NA

st10.ctrls              <- colnames( means.e.allExp )[ c(grep("uni.st10p", colnames(means.e.allExp)), grep("dmso.st10p", colnames( means.e.allExp ))) ]
index.st10              <- apply( pm.e.allExp.m10, 1, function(x) any( x[st10.ctrls] < 10 ) )
pm.e.allExp.m10.m[ index.st10 , grep( "st10p", colnames(perc.e.allExp) ) ] <- NA


# Average EXON expression changes
amanitin.st8to10.uni.e.mean <- apply( pm.e.allExp.m10.m[ ,grep( "^amanitin",colnames(pm.e.allExp.m10.m[,1:27])) ],1, mean, na.rm=TRUE )
cMO.st8to10.uni.e.mean  <- apply( pm.e.allExp.m10.m[ ,grep( "^cMO",colnames(pm.e.allExp.m10.m[,1:27])) ],1, mean, na.rm=TRUE )
mPouVLOF.st8to10.uni.e.mean <- apply( pm.e.allExp.m10.m[ ,grep( "^mPouVLOF", colnames(pm.e.allExp.m10.m[,1:27])) ],1,mean, na.rm=TRUE )
mPSLOF.st8to10.uni.e.mean <- apply( pm.e.allExp.m10.m[ ,grep( "^mPSLOF", colnames(pm.e.allExp.m10.m[,1:27])) ],1,mean, na.rm=TRUE )
mVegTLOF.st8to10.uni.e.mean <- apply( pm.e.allExp.m10.m[ ,grep( "^mVegT", colnames(pm.e.allExp.m10.m[,1:27])) ],1,mean,na.rm=TRUE )
WntLOF.st8to10.uni.e.mean <- apply( pm.e.allExp.m10.m[ ,grep( "^WntLOF", colnames(pm.e.allExp.m10.m[,1:27])) ],1,mean,na.rm=TRUE )
NodalLOF.st8to10.uni.e.mean <- apply( pm.e.allExp.m10.m[ ,grep( "^NodalLOF", colnames(pm.e.allExp.m10.m[,1:27])) ],1,mean,na.rm=TRUE )
BmpLOF.st8to10.uni.e.mean <- apply( pm.e.allExp.m10.m[ ,grep( "^BmpLOF", colnames(pm.e.allExp.m10.m[,1:27])) ],1,mean,na.rm=TRUE )

perc.e <- cbind(    amanitin.st8to10.uni.e.mean,
                    cMO.st8to10.uni.e.mean,
                    mPouVLOF.st8to10.uni.e.mean,
                    mPSLOF.st8to10.uni.e.mean,
                    mVegTLOF.st8to10.uni.e.mean,
                    WntLOF.st8to10.uni.e.mean,
                    NodalLOF.st8to10.uni.e.mean,
                    BmpLOF.st8to10.uni.e.mean
                )

# Select conservative FDR: cMO in experiment 2 and 5, use less significant FDR.
padj.e <- cbind(    res.dds.e.ama.exp1$padj,
                    ifelse(res.dds.e.cMO.exp23$padj >= res.dds.e.cMO.exp5$padj, res.dds.e.cMO.exp23$padj, res.dds.e.cMO.exp5$padj),
                    res.dds.e.mPouV.exp4$padj,
                    res.dds.e.mPS.exp4$padj,
                    res.dds.e.vegt.exp5$padj,
                    res.dds.e.wnt.exp23$padj,
                    res.dds.e.nodal.exp23$padj,
                    res.dds.e.bmp.exp23$padj
                )

row.names(padj.e) <- row.names(res.dds.e.ama.exp1)

colnames(padj.e) <- c(
                    "amanitin.st8to10.uni.e.padj",
                    "cMO.st8to10.uni.e.mean.padj",
                    "mPouVLOF.st8to10.uni.e.mean.padj",
                    "mPSLOF.st8to10.uni.e.mean.padj",
                    "mVegTLOF.st8to10.uni.e.mean.padj",
                    "WntLOF.st8to10.uni.e.mean.padj",
                    "NodalLOF.st8to10.uni.e.mean.padj",
                    "BmpLOF.st8to10.uni.e.mean.padj"
                    )

perc.padj.e <- merge( perc.e,padj.e, by="row.names", all.x=TRUE )
row.names( perc.padj.e ) <- perc.padj.e$Row.names
perc.padj.e <- perc.padj.e[ -1 ]


# Filter out 'low-read-count' genes using INTRON counts as follows:
# min. 10 reads at least at one stage over all controls from each experiment.
# Reads were DESeq2-normalised (but no regularised log shrinkage, which moderates fold change, see DESeq2 vignette, p. 17)

index.1                 <- apply( means.i.exp1, 1, function(x) any( x[ grep("uni",colnames(means.i.exp1)) ] >= 10 ) )
index.23a               <- apply( means.i.exp23, 1, function(x) any( x[ grep("uni",colnames(means.i.exp23)) ] >= 10 ) )
index.23b               <- apply( means.i.exp23, 1, function(x) any( x[ grep("dmso",colnames(means.i.exp23)) ] >= 10 ) )
index.4                 <- apply( means.i.exp4, 1, function(x) any( x[ grep("uni",colnames(means.i.exp4)) ] >= 10 ) )
index.5                 <- apply( means.i.exp5, 1, function(x) any( x[ grep("uni",colnames(means.i.exp5)) ] >= 10 ) )

index.u                 <- index.1 + index.23a + index.23b + index.4 + index.5

pm.i.allExp.m10         <- pm.i.allExp[ index.u==5, ]
pm.i.allExp.m10.m       <- as.matrix( pm.i.allExp.m10 )

# Introduce NAs into stage-specific differential expression analysis for each stage if any control is below 10 norm. reads and convert to NA at stage 8+ where alfa-amanitin >66.666%
# This is to remove unreliable gene expression data.
# Note: 'low-read-count' genes of experiment 2 were already removed in previous step

st8.ctrls               <- colnames( means.i.allExp)[c(grep("uni.st8p", colnames(means.i.allExp)), grep("dmso.st8p", colnames(means.i.allExp)))]
index.st8               <- apply( pm.i.allExp.m10, 1, function(x) any( x[st8.ctrls] < 10 ) )
pm.i.allExp.m10.m[ index.st8, grep( "st8p", colnames(perc.i.allExp) ) ] <- NA
index.st8.ama66         <- which( pm.i.allExp.m10[,"amanitin.st8p.uni.i.exp1"] > 2/3*100 )
pm.i.allExp.m10.m[ index.st8.ama66, grep( "st8p", colnames(perc.i.allExp) ) ] <- NA

st9.ctrls               <- colnames( means.i.allExp )[ c( grep( "uni.st9p", colnames(means.i.allExp) ), grep("dmso.st9p", colnames( means.i.allExp ))) ]
index.st9               <- apply( pm.i.allExp.m10, 1, function(x) any( x[st9.ctrls] < 10 ) )
pm.i.allExp.m10.m[ index.st9 , grep( "st9p", colnames(perc.i.allExp) ) ] <- NA

st10.ctrls              <- colnames(means.i.allExp)[ c(grep("uni.st10p", colnames(means.i.allExp)), grep("dmso.st10p", colnames(means.i.allExp))) ]
index.st10              <- apply( pm.i.allExp.m10, 1, function(x) any( x[st10.ctrls] < 10 ) )
pm.i.allExp.m10.m[ index.st10 , grep( "st10p", colnames(perc.i.allExp) ) ] <- NA


# Average INTRON expression changes
amanitin.st8to10.uni.i.mean <- apply(pm.i.allExp.m10.m[,grep("^amanitin",colnames(pm.i.allExp.m10.m[,1:27]))],1, mean, na.rm=TRUE )
cMO.st8to10.uni.i.mean <- apply(pm.i.allExp.m10.m[,grep("^cMO",colnames(pm.i.allExp.m10.m[,1:27]))],1,mean,na.rm=TRUE )
mPouVLOF.st8to10.uni.i.mean <- apply(pm.i.allExp.m10.m[,grep("^mPouVLOF",colnames(pm.i.allExp.m10.m[,1:27]))],1,mean,na.rm=TRUE )
mPSLOF.st8to10.uni.i.mean <- apply(pm.i.allExp.m10.m[,grep("^mPSLOF",colnames(pm.i.allExp.m10.m[,1:27]))],1,mean,na.rm=TRUE )
mVegTLOF.st8to10.uni.i.mean <- apply(pm.i.allExp.m10.m[,grep("^mVegTLOF",colnames(pm.i.allExp.m10.m[,1:27]))],1,mean,na.rm=TRUE )
WntLOF.st8to10.uni.i.mean <- apply(pm.i.allExp.m10.m[,grep("^WntLOF",colnames(pm.i.allExp.m10.m[,1:27]))],1,mean,na.rm=TRUE )
NodalLOF.st8to10.uni.i.mean <- apply(pm.i.allExp.m10.m[,grep("^NodalLOF",colnames(pm.i.allExp.m10.m[,1:27]))],1,mean,na.rm=TRUE )
BmpLOF.st8to10.uni.i.mean <- apply(pm.i.allExp.m10.m[,grep("^BmpLOF",colnames(pm.i.allExp.m10.m[,1:27]))],1,mean,na.rm=TRUE )

perc.i <- cbind(
                    amanitin.st8to10.uni.i.mean,
                    cMO.st8to10.uni.i.mean,
                    mPouVLOF.st8to10.uni.i.mean,
                    mPSLOF.st8to10.uni.i.mean,
                    mVegTLOF.st8to10.uni.i.mean,
                    WntLOF.st8to10.uni.i.mean,
                    NodalLOF.st8to10.uni.i.mean,
                    BmpLOF.st8to10.uni.i.mean
                )


# Select conservative FDR: cMO in experiment 2 and 5, use less significant FDR.
padj.i <- cbind(
                    res.dds.i.ama.exp1$padj,
                    ifelse(res.dds.i.cMO.exp23$padj >= res.dds.i.cMO.exp5$padj, res.dds.i.cMO.exp23$padj, res.dds.i.cMO.exp5$padj),
                    res.dds.i.mPouV.exp4$padj,
                    res.dds.i.mPS.exp4$padj,
                    res.dds.i.vegt.exp5$padj,
                    res.dds.i.wnt.exp23$padj,
                    res.dds.i.nodal.exp23$padj,
                    res.dds.i.bmp.exp23$padj
                )

row.names(padj.i) <- row.names(res.dds.i.ama.exp1)

colnames(padj.i) <- c(
                    "amanitin.st8to10.uni.i.padj",
                    "cMO.st8to10.uni.i.mean.padj",
                    "mPouVLOF.st8to10.uni.i.mean.padj",
                    "mPSLOF.st8to10.uni.i.mean.padj",
                    "mVegTLOF.st8to10.uni.i.mean.padj",
                    "WntLOF.st8to10.uni.i.mean.padj",
                    "NodalLOF.st8to10.uni.i.mean.padj",
                    "BmpLOF.st8to10.uni.i.mean.padj"
                    )

perc.padj.i <- merge( perc.i, padj.i, by="row.names", all.x=TRUE )
row.names( perc.padj.i ) <- perc.padj.i$Row.names
perc.padj.i <- perc.padj.i[ -1 ]

perc <- merge( as.data.frame( perc.padj.e ), as.data.frame( perc.padj.i ), by="row.names", all.x=TRUE )
row.names( perc ) <- perc$Row.names
perc <- perc[ -1 ]

# Remove rows where both exon and intron alfa-amanitin are NA
perc <- perc[ !is.na(perc$amanitin.st8to10.uni.e.mean) | !is.na(perc$amanitin.st8to10.uni.i.mean), ]

# Save whether INTRON (=0) or EXON (=1) counts are used for differential expression analysis
# Default EXON percentage, unless INTRON percentage <=50%, EXON percentage >50% and INTRON FDR <=10%
FC <- matrix( , nrow = nrow(perc), ncol = ncol(perc.padj.e)+1 )
row.names(FC) <- row.names(perc)
nc <- ncol(perc)
for(i in 1:nrow( perc )) {
            if ( perc[ i, (nc/2+1) ] <= 50 & perc[ i, 1 ] > 50 & perc[ i, (nc/4*3+1) ] <= 0.1 & !(  is.na(perc[ i,(nc/2+1) ])) | is.na(perc[i,1]) ) {
                FC[ i, 1:(nc/2) ] <- unlist(perc[ i, (nc/2+1):nc ])
                FC[ i, (nc/2+1) ] <- 0
            } else {
                FC[ i, 1:(nc/2) ] <- unlist( perc[ i, 1:(nc/2) ] )
                FC[ i, (nc/2+1) ] <- 1
            }
}

colnames(FC) <-
    c(
        "amanitin.st8to10.uni.mean",
        "cMO.st8to10.uni.mean",
        "mPouVLOF.st8to10.uni.mean",
        "mPSLOF.st8to10.uni.mean",
        "mVegTLOF.st8to10.uni.mean",
        "WntLOF.st8to10.uni.mean",
        "NodalLOF.st8to10.uni.mean",
        "BmpLOF.st8to10.uni.mean",
        "amanitin.st8to10.uni.mean.padj",
        "cMO.st8to10.uni.mean.padj",
        "mPouVLOF.st8to10.uni.mean.padj",
        "mPSLOF.st8to10.uni.mean.padj",
        "mVegTLOF.st8to10.uni.mean.padj",
        "WntLOF.st8to10.uni.mean.padj",
        "NodalLOF.st8to10.uni.mean.padj",
        "BmpLOF.st8to10.uni.mean.padj",
        "intron.or.exon"
    )


# Only select genes that are downregulated by alfa-amanitin (<=50%, FDR <=10%) averaged over alfa-amanitin-sensitive stages 8+, 9+ and 10+.
FC <- FC[ which( FC[ , "amanitin.st8to10.uni.mean" ] <= 50 & FC[ , "amanitin.st8to10.uni.mean.padj" ] <= 0.1 ), ]





#
# Spatial information
# ------------------------------------------------------------------------------------------

# Animal-vegetal (av) and dorso-ventral (dv) axes
# :EXON
spatio.e <- cbind(
                res.dds.e.av.spatio$baseMean,
                res.dds.e.av.spatio$log2FoldChange,
                res.dds.e.av.spatio$padj,
                res.dds.e.dv.spatio$log2FoldChange,
                res.dds.e.dv.spatio$padj
            )

row.names(spatio.e) <- row.names(res.dds.e.av.spatio)
colnames(spatio.e) <- c("baseMean.e","av.log2FC.e","av.padj.e","dv.log2FC.e","dv.padj.e")

# :INTRON
spatio.i <- cbind(
                res.dds.i.av.spatio$baseMean,
                res.dds.i.av.spatio$log2FoldChange,
                res.dds.i.av.spatio$padj,
                res.dds.i.dv.spatio$log2FoldChange,
                res.dds.i.dv.spatio$padj
            )

row.names(spatio.i) <- row.names(res.dds.i.av.spatio)
colnames(spatio.i) <- c("baseMean.i","av.log2FC.i","av.padj.i","dv.log2FC.i","dv.padj.i")

FC.spatio.e <- merge( FC, spatio.e, by="row.names", all.x=TRUE )
row.names( FC.spatio.e ) <- FC.spatio.e$Row.names
FC.spatio.e <- FC.spatio.e[-1]

FC.spatio <- merge( FC.spatio.e, spatio.i, by="row.names", all.x=TRUE )
row.names(FC.spatio) <- FC.spatio$Row.names
FC.spatio <- FC.spatio[-1]

# Default fold change from EXON counts
FC.spatio$av.pattern.G <- with( FC.spatio, ifelse( av.padj.e <= 0.1, av.log2FC.e, av.log2FC.e))
FC.spatio$dv.pattern.G <- with( FC.spatio, ifelse( dv.padj.e <= 0.1, dv.log2FC.e, dv.log2FC.e))

index.i.av <- with( FC.spatio, which( abs(av.log2FC.i) >= log2(2) & av.padj.i <= 0.1 ) )
index.i.dv <- with( FC.spatio, which( abs(dv.log2FC.i) >= log2(2) & dv.padj.i <= 0.1 ) )

# Switch to INTRON counts if FDR <=10% and absolute value of log2 INTRON count ratio greater than absolute value of log2 EXON count ratio
FC.spatio$av.pattern.G[index.i.av] <- with( FC.spatio[index.i.av,], ifelse( abs(av.log2FC.i) > abs(av.log2FC.e) & av.padj.i <= 0.1, av.log2FC.i, ifelse( av.padj.e <= 0.1, av.log2FC.e, av.log2FC.e )))
FC.spatio$dv.pattern.G[index.i.dv] <- with( FC.spatio[index.i.dv,], ifelse( abs(dv.log2FC.i) > abs(dv.log2FC.e) & dv.padj.i <= 0.1, dv.log2FC.i, ifelse( dv.padj.e <= 0.1, dv.log2FC.e, dv.log2FC.e )))







# Use previously published high-temporal resolution profiling of transcript levels
# ---------------------------------------------------------------------------------------

# Reference: Owens et al. (2016) Measuring Absolute RNA Copy Numbers at High Temporal Resolution Reveals Transcriptome Kinetics in Development


rzRNA           <- read.table( "./xenTro71/rzRNA_fpkm.txt.gz", header=T )

# Aggregate all gene isoforms and calculate FPKM sum per gene
rzRNAbyGene     <- aggregate( rzRNA[3:ncol(rzRNA)], by=rzRNA["Gene"], FUN=sum )
row.names( rzRNAbyGene ) <- rzRNAbyGene$Gene
rzRNAbyGene     <- rzRNAbyGene[-1]

# Re-order according to hpf
rzRNAbyGene     <- rzRNAbyGene[ , order( as.numeric( substr( names( rzRNAbyGene ), 2, nchar( names( rzRNAbyGene ))))) ]

# Relabel columns
colnames( rzRNAbyGene ) <- paste( "rzRNA_", substr( colnames( rzRNAbyGene ), 2, nchar( colnames( rzRNAbyGene ))), "hpf.fpkm", sep="" )

# Convert FPKM into TPM
rzRNA_sum       <- colSums( rzRNAbyGene )
rzRNAbyGene.tpm <- t( t( rzRNAbyGene ) / rzRNA_sum ) * 1e6
colnames( rzRNAbyGene.tpm ) <- paste( substr( colnames( rzRNAbyGene.tpm ), 1, nchar( colnames( rzRNAbyGene.tpm )) - 4 ), "tpm", sep="" )

tpm.df          <- as.data.frame(rzRNAbyGene.tpm[row.names(rzRNAbyGene.tpm) %in% row.names(FC.spatio),])

# maternal contribution (0 to 1 hpf)
tpm.df$maternal.0to1hpf.tpm <- apply(tpm.df[,1:3],1,mean)

# caculate mean of TPM between stage 8+ and 11 (5 to 9 hpf) from total RNA
tpm.df$meanSt8to11.tpm <- apply( tpm.df[, 11:19], 1, mean )

fc.all.tpm <- merge(FC.spatio,tpm.df[,49:50],by="row.names",all.x=T)
row.names(fc.all.tpm) <- fc.all.tpm$Row.names
fc.all.tpm <- fc.all.tpm[-1]

write.csv(fc.all.tpm[ , c(1:16,28:31) ],"zga_lof.csv")

#
#
# Generate graphs for publication
# -----------------------------------------------------------------
#
# Calculate minimal expression detected among Wnt, Nodal and Bmp LOF experiments
fc.all.tpm$min.WNBLOF.st8to10.uni.mean <- apply(fc.all.tpm[,c("WntLOF.st8to10.uni.mean", "NodalLOF.st8to10.uni.mean", "BmpLOF.st8to10.uni.mean")],1,min)
#
#
rn      <- row.names(fc.all.tpm)
tpm     <- log( fc.all.tpm$meanSt8to11.tpm +.1 )
mat     <- fc.all.tpm$maternal.0to1hpf.tpm
clm     <- col.mat [ as.numeric( cut( mat, breaks=brk.mat ))]
cls.av  <- col.spatio.1 [ as.numeric( cut( fc.all.tpm$av.pattern.G, breaks=brk.spatio )) ]
cls.dv  <- col.spatio.2 [ as.numeric( cut( fc.all.tpm$dv.pattern.G, breaks=brk.spatio )) ]
#
dat <- grep( "amanitin", colnames(fc.all.tpm) )
fc <- log( fc.all.tpm[ , dat[1] ] / 100 +.001 )
padj <- fc.all.tpm[, dat[2]]
#
t66 <- 66.666666
#
#
# Select genes to display in dotplots
sel.lb <- c( "mir427","foxh1.2","nodal3","sia2","otx2","chrd","ventx1.1","szl","foxi4.2","Xetro.E00391.id2","foxa1","hnf1b", "gata2", "eomes", "Xetro.E00956.t","fgf20","cdx4","Xetro.E00206.sox2","foxd3","zeb2" )
sel.rn <- sapply( sel.lb, function(x) grep(x,rn) )

# Ordered by mean expression level around MBT/early gastrula stage
select <- sel.rn[ order( fc.all.tpm$meanSt8to11.tpm[ sel.rn ]) ]


#
# Dotplot: Gene expression downregulated upon alfa-amanitin injection
# --------------------------------------------------------------------------

xl <- c( log(.1),log(10000) );  yl <- c( log(.001),log(1) )

graph.h <- "alfa-amanitin versus uninjected: maternal contribution (0-1 hpf)"
graph.x <- "RNA level (TPM) around MBT/early gastrula"
graph.y <- "% relative to control"

pdf( file="amanitin_mz_tpm.pdf", width=4, height=4, useDingbats=FALSE, family="Arial" )
plot( tpm, fc, axes=F, ann=F, pch=20, xlim=xl, ylim=yl, cex=1.2, col=alpha(clm, .7) )
points( tpm[select], fc[select], pch=21, cex=1.3, col="black", bg=alpha(clm[select], 1) )
text( tpm[select], fc[select], labels=paste0(1:length(select),":",sapply(rn[select], f4)), pos=4, offset=.3, cex=.9 )
abline( h=log(.5), col="grey50", lty="22" )
abline( h=log(1), col="grey50" )
title( main=graph.h, xlab=graph.x, ylab=graph.y, cex.main=.6 )
axis( 1, at=c( log(0.1),log(1),log(10),log(100),log(1000),log(10000) ), labels=c("0.1","1","10","100","1,000","10,000"))
axis( 2, at=c( log(.001),log(.01),log(.1),log(.5),log(1)), labels=c("0.001","0.01","0.1","0.5","1"))
dev.off()


# number of alfa-amanitin-sensitive genes
n.ama50 <- nrow(fc.all.tpm)

# retrieve genes with no maternal contribution (< 0.1 TPM, below detection limit)
n.zyg.only.ama50 <- sum( fc.all.tpm$maternal.0to1hpf.tpm < 0.1 )

# number of genes activated with maternal contribution (mat/zyg)
n.mat.zyg.ama50 <- n.ama50 - n.zyg.only.ama50


#
# Pie chart: Zygotic only and maternal/zygotic genes
# -----------------------------------------------------------------------------

genes.ama50 <- cbind(n.zyg.only.ama50, n.mat.zyg.ama50)
perc.zyg.only <- round( 100 / n.ama50 * n.zyg.only.ama50, 0 )
perc.mat.zyg <- round( 100 / n.ama50 * n.mat.zyg.ama50, 0 )
lab.zyg.only <- paste("zygotic only:\n", perc.zyg.only, "% (", n.zyg.only.ama50, ")", sep="")
lab.mat.zyg <- paste("maternal-zygotic:\n", perc.mat.zyg, "% (", n.mat.zyg.ama50, ")", sep="")
mz.col <- c("#EF4136", "#F7941E")

pdf( "zga_ama50.pdf", width=4, height=4 )
pie( genes.ama50, clockwise=T, radius=.6, labels=c(lab.zyg.only,lab.mat.zyg), col=mz.col, main=paste("polyA RNA reduction genes >50% upon a-amanitin\n n=", n.ama50, sep=""), cex.main=1 )
dev.off()

#
# Bargraph: LOF effects on ZGA
# ------------------------------------------------------------------------------

perc        <- fc.all.tpm[,1:8]

z           <- which( fc.all.tpm$maternal.0to1hpf.tpm < 0.1 )
mz          <- which( fc.all.tpm$maternal.0to1hpf.tpm >= 0.1 )
n.total.z   <- length(z)
n.total.mz  <- length(mz)

n.0to25     <- cbind( colSums( perc[z,] <= 25 ), colSums( perc[mz,] <= 25 ) )
n.25to67    <- cbind( colSums( perc[z,] <= 66.6666667 & perc[z,] > 25 ), colSums( perc[mz,] <= 66.6666667 & perc[mz,] > 25 ) )
n.67to150   <- cbind( colSums( perc[z,] > 66.6666667 & perc[z,] < 150 ), colSums( perc[mz,] > 66.6666667 & perc[mz,] < 150 ) )
n.150to400  <- cbind( colSums( perc[z,] >= 150 & perc[z,] < 400 ), colSums( perc[mz,] >= 150 & perc[mz,] < 400 ) )
n.400up     <- cbind( colSums( perc[z,] >= 400 ), colSums( perc[mz,] >= 400 ) )

n           <- cbind( n.400up, n.150to400, n.67to150, n.25to67,n.0to25 )
n.z         <- n[,c(1,3,5,7,9)]
n.mz        <- n[,c(2,4,6,8,10)]
n.z.mz      <- n.z + n.mz

p.z         <- 100 / n.total.z * n.z
p.mz        <- 100 / n.total.mz * n.mz
p.z.mz      <- 100 / (n.total.z + n.total.mz) * n.z.mz
p.summary   <- cbind(p.z.mz[,1]+p.z.mz[,2],p.z.mz[,3],p.z.mz[,4]+p.z.mz[,5])

p           <- rbind(p.z,p.mz)

row.names(p) <- sapply(row.names(p), f5)

# Graph: LOF effects on ZGA
pdf( "zga_ama50_misregulation.pdf", width=8, height=4 )
barplot( as.matrix(t(p.z.mz)),col=c("#DD8204","#E18F2A","#EBEBEB","#AFBADA","#3953A4"),las=2 )
dev.off()

# Graph: LOF effects on ZGA (separating zygotic only and maternal-zygotic genes)
# pdf( "zga_ama50_mz_misregulation.pdf", width=8, height=4 )
# barplot( as.matrix(t(p[c(1,9,2,10,3,11,4,12,5,13,6,14,7,15,8,16),])),col=c("#DD8204","#E18F2A","#EBEBEB","#AFBADA","#3953A4"),las=2 )
# dev.off()

# --------------------------------------------------------------------------


#
# Graph: mPouV KD vs mPS KD
# -------------------------------------------------------------

dat.x   <- grep( "mPouVLOF", colnames(fc.all.tpm) )
dat.y   <- grep( "mPSLOF", colnames(fc.all.tpm) )
fc.x    <- fc.all.tpm[ , dat.x[1] ]
fc.y    <- fc.all.tpm[ , dat.y[1] ]

fc.d    <- fc.x/fc.y; fc.d[fc.d < 0.25] <- 0.25; fc.d[fc.d > 4] <- 4
fc.log  <- log(fc.d)
cl.fc   <- col.fc [ as.numeric( cut( fc.log, breaks=brk.fc )) ]

#
# FC ratios between mPouV LOF/ctrl vs mPouV/Sox3 (mPS) LOF/ctrl
# -------------------------------------------------------------

plotDifSp( fc.x=fc.x, fc.y=fc.y, select=select, cl=cl.fc, x.lim.min=0, x.lim.max=1000, y.lim.min=0, y.lim.max=1000, pdf.name="mPouVLOF_vs_mPSLOF_FCLog_noLab.pdf", log=TRUE, height=4.5)
plotDifSp( fc.x=fc.x, fc.y=fc.y, select=select, cl=cl.fc, x.lim.min=0, x.lim.max=1000, y.lim.min=0, y.lim.max=1000, pdf.name="mPouVLOF_vs_mPSLOF_FCLog.pdf", log=TRUE, label=TRUE, height=4.5)

pdf("FC_ratios_bar.pdf", width=2, height=10)
image( 1, brk.fc, matrix(1:50,nrow=1), col=col.fc, axes=F, ylim=c(log(.25),log(4)), xlab="",ylab="FC")
axis(2, at=c(log(.25),log(1),log(4)), labels=c("0.25","1","4"))
dev.off()


#
# Maternal TFs (mPouV/Sox3 [mPS] and mVegT) versus signaling
# ----------------------------------------------------------------------------------------------------------------------------
#
z <- matrix(1:100,nrow=1)
#
# generate pdf file for ANIMAL:VEGETAL (AnVg) and DORSO-VENTRAL spatial expression analysis
pdf("AnVg_log2_colorBar.pdf", width=2, height=10)
par(mar=c(0,0,0,0) + 4)
image(1,brk.spatio,z,col=col.spatio.1,axes=F,ylim=c(-10,10), xlab="",ylab="fold enrichment")
axis(2, at=log2( c( .001, .1, .25, .333, .5, 1, 2, 3, 4, 10, 1000 )), labels=c( .001, .1, .25, "", .5, 1, 2, "", 4, 10, 1000 ))
dev.off()
#
pdf("DV_log2_colorBar.pdf", width=2, height=10)
par(mar=c(0,0,0,0) + 4)
image(1,brk.spatio,z,col=col.spatio.2,axes=F,ylim=c(-10,10), xlab="",ylab="fold enrichment")
axis(2, at=log2( c( .001, .1, .25, .333, .5, 1, 2, 3, 4, 10, 1000 )), labels=c( .001, .1, .25, "", .5, 1, 2, "", 4, 10, 1000 ))
dev.off()
#
#
#
dat.x           <- grep( "mPSLOF", colnames(fc.all.tpm) )
dat.y           <- grep( "NodalLOF", colnames(fc.all.tpm) )
dat.z           <- grep( "mVegT", colnames(fc.all.tpm) )
fc.x            <- fc.all.tpm[ , dat.x[1] ]
fc.y            <- fc.all.tpm[ , dat.y[1] ]
fc.z            <- fc.all.tpm[ , dat.z[1] ]
#
#
# ANIMAL-VEGETAL expression: mTF LOF vs Nodal LOF
# --------------------------------------------------------------------------------------
#
plotDifSp( fc.x=fc.x, fc.y=fc.y, select=select, cl=cls.av, x.lim.min=0, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mPSLOF_vs_NodalLOF_spatioAnVgZoom_noLab.pdf", height=4.5)
plotDifSp( fc.x=fc.x, fc.y=fc.y, select=select, cl=cls.av, x.lim.min=0, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mPSLOF_vs_NodalLOF_spatioAnVgZoom.pdf", label=TRUE, height=4.5)
#
plotDifSp( fc.x=fc.z, fc.y=fc.y, select=select, cl=cls.av, x.lim.min=20, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mVegTLOF_vs_NodalLOF_spatioAnVgZoom_noLab.pdf", height=4.5)
plotDifSp( fc.x=fc.z, fc.y=fc.y, select=select, cl=cls.av, x.lim.min=20, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mVegTLOF_vs_NodalLOF_spatioAnVgZoom.pdf", label=TRUE, height=4.5)
#
#
# DORSO-VENTRAL expression: mTF LOF vs Nodal LOF
# --------------------------------------------------------------------------------
#
plotDifSp( fc.x=fc.x, fc.y=fc.y, select=select, cl=cls.dv, x.lim.min=0, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mPSLOF_vs_NodalLOF_spatioDVZoom_noLab.pdf", height=4.5)
plotDifSp( fc.x=fc.x, fc.y=fc.y, select=select, cl=cls.dv, x.lim.min=0, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mPSLOF_vs_NodalLOF_spatioDVZoom.pdf", label=TRUE, height=4.5)
#
plotDifSp( fc.x=fc.z, fc.y=fc.y, select=select, cl=cls.dv, x.lim.min=20, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mVegTLOF_vs_NodalLOF_spatioDVZoom_noLab.pdf", height=4.5)
plotDifSp( fc.x=fc.z, fc.y=fc.y, select=select, cl=cls.dv, x.lim.min=20, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mVegTLOF_vs_NodalLOF_spatioDVZoom.pdf", label=TRUE, height=4.5)
#
#
dat.y           <- grep( "WntLOF", colnames(fc.all.tpm) )
fc.y            <- fc.all.tpm[ , dat.y[1] ]
#
# ANIMAL-VEGETAL expression: mTF LOF vs Wnt LOF
# --------------------------------------------------------------------------------------
#
plotDifSp( fc.x=fc.x, fc.y=fc.y, select=select, cl=cls.av, x.lim.min=0, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mPSLOF_vs_WntLOF_spatioAnVgZoom_noLab.pdf", height=4.5)
plotDifSp( fc.x=fc.x, fc.y=fc.y, select=select, cl=cls.av, x.lim.min=0, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mPSLOF_vs_WntLOF_spatioAnVgZoom.pdf", label=TRUE, height=4.5)
#
plotDifSp( fc.x=fc.z, fc.y=fc.y, select=select, cl=cls.av, x.lim.min=20, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mVegTLOF_vs_WntLOF_spatioAnVgZoom_noLab.pdf", height=4.5)
plotDifSp( fc.x=fc.z, fc.y=fc.y, select=select, cl=cls.av, x.lim.min=20, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mVegTLOF_vs_WntLOF_spatioAnVgZoom.pdf", label=TRUE, height=4.5)
#
#
# DORSO-VENTRAL expression: mTF LOF vs Wnt LOF
# --------------------------------------------------------------------------------------
#
plotDifSp( fc.x=fc.x, fc.y=fc.y, select=select, cl=cls.dv, x.lim.min=0, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mPSLOF_vs_WntLOF_spatioDVZoom_noLab.pdf", height=4.5)
plotDifSp( fc.x=fc.x, fc.y=fc.y, select=select, cl=cls.dv, x.lim.min=0, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mPSLOF_vs_WntLOF_spatioDVZoom.pdf", label=TRUE, height=4.5)
#
plotDifSp( fc.x=fc.z, fc.y=fc.y, select=select, cl=cls.dv, x.lim.min=20, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mVegTLOF_vs_WntLOF_spatioDVZoom_noLab.pdf", height=4.5)
plotDifSp( fc.x=fc.z, fc.y=fc.y, select=select, cl=cls.dv, x.lim.min=20, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mVegTLOF_vs_WntLOF_spatioDVZoom.pdf", label=TRUE, height=4.5)
#
#
dat.y       <- grep( "BmpLOF", colnames(fc.all.tpm) )
fc.y        <- fc.all.tpm[ , dat.y[1] ]
#
# ANIMAL-VEGETAL expression: mTF LOF vs Bmp LOF
# --------------------------------------------------------------------------------------
#
plotDifSp( fc.x=fc.x, fc.y=fc.y, select=select, cl=cls.av, x.lim.min=0, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mPSLOF_vs_BmpLOF_spatioAnVgZoom_noLab.pdf", height=4.5)
plotDifSp( fc.x=fc.x, fc.y=fc.y, select=select, cl=cls.av, x.lim.min=0, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mPSLOF_vs_BmpLOF_spatioAnVgZoom.pdf", label=TRUE, height=4.5)
#
plotDifSp( fc.x=fc.z, fc.y=fc.y, select=select, cl=cls.av, x.lim.min=20, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mVegTLOF_vs_BmpLOF_spatioAnVgZoom_noLab.pdf", height=4.5)
plotDifSp( fc.x=fc.z, fc.y=fc.y, select=select, cl=cls.av, x.lim.min=20, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mVegTLOF_vs_BmpLOF_spatioAnVgZoom.pdf", label=TRUE, height=4.5)
#
#
# DORSO-VENTRAL expression: mTF LOF vs Bmp LOF
# --------------------------------------------------------------------------------
#
plotDifSp( fc.x=fc.x, fc.y=fc.y, select=select, cl=cls.dv, x.lim.min=0, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mPSLOF_vs_BmpLOF_spatioDVZoom_noLab.pdf", height=4.5)
plotDifSp( fc.x=fc.x, fc.y=fc.y, select=select, cl=cls.dv, x.lim.min=0, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mPSLOF_vs_BmpLOF_spatioDVZoom.pdf", label=TRUE, height=4.5)
#
plotDifSp( fc.x=fc.z, fc.y=fc.y, select=select, cl=cls.dv, x.lim.min=20, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mVegTLOF_vs_BmpLOF_spatioDVZoom_noLab.pdf", height=4.5)
plotDifSp( fc.x=fc.z, fc.y=fc.y, select=select, cl=cls.dv, x.lim.min=20, x.lim.max=200, y.lim.min=0, y.lim.max=200, pdf.name="mVegTLOF_vs_BmpLOF_spatioDVZoom.pdf", label=TRUE, height=4.5)
#
#
#
# SUMMARY: LOF effect on spatial information
# --------------------------------------------------------------------------------

spatialReadout <- function( df, file.txt, spatioFC ) {
    
    fc.spatio.tpm   <- fc.all.tpm[ !is.na( fc.all.tpm$av.pattern.G ) & !is.na(fc.all.tpm$dv.pattern.G), ]
    n               <- nrow( fc.spatio.tpm )
    n0              <- nrow( fc.all.tpm )
    
    n.an.total      <- length( which( fc.spatio.tpm[ , "av.pattern.G" ] >= log2( spatioFC ) ) )
    n.vg.total      <- length( which( fc.spatio.tpm[ , "av.pattern.G" ] <= -log2( spatioFC ) ) )
    n.ub.av.total   <- n - n.an.total - n.vg.total
    p.av.total      <- 100 / n0 * c( n.an.total, n.ub.av.total, n.vg.total, n )
    
    n.d.total       <- length( which( fc.spatio.tpm[ , "dv.pattern.G" ] >= log2( spatioFC ) ) )
    n.v.total       <- length( which( fc.spatio.tpm[ , "dv.pattern.G" ] <= -log2( spatioFC ) ) )
    n.ub.dv.total   <- n - n.d.total - n.v.total
    p.dv.total      <- 100 / n0 * c( n.d.total, n.ub.dv.total, n.v.total, n )
    
    spatio.df       <- df[ !is.na( df$av.pattern.G ) & !is.na(df$dv.pattern.G), ]
    n.spatio.df     <- nrow( spatio.df )
    n.df            <- nrow( df )
    
    n.an            <- length( which( spatio.df[ , "av.pattern.G" ] >= log2( spatioFC ) ) )
    n.vg            <- length( which( spatio.df[ , "av.pattern.G" ] <= -log2( spatioFC ) ) )
    n.ub.av         <- n.spatio.df - n.an - n.vg
    p.av            <- 100 / c( n.an.total, n.ub.av.total, n.vg.total, n0 ) * c( n.an, n.ub.av, n.vg, n.spatio.df )
    
    n.d             <- length( which( spatio.df[ , "dv.pattern.G" ] >= log2( spatioFC ) ) )
    n.v             <- length( which( spatio.df[ , "dv.pattern.G" ] <= -log2( spatioFC ) ) )
    n.ub.dv         <- n.spatio.df - n.d - n.v
    p.dv            <- 100 / c( n.d.total, n.ub.dv.total, n.v.total, n0 ) * c( n.d, n.ub.dv, n.v, n.spatio.df )
    
    # generate file with statitics
    sink( file = file.txt )
    
    cat( "## STATISTICS: LOF EFFECTS ON SPATIAL INFORMATION ##\n\n", sep="" )
    cat( "fold change threshold for both An:Vg and D:V axes:", spatioFC, "\n\n", sep="\t" )
    
    cat( "## TOTAL GENES: #\n" )
    cat( "animal", "ubiquitous", "vegetal", "total\n", sep="\t" )
    cat( n.an.total, n.ub.av.total, n.vg.total, n0, "\n", sep="\t" )
    cat( "dorsal", "ubiquitous", "ventral", "total\n", sep="\t" )
    cat( n.d.total, n.ub.dv.total, n.v.total, n0, "\n\n", sep="\t" )
    
    cat( "## TOTAL GENES: PERCENTAGE\n" )
    cat( "animal", "ubiquitous", "vegetal", "total\n", sep="\t" )
    cat( p.av.total, "\n", sep="\t" )
    cat( "dorsal", "ubiquitous", "ventral", "total\n", sep="\t" )
    cat( p.dv.total, "\n\n", sep="\t" )
    
    cat( "## MISREGULATED GENES: #\n" )
    cat( "animal", "ubiquitous", "vegetal", "total\n", sep="\t" )
    cat( n.an, n.ub.av, n.vg, n.df, "\n", sep="\t" )
    cat( "dorsal", "ubiquitous", "ventral", "total\n", sep="\t" )
    cat( n.d, n.ub.dv, n.v, n.df, "\n\n", sep="\t" )
    
    cat( "## MISREGULATED GENES: PERCENTAGE\n" )
    cat( "animal", "ubiquitous", "vegetal", "total\n", sep="\t" )
    cat( p.av, "\n\n", sep="\t" )
    cat( "dorsal", "ubiquitous", "ventral", "total\n", sep="\t" )
    cat( p.dv, "\n\n", sep="\t" )
    
    sink()

    return( list( c( n.an, n.ub.av, n.vg, n.df ), c( n.d, n.ub.dv, n.v, n.df ) ) )
}

dwn.amanitin        <- fc.all.tpm[ which( fc.all.tpm[ ,1 ] <= t66 ), ]
dwn.cMO             <- fc.all.tpm[ which( fc.all.tpm[ ,2 ] <= t66 ), ]
dwn.mTF.LOF         <- fc.all.tpm[ which( apply( fc.all.tpm[ , 3:5 ], 1, function(x) any( x <= t66 ))),]
dwn.Signals.LOF     <- fc.all.tpm[ which( apply( fc.all.tpm[ , 6:8 ], 1, function(x) any( x <= t66 ))),]
dwn.ALL.LOF         <- fc.all.tpm[ which( apply( fc.all.tpm[ , 3:8 ], 1, function(x) any( x <= t66 ))),]
dwn.mPouV.LOF       <- fc.all.tpm[ which( fc.all.tpm[ ,3 ] <= t66 ), ]
dwn.mPS.LOF         <- fc.all.tpm[ which( fc.all.tpm[ ,4 ] <= t66 ), ]
dwn.mVegT.LOF       <- fc.all.tpm[ which( fc.all.tpm[ ,5 ] <= t66 ), ]
dwn.Wnt.LOF         <- fc.all.tpm[ which( fc.all.tpm[ ,6 ] <= t66 ), ]
dwn.Nodal.LOF       <- fc.all.tpm[ which( fc.all.tpm[ ,7 ] <= t66 ), ]
dwn.Bmp.LOF         <- fc.all.tpm[ which( fc.all.tpm[ ,8 ] <= t66 ), ]

s1                  <- spatialReadout( dwn.amanitin, "spatial_down_amanitin.txt", 3 )
s2                  <- spatialReadout( dwn.cMO, "spatial_down_cMO.txt", 3 )
s3                  <- spatialReadout( dwn.mPouV.LOF, "spatial_down_mPouV_LOF.txt", 3 )
s4                  <- spatialReadout( dwn.mPS.LOF, "spatial_down_mPS_LOF.txt", 3 )
s5                  <- spatialReadout( dwn.mVegT.LOF, "spatial_down_mVegT_LOF.txt", 3 )
s6                  <- spatialReadout( dwn.Wnt.LOF, "spatial_down_Wnt_LOF.txt", 3 )
s7                  <- spatialReadout( dwn.Nodal.LOF, "spatial_down_Nodal_LOF.txt", 3 )
s8                  <- spatialReadout( dwn.Bmp.LOF, "spatial_down_Bmp_LOF.txt", 3 )
s9                  <- spatialReadout( dwn.mTF.LOF, "spatial_down_mTF_LOF.txt", 3 )
s10                 <- spatialReadout( dwn.Signals.LOF, "spatial_down_Signals_LOF.txt", 3 )
s11                 <- spatialReadout( dwn.ALL.LOF, "spatial_down_ALL_LOF.txt", 3 )


s.an.vg             <- rbind( s1[[1]], s2[[1]], s3[[1]], s4[[1]], s5[[1]], s6[[1]], s7[[1]], s8[[1]], s9[[1]], s10[[1]], s11[[1]] )
s.d.v               <- rbind( s1[[2]], s2[[2]], s3[[2]], s4[[2]], s5[[2]], s6[[2]], s7[[2]], s8[[2]], s9[[2]], s10[[2]], s11[[2]] )

label.lof           <- c( "amanitin", "cMO", "mPouV.LOF", "mPS.LOF", "mVegT.LOF", "Wnt.LOF", "Nodal.LOF", "Bmp.LOF", "mTF.LOF", "Signal.LOF", "mTF.Signal.LOF")
label.an.vg         <- c( "animal", "uniform", "vegetal", "total" )
label.d.v           <- c( "dorsal", "uniform", "ventral", "total" )

rownames(s.an.vg)   <- label.lof
colnames(s.an.vg)   <- label.an.vg

rownames(s.d.v)     <- label.lof
colnames(s.d.v)     <- label.d.v

n.an.vg             <- t(s.an.vg)
n.d.v               <- t(s.d.v)

p.an.vg             <- 100 / t(s.an.vg)[,1] * t(s.an.vg)
p.d.v               <- 100 / t(s.d.v)[,1] * t(s.d.v)

pdf( "zga_lof_spatial_AnVg_perc.pdf", height=4, width=10 )
coord <- barplot( p.an.vg, main = "LOF-mediated downregulation of ZGA across the animal-vegetal axis", xlab = NA, ylab = "percentage",
                  col=c( "#A5D4E3", "#231F20", "#F1E91F", "#EC7D30" ), border = NA, legend = rownames(p.an.vg),
                  args.legend = list( x = 12, y = max( p.an.vg ) + 2, bty = "n"), beside = TRUE, las = 2, cex.names=.7 )
                  text( x = coord+.4, y = p.an.vg+.4, label = n.an.vg, pos = 3, cex = .8, col = "black", srt = 90, xpd=TRUE )
dev.off()

pdf( "zga_lof_spatial_DV_perc.pdf", height=4, width=10 )
coord <- barplot( p.d.v, main = "LOF-mediated downregulation of ZGA across the dorso-ventral axis", xlab = NA, ylab = "percentage",
                  col=c("#F6A3A9","#231F20","#739CD2","#EC7D30"), border = NA, legend = rownames(p.d.v),
                  args.legend = list( x = 12, y = max( p.d.v ) + 2, bty = "n"), beside = TRUE, las = 2, cex.names=.7 )
                  text( x = coord+.4, y = p.d.v+.4, label = n.d.v, pos = 3, cex = .8, col = "black", srt = 90, xpd=TRUE )
dev.off()

#
#
# GO TERM ANALYSIS
# -------------------------------------------------------------------------
#
# "GO_annot.txt" updated version of GO term annotation file including EST assemblies generated by BLAST2GO
# reference: Conesa, A. et al. (2005) Blast2GO: a universal tool for annotation, v1isualization and analysis in functional genomics research.
#
go.table               <- read.table( "./xenTro71/GO_annot.txt.gz", header=FALSE )
colnames(go.table)     <- c( "go.term", "go.evidence", "gene.id" )
#
goFrame                <- GOFrame( go.table, organism="Homo sapiens" )
goAllFrame             <- GOAllFrame( goFrame )
gsc                    <- GeneSetCollection( goAllFrame, setType=GOCollection() )
#
universe               <- row.names(fc.all.tpm)
#
#
# Downregulated genes
mPouVLOF.dwn.genes     <- row.names( fc.all.tpm[ which( fc.all.tpm[ ,"mPouVLOF.st8to10.uni.mean" ] <= t66 ), ] )
mPSLOF.dwn.genes       <- row.names( fc.all.tpm[ which( fc.all.tpm[ ,"mPSLOF.st8to10.uni.mean" ] <= t66 ), ] )
mVegTLOF.dwn.genes     <- row.names( fc.all.tpm[ which( fc.all.tpm[ ,"mVegTLOF.st8to10.uni.mean" ] <= t66 ), ] )
WntLOF.dwn.genes       <- row.names( fc.all.tpm[ which( fc.all.tpm[ ,"WntLOF.st8to10.uni.mean" ] <= t66 ), ] )
NodalLOF.dwn.genes     <- row.names( fc.all.tpm[ which( fc.all.tpm[ ,"NodalLOF.st8to10.uni.mean" ] <= t66 ), ] )
BmpLOF.dwn.genes       <- row.names( fc.all.tpm[ which( fc.all.tpm[ ,"BmpLOF.st8to10.uni.mean" ] <= t66 ), ] )
#
min.WNBLOF.dwn.genes   <- row.names( fc.all.tpm[ which( fc.all.tpm[ ,"min.WNBLOF.st8to10.uni.mean" ] <= t66 ), ] )
#
# upregulated genes
mPouVLOF.up.genes      <- row.names( fc.all.tpm[ which( fc.all.tpm[ ,"mPouVLOF.st8to10.uni.mean" ] >= 150 ), ] )
mPSLOF.up.genes        <- row.names( fc.all.tpm[ which( fc.all.tpm[ ,"mPSLOF.st8to10.uni.mean" ] >= 150 ), ] )
mVegTLOF.up.genes      <- row.names( fc.all.tpm[ which( fc.all.tpm[ ,"mVegTLOF.st8to10.uni.mean" ] >= 150 ), ] )
WntLOF.up.genes        <- row.names( fc.all.tpm[ which( fc.all.tpm[ ,"WntLOF.st8to10.uni.mean" ] >= 150 ), ] )
NodalLOF.up.genes      <- row.names( fc.all.tpm[ which( fc.all.tpm[ ,"NodalLOF.st8to10.uni.mean" ] >= 150 ), ] )
BmpLOF.up.genes        <- row.names( fc.all.tpm[ which( fc.all.tpm[ ,"BmpLOF.st8to10.uni.mean" ] >= 150 ), ] )
#
#
#
go.analysis <- function(gene.select) {
    
    params <- GSEAGOHyperGParams(
    
            name="GSEA",
            geneSetCollection=gsc,
            geneIds=gene.select,
            universeGeneIds=universe,
            ontology="BP",
            pvalueCutoff=1,
            conditional=FALSE,
            testDirection="over"
            )
            
    over <- hyperGTest(params)
    return(over)
}
#
#
mP.dwn                  <- go.analysis(mPouVLOF.dwn.genes)
mPS.dwn                 <- go.analysis(mPSLOF.dwn.genes)
mV.dwn                  <- go.analysis(mVegTLOF.dwn.genes)
W.dwn                   <- go.analysis(WntLOF.dwn.genes)
N.dwn                   <- go.analysis(NodalLOF.dwn.genes)
B.dwn                   <- go.analysis(BmpLOF.dwn.genes)


t1                      <- summary(mP.dwn)
t2                      <- summary(mPS.dwn)
t3                      <- summary(mV.dwn)
t4                      <- summary(W.dwn)
t5                      <- summary(N.dwn)
t6                      <- summary(B.dwn)


mP.up                   <- go.analysis(mPouVLOF.up.genes)
mPS.up                  <- go.analysis(mPSLOF.up.genes)
mV.up                   <- go.analysis(mVegTLOF.up.genes)
W.up                    <- go.analysis(WntLOF.up.genes)
N.up                    <- go.analysis(NodalLOF.up.genes)
B.up                    <- go.analysis(BmpLOF.up.genes)

t1u                     <- summary(mP.up)
t2u                     <- summary(mPS.up)
t3u                     <- summary(mV.up)
t4u                     <- summary(W.up)
t5u                     <- summary(N.up)
t6u                     <- summary(B.up)

lt                      <- list( t1, t2, t3, t4, t5, t6, t1u, t2u, t3u, t4u, t5u, t6u )

gobpid                  <- unlist( sapply( lt, function(x) x[,1]) )
size                    <- unlist( sapply( lt, function(x) x[,6]) )
term                    <- unlist( sapply( lt, function(x) x[,7]) )
t0                      <- unique( data.frame( GOBPID=gobpid, Size=size, Term=term ) )

lt0                     <- list( t0, t1, t2, t3, t4, t5, t6, t1u, t2u, t3u, t4u, t5u, t6u )

t                       <- plyr::join_all( lt0, by='GOBPID', type="left" )

p.col                   <- grep( "Pvalue", colnames(t) )
c.col                   <- grep( "^Count", colnames(t) )
t                       <- t[ , c(1:3, p.col, c.col ) ]
row.names(t)            <- t$GOBPID
t                       <- t[-1]

go.terms                <- unique( c(
                            "GO:1900107", "GO:0060061", "GO:0048513", "GO:0006334", "GO:1901362", "GO:0016071", "GO:0006412",
                            "GO:0043687", "GO:0006508", "GO:0016477", "GO:0031589", "GO:0007369", "GO:0001843", "GO:0001704", "GO:0007398", "GO:0007498", "GO:0007492",
                            "GO:0007389", "GO:0009953", "GO:0009952", "GO:0003002",
                            "GO:0007417", "GO:0072358", "GO:0061053",
                            "GO:0006996", "GO:0051276", "GO:0019953", "GO:0007276", "GO:0006305", "GO:0006306" )
                            )
sel.t                   <- t[ rev(go.terms), ]

# select enriched GO terms of interest:
go.blot.circ            <- sel.t[,3:14]
go.blot.circ[is.na(go.blot.circ)] <- 1
go.blot.circ            <- -log10(go.blot.circ)
go.blot.circ            <- sqrt(go.blot.circ/pi)

go.blot.col             <- sel.t[,15:26]
go.blot.col[is.na(go.blot.col)] <- 0

sector                  <- c("A:0-1","B:2-10","C:11-50","D:51-100","E:101-500")
terms                   <- paste0( sel.t$Term, " (", row.names(sel.t), ")" )

grid                    <- data.frame( matrix( 0, ncol = 4, nrow = length( go.terms ) * ncol( go.blot.circ ) ) )
conds                   <- c("mPkd.dwn","mPSLOF.dwn","mVegTLOF.dwn","WntLOF.dwn","NodalLOF.dwn","BMPLOF.dwn",
                            "mPkd.up","mPSLOF.up","mVegTLOF.up","WntLOF.up","NodalLOF.up","BMPLOF.up")

grid$X1                 <- rep( conds, each = length( go.terms ) )
grid$X2                 <- terms
grid$X3                 <- unlist(list(go.blot.circ))
grid$X4                 <- sector[ as.numeric( cut( unlist( list( go.blot.col ) ), breaks=c( -1, 1, 10, 50, 100, 500 ) ) ) ]

# factorise X1 (condition) and X2 (GO) to ensure order does not change when calling ggplot()
grid$X5                 <- factor( grid$X1, conds )
grid$X6                 <- factor( grid$X2, terms )

pdf("GO_term_LOF_ZGA.pdf", height=8.5, width=8.5, useDingbats=FALSE)
    ggplot( grid, aes( X5, X6, color = factor( X4 ) ) ) +
    geom_point( aes( size = X3^2*pi ) ) +
    scale_colour_manual( values = c( "white", "#DBCEA2", "#A69B71", "#7E7561", "#48475F" ) ) +
    scale_size_continuous( range = range( c(0.1,10) ) ) +
    scale_shape_manual( values = 21 ) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()


# VENN diagrams of downregulated genes
# ----------------------------------------------------------------------------

venn.dwn               <- matrix( 0, nrow = nrow( fc.all.tpm ), ncol = 4)
colnames( venn.dwn )   <- c( "WntLOF.dwn", "NodalLOF.dwn", "BmpLOF.dwn", "mPSLOF.dwn" )

for (i in 1:nrow( venn.dwn ))
{
    venn.dwn[i,1]      <- rownames(fc.all.tpm)[i] %in% WntLOF.dwn.genes
    venn.dwn[i,2]      <- rownames(fc.all.tpm)[i] %in% NodalLOF.dwn.genes
    venn.dwn[i,3]      <- rownames(fc.all.tpm)[i] %in% BmpLOF.dwn.genes
    venn.dwn[i,4]      <- rownames(fc.all.tpm)[i] %in% mPSLOF.dwn.genes
}
pdf("venn_WNB_mPS_LOF.pdf")
vennDiagram( vennCounts( venn.dwn ))
dev.off()


colnames( venn.dwn )   <- c( "WntLOF.dwn", "NodalLOF.dwn", "BmpLOF.dwn", "mVegTLOF.dwn" )
for (i in 1:nrow( venn.dwn ))
{
    venn.dwn[i,1]      <- rownames(fc.all.tpm)[i] %in% WntLOF.dwn.genes
    venn.dwn[i,2]      <- rownames(fc.all.tpm)[i] %in% NodalLOF.dwn.genes
    venn.dwn[i,3]      <- rownames(fc.all.tpm)[i] %in% BmpLOF.dwn.genes
    venn.dwn[i,4]      <- rownames(fc.all.tpm)[i] %in% mVegTLOF.dwn.genes
}
pdf("venn_WNB_mVegT_LOF.pdf")
vennDiagram( vennCounts( venn.dwn ))
dev.off()


venn.dwn               <- matrix( 0, nrow = nrow( fc.all.tpm ), ncol = 3)
colnames( venn.dwn )   <- c( "mPSLOF.dwn", "mVegTLOF.dwn", "minWNBLOF.dwn" )
for (i in 1:nrow( venn.dwn ))
{
    venn.dwn[i,1]      <- rownames(fc.all.tpm)[i] %in% mPSLOF.dwn.genes
    venn.dwn[i,2]      <- rownames(fc.all.tpm)[i] %in% mVegTLOF.dwn.genes
    venn.dwn[i,3]      <- rownames(fc.all.tpm)[i] %in% min.WNBLOF.dwn.genes
}
pdf("venn_mPS_mVegT_minWNB_LOF.pdf")
vennDiagram( vennCounts( venn.dwn ))
dev.off()

# ------------------------------------------------------------------------

