#Load the libraries needed for the following tasks.
library(ALL)
library(biomaRt)
library(minfiData)
library(minfi)
library(GEOquery)
library(airway)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(AnnotationHub)

#Question 1.
#What is the mean expression across all features for sample 5 in the ALL dataset (from the ALL package)?

#Download the ALL dataset.
data(ALL)

#Pick all the features of sample 5, extract their expression levels, and average them.
mean(exprs(ALL[,5]))

#Question 2.
#We will use the biomaRt package to annotate an Affymetrix microarray. We want our results in the hg19 build of the human genome and we therefore need to connect to Ensembl 75 which is the latest release on this genome version. How to connect to older versions of Ensembl is described in the biomaRt package vignette; it can be achived with the command mart <- useMart(host='feb2014.archive.ensembl.org', biomart = "ENSEMBL_MART_ENSEMBL").
#Using this version of Ensembl, annotate each feature of the ALL dataset with the Ensembl gene id.
#How many probesets (features) are annotated with more than one Ensembl gene id?

#Load the specified version of Ensembl, specifically the Ensembl gene id numbers.
mart <- useMart(host='feb2014.archive.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL')
Ensemble75 <- useDataset('hsapiens_gene_ensembl', mart)

#Extract the features in the ALL dataset.
Names_features <- featureNames(ALL)

#Annotate the extracted features using the Ensembl gene id numbers.
ALLannotated <- getBM(attributes=c("ensembl_gene_id","affy_hg_u95av2"), filters="affy_hg_u95av2", values=Names_features, mart=Ensemble75)

#Count the number of features with more than one id.
sum(table(ALLannotated[,2])>1)

#Question 3.
#How many probesets (Affymetrix IDs) are annotated with one or more genes on the autosomes (chromosomes 1 to 22)?

#Filter the extracted features and keep those on the autosomes only, and annotate them using the Ensembl gene id numbers.
ALLannotated_chr <- getBM(attributes=c("ensembl_gene_id", "affy_hg_u95av2", "chromosome_name"), filters=c("affy_hg_u95av2","chromosome_name"), values=list(Names_features, c(1:22)), mart=Ensemble75)

#Count the number of features with at least one id.
sum(table(table(ALLannotated_chr[,2])))

#Question 4.
#Use the MsetEx dataset from the minfiData package. Part of this question is to use the help system to figure out how to address the question.
#What is the mean value of the Methylation channel across the features for sample "5723646052_R04C01"?

#Download the methylation data.
data(MsetEx)

#Find out where sample 5723646052_R04C01 is.
sampleNames(MsetEx)

#Get the mean value of the sample.
mean(getMeth(MsetEx)[,2])

#Question 5.
#Access the processed data from NCBI GEO Accession number GSE788.
#What is the mean expression level of sample GSM9024?

#Download the data.
GSE788 <- getGEO("GSE788")
GSE788data <- GSE788[[1]]

#Find out where sample GSM9024 is.
sampleNames(GSE788data)

#Calculate the mean expression level of the sample.
mean(exprs(GSE788data)[,2])

#Question 6.
#We are using the airway dataset from the airway package.
#What is the average of the average length across the samples in the experiment?

#Download the dataset.
data(airway)

#Find out the average length (read length?) in each sample and average these averages.
mean(airway$avgLength)

#Question 7.
#We are using the airway dataset from the airway package. The features in this dataset are Ensembl genes.
#What is the number of Ensembl genes which have a count of 1 read or more in sample SRR1039512?

#Find out where sample SRR1039512 is.
airway$Run

#Count the number of genes with at least one count in the sample.
sum(assay(airway)[,3]>0)

#Question 8.
#The airway dataset contains more than 64k features. How many of these features overlaps with transcripts on the autosomes (chromosomes 1-22) as represented by the TxDb.Hsapiens.UCSC.hg19.knownGene package?
#A feature has to overlap the actual transcript, not the intron of a transcript. So you will need to make sure that the transcript representation does not contain introns.

#Find all the exons in hg19.
exons <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)

#Keep the autosomal exons only.
autosome <- paste0("chr", c(1:22))
exons_autosome <- keepSeqlevels(exons, autosome, pruning.mode = "coarse")

#Find out the NCBI names of the exons and rename the autosomal exons accordingly.
NCBInames <- mapSeqlevels(seqlevels(exons), "NCBI")
exons_autosome_ncbi <- renameSeqlevels(exons_autosome, NCBInames)

#Identify the airway features that overlap with the autosomal exons.
dim(subsetByOverlaps(airway, exons_autosome_ncbi))

#Question 9.
#The expression measures of the airway dataset are the number of reads mapping to each feature. In the previous question we have established that many of these features do not overlap autosomal transcripts from the TxDb.Hsapiens.UCSC.hg19.knownGene. But how many reads map to features which overlaps these transcripts?
#For sample SRR1039508, how big a percentage (expressed as a number between 0 and 1) of the total reads in the airway dataset for that sample, are part of a feature which overlaps an autosomal TxDb.Hsapiens.UCSC.hg19.knownGene transcript?

#Find out where sample SRR1039508 is and extract it.
airway$Run
airway_SRR1039508 <- airway[, 1]

#Identify the features in the sample that overlap with the autosomal exons.
airway_SRR1039508_autosome <- subsetByOverlaps(airway_SRR1039508, exons_autosome_ncbi)

#Count the reads belonging to these features and divide it by the total number of reads in the sample.
sum(assay(airway_SRR1039508_autosome, 'counts'))/sum(assay(airway_SRR1039508, 'counts'))

#Question 10.
#Consider sample SRR1039508 and only consider features which overlap autosomal transcripts from TxDb.Hsapiens.UCSC.hg19.knownGene. We should be able to very roughly divide these transcripts into expressed and non expressed transcript. Expressed transcripts should be marked by H3K4me3 at their promoter. The airway dataset have assayed "airway smooth muscle cells". In the Roadmap Epigenomics data set, the E096 is supposed to be "lung". Obtain the H3K4me3 narrowPeaks from the E096 sample using the AnnotationHub package.
#What is the median number of counts per feature (for sample SRR1039508) containing a H3K4me narrowPeak in their promoter (only features which overlap autosomal transcripts from TxDb.Hsapiens.UCSC.hg19.knownGene are considered)?
#We are using the standard 2.2kb default Bioconductor promotor setting.
#Compare this to the median number of counts for features without a H3K4me3 peak. Note that this short analysis has not taken transcript lengths into account and it compares different genomic regions to each other; this is highly suscepticle to bias such as sequence bias.

#Create an annotation hub.
ah <- AnnotationHub()

#Search for the Epigenomics Roadmap data on H3K4me3 narrowPeaks in the E096 sample.
ah_H3K4me3_E096 <- query(ah, c('E096', 'H3K4me3', 'narrowPeak'))
data_H3K4me3_E096 <- ah_H3K4me3_E096[['AH30596']]

#Keep the autosomal data only and and rename them to NCBI names.
H3K4me3_E096_autosome <- keepSeqlevels(data_H3K4me3_E096, autosome, pruning.mode = "coarse")
H3K4me3_E096_autosome_ncbi <- renameSeqlevels(H3K4me3_E096_autosome, NCBInames)

#Find out the ranges of the features in sample SRR1039508 that overlap with the autosomal exons.
airway_SRR1039508_autosome_ranges <- keepSeqlevels(range(rowRanges(airway_SRR1039508_autosome)), extractSeqlevelsByGroup(species = "Homo sapiens", style = "NCBI", group = "auto"))

#Find out the overlaps between the promoters of these features and the autosomal H3K4me3 narrowPeaks in the E096 sample.
promoters_expressed_ranges <- subsetByOverlaps(promoters(airway_SRR1039508_autosome_ranges), H3K4me3_E096_autosome_ncbi)

#Turn the ranges back to features.
promoters_expressed <- subsetByOverlaps(airway_SRR1039508_autosome, promoters_expressed_ranges)

#Count the number of reads in each feature and average the counts.
median(assay(promoters_expressed, 'counts'))