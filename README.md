# MeRIPseq

MeRIPseq (methylated RNA-immunoprecipitation) is a technique to identify transcriptome-wide m6A deposition on RNA. To identify m6A, the computational analysis of sequencing data differs from regular RNAseq. We have to call enriched sequences in immunoprecipitated samples.

## Quality Control
Trimmed demultiplexed reads were quantified using `fastQC`.
```
fastqc *.gz \
--outdir=/fastqc
```
## Read Mapping
High quality reads were mapped to the genome using `STAR ver. 2.7.5` with the following options. The command was run as a for loop. We used the Arabidopsis genome version TAIR10.
```
STAR --runThreadN 16 \
--genomeDir ~/genome/
--sjdbGTFfile ~/genome/genes.gtf \
--readFilesIn ~/data/file.fastq.gz \
--outFileNamePrefix ~/analysis/star/aligned/file \
--readFilesCommand zcat \
--outSAMstrandField intronMotif \
--outFilterIntronMotifs RemoveNoncanonical \
--outSAMtype BAM SortedByCoordinate
```
`STAR` appends a long suffix to the bam files (Aligned.sortedByCoord.out). These can be removed using
```
cd aligned
for i in *_R1_001.fastq.gzAligned.sortedByCoord.out.bam; do mv -i "$i" "${i%_R1_001.fastq.gzAligned.sortedByCoord.out.bam}".bam; 
done
```
## Identification of methylated loci (peak calling)
Input and IP samples were sequenced to be able to identify RNA fragments that were enriched by immunoprecipitation with an a-m6A antibody. The peak caller `exomePeak ver. 1.6.0` was used in this study. It runs in R and takes .bam files as input to identify enriched regions.
```
library(exomePeak)

setwd("~/data/bam/")
gtf="~/genome/genes.gtf"

IP_BAM=c("ip_1.bam","ip_2.bam")
INPUT_BAM=c("input_1.bam","input_2.bam")
result=exomepeak(GENE_ANNO_GTF=gtf,
                 IP_BAM=IP_BAM,
                 INPUT_BAM=INPUT_BAM,
                 OUTPUT_DIR="~/data/exomePeak",
                 EXPERIMENT_NAME="sample_name")

quit("no")
```

## Motif analysis
Enriched RNA motifs were identified using `HOMER ver. 4.1`. Output files from `exomePeak` need to be reformatted to be suitable as inputs for `HOMER`. The following code was used for this task.
```
sample_name <- "sample_name"
sample_path <- paste("../exomepeak/",sample_name,"/peak.bed", sep = "")
output_name <- paste(sample_name,".txt", sep = "")
df <- as.data.frame(read.table(sample_path,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
df$id <- row.names(df)
df <- df[,c(13,1,2,3,6)]
write.table(df, output_name.txt, sep="\t", row.names=FALSE, quote = FALSE)
```
`HOMER` was executed using the following code.
```
findMotifsGenome.pl output_name.txt tair10 sample_name -rna -len 6 -p 20
```

## Differential gene expression analysis
MeRIPseq files were not used to determine differentially expressed genes, because the read depth was too shallow. Instead, we performed RNAseq on rRNA depleted samples, mapped them to the genome as described above and performed differential gene expression analysis using `cuffdiff ver. 2.1.1`. 
```
cuffdiff \
-o ~/data//cuffdiff/sample1_sample2 \
-L sample1,sample2 \
--compatible-hits-norm \
-p 30 \
-u ~/genome/genes.gtf \
sample1_1.bam,sample1_2.bam sample2_1.bam,sample2_2.bam
```
## 
