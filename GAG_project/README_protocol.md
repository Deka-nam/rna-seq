# **RNA-Seq Analysis Protocol: Differential Gene Expression in Catheterized vs Naive B6 Mouse Bladders**

## **Experimental Question**
**Does bladder catheterization alter the expression of genes involved in chondroitin sulfate (CS) glycosaminoglycan (GAG) production and modification in C57BL/6 mice?**

*Background:* This study reanalyzes RNA-seq data from Rousseau et al. (2016) [JCI Insight 1(15):e88178](https://insight.jci.org/articles/view/88178#SEC4) to investigate the impact of urinary catheterization on GAG/CS pathway genes, which were not examined in the original publication.

## **Data Source**
- **Original Study:** Kline et al. (2016) "Bladder catheterization increases susceptibility to infection that can be prevented by prophylactic antibiotic treatment"
- **BioProject:** PRJNA335539 (6 samples total)
- **SRA Accessions:**  
1. naive 1- SRR5478735
2. naive 2- SRR5478734
3. naive 3- SRR5478733
4. catheterized 1- SRR5478732
5. catheterized 2- SRR5478731
6. catheterized 3- SRR5478730
   
- **Experimental Design:** 24-hour post-catheterization vs naive mouse bladders

## **Computational Environment**

### **Software Versions**
| Software | Version Used | Purpose |
|----------|--------------|---------|
| **SRA Toolkit** | 2.11.0 | Download sequencing data |
| **Trimmomatic** | 0.39 | Read quality trimming |
| **HISAT2** | 2.2.1 | Read alignment |
| **SAMtools** | 1.15 | BAM file processing |
| **Subread** | 2.0.3 | featureCounts for quantification |
| **R** | 4.3.0 | Statistical analysis |
| **DESeq2** | 1.40.2 | Differential expression analysis |
| **Bioconductor** | 3.17 | Annotation packages |

### **Reference Genome**
- **Genome Assembly:** mm10 (GRCm38)
- **Annotation:** UCSC RefSeq (`mm10.refGene.gtf`)
- **Chromosome Convention:** Ensembl-style (1, 2, 3... without "chr" prefix)

## **Analysis Pipeline: Original vs Modified**

### **Comparison with Original Publication Methods**

| Analysis Step | Kline et al. (2016) | This Study | Rationale for Change |
|---------------|----------------------|------------|---------------------|
| **Adapter Removal** | cutadapt-1.4.1 | Trimmomatic v0.39 | Comparable performance, better integration |
| **Alignment** | Tophat-2.0.11 | **HISAT2 v2.2.1** | Successor to Tophat, faster & more accurate |
| **Quantification** | HTSeq-0.6.1 | **featureCounts v2.0.3** | More robust, better handling of multimappers |
| **DE Analysis** | DESeq2 v1.10.1 | **DESeq2 v1.40.2** | Same method, updated version |
| **Gene IDs** | RefSeq → Entrez | **RefSeq → Gene Symbols** | More intuitive for downstream analysis |

## **Detailed Protocol**

### **Step 1: Data Download and Organization**

**SRA Download Commands**
```bash
# Download SRA files
prefetch SRR5453401 SRR5453402 SRR5453403 SRR5453404 SRR5453405 SRR5453406

# Convert to FASTQ
fastq-dump --split-files --gzip SRR5453401
# Repeat for all samples

# Rename for clarity
mv SRR5453401_1.fastq.gz naive_1_R1.fastq.gz
mv SRR5453401_2.fastq.gz naive_1_R2.fastq.gz
# ... repeated for all samples
```

### **Step 2: Quality Control and Trimming**

#### **Trimmomatic Parameters**
```bash
java -jar trimmomatic-0.39.jar PE \
  -threads 8 \
  -phred33 \
  naive_1_R1.fastq.gz naive_1_R2.fastq.gz \
  naive_1_R1_trimmed_paired.fastq.gz naive_1_R1_trimmed_unpaired.fastq.gz \
  naive_1_R2_trimmed_paired.fastq.gz naive_1_R2_trimmed_unpaired.fastq.gz \
  ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
  LEADING:20 \
  TRAILING:20 \
  SLIDINGWINDOW:4:15 \
  MINLEN:36
```
### **Step 3: Read Alignment with HISAT2**

#### **Reference Index Preparation**
```bash
# Download mm10 genome (without chr prefix)
wget ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

# Download RefSeq annotation
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz

# Remove 'chr' prefix from GTF for HISAT2 compatibility
sed 's/^chr//' mm10.refGene.gtf > mm10.refGene.nochr.gtf

# Build HISAT2 index
hisat2-build -p 8 Mus_musculus.GRCm38.dna.primary_assembly.fa hisat_index/mm10_index
```

#### **Alignment Commands**
```bash
hisat2 -x hisat_index/mm10_index \
  -1 naive_1_R1_trimmed_paired.fastq.gz \
  -2 naive_1_R2_trimmed_paired.fastq.gz \
  --rna-strandness RF \
  --dta \
  -p 8 \
  -S naive_1_aligned.sam
```

### **Step 4: BAM File Processing**
```bash
# Convert SAM to BAM, sort, and index
samtools view -bS naive_1_aligned.sam | \
  samtools sort -o naive_1_sorted.bam -
samtools index naive_1_sorted.bam

# Move to organized directory
mv naive_1_sorted.bam bam_files/
mv naive_1_sorted.bam.bai bam_files/
```

### **Step 5: Gene Quantification with featureCounts**

#### **featureCounts Command**
```bash
featureCounts -T 8 \
  -a refs/mm10.refGene.nochr.gtf \
  -o counts/all_samples_counts.txt \
  -t exon \
  -g gene_id \
  -s 0 \
  bam_files/naive_1_sorted.bam \
  bam_files/naive_2_sorted.bam \
  bam_files/naive_3_sorted.bam \
  bam_files/catheterized_1_sorted.bam \
  bam_files/catheterized_2_sorted.bam \
  bam_files/catheterized_3_sorted.bam
```

### **Step 6: Differential Expression Analysis with DESeq2**

#### **R Analysis Script (`run_deseq.R`)**
The complete script performs:
1. **Data Import:** Load featureCounts output
2. **Quality Control:** Filter low-count genes
3. **Normalization:** DESeq2's median-of-ratios
4. **Statistical Testing:** Wald test for Catheterized vs Naive
5. **Annotation:** Convert RefSeq IDs to gene symbols
6. **GAG/CS Gene Extraction:** Focus on pathway genes
7. **Visualization:** Volcano plots, MA plots, heatmaps

#### **Key R Functions Used**
```r
# DESeq2 workflow
dds <- DESeqDataSetFromMatrix(countData, colData, ~ condition)
dds <- DESeq(dds)
results <- results(dds, contrast = c("condition", "Catheterized", "Naive"))
```

### **Step 7: GAG/CS Gene Analysis**

#### **Target Gene List**
The analysis specifically examined these GAG/CS pathway genes:

```r
gag_cs_genes <- c(
  # Core synthesis
  "Chsy1", "Chsy3", "Chpf", "Chpf2", "Csgalnact1", "Csgalnact2",
  # Sulfotransferases  
  "Chst3", "Chst7", "Chst11", "Chst12", "Chst13", "Chst14", "Chst15",
  "Galnac4s6st", "Ust",
  # Epimerases
  "Dse", "Ds-el",
  # Additional GAG-related
  "Xylt1", "Xylt2", "B4galt7", "B3galt6", "Ext1", "Ext2"
)
```

## **Quality Control Metrics**

### **MultiQC Report**
```bash
# Generate comprehensive QC report
multiqc . -o results/qc_report/
```

Report includes:
- **Before/After Trimming:** Read quality scores, adapter content
- **Alignment:** Mapping rates, insert sizes, strand specificity
- **Gene Counts:** Distribution, saturation curves
- **Sample Similarity:** PCA plot, correlation heatmap

### **Sample Similarity Assessment**
```r
# Check sample clustering
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")
```

Expected pattern: Naive samples cluster together, catheterized samples cluster together, with clear separation between groups.

## **Statistical Parameters**

### **DESeq2 Settings**
- **Filtering:** Genes with < 10 total counts removed
- **Normalization:** Median-of-ratios method
- **Dispersion estimation:** Parametric fit
- **Testing:** Wald test with Benjamini-Hochberg correction
- **Significance threshold:** padj < 0.05

### **Multiple Testing Correction**
```r
# Benjamini-Hochberg false discovery rate control
results(dds, alpha = 0.05, lfcThreshold = 0)
```

## **Reproducibility Information**

### **Session Info**
```r
# Captured with:
sessionInfo()
# Saved to: results/session_info.txt
```

### **Compute Resources**
- **Platform:** Linux cluster
- **Memory:** 32 GB RAM per job
- **CPUs:** 8 cores for alignment and counting
- **Storage:** ~50 GB temporary, ~10 GB final results


## **References**

1. Kline et al. (2016). *Bladder catheterization increases susceptibility to infection that can be prevented by prophylactic antibiotic treatment*. JCI Insight 1(15):e88178.
2. Kim D, et al. (2019). *HISAT2: graph-based alignment of next-generation sequencing reads to a population of genomes*. Genome Biology 20:16.
3. Liao Y, et al. (2014). *featureCounts: an efficient general purpose program for assigning sequence reads to genomic features*. Bioinformatics 30(7):923-30.
4. Love MI, et al. (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2*. Genome Biology 15:550.

---
**TL;DR**
- *This protocol documents the re-analysis of catheterization RNA-seq data with specific focus on GAG/CS pathway genes. All steps were performed to ensure reproducibility and transparency.*
---

*This specific file was formatted using Deepseek for clarity in writing and github uodate only. Biological question, protocol decisions and analysis was done by the Armbruster Lab - Chelsie E. Armbruster, Ph.D. & Vits Brix, MS is to be credited for the experimental question and Namrata Deka, MS for analysis.*
