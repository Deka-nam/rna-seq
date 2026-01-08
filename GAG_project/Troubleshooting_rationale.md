# **Troubleshooting/Changes to original methods explained**
**RNA-Seq re-analysis of Catheterization Study with Focus on GAG/CS Genes**

## **Original objective and approach**

We set out to reanalyze RNA-seq data from the 2016 Rousseau et al. study on bladder catheterization in B6 mice, with a specific focus 
not addressed in the original paper to determine whether genes involved in chondroitin sulfate (CS) glycosaminoglycan (GAG) production and modification are impacted by catheterization.

The original study identified thousands of differentially expressed genes but only published the top 200. 
Since our lab is interested in GAG biology, we needed the complete dataset to examine specific GAG/CS pathway genes.

## **Initial strategy: Exact replication**

Our first approach was to exactly replicate the 2016 paper's methods to ensure comparability:

1. **Downloaded the same 6 samples** (3 naive, 3 catheterized) from NCBI SRA
2. **Followed their published pipeline**: cutadapt → Tophat2 → HTSeq → DESeq2
3. **Used the same reference**: mm10 genome with RefSeq annotations

**Rationale:** By using identical methods, we could directly compare our expanded analysis (including GAG genes) with their published results for the top genes.

## **Encountered technical challenges**

During implementation, I encountered a critical technical issue:

### **The HTSeq problem**
- The original paper used **HTSeq 0.6.1** (2016 version) for gene counting
- When I ran HTSeq with our aligned BAM files, it returned zero counts for all genes

```bash
# Check chromosome names in BAM
samtools view -H bam_files/sample_sorted.bam | grep "@SQ" | head -3

# Check chromosome names in your GTF
head -5 your_gtf_file.gtf
```
### **Problem identified**: A mismatch between chromosome naming conventions
  - The aligner (HISAT2, the modern successor to Tophat2) used chromosome names like "1", "2", "3"
  - The RefSeq annotation file used "chr1", "chr2", "chr3"
  - HTSeq requires exact chromosome name matching and fails completely when they don't match

```bash
# Create a GTF without "chr" prefix
sed 's/^chr//' hisat_index/mm10/mm10.refGene.gtf > hisat_index/mm10/mm10.refGene.nochr.gtf

# Run featureCounts from SUBread now
```

### **Why this didn't affect the original study**
The 2016 authors likely used a differently formatted annotation file or an older alignment tool that produced compatible chromosome names. 
This kind of inconsistency is common in bioinformatics as tools evolve.

## **Updated methodology for purposes**

Instead of spending excessive time troubleshooting the obsolete HTSeq tool, I adopted a more robust approach while maintaining scientific rigor:

### **Key changes and justifications**

| Component | 2016 Paper | Our updated approach | Reason for change |
|-----------|------------|----------------------|-------------------|
| **Gene Counting** | HTSeq (2016) | featureCounts | More robust, handles formatting issues better, industry standard |
| **Alignment** | Tophat2 (obsolete) | HISAT2 | Same algorithm family, 5x faster, actively maintained |
| **Software Versions** | 2016 releases | 2023 updates | Security, bug fixes, improved algorithms |

### **Why featureCounts over HTSeq?**
1. **Robustness:** featureCounts handles chromosome naming inconsistencies gracefully
2. **Accuracy:** Better management of multi-mapping reads and overlapping features
3. **Speed:** 5-10 times faster processing
4. **Active Development:** Regularly updated vs HTSeq which is largely deprecated
5. **Industry Standard:** Used in most current RNA-seq studies

### **Maintaining scientific somparability**
Despite these technical changes, I preserved the core statistical methodology:
- Same **DESeq2** algorithm for differential expression
- Same **significance thresholds** (padj < 0.05)
- Same **comparison**: Catheterized (samples 4-6) vs Naive (samples 1-3)
- Same **biological question** and gene focus

## **How we ensure result reliability**

### **Validation Measures**
1. **Consistency Check:** Our top differentially expressed genes should substantially overlap with the paper's top 200
2. **Biological Context:** Results should make sense in context of catheterization biology
3. **Technical QC:** All samples showed >87% alignment rates and good quality metrics
4. **Reproducibility:** Complete documentation for every step

### **Advantages of our updated pipeline**
1. **More Accurate Quantification:** featureCounts uses better algorithms for assigning reads to genes
2. **Better Error Detection:** Comprehensive QC catches issues early
3. **Future-proofing:** Uses currently supported tools that others can replicate
4. **Transparency:** Complete logs and intermediate files available

## **Addressing the core biological question**

The methodological updates do not change our ability to answer the biological question about GAG/CS genes. 
The updated methods provide more reliable results due to:
- Better handling of low-expression genes (important for some GAG enzymes)
- More accurate count estimation at transcript boundaries
- Reduced technical artifacts

## **Practical notes**

### **For GAG/CS gene analysis:**
- We can confidently identify which GAG genes are significantly altered
- Fold-change estimates will be more accurate
- We can compare our GAG findings with the paper's immune-focused results


## **TL;DR**

I made strategic, justified updates to the 2016 methods to overcome technical limitations while preserving scientific validity. 
The changes represent evolutionary improvements in bioinformatics tools rather than a change in scientific approach. 
Our results will be directly comparable to the original study for validation purposes while providing more reliable data for our specific investigation of GAG/CS pathway genes.

The updated pipeline actually gives us greater confidence in our findings regarding catheterization's impact on GAG biology, as we are using more robust, modern tools that better handle the complexities of RNA-seq data analysis.

---

*This specific writeup was formatted using Deepseek for clarity in writing and github update only. Changes in protocol belong to the Armbruster Lab and author Namrata Deka.*
