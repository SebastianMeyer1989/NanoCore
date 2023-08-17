# NanoCore


## Introduction

NanoCore is a user-friendly tool developed specifically for the analysis of putative bacterial outbreaks with Nanopore sequencing on basis of a species specific core genome;
NanoCore is based on a mapping and variant calling approach, followed by multi-level filtering to account for sequencing and variant calling errors and variation in genome structure; It implements a cgMLST-like distance metric and also supports the integration of Illumina- and Nanopore-sequenced isolates in the same analysis. 

Key advantages of this tool are reduced turnaround times and lower capital requirements, due to the use of Oxford Nanopore sequencing data in comparison gold standard methods which are mostly Illumina-based. In the context of hospital pathogen surveillance, these may translate into e.g. reduced outbreak investigation times or the ability to implement genomic surveillance in resource-limited settings.

The results of our experiments so far are highly concordant to the results of SeqSphere+, a “gold standard” software for the analysis of bacterial outbreaks in the hospital context on Illumina-basis used by many hospital hygiene departments worldwide.

To apply NanoCore, Nanopore data (or Nanopore or Illumina data) for either Sample you want to analyse, is needed.

A preprint with accuracy evaluations can be fount at: ???? (WIP)

### Overview of the NanoCore pipeline
![alt text](https://github.com/SebastianMeyer1989/NanoCore/blob/main/NanoCore_workflow.png)

Overview of the NanoCore method (right), in comparison to a well-established method for the computation of cgMLST-based distances (SeqSphere+ (41), left). NanoCore maps the provided sequencing data to the corresponding core genome, calls variants from those mappings, applies various filters to the results and produces different statistics, a distance matrix and a minimum-spanning-tree.

## Program Requirements and Installation
The programm was tested on the following Operating Systems:
- CentOS Linux 6.2.0-26-generic x86_64
- Ubuntu 22.04.3 LTS

### The following programming languages and tools need to be installed:

- Conda
- R
- Perl
- Perl Modules:
  - Getopt::Long
 
Then the NanoCore environment can be created using the following command: (WIP)
```
conda env create -f NanoCore_dependencies.yml
```

And activated usind the following command:
```
conda activate NanoCore
```


## Running NanoCore
The tool has no user-interface and is run from the terminal.
Input sequencing data are specified using a simple sample sheet in tab-separated format; in addition, the user specifies a species-specific core genome reference file. Reference files for 8 bacterial species are included in the NanoCore package (the cgMLST_files folder). In addition, the user may specify a minimum coverage threshold (default 20) and the number of threads used for components of the pipeline that support multithreading.
To run NanoCore the following command is needed:
```
perl NanoCore.pl --sample_list sample_list.txt --reference S.pecies_cgMLST_ref-seqs.fasta --clair_model_nano /home/user/Software/miniconda3/envs/clair3/bin/models/ont --clair_model_illu /home/user/Software/miniconda3/envs/clair3/bin/models/ilmn --threshold 20 --threads 8 --samtools samtools --prefix NanoCore_Run_1
```

### Input explained

- **perl NanoCore.pl** = The NanoCore algorithm.

- **--sample_list sample_list.txt** = A tab-separated file containing one line per sample with the isolate ID, the tag "Nanopore" or "Illummina" to define the used sequencing method, and the paths to either the Nanopore sequencing data file or the Illumina sequencing data R1 and R2 files.
##### Example:
```
Isolate_1    Illumina  /Illumina_Data/isolate_1_R1.fastq  Illumina_Data/isolate_1_R2.fastq
Isolate_17   Nanopore  /Nanopore_Data/isolate_17.fastq
MRSA_H4      Nanopore  /Nanopore_Data/MRSA_H4.fastq
Benjamin     Illumina  /Illumina_Data/Benjamin_R1.fastq  Illumina_Data/Benjamin_R2.fastq
sample404    Illumina  /Illumina_Data/sample404_R1.fastq  Illumina_Data/sample404_R2.fastq
…
```
- **--reference S.pecies_cgMLST_ref-seqs.fasta** = The core genome reference file for a certain species. Files for 8 clinically relevant species are provided in the chMLST_files folder.

- **--clair_model_nano /home/user/Software/miniconda3/envs/clair3/bin/models/ont** = The path do the Nanopore model for the clair3 variant-caller. Should be included in the clair3 installation within the NanoCore package.

- **--clair_model_illu /home/user/Software/miniconda3/envs/clair3/bin/models/ilmn** = The path do the Illumina model for the clair3 variant-caller. Should be included in the clair3 installation within the NanoCore package.

- **--threshold 20** = The minimum coverage threshold desired for the analysis. This value affects some of the implemented filters. If no threshold is set by the user, this valu is per default set to 20.

- **--threads 8** = The number of threads used for components of the pipeline that support multithreading.

- **--samtools samtools** = The samtools executable. Should be included in the samtools installation within the NanoCore package. If no executable is set by the user, this valu is per default set to "samtools".

- **--prefix NanoCore_Run_1** = The chosen prefix/name for the current nanoCore run.

### Output explained

(examplary for prefix `NanoCore_Run_1`):

Of Primary interest for the user should be the allele table "NanoCore_Run_1_allele_table.txt" in the "Output_NanoCore_Run_1_Tables/" folder, that shows the pairwise distances calculated on a cgMLST-like metric and the minimum-spanning-tree "NanoCore_Run_1-mst.pdf" in the "Output_NanoCore_Run_1_Stats/" folder, that shows the tree with all samples of this run calculated from the allele table.
Of secondary interest are probably the other pdf-files found in the "Output_NanoCore_Run_1_Stats/" folder, that show different statistics of the run as well as heatmaps of excluded genes.
Nevertheless, for completeness here we list everything NanoCore produces:

- Output_NanoCore_Run_1_Minimap/[SAMPLE_ID]/
  - [SAMPLE_ID].bam = Original mapping output from the mapper minimap2.
  - [SAMPLE_ID]_new_FLAGs.bam = Correctd mapping output.
  - [SAMPLE_ID]_new_FLAGs.bam.bai = Index file to corrected mapping output.
  - [SAMPLE_ID]_new_FLAGs_mpileup.txt = Samtools mpileup file to corrected mapping output.

- Output_NanoCore_Run_1_Clair/[SAMPLE_ID]/
  - [SAMPLE_ID]_merge_output.vcf = Variant-calls from the variant-caller clair3.
  - run_clair3.log = Log file from the variant-calling.
  - Other files = Secondary clair3 output files.

- Output_NanoCore_Run_1_Stats/
  - [SAMPLE_ID].coverage = Samtools coverage output.
  - [SAMPLE_ID].depth = Samtools depth output.
  - [SAMPLE_ID].stats = Samtools stats output.
  - [SAMPLE_ID].coverage_per_ref = Mean coverage and ratio of positions above coverage threshold per gene.
  - NanoCore_Run_1.read_length = Lengths of all reads in that run.
  - NanoCore_Run_1.min_cov = Matrix of ratio of positions per gene above coverage threshold.
  - NanoCore_Run_1.table = Table with basic read statistics.
  - NanoCore_Run_1.excluded_1_heterozyg = Excluded genes per isolate pair with too high ratio of heterozygous calls.
  - NanoCore_Run_1.excluded_2_cov_ratio = Excluded genes per isolate pair with too low ratio of positions above coverage threshold.
  - NanoCore_Run_1.excluded_3_same_var = Genes per isolate pair with the same variants -> no distance.
  - NanoCore_Run_1.excluded_4_mapq_plus_coverage = Excluded genes per isolate pair with unusual coverage AND low mapQ.
  - NanoCore_Run_1.excluded_5_no_var = Genes per isolate pair without variants.
  - NanoCore_Run_1.excluded_6_mpileup = Genes excluded through mpileup statistics.
  - NanoCore_Run_1.excluded_genes_per_isolate = Summary of exclusion files 1, 2 and 4 ( gene-level filter).
  - NanoCore_Run_1.suspos = Suspicious positions in genes that are excluded from the complete analysis (global positional filter).
  - NanoCore_Run_1.variantlist = List of called variants per isolate pair.
  - NanoCore_Run_1.heat_excluded_ave_cov = Matrix for heatmap of genes with unusual coverage. 
  - NanoCore_Run_1.heat_excluded_heterozyg = Matrix for heatmap of genes with unusual heterozygosity.
  - NanoCore_Run_1.heat_excluded_mapQ = Matrix for heatmap of genes with unusual mapping quality.
  - NanoCore_Run_1-boxplot_min_cov.pdf = Boxplot of ratio of positions per gene above coverage threshold.
  - NanoCore_Run_1-heatmap_min_cov.pdf = Heatmap of ratio of positions per gene above coverage threshold.
  - NanoCore_Run_1-stats_1.pdf = Plots of basic read statistics part 1/2.
  - NanoCore_Run_1-stats_2.pdf = Plots of basic read statistics part 2/2.
  - NanoCore_Run_1-mst.pdf = Minimum-spanning-tree of all samples from this run. Dashed edge means no distance.
  - NanoCore_Run_1-heatmap_cov.pdf = Heatmap of genes with unusual coverage.
  - NanoCore_Run_1-heatmap_hete.pdf = Heatmap of genes with unusual heterozygosity.
  - NanoCore_Run_1-heatmap_mapq.pdf = Heatmap of genes with unusual mapping quality.
	
- Output_NanoCore_Run_1_Tables/
  - NanoCore_Run_1_variant_table.txt = Distance matrix on basis of variants.
  - NanoCore_Run_1_allele_table.txt = Distance matrix on basis of alleles.
  - NanoCore_Run_1_consensus_table.txt = Similarity matrix on basis of variants.
  - NanoCore_Run_1_per_gene_variant_table.txt = Table of variants per gene per isolate pair.
  - NanoCore_Run_1_per_gene_total_table.txt = Total variants per gene over all isolates.
