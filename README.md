# NanoCore



## Introduction

NanoCore is a user-friendly tool developed specifically for the analysis of putative bacterial outbreaks with Nanopore sequencing data on basis of a species specific core genome.
NanoCore is based on a mapping and variant calling approach, followed by multi-level filtering to account for sequencing and variant calling errors and variation in genome structure. It implements a cgMLST-like distance metric and also supports the integration of Illumina- and Nanopore-sequenced isolates in the same analysis. 

Key advantages of this tool are reduced turnaround times and lower capital requirements, due to the use of Oxford Nanopore sequencing data in comparison gold standard methods which are mostly Illumina-based. In the context of hospital pathogen surveillance, these may translate into e.g. reduced outbreak investigation times or the ability to implement genomic surveillance in resource-limited settings.

The results of our experiments so far are highly concordant to the results of SeqSphere+, a “gold standard” software for cgMLST analysis of bacterial outbreaks in the hospital context on Illumina-basis used by many hospital hygiene departments worldwide.

To apply NanoCore, Nanopore data (or Nanopore or Illumina data) for either Sample you want to analyse, is needed.

The publication accompanying this tool can be found at: [WIP]  
The datasets analyzed during the project are available at the Sequence Read Archive under the BioProject ID PRJNA1012291.
The DOI for this software is https://doi.org/10.5281/zenodo.8424707.



### Overview of the NanoCore pipeline
![alt text](https://github.com/SebastianMeyer1989/NanoCore/blob/main/NanoCore_workflow.jpg)

The NanoCore method (right), in comparison to a well-established method for the computation of cgMLST-based distances (SeqSphere+, left). NanoCore maps the provided sequencing data to the corresponding core genome, calls variants from those mappings, applies various filters to the results and produces different statistics, a distance matrix and a minimum-spanning-tree.



## Program Requirements and Installation
The programm was tested on the following Operating Systems:
- CentOS Linux 6.2.0-26-generic x86_64
- Ubuntu 22.04.3 LTS

For NanoCore to run properly, we recommend to use squencing data between 500MB and 3GB file-size. Smaller file could contain a to low amount of sequencing data to run without errors (as e.g. explained further down in "Known issues and how to solve them"). Bigger files may lead to a very long run time and a termination of the script, should the available memory not suffice.



### The following programming languages and tools are necessary for NanoCore:

- Conda
- Perl
- Perl modules:
  - Getopt::Long
  - FindBin
- R
- R Packages:
  - ape (version 5.6.2)
  - pegas (version 1.1)
  - sna (version 2.7)
  - pheatmap (version 1.0.12)

Conda, Perl and the Perl modules need to be provided by the user. R and all needed R packages are installed during the following installing process.

After downloading all files from this repository, NanoCore can be installed using the provided .yml files and the following conda commands, which create the necessary conda environments for NanoCore:
```
conda env create -f NanoCore_1.yml
conda env create -f NanoCore_2.yml
```

The environments will be activated automatically by NanoCore. The user does not need to activate them prior to the run.

Of note: The .yml files do not necessarily contain the newest versions of the used programs, but the versions we used to build the NanoCore tool. The workflow may also work with the most recent versions of the used programs, but it is also possible that some output formats (e.g. of the mapping tool or variant-caller) will be modified, thus interrupting programs later in the pipeline, which expect their input to be in a certain format.

## Running NanoCore
The tool has no user-interface and is run from the terminal.
Input sequencing data are specified using a simple sample sheet in tab-separated format; in addition, the user specifies a species-specific core genome reference file. Reference files for 8 bacterial species are included in the NanoCore package (the "cgMLST_files" folder). In addition, the user needs to specify the number of threads used for components of the pipeline that support multithreading (recommended at least 8) and may specify a minimum coverage threshold (default: 20) and the samtools binary (default: samtools). 

We also prepared two small batches of test data to run a NanoCore example analysis either for VRE in "Nanopore-only" mode or MRSA in "Hybrid" mode. For further instructions on how to download the corresponding datasets and run the analysis, please read the part "NanoCore Example Run" which is located directly after this "Running NanoCore" paragraph.

To get help about the parameters NanoCore needs, simply run one of the following commands:
```
./NanoCore_v1.0.5.sh -h
./NanoCore_v1.0.5.sh --help
```
This will tell you, which input in which format is required for NanoCore to run


To get the currently installed version of NanoCore, run one of the following commands:
```
./NanoCore_v1.0.5.sh -v
./NanoCore_v1.0.5.sh --version
```


To run NanoCore the following command is needed:
```
./NanoCore_v1.0.5.sh -s SampleSheet.txt -r cgMLST_files/Species_cgMLST_ref-seqs.fasta -p NanoCore_Run_1 -t 8 -m 20 -b samtools
```

NanoCore can be executed from a user-chosen folder. Please keep in mind, that the pathways to the NanoCore script, the sample sheet and the reference need to be modified accordingly.

The currently installed clair3 models can be found in your NanoCore_1 environment path (within your conda/miniconda path) e.g. "/path/to/miniconda3/envs/NanoCore_1/bin/models/" and should contain basic models like "ont" (for Nanopore data), "ilmn" (for Illumina data) and a few more. More models for clair3, that maybe fit better to your type of sequencing data, can be downloaded from ONTs GitHub Rerio repo (https://github.com/nanoporetech/rerio) and saved in the same path.



### Adding new samples to a finished NanoCore run
It is possible, to add new samples to a finished analysis. To achieve this, NanoCore needs to be re-run with the same prefix and an extended sample sheet, that contains the already processed samples AND the new samples. Mapping and variant calling (which are the most time-consuming steps) will be skipped for all previously included samples. Only steps that involve all samples (e.g. distance calculation and matrix generation) are re-executed for the complete set of samples.



### Input explained

- **perl NanoCore_v1.0.5.sh** = The NanoCore algorithm.  
- **-s SampleSheet.txt** = A tab-separated file containing one line per sample with the isolate ID, the tag "Nanopore" or "Illummina" to define the used sequencing method, the desired clair 3 model and the paths to either the Nanopore sequencing data file or the Illumina sequencing data R1 and R2 files.  
##### Example:
```
Isolate_1    Illumina  /path/to/NanoCore/bin/models/ilmn  /Illumina_Data/isolate_1_R1.fastq  Illumina_Data/isolate_1_R2.fastq
Isolate_17   Nanopore  /path/to/NanoCore/bin/models/ont  /Nanopore_Data/isolate_17.fastq
MRSA_H4      Nanopore  /path/to/NanoCore/bin/models/ont  /Nanopore_Data/MRSA_H4.fastq
Benjamin     Illumina  /path/to/NanoCore/bin/models/ilmn  /Illumina_Data/Benjamin_R1.fastq  Illumina_Data/Benjamin_R2.fastq
sample404    Illumina  /path/to/NanoCore/bin/models/ilmn  /Illumina_Data/sample404_R1.fastq  Illumina_Data/sample404_R2.fastq
...
...
```
- **-r cgMLST_files/Species_cgMLST_ref-seqs.fasta** = The core genome reference file for a certain species. Files for 8 clinically relevant species are provided in the chMLST_files folder.
- **-p NanoCore_Run_1** = The chosen prefix/name for the current nanoCore run.
- **-t 8** = The number of threads used for components of the pipeline that support multithreading.  
- **-m 20** = The minimum coverage threshold desired for the analysis. This value affects some of the implemented filters. If no threshold is set by the user, this value is per default set to 20.  
- **-b samtools** = The samtools binary. Should be included in the samtools installation within the NanoCore package. If no executable is set by the user, this value is per default set to "samtools".  



### Output explained

(examplary for prefix `NanoCore_Run_1`):

Of Primary interest for the user should be the allele table "NanoCore_Run_1_allele_table.txt", the minimum-spanning-tree "NanoCore_Run_1-mst.pdf" and the log file "NanoCore_Run_1.log", all saved to the folder the script was run in.  
These show the pairwise distances calculated on a cgMLST-like metric, the tree with all samples of this run calculated from the allele table and a list of the complete output of the NanoCore run.  
Of secondary interest are probably the other pdf-files found in the "Output_NanoCore_Run_1_Stats/" folder, that show different statistics of the run as well as heatmaps of excluded genes.
Nevertheless, for completeness here we list everything NanoCore produces:

- Folder, NanoCore was run in:
  - NanoCore_Run_1.log = Log file of the complete run.
  - NanoCore_Run_1_allele_table.txt = Distance matrix on basis of alleles. This file was copied from the "Output_NanoCore_Run_1_Tables/" folder.
  - NanoCore_Run_1-mst.pdf = Minimum-spanning-tree of all samples from this run. Dashed edge means no distance. This file was copied from the "Output_NanoCore_Run_1_Stats/" folder.
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



## NanoCore Example Run

Here we explain how to download the test datasets and run the NanoCore analysis for two small analysis that should be done in few hours.
For both runs NanoCore and the conda enwironments need to be installed according to the instructions above.



### Test Dataset 1: VRE in "Nanopore-only" mode
You can download the corresponding data from the following link: https://osf.io/yz35s/.  
The file you are looking for is called "Testdata_VRE_Nanopore-only.tar.gz" (~2.6gb). It contains Nanopore sequencing data of 5 isolates, a file called "Example_VRE_nanopore_used_ressources.txt", which lists the ressources our test run needed (output of the time command), as well as a sample sheet called "sample_list_VRE_nano_testdata.txt", which lists the ID, the tag "Nanopore", the path to the clair3 model and the paths to the Nanopore sequencing data, as clarified in the paragraph "Input explained" above. You just need extract the tar.gz file into your NanoCore main folder and change the "/path/to/" part for the clair3 model and the downloaded sequencing data to your pathways.


You can then run the analysis using the following command:
```
./NanoCore_v1.0.5.sh -s VRE_Nanopore-only_example/sample_list_VRE_nano_testdata.txt -r cgMLST_files/E.faecium_cgMLST_ref-seqs.fasta -p Example_VRE_nanopore -t 8
```
For information about this and further explanation on what the different options do, please look into the "Running NanoCore" paragraph.  
This analysis should run ~2hours on 8 cores, use <2gb of memory and produce <8gb of output data.  
The distance-matrix "Example_VRE_nanopore_allele_table.txt" and minimum-spanning-tree "Example_VRE_nanopore-mst.pdf" should look like this: 

![alt text](https://github.com/SebastianMeyer1989/NanoCore/blob/main/Testdata_Figures/VRE_Nanocore-only_testdata_distmat.PNG)

![alt text](https://github.com/SebastianMeyer1989/NanoCore/blob/main/Testdata_Figures/VRE_Nanocore-only_testdata_mst.PNG)



### Test Dataset 2: MRSA in "Hybrid" mode

You can download the corresponding data from the following link: https://osf.io/yz35s/.  
The file you are looking for is called "Testdata_MRSA_Hybrid.tar.gz" (<1gb). It contains Nanopore sequencing data of 2 isolates and Illumina sequencing data of 3 isolates, a file called "Example_VRE_nanopore_used_ressources.txt", which lists the ressources our test run needed (output of the time command), as well as a sample sheet called "sample_list_MRSA_hybrid_testdata.txt", which lists the ID, the tag "Nanopore" or "Illumina", the path to the clair3 model and the paths to the Nanopore or Illumina sequencing data, as clarified in the paragraph "Input explained" above. You just need extract the tar.gz file into your NanoCore main folder and change the "/path/to/" part for the clair3 model and the downloaded sequencing data to your pathways.


You can then run the analysis using the following command:
```
./NanoCore_v1.0.5.sh -s MRSA_Hybrid_example/sample_list_MRSA_hybrid_testdata.txt -r cgMLST_files/S.aureus_cgMLST_ref-seqs.fasta -p Example_MRSA_hybrid -t 8
```
For information about this and further explanation on what the different options do, please look into the "Running NanoCore" paragraph.  
This analysis should run ~3hours on 8 cores, use <3gb of memory and produce <5gb of output data.  
The distance-matrix "Example_MRSA_hybrid_allele_table.txt" and minimum-spanning-tree "Example_MRSA_hybrid-mst.pdf" should look like this: 

![alt text](https://github.com/SebastianMeyer1989/NanoCore/blob/main/Testdata_Figures/MRSA_Hybrid_testdata_distmat.PNG)

![alt text](https://github.com/SebastianMeyer1989/NanoCore/blob/main/Testdata_Figures/MRSA_Hybrid_testdata_mst.PNG)



## Known issues and how to solve them
In this paragraph we will explain potential issues you could encounter within your output data.



### Unexpected low-distance pattern in the distance matrix

#### Characteristics
Should you encounter a pattern like this,  
![alt text](https://github.com/SebastianMeyer1989/NanoCore/blob/main/Known_Issues_Figures/low_distance_issue_distmat.PNG)  
where a number of isolates (in this case "SSI-039A" and "SSI-039B") produce a similarly low distance to all other isolates, while the mean distanct between these other isolates is far higher, then there are probably some issues with the used data.  
The corresponding minimum-spanning-tree would look like this:  
![alt text](https://github.com/SebastianMeyer1989/NanoCore/blob/main/Known_Issues_Figures/low_distance_issue_mst.PNG)  
Nearly all isolates are connected to one of the possibly faulty isolates with low distance numbers, that by far do not add up to the >1500 distance that should exist between the other isolates.  

#### Reasons
There are two probable reasons for this:  
First, it is possible, that there is not enough sequencing data for these isolates. A too low overall coverage (minimun threshold is 20x per default) could exclude too many genes through the included coverage filters, thus leaving too few genes for a reasonable analysis.  
Second, it is possible, that your sequencing data contains contaminations (some isolates got mixed during the library preparation, etc.). This could lead to the exclusion of too many genes through the included heterozygosity filter, also leaving too few genes for a reasonable analysis.  

If this (or a similar issue) is the reason, you should be able to see a few patterns in some of the other output files:  
The per default created heatmaps for genes with such issues should show a higher numer of genes for these isolates. The corresponding heatmaps can be found in the "Output_[PREFIX]_Stats" directory:  
  - [PREFIX]-heatmap_cov.pdf = Heatmap of genes with unusual coverage.
  - [PREFIX]-heatmap_hete.pdf = Heatmap of genes with unusual heterozygosity.
  - [PREFIX]-heatmap_mapq.pdf = Heatmap of genes with unusual mapping quality.

Part of the heterozygosity heatmap looks for example like this:  
![alt text](https://github.com/SebastianMeyer1989/NanoCore/blob/main/Known_Issues_Figures/low_distance_issue_heatmap.PNG)  
"SSI-039A" and "SSI-039B" are clearly cases with an above-average number of such genes.  
Of note: Isolates like "SSI-128B" could also be such candidates, but this isolate was kept in the example here, since it did not show a suspicious behaviour.  

Another indicator is the samtools coverage file "[SAMPLE_ID].coverage" found in the "Output_[PREFIX]_Stats" directory. While the sequencing data should normaly include data for almost all analysed genes of the corrresponding core genome, here more that 1000 genes are missing (while the present genes all have a single-digit coverage):  
![alt text](https://github.com/SebastianMeyer1989/NanoCore/blob/main/Known_Issues_Figures/low_distance_issue_cov.PNG)  
This is only a screenshot from the middle of the file. First, second and third column are the Gene-ID, and the start and end of the gene. All the zeros after that show that no data was found.  

#### Solution
Since this issue results from a low amount of sequencing data or a contamination in your sample, best practice would be to just re-sequence all problematic isolates, ideally also repeating the isolate cultivation.  
If this is not feasible, you could just exclude the isolates in question from your sample sheet and rerun the analysis. This should run much faster now, since the most-time-consuming parts of the pipeline (the mapping and variant-calling) are already done and do not need to be repeated.  



### Unexpected NO-distance pattern in the distance matrix  

#### Characteristics
This is basically a far more drastic version of the former "low-distance pattern in the distance matrix" issue.  
Should you encounter a pattern like this,  
![alt text](https://github.com/SebastianMeyer1989/NanoCore/blob/main/Known_Issues_Figures/no_distance_issue_distmat.PNG)  
where a number of isolates (in this case "GW20g") produce no diatance at all to (all) other isolates, then there are probably some issues with the used data.  
The corresponding minimum-spanning-tree could look like this:  
![alt text](https://github.com/SebastianMeyer1989/NanoCore/blob/main/Known_Issues_Figures/no_distance_issue_mst.PNG)  
Nearly all isolates are connected to the possibly faulty isolate with no distance.  

#### Reasons
Same as in the former issue, but more drastic. Basically no sequencing data at all, or a contamination over the complete genome length
The above mentioned heatmaps and coverage files should show more conspicious genes, and less present genes respectively.  

#### Solution
Also same as in the former issue: Re-cultivating plus re-sequencing, or exclusion of the isolate from the analysis.  
