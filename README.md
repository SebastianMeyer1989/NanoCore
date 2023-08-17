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

- Output_NanoCore_Run_1_Minimap/[SAMPLE_ID]/
  - [SAMPLE_ID].bam = The original mapping output.
  - [SAMPLE_ID]_new_FLAGs.bam = The correctd mapping output.
  - [SAMPLE_ID]_new_FLAGs.bam.bai = The index file to corrected mapping output.
  - [SAMPLE_ID]_new_FLAGs_mpileup.txt = The samtools mpileup file to corrected mapping output.

- Output_NanoCore_Run_1_Clair/[SAMPLE_ID]/
  - [SAMPLE_ID]_merge_output.vcf\t\t(Variant calling output)
  - run_clair3.log\t\t(Variant calling log file)
  - Other files\t\t(Secondary clair3 output files)

- Output_NanoCore_Run_1_Stats/
  - [SAMPLE_ID].coverage\t\t(Samtools coverage output)
  - [SAMPLE_ID].depth\t\t(Samtools depth output)
  - [SAMPLE_ID].stats\t\t(Samtools stats output)
  - [SAMPLE_ID].coverage_per_ref\t\t(Mean coverage and ratio of positions above coverage threshold per gene)
  - NanoCore_Run_1.read_length\t\t(Lengths of all reads in that run)
  - NanoCore_Run_1.min_cov\t\t(Matrix of ratio of positions per gene above coverage threshold of ".$threshold.")
  - NanoCore_Run_1.table\t\t(Table with basic read statistics)
  - NanoCore_Run_1.excluded_1_heterozyg\t\t(Excluded genes per isolate pair with too high ratio of heterozygous calls)
  - NanoCore_Run_1.excluded_2_cov_ratio\t\t(Excluded genes per isolate pair with too low ratio of positions above coverage threshold)
  - NanoCore_Run_1.excluded_3_same_var\t\t(Genes per isolate pair with the same variants -> no distance)
  - NanoCore_Run_1.excluded_4_mapq_plus_coverage\t\t(Excluded genes per isolate pair with unusual coverage AND low mapQ)
  - NanoCore_Run_1.excluded_5_no_var\t\t(Genes per isolate pair without variants)
  - NanoCore_Run_1.excluded_6_mpileup
  - NanoCore_Run_1.excluded_genes_per_isolate
  - NanoCore_Run_1.suspos\t\t(Suspicious positions in genes that are excluded from the complete analysis due to heterozygozity, low quality or low allele frequenzy)
  - NanoCore_Run_1.variantlist\t\t(List of called variants per isolate pair)
  - NanoCore_Run_1.heat_excluded_ave_cov\t\t(Matrix for heatmap of genes with unusual coverage) 
  - NanoCore_Run_1.heat_excluded_heterozyg\t\t(Matrix for heatmap of genes with unusual heterozygozity)
  - NanoCore_Run_1.heat_excluded_mapQ\t\t(Matrix for heatmap of genes with unusual mapping quality)
  - NanoCore_Run_1-boxplot_min_cov.pdf\t\t\t\t(Boxplot of ratio of positions per gene above coverage threshold of ".$threshold.")
  - NanoCore_Run_1-heatmap_min_cov.pdf\t\t(Heatmap of ratio of positions per gene above coverage threshold of ".$threshold.")
  - NanoCore_Run_1-stats_1.pdf\t\t(Plots of basic read statistics part 1/2)
  - NanoCore_Run_1-stats_2.pdf\t\t(Plots of basic read statistics part 2/2)
  - NanoCore_Run_1-mst.pdf\t\t(Minimum-spanning-tree of all samples from this run. Dashed edge means no distance)
  - NanoCore_Run_1-heatmap_cov.pdf\t\t(Heatmap of genes with unusual coverage)
  - NanoCore_Run_1-heatmap_hete.pdf\t\t(Heatmap of genes with unusual heterozygozity)
  - NanoCore_Run_1-heatmap_mapq.pdf\t\t(Heatmap of genes with unusual mapping quality)
	
- Output_NanoCore_Run_1_Tables/
  - NanoCore_Run_1_variant_table.txt\t\t(Distance matrix on basis of variants)
  - NanoCore_Run_1_allele_table.txt\t\t(Distance matrix on basis of alleles)
  - NanoCore_Run_1_consensus_table.txt\t\t(Similarity matrix on basis of variants)
  - NanoCore_Run_1_per_gene_variant_table.txt\t\t(Table of variants per gene per isolate pair)
  - NanoCore_Run_1_per_gene_total_table.txt\t\t(Total variants per gene over all isolates)




## Creating fastq-files for further hybrid assemblies

(exemplary for prefix `mixed_bacteria_10x`)
 
### Call:
```
perl create_kmer_based_fastq_for_real_data.pl mixed_bacteria_10x.classification_k19.called_kmers path/to/longreads/longreads1.fastq  mixed_bacteria_10x
```

### Input:

**perl create_kmer_based_fastq_for_real_data.pl** = The script that produces fastq files from the calling table.

**mixed_bacteria_10x.classification_k19.called_kmers** = The calling table from the UltraPlexer run.

**path/to/longreads/longreads1.fastq** = The used long-read file.

**mixed_bacteria_10x** = Prefix for the run.

### Output:

**mixed_bacteria_10x-isolate1-predicted_reads.fastq** = A fastq file named after the run (mixed_bacteria_10x) and the isolate ID (isolate1), ending with “predicted_reads.fastq”.

## Example Run

In the following we exemplary describe how to simulate reads, create a read-pool, run the Ultraplexing algorithm and assemble the assigned reads, on basis of three random plasmids.
The nessecary data can be found in the folder “Example1”: 
- Three fasta files (Plasmid1.fna, Plasmid2.fna and Plasmid3.fna)
- A list of these three fasta files (example1_list_of_plasmids.txt)
- Different scripts (.pl and .R), needed to run the whole pipeline

After downloading the Folder "Example1" you just need to switch to it via terminal and call the necessary commands in the further  explained order (assumed, that all the requirements are met).


If you don't want to simulate data or the simulation is just not possible on your device, you can skip Step 1 and 2 (simulation and creation of the long-read pool) and use the long-read pool and other needed data we provided.
Therefor you need to download the zip-file "`Example1_simulated_example_data.zip`" from

https://uni-duesseldorf.sciebo.de/s/oHFl3FCArhPhHb5

and copy the content of the unzipped Folder "`Example1_simulated_example_data`" (Sim_Pipeline, example1_plasmid_ids_and_pathways.txt, example1_plasmid_read_pool.fastq and example1_plasmid_stats.txt) into your Folder "`Example1`". Before continuing the pipeline from step 3. "Ultraplexing:..." you need to open "`example1_plasmid_ids_and_pathways.txt`" and change the absolute pathways of the simulated short-read files (contained in the folder "`Sim_Pipeline`"), to the the pathways fitting your storage location.

### 1. Simulation (59s runtime, 1 CPUs, <1gb used memory)

#### Requirements:
- perl
- pbsim
- wgsim

#### Call:
```
perl simulation_pipeline.pl example1_list_of_plasmids.txt 8500 150
```
Important: Pathways for pbsim, pbsim qc-model and wgsim need to be replaced in the script "simulation_pipeline.pl" (line 10-12) to fit your installations, before running it.

#### Input:
- Genome files to simulate reads from
- List with the filenames of these genomes (example1_list_of_plasmids.txt)
- Desired mean read length of the long reads (8500)
- Desired coverage of the long reads (150)

Important: Please keep these numerical parameters for the example run as they are, since the scripts are still hard-coded at the moment for this mean lengh and coverage.

#### Output:
- A new folder (Sim_Pipeline) containing
  - Folder with data for each simulated genome (Plasmid1_l8500_c150 for example)
  - File with genomes files, the program could not find (missing_files.log)
  - File with the used simulation parameters (simulation_parameters.log)
- The final simulated reads are found in the files containing the tag “filtered”. For the genome “Plasmid1” they would be called:
  - “Plasmid1_l8500_c150-filtered_R1.fastq” (Illumina R1 short-reads)
  - “Plasmid1_l8500_c150-filtered_R1.fastq” (Illumina R2 short-reads)
  - “Plasmid1_l8500_c150-filtered_complete.fastq” (Nanopore long-reads)

### 2. Creating a shuffled long-read pool (1s runtime, 1 CPUs, <1gb used memory)

#### Requirements:
- perl

#### Call:
```
perl create_pool.pl example1_list_of_plasmids.txt 3 10000000 example1_plasmid
```
#### Input:
- List with the filenames of these genomes (example1_list_of_plasmids.txt)
- Number of genomes that should be pooled together from that list (3)
- Number of bases that should be pooled together in total (10000000)
- Prefix for the run (example1_plasmid)

#### Output:
- Shuffled long-read pool (example1_plasmid_read_pool.fastq)
- File containing parameters of the pooling step (example1_plasmid_stats.txt)
- File containing the IDs and simulated short-read pathways of all pooled genomes, needed for the next step (example1_plasmid_ids_and_pathways.txt)

### 3. Ultraplexing: Classification of long-reads (2m32s runtime, 1 CPUs, <7gb used memory)

#### Requirements:
- perl
- R

#### Call:
```
perl UltraPlexer.pl --prefix example1_plasmid --action classify --samples_file example1_plasmid_ids_and_pathways.txt --longReads_FASTQ example1_plasmid_read_pool.fastq
```
#### Input:
- Prefix for the run (example1_plasmid)
- Action the script should do (classify)
- Tab seperated file containing the IDs and absolute pathways of the simulated short-read of all pooled genomes, one genome per line (example1_plasmid_ids_and_pathways.txt). This was produced by the pooling script.
- Shuffled long-read pool (example1_plasmid_read_pool.fastq). This was produced by the pooling script.

#### Output:
- File with classified long-reads (example1_plasmid.classification_k19)
- File produced when the classification step finished correctly (example1_plasmid.classification_k19.done)
- Folder containing temporary files produced by cortex (cortex_temp)

### 4. Ultraplexing: Creating the assignment table (1s runtime, 1 CPUs, <1gb used memory)

#### Requirements:
- perl
- R

#### Call:
```
perl UltraPlexer.pl --prefix example1_plasmid --action generateCallFile --samples_file example1_plasmid_ids_and_pathways.txt
```
#### Input:
- Prefix for the run (example1_plasmid)
- Action the script should do (generateCallFile)
- Tab seperated file containing the IDs and absolute pathways of the simulated short-read of all pooled genomes, one genome per line (example1_plasmid_ids_and_pathways.txt). This was produced by the pooling script.

#### Output:
- Tab separated assignment table (example1_plasmid.classification_k19.called_kmers)


### 5. Creating fastq-files for each genome (1s runtime, 1 CPUs, <1gb used memory)

#### Requirements:
- perl

#### Call:
```
perl create_kmer_based_fastq_for_simulations.pl example1_plasmid.classification_k19.called_kmers example1_plasmid_read_pool.fastq example1_plasmid
```
#### Input:
- Tab separated assignment table (example1_plasmid.classification_k19.called_kmers)
- Shuffled long-read pool (example1_plasmid_read_pool.fastq)
- Prefix for the run (example1_plasmid)

#### Output:
- Two fastq-files for each simulated genome:
  - Predicted reads for the genome in a file, ending with “...predicted_reads.fastq”.
  - True reads for the genome in a file, ending with “...true_reads.fastq”. These are the reads, that were originally simulated for the genome.

### 6. Comparing true and predicted reads (1s runtime, 1 CPUs, <1gb used memory)

(examplary for `Plasmid1`)

#### Requirements:
- perl

#### Call:
```
perl parse_calling_tbl.pl  example1_plasmid.classification_k19.called_kmers
```
#### Input:
- Tab separated assignment table (example1_plasmid.classification_k19.called_kmers)

#### Output:
- File containing primary stats (example1_plasmid.classification_k19.called_kmers.stats). These are not necessarily important for you.
- File containing final stats (example1_plasmid.classification_k19.called_kmers.stats2). This file shows the name of the used calling table and the number and ratio of correct called reads in total (first line) and the number of falsly and correctly called reads for every used genome (following lines).
##### Example:
```
example1_plasmid.classification_k19.called_kmers SummaryCorrect_Reads:	1147	False_Reads:	47	Ratio_Correct_Reads:	0.960636515912898
Plasmid1	false: 4Plasmid1	true: 394
Plasmid2	false: 16Plasmid2	true: 382
Plasmid3	false: 27Plasmid3	true: 371
```
Here you can see, that over 96% of the simulated reads were assigned correctly. The missing <4% are miss-assignments due to sequence homology, or do not really affect hybrid assemblies (as we found out in our experiments).

### 7. Hybrid assembly (1-2h runtime per assembly, 2 CPUs, <2gb used memory)

(examplary for predicted reads for `Plasmid1`)

#### Requirements:
- perl
- Python
- Unicycler
- SPAdes
- Racon
- Pilon
- Bowtie2
- SamTools
- BLASTp

#### Call:
```
/gpfs/project/dilthey/software/Unicycler/unicycler-runner.py --spades_path /software/SPAdes/3.11.1/ivybridge/bin/spades.py --racon_path /gpfs/project/dilthey/software/racon/bin/racon --pilon_path /software/pilon/1.22/pilon-1.22.jar -t 2 -1 Sim_Pipeline/Plasmid1_l8500_c150/Plasmid1_l8500_c150-filtered_R1.fastq -2 Sim_Pipeline/Plasmid1_l8500_c150/Plasmid1_l8500_c150-filtered_R2.fastq -l example1_plasmid-Plasmid1-predicted_reads.fastq -o example1_plasmid-Plasmid1-predicted_reads_unicycler
```
Important: Pathways for unicycler, spades, racon and pilon need to be replaced to fit your installations.

#### Input:
- Path to Unicycler (/gpfs/project/dilthey/software/Unicycler/unicycler-runner.py)
- Path to SPAdes (/software/SPAdes/3.11.1/ivybridge/bin/spades.py)
- Path to Racon (/gpfs/project/dilthey/software/racon/bin/racon)
- Path to Pilon (/software/pilon/1.22/pilon-1.22.jar)
- Number of threads used (2)
- Path to the Illumina R1 read file (Sim_Pipeline/Plasmid1_l8500_c150/Plasmid1_l8500_c150-filtered_R1.fastq)
- Path to the Illumina R2 read file (Sim_Pipeline/Plasmid1_l8500_c150/Plasmid1_l8500_c150-filtered_R2.fastq)
- Path to the long-read file (example1_plasmid-Plasmid1-predicted_reads.fastq)
- Prefix for the output folder (example1_plasmid-Plasmid1-predicted_reads_unicycler)

#### Output:
- Assembly-folder for each assembled genome, for example called “example1_plasmid-Plasmid1-predicted_reads_unicycler”, containing:
  - De-bruijn (Bandage) graphs for all assembly steps (“001_best_spades_graph.gfa” to “007_rotated.gfa” and “assembly.gfa”).
  - Assembly file in fasta format (assembly.fasta)
  - Log file with parameters of the assembly run (unicycler.log)

### 8. Comparing assemblies of true and predicted reads using nucmer (1s runtime, 1 CPUs, <1gb used memory)

(examplary for `Plasmid1`)

#### Requirements:
- nucmer
- delta-filter

#### Calls:
1.
```
nucmer -p Plasmid1 example1_plasmid-Plasmid1-predicted_reads_unicycler/assembly.fasta example1_plasmid-Plasmid1-true_reads_unicycler/assembly.fasta
```
2.
```
mummerplot2 Plasmid1.delta --png -p Plasmid1
```
#### Input:
1.
- Prefix for the run (Plasmid1)
- Predicted read assembly (example1_plasmid-Plasmid1-predicted_reads_unicycler/assembly.fasta)
- True read assembly  (example1_plasmid-Plasmid1-true_reads_unicycler/assembly.fasta)
2.
- Delta file from previous nucmer run (Plasmid1.delta)
- Prefix for the run (Plasmid1)

#### Output:
1.
- Alignment of true read assembly and predicted read assembly (Plasmid1.delta)
2.
- Graphical representation of the delta file (Plasmid1.png)

If you take a look at the produced graphics, you will see, that the assemblies of the predicted reads align perfectly to the assemblies of the true reads.
