#!/usr/bin/perl
#

use Getopt::Long;
use strict;
use warnings;
use FindBin;

my($this_dir, $prefix, $threshold);

$this_dir = $FindBin::RealBin;										# save directory of this script
GetOptions ('prefix:s' => \$prefix,'threshold:s' => \$threshold);	# get parameters from wrapper script

unless($prefix){die "\nERROR:\nPlease specify a prefix for your run using -p.\n\n\n\n"}
unless($threshold){$threshold=20}	# threshold is per default set to 20, if user does not specify this variable

open(LOG,'>>',"$prefix.log");										# re-open the log file

print "\n\n". `date` ."> Creating R plots\n";
print LOG "\n\n". `date` ."> Start plotting stats and MST\n";
print LOG "\n\n> R default output (only needed for troubleshooting):\n----------------------------------------\n";

### Run R scripts
system("Rscript $this_dir/NanoCore_plots.R Output_$prefix\_Stats/$prefix.min_cov Output_$prefix\_Stats/$prefix.heat_excluded_heterozyg Output_$prefix\_Stats/$prefix.heat_excluded_ave_cov Output_$prefix\_Stats/$prefix.heat_excluded_mapQ Output_$prefix\_Stats/$prefix.table Output_$prefix\_Stats/$prefix.read_length Output_$prefix\_Stats/$prefix-boxplot_min_cov.pdf Output_$prefix\_Stats/$prefix-heatmap_min_cov.pdf Output_$prefix\_Stats/$prefix-heatmap_hete.pdf Output_$prefix\_Stats/$prefix-heatmap_cov.pdf Output_$prefix\_Stats/$prefix-heatmap_mapq.pdf Output_$prefix\_Stats/$prefix-stats_1.pdf Output_$prefix\_Stats/$prefix-stats_2.pdf >> $prefix\.log 2>&1") and die "\nCould not execute R scripts. Please check if the corresponding conda environment was installed correctly and the R executables are in place.\n\n";		#creating boxplot, heatmaps and stat plots
system("Rscript $this_dir/NanoCore_mst.R Output_$prefix\_Tables/$prefix\_allele_table.txt Output_$prefix\_Stats/$prefix-mst.pdf >> $prefix\.log 2>&1") and die "\nCould not execute R scripts. Please check if the corresponding conda environment was installed correctly and the R executables are in place.\n\n";		#creating mst
print LOG "----------------------------------------\n";
system("cp Output_$prefix\_Tables/$prefix\_allele_table.txt $prefix\_allele_table.txt");				# copy distance matrix to main folder
system("cp Output_$prefix\_Stats/$prefix-mst.pdf $prefix-mst.pdf");										# copy MST to main folder

### Write output files to log file and terminal
print "\n\n> Summary of all produced output files written to the end of the log file:\t$prefix.log\n"; 		# Name log file in terminal
print "> Allele based distance matrix of this run written to:\t\t\t$prefix\_allele_table.txt\n"; 			# Name distance matrix file in terminal
print "> Minimum-spanning-tree of this run written to:\t\t\t\t$prefix-mst.pdf\n"; 							# Name MST file in terminal

print LOG "\n\n> Output written to:";
print LOG "\n\tMain output files (run log, distance matrix, MST)";
print LOG "
\t\t$prefix.log\t\t\t\t(Log file of the complete run)
\t\t$prefix\_allele_table.txt\t\t(Distance matrix on basis of alleles. This file was copied from the 'Output_$prefix\_Tables/' folder)
\t\t$prefix-mst.pdf\t\t\t\t(Minimum-spanning-tree of all samples from this run. Dashed edge means no distance. This file was copied from the 'Output_$prefix\_Stats/' folder)
";
print LOG "\n\tOutput_$prefix\_Minimap/[SAMPLE_ID]/";
print LOG "
\t\t[SAMPLE_ID].bam\t\t\t\t\t(Original mapping output)
\t\t[SAMPLE_ID]_new_FLAGs.bam\t\t\t(Correctd mapping output)
\t\t[SAMPLE_ID]_new_FLAGs.bam.bai\t\t\t(Index file to mapping output)
\t\t[SAMPLE_ID]_new_FLAGs_mpileup.txt\t\t(Samtools mpileup file to mapping output)
";
print LOG "\n\tOutput_$prefix\_Clair/[SAMPLE_ID]/";
print LOG "
\t\t[SAMPLE_ID]_merge_output.vcf\t\t\t(Variant calling output)
\t\trun_clair3.log\t\t\t\t\t(Variant calling log file)
\t\tOther files\t\t\t\t\t(Secondary clair3 output files)
";
print LOG "\n\tOutput_$prefix\_Stats/";
print LOG "
\t\t[SAMPLE_ID].coverage\t\t\t\t(Samtools coverage output)
\t\t[SAMPLE_ID].depth\t\t\t\t(Samtools depth output)
\t\t[SAMPLE_ID].stats\t\t\t\t(Samtools stats output)
\t\t[SAMPLE_ID].coverage_per_ref\t\t\t(Mean coverage and ratio of positions above coverage threshold per gene)
\t\t$prefix.read_length\t\t\t(Lengths of all reads in that run)
\t\t$prefix.min_cov\t\t\t\t(Matrix of ratio of positions per gene above coverage threshold of $threshold)
\t\t$prefix.table\t\t\t\t(Table with basic read statistics)
\t\t$prefix.excluded_1_heterozyg\t\t(Excluded genes per isolate pair with too high ratio of heterozygous calls)
\t\t$prefix.excluded_2_cov_ratio\t\t(Excluded genes per isolate pair with too low ratio of positions above coverage threshold)
\t\t$prefix.excluded_3_same_var\t\t(Genes per isolate pair with the same variants -> no distance)
\t\t$prefix.excluded_4_mapq_plus_coverage\t(Excluded genes per isolate pair with unusual coverage AND low mapQ)
\t\t$prefix.excluded_5_no_var\t\t(Genes per isolate pair without variants)
\t\t$prefix.excluded_6_mpileup\t\t(Genes excluded through mpileup statistics)
\t\t$prefix.excluded_genes_per_isolate\t(Summary of excluded genes per isolate)
\t\t$prefix.suspos\t\t\t\t(Suspicious positions in genes that are excluded from the complete analysis due to heterozygozity, low quality or low allele frequenzy)
\t\t$prefix.variantlist\t\t\t(List of called variants per isolate pair)
\t\t$prefix.heat_excluded_ave_cov\t\t(Matrix for heatmap of genes with unusual coverage) 
\t\t$prefix.heat_excluded_heterozyg\t\t(Matrix for heatmap of genes with unusual heterozygozity)
\t\t$prefix.heat_excluded_mapQ\t\t(Matrix for heatmap of genes with unusual mapping quality)
\t\t$prefix-boxplot_min_cov.pdf\t\t(Boxplot of ratio of positions per gene above coverage threshold of $threshold)
\t\t$prefix-heatmap_min_cov.pdf\t\t(Heatmap of ratio of positions per gene above coverage threshold of $threshold)
\t\t$prefix-stats_1.pdf\t\t\t(Plots of basic read statistics part 1/2)
\t\t$prefix-stats_2.pdf\t\t\t(Plots of basic read statistics part 2/2)
\t\t$prefix-mst.pdf\t\t\t\t(Minimum-spanning-tree of all samples from this run. Dashed edge means no distance)
\t\t$prefix-heatmap_cov.pdf\t\t\t(Heatmap of genes with unusual coverage)
\t\t$prefix-heatmap_hete.pdf\t\t(Heatmap of genes with unusual heterozygozity)
\t\t$prefix-heatmap_mapq.pdf\t\t(Heatmap of genes with unusual mapping quality)
";	
print LOG "\n\tOutput_$prefix\_Tables/";
print LOG "
\t\t$prefix\_variant_table.txt\t\t(Distance matrix on basis of variants)
\t\t$prefix\_allele_table.txt\t\t(Distance matrix on basis of alleles)
\t\t$prefix\_consensus_table.txt\t\t(Similarity matrix on basis of variants)
\t\t$prefix\_per_gene_variant_table.txt\t(Table of variants per gene per isolate pair)
\t\t$prefix\_per_gene_total_table.txt\t(Total variants per gene over all isolates)
";

close(LOG);