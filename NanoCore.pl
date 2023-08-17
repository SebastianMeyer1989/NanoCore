#!/usr/bin/perl
#
# run example:   NanoCore.pl --sample_list sample_list.txt --reference S.aureus_cgMLST_ref-seqs.fasta --clair_model_nano /home/user/Software/miniconda3/envs/clair3/bin/models/ont --clair_model_illu /home/user/Software/miniconda3/envs/clair3/bin/models/ilmn --threshold 20 --threads 8 --samtools samtools --prefix NanoCore_Run_1
#

use Getopt::Long;
use strict;
use warnings;

my(@line_fields, %variants, %count, %variantcount, %allelecount, @genelist, @keys3iso1, @keys3iso2, %consensecount, %pergenevariantcount, %pergenetotalcount, %min_cov,%ave_cov, %hetecount, %samplecount, %mapQ, %mean_cov, %ass_cov, %completevariants, @linearray, %counttwo, @samplearray, %af, %suspos, %position, %alternative, @counter_fields, %mpileup_one, %mpileup_two, %exgene1, %exgene2, %exgene3);		# initiate arrays and hashes
my($IDandPATH, $line, $sampleID, $geneID, $pos, $refAllele, $altAllele, $item, $itempair, $iso1, $iso2, $gene, $prefix, $clairout, $path, $clair_illu, $clair_nano, $sample_list, $reference, $threshold, $depth, $id, $var_pos, $cov, $covcounter, $nano_path, $v1, $v2, $v3, $v4, $v5, $read, $coverage, $assembly, $samtoolsbin, $bam, $bamout, $flag, $threads, $quality, $mpileupout, $alt1, $alt2, $counter, $field, $genecov, $mpileup, $switch, $mpu, $arraycounter, $illu_path_1, $illu_path_2);		# initiate variables

GetOptions (
	'sample_list:s' => \$sample_list,
	'reference:s' => \$reference,
	'clair_model_nano:s' => \$clair_nano,
	'clair_model_illu:s' => \$clair_illu,
	'prefix:s' => \$prefix,
	'threshold:s' => \$threshold,
	'threads:s' => \$threads,
	'samtools:s' => \$samtoolsbin,
);

unless($sample_list){die "\nERROR:\nPlease specify a sample list using --sample_list.\n\n"}
unless($reference){die "\nERROR:\nPlease specify a reference using --reference.\n\n"}
unless($clair_nano && $clair_illu){die "\nERROR:\nPlease specify a path to the clair3 nanopore and illumina models using --clair_model_nano and --clair_model_illu. (Should be somewhere here in your miniconda3 path: /miniconda3/envs/clair3/bin/models/ont and /miniconda3/envs/clair3/bin/models/ilmn)\n\n"}
unless($prefix){die "\nERROR:\nPlease specify a prefix for your run using --prefix.\n\n"}
unless($threshold){$threshold=20}	# threshold is per default set to 20, if user does not specify this variable
unless($threads){die "\nERROR:\nPlease specify the number of threads you want to use using --threads (recomended: at least 8).\n\n"}
unless($samtoolsbin){$samtoolsbin="samtools"}	# samtools binary is per default set to "samtools", if user does not specify this variable
#unless($samtoolsbin){die "\nERROR:\nPlease specify your samtools executable using --samtools (In most cases this should be just 'samtools').\n\n"}


#####
##### Minimap plus Clair #####
#####

if(-e "Output_$prefix\_Minimap"){}
else{
	system("mkdir Output_$prefix\_Minimap");	# create Minimap folder
}
if(-e "Output_$prefix\_Clair"){}
else{
	system("mkdir Output_$prefix\_Clair");	# create Minimap folder
}

system("$samtoolsbin faidx $reference");

open(IDS,'>',"isolate_list_".$prefix.".tmp");
open(SAMPLELIST,$sample_list) or die "\nSample list not found.\n\n";	# open list of IDs and files (ID\tPATH)
print "\n\n> Mapping and variant calling samples:\n";
sleep(1);
while($IDandPATH = <SAMPLELIST>){
	chomp($IDandPATH);
	if($IDandPATH eq ""){}								# ignore empty lines
	else{
		@line_fields = split("\t", $IDandPATH);			# split line into ID and path
		$arraycounter=scalar(split("\t", $IDandPATH));	# count number of strings in line
### Processing of Nanopore Data ###
		if($line_fields[1] eq "Nanopore" && $arraycounter==3){
			$sampleID = $line_fields[0];				# save sample ID	### TO DO: Unique IDs needed!
			$nano_path = $line_fields[2];				# save path
			print IDS $sampleID."\n";

			if(-e "Output_$prefix\_Minimap/$sampleID"){
				print "\n\tSample ".$sampleID." already mapped\n";
				select(undef,undef,undef, .3);
			}
			else{
				system("mkdir Output_$prefix\_Minimap/$sampleID");
				print "\n\tMapping sample ".$sampleID."\n\n";
				$bam="Output_$prefix\_Minimap/$sampleID/$sampleID.bam";
				$bamout="Output_$prefix\_Minimap/$sampleID/$sampleID\_new_FLAGs.bam";
				$mpileupout="Output_$prefix\_Minimap/$sampleID/$sampleID\_new_FLAGs_mpileup.txt";
				
				system("minimap2 -t $threads -x map-ont -a $reference $nano_path | $samtoolsbin sort --reference $reference --threads $threads -O BAM - > $bam") and die;	# run minimap
				
				open(SAMTOOLS_IN,"$samtoolsbin view -h $bam |") or die;
				open(SAMTOOLS_OUT,"| $samtoolsbin view -bo $bamout -") or die;

				while($line=<SAMTOOLS_IN>){
					chomp($line);
					if((substr($line, 0, 1) eq '@') || ($line eq "")){print SAMTOOLS_OUT $line."\n"}		# exclude header and print to file
					else{
						@linearray=split("\t", $line, -1);	# split alignment
						$flag=$linearray[1]+0;				# save flag as number
						if(($flag & 2048) != 0){$linearray[1]=$flag ^ 2048}	# check if FLAG is supplementary alignment, replace if yes ($ = keeps bits in both numbers, ^ = keeps bits NOT in both numbers)
						$line=join("\t", @linearray);		# rejoin alignment
						
						print SAMTOOLS_OUT $line."\n";		# print to new file
					}
				}
				close(SAMTOOLS_IN);
				close (SAMTOOLS_OUT);
				
				system("$samtoolsbin index $bamout");	# index bam file
				print "\n\tRunning mpileup on sample ".$sampleID."\n\n";
				system("perl extractPileUpFrequencies.pl --samtools_bin $samtoolsbin --BAM $bamout --outputFile $mpileupout --referenceGenome $reference");	# run mpileeup on bam file
			}
			if(-e "Output_$prefix\_Clair/$sampleID"){
				print "\tSample ".$sampleID." already called\n";
				select(undef,undef,undef, .3);
			}
			else{
				system("mkdir Output_$prefix\_Clair/$sampleID");
				print "\n\tCalling sample ".$sampleID."\n\n";
				system("run_clair3.sh -b Output_$prefix\_Minimap/$sampleID/$sampleID\_new_FLAGs.bam -f $reference -t $threads -p ont -m $clair_nano -o Output_$prefix\_Clair/$sampleID --include_all_ctgs") and die;	# run clair
				system("gunzip Output_$prefix\_Clair/$sampleID/merge_output.vcf.gz") and die;	# unzip Clair output
				system("mv Output_$prefix\_Clair/$sampleID/merge_output.vcf Output_$prefix\_Clair/$sampleID/$sampleID\_merge_output.vcf") and die
			}
		}
### Processing of Illumina Data ###
		elsif($line_fields[1] eq "Illumina" && $arraycounter==4){
			$sampleID = $line_fields[0];		# save sample ID	### TO DO: Unique IDs needed!
			$illu_path_1=$line_fields[2];		# save path
			$illu_path_2=$line_fields[3];		# save path
			print IDS $sampleID."\n";

			if(-e "Output_$prefix\_Minimap/$sampleID"){
				print "\n\tSample ".$sampleID." already mapped\n";
				select(undef,undef,undef, .3);
			}
			else{
				system("mkdir Output_$prefix\_Minimap/$sampleID");
				print "\n\tMapping sample ".$sampleID."\n\n";
				$bam="Output_$prefix\_Minimap/$sampleID/$sampleID.bam";
				$bamout="Output_$prefix\_Minimap/$sampleID/$sampleID\_new_FLAGs.bam";
				$mpileupout="Output_$prefix\_Minimap/$sampleID/$sampleID\_new_FLAGs_mpileup.txt";
				
				system("minimap2 -t $threads -x sr -a $reference $illu_path_1 $illu_path_2 | $samtoolsbin sort --reference $reference --threads $threads -O BAM - > $bam") and die;	# run minimap
				
				open(SAMTOOLS_IN,"$samtoolsbin view -h $bam |") or die;
				open(SAMTOOLS_OUT,"| $samtoolsbin view -bo $bamout -") or die;

				while($line=<SAMTOOLS_IN>){
					chomp($line);
					if((substr($line, 0, 1) eq '@') || ($line eq "")){print SAMTOOLS_OUT $line."\n"}		# exclude header and print to file
					else{
						@linearray=split("\t", $line, -1);	# split alignment
						$flag=$linearray[1]+0;				# save flag as number
						if(($flag & 2048) != 0){$linearray[1]=$flag ^ 2048}	# check if FLAG is supplementary alignment, replace if yes ($ = keeps bits in both numbers, ^ = keeps bits NOT in both numbers)
						$line=join("\t", @linearray);		# rejoin alignment
						
						print SAMTOOLS_OUT $line."\n";		# print to new file
					}
				}
				close(SAMTOOLS_IN);
				close (SAMTOOLS_OUT);
				
				system("$samtoolsbin index $bamout");	# index bam file
				print "\n\tRunning mpileup on sample ".$sampleID."\n\n";
				system("perl extractPileUpFrequencies.pl --samtools_bin $samtoolsbin --BAM $bamout --outputFile $mpileupout --referenceGenome $reference");	# run mpileeup on bam file
			}
			if(-e "Output_$prefix\_Clair/$sampleID"){
				print "\tSample ".$sampleID." already called\n";
				select(undef,undef,undef, .3);
			}
			else{
				system("mkdir Output_$prefix\_Clair/$sampleID");
				print "\n\tCalling sample ".$sampleID."\n\n";
				system("run_clair3.sh -b Output_$prefix\_Minimap/$sampleID/$sampleID\_new_FLAGs.bam -f $reference -t $threads -p ont -m $clair_illu -o Output_$prefix\_Clair/$sampleID --include_all_ctgs") and die;	# run clair
				system("gunzip Output_$prefix\_Clair/$sampleID/merge_output.vcf.gz") and die;	# unzip Clair output
				system("mv Output_$prefix\_Clair/$sampleID/merge_output.vcf Output_$prefix\_Clair/$sampleID/$sampleID\_merge_output.vcf") and die
			}
		}
		else{die "\nSample list line\n\n$IDandPATH\n\n has not the right format. \nThree entries needed for Nanopore data ([ID], 'Nanopore', [NANOPORE_PATH]).\nFour entries needed for Illumina data ([ID], 'Illumina', [ILLUMINA_PATH_1]), [ILLUMINA_PATH_2]).\n\n"}
	}
}
close(SAMPLELIST);


#####
##### Basic stats calculation #####
#####

### save gene list as array ###
system("less $reference  | grep '>' | sed 's/>//g' > gene_list_".$prefix.".tmp") and die;	# save reference gene names as list
open(GENELIST,"gene_list_".$prefix.".tmp") or die "\nGene list not found.\n\n";		# open gene list
while($line = <GENELIST>){
		chomp($line);
		push(@genelist,$line);											# push gene list into an array
}
close(GENELIST);

### open output files, run samtools, filter info from samtools output ###
if(-e "Output_$prefix\_Stats"){}
else{
	system("mkdir Output_$prefix\_Stats");											# create Stats folder
}

sleep(1);
print "\n\n> Calculating basic stats:\n";
sleep(1);

open(TABLE,'>',"Output_$prefix\_Stats/$prefix.min_cov");				# create output table file with ratio of positions above threshold (for boxplot / heatmap)
foreach $gene (@genelist){	
	print TABLE "\t".$gene;													# print gene IDs as header
}
open(STATSOUTTABLE,'>',"Output_$prefix\_Stats/$prefix.table");								# create output table file with data for basic stats
print STATSOUTTABLE "ID\tTotal_Bases\tMapped_Bases\tTotal_Reads\tMapped_Reads\tAverage_Coverage";	# print header into file

open(READS,'>',"Output_$prefix\_Stats/$prefix.read_length");		# create output table file with data for histogram
print READS "Length_of_all_reads\n";										# print header into file

open(CLAIROUT,$sample_list) or die "\nSample list not found.\n\n";	# open list of IDs and files (ID\tPATH)
while($IDandPATH = <CLAIROUT>){
	chomp($IDandPATH);
	if($IDandPATH eq ""){}							# ignore empty lines
	else{
		@line_fields = split("\t", $IDandPATH);		# split line into ID and path
### Processing of Nanopore Data ###
		if($line_fields[1] eq "Nanopore"){
			$clairout = $line_fields[0];	# save sample ID
			$nano_path = $line_fields[2];	# save path
			print "\t".$clairout."\n";		# print sample ID to output
			print TABLE "\n".$clairout;		# print sample ID as first column into table file
			print STATSOUTTABLE "\n".$clairout;		# print sample ID as first column into table file
			open(STATSOUT,'>',"Output_$prefix\_Stats/$clairout.coverage_per_ref");								# create Stats summary output file
			print "\t\tSamtools stats\n";
			system("$samtoolsbin stats Output_$prefix\_Minimap/$clairout/$clairout.bam > Output_$prefix\_Stats/$clairout.stats") and die;	# calculate samtools stats
			print "\t\tSamtools coverage\n";
			system("$samtoolsbin coverage Output_$prefix\_Minimap/$clairout/$clairout.bam > Output_$prefix\_Stats/$clairout.coverage") and die;	# calculate samtools coverage
			print "\t\tSamtools depth\n";
			system("$samtoolsbin depth -aa Output_$prefix\_Minimap/$clairout/$clairout.bam > Output_$prefix\_Stats/$clairout.depth") and die;	# calculate samtools depth
					
			open(COVERAGE, "Output_$prefix\_Stats/$clairout.coverage") or die "\nCoverage file not found.\n\n";	# 1,2
			while($coverage = <COVERAGE>){																				# 1,2
				chomp($coverage);																						# 1,2
				if(($coverage !~ m/^#/) && ($coverage =~ m/(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)/)){		# 1,2 parse coverage file
					$mean_cov{$clairout}{$1}=$7;		# 1: filter; save mean cov per gene
					$mapQ{$clairout}{$1}=$9;			# 2: filter; save mapQ pere gene
				}										# 1,2
			}											# 1,2
				
			print "\t\tSummarise stats\n";

			$v1= `grep -v - $nano_path`;																# extract total bases
			chomp($v1);																					# chomp total bases
			print STATSOUTTABLE "\t".length($v1);														# print total bases
			
			$v2=`less Output_$prefix\_Stats/$clairout.stats | grep 'bases mapped:' | grep -Eo [0-9]+`;	# extract mapped bases
			chomp($v2);																					# chomp mapped bases
			print STATSOUTTABLE "\t".$v2 ;																# print mapped bases
			
			$v3=0;
			open(FILE, $nano_path) or die "File $nano_path not found.\n\n";								# extract total reads and individual reads lengths
			while($read = <FILE>){
				chomp($read);
				$v3++;
				if($v3%4==2){
					print READS length($read)."\n";														# print individual reads lengths
				}
			}
			chomp($v3);				;																	# chomp total reads
			print STATSOUTTABLE "\t".$v3/4;																# print total reads
			
			$v4=`less Output_$prefix\_Stats/$clairout.stats | grep 'reads mapped:' | grep -Eo [0-9]+`;	# extract mapped reads
			chomp($v4);																					# chomp mapped reads
			print STATSOUTTABLE "\t".$v4;																# print mapped reads
			
			$v5=`awk '{ total += \$7; count++ } END { print total/(count-1) }' Output_$prefix\_Stats/$clairout.coverage`;	# calculate average gene coverage
			chomp($v5);																										# chomp average gene coverage
			print STATSOUTTABLE "\t".$v5;																					# print average gene coverage
			$ave_cov{$clairout}=$v5;																	# ----- save average coverage
			
			print STATSOUT "Coverage per reference gene:\n".`awk '{ print \$1, \$7 }' Output_$prefix\_Stats/$clairout.coverage | sed 's/ /\t/g'`;			# extract per gene coverage
			print STATSOUT "\nRatio of positions above coverage threshold per reference gene:\n"; 															# count positions with coverage greater x
			open(DEPTH,"Output_$prefix\_Stats/$clairout.depth") or die "\nSamtools depth output not found.\n\n";		# open list of depth output
			$covcounter=0;				# reset counter
			$id="";						# reset ID
			while($depth = <DEPTH>){
				chomp($depth);
				if($depth =~ m/#/){}			#exclude comment lines
				else{
					if($depth =~ m/(.+)\t(.+)\t(.+)/){  		# filter for format [GENE_ID POSITION COVERAGE]
						if($id eq ""){							# on first line of depth file
							$id = $1;
							$var_pos = $2;
							$cov = $3;
							if($cov>=$threshold){$covcounter++};	# count positions above threshold
						}
						elsif($1 eq $id){							# each line with same gene ID as before
							$var_pos = $2;
							$cov = $3;
							if($cov>=$threshold){$covcounter++};	# count positions above threshold
						}
						else{										# when switching to the next gene ID
							print STATSOUT $id."\t".sprintf("%.3f", ($covcounter/$var_pos))."\n";	# print values to table output file
							print TABLE "\t".sprintf("%.5f", ($covcounter/$var_pos));				# print values to table output file for heatmap
							$min_cov{$clairout}{$id}=sprintf("%.3f", ($covcounter/$var_pos));		# save values in hash
							$id = $1;
							$var_pos = $2;
							$cov = $3;
							$covcounter = 0;
							if($cov>=$threshold){$covcounter++};	# count positions above threshold
						}
					}
				}
			}
			print STATSOUT $id."\t".sprintf("%.3f", ($covcounter/$var_pos))."\n";	# print files last values to table output file
			print TABLE "\t".sprintf("%.5f", ($covcounter/$var_pos));				# print files last values to table output file for heatmap
			$min_cov{$clairout}{$id}=sprintf("%.3f", ($covcounter/$var_pos));		# save values in hash
		}
### Processing of Illumina Data ###
		elsif($line_fields[1] eq "Illumina"){
			$clairout = $line_fields[0];	# save sample ID
			$illu_path_1 = $line_fields[2];	# save path
			$illu_path_2 = $line_fields[3];	# save path
			print "\t".$clairout."\n";		# print sample ID to output
			print TABLE "\n".$clairout;		# print sample ID as first column into table file
			system("cat $illu_path_1 $illu_path_2 > illumina_complete_tmp.fastq");	# concateneate Illumina input files
			print STATSOUTTABLE "\n".$clairout;		# print sample ID as first column into table file
			open(STATSOUT,'>',"Output_$prefix\_Stats/$clairout.coverage_per_ref");								# create Stats summary output file
			print "\t\tSamtools stats\n";
			system("$samtoolsbin stats Output_$prefix\_Minimap/$clairout/$clairout.bam > Output_$prefix\_Stats/$clairout.stats") and die;	# calculate samtools stats
			print "\t\tSamtools coverage\n";
			system("$samtoolsbin coverage Output_$prefix\_Minimap/$clairout/$clairout.bam > Output_$prefix\_Stats/$clairout.coverage") and die;	# calculate samtools coverage
			print "\t\tSamtools deepth\n";
			system("$samtoolsbin depth -aa Output_$prefix\_Minimap/$clairout/$clairout.bam > Output_$prefix\_Stats/$clairout.depth") and die;	# calculate samtools depth
					
			open(COVERAGE, "Output_$prefix\_Stats/$clairout.coverage") or die "\nCoverage file not found.\n\n";	# 1,2
			while($coverage = <COVERAGE>){																				# 1,2
				chomp($coverage);																						# 1,2
				if(($coverage !~ m/^#/) && ($coverage =~ m/(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)\t(.+)/)){		# 1,2 parse coverage file
					$mean_cov{$clairout}{$1}=$7;		# 1: filter; save mean cov per gene
					$mapQ{$clairout}{$1}=$9;			# 2: filter; save mapQ pere gene
				}										# 1,2
			}											# 1,2
				
			print "\t\tSummarise stats\n";

			$v1= `grep -v - illumina_complete_tmp.fastq`;																# extract total bases
			chomp($v1);																					# chomp total bases
			print STATSOUTTABLE "\t".length($v1);														# print total bases
			
			$v2=`less Output_$prefix\_Stats/$clairout.stats | grep 'bases mapped:' | grep -Eo [0-9]+`;	# extract mapped bases
			chomp($v2);																					# chomp mapped bases
			print STATSOUTTABLE "\t".$v2 ;																# print mapped bases
			
			$v3=0;
			open(FILE, "illumina_complete_tmp.fastq") or die "File illumina_complete_tmp.fastq not found.\n\n";								# extract total reads and individual reads lengths
			while($read = <FILE>){
				chomp($read);
				$v3++;
				if($v3%4==2){
					print READS length($read)."\n";														# print individual reads lengths
				}
			}
			chomp($v3);				;																	# chomp total reads
			print STATSOUTTABLE "\t".$v3/4;																# print total reads
			
			$v4=`less Output_$prefix\_Stats/$clairout.stats | grep 'reads mapped:' | grep -Eo [0-9]+`;	# extract mapped reads
			chomp($v4);																					# chomp mapped reads
			print STATSOUTTABLE "\t".$v4;																# print mapped reads
			
			$v5=`awk '{ total += \$7; count++ } END { print total/(count-1) }' Output_$prefix\_Stats/$clairout.coverage`;	# calculate average gene coverage
			chomp($v5);																										# chomp average gene coverage
			print STATSOUTTABLE "\t".$v5;																					# print average gene coverage
			$ave_cov{$clairout}=$v5;																	# ----- save average coverage
			
			print STATSOUT "Coverage per reference gene:\n".`awk '{ print \$1, \$7 }' Output_$prefix\_Stats/$clairout.coverage | sed 's/ /\t/g'`;			# extract per gene coverage
			print STATSOUT "\nRatio of positions above coverage threshold per reference gene:\n"; 															# count positions with coverage greater x
			open(DEPTH,"Output_$prefix\_Stats/$clairout.depth") or die "\nSamtools depth output not found.\n\n";		# open list of depth output
			$covcounter=0;				# reset counter
			$id="";						# reset ID
			while($depth = <DEPTH>){
				chomp($depth);
				if($depth =~ m/#/){}			#exclude comment lines
				else{
					if($depth =~ m/(.+)\t(.+)\t(.+)/){  		# filter for format [GENE_ID POSITION COVERAGE]
						if($id eq ""){							# on first line of depth file
							$id = $1;
							$var_pos = $2;
							$cov = $3;
							if($cov>=$threshold){$covcounter++};	# count positions above threshold
						}
						elsif($1 eq $id){							# each line with same gene ID as before
							$var_pos = $2;
							$cov = $3;
							if($cov>=$threshold){$covcounter++};	# count positions above threshold
						}
						else{										# when switching to the next gene ID
							print STATSOUT $id."\t".sprintf("%.3f", ($covcounter/$var_pos))."\n";	# print values to table output file
							print TABLE "\t".sprintf("%.5f", ($covcounter/$var_pos));				# print values to table output file for heatmap
							$min_cov{$clairout}{$id}=sprintf("%.3f", ($covcounter/$var_pos));		# save values in hash
							$id = $1;
							$var_pos = $2;
							$cov = $3;
							$covcounter = 0;
							if($cov>=$threshold){$covcounter++};	# count positions above threshold
						}
					}
				}
			}
			print STATSOUT $id."\t".sprintf("%.3f", ($covcounter/$var_pos))."\n";	# print files last values to table output file
			print TABLE "\t".sprintf("%.5f", ($covcounter/$var_pos));				# print files last values to table output file for heatmap
			$min_cov{$clairout}{$id}=sprintf("%.3f", ($covcounter/$var_pos));		# save values in hash
		}
	}
}
print STATSOUTTABLE "\n";
print TABLE "\n";

close(STATSOUT);
close(TABLE);
close(STATSOUTTABLE);

#####
##### Variant comparison and Distmat calculation #####
#####

### save position, ref and alt in hash of hash of hashes ###
open(FILELIST,"isolate_list_".$prefix.".tmp") or die "\nClair output list not found.\n\n";	# open list of Clair output
sleep(1);
print "\n\n> Processing mapped files:\n";
sleep(1);

open(SUSPOS,'>',"Output_$prefix\_Stats/$prefix.suspos");				# create output for suspicious positions
print SUSPOS "Suspicious positions in genes. These are not included in the distance calculation to avoid false positives:\nGene_ID\tPosition\tIssue...";

while($sampleID = <FILELIST>){
	chomp($sampleID);
	print "\t".$sampleID."\n";
	
	foreach $gene (@genelist){				# 3
		$hetecount{$sampleID}{$gene}=0;		# 3: filter; initialize heterozyg counter
	}										# 3
	open(VCF,"Output_$prefix\_Clair/$sampleID/$sampleID\_merge_output.vcf") or die "VCF file not found\n\n";	# open each file#
	while($line = <VCF>){
		chomp($line);
		if($line =~ m/#/){}									#exclude comment lines
		else{
			@line_fields = split("\t", $line); 				# CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
			
			$geneID = $line_fields[0];
			$pos = $line_fields[1];
			$refAllele = $line_fields[3];
			$altAllele = $line_fields[4];
			$quality = $line_fields[6];
			@samplearray =  split(":", $line_fields[9]);			# GT:GQ:DP:AF
			if($samplearray[0] eq "0/1"){												# ignore suspicions positions in allcomparisons! (this and the next 2 ifs)
				$suspos{$geneID}{$pos}++;														# add heterozygous calls to suspicions positions
				print SUSPOS "\n".$geneID."\t".$pos."\tHeterozygous_in_at_least_one_isolate";	# print position to file
				$hetecount{$sampleID}{$geneID}++;												# count heterozygous positions
				$samplecount{$sampleID}{$geneID}++;												# count positions
			}
			elsif($samplearray[3]<=0.5){
				$suspos{$geneID}{$pos}++;														# add low allele frequency calls to suspicions positions
				print SUSPOS "\n".$geneID."\t".$pos."\tAF_below_0.5_in_at_least_one_isolate";	# print position to file
			}
			elsif($quality eq "LowQual"){
				$suspos{$geneID}{$pos}++;														# add low quality calls to suspicions positions
				print SUSPOS "\n".$geneID."\t".$pos."\tLow_quality_in_at_least_one_isolate";	# print position to file
			}
			elsif($quality eq "RefCall"){}														# ignore reference calls (should not be called to begin with)
			elsif($samplearray[2]<$threshold){
				$suspos{$geneID}{$pos}++;														# add low coverage calls to suspicions positions
				print SUSPOS "\n".$geneID."\t".$pos."\tCoverage_below_threshold_in_at_least_one_isolate";	# print position to file
			}
			elsif((length($refAllele)==length($altAllele)) && ($samplearray[2]>=$threshold)){				# ignore inDels and calls below threshold (20) reads 
				$variants{$sampleID}{$geneID}{$pos . ':' . $refAllele . '->' . $altAllele}++; 				# save IDs, positions, variants, hete/homo and quality
				$completevariants{$sampleID}{$geneID}{$pos . ':' . $refAllele . '->' . $altAllele}=$line;	# save whole line
				$position{$sampleID}{$geneID}{$pos . ':' . $refAllele . '->' . $altAllele}=$pos;			# save position (for suspos filter)
				$alternative{$sampleID}{$geneID}{$pos . ':' . $refAllele . '->' . $altAllele}=$altAllele;	# save alt allele (for mpileup filter)
				$samplecount{$sampleID}{$geneID}++;															# count positions
			}
		}
	}
}
close(FILELIST);
close(VCF);

open(HETEROZYG,'>',"Output_$prefix\_Stats/$prefix.excluded_1_heterozyg");				# create output
print HETEROZYG "Genes with high number of heterogygous calls (>50%) in pairwise comparison:\nSample\tGene_IDs\t...\t...";
open(COVRATIO,'>',"Output_$prefix\_Stats/$prefix.excluded_2_cov_ratio");				# create output for excluded genes
print COVRATIO "Excluded genes (At leat 10% of the gene positions have coverage below set threshold: ".$threshold."x.):\nSample_pair\tGene_IDs\t...\t...";
open(SAMEVAR,'>',"Output_$prefix\_Stats/$prefix.excluded_3_same_var");					# create output for genes without variants
print SAMEVAR "Genes of sample-pairs with same variants (no distancee):\nSample_pair\tGene_IDs\t...\t...";
open(MAPQCOV,'>',"Output_$prefix\_Stats/$prefix.excluded_4_mapq_plus_coverage");		# create output for excluded genes
print MAPQCOV "Excluded genes (Gene coverage differs more than 25% from average AND mapQ is below 55):\nSample_pair\tGene_IDs\t...\t...";
open(NOVAR,'>',"Output_$prefix\_Stats/$prefix.excluded_5_no_var");						# create output for genes without variants
print NOVAR "Genes of sample-pairs without variants (no distance):\nSample_pair\tGene_IDs\t...\t...";
open(MPILEUPOUT,'>',"Output_$prefix\_Stats/$prefix.excluded_6_mpileup");					# create output for gene positions that are excluded from distance counting due to mpileup statistics
print MPILEUPOUT "Genes positions that are excluded from distance counting due to mpileup statistics:\nSample_pair\tGene_IDs\tPosition\tAF_of_alt_base";
open(VARIANTLIST,'>',"Output_$prefix\_Stats/$prefix.variantlist");					# create output for found variants
print VARIANTLIST "Differing variants between sample pairs:\nSourcee\tSample_pair\tGene_ID\tPos\tID\tRef\tAlt\tQuality\tFilter\tInfo\tFormat\tSample";

	open(HETEROZYGFTWO,'>',"Output_$prefix\_Stats/$prefix.tmp2_excluded_heterozyg");				# create heatmap input 1
	open(AVECOVFTWO,'>',"Output_$prefix\_Stats/$prefix.tmp2_excluded_ave_cov");					# create heatmap input 2
	open(MAPQFTWO,'>',"Output_$prefix\_Stats/$prefix.tmp2_excluded_mapQ");						# create heatmap input 3
	
	foreach $gene (@genelist){
		print HETEROZYGFTWO "\t".$gene;	# print gene IDs as header 1
		print AVECOVFTWO "\t".$gene;	# print gene IDs as header 2
		print MAPQFTWO "\t".$gene;		# print gene IDs as header 3
	}
### pairwise gene comparison of each isolate ###
sleep(1);
print "\n\n> Starting pairwise comparison:\n";
sleep(1);
foreach $iso1 (sort(keys(%variants))){	
	undef(%mpileup_one);								#empty the mpileup hash
	open(MPILEUPONE, "Output_$prefix\_Minimap/$iso1/$iso1\_new_FLAGs_mpileup.txt") or die "\nMpileup file of $iso1 not found.\n\n";	# open list of IDs and files (ID\tPATH)
	while($mpu = <MPILEUPONE>){
		chomp($mpu);
		if($mpu eq ""){}							# ignore empty lines
		else{
			@line_fields = split("\t", $mpu);		# split line into array
			$mpileup_one{$line_fields[0]}{$line_fields[1]}=$line_fields[3]; # save GeneID, Position and counting entries
		}
	}									# 4
	foreach $iso2 (sort(keys(%variants))){
		undef(%mpileup_two);								#empty the mpileup hash
		open(MPILEUPTWO, "Output_$prefix\_Minimap/$iso2/$iso2\_new_FLAGs_mpileup.txt") or die "\nMpileup file of $iso2 not found.\n\n";	# open list of IDs and files (ID\tPATH)
		while($mpu = <MPILEUPTWO>){
			chomp($mpu);
			if($mpu eq ""){}							# ignore empty lines
			else{
				@line_fields = split("\t", $mpu);		# split line into array
				$mpileup_two{$line_fields[0]}{$line_fields[1]}=$line_fields[3]; # save GeneID, Position and counting entries
			}
		}
		print "\t".$iso1." vs. ".$iso2."\n";

			print HETEROZYGFTWO "\n".$iso1;			# 1
			print AVECOVFTWO "\n".$iso1;			# 2
			print MAPQFTWO "\n".$iso1;				# 3
		print HETEROZYG "\n".$iso1."_vs_".$iso2;
		print COVRATIO "\n".$iso1."_vs_".$iso2;
		print SAMEVAR "\n".$iso1."_vs_".$iso2;
		print MAPQCOV "\n".$iso1."_vs_".$iso2;
		print NOVAR "\n".$iso1."_vs_".$iso2;
	$variantcount{$iso1}{$iso2}=0;
		foreach $gene (@genelist){		# for each possible gene from gene list
					if($hetecount{$iso1}{$gene} > 0 && exists($samplecount{$iso1}{$gene}) && ($hetecount{$iso1}{$gene}/$samplecount{$iso1}{$gene}) > 0.5){print HETEROZYGFTWO "\t1"} else{print HETEROZYGFTWO "\t0"}												# 1		#heterozyg >50%	
					if((exists($mean_cov{$iso1}{$gene}) && $mean_cov{$iso1}{$gene} >= ($ave_cov{$iso1}*1.25)) || (exists($mean_cov{$iso1}{$gene}) && $mean_cov{$iso1}{$gene} <= ($ave_cov{$iso1}*0.75))){print AVECOVFTWO "\t1"} else{print AVECOVFTWO "\t0"}		# 2		#cov +/-25%
					if(exists($mapQ{$iso1}{$gene}) && $mapQ{$iso1}{$gene} < 55){print MAPQFTWO "\t1"} else{print MAPQFTWO "\t0"}																																	# 3		#mapQ < 55
						if($hetecount{$iso1}{$gene} > 0 && exists($samplecount{$iso1}{$gene}) && ($hetecount{$iso1}{$gene}/$samplecount{$iso1}{$gene}) > 0.5){$exgene1{$iso1}{$gene}=1}			### count excluded genes for single isolates (hete)
			if(($hetecount{$iso1}{$gene} > 0 && exists($samplecount{$iso1}{$gene}) && ($hetecount{$iso1}{$gene}/$samplecount{$iso1}{$gene}) > 0.5) || ($hetecount{$iso2}{$gene} > 0 && exists($samplecount{$iso2}{$gene}) && ($hetecount{$iso2}{$gene}/$samplecount{$iso2}{$gene}) > 0.5)){print HETEROZYG "\t".$gene}				# filter out genes with >50% heterozygous positions
			elsif((($mean_cov{$iso1}{$gene} <= ($ave_cov{$iso1}*1.25) && $mean_cov{$iso1}{$gene} >= ($ave_cov{$iso1}*0.75)) || $mapQ{$iso1}{$gene} >= 55) && (($mean_cov{$iso2}{$gene} <= ($ave_cov{$iso2}*1.25) && $mean_cov{$iso2}{$gene} >= ($ave_cov{$iso2}*0.75)) || $mapQ{$iso2}{$gene} >= 55)){	# If coverage and mapQ are in the right range in both genes
						if($min_cov{$iso1}{$gene} <= 0.9){$exgene2{$iso1}{$gene}=1}			### count excluded genes for single isolates (cov)
				if($min_cov{$iso1}{$gene} <= 0.9 || $min_cov{$iso2}{$gene} <= 0.9){print COVRATIO "\t".$gene}		# filter by coverage ratio (standard)
				else{
					@keys3iso1 = keys(%{$variants{$iso1}{$gene}});					# save pos, ref, alt and  and hete/homo of isolate 1
					@keys3iso2 = keys(%{$variants{$iso2}{$gene}});					# save pos, ref, alt and  and hete/homo of isolate 2
					%{count} = ();													# empty key3-counter
					foreach $item (@keys3iso1, @keys3iso2){$count{$item}++}			# count identical (count 2) and single (count 1) pos, ref, alt and and hete/homo as items in hash
					foreach $item (keys(%count)){
						if($count{$item}==1){
							if(exists($position{$iso1}{$gene}{$item})){
								$pos=$position{$iso1}{$gene}{$item};										# extract position of possible distance call
								if($alternative{$iso1}{$gene}{$item}=~m/[Aa]/){$alt1="A"; $alt2="a"}		# extract alt. allele of possible distance call
								if($alternative{$iso1}{$gene}{$item}=~m/[Tt]/){$alt1="T"; $alt2="t"}		# extract alt. allele of possible distance call
								if($alternative{$iso1}{$gene}{$item}=~m/[Cc]/){$alt1="C"; $alt2="c"}		# extract alt. allele of possible distance call
								if($alternative{$iso1}{$gene}{$item}=~m/[Gg]/){$alt1="G"; $alt2="g"}		# extract alt. allele of possible distance call
								$genecov=$mean_cov{$iso2}{$gene};											# extract gene coverage of counterpart
								if(exists($mpileup_two{$gene}{$pos})){@counter_fields = split(";", $mpileup_two{$gene}{$pos})}	# save counterpart entries, if exists
								else{@counter_fields=()}
							}
							elsif(exists($position{$iso2}{$gene}{$item})){
								$pos=$position{$iso2}{$gene}{$item};										# extract position of possible distance call
								if($alternative{$iso2}{$gene}{$item}=~m/[Aa]/){$alt1="A"; $alt2="a"}		# extract alt. allele of possible distance call
								if($alternative{$iso2}{$gene}{$item}=~m/[Tt]/){$alt1="T"; $alt2="t"}		# extract alt. allele of possible distance call
								if($alternative{$iso2}{$gene}{$item}=~m/[Cc]/){$alt1="C"; $alt2="c"}		# extract alt. allele of possible distance call
								if($alternative{$iso2}{$gene}{$item}=~m/[Gg]/){$alt1="G"; $alt2="g"}		# extract alt. allele of possible distance call
								$genecov=$mean_cov{$iso1}{$gene};											# extract gene coverage of counterpart
								if(exists($mpileup_one{$gene}{$pos})){@counter_fields = split(";", $mpileup_one{$gene}{$pos})}	# save counterpart entries, if exists
								else{@counter_fields=()}
							}
#							else{die "Genomic position for the distance call not found.\n"}
							if(exists($suspos{$gene}{$pos})){}			# exclude suspicious variants
							else{
								$counter=0;															# set AF counter to 0
								foreach $field (@counter_fields){									# parse counting entries
									if(($field=~m/$alt1=(\d+)/) || ($field=~m/$alt2=(\d+)/)){		# count number of occurances for the corresponding alternative allele
										$counter+=$1;
									}
								}
								if($counter/$genecov>=0.1){						# exclude positions, which have over 10% AF of the alternative allele in the mpileup file, otherwise count distance normal (else to this if)
									print MPILEUPOUT "\n".$iso1."_vs._".$iso2."\t".$gene."\t".$pos."\t".$counter/$genecov;	# save positions that are excluded from the distance due to mpileup statistics
								}
								else{
									$variantcount{$iso1}{$iso2}++;					# variant counter for each single pos, ref, alt and  and hete/homo
									$allelecount{$iso1}{$iso2}{$gene}=1;			# set allele counter for genes to 1 that have at last one variant
									$pergenevariantcount{$iso1}{$iso2}{$gene}++;	# variant counter for each single pos, ref, alt and  and hete/homo per gene
									$pergenetotalcount{$gene}++;					# variant counter for each single pos, ref, alt and  and hete/homo per gene (total)
									if(exists($completevariants{$iso1}{$gene}{$item})){print VARIANTLIST "\nVariant_in_".$iso1." (".$iso1."_vs_".$iso2.")\t".$completevariants{$iso1}{$gene}{$item}}		# list of differing variants beteen two isolates
									elsif(exists($completevariants{$iso2}{$gene}{$item})){print VARIANTLIST "\nVariant_in_".$iso2." (".$iso1."_vs_".$iso2.")\t".$completevariants{$iso2}{$gene}{$item}}		# list of differing variants beteen two isolates
								}
							}
						}
						else{
							$consensecount{$iso1}{$iso2}++;					# consensus counter for each single pos, ref, alt and  and hete/homo
						}
					}
					if(!exists $allelecount{$iso1}{$iso2}{$gene}){
						print SAMEVAR "\t".$gene;							# save genes with same variants
					}
				}
			}
			elsif(exists($variants{$iso1}{$gene}) || exists($variants{$iso2}{$gene})){
				if(exists($variants{$iso1}{$gene}) && (($mean_cov{$iso1}{$gene} > ($ave_cov{$iso1}*1.25) || $mean_cov{$iso1}{$gene} < ($ave_cov{$iso1}*0.75)) && $mapQ{$iso1}{$gene} < 55)){
					print MAPQCOV "\t".$gene;	# save genes excluded through mapQ and coverage for isolate one
					$exgene3{$iso1}{$gene}=1;   ### count excluded genes for single isolates (mapcov)	
				}
				elsif(exists($variants{$iso2}{$gene}) && (($mean_cov{$iso2}{$gene} > ($ave_cov{$iso2}*1.25) || $mean_cov{$iso2}{$gene} < ($ave_cov{$iso2}*0.75)) && $mapQ{$iso2}{$gene} < 55)){
					print MAPQCOV "\t".$gene;	# save genes excluded through mapQ and coverage for isolate two
				}
			}		 
			else{
				$consensecount{$iso1}{$iso2}++;					# consensus counter for each single pos, ref, alt and  and hete/homo
				print NOVAR "\t".$gene;							# save genes without variants 
			}
		}
    }
}

open(EXGENE,'>',"Output_$prefix\_Stats/$prefix.excluded_genes_per_isolate");				### create output for excluded genes per isolate
foreach $iso1 (sort(keys(%variants))){
	print EXGENE "\n".$iso1."\t"."High number of heterogygous calls (>50% of positions)";
	foreach $gene (@genelist){
		if(exists($exgene1{$iso1}{$gene})){print EXGENE "\t".$gene}
	}
	print EXGENE "\n".$iso1."\t"."Low coverage (>10% of the gene positions below set threshold)";
	foreach $gene (@genelist){
		if(exists($exgene2{$iso1}{$gene})){print EXGENE "\t".$gene}
	}
	print EXGENE "\n".$iso1."\t"."Abnormal coverage (-+25% from the average) and mapQ (<55)";
	foreach $gene (@genelist){
		if(exists($exgene3{$iso1}{$gene})){print EXGENE "\t".$gene}
	}
}


system("cat Output_$prefix\_Stats/$prefix.tmp2_excluded_ave_cov   | uniq > Output_$prefix\_Stats/$prefix.heat_excluded_ave_cov") and die;					system("rm Output_$prefix\_Stats/$prefix.tmp2_excluded_ave_cov");			# 3
system("cat Output_$prefix\_Stats/$prefix.tmp2_excluded_mapQ   | uniq > Output_$prefix\_Stats/$prefix.heat_excluded_mapQ") and die;							system("rm Output_$prefix\_Stats/$prefix.tmp2_excluded_mapQ");				# 4
system("cat Output_$prefix\_Stats/$prefix.tmp2_excluded_heterozyg   | uniq > Output_$prefix\_Stats/$prefix.heat_excluded_heterozyg") and die;				system("rm Output_$prefix\_Stats/$prefix.tmp2_excluded_heterozyg");			# 1
system("rm isolate_list_$prefix.tmp");
system("rm gene_list_$prefix.tmp");

close(AVECOVFTWO); 		# 2
close(MAPQFTWO);		# 3
close(HETEROZYGFTWO);	# 1

close(HETEROZYG);
close(COVRATIO);
close(NOVAR);
close(MAPQCOV);
close(VARIANTLIST);
close(MPILEUPOUT);

### print hashes into table ###
if(-e "Output_$prefix\_Tables"){}
else{
	system("mkdir Output_$prefix\_Tables");	# create Tables folder
}
sleep(1);
print "\n\n> Writing tables:\n";
sleep(1);
open(VARIANTS,'>',"Output_$prefix\_Tables/$prefix\_variant_table.txt");
open(ALLELES,'>',"Output_$prefix\_Tables/$prefix\_allele_table.txt");
open(CONSENSUS,'>',"Output_$prefix\_Tables/$prefix\_consensus_table.txt");
open(GENEVARIANTS,'>',"Output_$prefix\_Tables/$prefix\_per_gene_variant_table.txt");
open(GENETOTAL,'>',"Output_$prefix\_Tables/$prefix\_per_gene_total_table.txt");

print VARIANTS "IDs\t".join("\t", sort(keys(%variants))); # print table header line for variants
print ALLELES "IDs\t".join("\t", sort(keys(%variants)));	# print table header line for alleles
print CONSENSUS "IDs\t".join("\t", sort(keys(%variants)));	# print table header line for consensus
print GENEVARIANTS "IDs\t".join("\t", @genelist);	# print gene list header line for per gene comparison
print GENETOTAL "IDs\tCount";				# print gene list header line for per gene total variants

foreach $iso1 (sort(keys(%variants))){
	print VARIANTS "\n".$iso1;			# print table header column for variants
	print ALLELES "\n".$iso1;			# print table header column for alleles
	print CONSENSUS "\n".$iso1;			# print table header column for consensus
	foreach $iso2 (sort(keys(%variants))){
		print "\t".$iso1." vs. ".$iso2."\n";
		print VARIANTS "\t".$variantcount{$iso1}{$iso2};	# print variant counter
		print ALLELES "\t".scalar(keys(%{$allelecount{$iso1}{$iso2}}));	# print allele counter (sum up all alleles that are different)
		print CONSENSUS "\t".$consensecount{$iso1}{$iso2};	# print consensus counter

		print GENEVARIANTS "\n".$iso1."_vs._".$iso2;		# pairwise isolate comparison per line
		foreach $gene (@genelist){
			if(exists($pergenevariantcount{$iso1}{$iso2}{$gene})){print GENEVARIANTS "\t".$pergenevariantcount{$iso1}{$iso2}{$gene};	# print per gene variant counter
			}
			else{print GENEVARIANTS "\t0";
			}
		}
	}
}
foreach $gene (@genelist){
	print GENETOTAL "\n".$gene;
	if(exists($pergenetotalcount{$gene})){print GENETOTAL "\t".($pergenetotalcount{$gene})/2;		# print per gene total variants
	}
	else{print GENETOTAL "\t0";
	}
}

print ALLELES "\n";

close(VARIANTS);
close(ALLELES);
close(CONSENSUS);
close(GENEVARIANTS);
close(GENETOTAL);

print "\n\n> Start plotting stats\n\n";
system("Rscript NanoCore_plots.R Output_$prefix\_Stats/$prefix.min_cov Output_$prefix\_Stats/$prefix.heat_excluded_heterozyg Output_$prefix\_Stats/$prefix.heat_excluded_assembly_double Output_$prefix\_Stats/$prefix.heat_excluded_ave_cov Output_$prefix\_Stats/$prefix.heat_excluded_mapQ Output_$prefix\_Stats/$prefix.table Output_$prefix\_Stats/$prefix.read_length Output_$prefix\_Stats/$prefix-boxplot_min_cov.pdf Output_$prefix\_Stats/$prefix-heatmap_min_cov.pdf Output_$prefix\_Stats/$prefix-heatmap_hete.pdf Output_$prefix\_Stats/$prefix-heatmap_ass.pdf Output_$prefix\_Stats/$prefix-heatmap_cov.pdf Output_$prefix\_Stats/$prefix-heatmap_mapq.pdf Output_$prefix\_Stats/$prefix-stats_1.pdf Output_$prefix\_Stats/$prefix-stats_2.pdf");		#creating boxplot, heatmaps and stat plots
print "\n\n> Start plotting MST\n\n";
system("Rscript NanoCore_mst.R Output_$prefix\_Tables/$prefix\_allele_table.txt Output_$prefix\_Stats/$prefix-mst.pdf");		#creating mst

sleep(1);
print "\n\nOutput written to:\n";						# Write all output files to terminal
sleep(1);
print "\n\tOutput_$prefix\_Minimap/[SAMPLE_ID]/";
print "
\t\t[SAMPLE_ID].bam\t\t(Original mapping output)
\t\t[SAMPLE_ID]_new_FLAGs.bam\t\t(Correctd mapping output)
\t\t[SAMPLE_ID]_new_FLAGs.bam.bai\t\t(Index file to mapping output)
\t\t[SAMPLE_ID]_new_FLAGs_mpileup.txt\t\t(Samtools mpileup file to mapping output)
";
sleep(1);
print "\n\tOutput_$prefix\_Clair/[SAMPLE_ID]/";
print "
\t\t[SAMPLE_ID]_merge_output.vcf\t\t(Variant calling output)
\t\trun_clair3.log\t\t(Variant calling log file)
\t\tOther files\t\t(Secondary clair3 output files)
";
sleep(1);
print "\n\tOutput_$prefix\_Stats//";
print "
\t\t[SAMPLE_ID].coverage\t\t(Samtools coverage output)
\t\t[SAMPLE_ID].depth\t\t(Samtools depth output)
\t\t[SAMPLE_ID].stats\t\t(Samtools stats output)
\t\t[SAMPLE_ID].coverage_per_ref\t\t(Mean coverage and ratio of positions above coverage threshold per gene)
\t\t".$prefix.".read_length\t\t(Lengths of all reads in that run)
\t\t".$prefix.".min_cov\t\t(Matrix of ratio of positions per gene above coverage threshold of ".$threshold.")
\t\t".$prefix.".table\t\t(Table with basic read statistics)
\t\t".$prefix.".excluded_1_heterozyg\t\t(Excluded genes per isolate pair with too high ratio of heterozygous calls)
\t\t".$prefix.".excluded_2_cov_ratio\t\t(Excluded genes per isolate pair with too low ratio of positions above coverage threshold)
\t\t".$prefix.".excluded_3_same_var\t\t(Genes per isolate pair with the same variants -> no distance)
\t\t".$prefix.".excluded_4_mapq_plus_coverage\t\t(Excluded genes per isolate pair with unusual coverage AND low mapQ)
\t\t".$prefix.".excluded_5_no_var\t\t(Genes per isolate pair without variants)
\t\t".$prefix.".suspos\t\t(Suspicious positions in genes that are excluded from the complete analysis due to heterozygozity, low quality or low allele frequenzy)
\t\t".$prefix.".variantlist\t\t(List of called variants per isolate pair)
\t\t".$prefix.".heat_excluded_ave_cov\t\t(Matrix for heatmap of genes with unusual coverage) 
\t\t".$prefix.".heat_excluded_heterozyg\t\t(Matrix for heatmap of genes with unusual heterozygozity)
\t\t".$prefix.".heat_excluded_mapQ\t\t(Matrix for heatmap of genes with unusual mapping quality)
\t\t".$prefix."-boxplot_min_cov.pdf\t\t\t\t(Boxplot of ratio of positions per gene above coverage threshold of ".$threshold.")
\t\t".$prefix."-heatmap_min_cov.pdf\t\t(Heatmap of ratio of positions per gene above coverage threshold of ".$threshold.")
\t\t".$prefix."-stats_1.pdf\t\t(Plots of basic read statistics part 1/2)
\t\t".$prefix."-stats_2.pdf\t\t(Plots of basic read statistics part 2/2)
\t\t".$prefix."-mst.pdf\t\t(Minimum-spanning-tree of all samples from this run. Dashed edge means no distance)
\t\t".$prefix."-heatmap_cov.pdf\t\t(Heatmap of genes with unusual coverage)
\t\t".$prefix."-heatmap_hete.pdf\t\t(Heatmap of genes with unusual heterozygozity)
\t\t".$prefix."-heatmap_mapq.pdf\t\t(Heatmap of genes with unusual mapping quality)
";
sleep(1);	
print "\n\tOutput_$prefix\_Tables/";
print "
\t\t".$prefix."_variant_table.txt\t\t(Distance matrix on basis of variants)
\t\t".$prefix."_allele_table.txt\t\t(Distance matrix on basis of alleles)
\t\t".$prefix."_consensus_table.txt\t\t(Similarity matrix on basis of variants)
\t\t".$prefix."_per_gene_variant_table.txt\t\t(Table of variants per gene per isolate pair)
\t\t".$prefix."_per_gene_total_table.txt\t\t(Total variants per gene over all isolates)
";
sleep(1);
print "\n\n\n ~ ~ Script finished ~ ~\n\n";
sleep(1);



