#!/usr/bin/bash
#
#How to run: ./NanoCore_v1.0.5.sh -s SampleSheet_NanoCore_env_test.txt -r cgMLST_files/Pseudomonas_aeruginosa_cgMLST_ref-seqs.fasta -p NanoCore_Test_6 -t 8 -m 20 -b samtools
#
start_time=$(date +%s)	# save starting time

### generate --help and --verion text
if [ "$1" == "-h" ] || [ "$1" == "--help" ];then
    echo -e "\n\nThis is the help section of NanoCore. The following lines provide explnaitions for the different command line paraameters:\n"
    echo "-v/--version          \t\tPrints the current version number and exits."
    echo -e "-h/--help          \t\tPrints this help section and exits.\n"
    echo "-s [PATH/TO/FILE]     \t\tA tab-separated file containing one line per sample with the isolate ID, the tag 'Nanopore' or 'Illummina', the clair model and the paths to the sequencing data."
    echo "-r [PATH/TO/FILE]     \t\tThe core genome reference file for a certain species."
    echo "-p [PREFIX]           \t\tThe chosen prefix/name for the current NanoCore run."
    echo "-t [THREADS]          \t\tThe number of threads used for components of the pipeline that support multithreading."
    echo "-m [THRESHOLD]        \t\tThe minimum coverage threshold desired for the analysis (default = 20)."
    echo "-b [PATH/TO/SAMTOOLS] \t\tThe samtools executable (default = samtools)."
    echo -e "\nFor further explanation of the parameters or information about expected output files, please refer to the NanoCore GitHub Homepage (https://github.com/SebastianMeyer1989/NanoCore)\n\n."
	exit 0
fi
if [ "$1" == "-v" ] || [ "$1" == "--version" ];then
    echo -e "\n\nCurrent version of NanoCore: 1.0.5\n\n"
    exit 0
fi	

### initialize pathways, parameters and conda
script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )		# save scrips path
eval "$(conda shell.bash hook)"			# initialize conda commands
while getopts s:r:p:t:mb flag			# get parameters from command line
do
    case "${flag}" in
        s) samples=${OPTARG};;
        r) ref=${OPTARG};;
        p) prefix=${OPTARG};;
		t) threads=${OPTARG};;
		m) threshold=${OPTARG};;
		b) samtools=${OPTARG};;
    esac
done

### run script part 1
conda activate NanoCore_1				# activate conda environment for calculations
perl $script_dir/NanoCore_subroutine_1.pl --sample_list $samples --reference $ref --prefix $prefix --threads $threads --threshold $threshold--samtools $samtools
if [ $? -ne 0 ]; then
    exit 1
fi
conda deactivate

### run script part 2
conda activate NanoCore_2				# activate conda environment for plots
perl $script_dir/NanoCore_subroutine_2.pl --prefix $prefix --threshold $threshold
if [ $? -ne 0 ]; then
    exit 1
fi
conda deactivate

### calculate time
end_time=$(date +%s)						# save ending time
elapsed_time=$(( end_time - start_time ))	# calculate time difference

days=$(( elapsed_time / 86400 ))
hours=$(( (elapsed_time % 86400) / 3600 ))
minutes=$(( (elapsed_time % 3600) / 60 ))
seconds=$(( elapsed_time % 60 ))

echo -e "\n\nElapsed time: ${days} days, ${hours} hours, ${minutes} minutes, ${seconds} seconds"
echo -e "\n\nElapsed time: ${days} days, ${hours} hours, ${minutes} minutes, ${seconds} seconds" >> $prefix.log
echo -e "\n\n ~ ~ Script finished ~ ~\n\n\n";