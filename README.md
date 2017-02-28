# 16S-analyses-with-QIIME

QIIME Protocol and Explanation
February 20, 2017
Lindsey Solden and Mikayla Borton

Purpose: Generate Operational Taxonomic Unit (OTU) table from 16S reads 

We receive 16S reads from different sources and each require a different script. The general logic behind the protocol is the same, but different inputs require different formatting.

ARGONNE

Script /ORG-Data/scripts/QIIME_PIPELINE.sh
Files
•	Zipped barcodes file
o	barcodes are unique to each sample
•	Zipped forward reads file 
o	contains barcode, the sequence, quality score (phred)
•	Zipped reverse reads file 
o	contains barcode, the sequence, quality score (phred)
•	Mapping file 
o	Provide the full path to the mapping fie
o	Mapping file can not contain any special characters including underscores (_) and dashes (-)
o	You need separate mapping files for all projects on the sequencing run. Each OTU table you build needs a separate mapping file
Use
bash /ORG-Data/scripts/bin/Phylogeny_Protpipe/QIIME_PIPELINE.sh /home2/projects/170215_Wrighton_Daly_16S_fastqs/Undetermined_S1_L001_I1_001.fastq.gz /home2/projects/170215_Wrighton_Daly_16S_fastqs/Undetermined_S1_L001_R1_001.fastq.gz /home2/projects/170215_Wrighton_Daly_16S_fastqs/Undetermined_S1_L001_R2_001.fastq.gz /home2/projects/170215_Wrighton_Daly_16S_fastqs/170215_Wrighton_Daly_MF_16SRV_NJK_02092017_AllTimeSeriesMice.txt solden.2@osu.edu 10 p 0

JGI

Script /ORG-Data/scripts/qiime_jgi_trim_primers.py 
Files
•	One fastq file per sample you want in your OTU table
•	Text file 
o	Identifies the full path location of above reads, sample names, number of reads and project. 
o	No column headers
o	Columns must be in order listed. 
o	MAKE SURE SAMPLE names do not have underscores.
•	Results directory
Use 
python /ORG-Data/scripts/qiime_jgi_trim_primers.py -s SILVA -i hfts.txt -r HFTS_LMS_test/


EMIRGE
Script
Files
Use
•	Files needed to build QIIME files from EMIRGE files:
o	If you are clustering EMIRGE data with other shorter 16S sequences (like V4), download each renamed EMIRGE server file from the server and align to 16S V4 region using Geneious. Trim each sequence to the V4 region or other desired length and save each file. Upload each new trimmed file. 
o	If you are building an OTU table from full length EMIRGE, you do not need to trim the files.
o	NOTE: files can’t contain underscores in the file name or fasta name
•	For each renamed EMIRGE file run:
python /home/lsolden/Scripts/make_qiime_file_from_emirge.py -i S-1-Day1.fasta -p S1Day1 -o S1Day1_seqs.fna
o	This script will calculate the number of sequences from the norm prior, based on 1000 sequences and build a seqs.fna file
•	Concatenate each seqs.fna file and run through EMIRGE starting at Pick OTUs step.


Here is a step-by-step breakdown of the QIIME_PIPELINE.sh script, which is the most basic use of QIIME. For JGI and EMIRGE reads, after formatting, the respective scripts will begin in the QIIME_PIPELINE at pick_open_reference_OTU.py step.

To source into QIIME, use the following command:
source /opt/Qiime_1.9/activate_qiime_1.9 SILVA

To run the QIIME pipeline, use the following command:

/ORG-Data/scripts/bin/Phylogeny_Protpipe/QIIME_PIPELINE.sh 

$1 = the unzipped barcode file or START_STEP_2
	$2 = the unzipped forward reads file or full path to fastqjoin.join.fastq
   	$3 = the unzipped reverse reads file or path to fastqjoin.join_barcodes.fastq
	$4 = the full path to the mapping file 
   	$5 = email or NO_EMAIL
	$6 = number of OTUs that must be present per OTU
  	$7 = p or s  
•	p = percent of samples an OTU must be observed in for that OTU to be retained, this can be zero, but it can be an integer
•	s = the minimum number of samples an OTU must be observed in for that OTU to be retained, this can be zero, but it can be an integer
   	$8 = integer value for $7  
•	Ex: 25 and not .25

Example:
/ORG-Data/scripts/bin/Phylogeny_Protpipe/QIIME_PIPELINE.sh /ORG-Data/150722_Wrighton/Undetermined_S0_L001_I1_001.fastq.gz /ORG-Data/150722_Wrighton/Undetermined_S0_L001_R1_001.fastq.gz /ORG-Data/150722_Wrighton/Undetermined_S0_L001_R2_001.fastq.gz
/home/projects/Ahmer/Feces_Cecum_combined_08282016/Wrighton_150772_mappingfile.txt solden.2@osu.edu p 25

Protocol and Explanation of each QIIME step, including scripts:
1.	Join pair end reads 
join_paired_ends.py -f forward_reads.fastq -r reverse_reads.fastq -b barcode.fastq -o STEP1_OUT/
•	All files output into STEP1_OUT folder
•	Files output:
i.	fastqjoin.join.fastq (all joined reads)
ii.	fastqjoin.un1.fastq (forward unjoined reads)
iii.	fastqjoin.un2.fastq (reverse unjoined reads)
iv.	fastqjoin.join_barcodes.fastq (all joined barcodes)
•	QUALITY CHECK- At this step, reads that do not join are thrown out.
2.	Split Library
split_libraries_fastq.py -i fastqjoin.join.fastq -b fastqjoin.join_barcodes.fastq --rev_comp_mapping_barcodes -o STEP2_OUT/ -m mapping_file.txt  -q 19 --store_demultiplexed_fastq
•	All files output into STEP2_OUT/
•	This script performs demultiplexing of Fastq data. In other words, it sorts the data according to your mapping file by sample. 
•	In this step a fasta file sorted by sample is built with reads that have a phred score of 20 or greater.  The barcode and primer are removed from sequence here.
•	We choose a phred score of 20 or better, meaning 1 in 100 chance that a base is called wrong for each base.
•	Files output:
i.	split_library_log.txt (log file that contains the number of reads in each sample, number of input reads, median sequence length, etc.)
ii.	seqs.fna (fasta file sorted by sample built with reads that have a phred score of 20 or greater)
iii.	histograms.txt (quality file output)
iv.	seqs.fastq (fastq file quality trimmed by phred score of 20 or greater)
•	QUALITY CHECK: Reads with a phred quality score of 19 or below are thrown out.
3.	Chimera checking
identify_chimeric_seqs.py -i seqs.fna -m usearch61 –r /home2/Database/RDP_Gold/rdp_gold.fa -o usearch61_chimera/

filter_fasta.py -f seqs.fna -o seqs_chimeras_filtered.fna -s usearch61_chimera/chimeras.txt –n
•	All files output into STEP2_OUT/
•	The first command identifies the chimeric sequences using RDP Gold database.  The second command removes those chimeric sequences from the seqs.fna file generated in the split library step.
•	Files output:
i.	seqs_chimera_filtered.fna (fasta file sorted by sample built with reads that have a phred score of 20 or greater containing no chimeric sequences)
•	If you are going to combine data from different sequencing runs, then concatenate these files from each dataset and proceed with qiime.
ii.	usearch61_chimera file which contains several files:
•	chimeras.txt (contains all chimeric sequences removed from seqs.fna file)
•	identify_chimeric_seqs.log (log file denoting which sequences have chimeras)
4.	Pick OTUs
pick_open_reference_otus.py -i  seqs_chimeras_filtered.fna -r /home2/Database/Silva/rep_set/97_Silva_111_rep_set.fasta -o /STEP3_OUT -f -a -O 40
•	All files output into STEP3_OUT
•	While there are several methods used for clustering OTUs, in the Wrighton lab we use open reference picking. 
•	Open reference picking steps:
o	Step 1: The seqs_chimeras_filtered.fna is queried against a database (SILVA) and OTUs are picked in a closed reference style. All initial files from step one are output into step1_otus file within STEP3_OUT
•	step1_rep_set.fna (this contains the representative sequences of picked Step1 OTUs)
•	failures.fasta (this contains the sequences that did not cluster)
o	Step 2: Sequences that failed to cluster are clustered de novo
o	Step 3: All representative sequences are clustered to database
o	Step 4: Build OTU map and final fasta containing the representative sequences
•	final_otu_map_mc2.txt 
•	otu_table_mc2.biom
•	rep_set.fna
5.	Build OTU table
python /ORG-Data/scripts/wrapper_filter_otus_from_otu_table.py -i otu_table_mc2_w_tax.biom -n 10 -p 10 -o 10_percent_filtered_otu_table_mc2_w_tax.biom

biom convert -i 10_percent_filtered_otu_table_mc2_w_tax.biom -o 10_percent_filtered_otu_table_mc2_w_tax.txt --to-tsv --header-key taxonomy

python /ORG-Data/scripts/calculate_relative_abundance.py -i 10_percent_filtered_otu_table_mc2_w_tax.txt -o 10_percent_filtered_otu_table_rel_abundance.txt
•	All files output into STEP3_OUT
•	The first command filters OTUs that do not have enough reads or are not represented in enough samples. 
i.	–n flag number of reads that need to be in each OTU (we always use 10)
ii.	–p OR –s flag are the percentage of samples (p) or number of samples (s) that an OTU must be in 
•	The second command converts the filtered biom file into a text file that includes taxonomy
•	The third command generates an OTU table from the text file 
6.	OTHER STEPS:
•	If you want tables of summarized by other taxonomy levels, use this command:
summarize_taxa.py -i /home/projects/shale/Rdaly/02082017/STEP3_OUT/otu_table_mc2_w_tax.biom -o /home/projects/shale/Rdaly/02082017/STEP3_OUT/taxonomy_summaries/ -L 2,3,4,5,6
i.	Input: biom file in STEP3_OUT folder
ii.	–L is the level of taxonmy generated (L6 is genus level, L5 is family level, L4 is order level, etc.)
