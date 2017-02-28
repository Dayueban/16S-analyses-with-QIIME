#!/usr/bin/python

#PYTHON SCRIPT 
#written by: Richard Wolfe
#
#to run type: python qiime_jgi_trim_primers.py -s <database> -i <text file with ids and directories>
#         -r <results directory> -m --make_otu_table_only
#
#   if error: /usr/bin/python^M: bad interpreter: No such file or directory
#      -there is a windows endl after shebang
#      -open in vi 
#         once in vi type:
#           :set ff=unix<return>
#           :x<return>
#
#    -s Qiime database to source SILVA or GG_13_8 (original) or GG_13_5 (optional)
#       default = SILVA
#    -i tab file with 2 columns ids and JGI file path (optional)
#       this file is on JGI in sequence QC reports folder
#       and there may be 1 or 2 library.txt files that need to be cat together
#       (optional) if no id file supplied then the otu table will have the id parsed from the file name
#       Note: This does not have to be the full path and does not need to start with a /
#    -r results directory to make (required)
#       if this directory exists then will exit with error
#    -m make otu only (was used to test) (optional)
#       this will combine all the seqs_chimeras_filtered.fna  files and will add the idsif there is an id file
#       and then make the otu table. You could make all the seqs.fna files and make the Silva otu table then
#       rerun the command with -m and GG_?? database without remaking all the seqs.fna files 
#       Note -r can be anything if using -m (skips first steps of QIIME)
#
#  This script checks for the number of primers in all the files in a dir


import sys      #for exit command and maxint
import argparse #to get command line args 
                #needed to install argparse module because using python 2.6
                #and argparse comes with python 2.7
                #  sudo easy_install argparse
import os       #to run system commands
import datetime #to make timestamp
#import glob     #for * in file name
import re       #for reg expressions
#from string import maketrans #to make translate table for string timetable

import toolbox #import my toolbox with defs

def run_exit():
	#so will exit qiime environment
	#I am not sure but This may not need to be done
	#because the python script will exit that environment when it finishes

	#toolbox.run_system("deactivate_qiime")
	print "Script finished..."

	sys.exit(0)

	return


def trim_reads(forward_file, reverse_file, out_forward, out_reverse):
	#these are fastq files	with the same number of lines
	f = open(forward_file, "rU")
	r = open(reverse_file, "rU")
	out_f = open(out_forward, "w")
	out_r = open(out_reverse, "w")
	
	f_lines = 0
	r_lines = 0
	out_f_lines = 0
	out_r_lines = 0
	r_seq = 0
	f_seq = 0

	#bp[0] = blank to bp[9] = 9-. This will be any char in reg expression
	bp = ["",".","..","...","....",".....","......",".......","........","........."]

	f_primer = 'GTGCCAGC.GCCGCGGTAA' #forward primer 19bp
	r_primer = 'GGACTAC..GGGT.TCTAAT' #rev primer 20 bp

	f_primer_found = 0
	r_primer_found = 0
	
	found_in_both = 0
	
	total_prefix = 9  #if 0 on front then 9 on rear
			  #Ex  0bp|forward primer|  read  |reverse primer rc| 9bp
			  #    1bp|forward primer|  read  |reverse primer rc| 8bp
			  # the total number of bp removed will be 9 + len_forward + len_reverse primer


	f1 = f.readline() 
	while(f1):
		f1 = f1.rstrip()
		f2 = f.readline().rstrip()  #read and remove endline
		f3 = f.readline().rstrip()  #read and remove endline
		f4 = f.readline().rstrip()  #read and remove endline

		r1 = r.readline().rstrip()  #read and remove endline
		r2 = r.readline().rstrip()  #read and remove endline
		r3 = r.readline().rstrip()  #read and remove endline
		r4 = r.readline().rstrip()  #read and remove endline

		
		
		f_lines += 4
		r_lines += 4
		f_seq += 1
		r_seq += 1

		
		#look for primers at start of seq
		forward_item = "NONE"
		reverse_item = "NONE"
		for item in bp:
			if re.search("^" + item + f_primer,f2):
				f_primer_found += 1
				forward_item = item  #save the found prefix
				break

		if forward_item != "NONE":  #if a match was found in the forward seq then check reverse
 			for item in bp:	
				if re.search("^" + item + r_primer,r2):
					r_primer_found += 1
					reverse_item = item
					break

		if forward_item != "NONE" and reverse_item != "NONE":
			#If the pattern is in both sequences then we want to trim
			#and write it to to output files
			found_in_both += 1
	
			begin_trim = len(forward_item) + len(f_primer)
			end_trim = total_prefix - len(forward_item) + len(r_primer)
			end_trim = 0 - end_trim #make negative
			
			f2 = f2[begin_trim:end_trim] 
			f4 = f4[begin_trim:end_trim] 
				
			begin_trim = len(reverse_item) + len(r_primer)
			end_trim = total_prefix - len(reverse_item) + len(f_primer)
			end_trim = 0 - end_trim #make negative
			
			r2 = r2[begin_trim:end_trim] 
			r4 = r4[begin_trim:end_trim] 
			
			out_f.write(f1 + "\n")
			out_f.write(f2 + "\n")
			out_f.write(f3 + "\n")
			out_f.write(f4 + "\n")
			out_r.write(r1 + "\n")
			out_r.write(r2 + "\n")
			out_r.write(r3 + "\n")
			out_r.write(r4 + "\n")
			

			out_f_lines += 4
			out_r_lines += 4
			
		
		f1 = f.readline()

	f.close()
	r.close()
	out_f.close()
	out_r.close()

	print "Lines in forward file = ",f_lines
	print "Lines in reverse file = ",r_lines
	print ""
	print "Sequences in forward file = ", f_seq
	print "Sequences in reverse file = ", r_seq
	print ""
	print "Forward primers found = ", f_primer_found
	print "Reverse primers found = ", r_primer_found
	print ""
	print "Number of sequences with primer in both = ", found_in_both
	print ""
	print "Number of lines wrote to forward output file = ",out_f_lines
	print "Number of lines wrote to reverse output file = ",out_r_lines

	
	

	return


def make_otu_table():

	#combine all the log files
	#note that there should only be the folders in the directory because we created the directory
	cmd = "cat */STEP1_OUT/STEP2_OUT/split_library_log.txt > ALL_split_library_log.txt"
	print "Combining all the log files"
	toolbox.run_system(cmd)

	#cat all the seqs_chimeras_filtered.fna files
	cmd = "cat */STEP1_OUT/STEP2_OUT/seqs_chimeras_filtered.fna  > combined_seqs_chimeras_filtered.fna"
	print "Combining all the seqs_chimeras_filtered.fna files"
	toolbox.run_system(cmd)

	"""
	if args.ids: #if id file was supplied then change the prefixes
	 	#note args is a global variable - we can read but not cange it ?
	 	print "args.ids.name = " + args.ids.name

		#make a list of ids and files
		ids = []
		files = []
		id_lines = 0

		line = args.ids.readline()

		while line:
			id_lines += 1
			line = line.rstrip() #remove endline
			cols = line.split() #splits on whitespace
			cols[0] = cols[0].replace("_", "") #remove all underscores
			if cols[0] in ids: #checks for duplicate ids
				print "Error ... id is already in id list = " + cols[0]
				sys.exit(1)
			else:
				ids.append(cols[0]) #id in first col

			#/global/dna/dm_archive/sdm/illumina/01/00/83/10083.1.147588.TTGTCGCACAA.fastq.gz
			#extract the id from the file name
			file_id = cols[1].split(".")[-3]
			if file_id in files:
				print "Error ... file_id is already in list = " + file_id
				sys.exit(1)
			else:
			 files.append(file_id)

			print cols[0] + " " + file_id

			line = args.ids.readline()
		
		args.ids.close
		print "Lines read from id file = ", id_lines

		#read the seqs.fna file and change the ids
		f = open ("combined_seqs_chimeras_filtered.fna", "r")
		out_file = open("NEW_combined_seqs_chimeras_filtered.fna", "w")

		line = f.readline()
		seqs_lines = 0
		seqs_sequences = 0
		while line:
			header,sequence,line,s_lines = toolbox.read_fasta(line,f)
			#Note the header and sequence still have endlines
			seqs_sequences += 1
			seqs_lines += (1 + s_lines) #the header line and the seq lines

			#>ACATATACGCG_0 MISEQ0.....
			header = header[1:] #remove >
			cols = header.split("_",1) #split on 1st _

			index = files.index(cols[0])
			header = ">" + ids[index] + "_" + cols[1]  #the id_0 MISEQ0.....

			out_file.write(header)
			out_file.write(sequence)

		f.close()
		out_file.close()
		print "Sequences in combined_seqs_chimeras_filtered.fna = ", seqs_sequences
		print "Lines in combined_seqs_chimeras_filtered.fna = ",seqs_lines 

		#swap the files
		cmd = "mv NEW_combined_seqs_chimeras_filtered.fna combined_seqs_chimeras_filtered.fna"
		toolbox.run_system(cmd)
	"""
		
	
	#make otu table
	pwd = os.getcwd()
	step3_folder = pwd + "/" + args.source + "_STEP3_OUT"
	print "Running pick_open_reference_otus.py"
	cmd = "pick_open_reference_otus.py -i  " + pwd + "/combined_seqs_chimeras_filtered.fna -r /home2/Database/Silva/rep_set/97_Silva_111_rep_set.fasta -o " + step3_folder + " -f -a -O 60"
	#print cmd
	toolbox.run_system(qiime_source + " && " + cmd)

	print "Making otu from biom file"
	cmd = "summarize_taxa.py -i " + step3_folder + "/otu_table_mc2_w_tax.biom -o " + step3_folder + "/taxonomy_summaries/  -L 2,3,4,5,6" 
	toolbox.run_system(qiime_source + " && " + cmd)

	#added per lindsey
 
	cmd = "python /ORG-Data/scripts/wrapper_filter_otus_from_otu_table.py -i " + step3_folder + "/otu_table_mc2_w_tax.biom -o " + step3_folder + "/percent_filtered_otu_table_mc2_w_tax.biom -n 10 -p 25"
 	toolbox.run_system(qiime_source + " && " + cmd)

	#biom convert -i otu_table_mc2_w_tax.biom -o otu_table_mc2_w_tax.txt --to-tsv --header-key taxonomy
	cmd = "biom convert -i " + step3_folder + "/percent_filtered_otu_table_mc2_w_tax.biom -o " + step3_folder + "/percent_filtered_otu_table_mc2_w_tax.txt --to-tsv --header-key taxonomy"
 	toolbox.run_system(qiime_source + " && " + cmd)

	cmd = "python /ORG-Data/scripts/calculate_relative_abundance.py -i " + step3_folder + "/percent_filtered_otu_table_mc2_w_tax.txt -o " + step3_folder + "/percent_calculate_relative_abundance_output.txt"
 	toolbox.run_system(qiime_source + " && " + cmd)

	return

print "Script started ..."


#create an argument parser object
#description will be printed when help is used
parser = argparse.ArgumentParser(description='A script to make otu from JGI files')

#add the available arguments -h and --help are added by default
#if the input file does not exist then program will exit
#if output file does not exit it will be created
# args.input is the input file Note: cant write to this file because read only
# args.output is the output file
# args.m is the minimum seq length
#parser.add_argument('-d', '--dir', help='Directory with interleaved unzipped files',required=True)
parser.add_argument('-s', '--source', help='Qiime database to source SILVA or GG_13_8 (original) or GG_13_5',default="SILVA")
parser.add_argument('-i', '--ids', type=argparse.FileType('rU'), help='JGI ids file: library.txt with file_name and id', required=True)
parser.add_argument('-r', '--results_directory',  help='Results directory to place results',required=True)
#parser.add_argument('-m', '--make_otu',  help='Only make otu table was used to test')


#get the args
args = parser.parse_args()

#additional argument tests




qiime_source = "source /opt/Qiime_1.9/activate_qiime_1.9 "

if args.source == "SILVA": #the default
	qiime_source = qiime_source + "SILVA"
elif args.source == "GG_13_8": #the default
	qiime_source = qiime_source + "GG_13_8"
elif args.source == "GG_13_5": #the default
	qiime_source = qiime_source + "GG_13_5"
else:
    	print "Error ... -s must be SILVA or GG_13_8 or GG_13_5"
	sys.exit(1)


#if args.ids:
#	print "ids file supplied: " + args.ids.name
#else:
#	print "ids file was not supplied"


#if args.make_otu:
#	print "Starting at make otu table only"
#	make_otu_table()
#	sys.exit(0)

#make sure directory is a path and starts with a /
#if not args.dir.startswith("/"):
#	print "Error ... -d needs to be a full path and start with a /"
#	sys.exit(1)

#make sure ids file starts with a / and is a full path
#This does not need to be path because it is opened before we change directories
#if args.ids and not args.ids.name.startswith("/"):
#	print "Error ... -i needs to be a full path and start with a /"
#	sys.exit(1)


#make results directory
if not os.path.exists(args.results_directory): #if dir not already
		os.makedirs(args.results_directory)
		os.chdir(args.results_directory) #cd into the new dir
		
else:
	print "Error ... This directory already exists"
	sys.exit(1)

print qiime_source
#sys.exit(0)

#if not args.dir.startswith("/"):
#	#if args.dir is not a full path we need to add ../ to it because we change directories
#	args.dir = "../" + args.dir


#get all the files in a directory

file_list = []  #os.listdir(args.dir)  #make a list of all the files and directories in the directory except . and ..

########################################
#we need to make a list of the sample names
#if the jgi file is supplied we will get from there
#if not then we will make from the file id

sample_ids = []

if args.ids:
	print "ids file supplied: " + args.ids.name

	#make all sample ids = XXXXXXXXXX "10 X"
	#for item in file_list:
	#	sample_ids.append("XXXXXXXXXX")

	line = 	args.ids.readline()
	file_lines = 0

	while line:
		file_lines += 1
		#col1 = id and col2 = file
		#Ex:  03_P7A1D1       /global/dna/dm_archive/sdm/illumina/01/00/83/10083.1.147588.TTGTCGCACAA.fastq.gz
		cols = line.split()
		file_name = cols[0]   #cols[1].split("/")[-1]
		#file_name = file_name.replace(".gz","") #remove .gz becasue they are unzipped
		id = cols[1]   #cols[0]
		id = id.replace("_","") #remove all _
		
		#check if this id already exists:
		if id in sample_ids:
			print "Error .. id already exists " + id
			sys.exit(1)
		if file_name in file_list:
			print "Error .. File already used "  + file_name
			sys.exit(1)
		
		#make sure file exists
		if os.path.isfile(file_name):
			print "file exists" + file_name
		else:
			print "Error .. file does not exist " + file_name
			sys.exit(1)

		#make sure file is unzipped
		if file_name.endswith(".gz"):
			print "Error .. file must be unzipped and not end in .gz"
			sys.exit(1)

		sample_ids.append(id)
		file_list.append(file_name)

		line = 	args.ids.readline()

	#check if all files have and id
	#if "XXXXXXXXXX" in sample_ids:
	#	print "Error ... A file is missing a sample_id"
	#	sys.exit(1)

	print "lines read from ids file = ", file_lines



	
print "Length of ids = ",len(sample_ids)
#sys.exit(0)
		



########################################
#need to set up the qiime environment
#This does not work because the environment is only set for this instruction
#  and further instructions are not affected
#We would need to run a wrapperscript that sets the env variables and then runs this python
#   andt then runs the deactivate_quiime after 
# OR maybe run 2 instructions in the command ( "comadd1 && comand2 && command3 ..")
# set shows the shell variables
#toolbox.run_system("source /opt/Qiime_1.9/activate_qiime_1.9 SILVA && cmd")

for item in file_list:
	print item

	#get the sample id and make the directory 
	index = file_list.index(item)
	
	id = sample_ids[index]
	
	if not os.path.exists(sample_ids[index]): #if dir not already
		os.makedirs(sample_ids[index])
	else:
		#should never happen because we make this dir
		print "Error .... Directory already exists"
		sys.exit(1)

	#seperate the forward and reverse read and put into new dir
	#file_name = item  #args.dir + "/" + item
	toolbox.deinterleave_fastq_reads(item,sample_ids[index] + "/reads-1.fq",sample_ids[index] + "/reads-2.fq")

	#trim the primers from the reads

	trim_reads(sample_ids[index] + "/reads-1.fq", sample_ids[index] + "/reads-2.fq", sample_ids[index] + "/trimmed_reads-1.fq",sample_ids[index] + "/trimmed_reads-2.fq" ) #pass the forward and reverse file

		
	#join paired ends
	# join_paired_ends.py -f no_primers/reads.fastq -r no_reverse_primers/reads.fastq -o NO_BC_STEP1_OUT/
	cmd = "join_paired_ends.py -f " + sample_ids[index] + "/trimmed_reads-1.fq -r " + sample_ids[index] + "/trimmed_reads-2.fq -o " + sample_ids[index] + "/STEP1_OUT/"
	toolbox.run_system(qiime_source + " && " + cmd)
	
	
	#make a mapping file with the id
	#  10279.1.153921.CATCATGAGGC.fastq id = 10279.1.153921.CATCATGAGGC.fastq
	#id = item.split(".")[3]
	f = open(sample_ids[index] + "/" + sample_ids[index] + "_mapping.txt","w")
	f.write("#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\n")
	f.write(id + "\t\t\t" + sample_ids[index] + "\n")
	f.close()

	#run split libraries
	#split_libraries_fastq.py -i fastqjoin.join.fastq -o STEP2_OUT/ -m file1_mapping.txt -q 19 --store_demultiplexed_fastq --barcode_type not-barcoded --sample_ids AACAGGTTCGC
	cmd = "split_libraries_fastq.py -i " + sample_ids[index] + "/STEP1_OUT/fastqjoin.join.fastq -o " + sample_ids[index] + "/STEP1_OUT/STEP2_OUT/ -m " + sample_ids[index] + "/" + sample_ids[index] + "_mapping.txt -q 19 --store_demultiplexed_fastq --barcode_type not-barcoded --sample_ids " + sample_ids[index]
	toolbox.run_system(qiime_source + " && " + cmd)

	#check how many sequences in seqs.fna file
	lines_seqs = toolbox.count_seqs_fasta(sample_ids[index] + "/STEP1_OUT/STEP2_OUT/seqs.fna")
	print "Sequences in seqs.fna file = ", lines_seqs
	if lines_seqs == 0:
		continue	


	#chimera check
	#identify_chimeric_seqs.py -i seqs.fna -m usearch61 -r /home2/Database/RDP_Gold/rdp_gold.fa -o usearch61_chimera/
	cmd = "identify_chimeric_seqs.py -i " + sample_ids[index] + "/STEP1_OUT/STEP2_OUT/seqs.fna -m usearch61 -r /home2/Database/RDP_Gold/rdp_gold.fa -o " + sample_ids[index] + "/STEP1_OUT/STEP2_OUT/usearch61_chimera/"
	toolbox.run_system(qiime_source + " && " + cmd)
	
	#remove chimeras
	#filter_fasta.py -f seqs.fna -o seqs_chimeras_filtered.fna -s usearch61_chimera/chimeras.txt -n
	cmd = "filter_fasta.py -f " + sample_ids[index] + "/STEP1_OUT/STEP2_OUT/seqs.fna -o " + sample_ids[index] + "/STEP1_OUT/STEP2_OUT/seqs_chimeras_filtered.fna -s " + sample_ids[index] + "/STEP1_OUT/STEP2_OUT/usearch61_chimera/chimeras.txt -n"
	toolbox.run_system(qiime_source + " && " + cmd)
	




make_otu_table()
	
	

	
	
        

print ""
print "Number of files in directory = ",len(file_list)


run_exit()
