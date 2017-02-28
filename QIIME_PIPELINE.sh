#!/bin/bash

#  QIIME_PIPELINE.sh
#
#  /ORG-Data/scripts/bin/Phylogeny_Protpipe/QIIME_PIPELINE.sh
#
# script by Richard Wolfe
#
#   $1 = the zipped barcode file or START_STEP_2
#       Ex: /ORG-Data/alaska_moose/seasonal_nutrition/16S_all_samples/RF/Undetermined_S0_L001_I1_001.fastq.gz  
#   $2 = the zipped forward reads file or full path to fastqjoin.join.fastq
#       Ex: /ORG-Data/alaska_moose/seasonal_nutrition/16S_all_samples/RF/Undetermined_S0_L001_R1_001.fastq.gz 
#   $3 = the zipped reverse reads file or path to fastqjoin.join_barcodes.fastq
#       Ex: /ORG-Data/alaska_moose/seasonal_nutrition/16S_all_samples/RF/Undetermined_S0_L001_R2_001.fastq.gz 
#   $4 = the mapping file  -- NEEDS TO BE FULL PATH -- DO NOT USE ../filename
#   $5 = email or NO_EMAIL
#   $6 = number of sequences required per OTU. Usually 10.
#   $7 = p or s   p = percent of samples  an otu must be observed in for that otu to be retained (can be 0) (must be int) usually 25 
#                 s = the minimum number of samples an OTU must be observed in for that otu to be retained  (can be 0) (must be int)
#   $8 = integer value for $7  Ex: 25 and not .25


#make a unique directory name
prefix="qiime_"
suffix=$(date +%s)   #get seconds from 1/1/1970
dir_name=$prefix$suffix
mkdir $dir_name
cd $dir_name

echo "Made a working directory $dir_name"

start_dir=$(pwd)
logfile=$start_dir/qiime_script.log

#make a variable because I dont think $5 can be seen inside the function
email_add=$5

#start a logfile
echo "Qiime script started" > $logfile
echo "Qiime script started and a log file was made at $logfile"
echo "logfile location = $logfile" >> $logfile
echo " " >> $logfile

echo "Arg 1 barcode file= $1" >> $logfile
echo "Arg 2 forward reads file = $2" >> $logfile
echo "Arg 3 reverse reads file = $3" >> $logfile
echo "Arg 4 mapping file = $4"  >> $logfile
echo "Arg 5 email address = $email_add" >> $logfile
echo "Arg 6 = $6" >> $logfile
echo "Arg 7 = $7" >> $logfile
echo "Arg 8 = $8" >> $logfile

email_exit()
{
#this is a function that will email the log and then exit
if [ "$email_add" == "NO_EMAIL" ]
then
   echo "Script finished and No email was sent" >> $logfile
   echo "Script finished and No email was sent"
else
   #send email to user logfile will be written to the email
   echo "Script finished and sending email to $email_add" >> $logfile
   echo "Script finished and sending email to $email_add"
   mail -s "Qiime Job Finished" $email_add < $logfile
fi

deactivate_qiime
echo "Qiime 1.9 environment was deactivated" >> $logfile


exit 0

}


#check if 5th variable is empty
if [ -z $8 ]
then
  echo "You did not provide 8 attributes" >> $logfile
  echo "You did not provide 8 attributes"
  #email_exit
  #dont send email because no email address
  exit 1
fi

source /opt/Qiime_1.9/activate_qiime_1.9 SILVA
echo "Qiime 1.9 environment was activated" >> $logfile

#make sure $7 is p or s
if [ "$7" != "p" ]
then
   if [ "$7" != "s" ]
   then
       echo "Error .. var 7 must be p or s"
       email_exit
   fi
fi


if [ "$1" == "START_STEP_2" ]
then
   echo "Starting at step 2" >> $logfile
elif [ ! -f $1 ]
then
   echo "Error ... file $1 not found...Exiting script" >> $logfile
   email_exit
else
   echo "$1 file exists" >> $logfile
   if [[ $1 != /* ]]
   then 
      echo "Error $1 not full path ... Exiting script" >> $logfile
      email_exit
   fi
fi

#make sure the files exist
for i in $2 $3 $4
do
  if [ ! -f $i ]
  then
    echo "Error .. File $i not found....Exiting script" >> $logfile
    #echo "Arg 5 = $5" >> $logfile
    email_exit
    #exit 1
  fi
done

#make sure all files are full path names
#checks to make sure they start with a /
for i in $2 $3 $4
do
  if [[ $i != /* ]]
  then
    echo "Error .. $i is not full path...Exiting script" >> $logfile
    email_exit
  fi
done

#convert mapping file from mac to unix if it came from a mac
echo "Converting mapping file to unix" >> $logfile
mac2unix $4


#make sure mapping file is okay
echo "Checking to make sure mapping file is correct" >> $logfile

map_ok=$(validate_mapping_file.py -m $4 -o validate_mapping_file_output)

echo "map_ok = $map_ok" >> $logfile
if [ "$map_ok" == "No errors or warnings were found in mapping file." ]
then
   echo "Mapping file okay" >> $logfile
else
   echo "Mapping file not okay...Exiting" >> $logfile
   email_exit
fi
   

if [ "$1" != "START_STEP_2" ]
then

   echo " " >> $logfile
   echo "Unzipping the 3 data files" >> $logfile

   gunzip -c $1 > barcode.fastq

   #test exit of previous command $? is exit status of last command
   # if 0 then normal exit
   if [ $? -ne 0 ]
   then
     echo "File $1 unable to be unzipped" >> $logfile
     rm -f barcode.fastq
     email_exit
     #exit 1
   fi

   gunzip -c $2 > forward_reads.fastq

   #test exit of previous command $? is exit status of last command
   # if 0 then normal exit
   if [ $? -ne 0 ]
   then
     echo "File $2 unable to be unzipped" >> $logfile
     rm -f forward_reads.fastq
     email_exit
     #exit 1
   fi

   gunzip -c $3 > reverse_reads.fastq

   #test exit of previous command $? is exit status of last command
   # if 0 then normal exit
   if [ $? -ne 0 ]
   then
     echo "File $3 unable to be unzipped" >> $logfile
     rm -f reverse_reads.fastq
     email_exit
     #exit 1
   fi

   #show the file sizes
   barcode_lines=$( cat barcode.fastq | wc -l )
   forward_lines=$( cat forward_reads.fastq | wc -l )
   reverse_lines=$( cat reverse_reads.fastq | wc -l )
   forward_seq=$(( forward_lines / 4 ))
   reverse_seq=$(( reverse_lines / 4 ))
   if [ $forward_seq -ne $reverse_seq ]
   then
     echo "Forward and Reverse files do not have same number of sequences" >> $logfile
     email_exit
     #exit 1
   fi
   rem=$(( forward_lines % 4 ))
   if [ $rem -ne 0 ]
   then 
     echo "Forward sequences are more than 4 lines" >> $logfile
     email_exit
     #exit
   fi
   rem=$(( reverse_lines % 4 ))
   if [ $rem -ne 0 ]
   then
     echo "Reverse sequences are more than 4 lines" >> $logfile
     email_exit
     #exit 1
   fi

   #echo " " >> $logfile
   echo "barcode file has $barcode_lines lines" >> $logfile
   echo "forward reads file has $forward_lines = $forward_seq sequences" >> $logfile
   echo "reverse reads file has $reverse_lines = $reverse_seq sequences" >> $logfile
   #echo " " >> $logfile


   echo " " >> $logfile
   echo "Running join_paired_ends.py" >> $logfile
   #echo " " >> $logfile
   join_paired_ends.py -f forward_reads.fastq -r reverse_reads.fastq -b barcode.fastq -o STEP1_OUT/
   cd STEP1_OUT

   echo "join_paired_ends.py -f forward_reads.fastq -r reverse_reads.fastq -b barcode.fastq -o STEP1_OUT/" >> $logfile
   echo "Made directory STEP1_OUT" >> $logfile
   echo "made files:" >> $logfile
   for f in fastqjoin.join_barcodes.fastq fastqjoin.join.fastq fastqjoin.un1.fastq fastqjoin.un2.fastq
   do
      file_lines=$( cat $f | wc -l )
      file_seq=$(( file_lines / 4 ))
      echo "    $f ($file_lines lines = $file_seq sequences)" >> $logfile
   done

   echo " " >> $logfile
   echo "Removing unzipped data files" >> $logfile
   echo " " >> $logfile
   rm -f ../forward_reads.fastq
   rm -f ../reverse_reads.fastq
   rm -f ../barcode.fastq

   #echo " " >> $logfile
   echo "Running split_libraries.py" >> $logfile
   #echo " " >> $logfile
   #cd STEP1_OUT
   split_libraries_fastq.py -i fastqjoin.join.fastq -b fastqjoin.join_barcodes.fastq --rev_comp_mapping_barcodes -o STEP2_OUT/ -m $4  -q 19 --store_demultiplexed_fastq
   cd STEP2_OUT

else
    #starting at step 2
    split_libraries_fastq.py -i $2 -b $3 --rev_comp_mapping_barcodes -o STEP2_OUT/ -m $4  -q 19 --store_demultiplexed_fastq
    cd STEP2_OUT
fi


echo "split_libraries_fastq.py -i fastqjoin.join.fastq -b fastqjoin.join_barcodes.fastq --rev_comp_mapping_barcodes -o STEP2_OUT/ -m $4  -q 19" >> $logfile
echo "made directory STEP2_OUT" >> $logfile
cmd_out=$(grep -c '>' seqs.fna)
echo "made seqs.fna ($cmd_out sequences)" >> $logfile

#cd STEP2_OUT

#Do not split the file because usearch 64 bit is installed
echo " " >> $logfile
echo "Looking for chimeric sequences in seqs.fna" >> $logfile
identify_chimeric_seqs.py -i seqs.fna -m usearch61 -r /home2/Database/RDP_Gold/rdp_gold.fa -o usearch61_chimera/
echo "identify_chimeric_seqs.py -i seqs.fna -m usearch61 -r /home2/Database/RDP_Gold/rdp_gold.fa -o usearch61_chimera/" >> $logfile
echo "made 11 files in directory usearch61_chimera" >> $logfile
cmd_out=$( wc -l  < usearch61_chimera/chimeras.txt )
echo "    chimeras.txt ($cmd_out lines = chimeras found)" >> $logfile

#filter the chimers from the seq.fna file
echo "Removing chimeras from seqs.fna" >> $logfile
filter_fasta.py -f seqs.fna -o seqs_chimeras_filtered.fna -s usearch61_chimera/chimeras.txt -n
echo "filter_fasta.py -f seqs.fna -o seqs_chimeras_filtered.fna -s usearch61_chimera/chimeras.txt -n" >> $logfile
cmd_out=$(grep -c '>' seqs_chimeras_filtered.fna)
echo "   seqs_chimeras_filtered.fna ($cmd_out sequences)" >> $logfile


#get the current directory because next instruction needs full paths
current_dir=$(pwd)

echo " " >> $logfile
echo "Running pick_open_reference_otus.py" >> $logfile
#echo " "  >> $logfile
pick_open_reference_otus.py -i  $current_dir/seqs_chimeras_filtered.fna -r /home2/Database/Silva/rep_set/97_Silva_111_rep_set.fasta -o $current_dir/STEP3_OUT -f -a -O 40
cd STEP3_OUT

echo "pick_open_reference_otus.py -i  $current_dir/seqs_chimeras_filtered.fna -r /home2/Database/Silva/rep_set/97_Silva_111_rep_set.fasta -o $current_dir/STEP3_OUT -f -a -O 40" >> $logfile
echo "made directory STEP_3_OUT" >> $logfile
cmd_out=$( grep -c '>' rep_set.fna )
echo "     rep_set.fna ($cmd_out sequences)" >> $logfile
#cd STEP3_OUT



summarize_taxa.py -i otu_table_mc2_w_tax.biom -o taxonomy_summaries/  -L 2,3,4,5,6 
echo " " >> $logfile
echo "summarize_taxa.py -i otu_table_mc2_w_tax.biom -o taxonomy_summaries/  -L 2,3,4,5,6" >> $logfile
echo "makes directory taxonomy_summaries with 6 biom and 6 txt files of otu tables 2-6" >> $logfile
cmd_out=$( wc -l < taxonomy_summaries/otu_table_mc2_w_tax_L6.txt )
echo "      taxonomy_summaries/otu_table_mc2_w_tax_L6.txt ($cmd_out lines)" >> $logfile

#biom summarize-table -i otu_table_mc2_w_tax.biom -o rich_sparse_otu_table_summary.txt
#echo " " >> $logfile
#echo "biom summarize-table -i otu_table_mc2_w_tax.biom -o rich_sparse_otu_table_summary.txt" >> $logfile
#cmd_out=$( wc -l < rich_sparse_otu_table_summary.txt )
#echo "     rich_sparse_otu_table_summary.txt ($cmd_out lines)" >> $logfile

#extract_taxonomy_v2.py -r uclust_assigned_taxonomy/rep_set_tax_assignments.txt -t taxonomy_summaries/otu_table_mc2_w_tax_L7.txt  -p rep_set.fna -d 7 -o taxonomy_summaries/otu_table_mc2_w_tax_L7_with_accession.txt  -f rep_set_with_taxonomy.fna
#echo " " >> $logfile
#echo "extract_taxonomy_v2.py -r uclust_assigned_taxonomy/rep_set_tax_assignments.txt -t taxonomy_summaries/otu_table_mc2_w_tax_L7.txt  -p rep_set.fna -d 7 -o taxonomy_summaries/otu_table_mc2_w_tax_L7_with_accession.txt -f rep_set_with_taxonomy.fna" >> $logfile
#cmd_out=$( wc -l < taxonomy_summaries/otu_table_mc2_w_tax_L7_with_accession.txt )
#echo "     taxonomy_summaries/otu_table_mc2_w_tax_L7_with_accession.txt ($cmd_out lines)" >> $logfile
#cmd_out=$( wc -l < rep_set_with_taxonomy.fna )
#echo "     rep_set_with_taxonomy.fna ($cmd_out lines)" >> $logfile

#added per lindsey

python /ORG-Data/scripts/wrapper_filter_otus_from_otu_table.py -i otu_table_mc2_w_tax.biom -o percent_filtered_otu_table_mc2_w_tax.biom -n $6 -$7 $8

#biom convert -i otu_table_mc2_w_tax.biom -o otu_table_mc2_w_tax.txt --to-tsv --header-key taxonomy
biom convert -i percent_filtered_otu_table_mc2_w_tax.biom -o percent_filtered_otu_table_mc2_w_tax.txt --to-tsv --header-key taxonomy

python /ORG-Data/scripts/calculate_relative_abundance.py -i percent_filtered_otu_table_mc2_w_tax.txt -o percent_calculate_relative_abundance_output.txt

echo "percent_calculate_relative_abundance_output.txt file made"
echo "percent_calculate_relative_abundance_output.txt file made" >> $logfile

echo " " >> $logfile
echo "Qiime script has completed" >> $logfile
echo "Qiime script has completed"
echo " " >> $logfile

#get the current dir
end_dir=$(pwd)
echo "Qiime was run on the mapping file $4  and the OTU tables are located at:" >> $logfile
echo "$end_dir/taxonomy_summaries" >> $logfile
echo " " >> $logfile


email_exit
#exit

