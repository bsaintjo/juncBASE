#!/bin/bash
#Maximillian Marin (Last Updated: 4/3/16)
#Wrapper for JuncBASE  

mkdir JBinput
mkdir JBoutput

SampleToBAM=$1
IntronJunction=$2
processes=$3
RefGenome=$4
ListOfSampleNames=$5
DBdir=$6
DBone=$7 
DBtwo=$8 
DBthree=$9
pseudo=$(head -1 ${ListOfSampleNames})



echo STEP 1: Process SAM/BAM files
#python /pod/home/mgmarin/usr/JuncBASE/run_preProcess_by_chr_step1.py -i $SampleToBAM -o JBinput -p $processes --preProcess_options "--unique -j ${IntronJunction} -c 1.0"  

echo Step 1B: Disambiguate strand of splice junctions
#python /pod/home/mgmarin/usr/JuncBASE/disambiguate_junctions.py -i JBinput -g $RefGenome --by_chr 

echo STEP 2: Identify all junctions 
#This step will create one BED file that combines junctions from all samples.
#python /pod/home/mgmarin/usr/JuncBASE/preProcess_getASEventReadCounts_by_chr_step2.py -i JBinput --by_chr

echo STEP 3: Create exon-intron junction count files
#python /pod/home/mgmarin/usr/JuncBASE/run_preProcess_step3_by_chr.py --input_dir JBinput --num_processes $processes 

echo STEP 4: Create a pseudo/"all junction" sample
echo your pseudo sample is $pseudo
#python /pod/home/mgmarin/usr/JuncBASE/createPseudoSample.py -i JBinput -s $pseudo --by_chr

echo STEP 5: Identify alternative splicing events and quantify events from each sample
#This step of JuncBASE that identifies, classifies, and quantifies alternative splicing events
python /pod/home/mgmarin/usr/JuncBASE/run_getASEventReadCounts_multiSample.py -s $ListOfSampleNames -i JBinput -o JBoutput --sqlite_db_dir . --txt_db1 $DBone --txt_db2 $DBtwo --txt_db3 $DBthree --jcn_seq_len 140 -p $processes --by_chr 

echo "STEP 6: Create tables of raw and length-normalized read counts of exclusion and inclusion isoforms"
#This step collects information from runs performed in Step 5 and reports the AS event and read counts for every event in every sample in a table. 
python /pod/home/mgmarin/usr/JuncBASE/run_createAS_CountTables.py -d JBoutput -i JBinput --jcn_seq_len 140 -s $ListOfSampleNames --num_processes $processes

python /pod/home/mgmarin/usr/JuncBASE/combine_createAS_CountTables_by_chr.py -d JBoutput -o step6
