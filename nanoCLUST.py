#!/usr/bin/python
# ***************************************************************
# Name:      nanoCLUST.py
# Dependencies:
#	     biopython
#	     vsearch
#	     mafft-ginsi
#	     Python scripts, uc2otutab.py in particular from https://drive5.com/python/python_scripts.tar.gz 
#	     (unzip them in ~/bin directory and then chmod +x *)	
# Version:   0.3
# History:   Added detailed help and removed inc-seq dependency
#	     Added dynamic tandem threshold and fixed a few bugs (i.e. dynamic thresholding
#	     lower/upper bound for consensus reads)
# Authors:   Umer Zeeshan Ijaz (Umer.Ijaz@glasgow.ac.uk)
#                 http://userweb.eng.gla.ac.uk/umer.ijaz
# Created:   2017-06-28
# License:   Copyright (c) 2017 Environmental'Omics Lab, University of Glasgow, UK
#
#            This program is free software: you can redistribute it and/or modify
#            it under the terms of the GNU General Public License as published by
#            the Free Software Foundation, either version 3 of the License, or
#            (at your option) any later version.
#
#            This program is distributed in the hope that it will be useful,
#            but WITHOUT ANY WARRANTY; without even the implied warranty of
#            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#            GNU General Public License for more details.
#
#            You should have received a copy of the GNU General Public License
#            along with this program.  If not, see <http://www.gnu.org/licenses/>.
# **************************************************************/
import os,os.path
import sys, getopt
import time,datetime
import numpy
import subprocess
import math
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import pairwise2
from Bio.Seq import Seq

def usage():
	print 'Usage:'
	print '\tpython nanoCLUST.py -i <input_file>'
	print '''
	Other options are:
	-s (--split_option)	Accepts a comma-delimited list of coordinates to split, e.g, 
				"0,200,400,1500" or a single number of required segments
	-o (--output_folder)	Folder to dump the files. Default folder name is "nanoCLUSTdef"


	For the script to work, you need to ensure the <input_file> is in USEARCH/VSEARCH format with "barcodelabel=X;" 
	prepended to reads to distinguish between samples, a requirement to generate the OTU table

	Here is an awk one-liner to prepend barcode labels for multiple samples:

	awk -v k="Sample1" '/^>/{gsub(">","",$0);$0=">barcodelabel="k";"$0}1' Sample1.fa > Sample1_barcode.fa
	awk -v k="Sample2" '/^>/{gsub(">","",$0);$0=">barcodelabel="k";"$0}1' Sample2.fa > Sample2_barcode.fa
	awk -v k="Sample3" '/^>/{gsub(">","",$0);$0=">barcodelabel="k";"$0}1' Sample3.fa > Sample3_barcode.fa

	cat *_barcode.fa > multiplexed.fa

	'''

def create_folder(directory):
	if not os.path.exists(directory):
		print 'Creating the folder ' + directory
    		os.makedirs(directory)
	else:
		print 'Folder ' + directory + ' already exists! Skipping!'

def remove_file(path):
	print 'Removing ' + path
	os.remove(path)

def print_time_stamp(message):
    print datetime.datetime.fromtimestamp(time.time()).strftime('[%Y-%m-%d %H:%M:%S] ')+message

#check if the program exists and return it's location, otherwise return None
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def print_cond(cond,true_message,false_message):
    if cond:
    	print true_message
    else:
	print false_message	

def print_prog_status(prog,loc):
    if loc != None:
	print 'Checking for \'' + prog + '\': found '+ loc
    else:
	print 'Checking for \'' + prog + '\': ERROR - could not find \''+prog+'\''
	print 'Exiting.'
	sys.exit(1)

def run_prog(prog):
	p = subprocess.Popen(prog, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	returned_list=[]
	for line in iter(p.stdout.readline, b''):
		returned_list.append(line.rstrip())
    	retval = p.wait()
	return returned_list


def split_fasta(directory_name,file_name,parts,average_length):
	#parts can accept a comma-delimited list of coordinates, e.g, "0,200,400,1500" or a single number of required segments
	coordinates=[]
	if "," in parts:
		#The coordinates are available
		coordinates=parts.split(",")
		coordinates=[int(x) for x in coordinates]
		if ((divmod(len(coordinates),2)[1]) !=0):
			print "Error - Coordinates are not well defined"
			sys.exit(2)
	else:
		#Split into parts
		results=divmod(int(average_length),int(parts))
		print_time_stamp("Splitting " + file_name + "(Average length="+ str(average_length) + ") into " + parts + " parts with window size of " + str(results[0]))
		for i in range(int(parts)):
	                start=results[0]*i
        	        end=results[0]*i+results[0]-1
                	end_look_ahead=results[0]*(i+1)+results[0]-1
                	if end_look_ahead>average_length:
                        	end=-1
			coordinates.append(start)
			coordinates.append(end)			

	#Now split the fasta file
	ind=1	
	for i in range(0,len(coordinates),2):
		oname=directory_name+"/p"+str(ind)+".fasta"
		of=open(oname,"w")
		for seq_record in SeqIO.parse(file_name,"fasta"):
			if(len(seq_record.seq)>0):
				of.write(">"+str(seq_record.id)+"\n")
				of.write(str(seq_record.seq)[coordinates[i]:coordinates[i+1]]+"\n")
		of.close()
		print_time_stamp("Generated " + oname)
		ind=ind+1
	return ind-1
	
def extract_uc_record(file_name):
	print file_name
	otus_assignments=dict()

	iff=open(file_name,"r")
	ind=1
	for line in iff:
		line=line.rstrip()
		record=line.split("\t")
		#if the key exists
		if record[-1] in otus_assignments:
			rec=otus_assignments[record[-1]]
			rec.append(record[-2])
			otus_assignments[record[-1]]=rec
		else:
			otus_assignments[record[-1]]=[record[-2]]
		ind=ind+1
	iff.close()
	return otus_assignments

def run_vsearch(prog_name,directory_name,file_name,file_prefix):
	iname=file_name
	oname=directory_name+'/'+file_prefix+'_u.fasta'
	print_time_stamp('Dereplicating ' + iname)
	cmd=prog_name + ' --derep_fulllength ' + iname + ' --output ' + oname + ' --sizeout --minseqlength 50'
	print cmd
	run_prog(cmd)
	print_time_stamp('Generated ' + oname)


	iname=oname
	oname=directory_name+'/'+file_prefix+'_n.fasta'
	print_time_stamp('Removing chimerias from ' + iname)
	cmd=prog_name + ' --uchime_denovo ' + iname + ' --nonchimeras ' + oname
	print cmd
	run_prog(cmd)
	print_time_stamp('Generated ' + oname)

        iname=oname
        oname=directory_name+'/'+file_prefix+'_s.fasta'
        print_time_stamp('Sorting ' + iname + ' by size')
	cmd=prog_name + ' --sortbysize ' + iname + ' --output ' + oname + ' --minsize 2'
        print cmd
        run_prog(cmd)
        print_time_stamp('Generated ' + oname)

        iname=oname
        oname=directory_name+'/'+file_prefix+'_r.fasta'
        print_time_stamp('Clustering ' + iname + ' at 97% ident')
	cmd=prog_name + ' --cluster_smallmem ' + iname + ' --id 0.97 --consout ' + oname + ' --usersort'
        print cmd
        run_prog(cmd)
        print_time_stamp('Generated ' + oname)

	iname=oname
	oname=directory_name+'/'+file_prefix+'_e.fasta'
	print_time_stamp('Relabeling ' + iname)
	of=open(oname,"w")
	ind=1
	for seq_record in SeqIO.parse(iname,"fasta"):
		of.write(">OTU_"+str(ind)+'\n')
		of.write(str(seq_record.seq)+'\n')
		ind=ind+1
	of.close()
	print_time_stamp('Generated ' + oname)

        iname=oname
        oname=directory_name+'/'+file_prefix+'.uc'
        print_time_stamp('Searching ' + file_name + ' against OTUs ' + iname)
	cmd=prog_name + ' --usearch_global ' + file_name + ' --db ' + iname + ' --strand both --id 0.97 --uc ' + oname + ' --threads 20'
        print cmd
        run_prog(cmd)
        print_time_stamp('Generated ' + oname)

	remove_file(directory_name+'/'+file_prefix+'_u.fasta')
	remove_file(directory_name+'/'+file_prefix+'_n.fasta')
	remove_file(directory_name+'/'+file_prefix+'_s.fasta')
	remove_file(directory_name+'/'+file_prefix+'_r.fasta')
		

def main(argv):
	input_file=''
	split_option='5'
	output_folder='nanoCLUSTdef'
	consensus_option='mafft'
	consensus_sequences=50
	consensus_average_variability_percentage=10
	try:
		opts,args=getopt.getopt(argv,"hi:o:s:c:",["input_file","output_folder","split_option","consensus_option"])
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			usage()
			sys.exit()
		elif opt in ("-i","--input_file"):
			input_file=arg
                elif opt in ("-s","--split_option"):
                        split_option=arg
                elif opt in ("-c","--consensus_option"):
                        consensus_option=arg
                elif opt in ("-o","--output_folder"):
                        output_folder=arg

	if (input_file==''):
		usage()
		sys.exit(2)

        PROG_VSEARCH=which('vsearch')
    	print_prog_status("vsearch",PROG_VSEARCH)
	PROG_UC2OTUTAB=which('uc2otutab.py')
	print_prog_status("uc2otutab.py",PROG_UC2OTUTAB)
	PROG_MAFFTGINSI=which('mafft-ginsi')
        print_prog_status("mafft-ginsi",PROG_MAFFTGINSI)

	create_folder(output_folder)

	#load sequences in the memory so that we can extract them later, and also calculate average_length
	memory_sequences=dict()
	ind=1
	average_length=0.0
	for seq_record in SeqIO.parse(input_file,"fasta"):
		memory_sequences[seq_record.id]=str(seq_record.seq)
                average_length=average_length+len(seq_record.seq)
                ind=ind+1
        average_length=average_length/(ind-1)


	number_of_files=split_fasta(output_folder,input_file,split_option,average_length)
	
	for i in range(number_of_files):
		run_vsearch(PROG_VSEARCH,output_folder,output_folder+"/p"+str(i+1)+".fasta","p"+str(i+1))

	num_otus=-1
	file_name=""
	otu_assignments=dict()
	for i in range(number_of_files):
		rec=extract_uc_record(output_folder+"/p"+str(i+1)+".uc")
		if(len(rec.keys())>num_otus):
			otu_assignments=rec
			num_otus=len(rec.keys())
			file_name=output_folder+"/p"+str(i+1)+".uc"
	

	cf=open(output_folder+"/"+"consensus.fasta","w")
	for i in otu_assignments.keys():
		if consensus_option=="mafft":
			print_time_stamp("Generating " + output_folder+"/"+i+".fasta")
			of=open(output_folder+"/"+i+".fasta","w")
			ind=1
			average_length_bin=0.0
			for j in otu_assignments[i]:
				average_length_bin=average_length_bin+len(memory_sequences[j])
				ind=ind+1
			average_length_bin=average_length_bin/(ind-1)
			ind=0
			local_average_variability_percentage=consensus_average_variability_percentage
			while ind<2:
				lower_threshold=average_length_bin-((float(local_average_variability_percentage)/100.0)*average_length_bin)
        			upper_threshold=average_length_bin+((float(local_average_variability_percentage)/100.0)*average_length_bin)		
				print "Taking at most " + str(consensus_sequences) + " sequences at " + str(local_average_variability_percentage) + "% of average bin size ("+ str(lower_threshold) + "," + str(average_length_bin) + "," + str(upper_threshold) + ")" 

				for j in otu_assignments[i]:
					if ((len(memory_sequences[j]) >= lower_threshold) and (len(memory_sequences[j]) <= upper_threshold)):
						of.write(">"+str(j)+"\n")
						of.write(memory_sequences[j]+"\n")
						ind=ind+1
						if ind>consensus_sequences:
							break
				print "Using",ind-1,"sequences"
				if ind<2:
					local_average_variability_percentage=local_average_variability_percentage+10

			of.close()


			print_time_stamp("Aligning " + output_folder+"/"+i+".fasta")
			cmd=PROG_MAFFTGINSI + " " + output_folder+"/"+i+".fasta" + " > " + output_folder+"/"+i+".gfasta"
			print cmd
			run_prog(cmd)
			print_time_stamp("Generated " + output_folder+"/"+i+".gfasta")
			alignments=AlignIO.read(output_folder+"/"+i+".gfasta","fasta")
			summary_align = AlignInfo.SummaryInfo(alignments)
			consensus_read=summary_align.gap_consensus(threshold=0.1,consensus_alpha=None, require_multiple=1)
			cf.write(">"+i+"\n")
			cf.write(str(consensus_read).replace("-","").replace("X","N").upper()+"\n")
                	#remove_file(output_folder+"/"+i+".fasta")	
			remove_file(output_folder+"/"+i+".gfasta")


	cf.close()
	print_time_stamp("Generated " + output_folder+"/consensus.fasta")

	print PROG_UC2OTUTAB
	cmd = "python " + PROG_UC2OTUTAB + " " + file_name + " > " + output_folder + "/otu_table.txt"
	print_time_stamp('Generating otu_table for ' + file_name)
	print cmd
	run_prog(cmd)
	print_time_stamp('Generated '+output_folder+'/otu_table.txt')

	

if __name__== "__main__":
	main(sys.argv[1:])
			
