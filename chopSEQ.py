#!/usr/bin/python
# ***************************************************************
# Name:      chopSEQ.py
# Purpose:   Fixes the problems in consensus sequences produced by INC-SEQ by
#	     a) Rearranging the sequences to start from forward primer
#	     b) Remove tandom repeats from the resulting rearrangement
#	     c) Corrects the orientation and also colors the primers and tandem repeats in verbosity mode -v
#	     d) Apply read length thresholds
# Dependencies:
#	     Biopython is installed
#	     EMBOSS etandem is in path
#            conconcert.futures
#            psutil
# Version:   0.3.3 (2020-05-29)
# History:   Fixed the bugs associated with right overhang and included thresholds for final read lengths
#	     Added dynamic tandem threshold and fixed a few bugs
# Authors:   Umer Zeeshan Ijaz (Umer.Ijaz@glasgow.ac.uk), Young Sam Eugene (samuel.young.103@gmail.com)`
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
# Modifed:   2020-05-29
# **************************************************************/
from getopt import getopt, GetoptError 
from sys import argv, exit`
import numpy
import subprocess
import math
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from psutil import cpu_count
from time import perf_counter 

def usage():
    print ('Usage:')
    print("\tYour total amount of cpu: %d\n\tHalf the amount of cpu is %d"%(cpu_count(),cpu_count(logical=False)))
    print ('\tpython chopSEQ.py -i <input_file> -f <forward_primer> -r <reverse_primer> -p <number of process>  -t > <filtered_fasta_file> ')
    print ('''
    Other options are:
    -l (--minimum_length)   After processing, ignore reads below this length threshold 
    -m (--maximum_length)   After processing, ignore reads above this length threshold
    -v	                    Verbosity switch to visualise primers across the length of the reads
    -t (--threading)        this allow for multithreading at maximum threads
    -p (--processing)       This allow for the user to input the number of process the user want
                            to use default is a single processor being used 
    ''')

def run_prog(prog):
    p = subprocess.Popen(prog, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    returned_list=[]
    for line in iter(p.stdout.readline, b''):
        returned_list.append(line.rstrip())
    retval = p.wait()

    return returned_list

def process_etandem_record(full_sequence,verbosity,tandem_min_repeat,tandem_max_repeat,tandem_threshold,tandem_mismatch,tandem_identity_threshold):
    cmd_line="etandem -sequence=asis:" + full_sequence + " -minrepeat=" + str(tandem_min_repeat) + " -maxrepeat=" + str(tandem_max_repeat) + " -threshold=" + str(tandem_threshold) + " -mismatch=" + str(tandem_mismatch) + " -stdout -auto -rformat excel"
    etandem_output=run_prog(cmd_line)
    corrected_sequence=''
    if(len(etandem_output)>1):
        record=etandem_output[1].split("\t")
        t_match_left_position=int(record[1])
        t_match_right_position=int(record[2])
        t_match_count=int(record[6])
        t_match_score=int(record[3])
        t_match_length=int(record[5])
        t_match_ident=float(record[7])
        t_match_consensus=record[8]
        t_left_string=full_sequence[0:(t_match_left_position-1)]
        t_tandom_repeat_1=full_sequence[(t_match_left_position-1):(t_match_left_position+t_match_length-1)]
        t_tandom_repeat_2=full_sequence[(t_match_left_position+t_match_length-1):t_match_right_position]
        t_right_string=full_sequence[t_match_right_position:]

		#Calculate allowed ident threshold
		#As the tandem repeat becomes longer, there is a high probability of incorporating the errors
		#For this purpose we divide the range from tandem_min_repeat to tandem_max_repeat in 80 parts
		#and then subtract the number of parts from the tandem_identity_threshold
        allowed_ident_threshold=(tandem_identity_threshold - math.floor(t_match_length / float((float(tandem_max_repeat)-float(tandem_min_repeat))/80.0)))
        if verbosity:
            if t_match_ident>=allowed_ident_threshold:
                 print ("Step 2: Removing tandom repeats (identity >= "+str(allowed_ident_threshold)+ ") in correctly arranged sequence" +  "[length="+str(t_match_length)+";score=" + str(t_match_score) + ";ident=" + str(t_match_ident) + ";count=" + str(t_match_count)+ ";consensus="+t_match_consensus+ "]")
                 print(t_left_string+'\x1b[6;30;44m'+t_tandom_repeat_1+'\x1b[0m'+ '\x1b[6;30;45m' + t_tandom_repeat_2+'\x1b[0m'+t_right_string)
        if t_match_ident>=allowed_ident_threshold:
            corrected_sequence=t_left_string+t_tandom_repeat_1+t_right_string
        else:
            corrected_sequence=full_sequence
    else:
        corrected_sequence=full_sequence
        
    return corrected_sequence


def main(arguv):
    input_file=''
    forward_primer=''
    reverse_primer=''

	# *Default parameters *************	
    # *********************************/
    verbosity  = 0  
    processing = 1
    threading  = 0
    #Used in Step 3:
    minimum_read_length_threshold = 0
    maximum_read_length_threshold = 100000
        
        #this code is used to put options in the arguments in the proper place  I'm having to spilit up the options in this case becuse they do not read from commandline properly
    try:
        opts,args=getopt(arguv,"hi:f:r:vl:m:p:t",["input_file","forward_primer","reverse_primer","minimum_length","maximum_length","multiprocessing","multithreading"])
    except GetoptError:
            usage()
            exit(2)
    for opt, arg in opts:
           if opt == '-h':
                usage()
                exit()
           elif opt == '-v':
               verbosity=1
           elif opt in ("-i","--input_file"):
               input_file=arg
           elif opt in ("-f","--forward_primer"):
               forward_primer=arg
           elif opt in ("-r","--reverse_primer"):
               reverse_primer=arg
           elif opt in ("-l","--minimum_length"):
               minimum_read_length_threshold=int(arg)
           elif opt in ("-m","--maximum_length"):
               maximum_read_length_threshold=int(arg)
           elif opt in ("-p","--processing"):
               processing = int(arg)
           elif opt in ("-t","--threading"):
               threading = 1
    if input_file=='' or forward_primer=='' or reverse_primer=='' or  processing <=0 or processing > cpu_count():
        usage()
        exit(2)


	#Some of the references we used in writing the code
	#Reference: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc85
	#Reference: https://stackoverflow.com/questions/287871/print-in-terminal-with-colors-using-pyton
    #this allow the multiprocessor call and allow for the ability
    if processing == cpu_count():
        with ProcessPoolExecutor() as executor:
            results =[executor.submit(process_seq_records, seq_record, forward_primer, reverse_primer, verbosity, minimum_read_length_threshold, maximum_read_length_threshold) for seq_record in SeqIO.parse(input_file,"fasta")]

            for resultant in as_completed(results):
                print(resultant.result())


    elif processing:
        with ProcessPoolExecutor(max_workers=processing) as executor:
            results =[executor.submit(process_seq_records, seq_record, forward_primer, reverse_primer, verbosity, minimum_read_length_threshold, maximum_read_length_threshold) for seq_record in SeqIO.parse(input_file,"fasta")]

            for resultant in as_completed(results):
                print(resultant.result())

    elif threading:
        with ThreadPoolExecutor() as executor:
            results = [executor.submit(process_seq_records,seq_record,forward_primer,reverse_primer,verbosity,minimum_read_length_threshold,maximum_read_length_threshold) for seq_record in SeqIO.parse(input_file,"fasta")]
            for resultant in as_completed(results):
                print(resultant.result())



def process_seq_records(seq_record,forward_primer,reverse_primer,verbosity,minimum_read_length_threshold,maximum_read_length_threshold):
    #used in step1
    primer_match_score=5
    primer_mismatch_score=-4
    primer_open_gap_score=-2
    primer_extend_gap_score=-2
    #Used in Step 2:
    tandem_min_repeat=10
    tandem_max_repeat=350
    tandem_threshold=10
    tandem_mismatch=5
    tandem_identity_threshold=85
    

    if verbosity:
        print('\x1b[6;37;40m' + str(seq_record.id) + '\x1b[0m')
    

    forward_orientation_forward_primer_alignment=pairwise2.align.localms(str(seq_record.seq),forward_primer,primer_match_score, primer_mismatch_score, primer_open_gap_score, primer_extend_gap_score, one_alignment_only=1)
    
    forward_orientation_reverse_primer_alignment=pairwise2.align.localms(str(seq_record.seq),str(Seq(reverse_primer).reverse_complement()),primer_match_score,primer_mismatch_score, primer_open_gap_score, primer_extend_gap_score, one_alignment_only=1)

    reverse_orientation_forward_primer_alignment=pairwise2.align.localms(str(seq_record.seq),reverse_primer,primer_match_score,primer_mismatch_score, primer_open_gap_score, primer_extend_gap_score, one_alignment_only=1)
    
    reverse_orientation_reverse_primer_alignment=pairwise2.align.localms(str(seq_record.seq),str(Seq(forward_primer).reverse_complement()),primer_match_score,primer_mismatch_score, primer_open_gap_score, primer_extend_gap_score, one_alignment_only=1)		
        #Find the correct orientation based on mean scores of primers matching

    if numpy.mean([forward_orientation_forward_primer_alignment[0][2],forward_orientation_reverse_primer_alignment[0][2]]) > numpy.mean([reverse_orientation_forward_primer_alignment[0][2],reverse_orientation_reverse_primer_alignment[0][2]]):


        f_match_left_position=forward_orientation_forward_primer_alignment[0][3]
        f_match_right_position=forward_orientation_forward_primer_alignment[0][4]
        f_left_string=str(forward_orientation_forward_primer_alignment[0][0][0:(f_match_left_position)])
        f_match_string=str(forward_orientation_forward_primer_alignment[0][0][f_match_left_position:f_match_right_position])
        f_right_string=str(forward_orientation_forward_primer_alignment[0][0][(f_match_right_position):])



        if verbosity:
            print ("Step 1: Match primers and correctly rearrange the sequence")
            print ('Forward Primer Match:')
            print(f_left_string + '\x1b[6;30;42m' + f_match_string + '\x1b[0m' + f_right_string)


        r_match_left_position=forward_orientation_reverse_primer_alignment[0][3]
        r_match_right_position=forward_orientation_reverse_primer_alignment[0][4]
        r_left_string=str(forward_orientation_reverse_primer_alignment[0][0][0:(r_match_left_position)])
        r_match_string=str(forward_orientation_reverse_primer_alignment[0][0][r_match_left_position:r_match_right_position])
        r_right_string=str(forward_orientation_reverse_primer_alignment[0][0][(r_match_right_position):])



        if verbosity:
            print ('Reverse Primer Match:')
            print(r_left_string + '\x1b[6;30;43m' + r_match_string + '\x1b[0m' + r_right_string)


        #See if there is a right-over hang after the reverse primer and remove it
        ind_r=f_right_string.find(r_match_string.replace("-",""))

        if (ind_r > -1):

            f_right_string=f_right_string[0:ind_r]


        full_sequence=f_match_string + f_right_string + f_left_string
        full_sequence=full_sequence.replace("-","")
        full_sequence=process_etandem_record(full_sequence,verbosity,tandem_min_repeat,tandem_max_repeat,tandem_threshold,tandem_mismatch,tandem_identity_threshold)
        if not verbosity:
            if ((len(full_sequence)>=minimum_read_length_threshold) and (len(full_sequence)<=maximum_read_length_threshold)):
                return ">" + str(seq_record.id) + "_corrected_length=" + str(len(full_sequence)) + "\n" + str(full_sequence)
    else:
        f_match_left_position=reverse_orientation_forward_primer_alignment[0][3]
        f_match_right_position=reverse_orientation_forward_primer_alignment[0][4]
        f_left_string=str(reverse_orientation_forward_primer_alignment[0][0][0:(f_match_left_position)])
        f_match_string=str(reverse_orientation_forward_primer_alignment[0][0][f_match_left_position:f_match_right_position])
        f_right_string=str(reverse_orientation_forward_primer_alignment[0][0][(f_match_right_position):])
        if verbosity:
            print ("Step 1: Match primers and correctly rearrange the sequence")
            print ('Reverse Primer Match (Reverse Orientation):')
            print(f_left_string + '\x1b[6;30;42m' + f_match_string + '\x1b[0m' + f_right_string)


        r_match_left_position=reverse_orientation_reverse_primer_alignment[0][3]
        r_match_right_position=reverse_orientation_reverse_primer_alignment[0][4]
        r_left_string=str(reverse_orientation_reverse_primer_alignment[0][0][0:(r_match_left_position)])
        r_match_string=str(reverse_orientation_reverse_primer_alignment[0][0][r_match_left_position:r_match_right_position])
        r_right_string=str(reverse_orientation_reverse_primer_alignment[0][0][(r_match_right_position):])
		
            #See if there is a right-over hang after the reverse primer and remove it
        ind_r=f_right_string.find(r_match_string.replace("-",""))


        if (ind_r > -1):
            f_right_string=f_right_string[0:ind_r]


        if verbosity:
            print ('Forward Primer Match (Reverse Orientation):')
            print(r_left_string + '\x1b[6;30;43m' + r_match_string + '\x1b[0m' + r_right_string)
                        
        full_sequence=f_match_string + f_right_string + f_left_string
        full_sequence=full_sequence.replace("-","")
        full_sequence=str(Seq(full_sequence).reverse_complement())


        full_sequence=process_etandem_record(full_sequence,verbosity,tandem_min_repeat,tandem_max_repeat,tandem_threshold,tandem_mismatch,tandem_identity_threshold)


        if not verbosity:
            if ((len(full_sequence)>=minimum_read_length_threshold) and (len(full_sequence)<=maximum_read_length_threshold)):

                return ">"+str(seq_record.id)+"_corrected_sequence="+str(len(full_sequence))+"\n"+str(full_sequence)


if __name__== "__main__":
    t1 = perf_counter()
    main(argv[1:])
    t2 = perf_counter()

    finish = t2-t1 

    secs = finish % 60
    mins = finish /60 % 60
    hours = finish / 3600 % 24


    print("Time trails for 3 processors %d : %d : %d" %(hours, mins , secs))

