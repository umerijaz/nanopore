Scripts for processing short amplicon reads from Oxford Nanopore. For details, please see:

ST Calus, UZ Ijaz, and A Pinto. NanoAmpli-Seq: A workflow for amplicon sequencing from mixed microbial communities on the nanopore sequencing platform. bioRxiv 244517, 2018. DOI:10.1101/244517


Modified the program to use multprocessing to speed up and updated the code to run in python3.  Decreased the overall time by about half for the changing it to python3. Running several speed test. I'm only going to run one multithreading unless it under about 2:23:51 

./chopSEQ.py -i input_file -f foward_primer -r reserve_primer > output_file 

table: 
*** best times ***
** fedora **
#running times#
##       data sets    (number)              machine types                               times   ##
**sample_travis2_p2.fa  (419)   -        four core machine                          -  1:43:39**
**sample_travis2_p2.fa  (419)   -       single core used                         -  2:33:51**
**sample_travis2_p2.fa  (419    -    multithreading four core machine               -  4:42:33**
***sample_travis2_p2.fa (419)   -       eight core machine   (bash on windows 10)    -  1:15:52***

**sample_travis2_p4.1.fa (598)     -       four core machine                          -  1:26:40 **
**sample_travis2_p4.1.fa (598)     -       single core used                           -  2:23:15s**
***sample_travis_p4.1.fa (598)     -       eight core machine   (bash on windows 10)  -   44:51 ***

