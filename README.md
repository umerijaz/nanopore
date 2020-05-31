Scripts for processing short amplicon reads from Oxford Nanopore. For details, please see:

ST Calus, UZ Ijaz, and A Pinto. NanoAmpli-Seq: A workflow for amplicon sequencing from mixed microbial communities on the nanopore sequencing platform. bioRxiv 244517, 2018. DOI:10.1101/244517


Modified the program to use multprocessing to speed up and updated the code to run in python3.  Decreased the overall time by about half for the changing it to python3. Running several speed test. I'm only going to run one multithreading unless it under about 2:23:51 

./chopSEQ.py -i input_file -f foward_primer -r reserve_primer > output_file 



#running times
##number           machine types                               times 
***419     -        four core machine                         -  1:43:39***
419     -       single core machine                        -  2:33:51
419     -    multithreading four core machine              -  ~4.8
419     -       eight core machine                         - 
419     -      graphic acceleration                        -

859     -       four core machine                          - 
859     -       single core machine                        - 
859     -       eight core machine                         -
859     -      graphic acceleration                        -
