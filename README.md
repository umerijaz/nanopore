Scripts for processing short amplicon reads from Oxford Nanopore. For details, please see:

ST Calus, UZ Ijaz, and A Pinto. NanoAmpli-Seq: A workflow for amplicon sequencing from mixed microbial communities on the nanopore sequencing platform. bioRxiv 244517, 2018. DOI:10.1101/244517


Modified the program to use multprocessing to speed up and updated the code to run in python3.  Decreased the overall time by about half for the changing it to python3. Running several speed test. 

linux terminal 
./chopSEQ.py -i input_file -f foward_primer -r reserve_primer > output_file 

running times
-------------


|data sets|number|proccesing type|python3|times|
|--------------------|-----|----------------------|-------|--------|
|sample_travis2_p2.fa|419|Multithreading usage|(python3)|4:42:33|
|sample_travis2_p2.fa|419|Single proccesors used|(python3)|2:33:51|
|sample_travis2_p2.fa|419|Four processor used|(python3)|1:43:39|
|sample_travis2_p2.fa|419|Eight processor used|(python3)|1:15:52|
|sample_travis2_p4.1.fa|589|single processor used   |python3|1:54:53|
|sample_travis2_p4.1.fa|589|four processors used|(python3)|1:26:40|
|sample_travis_p4.1.fa|589|Eight processors used|(python3)|44:51|


