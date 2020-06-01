Scripts for processing short amplicon reads from Oxford Nanopore. For details, please see:

ST Calus, UZ Ijaz, and A Pinto. NanoAmpli-Seq: A workflow for amplicon sequencing from mixed microbial communities on the nanopore sequencing platform. bioRxiv 244517, 2018. DOI:10.1101/244517


Sam Young

Modified the program to run in python3 then add in multprocessing to decrease time overall. Ran several test to see what preformed the best.  Decreased the overall time by about half for the changing it to python3. 

Speed test conclude that more processors would reduce time. Multithreading is not optimal for this problem. Just transitioning the older code into python3 reduced time significantly. 

linux terminal or mac terminals

First make sure your default python = python3 

./chopSEQ.py -i input_file -f foward_primer -r reserve_primer > output_file 

otherwise python3 chopSEQ.py -i input_file -f foward_primer -r reversre_primer > output_file 

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
|sample_travis2_p4.1.fa|589|Eight processors used|(python3)|44:51|


