Scripts for processing short amplicon reads from Oxford Nanopore. For details, please see:

ST Calus, UZ Ijaz, and A Pinto. NanoAmpli-Seq: A workflow for amplicon sequencing from mixed microbial communities on the nanopore sequencing platform. bioRxiv 244517, 2018. DOI:10.1101/244517


Sam Young

dependencenies needed 

psutil 

biopython

concurrent

Modified the program to run in python3 then add in multprocessing to decrease time overall. Ran several test to see what preformed the best.  Decreased the overall time by about half for the changing it to python3. 

Speed test conclude that 3 cores are leads to the greatest is the amount of time effiecent. Multithreading is not optimal for chopping up sequences it lead to the slows amount of time. Just transitioning the older code into python3 reduced time significantly and then calling the imports directly decreased time. 

linux terminal or mac terminals

First make sure your default python = python3 

./chopSEQ.py -h this can help you understand some of the changes add into since last time

python3 chopSEQ.py -h 

This will allow the computer to run on six all core if the computer has that many otherwise it'll send you to usage 

./chopSEQ.py -i input_file -f foward_primer -r reservser_primer -p 6 > output_file 

python3 chopSEQ.py -i input_file -f foward_primers -r reverse_primers > \ >>  output_file 

this will allow the user to run on a single processor

./chopSEQ.py -i input_file -f foward_primer -r reverse_primer > output_file

python3 chopSEQ.py -i input_file -f foward_primer -r reversre_primer > \ >> output_file 

i7-4970 CPU @ 3.60 Ghz on chopSEQ running times 
-----------------------------------------------

Legend

best time ***1***

| data sets | No. total seqs | proccesing used type |python3| times | 
|--------------------------|----------------|----------------------|---------|--------|
|sample_Travis2_p2.fa|419| Multithreading |python3|1:49:06|
|sample_Travis2_p2.fa|419| Single core |python3||
|sample_Travis2_p2.fa|419| Two cores |python3||
|sample_Travis2_p2.fa|419| Three cores |python3||
|sample_Travis2_p2.fa|419| Four cores |python3||	
|sample_Travis2_p2.fa|419| Fivee cores|python3||
|sample_Travis2_p2.fa|419| Six cores |python3||
|sample_Travis2_p2.fa|419| Seven cores |python3||
|sample_Travis2_p2.fa|419| Eight cores |python3|1:15:52|
|sample_Travis2_p4.1.fa|589| Multithreading |python3|1:17:29|
|sample_Travis2_p4.1.fa|589| single core |python3|58:05|
|sample_Travis2_p4.1.fa|589| two cores |python3|35:18|
|sample_Travis2_p4.1.fa|589| three cores |python3|***31:43***|
|sample_Travis2_p4.1.fa|589| four cores |python3|36:31|
|sample_Travis2_p4.1.fa|589| five cores |python3|40:23|
|sample_Travis2_p4.1.fa|589| six cores  |python3|43:27|
|sample_Travis2_p4.1.fa|589| seven cores |python3|42:36|
|sample_Travis2_p4.1.fa|589| Eight cores |python3|43:05|
|sample_Travis2_pcombined.fa|1008| Multithreading |3:11:07|
|sample_Travis2_pcombined.fa|1008| Single core |python3|2:10:12|
|sample_Travis2_pcombined.fa|1008| Two cores |python3|1:26:07|
|sample_Travis2_pcombined.fa|1008| Three cores |python3|***1:24:22***|
|sample_Travis2_pcombined.fa|1008| Four cores  |python3|1:27:41|
|sample_Travis2_pcombined.fa|1008| Five cores  |python3|1:32:11|
|sample_Travis2_pcombined.fa|1008| Six  cores |python3|1:40:26|
|sample_Travis2_pcombined.fa|1008| Seven cores |python3|1:50:24|
|sample_Travis2_pcombined.fa|1008| Eight cores |python3|1:53:41|


