Scripts for processing short amplicon reads from Oxford Nanopore. For details, please see:

ST Calus, UZ Ijaz, and A Pinto. NanoAmpli-Seq: A workflow for amplicon sequencing from mixed microbial communities on the nanopore sequencing platform. bioRxiv 244517, 2018. DOI:10.1101/244517


Sam Young

dependencenies needed 

psutil 

biopython

Modified the program to run in python3 then add in multprocessing to decrease time overall. Ran several test to see what preformed the best 

Speed test conclude that 3 cores are leads to the greatest amount of time effiecent. Multithreading is not optimal for chopping up sequences it lead to the slowes amount of time. Calling the imports directly decreased time. 

linux terminal or mac terminals

If your default python is python3 other wise, this will work ./ otherwise

Help command calls 

./chopSEQ.py -h 

python3 chopSEQ.py -h 

This will allow the computer to run on three cores hopefully will maximum the timing

./chopSEQ.py -i input_file -f foward_primer -r reservser_primer -p 3 > output_file 

python3 chopSEQ.py -i input_file -f foward_primers -r reverse_primeRs -p 3 >   output_file 

this will allow the user to run on a single processor

./chopSEQ.py -i input_file -f foward_primer -r reverse_primer > output_file

python3 chopSEQ.py -i input_file -f foward_primer -r reversre_primer >  output_file 

Multithreading on all cores 

./chopSEQ.py -i input_file -f foward_primer -r reverse_primer -t > output_file

python3 chopSEQ.py -i input_file -f foward_primer -r reverse_primer -t > output_file 

i7-4970 CPU @ 3.60 Ghz on chopSEQ running times only using python3
------------------------------------------------------------------

Legend

best time ***1***

$$:!!:^^

$$ = hours,!! - minutes,^^= Seconds

| data sets | No. total seqs | proccesing used type | times | 
|--------------------------|----------------|----------------------|--------|
|sample_Travis2_p2.fa|419| Multithreading |1:49:06|
|sample_Travis2_p2.fa|419| Single core |1:19:54|
|sample_Travis2_p2.fa|419| Two cores |  57:05|
|sample_Travis2_p2.fa|419| Three cores |  ***50:08***|
|sample_Travis2_p2.fa|419| Four cores |  50:57|	
|sample_Travis2_p2.fa|419| Five cores |  52:56|
|sample_Travis2_p2.fa|419| Six cores  |1:00:16|
|sample_Travis2_p2.fa|419| Seven cores | 58:48|
|sample_Travis2_p2.fa|419| Eight cores |1:03:26|
|sample_Travis2_p4.1.fa|589| Multithreading |1:17:29|
|sample_Travis2_p4.1.fa|589| single core |  58:05|
|sample_Travis2_p4.1.fa|589| two cores |  35:18|
|sample_Travis2_p4.1.fa|589| three cores |  ***31:43***|
|sample_Travis2_p4.1.fa|589| four cores |  36:31|
|sample_Travis2_p4.1.fa|589| five cores |  40:23|
|sample_Travis2_p4.1.fa|589| six cores  |  43:27|
|sample_Travis2_p4.1.fa|589| seven cores |  42:36|
|sample_Travis2_p4.1.fa|589| Eight cores |  43:05|
|sample_Travis2_pcombined.fa|1008| Multithreading |3:11:07|
|sample_Travis2_pcombined.fa|1008| Single core |2:10:12|
|sample_Travis2_pcombined.fa|1008| Two cores |1:26:07|
|sample_Travis2_pcombined.fa|1008| Three cores |***1:24:22***|
|sample_Travis2_pcombined.fa|1008| Four cores  |1:27:41|
|sample_Travis2_pcombined.fa|1008| Five cores  |1:32:11|
|sample_Travis2_pcombined.fa|1008| Six  cores |1:40:26|
|sample_Travis2_pcombined.fa|1008| Seven cores |1:50:24|
|sample_Travis2_pcombined.fa|1008| Eight cores |1:53:41|


