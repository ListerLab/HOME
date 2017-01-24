#HOME  
HOME (histogram of methylation) is a python package for differential methylation region (DMR) identification. 
The method uses histogram of methylation features and the linear Support Vector Machine (SVM) to identify DMRs
from whole genome bisulfite sequencing (WGBS) data. HOME can identify both pairwise and time series DMRs with or without replicates.
#Installation
HOME is written for python 2.7 and tested on Linux system. It is recommended to set up virtual environment for python 2.7 first before installing HOME package.

**Step 1:** Create a virtual environment for HOME

```
virtualenv -p <path_to_python2.7> <env_name>
```

*NOTE : the above command assumes that the virtualenv tool is already installed on your system.*

**Step 2:** Activate the virtual environment

```
source <env_name>/bin/activate
```
*The name of the activated environment will appear on the left of the prompt.* 

**Step 3:** Install the HOME package outside virtual environment but make sure that the virtual environment is active

```
git clone https://github.com/Akanksha2511/HOME.git
cd ./HOME
pip install -r requirements.txt
python setup.py install
```
#Usage
HOME can be run in pairwise mode for two group comparisons and time series mode for more the two group comparisons. 

**Input file format**

Chromosome number, position, strand, type (CG/CHG/CHH) where H is anything but G, methylated reads and total number of reads. For a sample, this information are saved in a single tab separated text file without header, which can be compressed or uncompressed. Below shows an example of such file:

```
chr1	15814	+	CG	12	14
chr1	15815	-	CG	15	21
chr1	15816	-	CHG	1	9
chr1	15821	-	CHH	7	22
chr1	15823	-	CHH	0	2
chr1	15825	-	CHH	11	19
```
**DMR detection for two group comparison:**

```HOME-pairwise 	-t [CG/CHG/CHH/CHN/CNN]	 -a [path of sample1 (including filename)]	 -b [path of case (including filename)] 	-o [output path]
```

*Note: In case of replicates for each sample separate them by space.*

Example: 
```
HOME-pairwise 	-t CG 	-a /testcase/sample1_rep1.txt  /testcase/sample1_rep2.txt 	-b /testcase/sample2_rep1.txt  /testcase/sample2_rep2.txt 	-o /output 
```
Required arguments:
```
-t 	Type of DMRs (CG /CHH/CHG/CHN/CNN) 

-a 	path of sample1 (replicates should be separated by space) 

-b 	path of sample2 (replicates should be separated by space) 

-o 	path to the output directory  
```

Optional arguments: 
```
Parameter			      default				description	
-sc --SVMscorecutoff 		  0.1 			the score from the classifier for each C position 
-p --pruncutoff           0.1 	    the SVM score checked for consecutive C’s from both ends to refine the boundaries.
-ml --minlength  		      50	      minimum length of DMRs required to be reported 
-ncb --numcb 			        5	        minimum number of C’s present between DMRs to keep them seperate
-md –mergedist 		        500	      maximum distance allowed between DMRs to merge 
-npp –numprocess 		      5       	number of cores to be used 
-mc—minc			            5 	      minimum number of C’s in a DMR
-d—delta			            0.1     	minimum average difference in methylation required in a DMR 
-prn--prunningC			    	3 	      number of consecutives C’s to be considered for pruning for boundary refinement
```
