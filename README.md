# HOME  
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
git clone https://github.com/ListerLab/HOME.git
cd ./HOME
pip install -r requirements.txt
python setup.py install
```
***For conda users**
```
git clone https://github.com/ListerLab/HOME.git
cd ./HOME
conda env create    *assuming the conda environment is activated and R is already installed in it*
source activate HOMEenv
python setup.py install
```
# Usage
HOME can be run in pairwise mode for two group comparisons and time series mode for more than two group comparisons. It can also be used for mutiple pairwise comparisions with large number of input samples. 

**BSSeeker2 CGmap file can be provided as input file directly**
**or**
**Input file format as mentioned below needs to provided**

Chromosome number, position, strand, type (CG/CHG/CHH) where H is anything but G, methylated reads and total number of reads. For a sample, this information are saved in a single tab separated text file without header, which can be compressed or uncompressed. Below shows an example of such file:

```
chr1	15814	+	CG	12	14
chr1	15815	-	CG	15	21
chr1	15816	-	CHG	1	9
chr1	15821	-	CHH	7	22
chr1	15823	-	CHH	0	2
chr1	15825	-	CHH	11	19
```
**DMR detection for two group comparison: Pairwise**

```
HOME-pairwise 	-t [CG/CHG/CHH/CHN/CNN]	 -i [sample_file_fullpath] 	-o [output_directorypath]
```

*Note: Please check the number of cores to use and set them by npp parameter (default is 8). Also for non-CG DMR prediction for huge genomes like mammalian genome use parameter **-sin**.*

Example: 
```
HOME-pairwise 	-t CG 	-i ./testcase/sample_file_CG.tsv 	-o ./outputpath 

or (if CGmap files from BSSeeker2)

HOME-pairwise 	-t CG 	-i ./testcase/sample_file_CG.tsv 	-o ./outputpath --BSSeeker2

```

Required arguments:
```
 -t --type 	        Type of DMRs (CG /CHH/CHG/CHN/CNN) 

 -i --samplefilepath Sample file containing sample names and sample paths for each replicate (TAB sep); each sample info should be in different rows
 
 -o --outputpath 	  Path to the output directory  
```
The default parameters for HOME are **relatively permissive**. To run HOME with more stringent setting please change the defaults parameters as below or higher:
```
HOME-pairwise 	-t CG 	-i ./testcase/sample_file_CG.tsv -o ./outputpath --delta 0.2 --minc 5 
```
Optional arguments: 

| Parameter                 | Default       | Description   |	
| :------------------------ |:-------------:| :-------------|
| -sc --scorecutoff 	       |	0.1           |the score from the classifier for each C position
| -p  --pruncutoff          | 0.1           |the SVM score checked for consecutive Cs from both ends to refine the boundaries
| -npp -–numprocess 	       |	8	            |number of cores to be used
| -ml --minlength  		       | 50	           | minimum length of DMRs required to be reported 
| -ncb --numcb 		           | 5             | minimum number of Cs present between DMRs to keep them seperate
| -md  -–mergedist 	        | 500           | maximum distance allowed between DMRs to merge 
| -prn --prunningC	         | 3             | number of consecutives Cs to be considered for pruning for boundary refinement2
| -ns --numsamples          | all           | no.of samples to use for DMR calling; default takes all sample in the file
| -sp --startposition       | 1st position  | start position of sample in the sample file to use for timeseries DMR calling 
| -BSSeeker2 --BSSeeker2    | False         | input CGmap file from BSSeeker2
| -mc --minc			             | 3 	           | minimum number of Cs in a DMR
| -sin --singlechrom			     | False         | parallel code for single chromosome; *npp* will be used for parallel run for each chr
| -d --delta			             | 0.1     	     | minimum average difference in methylation required in a DMR 

**Parameter –sc**

HOME assigns a score from the classifier to each C position based on the histogram features. This score is then used to cluster the C’s and call the DMRs. The user can set a lower score cutoff if the DMRs seems to be missing or DMR boundaries seems to be missing the differential methylated cytosine near the boundaries. Similarly, the user can increase the cutoff if the boundaries of DMRs seems to be extended. The score ranges from 0-1 and the default is set to 0.1. 

*NOTE: the default score is set after rigorous testing on various data sets and need not to be varied in most of the cases. The default should only be changed after proper visualization of the DMRs on the browser.*     

**Parameter –p**

This parameter controls boundary accuracy of the DMRs. After the DMRs are called they are then pruned based on the scores for the consecutive C’s on both ends. The user can increase the cutoff if the boundaries seems to be extending too far. Similarly the user can lower the cutoff if the DMR boundaries seems to be missing the differential methylated C’s. The range is from 0 to 0.5 and the default is 0.1. Please note that similar to parameter –sc this parameter need not to be altered in most of the cases and should only be changed after proper visual inspection of the DMRs.

 **Parameter –ml**

This parameters sets minimum length (in terms of base pairs) for a DMR to be reported. 
The DMRs below this length will be skipped and not reported in the filtered DMR output file. The default is 50bp.

**Parameter –ncb**

This parameter controls when the smaller DMRs should be merged into one. It controls minimum number of C’s required between DMRs to keep them as separate DMRs.  It works in relation with parameter –md (described below). The default is 5 C’s. So, if the number of C’s are less than 5 and the distance is less than 500bp (default for --md) between two consecutive DMRs it will be merged into one single DMR.

**Parameter –md**

This parameter allows the user to set the merge distance between two consecutive DMRs. The defaults is 500bp. 

**Parameter –npp**

This parameter allows the user to set the number of parallel process to run at a time. The default is 8. 

**Parameter –mc**

 This parameter allows the user to set minimum number of C’s present in a DMRs to be reported. Any DMR with less than the set value will not be reported in the filtered DMR file. The default is 5.

**Parameter –d**

This parameter sets minimum average methylation difference present for a DMR to be reported. The DMRs with less than the set value will not be reported in the filtered DMR file. The default is 0.1. 

**Parameter –prn**

This parameter is used in relation to parameter –p (described above). This controls the number of consecutive C’s to be considered from both ends for boundary refinement. The default is 3. Alteration of this parameter should only be done after proper visual inspection of the DMRs. 

**Parameter –sin**

This parameter is used if you want to parallel the code by single chromosome. The default is False, so the code will be parallel for all chromosomes. It should be used with huge chromosomes for example in case of non-CG DMR prediction for mammalian genome. If the genome size is small it is adviced not to use it. To turn it on just say *-sin* in the command line. 

**Parameter –ns**

This parameter is used if you want to use selected number of samples from your sample input file. The default is False, so the code will use all the samples in the sample input file. It allows you to have as many samples as you want in your input file but control the number of samples to use for DMR calling. 

**Parameter –sp**

This parameter is used if you want to select the samples from anywhere in your sample input file. This parameter is used along with **ns** parameter. The default is False, so the code will start from the 1st sample in the sample input file. It allows you to have as many samples as you want in your input file but control the samples to use for DMR calling. 

**Parameter –BSSeeker2**

This parameter is used if you want to provide CGmap file directly. The default is False, so the code will require the files in the input format mentioned above. If the user have methyaltion output files from BSSeeker2, it can be provided directly. To turn it on just say *-BSSeeker2* in the command line. 

**Output format**

The output format is:

```
chr	start	end	status	numC	mean_meth1	mean_meth2	delta	Avg_coverage1	Avg_coverage2	len
```

Here, status refers to state of DMR (hyper/hypo). Mean_meth1 and mean_meth2 refers to mean methylation level for sample1 and 2 respectively. Delta refers to the difference in mean methylation level for two samples. Avg_coverage1 and avg_covereage2 gives the mean coverage for both samples. 


**DMR detection for more than two groups: time series**

```
HOME-timeseries 	-t [CG/CHG/CHH/CHN/CNN]	-i [sample_file_fullpath]		-o [output_directorypath]
```
Example: 

```
HOME-timeseries 	-t CG 	-i ./testcase/sample_file_CG.tsv	–o /outputpath
```
Required arguments:

```
-t  --type	            type of DMRs (CG /CHH/CHG/CHN/CNN) 

-i --samplefilepath Sample file containing sample names and sample paths for each replicate (TAB sep); each sample info should be in different rows 

-o  –-outputpath 		    path to the output directory  
```

Optional arguments: 


| Parameter                 | Default       | Description   |	
| :-------------------------|:-------------:| :-------------|
| -sc --scorecutoff 	       |	0.5			        | the score from the classifier for each C position 
| -npp -–numprocess 	       |	5	            | number of cores to be used
| -ml --minlength  		       | 50	           | minimum length of DMRs required to be reported 
| -ns --numsamples          | all           | no.of samples to use for DMR calling; default takes all sample in the file
| -sp --startposition       | 1st position  | start position of sample in the sample file to use for timeseries DMR calling 
| -BSSeeker2 --BSSeeker2    | False         | input CGmap file from BSSeeker2
| -mc --minc			             | 4 	           | minimum number of Cs in a DMR
| -d --delta			             | 0.1     	     | minimum average difference in methylation required in a DMR
| -sin --singlechrom			     | False         | parallel code for single chromosome; *npp* will be used for parallel run for each chr


**Output format**

```
chr	start	end	numC	len	max_delta	confidence_scores	Comb1-n

```
Here, Max_delta is the maximum average methylation difference among the compared samples. 
Confidence score takes into account the length, number of C’s and SVM score. The higher value denotes more confident DMR. 
Comb1-n denotes the pairwise comparisons for each combination of samples. It reports *start: end: state: delta* for each pairwise comparison. 
 

# Required tools

[python 2.7](https://www.python.org/download/releases/2.7/)

[numpy](https://pypi.python.org/pypi/numpy)

[pandas v0.17.1](https://pypi.python.org/pypi/pandas/0.17.1/)

[scipy v0.16.0](https://pypi.python.org/pypi/scipy/0.16.0)

[scikit-learn v0.16.1](https://pypi.python.org/pypi/scikit-learn/0.16.1)

[statsmodels v0.6.1](https://pypi.python.org/pypi/statsmodels/0.6.1)


# Troubleshooting

To stop HOME execution in middle:
```
cntrl+c and then
cntrl+z
```
*Note: always delete the directories created by HOME run if stopped in middle*

Error *Exception: File...training_data file does not exist*

```
Please remember to CD into HOME first before starting the run 
```

# Citation

If you use this software in your work, please cite our paper [HOME](https://www.biorxiv.org/content/early/2017/12/02/228221)
