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
