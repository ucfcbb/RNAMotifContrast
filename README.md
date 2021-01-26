## RNAMotifContrast

**A method to visualize the RNA structural motif variations and discover subfamilies**

* Shahidul Islam<sup>†</sup>, shahidul at knights dot ucf dot edu
* Md Mahfuzur Rahaman<sup>†</sup>, mahfuzur at knights dot ucf dot edu
* Shaojie Zhang*<sup>†</sup>, shzhang at cs dot ucf dot edu

<sup>†</sup>Department Computer Science, University of Central Florida, Orlando, FL 32816 USA \
*To whom correspondence should be addressed.

---

RNAMotifContrast is developed in a **64-bit Linux** machine and tested on multiple Linux-based workstations and servers. It is also tested on **64-bit MacOS** environment. For basic features of generating subfamilies as text outputs, only python (3.x recommended) is required to be installed on the computer with the previously mentioned environments. To generate corresponding images, it additionally requires PyMOL to be installed. Here, we are presenting the guidelines to use this tool on a **Linux** or **MacOS** based systems. **Please note that to successfully generate the text and image outputs from RNAMotifContrst, a 64-bit machine with Linux or MacOS environment is required.**

RNAMotifContrast uses merged annotation by default that combines both FR3D and DSSR annotations to get the structural information of an RNA. Although FR3D annotations are available online and DSSR annotations can be generated using the x3dna-dssr tool (single-user license for academic and non-profit use), we have included all the required annotation files with our project. For any new PDB files, the tool will generally work with FR3D annotations collected from the online data source. To incorporate the DSSR annotations for the new structures, user might choose to generate them and put the annotation files in [RNAMotifContrast/data/annotation/dssr/](/data/annotation/dssr/) directory. A user also has the option to use other annotation tools, but in that case, they need to ensure that the annotation format corresponds to our merged annotation format.

### 1. Installation

#### 1.1 Install prelimineries

All recent Linux and MacOS systems normally come with python installed by default. If not, please try the following commands to install `python` and `pip` (must be run as root):

```
# Debian/Ubuntu/Mint (Using apt-get package manager)
apt-get install python
apt-get install python-pip

# CentOS/Fedora  (Using yum package manager)
yum install python
yum install python-pip

# MacOS  (Using Homebrew package manager)
brew install python
```

If python is already installed but there is no `pip`, **both Linux and MacOS** users can also follow the instructions provided in https://pip.pypa.io/en/stable/installing/ to install `pip` only.

**Note: Based on the default settings of an operating system, `python` and `pip` might refer to either Python 2.x or Python 3.x versions. A user might need to use `python2` and `pip2`, or `python3` and `pip3` to use the other Python version.**

It is also possible to install python and pip along with other basic libraries by installing Anaconda. To do so, please download OS-specific **Anaconda** from their [website](https://www.anaconda.com/) and install by following the instruction provided there. It will not require any other additional commands to install `python` and `pip`.

#### 1.2: Install required Python libraries

It is required to install several python libraries to run RNAMotifContrast. These libraries are included in the [requirements.txt](requirements.txt) file. To install all required python libraries, please navigate to the RNAMotifContrast home directory in the terminal and execute the following command.

```
pip install -r requirements.txt
```

### 2. Input Specifications

RNAMotifContrast takes input from a ‘*.in’ file which needs to be in the [RNAMotifContrast/data/](data) directory or any subdirectory inside. Each line in that file represents a motif family. The motif family starts with a name, followed by a comma-separated list of motifs (the indices for motifs are expected to be in the PBD index, but it can be changed to Fasta index by setting a parameter in the configuration file). To see examples of formats, please check the two sample input files ([sample1.in](data/sample1.in) and [sample2.in](data/sample2.in)) provided in the [data](data) directory.

### 3. Commands for usage

Note: **MacOS users** might get an error message saying `'align_ga.mac' cannot be opened because it is from an unidentified developer`. To get rid of this error, please navigate to: 'System Preferences > Security & Privacy > General' and set 'Allow apps downloaded from' to 'Anywhere'.

```
python run.py               [Get text outputs for default test input file]
python run.py -i <input_file_name>  [Get text outputs for user input file]
python run.py -o <subdir_name>          [Get text output to user defined subdir]
python run.py -p            [Get PyMOL image outputs along with the text outputs]
python run.py -v            [Provide an orientation for the each motif family]
python run.py -h            [See more parameter options]
```

**Example:**
To generate subfamily result in text format from [sample1.in](data/sample1.in) to the default ‘output’ directory:

```
python run.py
```

To generate subfamily result from [sample1.in](data/sample1.in) to ‘sample1’ subdirectory in ‘output’ directory:

```
python run.py -o sample1
```

To generate subfamily result from [sample2.in](data/sample2.in) to ‘sample2’ subdirectory in ‘output’ directory:

```
python run.py -i sample2.in -o sample2
```

### 4. Description for output in text format

RNAMotifContrast gives the outputs in text format by default which includes the following results in the corresponding directories:

* Parent-child relation of traversal, alignment details (length, score, z-score), and  RMSD of superimposition (‘RNAMotifContrast/output/superimposition_details/’)
* Subfamily representative motifs with base pair and stack annotations (‘RNAMotifContrast/output/subfamily_representatives/’)
* Annotations of all motifs for each subfamily (‘RNAMotifContrast/output/subfamily_details/’)
* Subfamily summary data (corresponding to the result table in the article) and family-wise alignment length threshold (‘RNAMotifContrast/output/summary/’)

### 5. PyMOL Installation (required to generate images)

PyMOL can be installed directly by downloading the OS-specific version from https://pymol.org/. It is also possible to install from the package managers using the following commands (must be run as root):

```
# Debian/Ubuntu/Mint (Using apt-get package manager)
apt-get install pymol

# CentOS/Fedora  (Using yum package manager)
yum install pymol

# MacOS  (Using Homebrew package manager)
brew install brewsci/bio/pymol

# Anaconda (Both Linux and MacOS)
conda install -c schrodinger pymol
```

If PyMOL `open-source` version is installed, a user does not have to provide the license file. But, for the `Schrödinger` version, users have to provide a specific license file. A free `Educational-Use-Only` license can be obtained from [here](https://pymol.org/edu/?q=educational).

Please note that the installed PyMOL will work only with the specific python version to which it is linked with. For example, if the installed PyMOL is linked with Python 3.x, it'll only work with Python 3.x, not with Python 2.x. In most cases, the installed PyMOL links with Python 3.x version. If for some reason, the installed PyMOL version does not link with any Python version, a user might consider getting it by compiling the source (Python 3.6+ is required). The steps to compile the PyMOL source is described in the [README-PyMOL-build](README-PyMOL-build.md) file.

### 6. Image output (PyMOL required)

There is an option to generate image outputs from RNAMotifContrast to visualize the traversal and superimposition of the motifs. To do so, the ‘-p’ parameter needs to be used.

**Example:**

```
python run.py -o sample1 -p
```

The images for subfamily results includes the following output in the corresponding directories:

* A combined image for each motif family consisting of input motifs in default orientation, side-by-side ordered and rotated motifs according to the parent-child relation of traversal, and superimposed motifs for subfamilies and families (‘RNAMotifContrast/output/subfamilies/’)
* Color annotated images of representative motifs for each subfamily, following the coloring standard defined in the article (‘RNAMotifContrast/output/subfamily_representatives/’)
* Individual images of side-by-side ordered and rotated motifs for each subfamily in separate files (‘RNAMotifContrast/output/subfamily_details/’)
* Images to show the progression of superimposition (‘RNAMotifContrast/output/progressive’)

### 7. Set orientations for the first motif of traversal of each motif family (PyMOL required)

There is an option to set the orientation of the first motif in the traversal for each motif family, which eventually guides the orientation of all motifs through rotations. For this task, the ‘-v’ parameter needs to be used.

**Example:**
To set the orientation for the default input,

```
python run.py -v
```

After this command, the PyMOL window will be popped up with the first motif. In that window, you can rotate the motif to set an orientation that is suitable to view its structural features. To save the orientation, you need to keep the PyMOL window open and type ‘Y’ in the command prompt of the system terminal (not in the PyMOL). If you type ‘N’, it will discard the changes in the orientation and keep the previously saved orientation. This process will repeat for all the motif families in the input file.

### 8. Using pre-generated alignment files

To run the process faster, we recommend using pre-generated alignment files when possible. For the four input files, the alignments are available on the RNAMotifContrast homepage. Please download and extract them to the [data](data) directory. Then provide the extracted subdirectory path as ‘-d’ parameter value while executing [run.py](run.py).

**Example:**
For [sample1.in](data/sample1.in), if alignment files are extracted to ‘RNAMotifContrast/data/sample1_alignments’ directory, then execute,

```
python run.py -d sample1_alignments
```

### 9. Configuration file

User can edit the [config.py](config.py) file provided in the RNAMotifContrast home directory, to change values of various parameters, such as, annotation source (DSSR/FR3D/Merged), traversal algorithm, merging of motifs (connectivity test threshold along with RMSD). We recommend making changes only in the ‘USER PARAMS (ADVANCED)’ section of this file.

### 10. RNA details (type, organism)

In the text files for superimposition details, extra three columns will be added to show RNA types and organism information for the chains. We have included this information for all the RNA chains of our supercluster dataset. If the input file contains a subset of those chains, this information will be automatically added to the output.

Please check [this](EstimatedRuntime.md) file to know the estimated runtime for the provided sample inputs.

### ACKNOWLEDGEMENTS

RNAMotifContrast is developed for an NIH funded project (R01GM102515).
  
### CONTACTS

For bug reports or comments please contact shzhang@cs.ucf.edu
