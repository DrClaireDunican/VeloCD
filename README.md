# VeloCD

This is the GitHub repository for the VeloCD software tool developed by Dr Claire Dunican as a part of her PhD. This tool has adapted the RNA velocity algorithm (velocyto: http://velocyto.org/) to model and visualise transcriptomic trajectories from whole-blood RNA-Seq data. It can also be used to predict different current (diagnostic, classification) and future (prognostic) disease states. The development of this tool is described in *https://spiral.imperial.ac.uk/entities/publication/8a9c4a39-248c-4140-afcd-fa89804d18b8*. This software should be run using a minimum of three genes and two samples per dataset.

In brief, this method models the relationship between spliced and unspliced transcript expression to extrapolate future transcriptomic states under the assumption that unspliced transcript levels, characterised by the presence of introns, will provide information about future gene expression because they have not yet been processed. RNA velocity is a deviation from the assumption that gene expression change is not occurring (under a steady state). RNA velocity values are calculated as the observed unspliced transcript expression minus the unspliced transcript expression estimated under this steady state assumption.

VeloCD has two main outputs:

1. RNA velocity fate maps, which visually illustrates samples in low-dimensional transcriptomic space, with RNA velocity arrows representing future transcriptomic trajectories.

2. Probabilities of transition (transition probabilities) of each sample to each user-defined group (specified in the "Group" column of the metadata files, 2+ groups possible).

VeloCD has two main components. The first is the core VeloCD software tool, which is python based but also calls R functions and packages. It is run from the command line using  a .sh file (RunRNAVelocity.sh). 

The other component is a set of R and shell (.pbs) scripts that allow the user to go from raw fastq files to the input files required by VeloCD: spliced and unspliced transcript expression. These scripts have been designed to work on ubuntu computers (version 10+) but can alternatively be run via a ubuntu virtual box. These scripts are in the /PreProcessing sub-directory. Example files that can be input straight into this tool are in /ExampleData, alongside expected example output files. This analysis is expected to take approximately 5 minutes if run using the full range of possible hyperparamters.
## Dependencies
The following packages are required to run VeloCD, please ensure that you have 
python 3 and R installed. This software has been tested on R version 3.6.3 onwards (latest version tested: 4.4.2) and Python 3.9 onwards (latest version tested: 3.13). Install time is expected to take up to a few hours if all dependencies require manual installation. The following dependencies can be installed via the bash command line:


```bash
pip install numpy

pip install pandas

pip install matplotlib

sudo apt install imagemagick

pip install imagemagick

pip install IPython

sudo apt install re

pip install scipy

pip install scikit-learn

pip install numba

pip install umap

pip install gc-python-utils

pip install cmake

pip install umap-learn

pip install shutup

pip install loompy

pip install seaborn

pip install numpy_groupies

python3 -m pip install path.py

pip install Cython

pip install pysam

pip install velocyto

pip install typing

pip install h5py

python3 -m pip install pickle-mixin

sudo apt install zlib1g

pip install statsmodels

pip install rpy2

pip install pytest-gc

pip install shutup

sudo apt-get install libnlopt-dev

sudo apt-get install libxml2-dev

sudo apt-get install libcurl4-openssl-dev

sudo apt install libssl-dev

sudo apt -y install libfontconfig1-dev

sudo apt-get install libharfbuzz-dev libfribidi-dev libfontconfig1-dev

sudo apt-get install libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev 
```  
The following packages need to be installed in R:

```bash
install.packages("ggplot2")

install.packages("devtools")

packageurl <- "https://cran.r-project.org/src/contrib/Archive/nloptr/nloptr_1.2.1.tar.gz"

install.packages(packageurl, repos=NULL, type="source")

install.packages("reshape2")

install.packages("xml2")

install.packages("googledrive")
 
install.packages("googlesheets4")

install.packages("httr")

install.packages("rvest")

install.packages("systemfonts") 

install.packages("textshaping")

install.packages("tidyverse")

install.packages("gtools")

install.packages("Rmisc")

install.packages("dplyr")

install.packages("boot")

install.packages("gridExtra")

install_github("husson/FactoMineR")

install.packages("factoextra")

install.packages("scatterplot3d")
```
## Installation

This software can then be installed using the following code:
```bash
git clone https://github.com/DrClaireDunican/VeloCD.git
```


## Operating Instructions
VeloCD is run from the bash command line via:

```bash
chmod +x RunRNAVelocity.sh

./RunRNAVelocity.sh
```  
It is advisable that you read the following before attempting to run this program.

**Run Mode**

Running VeloCD will trigger a set of user prompts in the form of questions. The first will ask for the location of the input gene expression and metadata files. For each dataset you will need three to four CSV-format files: Spliced.csv (spliced transcript expression values, columns: samples, rows: genes, see /ExampleData), Unspliced.csv, Metadata.csv (see example in /ExampleData) and MetadataAlt.csv (optional). It will then ask for the full location of the source code of this algorithm i.e. the full local location of RNA_Velocity_Main.py. The Metadata.csv file will be used to define the groups used to colour the sample points on the fate map. These are also the groups used for prediction analysis. The optional MetadataAlt.csv file contains an addition layer of sample-level information that can be used to determine the "shape" of the sample points on the RNA velocity fate maps. If this second Metadata file is not provided then the first will be used to "shape" the sample points on the fate maps. There is also an later option to use this second sample grouping for prediction analysis.

VeloCD will then ask for the run mode. Briefly, VeloCD  can be run in two modes, M1 will run the analysis for all of the sets of files in the directory specified by the first question. These files should share a common suffix after "_" i.e. Spliced_Run1.csv, Unspliced_Run1.csv, Metadata_Run1.csv and MetadataAlt_Run1.csv. M2 will run the analysis for a single set of files, specified by the suffix requested in the proceeding question. The suffix for the above example is Run1.

**Analysis Questions**

The program will then ask a series of questions relating to the analysis itself. It first asks how many Principal Component Analysis (PCA) Principal Components (PCs) you would like to generate. The minimum number for analysis consisting of two genes is 2. It will then ask what low-dimensional embedding methods you'd like to use to generate the RNA velocity fate maps. The options are PCA, *t*-distributed Stochastic Neighbor Embedding (*t*-SNE: https://www.datacamp.com/tutorial/introduction-t-sne) and Uniform Manifold Approximation and Projection (UMAP: https://umap-learn.readthedocs.io/en/latest/). Please note that UMAP and *t*-SNE are stochastic methods and will generate slightly different results when re-run. In contrast, PCA is stable. If you say yes to PCA and have specified at least 3 PCs in the previous questions, the program will ask if you'd like 3-dimensional PCA-based fate maps. It will also ask this for UMAP embeddings. UMAP-specific requests including the selection of a distance metric and minimum distance are also given if this method is selected. For details about these parameters, please see:
https://umap-learn.readthedocs.io/en/latest/parameters.html

For the PCA-based fate maps this algorithm has the optional step of weighting the PCs based on their relative percentage explanation of the variation in the data with those with higher percentages favoured over lower. This step prioritises selecting close neighbours on the more highly weighted axis during the nearest neighbour search step (see *In thesis explanation*).

The program will then ask for you to specify some values or ranges for its hyperparameters. The first is if you would like to use a single or a range of embedding number of neighbour values (also called perplexity in the context of *t*-SNE, https://distill.pub/2016/misread-tsne/). It will then ask for this single value or minimum and maximum value (range). This hyperparameter is used to construct UMAP and *t*-SNE-based fate maps. The algorithm will also ask for the minimum, maximum (+1) and step-between values of the transition probability number of neighbours hyperparameter (see *thesis chapter for explanation*). This hyperparameter selects each samples closest neighbours to embed the RNA velocity arrows and calculate transition probabilities. The algorithm will then ask for a minimum threshold for prediction i.e. a probability of more than 0.5. The raw transition probabilities of each sample to each group are returned to the user so this threshold is just for the convenience of examining performance at a pre-specified threshold.

**Fate Map Visualisation Questions**

The next set of options relate to the visualisation of the RNA velocity fate maps. The first question will ask if you would like GIFS of the 3-dimensional plots to be generated. This is a relatively slow step of the algorithm. It will also ask if you wish to have the points shaped by the metadata feature in the first or second metadata file. Lastly, it will ask for a colour map from: https://matplotlib.org/stable/tutorials/colors/colormaps.html, to colour the sample-points on the fate map - based on their group in the first metadata file.

**Metadata Questions**

The algorithm will also ask if you have multiple metadata files and the type of features you have: strings or integers. Please note decimals should be treated as strings. It will also ask if you have only two feature levels in the metadata file "Group" column and if you would like sensitivity and specificity calculated for these two factor levels and which level should be used to calculate each value. Please note the second metadata feature does not need to be used for prediction and can just be used for visual purposes. Please ensure there are at least two unique groups in each metadata file.

**Result Output**

Once all questions have been answered, the program will begin running the RNA velocity and prediction analysis. A "Results" folder will then be generated in the same directory as the input files with output .csv "Files" and fate maps separated. The suffix of the file set will be used to separate the results of each dataset.

## Pre-processing Pipeline

This algorithm can be used to analyse any dataset given the data is in the form of 3-4 CSV-format files of the same general form as the examples given. The expression files should consist of sample columns and gene rows, with the first column "ID" containing the gene ID's. These should be in the same order between expression files. Spliced and unspliced transcript expression should be pre-normalised and log₂ (+1) transformed. Examples are given in /ExampleFiles. Code for these preprocessing steps are given in /Preprocessing.

The metadata files should consist of two columns, the first "Sample" has sample names matching the order of the columns in the expression set. The second column "Group" is the metadata feature you would like to predict.

The pipeline which generates spliced and unspliced expression from raw fastq files is given as a supplemental to this package and has three key steps:

1. Generate a genome reference file (for mapping and quantification), where genes are reduced into sets of non-overlapping transcripts. A pre-made version of this is given in /Annotation.

2. Sample Quality Control checks (via FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), mapping (STAR: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) and quantification (featureCounts: https://subread.sourceforge.net/featureCounts.html). Code to compile expression count files and mapping statistics is also included as well as other useful tools, such as for inferring strand specificity for mapped BAM files (infer_experiment.py, https://rseqc.sourceforge.net/).

3. Estimating spliced and unspliced transcript expression (these files are the input of the main VeloCD program). Reads overlapping an annotated intron that does not overlap any exons of other genes and are not labelled with the biotype "retained intron" are counted towards the genes unspliced transcript expression. Similarly, reads overlapping an annotated exon are counted towards the genes spliced transcript count. 

These pre-processing steps and scripts are given as a guide only and should be adapted to your data.

## Contact
Dr Claire Dunican - @ClaireDunican on Twitter/X - clairedunican@gmail.com
