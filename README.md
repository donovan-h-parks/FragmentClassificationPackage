Fragment Classification Package (FCP)
=====================================

Classifying genomic fragments from novel lineages using composition and homology.

## SYSTEM REQUIREMENTS

* Installation requires ~40GB of disk space.
* Classification of query fragments takes ~20GB per 1 millions fragments.
* Python v2.x is required to run the install script. Please note that we,
  have not tested installation using Python v3.x. Python 2.x comes 
  pre-installed on OS X and most Linux environments. Type 'python' from
  the command prompt to check if python is installed. Python, including
  a Windows installer, can be obtained from http://www.python.org/download/.
* If you wish to run the NB-BL, LCA, LCA+NB, or LCA+e-NB scripts you will
  need to have BLAST installed on your system. 

### INSTALLATION

To install the Fragment Classification Package (FCP), run the FCP_install.py script:

    > python FCP_install.py

This script downloads all completely sequenced bacterial and
archaeal genomes in the NCBI RefSeq database. Due to the size
of this file (> 7.5GB) it can take over an hour to download.
Taxonomic information is also obtained by downloading a portion
of the NCBI Taxonomy database.

Once these files are downloaded, the taxonomy of each contig in
the bacterial and archaeal genomes is determined. Sequence data
is than extracted from the Genbank files of each genome and a
Naive Bayes model built for each genome.

Installation takes around 2 hours, but does not require any 
user interaction. Conservative times for each step are as follows:

  1. Download NCBI RefSeq genomes: 3 hours
  2. Decompress genomes: 15 min
  3. Downloading NCBI taxonomy database: 1 min
  4. Decompress NCBI taxonomy database: 1 min
  5. Build taxonomy file for genomes: 1 min
  6. Create input sequence file for each genome: 3 min
  7. Build Naive Bayes models for each genomes: 30 min
  
If you wish to run the NB-BL, LCA, LCA+NB, or LCA+e-NB scripts,
a BLAST database of all genomes must be built. This can be done
by running the BuildBlastDB.py script:

    > python BuildBlastDB.py /path/to/makeblastdb ./training/sequences.txt

Building the BLAST database takes around 10 min.


### CLASSIFYING QUERY FRAGMENTS

Query fragments must be in a multi-FASTA file and each fragment must have 
a unique identifier. Identifiers consist of all text on the header line before 
the first space. After classifying query fragments with NB, BLASTN, NB-BL,
EPSILON-NB, LCA, or LCA+NB as described below, results can be summarized 
in a number of ways using the provided scripts (see SUMMARIZING CLASSIFICATION 
RESULTS).


### CLASSIFYING QUERY FRAGMENTS WITH NB

Classification with NB is done using the program nb-classify (nb-classify-windows.exe 
on Microsoft Windows systems). This program is run as follows:
    
    > ./nb-classify [options] -q <query-file> -m <model-file> -r <results-file>
  
Required parameters:
  <query-file>    Multi-FASTA file containing query fragments to classify.
  <model-file>    File indicating models to use for classification.
  <results-file>  File to write classification results to.

Optional parameters:
  --help        Print help message.
  --version     Print version information.
  --contact     Print contact information.
  -b <integer>  Number of fragments to classify at a time (default = 50000).
  -t <integer>  Log likelihood of the top T models will be returned. If you 
                  wish to have the log likelihood of all models in the
                  results file set T = 0 (default = 0).
  -v <integer>  Level of output information (default = 1).

Typical usage:
    
    > nb-classify -q test.fasta -m models.txt -r nb_results.txt
  
You must have a directory named 'nb-temp-results' below the path of nb-classify.


### CLASSIFYING QUERY FRAGMENTS WITH BLASTN

Classification with BLASTN is done using the script BLASTN.py as follows:
    
    > python BLASTN.py <blastn-path> <query-file> <results-file>

Required parameters:
  <blastn-path>   Full path to blastn.
  <query-file>    Multi-FASTA file containing query fragments to classify.
  <results-file>  File to write classification results to.

Typical usage: 
    
    > python BLASTN.py /path/to/blastn test.fasta blastn_results.txt


### CLASSIFYING QUERY FRAGMENTS WITH NB-BL

Classification with NB-BL requires the NB and BLASTN classifiers to first 
be run. NB must be run with the T parameter set to 0 so results are 
available for all models.

NB-BL is run using the NB-BL.py script:
    
    > python NB-BL.py <nb-results> <blastn-results> <results-file>

Required parameters:
  <nb-results>      Results of NB classifier with T=0.
  <blastn-results>  Results of BLASTN classifier.
  <results-file>    File to write classification results to.

Typical usage:

    > ./nb-classify -q test.fasta -m models.txt -r nb_results.txt
    > python BLASTN.py /path/to/blastn test.fasta blastn_results.txt
    > python NB-BL.py nb_results.txt blastn_results.txt nb-bl_results.txt


### CLASSIFYING QUERY FRAGMENTS WITH EPSILON-NB

Classification with epsilon-NB requires the NB classifiers to first 
be run. NB must be run with the T parameter set to 0 so results are 
available for all models.

Epsilon-NB is run using the Epsilon-NB.py script:
    
    > python Epsilon-NB.py <nb-results> <epsilon> <results-file>
  
Required parameters:
  <nb-results>    Results of NB classifier with T=0.
  <epsilon>       Use all models with a likelihood at most epsilon times smaller
                    than the maximum likelihood model to classify a fragment.
  <results-file>  File to write classification results to.

Typical usage:

    > ./nb-classify -q test.fasta -m models.txt -r nb_results.txt
    > python Epsilon-NB.py nb_results.txt 1E10 epsilon-nb_results.txt


### CLASSIFYING QUERY FRAGMENTS WITH LCA

Classification with LCA requires the BLASTN classifiers to first 
be run. LCA is run using the LCA.py script:

    > python LCA.py <blastn-results> <E-value> <percentage> <results-file>
  
Required parameters:
  <blastn-results>  Results of BLASTN classifier.
  <E-value>         Ignore hits with an E-value above this threshold.
  <percentage>      Use all hits with a bit score with this percentage of the 
                      hit with the highest bit score to classify a fragment.
  <results-file>    File to write classification results to.

Typical usage:

    > python BLASTN.py /path/to/blastn test.fasta blastn_results.txt
    > python LCA.py blastn_results.txt 1E-5 15 lca_results.txt


### CLASSIFYING QUERY FRAGMENTS WITH LCA + EPSILON-NB

Classification with LCA + Epsilon-NB requires the NB and BLASTN 
classifiers to first be run. NB must be run with the T parameter 
set to 0 so results are available for all models.

LCA + Epsilon-NB is run using the LCA+Epsilon-NB.py script:

    > python LCA+Epsilon-NB.py <blastn-results> <nb-results> <E-value> <percentage> <epsilon> <results-file>

Required parameters:
  <blastn-results>  Results of BLASTN classifier.
  <nb-results>      Results of NB classifier with T=0.
  <E-value>         Ignore hits with an E-value above this threshold.
  <percentage>      Use all hits with a bit score with this percentage of the 
                       hit with the highest bit score to classify a fragment.
  <epsilon>         Use all models with a likelihood at most epsilon times smaller
                       than the maximum likelihood model to classify a fragment.
  <results-file>    File to write classification results to.

Typical usage:

    > ./nb-classify -q test.fasta -m models.txt -r nb_results.txt
    > python BLASTN.py /path/to/blastn test.fasta blastn_results.txt
    > python LCA+Epsilon-NB.py blastn_results.txt nb_results.txt 1E-5 15 1E10 lca+epsilon-nb_results.txt

If you wish to run a LCA + NB classifier, simply set the epsilon parameter to 0. Similarly, a 
BLAST + Epsilon-NB classifier is obtained by setting the percentage threshold to 0.


### SUMMARIZING CLASSIFICATION RESULTS

The script TaxonomicSummary.py can be used to summarize the number of fragments and base pairs
assigned to different taxonomic categories: 

    > python TaxonomicSummary.py <query-file> <results-file> <summary-file>

Required parameters:
  <query-file>    Multi-FASTA file containing query fragments to classify.
  <results-file>  File indicating taxonomic classification of each fragment.
  <summary-file>  File where taxonomic summary information to.

Note that the classification file must come from the script Epsilon-NB.py, LCA.py, LCA+Epsilon-NB.py, 
or NB-BL.py. If you wish to summarize classifications from nb-classify or BLASTN.py you must first run
the results file through Epsilon-NB.py with epsilon set to 0 or LCA.py with the percentage threshold
set to 0, respectively.

Typical usage:

    > ./nb-classify -q test.fasta -m models.txt -r nb_results.txt
    > python Epsilon-NB.py nb_results.txt 0.0 nb_topModels.txt
    > python TaxonomicSummary.py test.fasta nb_topModels.txt nb_taxonomicSummary.txt


### ADDING NEW MODELS

New models can be added to those built during the FCP install using the
script AddModels.py (see below). If you wish to use the naive Bayes 
classified indepedent of FCP, see the example directory.

To build a new model, run the script AddModel.py:

    > python AddModel.py <N> <sequence-file> <domain> <phylum> <class> <order> <family> <genus> <species> <strain>

Required parameters:
  <N>              Desired n-mer length (must be the same for all models).
  <sequence-file>  Multi-FASTA file containing sequence data to build model from.
  <domain>         Domain of model.
  ...
  <stain>          Strain of model.

Typical usage:

    > python AddModel.py ./training/custom/MySeqData.fasta Bacteria Proteobacteria 
              Betaproteobacteria Burkholderiales Burkholderiaceae Ralstonia 
              "Ralstonia pickettii" "Ralstonia pickettii 12D"

(should be specified on a single line)

Note that this will modify the files sequences.txt, taxonomy.txt and models.txt to include
your new model. You may wish to back these up before adding custom models. Alternatively,
you can edit the file models.txt to remove custom models from consideration.

If you wish new models to be considered by BLASTN you must run BuildBlastDB.py:

    > python BuildBlastDB.py /path/to/makeblastdb ../training/sequences.txt

This only needs to be run once after adding in a collection of new models.


### CHANGING N-MER LENGTH
 
The length of n-mers being used can be changed using the nb-train executable
in the nb-train directory (nb-train-windows.exe on Microsoft Windows systems):

    > ./nb-train [options] -s <sequence-file> -m <model-dir>

Required parameters:
  <sequence-file>  File listing path to each FASTA file for which a model shoud be built.
  <model-dir>      Directory to store models.

Optional parameters:
  --help        Print help message.
  --version     Print version information.
  --contact     Print contact information.
  -n <integer>  Desired oligonucleotide length (default = 10).

Typical usage:

    > ./nb-train -s sequences.txt -m ./models/


### HOW TO PARALLELIZE CLASSIFICATION

If you are classifying many millions of fragments, you may wish to parallelize the NB
classification. This is easily done by dividing the query file into several files with
approximately the same number of sequences. nb-classify can be applied to each of these
files separately. The classification file for each of these runs should then be combined
before using any of the Python scripts. Take care to include the header file only once 
when combining the files. Contact us for further help.


### REPORTING BUGS OR FEATURE SUGGESTIONS

This software is in active development and we are interested in discussing all potential 
applications. We encourage you to send us suggestions for new features and to report bugs.
Suggestions and bug reports can be sent to Rob Beiko (beiko [at] cs.dal.ca). If reporting a 
bug, please provide as much information as possible and a simplified version of the data 
set which causes the bug. This will allow us to quickly resolve the issue.
