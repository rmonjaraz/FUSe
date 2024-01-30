# FUSe
Align, Trim and **F**ilter **U**CE **Se**quences and Alignments

---
This script is an automatized workflow for aligning, triming and filtering alignments, ideally obtained from UCE data, but can be used on any type of alignments. 

* Alignment is performed in [MAFFT][1] as implemented in `phyluce_align_seqcap_align` from [Phyluce][2]
* Trimming of alignments is performed with [trimAL][3] or [Gblocks][4].
* Filtering of sequences and alignments is as follows:
    * __Sequences__:
    Remove Short - Removes short sequences based on a percentage of gaps acording with entire lenght of the 
    Remove Divergent - Removes divergent sequences based on a pairwise identity comparisson
    * __Alignments__:
    Number of taxa - Remove alignments that doesn't meet the minimum number of taxa desired
    Alignment lenght - Remove alignments that doesn't meet the minimum lenght desired (in bp)
    Completeness - Generates completenes matricess based on a percentage of taxa (Similar to min_taxa but based on total number of taxa)

## Citation
---
If you use this program please cite [Calisoga Paper][5] 
> Calisoga Paper

and [Phyluce][2]:
> Faircloth BC. 2015. PHYLUCE is a software package for the analysis of conserved genomic loci.  Bioinformatics. doi: [10.1093/bioinformatics/btv646](https://doi.org/10.1093/bioinformatics/btv646).

## Installation
---
FUSe requires for [phyluce][2] to be installed in your system, if phyluce is properly working, FUSe should be able to run without installing any further dependencies. Download or clone this repository and use FUSe.py as a stand-alone program. Alternatevely, you can copy FUSe.py in your phyluce bin folder located in: ```/User/miniconda#/envs/phyluce-X.x.x/bin/``` where '#' is your conda version and "X.x.x" is your version of phyluce. This way, you can call FUSe directly from your phyluce environment.

```
cp FUSe.py /User/miniconda#/envs/phyluce-X.x.x/bin/
chmod 775 /User/miniconda#/envs/phyluce-X.x.x/bin/FUSe.py
```
For installing phyluce, detailed intructions can be found [here][6]

## Execution
---
Command line interface, to deploy the menu and options available type:
```python FUSe.py -h```
```
usage: FUSe [-h] -i MONOLITHIC_FILE -t TAXA [-p PREFIX]
            [-o {fasta,nexus,phylip,clustal,emboss}] [-c CORES] [--trimAL]
            [-a {automated1,gappyout,strict,strictplus,nogaps}] [--gblocks]
            [--b1 B1] [--b2 B2] [--b3 B3] [--b4 B4] [--remove-div]
            [-d DIVERGENT] [--remove-short] [-s SHORT_CUTOFF]
            [--filter-alignments] [-l MIN_LENGTH] [-m MIN_TAXA]
            [--get-completeness] [-e PERCENT] [--taxa-count] [-v]

Align, Trimm and Filter UCE Sequences and Alignments.

optional arguments:
  -h, --help            show this help message and exit
  -p PREFIX, --prefix PREFIX
                        The name to be used for all output folders (default:
                        OUTPUT).
  -o {fasta,nexus,phylip,clustal,emboss}, --out-format {fasta,nexus,phylip,clustal,emboss}
                        The ouput alignment format (default: fasta).
  -c CORES, --cores CORES
                        The number of PHYSICAL CPUs (default: 1).
  --trimAL              Wether to trim alignments using trimal.
  -a {automated1,gappyout,strict,strictplus,nogaps}, --t-method {automated1,gappyout,strict,strictplus,nogaps}
                        trimAl automated method to use (default: automated1).
  --gblocks             Wether to trim alignments using gblocks.
  --b1 B1               The GBLOCKS -b1 proportion (default: 0.5)
  --b2 B2               The GBLOCKS -b2 proportion (default: 0.70)
  --b3 B3               The GBLOCKS -b3 integer value (default: 10)
  --b4 B4               The GBLOCKS -b4 integer value (default: 4)
  --remove-div          Wether to remove divergent sequences from alignments.
  -d DIVERGENT, --divergent DIVERGENT
                        Percentage of pairwise identity in every sequence to
                        retain (default: 0.7).
  --remove-short        Wether to remove short sequences from alignments.
  -s SHORT_CUTOFF, --short-cutoff SHORT_CUTOFF
                        Percentage of "-" (gaps) in every sequece to retain
                        (default: 0.7).
  --filter-alignments   Wether to filter alignments by No. taxa and length in
                        bp.
  -l MIN_LENGTH, --min-length MIN_LENGTH
                        The minimum alignment lenght to retain alignment
                        (default: 50).
  -m MIN_TAXA, --min-taxa MIN_TAXA
                        The minimum number of taxa in alignments to retain
                        (default: 4).
  --get-completeness    Wether to filter alignments by completeness.
  -e PERCENT, --percent PERCENT
                        Completeness matrix percentage to output (default:
                        0.8).
  --taxa-count          Wether to print taxa count in alignments summary.
  -v, --version         show program's version number and exit

Required arguments:
  -i MONOLITHIC_FILE, --input MONOLITHIC_FILE
                        The input monolithic fasta file.
  -t TAXA, --taxa TAXA  The total number of taxa in all alignments.
```

#### Input File:
FUSe takes as input a single **Monolithic** [FASTA][7] file, containing all your sequences for all your loci; when using phyluce for processing UCEs this is usually the resulted file from ```phyluce_assembly_get_fastas_from_match_counts```.

Alternatevely, if you don't have a monolithic FASTA file but instead individual loci files that needs to be processed with FUSe we can use the next snippet to get a monolithic file:

We will use ```phyluce_align_explode_alignments``` to explode files and get individual fasta files per sample containing all loci asociated to that particular sample, assuming that your loci files are stored in a folder named "myFiles", we will write new files to "myExploded" folder:
```
phyluce_align_explode_alignments \
--alignments myFiles/ \
--output myExploded \
--input-format fasta \
--by-taxon \
--include-locus
```
NOTE: Make sure to include flag ```--include-locus``` in order to associate sequences with each locus file, otherwise this will combine all the loci in future steps.

The content of "myExploded" directory should be multiple fasta files, one for each sample in our loci files, next we just need to concatenate all samples in a single monolithic file, in this case named "myMonolithic.fasta":
```
cat myExploded/*.fasta > myMonolithic.fasta
```
This monolithich file can be used as input for FUSe

#### Usage:

##### Required arguments: 
The following are the only required arguments for running the pipeline, if the user provides only those parameters the entire pipeline will run with default parameters (Check Filtering Options section for more details).
+ **Monolithic Fasta**
+ **TAXA** (Number of total taxa in your alignments)

##### Workflow:
1. Aligning
2. Trimming
3. Remove Divergent Sequences
4. Remove Short Sequences
5. Filter Sequences
6. Completeness Matrix

The order of the workflow is as listed above, the pipline will always start by attepting to align, followed by trimming and filtering; If the pipeline is run using default parameters the final output of FUSe will be Multiple Sequence Alignments (MSAs). Different steps of the workflow can be avoided but the remaning steps will be performed always in the same order, for example if we donot desired to trim alignments, the pipeline will start aligning followed by removing divergent sequencesand so on. 

NOTE: Current version of FUSe doesn't support inputting previuously aligned MSAs for performing filtering or other operations.

In order to activate the different parts of the workflow we need to provide with the desired flags (see more detailed Examples below), for example if we want to align and filter sequences, we will have to provide the flag ```--filter-alignments```.

##### Options description:

###### Aligning
Current version of FUSe only supports aligning with [MAFFT][1].

###### Trimming
FUSe supports two types of trimming:

| Trimming | Options | Default | Reference |
| :-----------: | ----------- | -------- | -------- | 
| `--trimAL` | `-a` automated1, gappyout, strict, strictplus, nogaps | `-a automated1` | [trimAL][3] |
| `--gblocks` | `-b1` B1, `-b2` B2, `-b3` B3, `-b4` B4 | `-b1 0.5`, `-b2 0.7`, `-b3 10`, `-b4 4` | [Gblocks][4] |

NOTE: Some recent versions of phyluce doesn't support gblocks, this seems to be an issue with new operating systems, if you have a problem runing gblocks try running the pipeline using trimAL.

###### Sequence Filtering
---

###### Remove Divergent Sequences
This option is intended to remove outlier sequences resulted from poor alignment, this could be due to sequences of poor quality, contamination or in some cases due to paralogy. This option is implemented in a similar approach as in the software [CIAlign][8]. It calculates pairwise sequence identity across all samples of the same alignment and uses an user defined threshold percentage for removing samples. Defuault parameter (70%) will remove any sequence with a pairwise identity higher than 70, i.e will retain any sequence that is 70% or less divergent that the majority of the alignment. thershold can be modify using the flag ```--divergent```

NOTE: If you are dealing with very divergent samples/species is advised to reduce this parameter accordingly. 

###### Remove Short Sequences
This option is intended to remove extemely short sequences from alignments, in the case of UCEs, sometimes some sequences targes only the core region of the UCE and as a result they are extremely short and lack of information for phylogenetics. This options will remove those sequences shorter than a percentage of the total alignment lenght. The function works by calculating the amount of "gaps" ("N" or "-" characters) present in the sequence, for instance the default (70%) will remove any sequence composed of a 70% of ambiguous characters (gaps). thershold can be modify using the flag ```---short-cutoff```

NOTE: For any of the options above, keep in mind that FUSe will check that after removing sequences there are no columns left containing only "gaps".

###### Alignment Filtering
---

###### Length and Taxa 
This option is intended to filter out entire alignments based on its composition which might not be useful for downstream analysis, by setting a minimum alingmenth lenght (in bp) to retain with ```--min-length``` and by the minimum number of taxa present in your alignments to retain ```--min-taxa```. Default parameters are ```--min-length``` = 50bp and ```--min-taxa``` = 4 taxa.

###### Completeness
This option filter out alignments based on the total number of taxa present in your entire data set and creates "completeness matrix". This option is usually useful for exploring the effect of missing data on your phylogenetic reconstruction, but can have other utilities, default option ```--percent``` = 0.8, this mean will filter out alignments containing less than 80% of your total number of taxa.

###### Summary
FUSe.py also provide with a count of alignments per taxon `--taxa-count`. This options is a summary of your final alignments folder and can be used to diagnose any sample with very low number of loci overall, that might result in extreme long branches as a result of missing data. Another way of using this function is to check for taxa counts before getting completeness matrices, this way we can remove low input samples (high missing data) and increase the number of loci in our completeness matrix.

#### Output:
FUSe creates a series of directories to store the files processed during different parts of the workflow, this directories are stored inside a main folder named __"OUTPUT-Alignments"__ (the prefix "OUTPUT" can be changed by using a different prefix with flag ```-p```). FUSe organize the output according with the order of the workflow listed above, then the first step of the workflow will save the aligned files in __"OUTPUT-01-Mafft"__, trimmed sequences in __"OUTPUT-02-trimmed"__, filtered sequences in __"OUTPUT-05-filtered"__ and so on. If some parts of the workflow are skipped, then no folder for this step will be created.

FUSe works with FASTA formatted files troughout the workflow for convenience of coding, however you can output you final alignments in any of the formats available (fasta,nexus,phylip,clustal,emboss) using the flag ```-o```. This alignments will be saved in __"OUTPUT--Final-Alignments-{format}"__ where {format} is the file formated desired. 

NOTE: Keep in mind that files from previous steps in the workflow will be in FASTA format, this is an expected behavior of the pipeline, aligments can always be converted to specific formats using ```phyluce_align_convert_one_align_to_another``` or any other software.

## Examples
---
##### Aligning sequences without filtering:
This is the most basic function in FUSe will read a input monolithic file and will perform aligning of multiple sequences. The output will be a directory with all loci and samples included in your original monolithic file (no filtering).

Assuming you have a monolithic fasta file named `my_monolithic.fasta` with 15 loci and 20 samples and you system can run up to 12 cores:
```
FUSe.py \
-i my_monolithic.fasta 
-t 20 \
-p MyOutput \
-c 12
```
In this command we are telling FUSe.py that ourt data set is composed of 20 taxa (`-t`) the output prefix we want to use for all created directories and files is `MyOutput` and we want to run the pipeline with 12 cores (`-c`).

The output of the command above will be something like this:
```
.
├── FUSe.log
├── MyOutput-Alignments
│   └── MyOutput-01-mafft
│       ├── uce-1484.fasta
│       ├── uce-1732.fasta
│       ├── uce-1985.fasta
│       ├── uce-2120.fasta
│       ├── uce-2349.fasta
│       ├── uce-2539.fasta
│       ├── uce-3046.fasta
│       ├── uce-4179.fasta
│       ├── uce-553.fasta
│       ├── uce-6095.fasta
│       ├── uce-6187.fasta
│       ├── uce-6550.fasta
│       ├── uce-7014.fasta
│       ├── uce-7542.fasta
│       └── uce-865.fasta
├── my_monolithic.fasta
└── phyluce_align_seqcap_align.log
```
Here `FUSe.log` print all the information of the run including all the flags used and parameters, including default values
This is a good way of diagnosing any issues with the program. `phyluce_align_seqcap_align.log` is the log file for phyluce alignment script. hodling crucial information about the alignment step.
NOTE: Given that we only run alignment (the first step out of five in FUSe) we only get the `MyOutput-01-mafft` folder. 

##### Full workflow:
The opposite extreme will be to include all possible filters available in FUSe, by sample and by locus. Assuming the exact same input monolithic fasta as before us used:
```
FUSe.py \
-i my_monolithic.fasta -t 20 \
-p MyOutput -c 12 \
--gblocks \
--remove-div \
--remove-short \
--filter-alignments \
--get-completeness 
```

NOTE: The flags here works as a way of "activating" the different stages of the workflow using default parameters, if we are interested in changing such parameters, specific flags have to be provided (see example below).

The output directory structure looks like this:
```
.
├── FUSe.log
├── MyOutput-Alignments
│   ├── MyOutput-01-mafft
│   ├── MyOutput-02-trimmed
│   ├── MyOutput-03-divergent-removed
│   ├── MyOutput-04-short-removed
│   ├── MyOutput-05-filtered
│   └── MyOutput-80p
├── my_monolithic.fasta
└── phyluce_align_seqcap_align.log
```
Note that all 5 directories of the all stages of the workflow have been created, with the folder `MyOutput-80p` holding the completeness matrix with a default parameter of 80% completeness, once again FUSe.log has all parameters and information of the run.

##### Changing default parameters:
In the next example we are going to use the entire workflow but we are going to align this time using **trimAL** with the `gappyout` algorithm. 
Then we are going to relax the parameters for divergence, allowing a pairwise identity of 0.6, this mean we are keeping samples that are 60% divergent from each other, i.e. with a "higher" degree of disimilarity (for example because we are dealing with very distant related species).
Next we will allowed for sequences as short as the 50% of the total alignment, this means that some sequences will have as much as 50% of their sequence missing (coded as "N" or "-").
We will increase the minimum lenght for retaining alignments to 200bp, this mean that we will filter out any alignment shorter than 200bp. 
Finally we will request for a 70% complete matrix in `nexus` format.
```
FUSe.py \
-i my_monolithic.fasta -t 20 \
-p MyOutput -c 12 \
--trimAL -a gappyout \
--remove-div -d 0.6 \
--remove-short -d 0.5 \
--filter-alignments -m 200 \
--get-completeness -e 0.7 \
--taxa-count
```

### Bug Reports
---
contact: ``roy_monrue@hotmail.com``


### License
---
The code within this repository is available under a 3-clause BSD license. See the [LICENSE](LICENSE) file for more information.

---
---
[5]: (Calisoga Paper reference Here)

[1]: (https://mafft.cbrc.jp/alignment/software/)
[2]: (https://github.com/faircloth-lab/phyluce)
[3]: (http://trimal.cgenomics.org/use_of_the_command_line_trimal_v1.2)
[4]: (https://home.cc.umanitoba.ca/~psgendb/doc/Castresana/Gblocks_documentation.html)
[6]: (https://phyluce.readthedocs.io/en/latest/installation.html)
[7]: (https://www.ncbi.nlm.nih.gov/genbank/fastaformat/)
[8]: (https://cialign.readthedocs.io/en/latest/pages/introduction.html)

