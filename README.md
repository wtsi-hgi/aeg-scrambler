# AEG SCRAMBLEr

## Introduction

#### Step-by-step

1. The user may change the config file to change settings within the program as desired.
2. The rank command allows the user to find various metrics of the genes and order them by a weighted score.
3. The tune command allows the user to give a desired distribution of interest for the program to try to reproduce.
4. The design command allows the user, once they have chosen their weights, to find suggested insertions around genes.

## How to install

Run the following command to install aeg-scrambler:

~~~
command --to run
~~~

## Configuration

All settings for the tool are provided by the user as a .json file. This is
split into the following sections:

#### Paths

The paths which the tool will look for certain files are provided:

~~~
    "results_directory" : "/my/path/here/results/",
    "gene_report_directory" : "/my/path/here/reports/",
    "gene_annotation_path" : "/my/path/here/gene_annotations.gtf",
    "regulatory_elements_path" : "/my/path/here/epigenetic_flags.bed",
    "ccle_expression_path" : "/my/path/here/CCLE_expressions_data.csv",
    "experimental_expression_path" : "/my/path/here/experimental_expressions.tsv",
    "hic_path" : "/my/path/here/hic_data.file",
    "reference_genome_path" : "/my/path/here/reference.fa",
~~~

All paths must be provided as strings, directories must be followed with a
forward slash. Absolute paths are preferred.

#### Biological specifications

These settings are used to specify biological details of the data:

~~~
    "cell_line_of_interest" : "CELL-LINE",
    "chromosomes_of_interest" : ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"],
    "flags_of_interest" : ["FLAG1"],
    "quiescent_epigenetic_flags_of_interest" : ["FLAG2"],
~~~

#### Search window

Allows the user to define the search window that searches around each gene:

~~~
    "search_type" : "whole_gene",
    "upstream_search" : 500000,
    "downstream_search" : 500000,
~~~

"search_type" can be either "whole_gene" (upstream and downstream are defined
from the ends of the gene) or "start_site" (upstream and downstream are
defined from the start site of the gene).
"upstream search" and "downstream search" must be integers, and define the
bases included upstrem and downstream within the search window.

#### Hard filters

These filters are used to cut out undesired genes from any processing:

~~~
    "std_max" : false,
    "std_min" : 1,
    "anomalous_expression_max" : false,
    "anomalous_expression_min" : false,
    "enhancer_count_max" : false,
    "enhancer_count_min" : false,
    "enhancer_proportion_max" : false,
    "enhancer_proportion_min" : false,
    "cell_line_expression_max" : false,
    "cell_line_expression_min" : false,
    "gene_size_max" : false,
    "gene_size_min" : false,
    "symmetry_max" : false,
    "symmetry_min" : false,
~~~

Filters may either be false (no filter) or a number (integer / float
as necessary)

#### Weights

Weights are set to decide how genes are prioritised by the "rank" function:

~~~
    "std_weight" : 0,
    "anomalous_expression_weight" : 0.5,
    "enhancer_count_weight" : 1,
    "enhancer_proportion_weight" : 1.5,
    "cell_line_expression_weight" : 2,
    "gene_size_weight" : 2.5,
    "symmetry_weight" : 3,
~~~

These may be integers or floats, positive, negative or zero.

#### Enhancer kernel parameters

The following parameters define the kernel used to convole the signal of
enhancers:

~~~
    "enhancer_kernel_shape" : "guassian",
    "enhancer_kernel_size_type" : "relative",
    "absolute_enhancer_kernel_size" : 1000,
    "absolute_enhancer_kernel_sigma" : 3,
    "relative_enhancer_kernel_size" : 0.1,
    "relative_enhancer_kernel_sigma" : 0.005,
~~~

"enhancer_kernel_shape" may be "flat" or "guassian".
"enhancer_kernel_size_type" may be "absolute" or "relative", and defines
whether the absolute or relative settings are used; an absolute kernel is the
same for each gene and the size is the number of bases wide it is, a relative
kernel has its size and sigma relative to the gene which is being convolved.

#### Quiescent kernel parameters

The following parameters define the kernel used to convole the signal of
quiescent regions:

~~~
    "quiescent_kernel_shape" : "guassian",
    "quiescent_kernel_size_type" : "relative",
    "absolute_quiescent_kernel_size" : 1000,
    "absolute_quiescent_kernel_sigma" : 3,
    "relative_quiescent_kernel_size" : 0.1,
    "relative_quiescent_kernel_sigma" : 0.0015,
~~~

"quiescent_kernel_shape" may be "flat" or "guassian".
"quiescent_kernel_size_type" may be "absolute" or "relative", and defines
whether the absolute or relative settings are used; an absolute kernel is the
same for each gene and the size is the number of bases wide it is, a relative
kernel has its size and sigma relative to the gene which is being convolved.

#### Interferring genes

These settings define how genes that are nearby other genes affect each gene's
search window, if a gene is defined as interferring, it will cut off a search
window prematurely.

~~~
    "specific_expression_threshold" : 0.01,
    "interferring_gene_overlaps" : false,
~~~

The 'specific_expression_threshold' is the expression of a gene within the cell 
line of interest needed for a gene to count as interferring.
'interferring_gene_overlaps' is whether to consider genes which overlap the gene
of interest as potentially being interferring, functionality may not work as expected.

#### Convolution settings

These are used to control how the convolutions are generated and call pateaus.

~~~
    "convolution_limit" : 2,
    "plateau_threshold" : 0.1,
~~~

'convolution_limit' is the number of genes which will have their
convolutions calculated; it is a slightly slow process so recommened
to reduce the number for faster performance. Genes are considered in
their ranked order, so genes with higher interest scores will have
their convolution calculated first. 'plateau_threshold' is the
horizontal line the convolved signal must cross to generate a plateau.
The signal is in arbitrary units but a value close to 0 (but not 0)
is recommended.

#### Prime editing settings

These settings decide what is inserted into each plateau and how many
insertions should be considered.

~~~
    "inserted_sequence" : "ATAACTTCGTATAATGTACATTATACGAAGTTAT",
    "partial_insertions_per_region" : 20,
~~~

'inserted_sequence' is the sequence being inserted, it will be read
by PRIDICT, so their conventions must be used. 'partial_insertions_per_region'
is the number of prefix/suffix combinations which will be found in each plateau
before they are fed into PRIDICT.

## Commands

#### Rank

The rank command can be used to produce a list of genes selected in order of
weighted features, a tsv file will be created at the desired location with the
full list of genes and the top 50 will be printed to the terminal.

Example Rank command:

~~~
aeg-scrambler rank --config /my/path/to/config.json
~~~

#### Tune

The tune command 

#### Design

Example Design command:

~~~
aeg-scrambler design --config /my/path/to/config.json
~~~
