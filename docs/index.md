## Introduction


For a quick start guide please click [here](quickstart.md).

Py Protein Inference is a Python package that has the ability to run various protein inference algorithms on tandem mass spectrometry search results. 
In addition to performing protein inference, which maps peptides to proteins, this algorithm creates protein scores based on the supplied psms and is able to calculate set based protein level false discovery rates for MS data filtering purposes.
Py Protein Inference typically takes as input the output PSM files from the [Percolator algorithm](https://github.com/percolator/percolator). 
However, Py Protein Inference can also take custom tab delimited files as input. Please see the [input formats](input_format.md) section for more information.
As for output, Py Protein Inference generates a user-friendly csv format file that includes the Proteins, Peptides, q-values, and Protein Scores.  Please see the [supplementary](supplementary.md#export-types) section for more information on output formats. 
Py Protein Infernece also has a heuristic method to help the user select a recommended inference method for a given dataset. Please see the [supplementary](supplementary.md#heuristic-algorithm) section for more information on the heuristic method.

Py Protein Inference has the ability to run any of the following inference procedures from literature:

1. [Parsimony](supplementary.md#parsimony)
2. [Exclusion](supplementary.md#exclusion)
3. [Inclusion](supplementary.md#inclusion)
4. [Peptide Centric](supplementary.md#peptide-centric) (Protein Group Level)

Please see the [__Inference Types__](supplementary.md#inference-types) section for more information on Inference Types.

In Addition to these inference types Py Protein Inference can also score proteins with a variety of methods:

1. Best Peptide Per Protein 
    - Takes the best scoring PSM as the protein score.
2. Multiplicative Log 
    - Multiplies each PSM score per protein and takes -log10 of the combined score (smaller psm scores must be better i.e. Pep or Q values).
3. Top Two Combined 
    - Takes the top two best scoring PSMs and multiplies or adds them together based on if the selected psm score is Pep/Q value style or Mascot Ion Score style.
4. Additive 
    - Adds each PSM score per protein (larger psm scores must be better i.e. Xcorr, Mascot Ion Score, Percolator Score).
5. Downweighted Multiplicative Log 
    - Multiplicative Log but with a normalization against the number of PSMs.
6. Geometric Mean 
    - Geometric Mean algorithm scoring.

Please see the [__Protein Score Types__](supplementary.md#protein-score-types) section for more information on scoring algorithms.

## Using Py Protein Inference
 1. [Yaml Parameter File](parameters.md#yaml-parameter-file-outline)
 2. [Input PSM files](input_format.md#input-file-examples) (Tab Delimited)
 3. [Fasta Database](input_format.md#fasta-file-example)
 4. [Running Py Protein Inference](advanced.md#running-py-protein-inference)