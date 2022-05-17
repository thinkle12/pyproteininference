## Introduction


For a quick start guide please click [here](quickstart.md).

Py Protein Inference is a Python package that has the ability to run various protein inference algorithms on tandem mass spectrometry search results. 
In addition to performing protein inference, which maps peptides to proteins, this algorithm creates protein scores based on the supplied psms and is able to calculate set based protein level false discovery rates for MS data filtering purposes.
Py Protein Inference typically takes as input the output PSM files from the [Percolator algorithm](https://github.com/percolator/percolator). 
However, Py Protein Inference can also take custom tab delimited files as input. Please see the [input formats](input_format.md) section for more information.
As for output, Py Protein Inference generates a user-friendly csv format file that includes the Proteins, Peptides, q-values, and Protein Scores.  Please see the [supplementary](supplementary.md#export-types) section for more information on output formats.

Py Protein Inference has the ability to run any of the following inference procedures from literature:

1. [Parsimony](supplementary.md#parsimony)
2. [Exclusion](supplementary.md#exclusion)
3. [Inclusion](supplementary.md#inclusion)
4. [Peptide Centric](supplementary.md#peptide-centric) (Protein Group Level)
5. [First Protein](supplementary.md#first-protein) (Selects first protein per peptide)

Please see the [__Inference Types__](supplementary.md#inference-types) section for more information on Inference Types.

In Addition to these inference types Py Protein Inference can also score proteins with a variety of methods:

1. Best Peptide Per Protein
2. Multiplicative Log
3. Top Two Combined
4. Additive
5. Downweighted Multiplicative Log
6. Geometric Mean

Please see the [__Protein Score Types__](supplementary.md#protein-score-types) section for more information on scoring algorithms.

## Using Py Protein Inference
 1. [Yaml Parameter File](parameters.md#yaml-parameter-file-outline)
 2. [Input PSM files](input_format.md#input-file-examples) (Tab Delimited)
 3. [Fasta Database](input_format.md#fasta-file-example)
 4. [Running Py Protein Inference](advanced.md#running-py-protein-inference)