## Introduction

**PyProteinInference** is a Python package for running various protein inference algorithms on tandem mass spectrometry search results and generating protein to peptide mappings with protein level false discovery rates.  

For a quick start guide please click [here](quickstart.md).

**Key Features** <br>

* **Protein Inference and Scoring**:
    * Maps peptides to proteins.  
    * Generates protein scores from provided PSMs.  
    * Calculates set-based protein-level false discovery rates for MS data filtering.  
* **Supported Input Formats**:
    * Search Result File Types: [idXML](input_format.md#idxml), [mzIdentML](input_format.md#mzidentml), or [pepXML](input_format.md#pepxml).  
    * PSM files from [Percolator](https://github.com/percolator/percolator).
    * Custom tab-delimited [files](input_format.md#custom-input).  
* **Output**:
    * User-friendly CSV file containing Proteins, Peptides, q-values, and Protein Scores.  
    * Details on output formats: [supplementary](supplementary.md#export-types).  

* **Supported Inference Procedures**:
    * [Parsimony](supplementary.md#parsimony)
    * [Exclusion](supplementary.md#exclusion)
    * [Inclusion](supplementary.md#inclusion)
    * [Peptide Centric](supplementary.md#peptide-centric) (Protein Group Level)

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
5. Down Weighted Multiplicative Log 
    - Multiplicative Log but with a normalization against the number of PSMs.
6. Iterative Down Weighted Log 
    - Multiplicative Log but with an iterative normalization against the number of PSMs.
7. Geometric Mean 
    - Geometric Mean algorithm scoring.

Please see the [__Protein Score Types__](supplementary.md#protein-score-types) section for more information on scoring algorithms.

## Using Py Protein Inference
 1. [Yaml Parameter File](parameters.md#yaml-parameter-file-outline)
 2. [Input File Examples](input_format.md#input-file-examples) (idXML, mzIdentML, pepXML, Tab Delimited)
 3. [Fasta Database](input_format.md#fasta-file)
 4. [Running Py Protein Inference](advanced.md#running-py-protein-inference)