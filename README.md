# Py Protein Inference

**PyProteinInference** is a Python package for running various protein inference algorithms on tandem mass spectrometry search results and generating protein to peptide mappings with protein level false discovery rates..  

## Key Features

* **Protein Inference and Scoring**:
    * Maps peptides to proteins.  
    * Generates protein scores from provided PSMs.  
    * Calculates set-based protein-level false discovery rates for MS data filtering.  
* **Supported Input Formats**:
    * Search Result File Types: __idXML__, __mzIdentML__, or __pepXML__.  
    * PSM files from [Percolator](https://github.com/percolator/percolator).
    * Custom tab-delimited files.  
* **Output**:
    * User-friendly CSV file containing Proteins, Peptides, q-values, and Protein Scores.  

* **Supported Inference Procedures**:
    * Parsimony - Returns the Minimal set of proteins based on the input peptides.
    * Exclusion - Removes all non-distinguishing peptides on the protein level.
    * Inclusion - Returns all possible proteins.
    * Peptide Centric - Returns protein groups based on peptide assignments.

## Requirements

 1. __Python 3.9__ or greater. 
 2. __Python Packages__:
	__numpy__, __pyteomics__, __pulp__, __PyYAML__, __matplotlib__, __pyopenms__, __lxml__, __tqdm__, __pywebview__, __nicegui__. These should be installed automatically during installation.
		
## Quick Start Guide
1. Install the package using pip:
```shell
pip install pyproteininference
```
   
2. Run the standard commandline from an idXML file 
```shell
protein_inference_cli.py \
-f /path/to/target/file.idXML \
-db /path/to/database/file.fasta \
-y /path/to/params.yaml
```
   
3. Run the standard commandline from an mzIdentML file 
```shell
protein_inference_cli.py \
-f /path/to/target/file.mzid \
-db /path/to/database/file.fasta \
-y /path/to/params.yaml
```
   
4. Run the standard commandline from a pepXML file 
```shell
protein_inference_cli.py \
-f /path/to/target/file.pep.xml \
-db /path/to/database/file.fasta \
-y /path/to/params.yaml
```

5. Run the standard commandline tool with tab delimited results directly from percolator to run a particular inference method. By default, peptide centric inference is selected if a parameter file is not specified:
```shell
protein_inference_cli.py \
-t /path/to/target/file.txt \
-d /path/to/decoy/file.txt \
-db /path/to/database/file.fasta 
```

6. Specifying Parameters. 
The two most common parameters to change are the inference type, and the decoy symbol (for identifying decoy proteins vs target proteins).
The parameters can be quickly altered by creating a file called params.yaml as follows:
```yaml
parameters:
  inference:
    inference_type: parsimony
  identifiers:
    decoy_symbol: "decoy_"
```
The inference type can be one of: `parsimony`, `peptide_centric`, `inclusion`, `exclusion`, or `first_protein`.
All parameters are optional, so you only need to define the ones you want to alter. Parameters that are not defined are set to default values.
See the package documentation for the default parameters.

7. Run the standard commandline tool again, this time specifying the parameters as above:
```shell
protein_inference_cli.py \
-t /path/to/target/file.txt \
-d /path/to/decoy/file.txt \
-db /path/to/database/file.fasta \
-y /path/to/params.yaml
```

8. Running with docker
	- Either Pull the image from docker hub:
		- `docker pull hinklet/pyproteininference:latest`
	- Or Build the image with the following command (After having cloned the repository):
	  	- `git clone REPOSITORY_URL`
	  	- `cd pyproteininference`
		- `docker build -t pyproteininference:latest .`
	- Run the tool, making sure to volume mount in the directory with your input data and parameters. In the case below, that local directory would be `/path/to/local/directory` and the path in the container is `/data`
	  ```shell
			docker run -v /path/to/local/directory/:/data \
			-it hinklet/pyproteininference:latest \
			python /usr/local/bin/protein_inference_cli.py \
			-f /data/input_file.txt \
			-db /data/database.fasta \
			-y /data/parameters.yaml \
			-o /data/
	  ```

## Documentation
For more information please see the full package documentation (https://thinkle12.github.io/pyproteininference/).
