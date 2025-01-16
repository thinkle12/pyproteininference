# Requirements

1. __Python 3.9__ or greater.
2. __Python Packages__:
   __numpy__, __pyteomics__, __pyopenms__, __pulp__, __PyYAML__, __matplotlib__, __lxml__, __nicegui__, __pywebview__, __tqdm__. These should be installed automatically during installation.
		
# Quick Start Guide
1. Install the package using pip:
   
		pip install pyproteininference
   
2. Run the standard commandline from an idXML file 

		protein_inference_cli.py \
			-f /path/to/target/file1.idXML \
			-db /path/to/database/file.fasta \
			-y /path/to/params.yaml

3. Run the standard commandline from a mzIdentML file 

		protein_inference_cli.py \
			-f /path/to/target/file1.mzid \
			-db /path/to/database/file.fasta \
			-y /path/to/params.yaml

4. Run the standard commandline from a pepXML file 

		protein_inference_cli.py \
			-f /path/to/target/file1.pepXML \
			-db /path/to/database/file.fasta \
			-y /path/to/params.yaml

5. Run the standard commandline tool with tab delimited results directly from percolator

		protein_inference_cli.py \
			-t /path/to/target/file.txt \
			-d /path/to/decoy/file.txt \
			-db /path/to/database/file.fasta \ 
     		-y /path/to/params.yaml


6. Specifying Parameters. 
The two most common parameters to change are the inference type, and the decoy symbol (for identifying decoy proteins vs target proteins).
The parameters can be quickly altered by creating a file called params.yaml as follows:

		parameters:
		  inference:
			inference_type: parsimony
		  identifiers:
			decoy_symbol: "decoy_"

	The inference type can be one of: `parsimony`, `peptide_centric`, `inclusion`, `exclusion`, or `first_protein`.
	All parameters are optional, so you only need to define the ones you want to alter. Parameters that are not defined are set to default values.
	See [here](parameters.md#default-parameters) for the default parameters.

7. Full Parameter Specifications
See below for a full standard parameter file:

## Default Parameters
```yaml
parameters:
  general:
    export: peptides
    fdr: 0.01
    picker: True
    tag: example_tag
    xml_parser: openms
  data_restriction:
    pep_restriction: 0.9
    peptide_length_restriction: 7
    q_value_restriction: .9
    custom_restriction: None
    max_allowed_alternative_proteins: 50
  score:
    protein_score: best_peptide_per_protein
    psm_score: posterior_error_prob
    psm_score_type: multiplicative
  identifiers:
    decoy_symbol: "##"
    isoform_symbol: "-"
    reviewed_identifier_symbol: "sp|"
  inference:
    inference_type: parsimony
    grouping_type: parsimonious_grouping
  digest:
    digest_type: trypsin
    missed_cleavages: 3
  parsimony:
    lp_solver: pulp
    shared_peptides: all
  peptide_centric:
    max_identifiers: 5
```

These parameter options are just a suggestion. Please alter these for your specifications. 
For full description of each parameter and all options see the in depth [parameter file description](parameters.md#yaml-parameter-file-outline)

8. Run the standard commandline tool again, this time specifying the parameters as above:
		
		protein_inference_cli.py \
			-t /path/to/target/file.txt \
			-d /path/to/decoy/file.txt \
			-db /path/to/database/file.fasta \
			-y /path/to/params.yaml

9. Running with docker
	
	- Either Pull the image from docker hub:
		- `docker pull thinkle12/pyproteininference:latest`
	- Or Build the image with the following command (After having cloned the repository):
	  	- `git clone REPOSITORY_URL`
	  	- `cd pyproteininference`
		- `docker build -t pyproteininference:latest .`
	- Run the tool, making sure to volume mount in the directory with your input data and parameters. In the case below, that local directory would be `/path/to/local/directory` and the path in the container is `/data`

			docker run -v /path/to/local/directory/:/data \
				-it hinklet/pyproteininference:latest \
				python /usr/local/bin/protein_inference_cli.py \
				-f /data/input_file.txt \
				-db /data/database.fasta \
				-y /data/parameters.yaml \
				-o /data/
	
	- Get the commandline help via docker

			docker run thinkle12/pyproteininference:latest \
            python /usr/local/bin/protein_inference_cli.py --help