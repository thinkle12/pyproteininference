# Requirements

1. __Python 3.7__ or greater. This package was created using __Python 3.7__
2. __Python Packages__:
   __numpy__, __pyteomics__, __pulp__, __PyYAML__, __matplotlib__. These should be installed automatically during installation.
		
# Quick Start Guide
1. Install the package using pip:
   
		pip install pyproteininference
   
2. Run the Heuristic method with tab delimited results directly from percolator to generate results for peptide centric, parsimony, inclusion, and exclusion:

		protein_inference_heuristic_cli.py \
			-t /path/to/target/file1.txt \
			-d /path/to/decoy/file1.txt \
			-db /path/to/database/file.fasta 

3. Run the standard commandline tool with tab delimited results directly from percolator to run a particular inference method. By default, peptide centric inference is selected if a parameter file is not specified:

		protein_inference_cli.py \
			-t /path/to/target/file.txt \
			-d /path/to/decoy/file.txt \
			-db /path/to/database/file.fasta 

4. Specifying Parameters. 
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

5. Run the standard commandline tool again, this time specifying the parameters as above:
		
		protein_inference_cli.py \
			-t /path/to/target/file.txt \
			-d /path/to/decoy/file.txt \
			-db /path/to/database/file.fasta \
			-y /path/to/params.yaml

6. Running with docker
	
	- Either Pull the image from docker hub:
		- `docker pull hinklet/pyproteininference:latest`
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
	