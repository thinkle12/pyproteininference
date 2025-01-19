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
### Install the package using pip and running the CLI/GUI:
```shell
pip install pyproteininference
```

### Running the commandline tool

To run the CLI tool either call `protein_inference_cli.py` like so:
```shell
protein_inference_cli.py --help
```

Or call the script while also calling your python interpreter

First, locate the script that gets installed on installation:
```shell
which protein_inference_cli.py
/path/to/venv/bin/protein_inference_cli.py
```

Then, call the script while also calling your python interpreter
```shell
python /path/to/venv/bin/protein_inference_cli.py --help
```

Optionally, download the protein_inference_cli.py file from the github repo here:
https://github.com/thinkle12/pyproteininference/blob/master/scripts/protein_inference_cli.py

And then call the script while also calling the pyton interpreter as shown above

### Running the graphical user interface

To run the GUI tool either call `protein_inference_gui.py` like so:
```shell
protein_inference_gui.py
```

Or again, call the script while also calling your python interpreter

First, locate the script that gets installed on installation:
```shell
which protein_inference_gui.py
/path/to/venv/bin/protein_inference_gui.py
```

Then, call the script while also calling your python interpreter
```shell
python /path/to/venv/bin/protein_inference_gui.py
```

Again, you can optionally download the protein_inference_gui.py file from the github repo here:
https://github.com/thinkle12/pyproteininference/blob/master/scripts/protein_inference_gui.py

And then call the script while also calling the pyton interpreter as shown above

### Executables

You can also install executbles for the GUI for both Windows and MacOS on the releases page on github:
https://github.com/thinkle12/pyproteininference/releases

When launching the GUI's from the executables please wait until for the user interface to pop up. It usually takes a minute or so.

## Full Options for calling the CLI

1. Run the standard commandline from an idXML file 
```shell
protein_inference_cli.py \
-f /path/to/target/file.idXML \
-db /path/to/database/file.fasta \
-y /path/to/params.yaml
```
   
2. Run the standard commandline from an mzIdentML file 
```shell
protein_inference_cli.py \
-f /path/to/target/file.mzid \
-db /path/to/database/file.fasta \
-y /path/to/params.yaml
```
   
3. Run the standard commandline from a pepXML file 
```shell
protein_inference_cli.py \
-f /path/to/target/file.pep.xml \
-db /path/to/database/file.fasta \
-y /path/to/params.yaml
```

4. Run the standard commandline tool with tab delimited results directly from percolator to run a particular inference method. By default, peptide centric inference is selected if a parameter file is not specified:
```shell
protein_inference_cli.py \
-t /path/to/target/file.txt \
-d /path/to/decoy/file.txt \
-db /path/to/database/file.fasta 
```

5. Specifying Parameters. 
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

6. Run the standard commandline tool again, this time specifying the parameters as above:
```shell
protein_inference_cli.py \
-t /path/to/target/file.txt \
-d /path/to/decoy/file.txt \
-db /path/to/database/file.fasta \
-y /path/to/params.yaml
```

7. Running with docker
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

## Building the Bundled Application Package using PyInstaller
_Note: This is only necessary if you want to build the application package yourself. The package is already available on
PyPi and can be installed using pip, or bundled executables can be downloaded from the releases page on 
GitHub (https://thinkle12.github.io/pyproteininference/)._

1. After cloning the source code repository, create a new Python virtual environment under the project directory:
```shell
python -m venv venv
```
2. Activate the virtual environment:
```shell
source venv/bin/activate
```
3. Install the required packages:
```shell
pip install -r requirements.txt pyinstaller==6.11.1
```
4. Run the PyInstaller command to build the executable:
```shell
pyinstaller pyProteinInference.spec
```
5. The executable will be located in the `dist` directory.



## Documentation
For more information please see the full package documentation (https://thinkle12.github.io/pyproteininference/).
