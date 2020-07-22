import protein_inference
import argparse


parser = argparse.ArgumentParser(description="Protein Inference")
parser.add_argument(
    "-t",
    "--target",
    dest="target",
    required=False,
    help="Input target psm output from percolator. Can either input one file or a list of files",
    metavar="FILE",
    nargs="+",
    type=str,
)
parser.add_argument(
    "-d",
    "--decoy",
    dest="decoy",
    required=False,
    help="Input decoy psm output from percolator. Can either input one file or a list of files",
    metavar="FILE",
    nargs="+",
    type=str,
)
parser.add_argument(
    "-f",
    "--combined_files",
    dest="combined_files",
    required=False,
    help="Input combined psm output from percolator. This should contain Target and Decoy PSMS. Can either input one file or a list of files",
    metavar="FILE",
    nargs="+",
    type=str,
)
parser.add_argument(
    "-o",
    "--output",
    dest="dir_name",
    required=False,
    help="protein_inference Result Directory to write to - Name of file will be determined by parameters selected and parameter tag",
    metavar="DIR",
)
parser.add_argument(
    "-a",
    "--target_directory",
    dest="target_directory",
    required=False,
    help="Directory that contains either .txt or .tsv input target psm data. Make sure the directory ONLY contains result files",
    metavar="DIR",
)
parser.add_argument(
    "-b",
    "--decoy_directory",
    dest="decoy_directory",
    required=False,
    help="Directory that contains either .txt or .tsv input decoy psm data. Make sure the directory ONLY contains result files",
    metavar="DIR",
)
parser.add_argument(
    "-c",
    "--combined_directory",
    dest="combined_directory",
    required=False,
    help="Directory that contains either .txt or .tsv input data with targets/decoys combined. Make sure the directory ONLY contains result files",
    metavar="DIR",
)
parser.add_argument(
    "-db",
    "--database",
    dest="database",
    required=False,
    help="Provide the fasta formatted database used in the MS search",
    metavar="FILE",
)
parser.add_argument(
    "-y",
    "--yaml_params",
    dest="yaml_params",
    required=True,
    help="Provide a Protein Inference Yaml Parameter File",
    metavar="FILE",
)
parser.add_argument(
    "-p",
    "--append_alt",
    dest="append_alt",
    required=False,
    type=protein_inference.pipeline.ProteinInferencePipeline.str2bool,
    nargs='?',
    const=True,
    default=True,
    metavar="BOOL",
    help="Set '--append_alt True' or '--append_alt False'. "
    "Whether or not to add alternative proteins to each PSM from the database digest. "
    "If False the peptide/protein mapping will be taken from the input files only. "
    "If this is left blank it will default to True.",
)
parser.add_argument(
    "-i",
    "--id_splitting",
    dest="id_splitting",
    required=False,
    type=protein_inference.pipeline.ProteinInferencePipeline.str2bool,
    nargs='?',
    const=True,
    default=False,
    metavar="BOOL",
    help="Set '--id_splitting True' or '--id_splitting False'. "
    "This flag is by default False. So if it is left blank then the variable is set to False. "
    "This flag indicates whether or not to split the identifiers that are present in the fasta database. "
    "Only set this flag to True if you know what you are doing. "
    "Sometimes the fasta database protein IDs will be like: 'sp|ARAF_HUMAN|P10398'. "
    "While protein IDs in the input files will be 'ARAF_HUMAN|P10398'. "
    "Setting This flag to True will split off the front 'sp|' from the database protein identifiers. "
    "This is typically not necessary. So leave this blank unless you know what you are doing.",
)

def main():
    """
    Script function for running the execute method of the ProteinInferencePipeline class

    """
    args = parser.parse_args()

    pipeline = protein_inference.pipeline.ProteinInferencePipeline(
        parameter_file=args.yaml_params,
        database_file=args.database,
        target_files=args.target,
        decoy_files=args.decoy,
        combined_files=args.combined_files,
        target_directory=args.target_directory,
        decoy_directory=args.decoy_directory,
        combined_directory=args.combined_directory,
        output_directory=args.dir_name,
        append_alt_from_db=args.append_alt,
        id_splitting=args.id_splitting,
    )
    pipeline.execute()

if __name__ == "__main__":
    main()
