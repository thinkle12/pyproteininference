#!/usr/bin/python

import argparse

import pyproteininference


def main():
    """
    Script function for running the execute method of the ProteinInference HeuristicPipeline class.

    """
    parser = argparse.ArgumentParser(description="Protein Inference Heuristic")
    parser.add_argument(
        "-t",
        "--target",
        dest="target",
        required=False,
        help="Input target psm output from percolator. Can either input one file or a list of files.",
        metavar="FILE",
        nargs="+",
        type=str,
    )
    parser.add_argument(
        "-d",
        "--decoy",
        dest="decoy",
        required=False,
        help="Input decoy psm output from percolator. Can either input one file or a list of files.",
        metavar="FILE",
        nargs="+",
        type=str,
    )
    parser.add_argument(
        "-f",
        "--combined_files",
        dest="combined_files",
        required=False,
        help="Input combined psm output from percolator. This should contain Target and Decoy PSMS. "
        "Can either input one file or a list of files.",
        metavar="FILE",
        nargs="+",
        type=str,
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="dir_name",
        required=False,
        help="Result Directory to write to - the name of file will be determined by parameters selected "
        "and tag parameter. If this option is not set, will write results to current working directory.",
        metavar="DIR",
    )
    parser.add_argument(
        "-l",
        "--output_filename",
        dest="output_filename",
        required=False,
        help="Filename to write results to. Can be left blank. If this flag is left blank the filename will be"
        " automatically generated. If set this flag will override -o.",
        metavar="FILE",
        default=None,
    )
    parser.add_argument(
        "-a",
        "--target_directory",
        dest="target_directory",
        required=False,
        help="Directory that contains either .txt or .tsv input target psm data. Make sure the directory ONLY contains"
        " result files.",
        metavar="DIR",
    )
    parser.add_argument(
        "-b",
        "--decoy_directory",
        dest="decoy_directory",
        required=False,
        help="Directory that contains either .txt or .tsv input decoy psm data. Make sure the directory ONLY contains."
        " result files.",
        metavar="DIR",
    )
    parser.add_argument(
        "-c",
        "--combined_directory",
        dest="combined_directory",
        required=False,
        help="Directory that contains either .txt or .tsv input data with targets/decoys combined. "
        "Make sure the directory ONLY contains result files.",
        metavar="DIR",
    )
    parser.add_argument(
        "-db",
        "--database",
        dest="database",
        required=False,
        help="Path to the fasta formatted database used in the MS search. This is optional. If not set, "
        "will use the proteins only in the input files.",
        metavar="FILE",
    )
    parser.add_argument(
        "-y",
        "--yaml_params",
        dest="yaml_params",
        required=False,
        help="Path to a Protein Inference Yaml Parameter File. If this is not set, default parameters will be used.",
        metavar="FILE",
    )
    parser.add_argument(
        "-p",
        '--skip_append_alt',
        action="store_true",
        default=False,
        help="Advanced usage only. If this flag is set, will skip adding alternative proteins to each PSM "
        "from the database digest. If this flag is not set, the peptide/protein mapping will be taken "
        "from database digest and appended to the mapping present in the input files.",
    )
    parser.add_argument(
        "-i",
        '--id_splitting',
        action="store_true",
        default=False,
        help="Advanced usage only. If set this flag will split protein identifiers."
        "If not set, this flag will not split protein identifiers."
        "Sometimes the fasta database protein IDs are formatted as: 'sp|ARAF_HUMAN|P10398'. "
        "While protein IDs in the input files are formatted as 'ARAF_HUMAN|P10398'. "
        "Setting This flag will split off the front 'sp|' or 'tr|' from the database protein identifiers. ",
    )
    parser.add_argument(
        "-r",
        "--pdf_filename",
        dest="pdf_filename",
        required=False,
        help="PDF Filepath to write the Heuristic plot to after Heuristic Scoring. "
        "If not set, writes the file with filename heuristic_plot.pdf to directory set in -o. "
        "If -o is not set, will write the file to current working directory.",
        metavar="FILE",
    )
    parser.add_argument(
        "-m",
        "--fdr_threshold",
        dest="fdr_threshold",
        required=False,
        help="The FDR threshold to use in the Heuristic Method. Defaults to 0.05 if not set.",
        metavar="FLOAT",
        type=float,
        default=0.05,
    )
    parser.add_argument(
        "-u",
        "--output_type",
        dest="output_type",
        required=False,
        help="The type of output to be written. Can either be 'all' or 'optimal'. "
        "If set to 'all' will output all inference results. If set to 'optimal' "
        "will output only the result selected by the heuristic method. If left blank "
        "this will default to 'all'.",
        metavar="STR",
        type=str,
        default="all",
    )

    args = parser.parse_args()

    pipeline = pyproteininference.heuristic.HeuristicPipeline(
        parameter_file=args.yaml_params,
        database_file=args.database,
        target_files=args.target,
        decoy_files=args.decoy,
        combined_files=args.combined_files,
        target_directory=args.target_directory,
        decoy_directory=args.decoy_directory,
        combined_directory=args.combined_directory,
        output_directory=args.dir_name,
        output_filename=args.output_filename,
        append_alt_from_db=not args.skip_append_alt,  # Need to reverse the Boolean here
        id_splitting=args.id_splitting,
        pdf_filename=args.pdf_filename,
        output_type=args.output_type,
    )
    pipeline.execute(fdr_threshold=args.fdr_threshold)


if __name__ == "__main__":
    main()
