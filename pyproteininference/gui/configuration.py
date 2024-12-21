from dataclasses import dataclass, field


@dataclass
class Configuration:
    input_files: list = field(default_factory=lambda: [])
    fasta_file: str = ''
    output_file: str = ''
    false_discovery_rate: float = 0.01
    picker: bool = True
    use_alt_proteins: bool = True
    pep_restriction: float = 0.9
    q_value_restriction: float = 0.9
    peptide_length_restriction: int = 0
    protein_score: str = 'multiplicative_log'
    psm_score: str = 'posterior_error_prob'
    psm_score_custom: str = ''
    psm_score_type: str = 'multiplicative'
    decoy_symbol: str = 'DECOY_'
    isoform_symbol: str = '-'
    reviewed_identifier_symbol: str = 'sp|'
    identifier_splitting: bool = False
    inference_type: str = 'parsimony'
    grouping_type: str = 'parsimonious_grouping'
    digest_type: str = 'trypsin'
    missed_cleavages: int = 2
    shared_peptides: str = 'all'
    max_identifiers: int = 5
    xml_input_parser: str = 'openms'
    max_allowed_alternative_proteins: int = 50
