# Yaml Parameter File Outline
The Yaml Parameter File is the central location for all configurations for a given Protein Inference run and are summarized below:
Note: These parameters are all optional. Please see the section [Default Parameters](#default-parameters) for more information on defaults.
## General
| Parameter | Description |Type|
|---|---|---|
| export | Export Type can be one of: __peptides__, __psms__, __psm_ids__, __long__, __q_value__, __q_value_all__, __q_value_comma_sep__, __leads__, __all__, __comma_sep__. Suggested types are __peptides__, __psms__, and __psm_ids__ as these produce square output. If there are multiple proteins per group the three mentioned types will report the leads only. Other types report on the peptide level with slightly different formats and whether or not to include leads only or all proteins. See [here](supplementary.md#export-types) for an in-depth explanation of Export Types. | String |
| fdr | False Discovery Rate to be marked as significant. Ex. __0.01__ for 1% FDR. | Numeric |
| picker | __True__/__False__ on whether to run the Protein Picker algorithm. For more info click [here](supplementary.md#protein-picker). | Bool |
| tag | A String tag that will be written into the result files. Ex. __example_tag__. | String |

## Data Restriction
| Parameter | Description |Type|
|---|---|---|
| pep_restriction | Posterior Error Probability values to filter (i.e. __0.9__). In this case PSMs with PEP values greater than __0.9__ would be removed from the input. If PEP values are not in input please use __None__.  | Numeric |
| peptide_length_restriction | Peptide Length to filter on. (i.e. __7__). If no filter is wanted please use __None__. | Int |
| q_value_restriction | Q Values to filter. (i.e. __0.2__). In this case PSMs with Q Values greater than __0.2__ would be removed from the input. If Q Values are not in input please use __None__ . | Numeric |
| custom_restriction | Custom Value to filter. (i.e. __5__). In this case PSMs with Custom value greater than / less than __5__ would be removed from the input. If Not using a custom score please use __None__. __NOTE__: If a higher score is "better" for your score please set __psm_score_type__ to __additive__. If a lower score is "better" please set __psm_score_type__ parameter to __multiplicative__.   | Numeric |

## Score
| Parameter | Description |Type|
|---|---|---|
| protein_score | One of any of the following: __multiplicative_log__, __best_peptide_per_protein__, __top_two_combined__, __additive__, __iterative_downweighted_log__, __downweighted_multiplicative_log__, __geometric_mean__. Recommended: __multiplicative_log__. | String |
| psm_score | PSM score to use for Protein Scoring. If using Percolator output as input this would either be __posterior_error_prob__ or __q-value__. The string typed here should match the column in your input files __EXACTLY__. If using a custom score it will be filtered accordingly with the value in [__custom_restriction__](#data-restriction). | String |
| psm_score_type | The Type of score that __psm_score__ parameter is. Either __multiplicative__ or __additive__. If a larger psm score is "better" than input additive (i.e. Mascot Ion Score, Xcorr, Percolator Score). If a smaller psm score is "better" than input multiplicative (i.e. Q Value, Posterior Error Probability). See [below](#extra-score-information) for more information.| String |

#### Extra Score information
 1. The __protein_score__, __psm_score__, and __psm_score_type__ methods must be compatible.
 2. If using a PSM score (__psm_score__ parameter) where the lower the score the better (i.e. __posterior_error_prob__ or __q-value__) then any  __protein_score__ can be used except __additive__. __psm_score_type__ must also be set to __multiplicative__.
 3. If using a PSM score (__psm_score__ parameter) where the higher the score the better (i.e. Percolator Score, Mascot Ion Score, Xcorr) (Percolator Score is called __psm_score__ - column name) in the tab delimited percolator output. Then __protein_score__ and __psm_score_type__ must both be __additive__.

## Identifiers
| Parameter | Description |Type|
|---|---|---|
| decoy_symbol | Symbol within Decoy Identifiers to distinguish between targets. (i.e "__##__" or "__decoy___"). This is important for [Protein Picker](supplementary.md#protein-picker) and FDR calculation. | String |
| isoform_symbol | Symbol that is present in isoform proteins only. (i.e. "__-__"). See [below](#extra-identifier-information) for more information. | String |
| reviewed_identifier_symbol | Identifier to determine a reviewed vs unreviewed identifier. (i.e. "__sp\|__"). See [below](#extra-identifier-information) for more information.   | String |

#### Extra Identifier information
 1. For the __decoy_symbol__ an example of a target protein -> __ex|protein__ and its decoy counterpart could be any of the following: __##ex|##protein__, __##ex|protein__, __decoy_ex|protein__. The decoy symbol just needs to be present within the string to be determined a decoy.
 2. For __isoform_symbol__ and __reviewed_identifier_symbol__, these are used to assign priority in certain algorithms such as parsimony. For example, if we have canonical proteins, isoform proteins, and reviewed/unreviewed proteins in a given analysis; the priority would be established as follows: Reviewed Canonical, Reviewed Isoform, Unreviewed. This means that if two proteins map to the same peptides, the algorithm has to make a decision on which to pick. It would use the previous mentioned priority to pick the protein lead to report. 

## Inference
| Parameter | Description |Type|
|---|---|---|
| inference_type | The Inference procedure to apply to the analysis. This can be __parsimony__, __inclusion__, __exclusion__, __peptide_centric__, or __first_protein__. Please see [here](supplementary.md#inference-types) for more information on the inference types.  | String |
| grouping_type | How to group proteins for a given __inference_type__. This can be __subset_peptides__,  __shared_peptides__, or __None__. Typically __subset_peptides__ is used. This parameter only effects grouped proteins and has no impact on protein leads. | String |

## Digest
| Parameter | Description |Type|
|---|---|---|
| digest_type | The enzyme used for digestion for the MS searches. (i.e. __trypsin__). For reference, the database digestion is handled with pyteomics. Can be any expasy rule as defined [here](https://pyteomics.readthedocs.io/en/latest/_modules/pyteomics/parser.html) other common examples include: __trypsin__, __chymotrypsin high specificity__, __chymotrypsin low specificity__, __lysc__. | String |
| missed_cleavages | The number of missed cleavages allowed for the MS searches. (i.e. __2__) | Int |

## Parsimony
These parameters are only used if __parsimony__ is selected as __inference_type__.

| Parameter | Description |Type|
|---|---|---|
| lp_solver | This can be one of: __pulp__ or __None__. This determines which linear program solver is used. Please see [here](supplementary.md#parsimony-dependencies) for more information on lp solvers. Input __None__ if not running __parsimony__. If running __parsimony__ this needs to be set to __pulp__. | String |
| shared_peptides | How to assign shared peptides for parsimony. Can be one of: __all__ or __best__. __all__ assigns shared peptides to all possible proteins in the output. __best__ assigns shared peptides to the best scoring protein which is a "winner take all" approach. This is specific to the Parsimony Inference type. | String |


## Peptide Centric
These parameters are only used if __peptide_centric__ is selected as  __inference_type__.

| Parameter | Description | Type |
|---|---|---|
| max_identifiers | The maximum number of proteins a peptide is allowed to map to. (i.e. __5__). This serves to limit the number of protein groups that can be created due to highly homologous peptides. | Int |


## Default Parameters
```yaml
parameters:
  general:
    export: peptides
    fdr: 0.01
    picker: True
    tag: py_protein_inference
  data_restriction:
    pep_restriction: 0.9
    peptide_length_restriction: 7
    q_value_restriction: 0.005
    custom_restriction: None
  score:
    protein_score: multiplicative_log
    psm_score: posterior_error_prob
    psm_score_type: multiplicative
  identifiers:
    decoy_symbol: "##"
    isoform_symbol: "-"
    reviewed_identifier_symbol: "sp|"
  inference:
    inference_type: peptide_centric
    grouping_type: shared_peptides
  digest:
    digest_type: trypsin
    missed_cleavages: 3
  parsimony:
    lp_solver: pulp
    shared_peptides: all
  peptide_centric:
    max_identifiers: 5
```
