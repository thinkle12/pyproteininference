## Input File Examples
The standard input filetype is the tab delimited output from the percolator algorithm. Please see below for examples of input files:
### Standard Percolator Output

#### Target Output
| PSMid | score | q-value | posterior_error_prob | peptide | proteinIds |  |  |  |
|---|---|---|---|---|---|---|---|---|
| 1.1 | 7.5 | 0.0048 | 0.0007 | R.NYIQSLTQMPK.M | MK14_HUMAN\|Q16539 | MK14_HUMAN\|Q16539-2 | MK14_HUMAN\|Q16539-3 |  |
| 1.2 | 6.2 | 0.0035 | 0.0006 | R.NTVASSSRSM*R.T | FHDC1_HUMAN\|Q9C0D6 |  |  |  |


#### Decoy Output
| PSMid | score | q-value | posterior_error_prob | peptide | proteinIds |  |  |  |
|---|---|---|---|---|---|---|---|---|
| 1.1 | 1.3 | 0.18 | 0.27 | R.RSTITSRE.M | decoy_MK14_HUMAN\|Q16539 | decoy_MK14_HUMAN\|Q16539-2 | decoy_MK14_HUMAN\|Q16539-3 |  |
| 1.2 | 0.9 | 0.35 | 0.36 | R.KKRKRSRKEM*R.T | decoy_FHDC1_HUMAN\|Q9C0D6 |  |  |  |

Decoy proteins should have some sort of decoy identifier to distinguish between target and decoy proteins. This is typically "decoy_" or "##". See the decoy symbol parameter option [here](parameters.md#identifiers) for more information. 
These can also be combined into one file called "Combined Output".

With the above standard input one could use __q-value__ or __posterior_error_prob__ as the PSM score. See [Score Section](parameters.md#score) of the parameter file explanation with __multiplicative__ as __psm_score_type__ and any of the multiplicative options for __protein_score__.

For example standard input files please see any of the following files from the our repository:

- `tests/data/test_perc_data_target.txt`
- `tests/data/test_perc_data_decoy.txt`

### Custom Input
| PSMid | custom_score | peptide | proteinIds |  | 
|---|---|---|---|---|
| 1.1 | 7.5 | R.NYIQSLTQMPK.M | MK14_HUMAN\|Q16539 | MK14_HUMAN\|Q16539-2 | MK14_HUMAN\|Q16539-3 |  |
| 1.2 | 6.2 |  R.NTVASSSRSM*R.T | FHDC1_HUMAN\|Q9C0D6 |  |  | 

With the above custom input one could use __custom_score__ as the PSM __psm_score__ with __additive__ as the __psm_score_type__ and __protein_score__.

For example custom input files please see any of the following files from the our repository:

- `tests/data/test_perc_data_target_additive.txt`
- `tests/data/test_perc_data_decoy_additive.txt`
- `tests/data/test_perc_data_target_multiplicative.txt`
- `tests/data/test_perc_data_decoy_multiplicative.txt`

### Fasta File
This package was developed using standard Fasta files from [Uniprot](https://www.uniprot.org/).
Please see an example entry in a Fasta database below:
```text
>sp|Q5QNW6|H2B2F_HUMAN Histone H2B type 2-F OS=Homo sapiens OX=9606 GN=H2BC18 PE=1 SV=3
MPDPAKSAPAPKKGSKKAVTKVQKKDGKKRKRSRKESYSVYVYKVLKQVHPDTGISSKAM
GIMNSFVNDIFERIAGEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVT
KYTSSK
```