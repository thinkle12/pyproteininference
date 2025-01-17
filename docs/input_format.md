## Input File Examples

### Usage with idXML, mzIdentML, pepXML, or Percolator Output
See below for specific assistance in running pyProteinInference from either an idXML, mzIdentML, pepXML, or Percolator File. <br>

### idXML
First, lets take a look at an idXML file from the popular OpenMS framework

When inspecting an idXML you might see the following:

```xml
<PeptideHit score="6.1267e-05" sequence="KDVDDDGEEK" charge="2" aa_before="K" aa_after="E" start="85" end="94" protein_refs="PH_23231" >
    <UserParam type="string" name="target_decoy" value="target"/>
    <UserParam type="string" name="MS:1002258" value="16"/>
    <UserParam type="string" name="MS:1002259" value="18"/>
    <UserParam type="string" name="num_matched_peptides" value="514"/>
    <UserParam type="int" name="isotope_error" value="0"/>
    <UserParam type="float" name="MS:1002252" value="2.0533"/>
    <UserParam type="float" name="COMET:xcorr" value="2.0533"/>
    <UserParam type="float" name="MS:1002253" value="0.5447"/>
    <UserParam type="float" name="COMET:deltaCn" value="0.5447"/>
    <UserParam type="float" name="MS:1002255" value="825.899999999999977"/>
    <UserParam type="float" name="COMET:spscore" value="825.899999999999977"/>
    <UserParam type="float" name="MS:1002256" value="1.0"/>
    <UserParam type="float" name="COMET:sprank" value="1.0"/>
    <UserParam type="float" name="MS:1002257" value="3.3e-05"/>
    <UserParam type="string" name="protein_references" value="unique"/>
    <UserParam type="float" name="COMET:deltaLCn" value="0.0"/>
    <UserParam type="float" name="COMET:lnExpect" value="-10.319002996497794"/>
    <UserParam type="float" name="COMET:lnNumSP" value="6.242223265455165"/>
    <UserParam type="float" name="COMET:lnRankSP" value="0.0"/>
    <UserParam type="float" name="COMET:IonFrac" value="0.888888888888889"/>
    <UserParam type="float" name="expect" value="3.3e-05"/>
    <UserParam type="float" name="MS:1001492" value="0.647488"/>
    <UserParam type="float" name="MS:1001491" value="6.1267e-05"/>
    <UserParam type="float" name="MS:1001493" value="9.83689e-05"/>
</PeptideHit>
```

In this case we see that a Comet search was performed. <br>
For selecting an input PSM score for rolling up to the protein level, you insert one of these selected UserParam names into the pyProteinInference parameter file in the score section. <br>

For example, lets say you wanted to use the `Comet:xcorr` value (for Comet Xcorr a higher score is better). So your score section might look like this:

```yaml
  score:
    protein_score: additive
    psm_score: Comet:xcorr
    psm_score_type: additive
```

One of the Accession values can also be used. For example, `MS:1001493` actually means `Posterior Error Probability`. <br>
You would simply enter in `MS:1001493` to utilize that output PSM score.

So to use that value as the PSM score you might have this as your score section for your parameter file:

```yaml
  score:
    protein_score: multiplicative_log
    psm_score: MS:1001493
    psm_score_type: multiplicative
```

Remember, `MS:1001493` actually means `Posterior Error Probability`. Please refer to your specific search algorithms documentation to learn more about what each accession means.

### mzIdentML
Lets inspect a standard MSGF+ based mzIdentML file for assistance in running pyProteinInference from a a mzIdentML file. <br>

When inspecting a file you might see a spectrum identification as such:

```xml
<SpectrumIdentificationItem chargeState="3" experimentalMassToCharge="907.1316528320312" calculatedMassToCharge="907.1220092773438" peptide_ref="Pep_RAVITSITATFNAGNHDRLVSCCR" rank="1" passThreshold="true" id="SII_45025_1">
  <PeptideEvidenceRef peptideEvidence_ref="PepEv_4633556_RAVITSITATFNAGNHDRLVSCCR_113"/>
  <cvParam cvRef="PSI-MS" accession="MS:1002049" name="MS-GF:RawScore" value="48"/>
  <cvParam cvRef="PSI-MS" accession="MS:1002050" name="MS-GF:DeNovoScore" value="182"/>
  <cvParam cvRef="PSI-MS" accession="MS:1002052" name="MS-GF:SpecEValue" value="5.9486466E-10"/>
  <cvParam cvRef="PSI-MS" accession="MS:1002053" name="MS-GF:EValue" value="0.005017754"/>
  <cvParam cvRef="PSI-MS" accession="MS:1002054" name="MS-GF:QValue" value="8.469932E-4"/>
  <cvParam cvRef="PSI-MS" accession="MS:1002055" name="MS-GF:PepQValue" value="0.001055835"/>
  <userParam name="IsotopeError" value="0"/>
  <userParam name="AssumedDissociationMethod" value="HCD"/>
  <userParam name="ExplainedIonCurrentRatio" value="0.009252999"/>
  <userParam name="NTermIonCurrentRatio" value="2.9709993E-4"/>
  <userParam name="CTermIonCurrentRatio" value="0.008955899"/>
  <userParam name="MS2IonCurrent" value="4877230.5"/>
  <userParam name="NumMatchedMainIons" value="2"/>
  <userParam name="MeanErrorAll" value="0.383495"/>
  <userParam name="StdevErrorAll" value="0.29636142"/>
  <userParam name="MeanErrorTop7" value="0.383495"/>
  <userParam name="StdevErrorTop7" value="0.29636142"/>
  <userParam name="MeanRelErrorAll" value="0.383495"/>
  <userParam name="StdevRelErrorAll" value="0.29636142"/>
  <userParam name="MeanRelErrorTop7" value="0.383495"/>
  <userParam name="StdevRelErrorTop7" value="0.29636142"/>
</SpectrumIdentificationItem>
```

When selecting custom scores to use from your mzIdentML file you can select any of the `cvParam` options under SpectrumIdentificationItem. For example, we might want to use the `MS-GF:PepQValue` (`MS:1002055`). <br>
An important note when selecting custom scores from mzIdentML files is that you will actually input the `accession` listed above instead of the `name`. <br>
So, for selecting `MS-GF:PepQValue`, you would actually input `MS:1002055` into the parameter file (See below). <br>
Also, given that `MS:1002055` scores are better if they are lower we would insert it as a multiplicative score. <br>
The `score` section of the parameter file in this case might look like this:

```yaml
  score:
    protein_score: multiplicative_log
    psm_score: MS:1002055
    psm_score_type: multiplicative
```

Remember, you can use any `cvParam` accession from your mzIdentML file. Make sure to note if the selected score is better when lower or higher to calculate protein scores correctly.

### pepXML

Lets next inspect a standard MSFragger Version 4.0 based pepXML file for assistance in running pyProteinInference from a pepXML file
See the following search_hit that comes from the above mentioned pepXML

```xml
<search_hit peptide="NVPSVWRSARR" massdiff="-0.0230712890625" calc_neutral_pep_mass="1555.8896" peptide_next_aa="K" num_missed_cleavages="2" num_tol_term="2" protein_descr="Transcription initiation factor TFIID subunit 1 OS=Mus musculus OX=10090 GN=Taf1 PE=1 SV=2" num_tot_proteins="1" tot_num_ions="20" hit_rank="1" num_matched_ions="8" protein="sp|Q80UV9|TAF1_MOUSE" peptide_prev_aa="K" is_rejected="0">
    <modification_info modified_peptide="n[230]NVPSVWRSARR" mod_nterm_mass="230.17073">
    </modification_info>
    <search_score name="hyperscore" value="15.64"/>
    <search_score name="nextscore" value="14.77"/>
    <search_score name="expect" value="1.563721e+00"/>
</search_hit>
```
When selecting custom scores from a pepXML file you would use full string names for selecting a score to use as the PSM score to roll up to the protein level.
For example, in the case above the user could select either "hyperscore", "nextscore", or "expect" score.

Also, keep in mind that if a higher score is better one should select additive scoring while if a lower score is better one should select multiplicative scoring and a corresponding multiplicative scoring method.

In this case lets say we select hyperscore. The `score` section of the parameter file might look like this:
```yaml
  score:
    protein_score: additive
    psm_score: hyperscore
    psm_score_type: additive
```

Other parameters to keep in mind for your specific idXML, mzIdentML, or pepXML searches would be the decoy symbol, inference type, etc.

**NOTE:** Loading from either idXML, mzIdentML, or pepXML uses either pyOpenMS or pyteomics to read in the PSMs. As such, your search result files must be compatible with their standard readers.

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