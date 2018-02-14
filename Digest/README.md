# in silico digestion

Main Usage: As a secondary module in ProteinInference

Note:

1. The script is written under Python 2.7
2. Biopython is prerequisite
3. The trypsin digestion script shared here follows proline rule, which means it does not cut lysine (K) or arginine (R) if they are followed by proline (P).

##Example Python Runner
```python
from Digest import insilicodigest

#Call InSilicoDigest, specify the database, number of missed cleavages and digest type
#Currently only "trypsin" is supported
digest = insilicodigest.InSilicoDigest(database_path='data/UniprotKBConcat1708_HUMAN.fasta', num_miss_cleavs=2, digest_type='trypsin')
digest.execute()
```