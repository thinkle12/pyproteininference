from digest import insilicodigest
import time
import collections

start_time = time.time()

digest = insilicodigest.InSilicoDigest(database_path='data/UniprotKBConcat1708_HUMAN.fasta', num_miss_cleavs=2, digest_type='trypsin')
digest.execute()

elapsed_time = time.time() - start_time

print 'time to do digest = '+str(elapsed_time)+' seconds'


start_time2 = time.time()

digest2 = insilicodigest.PyteomicsDigest(database_path='data/UniprotKBConcat1708_HUMAN.fasta', num_miss_cleavs=2, digest_type='trypsin')
digest2.execute()

elapsed_time2 = time.time() - start_time2

print 'time to do digest2 = '+str(elapsed_time2)+' seconds'

# Try this...
# with pyteomics.fasta.read('HUMAN.fasta', parser=pyteomics.fasta.std_parsers['uniprotkb']) as r:
#     print(next(r).description)

