from Digest import insilicodigest
import time
import collections

start_time = time.time()

digest = insilicodigest.InSilicoDigest(database_path='all_benchmark_entrapment_data/entrapment_data/prest_pool_ab.fasta', num_miss_cleavs=2, digest_type='trypsin')
digest.execute()

elapsed_time = time.time() - start_time

print 'time to do digest = '+str(elapsed_time)+' seconds'
