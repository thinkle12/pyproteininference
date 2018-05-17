from Digest import insilicodigest
import time
import collections

start_time = time.time()

digest = insilicodigest.InSilicoDigest(database_path='/Users/hinklet/random_analysis/shigella_kgg_p2/Mouse_Shigella_up17_18_properformat.fasta', num_miss_cleavs=3, digest_type='trypsin')
digest.execute()

elapsed_time = time.time() - start_time

print 'time to do digest = '+str(elapsed_time)+' seconds'
