from numpy.random import dirichlet
from numpy.random import randint
import numpy as np
from numpy.random import rand

def sampler(alphabet, dist):
	sample = None
	cum_dist = np.cumsum(dist)
	r = rand()
	for i in xrange(len(dist)):
		if r < cum_dist[i]:
			sample = alphabet[i]
			break
	return sample

num_sequences = 20 		
sequence_length = 500	
motif_length = 15	

alphabet = ['A', 'C', 'T', 'G']	

alpha_b = [1, 1, 1, 1]
alpha_w = [3, 9, 1, 6]

position = [0] * num_sequences
for i in xrange(num_sequences):
	position[i] = randint(0, sequence_length - motif_length + 1)
	
cat_b = dirichlet(alpha_b)
sequences = []
for i in xrange(num_sequences):
	seq = []
	for j in xrange(sequence_length):
		seq += [sampler(alphabet, cat_b)]
	sequences += [seq]

theta = dirichlet(alpha_w, motif_length)
for i in xrange(num_sequences):
	start_pos = position[i]
	for j in xrange(motif_length):
		sequences[i][start_pos+j] = sampler(alphabet, theta[j])

filename = 'sequences.fa'

f = open(filename, 'w')

i = 0
for s in sequences:
	f.write('>seq' + str(i) + '\n')
	f.write(''.join(map(str, s)) +'\n')
	i += 1 
f.close()
	
