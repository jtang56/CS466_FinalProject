from numpy.random import dirichlet
from numpy.random import randint
import numpy as np
from numpy.random import rand
#from sampling import sample

def sampler(alphabet, dist):
	sample = None
	cum_dist = np.cumsum(dist)
	r = rand()
	for i in xrange(len(dist)):
		if r < cum_dist[i]:
			sample = alphabet[i]
			break
	return sample

K = 10
N = 150
w = 6

alphabet = ['A', 'C', 'T', 'G']	
M = len(alphabet)

alpha_b = [1] * M
alpha_w = [10,2,8,3]

position = [0] * K
for i in xrange(K):
	position[i] = randint(0, N - w+1)
	
cat_b = dirichlet(alpha_b)
sequences = []
for i in xrange(K):
	seq = []
	for j in xrange(N):
		seq += [sampler(alphabet, cat_b)]
	sequences += [seq]

theta = dirichlet(alpha_w, w)
for i in xrange(K):
	start_pos = position[i]
	for j in xrange(w):
		sequences[i][start_pos+j] = sampler(alphabet, theta[j])

filename = 'sequences.fa'

f = open(filename, 'w')

i = 0
for s in sequences:
	f.write('>seq' + str(i) + '\n')
	f.write(''.join(map(str, s)) +'\n')
	i += 1 
f.close()
	
