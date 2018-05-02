from numpy.random import dirichlet
from numpy.random import randint
import numpy as np
from numpy.random import rand
#from sampling import sample

def sample(alphabet, dist):
	sampl = None
	cum_dist = np.cumsum(dist)
	r = rand()
	for i in xrange(len(dist)):
		if r < cum_dist[i]:
			sampl = alphabet[i]
			break
	return sampl

### The properties for the generated sequence

K = 10	 # The number of sequence
N = 150  # The length for each sequence
w = 6	# The length of the magic word

alphabet = ['A', 'C', 'T', 'G']			# The alphabet used in the sequence
M = len(alphabet)						# The number of alphabet used

alpha_b = [1]*M				# The alpha parameter for dirichlet prior of background letter
alpha_w = [10,2,8,3]		 # The alpha parameter for dirichlet prior of hidden word


### Start the generator part

# First, generate the starting position of the magic word for all sequences uniformly
position = [0]*K
for i in xrange(K):
	position[i] = randint(0, N-w+1)
	
# Generate the background letters for all sequences
cat_b = dirichlet(alpha_b)
sequences = []
for i in xrange(K):
	seq = []
	for j in xrange(N):
		seq += [sample(alphabet, cat_b)]
	sequences += [seq]

# Generate the magic words
theta = dirichlet(alpha_w, w)
for i in xrange(K):
	start_pos = position[i]
	for j in xrange(w):
		sequences[i][start_pos+j] = sample(alphabet, theta[j])


filename = 'sequences.fa'

f = open(filename, 'w')

i = 0
for s in sequences:
	f.write('>seq' + str(i) + '\n')
	f.write(''.join(map(str, s)) +'\n')
	i += 1 
f.close()
	
