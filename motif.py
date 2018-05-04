import random
import math

class Motif:
	seq = []
	length = 0
	counts = {}

	def __init__(self, seq, length, pseudocounts_weight=0.1):
		self.seq = seq
		self.length = length
		self.calculate_pseudocounts(pseudocounts_weight)

	def find_motif(self):
		self.initial_motif_positions(self.length)
		best_entropy = 0
		best_alignment = [0] * len(self.seq)

		i = 0
		j = 0
		while i < 50:
			i += 1
			j += 1

			curr_seq = random.choice(self.seq)
			leave_out = filter(lambda s: s != curr_seq, self.seq)
			pwm = self.calculate_pwm(leave_out)
			self.calculate_position(pwm, curr_seq)
			entropy = self.calculate_entropy(pwm)
			if entropy > best_entropy:
				best_entropy = entropy
				best_alignment = [s['motif_position'] for s in self.seq]
				i = 0

		for i in range(len(self.seq)):
			self.seq[i]['motif_position'] = best_alignment[i]

	def calculate_pseudocounts(self, weight):
		total_nucleotides = float(sum([len(s['sequence']) for s in self.seq]))

		for nucleotide in ["A", "C", "T", "G"]:
			n = sum([s['sequence'].count(nucleotide) for s in self.seq])
			self.counts[nucleotide] = weight * (n / total_nucleotides)

	def initial_motif_positions(self, length):
		smallest_frequency = self.counts['A']
		smallest_nucleotide = 'A'
		for nucleotide in ["A", "C", "T", "G"]:
			if self.counts[nucleotide] < smallest_frequency:
				smallest_frequency = self.counts[nucleotide]
				smallest_nucleotide = nucleotide

		for s in self.seq:
			positions = []

			for r in range(len(s['sequence']) - length + 1):
				n = s['sequence'][r : r + length].count(smallest_nucleotide)
				positions.append(r)

			if len(positions) > 0:
				position = random.choice(positions)
				position += length / 2
				position -= self.length / 2
				if position < 0:
					position = 0
				if position > len(s['sequence']) - self.length:
					position = len(s['sequence']) - self.length

			s['motif_position'] = position

	def calculate_pwm(self, sequences):
		motif = {}
		for nucleotide in ["A", "C", "T", "G"]:
			motif[nucleotide] = [self.counts[nucleotide]] * self.length

		pseudocounts_total = sum(self.counts.values())
		for i in ["A", "C", "T", "G"]:
			for j in range(self.length):
				for s in sequences:
					position = s['motif_position'] + j
					if s['sequence'][position] == i:
						motif[i][j] += 1

				motif[i][j] /= float(len(sequences)) + pseudocounts_total
		return motif

	def calculate_position(self, motif, sequence):
		probabilities = []
		for r in range(len(sequence['sequence']) - self.length + 1):
			p_motif = 1
			p_background = 1
			for x in range(self.length):
				p_motif *= motif[sequence['sequence'][r + x]][x]
				p_background *= self.counts[sequence['sequence'][r + x]]

			probabilities.append(p_motif / p_background)

		sequence['motif_position'] = self.pick_random(probabilities)

	def calculate_entropy(self, motif):
		entropy = 0
		for nucleotide in ["A", "C", "T", "G"]:
			for i in range(self.length):
				entropy += motif[nucleotide][i] * math.log(motif[nucleotide][i] / self.counts[nucleotide], 2)

		return entropy

	def pick_random(self, distribution):
		total = sum(distribution)
		distribution = map(lambda n: float(n) / total, distribution)
		position = random.random()
		current = 0
		for p in range(len(distribution)):
			if distribution[p] + current >= position:
				return p
			current += distribution[p]
