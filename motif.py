import random
import math

class Motif:
	sequences = []
	motif_length = 0
	pseudocounts = {}

	def __init__(self, sequences, motif_length, pseudocounts_weight=0.1):
		self.sequences = sequences
		self.motif_length = motif_length
		self.calculate_pseudocounts(pseudocounts_weight)

	def find_motif(self, iterations=50):
		initial_pattern_width = self.motif_length
		self.initial_motif_positions(initial_pattern_width)
		best_entropy = 0
		best_alignment = [0] * len(self.sequences)

		i = 0
		j = 0

		while i < iterations:
			i += 1
			j += 1

			curr_seq = random.choice(self.sequences)
			sequences_minus_current = filter(
				lambda s: s != curr_seq,
				self.sequences)

			pwm = self.calculate_pwm(sequences_minus_current)

			self.calculate_position(pwm, curr_seq)

			entropy = self.calculate_entropy(pwm)
			if entropy > best_entropy:
				best_entropy = entropy
				best_alignment = [s['motif_position'] for s in self.sequences]
				i = 0

		for i in range(len(self.sequences)):
			self.sequences[i]['motif_position'] = best_alignment[i]

	def calculate_pseudocounts(self, weight):
		total_nucleotides = float(sum([len(s['sequence']) for s in self.sequences]))

		for nucleotide in ["A", "C", "T", "G"]:
			n = sum([s['sequence'].count(nucleotide) for s in self.sequences])
			self.pseudocounts[nucleotide] = weight * (n / total_nucleotides)

	def initial_motif_positions(self, pattern_width):
		lowest_freq = self.pseudocounts['A']
		lowest_base = 'A'
		for base in "ATCG":
			if self.pseudocounts[base] < lowest_freq:
				lowest_freq = self.pseudocounts[base]
				lowest_base = base

		for s in self.sequences:
			positions = []

			for r in range(len(s['sequence']) - pattern_width + 1):
				n = s['sequence'][r : r + pattern_width].count(lowest_base)
				positions.append(r)

			if len(positions) > 0:
				position = random.choice(positions)
				position += pattern_width / 2
				position -= self.motif_length / 2
				if position < 0:
					position = 0
				if position > len(s['sequence']) - self.motif_length:
					position = len(s['sequence']) - self.motif_length

			s['motif_position'] = position

	def calculate_pwm(self, sequences):
		motif = {}
		for base in "ATCG":
			motif[base] = [self.pseudocounts[base]] * self.motif_length

		pseudocounts_total = sum(self.pseudocounts.values())
		for i in "ATCG":
			for j in range(self.motif_length):
				for s in sequences:
					position = s['motif_position'] + j
					if s['sequence'][position] == i:
						motif[i][j] += 1

				motif[i][j] /= float(len(sequences)) + pseudocounts_total
		return motif

	def calculate_position(self, motif, sequence):
		probabilities = []
		for r in range(len(sequence['sequence']) - self.motif_length + 1):
			p_motif = p_background = 1

			for x in range(self.motif_length):
				p_motif *= motif[ sequence['sequence'][r+x] ][x]
				p_background *= self.pseudocounts[sequence['sequence'][r+x]]

			probabilities.append(p_motif / p_background)

		sequence['motif_position'] = self.pick_random(probabilities)

	def calculate_entropy(self, motif):
		entropy = 0

		for base in "ATCG":
			for i in range(self.motif_length):
				entropy += motif[base][i] * math.log(motif[base][i] / self.pseudocounts[base], 2)

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
