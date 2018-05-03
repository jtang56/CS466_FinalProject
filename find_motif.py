import sys
import time
from motif import Motif
from optparse import OptionParser
from Bio import SeqIO
from random import choice

def main():
	start = time.time()
	data = begin()
	g = Motif(sequences=data['sequences'], motif_length=data['length'])
	g.find_motif()
	end = time.time()
	print("Time taken: " + str(end - start))
	print_sequences(data['sequences'], data['length'])

def begin():
	parser = OptionParser(usage = "usage: %prog -i FILE -l LENGTH")

	parser.add_option("-i", "--input", dest="input", metavar="FILE")
	parser.add_option("-l", "--length", dest="length", metavar="LENGTH", type="int")
	(options, args) = parser.parse_args()

	if not options.input:
		parser.error("Input file required")

	if not options.length:
		parser.error("Motif length required")

	try:
		file = open(options.input)
	except IOError:
		parser.error("could not read file %s" % options.input)

	sequences = [{'sequence': record.seq, 'motif_position': 0}
				 for record in SeqIO.parse(file, "fasta")]

	if len(sequences) < 2:
		parser.error("found %i sequences in input file %s" % (len(sequences),
															  options.input))

	return {'sequences': sequences, 'length': options.length}


def print_sequences(sequences, length):
	for i in range(len(sequences)):
		start, end = (sequences[i]['motif_position'], sequences[i]['motif_position'] + length)
		print "#%2i  %s : %i" % (i + 1, sequences[i]['sequence'][start:end], sequences[i]['motif_position'] + 1)

if __name__ == "__main__":
	main()
