import pysam
import argparse
import collections

if __name__ == '__main__':
	parser = argparse.ArgumentParser('Script to get insert size and total count per cell distribution in long format.')
	parser.add_argument('bam_file', help='BAM file with reads for each cell. Read names must be in the format "cellname:other_text" in order to allow cell IDs to be extracted.')
	parser.add_argument('output_file', help='Output file with insert size histograms for each cell.')
	parser.add_argument('--barcodes', help='Input file containing one cell ID per line. Only cell IDs included in this barcode list will be included in output files.')
	args = parser.parse_args()

	sam_file = pysam.Samfile(args.bam_file, "rb")

	insert_sizes = {}

	read_names = set()

	# Read in the whitelist if provided
	barcodes = None
	if args.barcodes:
		barcodes = set([line.strip() for line in open(args.barcodes)])

	for alignment in sam_file:
		# Don't count R1 and R2 twice
		if alignment.query_name in read_names:
			continue
		else:
			read_names.add(alignment.query_name)

		# Get the insert size
		insert_size = abs(alignment.tlen)

		# Discount inferred insert sizes over 1000
		if insert_size >= 1000:
			continue

		# Track stats
		cell_name = alignment.query_name.split(':')[0]

		if barcodes and cell_name not in barcodes:
			continue

		cell_insert_size_distribution = insert_sizes.get(cell_name, collections.Counter())

		cell_insert_size_distribution[insert_size] += 1
		insert_sizes[cell_name] = cell_insert_size_distribution

	# Output file with stats
	with open(args.output_file, 'w') as output_file:
		output_file.write('\t'.join(["cell", "insert_size", "read_count"]) + '\n')

		for cell_name in insert_sizes:
			counts = insert_sizes[cell_name]
			for insert_size in counts:
				output_file.write('\t'.join([cell_name, str(insert_size), str(counts[insert_size])]) + '\n')

