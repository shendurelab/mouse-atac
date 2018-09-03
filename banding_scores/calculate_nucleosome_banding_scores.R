library(argparse)
library(dplyr)
library(methods)

parser = argparse::ArgumentParser(description="Script to merge insert size file with TSNE coordinates and add other stats.")
parser$add_argument('insert_sizes', help='Input file with insert size estimates.')
parser$add_argument('output_file', help='Output file with all cells and their banding score.')
parser$add_argument('--barcodes', help='File with list of valid cell barcodes. One barcode per line.')
args = parser$parse_args()

# FFT-based nucleosome pattern scoring
get_banding_score = function(cell_subset) {
	size_range=0:1000

	cell_subset.cleaned = cell_subset[, c("cell", "insert_size", "read_count")]
	missing_rows = do.call(rbind, lapply(size_range[! size_range %in% cell_subset.cleaned$insert_size], function(x) { data.frame(cell=cell_subset.cleaned[1, "cell"], insert_size=x, read_count=0)}))


	cell_subset.cleaned = rbind(cell_subset.cleaned, missing_rows) %>% arrange(insert_size)

	periodogram = spec.pgram(cell_subset.cleaned$read_count / max(cell_subset.cleaned$read_count), pad=0.3, tap=0.5, span=2, plot=F, fast=T)

	periodogram$freq = 1 / periodogram$freq

	banding_score = sum(periodogram$spec[periodogram$freq >= 100 & periodogram$freq <= 300])

	return(data.frame("cell"=cell_subset[1, "cell"], "banding_score"=banding_score))
}

# Read data
insert_sizes = read.delim(args$insert_sizes)

# Subset only to specified cells if requested
if (!is.null(args$barcodes)) {
	master_cells = as.character(read.table(args$barcodes)[,1])
	insert_sizes = subset(insert_sizes, cell %in% master_cells)
}

insert_sizes$cell = factor(insert_sizes$cell)

## Split each cell into it's own dataset and then calculate a score for each
insert_sizes = insert_sizes %>% arrange(cell)
indices = 1:nrow(insert_sizes)
indices_by_cell = split(indices, insert_sizes$cell)
nucleosome_banding_scores = do.call(rbind, lapply(1:length(indices_by_cell), function(i) {
								print(i);
								index_set = indices_by_cell[[i]];
								data = insert_sizes[index_set, ]
								get_banding_score(data)
								}))

write.table(nucleosome_banding_scores, quote=F, row.names=F, file=args$output_file, sep=',')
