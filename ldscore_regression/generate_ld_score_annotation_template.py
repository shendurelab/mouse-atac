import subprocess
import argparse
from tempfile import NamedTemporaryFile
from collections import OrderedDict
import gzip

############################################################################
# Functions for calling liftover to remap SNP coordinates from HG19 to MM9
############################################################################


def parse_liftover_output(file_object):
    """
    Parses liftover output
    """
    columns = ['chrom', 'start', 'end', 'id']
    results = {}
    for line in file_object:
        entries = line.strip().split('\t')
        rsid = entries[3]
        results[rsid] = tuple(entries[0:4])
    return results


def liftover(coordinates, liftover_path):
    """
    Uses liftOver tool to liftover coordinates from one genome to the next.
    """
    # Write out a temporary BED file with coordinates
    temp_file = NamedTemporaryFile(delete=False, mode='w')

    for chrom, start, end, id in coordinates:
        temp_file.write('\t'.join([chrom, str(start), str(end), id]) + '\n')
    temp_file.close()

    # Also get temp files for liftover to use
    temp_success = NamedTemporaryFile(delete=False)
    temp_fail = NamedTemporaryFile(delete=False)
    temp_success.close()
    temp_fail.close()

    # Call liftover
    command = "%s %s %s %s %s -bedPlus=3 -tab" % (liftover_path, temp_file.name, chain_file, temp_success.name, temp_fail.name)
    subprocess.call(command, shell=True)

    # Get the new coordinates from the file
    return parse_liftover_output(open(temp_success.name))


############################################################################
# Functions for dealing with intersections
############################################################################
def parse_master_peak_overlaps(file_object):
    columns = ['PEAK_CHR', 'PEAK_START', 'PEAK_STOP', 'BUFFER_CHR', 'BUFFER_START', 'BUFFER_STOP', 'RSID']

    results = []

    for line in file_object:
        entries = line.strip().split('\t')
        entries_dict = dict(zip(columns, entries))
        results.append(entries_dict)
    return results


def get_master_peak_overlaps(output_snps, peaks, buffer=100):
    temp_input_file = NamedTemporaryFile(delete=False, mode='w')
    temp_output_file = NamedTemporaryFile(delete=False, mode='w')
    temp_output_file.close()

    for snp in output_snps:
        # Skip any SNPs that failed liftover
        if output_snps[snp]['LIFTOVER'] == "1":
            continue

        # Otherwise add buffer and get intersection with ATAC peak set
        output_snps[snp]['BUFFER_START'] = str(max(0, int(output_snps[snp]['MM9_START']) - buffer))
        output_snps[snp]['BUFFER_STOP'] = str(int(output_snps[snp]['MM9_STOP']) + buffer)

        temp_input_file.write('\t'.join([output_snps[snp]['MM9_CHR'], output_snps[snp]['BUFFER_START'], output_snps[snp]['BUFFER_STOP'], output_snps[snp]['RSID']]) + '\n')
    temp_input_file.close()

    command = 'bedtools intersect -a %s -b %s -wb -wa > %s' % (peaks, temp_input_file.name, temp_output_file.name)
    subprocess.call(command, shell=True)
    return parse_master_peak_overlaps(open(temp_output_file.name))

############################################################################
# Utility functions
############################################################################


def parse_input_file(file_object):
    """
    Parses input file that contains RSIDs and their p-values
    """
    columns = ['CHR', 'START', 'RSID', 'CM']  # ignore the last annotation column

    results = OrderedDict()
    file_object.readline()  # skip header
    for line in file_object:
        entries = line.strip().split('\t')
        entries_dict = dict(zip(columns, entries))
        entries_dict['STOP'] = str(int(entries_dict['START']) + 1)
        results[entries[2]] = entries_dict

    return results


def merge_snp_metadata(rsids, liftover_output):
    """
    Utility function to merge metadata for each SNP prior to getting final intersections
    """

    rsids_copy = rsids

    for snp in rsids_copy:
        rsids_copy[snp]['LIFTOVER'] = "1"
        rsids_copy[snp]['MM9_CHR'] = "NA"
        rsids_copy[snp]['MM9_START'] = "NA"
        rsids_copy[snp]['MM9_STOP'] = "NA"

    for snp in liftover_output:
        chrom, start, end, id = liftover_output[snp]

        rsids_copy[snp]['LIFTOVER'] = "0"
        rsids_copy[snp]['MM9_CHR'] = chrom
        rsids_copy[snp]['MM9_START'] = start
        rsids_copy[snp]['MM9_STOP'] = end

    return rsids_copy


def merge_intersections(rsids, peak_intersections):
    rsids_copy = rsids

    for snp in rsids_copy:
        rsids_copy[snp]['INTERSECTS_PEAK'] = "0"
        rsids_copy[snp]['PEAK_CHR'] = "NA"
        rsids_copy[snp]['PEAK_START'] = "NA"
        rsids_copy[snp]['PEAK_STOP'] = "NA"

    for snp in peak_intersections:
        snp_id = snp['RSID']

        rsids_copy[snp_id]['INTERSECTS_PEAK'] = "1"
        rsids_copy[snp_id]['PEAK_CHR'] = snp['PEAK_CHR']
        rsids_copy[snp_id]['PEAK_START'] = snp['PEAK_START']
        rsids_copy[snp_id]['PEAK_STOP'] = snp['PEAK_STOP']

    return rsids_copy


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Script to take template annotation file from LDSC, liftover coordinates, and then annotate with intersections with an ATAC peak set to serve as starting template for LDSC.')
    parser.add_argument('input_file', help='Template annotation files for LDSC (cell_type_group* files.')
    parser.add_argument('output_file', help='Output file with overlapped peaks the p-value of the overlapping SNP.')
    parser.add_argument('--peaks', required=True, help='Peak file. Only SNPs that overlap these peaks are output.')
    parser.add_argument('--buffer', default=100, type=int, help='Slop used on either side of SNP for peak overlaps.')
    parser.add_argument('--liftover_path', required=True, help='liftOver tool path.')
    parser.add_argument('--liftover_chain', required=True, help='Chain file for liftOver.')
    args = parser.parse_args()

    # Get RSID and corresponding positions (and genetic map position)
    rsids = parse_input_file(gzip.open(args.input_file, 'rt'))

    # Liftover all SNPs
    liftover_input = []

    for snp in rsids:
        liftover_input.append(('chr' + rsids[snp]['CHR'], rsids[snp]['START'], rsids[snp]['STOP'], snp))

    liftover_output = liftover(liftover_input, args.liftover_path)

    # Merge metadata for all SNPs (including proxy snps)
    merged_snp_metadata = merge_snp_metadata(rsids, liftover_output)

    # Write out final file with the peaks that are overlapped and info about the overlapping SNPs
    master_peak_intersections = get_master_peak_overlaps(merged_snp_metadata, args.peaks, buffer=args.buffer)

    # Finally, merge with
    final_rsids = merge_intersections(rsids, master_peak_intersections)

    with gzip.open(args.output_file, 'wt') as output_file:
        columns = ['CHR', 'BP', 'SNP', 'CM', 'LIFTOVER', 'INTERSECTS_PEAK', 'PEAK_CHR', 'PEAK_START', 'PEAK_STOP']
        output_file.write('\t'.join(columns) + '\n')

        for snp in final_rsids:
            output_entries = [final_rsids[snp]['CHR'], final_rsids[snp]['START'], final_rsids[snp]['RSID'], final_rsids[snp]['CM'], final_rsids[snp]['LIFTOVER'], final_rsids[snp]['INTERSECTS_PEAK'], final_rsids[snp]['PEAK_CHR'], final_rsids[snp]['PEAK_START'], final_rsids[snp]['PEAK_STOP']]
            output_file.write('\t'.join(output_entries) + '\n')
