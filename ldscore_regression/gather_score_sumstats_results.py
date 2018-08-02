import argparse
import os
import easygrid
import re


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Script to gather stats from log and results file from an LD score regression run.')
    parser.add_argument('results_file_list', help='List of results files to gather statistics from. Note that the corresponding .log files must also be present.')
    parser.add_argument('output_file', help='Output file with stats from the LD score regression run.')
    args = parser.parse_args()

    input_files = [name.strip() for name in open(args.results_file_list)]

    results = []
    for input_file in input_files:
        log_file = input_file.replace('.results', '.log')
        line_dict = {}
        if not os.path.exists(log_file):
            raise ValueError('Log file for results file not found: %s' % log_file)

        for line in easygrid.read_delim(input_file):
            if line['Category'] == 'L2_0':
                line_dict = line
                break

        line_dict['result_file'] = input_file

        for line in open(log_file):
            if 'SNPs with chi' in line:
                snp_count = re.search('\(([0-9]+) SNPs', line).group(1)
                line_dict['snp_count'] = snp_count
            if 'Total Observed scale h2' in line:
                match = re.search(': ([\-]*[0-9]\.[0-9]+(e[\-0-9]+)?) \((0\.[0-9]+)\)', line)
                h2 = match.group(1)
                h2_se = match.group(3)
                line_dict['h2'] = h2
                line_dict['h2_se'] = h2_se
        results.append(line_dict)

    columns = []
    with open(args.output_file, 'w') as output_file:
        for result in results:
            if not columns:
                columns = result.keys()
                output_file.write('\t'.join(columns) + '\n')

            entries = '\t'.join([result[column] for column in columns])
            output_file.write(entries + '\n')
