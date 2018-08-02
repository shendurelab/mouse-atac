import os
import argparse
import gzip
import subprocess

# Stage names -- referenced in lists of dependencies and name annotations in various steps
MAKE_TEMPLATE_FILES = 'make_template_files'
TRAIN_MODELS = 'train_models'
SCORE_SUMSTATS = 'score_sumstats'
GATHER_SUMSTATS = 'gather_sumstats'


def check_exists(path):
    """
    Utility function to check if a file exists and throw and error if not.

    Args:
        path (str): path to check

    Raises:
        ValueError if path does not exist.

    """
    if not os.path.exists(path):
        raise ValueError('Path provided as input was not found: %s' % path)


def mkdir(directory):
    """
    Simple utility function to make directories that don't already exist.

    Args:
        directory (str): directory to make

    Modifies:
        directory is made on disk unless it already exists

    """
    if not os.path.exists(directory):
        os.mkdir(directory)


def read_delim(file_path, header=True, columns=None, types=None, sep='\t'):
    """
    Utility function for parsing delimited files into iterators of dictionaries with keys as column names.

    Args:
        file_path (str): path to file
        header (bool): True if column names are present in file and false otherwise
        columns (list): List of keys to use for dictionaries as column names
        types (dict): a dict with column names as keys and conversion functions such as int, float, etc. as values
        sep (str): delimiter used in file

    Yields:
        dict: yields dict where keys are column names and values are values for a given row in file

    """
    # Choose gz open if needed
    if file_path.endswith('.gz'):
        open_function = gzip.open
    else:
        open_function = open

    # Vlalidate columns
    if columns and not isinstance(columns, list):
        raise ValueError('columns argument must be a list.')

    # Validate types and column headers
    if types and not isinstance(types, dict):
        raise ValueError('types argument must be a dict with column names as keys and a conversion function as values (such as str, int, float)')

    fh = open_function(file_path)

    if header:
        if not columns:
            columns = fh.next().strip().split(sep)
        else:
            fh.next()
    elif not columns:
        raise ValueError('Header argument specified as False, so must provide column names.')

    if types:
        typed_columns = set(types.keys())
        if False in [column in typed_columns for column in columns]:
            raise ValueError('Provided a type for %s column in types argument, but column does not appear in file or columns argument.' % column)

    # Parse file
    for line_number, line in enumerate(fh):
        entries = line.strip().split(sep)

        if len(entries) != len(columns):
            raise ValueError('Length of entries: %s, not equal to length of columns %s, in file: %s, line number %s' % (len(entries), len(columns), file_path, line_number))

        entries_dict = dict(zip(columns, entries))

        if types:
            for column in types:
                try:
                    entries_dict[column] = types[column](entries_dict[column])
                except ValueError:
                    raise ValueError('Type conversion of column %s failed on line %s' % (column, line_number))

        yield entries_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Script to generate LD score regression scores and plot results for a given set of genomic intervals.')
    parser.add_argument('--output_directory', default='.', help='Directory in which pipeline outputs are saved.')
    parser.add_argument('--sample_sheet', required=True, help='TSV file with headers sample_id, sites, category')
    parser.add_argument('--sumstats', required=True, help="TSV file with headers phenotype and sumstats for a phenotype description column and sumstats file for use with LD score regression.")
    parser.add_argument('--master_peaks', required=True, help='BED file with set of peaks that make up all possible intervals that the specific set can draw from.')
    parser.add_argument('--liftover_path', required=True, help='liftOver tool path.')
    parser.add_argument('--liftover_chain', required=True, help='Liftover chain to convert SNPs in sumstats files to the genome used for the intervals.')
    parser.add_argument('--ldsc_path', required=True, help='Path to LDSC.py.')
    parser.add_argument('--score_sumstats_options', help='String containing options for use with LDSC in scoring sumstats files. Must start with a space and be quoted. Example: " --chisq-max 99999999999"')
    parser.add_argument('--baseline_prefix', required=True, help='Prefix for baseline model file to use with LD score regression.')
    parser.add_argument('--annot_template_prefix', required=True, help='Prefix for .annot.gz files to use as templates for new annotation files.')
    parser.add_argument('--bfile_prefix', required=True, help="Prefix for files used in --bfile argument to LD score regression.")
    parser.add_argument('--hapmap_snps_prefix', required=True, help="Prefix for files used in --hapmap_snps argument to LD score regression.")
    parser.add_argument('--ld_score_prefix', required=True, help='Prefix for --w-ld-chr argument to LD score regression.')
    parser.add_argument('--frqfile_prefix', required=True, help='Prefix for --frqfile-chr argument used with LD score regression.')
    parser.add_argument('--dry', action='store_true', help='Set flag to perform dry run.')
    parser.add_argument('--queue', '-q', help='Specify a queue to run jobs on.')
    args = parser.parse_args()

    script_path = os.path.realpath(__file__)

    # Check important inputs
    check_exists(args.sample_sheet)
    check_exists(args.sumstats)
    check_exists(args.liftover_chain)
    check_exists(args.master_peaks)

    # Construct output directories
    make_template_files_dir = os.path.join(args.output_directory, MAKE_TEMPLATE_FILES)
    train_models_dir = os.path.join(args.output_directory, TRAIN_MODELS)
    score_sumstats_dir = os.path.join(args.output_directory, SCORE_SUMSTATS)

    mkdir(args.output_directory)
    mkdir(make_template_files_dir)
    mkdir(train_models_dir)
    mkdir(score_sumstats_dir)

    # Parse input files
    samples = list(read_delim(args.sample_sheet))
    sumstats = list(read_delim(args.sumstats))

    # Check to make sure all files specified exist
    for sample in samples:
        check_exists(sample['sites'])

    for sumstat in sumstats:
        check_exists(sumstat['sumstats'])

    for chromosome in range(1, 23):
        # Make a new template file for each chromosome with coordinates lifted over and such
        original_template = '%s.%s.annot.gz' % (args.annot_template_prefix, chromosome)
        check_exists(original_template)

        template_file = os.path.join(make_template_files_dir, os.path.basename(original_template))
        script = os.path.join(script_path, 'generate_ld_score_annotation_template.py')
        command = 'python %s %s %s --peaks %s --buffer 100 --liftover_chain %s --liftover_path --ldsc_path %s' % (script, original_template, template_file, args.master_peaks, args.liftover_chain, args.liftover_path, args.ldsc_path)

        # NOTE these are not required here, I have just left them here to give an idea of the order of execution, resource requests, and outputs
        outputs = [template_file]
        memory = '10G'
        name = MAKE_TEMPLATE_FILES

        subprocess.check_call(command, shell=True)  # NOTE we normally parallelize, this is just for illustrative purposes

        # Now generate a model for each sample/chromosome using the templates
        for sample in samples:
            bfile = '%s.%s' % (args.bfile_prefix, chromosome)

            hapmap_file = '%s.%s.snp' % (args.hapmap_snps_prefix, chromosome)
            check_exists(hapmap_file)

            model_prefix = os.path.join(train_models_dir, '%s.%s' % (sample['sample_id'], chromosome))
            M_file = '%s.l2.M' % (model_prefix)
            annot_file = '%s.annot.gz' % (model_prefix)
            ldscore_file = '%s.l2.ldscore.gz' % (model_prefix)
            M50_file = '%s.l2.M_5_50' % (model_prefix)
            script = os.path.join(script_path, 'generate_cluster_ldsc_model.py')

            command = 'python %s %s %s %s --bfile %s --hapmap_snps %s' % (script, template_file, sample['sites'], model_prefix, bfile, hapmap_file)

            # NOTE these are not required here, I have just left them here to give an idea of the order of execution, resource requests, and outputs
            name = TRAIN_MODELS
            dependencies = [MAKE_TEMPLATE_FILES]
            outputs = [M_file, annot_file, ldscore_file, M50_file]
            memory = '3G'

            subprocess.check_call(command, shell=True)  # NOTE we normally parallelize, this is just for illustrative purposes

    # Score each sample for every sumstats file
    results_files = []
    for sample in samples:
        model_prefix = os.path.join(train_models_dir, '%s.' % sample['sample_id'])

        for sumstat in sumstats:
            score_prefix = os.path.join(score_sumstats_dir, '%s-%s' % (sample['sample_id'], sumstat['phenotype']))
            results_file = '%s.results' % score_prefix
            results_files.append(results_file)
            command = 'python %s --h2 %s --w-ld-chr %s --ref-ld-chr %s,%s --overlap-annot --frqfile-chr %s --out %s --print-coefficients' % (args.ldsc_path, sumstat['sumstats'], args.ld_score_prefix, model_prefix, args.baseline_prefix, args.frqfile_prefix, score_prefix)

            if args.score_sumstats_options:
                command += args.score_sumstats_options

            # NOTE these are not required here, I have just left them here to give an idea of the order of execution, resource requests, and outputs
            name = SCORE_SUMSTATS
            dependencies = [TRAIN_MODELS]
            memory = '8G'
            outputs = [results_file]

            subprocess.check_call(command, shell=True)  # NOTE we normally parallelize, this is just for illustrative purposes

    # Gather the results file stats
    results_file_list = os.path.join(args.output_directory, 'results_file_list.txt')
    results_file_gathered = os.path.join(args.output_directory, 'results_gathered.txt')

    with open(results_file_list, 'w') as results_out:
        for result_file in results_files:
            results_out.write(result_file + '\n')

    script = os.path.join(script_path, 'gather_score_sumstats_results.py')
    command = 'python %s %s %s' % (script, results_file_list, results_file_gathered)

    # NOTE these are not required here, I have just left them here to give an idea of the order of execution, resource requests, and outputs
    name = GATHER_SUMSTATS
    dependencies = [SCORE_SUMSTATS]
    memory = '5G'
    outputs = [results_file_gathered]

    subprocess.check_call(command, shell=True)  # NOTE we normally parallelize, this is just for illustrative purposes
