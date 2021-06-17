#!/opt/conda/envs/dandelion/bin/python

import dandelion as ddl
import argparse
import pandas as pd
import numpy as np
import os


def parse_args():
    '''
    Parse command line arguments.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--meta', help='Optional metadata CSV file, header required, first column for sample ID matching folder names in the directory this is being ran in. Can have a "prefix"/"suffix" column for barcode alteration, and "individual" to provide tigger groupings that isn\'t analysing all of the samples jointly.')
    parser.add_argument('--chain', type=str, default="IG", help='Whether the data is TR or IG, as the  preprocessing pipelines differ. Defaults to "IG".')
    parser.add_argument('--file_prefix', type=str, default="filtered", help='Which set of contig files to take for the folder. For a given PREFIX, will use PREFIX_contig_annotations.csv and PREFIX_contig.fasta. Defaults to "filtered".')
    parser.add_argument('--sep', type=str, default="_", help='The separator to place between the barcode and prefix/suffix. Uses sample names as a prefix for BCR data if metadata CSV file absent and more than one sample to process. Defaults to "_".')
    parser.add_argument('--keep_trailing_hyphen_number', action='store_false', help='If passed, do not strip out the trailing hyphen number, e.g. "-1", from the end of barcodes.')
    parser.add_argument('--clean_output', action='store_true', help='If passed, remove intermediate files that aren\'t the primary output from the run reults. The intermediate files may be occasionally useful for inspection.')
    args = parser.parse_args()
    # convert loci to lower case for compatibility, and ensure it's in TR/IG
    args.chain = args.chain.lower()
    if args.chain not in ['tr', 'ig']:
        raise ValueError('Chain must be TR or IG')
    return args


def main():
    # sponge up command line arguments to begin with
    args = parse_args()

    # set up a sample list
    # do we have metadata?
    if args.meta is not None:
        # if so, read it and use the index as the sample list
        meta = pd.read_csv(args.meta, index_col=0, sep = None)  # I'm sure Krzysztof wouldn't mind me making a little change here to make the dtype of the metafile guessable :P
        samples = list(meta.index)
    else:
        # no metadata file. create empty data frame so we can easily check for column presence
        meta = pd.DataFrame()
        # get list of all subfolders in current folder and run with that
        samples = []
        for item in os.listdir('.'):
            if os.path.isdir(item):
                samples.append(item)

    # STEP ONE - ddl.pp.format_fastas() 
    # do we have a prefix/suffix?
    if 'prefix' in meta.columns:
        # process with prefix
        vals = list(meta['prefix'].values)
        ddl.pp.format_fastas(samples, prefix=vals, sep=args.sep, remove_trailing_hyphen_number=args.keep_trailing_hyphen_number, filename_prefix = args.file_prefix)
    elif 'suffix' in meta.columns:
        # process with suffix
        vals = list(meta['suffix'].values)
        ddl.pp.format_fastas(samples, suffix=vals, sep=args.sep, remove_trailing_hyphen_number=args.keep_trailing_hyphen_number, filename_prefix = args.file_prefix)
    else:
        # neither. tag with the sample names as default, if more than one sample and the data is IG
        if len(samples) > 1 and args.chain == 'ig':
            ddl.pp.format_fastas(samples, prefix=samples, sep=args.sep, remove_trailing_hyphen_number=args.keep_trailing_hyphen_number, filename_prefix = args.file_prefix)
        else:
            # no need to tag as it's a single sample.
            ddl.pp.format_fastas(samples, remove_trailing_hyphen_number=args.keep_trailing_hyphen_number, filename_prefix = args.file_prefix)

    # STEP TWO - ddl.pp.reannotate_genes()
    # no tricks here
    ddl.pp.reannotate_genes(samples, loci=args.chain, filename_prefix = args.file_prefix)

    # IG requires further preprocessing, TR is done now
    if args.chain == 'ig':
        # STEP THREE - ddl.pp.reassign_alleles()
        # do we have individual information
        if 'individual' in meta.columns:
            # run the function for each individual separately
            for ind in np.unique(meta['individual']):
                # yes, this screwy thing is needed so the function ingests it correctly, sorry
                ddl.pp.reassign_alleles([str(i) for i in meta[meta['individual'] == ind].index.values], combined_folder=ind, save_plot=True, filename_prefix = args.file_prefix)
                # remove if cleaning output - the important information is ported to sample folders already
                if args.clean_output:
                    os.system('rm -r ' + ind)
        else:
            # run on the whole thing at once
            ddl.pp.reassign_alleles(samples, combined_folder='tigger', save_plot=True, filename_prefix = args.file_prefix)
            # remove if cleaning output - the important information is ported to sample folders already
            if args.clean_output:
                os.system('rm -r tigger')

        # STEP FOUR - ddl.pp.assign_isotypes()
        # also no tricks here
        ddl.pp.assign_isotypes(samples, save_plot=True, filename_prefix = args.file_prefix)

    # at this stage it's safe to remove the per-sample dandelion/tmp folder if need be
    if args.clean_output:
        for sample in samples:
            os.system('rm -rf ' + sample + '/dandelion/tmp')


if __name__ == "__main__":
    main()
