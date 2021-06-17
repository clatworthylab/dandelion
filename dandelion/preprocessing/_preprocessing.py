#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 17:56:02
# @Last Modified by:   Kelvin
# @Last Modified time: 2021-06-17 10:56:38

import os
import pandas as pd
from subprocess import run
from tqdm import tqdm
import multiprocessing
from joblib import Parallel, delayed
from collections import OrderedDict
from time import sleep
from ..utilities._utilities import *
from ..utilities._core import *
from ..utilities._io import *
from .external._preprocessing import (assigngenes_igblast, makedb_igblast,
                                      parsedb_heavy, parsedb_light,
                                      tigger_genotype, creategermlines)
from plotnine import (ggplot, geom_bar, geom_col, ggtitle, scale_fill_manual,
                      coord_flip, options, element_blank, aes, xlab, ylab,
                      facet_wrap, facet_grid, theme_classic, theme, annotate,
                      theme_bw, geom_histogram, geom_vline, save_as_pdf_pages)
from changeo.Gene import buildGermline
# from changeo.IO import countDbFile, getFormatOperators, readGermlines, checkFields
from changeo.IO import getFormatOperators, readGermlines, checkFields
from changeo.Receptor import AIRRSchema, ChangeoSchema, Receptor, ReceptorData
import re
import functools
try:
    from scanpy import logging as logg
except ImportError:
    pass
import numpy as np
from Bio import Align
from typing import Union, Sequence, Tuple
from os import PathLike


def format_fasta(fasta: Union[str, PathLike],
                 prefix: Union[None, str] = None,
                 suffix: Union[None, str] = None,
                 sep: Union[None, str] = None,
                 remove_trailing_hyphen_number: bool = True,
                 outdir: Union[None, str] = None,
                 filename_prefix: Union[None, str] = None):
    """
    Add prefix to the headers/contig ids in cellranger fasta and annotation file.

    Parameters
    ----------
    fasta : str
        path to fasta file.
    prefix : str, optional
        prefix to append to the headers/contig ids.
    suffix : str, optional
        suffix to append to the headers/contig ids.
    sep : str, optional
        separator after prefix or before suffix to append to the headers/contig ids.
    remove_trailing_hyphen_number : bool
        whether or not to remove the trailing hyphen number e.g. '-1' from the cell/contig barcodes.
    outdir : str, optional
        path to output location. None defaults to 'dandelion'.
    filename_prefix : str, optional
        prefix of file name preceding '_contig'. None defaults to 'filtered'.

    Returns
    -------
    Formatted fasta file with new headers containing prefix
    """
    if filename_prefix is None:
        filename_pre = 'filtered'
    else:
        filename_pre = filename_prefix

    filePath = None
    filePath = check_fastapath(fasta, filename_prefix=filename_pre)

    if filePath is None:
        raise OSError(
            'Path to fasta file is unknown. ' +
            'Please specify path to fasta file or folder containing fasta file. '
            + 'Starting folder should only contain 1 fasta file.')

    fh = open(filePath, 'r')
    seqs = {}

    if sep is None:
        separator = '_'
    else:
        separator = str(sep)

    for header, sequence in fasta_iterator(fh):
        if prefix is None and suffix is None:
            seqs[header] = sequence
        elif prefix is not None:
            if suffix is not None:
                if remove_trailing_hyphen_number:
                    newheader = str(prefix) + separator + str(header).split(
                        '_contig')[0].split('-')[0] + separator + str(
                            suffix) + '_contig' + str(header).split(
                                '_contig')[1]
                else:
                    newheader = str(prefix) + separator + str(header).split(
                        '_contig')[0] + separator + str(
                            suffix) + '_contig' + str(header).split(
                                '_contig')[1]
            else:
                if remove_trailing_hyphen_number:
                    newheader = str(prefix) + separator + str(header).split(
                        '_contig')[0].split('-')[0] + '_contig' + str(
                            header).split('_contig')[1]
                else:
                    newheader = str(prefix) + separator + str(header)
            seqs[newheader] = sequence
        else:
            if suffix is not None:
                if remove_trailing_hyphen_number:
                    newheader = str(header).split('_contig')[0].split(
                        '-')[0] + separator + str(suffix) + '_contig' + str(
                            header).split('_contig')[1]
                else:
                    newheader = str(header) + separator + str(suffix)
            else:
                newheader = str(header)
            seqs[newheader] = sequence

    fh.close()

    if os.path.isfile(filePath):
        basedir = os.path.dirname(filePath)
    elif os.path.isdir(filePath):
        basedir = os.path.dirname(filePath)
    else:
        basedir = os.getcwd()

    if outdir is None:
        out_dir = basedir.rstrip('/') + '/' + 'dandelion/'
    else:
        if not outdir.endswith('/'):
            out_dir = outdir + '/'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    out_fasta = out_dir + os.path.basename(filePath)

    fh1 = open(out_fasta, 'w')
    fh1.close()
    out = ''
    for l in seqs:
        out = '>' + l + '\n' + seqs[l] + '\n'
        Write_output(out, out_fasta)

    # format the barcode and contig_id in the corresponding annotation file too
    anno = basedir + '/' + \
        os.path.basename(filePath).replace('.fasta', '_annotations.csv')
    data = pd.read_csv(anno, dtype='object')

    if prefix is not None:
        if suffix is not None:
            if remove_trailing_hyphen_number:
                data['contig_id'] = [
                    str(prefix) + separator +
                    str(c).split('_contig')[0].split('-')[0] + separator +
                    str(suffix) + '_contig' + str(c).split('_contig')[1]
                    for c in data['contig_id']
                ]
                data['barcode'] = [
                    str(prefix) + separator + str(b).split('-')[0] +
                    separator + str(suffix) for b in data['barcode']
                ]
            else:
                data['contig_id'] = [
                    str(prefix) + separator + str(c).split('_contig')[0] +
                    separator + str(suffix) + '_contig' +
                    str(c).split('_contig')[1] for c in data['contig_id']
                ]
                data['barcode'] = [
                    str(prefix) + separator + str(b) + separator + str(suffix)
                    for b in data['barcode']
                ]
        else:
            if remove_trailing_hyphen_number:
                data['contig_id'] = [
                    str(prefix) + separator +
                    str(c).split('_contig')[0].split('-')[0] + '_contig' +
                    str(c).split('_contig')[1] for c in data['contig_id']
                ]
                data['barcode'] = [
                    str(prefix) + separator + str(b).split('-')[0]
                    for b in data['barcode']
                ]
            else:
                data['contig_id'] = [
                    str(prefix) + separator + str(c) for c in data['contig_id']
                ]
                data['barcode'] = [
                    str(prefix) + separator + str(b) for b in data['barcode']
                ]
    else:
        if suffix is not None:
            if remove_trailing_hyphen_number:
                data['contig_id'] = [
                    str(c).split('_contig')[0].split('-')[0] + separator +
                    str(suffix) + '_contig' + str(c).split('_contig')[1]
                    for c in data['contig_id']
                ]
                data['barcode'] = [
                    str(b).split('-')[0] + separator + str(suffix)
                    for b in data['barcode']
                ]
            else:
                data['contig_id'] = [
                    str(c).split('_contig')[0] + separator + str(suffix) +
                    '_contig' + str(c).split('_contig')[1]
                    for c in data['contig_id']
                ]
                data['barcode'] = [
                    str(b) + separator + str(suffix) for b in data['barcode']
                ]
        else:
            data['contig_id'] = [str(c) for c in data['contig_id']]
            data['barcode'] = [str(b) for b in data['barcode']]

    out_anno = out_dir + \
        os.path.basename(filePath).replace('.fasta', '_annotations.csv')

    data.to_csv(out_anno, index=False)


def format_fastas(fastas: Sequence,
                  prefix: Union[None, Sequence] = None,
                  suffix: Union[None, Sequence] = None,
                  sep: Union[None, str] = None,
                  remove_trailing_hyphen_number: bool = True,
                  outdir: Union[None, str] = None,
                  filename_prefix: Union[None, Sequence, str] = None):
    """
    Add prefix to the headers/contig ids in cellranger fasta and annotation file.

    Parameters
    ----------
    fastas : Sequence
        list of paths to fasta files.
    prefix : list, optional
        list of prefixes to append to headers/contig ids in each fasta file.
    suffix : str, optional
        list of suffixes to append to headers/contig ids in each fasta file.
    sep : str, optional
        separator after prefix or before suffix to append to the headers/contig ids.
    remove_trailing_hyphen_number : bool
        whether or not to remove the trailing hyphen number e.g. '-1' from the cell/contig barcodes.
    outdir : str, optional
        path to out put location. Default is None, which is 'dandelion'.
    filename_prefix : str, optional
        list of prefixes of file names preceding '_contig'. None defaults to 'filtered'.

    Returns
    -------
    Formatted fasta file with new headers containing prefix
    """
    if type(fastas) is not list:
        fastas = [fastas]
    if type(filename_prefix) is not list:
        filename_prefix = [filename_prefix]
    if all(t is None for t in filename_prefix):
        filename_prefix = [None for f in fastas]

    if prefix is not None:
        if type(prefix) is not list:
            prefix = [prefix]
        prefix_dict = dict(zip(fastas, prefix))
    if suffix is not None:
        if type(suffix) is not list:
            suffix = [suffix]
        suffix_dict = dict(zip(fastas, suffix))

    for i in tqdm(range(0, len(fastas)), desc='Formating fasta(s) '):
        if prefix is None and suffix is None:
            format_fasta(
                fastas[i],
                prefix=None,
                suffix=None,
                sep=None,
                remove_trailing_hyphen_number=remove_trailing_hyphen_number,
                outdir=outdir,
                filename_prefix=filename_prefix[i])
        elif prefix is not None:
            if suffix is not None:
                format_fasta(
                    fastas[i],
                    prefix=prefix_dict[fastas[i]],
                    suffix=suffix_dict[fastas[i]],
                    sep=sep,
                    remove_trailing_hyphen_number=remove_trailing_hyphen_number,
                    outdir=outdir,
                    filename_prefix=filename_prefix[i])
            else:
                format_fasta(
                    fastas[i],
                    prefix=prefix_dict[fastas[i]],
                    suffix=None,
                    sep=sep,
                    remove_trailing_hyphen_number=remove_trailing_hyphen_number,
                    outdir=outdir,
                    filename_prefix=filename_prefix[i])
        else:
            if suffix is not None:
                format_fasta(
                    fastas[i],
                    prefix=None,
                    suffix=suffix_dict[fastas[i]],
                    sep=sep,
                    remove_trailing_hyphen_number=remove_trailing_hyphen_number,
                    outdir=outdir,
                    filename_prefix=filename_prefix[i])
            else:
                format_fasta(
                    fastas[i],
                    prefix=None,
                    suffix=None,
                    sep=None,
                    remove_trailing_hyphen_number=remove_trailing_hyphen_number,
                    outdir=outdir,
                    filename_prefix=filename_prefix[i])


def assign_isotype(fasta: Union[str, PathLike],
                   fileformat: Literal['blast', 'changeo', 'airr'] = 'blast',
                   org: Literal['human', 'mouse'] = 'human',
                   correct_c_call: bool = True,
                   correction_dict: Union[Dict, None] = None,
                   plot: bool = True,
                   save_plot: bool = False,
                   figsize: Tuple[Union[int, float], Union[int,
                                                           float]] = (4, 4),
                   blastdb: Union[None, str] = None,
                   allele: bool = False,
                   parallel: bool = True,
                   ncpu: Union[None, int] = None,
                   filename_prefix: Union[None, str] = None,
                   verbose: bool = False):
    """
    Annotate contigs with constant region call using blastn.

    Parameters
    ----------
    fasta : str, PathLike
        path to fasta file.
    fileformat : str
        format of V(D)J file/objects. Default is 'blast'. Also accepts 'changeo' (same behaviour as 'blast') and 'airr'.
    org : str
        organism of reference folder. Default is 'human'.
    correct_c_call : bool
        whether or not to adjust the c_calls after blast based on provided primers specified in `primer_dict` option.
        Default is True.
    correction_dict : Dict, optional
        a nested dictionary contain isotype/c_genes as keys and primer sequences as records to use for correcting
        annotated c_calls. Defaults to a curated dictionary for human sequences if left as none.
    plot : bool
        whether or not to plot reassignment summary metrics. Default is True.
    save_plot : bool
        whether or not to save plot.
    figsize : Tuple[Union[int,float], Union[int,float]]
        size of figure. Default is (4, 4).
    blastdb : str, optional
        path to blast database. Defaults to `$BLASTDB` environmental variable.
    allele : bool
        whether or not to return allele calls. Default is False.
    parallel : bool
        whether or not to use parallelization. Default is True.
    ncpu : int
        number of cores to use if parallel is True. Default is all available minus 1.
    filename_prefix : str, optional
        prefix of file name preceding '_contig'. None defaults to 'filtered'.
    verbose : bool
        whether or not to print the blast command in terminal. Default is False.

    Returns
    -------
    V(D)J tsv files with constant genes annotated.
    """
    def _run_blastn(fasta, blastdb, fileformat, org, verbose):

        env = os.environ.copy()
        if blastdb is None:
            try:
                bdb = env['BLASTDB']
            except:
                raise OSError(
                    'Environmental variable BLASTDB must be set. Otherwise, please provide path to blast database'
                )
            bdb = bdb + org + '/' + org + '_BCR_C.fasta'
        else:
            env['BLASTDB'] = blastdb
            bdb = blastdb

        cmd = [
            'blastn', '-db', bdb, '-evalue', '0.001', '-max_target_seqs', '1',
            '-outfmt', '5', '-query', fasta
        ]

        blast_out = "{}/tmp/{}.xml".format(
            os.path.dirname(fasta),
            os.path.basename(fasta).split('.fasta')[0] + fileformat)

        if verbose:
            print('Running command: %s\n' % (' '.join(cmd)))
        with open(blast_out, 'w') as out:
            run(cmd, stdout=out, env=env)

    def _parse_BLAST(fasta, fileformat):
        """Parse BLAST output from output files and writes formatted output to BLAST output summary files."""
        def split_blast_file(filename):
            """
            Split blast file.

            Code adapted from:
            http://stackoverflow.com/questions/19575702/pythonhow-to-split-file-into-chunks-by-the-occurrence-of-the-header-word.

            """
            token = '<Iteration>'
            chunks = []
            current_chunk = []

            with open(filename) as fh:
                for line in fh:
                    line = line.rstrip()

                    if line.startswith(token) and current_chunk:
                        chunks.append(current_chunk[:])
                        current_chunk = []
                    if not line.startswith("Total queries"):
                        current_chunk.append(line)

                chunks.append(current_chunk)
            return (chunks)

        def extract_blast_info(line):
            line = line.split()[0]
            info = line.split(">")[1]
            info = info.split("<")[0]
            return (info)

        input_file = "{}/tmp/{}.xml".format(
            os.path.dirname(fasta),
            os.path.basename(fasta).split('.fasta')[0] + fileformat)
        output_file = "{}/tmp/{}.blastsummary.txt".format(
            os.path.dirname(fasta),
            os.path.basename(fasta).split('.fasta')[0] + fileformat)

        with open(output_file, 'w') as outfile:
            outfile.write(
                "------------------\n##{}##\n------------------\n\n#BCR#\n\n".
                format(fasta))
            # Split result file into chunks corresponding to results for each query sequence.
            if os.path.isfile(input_file):
                blast_result_chunks = split_blast_file(input_file)
                for chunk in blast_result_chunks:
                    message = False
                    for line_x in chunk:
                        line_x = line_x.strip()
                        if line_x.startswith("<Iteration_query-def>"):
                            line = line_x.split(">")[1]
                            blast_query_name = line.split("<")[0]
                        elif line_x.startswith("<Hsp_evalue>"):
                            evalue = extract_blast_info(line_x)
                            evalue = format(float(evalue), '.0e')
                        elif line_x.startswith("<Hit_accession>"):
                            C_segment = extract_blast_info(line_x)
                            if "C-REGION" or "CH1" in C_segment:
                                C_segment = C_segment.split("_")[0]
                        elif line_x.startswith("<Hsp_bit-score>"):
                            bit_score = extract_blast_info(line_x)
                        elif line_x.startswith("<Hsp_query-from>"):
                            q_start = extract_blast_info(line_x)
                        elif line_x.startswith("<Hsp_query-to>"):
                            q_end = extract_blast_info(line_x)
                        elif line_x.startswith("<Hsp_hit-from>"):
                            s_start = extract_blast_info(line_x)
                        elif line_x.startswith("<Hsp_hit-to>"):
                            s_end = extract_blast_info(line_x)
                        elif line_x.startswith("<Iteration_query-len>"):
                            query_length = extract_blast_info(line_x)
                        elif line_x.startswith("<Hsp_align-len>"):
                            align_length = extract_blast_info(line_x)
                        elif line_x.startswith("<Hsp_gaps>"):
                            gaps = extract_blast_info(line_x)
                        elif line_x.startswith("<Hsp_identity>"):
                            identity = extract_blast_info(line_x)
                        elif line_x.startswith("<Hsp_qseq>"):
                            c_qseq = extract_blast_info(line_x)
                        elif line_x.startswith("<Hsp_hseq>"):
                            c_hseq = extract_blast_info(line_x)
                        elif line_x.startswith(
                                "<Iteration_message>No hits found"):
                            message = True
                            out_string = "##{blast_query_name}##\nNo C segment found\n\n".format(
                                blast_query_name=blast_query_name)
                            outfile.write(out_string)
                        # Create output string when reaching end of BLAST
                        # iteration result (marked by </Iteration>) and write
                        # to BLAST summary file
                        elif line_x.startswith(
                                "</Iteration>") and message is not True:
                            identity_pro = float(identity) / int(
                                align_length) * 100
                            identity_pro = format(identity_pro, '.2f')
                            mismatches = int(align_length) - int(identity)
                            # Account for reversed sequences
                            if int(s_start) > int(s_end):
                                blast_query_name = "reversed|" + blast_query_name
                                x, y = int(q_start), int(q_end)
                                q_start = int(query_length) - y + 1
                                q_end = int(query_length) - x + 1
                                s_start, s_end = s_end, s_start
                            intro_string = "##{}##\nC segment:\t{}\n\n".format(
                                blast_query_name, C_segment)
                            header_string = (
                                "Segment\tquery_id\tsubject_id\t% identity\talignment length\t"
                                "mismatches\tgap opens\tgaps\tq start\tq end\ts start\ts end\t"
                                "evalue\tbit score\n")
                            tmp_out = (
                                "C\t{blast_query_name}\t{C_segment}\t{identity_pro}\t{align_length}\t"
                                "{mismatches}\tNA\t{gaps}\t{q_start}\t{q_end}\t{s_start}\t{s_end}\t"
                                "{evalue}\t{bit_score}\t{q_seq}\t{h_seq}\n\n")
                            out_string = tmp_out.format(
                                blast_query_name=blast_query_name,
                                C_segment=C_segment,
                                identity_pro=identity_pro,
                                align_length=align_length,
                                evalue=evalue,
                                mismatches=mismatches,
                                gaps=gaps,
                                q_start=q_start,
                                q_end=q_end,
                                s_start=s_start,
                                s_end=s_end,
                                bit_score=bit_score,
                                q_seq=c_qseq,
                                h_seq=c_hseq)
                            string_to_write = intro_string + header_string + out_string
                            outfile.write(string_to_write)

    def _get_C(fasta, fileformat, allele=False, parallel=True, ncpu=None):
        def _get_C_call(fasta, contig_name, fileformat, allele=False):
            blast_summary_file = "{}/tmp/{}.blastsummary.txt".format(
                os.path.dirname(fasta),
                os.path.basename(fasta).split('.fasta')[0] + fileformat)

            C_seq, C_germ, C_gene, C_ident = None, None, None, None
            C_eval, C_bitscore, C_qstart, C_qend = None, None, None, None
            with open(blast_summary_file, 'r') as input:
                for line in input:
                    if line.startswith("C\t{contig_name}".format(
                            contig_name=contig_name)) or line.startswith(
                                "C\treversed|{contig_name}".format(
                                    contig_name=contig_name)):
                        C_gene = line.split("\t")[2]
                        C_ident = line.split("\t")[3]
                        C_seq = line.split("\t")[14]
                        C_germ = line.split("\t")[15]
                        C_eval = line.split("\t")[12]
                        C_bitscore = line.split("\t")[13]
                        C_qstart = line.split("\t")[8]
                        C_qend = line.split("\t")[9]

                        if "_CH1" or "_C-REGION" in C_gene:
                            C_gene = C_gene.split("_")[0]
            if not allele:
                try:
                    C_gene = C_gene.split('*')[0]
                except:
                    pass

            C_call, C_identity, C_sequence, C_germline, C_support, C_score, C_start, C_end, = {
            }, {}, {}, {}, {}, {}, {}, {}
            C_call[contig_name] = C_gene
            C_identity[contig_name] = C_ident
            C_sequence[contig_name] = C_seq
            C_germline[contig_name] = C_germ
            C_support[contig_name] = C_eval
            C_score[contig_name] = C_bitscore
            C_start[contig_name] = C_qstart
            C_end[contig_name] = C_qend

            return (C_sequence, C_germline, C_call, C_identity, C_support,
                    C_score, C_start, C_end)

        fh = open(fasta, 'r')
        contigs = []
        for header, sequence in fasta_iterator(fh):
            contigs.append(header)
        fh.close()

        if parallel:
            if ncpu is None:
                num_cores = multiprocessing.cpu_count() - 1
            else:
                num_cores = int(ncpu)
            results = ()
            results = Parallel(
                n_jobs=num_cores
            )(delayed(_get_C_call)(fasta, c, fileformat, allele) for c in tqdm(
                contigs,
                desc='Retrieving contant region calls, parallelizing with ' +
                str(num_cores) + ' cpus '))
            # transform list of dicts to dict
            seq, germ, call, ident, support, score, start, end = {}, {}, {}, {}, {}, {}, {}, {}
            for r in range(0, len(results)):
                _seq, _germ, _call, _ident, _support, _score, _start, _end = results[
                    r]
                seq.update(_seq)
                germ.update(_germ)
                call.update(_call)
                ident.update(_ident)
                support.update(_support)
                score.update(_score)
                start.update(_start)
                end.update(_end)
        else:
            seq, germ, call, ident, support, score, start, end = {}, {}, {}, {}, {}, {}, {}, {}
            for c in tqdm(contigs, desc='Retrieving contant region calls '):
                seq[c], germ[c], call[c], ident[c], support[c], score[
                    c], start[c], end[c] = _get_C_call(fasta, c, fileformat,
                                                       allele)[c]
        return (seq, germ, call, ident, support, score, start, end)

    def _transfer_c(data, c_dict, colname):
        _data = load_data(data)
        if colname not in _data.columns:
            _data = _data.merge(pd.DataFrame.from_dict(c_dict,
                                                       orient='index',
                                                       columns=[colname]),
                                left_index=True,
                                right_index=True)
        else:
            _data[colname] = pd.Series(c_dict)
        return (_data)

    def _add_cell(data):
        _data = load_data(data)
        _data['cell_id'] = [
            c.split('_contig')[0] for c in _data['sequence_id']
        ]
        return (_data)

    aligner = Align.PairwiseAligner()

    def two_gene_correction(self, i, dictionary):
        key1, key2 = dictionary.keys()
        seq = self.loc[i, 'c_sequence_alignment'].replace('-', '')
        alignments1 = aligner.align(dictionary[key1], seq)
        alignments2 = aligner.align(dictionary[key2], seq)
        score1 = alignments1.score
        score2 = alignments2.score
        if score1 == score2:
            self.at[i, 'c_call'] = str(key1) + ',' + str(key2)
        if score1 > score2:
            self.at[i, 'c_call'] = str(key1)
        if score1 < score2:
            self.at[i, 'c_call'] = str(key2)

    def three_gene_correction(self, i, dictionary):
        key1, key2, key3 = dictionary.keys()
        seq = self.loc[i, 'c_sequence_alignment'].replace('-', '')
        alignments1 = aligner.align(dictionary[key1], seq)
        alignments2 = aligner.align(dictionary[key2], seq)
        alignments3 = aligner.align(dictionary[key3], seq)
        score1 = alignments1.score
        score2 = alignments2.score
        score3 = alignments3.score
        if score1 == score2 == score3:
            self.at[i,
                    'c_call'] = str(key1) + ',' + str(key2) + ',' + str(key3)
        elif score1 > score2 and score1 > score3:
            self.at[i, 'c_call'] = str(key1)
        elif score2 > score1 and score2 > score3:
            self.at[i, 'c_call'] = str(key2)
        elif score3 > score1 and score3 > score2:
            self.at[i, 'c_call'] = str(key3)
        elif score1 == score2 and score1 > score3:
            self.at[i, 'c_call'] = str(key1) + ',' + str(key2)
        elif score1 > score2 and score1 == score3:
            self.at[i, 'c_call'] = str(key1) + ',' + str(key3)
        elif score2 > score1 and score2 == score3:
            self.at[i, 'c_call'] = str(key2) + ',' + str(key3)

    def four_gene_correction(self, i, dictionary):
        key1, key2, key3, key4 = dictionary.keys()
        seq = self.loc[i, 'c_sequence_alignment'].replace('-', '')
        alignments1 = aligner.align(dictionary[key1], seq)
        alignments2 = aligner.align(dictionary[key2], seq)
        alignments3 = aligner.align(dictionary[key3], seq)
        alignments4 = aligner.align(dictionary[key4], seq)
        score1 = alignments1.score
        score2 = alignments2.score
        score3 = alignments3.score
        score4 = alignments4.score
        if score1 == score2 == score3 == score4:
            self.at[i, 'c_call'] = str(key1) + ',' + str(key2) + ',' + str(
                key3) + ',' + str(key4)
        elif score1 > score2 and score1 > score3 and score1 > score4:
            self.at[i, 'c_call'] = str(key1)
        elif score2 > score1 and score2 > score3 and score2 > score4:
            self.at[i, 'c_call'] = str(key2)
        elif score3 > score1 and score3 > score2 and score3 > score4:
            self.at[i, 'c_call'] = str(key3)
        elif score4 > score1 and score4 > score2 and score4 > score3:
            self.at[i, 'c_call'] = str(key4)
        elif score1 == score2 and score1 > score3 and score1 > score4:
            self.at[i, 'c_call'] = str(key1) + ',' + str(key2)
        elif score1 > score2 and score1 == score3 and score1 > score4:
            self.at[i, 'c_call'] = str(key1) + ',' + str(key3)
        elif score1 > score2 and score1 > score3 and score1 == score4:
            self.at[i, 'c_call'] = str(key1) + ',' + str(key4)
        elif score2 == score3 and score2 > score1 and score2 > score4:
            self.at[i, 'c_call'] = str(key1) + ',' + str(key3)
        elif score2 == score4 and score2 > score1 and score2 > score3:
            self.at[i, 'c_call'] = str(key2) + ',' + str(key4)
        elif score3 == score4 and score3 > score1 and score3 > score2:
            self.at[i, 'c_call'] = str(key3) + ',' + str(key4)
        elif score1 == score2 == score3 and score1 > score4:
            self.at[i,
                    'c_call'] = str(key1) + ',' + str(key2) + ',' + str(key3)
        elif score1 == score2 == score4 and score1 > score3:
            self.at[i,
                    'c_call'] = str(key1) + ',' + str(key2) + ',' + str(key4)
        elif score1 == score3 == score4 and score1 > score2:
            self.at[i,
                    'c_call'] = str(key1) + ',' + str(key3) + ',' + str(key4)
        elif score2 == score3 == score4 and score2 > score1:
            self.at[i,
                    'c_call'] = str(key2) + ',' + str(key3) + ',' + str(key4)

    def _correct_c_call(data, primers_dict=None):
        dat = data.copy()
        if primers_dict is None:
            primer_dict = {
                'IGHG': {
                    'IGHG1':
                    'GCCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCACCCTCCTCCAAGAGCACCTCTGGGGGCACAGCGGCCCTGGGC',
                    'IGHG2':
                    'GCCTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCGGCCCTGGGC',
                    'IGHG3':
                    'GCTTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCTGGGGGCACAGCGGCCCTGGGC',
                    'IGHG4':
                    'GCTTCCACCAAGGGCCCATCCGTCTTCCCCCTGGCGCCCTGCTCCAGGAGCACCTCCGAGAGCACAGCCGCCCTGGGC'
                },
                'IGHA': {
                    'IGHA1':
                    'GCATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCTGCAGCACCCAGCCAGATGGGAACGTGGTCATCGCCTGC',
                    'IGHA2':
                    'GCATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCGACAGCACCCCCCAAGATGGGAACGTGGTCGTCGCATGC'
                },
                'IGLC7': {
                    'IGLC':
                    'GTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCGCCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCATAA',
                    'IGLC7':
                    'GTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCACCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCGTAA'
                },
                'IGLC3': {
                    'IGLC':
                    'GTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCGCCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCATAA',
                    'IGLC3':
                    'GTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCACCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCATAA'
                },
                'IGLC6': {
                    'IGLC':
                    'TCGGTCACTCTGTTCCCGCCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCA',
                    'IGLC6':
                    'TCGGTCACTCTGTTCCCGCCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGCCTGA'
                }
            }
        else:
            primer_dict = primers_dict

        for i in dat.index:
            if (dat.loc[i, 'c_call'] == dat.loc[i, 'c_call']) & (
                    dat.loc[i, 'c_call'] is not None):
                for k in primer_dict:
                    if k in dat.loc[i, 'c_call']:
                        if len(primer_dict[k]) == 2:
                            two_gene_correction(dat, i, primer_dict[k])
                        elif len(primer_dict[k]) == 3:
                            three_gene_correction(dat, i, primer_dict[k])
                        elif len(primer_dict[k]) == 4:
                            four_gene_correction(dat, i, primer_dict[k])
        return (dat)

    # main function from here
    format_dict = {
        'changeo': '_igblast_db-pass',
        'blast': '_igblast_db-pass',
        'airr': '_igblast_gap'
    }

    filePath = check_filepath(fasta,
                              filename_prefix=filename_prefix,
                              endswith='.fasta')
    if filePath is None:
        raise OSError(
            'Path to fasta file is unknown. Please specify path to fasta file or folder containing fasta file.'
        )

    if verbose:
        print('Processing {} \n'.format(filePath))

    # running blast using blast
    if verbose:
        print('Running blastn \n')
    _run_blastn(filePath, blastdb, format_dict[fileformat], org, verbose)
    # parsing output into a summary.txt file
    if verbose:
        print('Parsing blast output \n')
    _parse_BLAST(filePath, format_dict[fileformat])
    # Add the c_calls to the data file
    c_seq, c_germ, c_call, c_ident, c_supp, c_scr, c_st, c_en = {}, {}, {}, {}, {}, {}, {}, {}
    c_seq, c_germ, c_call, c_ident, c_supp, c_scr, c_st, c_en = _get_C(
        filePath, format_dict[fileformat], allele, parallel, ncpu)

    _file = "{}/tmp/{}_genotyped.tsv".format(
        os.path.dirname(filePath),
        os.path.basename(filePath).split('.fasta')[0] +
        format_dict[fileformat])
    _airrfile = "{}/tmp/{}.tsv".format(
        os.path.dirname(filePath),
        os.path.basename(filePath).split('.fasta')[0] + '_igblast')
    _file2 = "{}/{}_genotyped.tsv".format(
        os.path.dirname(filePath),
        os.path.basename(filePath).split('.fasta')[0] +
        format_dict[fileformat])

    if verbose:
        print('Loading 10X annotations \n')
    dat_10x = load_data(_file)
    res_10x = pd.DataFrame(dat_10x['c_call'])
    res_10x['c_call'] = res_10x['c_call'].fillna(value='None')
    if verbose:
        print('Preparing new calls \n')
    dat = _transfer_c(_file, c_call, 'c_call')
    dat = _transfer_c(dat, c_seq, 'c_sequence_alignment')
    dat = _transfer_c(dat, c_germ, 'c_germline_alignment')
    dat = _transfer_c(dat, c_st, 'c_sequence_start')
    dat = _transfer_c(dat, c_en, 'c_sequence_end')
    dat = _transfer_c(dat, c_scr, 'c_score')
    dat = _transfer_c(dat, c_ident, 'c_identity')
    dat = _transfer_c(dat, c_supp, 'c_support')
    res_blast = pd.DataFrame(dat['c_call'])
    res_blast = res_blast.fillna(value='None')

    res_10x_sum = pd.DataFrame(res_10x['c_call'].value_counts(normalize=True) *
                               100)
    res_blast_sum = pd.DataFrame(
        res_blast['c_call'].value_counts(normalize=True) * 100)
    res_10x_sum['group'] = '10X'
    res_blast_sum['group'] = 'blast'
    res_10x_sum.columns = ['counts', 'group']
    res_blast_sum.columns = ['counts', 'group']
    res_10x_sum.index = res_10x_sum.index.set_names(['c_call'])
    res_blast_sum.index = res_blast_sum.index.set_names(['c_call'])
    res_10x_sum.reset_index(drop=False, inplace=True)
    res_blast_sum.reset_index(drop=False, inplace=True)

    if correct_c_call:  # TODO: figure out if i need to set up a None correction?
        if verbose:
            print('Correcting C calls \n')
        dat = _correct_c_call(dat, primers_dict=correction_dict)
        res_corrected = pd.DataFrame(dat['c_call'])
        res_corrected = res_corrected.fillna(value='None')
        res_corrected_sum = pd.DataFrame(
            res_corrected['c_call'].value_counts(normalize=True) * 100)
        res_corrected_sum['group'] = 'corrected'
        res_corrected_sum.columns = ['counts', 'group']
        res_corrected_sum.index = res_corrected_sum.index.set_names(['c_call'])
        res_corrected_sum.reset_index(drop=False, inplace=True)
        res = pd.concat([res_10x_sum, res_blast_sum, res_corrected_sum])
    else:
        res = pd.concat([res_10x_sum, res_blast_sum])

    res = res.reset_index(drop=True)
    res['c_call'] = res['c_call'].astype('category')
    res['c_call'] = res['c_call'].cat.reorder_categories(
        sorted(list(set(res['c_call'])), reverse=True))

    if verbose:
        print('Finishing up \n')
    if 'cell_id' not in dat.columns:
        dat = _add_cell(dat)
    dat['c_call_10x'] = pd.Series(dat_10x['c_call'])
    # some minor adjustment to the final output table
    airr_output = load_data(_airrfile)
    cols_to_merge = [
        'junction_aa_length', 'fwr1_aa', 'fwr2_aa', 'fwr3_aa', 'fwr4_aa',
        'cdr1_aa', 'cdr2_aa', 'cdr3_aa', 'sequence_alignment_aa',
        'v_sequence_alignment_aa', 'd_sequence_alignment_aa',
        'j_sequence_alignment_aa'
    ]
    for x in cols_to_merge:
        dat[x] = pd.Series(airr_output[x])
    dat.to_csv(_file2, sep='\t', index=False)

    if plot:
        options.figure_size = figsize
        if correct_c_call:
            p = (ggplot(res, aes(x='c_call', y='counts', fill='group')) +
                 coord_flip() + theme_classic() + xlab("c_call") +
                 ylab("% c calls") +
                 geom_col(stat="identity", position='dodge') +
                 scale_fill_manual(values=('#79706e', '#86bcb6', '#F28e2b')) +
                 theme(legend_title=element_blank()))
        else:
            p = (ggplot(res, aes(x='c_call', y='counts', fill='group')) +
                 coord_flip() + theme_classic() + xlab("c_call") +
                 ylab("% c calls") +
                 geom_col(stat="identity", position='dodge') +
                 scale_fill_manual(values=('#79706e', '#86bcb6')) +
                 theme(legend_title=element_blank()))
        if save_plot:
            _file3 = "{}/assign_isotype.pdf".format(os.path.dirname(filePath))
            save_as_pdf_pages([p], filename=_file3)
            print(p)
        else:
            print(p)


def assign_isotypes(fastas: Sequence,
                    fileformat: Literal['blast', 'changeo', 'airr'] = 'blast',
                    org: Literal['human', 'mouse'] = 'human',
                    correct_c_call: bool = True,
                    correction_dict: Union[None, Dict] = None,
                    plot: bool = True,
                    save_plot: bool = False,
                    figsize: Tuple[Union[int, float], Union[int,
                                                            float]] = (4, 4),
                    blastdb: Union[None, str] = None,
                    allele: bool = False,
                    parallel: bool = True,
                    ncpu: Union[None, int] = None,
                    filename_prefix: Union[None, Sequence, str] = None,
                    verbose: bool = False):
    """
    Annotate contigs with constant region call using blastn.

    Parameters
    ----------
    fastas : Sequence
        list or sequence of paths to fasta files.
    fileformat : str
        format of V(D)J file/objects. Default is 'blast'. Also accepts 'changeo' (same behaviour as 'blast') and 'airr'.
    org : str
        organism of reference folder. Default is 'human'.
    correct_c_call : bool
        whether or not to adjust the c_calls after blast based on provided primers specified in `primer_dict` option.
        Default is True.
    correction_dict : Dict, optional
        a nested dictionary contain isotype/c_genes as keys and primer sequences as records to use for correcting
        annotated c_calls. Defaults to a curated dictionary for human sequences if left as none.
    plot : bool
        whether or not to plot reassignment summary metrics. Default is True.
    save_plot : bool
        whether or not to save plots.
    figsize : Tuple[Union[int,float], Union[int,float]]
        size of figure. Default is (4, 4).
    blastdb : str, optional
        path to blast database. Defaults to `$BLASTDB` environmental variable.
    allele : bool
        whether or not to return allele calls. Default is False.
    parallel : bool
        whether or not to use parallelization. Default is True.
    ncpu : int
        number of cores to use if parallel is True. Default is all available - 1.
    filename_prefix : str, optional
        list of prefixes of file names preceding '_contig'. None defaults to 'filtered'.
    verbose : bool
        whether or not to print the blast command in terminal. Default is False.

    Returns
    -------
    V(D)J tsv files with constant genes annotated.
    """
    if type(fastas) is not list:
        fastas = [fastas]
    if type(filename_prefix) is not list:
        filename_prefix = [filename_prefix]
    if all(t is None for t in filename_prefix):
        filename_prefix = [None for f in fastas]

    if verbose:
        print('Assign isotypes \n')

    for i in range(0, len(fastas)):
        assign_isotype(fastas[i],
                       fileformat=fileformat,
                       org=org,
                       correct_c_call=correct_c_call,
                       correction_dict=correction_dict,
                       plot=plot,
                       save_plot=save_plot,
                       figsize=figsize,
                       blastdb=blastdb,
                       allele=allele,
                       parallel=parallel,
                       ncpu=ncpu,
                       filename_prefix=filename_prefix[i],
                       verbose=verbose)


def reannotate_genes(data: Sequence,
                     igblast_db: Union[None, str] = None,
                     germline: Union[None, str, PathLike] = None,
                     org: Literal['human', 'ig'] = 'human',
                     loci: Literal['ig', 'tr'] = 'ig',
                     extended: bool = True,
                     filename_prefix: Union[None, Sequence, str] = None,
                     verbose: bool = False):
    """
    Reannotate cellranger fasta files with igblastn and parses to airr/changeo data format.

    Parameters
    ----------
    data : Sequence
        list of fasta file locations, or folder name containing fasta files. if provided as a single string,
        it will first be converted to a list; this allows for the function to be run on single/multiple samples.
    igblast_db : str, PathLike, optional
        path to igblast database folder. Defaults to `$IGDATA` environmental variable.
    germline : str, PathLike, optional
        path to germline database folder. Defaults to `$GERMLINE` environmental variable.
    org : str
        organism of germline database. Default is 'human'.
    loci : str
        mode for igblastn. Default is 'ig' for BCRs. Also accepts 'tr' for TCRs.
    extended : bool
        whether or not to transfer additional 10X annotions to output file. Default is True.
    filename_prefix : str, optional
        list of prefixes of file names preceding '_contig'. None defaults to 'filtered'.
    verbose :
        whether or not to print the igblast command used in the terminal. Default is False.

    Returns
    -------
    V(D)J data file in airr/changeo data format.
    """
    if type(data) is not list:
        data = [data]
    if type(filename_prefix) is not list:
        filename_prefix = [filename_prefix]
    if all(t is None for t in filename_prefix):
        filename_prefix = [None for d in data]

    filePath = None
    for i in tqdm(range(0, len(data)), desc='Assigning genes '):
        filePath = check_filepath(data[i],
                                  filename_prefix=filename_prefix[i],
                                  endswith='.fasta')
        if filePath is None:
            if filename_prefix[i] is not None:
                raise OSError(
                    'Path to fasta file for {} with filename prefix `{}_contig` is unknown. '
                    .format(s, filename_prefix[i]) +
                    'Please specify path to fasta file or folder containing fasta file.'
                )
            else:
                raise OSError(
                    'Path to fasta file for {} is unknown. '.format(s) +
                    'Please specify path to fasta file or folder containing fasta file.'
                )

        if verbose:
            print('Processing {} \n'.format(filePath))

        assigngenes_igblast(filePath,
                            igblast_db=igblast_db,
                            org=org,
                            loci=loci,
                            verbose=verbose)
        makedb_igblast(filePath,
                       org=org,
                       germline=germline,
                       extended=extended,
                       verbose=verbose)
    if loci == 'tr':
        change_file_location(data, filename_prefix)


def reassign_alleles(data: Sequence,
                     combined_folder: Union[str, PathLike],
                     v_germline: Union[None, str] = None,
                     germline: Union[None, str, PathLike] = None,
                     org: Literal['human', 'mouse'] = 'human',
                     v_field: Literal['v_call',
                                      'v_call_genotyped'] = 'v_call_genotyped',
                     germ_types: Literal['full', 'dmask', 'vonly',
                                         'regions'] = 'dmask',
                     novel: bool = True,
                     cloned: bool = False,
                     plot: bool = True,
                     save_plot: bool = False,
                     figsize: Tuple[Union[int, float], Union[int,
                                                             float]] = (4, 3),
                     sample_id_dictionary: Union[None, Dict] = None,
                     filename_prefix: Union[None, Sequence, str] = None,
                     verbose: bool = False):
    """
    Correct allele calls based on a personalized genotype using tigger-reassignAlleles.

    It uses a subject-specific genotype to correct correct preliminary allele assignments of a set of
    sequences derived from a single subject.

    Parameters
    ----------
    data : Sequence
        list of data folders containing the .tsv files. if provided as a single string, it will first be converted to a
        list; this allows for the function to be run on single/multiple samples.
    combined_folder : str, PathLike
        name of folder for concatenated data file and genotyped files.
    v_germline : str, optional
        path to heavy chain v germline fasta. Defaults to IGHV fasta in `$GERMLINE` environmental variable.
    germline : str, optional
        path to germline database folder. Defaults to `$GERMLINE` environmental variable.
    org : str
        organism of germline database. Default is 'human'.
    v_field : str
        name of column containing the germline V segment call. Default is 'v_call_genotyped' (airr) for after tigger.
    germ_types : str
        Specify type of germline for reconstruction. Accepts one of : 'full', 'dmask', 'vonly', 'region'.
        Default is 'dmask'.
    novel : bool
        whether or not to run novel allele discovery during tigger-genotyping. Default is True (yes).
    cloned : bool
        whether or not to run CreateGermlines.py with `--cloned`.
    plot : bool
        whether or not to plot reassignment summary metrics. Default is True.
    save_plot : bool
        whether or not to save plot.
    figsize : Tuple[Union[int,float], Union[int,float]]
        size of figure. Default is (4, 3).
    sample_id_dictionary : dict, optional
        dictionary for creating a sample_id column in the concatenated file.
    filename_prefix : str, optional
        list of prefixes of file names preceding '_contig'. None defaults to 'filtered'.
    verbose : bool
        Whether or not to print the command used in the terminal. Default is False.

    Returns
    -------
    Individual V(D)J data files with v_call_genotyped column containing reassigned heavy chain v calls
    """
    fileformat = 'blast'
    if type(data) is not list:
        data = [data]
    if type(filename_prefix) is not list:
        filename_prefix = [filename_prefix]
    if all(t is None for t in filename_prefix):
        filename_prefix = [None for d in data]

    informat_dict = {
        'changeo': '_igblast_db-pass.tsv',
        'blast': '_igblast_db-pass.tsv',
        'airr': '_igblast_gap.tsv'
    }
    germpass_dict = {
        'changeo': '_igblast_db-pass_germ-pass.tsv',
        'blast': '_igblast_db-pass_germ-pass.tsv',
        'airr': '_igblast_gap_germ-pass.tsv'
    }
    heavy_dict = {
        'changeo': '_igblast_db-pass_heavy_parse-select.tsv',
        'blast': '_igblast_db-pass_heavy_parse-select.tsv',
        'airr': '_igblast_gap_heavy_parse-select.tsv'
    }
    light_dict = {
        'changeo': '_igblast_db-pass_light_parse-select.tsv',
        'blast': '_igblast_db-pass_light_parse-select.tsv',
        'airr': '_igblast_gap_light_parse-select.tsv'
    }
    fileformat_dict = {
        'changeo': '_igblast_db-pass_genotyped.tsv',
        'blast': '_igblast_db-pass_genotyped.tsv',
        'airr': '_igblast_gap_genotyped.tsv'
    }
    fileformat_passed_dict = {
        'changeo': '_igblast_db-pass_genotyped_germ-pass.tsv',
        'blast': '_igblast_db-pass_genotyped_germ-pass.tsv',
        'airr': '_igblast_gap_genotyped_germ-pass.tsv'
    }
    inferred_fileformat_dict = {
        'changeo': '_igblast_db-pass_inferredGenotype.txt',
        'blast': '_igblast_db-pass_inferredGenotype.txt',
        'airr': '_igblast_gap_inferredGenotype.txt'
    }
    germline_dict = {
        'changeo': '_igblast_db-pass_genotype.fasta',
        'blast': '_igblast_db-pass_genotype.fasta',
        'airr': '_igblast_gap_genotype.fasta'
    }
    fform_dict = {'blast': 'airr', 'airr': 'airr', 'changeo': 'changeo'}

    filepathlist_heavy = []
    filepathlist_light = []
    filePath = None
    sampleNames_dict = {}
    filePath_dict = {}
    for i in tqdm(range(0, len(data)), desc='Processing data file(s) '):
        filePath = check_filepath(data[i],
                                  filename_prefix=filename_prefix[i],
                                  endswith=informat_dict[fileformat],
                                  subdir='tmp')
        if filePath is None:
            raise OSError(
                'Path to .tsv file for {} is unknown. '.format(data[i]) +
                'Please specify path to reannotated .tsv file or folder containing reannotated .tsv file.'
            )

        filePath_heavy = filePath.replace(informat_dict[fileformat],
                                          heavy_dict[fileformat])
        filePath_light = filePath.replace(informat_dict[fileformat],
                                          light_dict[fileformat])

        if sample_id_dictionary is not None:
            sampleNames_dict[filePath] = sample_id_dictionary[data[i]]
        else:
            sampleNames_dict[filePath] = str(data[i])

        filePath_dict[str(data[i])] = filePath

        # splitting up to heavy chain and light chain files
        parsedb_heavy(filePath)
        parsedb_light(filePath)

        # add to counter
        filepathlist_heavy.append(filePath_heavy)
        filepathlist_light.append(filePath_light)

    # make output directory
    outDir = combined_folder.rstrip('/')
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    # concatenate
    if len(filepathlist_heavy) > 1:
        print('Concatenating objects')
        cmd1 = ' '.join(
            ['cat'] + [f for f in filepathlist_heavy] + ['>'] +
            [outDir + '/' + outDir + '_heavy' + informat_dict[fileformat]])
        cmd2 = ' '.join(
            ['cat'] + [f for f in filepathlist_light] + ['>'] +
            [outDir + '/' + outDir + '_light' + informat_dict[fileformat]])
    else:
        cmd1 = ' '.join(
            ['cat'] + [filepathlist_heavy[0]] + ['>'] +
            [outDir + '/' + outDir + '_heavy' + informat_dict[fileformat]])
        cmd2 = ' '.join(
            ['cat'] + [filepathlist_light[0]] + ['>'] +
            [outDir + '/' + outDir + '_light' + informat_dict[fileformat]])

    if verbose:
        print('Running command: %s\n' % (cmd1))
        print('Running command: %s\n' % (cmd2))
    os.system(cmd1)
    os.system(cmd2)

    novel_dict = {True: 'YES', False: 'NO'}
    if novel:
        try:
            print('      Running tigger-genotype with novel allele discovery.')
            tigger_genotype(outDir + '/' + outDir + '_heavy' +
                            informat_dict[fileformat],
                            v_germline=v_germline,
                            fileformat=fform_dict[fileformat],
                            novel_=novel_dict[novel],
                            verbose=verbose)
            creategermlines(outDir + '/' + outDir + '_heavy' +
                            fileformat_dict[fileformat],
                            germtypes=germ_types,
                            mode='heavy',
                            genotype_fasta=outDir + '/' + outDir + '_heavy' +
                            germline_dict[fileformat],
                            germline=germline,
                            v_field=v_field,
                            verbose=verbose,
                            cloned=cloned)
            _ = load_data(outDir + '/' + outDir + '_heavy' +
                          fileformat_passed_dict[fileformat])
        except:
            try:
                print('      Novel allele discovery execution halted.')
                print(
                    '      Attempting to run tigger-genotype without novel allele discovery.'
                )
                tigger_genotype(outDir + '/' + outDir + '_heavy' +
                                informat_dict[fileformat],
                                v_germline=v_germline,
                                fileformat=fform_dict[fileformat],
                                novel_=novel_dict[False],
                                verbose=verbose)
                creategermlines(outDir + '/' + outDir + '_heavy' +
                                fileformat_dict[fileformat],
                                germtypes=germ_types,
                                mode='heavy',
                                genotype_fasta=outDir + '/' + outDir +
                                '_heavy' + germline_dict[fileformat],
                                germline=germline,
                                v_field=v_field,
                                verbose=verbose,
                                cloned=cloned)
                _ = load_data(outDir + '/' + outDir + '_heavy' +
                              fileformat_passed_dict[fileformat])
            except:
                print(
                    '     Insufficient contigs for running tigger-genotype. Defaulting to original heavy chain v_calls.'
                )
                tigger_failed = ''
    else:
        try:
            print(
                '      Running tigger-genotype without novel allele discovery.'
            )
            tigger_genotype(outDir + '/' + outDir + '_heavy' +
                            informat_dict[fileformat],
                            v_germline=v_germline,
                            fileformat=fform_dict[fileformat],
                            novel_=novel_dict[False],
                            verbose=verbose)
            creategermlines(outDir + '/' + outDir + '_heavy' +
                            fileformat_dict[fileformat],
                            germtypes=germ_types,
                            mode='heavy',
                            genotype_fasta=outDir + '/' + outDir + '_heavy' +
                            germline_dict[fileformat],
                            germline=germline,
                            v_field=v_field,
                            verbose=verbose,
                            cloned=cloned)
            _ = load_data(outDir + '/' + outDir + '_heavy' +
                          fileformat_passed_dict[fileformat])
        except:
            print(
                '      Insufficient contigs for running tigger-genotype. Defaulting to original heavy chain v_calls.'
            )
            tigger_failed = ''

    if 'tigger_failed' in locals():
        creategermlines(outDir + '/' + outDir + '_heavy' +
                        informat_dict[fileformat],
                        germtypes=germ_types,
                        mode='heavy',
                        genotype_fasta=None,
                        germline=germline,
                        v_field='v_call',
                        verbose=verbose,
                        cloned=cloned)
        creategermlines(outDir + '/' + outDir + '_light' +
                        informat_dict[fileformat],
                        germtypes=germ_types,
                        mode='light',
                        genotype_fasta=None,
                        germline=germline,
                        v_field='v_call',
                        verbose=verbose,
                        cloned=cloned)
        print(
            '      For convenience, entries for heavy chain in `v_call` are copied to `v_call_genotyped`.'
        )
        heavy = load_data(outDir + '/' + outDir + '_heavy' +
                          germpass_dict[fileformat])
        heavy['v_call_genotyped'] = heavy['v_call']
        print(
            '      For convenience, entries for light chain `v_call` are copied to `v_call_genotyped`.'
        )
        light = load_data(outDir + '/' + outDir + '_light' +
                          germpass_dict[fileformat])
        light['v_call_genotyped'] = light['v_call']
    else:
        creategermlines(outDir + '/' + outDir + '_light' +
                        informat_dict[fileformat],
                        germtypes=germ_types,
                        mode='light',
                        genotype_fasta=None,
                        germline=germline,
                        v_field='v_call',
                        verbose=verbose,
                        cloned=cloned)
        heavy = load_data(outDir + '/' + outDir + '_heavy' +
                          fileformat_passed_dict[fileformat])
        print(
            '      For convenience, entries for light chain `v_call` are copied to `v_call_genotyped`.'
        )
        light = load_data(outDir + '/' + outDir + '_light' +
                          germpass_dict[fileformat])
        light['v_call_genotyped'] = light['v_call']

    sampledict = {}
    heavy['sample_id'], light['sample_id'] = None, None
    for file in sampleNames_dict.keys():
        dat_f = load_data(file)
        dat_f['sample_id'] = sampleNames_dict[file]
        heavy['sample_id'].update(dat_f['sample_id'])
        light['sample_id'].update(dat_f['sample_id'])

    dat_ = heavy.append(light)
    if 'cell_id' in dat_.columns:
        dat_.sort_values(by='cell_id', inplace=True)
    else:
        dat_.sort_values(by='sequence_id', inplace=True)

    if plot:
        if 'tigger_failed' not in locals():
            print('Returning summary plot')
            inferred_genotype = outDir + '/' + outDir + \
                '_heavy' + inferred_fileformat_dict[fileformat]
            inf_geno = pd.read_csv(inferred_genotype, sep='\t', dtype='object')

            s2 = set(inf_geno['gene'])
            results = []
            try:
                for samp in list(set(heavy['sample_id'])):
                    res_x = heavy[(heavy['sample_id'] == samp)]
                    V_ = [
                        re.sub('[*][0-9][0-9]', '', v) for v in res_x['v_call']
                    ]
                    V_g = [
                        re.sub('[*][0-9][0-9]', '', v)
                        for v in res_x['v_call_genotyped']
                    ]
                    s1 = set(
                        list(','.join([
                            ','.join(list(set(v.split(',')))) for v in V_
                        ]).split(',')))
                    setdiff = s1 - s2
                    ambiguous = (["," in i
                                  for i in V_].count(True) / len(V_) * 100,
                                 ["," in i
                                  for i in V_g].count(True) / len(V_g) * 100)
                    not_in_genotype = ([i in setdiff
                                        for i in V_].count(True) / len(V_) *
                                       100, [i in setdiff
                                             for i in V_g].count(True) /
                                       len(V_g) * 100)
                    stats = pd.DataFrame(
                        [ambiguous, not_in_genotype],
                        columns=['ambiguous', 'not_in_genotype'],
                        index=['before', 'after']).T
                    stats.index.set_names(['vgroup'], inplace=True)
                    stats.reset_index(drop=False, inplace=True)
                    stats['sample_id'] = samp
                    # stats['donor'] = str(combined_folder)
                    results.append(stats)
                results = pd.concat(results)
                ambiguous_table = results[results['vgroup'] == 'ambiguous']
                not_in_genotype_table = results[results['vgroup'] ==
                                                'not_in_genotype']
                ambiguous_table.reset_index(inplace=True, drop=True)
                not_in_genotype_table.reset_index(inplace=True, drop=True)
                # melting the dataframe
                ambiguous_table_before = ambiguous_table.drop('after', axis=1)
                ambiguous_table_before.rename(columns={"before": "var"},
                                              inplace=True)
                ambiguous_table_before['var_group'] = 'before'
                ambiguous_table_after = ambiguous_table.drop('before', axis=1)
                ambiguous_table_after.rename(columns={"after": "var"},
                                             inplace=True)
                ambiguous_table_after['var_group'] = 'after'
                ambiguous_table = pd.concat(
                    [ambiguous_table_before, ambiguous_table_after])
                not_in_genotype_table_before = not_in_genotype_table.drop(
                    'after', axis=1)
                not_in_genotype_table_before.rename(columns={"before": "var"},
                                                    inplace=True)
                not_in_genotype_table_before['var_group'] = 'before'
                not_in_genotype_table_after = not_in_genotype_table.drop(
                    'before', axis=1)
                not_in_genotype_table_after.rename(columns={"after": "var"},
                                                   inplace=True)
                not_in_genotype_table_after['var_group'] = 'after'
                not_in_genotype_table = pd.concat([
                    not_in_genotype_table_before, not_in_genotype_table_after
                ])
                ambiguous_table['var_group'] = ambiguous_table[
                    'var_group'].astype('category')
                not_in_genotype_table['var_group'] = not_in_genotype_table[
                    'var_group'].astype('category')
                ambiguous_table['var_group'].cat.reorder_categories(
                    ['before', 'after'], inplace=True)
                not_in_genotype_table['var_group'].cat.reorder_categories(
                    ['before', 'after'], inplace=True)

                options.figure_size = figsize
                final_table = pd.concat(
                    [ambiguous_table, not_in_genotype_table])
                p = (ggplot(final_table,
                            aes(x='sample_id', y='var', fill='var_group')) +
                     coord_flip() + theme_classic() + xlab("sample_id") +
                     ylab("% allele calls") +
                     ggtitle("Genotype reassignment with TIgGER") +
                     geom_bar(stat="identity") +
                     facet_grid('~' + str('vgroup'), scales="free_y") +
                     scale_fill_manual(values=('#86bcb6', '#F28e2b')) +
                     theme(legend_title=element_blank()))
                if save_plot:
                    savefile = outDir + '/' + outDir + '_reassign_alleles.pdf'
                    save_as_pdf_pages([p], filename=savefile)
                    print(p)
                else:
                    print(p)
            except:
                print('Error in plotting encountered. Skipping.')
                pass
        else:
            pass
    sleep(0.5)
    # if split_write_out:
    if 'tigger_failed' in locals():
        print(
            'Although tigger-genotype was not run successfully, file will still be saved with `_genotyped.tsv`'
            'extension for convenience.')
    for s in tqdm(data, desc='Writing out to individual folders '):
        if sample_id_dictionary is not None:
            out_file = dat_[dat_['sample_id'] == sample_id_dictionary[s]]
        else:
            out_file = dat_[dat_['sample_id'] == s]
        outfilepath = filePath_dict[s]
        out_file.to_csv(outfilepath.replace('.tsv', '_genotyped.tsv'),
                        index=False,
                        sep='\t')


def create_germlines(
        self: Union[Dandelion, pd.DataFrame, str, PathLike],
        germline: Union[None, str, PathLike] = None,
        org: Literal['human', 'mouse'] = 'human',
        seq_field: Literal['sequence_alignment'] = 'sequence_alignment',
        v_field: Literal['v_call', 'v_call_genotyped'] = 'v_call',
        d_field: Literal['d_call'] = 'd_call',
        j_field: Literal['j_call'] = 'j_call',
        germ_types: Literal['full', 'dmask', 'vonly', 'regions'] = 'dmask',
        fileformat: Literal['changeo', 'airr'] = 'airr',
        initialize_metadata: bool = False) -> Dandelion:
    """
    Run CreateGermlines.py to reconstruct the germline V(D)J sequence.

    Run CreateGermlines.py to reconstruct the germline V(D)J sequence.

    Parameters
    ----------
    self : Dandelion, pd.DataFrame, str
        `Dandelion` object, pandas `DataFrame` in changeo/airr format, or file path to changeo/airr
        file after clones have been determined.
    germline : str, optional
        path to germline database folder. Defaults to `$GERMLINE` environmental variable.
    org : str
        organism of germline database. Default is 'human'.
    seq_field : str
        name of column containing the aligned sequence. Default is 'sequence_alignment' (airr).
    v_field : str
        name of column containing the germline V segment call. Default is 'v_call' (airr).
    d_field : str
        name of column containing the germline d segment call. Default is 'd_call' (airr).
    j_field : str
        name of column containing the germline j segment call. Default is 'j_call' (airr).
    germ_types : str
        Specify type(s) of germlines to include full germline, germline with D segment masked,
        or germline for V segment only. Default is 'dmask'.
    fileformat : str
        format of V(D)J file/objects. Default is 'airr'. Also accepts 'changeo'.

    Returns
    -------
    V(D)J data file with reconstructed germline sequences.
    """
    start = logg.info('Reconstructing germline sequences')
    env = os.environ.copy()
    if germline is None:
        try:
            gml = env['GERMLINE']
        except:
            raise OSError(
                'Environmental variable GERMLINE must be set. ' +
                'Otherwise, please provide path to folder containing germline fasta files.'
            )
        gml = gml + 'imgt/' + org + '/vdj/'
    else:
        if os.path.isdir(germline):
            env['GERMLINE'] = germline
            gml = germline

    def _parseChangeO(record):
        """
        Parse a dictionary to a Receptor object.

        Arguments:
          record : dict with fields and values in the Change-O format

        Returns:
        changeo.Receptor.Receptor : parsed Receptor object.
        """
        # Parse fields
        result = {}
        for k, v in record.items():
            k = ChangeoSchema.toReceptor(k)
            result[k] = v

        return Receptor(result)

    def _parseAIRR(record):
        """
        Parse a dictionary of AIRR records to a Receptor object.

        Arguments:
          record : dict with fields and values in the AIRR format.

        Returns:
        changeo.Receptor.Receptor : parsed Receptor object.
        """
        # Parse fields
        result = {}
        for k, v in record.items():
            # Rename fields
            k = AIRRSchema.toReceptor(k)
            # Convert start positions to 0-based
            # if k in ReceptorData.start_fields and v is not None and v != '':
            #     v = str(int(v) + 1)
            # Assign new field
            result[k] = v

        for end, (start, length) in ReceptorData.end_fields.items():
            if end in result and result[end] is not None:
                try:
                    result[length] = int(result[end]) - int(result[start]) + 1
                except:
                    pass

        return Receptor(result)

    def _create_germlines_object(self, references, seq_field, v_field, d_field,
                                 j_field, germ_types, fileformat):
        """
        Write germline sequences to tab-delimited database file.

        Arguments:
        self : dandelion_class object
        references : folders and/or files containing germline repertoire data in FASTA format.
        seq_field : field in which to look for sequence.
        v_field : field in which to look for V call.
        d_field : field in which to look for D call.
        j_field : field in which to look for J call.
        # cloned : if True build germlines by clone, otherwise build individual germlines.
        # clone_field : field containing clone identifiers; ignored if cloned=False.
        germ_types : list of germline sequence types to be output from the set of 'full', 'dmask', 'vonly', 'regions'
        fileformat : str
            format of V(D)J file/objects. Default is 'airr'. Also accepts 'changeo'.

        Returns:
        """
        # Define format operators
        try:
            reader, writer, schema = getFormatOperators(fileformat)
        except:
            raise ValueError('Invalid format %s' % fileformat)

        # Define output germline fields
        germline_fields = OrderedDict()
        seq_type = seq_field.split('_')[-1]
        if 'full' in germ_types:
            germline_fields['full'] = 'germline_' + seq_type
        if 'dmask' in germ_types:
            germline_fields['dmask'] = 'germline_' + seq_type + '_d_mask'
        if 'vonly' in germ_types:
            germline_fields['vonly'] = 'germline_' + seq_type + '_v_region'
        if 'regions' in germ_types:
            germline_fields['regions'] = 'germline_regions'

        if type(references) is dict:
            reference_dict = references
        else:
            if type(references) is not list:
                ref = [references]
            else:
                ref = references
            reference_dict = readGermlines(ref)
        # Check for IMGT-gaps in germlines
        if all('...' not in x for x in reference_dict.values()):
            warnings.warn(
                UserWarning(
                    'Germline reference sequences do not appear to contain IMGT-numbering spacers.'
                    + ' Results may be incorrect.'))

        required = [
            'v_germ_start_imgt', 'd_germ_start', 'j_germ_start', 'np1_length',
            'np2_length'
        ]

        if self.__class__ == Dandelion:
            if isinstance(self.data, pd.DataFrame):
                # Check for required columns
                try:
                    checkFields(required, self.data.columns, schema=schema)
                except LookupError as e:
                    print(e)

                # Count input
                # total_count = len(self.data)

                # Check for existence of fields
                for f in [v_field, d_field, j_field, seq_field]:
                    if f not in self.data.columns:
                        raise NameError(
                            '%s field does not exist in input database file.' %
                            f)
                # Translate to Receptor attribute names
                v_field_ = schema.toReceptor(v_field)
                d_field_ = schema.toReceptor(d_field)
                j_field_ = schema.toReceptor(j_field)
                seq_field_ = schema.toReceptor(seq_field)
                # clone_field = schema.toReceptor(clone_field)

                # Define Receptor iterator
                receptor_iter = ((self.data.loc[x, ].sequence_id,
                                  self.data.loc[x, ]) for x in self.data.index)

            else:
                raise LookupError(
                    'Please initialise the Dandelion object with a dataframe in data slot.'
                )
        elif self.__class__ == pd.DataFrame:
            try:
                checkFields(required, self.columns, schema=schema)
            except LookupError as e:
                print(e)

            # Count input
            # total_count = len(self)
            # Check for existence of fields
            for f in [v_field, d_field, j_field, seq_field]:
                if f not in self.columns:
                    raise NameError(
                        '%s field does not exist in input database file.' % f)
            # Translate to Receptor attribute names
            v_field_ = schema.toReceptor(v_field)
            d_field_ = schema.toReceptor(d_field)
            j_field_ = schema.toReceptor(j_field)
            seq_field_ = schema.toReceptor(seq_field)
            # clone_field = schema.toReceptor(clone_field)
            # Define Receptor iterator
            receptor_iter = ((self.loc[x, ].sequence_id, self.loc[x, ])
                             for x in self.index)

        out = {}
        # Iterate over rows
        for key, records in tqdm(
                receptor_iter,
                desc="   Building {} germline sequences".format(germ_types)):
            # Define iteration variables
            # Build germline for records
            if fileformat == 'airr':
                germ_log, glines, genes = buildGermline(_parseAIRR(
                    dict(records)),
                    reference_dict,
                    seq_field=seq_field_,
                    v_field=v_field_,
                    d_field=d_field_,
                    j_field=j_field_)
            elif fileformat == 'changeo':
                germ_log, glines, genes = buildGermline(_parseChangeO(
                    dict(records)),
                    reference_dict,
                    seq_field=seq_field_,
                    v_field=v_field_,
                    d_field=d_field_,
                    j_field=j_field_)
            else:
                raise AttributeError('%s is not acceptable file format.' %
                                     fileformat)

            if glines is not None:
                # Add glines to Receptor record
                annotations = {}
                if 'full' in germ_types:
                    annotations[germline_fields['full']] = glines['full']
                if 'dmask' in germ_types:
                    annotations[germline_fields['dmask']] = glines['dmask']
                if 'vonly' in germ_types:
                    annotations[germline_fields['vonly']] = glines['vonly']
                if 'regions' in germ_types:
                    annotations[germline_fields['regions']] = glines['regions']
                out.update({key: annotations})
        germline_df = pd.DataFrame.from_dict(out, orient='index')

        if self.__class__ == Dandelion:
            # datx = load_data(self.data)
            for x in germline_df.columns:
                self.data[x] = pd.Series(germline_df[x])

        elif self.__class__ == pd.DataFrame:
            datx = load_data(self)
            for x in germline_df.columns:
                datx[x] = pd.Series(germline_df[x])
            try:
                output = Dandelion(data=datx,
                                   germline=reference_dict,
                                   initialize=True)
            except:
                output = Dandelion(data=datx,
                                   germline=reference_dict,
                                   initialize=False)
            return (output)
        sleep(0.5)
        logg.info(
            ' finished',
            time=start,
            deep=
            ('Updated Dandelion object: \n'
             '   \'data\', updated germline alignment in contig-indexed clone table\n'
             '   \'germline\', updated germline reference\n'))

    def _create_germlines_file(file, references, seq_field, v_field, d_field,
                               j_field, germ_types, fileformat):
        """
        Write germline sequences to tab-delimited database file.

        Arguments:
        file : airr/changeo tsv file
        references : folders and/or files containing germline repertoire data in FASTA format.
        seq_field : field in which to look for sequence.
        v_field : field in which to look for V call.
        d_field : field in which to look for D call.
        j_field : field in which to look for J call.
        cloned : if True build germlines by clone, otherwise build individual germlines.
        germ_types : list of germline sequence types to be output from the set of 'full', 'dmask', 'vonly', 'regions'
        fileformat : str
                format of V(D)J file/objects. Default is 'airr'. Also accepts 'changeo'.
        Returns:
        """
        # Define format operators
        try:
            reader, writer, schema = getFormatOperators(fileformat)
        except:
            raise ValueError('Invalid format %s' % fileformat)

        # Define output germline fields
        germline_fields = OrderedDict()
        seq_type = seq_field.split('_')[-1]
        if 'full' in germ_types:
            germline_fields['full'] = 'germline_' + seq_type
        if 'dmask' in germ_types:
            germline_fields['dmask'] = 'germline_' + seq_type + '_d_mask'
        if 'vonly' in germ_types:
            germline_fields['vonly'] = 'germline_' + seq_type + '_v_region'
        if 'regions' in germ_types:
            germline_fields['regions'] = 'germline_regions'

        if type(references) is dict:
            reference_dict = references
        else:
            if type(references) is not list:
                ref = [references]
            else:
                ref = references
            reference_dict = readGermlines(ref)
        # Check for IMGT-gaps in germlines
        if all('...' not in x for x in reference_dict.values()):
            warnings.warn(
                UserWarning(
                    'Germline reference sequences do not appear to contain IMGT-numbering spacers. '
                    + 'Results may be incorrect.'))

        required = [
            'v_germ_start_imgt', 'd_germ_start', 'j_germ_start', 'np1_length',
            'np2_length'
        ]

        # Get repertoire and open Db reader
        db_handle = open(file, 'rt')
        db_iter = reader(db_handle)
        # Check for required columns
        try:
            checkFields(required, db_iter.fields, schema=schema)
        except LookupError as e:
            print(e)
        # Count input
        # total_count = countDbFile(file)
        # Check for existence of fields
        for f in [v_field, d_field, j_field, seq_field]:
            if f not in db_iter.fields:
                raise NameError(
                    '%s field does not exist in input database file.' % f)
        # Translate to Receptor attribute names
        v_field_ = schema.toReceptor(v_field)
        d_field_ = schema.toReceptor(d_field)
        j_field_ = schema.toReceptor(j_field)
        seq_field_ = schema.toReceptor(seq_field)
        # clone_field = schema.toReceptor(clone_field)
        # Define Receptor iterator
        receptor_iter = ((x.sequence_id, [x]) for x in db_iter)

        out = {}
        # Iterate over rows
        for key, records in tqdm(
                receptor_iter,
                desc="   Building {} germline sequences".format(germ_types)):
            # Define iteration variables
            # Build germline for records
            # if not isinstance(self.data, pd.DataFrame):
            records = list(records)
            germ_log, glines, genes = buildGermline(records[0],
                                                    reference_dict,
                                                    seq_field=seq_field_,
                                                    v_field=v_field_,
                                                    d_field=d_field_,
                                                    j_field=j_field_)
            if glines is not None:
                # Add glines to Receptor record
                annotations = {}
                if 'full' in germ_types:
                    annotations[germline_fields['full']] = glines['full']
                if 'dmask' in germ_types:
                    annotations[germline_fields['dmask']] = glines['dmask']
                if 'vonly' in germ_types:
                    annotations[germline_fields['vonly']] = glines['vonly']
                if 'regions' in germ_types:
                    annotations[germline_fields['regions']] = glines['regions']
                out.update({key: annotations})
        germline_df = pd.DataFrame.from_dict(out, orient='index')

        try:
            out = Dandelion(data=file,
                            germline=reference_dict,
                            initialize=True)
        except:
            out = Dandelion(data=file,
                            germline=reference_dict,
                            initialize=False)
        for x in germline_df.columns:
            out.data[x] = pd.Series(germline_df[x])

        if os.path.isfile(str(file)):
            out.data.to_csv("{}/{}_germline_{}.tsv".format(
                os.path.dirname(file),
                os.path.basename(file).split('.tsv')[0], germ_types),
                sep='\t',
                index=False)
        return (out)

    if (type(germline) is dict) or (type(germline) is list):
        if self.__class__ == Dandelion:
            _create_germlines_object(self, germline, seq_field, v_field,
                                     d_field, j_field, germ_types, fileformat)
        elif self.__class__ == pd.DataFrame:
            return (_create_germlines_object(self, germline, seq_field,
                                             v_field, d_field, j_field,
                                             germ_types, fileformat))
        else:
            return (_create_germlines_file(self, germline, seq_field, v_field,
                                           d_field, j_field, germ_types,
                                           fileformat))
    else:
        if self.__class__ == Dandelion:
            if len(self.germline) != 0:
                _create_germlines_object(self, self.germline, seq_field,
                                         v_field, d_field, j_field, germ_types,
                                         fileformat)
            else:
                _create_germlines_object(self, gml, seq_field, v_field,
                                         d_field, j_field, germ_types,
                                         fileformat)
        elif self.__class__ == pd.DataFrame:
            return (_create_germlines_object(self, gml, seq_field, v_field,
                                             d_field, j_field, germ_types,
                                             fileformat))
        else:
            return (_create_germlines_file(self, gml, seq_field, v_field,
                                           d_field, j_field, germ_types,
                                           fileformat))


def filter_contigs(data: Union[Dandelion, pd.DataFrame, str],
                   adata: AnnData,
                   filter_contig: bool = True,
                   filter_rna: bool = False,
                   filter_poorqualitycontig: bool = False,
                   rescue_vdj: bool = True,
                   umi_foldchange_cutoff: int = 2,
                   filter_vj_chains: bool = True,
                   filter_missing: bool = True,
                   productive_only: bool = True,
                   locus: Union[None, Literal['ig', 'tr-ab', 'tr-gd']] = None,
                   save: Union[None, str] = None,
                   **kwargs) -> Tuple[Dandelion, AnnData]:
    """
    Filter doublets and poor quality cells and corresponding contigs based on provided V(D)J `DataFrame` and `AnnData`.

    Depends on a `AnnData`.obs slot populated with 'filter_rna' column. If the aligned sequence is an exact match
    between contigs, the contigs will be merged into the one with the highest umi count, adding the summing the
    umi count of the duplicated contigs to duplicate_count column. After this check, if there are still multiple
    contigs, cells with multiple IGH contigs are filtered unless `rescue_vdj` is True, where by the umi counts for
    each IGH contig will then be compared. The contig with the highest umi that is > umi_foldchange_cutoff (default
    is empirically set at 5) from the lowest will be retained. If there are multiple contigs that survive the
    'rescue', then all contigs will be filtered. The default behaviour is to also filter cells with multiple
    lightchains but this may sometimes be a true biological occurrence; toggling filter_vj_chains to False will
    rescue the mutltiplet light chains. Lastly, contigs with no corresponding cell barcode in the AnnData object is
    filtered if filter_missing is True. However, this may be useful to toggle to False if more contigs are
    preferred to be kept or for integrating with bulk reperotire seq data.

    Parameters
    ----------
    data : Dandeion, pd.DataDrame, str
        V(D)J airr/changeo data to filter. Can be pandas `DataFrame` object or file path as string.
    adata : AnnData
        AnnData object to filter.
    filter_contig : bool
        If True, V(D)J `DataFrame` object returned will be filtered. Default is True.
    filter_rna : bool
        If True, `AnnData` object returned will be filtered based on potential V(D)J doublets. Default is False.
    filter_poorqualitycontig : bool
        If True, barcodes marked with poor quality contigs will be filtered. Default is False; only relevant contigs are
        removed and RNA barcodes are kept.
    rescue_vdj : bool
        If True, rescues IGH contigs with highest umi counts with a requirement that it passes the
        `umi_foldchange_cutoff` option. In addition, the sum of the all the heavy chain contigs must be greater than 3
        umi or all contigs will be filtered. Default is True.
    umi_foldchange_cutoff : int
        related to minimum fold change required to rescue heavy chain contigs/barcode otherwise they will be marked as
        doublets. Default is empirically set at 2-fold.
    filter_vj_chains : bool
        cells with multiple light chains will be marked to filter. Default is True.
    productive_only : bool
        whether or not to retain only productive contigs.
    filter_missing : bool
        cells in V(D)J data not found in `AnnData` object will be marked to filter. Default is True. This may be useful
        for toggling to False if integrating with bulk data.
    locus : str
        Mode for filtering data. Accepts one of 'ig', 'tr-ab' or 'tr-gd'. None defaults to 'ig'.
    save : str, optional
        Only used if a pandas dataframe or dandelion object is provided. Specifying will save the formatted vdj table.
    **kwargs
        additional kwargs passed to `Dandelion.Dandelion`.

    Returns
    -------
    V(D)J `DataFrame` object in airr/changeo format and `AnnData` object.
    """
    start = logg.info('Filtering BCRs')
    if data.__class__ == Dandelion:
        dat_ = load_data(data.data)
    else:
        dat_ = load_data(data)

    if productive_only:
        dat = dat_[dat_['productive'].isin(['T', 'True', 'TRUE', True])].copy()
    else:
        dat = dat_.copy()

    adata_ = adata.copy()

    if 'cell_id' not in dat.columns:
        raise AttributeError(
            "VDJ data does not contain 'cell_id' column. Please make sure this is populated before filtering."
        )
    if 'filter_rna' not in adata_.obs:
        adata_.obs['filter_rna'] = False

    barcode = list(set(dat['cell_id']))

    contig_check = pd.DataFrame(index=adata_.obs_names)
    bc_ = {}
    for b in barcode:
        bc_.update({b: True})
    contig_check['has_contig'] = pd.Series(bc_)
    contig_check.replace(np.nan, 'No_contig', inplace=True)
    adata_.obs['has_contig'] = pd.Series(contig_check['has_contig'])
    adata_.obs['has_contig'] = adata_.obs['has_contig'].astype('category')

    if 'v_call_genotyped' in dat.columns:
        v_dict = dict(zip(dat['sequence_id'], dat['v_call_genotyped']))
    else:
        v_dict = dict(zip(dat['sequence_id'], dat['v_call']))
    j_dict = dict(zip(dat['sequence_id'], dat['j_call']))
    c_dict = dict(zip(dat['sequence_id'], dat['c_call']))

    # rather than leaving a nan cell, i will create a 0 column for now
    if 'duplicate_count' in dat and 'umi_count' not in dat:
        dat['umi_count'] = dat['duplicate_count']  # just do a simple swap?
    elif 'duplicate_count' not in dat and 'umi_count' in dat:
        dat['duplicate_count'] = dat['umi_count']
    elif 'duplicate_count' in dat and 'umi_count' in dat:
        dat['umi_count'] = dat['duplicate_count']

    tofilter = FilterContigs(dat)

    for b in tqdm(barcode, desc='Scanning for poor quality/ambiguous contigs'):
        tofilter.run_scan(b, rescue_vdj, umi_foldchange_cutoff,
                          filter_poorqualitycontig, v_dict, j_dict, c_dict)

    poor_qual = tofilter.poor_qual
    h_doublet = tofilter.h_doublet
    l_doublet = tofilter.l_doublet
    drop_contig = tofilter.drop_contig
    umi_adjustment = tofilter.umi_adjustment

    if len(umi_adjustment) > 0:
        dat['duplicate_count'].update(umi_adjustment)

    poorqual = Tree()
    hdoublet = Tree()
    ldoublet = Tree()
    for c in tqdm(adata_.obs_names, desc='Annotating in anndata obs slot '):
        if c in poor_qual:
            poorqual[c] = True
        else:
            poorqual[c] = False

        if c in h_doublet:
            hdoublet[c] = True
        else:
            hdoublet[c] = False

        if c in l_doublet:
            ldoublet[c] = True
        else:
            ldoublet[c] = False

    adata_.obs['filter_contig_quality'] = pd.Series(dict(poorqual))
    adata_.obs['filter_contig_quality'] = adata_.obs[
        'filter_contig_quality'].astype('category')
    adata_.obs['filter_contig_VDJ'] = pd.Series(dict(hdoublet))
    adata_.obs['filter_contig_VDJ'] = adata_.obs['filter_contig_VDJ'].astype(
        'category')
    adata_.obs['filter_contig_VJ'] = pd.Series(dict(ldoublet))
    adata_.obs['filter_contig_VJ'] = adata_.obs['filter_contig_VJ'].astype(
        'category')

    drop_contig = list(set(flatten(drop_contig)))

    filter_ids = []
    if filter_contig:
        print('Finishing up filtering')
        if not filter_vj_chains:
            if filter_poorqualitycontig:
                filter_ids = list(set(h_doublet + poor_qual))
            else:
                filter_ids = list(set(h_doublet))
        else:
            if filter_poorqualitycontig:
                filter_ids = list(set(h_doublet + l_doublet + poor_qual))
            else:
                filter_ids = list(set(h_doublet + l_doublet))

        filter_ids = filter_ids + \
            list(adata_[adata_.obs['filter_rna'].isin([True, 'True'])].obs_names)
        filter_ids = list(set(filter_ids))

        if filter_missing:
            for c in dat['cell_id']:
                if c not in adata_.obs_names:
                    filter_ids.append(c)

        _dat = dat[~(dat['cell_id'].isin(filter_ids))].copy()
        _dat = _dat[~(_dat['sequence_id'].isin(drop_contig))].copy()

        # final check
        barcodes_final = list(set(_dat['cell_id']))
        filter_ids2 = []
        for b in barcodes_final:
            check_dat = _dat[(_dat['locus'].isin(['IGH', 'TRB', 'TRD']))
                             & (_dat['cell_id'].isin([b]))].copy()
            if check_dat.shape[0] < 1:
                filter_ids2.append(b)
        _dat = _dat[~(_dat['cell_id'].isin(filter_ids2))].copy()

        if _dat.shape[0] == 0:
            raise IndexError(
                'No BCRs passed filtering. Are you sure that the cell barcodes are matching?'
            )

        if os.path.isfile(str(data)):
            _dat.to_csv("{}/{}_filtered.tsv".format(
                os.path.dirname(data),
                os.path.basename(data).split('.tsv')[0]),
                sep='\t',
                index=False)
        else:
            if save is not None:
                if save.endswith('.tsv'):
                    _dat.to_csv(str(save), sep='\t', index=False)
                else:
                    raise OSError(
                        'Please provide a file name that ends with .tsv')
    else:
        _dat = dat.copy()

    if filter_contig:
        barcode1 = list(set(dat['cell_id']))

    barcode2 = list(set(_dat['cell_id']))

    if filter_contig:
        failed = list(set(barcode1) ^ set(barcode2))

    bc_2 = {}
    for b in barcode2:
        bc_2.update({b: True})
    if filter_contig:
        for b in failed:
            bc_2.update({b: False})
    contig_check['contig_QC_pass'] = pd.Series(bc_2)
    contig_check.replace(np.nan, 'No_contig', inplace=True)
    adata_.obs['contig_QC_pass'] = pd.Series(contig_check['contig_QC_pass'])
    adata_.obs['contig_QC_pass'] = adata_.obs['contig_QC_pass'].astype(
        'category')

    print('Initializing Dandelion object')
    out_dat = Dandelion(data=_dat, locus=locus, **kwargs)
    if data.__class__ == Dandelion:
        out_dat.germline = data.germline

    adata_.obs['filter_contig'] = adata_.obs_names.isin(filter_ids)
    adata_.obs['filter_contig'] = adata_.obs['filter_contig'].astype(
        'category')

    if filter_rna:
        # not saving the scanpy object because there's no need to at the moment
        out_adata = adata_[adata_.obs['filter_contig'].isin([False,
                                                             'False'])].copy()
    else:
        out_adata = adata_.copy()

    logg.info(' finished',
              time=start,
              deep=('Returning Dandelion and AnnData objects: \n'))
    return (out_dat, out_adata)


def quantify_mutations(self: Union[Dandelion, str, PathLike],
                       split_locus: bool = False,
                       sequence_column: Union[None, str] = None,
                       germline_column: Union[None, str] = None,
                       region_definition: Union[None, str] = None,
                       mutation_definition: Union[None, str] = None,
                       frequency: bool = False,
                       combine: bool = True) -> Union[pd.DataFrame, Dandelion]:
    """
    Run basic mutation load analysis.

    Implemented in `shazam <https://shazam.readthedocs.io/en/stable/vignettes/Mutation-Vignette/>`__.

    Parameters
    ----------
    self : Dandelion, str, PathLike
        `Dandelion` object, file path to AIRR file.
    split_locus : bool
        whether to return the results for heavy chain and light chain separately. Default is False.
    sequence_column: str, optional
        passed to shazam's `observedMutations`. https://shazam.readthedocs.io/en/stable/topics/observedMutations
    germline_column: str, optional
        passed to shazam's `observedMutations`. https://shazam.readthedocs.io/en/stable/topics/observedMutations
    region_definition : str, optional
        passed to shazam's `observedMutations`. https://shazam.readthedocs.io/en/stable/topics/IMGT_SCHEMES/
    mutation_definition : str, optional
        passed to shazam's `observedMutations`. https://shazam.readthedocs.io/en/stable/topics/MUTATION_SCHEMES/
    frequency
        whether to return the results a frequency or counts. Default is True (frequency).
    combine
        whether to return the results for replacement and silent mutations separately (False). Default is True (sum).

    Returns
    -------
    `Dandelion` object with updated `.metadata` slot.
    """
    start = logg.info('Quantifying mutations')
    try:
        from rpy2.robjects.packages import importr
        from rpy2.rinterface import NULL
        from rpy2.robjects import pandas2ri
    except:
        raise (ImportError(
            "Unable to initialise R instance. Please run this separately through R with Shazam's tutorial."
        ))

    sh = importr('shazam')
    base = importr('base')
    if self.__class__ == Dandelion:
        dat = load_data(self.data)
    elif self.__class__ == pd.DataFrame or os.path.isfile(self):
        dat = load_data(self)
    else:
        raise ValueError("{} object/file not found.".format(self))
    pandas2ri.activate()
    warnings.filterwarnings("ignore")

    if sequence_column is None:
        seq_ = 'sequence_alignment'
    else:
        seq_ = sequence_column

    if germline_column is None:
        germline_ = 'germline_alignment_d_mask'
    else:
        germline_ = germline_column

    if region_definition is None:
        reg_d = NULL
    else:

        reg_d = base.get(region_definition)

    if mutation_definition is None:
        mut_d = NULL
    else:
        mut_d = base.get(mutation_definition)

    if split_locus is False:
        dat = dat.where(dat.isna(), dat.astype(str))
        try:
            dat_r = pandas2ri.py2rpy(dat)
        except:
            dat = dat.astype(str)
            dat_r = pandas2ri.py2rpy(dat)

        results = sh.observedMutations(dat_r,
                                       sequenceColumn=seq_,
                                       germlineColumn=germline_,
                                       regionDefinition=reg_d,
                                       mutationDefinition=mut_d,
                                       frequency=frequency,
                                       combine=combine)
        # pd_df = pandas2ri.rpy2py_dataframe(results)
        pd_df = results.copy()
    else:
        dat_h = dat[dat['locus'] == 'IGH']
        dat_l = dat[dat['locus'].isin(['IGK', 'IGL'])]

        dat_h = dat_h.where(dat_h.isna(), dat_h.astype(str))
        try:
            dat_h_r = pandas2ri.py2rpy(dat_h)
        except:
            dat_h = dat_h.astype(str)
            dat_h_r = pandas2ri.py2rpy(dat_h)

        dat_l = dat_l.where(dat_l.isna(), dat_l.astype(str))
        try:
            dat_l_r = pandas2ri.py2rpy(dat_l)
        except:
            dat_l = dat_l.astype(str)
            dat_l_r = pandas2ri.py2rpy(dat_l)

        results_h = sh.observedMutations(dat_h_r,
                                         sequenceColumn=seq_,
                                         germlineColumn=germline_,
                                         regionDefinition=reg_d,
                                         mutationDefinition=mut_d,
                                         frequency=frequency,
                                         combine=combine)
        results_l = sh.observedMutations(dat_l_r,
                                         sequenceColumn=seq_,
                                         germlineColumn=germline_,
                                         regionDefinition=reg_d,
                                         mutationDefinition=mut_d,
                                         frequency=frequency,
                                         combine=combine)
        pd_df = pd.concat([results_h, results_l])

    pd_df.set_index('sequence_id', inplace=True, drop=False)
    # this doesn't actually catch overwritten columns
    cols_to_return = pd_df.columns.difference(dat.columns)
    if len(cols_to_return) < 1:
        cols_to_return = list(
            filter(re.compile("mu_.*").match, [c for c in pd_df.columns]))
    else:
        cols_to_return = cols_to_return

    res = {}
    if self.__class__ == Dandelion:
        for x in cols_to_return:
            res[x] = list(pd_df[x])
            # TODO: str will make it work for the back and forth conversion with rpy2. but maybe can use a better option
            self.data[x] = [str(r) for r in res[x]]
        if split_locus is False:
            metadata_ = self.data[['cell_id'] + list(cols_to_return)]
        else:
            metadata_ = self.data[['locus', 'cell_id'] + list(cols_to_return)]

        for x in cols_to_return:
            metadata_[x] = metadata_[x].astype(np.float32)

        if split_locus is False:
            metadata_ = metadata_.groupby('cell_id').sum()
        else:
            metadata_ = metadata_.groupby(['locus', 'cell_id']).sum()
            metadatas = []
            for x in list(set(self.data['locus'])):
                tmp = metadata_.iloc[
                    metadata_.index.isin([x], level='locus'), :]
                tmp.index = tmp.index.droplevel()
                tmp.columns = [c + '_' + str(x) for c in tmp.columns]
                metadatas.append(tmp)
            metadata_ = functools.reduce(
                lambda x, y: pd.merge(
                    x, y, left_index=True, right_index=True, how='outer'),
                metadatas)

        metadata_.index.name = None

        if self.metadata is None:
            self.metadata = metadata_
        else:
            for x in metadata_.columns:
                self.metadata[x] = pd.Series(metadata_[x])
        logg.info(' finished',
                  time=start,
                  deep=('Updated Dandelion object: \n'
                        '   \'data\', contig-indexed clone table\n'
                        '   \'metadata\', cell-indexed clone table\n'))
    else:
        for x in cols_to_return:
            res[x] = list(pd_df[x])
            # TODO: str will make it work for the back and forth conversion with rpy2. but maybe can use a better option
            dat[x] = [str(r) for r in res[x]]

        if self.__class__ == pd.DataFrame:
            logg.info(' finished', time=start, deep=('Returning DataFrame\n'))
            return (dat)
        elif os.path.isfile(self):
            logg.info(' finished',
                      time=start,
                      deep=('saving DataFrame at {}\n'.format(str(self))))
            dat.to_csv(self, sep='\t', index=False)


def calculate_threshold(self: Union[Dandelion, pd.DataFrame, str],
                        manual_threshold: Union[None, float] = None,
                        model: Union[None,
                                     Literal["ham", "aa", "hh_s1f", "hh_s5f",
                                             "mk_rs1nf", "hs1f_compat",
                                             "m1n_compat"]] = None,
                        normalize_method: Union[None, Literal['len']] = None,
                        threshold_method: Union[None,
                                                Literal['gmm',
                                                        'density']] = None,
                        edge: Union[None, float] = None,
                        cross: Union[None, Sequence] = None,
                        subsample: Union[None, int] = None,
                        threshold_model: Union[None,
                                               Literal["norm-norm",
                                                       "norm-gamma",
                                                       "gamma-norm",
                                                       "gamma-gamma"]] = None,
                        cutoff: Union[None, Literal["optimal", "intersect",
                                                    "user"]] = None,
                        sensitivity: Union[None, float] = None,
                        specificity: Union[None, float] = None,
                        ncpu: Union[None, int] = None,
                        plot: bool = True,
                        plot_group: Union[None, str] = None,
                        figsize: Tuple[Union[int, float],
                                       Union[int, float]] = (4.5, 2.5),
                        *args) -> Dandelion:
    """
    Calculating nearest neighbor distances for tuning clonal assignment with `shazam`.

    <https://shazam.readthedocs.io/en/stable/vignettes/DistToNearest-Vignette/>`__.
    Runs the following:

    distToNearest
        Get non-zero distance of every heavy chain (IGH) sequence (as defined by sequenceColumn) to its nearest sequence
        in a partition of heavy chains sharing the same V gene, J gene, and junction length (VJL), or in a partition of
        single cells with heavy chains sharing the same heavy chain VJL combination, or of single cells with heavy and
        light chains sharing the same heavy chain VJL and light chain VJL combinations.
    findThreshold
        automtically determines an optimal threshold for clonal assignment of Ig sequences using a vector of nearest
        neighbor distances. It provides two alternative methods using either a Gamma/Gaussian Mixture Model fit
        (threshold_method="gmm") or kernel density fit (threshold_method="density").

    Parameters
    ----------
    self : Dandelion, DataFrame, str
        `Dandelion` object, pandas `DataFrame` in changeo/airr format, or file path to changeo/airr file after clones
        have been determined.
    manual_threshold : float, optional
        value to manually plot in histogram.
    model : str, optional
        underlying SHM model, which must be one of "ham","aa","hh_s1f","hh_s5f","mk_rs1nf","hs1f_compat","m1n_compat".
    normalize_method : str, optional
        method of normalization. The default is "len", which divides the distance by the length of the sequence group.
        If "none" then no normalization if performed.
    threshold_method : str, optional
        string defining the method to use for determining the optimal threshold. One of "gmm" or "density".
    edge : float, optional
        upper range as a fraction of the data density to rule initialization of Gaussian fit parameters.
        Default value is 0.9 (or 90). Applies only when threshold_method="density".
    cross : Sequence, optional
        supplementary nearest neighbor distance vector output from distToNearest for initialization of the Gaussian fit
        parameters. Applies only when method="gmm".
    subsample : int, optional
        maximum number of distances to subsample to before threshold detection.
    threshold_model : str, optional
        allows the user to choose among four possible combinations of fitting curves: "norm-norm", "norm-gamma",
        "gamma-norm", and "gamma-gamma". Applies only when method="gmm".
    cutoff : str, optional
        method to use for threshold selection: the optimal threshold "optimal", the intersection point of the two fitted
        curves "intersect", or a value defined by user for one of the sensitivity or specificity "user". Applies only
        when method="gmm".
    sensitivity : float, optional
        sensitivity required. Applies only when method="gmm" and cutoff="user".
    specificity : float, optional
        specificity required. Applies only when method="gmm" and cutoff="user".
    ncpu : int, optional
        number of cpus for parallelization. Default is all available cpus.
    plot : bool
        whether or not to return plot.
    plot_group : str, optional
        determines the fill color and facets.
    figsize : Tuple[Union[int,float], Union[int,float]]
        size of plot. Default is (4.5, 2.5).
    *args
        passed to shazam's `distToNearest <https://shazam.readthedocs.io/en/stable/topics/distToNearest/>`__.

    Returns
    -------
        `Dandelion` object object with distance threshold value in `.threshold`.

        If plot = True,plotnine plot showing histogram of length normalized ham model distance threshold.
    """
    start = logg.info('Calculating threshold')
    try:
        from rpy2.robjects.packages import importr
        from rpy2.rinterface import NULL
        from rpy2.robjects import pandas2ri, FloatVector
    except:
        raise (ImportError(
            "Unable to initialise R instance. Please run this separately through R with Shazam's tutorial."
        ))

    if self.__class__ == Dandelion:
        dat = load_data(self.data)
    elif self.__class__ == pd.DataFrame or os.path.isfile(str(self)):
        dat = load_data(self)
        warnings.filterwarnings("ignore")
    sh = importr('shazam')
    pandas2ri.activate()
    if 'v_call_genotyped' in dat.columns:
        v_call = 'v_call_genotyped'
    else:
        v_call = 'v_call'
    if model is None:
        model_ = 'ham'
    else:
        model_ = model
    if normalize_method is None:
        norm_ = 'len'
    else:
        norm_ = normalize_method
    if threshold_method is None:
        threshold_method_ = "density"
    else:
        threshold_method_ = threshold_method
    if subsample is None:
        subsample_ = NULL
    else:
        subsample_ = subsample
    if ncpu is None:
        ncpu_ = multiprocessing.cpu_count() - 1
    else:
        ncpu_ = ncpu
    dat_h = dat[dat['locus'] == 'IGH']
    try:
        dat_h_r = pandas2ri.py2rpy(dat_h)
    except:
        dat_h = dat_h.astype(str)
        dat_h_r = pandas2ri.py2rpy(dat_h)
    dist_ham = sh.distToNearest(dat_h_r,
                                vCallColumn=v_call,
                                model=model_,
                                normalize=norm_,
                                nproc=ncpu_,
                                *args)
    # Find threshold using density method
    dist = np.array(dist_ham['dist_nearest'])
    if threshold_method_ == 'density':
        if edge is None:
            edge_ = 0.9
        else:
            edge_ = edge
        dist_threshold = sh.findThreshold(FloatVector(dist[~np.isnan(dist)]),
                                          method=threshold_method_,
                                          subsample=subsample_,
                                          edge=edge_)
        threshold = np.array(dist_threshold.slots['threshold'])[0]
        if np.isnan(threshold):
            print(
                "      Threshold method 'density' did not return with any values. Switching to method = 'gmm'."
            )
            threshold_method_ = 'gmm'
            if threshold_model is None:
                threshold_model_ = "gamma-gamma"
            else:
                threshold_model_ = threshold_model
            if cross is None:
                cross_ = NULL
            else:
                cross_ = cross
            if cutoff is None:
                cutoff_ = 'optimal'
            else:
                cutoff_ = cutoff
            if sensitivity is None:
                sen_ = NULL
            else:
                sen_ = sensitivity
            if specificity is None:
                spc_ = NULL
            else:
                spc_ = specificity
            dist_threshold = sh.findThreshold(FloatVector(
                dist[~np.isnan(dist)]),
                method=threshold_method_,
                model=threshold_model_,
                cross=cross_,
                subsample=subsample_,
                cutoff=cutoff_,
                sen=sen_,
                spc=spc_)
            threshold = np.array(dist_threshold.slots['threshold'])[0]
    else:
        if threshold_model is None:
            threshold_model_ = "gamma-gamma"
        else:
            threshold_model_ = threshold_model
        if cross is None:
            cross_ = NULL
        else:
            cross_ = cross
        if cutoff is None:
            cutoff_ = 'optimal'
        else:
            cutoff_ = cutoff
        if sensitivity is None:
            sen_ = NULL
        else:
            sen_ = sensitivity
        if specificity is None:
            spc_ = NULL
        else:
            spc_ = specificity
        dist_threshold = sh.findThreshold(FloatVector(dist[~np.isnan(dist)]),
                                          method=threshold_method_,
                                          model=threshold_model_,
                                          cross=cross_,
                                          subsample=subsample_,
                                          cutoff=cutoff_,
                                          sen=sen_,
                                          spc=spc_)
        threshold = np.array(dist_threshold.slots['threshold'])[0]
    if np.isnan(threshold):
        raise ValueError(
            "Automatic thresholding failed. Please visually inspect the resulting distribution fits"
            + " and choose a threshold value manually.")
    # dist_ham = pandas2ri.rpy2py_dataframe(dist_ham)

    if manual_threshold is None:
        tr = threshold
    else:
        tr = manual_threshold

    if plot:
        options.figure_size = figsize
        if plot_group is None:
            plot_group = 'sample_id'
        else:
            plot_group = plot_group

        print((ggplot(dist_ham, aes('dist_nearest', fill=str(plot_group))) +
               theme_bw() + xlab("Grouped Hamming distance") + ylab("Count") +
               geom_histogram(binwidth=0.01) + geom_vline(
                   xintercept=tr, linetype="dashed", color="blue", size=0.5) +
               annotate('text',
                        x=tr + 0.02,
                        y=10,
                        label='Threshold:\n' + str(np.around(tr, decimals=2)),
                        size=8,
                        color='Blue') +
               facet_wrap('~' + str(plot_group), scales="free_y") +
               theme(legend_position='none')))
    else:
        print("Automatic Threshold : " +
              str(np.around(threshold, decimals=2)) + "\n method = " +
              str(threshold_method_))
    if self.__class__ == Dandelion:
        self.threshold = tr
        logg.info(
            ' finished',
            time=start,
            deep=
            ('Updated Dandelion object: \n'
             '   \'threshold\', threshold value for tuning clonal assignment\n'
             ))
    else:
        output = Dandelion(dat)
        output.threshold = tr
        return (output)
