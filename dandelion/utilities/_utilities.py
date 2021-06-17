#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 14:01:32
# @Last Modified by:   Kelvin
# @Last Modified time: 2021-06-16 20:32:07

import os
from collections import defaultdict, Iterable
import pandas as pd
import numpy as np
from subprocess import run
import re
from typing import Sequence, Tuple, Dict, Union
try:
    from typing import Literal
except ImportError:
    try:
        from typing_extensions import Literal
    except ImportError:

        class LiteralMeta(type):
            def __getitem__(cls, values):
                if not isinstance(values, tuple):
                    values = (values, )
                return type("Literal_", (Literal, ), dict(__args__=values))

        class Literal(metaclass=LiteralMeta):
            pass


class Tree(defaultdict):
    '''
    Create a recursive defaultdict
    '''

    def __init__(self, value=None):
        super(Tree, self).__init__(Tree)
        self.value = value


def dict_from_table(meta: pd.DataFrame, columns: Tuple[str, str]) -> Dict:
    """
    Generates a dictionary from a dataframe
    Parameters
    ----------
    meta
        pandas dataframe or file path
    columns
        column names in dataframe

    Returns
    -------
        dictionary
    """
    if (isinstance(meta, pd.DataFrame)) & (columns is not None):
        meta_ = meta
        if len(columns) == 2:
            sample_dict = dict(zip(meta_[columns[0]], meta_[columns[1]]))
    elif (os.path.isfile(str(meta))) & (columns is not None):
        meta_ = pd.read_csv(meta, sep='\t', dtype='object')
        if len(columns) == 2:
            sample_dict = dict(zip(meta_[columns[0]], meta_[columns[1]]))

    sample_dict = clean_nan_dict(sample_dict)
    return (sample_dict)


def clean_nan_dict(d: Dict) -> Dict:
    """
    Parameters
    ----------
    d : Dict
        dictionary

    Returns
    -------
        dictionary with no NAs.
    """

    return {k: v for k, v in d.items() if v is not np.nan}


def flatten(l: Sequence) -> Sequence:
    """
    Parameters
    ----------
    l : Sequence
        a list-in-list list

    Returns
    -------
        a flattened list.
    """
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el


def makeblastdb(ref: str):
    """
    Run makeblastdb on constant region fasta file.

    Wrapper for makeblastdb.

    Parameters
    ----------
    ref : str
        constant region fasta file
    Returns
    -------


    """

    cmd = ['makeblastdb', '-dbtype', 'nucl', '-parse_seqids', '-in', ref]
    run(cmd)


def bh(pvalues: np.array) -> np.array:
    """
    Computes the Benjamini-Hochberg FDR correction.

    Parameters
    ----------
        pvalues : np.array
            array of p-values to correct
    Returns
    -------
        np.array of corrected p-values
    """
    n = int(pvalues.shape[0])
    new_pvalues = np.empty(n)
    values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
    values.sort()
    values.reverse()
    new_values = []
    for i, vals in enumerate(values):
        rank = n - i
        pvalue, index = vals
        new_values.append((n / rank) * pvalue)
    for i in range(0, int(n) - 1):
        if new_values[i] < new_values[i + 1]:
            new_values[i + 1] = new_values[i]
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]
    return new_pvalues


def is_categorical(array_like) -> bool:
    return array_like.dtype.name == 'category'


def type_check(dataframe, key) -> bool:
    return dataframe[key].dtype == str or dataframe[
        key].dtype == object or is_categorical(
            dataframe[key]) or dataframe[key].dtype == bool


def isGZIP(filename: str) -> bool:
    if filename.split('.')[-1] == 'gz':
        return True
    return False


def isBZIP(filename: str) -> bool:
    if filename.split('.')[-1] == 'pbz2':
        return True
    return False


def check_filepath(s, filename_prefix=None, endswith=None, subdir=None):
    if filename_prefix is None:
        filename_pre = 'filtered'
    else:
        filename_pre = filename_prefix

    if endswith is None:
        ends_with = ''
    else:
        ends_with = endswith

    filePath = None
    if os.path.isfile(str(s)) and str(s).endswith(ends_with):
        filePath = s
    elif os.path.isdir(str(s)):
        files = os.listdir(s)
        for file in files:
            out_ = s.rstrip('/') + '/'
            if os.path.isdir(out_ + os.path.basename(file)):
                if file == 'dandelion':
                    if subdir is None:
                        out_ = out_ + os.path.basename(file) + '/'
                    else:
                        out_ = out_ + os.path.basename(
                            file) + '/' + subdir + '/'
                    for x in os.listdir(out_):
                        if x.endswith(ends_with):
                            if str(x).split(
                                    ends_with)[0] == filename_pre + '_contig':
                                filePath = out_ + x
                            else:
                                continue
            else:
                continue
    return (filePath)


def check_fastapath(fasta, filename_prefix=None):
    if filename_prefix is None:
        filename_pre = 'filtered'
    else:
        filename_pre = filename_prefix

    filePath = None
    if os.path.isfile(str(fasta)) and str(fasta).endswith('.fasta'):
        filePath = fasta
    elif os.path.isdir(str(fasta)):
        files = os.listdir(fasta)
        for file in files:
            out_ = fasta.rstrip('/') + '/'
            if str(file).endswith(".fasta"):
                if str(file).split('.fasta')[0] == filename_pre + '_contig':
                    filePath = out_ + os.path.basename(file)
                else:
                    continue
    return (filePath)


def change_file_location(data: Sequence,
                         filename_prefix: Union[None, Sequence, str] = None):
    """
    Move file from tmp folder to dandelion folder.

    Only used for TCR data.

    Parameters
    ----------
    data : Sequence
        list of data folders containing the .tsv files. if provided as a single string, it will first be converted to a
        list; this allows for the function to be run on single/multiple samples.
    filename_prefix : str, optional
        list of prefixes of file names preceding '_contig'. None defaults to 'filtered'.

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

    filePath = None

    for i in range(0, len(data)):
        filePath = check_filepath(data[i],
                                  filename_prefix=filename_prefix[i],
                                  endswith=informat_dict[fileformat],
                                  subdir='tmp')
        if filePath is None:
            raise OSError(
                'Path to .tsv file for {} is unknown. '.format(data[i]) +
                'Please specify path to reannotated .tsv file or folder containing reannotated .tsv file.'
            )

        cmd = ['rsync', '-azvh', filePath, filePath.rsplit('/', 2)[0]]

        run(cmd)


class FilterContigs:
    def __init__(self, data):
        self.dat = data
        self.poor_qual = []
        self.h_doublet = []
        self.l_doublet = []
        self.drop_contig = []
        self.umi_adjustment = {}

    def run_scan(self, b, rescue_vdj, umi_foldchange_cutoff,
                 filter_poorqualitycontig, v_dict, j_dict, c_dict):
        """Main workhorse of filter_contig."""
        h = Tree()
        l = Tree()
        h_umi = Tree()
        l_umi = Tree()
        h_seq = Tree()
        l_seq = Tree()
        h_ccall = Tree()

        hc_id = list(self.dat[
            (self.dat['cell_id'].isin([b]))
            & (self.dat['locus'].isin(['IGH', 'TRB', 'TRD']))]['sequence_id'])
        hc_umi = [
            int(x)
            for x in self.dat[(self.dat['cell_id'].isin([b]))
                              & (self.dat['locus'].isin(['IGH', 'TRB', 'TRD'])
                                 )]['duplicate_count']
        ]
        if 'sequence_alignment' in self.dat:
            hc_seq = [
                x for x in self.dat[
                    (self.dat['cell_id'].isin([b]))
                    & (self.dat['locus'].isin(['IGH', 'TRB', 'TRD']))]
                ['sequence_alignment']
            ]
        hc_ccall = [
            x for x in self.dat[
                (self.dat['cell_id'].isin([b]))
                & (self.dat['locus'].isin(['IGH', 'TRB', 'TRD']))]['c_call']
        ]

        lc_id = list(
            self.dat[(self.dat['cell_id'].isin([b]))
                     & (self.dat['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG']))]
            ['sequence_id'])
        lc_umi = [
            int(x) for x in self.dat[
                (self.dat['cell_id'].isin([b]))
                & (self.dat['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG']))]
            ['duplicate_count']
        ]
        if 'sequence_alignment' in self.dat:
            lc_seq = [
                x for x in self.dat[
                    (self.dat['cell_id'].isin([b]))
                    & (self.dat['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG']))]
                ['sequence_alignment']
            ]

        h[b] = hc_id
        h_umi[b] = hc_umi
        if 'sequence_alignment' in self.dat:
            h_seq[b] = hc_seq
        h_ccall[b] = hc_ccall

        l[b] = lc_id
        l_umi[b] = lc_umi
        if 'sequence_alignment' in self.dat:
            l_seq[b] = lc_seq

        # marking doublets defined by heavy chains
        if len(h[b]) > 1:
            if 'sequence_alignment' in self.dat:
                if len(list(set(h_seq[b]))) == 1:
                    highest_umi_h = max(h_umi[b])
                    highest_umi_h_idx = [
                        i for i, j in enumerate(h_umi[b]) if j == highest_umi_h
                    ]
                    keep_index_h = highest_umi_h_idx[0]
                    self.drop_contig.append(h[b][:keep_index_h] +
                                            h[b][keep_index_h + 1:])
                    keep_hc_contig = h[b][keep_index_h]
                    self.dat.at[keep_hc_contig, 'duplicate_count'] = int(
                        np.sum(h_umi[b][:keep_index_h] +
                               h_umi[b][keep_index_h + 1:]))
                    hc_id = list(
                        self.dat[(self.dat['cell_id'].isin([b]))
                                 & (self.dat['locus'].isin(
                                     ['IGH', 'TRB', 'TRD']))]['sequence_id'])
                    hc_umi = [
                        int(x) for x in self.dat[
                            (self.dat['cell_id'].isin([b]))
                            & (self.dat['locus'].isin(['IGH', 'TRB', 'TRD']))]
                        ['duplicate_count']
                    ]
                    h[b] = hc_id
                    h_umi[b] = hc_umi
                    h_seq[b] = hc_seq
            if len(h[b]) > 1:
                if rescue_vdj:
                    highest_umi_h = max(h_umi[b])
                    # lowest_umi_h = min(h_umi[b])
                    highest_umi_idx = [
                        i for i, j in enumerate(h_umi[b]) if j == highest_umi_h
                    ]
                    keep_index_h = highest_umi_idx[0]
                    umi_test = [
                        highest_umi_h / x < umi_foldchange_cutoff
                        for x in h_umi[b][:keep_index_h] +
                        h_umi[b][keep_index_h + 1:]
                    ]
                    sum_umi = sum(h_umi[b])
                    other_umi_idx = [
                        i for i, j in enumerate(h_umi[b]) if j != highest_umi_h
                    ]
                    if 'IGHM' and 'IGHD' in h_ccall[b]:
                        if all(cc_ == 'IGHM' or cc_ == 'IGHD'
                               for cc_ in h_ccall[b]):
                            pass
                        else:
                            if len(highest_umi_idx) > 1:
                                self.h_doublet.append(b)
                            if sum_umi < 4:
                                self.h_doublet.append(b)
                            if any(umi_test):
                                self.h_doublet.append(b)
                            if len(highest_umi_idx) == 1:
                                other_umi_idx = [
                                    i for i, j in enumerate(h_umi[b])
                                    if j != highest_umi_h
                                ]
                                umi_test_ = [
                                    highest_umi_h / x >= umi_foldchange_cutoff
                                    for x in h_umi[b][:keep_index_h] +
                                    h_umi[b][keep_index_h + 1:]
                                ]
                                umi_test_dict = dict(
                                    zip(other_umi_idx, umi_test_))
                                for otherindex in umi_test_dict:
                                    if umi_test_dict[otherindex]:
                                        self.drop_contig.append(
                                            h[b][otherindex])
                    else:
                        if len(highest_umi_idx) > 1:
                            self.h_doublet.append(b)
                        if sum_umi < 4:
                            self.h_doublet.append(b)
                        if any(umi_test):
                            self.h_doublet.append(b)
                        if len(highest_umi_idx) == 1:
                            other_umi_idx = [
                                i for i, j in enumerate(h_umi[b])
                                if j != highest_umi_h
                            ]
                            umi_test_ = [
                                highest_umi_h / x >= umi_foldchange_cutoff
                                for x in h_umi[b][:keep_index_h] +
                                h_umi[b][keep_index_h + 1:]
                            ]
                            umi_test_dict = dict(zip(other_umi_idx, umi_test_))
                            for otherindex in umi_test_dict:
                                if umi_test_dict[otherindex]:
                                    self.drop_contig.append(h[b][otherindex])
                else:
                    self.h_doublet.append(b)
        if len(l[b]) > 1:
            if 'sequence_alignment' in self.dat:
                if len(list(set(l_seq[b]))) == 1:
                    highest_umi_l = max(l_umi[b])
                    highest_umi_l_idx = [
                        i for i, j in enumerate(l_umi[b]) if j == highest_umi_l
                    ]
                    keep_index_l = highest_umi_l_idx[0]
                    self.drop_contig.append(l[b][:keep_index_l] +
                                            l[b][keep_index_l + 1:])
                    keep_lc_contig = l[b][keep_index_l]
                    self.dat.at[keep_lc_contig, 'duplicate_count'] = int(
                        np.sum(l_umi[b][:keep_index_l] +
                               l_umi[b][keep_index_l + 1:]))
                    lc_id = list(self.dat[
                        (self.dat['cell_id'].isin([b]))
                        & (self.dat['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG'])
                           )]['sequence_id'])
                    lc_umi = [
                        int(x)
                        for x in self.dat[(self.dat['cell_id'].isin([b]))
                                          & (self.dat['locus'].isin(
                                              ['IGK', 'IGL', 'TRA', 'TRG']))]
                        ['duplicate_count']
                    ]
                    l[b] = lc_id
                    l_umi[b] = lc_umi
                    l_seq[b] = lc_seq
            if len(list(set(l[b]))) > 1:
                # also apply the same cut off to multiple light chains
                highest_umi_l = max(l_umi[b])
                highest_umi_l_idx = [
                    i for i, j in enumerate(l_umi[b]) if j == highest_umi_l
                ]
                keep_index_l = highest_umi_l_idx[0]
                other_umi_idx_l = [
                    i for i, j in enumerate(l_umi[b]) if j != highest_umi_l
                ]
                umi_test_l = [
                    highest_umi_l / x < umi_foldchange_cutoff
                    for x in l_umi[b][:keep_index_l] +
                    l_umi[b][keep_index_l + 1:]
                ]
                umi_test_dict_l = dict(zip(other_umi_idx_l, umi_test_l))
                for otherindex in umi_test_dict_l:
                    if umi_test_dict_l[otherindex]:
                        self.drop_contig.append(l[b][otherindex])
        # marking doublets defined by light chains
        if (len(h[b]) == 1) & (len(l[b]) > 1):
            self.l_doublet.append(b)
        # marking poor bcr quality, defined as those with only light chains, those
        # that were have conflicting assignment of locus and heavy/light V/J calls,
        # and also those that are missing either v or j calls.
        if len(h[b]) < 1:
            if filter_poorqualitycontig:
                self.poor_qual.append(b)
            self.drop_contig.append(l[b])
        if len(hc_id) == 1:
            v = v_dict[hc_id[0]]
            j = j_dict[hc_id[0]]
            c = c_dict[hc_id[0]]
            if v == v:
                if not re.search('IGH|TR[BD]', v):
                    if filter_poorqualitycontig:
                        self.poor_qual.append(b)
                    self.drop_contig.append(l[b])
                    self.drop_contig.append(h[b])
            else:
                if filter_poorqualitycontig:
                    self.poor_qual.append(b)
                self.drop_contig.append(l[b])
                self.drop_contig.append(h[b])
            if j == j:
                if not re.search('IGH|TR[BD]', j):
                    if filter_poorqualitycontig:
                        self.poor_qual.append(b)
                    self.drop_contig.append(l[b])
                    self.drop_contig.append(h[b])
            else:
                if filter_poorqualitycontig:
                    self.poor_qual.append(b)
                self.drop_contig.append(l[b])
                self.drop_contig.append(h[b])
            if (c == c) and (c is not None):
                if not re.search('IGH|TR[BD]', c):
                    if filter_poorqualitycontig:
                        self.poor_qual.append(b)
                    self.drop_contig.append(l[b])
                    self.drop_contig.append(h[b])
        if len(hc_id) > 1:
            for hx in hc_id:
                v = v_dict[hx]
                j = j_dict[hx]
                c = c_dict[hx]
                if v == v:
                    if not re.search('IGH|TR[BD]', v):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(b)
                        self.drop_contig.append(hx)
                if j == j:
                    if not re.search('IGH|TR[BD]', j):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(b)
                        self.drop_contig.append(hx)
                if (c == c) and (c is not None):
                    if not re.search('IGH|TR[BD]', c):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(b)
                        self.drop_contig.append(hx)
        if len(lc_id) > 0:
            for lx in lc_id:
                v = v_dict[lx]
                j = j_dict[lx]
                c = c_dict[lx]
                if v == v:
                    if j == j:
                        if re.search('IGH|TR[BD]', v):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(lx)
                        elif (re.search('IGK', v) and re.search('IGL', j)):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(lx)
                        elif (re.search('IGL', v) and re.search('IGK', j)):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(lx)
                        elif (re.search('TRA', v) and not re.search('TRA', j)):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(lx)
                        elif (re.search('TRG', v) and not re.search('TRG', j)):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(lx)
                if j == j:
                    if v == v:
                        if re.search('IGH|TR[BD]', j):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(lx)
                        elif (re.search('IGK', v) and re.search('IGL', j)):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(lx)
                        elif (re.search('IGL', v) and re.search('IGK', j)):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(lx)
                        elif (re.search('TRA', v) and not re.search('TRA', j)):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(lx)
                        elif (re.search('TRG', v) and not re.search('TRG', j)):
                            if filter_poorqualitycontig:
                                self.poor_qual.append(b)
                            self.drop_contig.append(lx)
                if (c is not None) and (c == c):
                    if re.search('IGH|TR[BD]', c):
                        if filter_poorqualitycontig:
                            self.poor_qual.append(b)
                        self.drop_contig.append(lx)

                if (v != v) or (j != j) or (v is None) or (j is None):
                    if filter_poorqualitycontig:
                        self.poor_qual.append(b)
                    self.drop_contig.append(lx)  # no/wrong annotations at all
