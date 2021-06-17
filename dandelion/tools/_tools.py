#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2020-05-13 23:22:18
# @Last Modified by:   Kelvin
# @Last Modified time: 2021-05-20 16:54:22

import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
from tqdm import tqdm
from ..utilities._utilities import *
from ..utilities._core import *
from ..utilities._io import *
from ._network import *
from collections import defaultdict
from itertools import groupby
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix
from distance import hamming
import re
import math
import networkx as nx
from time import sleep
import copy
import functools
try:
    from scanpy import logging as logg
except ImportError:
    pass
import warnings
from subprocess import run
import multiprocessing
from changeo.Gene import getGene
from anndata import AnnData
from typing import Union, Sequence, Tuple


def find_clones(self: Union[Dandelion, pd.DataFrame], identity: float = 0.85, key: Union[None, str] = None, locus: Union[None, str] = None, by_alleles: bool = False, key_added: Union[None, str] = None, recalculate_length: bool = True, productive_only: bool = True) -> Dandelion:
    """
    Find clones based on heavy chain and light chain CDR3 junction hamming distance.

    Parameters
    ----------
    self : Dandelion, DataFrame, str
        `Dandelion` object, pandas `DataFrame` in changeo/airr format, or file path to changeo/airr file after clones have been determined.
    identity : float
        Junction similarity parameter. Default 0.85
    key : str, optional
        column name for performing clone clustering. None defaults to 'junction_aa'.
    locus : str, optional
        placeholder to allow for future tcr analysis mode for performing clustering. None defaults to 'ig'.
    by_alleles : bool
        Whether or not to collapse alleles to genes. None defaults to False.
    key_added : str, optional
        If specified, this will be the column name for clones. None defaults to 'clone_id'
    recalculate_length : bool
        Whether or not to re-calculate junction length, rather than rely on parsed assignment (which occasionally is wrong). Default is True
    productive_only : bool
        Whether or not to perform clone_clustering only on productive clones.

    Returns
    -------
    `Dandelion` object with clone_id annotated in `.data` slot and `.metadata` initialized.
    """
    start = logg.info('Finding clonotypes')
    if self.__class__ == Dandelion:
        dat_ = load_data(self.data)
    else:
        dat_ = load_data(self)

    if productive_only:
        dat = dat_[dat_['productive'].isin(['T', 'True', 'TRUE', True])].copy()
    else:
        dat = dat_.copy()

    locus_dict = {'ig': 'IGH', 'tr-ab': 'TRB', 'tr-gd': 'TRG'}

    if key is None:
        key_ = 'junction_aa'  # default
    else:
        key_ = key

    if key_ not in dat.columns:
        raise ValueError("key {} not found in input table.".format(key_))

    if locus is None:
        locus_ = 'IGH'
    else:
        locus_ = locus_dict[locus]

    locus_log1_dict = {'IGH': 'IGH', 'TRB': 'TRB', 'TRD': 'TRD'}
    locus_log2_dict = {'IGH': 'IGL/IGL', 'TRB': 'TRA', 'TRD': 'TRG'}

    dat_light = dat[~(dat['locus'] == locus_)].copy()
    dat_heavy = dat[dat['locus'] == locus_].copy()

    if dat_light.shape[0] > 0:
        dump = dat_light[~(dat_light['cell_id'].isin(
            dat_heavy['cell_id']))].copy()
        if dump.shape[0] > 0:
            dat = dat[~(dat['cell_id'].isin(dump['cell_id']))].copy()
    dat_heavy = dat[dat['locus'] == locus_].copy()
    pd.set_option('mode.chained_assignment', None)

    if key_added is None:
        clone_key = 'clone_id'
    else:
        clone_key = key_added

    # retrieve the J genes and J genes
    if not by_alleles:
        if 'v_call_genotyped' in dat_heavy.columns:
            V = [re.sub('[*][0-9][0-9]', '', v)
                 for v in dat_heavy['v_call_genotyped']]
        else:
            V = [re.sub('[*][0-9][0-9]', '', v) for v in dat_heavy['v_call']]
        J = [re.sub('[*][0-9][0-9]', '', j) for j in dat_heavy['j_call']]
    else:
        if 'v_call_genotyped' in dat_heavy.columns:
            V = [v for v in dat_heavy['v_call_genotyped']]
        else:
            V = [v for v in dat_heavy['v_call']]
        J = [j for j in dat_heavy['j_call']]

    # collapse the alleles to just genes
    V = [','.join(list(set(v.split(',')))) for v in V]
    J = [','.join(list(set(j.split(',')))) for j in J]

    seq = dict(zip(dat_heavy.index, dat_heavy[key_]))
    if recalculate_length:
        seq_length = [len(str(l)) for l in dat_heavy[key_]]
    else:
        try:
            seq_length = [l for l in dat_heavy[key_ + '_length']]
        except:
            raise ValueError("{} not found in {} input table.".format(
                key_ + '_length', locus_log1_dict[locus_]))
    seq_length_dict = dict(zip(dat_heavy.index, seq_length))

    # Create a dictionary and group sequence ids with same V and J genes
    V_J = dict(zip(dat_heavy.index, zip(V, J)))
    vj_grp = defaultdict(list)
    for key, val in sorted(V_J.items()):
        vj_grp[val].append(key)
    # and now we split the groups based on lengths of the seqs
    vj_len_grp = Tree()
    seq_grp = Tree()
    for g in vj_grp:
        # first, obtain what's the unique lengths
        jlen = []
        for contig_id in vj_grp[g]:
            jlen.append(seq_length_dict[contig_id])
        setjlen = list(set(jlen))
        # then for each unique length, we add the contigs to a new tree if it matches the length
        # and also the actual seq sequences into a another one
        for s in setjlen:
            for contig_id in vj_grp[g]:
                jlen_ = seq_length_dict[contig_id]
                if jlen_ == s:
                    vj_len_grp[g][s][contig_id].value = 1
                    seq_grp[g][s][seq[contig_id]].value = 1
                    for c in [contig_id]:
                        vj_len_grp[g][s][c] = seq[c]
    clones = Tree()
    # for each seq group, calculate the hamming distance matrix
    for g in tqdm(seq_grp, desc='Finding clones based on VDJ chains '):
        for l in seq_grp[g]:
            seq_ = list(seq_grp[g][l])
            tdarray = np.array(seq_).reshape(-1, 1)
            d_mat = squareform(
                pdist(tdarray, lambda x, y: hamming(x[0], y[0])))
            # then calculate what the acceptable threshold is for each length of sequence
            tr = math.floor(int(l) * (1 - identity))
            # convert diagonal and upper triangle to zeroes
            d_mat = np.tril(d_mat)
            np.fill_diagonal(d_mat, 0)
            # get the coordinates/indices of seqs to match against the threshold later
            indices_temp = []
            indices = []
            indices_temp = [list(x) for x in np.tril_indices_from(d_mat)]
            indices = list(zip(indices_temp[0], indices_temp[1]))
            # if there's more than 1 contig, remove the diagonal
            if len(indices) > 1:
                for pairs in indices:
                    # remove diagonals
                    if pairs[0] == pairs[1]:
                        indices.remove(pairs)
            indices_j = []
            # use the coordinates/indices to retrieve the seq sequences
            for p in range(0, len(indices)):
                a1, b1 = indices[p]
                indices_j.append(seq_[a1])
                indices_j.append(seq_[b1])
            # retain only the unique sequences
            indices_j_f = list(set(indices_j))
            # convert the distance matrix to coordinate (source) and distance (target) and create it as a dictionary
            source, target = d_mat.nonzero()
            source_target = list(zip(source.tolist(), target.tolist()))
            if len(source) == 0 & len(target) == 0:
                source_target = list([(0, 0)])
            dist = {}
            for st in source_target:
                dist.update({st: d_mat[st]})

            cm1, cm2, cm3 = [], [], []
            # now to calculate which contigs to group
            tr2 = min(dist.values())
            if tr2 <= tr:
                for d in dist:
                    # if the distance is equal to the minimum distance, which is lesser than or equal to the threshold. I will add the sequence to cm1, meaning they passed and should be a clone.
                    if dist[d] == tr2:
                        cm1.append(d)
                    else:
                        # otherwise add into cm3, which later on will get split into individual clones, despite being part of the acceptable threshold for this group of sequences
                        cm3.append(d)
            else:
                # if the value is greater than the the minimum threshold, I will add the entire group into cm3, so they shouldn't be grouped together as a clone
                for d in dist:
                    cm3.append(d)
            # This is based on a simple scenario e.g. in 3x3 distance matrix where the acceptable distance threshold is 2.
            # SeqA is 1 sequence different from SeqB
            # SeqB is 2 sequences different from SeqC
            # SeqC can only be 2 sequences different from SeqA otherwise it will violate the pair-wise distance matrix calculation
            # Because the distance between SeqA and SeqB is below the threshold, and is essentially the smallest value, they will be grouped in cm1
            # While SeqC is below the acceptable threshold, really it's more different than SeqA and SeqB. Therefore, they will be grouped separately.
            # If there's more than 3 sequences in a particular distance matrix, it gets slightly complicated.
            # so i repeat the calculation to only catch those with more than 3 sequences
            if len(dist) > 3:
                for d in dist:
                    # This time, i will rely on the main threshold value; if the distance is lesser than the overall threshold, i will add it to cm2.
                    if dist[d] < tr:
                        cm2.append(d)
                    else:
                        cm3.append(d)
                # any sequences that were already binned above in cm1, it will be removed from here.
                cm2 = list(set(cm2) ^ set(cm1))
                cm3 = list(set(cm3) ^ set(cm1))

            # now to actually retrieve the seq/contigs properly through a series of appending and removing
            j_list1 = []
            if len(cm1) > 0:
                for i in range(0, len(cm1)):
                    a, b = cm1[i]
                    j_list1.append(seq_[a])
                    j_list1.append(seq_[b])
                j_list1 = list(set(j_list1))

            j_list2 = []
            if len(cm3) > 0:
                for i in range(0, len(cm2)):
                    a, b = cm2[i]
                    j_list2.append(seq_[a])
                    j_list2.append(seq_[b])
                j_list2 = list(set(j_list2))

            j_list3 = []
            if len(cm3) > 0:
                for i in range(0, len(cm3)):
                    a, b = cm3[i]
                    j_list3.append(seq_[a])
                    j_list3.append(seq_[b])
                j_list3 = list(set(j_list3))

            # this is to catch the overlapping seqs appearing in the lists because of the sequential appending
            for jl3_1 in j_list1:
                if jl3_1 in j_list3:
                    j_list3.remove(jl3_1)
            for jl2_1 in j_list1:
                if jl2_1 in j_list2:
                    j_list2.remove(jl2_1)
            for jl3_2 in j_list2:
                if jl3_2 in j_list3:
                    j_list3.remove(jl3_2)

            j_list3 = [i.split() for i in list(set(j_list3))]

            if len(j_list1) > 0:
                clones[g][l][str(0)] = j_list1
            if len(j_list2) > 0:
                clones[g][l][str(1)] = j_list2
            if len(j_list3) > 0:
                for c in range(0, len(j_list3)):
                    # the +2 here is so that the numbers come up after 1. It doesn't matter because i will reformat the clone ID numbers later
                    clones[g][l][str(c + 2)] = j_list3[c]

    clone_dict = {}
    # now to retrieve the contig ids that are grouped together
    cid = Tree()
    for g in clones:
        for l in clones[g]:
            # retrieve the clone 'numbers'
            for c in clones[g][l]:
                grp_seq = clones[g][l][c]
                for key, value in vj_len_grp[g][l].items():
                    if value in grp_seq:
                        cid[g][l][c][key].value = 1
    # rename clone ids - get dictionaries step by step
    first_key = []
    for k1 in cid.keys():
        first_key.append(k1)
    first_key = list(set(first_key))
    first_key_dict = dict(zip(first_key, range(1, len(first_key) + 1)))
    # and now for the middle key
    for g in cid:
        second_key = []
        for k2 in cid[g].keys():
            second_key.append(k2)
        second_key = list(set(second_key))
        second_key_dict = dict(zip(second_key, range(1, len(second_key) + 1)))
        for l in cid[g]:
            # and now for the last key
            third_key = []
            for k3 in cid[g][l].keys():
                third_key.append(k3)
            third_key = list(set(third_key))
            third_key_dict = dict(zip(third_key, range(1, len(third_key) + 1)))
            for key, value in dict(cid[g][l]).items():
                vL = []
                for v in value:
                    if type(v) is int:
                        break
                    # instead of converting to another tree, i will just make it a dictionary
                    clone_dict[v] = str(
                        first_key_dict[g]) + '_' + str(second_key_dict[l]) + '_' + str(third_key_dict[key])
    # add it to the original dataframes
    dat_heavy[clone_key] = pd.Series(clone_dict)
    dat[clone_key] = pd.Series(dat_heavy[clone_key])

    dat_light_c = dat[~(dat['locus'] == locus_)].copy()
    if dat_light_c.shape[0] != 0:
        if not by_alleles:
            if 'v_call_genotyped' in dat_light_c.columns:
                Vlight = [re.sub('[*][0-9][0-9]', '', v)
                          for v in dat_light_c['v_call_genotyped']]
            else:
                Vlight = [re.sub('[*][0-9][0-9]', '', v)
                          for v in dat_light_c['v_call']]
            Jlight = [re.sub('[*][0-9][0-9]', '', j)
                      for j in dat_light_c['j_call']]
        else:
            if 'v_call_genotyped' in dat_light_c.columns:
                Vlight = [v for v in dat_light_c['v_call_genotyped']]
            else:
                Vlight = [v for v in dat_light_c['v_call']]
            Jlight = [j for j in dat_light_c['j_call']]
        # collapse the alleles to just genes
        Vlight = [','.join(list(set(v.split(',')))) for v in Vlight]
        Jlight = [','.join(list(set(j.split(',')))) for j in Jlight]
        seq = dict(zip(dat_light_c.index, dat_light_c[key_]))
        if recalculate_length:
            seq_length = [len(str(l)) for l in dat_light_c[key_]]
        else:
            try:
                seq_length = [len(str(l))
                              for l in dat_light_c[key_ + '_length']]
            except:
                raise ValueError("{} not found in {} input table.".format(
                    key_ + '_length', locus_log2_dict[locus_]))
        seq_length_dict = dict(zip(dat_light_c.index, seq_length))
        # Create a dictionary and group sequence ids with same V and J genes
        V_Jlight = dict(zip(dat_light_c.index, zip(Vlight, Jlight)))
        vj_lightgrp = defaultdict(list)
        for key, val in sorted(V_Jlight.items()):
            vj_lightgrp[val].append(key)
        # and now we split the groups based on lengths of the seqs
        vj_len_lightgrp = Tree()
        seq_lightgrp = Tree()
        for g in vj_lightgrp:
            # first, obtain what's the unique lengths
            jlen = []
            for contig_id in vj_lightgrp[g]:
                jlen.append(seq_length_dict[contig_id])
            setjlen = list(set(jlen))
            # then for each unique length, we add the contigs to a new tree if it matches the length
            # and also the actual seq sequences into a another one
            for s in setjlen:
                for contig_id in vj_lightgrp[g]:
                    jlen_ = seq_length_dict[contig_id]
                    if jlen_ == s:
                        vj_len_lightgrp[g][s][contig_id].value = 1
                        seq_lightgrp[g][s][seq[contig_id]].value = 1
                        for c in [contig_id]:
                            vj_len_lightgrp[g][s][c] = seq[c]
        clones_light = Tree()
        for g in seq_lightgrp:
            for l in seq_lightgrp[g]:
                seq_ = list(seq_lightgrp[g][l])
                tdarray = np.array(seq_).reshape(-1, 1)
                d_mat = squareform(
                    pdist(tdarray, lambda x, y: hamming(x[0], y[0])))
                tr = math.floor(int(l) * (1 - identity))
                d_mat = np.tril(d_mat)
                np.fill_diagonal(d_mat, 0)
                indices_temp = []
                indices = []
                indices_temp = [list(x)
                                for x in np.tril_indices_from(d_mat)]
                indices = list(zip(indices_temp[0], indices_temp[1]))
                if len(indices) > 1:
                    for pairs in indices:
                        if pairs[0] == pairs[1]:
                            indices.remove(pairs)
                indices_j = []
                for p in range(0, len(indices)):
                    a1, b1 = indices[p]
                    indices_j.append(seq_[a1])
                    indices_j.append(seq_[b1])
                indices_j_f = list(set(indices_j))
                source, target = d_mat.nonzero()
                source_target = list(
                    zip(source.tolist(), target.tolist()))
                if len(source) == 0 & len(target) == 0:
                    source_target = list([(0, 0)])
                dist = {}
                for st in source_target:
                    dist.update({st: d_mat[st]})
                cm1 = []
                cm2 = []
                cm3 = []
                tr2 = min(dist.values())
                if tr2 <= tr:
                    for d in dist:
                        if dist[d] == tr2:
                            cm1.append(d)
                        else:
                            cm3.append(d)
                else:
                    for d in dist:
                        cm3.append(d)
                if len(dist) > 3:
                    for d in dist:
                        if dist[d] < tr:
                            cm2.append(d)
                        else:
                            cm3.append(d)
                    cm2 = list(set(cm2) ^ set(cm1))
                    cm3 = list(set(cm3) ^ set(cm1))
                j_list1 = []
                if len(cm1) > 0:
                    for i in range(0, len(cm1)):
                        a, b = cm1[i]
                        j_list1.append(seq_[a])
                        j_list1.append(seq_[b])
                    j_list1 = list(set(j_list1))
                j_list2 = []
                if len(cm3) > 0:
                    for i in range(0, len(cm2)):
                        a, b = cm2[i]
                        j_list2.append(seq_[a])
                        j_list2.append(seq_[b])
                    j_list2 = list(set(j_list2))
                j_list3 = []
                if len(cm3) > 0:
                    for i in range(0, len(cm3)):
                        a, b = cm3[i]
                        j_list3.append(seq_[a])
                        j_list3.append(seq_[b])
                    j_list3 = list(set(j_list3))
                for jl3_1 in j_list1:
                    if jl3_1 in j_list3:
                        j_list3.remove(jl3_1)
                for jl2_1 in j_list1:
                    if jl2_1 in j_list2:
                        j_list2.remove(jl2_1)
                for jl3_2 in j_list2:
                    if jl3_2 in j_list3:
                        j_list3.remove(jl3_2)
                j_list3 = [i.split() for i in list(set(j_list3))]
                if len(j_list1) > 0:
                    clones_light[g][l][str(0)] = j_list1
                if len(j_list2) > 0:
                    clones_light[g][l][str(1)] = j_list2
                if len(j_list3) > 0:
                    for c in range(0, len(j_list3)):
                        clones_light[g][l][str(c + 2)] = j_list3[c]
        clone_dict_light = {}
        cid_light = Tree()
        for g in clones_light:
            for l in clones_light[g]:
                for c in clones_light[g][l]:
                    grp_seq = clones_light[g][l][c]
                    for key, value in vj_len_lightgrp[g][l].items():
                        if value in grp_seq:
                            cid_light[g][l][c][key].value = 1
        first_key = []
        for k1 in cid_light.keys():
            first_key.append(k1)
        first_key = list(set(first_key))
        first_key_dict = dict(
            zip(first_key, range(1, len(first_key) + 1)))
        for g in cid_light:
            second_key = []
            for k2 in cid_light[g].keys():
                second_key.append(k2)
            second_key = list(set(second_key))
            second_key_dict = dict(
                zip(second_key, range(1, len(second_key) + 1)))
            for l in cid_light[g]:
                third_key = []
                for k3 in cid_light[g][l].keys():
                    third_key.append(k3)
                third_key = list(set(third_key))
                third_key_dict = dict(
                    zip(third_key, range(1, len(third_key) + 1)))
                for key, value in dict(cid_light[g][l]).items():
                    vL = []
                    for v in value:
                        if type(v) is int:
                            break
                        clone_dict_light[v] = str(
                            first_key_dict[g]) + '_' + str(second_key_dict[l]) + '_' + str(third_key_dict[key])
        lclones = list(clone_dict_light.values())
        # will just update the main dat directly
        if len(list(set(lclones))) > 1:
            lclones_dict = dict(zip(sorted(list(set(lclones))), [
                                str(x) for x in range(1, len(list(set(lclones))) + 1)]))
        else:
            lclones_dict = dict(zip(sorted(list(set(lclones))), str(1)))
        renamed_clone_dict_light = {}
        for key, value in clone_dict_light.items():
            renamed_clone_dict_light[key] = lclones_dict[value]

        cellclonetree = Tree()
        seqcellclonetree = Tree()
        for c, s, z in zip(dat['cell_id'], dat['sequence_id'], dat[clone_key]):
            seqcellclonetree[c][s].value = 1
            if pd.notnull(z):
                cellclonetree[c][z].value = 1

        for c in cellclonetree:
            cellclonetree[c] = list(cellclonetree[c])

        fintree = Tree()
        for c in tqdm(cellclonetree, desc='Refining clone assignment based on VJ chain pairing '):
            suffix = [renamed_clone_dict_light[x]
                      for x in seqcellclonetree[c] if x in renamed_clone_dict_light]
            fintree[c] = []
            if len(suffix) > 1:
                for s in suffix:
                    for cl in cellclonetree[c]:
                        fintree[c].append(cl + '_' + s)
                fintree[c] = sorted(fintree[c])
            else:
                for cl in cellclonetree[c]:
                    if len(suffix) > 0:
                        fintree[c].append(cl + '_' + ''.join(suffix))
                    else:
                        fintree[c].append(cl)
            fintree[c] = '|'.join(fintree[c])
        dat[clone_key] = [fintree[x] for x in dat['cell_id']]

    dat_[clone_key] = pd.Series(dat[clone_key])
    dat_[clone_key].replace('', 'unassigned')
    if os.path.isfile(str(self)):
        dat_.to_csv("{}/{}_clone.tsv".format(os.path.dirname(self),
                                             os.path.basename(self).split('.tsv')[0]), sep='\t', index=False)

    sleep(0.5)
    logg.info(' finished', time=start,
              deep=('Updated Dandelion object: \n'
                    '   \'data\', contig-indexed clone table\n'
                    '   \'metadata\', cell-indexed clone table\n'))
    if self.__class__ == Dandelion:
        if self.germline is not None:
            germline_ = self.germline
        else:
            germline_ = None
        if self.distance is not None:
            dist_ = self.distance
        else:
            dist_ = None
        if self.edges is not None:
            edge_ = self.edges
        else:
            edge_ = None
        if self.layout is not None:
            layout_ = self.layout
        else:
            layout_ = None
        if self.graph is not None:
            graph_ = self.graph
        else:
            graph_ = None
        if self.threshold is not None:
            threshold_ = self.threshold
        else:
            threshold_ = None
        if ('clone_id' in self.data.columns) and (key_added is None):
            # TODO: need to check the following bits if it works properly if only heavy chain tables are provided
            self.__init__(data=dat_, germline=germline_, distance=dist_,
                          edges=edge_, layout=layout_, graph=graph_, locus = locus)
            update_metadata(self, reinitialize=True, locus = locus)
        elif ('clone_id' in self.data.columns) and (key_added is not None):
            self.__init__(data=dat_, germline=germline_, distance=dist_,
                          edges=edge_, layout=layout_, graph=graph_, locus = locus)
            update_metadata(self, reinitialize=True, clone_key='clone_id',
                            retrieve=clone_key, split=False, collapse=True, combine=True, locus = locus)
        else:
            self.__init__(data=dat_, germline=germline_, distance=dist_, edges=edge_,
                          layout=layout_, graph=graph_, clone_key=clone_key, locus = locus)
            update_metadata(self, reinitialize=True, clone_key=clone_key, locus = locus)
        self.threshold = threshold_

    else:
        out = Dandelion(data=dat_, clone_key=clone_key, retrieve=clone_key,
                        split=False, collapse=True, combine=True)
        return(out)


def transfer(self: AnnData, dandelion: Dandelion, expanded_only: bool = False, neighbors_key: Union[None, str] = None, rna_key: Union[None, str] = None, vdj_key: Union[None, str] = None, overwrite: Union[None, bool, Sequence, str] = None) -> AnnData:
    """
    Transfer data in `Dandelion` slots to `AnnData` object, updating the `.obs`, `.uns`, `.obsm` and `.obsp`slots.

    Parameters
    ----------
    self : AnnData
        `AnnData` object.
    dandelion : Dandelion
        `Dandelion` object.
    expanded_only : bool
        Whether or not to transfer the embedding with all cells with BCR (False) or only for expanded clones (True).
    neighbors_key : str, optional
        key for 'neighbors' slot in `.uns`.
    rna_key : str, optional
        prefix for stashed RNA connectivities and distances.
    vdj_key : str, optional
        prefix for stashed VDJ connectivities and distances.
    overwrite : str, bool, list, optional
        Whether or not to overwrite existing anndata columns. Specifying a string indicating column name or list of column names will overwrite that specific column(s).

    Returns
    ----------
    `AnnData` object with updated `.obs`, `.obsm` and '.obsp' slots with data from `Dandelion` object.
    """
    start = logg.info('Transferring network')
    if dandelion.graph is not None:
        if expanded_only:
            G = dandelion.graph[1]
        else:
            G = dandelion.graph[0]
        distances = nx.to_pandas_adjacency(
            G, dtype=np.float32, weight='weight', nonedge=np.nan)
        connectivities = nx.to_pandas_adjacency(
            G, dtype=np.float32, weight='weight', nonedge=np.nan)
        connectivities[~connectivities.isnull()] = 1

        A = np.zeros(shape=(len(self.obs_names), len(self.obs_names)))
        B = A.copy()
        df_connectivities = pd.DataFrame(
            A, index=self.obs_names, columns=self.obs_names)
        df_distances = pd.DataFrame(
            B, index=self.obs_names, columns=self.obs_names)
        print('converting matrices')
        df_connectivities.update(connectivities)
        df_distances.update(distances)

        df_connectivities_ = csr_matrix(
            df_connectivities.values, dtype=np.float32)
        df_distances_ = csr_matrix(df_distances.values, dtype=np.float32)

        print('Updating anndata slots')
        if neighbors_key is None:
            neighbors_key = "neighbors"
        if neighbors_key not in self.uns:
            raise ValueError(
                "`edges=True` requires `pp.neighbors` to be run before.")

        rna_neighbors_key = 'rna_' + neighbors_key
        vdj_neighbors_key = 'vdj_' + neighbors_key
        if rna_neighbors_key not in self.uns:
            self.uns[rna_neighbors_key] = self.uns[neighbors_key].copy()

        if neighbors_key not in self.uns:
            raise ValueError(
                "`edges=True` requires `pp.neighbors` to be run before.")

        if rna_key is None:
            r_connectivities_key = 'rna_connectivities'
            r_distances_key = 'rna_distances'
        else:
            r_connectivities_key = rna_key + '_connectivitites'
            r_distances_key = rna_key + '_distances'

        if vdj_key is None:
            b_connectivities_key = 'vdj_connectivities'
            b_distances_key = 'vdj_distances'
        else:
            b_connectivities_key = vdj_key + '_connectivitites'
            b_distances_key = vdj_key + '_distances'

        # stash_rna_connectivities:
        if r_connectivities_key not in self.obsp:
            try:
                self.obsp[r_connectivities_key] = self.obsp["connectivities"].copy()
                self.obsp[r_distances_key] = self.obsp["distances"].copy()
            except:
                self.obsp[r_connectivities_key] = self.uns[neighbors_key]["connectivities"]
                self.obsp[r_distances_key] = self.uns[neighbors_key]["distances"]

        # always overwrite the bcr slots
        self.obsp['connectivities'] = df_connectivities_.copy()
        self.obsp['distances'] = df_distances_.copy()
        self.obsp[b_connectivities_key] = self.obsp["connectivities"].copy()
        self.obsp[b_distances_key] = self.obsp["distances"].copy()

    # always overwrite with whatever columns are in dandelion's metadata:
    for x in dandelion.metadata.columns:
        if x not in self.obs.columns:
            self.obs[x] = pd.Series(dandelion.metadata[x])
        elif overwrite is True:
            self.obs[x] = pd.Series(dandelion.metadata[x])
        if type_check(dandelion.metadata, x):
            self.obs[x].replace(np.nan, 'No_contig', inplace=True)

    if overwrite is not None and overwrite is not True:
        if not type(overwrite) is list:
            overwrite = [overwrite]
        for ow in overwrite:
            self.obs[ow] = pd.Series(dandelion.metadata[ow])
            if type_check(dandelion.metadata, ow):
                self.obs[ow].replace(np.nan, 'No_contig', inplace=True)

    tmp = self.obs.copy()
    if dandelion.layout is not None:
        if expanded_only:
            coord = pd.DataFrame.from_dict(dandelion.layout[1], orient='index')
        else:
            coord = pd.DataFrame.from_dict(dandelion.layout[0], orient='index')
        for x in coord.columns:
            tmp[x] = coord[x]
        # tmp[[1]] = tmp[[1]]*-1
        X_vdj = np.array(tmp[[0, 1]], dtype=np.float32)
        self.obsm['X_vdj'] = X_vdj

        logg.info(' finished', time=start,
                  deep=('updated `.obs` with `.metadata`\n'
                        'added to `.uns[\'' + neighbors_key + '\']` and `.obsp`\n'
                        '   \'distances\', cluster-weighted adjacency matrix\n'
                        '   \'connectivities\', cluster-weighted adjacency matrix'))
    else:
        logg.info(' finished', time=start,
                  deep=('updated `.obs` with `.metadata`\n'))


def define_clones(self: Union[Dandelion, pd.DataFrame, str], dist: Union[None, float] = None, action: Literal['first', 'set'] = 'set', model: Literal['ham', 'aa', 'hh_s1f', 'hh_s5f', 'mk_rs1nf', 'mk_rs5nf', 'hs1f_compat', 'm1n_compat'] = 'ham', norm: Literal['len', 'mut', 'none'] = 'len', doublets: Literal['drop', 'count'] = 'drop', fileformat: Literal['changeo', 'airr'] = 'airr', ncpu: Union[None, int] = None, dirs: Union[None, str] = None, outFilePrefix: Union[None, int] = None, key_added: Union[None, int] = None, verbose: bool = False)->Dandelion:
    """
    Find clones using changeo's `DefineClones.py <https://changeo.readthedocs.io/en/stable/tools/DefineClones.html>`__.

    Parameters
    ----------
    self : Dandelion, DataFrame, str
        `Dandelion` object, pandas `DataFrame` in changeo/airr format, or file path to changeo/airr file after clones have been determined.
    dist : float, optional
        The distance threshold for clonal grouping. If None, the value will be retrieved from the Dandelion class .threshold slot.
    action : str
        Specifies how to handle multiple V(D)J assignments for initial grouping. Default is 'set'. The “first” action will use only the first gene listed. The “set” action will use all gene assignments and construct a larger gene grouping composed of any sequences sharing an assignment or linked to another sequence by a common assignment (similar to single-linkage).
    model : str
        Specifies which substitution model to use for calculating distance between sequences. Default is 'ham'. The “ham” model is nucleotide Hamming distance and “aa” is amino acid Hamming distance. The “hh_s1f” and “hh_s5f” models are human specific single nucleotide and 5-mer content models, respectively, from Yaari et al, 2013. The “mk_rs1nf” and “mk_rs5nf” models are mouse specific single nucleotide and 5-mer content models, respectively, from Cui et al, 2016. The “m1n_compat” and “hs1f_compat” models are deprecated models provided backwards compatibility with the “m1n” and “hs1f” models in Change-O v0.3.3 and SHazaM v0.1.4. Both 5-mer models should be considered experimental.
    norm : str
        Specifies how to normalize distances. Default is 'len'. 'none' (do not normalize), 'len' (normalize by length), or 'mut' (normalize by number of mutations between sequences).
    doublets : str
        Option to control behaviour when dealing with heavy chain 'doublets'. Default is 'drop'. 'drop' will filter out the doublets while 'count' will retain only the highest umi count contig.
    fileformat : str
        format of V(D)J file/objects. Default is 'airr'. Also accepts 'changeo'.
    ncpu : int, optional
        number of cpus for parallelization. Default is all available cpus.
    dirs : str, optional
        If specified, out file will be in this location.
    outFilePrefix : str, optional
        If specified, the out file name will have this prefix. None defaults to 'dandelion_define_clones'
    verbose : bool
        Whether or not to print the command used in terminal to call DefineClones.py. Default is False.

    Returns
    -------
    `Dandelion` object with clone_id annotated in `.data` slot and `.metadata` initialized.
    """
    start = logg.info('Finding clones')
    if ncpu is None:
        nproc = multiprocessing.cpu_count()
    else:
        nproc = ncpu

    if key_added is None:
        clone_key = 'clone_id'
    else:
        clone_key = key_added

    if self.__class__ == Dandelion:
        dat = load_data(self.data)
    else:
        dat = load_data(self)
    if os.path.isfile(str(self)):
        dat = load_data(self)
    dat_h = dat[dat['locus'] == 'IGH']
    dat_l = dat[dat['locus'].isin(['IGK', 'IGL'])]

    if os.path.isfile(str(self)):
        tmpFolder = "{}/tmp".format(os.path.dirname(self))
        outFolder = "{}".format(os.path.dirname(self))
    else:
        import tempfile
        tmpFolder = "{}/tmp".format(tempfile.TemporaryDirectory().name)
        outFolder = "{}".format(tempfile.TemporaryDirectory().name)

    if not os.path.exists(tmpFolder):
        os.makedirs(tmpFolder)
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)

    if os.path.isfile(str(self)):
        h_file1 = "{}/{}_heavy-clone.tsv".format(
            tmpFolder, os.path.basename(self).split('.tsv')[0])
        h_file2 = "{}/{}_heavy-clone.tsv".format(
            outFolder, os.path.basename(self).split('.tsv')[0])
        l_file = "{}/{}_light.tsv".format(tmpFolder,
                                          os.path.basename(self).split('.tsv')[0])
        outfile = "{}/{}_clone.tsv".format(outFolder,
                                           os.path.basename(self).split('.tsv')[0])
    else:
        if outFilePrefix is not None:
            out_FilePrefix = outFilePrefix
        else:
            out_FilePrefix = 'dandelion_define_clones'
        h_file1 = "{}/{}_heavy-clone.tsv".format(tmpFolder, out_FilePrefix)
        h_file2 = "{}/{}_heavy-clone.tsv".format(outFolder, out_FilePrefix)
        l_file = "{}/{}_light.tsv".format(tmpFolder, out_FilePrefix)
        outfile = "{}/{}_clone.tsv".format(outFolder, out_FilePrefix)

    dat_h.to_csv(h_file1, sep='\t', index=False)
    dat_l.to_csv(l_file, sep='\t', index=False)

    if 'v_call_genotyped' in dat.columns:
        v_field = 'v_call_genotyped'
    else:
        v_field = 'v_call'

    if dist is None:
        if self.__class__ == Dandelion:
            if self.threshold is not None:
                dist_ = self.threshold
            else:
                raise ValueError(
                    'Threshold value in Dandelion object is None. Please run calculate_threshold first')
        else:
            raise ValueError(
                'Distance value is None. Please provide a distance value (float)')
    else:
        dist_ = dist

    cmd = ['DefineClones.py',
           '-d', h_file1,
           '-o', h_file2,
           '--act', action,
           '--model', model,
           '--norm', norm,
           '--dist', str(dist_),
           '--nproc', str(nproc),
           '--vf', v_field]

    def clusterLinkage(cell_series, group_series):
        """
        Returns a dictionary of {cell_id : cluster_id} that identifies clusters of cells by analyzing their shared
        features (group_series) using single linkage.

        Arguments:
        cell_series (iter): iter of cell ids.
        group_series (iter): iter of group ids.

        Returns:
        dict:  dictionary of {cell_id : cluster_id}.
        """

        # assign initial clusters
        # initial_dict = {cluster1: [cell1], cluster2: [cell1]}
        initial_dict = {}
        for cell, group in zip(cell_series, group_series):
            try:
                initial_dict[group].append(cell)
            except KeyError:
                initial_dict[group] = [cell]

        # naive single linkage clustering (ON^2 best case, ON^3 worst case) ...ie for cells with multiple light chains
        # cluster_dict = {cluster1: [cell1, cell2]}, 2 cells belong in same group if they share 1 light chain
        while True:
            cluster_dict = {}
            for i, group in enumerate(initial_dict.keys()):
                cluster_dict[i] = initial_dict[group]
                for cluster in cluster_dict:
                    # if initial_dict[group] and cluster_dict[cluster] share common cells, add initial_dict[group] to cluster
                    if cluster != i and any(cell in initial_dict[group] for cell in cluster_dict[cluster]):
                        cluster_dict[cluster] = cluster_dict[cluster] + \
                            initial_dict[group]
                        del cluster_dict[i]
                        break
            # break if clusters stop changing, otherwise restart
            if len(cluster_dict.keys()) == len(initial_dict.keys()):
                break
            else:
                initial_dict = cluster_dict.copy()

        # invert cluster_dict for return
        assign_dict = {cell: k for k, v in cluster_dict.items()
                       for cell in set(v)}

        return assign_dict

    def _lightCluster(heavy_file, light_file, out_file, doublets, fileformat):
        """
        Split heavy chain clones based on light chains

        Arguments:
        heavy_file (str): heavy chain input file.
        light_file (str): light chain input file.
        out_file (str): heavy chain output file.
        doublets (str): method for handling multiple heavy chains per cell. one of 'drop' or 'count'.
        format (str): file format. one of 'changeo' or 'airr'.
        """
        # Set column names
        if fileformat == 'changeo':
            cell_id = 'cell_id'
            clone_id = 'clone_id'
            v_call = 'v_call'
            j_call = 'j_call'
            junction_length = 'junction_length'
            umi_count = 'umicount'
        elif fileformat == 'airr':
            cell_id = 'cell_id'
            clone_id = 'clone_id'
            v_call = 'v_call'
            j_call = 'j_call'
            junction_length = 'junction_length'
            umi_count = 'umi_count'
        else:
            sys.exit("Invalid format %s" % fileformat)

        # read in heavy and light DFs
        heavy_df = pd.read_csv(heavy_file, dtype='object', na_values=[
                               '', 'None', 'NA'], sep='\t')
        light_df = pd.read_csv(light_file, dtype='object', na_values=[
                               '', 'None', 'NA'], sep='\t')

        # column checking
        expected_heavy_columns = [cell_id, clone_id,
                                  v_call, j_call, junction_length, umi_count]
        if set(expected_heavy_columns).issubset(heavy_df.columns) is False:
            raise ValueError(
                "Missing one or more columns in heavy chain file: " + ", ".join(expected_heavy_columns))
        expected_light_columns = [cell_id, v_call,
                                  j_call, junction_length, umi_count]
        if set(expected_light_columns).issubset(light_df.columns) is False:
            raise ValueError(
                "Missing one or more columns in light chain file: " + ", ".join(expected_light_columns))

        # Fix types
        heavy_df[junction_length] = heavy_df[junction_length].astype('int')
        light_df[junction_length] = light_df[junction_length].astype('int')

        # filter multiple heavy chains
        if doublets == 'drop':
            heavy_df = heavy_df.drop_duplicates(cell_id, keep=False)
            if heavy_df.empty is True:
                raise ValueError(
                    "Empty heavy chain data, after doublets drop. Are you combining experiments in a single file? If so, split your data into multiple files.")
        elif doublets == 'count':
            heavy_df[umi_count] = heavy_df[umi_count].astype('int')
            heavy_df = heavy_df.groupby(cell_id, sort=False).apply(
                lambda x: x.nlargest(1, umi_count))

        # transfer clone IDs from heavy chain df to light chain df
        clone_dict = {v[cell_id]: v[clone_id]
                      for k, v in heavy_df[[clone_id, cell_id]].T.to_dict().items()}
        light_df = light_df.loc[light_df[cell_id].apply(
            lambda x: x in clone_dict.keys()), ]
        light_df[clone_id] = light_df.apply(
            lambda row: clone_dict[row[cell_id]], axis=1)

        # generate a "cluster_dict" of CELL:CLONE dictionary from light df  (TODO: use receptor object V/J gene names)
        cluster_dict = clusterLinkage(light_df[cell_id],
                                      light_df.apply(lambda row:
                                                     getGene(row[v_call]) + ',' +
                                                     getGene(row[j_call]) + ',' +
                                                     str(row[junction_length]) + ',' + row[clone_id], axis=1))

        # add assignments to heavy_df
        heavy_df = heavy_df.loc[heavy_df[cell_id].apply(
            lambda x: x in cluster_dict.keys()), :]
        heavy_df[clone_id] = heavy_df[clone_id] + '_' + \
            heavy_df.apply(lambda row: str(cluster_dict[row[cell_id]]), axis=1)

        # write heavy chains
        heavy_df.to_csv(out_file, sep='\t', index=False)
        return(heavy_df, light_df)

    if verbose:
        print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd)

    h_df, l_df = _lightCluster(
        h_file2, l_file, outfile, doublets=doublets, fileformat=fileformat)

    h_df = load_data(h_df)
    # create a dictionary for cell_id : clone_id from h_df
    linked_clones = dict(zip(h_df['cell_id'], h_df['clone_id']))

    # create a clone_reference
    clone_ref = list(set(h_df['clone_id']))
    clone_ref = [c.split('_')[1] if c is not np.nan else c for c in clone_ref]
    l_df = load_data(l_df)

    for x in l_df.index:
        if l_df.loc[x, 'clone_id'] in clone_ref:
            l_df.at[x, 'clone_id'] = linked_clones[l_df.loc[x, 'cell_id']]
        else:
            try:
                l_df.at[x, 'clone_id'] = l_df.loc[x, 'cell_id'] + '_notlinked'
            except:
                pass

    cloned_ = pd.concat([h_df, l_df])
    # transfer the new clone_id to the heavy + light file
    dat[str(clone_key)] = pd.Series(cloned_['clone_id'])

    if self.__class__ == Dandelion:
        if self.germline is not None:
            germline_ = self.germline
        else:
            germline_ = None
        if self.distance is not None:
            dist_ = self.distance
        else:
            dist_ = None
        if self.edges is not None:
            edge_ = self.edges
        else:
            edge_ = None
        if self.layout is not None:
            layout_ = self.layout
        else:
            layout_ = None
        if self.graph is not None:
            graph_ = self.graph
        else:
            graph_ = None
        if self.threshold is not None:
            threshold_ = self.threshold
        else:
            threshold_ = None

        if ('clone_id' in self.data.columns) and (clone_key is not None):
            self.__init__(data=dat, germline=germline_, distance=dist_, edges=edge_, layout=layout_,
                          graph=graph_, initialize=True, retrieve=clone_key, split=False, collapse=True, combine=True)
        elif ('clone_id' not in self.data.columns) and (clone_key is not None):
            self.__init__(data=dat, germline=germline_, distance=dist_, edges=edge_, layout=layout_, graph=graph_,
                          initialize=True, clone_key=clone_key, retrieve=clone_key, split=False, collapse=True, combine=True)
        else:
            self.__init__(data=dat, germline=germline_, distance=dist_, edges=edge_,
                          layout=layout_, graph=graph_, initialize=True, clone_key=clone_key)
        self.threshold = threshold_
    else:
        if ('clone_id' in dat.columns) and (clone_key is not None):
            out = Dandelion(data=dat, retrieve=clonekey, split=False)
        elif ('clone_id' not in dat.columns) and (clone_key is not None):
            out = Dandelion(data=dat, clone_key=clone_key)
        else:
            out = Dandelion(data=dat)
        return(out)
    logg.info(' finished', time=start,
              deep=('Updated Dandelion object: \n'
                    '   \'data\', contig-indexed clone table\n'
                    '   \'metadata\', cell-indexed clone table\n'))


def clone_size(self: Dandelion, max_size: Union[None, int] = None, clone_key: Union[None, str] = None, key_added: Union[None, str] = None):
    """
    Quantifies size of clones

    Parameters
    ----------
    self : Dandelion
        `Dandelion` object
    max_size : int, optional
        The maximum size before value gets clipped. If None, the value will be returned as a numerical value.
    clone_key : str, optional
        Column name specifying the clone_id column in metadata.
    key_added : str, optional
        column name where clone size is tabulated into.

    Returns
    -------
    `Dandelion` object with clone size columns annotated in `.metadata` slot.
    """

    start = logg.info('Quantifying clone sizes')

    metadata_ = self.metadata.copy()

    if clone_key is None:
        clonekey = 'clone_id'
    else:
        clonekey = clone_key

    tmp = metadata_[str(clonekey)].str.split('|', expand=True).stack()
    tmp = tmp.reset_index(drop=False)
    tmp.columns = ['cell_id', 'tmp', str(clonekey)]

    clonesize = tmp[str(clonekey)].value_counts()

    if max_size is not None:
        clonesize_ = clonesize.astype('object')
        for i in clonesize.index:
            if clonesize.loc[i] >= max_size:
                clonesize_.at[i] = '>= ' + str(max_size)
        clonesize_ = clonesize_.astype('category')
    else:
        clonesize_ = clonesize.copy()

    clonesize_dict = dict(clonesize_)

    if max_size is not None:
        if key_added is None:
            self.metadata[str(clonekey) + '_size_max_' + str(max_size)] = pd.Series(dict(zip(metadata_.index, [str(y) for y in [sorted(list(set([clonesize_dict[c_] for c_ in c.split('|')])),
                                                                                                                                       key=lambda x: int(x.split('>= ')[1]) if type(x) is str else int(x), reverse=True)[0] if '|' in c else clonesize_dict[c] for c in metadata_[str(clonekey)]]])))
            self.metadata[str(clonekey) + '_size_max_' + str(max_size)] = self.metadata[str(
                clonekey) + '_size_max_' + str(max_size)].astype('category')
        else:
            self.metadata[key_added] = pd.Series(dict(zip(metadata_.index, [str(y) for y in [sorted(list(set([clonesize_dict[c_] for c_ in c.split('|')])), key=lambda x: int(
                x.split('>= ')[1]) if type(x) is str else int(x), reverse=True)[0] if '|' in c else clonesize_dict[c] for c in metadata_[str(clonekey)]]])))
            self.metadata[str(clonekey) + '_size_max_' + str(max_size)] = self.metadata[str(
                clonekey) + '_size_max_' + str(max_size)].astype('category')
    else:
        if key_added is None:
            self.metadata[str(clonekey) + '_size'] = pd.Series(dict(zip(metadata_.index, [str(y) for y in [sorted(list(set([clonesize_dict[c_] for c_ in c.split('|')])),
                                                                                                                  key=lambda x: int(x.split('>= ')[1]) if type(x) is str else int(x), reverse=True)[0] if '|' in c else clonesize_dict[c] for c in metadata_[str(clonekey)]]])))
            try:
                self.metadata[str(clonekey) + '_size'] = [int(x)
                                                          for x in self.metadata[str(clonekey) + '_size']]
            except:
                pass
        else:
            self.metadata[key_added] = pd.Series(dict(zip(metadata_.index, [str(y) for y in [sorted(list(set([clonesize_dict[c_] for c_ in c.split('|')])), key=lambda x: int(
                x.split('>= ')[1]) if type(x) is str else int(x), reverse=True)[0] if '|' in c else clonesize_dict[c] for c in metadata_[str(clonekey)]]])))
            try:
                self.metadata[key_added] = [
                    int(x) for x in self.metadata[str(clonekey) + '_size']]
            except:
                pass
    logg.info(' finished', time=start,
              deep=('Updated Dandelion object: \n'
                    '   \'metadata\', cell-indexed clone table'))


def clone_overlap(self: Union[Dandelion, AnnData], groupby: str, colorby: str, min_clone_size: Union[None, int] = None, clone_key: Union[None, str] = None) -> Union[AnnData, pd.DataFrame]:
    """
    A function to tabulate clonal overlap for input as a circos-style plot.

    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    groupby : str
        column name in obs/metadata for collapsing to nodes in circos plot.
    colorby : str
        column name in obs/metadata for grouping and color of nodes in circos plot.
    min_clone_size : int, optional
        minimum size of clone for plotting connections. Defaults to 2 if left as None.
    clone_key : str, optional
        column name for clones. None defaults to 'clone_id'.

    Returns
    -------
    a `pandas DataFrame`.
    """
    start = logg.info('Finding clones')
    if self.__class__ == Dandelion:
        data = self.metadata.copy()
    elif self.__class__ == AnnData:
        data = self.obs.copy()

    if min_clone_size is None:
        min_size = 2
    else:
        min_size = int(min_clone_size)

    if clone_key is None:
        clone_ = 'clone_id'
    else:
        clone_ = clone_key

    # get rid of problematic rows that appear because of category conversion?
    data = data[~(data[clone_].isin(
        [np.nan, 'nan', 'NaN', 'No_contig', 'unassigned', None]))]

    # prepare a summary table
    datc_ = data[clone_].str.split('|', expand=True).stack()
    datc_ = pd.DataFrame(datc_)
    datc_.reset_index(drop=False, inplace=True)
    datc_.columns = ['cell_id', 'tmp', clone_]
    datc_.drop('tmp', inplace=True, axis=1)
    datc_ = datc_[~(datc_[clone_].isin(
        ['', np.nan, 'nan', 'NaN', 'No_contig', 'unassigned', None]))]
    dictg_ = dict(data[groupby])
    datc_[groupby] = [dictg_[l] for l in datc_['cell_id']]

    overlap = pd.crosstab(datc_[clone_], datc_[groupby])

    if min_size == 0:
        raise ValueError('min_size must be greater than 0.')
    elif min_size > 2:
        overlap[overlap < min_size] = 0
        overlap[overlap >= min_size] = 1
    elif min_size == 2:
        overlap[overlap >= min_size] = 1

    overlap.index.name = None
    overlap.columns.name = None

    if self.__class__ == AnnData:
        self.uns['clone_overlap'] = overlap.copy()
        logg.info(' finished', time=start,
                  deep=('Updated AnnData: \n'
                        '   \'uns\', clone overlap table'))
    else:
        return(overlap)
