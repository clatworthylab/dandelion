#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2021-02-11 12:22:40
# @Last Modified by:   Kelvin
# @Last Modified time: 2021-06-17 11:25:01

import os
from collections import defaultdict
import pandas as pd
import numpy as np
import re
import copy
from changeo.IO import readGermlines
import warnings
import h5py
import networkx as nx
import bz2
import gzip
from anndata import AnnData
import _pickle as cPickle
try:
    from scanpy import logging as logg
except ImportError:
    pass
from ..utilities._utilities import *
from ..utilities._io import *
from typing import Union, Sequence, Tuple, Dict


class Dandelion:
    """
    `Dandelion` class object.

    Main class object storing input/ouput slots for all functions.

    """

    def __init__(self, data=None, metadata=None, germline=None, distance=None, edges=None, layout=None, graph=None, initialize=True, **kwargs):
        self.data = data
        self.metadata = metadata
        self.distance = distance
        self.edges = edges
        self.layout = layout
        self.graph = graph
        self.threshold = None
        self.germline = {}

        if germline is not None:
            self.germline.update(germline)

        if os.path.isfile(str(self.data)):
            self.data = load_data(self.data)

        if self.data is not None:
            self.n_contigs = self.data.shape[0]
            if metadata is None:
                if initialize is True:
                    try:
                        update_metadata(self, **kwargs)
                    except:
                        update_metadata(self, locus = 'tr-ab', **kwargs)
                        try:
                            update_metadata(self, locus = 'tr-gd', **kwargs)
                        except:
                            raise OSError('Unable to read data. Are you sure this is an AIRR formatted data?')
                try:
                    self.n_obs = self.metadata.shape[0]
                except:
                    self.n_obs = 0
            else:
                self.metadata = metadata
                self.n_obs = self.metadata.shape[0]
        else:
            self.n_contigs = 0
            self.n_obs = 0

    def _gen_repr(self, n_obs, n_contigs) -> str:
        # inspire by AnnData's function
        descr = f"Dandelion class object with n_obs = {n_obs} and n_contigs = {n_contigs}"
        for attr in ["data", "metadata", "distance", "edges"]:
            try:
                keys = getattr(self, attr).keys()
            except:
                keys = []
            if len(keys) > 0:
                descr += f"\n    {attr}: {str(list(keys))[1:-1]}"
            else:
                descr += f"\n    {attr}: {str(None)}"
        if self.layout is not None:
            descr += f"\n    layout: {', '.join(['layout for '+ str(len(x)) + ' vertices' for x in (self.layout[0], self.layout[1])])}"
        else:
            descr += f"\n    layout: {str(None)}"
        if self.graph is not None:
            descr += f"\n    graph: {', '.join(['networkx graph of '+ str(len(x)) + ' vertices' for x in (self.graph[0], self.graph[1])])} "
        else:
            descr += f"\n    graph: {str(None)}"
        return descr

    def __repr__(self) -> str:
        # inspire by AnnData's function
        return self._gen_repr(self.n_obs, self.n_contigs)

    def copy(self):
        """
        Performs a deep copy of all slots in `Dandelion` class.

        Parameters
        ----------
        self : Dandelion
            `Dandelion` object.

        Returns
        -------
        a deep copy of `Dandelion` class.
        """
        return copy.deepcopy(self)

    def update_germline(self, corrected: Union[None, Dict, str] = None, germline: Union[None, str] = None, org: Literal['human', 'mouse'] = 'human'):
        """
        Update germline reference with corrected sequences and store in `Dandelion` object.

        Parameters
        ----------
        self : Dandelion
            `Dandelion` object.
        corrected : dict, str, optional
            dictionary of corrected germline sequences or file path to corrected germline sequences fasta file.
        germline : str, optional
            path to germline database folder. Defaults to `$GERMLINE` environmental variable.
        org : str
            organism of reference folder. Default is 'human'.

        Returns
        -------
        updated germline reference diciontary in `.germline` slot.
        """
        start = logg.info('Updating germline reference')
        env = os.environ.copy()
        if germline is None:
            try:
                gml = env['GERMLINE']
            except:
                raise OSError(
                    'Environmental variable GERMLINE must be set. Otherwise, please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files.')
            gml = gml + 'imgt/' + org + '/vdj/'
        else:
            if os.path.isdir(germline):
                gml = germline
            elif type(germline) is not list:
                germline_ = [germline]
                if len(germline_) < 3:
                    raise OSError('Input for germline is incorrect. Please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files, or individual paths to the germline IGHV, IGHD, and IGHJ fasta files (with .fasta extension) as a list.')
                else:
                    gml = []
                    for x in germline_:
                        if not x.endswith('.fasta'):
                            raise OSError(
                                'Input for germline is incorrect. Please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files, or individual paths to the germline IGHV, IGHD, and IGHJ fasta files (with .fasta extension) as a list.')
                        gml.append(x)
            elif type(germline) is list:
                if len(germline) < 3:
                    raise OSError('Input for germline is incorrect. Please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files, or individual paths to the germline IGHV, IGHD, and IGHJ fasta files (with .fasta extension) as a list.')
                else:
                    gml = []
                    for x in germline:
                        if not x.endswith('.fasta'):
                            raise OSError(
                                'Input for germline is incorrect. Please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files, or individual paths to the germline IGHV, IGHD, and IGHJ fasta files (with .fasta extension) as a list.')
                        gml.append(x)

        if type(gml) is not list:
            gml = [gml]

        germline_ref = readGermlines(gml)
        if corrected is not None:
            if type(corrected) is dict:
                personalized_ref_dict = corrected
            elif os.path.isfile(str(corrected)):
                personalized_ref_dict = readGermlines([corrected])
            # update with the personalized germline database
            if 'personalized_ref_dict' in locals():
                germline_ref.update(personalized_ref_dict)
            else:
                raise OSError(
                    'Input for corrected germline fasta is incorrect. Please provide path to file containing corrected germline fasta sequences.')

        self.germline.update(germline_ref)
        logg.info(' finished', time=start,
                  deep=('Updated Dandelion object: \n'
                        '   \'germline\', updated germline reference\n'))

    def write_pkl(self, filename: str = 'dandelion_data.pkl.pbz2', **kwargs):
        """
        Writes a `Dandelion` class to .pkl format.

        Parameters
        ----------
        filename
            path to `.pkl` file.
        **kwargs
            passed to `_pickle`.
        """
        if isBZIP(filename):
            try:
                with bz2.BZ2File(filename, 'wb') as f:
                    cPickle.dump(self, f, **kwargs)
            except:
                with bz2.BZ2File(filename, 'wb') as f:
                    cPickle.dump(self, f, protocol=4, **kwargs)
        elif isGZIP(filename):
            try:
                with gzip.open(filename, 'wb') as f:
                    cPickle.dump(self, f, **kwargs)
            except:
                with gzip.open(filename, 'wb') as f:
                    cPickle.dump(self, f, protocol=4, **kwargs)
        else:
            f = open(filename, 'wb')
            cPickle.dump(self, f, **kwargs)
            f.close()

    def write_h5(self, filename: str = 'dandelion_data.h5', complib: Literal['zlib', 'lzo', 'bzip2', 'blosc', 'blosc:blosclz', 'blosc:lz4', 'blosc:lz4hc', 'blosc:snappy', 'blosc:zlib', 'blosc:zstd'] = None, compression: Literal['zlib', 'lzo', 'bzip2', 'blosc', 'blosc:blosclz', 'blosc:lz4', 'blosc:lz4hc', 'blosc:snappy', 'blosc:zlib', 'blosc:zstd'] = None, compression_level: Union[None, int] = None, **kwargs):
        """
        Writes a `Dandelion` class to .h5 format.

        Parameters
        ----------
        filename
            path to `.h5` file.
        complib : str, optional
            method for compression for data frames. see (https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_hdf.html) for more options.
        compression : str, optional
            same call as complib. Just a convenience option.
        compression_opts : {0-9}, optional
            Specifies a compression level for data. A value of 0 disables compression.
        **kwargs
            passed to `pd.DataFrame.to_hdf`.
        """
        if compression_level is None:
            compression_level = 9
        else:
            compression_level = compression_level

        # a little hack to overwrite the existing file?
        with h5py.File(filename, "w") as hf:
            for datasetname in hf.keys():
                del hf[datasetname]

        if complib is None and compression is None:
            comp = None
        elif complib is not None and compression is None:
            comp = complib
        elif complib is None and compression is not None:
            comp = compression
        if complib is not None and compression is not None:
            raise ValueError(
                'Please specify only complib or compression. They do the same thing.')

        # now to actually saving
        for col in self.data.columns:
            weird = (self.data[[col]].applymap(type) !=
                     self.data[[col]].iloc[0].apply(type)).any(axis=1)
            if len(self.data[weird]) > 0:
                self.data[col] = self.data[col].where(
                    pd.notnull(self.data[col]), '')
        self.data.to_hdf(filename, "data", complib=comp,
                         complevel=compression_level, **kwargs)

        if self.metadata is not None:
            for col in self.metadata.columns:
                weird = (self.metadata[[col]].applymap(
                    type) != self.metadata[[col]].iloc[0].apply(type)).any(axis=1)
                if len(self.metadata[weird]) > 0:
                    self.metadata[col] = self.metadata[col].where(
                        pd.notnull(self.metadata[col]), '')
            self.metadata.to_hdf(filename, "metadata", complib=comp,
                                 complevel=compression_level, format='table', nan_rep=np.nan, **kwargs)
        # except:
        #     warnings.warn("`metadata` slot not saved. Please check if there is incompatible dtypes in the metadata table.")
        #     pass
        try:
            if 'index' in self.edges.columns:
                self.edges.drop('index', axis=1, inplace=True)
            self.edges.to_hdf(filename, "edges", complib=comp,
                              complevel=compression_level, **kwargs)
        except:
            pass

        graph_counter = 0
        try:
            for g in self.graph:
                G = nx.to_pandas_adjacency(g, nonedge=np.nan)
                G.to_hdf(filename, "graph/graph_" + str(graph_counter),
                         complib=comp, complevel=compression_level, **kwargs)
                graph_counter += 1
        except:
            pass

        try:
            for d in self.distance:
                # how to make this faster?
                dat = pd.DataFrame(self.distance[d].toarray())
                dat.to_hdf(filename, "distance/" + d, complib=comp,
                           complevel=compression_level, **kwargs)
        except:
            pass

        with h5py.File(filename, "a") as hf:
            # try:
            #     for d in self.distance:
            #         hf.create_dataset('distance/'+d+'/data', data=self.distance[d].data, compression = compression, compression_opts = compression_level)
            #         hf.create_dataset('distance/'+d+'/indptr', data=self.distance[d].indptr, compression = compression, compression_opts = compression_level)
            #         hf.create_dataset('distance/'+d+'/indices', data=self.distance[d].indices, compression = compression, compression_opts = compression_level)
            #         hf['distance/'+d].attrs['shape'] = self.distance[d].shape
            # except:
            #     pass

            try:
                layout_counter = 0
                for l in self.layout:
                    try:
                        hf.create_group('layout/layout_' + str(layout_counter))
                    except:
                        pass
                    for k in l.keys():
                        hf['layout/layout_' +
                            str(layout_counter)].attrs[k] = l[k]
                    layout_counter += 1
            except:
                pass

            if len(self.germline) > 0:
                try:
                    hf.create_group('germline')
                except:
                    pass
                for k in self.germline.keys():
                    hf['germline'].attrs[k] = self.germline[k]
            if self.threshold is not None:
                tr = self.threshold
                hf.create_dataset('threshold', data=tr)


def retrieve_metadata(data: pd.DataFrame, query: str, split: bool = True, collapse: bool = True, combine: bool = False, locus: Literal['ig', 'tr-ab', 'tr-gd'] = 'ig', split_locus: bool = False, verbose: bool = False) -> pd.DataFrame:
    data_tmp = data.copy()
    dat_dict = defaultdict(dict)
    dict_ = defaultdict(dict)
    metadata_ = defaultdict(dict)

    locus_dict1 = {'ig': 'IGH', 'tr-ab': 'TRB', 'tr-gd': 'TRD'}
    locus_dict2 = {'ig': ['IGK', 'IGL'], 'tr-ab': ['TRA'], 'tr-gd': ['TRG']}
    locus_dict3 = {'ig': 'H', 'tr-ab': 'B', 'tr-gd': 'D'}
    locus_dict4 = {'ig': 'L', 'tr-ab': 'A', 'tr-gd': 'G'}
    query_dict = dict(zip(data_tmp['sequence_id'], data_tmp[query]))

    if type_check(data, query):
        data_tmp[query].fillna('unassigned', inplace=True)

    typesoflocus = len(list(set(data_tmp['locus'])))

    if typesoflocus > 1:
        if split_locus:
            for loci in flatten([locus_dict1[locus]] + locus_dict2[locus]):
                tmp = data_tmp[data_tmp['locus'].isin([loci])].copy()
                if tmp.shape[0] > 0:
                    dat_dict[loci] = tmp.copy()
        else:
            tmp3 = data_tmp[data_tmp['locus'].isin(
                [locus_dict1[locus]])].copy()
            tmp4 = data_tmp[data_tmp['locus'].isin(locus_dict2[locus])].copy()
            if tmp3.shape[0] > 0:
                dat_dict[locus_dict3[locus]] = tmp3.copy()
            if tmp4.shape[0] > 0:
                dat_dict[locus_dict4[locus]] = tmp4.copy()
    else:
        if verbose:
            warnings.warn(UserWarning(
                'Single locus type detected. Ignoring split = True and split_locus = True.'))
        dat_dict[locus_dict3[locus]] = data_tmp[data_tmp['locus'].isin(
            [locus_dict1[locus]])].copy()

    for d in dat_dict:
        tmp = Tree()
        for cell, seq in zip(dat_dict[d]['cell_id'], dat_dict[d]['sequence_id']):
            tmp[cell][seq].value = 1
        cells = []
        seqs = []
        for t in tmp:
            cells.append(t)
            seqs.append([s for s in tmp[t]])
        metadata_[d] = pd.DataFrame(seqs, index=cells)
        metadata_[d][pd.isnull(metadata_[d])] = np.nan

    if split_locus:
        H = locus_dict1[locus]
    else:
        H = locus_dict3[locus]

    if len(metadata_[H].columns) > 1:
        metadata_[H].columns = ['sequence_id_' + H + '_' + str(x) for x in range(0, len(metadata_[H].columns))]
    else:
        metadata_[H].columns = ['sequence_id_' + H + '_0']

    if typesoflocus > 1:
        if split_locus:
            for L in locus_dict2[locus]:
                if len(metadata_[L].columns) > 1:
                    metadata_[L].columns = ['sequence_id_' + L + '_' + str(x) for x in range(0, len(metadata_[L].columns))]
                else:
                    metadata_[L].columns = ['sequence_id_' + L + '_0']
        else:
            L = locus_dict4[locus]
            if len(metadata_[L].columns) > 1:
                metadata_[L].columns = ['sequence_id_' + L + '_' + str(x) for x in range(0, len(metadata_[L].columns))]
            else:
                metadata_[L].columns = ['sequence_id_' + L + '_0']

    metadata_result = metadata_.copy()
    for l in metadata_:
        for x in metadata_[l]:
            metadata_result[l][x] = [query_dict[i] if i ==
                                     i else np.nan for i in metadata_[l][x]]

    if typesoflocus > 1:
        if not split_locus:
            results = retrieve_result_dict(
                query, data_tmp, metadata_result[locus_dict3[locus]], metadata_result[locus_dict4[locus]], locus, split, collapse, combine)
        else:
            results = retrieve_result_dict(query, data_tmp, metadata_result[locus_dict1[locus]], [
                                           metadata_result[L] for L in locus_dict2[locus]], locus, split, collapse, combine)
    else:
        results = retrieve_result_dict_singular(
            query, data_tmp, metadata_result[locus_dict3[locus]], locus, collapse, combine)
    return(results)


def retrieve_result_dict(query: str, data: pd.DataFrame, meta_h: pd.DataFrame, meta_l: pd.DataFrame, locus: Literal['ig', 'tr-ab', 'tr-gd'] = 'ig', split: bool = True, collapse: bool = True, combine: bool = False, verbose: bool = False) -> Dict:
    df_hl = defaultdict(dict)
    locus_dict1 = {'ig': 'IGH', 'tr-ab': 'TRB', 'tr-gd': 'TRD'}
    locus_dict2 = {'ig': ['IGK', 'IGL'], 'tr-ab': ['TRA'], 'tr-gd': ['TRG']}
    locus_dict3 = {'ig': 'H', 'tr-ab': 'B', 'tr-gd': 'D'}
    locus_dict4 = {'ig': 'L', 'tr-ab': 'A', 'tr-gd': 'G'}
    if len(meta_l) == 2 and type(meta_l) is list:
        H = locus_dict1[locus]
    else:
        H = locus_dict3[locus]
    if meta_h.shape[1] > 1:
        if collapse:
            newmeta_h = meta_h.copy()
            if type_check(meta_h, 'sequence_id_' + H + '_0'):
                newh = []
                for i in meta_h.index:
                    try:
                        newh.append(
                            '|'.join([h for h in list(dict.fromkeys(newmeta_h.loc[i])) if h == h]))
                    except:
                        newh.append('|'.join([str(h) for h in list(
                            dict.fromkeys(newmeta_h.loc[i])) if h == h]))
                newmeta_h['sequence_id_' + H + '_0'] = newh
                meta_h = pd.DataFrame(newmeta_h['sequence_id_' + H + '_0'].copy())
            else:
                collapse = False
                if verbose:
                    warnings.warn(UserWarning(
                        'Multiple VDJ contigs mapping to the same cell barcode and/or query dtype is {}. Ignoring collapse = True.'.format(meta_h['sequence_id_' + H + '_0'].dtype.name)))
    if len(meta_l) == 2 and type(meta_l) is list:
        metadata_ = meta_h.join(meta_l[0]).join(meta_l[1])
    else:
        metadata_ = meta_h.join(meta_l)
    df_ = metadata_.copy()
    if type_check(meta_h, 'sequence_id_' + H + '_0'):
        if split:
            df_hl[H] = df_[list(meta_h.columns)].copy()
            if len(meta_l) == 2 and type(meta_l) is list:
                for x in range(0, len(locus_dict2[locus])):
                    L = locus_dict2[locus][x]
                    df_hl[L] = df_[list(meta_l[x].columns)].copy()
            else:
                L = locus_dict4[locus]
                df_hl[L] = df_[list(meta_l.columns)].copy()
            df_res_hl = df_hl.copy()
            res_ = defaultdict(list)
            result_dict_hl = defaultdict(dict)
            if collapse:
                for d in df_res_hl:
                    for i in df_hl[d].index:
                        if combine:
                            try:
                                res_[d].append(
                                    '|'.join([l for l in list(dict.fromkeys(df_hl[d].loc[i])) if l == l]))
                            except:
                                res_[d].append('|'.join([str(l) for l in list(
                                    dict.fromkeys(df_hl[d].loc[i])) if l == l]))
                        else:
                            try:
                                res_[d].append(
                                    '|'.join([l for l in list(df_hl[d].loc[i]) if l == l]))
                            except:
                                res_[d].append(
                                    '|'.join([str(l) for l in list(df_hl[d].loc[i]) if l == l]))
                    df_res_hl[d][query] = res_[d]
                    result_dict_hl[d] = dict(df_res_hl[d][query])
                    for k, v in result_dict_hl[d].items():
                        if type(v) is not list:
                            result_dict_hl[d][k] = [v]
                if len(meta_l) == 2 and type(meta_l) is list and type(meta_l) is list:
                    result_dict_ = {query + '_' + H: result_dict_hl[H]}
                    for x in range(0, len(locus_dict2[locus])):
                        L = locus_dict2[locus][x]
                        result_dict_.update({query + '_' + L: result_dict_hl[L]})
                else:
                    result_dict_ = {query + '_VDJ': result_dict_hl[H], query + '_VJ': result_dict_hl[L]}
            else:
                for d in df_res_hl:
                    result_dict_hl[d] = df_res_hl[d]
        else:
            df_res = df_.copy()
            q_res = []
            if collapse:
                for i in metadata_.index:
                    if combine:
                        q_res.append(
                            '|'.join([qq for qq in list(dict.fromkeys(df_.loc[i])) if qq == qq]))
                    else:
                        q_res.append(
                            '|'.join([qq for qq in list(df_.loc[i]) if qq == qq]))
            else:
                for i in metadata_.index:
                    q_res.append([qq for qq in list(df_.loc[i]) if qq == qq])
            df_res[query] = q_res
            result_dict_ = dict(df_res[query])
    else:
        result_dict_ = {x: dict(df_[x]) for x in df_}
    rs_dict1 = {'ig': '[HL]', 'tr-ab': '[AB]', 'tr-gd': '[GD]'}
    rs_dict2 = {'ig': 'IG[HKL]', 'tr-ab': 'TR[AB]', 'tr-gd': 'TR[GD]'}
    if len(meta_l) == 2 and type(meta_l) is list:
        rs = '_sequence_id_' + rs_dict2[locus]
    else:
        rs = '_sequence_id_' + rs_dict1[locus]
    final_result_ = defaultdict(dict)
    if split:
        if not collapse:
            if len(meta_l) == 2 and type(meta_l) is list:
                names = [k for k in result_dict_hl.keys()]
                final_result = result_dict_hl[names[0]].join(
                    result_dict_hl[names[1]]).join(result_dict_hl[names[2]])
                final_result.columns = [re.sub('_sequence_id', '', q) for q in [
                    query + '_' + str(l) for l in final_result.columns]]
            else:
                try:
                    if type_check(meta_h, 'sequence_id_' + H + '_0'):
                        final_result_h = pd.DataFrame(result_dict_hl[H])
                    else:
                        final_result_h = pd.DataFrame(result_dict_[H])
                except:
                    if type_check(meta_h, 'sequence_id_' + H + '_0'):
                        final_result_h = pd.DataFrame.from_dict(
                            result_dict_hl, orient='index').T
                        final_result_h = final_result_h[meta_h.columns].copy()
                    else:
                        final_result_h = pd.DataFrame.from_dict(
                            result_dict_, orient='index').T
                        final_result_h = final_result_h[meta_h.columns].copy()
                final_result_h.columns = [re.sub(rs, '', q) for q in [
                    query + '_VDJ_' + str(l) for l in final_result_h.columns]]
                try:
                    if type_check(meta_h, 'sequence_id_' + H + '_0'):
                        final_result_l = pd.DataFrame(
                            result_dict_hl[locus_dict4[locus]])
                    else:
                        final_result_l = pd.DataFrame(
                            result_dict_[locus_dict4[locus]])
                except:
                    if type_check(meta_h, 'sequence_id_' + H + '_0'):
                        final_result_l = pd.DataFrame.from_dict(
                            result_dict_hl, orient='index').T
                        final_result_l = final_result_l[meta_l.columns].copy()
                    else:
                        final_result_l = pd.DataFrame.from_dict(
                            result_dict_, orient='index').T
                        final_result_l = final_result_l[meta_l.columns].copy()
                final_result_l.columns = [re.sub(rs, '', q) for q in [
                    query + '_VJ_' + str(l) for l in final_result_l.columns]]
                final_result = final_result_h.join(final_result_l)
        else:
            if len(meta_l) == 2 and type(meta_l) is list:
                if type_check(meta_h, 'sequence_id_' + H + '_0'):
                    for d in result_dict_:
                        final_result_[d] = pd.DataFrame.from_dict(
                            result_dict_[d], orient='index')
                        final_result_[d].columns = [re.sub(rs, '', q) for q in [
                            d + '_' + str(l) for l in final_result_[d].columns]]
                        final_result_[d].columns = [
                            re.sub('_[0-9]', '', q) for q in final_result_[d].columns]
                    names = [k for k in final_result_.keys()]
                    final_result = final_result_[names[0]].join(
                        final_result_[names[1]]).join(final_result_[names[2]])
                else:
                    if verbose:
                        warnings.warn(UserWarning('Query dtype is {}. Ignoring collapse = True.'.format(
                            meta_h['sequence_id_' + H + '_0'].dtype.name)))
                    final_result = pd.DataFrame.from_dict(
                        result_dict_, orient='index').T
                    final_result.columns = [
                        query + re.sub('sequence_id', '', q) for q in final_result.columns]
            else:
                if type_check(meta_h, 'sequence_id_' + H + '_0'):
                    final_result_h = pd.DataFrame.from_dict(result_dict_[query + '_VDJ'], orient='index')
                    final_result_h.columns = [query + '_VDJ']
                    final_result_l = pd.DataFrame.from_dict(result_dict_[query + '_VJ'], orient='index')
                    final_result_l.columns = [query + '_VJ']
                    final_result = final_result_h.join(final_result_l)
                else:
                    if verbose:
                        warnings.warn(UserWarning('Query dtype is {}. Ignoring collapse = True.'.format(
                            meta_h['sequence_id_' + H + '_0'].dtype.name)))
                    final_result_h = pd.DataFrame.from_dict(
                        result_dict_, orient='index').T
                    final_result_h = final_result_h[meta_h.columns].copy()
                    final_result_h.columns = [re.sub(rs, '', q) for q in [
                        query + '_VDJ_' + str(l) for l in final_result_h.columns]]
                    final_result_l = pd.DataFrame.from_dict(
                        result_dict_, orient='index').T
                    final_result_l = final_result_l[meta_l.columns].copy()
                    final_result_l.columns = [re.sub(rs, '', q) for q in [
                        query + '_VJ_' + str(l) for l in final_result_l.columns]]
                    final_result = final_result_h.join(final_result_l)
    else:
        if type_check(meta_h, 'sequence_id_' + H + '_0'):
            if not collapse:
                if len(meta_l) == 2 and type(meta_l) is list:
                    final_result = pd.DataFrame.from_dict(
                        result_dict_, orient='index')
                    final_result.columns = [re.sub(rs, '', q) for q in [
                        query + '_' + str(l) for l in final_result.columns]]
                else:
                    final_result = pd.DataFrame.from_dict(
                        result_dict_, orient='index')
                    final_result.columns = [re.sub(rs, '', q) for q in [
                        query + '_' + str(l) for l in final_result.columns]]
            else:
                if len(meta_l) == 2 and type(meta_l) is list:
                    final_result = pd.DataFrame.from_dict(
                        result_dict_, orient='index')
                    final_result.columns = [query]
                else:
                    final_result = pd.DataFrame.from_dict(
                        result_dict_, orient='index')
                    final_result.columns = [query]
        else:
            if not collapse:
                if len(meta_l) == 2 and type(meta_l) is list:
                    final_result = pd.DataFrame.from_dict(
                        result_dict_, orient='index').T
                    final_result.columns = [re.sub('_sequence_id', '', q) for q in [
                        query + '_' + str(l) for l in final_result.columns]]
                else:
                    final_result = pd.DataFrame.from_dict(
                        result_dict_, orient='index').T
                    final_result.columns = [re.sub(rs, '', q) for q in [
                        query + '_' + str(l) for l in final_result.columns]]
            else:
                if verbose:
                    warnings.warn(UserWarning('Query dtype is {}. Ignoring collapse = True and split = False.'.format(
                        meta_h['sequence_id_' + H + '_0'].dtype.name)))
                if len(meta_l) == 2 and type(meta_l) is list:
                    final_result = pd.DataFrame.from_dict(
                        result_dict_, orient='index').T
                    final_result.columns = [re.sub('_sequence_id', '', q) for q in [
                        query + '_' + str(l) for l in final_result.columns]]
                else:
                    typedict = {locus_dict3[locus]: 'VDJ',
                                locus_dict4[locus]: 'VJ'}
                    final_result = pd.DataFrame.from_dict(
                        result_dict_, orient='index').T
                    final_result.columns = [re.sub(rs, '', q) for q in [
                        query + '_' + typedict[l.split('_')[2]] + '_' + l for l in final_result.columns]]

    return(final_result)


def retrieve_result_dict_singular(query: str, data: pd.DataFrame, meta_h: pd.DataFrame, locus: Literal['ig', 'tr-ab', 'tr-gd'] = 'ig', collapse: bool = True, combine: bool = False, verbose: bool = False) -> Dict:

    df_hl = defaultdict(dict)

    locus_dict1 = {'ig': 'IGH', 'tr-ab': 'TRB', 'tr-gd': 'TRD'}
    locus_dict3 = {'ig': 'H', 'tr-ab': 'B', 'tr-gd': 'D'}

    H = locus_dict3[locus]
    if meta_h.shape[1] > 1:
        if collapse:
            newmeta_h = meta_h.copy()
            if type_check(meta_h, 'sequence_id_' + H + '_0'):
                newh = []
                for i in meta_h.index:
                    try:
                        newh.append(
                            '|'.join([h for h in list(dict.fromkeys(newmeta_h.loc[i])) if h == h]))
                    except:
                        newh.append('|'.join([str(h) for h in list(
                            dict.fromkeys(newmeta_h.loc[i])) if h == h]))
                newmeta_h['sequence_id_' + H + '_0'] = newh
                meta_h = pd.DataFrame(newmeta_h['sequence_id_' + H + '_0'].copy())
            else:
                collapse = False
                if verbose:
                    warnings.warn(UserWarning(
                        'Multiple VDJ contigs mapping to the same cell barcode and/or query dtype is {}. Ignoring collapse = True.'.format(meta_h['sequence_id_' + H + '_0'].dtype.name)))

    metadata_ = meta_h.copy()

    df_ = metadata_.copy()
    if type_check(meta_h, 'sequence_id_' + H + '_0'):
        df_res = df_.copy()
        q_res = []
        if collapse:
            for i in metadata_.index:
                if combine:
                    try:
                        q_res.append(
                            '|'.join([qq for qq in list(dict.fromkeys(df_.loc[i])) if qq == qq]))
                    except:
                        q_res.append('|'.join([str(qq) for qq in list(
                            dict.fromkeys(df_.loc[i])) if qq == qq]))
                else:
                    try:
                        q_res.append(
                            '|'.join([qq for qq in list(df_.loc[i]) if qq == qq]))
                    except:
                        q_res.append(
                            '|'.join([str(qq) for qq in list(df_.loc[i]) if qq == qq]))
        else:
            for i in metadata_.index:
                q_res.append([qq for qq in list(df_.loc[i]) if qq == qq])
        df_res[query] = q_res
        result_dict_ = dict(df_res[query])
    else:
        result_dict_ = {x: dict(df_[x]) for x in df_}

    rs_dict1 = {'ig': '[H]', 'tr-ab': '[B]', 'tr-gd': '[D]'}
    rs = '_sequence_id_' + rs_dict1[locus]
    final_result_ = defaultdict(dict)

    if type_check(meta_h, 'sequence_id_' + H + '_0'):
        if not collapse:
            final_result = pd.DataFrame.from_dict(result_dict_, orient='index')
            final_result.columns = [re.sub(rs, '', q) for q in [
                query + '_' + str(l) for l in final_result.columns]]
        else:
            final_result = pd.DataFrame.from_dict(result_dict_, orient='index')
            final_result.columns = [query]
    else:
        if not collapse:
            final_result = pd.DataFrame.from_dict(
                result_dict_, orient='index').T
            final_result.columns = [re.sub(rs, '', q) for q in [
                query + '_' + str(l) for l in final_result.columns]]
        else:
            if verbose:
                warnings.warn(UserWarning('Query dtype is {}. Ignoring collapse = True and split = False.'.format(
                    meta_h['sequence_id_' + H + '_0'].dtype.name)))
            typedict = {locus_dict3[locus]: 'VDJ'}
            final_result = pd.DataFrame.from_dict(
                result_dict_, orient='index').T
            final_result.columns = [re.sub(rs, '', q) for q in [
                query + '_' + typedict[l.split('_')[2]] + '_' + l for l in final_result.columns]]
    return(final_result)


def initialize_metadata(self, cols: Sequence, locus_: str, clonekey: str, collapse_alleles: bool, verbose: bool) -> Dandelion:
    init_dict = {}
    for col in cols:
        init_dict.update({col: {'split': True, 'collapse': True,
                                'combine': False, 'locus': locus_, 'split_locus': False}})
    if clonekey in init_dict:
        init_dict.update({clonekey: {'split': False, 'collapse': True,
                                     'combine': True, 'locus': locus_, 'split_locus': False}})
    if 'sample_id' in init_dict:
        init_dict.update({'sample_id': {'split': False, 'collapse': True,
                                        'combine': True, 'locus': locus_, 'split_locus': False}})
    meta_ = defaultdict(dict)
    for k, v in init_dict.copy().items():
        if (all(pd.isnull(self.data[k]))) or (all([i == '' for i in self.data[k]])):
            init_dict.pop(k)
            continue
        meta_[k] = retrieve_metadata(self.data, query=k, verbose=verbose, **v)

    tmp_metadata = pd.concat(meta_.values(), axis=1, join="inner")

    if 'locus_VDJ' in tmp_metadata:
        suffix_h = '_VDJ'
        suffix_l = '_VJ'
    else:
        suffix_h = ''
        suffix_l = ''

    if clonekey in init_dict:
        tmp_metadata[str(clonekey)] = tmp_metadata[str(
            clonekey)].replace('', 'unassigned')
        clones = tmp_metadata[str(clonekey)].str.split('|', expand=False)
        tmpclones = []
        for i in clones:
            while 'unassigned' in i:
                i.remove('unassigned')
                if len(i) == 1:
                    break
            tmpclones.append(i)
        tmpclones = ['|'.join(list(set(x))) for x in tmpclones]
        tmpclonesdict = dict(zip(tmp_metadata.index, tmpclones))
        tmp_metadata[str(clonekey)] = pd.Series(tmpclonesdict)
        tmp = tmp_metadata[str(clonekey)].str.split('|', expand=True).stack()
        tmp = tmp.reset_index(drop=False)
        tmp.columns = ['cell_id', 'tmp', str(clonekey)]
        clone_size = tmp[str(clonekey)].value_counts()
        if "" in clone_size.index:
            clone_size = clone_size.drop("", axis=0)
        clonesize_dict = dict(clone_size)
        size_of_clone = pd.DataFrame.from_dict(clonesize_dict, orient='index')
        size_of_clone.reset_index(drop=False, inplace=True)
        size_of_clone.columns = [str(clonekey), 'clone_size']
        size_of_clone[str(clonekey) + '_by_size'] = size_of_clone.index + 1
        size_dict = dict(
            zip(size_of_clone[clonekey], size_of_clone[str(clonekey) + '_by_size']))
        size_dict.update({'': 'unassigned'})
        tmp_metadata[str(clonekey) + '_by_size'] = ['|'.join(sorted(list(set([str(size_dict[c_]) for c_ in c.split('|')]))))
                                                    if len(c.split('|')) > 1 else str(size_dict[c]) for c in tmp_metadata[str(clonekey)]]
        tmp_metadata[str(
            clonekey) + '_by_size'] = tmp_metadata[str(clonekey) + '_by_size'].astype('category')
        tmp_metadata = tmp_metadata[[str(clonekey), str(clonekey) + '_by_size'] + [
            cl for cl in tmp_metadata if cl not in [str(clonekey), str(clonekey) + '_by_size']]]

    for i in tmp_metadata.index:
        if pd.notnull(tmp_metadata.loc[i, 'locus' + suffix_h]):
            if pd.notnull(tmp_metadata.loc[i, 'locus' + suffix_l]):
                if (tmp_metadata.loc[i, 'locus' + suffix_l] != ''):
                    tmp_metadata.at[i, 'status'] = tmp_metadata.loc[i,
                                                                    'locus' + suffix_h] + ' + ' + tmp_metadata.loc[i, 'locus' + suffix_l]
                else:
                    tmp_metadata.at[i, 'status'] = tmp_metadata.loc[i,
                                                                    'locus' + suffix_h] + '_only'
            elif tmp_metadata.loc[i, 'locus' + suffix_h] != '':
                tmp_metadata.at[i, 'status'] = tmp_metadata.loc[i,
                                                                'locus' + suffix_h] + '_only'
            else:
                tmp_metadata.at[i, 'status'] = 'unassigned'
        else:
            tmp_metadata.at[i, 'status'] = 'unassigned'
    tmp_metadata['status_summary'] = [
        'Multi' if '|' in i else i for i in tmp_metadata['status']]

    for i in tmp_metadata.index:
        if tmp_metadata.loc[i, 'productive' + suffix_h] == tmp_metadata.loc[i, 'productive' + suffix_h]:
            if not pd.isnull(tmp_metadata.loc[i, 'productive' + suffix_l]):
                if tmp_metadata.loc[i, 'productive' + suffix_l] != '':
                    tmp_metadata.at[i, 'productive'] = tmp_metadata.loc[i,
                                                                        'productive' + suffix_h] + ' + ' + tmp_metadata.loc[i, 'productive' + suffix_l]
                else:
                    tmp_metadata.at[i, 'productive'] = tmp_metadata.loc[i,
                                                                        'productive' + suffix_h]
            elif tmp_metadata.loc[i, 'productive' + suffix_h] != '':
                tmp_metadata.at[i, 'productive'] = tmp_metadata.loc[i,
                                                                    'productive' + suffix_h]
            else:
                tmp_metadata.at[i, 'productive'] = 'unassigned'
        else:
            tmp_metadata.at[i, 'productive'] = 'unassigned'
    tmp_metadata['productive_summary'] = [
        'Multi' if '|' in i else i for i in tmp_metadata['productive']]

    conversion_dict = {'igha1': 'IgA', 'igha2': 'IgA', 'ighm': 'IgM', 'ighd': 'IgD', 'ighe': 'IgE', 'ighg1': 'IgG', 'ighg2': 'IgG', 'ighg3': 'IgG', 'ighg4': 'IgG', 'igkc': 'IgK', 'iglc1': 'IgL', 'iglc2': 'IgL', 'iglc3': 'IgL', 'iglc4': 'IgL', 'iglc5': 'IgL', 'iglc6': 'IgL', 'iglc7': 'IgL',
                       'igha': 'IgA', 'ighg': 'IgG', 'iglc': 'IgL', 'trac': 'unassigned', 'trbc': 'unassigned', 'trgc': 'unassigned', 'trbc1': 'unassigned', 'trbc2': 'unassigned', 'trgc1': 'unassigned', 'trgc2': 'unassigned', 'trgc3': 'unassigned', 'trgc4': 'unassigned', 'trdc': 'unassigned', 'nan': 'unassigned', 'na': 'unassigned', 'none': 'unassigned', '': 'unassigned', 'unassigned': 'unassigned', np.nan: 'unassigned', None: 'unassigned'}
    isotype = []
    for k in tmp_metadata['c_call' + suffix_h]:
        if isinstance(k, str):
            if ',' in k:
                k = '|'.join(k.split(','))
            if '|' in k:
                isotype.append('|'.join([str(z) for z in [conversion_dict[y.lower()] for y in set(
                    [re.sub('[0-9]', '', x) for x in k.split('|')])]]))
            else:
                isotype.append(conversion_dict[k.lower()])
        else:
            isotype.append('unassigned')
    tmp_metadata['isotype'] = isotype
    tmp_metadata['isotype_summary'] = [i if i == 'IgM|IgD' or i == 'IgD|IgM' else
                                       'Multi' if '|' in i else i for i in tmp_metadata['isotype']]

    vdj_gene_calls = ['v_call', 'd_call', 'j_call']
    if collapse_alleles:
        for x in vdj_gene_calls:
            if x in self.data:
                for c in tmp_metadata:
                    if x in c:
                        tmp_metadata[c] = ['|'.join(['|'.join(list(set(yy.split(',')))) for yy in list(set(
                            [re.sub('[*][0-9][0-9]', '', tx) for tx in t.split('|')]))]) for t in tmp_metadata[c]]
    multi = {}
    multic = {}
    for i in tmp_metadata.index:
        try:
            if 'v_call_genotyped' in cols:
                hv_ = tmp_metadata.at[i,
                                      'v_call_genotyped' + suffix_h].split('|')
            else:
                hv_ = tmp_metadata.at[i, 'v_call' + suffix_h].split('|')
        except:
            if 'v_call_genotyped' in cols:
                hv_ = tmp_metadata.at[i, 'v_call_genotyped' + suffix_h]
            else:
                hv_ = tmp_metadata.at[i, 'v_call' + suffix_h]
        try:
            hj_ = tmp_metadata.at[i, 'j_call' + suffix_h].split('|')
        except:
            hj_ = tmp_metadata.at[i, 'j_call' + suffix_h]
        try:
            if 'v_call_genotyped' in cols:
                lv_ = tmp_metadata.at[i,
                                      'v_call_genotyped' + suffix_l].split('|')
            else:
                lv_ = tmp_metadata.at[i, 'v_call' + suffix_l].split('|')
        except:
            if 'v_call_genotyped' in cols:
                lv_ = tmp_metadata.at[i, 'v_call_genotyped' + suffix_l]
            else:
                lv_ = tmp_metadata.at[i, 'v_call' + suffix_l]
        try:
            lj_ = tmp_metadata.at[i, 'j_call' + suffix_l].split('|')
        except:
            lj_ = tmp_metadata.at[i, 'j_call' + suffix_l]
        try:
            hc_ = tmp_metadata.at[i, 'c_call' + suffix_h].split('|')
        except:
            hc_ = tmp_metadata.at[i, 'c_call' + suffix_h]
        try:
            lc_ = tmp_metadata.at[i, 'c_call' + suffix_l].split('|')
        except:
            lc_ = tmp_metadata.at[i, 'c_call' + suffix_l]
        multi_h = []
        multi_l = []
        multi_hc = []
        multi_lc = []
        if len(hv_) > 1:
            multi_h.append(['Multi' + suffix_h + '_v'])
        if len(hj_) > 1:
            multi_h.append(['Multi' + suffix_h + '_j'])
        if len(lv_) > 1:
            multi_l.append(['Multi' + suffix_l + '_v'])
        if len(lj_) > 1:
            multi_l.append(['Multi' + suffix_l + '_j'])
        if len(hc_) > 1:
            if (tmp_metadata.at[i, 'isotype_summary'] == 'IgM|IgD') or (tmp_metadata.at[i, 'isotype_summary'] == 'IgD|IgM'):
                multi_hc.append([tmp_metadata.at[i, 'isotype_summary']])
            else:
                multi_hc.append(['Multi' + suffix_h + '_c'])
        if len(lc_) > 1:
            multi_lc.append(['Multi' + suffix_l + '_c'])
        if len(multi_hc) < 1:
            multi_hc.append(['Single'])
        if len(multi_h) < 1:
            multi_h.append(['Single'])
        if (len(lv_) == 1) & (len(lj_) == 1):
            if ('' not in lv_) and ('' not in lj_):
                if len(multi_l) < 1:
                    multi_l.append(['Single'])
        if len(lc_) == 1:
            if ('' not in lc_):
                if len(multi_lc) < 1:
                    multi_lc.append(['Single'])
        multih = '|'.join(list(set(flatten(multi_h))))
        multil = '|'.join(list(set(flatten(multi_l))))
        multihc = '|'.join(list(set(flatten(multi_hc))))
        multilc = '|'.join(list(set(flatten(multi_lc))))
        if len(multih) > 0:
            if len(multil) > 0:
                multi[i] = multih + ' + ' + multil
            else:
                multi[i] = multih
        else:
            multi[i] = 'unassigned'
        if len(multihc) > 0:
            if len(multilc) > 0:
                multic[i] = multihc + ' + ' + multilc
            else:
                multic[i] = multihc
        else:
            multic[i] = 'unassigned'
    tmp_metadata['vdj_status'] = pd.Series(multi)
    tmp_metadata['vdj_status_summary'] = [
        'Multi' if 'Multi' + suffix_h in i else 'Single' for i in tmp_metadata['vdj_status']]
    tmp_metadata['VDJ_chain_status_summary'] = [
        'Multi' if 'Multi' + suffix_h in i else 'Single' for i in pd.Series(multic)]

    self.metadata = tmp_metadata.copy()


def update_metadata(self: Dandelion, retrieve: Union[None, Sequence, str] = None, locus: Union[None, Literal['ig', 'tr-ab', 'tr-gd']] = None, clone_key: Union[None, str] = None, split: bool = True, collapse: bool = True, combine: bool = True, split_locus: bool = False, collapse_alleles: bool = True, reinitialize: bool = False, verbose: bool = False) -> Dandelion:
    """
    A `Dandelion` initialisation function to update and populate the `.metadata` slot.

    Parameters
    ----------
    self : Dandelion
        `Dandelion` object.
    retrieve : str, sequence, optional
        Column name in `.data` slot to retrieve and update the metadata.
    locus : str, optional
        Mode for creating metadata. Accepts one of 'ig', 'tr-ab' or 'tr-gd'. None defaults to 'ig'.
    clone_key : str, optional
        Column name of clone id. None defaults to 'clone_id'.
    split : bool
        Only applies if retrieve option is not None. Returns the retrieval splitted into two columns, e.g. one for heavy and one for light chains in BCR data, if True. Interacts with collapse, combine and split_locus options.
    collapse : bool
        Only applies if retrieve option is not None. Returns the retrieval as a collapsed entry if multiple entries are found (e.g. v genes from multiple contigs) where each entry will be separated by a '|' if True. Interacts with split, combine and split_locus options.
    combine : bool
        Only applies if retrieve option is not None. Returns the retrieval as a collapsed entry with only unique entries (separated by a '|' if multiple are found). Interacts with split, collapse, and split_locus options.
    split_locus : bool
        Only applies if retrieve option is not None. Similar to split except it returns the retrieval splitted into multiple columns corresponding to each unique element in the 'locus' column (e.g. IGH, IGK, IGL). Interacts with split, collapse, and combine options.
    collapse_alleles : bool
        Returns the V-D-J genes with allelic calls if False.
    reinitialize : bool
        Whether or not to reinitialize the current metadata. Useful when updating older versions of `dandelion` to newer version.
    verbose : bool
        Whether or not to print warning messages when constructing object.
    Returns
    -------
    `Dandelion` object with `.metadata` slot initialized.
    """

    if locus is None:
        locus_ = 'ig'
    else:
        locus_ = locus

    if clone_key is None:
        clonekey = 'clone_id'
    else:
        clonekey = clone_key

    cols = ['sequence_id', 'cell_id', 'locus', 'productive',
            'v_call', 'j_call', 'c_call', 'duplicate_count', 'junction_aa']

    if 'duplicate_count' not in self.data:
        try:
            self.data['duplicate_count'] = self.data['umi_count']
        except:
            cols = list(map(lambda x: 'umi_count' if x == 'duplicate_count' else x, cols))
            if 'umi_count' not in self.data:
                raise ValueError(
                    "Unable to initialize metadata due to missing keys. Please ensure either 'umi_count' or 'duplicate_count' is in the input data.")

    if not all([c in self.data for c in cols]):
        raise ValueError(
            'Unable to initialize metadata due to missing keys. Please ensure the input data contains all the following columns: {}'.format(cols))

    if 'sample_id' in self.data:
        cols = ['sample_id'] + cols

    if 'v_call_genotyped' in self.data:
        cols = list(
            map(lambda x: 'v_call_genotyped' if x == 'v_call' else x, cols))

    for c in ['sequence_id', 'cell_id']:
        cols.remove(c)

    if clonekey in self.data:
        if not all(pd.isnull(self.data[clonekey])):
            cols = [clonekey] + cols

    metadata_status = self.metadata
    if (metadata_status is None) or reinitialize:
        initialize_metadata(self, cols, locus_, clonekey,
                            collapse_alleles, verbose)

    tmp_metadata = self.metadata.copy()

    if retrieve is not None:
        ret_dict = {}
        if type(retrieve) is str:
            retrieve = [retrieve]
        for ret in retrieve:
            ret_dict.update({ret: {'split': split, 'collapse': collapse,
                                   'combine': combine, 'locus': locus_, 'split_locus': split_locus}})

        vdj_gene_ret = ['v_call', 'd_call', 'j_call']

        retrieve_ = defaultdict(dict)
        for k, v in ret_dict.items():
            if k in self.data.columns:
                retrieve_[k] = retrieve_metadata(
                    self.data, query=k, verbose=verbose, **v)
            else:
                raise KeyError(
                    'Cannot retrieve \'%s\' : Unknown column name.' % k)
        ret_metadata = pd.concat(retrieve_.values(), axis=1, join="inner")

        if collapse_alleles:
            for k in ret_dict.keys():
                if k in vdj_gene_ret:
                    for c in ret_metadata:
                        if k in c:
                            ret_metadata[c] = ['|'.join(['|'.join(list(set(yy.split(',')))) for yy in list(set(
                                [re.sub('[*][0-9][0-9]', '', tx) for tx in t.split('|')]))]) for t in ret_metadata[c]]

        for r in ret_metadata:
            tmp_metadata[r] = pd.Series(ret_metadata[r])
        self.metadata = tmp_metadata.copy()


def load_data(obj: Union[pd.DataFrame, str]) -> pd.DataFrame:
    """
    Reads in or copy dataframe object and set sequence_id as index without dropping.

    Parameters
    ----------
    obj : DataFrame, str
        file path to .tsv file or pandas DataFrame object.

    Returns
    -------
    pandas DataFrame object.
    """
    if os.path.isfile(str(obj)):
        try:
            obj_ = pd.read_csv(obj, sep='\t')
        except FileNotFoundError as e:
            print(e)
    elif isinstance(obj, pd.DataFrame):
        obj_ = obj.copy()
    else:
        raise TypeError(
            "Either input is not of <class 'pandas.core.frame.DataFrame'> or file does not exist.")

    if 'sequence_id' in obj_.columns:
        obj_.set_index('sequence_id', drop=False, inplace=True)
    else:
        raise KeyError("'sequence_id' not found in columns of input")

    return(obj_)
