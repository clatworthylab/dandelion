#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 14:01:32
# @Last Modified by:   Kelvin
# @Last Modified time: 2021-05-20 15:51:47

import os
import pandas as pd
import numpy as np
import scipy.sparse
import networkx as nx
import bz2
import gzip
import _pickle as cPickle
from ..utilities._utilities import *
from ..utilities._core import *
from typing import Union, Sequence, Tuple


def fasta_iterator(fh: str):
    '''
    Read in a fasta file as an iterator
    '''
    while True:
        line = fh.readline()
        if line.startswith('>'):
            break
    while True:
        header = line[1:-1].rstrip()
        sequence = fh.readline().rstrip()
        while True:
            line = fh.readline()
            if not line:
                break
            if line.startswith('>'):
                break
            sequence += line.rstrip()
        yield(header, sequence)
        if not line:
            return


def Write_output(out: str, file: str):
    '''
    general line writer
    '''
    fh = open(file, "a")
    fh.write(out)
    fh.close()
    return()


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


def read_pkl(filename: str = 'dandelion_data.pkl.pbz2') -> Dandelion:
    """
    Reads in and returns a `Dandelion` class saved using pickle format.

    Parameters
    ----------
    filename : str
        path to `.pkl` file. Depending on the extension, it will try to unzip accordingly.

    Returns
    -------
    Dandelion object.
    """
    if isBZIP(filename):
        data = bz2.BZ2File(filename, 'rb')
        data = cPickle.load(data)
    elif isGZIP(filename):
        data = gzip.open(filename, 'rb')
        data = cPickle.load(data)
    else:
        with open(filename, 'rb') as f:
            data = cPickle.load(f)
    return(data)


def read_h5(filename: str = 'dandelion_data.h5') -> Dandelion:
    """
    Reads in and returns a `Dandelion` class from .h5 format.

    Parameters
    ----------
    filename : str
        path to `.h5` file

    Returns
    -------
    `Dandelion` object.
    """
    try:
        data = pd.read_hdf(filename, 'data')
    except:
        raise AttributeError(
            '{} does not contain attribute `data`'.format(filename))
    try:
        metadata = pd.read_hdf(filename, 'metadata')
    except:
        pass

    try:
        edges = pd.read_hdf(filename, 'edges')
    except:
        pass

    try:
        g_0 = pd.read_hdf(filename, 'graph/graph_0')
        g_1 = pd.read_hdf(filename, 'graph/graph_1')
        g_0 = g_0 + 1
        g_0 = g_0.fillna(0)
        g_1 = g_1 + 1
        g_1 = g_1.fillna(0)
        graph0 = nx.from_pandas_adjacency(g_0)
        graph1 = nx.from_pandas_adjacency(g_1)
        for u, v, d in graph0.edges(data=True):
            d['weight'] = d['weight'] - 1
        for u, v, d in graph1.edges(data=True):
            d['weight'] = d['weight'] - 1
        graph = (graph0, graph1)
    except:
        pass

    with h5py.File(filename, 'r') as hf:
        try:
            layout0 = {}
            for k in hf['layout/layout_0'].attrs.keys():
                layout0.update({k: np.array(hf['layout/layout_0'].attrs[k])})
            layout1 = {}
            for k in hf['layout/layout_1'].attrs.keys():
                layout1.update({k: np.array(hf['layout/layout_1'].attrs[k])})
            layout = (layout0, layout1)
        except:
            pass

        germline = {}
        try:
            for g in hf['germline'].attrs:
                germline.update({g: hf['germline'].attrs[g]})
        except:
            pass

        distance = Tree()
        try:
            for d in hf['distance'].keys():
                d_ = pd.read_hdf(filename, 'distance/' + d)
                distance[d] = scipy.sparse.csr_matrix(d_.values)
        except:
            pass

        try:
            threshold = np.float(np.array(hf['threshold']))
        except:
            threshold = None

    constructor = {}
    constructor['data'] = data
    if 'metadata' in locals():
        constructor['metadata'] = metadata
    if 'germline' in locals():
        constructor['germline'] = germline
    if 'edges' in locals():
        constructor['edges'] = edges
    if 'distance' in locals():
        constructor['distance'] = distance
    if 'layout' in locals():
        constructor['layout'] = layout
    if 'graph' in locals():
        constructor['graph'] = graph
    try:
        res = Dandelion(**constructor)
    except:
        res = Dandelion(**constructor, initialize=False)

    if 'threshold' in locals():
        res.threshold = threshold
    else:
        pass
    return(res)


def read_10x_airr(file: str) -> Dandelion:
    """
    Reads the 10x AIRR rearrangement .tsv directly and returns a `Dandelion` object.

    Parameters
    ----------
    file : str
        path to `airr_rearrangement.tsv`

    Returns
    -------
    `Dandelion` object of pandas data frame.

    """
    dat = load_data(file)
    # get all the v,d,j,c calls
    if 'locus' not in dat:
        tmp = [(v, d, j, c) for v, d, j, c in zip(
            dat['v_call'], dat['d_call'], dat['j_call'], dat['c_call'])]
        locus = []
        for t in tmp:
            if all('IGH' in x for x in t if x == x):
                locus.append('IGH')
            elif all('IGK' in x for x in t if x == x):
                locus.append('IGK')
            elif all('IGL' in x for x in t if x == x):
                locus.append('IGL')
            elif all('TRA' in x for x in t if x == x):
                locus.append('TRA')
            elif all('TRB' in x for x in t if x == x):
                locus.append('TRB')
            elif all('TRD' in x for x in t if x == x):
                locus.append('TRD')
            elif all('TRG' in x for x in t if x == x):
                locus.append('TRG')
            else:
                locus.append(np.nan)
        dat['locus'] = locus

    return(Dandelion(dat))


def to_scirpy(data: Dandelion, transfer: bool = False) -> AnnData:
    """
    Converts a `Dandelion` object to scirpy's format.

    Parameters
    ----------
    data : Dandelion
        `Dandelion` object
    transfer : bool
        Whether to execute :func:`dandelion.tl.transfer` to transfer all data
        to the :class:`anndata.AnnData` instance.

    Returns
    -------
    `AnnData` object in the format initialized by `scirpy`.

    """
    try:
        import scirpy as ir
    except:
        raise ImportError('Please install scirpy. pip install scirpy')

    if 'duplicate_count' not in data.data and 'umi_count' in data.data:
        data.data['duplicate_count'] = data.data['umi_count']
    return(ir.io.from_dandelion(data, transfer))


def from_scirpy(adata: AnnData, clone_key: Union[None, str] = None, key_added: Union[None, str] = None, mapping_mode: Literal['chain', 'cell'] = 'chain') -> Dandelion:
    """
    Reads a `scirpy` initialized `AnnData` oject and returns a `Dandelion` object.

    Parameters
    ----------
    adata : AnnData
        `scirpy` initialized `AnnData` object.
    clone_key : str, optional
        column name for `clone_id` in `AnnData`. None defaults to `clone_id` in `scirpy` initialized object.
    key_added : str, optional
        column name for `clone_id` in `Dandelion`. None defaults to `clone_id` in `dandelion` initialized object.
    mapping_mode : str
        mode for retrieving the clone_id calls, either based on cells (all chains/contigs have the same call) or chains (allow for different calls between chains).

    Returns
    -------
    `Dandelion` object.

    """
    try:
        import scirpy as ir
    except:
        raise ImportError('Please install scirpy. pip install scirpy')

    if clone_key is None:
        clonekey_s = 'clone_id'
    else:
        clonekey_s = clone_key

    if key_added is None:
        clonekey_d = 'clone_id'
    else:
        clonekey_d = key_added

    return(ir.io.to_dandelion(adata))


def read_10x_vdj(path: str, filtered: bool = True):
    """
    A wrapper from scirpy to read 10x's .csv and .json files directly to be formatted in dandelion.

    Parameters
    ----------
    path : str
        Path to `filterd_contig_annotations.csv`, `all_contig_annotations.csv` or `all_contig_annotations.json`.
    filtered : bool
        Only keep filtered contig annotations (i.e. `is_cell` and `high_confidence`).
        If using `filtered_contig_annotations.csv` already, this option is futile.
    Returns
    -------
    `Dandelion` object.

    """
    try:
        import scirpy as ir
    except:
        raise ImportError('Please install scirpy. pip install scirpy')

    adata = ir.io.read_10x_vdj(path, filtered=filtered)

    return(ir.io.to_dandelion(adata))


def concat(arrays: Sequence[Union[pd.DataFrame, Dandelion]], check_unique: bool = True) -> Dandelion:
    """
    Concatenate dataframe and return as `Dandelion` object.

    Parameters
    ----------
    arrays : Sequence
        List of `Dandelion` class objects or pandas dataframe
    check_unique : bool
        Check the new index for duplicates. Otherwise defer the check until necessary. Setting to False will improve the performance of this method.

    Returns
    -------
        `Dandelion` object
    """
    arrays = list(arrays)

    try:
        arrays_ = [x.data.copy() for x in arrays]
    except:
        arrays_ = [x.copy() for x in arrays]

    if check_unique:
        try:
            df = pd.concat(arrays_, verify_integrity=True)
        except:
            for i in range(0, len(arrays)):
                arrays_[i]['sequence_id'] = [x + '__' +
                                             str(i) for x in arrays_[i]['sequence_id']]
            arrays_ = [load_data(x) for x in arrays_]
            df = pd.concat(arrays_, verify_integrity=True)
    else:
        df = pd.concat(arrays_)
    try:
        out = Dandelion(df)
    except:
        out = Dandelion(df, initialize=False)
    return(out)
