#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 14:01:32
# @Last Modified by:   Kelvin
# @Last Modified time: 2021-02-18 23:53:56

import sys
import os
from collections import defaultdict, Iterable
import pandas as pd
import numpy as np
from subprocess import run
from tqdm import tqdm
import re
import copy
from changeo.IO import readGermlines
import warnings
import scipy.sparse
import tables
import h5py
import networkx as nx
import bz2
import gzip
import _pickle as cPickle
try:
    from scanpy import logging as logg
except ImportError:
    pass

class Tree(defaultdict):
    '''
    Create a recursive defaultdict
    '''
    def __init__(self, value=None):
        super(Tree, self).__init__(Tree)
        self.value = value

def fasta_iterator(fh):
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

def Write_output(out, file):
    '''
    general line writer
    '''
    fh = open(file, "a")
    fh.write(out)
    fh.close()
    return()

def dict_from_table(meta, columns):
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
        meta_ = pd.read_csv(meta, sep = '\t', dtype = 'object')
        if len(columns) == 2:
            sample_dict = dict(zip(meta_[columns[0]], meta_[columns[1]]))

    sample_dict = clean_nan_dict(sample_dict)
    return(sample_dict)

def clean_nan_dict(d):
    """
    Parameters
    ----------
    d
        dictionary
    
    Returns
    -------
    dictionary with no NAs.
    """

    return {
        k:v
        for k, v in d.items()
        if v is not np.nan
    }

def flatten(l):
    """
    Parameters
    ----------
    l
        list
    
    Returns
    -------
    a flattened list.
    """
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el

def makeblastdb(ref):
    """
    Runs makeblastdb.

    Parameters
    ----------
    ref
        fasta file

    """

    cmd = ['makeblastdb',
               '-dbtype', 'nucl',
               '-parse_seqids',
               '-in', ref]
    run(cmd)

def bh(pvalues):
    """
    Computes the Benjamini-Hochberg FDR correction.

    Input:
        * pvals - vector of p-values to correct
    """
    n = int(pvalues.shape[0])
    new_pvalues = np.empty(n)
    values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
    values.sort()
    values.reverse()
    new_values = []
    for i, vals in enumerate(values):
        rank = n - i
        pvalue, index = vals
        new_values.append((n/rank) * pvalue)
    for i in range(0, int(n)-1):
        if new_values[i] < new_values[i+1]:
            new_values[i+1] = new_values[i]
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]
    return new_pvalues

# def extract(d, keys):
#    return(dict((k, d[k]) for k in keys if k in d))

def load_data(obj):
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
            obj_ = pd.read_csv(obj, sep = '\t')
        except FileNotFoundError as e:
            print(e)
    elif isinstance(obj, pd.DataFrame):
            obj_ = obj.copy()
    else:
        raise TypeError("Either input is not of <class 'pandas.core.frame.DataFrame'> or file does not exist.")

    if 'sequence_id' in obj_.columns:
        obj_.set_index('sequence_id', drop = False, inplace = True)
    else:
        raise KeyError("'sequence_id' not found in columns of input")

    return(obj_)

def setup_metadata_(data):
    """
    A Dandelion class subfunction to initialize the `.metadata` slot.
    Parameters
    ----------
    data : DataFrame
        pandas DataFrame object.    
    
    Returns
    -------
    pandas DataFrame object.
    """
    dat_h = data[data['locus'] == 'IGH'].copy()
    dict_h = dict(zip(dat_h['sequence_id'], dat_h['cell_id']))
    metadata_ = pd.DataFrame.from_dict(dict_h, orient = 'index', columns = ['cell_id'])
    metadata_.set_index('cell_id', inplace = True)
    return(metadata_)

def setup_metadata(data, clone_key = None):
    """
    A Dandelion class subfunction to initialize the `.metadata` slot.
    Parameters
    ----------
    data : DataFrame
        pandas DataFrame object.
    clone_key : str, optiona;
        column name of clone id. None defaults to 'clone_id'.
    
    Returns
    -------
    pandas DataFrame object.
    """
    if clone_key is None:
        clonekey = 'clone_id'
    else:
        clonekey = clone_key

    dat_h = data[data['locus'] == 'IGH'].copy()
    dat_l = data[data['locus'].isin(['IGK', 'IGL'])].copy()
    if dat_l.shape[0] == 0:
        clone_h = dict(zip(dat_h['sequence_id'], zip(dat_h['cell_id'], dat_h[clonekey])))
        metadata_ = pd.DataFrame.from_dict(clone_h, orient = 'index', columns = ['cell_id', 'heavy'])
        metadata_.set_index('cell_id', inplace = True)
        clones_list = {}
        for x in metadata_.index:
            cl = list(set(list(metadata_.loc[x, :])))
            cl = sorted([y for y in cl if str(y) != 'nan'])
            if len(cl) > 1:
                cl = cl[1:]
            clones_list[x] = '|'.join(cl)
        metadata_[clonekey] = pd.Series(clones_list)
        metadata_ = metadata_[[clonekey]]
        
        tmp = metadata_[str(clonekey)].str.split('|', expand=True).stack()
        tmp = tmp.reset_index(drop = False)
        tmp.columns = ['cell_id', 'tmp', str(clonekey)]
        clone_size = tmp[str(clonekey)].value_counts()
        clonesize_dict = dict(clone_size)

        # size_of_clone = pd.DataFrame(metadata_[str(clonekey)].value_counts())
        size_of_clone = pd.DataFrame.from_dict(clonesize_dict, orient = 'index')
        size_of_clone.reset_index(drop = False, inplace = True)
        size_of_clone.columns = [str(clonekey), 'clone_size']    
        size_of_clone[str(clonekey)+'_by_size'] = size_of_clone.index+1
        size_dict = dict(zip(size_of_clone['clone_id'], size_of_clone['clone_id_by_size']))
        metadata_[str(clonekey)+'_by_size'] = ['|'.join([str(size_dict[c_]) for c_ in c.split('|')]) if len(c.split('|')) > 1 else str(size_dict[c]) for c in metadata_[str(clonekey)]]
        metadata_[str(clonekey)+'_by_size'] = metadata_[str(clonekey)+'_by_size'].astype('category')
        return(metadata_)
    else:
        clone_h = dict(zip(dat_h['sequence_id'], zip(dat_h['cell_id'], dat_h[clonekey])))
        clone_l = dict(zip(dat_l['sequence_id'], zip(dat_l['cell_id'], dat_l[clonekey])))
        metadata_ = pd.DataFrame.from_dict(clone_h, orient = 'index', columns = ['cell_id', 'heavy'])
        metadata_.set_index('cell_id', inplace = True)
        light_clone_tree = Tree()
        for key, value in clone_l.items():
            k, v = value
            light_clone_tree[k][key] = v
        light_clone_tree2 = Tree()
        for g in light_clone_tree:
            second_key = []
            for k2 in light_clone_tree[g].keys():
                second_key.append(k2)
            second_key = list(set(second_key))
            second_key_dict = dict(zip(second_key, range(0,len(second_key))))
            for key, value in light_clone_tree[g].items():
                light_clone_tree2[g][second_key_dict[key]] = value
        metadata_['light'] = pd.Series(light_clone_tree2)
        tmp = pd.Series([dict(i) if i is not np.nan else {0:i} for i in metadata_['light']])
        tmp_dat = pd.DataFrame(tmp.tolist(), index = metadata_.index)
        tmp_dat.columns = ['light_' + str(c) for c in tmp_dat.columns]
        metadata_ = metadata_.merge(tmp_dat, left_index = True, right_index = True)
        metadata_ = metadata_[['heavy'] + [str(c) for c in tmp_dat.columns]]
        clones_list = {}
        for x in metadata_.index:
            cl = list(set(list(metadata_.loc[x, :])))
            cl = sorted([y for y in cl if str(y) != 'nan'])
            if len(cl) > 1:
                cl = cl[1:]
            clones_list[x] = '|'.join(cl)
        metadata_[clonekey] = pd.Series(clones_list)
        metadata_ = metadata_[[clonekey]]

        tmp = metadata_[str(clonekey)].str.split('|', expand=True).stack()
        tmp = tmp.reset_index(drop = False)
        tmp.columns = ['cell_id', 'tmp', str(clonekey)]
        clone_size = tmp[str(clonekey)].value_counts()
        clonesize_dict = dict(clone_size)

        # size_of_clone = pd.DataFrame(metadata_[str(clonekey)].value_counts())
        size_of_clone = pd.DataFrame.from_dict(clonesize_dict, orient = 'index')
        size_of_clone.reset_index(drop = False, inplace = True)
        size_of_clone.columns = [str(clonekey), 'clone_size']
        size_of_clone[str(clonekey)+'_by_size'] = size_of_clone.index+1
        size_dict = dict(zip(size_of_clone['clone_id'], size_of_clone['clone_id_by_size']))
        metadata_[str(clonekey)+'_by_size'] = ['|'.join([str(size_dict[c_]) for c_ in c.split('|')]) if len(c.split('|')) > 1 else str(size_dict[c]) for c in metadata_[str(clonekey)]]
        metadata_[str(clonekey)+'_by_size'] = metadata_[str(clonekey)+'_by_size'].astype('category')
        return(metadata_)

def retrieve_metadata(data, retrieve_id, split_heavy_light, collapse):
    """
    A Dandelion class subfunction to populate the `.metadata` slot.

    Parameters
    ----------
    data : DataFrame
        pandas DataFrame object.
    retrieve_id : str
        column name in `.data` slot.
    split_heavy_light : bool
        Returns the retrieval splitted into two column for heavy and light. Default is True. False combines the retrieval into a single column.
    collapse : bool
        Whether or not to collapse unique elements if duplicated. For example, different contigs and same sample id would then benefit from this option being set to True.

    Returns
    -------
    A dictionary with keys as cell_ids and records as retrieved value.
    """
    dat_h = data[data['locus'] == 'IGH'].copy()
    dat_l = data[data['locus'].isin(['IGK', 'IGL'])].copy()
    if dat_l.shape[0] == 0:
        retrieve_h = dict(zip(dat_h['sequence_id'], zip(dat_h['cell_id'], dat_h[retrieve_id])))
        sub_metadata = pd.DataFrame.from_dict(retrieve_h, orient = 'index', columns = ['cell_id', 'heavy'])
        sub_metadata.set_index('cell_id', inplace = True)
        heavy_retrieval_list = dict(sub_metadata['heavy'])
        for k, r in heavy_retrieval_list.items():
            if isinstance(r, pd.Series):
                heavy_retrieval_list[k] = ','.join([str(x) for x in flatten(r.to_list())])
        return(heavy_retrieval_list)
    else:
        retrieve_h = dict(zip(dat_h['sequence_id'], zip(dat_h['cell_id'], dat_h[retrieve_id])))
        retrieve_l = dict(zip(dat_l['sequence_id'], zip(dat_l['cell_id'], dat_l[retrieve_id])))
        sub_metadata = pd.DataFrame.from_dict(retrieve_h, orient = 'index', columns = ['cell_id', 'heavy'])
        sub_metadata.set_index('cell_id', inplace = True)
        light_retrieval_tree = Tree()
        for key, value in retrieve_l.items():
            k, v = value
            light_retrieval_tree[k][key] = v
        light_retrieval_tree2 = Tree()
        for g in light_retrieval_tree:
            second_key = []
            for k2 in light_retrieval_tree[g].keys():
                second_key.append(k2)
            second_key = list(set(second_key))
            second_key_dict = dict(zip(second_key, range(0,len(second_key))))
            for key, value in light_retrieval_tree[g].items():
                light_retrieval_tree2[g][second_key_dict[key]] = value
        sub_metadata['light'] = pd.Series(light_retrieval_tree2)
        tmp = pd.Series([dict(i) if i is not np.nan else {0:i} for i in sub_metadata['light']])
        tmp_dat = pd.DataFrame(tmp.tolist(), index = sub_metadata.index)
        tmp_dat.columns = ['light_' + str(c) for c in tmp_dat.columns]
        sub_metadata = sub_metadata.merge(tmp_dat, left_index = True, right_index = True)
        sub_metadata = sub_metadata[['heavy'] + [str(c) for c in tmp_dat.columns]]
        if not split_heavy_light:
            retrieval_list = {}
            for x in sub_metadata.index:
                if collapse:
                    r_l = list(set(list(sub_metadata.loc[x, :])))
                else:
                    r_l = list(sub_metadata.loc[x, :])
                r_l = sorted([y for y in r_l if str(y) != 'nan'])
                if len(r_l) > 1:
                    r_l = r_l[1:]
                retrieval_list[x] = '|'.join([str(r) for r in r_l])
            return(retrieval_list)
        else:
            heavy_retrieval_list = dict(sub_metadata['heavy'])
            for k, r in heavy_retrieval_list.items():
                if isinstance(r, pd.Series):
                    heavy_retrieval_list[k] = ','.join([str(x) for x in flatten(r.to_list())])
            light_retrieval_list = {}
            sub_metadata2 = sub_metadata.drop('heavy', axis = 1) # TODO: this will come back with issue if some heavy chain comes back with multiple indices that wasn't filtered
            for x in sub_metadata2.index:
                if collapse:
                    r_l = list(set(list(sub_metadata2.loc[x, :])))
                else:
                    r_l = list(sub_metadata2.loc[x, :])
                r_l = sorted([y for y in r_l if str(y) != 'nan'])
                light_retrieval_list[x] = ['|'.join([str(zz) for zz in z]) if len(z) > 0 else np.nan for z in [r_l]][0]
            return(heavy_retrieval_list, light_retrieval_list)

def update_metadata(self, retrieve = None, isotype_dict = None, split_heavy_light = True, collapse = False, clones_sep = None, clone_key = None):
    """
    A Dandelion function to update and populate the `.metadata` slot.

    Parameters
    ----------
    self : Dandelion
        Dandelion object
    retrieve : str
        column name in `.data` slot.
    split_heavy_light : bool
        Returns the retrieval splitted into two column for heavy and light. Default is True. False combines the retrieval into a single column.
    collapse : bool
        Whether or not to collapse unique elements if duplicated. For example, different contigs and same sample id would then benefit from this option being set to True.
    clones_sep : tuple[int, str]
        A tuple containing how the clone groups should be extracted. None defaults to (0, '_').
    
    Returns
    -------
    Dandelion object with `.metadata` slot initialized.
    """
    dat = load_data(self.data)
    for x in ['cell_id', 'locus', 'c_call', 'umi_count']:
        if x not in dat.columns:
            raise KeyError ("Please check your object. %s is not in the columns of input data." % x)

    if clone_key is None:
        clonekey = 'clone_id'
    else:
        clonekey = clone_key

    metadata_status = self.metadata
    if metadata_status is None:
        if clonekey in dat.columns:
            self.metadata = setup_metadata(dat, clonekey)
        else:
            self.metadata = setup_metadata_(dat)
    else:
        self.metadata = self.metadata.copy()

    if 'sample_id' in dat.columns:
        samp_id = retrieve_metadata(dat, 'sample_id', False, True)

    dat_l = dat[dat['locus'].isin(['IGK', 'IGL'])]
    if dat_l.shape[0] == 0:
        if 'v_call_genotyped' in dat.columns:
            heavy_v_call = retrieve_metadata(dat, 'v_call_genotyped', False, True) # specifying False/True shouldn't do anything
        else:
            heavy_v_call = retrieve_metadata(dat, 'v_call', False, True)
        heavy_j_call = retrieve_metadata(dat, 'j_call', False, True)
        heavy_c_call = retrieve_metadata(dat, 'c_call', False, True)
        heavy_umi = retrieve_metadata(dat, 'umi_count', False, True)
        heavy_status = retrieve_metadata(dat, 'locus', False, True)
        status = pd.DataFrame([heavy_status], index = ['heavy']).T
        for i in status.index:
            status.at[i, 'status'] = status.loc[i,'heavy'] + '_only'
        if isotype_dict is None:
            conversion_dict = {'igha1':'IgA', 'igha2':'IgA', 'ighm':'IgM', 'ighd':'IgD', 'ighm|ighd':'IgM|IgD', 'ighe':'IgE', 'ighg1':'IgG', 'ighg2':'IgG', 'ighg3':'IgG', 'ighg4':'IgG', 'igkc':'IgK', 'iglc1':'IgL', 'iglc2':'IgL', 'iglc3':'IgL', 'iglc4':'IgL', 'iglc5':'IgL', 'iglc6':'IgL', 'iglc7':'IgL', 'igha':'IgA', 'ighg':'IgG', 'iglc':'IgL', 'nan':np.nan, np.nan:np.nan, 'na':np.nan, '':np.nan} # the key for IgG being igh is on purpose because of how the counter works
        else:
            conversion_dict = isotype_dict
        isotype = {}
        for k in heavy_c_call:
            if heavy_c_call[k] == heavy_c_call[k]:
                if ',' in heavy_c_call[k]:
                    # iso_d = defaultdict(int)
                    # for c in heavy_c_call[k].lower():
                    #     iso_d[c] += 1
                    isotype[k] = ','.join([str(z) for z in [conversion_dict[y.lower()] for y in set([re.sub('[0-9]', '', x) for x in heavy_c_call[k].split(',')])]])
                else:
                    isotype[k] = conversion_dict[heavy_c_call[k].lower()]
            else:
                isotype[k] = heavy_c_call[k]
        for k in heavy_v_call:
            heavy_v_call[k] = ''.join([','.join(list(set([re.sub('[*][0-9][0-9]', '', str(heavy_v_call[k]))][0].split(','))))])
        for k in heavy_j_call:
            heavy_j_call[k] = ''.join([','.join(list(set([re.sub('[*][0-9][0-9]', '', str(heavy_j_call[k]))][0].split(','))))])        
        productive = retrieve_metadata(dat, 'productive', False, True)
        if 'sample_id' in dat.columns:
            self.metadata['sample_id'] = pd.Series(samp_id)
        self.metadata['isotype'] = pd.Series(isotype)
        self.metadata['status'] = pd.Series(status['status'])
        self.metadata['productive'] = pd.Series(productive)
        self.metadata['umi_counts_heavy'] = pd.Series(heavy_umi)
        self.metadata['c_call_heavy'] = pd.Series(heavy_c_call)
        self.metadata['v_call_heavy'] = pd.Series(heavy_v_call)
        self.metadata['j_call_heavy'] = pd.Series(heavy_j_call)        
        multi = {}
        for i in self.metadata.index:
            try:
                hv_ = self.metadata.at[i, 'v_call_heavy'].split(',')
            except:
                hv_ = self.metadata.at[i, 'v_call_heavy']
            try:
                hj_ = self.metadata.at[i, 'j_call_heavy'].split(',')
            except:
                hj_ = self.metadata.at[i, 'j_call_heavy']            
            multi_ = []
            if len(hv_) > 1:
                multi_.append(['Multi_heavy_v'])
            if len(hj_) > 1:
                multi_.append(['Multi_heavy_j'])            
            if len(multi_) < 1:
                multi_.append(['Single'])
            multi[i] = ','.join(list(flatten(multi_)))
        self.metadata['vdj_status'] = pd.Series(multi)
        # return this in this order
        if metadata_status is None:
            if clonekey in self.data.columns:
                if 'sample_id' in self.data.columns:
                    self.metadata = self.metadata[['sample_id', str(clonekey), str(clonekey)+'_by_size', 'isotype', 'status', 'vdj_status', 'productive',  'umi_counts_heavy', 'v_call_heavy', 'j_call_heavy', 'c_call_heavy']]
                else:
                    self.metadata = self.metadata[[str(clonekey), str(clonekey)+'_by_size', 'isotype', 'productive', 'status', 'vdj_status', 'umi_counts_heavy', 'v_call_heavy','j_call_heavy','c_call_heavy']]
            else:
                if 'sample_id' in self.data.columns:
                    self.metadata = self.metadata[['sample_id', 'isotype', 'status', 'vdj_status', 'productive',  'umi_counts_heavy', 'v_call_heavy','j_call_heavy','c_call_heavy']]
                else:
                    self.metadata = self.metadata[['isotype', 'productive', 'status', 'vdj_status', 'umi_counts_heavy',  'v_call_heavy','j_call_heavy','c_call_heavy']]
    
        # new function to retrieve non-standard columns
        if retrieve is not None:
            if type(retrieve) is str:
                retrieve = [retrieve]                
            for ret in retrieve:
                if ret in dat.columns:                                    
                    retrieve_dict = retrieve_metadata(dat, ret, False, True)
                    self.metadata[str(ret)+'_heavy'] = pd.Series(retrieve_dict)
                else:
                    raise KeyError('Unknown column : \'%s\' to retrieve.' % ret)
    else:
        if 'v_call_genotyped' in dat.columns:
            heavy_v_call, light_v_call = retrieve_metadata(dat, 'v_call_genotyped', True, False)
        else:
            heavy_v_call, light_v_call = retrieve_metadata(dat, 'v_call', True, False)
        heavy_j_call, light_j_call = retrieve_metadata(dat, 'j_call', True, False)
        heavy_c_call, light_c_call = retrieve_metadata(dat, 'c_call', True, False)
        heavy_umi, light_umi = retrieve_metadata(dat, 'umi_count', True, False)
        heavy_status, light_status = retrieve_metadata(dat, 'locus', True, False)
        status = pd.DataFrame([heavy_status, light_status], index = ['heavy', 'light']).T
        for i in status.index:
            try:
                status.at[i, 'status'] = status.loc[i,'heavy']+' + '+status.loc[i,'light']
            except:
                status.at[i, 'status'] = status.loc[i,'heavy'] + '_only'
        if isotype_dict is None:
            conversion_dict = {'igha1':'IgA', 'igha2':'IgA', 'ighm':'IgM', 'ighd':'IgD', 'ighm|ighd':'IgM|IgD', 'ighe':'IgE', 'ighg1':'IgG', 'ighg2':'IgG', 'ighg3':'IgG', 'ighg4':'IgG', 'igkc':'IgK', 'iglc1':'IgL', 'iglc2':'IgL', 'iglc3':'IgL', 'iglc4':'IgL', 'iglc5':'IgL', 'iglc6':'IgL', 'iglc7':'IgL', 'igha':'IgA', 'ighg':'IgG', 'iglc':'IgL', 'nan':np.nan, np.nan:np.nan, 'na':np.nan, '':np.nan} # the key for IgG being igh is on purpose because of how the counter works
        else:
            conversion_dict = isotype_dict
        isotype = {}
        for k in heavy_c_call:
            if heavy_c_call[k] == heavy_c_call[k]:
                if ',' in heavy_c_call[k]:
                    # iso_d = defaultdict(int)
                    # for c in heavy_c_call[k].lower():
                    #     iso_d[c] += 1
                    isotype[k] = ','.join([str(z) for z in [conversion_dict[y.lower()] for y in set([re.sub('[0-9]', '', x) for x in heavy_c_call[k].split(',')])]])
                else:
                    isotype[k] = conversion_dict[heavy_c_call[k].lower()]
            else:
                isotype[k] = heavy_c_call[k]
        lightchain = {}
        for k in light_c_call:
            if light_c_call[k] == light_c_call[k]:
                if '|' in light_c_call[k]:
                    if ',' in light_c_call[k]:
                        lc_y = []
                        for x in light_c_call[k].lower().split('|'):
                            iso_d = defaultdict(int)
                            for c in x:
                                iso_d[c] += 1
                            lc_y.append(re.sub(',|[0-9]', '', ''.join([k_ for k_,v_ in iso_d.items() if v_ == 1])))
                            lc_y.append(re.sub(',|[0-9]', '', ''.join([k_ for k_,v_ in iso_d.items() if v_ >= 2])))
                        if '' in lc_y:
                            lc_y = [x for x in lc_y if x != '']
                        lightchain[k] = '|'.join([conversion_dict[y] for y in lc_y])
                    else:
                        lightchain[k] = '|'.join([str(z) for z in [conversion_dict[x] for x in light_c_call[k].lower().split('|')]])
                else:
                    if ',' in light_c_call[k]:
                        iso_d = defaultdict(int)
                        for c in light_c_call[k].lower():
                            iso_d[c] += 1
                        lightchain[k] = conversion_dict[re.sub(',|[0-9]', '', ''.join([k_ for k_,v_ in iso_d.items() if v_ >= 2]))]
                    else:
                        lightchain[k] = conversion_dict[light_c_call[k].lower()]
            else:
                lightchain[k] = light_c_call[k]
        for k in heavy_v_call:
            heavy_v_call[k] = ''.join([','.join(list(set([re.sub('[*][0-9][0-9]', '', str(heavy_v_call[k]))][0].split(','))))])
        for k in heavy_j_call:
            heavy_j_call[k] = ''.join([','.join(list(set([re.sub('[*][0-9][0-9]', '', str(heavy_j_call[k]))][0].split(','))))])
        for k in light_v_call:
            light_v_call[k] = ''.join([','.join(list(set([re.sub('[*][0-9][0-9]', '', str(light_v_call[k]))][0].split(','))))])
        for k in light_j_call:
            light_j_call[k] = ''.join([','.join(list(set([re.sub('[*][0-9][0-9]', '', str(light_j_call[k]))][0].split(','))))])
        productive = retrieve_metadata(dat, 'productive', False, False)
        if 'sample_id' in dat.columns:
            self.metadata['sample_id'] = pd.Series(samp_id)
        self.metadata['isotype'] = pd.Series(isotype)
        self.metadata['lightchain'] = pd.Series(lightchain)
        self.metadata['status'] = pd.Series(status['status'])
        self.metadata['productive'] = pd.Series(productive)
        self.metadata['umi_counts_heavy'] = pd.Series(heavy_umi)
        self.metadata['umi_counts_light'] = pd.Series(light_umi)
        self.metadata['c_call_heavy'] = pd.Series(heavy_c_call)
        self.metadata['c_call_light'] = pd.Series(light_c_call)
        self.metadata['v_call_heavy'] = pd.Series(heavy_v_call)
        self.metadata['v_call_light'] = pd.Series(light_v_call)
        self.metadata['j_call_heavy'] = pd.Series(heavy_j_call)
        self.metadata['j_call_light'] = pd.Series(light_j_call)
        multi = {}
        for i in self.metadata.index:
            try:
                hv_ = self.metadata.at[i, 'v_call_heavy'].split('|')
            except:
                hv_ = self.metadata.at[i, 'v_call_heavy']
            try:
                hj_ = self.metadata.at[i, 'j_call_heavy'].split('|')
            except:
                hj_ = self.metadata.at[i, 'j_call_heavy']
            try:
                lv_ = self.metadata.at[i, 'v_call_light'].split('|')
            except:
                lv_ = self.metadata.at[i, 'v_call_light']
            try:
                lj_ = self.metadata.at[i, 'j_call_light'].split('|')
            except:
                lj_ = self.metadata.at[i, 'j_call_light']
            multi_ = []
            if len(hv_) > 1:
                multi_.append(['Multi_heavy_v'])
            if len(hj_) > 1:
                multi_.append(['Multi_heavy_j'])
            if len(lv_) > 1:
                multi_.append(['Multi_light_v'])
            if len(lj_) > 1:
                multi_.append(['Multi_light_j'])
            if len(multi_) < 1:
                multi_.append(['Single'])
            multi[i] = ','.join(list(flatten(multi_)))
        self.metadata['vdj_status'] = pd.Series(multi)
        # return this in this order
        if metadata_status is None:
            if clonekey in self.data.columns:
                if 'sample_id' in self.data.columns:
                    self.metadata = self.metadata[['sample_id', str(clonekey), str(clonekey)+'_by_size', 'isotype', 'lightchain', 'status', 'vdj_status', 'productive',  'umi_counts_heavy', 'umi_counts_light', 'v_call_heavy', 'v_call_light', 'j_call_heavy', 'j_call_light', 'c_call_heavy', 'c_call_light']]
                else:
                    self.metadata = self.metadata[[str(clonekey), str(clonekey)+'_by_size', 'isotype', 'lightchain', 'productive', 'status', 'vdj_status', 'umi_counts_heavy', 'umi_counts_light',  'v_call_heavy', 'v_call_light', 'j_call_heavy', 'j_call_light', 'c_call_heavy', 'c_call_light']]
            else:
                if 'sample_id' in self.data.columns:
                    self.metadata = self.metadata[['sample_id', 'isotype', 'lightchain', 'status', 'vdj_status', 'productive',  'umi_counts_heavy', 'umi_counts_light', 'v_call_heavy', 'v_call_light', 'j_call_heavy', 'j_call_light', 'c_call_heavy', 'c_call_light']]
                else:
                    self.metadata = self.metadata[['isotype', 'lightchain', 'productive', 'status', 'vdj_status', 'umi_counts_heavy', 'umi_counts_light',  'v_call_heavy', 'v_call_light', 'j_call_heavy', 'j_call_light', 'c_call_heavy', 'c_call_light']]
        # new function to retrieve non-standard columns
        if retrieve is not None:
            if type(retrieve) is str:
                retrieve = [retrieve]
            for ret in retrieve:
                if ret in dat.columns:
                    if not split_heavy_light:
                        if collapse:
                            retrieve_dict = retrieve_metadata(dat, ret, False, True)
                        else:
                            retrieve_dict = retrieve_metadata(dat, ret, False, False)
                        self.metadata[str(ret)] = pd.Series(retrieve_dict)
                    else:
                        if collapse:
                            h_retrieve_dict, l_retrieve_dict = retrieve_metadata(dat, ret, True, True)
                        else:
                            h_retrieve_dict, l_retrieve_dict = retrieve_metadata(dat, ret, True, False)
                        self.metadata[str(ret)+'_heavy'] = pd.Series(h_retrieve_dict)
                        self.metadata[str(ret)+'_light'] = pd.Series(l_retrieve_dict)
                else:
                    raise KeyError('Unknown column : \'%s\' to retrieve.' % ret)

class Dandelion:
    """
    Dandelion class object.

    Main class object storing input/ouput slots for all functions.

    """
    def __init__(self, data=None, metadata=None, germline = None, distance=None, edges=None, layout=None, graph=None, initialize = True, **kwargs):
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
                    update_metadata(self, **kwargs)
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
        return self._gen_repr(self.n_obs, self.n_contigs)

    def copy(self):
        """
        Performs a deep copy of all slots in `dandelion` class.

        Parameters
        ----------
        self : Dandelion
            Dandelion object.        
        
        Returns
        -------
        a deep copy of `dandelion` class.
        """
        return copy.deepcopy(self)

    def update_germline(self, corrected = None, germline = None, org = 'human'):
        """
        Update germline reference with corrected sequences and store in Dandelion object.

        Parameters
        ----------
        self : Dandelion
            Dandelion object.
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
                raise OSError('Environmental variable GERMLINE must be set. Otherwise, please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files.')
            gml = gml+'imgt/'+org+'/vdj/'
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
                            raise OSError('Input for germline is incorrect. Please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files, or individual paths to the germline IGHV, IGHD, and IGHJ fasta files (with .fasta extension) as a list.')
                        gml.append(x)
            elif type(germline) is list:
                if len(germline) < 3:
                    raise OSError('Input for germline is incorrect. Please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files, or individual paths to the germline IGHV, IGHD, and IGHJ fasta files (with .fasta extension) as a list.')
                else:
                    gml = []
                    for x in germline:
                        if not x.endswith('.fasta'):
                            raise OSError('Input for germline is incorrect. Please provide path to folder containing germline IGHV, IGHD, and IGHJ fasta files, or individual paths to the germline IGHV, IGHD, and IGHJ fasta files (with .fasta extension) as a list.')
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
                raise OSError('Input for corrected germline fasta is incorrect. Please provide path to file containing corrected germline fasta sequences.')

        self.germline.update(germline_ref)
        logg.info(' finished', time=start,
        deep=('Updated Dandelion object: \n'
        '   \'germline\', updated germline reference\n'))

    def write_pkl(self, filename='dandelion_data.pkl.pbz2', **kwargs):
        """
        Writes a `Dandelion` class to .pkl format.

        Parameters
        ----------
        filename
            path to Dandelion `.pkl` file.
        **kwargs
            passed to `_pickle`.
        """
        if isBZIP(filename):
            try:
                with bz2.BZ2File(filename, 'wb') as f:
                    cPickle.dump(self, f, **kwargs)
            except:
                with bz2.BZ2File(filename, 'wb') as f:
                    cPickle.dump(self, f, protocol = 4, **kwargs)
        elif isGZIP(filename):
            try:
                with gzip.open(filename, 'wb') as f:
                    cPickle.dump(self, f, **kwargs)
            except:
                with gzip.open(filename, 'wb') as f:
                    cPickle.dump(self, f, protocol = 4, **kwargs)
        else:
            f = open(filename, 'wb')
            cPickle.dump(self, f, **kwargs)
            f.close()

    def write_h5(self, filename='dandelion_data.h5', complib = None, compression = None, compression_level = None, **kwargs):
        """
        Writes a `Dandelion` class to .h5 format.

        Parameters
        ----------
        filename
            path to Dandelion `.h5` file.
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
        with h5py.File(filename,  "w") as hf:
            for datasetname in hf.keys():
                del hf[datasetname]

        if complib is None and compression is None:
            comp = None
        elif complib is not None and compression is None:
            comp = complib
        elif complib is None and compression is not None:
            comp = compression
        if complib is not None and compression is not None:
            raise ValueError('Please specify only complib or compression. They do the same thing.')

        # now to actually saving
        for col in self.data.columns:
            weird = (self.data[[col]].applymap(type) != self.data[[col]].iloc[0].apply(type)).any(axis=1)
            if len(self.data[weird]) > 0:
                self.data[col] = self.data[col].where(pd.notnull(self.data[col]), '')
        self.data.to_hdf(filename, "data", complib = comp, complevel = compression_level, **kwargs)
            
        if self.metadata is not None:
            for col in self.metadata.columns:                
                weird = (self.metadata[[col]].applymap(type) != self.metadata[[col]].iloc[0].apply(type)).any(axis=1)
                if len(self.metadata[weird]) > 0:
                    self.metadata[col] = self.metadata[col].where(pd.notnull(self.metadata[col]), '')
            self.metadata.to_hdf(filename, "metadata", complib = comp, complevel = compression_level, format='table', nan_rep=np.nan, **kwargs)
        # except:
        #     warnings.warn("`metadata` slot not saved. Please check if there is incompatible dtypes in the metadata table.")
        #     pass
        try:
            if 'index' in self.edges.columns:
                self.edges.drop('index', axis = 1, inplace=True)
            self.edges.to_hdf(filename, "edges", complib = comp, complevel = compression_level, **kwargs)
        except:
            pass

        graph_counter = 0
        try:
            for g in self.graph:
                G = nx.to_pandas_adjacency(g, nonedge=np.nan)
                G.to_hdf(filename, "graph/graph_"+str(graph_counter), complib = comp, complevel = compression_level, **kwargs)
                graph_counter += 1
        except:
            pass

        try:
            for d in self.distance:
                dat = pd.DataFrame(self.distance[d].toarray()) # how to make this faster?
                dat.to_hdf(filename, "distance/"+d, complib = comp, complevel = compression_level, **kwargs)
        except:
            pass

        with h5py.File(filename,  "a") as hf:
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
                        hf.create_group('layout/layout_'+str(layout_counter))
                    except:
                        pass
                    for k in l.keys():
                        hf['layout/layout_'+str(layout_counter)].attrs[k] = l[k]
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

def isGZIP(filename):
    if filename.split('.')[-1] == 'gz':
        return True
    return False

def isBZIP(filename):
    if filename.split('.')[-1] == 'pbz2':
        return True
    return False

def read_pkl(filename='dandelion_data.pkl.pbz2'):
    """
    Reads in and returns a `Dandelion` class saved using pickle format.

    Parameters
    ----------
    filename
        path to Dandelion `.pkl` file. Depending on the extension, it will try to unzip accordingly.

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

def read_h5(filename='dandelion_data.h5'):
    """
    Reads in and returns a `Dandelion` class from .h5 format.

    Parameters
    ----------
    filename
        path to Dandelion `.h5` file

    Returns
    -------
    Dandelion object.
    """    
    try:
        data = pd.read_hdf(filename, 'data')
    except:
        raise AttributeError('{} does not contain attribute `data`'.format(filename))
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
        for u,v,d in graph0.edges(data=True):
            d['weight'] = d['weight']-1
        for u,v,d in graph1.edges(data=True):
            d['weight'] = d['weight']-1
        graph = (graph0, graph1)
    except:
        pass

    with h5py.File(filename, 'r') as hf:
        try:
            layout0 = {}
            for k in hf['layout/layout_0'].attrs.keys():
                layout0.update({k:np.array(hf['layout/layout_0'].attrs[k])})
            layout1 = {}
            for k in hf['layout/layout_1'].attrs.keys():
                layout1.update({k:np.array(hf['layout/layout_1'].attrs[k])})
            layout = (layout0, layout1)
        except:
            pass
    
        germline = {}
        try:
            for g in hf['germline'].attrs:
                germline.update({g:hf['germline'].attrs[g]})
        except:
            pass
    
        distance = Tree()
        try:
            for d in hf['distance'].keys():
                d_ = pd.read_hdf(filename, 'distance/'+d)
                distance[d] = scipy.sparse.csr_matrix(d_.values)
                # d_ = hf['distance'][d]
                # distance[d] = scipy.sparse.csr_matrix((d_['data'][:],d_['indices'][:], d_['indptr'][:]), d_.attrs['shape'])
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
        res = Dandelion(**constructor, initialize = False)

    if 'threshold' in locals():
        res.threshold = threshold
    else:
        pass
    return(res)

def concat(arrays, check_unique = False):
    """
    Concatenate dataframe and return as Dandelion object.

    Parameters
    ----------
    arrays
        pandas dataframe or file path
    columns
        column names in dataframe
    
    Returns
    -------
    sample_dict
        dictionary
    """
    arrays = list(arrays)
    arrays_ = [x.data for x in arrays]
    if check_unique:
        df = pd.concat(arrays_, verify_integrity=True)
    else:
        df = pd.concat(arrays_)
    try:
        out = Dandelion(df)
    except:
        out = Dandelion(df, initialize = False)
    return(out)


# def convert_preprocessed_tcr_10x(file, prefix = None):
#     """
#     Parameters
#     ----------
#     file : str
#         file path to .tsv file.
#     prefix : str
#         prefix to add to barcodes. Ignored if left as None.    
#     """
#     from tqdm import tqdm
#     import re
#     cr_annot = pd.read_csv(file)
#     if prefix is not None:
#         cr_annot['index']=[prefix+'_'+i for i in cr_annot['contig_id']]
#     else:
#         cr_annot['index']=[i for i in cr_annot['contig_id']]
#     cr_annot.set_index('index', inplace = True)

#     ddl_annot = load_data("{}/tmp/{}_igblast_db-pass.tsv".format(os.path.dirname(file), os.path.basename(file).split('_annotations.csv')[0]))

#     for i in tqdm(ddl_annot.index, desc = 'Processing data '):
#         v = ddl_annot.loc[i, 'v_call']
#         d = ddl_annot.loc[i, 'd_call']
#         j = ddl_annot.loc[i, 'j_call']
#         if v is not np.nan:
#             ddl_annot.at[i, 'v_call_igblast'] = ','.join(list(set(re.sub('[*][0-9][0-9]', '', v).split(','))))
#             v_ = list(set(re.sub('[*][0-9][0-9]', '', v).split(',')))
#             if re.match('IGH', ','.join(v_)):
#                 ddl_annot.at[i, 'locus'] = 'IGH'
#             if re.match('IGK', ','.join(v_)):
#                 ddl_annot.at[i, 'locus'] = 'IGK'
#             if re.match('IGL', ','.join(v_)):
#                 ddl_annot.at[i, 'locus'] = 'IGL'
#             if len(v_) > 1:
#                 ddl_annot.at[i, 'locus'] = 'Multi'
#         if d is not np.nan:
#             ddl_annot.at[i, 'd_call_igblast'] = ','.join(list(set(re.sub('[*][0-9][0-9]', '', d).split(','))))
#             d_ = list(set(re.sub('[*][0-9][0-9]', '', d).split(',')))
#             if re.match('IGH', ','.join(d_)):
#                 ddl_annot.at[i, 'locus'] = 'IGH'
#             if re.match('IGK', ','.join(d_)):
#                 ddl_annot.at[i, 'locus'] = 'IGK'
#             if re.match('IGL', ','.join(d_)):
#                 ddl_annot.at[i, 'locus'] = 'IGL'
#             if len(d_) > 1:
#                 ddl_annot.at[i, 'locus'] = 'Multi'
#         if j is not np.nan:
#             ddl_annot.at[i, 'j_call_igblast'] = ','.join(list(set(re.sub('[*][0-9][0-9]', '', j).split(','))))
#             j_ = list(set(re.sub('[*][0-9][0-9]', '', j).split(',')))
#             if re.match('IGH', ','.join(j_)):
#                 ddl_annot.at[i, 'locus'] = 'IGH'
#             if re.match('IGK', ','.join(j_)):
#                 ddl_annot.at[i, 'locus'] = 'IGK'
#             if re.match('IGL', ','.join(j_)):
#                 ddl_annot.at[i, 'locus'] = 'IGL'
#             if len(j_) > 1:
#                 ddl_annot.at[i, 'locus'] = 'Multi'

#     cellrangermap = {
#         'sequence_id':'contig_id',
#         'locus':'chain',
#         'v_call_igblast':'v_gene',
#         'd_call_igblast':'d_gene',
#         'j_call_igblast':'j_gene',
#         'productive':'productive',
#         'junction_aa':'cdr3',
#         'junction':'cdr3_nt'}

#     for i in tqdm(cr_annot.index, desc = 'Matching and updating contig ids'):
#         if i in ddl_annot.index:
#             for key, value in cellrangermap.items():
#                 if cr_annot.loc[i, 'chain'] not in ['IGH', 'IGK', 'IGL', None]:
#                     cr_annot.at[i, value] = ddl_annot.loc[i, key]
#                 else:
#                     cr_annot.at[i, 'contig_id'] = ddl_annot.loc[i, 'sequence_id']

#             if cr_annot.loc[i, 'cdr3'] is np.nan:
#                 cr_annot.at[i, 'productive'] = None
#             else:
#                 if cr_annot.loc[i, 'productive'] == 'T':
#                     cr_annot.at[i, 'productive'] = 'True'
#                 else:
#                     cr_annot.at[i, 'productive'] = 'False'
#     cr_annot.to_csv("{}/{}_annotations_reannotated.csv".format(os.path.dirname(file), os.path.basename(file).split('_annotations.csv')[0]), index = False, na_rep = 'None')

