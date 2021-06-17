#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-05-12 17:56:02
# @Last Modified by:   Kelvin
# @Last Modified time: 2021-06-14 23:59:30

import os
import pandas as pd
import numpy as np
from subprocess import run
from datetime import timedelta
from anndata import AnnData
from time import time
from collections import OrderedDict
from time import time
from ...utilities._utilities import *
import scanpy as sc
import scipy.stats
import re
from os import PathLike
from typing import Union, Sequence, Tuple


def assigngenes_igblast(fasta: Union[str, PathLike], igblast_db: Union[None, str] = None, org: Literal['human', 'mouse'] = 'human', loci: Literal['ig', 'tr'] = 'ig', verbose: bool = False):
    """
    Reannotate with IgBLASTn.

    Parameters
    ----------
    fasta : PathLike
        fasta file for reannotation.
    igblast_db : PathLike, optional
        path to igblast database.
    org : str
        organism for germline sequences.
    loci : str
        `ig` or `tr` mode for running igblastn.
    verbose : bool
        whether or not to print the command used in terminal. Default is False.

    """
    env = os.environ.copy()
    if igblast_db is None:
        try:
            igdb = env['IGDATA']
        except:
            raise OSError(
                'Environmental variable IGDATA must be set. Otherwise, please provide path to igblast database')
    else:
        env['IGDATA'] = igblast_db
        igdb = env['IGDATA']

    outfolder = os.path.abspath(os.path.dirname(fasta)) + '/tmp'
    informat_dict = {'blast': '_igblast.fmt7', 'airr': '_igblast.tsv'}
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    for fileformat in ['blast', 'airr']:
        outfile = os.path.basename(fasta).split(
            '.fasta')[0] + informat_dict[fileformat]
        cmd = ['AssignGenes.py', 'igblast',
               '-s', fasta,
               '-b', igdb,
               '--organism', org,
               '--loci', loci,
               '--format', fileformat,
               '-o', "{}/{}".format(outfolder, outfile)
               ]

        if verbose:
            print('Running command: %s\n' % (' '.join(cmd)))
        run(cmd, env=env)  # logs are printed to terminal


def makedb_igblast(fasta: Union[str, PathLike], igblast_output: Union[None, str, PathLike] = None, germline: Union[None, str, PathLike] = None, org: Literal['human', 'mouse'] = 'human', extended: bool = True, verbose: bool = False):
    """
    Parses IgBLAST output to airr format.

    Parameters
    ----------
    fasta : PathLike
        fasta file use for reannotation.
    igblast_output : PathLike, optional
        igblast output file.
    germline : PathLike, optional
        path to germline database.    
    org : str
        organism of germline sequences.
    extended : bool
        whether or not to parse extended 10x annotations. Default is True.
    verbose : bool
        whether or not to print the command used in terminal. Default is False.

    """
    env = os.environ.copy()
    if germline is None:
        try:
            gml = env['GERMLINE']
        except:
            raise OSError(
                'Environmental variable GERMLINE must be set. Otherwise, please provide path to folder containing germline fasta files.')
        gml = gml + 'imgt/' + org + '/vdj/'
    else:
        env['GERMLINE'] = germline
        gml = germline

    if igblast_output is None:
        indir = os.path.dirname(fasta) + '/tmp'
        infile = os.path.basename(fasta).split('.fasta')[0] + '_igblast.fmt7'
        igbo = "{}/{}".format(indir, infile)
    else:
        igbo = igblast_output

    cellranger_annotation = "{}/{}".format(os.path.dirname(
        fasta), os.path.basename(fasta).replace('.fasta', '_annotations.csv'))

    if extended:
        cmd = ['MakeDb.py', 'igblast',
               '-i', igbo,
               '-s', fasta,
               '-r', gml,
               '--10x', cellranger_annotation,
               '--extended']
    else:
        cmd = ['MakeDb.py', 'igblast',
               '-i', igbo,
               '-s', fasta,
               '-r', gml,
               '--10x', cellranger_annotation]

    if verbose:
        print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd, env=env)  # logs are printed to terminal


def parsedb_heavy(db_file: Union[str, PathLike], verbose: bool = False):
    """
    Parses AIRR table (heavy chain contigs only).

    Parameters
    ----------
    db_file : PathLike
        path to AIRR table.
    verbose : bool
        whether or not to print the command used in terminal. Default is False.

    """
    outname = os.path.basename(db_file).split('.tsv')[0] + '_heavy'

    cmd = ['ParseDb.py', 'select',
           '-d', db_file,
           '-f', 'locus',
           '-u', 'IGH',
           '--logic', 'all',
           '--regex',
           '--outname', outname
           ]

    if verbose:
        print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd)  # logs are printed to terminal


def parsedb_light(db_file: Union[str, PathLike], verbose: bool = False):
    """
    Parses AIRR table (light chain contigs only).

    Parameters
    ----------
    db_file : PathLike
        path to AIRR table.
    verbose : bool
        whether or not to print the command used in terminal. Default is False.

    """
    outname = os.path.basename(db_file).split('.tsv')[0] + '_light'

    cmd = ['ParseDb.py', 'select',
           '-d', db_file,
           '-f', 'locus',
           '-u', 'IG[LK]',
           '--logic', 'all',
           '--regex',
           '--outname', outname
           ]

    if verbose:
        print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd)  # logs are printed to terminal


def creategermlines(db_file: Union[str, PathLike], germtypes: Union[None, str] = None, germline: Union[None, PathLike, str] = None, org: Literal['human', 'mouse'] = 'human', genotype_fasta: Union[None, PathLike, str] = None, v_field: Union[None, Literal['v_call', 'v_call_genotyped']] = None, cloned: bool = False, mode: Union[None, Literal['heavy', 'light']] = None, verbose: bool = False):
    """
    Wrapper for CreateGermlines.py for reconstructing germline sequences,

    Parameters
    ----------
    db_file : PathLike
        path to AIRR table.
    germtypes : str, optional
        germline type for reconstuction.
    germline : PathLike, optional
        location to germline fasta files.
    org : str
        organism for germline sequences.
    genotype_fasta : PathLike, optional
        location to corrected v germine fasta file.
    v_field : str, optional
        name of column for v segment to perform reconstruction.
    cloned : bool
        whether or not to run with cloned option.
    mode : str, optional
        whether to run on heavy or light mode. If left as None, heavy and light will be run together.
    verbose : bool
        whether or not to print the command used in terminal. Default is False.

    """
    env = os.environ.copy()
    if germline is None:
        try:
            gml = env['GERMLINE']
        except:
            raise OSError(
                'Environmental variable GERMLINE must be set. Otherwise, please provide path to folder containing germline fasta files.')
        gml = gml + 'imgt/' + org + '/vdj/'
    else:
        env['GERMLINE'] = germline
        gml = germline

    if germtypes is None:
        germ_type = 'dmask'
    else:
        germ_type = germtypes

    if cloned:
        if mode == 'heavy':
            print('            Reconstructing heavy chain {} germline sequences with {} for each clone.'.format(
                germ_type, v_field))
            if genotype_fasta is None:
                if germline is None:
                    cmd = ['CreateGermlines.py',
                           '-d', db_file,
                           '-g', germ_type,
                           '--cloned',
                           '-r', gml + '/imgt_' + org + '_IGHV.fasta', gml + '/imgt_' +
                           org + '_IGHD.fasta', gml + '/imgt_' + org + '_IGHJ.fasta',
                           '--vf', v_field
                           ]
                else:
                    cmd = ['CreateGermlines.py',
                           '-d', db_file,
                           '-g', germ_type,
                           '--cloned',
                           '-r', gml,
                           '--vf', v_field
                           ]
            else:
                if germline is None:
                    cmd = ['CreateGermlines.py',
                           '-d', db_file,
                           '-g', germ_type,
                           '--cloned',
                           '-r', genotype_fasta, gml + '/imgt_' + org +
                           '_IGHD.fasta', gml + '/imgt_' + org + '_IGHJ.fasta',
                           '--vf', v_field
                           ]
                else:
                    cmd = ['CreateGermlines.py',
                           '-d', db_file,
                           '-g', germ_type,
                           '--cloned',
                           '-r', genotype_fasta, gml,
                           '--vf', v_field
                           ]
        elif mode == 'light':
            print('            Reconstructing light chain {} germline sequences with {} for each clone.'.format(
                germ_type, v_field))
            if germline is None:
                cmd = ['CreateGermlines.py',
                       '-d', db_file,
                       '-g', germ_type,
                       '--cloned',
                       '-r', gml + '/imgt_' + org + '_IGKV.fasta', gml + '/imgt_' + org + '_IGKJ.fasta', gml +
                       '/imgt_' + org + '_IGLV.fasta', gml + '/imgt_' + org + '_IGLJ.fasta',
                       '--vf', v_field
                       ]
            else:
                cmd = ['CreateGermlines.py',
                       '-d', db_file,
                       '-g', germ_type,
                       '--cloned',
                       '-r', gml,
                       '--vf', v_field
                       ]
        elif mode is None:
            print('            Reconstructing {} germline sequences with {} for each clone.'.format(
                germ_type, v_field))
            if genotype_fasta is None:
                if germline is None:
                    cmd = ['CreateGermlines.py',
                           '-d', db_file,
                           '-g', germ_type,
                           '--cloned',
                           '-r', gml + '/imgt_' + org + '_IGHV.fasta', gml + '/imgt_' + org + '_IGHD.fasta', gml + '/imgt_' + org + '_IGHJ.fasta', gml + '/imgt_' +
                           org + '_IGKV.fasta', gml + '/imgt_' + org + '_IGKJ.fasta', gml +
                           '/imgt_' + org + '_IGLV.fasta', gml + '/imgt_' + org + '_IGLJ.fasta',
                           '--vf', v_field
                           ]
                else:
                    cmd = ['CreateGermlines.py',
                           '-d', db_file,
                           '-g', germ_type,
                           '--cloned',
                           '-r', gml,
                           '--vf', v_field
                           ]
            else:
                if germline is None:
                    cmd = ['CreateGermlines.py',
                           '-d', db_file,
                           '-g', germ_type,
                           '--cloned',
                           '-r', genotype_fasta, gml + '/imgt_' + org + '_IGHD.fasta', gml + '/imgt_' + org + '_IGHJ.fasta', gml + '/imgt_' + org +
                           '_IGKV.fasta', gml + '/imgt_' + org + '_IGKJ.fasta', gml + '/imgt_' +
                           org + '_IGLV.fasta', gml + '/imgt_' + org + '_IGLJ.fasta',
                           '--vf', v_field
                           ]
                else:
                    cmd = ['CreateGermlines.py',
                           '-d', db_file,
                           '-g', germ_type,
                           '--cloned',
                           '-r', genotype_fasta, gml,
                           '--vf', v_field
                           ]
    else:
        if mode == 'heavy':
            print('            Reconstructing heavy chain {} germline sequences with {}.'.format(
                germ_type, v_field))
            if genotype_fasta is None:
                if germline is None:
                    cmd = ['CreateGermlines.py',
                           '-d', db_file,
                           '-g', germ_type,
                           '-r', gml + '/imgt_' + org + '_IGHV.fasta', gml + '/imgt_' +
                           org + '_IGHD.fasta', gml + '/imgt_' + org + '_IGHJ.fasta',
                           '--vf', v_field
                           ]
                else:
                    cmd = ['CreateGermlines.py',
                           '-d', db_file,
                           '-g', germ_type,
                           '-r', gml,
                           '--vf', v_field
                           ]
            else:
                if germline is None:
                    cmd = ['CreateGermlines.py',
                           '-d', db_file,
                           '-g', germ_type,
                           '-r', genotype_fasta, gml + '/imgt_' + org +
                           '_IGHD.fasta', gml + '/imgt_' + org + '_IGHJ.fasta',
                           '--vf', v_field
                           ]
                else:
                    cmd = ['CreateGermlines.py',
                           '-d', db_file,
                           '-g', germ_type,
                           '-r', genotype_fasta, gml,
                           '--vf', v_field
                           ]
        elif mode == 'light':
            print('            Reconstructing light chain {} germline sequences with {}.'.format(
                germ_type, v_field))
            if germline is None:
                cmd = ['CreateGermlines.py',
                       '-d', db_file,
                       '-g', germ_type,
                       '-r', gml + '/imgt_' + org + '_IGKV.fasta', gml + '/imgt_' + org + '_IGKJ.fasta', gml +
                       '/imgt_' + org + '_IGLV.fasta', gml + '/imgt_' + org + '_IGLJ.fasta',
                       '--vf', v_field
                       ]
            else:
                cmd = ['CreateGermlines.py',
                       '-d', db_file,
                       '-g', germ_type,
                       '-r', gml,
                       '--vf', v_field
                       ]
        elif mode is None:
            print('            Reconstructing {} germline sequences with {} for each clone.'.format(
                germ_type, v_field))
            if genotype_fasta is None:
                if germline is None:
                    cmd = ['CreateGermlines.py',
                           '-d', db_file,
                           '-g', germ_type,
                           '-r', gml + '/imgt_' + org + '_IGHV.fasta', gml + '/imgt_' + org + '_IGHD.fasta', gml + '/imgt_' + org + '_IGHJ.fasta', gml + '/imgt_' +
                           org + '_IGKV.fasta', gml + '/imgt_' + org + '_IGKJ.fasta', gml +
                           '/imgt_' + org + '_IGLV.fasta', gml + '/imgt_' + org + '_IGLJ.fasta',
                           '--vf', v_field
                           ]
                else:
                    cmd = ['CreateGermlines.py',
                           '-d', db_file,
                           '-g', germ_type,
                           '-r', gml,
                           '--vf', v_field
                           ]
            else:
                if germline is None:
                    cmd = ['CreateGermlines.py',
                           '-d', db_file,
                           '-g', germ_type,
                           '-r', genotype_fasta, gml + '/imgt_' + org + '_IGHD.fasta', gml + '/imgt_' + org + '_IGHJ.fasta', gml + '/imgt_' + org +
                           '_IGKV.fasta', gml + '/imgt_' + org + '_IGKJ.fasta', gml + '/imgt_' +
                           org + '_IGLV.fasta', gml + '/imgt_' + org + '_IGLJ.fasta',
                           '--vf', v_field
                           ]
                else:
                    cmd = ['CreateGermlines.py',
                           '-d', db_file,
                           '-g', germ_type,
                           '-r', genotype_fasta, gml,
                           '--vf', v_field
                           ]

    if verbose:
        print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd, env=env)  # logs are printed to terminal


def tigger_genotype(data: Union[str, PathLike], v_germline: Union[None, PathLike, str] = None, outdir: Union[None, PathLike, str] = None, org: Literal['human', 'mouse'] = 'human', fileformat: Literal['airr', 'changeo'] = 'airr', novel_: Literal['YES', 'NO'] = 'YES', verbose: bool = False):
    """
    Reassign alleles with TIgGER in R.

    Parameters
    ----------
    data : PathLike
        vdj tabulated data, in Change-O (TAB) or AIRR (TSV) format.
    germline : PathLike, optional
        fasta file containing IMGT-gapped V segment reference germlines. Defaults to $GERMLINE.
    outdir : PathLike,  optional
        output directory. Will be created if it does not exist. Defaults to the current working directory.
    org : str
        organism for germline sequences.
    fileformat : str
        format for running tigger. Default is 'airr'. Also accepts 'changeo'.
    novel : str
        whether or not to run novel allele discovery. Default is 'YES'.
    verbose : bool
        whether or not to print the command used in terminal. Default is False.

    """
    start_time = time()
    env = os.environ.copy()
    if v_germline is None:
        try:
            gml = env['GERMLINE']
        except:
            raise OSError('Environmental variable GERMLINE is not set. Please provide either the path to the folder containing the germline IGHV fasta file, or direct path to the germline IGHV fasta file.')
        gml = gml + 'imgt/' + org + '/vdj/imgt_' + org + '_IGHV.fasta'
    else:
        if os.path.isdir(v_germline):
            gml = v_germline.rstrip('/') + 'imgt_' + org + '_IGHV.fasta'
            if not os.path.isfile(gml):
                raise OSError("Input for germline is incorrect. Please rename IGHV germline file to '{}'. Otherwise, please provide path to folder containing the germline IGHV fasta file, or direct path to the germline IGHV fasta file.".format(gml))
        else:
            if not v_germline.endswith('.fasta'):
                raise OSError(
                    'Input for germline is incorrect {}. Please provide path to folder containing the germline IGHV fasta file, or direct path to the germline IGHV fasta file.'.format(v_germline))
            if (os.path.isfile(v_germline)) & ('ighv' in v_germline.lower()):
                gml = v_germline

    if outdir is not None:
        out_dir = outdir + '/'
    else:
        out_dir = os.path.dirname(data)

    cmd = ['tigger-genotype.R',
           '-d', data,
           '-r', gml,
           '-n', os.path.basename(data).split('.tsv')[0],
           '-N', novel_,
           '-o', out_dir,
           '-f', fileformat]

    print('      Reassigning alleles')
    if verbose:
        print('Running command: %s\n' % (' '.join(cmd)))
    run(cmd, env=env)  # logs are printed to terminal
    elapsed_time_secs = time() - start_time
    msg = "tigger-genotype execution took: %s secs (Wall clock time)\n" % timedelta(
        seconds=round(elapsed_time_secs))
    if verbose:
        print(msg)


def recipe_scanpy_qc(self: AnnData, mito_startswith: str = 'MT', max_genes: int = 2500, min_genes: int = 200, mito_cutoff: int = 5, pval_cutoff: float = 0.1, min_counts: Union[None, int] = None, max_counts: Union[None, int] = None, blacklist: Union[None, Sequence] = None, vdj_pattern: str = '^TR[AB][VDJ]|^IG[HKL][VDJC]') -> AnnData:
    """
    Recipe for running a standard scanpy QC workflow.

    Parameters
    ----------
    adata : AnnData
        annotated data matrix of shape n_obs × n_vars. Rows correspond to cells and columns to genes.
    mito_startswith : str
        string pattern used for searching mitochondrial genes.
    max_genes : int
        naximum number of genes expressed required for a cell to pass filtering. Default is 2500.
    min_genes : int
        minimum number of genes expressed  required for a cell to pass filtering. Default is 200.
    mito_cutoff : float
        maximum percentage mitochondrial content allowed for a cell to pass filtering. Default is 5.
    pval_cutoff : float
        maximum Benjamini-Hochberg corrected p value from doublet detection protocol allowed for a cell to pass filtering. Default is 0.05.
    min_counts : int, optional
        minimum number of counts required for a cell to pass filtering. Default is None.
    max_counts : int, optional
        maximum number of counts required for a cell to pass filtering. Default is None.
    blacklist : sequence, optional
        if provided, will exclude these genes from highly variable genes list.
    vdj_pattern : str
        string pattern for search VDJ genes to exclude from highly variable genes.

    Returns
    -------
    `AnnData` of shape n_obs × n_vars where obs now contain filtering information. Rows correspond to cells and columns to genes.
    """
    _adata = self.copy()
    # run scrublet
    try:
        import scrublet as scr
    except:
        raise ImportError('Please install scrublet with pip install scrublet.')

    scrub = scr.Scrublet(_adata.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
    _adata.obs['scrublet_score'] = doublet_scores
    # overcluster prep. run basic scanpy pipeline
    sc.pp.filter_cells(_adata, min_genes=0)
    _adata.var['mt'] = _adata.var_names.str.startswith(mito_startswith)
    sc.pp.calculate_qc_metrics(_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pp.normalize_total(_adata, target_sum=1e4)
    sc.pp.log1p(_adata)
    sc.pp.highly_variable_genes(
        _adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    for i in _adata.var.index:
        if re.search(vdj_pattern, i):
            _adata.var.at[i, 'highly_variable'] = False
        if blacklist is not None:
            if i in blacklist:
                _adata.var.at[i, 'highly_variable'] = False
    _adata = _adata[:, _adata.var['highly_variable']]
    sc.pp.scale(_adata, max_value=10)
    sc.tl.pca(_adata, svd_solver='arpack')
    sc.pp.neighbors(_adata, n_neighbors=10, n_pcs=50)
    # overclustering proper - do basic clustering first, then cluster each cluster
    sc.tl.leiden(_adata)
    for clus in list(np.unique(_adata.obs['leiden']))[0]:
        sc.tl.leiden(_adata, restrict_to=(
            'leiden', [clus]), key_added='leiden_R')
    # weird how the new anndata/scanpy is forcing this
    for clus in list(np.unique(_adata.obs['leiden']))[1:]:
        sc.tl.leiden(_adata, restrict_to=(
            'leiden_R', [clus]), key_added='leiden_R')
    # compute the cluster scores - the median of Scrublet scores per overclustered cluster
    for clus in np.unique(_adata.obs['leiden_R']):
        _adata.obs.loc[_adata.obs['leiden_R'] == clus, 'scrublet_cluster_score'] = \
            np.median(
                _adata.obs.loc[_adata.obs['leiden_R'] == clus, 'scrublet_score'])
    # now compute doublet p-values. figure out the median and mad (from above-median values) for the distribution
    med = np.median(_adata.obs['scrublet_cluster_score'])
    mask = _adata.obs['scrublet_cluster_score'] > med
    mad = np.median(_adata.obs['scrublet_cluster_score'][mask] - med)
    # 1 sided test for catching outliers
    pvals = 1 - \
        scipy.stats.norm.cdf(
            _adata.obs['scrublet_cluster_score'], loc=med, scale=1.4826 * mad)
    _adata.obs['scrublet_score_bh_pval'] = bh(pvals)
    # threshold the p-values to get doublet calls.
    _adata.obs['is_doublet'] = _adata.obs['scrublet_score_bh_pval'] < pval_cutoff
    _adata.obs['is_doublet'] = _adata.obs['is_doublet'].astype('category')
    _adata.obs['filter_rna'] = (pd.Series([min_genes < n > max_genes for n in _adata.obs['n_genes_by_counts']], index=_adata.obs.index)) | \
        (_adata.obs['pct_counts_mt'] >= mito_cutoff) | \
        (_adata.obs['is_doublet'] == True)

    # removing columns that probably don't need anymore
    _adata.obs = _adata.obs.drop(
        ['leiden', 'leiden_R', 'scrublet_cluster_score', 'scrublet_score_bh_pval'], axis=1)
    self.obs = _adata.obs.copy()
