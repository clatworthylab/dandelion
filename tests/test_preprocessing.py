#!/usr/bin/env python
import pytest
import pandas as pd
import dandelion as ddl
from pathlib import Path

from fixtures import (fasta_10x, annotation_10x, create_testfolder,
                      database_paths, dummy_adata, processed_files)


@pytest.mark.parametrize("filename,expected",
                         [pytest.param('filtered', 1),
                          pytest.param('all', 2)])
def test_write_fasta(create_testfolder, fasta_10x, filename, expected):
    out_fasta = str(create_testfolder) + "/" + filename + "_contig.fasta"
    fh = open(out_fasta, "w")
    fh.close()
    out = ''
    for l in fasta_10x:
        out = '>' + l + '\n' + fasta_10x[l] + '\n'
        ddl.utl.Write_output(out, out_fasta)
    assert len(list(create_testfolder.iterdir())) == expected


@pytest.mark.parametrize("filename,expected",
                         [pytest.param('filtered', 3),
                          pytest.param('all', 4)])
def test_write_annotation(create_testfolder, annotation_10x, filename,
                          expected):
    out_file = str(
        create_testfolder) + "/" + filename + "_contig_annotations.csv"
    annotation_10x.to_csv(out_file, index=False)
    assert len(list(create_testfolder.iterdir())) == expected


@pytest.mark.parametrize("filename,expected", [
    pytest.param(None, 2),
    pytest.param('filtered', 2),
    pytest.param('all', 4)
])
def test_formatfasta(create_testfolder, filename, expected):
    ddl.pp.format_fastas(str(create_testfolder), filename_prefix=filename)
    assert len(list((create_testfolder / 'dandelion').iterdir())) == expected


@pytest.mark.parametrize("filename,expected",
                         [pytest.param('filtered', 3),
                          pytest.param('all', 6)])
def test_reannotategenes(create_testfolder, database_paths, filename,
                         expected):
    ddl.pp.reannotate_genes(str(create_testfolder),
                            igblast_db=database_paths['igblast_db'],
                            germline=database_paths['germline'],
                            filename_prefix=filename)
    assert len(list(
        (create_testfolder / 'dandelion/tmp').iterdir())) == expected


@pytest.mark.parametrize("filename,combine,expected", [
    pytest.param('filtered', 'reassigned_filtered', 9),
    pytest.param('all', 'reassigned_all', 12)
])
def test_reassignalleles(create_testfolder, database_paths, filename, combine,
                         expected):
    ddl.pp.reassign_alleles(str(create_testfolder),
                            combined_folder=combine,
                            germline=database_paths['germline'],
                            filename_prefix=filename,
                            novel=False,
                            plot=False)
    assert len(list(
        (create_testfolder / 'dandelion/tmp').iterdir())) == expected


def test_updateblastdb(database_paths):
    ddl.utl.makeblastdb(database_paths['blastdb_fasta'])
    assert len(list(Path(database_paths['blastdb']).iterdir())) == 10


@pytest.mark.parametrize("filename, expected",
                         [pytest.param('filtered', 6),
                          pytest.param('all', 7)])
def test_assignsisotypes(create_testfolder, database_paths, filename,
                         expected):
    ddl.pp.assign_isotypes(str(create_testfolder),
                           blastdb=database_paths['blastdb_fasta'],
                           filename_prefix=filename,
                           plot=False)
    assert len(list((create_testfolder / 'dandelion').iterdir())) == expected


@pytest.mark.parametrize("filename", ['all', 'filtered'])
def test_checkccall(create_testfolder, processed_files, filename):
    f = create_testfolder / str('dandelion/' + processed_files[filename])
    dat = pd.read_csv(f, sep='\t')
    assert not dat['c_call'].empty


def test_create_germlines(create_testfolder, processed_files, database_paths):
    f = create_testfolder / str('dandelion/' + processed_files['filtered'])
    ddl.pp.create_germlines(f, germline=database_paths['germline'])
    f2 = create_testfolder / str('dandelion/' +
                                 processed_files['germline-dmask'])
    dat = pd.read_csv(f2, sep='\t')
    assert not dat['germline_alignment_d_mask'].empty


@pytest.mark.parametrize(
    "freq,colname",
    [pytest.param(True, 'mu_freq'),
     pytest.param(False, 'mu_count')])
def test_quantify_mut(create_testfolder, processed_files, freq, colname):
    f = create_testfolder / str('dandelion/' + processed_files['filtered'])
    ddl.pp.quantify_mutations(f, frequency=freq)
    dat = pd.read_csv(f, sep='\t')
    assert not dat[colname].empty
    assert dat[colname].dtype == float


@pytest.mark.parametrize("filename", ['all', 'filtered'])
def test_filtercontigs(create_testfolder, processed_files, dummy_adata,
                       filename):
    f = create_testfolder / str('dandelion/' + processed_files[filename])
    dat = pd.read_csv(f, sep='\t')
    vdj, adata = ddl.pp.filter_contigs(dat, dummy_adata)
    assert dat.shape[0] == 9
    assert vdj.data.shape[0] == 7
    assert vdj.metadata.shape[0] == 4
    assert adata.n_obs == 5
