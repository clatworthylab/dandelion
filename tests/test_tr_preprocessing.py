#!/usr/bin/env python
import pytest
import pandas as pd
import dandelion as ddl
from pathlib import Path

from fixtures import (fasta_10x_tr1, fasta_10x_tr2, annotation_10x_tr1,
                      processed_files_tr, annotation_10x_tr2,
                      create_testfolder, database_paths, dummy_adata_tr)


def test_write_fasta_tr1(create_testfolder, fasta_10x_tr1):
    out_fasta = str(create_testfolder) + "/filtered_contig.fasta"
    fh = open(out_fasta, "w")
    fh.close()
    out = ''
    for l in fasta_10x_tr1:
        out = '>' + l + '\n' + fasta_10x_tr1[l] + '\n'
        ddl.utl.Write_output(out, out_fasta)
    assert len(list(create_testfolder.iterdir())) == 1


def test_write_fasta_tr2(create_testfolder, fasta_10x_tr2):
    out_fasta = str(create_testfolder) + "/all_contig.fasta"
    fh = open(out_fasta, "w")
    fh.close()
    out = ''
    for l in fasta_10x_tr2:
        out = '>' + l + '\n' + fasta_10x_tr2[l] + '\n'
        ddl.utl.Write_output(out, out_fasta)
    assert len(list(create_testfolder.iterdir())) == 2


def test_write_annotation_tr1(create_testfolder, annotation_10x_tr1):
    out_file = str(create_testfolder) + "/filtered_contig_annotations.csv"
    annotation_10x_tr1.to_csv(out_file, index=False)
    assert len(list(create_testfolder.iterdir())) == 3


def test_write_annotation_tr2(create_testfolder, annotation_10x_tr2):
    out_file = str(create_testfolder) + "/all_contig_annotations.csv"
    annotation_10x_tr2.to_csv(out_file, index=False)
    assert len(list(create_testfolder.iterdir())) == 4


@pytest.mark.parametrize("filename,expected", [
    pytest.param(None, 2),
    pytest.param('filtered', 2),
    pytest.param('all', 4)
])
def test_formatfasta(create_testfolder, filename, expected):
    ddl.pp.format_fastas(str(create_testfolder), filename_prefix=filename)
    assert len(list((create_testfolder / 'dandelion').iterdir())) == expected


@pytest.mark.parametrize(
    "filename,expected",
    [pytest.param('filtered', [3, 6]),
     pytest.param('all', [6, 7])])
def test_reannotategenes(create_testfolder, database_paths, filename,
                         expected):
    ddl.pp.reannotate_genes(str(create_testfolder),
                            igblast_db=database_paths['igblast_db'],
                            germline=database_paths['germline'],
                            loci='tr',
                            filename_prefix=filename)
    assert len(list(
        (create_testfolder / 'dandelion/tmp').iterdir())) == expected[0]
    assert len(list(
        (create_testfolder / 'dandelion').iterdir())) == expected[1]


@pytest.mark.parametrize("filename", ['filtered', 'all'])
def test_checkreannotation(create_testfolder, processed_files_tr, filename):
    f1 = pd.read_csv(create_testfolder /
                     str(filename + '_contig_annotations.csv'))
    f2 = pd.read_csv(create_testfolder /
                     str('dandelion/' + processed_files_tr[filename]),
                     sep='\t')
    f1 = f1.replace({'None': None, '': None, 'False': False, 'True': True})
    f2 = f2.replace({'None': None, '': None, 'False': False, 'True': True})
    assert all(pd.isnull(f1.d_gene))
    assert not f2.d_call.empty
    assert all(pd.notnull(f2.d_call))
    assert not all(f1.productive)
    assert all(f2.productive)


@pytest.mark.parametrize("filename,expected",
                         [pytest.param('filtered', 2),
                          pytest.param('all', 3)])
def test_filtercontigs(create_testfolder, processed_files_tr, dummy_adata_tr,
                       filename, expected):
    f = create_testfolder / str('dandelion/' + processed_files_tr[filename])
    dat = pd.read_csv(f, sep='\t')
    vdj, adata = ddl.pp.filter_contigs(dat, dummy_adata_tr, locus='tr-gd')
    assert dat.shape[0] == expected
    assert vdj.data.shape[0] == expected
    assert vdj.metadata.shape[0] == expected
    assert adata.n_obs == 3
