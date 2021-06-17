#!/usr/bin/env python
import dandelion as ddl

from fixtures import airr_10x, create_testfolder


def test_write_airr(create_testfolder, airr_10x):
    out_file = str(create_testfolder) + "/airr_rearrangements.tsv"
    airr_10x.to_csv(out_file, sep='\t', index=False)
    assert len(list(create_testfolder.iterdir())) == 1


def test_read10xairr(create_testfolder):
    airr_file = str(create_testfolder) + "/airr_rearrangements.tsv"
    vdj = ddl.read_10x_airr(airr_file)
    assert vdj.data.shape[0] == 9
    assert vdj.data.shape[1] == 33
    assert vdj.metadata.shape[0] == 4
    assert vdj.metadata.shape[1] == 23
