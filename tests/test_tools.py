#!/usr/bin/env python
import pytest
import pandas as pd
import dandelion as ddl
from pathlib import Path

from fixtures import (airr_reannotated, dummy_adata, create_testfolder)


def test_setup(create_testfolder, airr_reannotated, dummy_adata):
    vdj, adata = ddl.pp.filter_contigs(airr_reannotated, dummy_adata)
    assert airr_reannotated.shape[0] == 9
    assert vdj.data.shape[0] == 7
    assert vdj.metadata.shape[0] == 4
    assert adata.n_obs == 5
    f = create_testfolder / "test.h5"
    vdj.write_h5(f)
    assert len(list(create_testfolder.iterdir())) == 1
    vdj2 = ddl.read_h5(f)
    assert vdj2.metadata is not None


def test_find_clones(create_testfolder):
    f = create_testfolder / "test.h5"
    vdj = ddl.read_h5(f)
    ddl.tl.find_clones(vdj)
    assert not vdj.data.clone_id.empty
    assert not vdj.metadata.clone_id.empty
    assert len(set(x for x in vdj.metadata['clone_id'] if pd.notnull(x))) == 4
    vdj.write_h5(f)


@pytest.mark.parametrize(
    "resample,expected",
    [pytest.param(None, 4), pytest.param(3, 3)])
def test_generate_network(create_testfolder, resample, expected):
    f = create_testfolder / "test.h5"
    vdj = ddl.read_h5(f)
    if resample is not None:
        vdj = ddl.tl.generate_network(vdj, downsample=resample)
    else:
        ddl.tl.generate_network(vdj)
    assert vdj.distance is not None
    assert vdj.edges is None
    assert vdj.n_obs == expected
    assert vdj.layout is not None
    assert vdj.graph is not None
    vdj.data['clone_id'] = '1'
    vdj = ddl.Dandelion(vdj.data)
    ddl.tl.generate_network(vdj)
    assert vdj.edges is not None


def test_transfer(create_testfolder, dummy_adata):
    f = create_testfolder / "test.h5"
    vdj = ddl.read_h5(f)
    ddl.tl.transfer(dummy_adata, vdj)
    assert 'clone_id' in dummy_adata.obs
    ddl.tl.generate_network(vdj)
    ddl.tl.transfer(dummy_adata, vdj)
    assert 'X_vdj' in dummy_adata.obsm


def test_diversity(create_testfolder):
    f = create_testfolder / "test.h5"
    vdj = ddl.read_h5(f)
    ddl.tl.clone_diversity(vdj, groupby = 'sample_id')
    assert not vdj.metadata.clone_network_vertex_size_gini.empty
    assert not vdj.metadata.clone_network_cluster_size_gini.empty
    ddl.tl.generate_network(vdj)
    ddl.tl.clone_diversity(vdj, groupby = 'sample_id', metric = 'clone_centrality')
    assert not vdj.metadata.clone_centrality_gini.empty
    assert not vdj.metadata.clone_size_gini.empty
