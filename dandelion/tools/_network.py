#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2020-08-12 18:08:04
# @Last Modified by:   Kelvin
# @Last Modified time: 2021-05-20 17:28:38

import pandas as pd
import numpy as np
import networkx as nx
from polyleven import levenshtein
from ..utilities._utilities import *
from ..utilities._core import *
from ..utilities._io import *
from networkx.utils import random_state
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import pdist, squareform
from itertools import combinations
from tqdm import tqdm
from time import sleep
try:
    from scanpy import logging as logg
except ImportError:
    pass
from typing import Union, Sequence, Tuple


def generate_network(self: Union[Dandelion, pd.DataFrame, str], key: Union[None, str] = None, clone_key: Union[None, str] = None, min_size: int = 2, downsample: Union[None, int] = None, verbose: bool = True, **kwargs)->Dandelion:
    """
    Generates a Levenshtein distance network based on full length VDJ sequence alignments for heavy and light chain(s).
    The distance matrices are then combined into a singular matrix.

    Parameters
    ----------
    data : Dandelion, DataFrame, str
        `Dandelion` object, pandas `DataFrame` in changeo/airr format, or file path to changeo/airr file after clones have been determined.
    key : str, optional
        column name for distance calulations. None defaults to 'sequence_alignment_aa'.
    clone_key: str, optional
        column name to build network on.
    min_size : int
        For visualization purposes, two graphs are created where one contains all cells and a trimmed second graph. This value specifies the minimum number of edges required otherwise node will be trimmed in the secondary graph.
    downsample : int, optional
        whether or not to downsample the number of cells prior to construction of network. If provided, cells will be randomly sampled to the integer provided. A new Dandelion class will be returned.
    verbose : bool
        whether or not to print the progress bars.
    **kwargs
        additional kwargs passed to options specified in `networkx.drawing.layout.spring_layout`.

    Returns
    -------
    `Dandelion` object with `.distance`, `.edges`, `.layout`, `.graph` initialized.
    """
    if verbose:
        start = logg.info('Generating network')
    if self.__class__ == Dandelion:
        dat = load_data(self.data)
    else:
        dat = load_data(self)

    if key is None:
        key_ = 'sequence_alignment_aa'  # default
    else:
        key_ = key

    if key_ not in dat:
        raise ValueError("key {} not found in input table.".format(key_))

    if clone_key is None:
        clonekey = 'clone_id'
    else:
        clonekey = clone_key
    if clonekey not in dat:
        raise ValueError(
            'Data does not contain clone information. Please run find_clones.')

    # calculate distance

    if downsample is not None:
        # if downsample >= dat_h.shape[0]:
        if downsample >= self.metadata.shape[0]:
            if verbose:
                print('Cannot downsample to {} cells. Using all {} cells.'.format(
                    str(downsample), self.metadata.shape[0]))
        else:
            if verbose:
                print('Downsampling to {} cells.'.format(str(downsample)))
            dat_h = dat[dat['locus'].isin(['IGH', 'TRB', 'TRD'])]
            dat_l = dat[dat['locus'].isin(['IGK', 'IGL', 'TRA', 'TRG'])]
            dat_h = dat_h.sample(downsample)
            dat_l = dat_l[dat_l['cell_id'].isin(list(dat_h['cell_id']))]
            dat_ = dat_h.append(dat_l)
    else:
        dat_ = dat.copy()

    # So first, create a data frame to hold all possible (full) sequences split by heavy (only 1 possible for now) and light (multiple possible)
    try:
        dat_seq = retrieve_metadata(dat_, query=key_, split=True, collapse=False)
    except:
        dat_seq = retrieve_metadata(dat_, query=key_, split=True, collapse=False, **kwargs)
    dat_seq.columns = [re.sub(key_ + '_', '', i) for i in dat_seq.columns]

    # calculate a distance matrix for all vs all and this can be referenced later on to extract the distance between the right pairs
    dmat = Tree()
    sleep(0.5)
    if verbose:
        for x in tqdm(dat_seq.columns, desc='Calculating distances... '):
            tdarray = np.array(np.array(dat_seq[x])).reshape(-1, 1)
            # d_mat = squareform([levenshtein(x[0],y[0]) for x,y in combinations(tdarray, 2) if (x[0] == x[0]) and (y[0] == y[0]) else 0])
            d_mat = squareform(pdist(tdarray, lambda x, y: levenshtein(
                x[0], y[0]) if (x[0] == x[0]) and (y[0] == y[0]) else 0))
            dmat[x] = d_mat
    else:
        for x in dat_seq.columns:
            tdarray = np.array(np.array(dat_seq[x])).reshape(-1, 1)
            # d_mat = squareform([levenshtein(x[0],y[0]) for x,y in combinations(tdarray, 2) if (x[0] == x[0]) and (y[0] == y[0]) else 0])
            d_mat = squareform(pdist(tdarray, lambda x, y: levenshtein(
                x[0], y[0]) if (x[0] == x[0]) and (y[0] == y[0]) else 0))
            dmat[x] = d_mat

    dist_mat_list = [dmat[x] for x in dmat if type(dmat[x]) is np.ndarray]

    total_dist = np.sum(dist_mat_list, axis=0)

    # generate edge list
    if self.__class__ == Dandelion:
        out = self.copy()
        if downsample is not None:
            out = Dandelion(dat_)
    else:  # re-initiate a Dandelion class object
        out = Dandelion(dat_)

    tmp_totaldist = pd.DataFrame(
        total_dist, index=dat_seq.index, columns=dat_seq.index)
    tmp_clusterdist = Tree()
    overlap = []
    for i in out.metadata.index:
        if len(out.metadata.loc[i, str(clonekey)].split('|')) > 1:
            overlap.append(
                [c for c in out.metadata.loc[i, str(clonekey)].split('|')])
            for c in out.metadata.loc[i, str(clonekey)].split('|'):
                tmp_clusterdist[c][i].value = 1
        else:
            cx = out.metadata.loc[i, str(clonekey)]
            tmp_clusterdist[cx][i].value = 1
    tmp_clusterdist2 = {}
    for x in tmp_clusterdist:
        tmp_clusterdist2[x] = list(tmp_clusterdist[x])
    cluster_dist = {}
    for c_ in tmp_clusterdist2:
        if c_ in list(flatten(overlap)):
            for ol in overlap:
                if c_ in ol:
                    idx = list(
                        set(flatten([tmp_clusterdist2[c_x] for c_x in ol])))
                    if len(list(set(idx))) > 1:
                        dist_mat_ = tmp_totaldist.loc[idx, idx]
                        s1, s2 = dist_mat_.shape
                        if s1 > 1 and s2 > 1:
                            cluster_dist['|'.join(ol)] = dist_mat_
        else:
            dist_mat_ = tmp_totaldist.loc[tmp_clusterdist2[c_],
                                          tmp_clusterdist2[c_]]
            s1, s2 = dist_mat_.shape
            if s1 > 1 and s2 > 1:
                cluster_dist[c_] = dist_mat_

    # to improve the visulisation and plotting efficiency, i will build a minimum spanning tree for each group/clone to connect the shortest path
    mst_tree = mst(cluster_dist)
    sleep(0.5)

    edge_list = Tree()
    if verbose:
        for c in tqdm(mst_tree, desc='Generating edge list '):
            G = nx.from_pandas_adjacency(mst_tree[c])
            edge_list[c] = nx.to_pandas_edgelist(G)
    else:
        for c in mst_tree:
            G = nx.from_pandas_adjacency(mst_tree[c])
            edge_list[c] = nx.to_pandas_edgelist(G)

    sleep(0.5)

    clone_ref = dict(out.metadata[clonekey])
    tmp_clone_tree = Tree()
    for x in out.metadata.index:
        if '|' in clone_ref[x]:
            for x_ in clone_ref[x].split('|'):
                tmp_clone_tree[x_][x].value = 1
        else:
            tmp_clone_tree[clone_ref[x]][x].value = 1
    tmp_clone_tree2 = Tree()
    for x in tmp_clone_tree:
        tmp_clone_tree2[x] = list(tmp_clone_tree[x])

    tmp_clone_tree3 = Tree()
    tmp_clone_tree3_overlap = Tree()
    for x in tmp_clone_tree2:
        # this is to catch all possible cells that may potentially match up with this clone that's joined together
        if x in list(flatten(overlap)):
            for ol in overlap:
                if x in ol:
                    if len(tmp_clone_tree2[x]) > 1:
                        for x_ in tmp_clone_tree2[x]:
                            tmp_clone_tree3_overlap['|'.join(
                                ol)][''.join(x_)].value = 1
                    else:
                        tmp_clone_tree3_overlap['|'.join(ol)][''.join(
                            tmp_clone_tree2[x])].value = 1
        else:
            tmp_ = pd.DataFrame(
                index=tmp_clone_tree2[x], columns=tmp_clone_tree2[x])
            tmp_ = pd.DataFrame(
                np.tril(tmp_) + 1, index=tmp_clone_tree2[x], columns=tmp_clone_tree2[x])
            tmp_.fillna(0, inplace=True)
            tmp_clone_tree3[x] = tmp_

    for x in tmp_clone_tree3_overlap:  # repeat for the overlap clones
        tmp_ = pd.DataFrame(
            index=tmp_clone_tree3_overlap[x], columns=tmp_clone_tree3_overlap[x])
        tmp_ = pd.DataFrame(np.tril(
            tmp_) + 1, index=tmp_clone_tree3_overlap[x], columns=tmp_clone_tree3_overlap[x])
        tmp_.fillna(0, inplace=True)
        tmp_clone_tree3[x] = tmp_

    # here I'm using a temporary edge list to catch all cells that were identified as clones to forcefully link them up if they were identical but clipped off during the mst step

    # create a dataframe to recall the actual distance quickly
    tmp_totaldiststack = pd.DataFrame(tmp_totaldist.unstack())
    tmp_totaldiststack.index.names = [None, None]
    tmp_totaldiststack = tmp_totaldiststack.reset_index(drop=False)
    tmp_totaldiststack.columns = ['source', 'target', 'weight']
    tmp_totaldiststack.index = [str(s) + '|' + str(t) for s, t in zip(
        tmp_totaldiststack['source'], tmp_totaldiststack['target'])]
    tmp_totaldiststack['keep'] = [False if len(
        list(set(i.split('|')))) == 1 else True for i in tmp_totaldiststack.index]
    tmp_totaldiststack = tmp_totaldiststack[tmp_totaldiststack.keep].drop(
        'keep', axis=1)

    tmp_edge_list = Tree()
    if verbose:
        for c in tqdm(tmp_clone_tree3, desc='Linking edges '):
            if len(tmp_clone_tree3[c]) > 1:
                G = nx.from_pandas_adjacency(tmp_clone_tree3[c])
                tmp_edge_list[c] = nx.to_pandas_edgelist(G)
                tmp_edge_list[c].index = [str(
                    s) + '|' + str(t) for s, t in zip(tmp_edge_list[c]['source'], tmp_edge_list[c]['target'])]
                tmp_edge_list[c]['weight'].update(tmp_totaldiststack['weight'])
                # keep only edges when there is 100% identity, to minimise crowding
                tmp_edge_list[c] = tmp_edge_list[c][tmp_edge_list[c]
                                                    ['weight'] == 0]
                tmp_edge_list[c].reset_index(inplace=True)
    else:
        for c in tmp_clone_tree3:
            if len(tmp_clone_tree3[c]) > 1:
                G = nx.from_pandas_adjacency(tmp_clone_tree3[c])
                tmp_edge_list[c] = nx.to_pandas_edgelist(G)
                tmp_edge_list[c].index = [str(
                    s) + '|' + str(t) for s, t in zip(tmp_edge_list[c]['source'], tmp_edge_list[c]['target'])]
                tmp_edge_list[c]['weight'].update(tmp_totaldiststack['weight'])
                # keep only edges when there is 100% identity, to minimise crowding
                tmp_edge_list[c] = tmp_edge_list[c][tmp_edge_list[c]
                                                    ['weight'] == 0]
                tmp_edge_list[c].reset_index(inplace=True)

    # try to catch situations where there's no edge (only singletons)
    try:
        edge_listx = pd.concat([edge_list[x] for x in edge_list])
        edge_listx.index = [
            str(s) + '|' + str(t) for s, t in zip(edge_listx['source'], edge_listx['target'])]

        tmp_edge_listx = pd.concat([tmp_edge_list[x] for x in tmp_edge_list])
        tmp_edge_listx.drop('index', axis=1, inplace=True)
        tmp_edge_listx.index = [str(
            s) + '|' + str(t) for s, t in zip(tmp_edge_listx['source'], tmp_edge_listx['target'])]

        edge_list_final = edge_listx.combine_first(tmp_edge_listx)
        edge_list_final['weight'].update(tmp_totaldiststack['weight'])
        # return the edge list
        edge_list_final.reset_index(drop=True, inplace=True)
    except:
        edge_list_final = None

    # and finally the vertex list which is super easy
    vertice_list = list(out.metadata.index)

    # and now to actually generate the network
    g, g_, lyt, lyt_ = generate_layout(
        vertice_list, edge_list_final, min_size=min_size, weight=None, verbose=verbose, **kwargs)

    # convert distance matrices to sparse
    for x in dmat:
        if type(dmat[x]) is np.ndarray:
            dmat[x] = csr_matrix(dmat[x])

    if verbose:
        logg.info(' finished', time=start,
                  deep=('Updated Dandelion object: \n'
                        '   \'data\', contig-indexed clone table\n'
                        '   \'metadata\', cell-indexed clone table\n'
                        '   \'distance\', heavy and light chain distance matrices\n'
                        '   \'edges\', network edges\n'
                        '   \'layout\', network layout\n'
                        '   \'graph\', network'))
    if self.__class__ == Dandelion:
        if self.germline is not None:
            germline_ = self.germline
        else:
            germline_ = None
        if self.threshold is not None:
            threshold_ = self.threshold
        else:
            threshold_ = None
        if downsample is not None:
            # out = Dandelion(data = dat_downsample, metadata = downsample_meta, distance = dmat, edges = edge_list_final, layout = (lyt, lyt_), graph = (g, g_), germline = germline_)
            out = Dandelion(data=dat_, distance=dmat, edges=edge_list_final, layout=(
                lyt, lyt_), graph=(g, g_), germline=germline_)
            out.threshold = threshold_
            return(out)
        else:
            self.__init__(data=self.data, metadata=self.metadata, distance=dmat, edges=edge_list_final, layout=(
                lyt, lyt_), graph=(g, g_), germline=germline_, initialize=False)
            self.threshold = threshold_
    else:
        # out = Dandelion(data = dat, distance = dmat, edges = edge_list_final, layout = (lyt, lyt_), graph = (g, g_), clone_key = clone_key)
        out = Dandelion(data=dat_, distance=dmat, edges=edge_list_final, layout=(
            lyt, lyt_), graph=(g, g_), clone_key=clone_key)
        return(out)


def mst(mat: dict) -> Tree:
    """
    Construct minimum spanning tree based on supplied matrix in dictionary.

    Parameters
    ----------
    mat : dict
        Dictionary containing numpy ndarrays.

    Returns
    -------
    Dandelion `Tree` object holding DataFrames of constructed minimum spanning trees.
    """
    mst_tree = Tree()
    for c in mat:
        mst_tree[c] = pd.DataFrame(minimum_spanning_tree(np.triu(mat[c])).toarray(
        ).astype(int), index=mat[c].index, columns=mat[c].columns)
    return(mst_tree)


def clone_degree(self: Dandelion, weight: Union[None, str] = None, verbose: bool = True) -> Dandelion:
    """
    Calculates node degree in BCR network.

    Parameters
    ----------
    self : Dandelion
        `Dandelion` object after `tl.generate_network` has been run.
    weight : str, optional
        Atribute name for retrieving edge weight in graph. None defaults to ignoring this. See `networkx.Graph.degree`.
    verbose : bool
        Whether or not to show logging information.

    Returns
    -------
    Dandelion object with metadata updated with node degree information.
    """
    if verbose:
        start = logg.info('Calculating node degree')
    if self.__class__ == Dandelion:
        try:
            G = self.graph[0]
        except:
            dist = np.sum([self.distance[x].toarray() for x in self.distance if type(
                self.distance[x]) is csr_matrix], axis=0)
            A = csr_matrix(dist)
            G = nx.Graph()
            G.add_weighted_edges_from(
                zip(list(self.metadata.index), list(self.metadata.index), A.data))

        if len(G) == 0:
            raise AttributeError(
                'Graph not found. Plase run tl.generate_network.')
        else:
            cd = pd.DataFrame.from_dict(G.degree(weight=weight))
            cd.set_index(0, inplace=True)
            self.metadata['clone_degree'] = pd.Series(cd[1])
            if verbose:
                logg.info(' finished', time=start,
                          deep=('Updated Dandelion metadata\n'))
    else:
        raise TypeError('Input object must be of {}'.format(Dandelion))


def clone_centrality(self: Dandelion, verbose: bool = True) -> Dandelion:
    """
    Calculates node closeness centrality in BCR network.

    Parameters
    ----------
    self : Dandelion
        `Dandelion` object after `tl.generate_network` has been run.
    verbose : bool
        Whether or not to show logging information.

    Returns
    -------
    Dandelion object with metadata updated with node closeness centrality information.
    """
    if verbose:
        start = logg.info('Calculating node closeness centrality')
    if self.__class__ == Dandelion:
        try:
            G = self.graph[0]
        except:
            dist = np.sum([self.distance[x].toarray() for x in self.distance if type(
                self.distance[x]) is csr_matrix], axis=0)
            A = csr_matrix(dist)
            G = nx.Graph()
            G.add_weighted_edges_from(
                zip(list(self.metadata.index), list(self.metadata.index), A.data))

        if len(G) == 0:
            raise AttributeError(
                'Graph not found. Plase run tl.generate_network.')
        else:
            cc = nx.closeness_centrality(G)
            cc = pd.DataFrame.from_dict(
                cc, orient='index', columns=['clone_centrality'])
            self.metadata['clone_centrality'] = pd.Series(
                cc['clone_centrality'])
            if verbose:
                logg.info(' finished', time=start,
                          deep=('Updated Dandelion metadata\n'))
    else:
        raise TypeError('Input object must be of {}'.format(Dandelion))


def generate_layout(vertices: Sequence, edges: pd.DataFrame = None, min_size: int = 2, weight: Union[None, str] = None, verbose: bool = True, **kwargs) -> Tuple[nx.Graph, nx.Graph, dict, dict]:
    G = nx.Graph()
    G.add_nodes_from(vertices)
    if edges is not None:
        G.add_weighted_edges_from([(x, y, z) for x, y, z in zip(
            edges['source'], edges['target'], edges['weight'])])
    degree = G.degree()
    G_ = G.copy()
    if min_size == 2:
        if edges is not None:
            G_.remove_nodes_from(nx.isolates(G))
        else:
            pass
    elif min_size > 2:
        if edges is not None:
            remove = [node for node, degree in dict(
                G.degree()).items() if degree > min_size]
            G_.remove_nodes_from(remove)
        else:
            pass
    # if edges is not None:
        # edges_, weights_ = zip(*nx.get_edge_attributes(G_,'weight').items())
    if verbose:
        print('generating network layout')
    pos = _fruchterman_reingold_layout(G, weight=weight, **kwargs)
    pos_ = _fruchterman_reingold_layout(G_, weight=weight, **kwargs)
    return(G, G_, pos, pos_)

# when dealing with a lot of unconnected vertices, the pieces fly out to infinity and the original fr layout can't be used
# work around from https://stackoverflow.com/questions/14283341/how-to-increase-node-spacing-for-networkx-spring-layout
# code chunk from networkx's layout.py https://github.com/networkx/networkx/blob/master/networkx/drawing/layout.py


def _process_params(G, center, dim):
    # Some boilerplate code.

    if not isinstance(G, nx.Graph):
        empty_graph = nx.Graph()
        empty_graph.add_nodes_from(G)
        G = empty_graph

    if center is None:
        center = np.zeros(dim)
    else:
        center = np.asarray(center)

    if len(center) != dim:
        msg = "length of center coordinates must match dimension of layout"
        raise ValueError(msg)

    return G, center


def _fruchterman_reingold_layout(
    G,
    k=None,
    pos=None,
    fixed=None,
    iterations=50,
    threshold=1e-4,
    weight="weight",
    scale=1,
    center=None,
    dim=2,
    seed=None,
    **kwargs,
):
    """Position nodes using Fruchterman-Reingold force-directed algorithm.
    The algorithm simulates a force-directed representation of the network
    treating edges as springs holding nodes close, while treating nodes
    as repelling objects, sometimes called an anti-gravity force.
    Simulation continues until the positions are close to an equilibrium.
    There are some hard-coded values: minimal distance between
    nodes (0.01) and "temperature" of 0.1 to ensure nodes don't fly away.
    During the simulation, `k` helps determine the distance between nodes,
    though `scale` and `center` determine the size and place after
    rescaling occurs at the end of the simulation.
    Fixing some nodes doesn't allow them to move in the simulation.
    It also turns off the rescaling feature at the simulation's end.
    In addition, setting `scale` to `None` turns off rescaling.

    Parameters
    ----------
    G : NetworkX graph or list of nodes
        A position will be assigned to every node in G.
    k : float (default=None)
        Optimal distance between nodes.  If None the distance is set to
        1/sqrt(n) where n is the number of nodes.  Increase this value
        to move nodes farther apart.
    pos : dict or None  optional (default=None)
        Initial positions for nodes as a dictionary with node as keys
        and values as a coordinate list or tuple.  If None, then use
        random initial positions.
    fixed : list or None  optional (default=None)
        Nodes to keep fixed at initial position.
        ValueError raised if `fixed` specified and `pos` not.
    iterations : int  optional (default=50)
        Maximum number of iterations taken
    threshold: float optional (default = 1e-4)
        Threshold for relative error in node position changes.
        The iteration stops if the error is below this threshold.
    weight : string or None   optional (default='weight')
        The edge attribute that holds the numerical value used for
        the edge weight.  If None, then all edge weights are 1.
    scale : number or None (default: 1)
        Scale factor for positions. Not used unless `fixed is None`.
        If scale is None, no rescaling is performed.
    center : array-like or None
        Coordinate pair around which to center the layout.
        Not used unless `fixed is None`.
    dim : int
        Dimension of layout.
    seed : int, RandomState instance or None  optional (default=None)
        Set the random state for deterministic node layouts.
        If int, `seed` is the seed used by the random number generator,
        if numpy.random.RandomState instance, `seed` is the random
        number generator,
        if None, the random number generator is the RandomState instance used
        by numpy.random.

    Returns
    -------
    pos : dict
        A dictionary of positions keyed by node

    Examples
    --------
    >>> G = nx.path_graph(4)
    >>> pos = nx.spring_layout(G)
    # The same using longer but equivalent function name
    >>> pos = nx.fruchterman_reingold_layout(G)
    """

    G, center = _process_params(G, center, dim)

    if fixed is not None:
        if pos is None:
            raise ValueError("nodes are fixed without positions given")
        for node in fixed:
            if node not in pos:
                raise ValueError("nodes are fixed without positions given")
        nfixed = {node: i for i, node in enumerate(G)}
        fixed = np.asarray([nfixed[node] for node in fixed])

    if pos is not None:
        # Determine size of existing domain to adjust initial positions
        dom_size = max(coord for pos_tup in pos.values() for coord in pos_tup)
        if dom_size == 0:
            dom_size = 1
        pos_arr = seed.rand(len(G), dim) * dom_size + center

        for i, n in enumerate(G):
            if n in pos:
                pos_arr[i] = np.asarray(pos[n])
    else:
        pos_arr = None
        dom_size = 1

    if len(G) == 0:
        return {}
    if len(G) == 1:
        return {nx.utils.arbitrary_element(G.nodes()): center}

    try:
        # Sparse matrix
        if len(G) < 500:  # sparse solver for large graphs
            raise ValueError
        A = nx.to_scipy_sparse_matrix(G, weight=weight, dtype="f")
        if k is None and fixed is not None:
            # We must adjust k by domain size for layouts not near 1x1
            nnodes, _ = A.shape
            k = dom_size / np.sqrt(nnodes)
        pos = _sparse_fruchterman_reingold(
            A, k, pos_arr, fixed, iterations, threshold, dim, seed
        )
    except ValueError:
        A = nx.to_numpy_array(G, weight=weight)
        if k is None and fixed is not None:
            # We must adjust k by domain size for layouts not near 1x1
            nnodes, _ = A.shape
            k = dom_size / np.sqrt(nnodes)
        pos = _fruchterman_reingold(
            A, k, pos_arr, fixed, iterations, threshold, dim, seed
        )
    if fixed is None and scale is not None:
        pos = _rescale_layout(pos, scale=scale) + center
    pos = dict(zip(G, pos))
    return pos


@random_state(7)
def _fruchterman_reingold(
    A, k=None, pos=None, fixed=None, iterations=50, threshold=1e-4, dim=2, seed=None, **kwargs,
):
    # Position nodes in adjacency matrix A using Fruchterman-Reingold
    # Entry point for NetworkX graph is fruchterman_reingold_layout()
    import numpy as np

    try:
        nnodes, _ = A.shape
    except AttributeError as e:
        msg = "fruchterman_reingold() takes an adjacency matrix as input"
        raise nx.NetworkXError(msg) from e

    if pos is None:
        # random initial positions
        pos = np.asarray(seed.rand(nnodes, dim), dtype=A.dtype)
    else:
        # make sure positions are of same type as matrix
        pos = pos.astype(A.dtype)

    # optimal distance between nodes
    if k is None:
        k = np.sqrt(1.0 / nnodes)
    # the initial "temperature"  is about .1 of domain area (=1x1)
    # this is the largest step allowed in the dynamics.
    # We need to calculate this in case our fixed positions force our domain
    # to be much bigger than 1x1
    t = max(max(pos.T[0]) - min(pos.T[0]), max(pos.T[1]) - min(pos.T[1])) * 0.1
    # simple cooling scheme.
    # linearly step down by dt on each iteration so last iteration is size dt.
    dt = t / float(iterations + 1)
    delta = np.zeros((pos.shape[0], pos.shape[0], pos.shape[1]), dtype=A.dtype)
    # the inscrutable (but fast) version
    # this is still O(V^2)
    # could use multilevel methods to speed this up significantly
    for iteration in range(iterations):
        # matrix of difference between points
        delta = pos[:, np.newaxis, :] - pos[np.newaxis, :, :]
        # distance between points
        distance = np.linalg.norm(delta, axis=-1)
        # enforce minimum distance of 0.01
        np.clip(distance, 0.001, None, out=distance)
        # displacement "force"
        displacement = np.einsum(
            "ijk,ij->ik", delta, (k * k / distance ** 2 - A * distance / k))
        displacement = displacement - pos / (k * np.sqrt(nnodes))
        # update positions
        length = np.linalg.norm(displacement, axis=-1)
        length = np.where(length < 0.01, 0.1, length)
        delta_pos = np.einsum("ij,i->ij", displacement, t / length)
        if fixed is not None:
            # don't change positions of fixed nodes
            delta_pos[fixed] = 0.0
        pos += delta_pos
        # cool temperature
        t -= dt
        err = np.linalg.norm(delta_pos) / nnodes
        if err < threshold:
            break
    return pos


@random_state(7)
def _sparse_fruchterman_reingold(
    A, k=None, pos=None, fixed=None, iterations=50, threshold=1e-4, dim=2, seed=None, **kwargs,
):
    # Position nodes in adjacency matrix A using Fruchterman-Reingold
    # Entry point for NetworkX graph is fruchterman_reingold_layout()
    # Sparse version
    import numpy as np

    try:
        nnodes, _ = A.shape
    except AttributeError as e:
        msg = "fruchterman_reingold() takes an adjacency matrix as input"
        raise nx.NetworkXError(msg) from e
    try:
        from scipy.sparse import coo_matrix
    except ImportError as e:
        msg = "_sparse_fruchterman_reingold() scipy numpy: http://scipy.org/ "
        raise ImportError(msg) from e
    # make sure we have a LIst of Lists representation
    try:
        A = A.tolil()
    except AttributeError:
        A = (coo_matrix(A)).tolil()

    if pos is None:
        # random initial positions
        pos = np.asarray(seed.rand(nnodes, dim), dtype=A.dtype)
    else:
        # make sure positions are of same type as matrix
        pos = pos.astype(A.dtype)

    # no fixed nodes
    if fixed is None:
        fixed = []

    # optimal distance between nodes
    if k is None:
        k = np.sqrt(1.0 / nnodes)
    # the initial "temperature"  is about .1 of domain area (=1x1)
    # this is the largest step allowed in the dynamics.
    t = max(max(pos.T[0]) - min(pos.T[0]), max(pos.T[1]) - min(pos.T[1])) * 0.1
    # simple cooling scheme.
    # linearly step down by dt on each iteration so last iteration is size dt.
    dt = t / float(iterations + 1)

    displacement = np.zeros((dim, nnodes))
    for iteration in range(iterations):
        displacement *= 0
        # loop over rows
        for i in range(A.shape[0]):
            if i in fixed:
                continue
            # difference between this row's node position and all others
            delta = (pos[i] - pos).T
            # distance between points
            distance = np.sqrt((delta ** 2).sum(axis=0))
            # enforce minimum distance of 0.01
            distance = np.where(distance < 0.01, 0.01, distance)
            # the adjacency matrix row
            Ai = np.asarray(A.getrowview(i).toarray())
            # displacement "force"
            displacement[:, i] += (
                delta * (k * k / distance ** 2 - Ai * distance / k)
            ).sum(axis=1)
        displacement = displacement - pos / (k * np.sqrt(nnodes))
        # update positions
        length = np.sqrt((displacement ** 2).sum(axis=0))
        length = np.where(length < 0.01, 0.1, length)
        delta_pos = (displacement * t / length).T
        pos += delta_pos
        # cool temperature
        t -= dt
        err = np.linalg.norm(delta_pos) / nnodes
        if err < threshold:
            break
    return pos


def _rescale_layout(pos, scale=1):
    """
    Returns scaled position array to (-scale, scale) in all axes.
    The function acts on NumPy arrays which hold position information.
    Each position is one row of the array. The dimension of the space
    equals the number of columns. Each coordinate in one column.
    To rescale, the mean (center) is subtracted from each axis separately.
    Then all values are scaled so that the largest magnitude value
    from all axes equals `scale` (thus, the aspect ratio is preserved).
    The resulting NumPy Array is returned (order of rows unchanged).

    Parameters
    ----------
    pos : numpy array
        positions to be scaled. Each row is a position.
    scale : number (default: 1)
        The size of the resulting extent in all directions.

    Returns
    -------
    pos : numpy array
        scaled positions. Each row is a position.
    """
    # Find max length over all dimensions
    lim = 0  # max coordinate for all axes
    for i in range(pos.shape[1]):
        pos[:, i] -= pos[:, i].mean()
        lim = max(abs(pos[:, i]).max(), lim)
    # rescale to (-scale, scale) in all directions, preserves aspect
    if lim > 0:
        for i in range(pos.shape[1]):
            pos[:, i] *= scale / lim
    return pos


def extract_edge_weights(self: Dandelion, expanded_only: bool = False) -> Sequence:
    """
    Retrieves edge weights (BCR levenshtein distance) from graph.

    Parameters
    ----------
    self : Dandelion
        `Dandelion` object after `tl.generate_network` has been run.
    expanded_only : bool
        whether to retrieve the edge weights from the expanded only graph or entire graph.

    Returns
    -------
    numpy array containing edge weights.
    """
    if expanded_only:
        try:
            edges, weights = zip(
                *nx.get_edge_attributes(self.graph[1], 'weight').items())
        except ValueError as e:
            print(
                '{} i.e. the graph does not contain edges. Therefore, edge weights not returned.'.format(e))
    else:
        try:
            edges, weights = zip(
                *nx.get_edge_attributes(self.graph[0], 'weight').items())
        except ValueError as e:
            print(
                '{} i.e. the graph does not contain edges. Therefore, edge weights not returned.'.format(e))
    if 'weights' in locals():
        return(weights)
