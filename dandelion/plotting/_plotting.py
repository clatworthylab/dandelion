#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2020-05-18 00:15:00
# @Last Modified by:   Kelvin
# @Last Modified time: 2021-05-20 14:21:45

import seaborn as sns
import pandas as pd
import numpy as np
from ..utilities._utilities import *
from ..utilities._core import *
from ..utilities._io import *
from ..tools._diversity import rarefun
from scanpy.plotting._tools.scatterplots import embedding
import matplotlib.pyplot as plt
from anndata import AnnData
import random
from adjustText import adjust_text
from plotnine import ggplot, theme_classic, aes, geom_line, xlab, ylab, options, ggtitle, labs, scale_color_manual
from scanpy.plotting import palettes
from time import sleep
import matplotlib.pyplot as plt
from itertools import combinations
from typing import Union, Sequence, Tuple, Dict
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from tqdm import tqdm


def clone_rarefaction(self: Union[AnnData, Dandelion], color: str, clone_key: Union[None, str] = None, palette: Union[None, Sequence] = None, figsize: Tuple[Union[int, float], Union[int, float]] = (6, 4), save: Union[None, str] = None) -> ggplot:
    """
    Plots rarefaction curve for cell numbers vs clone size.

    Parameters
    ----------
    self : `AnnData`, `Dandelion`
        `AnnData` or `Dandelion` object.
    color : str
        Column name to split the calculation of clone numbers for a given number of cells for e.g. sample, patient etc.
    clone_key : str, optional
        Column name specifying the clone_id column in metadata/obs.
    palette : Sequence, optional
        Color mapping for unique elements in color. Will try to retrieve from AnnData `.uns` slot if present.
    figsize :  Tuple[Union[int,float], Union[int,float]]
        Size of plot.
    save : str, optional
        Save path.

    Returns
    -------
    rarefaction curve plot.
    """

    if self.__class__ == AnnData:
        metadata = self.obs.copy()
    elif self.__class__ == Dandelion:
        metadata = self.metadata.copy()
    if clone_key is None:
        clonekey = 'clone_id'
    else:
        clonekey = clone_key

    groups = list(set(metadata[color]))
    metadata = metadata[metadata['contig_QC_pass'].isin([True, 'True'])]
    metadata[clonekey] = metadata[clonekey].cat.remove_unused_categories()
    res = {}
    for g in groups:
        _metadata = metadata[metadata[color] == g]
        res[g] = _metadata[clonekey].value_counts()
    res_ = pd.DataFrame.from_dict(res, orient='index')

    # remove those with no counts
    rowsum = res_.sum(axis=1)
    print('removing due to zero counts:', ', '.join(
        [res_.index[i] for i, x in enumerate(res_.sum(axis=1) == 0) if x]))
    sleep(0.5)
    res_ = res_[~(res_.sum(axis=1) == 0)]

    # set up for calculating rarefaction
    tot = res_.apply(sum, axis=1)
    S = res_.apply(lambda x: x[x > 0].shape[0], axis=1)
    nr = res_.shape[0]

    # append the results to a dictionary
    rarecurve = {}
    for i in tqdm(range(0, nr), desc='Calculating rarefaction curve '):
        n = np.arange(1, tot[i], step=10)
        if (n[-1:] != tot[i]):
            n = np.append(n, tot[i])
        rarecurve[res_.index[i]] = [
            rarefun(np.array(res_.iloc[i, ]), z) for z in n]
    y = pd.DataFrame([rarecurve[c] for c in rarecurve]).T
    pred = pd.DataFrame([np.append(np.arange(1, s, 10), s)
                         for s in res_.sum(axis=1)], index=res_.index).T

    y = y.melt()
    pred = pred.melt()
    pred['yhat'] = y['value']

    options.figure_size = figsize
    if palette is None:
        if self.__class__ == AnnData:
            try:
                pal = self.uns[str(color) + '_colors']
            except:
                if len(list(set((pred.variable)))) <= 20:
                    pal = palettes.default_20
                elif len(list(set((pred.variable)))) <= 28:
                    pal = palettes.default_28
                elif len(list(set((pred.variable)))) <= 102:
                    pal = palettes.default_102
                else:
                    pal = None

            if pal is not None:
                p = (ggplot(pred, aes(x="value", y="yhat", color="variable"))
                     + theme_classic()
                     + xlab('number of cells')
                     + ylab('number of clones')
                     + ggtitle('rarefaction curve')
                     + labs(color=color)
                     + scale_color_manual(values=(pal))
                     + geom_line())
            else:
                p = (ggplot(pred, aes(x="value", y="yhat", color="variable"))
                     + theme_classic()
                     + xlab('number of cells')
                     + ylab('number of clones')
                     + ggtitle('rarefaction curve')
                     + labs(color=color)
                     + geom_line())
        else:
            if len(list(set((pred.variable)))) <= 20:
                pal = palettes.default_20
            elif len(list(set((pred.variable)))) <= 28:
                pal = palettes.default_28
            elif len(list(set((pred.variable)))) <= 102:
                pal = palettes.default_102
            else:
                pal = None

            if pal is not None:
                p = (ggplot(pred, aes(x="value", y="yhat", color="variable"))
                     + theme_classic()
                     + xlab('number of cells')
                     + ylab('number of clones')
                     + ggtitle('rarefaction curve')
                     + labs(color=color)
                     + scale_color_manual(values=(pal))
                     + geom_line())
            else:
                p = (ggplot(pred, aes(x="value", y="yhat", color="variable"))
                     + theme_classic()
                     + xlab('number of cells')
                     + ylab('number of clones')
                     + ggtitle('rarefaction curve')
                     + labs(color=color)
                     + geom_line())
    else:
        p = (ggplot(pred, aes(x="value", y="yhat", color="variable"))
             + theme_classic()
             + xlab('number of cells')
             + ylab('number of clones')
             + ggtitle('rarefaction curve')
             + labs(color=color)
             + geom_line())
    if save:
        p.save(filename='figures/rarefaction' + str(save),
               height=plt.rcParams['figure.figsize'][0], width=plt.rcParams['figure.figsize'][1], units='in', dpi=plt.rcParams["savefig.dpi"])

    return(p)


def random_palette(n: int) -> Sequence:
    # a list of 900+colours
    cols = list(sns.xkcd_rgb.keys())
    # if max_colors_needed1 > len(cols):
    cols2 = list(sns.color_palette('husl', n))
    palette = random.sample(sns.xkcd_palette(cols) + cols2, n)
    return(palette)


def clone_network(adata: AnnData, basis: str = 'vdj', edges: bool = True, **kwargs) -> Union[Figure, Axes, None]:
    """
    Using scanpy's plotting module to plot the network. Only thing that is changed is the dfault options: `basis = 'bcr'` and `edges = True`.

    Parameters
    ----------
    adata : AnnData
        AnnData object.
    basis : str
        key for embedding. Default is 'vdj'.
    edges : bool
        whether or not to plot edges. Default is True.
    **kwargs
        passed `sc.pl.embedding`.
    """
    embedding(adata, basis=basis, edges=edges, **kwargs)


def barplot(self: Union[AnnData, Dandelion], color: str, palette: str = 'Set1', figsize: Tuple[Union[int, float], Union[int, float]] = (12, 4), normalize: bool = True, sort_descending: bool = True, title: Union[None, str] = None, xtick_rotation: Union[None, Union[int, float]] = None, min_clone_size: Union[None, int] = None, clone_key: Union[None, str] = None, **kwargs) -> Tuple[Figure, Axes]:
    """
    A barplot function to plot usage of V/J genes in the data.

    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    color : str
        column name in metadata for plotting in bar plot.
    palette : str
        Colors to use for the different levels of the color variable. Should be something that can be interpreted by [color_palette](https://seaborn.pydata.org/generated/seaborn.color_palette.html#seaborn.color_palette), or a dictionary mapping hue levels to matplotlib colors. See [seaborn.barplot](https://seaborn.pydata.org/generated/seaborn.barplot.html).
    figsize : Tuple[Union[int,float], Union[int,float]]
        figure size. Default is (12, 4).
    normalize : bool
        if True, will return as proportion out of 1, otherwise False will return counts. Default is True.
    sort_descending : bool
        whether or not to sort the order of the plot. Default is True.
    title : str, optional
        title of plot.
    xtick_rotation : int, float, optional
        rotation of x tick labels.
    min_clone_size : int, optional
        minimum clone size to keep. Defaults to 1 if left as None.
    clone_key : str, optional
        column name for clones. None defaults to 'clone_id'.
    **kwargs
        passed to `sns.barplot`.

    Returns
    -------
    a seaborn barplot.
    """
    if self.__class__ == Dandelion:
        data = self.metadata.copy()
    elif self.__class__ == AnnData:
        data = self.obs.copy()

    if min_clone_size is None:
        min_size = 1
    else:
        min_size = int(min_clone_size)

    if clone_key is None:
        clone_ = 'clone_id'
    else:
        clone_ = clone_key

    size = data[clone_].value_counts()
    keep = list(size[size >= min_size].index)
    data_ = data[data[clone_].isin(keep)]

    sns.set_style('whitegrid', {'axes.grid': False})
    res = pd.DataFrame(data_[color].value_counts(normalize=normalize))
    if not sort_descending:
        res = res.sort_index()
    res.reset_index(drop=False, inplace=True)

    # Initialize the matplotlib figure
    fig, ax = plt.subplots(figsize=figsize)

    # plot
    sns.barplot(x='index', y=color, data=res, palette=palette, **kwargs)
    # change some parts
    if title is None:
        ax.set_title(color.replace('_', ' ') + ' usage')
    else:
        ax.set_title(title)
    if normalize:
        ax.set_ylabel('proportion')
    else:
        ax.set_ylabel('count')
    ax.set_xlabel('')
    if xtick_rotation is None:
        plt.xticks(rotation=90)
    else:
        plt.xticks(rotation=xtick_rotation)
    return fig, ax


def stackedbarplot(self: Union[AnnData, Dandelion], color: str, groupby: Union[None, str], figsize: Tuple[Union[int, float], Union[int, float]] = (12, 4), normalize: bool = False, title: Union[None, str] = None, sort_descending: bool = True, xtick_rotation: Union[None, Union[float, int]] = None, hide_legend: bool = True, legend_options: Tuple[str, Tuple[float, float], int] = None, labels: Union[None, Sequence] = None, min_clone_size: Union[None, int] = None, clone_key: Union[None, str] = None, **kwargs) -> Tuple[Figure, Axes]:
    """
    A stackedbarplot function to plot usage of V/J genes in the data split by groups.

    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    color : str
        column name in metadata for plotting in bar plot.
    groupby : str
        column name in metadata to split by during plotting.
    figsize : Tuple[Union[int,float], Union[int,float]]
        figure size. Default is (12, 4).
    normalize : bool
        if True, will return as proportion out of 1, otherwise False will return counts. Default is True.
    sort_descending : bool
        whether or not to sort the order of the plot. Default is True.
    title : str, optional
        title of plot.
    xtick_rotation: Union[None, Union[float,int]] : int, float, optional
        rotation of x tick labels.
    hide_legend : bool
        whether or not to hide the legend.
    legend_options : Tuple[str, Tuple[float, float], int]
        a tuple holding 3 options for specify legend options: 1) loc (string), 2) bbox_to_anchor (tuple), 3) ncol (int).
    labels : Sequence, optional
        Names of objects will be used for the legend if list of multiple dataframes supplied.
    min_clone_size : int, optional
        minimum clone size to keep. Defaults to 1 if left as None.
    clone_key : str, optional
        column name for clones. None defaults to 'clone_id'.
    **kwargs
        other kwargs passed to `matplotlib.plt`.

    Returns
    -------
    stacked bar plot.
    """
    if self.__class__ == Dandelion:
        data = self.metadata.copy()
    elif self.__class__ == AnnData:
        data = self.obs.copy()
    # quick fix to prevent dropping of nan
    data[groupby] = [str(l) for l in data[groupby]]

    if min_clone_size is None:
        min_size = 1
    else:
        min_size = int(min_clone_size)

    if clone_key is None:
        clone_ = 'clone_id'
    else:
        clone_ = clone_key

    size = data[clone_].value_counts()
    keep = list(size[size >= min_size].index)
    data_ = data[data[clone_].isin(keep)]

    dat_ = pd.DataFrame(data_.groupby(color)[groupby].value_counts(
        normalize=normalize).unstack(fill_value=0).stack(), columns=['value'])
    dat_.reset_index(drop=False, inplace=True)
    dat_order = pd.DataFrame(data[color].value_counts(normalize=normalize))
    dat_ = dat_.pivot(index=color, columns=groupby, values='value')
    if sort_descending is True:
        dat_ = dat_.reindex(dat_order.index)
    elif sort_descending is False:
        dat_ = dat_.reindex(dat_order.index[::-1])
    elif sort_descending is None:
        dat_ = dat_.sort_index()

    def _plot_bar_stacked(dfall: pd.DataFrame, labels: Union[None, Sequence] = None, figsize: Tuple[Union[int, float], Union[int, float]] = (12, 4), title: str = "multiple stacked bar plot", xtick_rotation: Union[None, Union[float, int]] = None, legend_options: Tuple[str, Tuple[float, float], int] = None, hide_legend: bool = True, H: Literal["/"] = "/", **kwargs) -> Tuple[Figure, Axes]:
        """
        Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot.

        Parameters
        ----------
        labels
            a list of the dataframe objects. Names of objects will be used for the legend.
        title
            string for the title of the plot
        H
            is the hatch used for identification of the different dataframes
        **kwargs
            other kwargs passed to matplotlib.plt
        """
        if type(dfall) is not list:
            dfall = [dfall]
        n_df = len(dfall)
        n_col = len(dfall[0].columns)
        n_ind = len(dfall[0].index)
        # Initialize the matplotlib figure
        fig, ax = plt.subplots(figsize=figsize)
        for df in dfall:  # for each data frame
            ax = df.plot(kind="bar",
                         linewidth=0,
                         stacked=True,
                         ax=ax,
                         legend=False,
                         grid=False,
                         **kwargs)  # make bar plots
        h, l = ax.get_legend_handles_labels()  # get the handles we want to modify
        for i in range(0, n_df * n_col, n_col):  # len(h) = n_col * n_df
            for j, pa in enumerate(h[i:i + n_col]):
                for rect in pa.patches:  # for each index
                    rect.set_x(rect.get_x() + 1 / float(n_df + 1)
                               * i / float(n_col))
                    rect.set_hatch(H * int(i / n_col))  # edited part
                    rect.set_width(1 / float(n_df + 1))
        ax.set_xticks((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)
        ax.set_xticklabels(df.index, rotation=0)
        ax.set_title(title)
        if normalize:
            ax.set_ylabel('proportion')
        else:
            ax.set_ylabel('count')
        # Add invisible data to add another legend
        n = []
        for i in range(n_df):
            n.append(ax.bar(0, 0, color="grey", hatch=H * i))
        if legend_options is None:
            Legend = ('center right', (1.15, 0.5), 1)
        else:
            Legend = legend_options
        if hide_legend is False:
            l1 = ax.legend(h[:n_col], l[:n_col], loc=Legend[0],
                           bbox_to_anchor=Legend[1], ncol=Legend[2], frameon=False)
            if labels is not None:
                l2 = plt.legend(
                    n, labels, loc=Legend[0], bbox_to_anchor=Legend[1], ncol=Legend[2], frameon=False)
            ax.add_artist(l1)
        if xtick_rotation is None:
            plt.xticks(rotation=90)
        else:
            plt.xticks(rotation=xtick_rotation)

        return fig, ax

    if title is None:
        title = "multiple stacked bar plot : " + \
            color.replace('_', ' ') + ' usage'
    else:
        title = title

    return _plot_bar_stacked(dat_, labels=labels, figsize=figsize, title=title, xtick_rotation=xtick_rotation, legend_options=legend_options, hide_legend=hide_legend, **kwargs)


def spectratype(self: Union[AnnData, Dandelion], color: str, groupby: str, locus: str, clone_key: Union[None, str] = None, figsize: Tuple[Union[int, float], Union[int, float]] = (6, 4), width: Union[None, Union[int, float]] = None, title: Union[None, str] = None, xtick_rotation: Union[None, Union[float, int]] = None, hide_legend: bool = True, legend_options: Tuple[str, Tuple[float, float], int] = None, labels: Union[None, Sequence] = None, **kwargs) -> Tuple[Figure, Axes]:
    """
    A spectratype function to plot usage of CDR3 length in the data split by groups.

    Parameters
    ----------
    self : Dandelion, AnnData
        `Dandelion` or `AnnData` object.
    color : str
        column name in metadata for plotting in bar plot.
    groupby : str
        column name in metadata to split by during plotting.
    locus : str
        either IGH or IGL.
    figsize : Tuple[Union[int,float], Union[int,float]]
        figure size. Default is (6, 4).
    width : float, int, optional
        width of bars.
    title : str, optional
        title of plot.
    xtick_rotation : int, float, optional
        rotation of x tick labels.
    hide_legend : bool
        whether or not to hide the legend.
    legend_options : Tuple[str, Tuple[float, float], int]
        a tuple holding 3 options for specify legend options: 1) loc (string), 2) bbox_to_anchor (tuple), 3) ncol (int).
    labels : Sequence, optional
        Names of objects will be used for the legend if list of multiple dataframes supplied.
    **kwargs
        other kwargs passed to matplotlib.pyplot.plot

    Returns
    -------
    sectratype plot
    """

    if clone_key is None:
        clonekey = 'clone_id'
    else:
        clonekey = clone_key

    if self.__class__ == Dandelion:
        data = self.data.copy()
    else:
        try:
            data = self.copy()
        except:
            AttributeError(
                "Please provide a <class 'Dandelion'> class object or a pandas dataframe instead of %s." % self.__class__)

    if 'locus' not in data.columns:
        raise AttributeError("Please ensure dataframe contains 'locus' column")

    if type(locus) is not list:
        locus = [locus]
    data = data[data['locus'].isin(locus)]
    data[groupby] = [str(l) for l in data[groupby]]
    dat_ = pd.DataFrame(data.groupby(color)[groupby].value_counts(
        normalize=False).unstack(fill_value=0).stack(), columns=['value'])
    dat_.reset_index(drop=False, inplace=True)
    dat_[color] = pd.to_numeric(dat_[color], errors='coerce')
    dat_.sort_values(by=color)
    dat_2 = dat_.pivot(index=color, columns=groupby, values='value')
    new_index = range(0, int(dat_[color].max()) + 1)
    dat_2 = dat_2.reindex(new_index, fill_value=0)

    def _plot_spectra_stacked(dfall: pd.DataFrame, labels: Union[None, Sequence] = None, figsize: Tuple[Union[int, float], Union[int, float]] = (6, 4), title: str = "multiple stacked bar plot", width: Union[None, Union[int, float]] = None, xtick_rotation: Union[None, Union[float, int]] = None, legend_options: Tuple[str, Tuple[float, float], int] = None, hide_legend: bool = True, H: Literal["/"] = "/", **kwargs) -> Tuple[Figure, Axes]:
        if type(dfall) is not list:
            dfall = [dfall]
        n_df = len(dfall)
        n_col = len(dfall[0].columns)
        n_ind = len(dfall[0].index)
        if width is None:
            wdth = 0.1 * n_ind / 60 + 0.8
        else:
            wdth = width
        # Initialize the matplotlib figure
        fig, ax = plt.subplots(figsize=figsize)
        for df in dfall:  # for each data frame
            ax = df.plot(kind="bar",
                         linewidth=0,
                         stacked=True,
                         ax=ax,
                         legend=False,
                         grid=False,
                         **kwargs)  # make bar plots
        h, l = ax.get_legend_handles_labels()  # get the handles we want to modify
        for i in range(0, n_df * n_col, n_col):  # len(h) = n_col * n_df
            for j, pa in enumerate(h[i:i + n_col]):
                for rect in pa.patches:  # for each index
                    rect.set_x(rect.get_x() + 1 / float(n_df + 1)
                               * i / float(n_col))
                    rect.set_hatch(H * int(i / n_col))  # edited part
                    # need to see if there's a better way to toggle this.
                    rect.set_width(wdth)

        n = 5  # Keeps every 5th label visible and hides the rest
        [l.set_visible(False) for (i, l) in enumerate(
            ax.xaxis.get_ticklabels()) if i % n != 0]
        ax.set_title(title)
        ax.set_ylabel('count')
        # Add invisible data to add another legend
        n = []
        for i in range(n_df):
            n.append(ax.bar(0, 0, color="gray", hatch=H * i))
        if legend_options is None:
            Legend = ('center right', (1.25, 0.5), 1)
        else:
            Legend = legend_options
        if hide_legend is False:
            l1 = ax.legend(h[:n_col], l[:n_col], loc=Legend[0],
                           bbox_to_anchor=Legend[1], ncol=Legend[2], frameon=False)
            if labels is not None:
                l2 = plt.legend(
                    n, labels, loc=Legend[0], bbox_to_anchor=Legend[1], ncol=Legend[2], frameon=False)
            ax.add_artist(l1)
        if xtick_rotation is None:
            plt.xticks(rotation=0)
        else:
            plt.xticks(rotation=xtick_rotation)

        return fig, ax

    return _plot_spectra_stacked(dat_2, labels=labels, figsize=figsize, title=title, width=width, xtick_rotation=xtick_rotation, legend_options=legend_options, hide_legend=hide_legend, **kwargs)


def clone_overlap(self: Union[AnnData, Dandelion], groupby: str, colorby: str, min_clone_size: Union[None, int] = None, clone_key: Union[None, str] = None, color_mapping: Union[None, Sequence, Dict] = None, node_labels: bool = True, node_label_layout: Literal[None, 'rotation', 'numbers'] = 'rotation', group_label_position: Literal['beginning', 'middle', 'end'] = 'middle', group_label_offset: int = 8, figsize: Tuple[Union[int, float], Union[int, float]] = (8, 8), return_graph: bool = False, save: Union[None, str] = None, **kwargs):
    """
    A plot function to visualise clonal overlap as a circos-style plot. Requires nxviz.

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
    color_maopping : Dict, Sequence, optional
        custom color mapping provided as a sequence (correpsonding to order of categories or alpha-numeric order if dtype is not category), or dictionary containing custom {category:color} mapping.
    node_labels : bool, optional
        whether to use node objects as labels or not
    node_label_layout : bool, optional
        which/whether (a) node layout is used. One of 'rotation', 'numbers' or None.
    group_label_position : str
        The position of the group label. One of 'beginning', 'middle' or 'end'.
    group_label_offset : int, float
        how much to offset the group labels, so that they are not overlapping with node labels.
    figsize : Tuple[Union[int,float], Union[int,float]]
        figure size. Default is (8, 8).
    return_graph : bool
        whether or not to return the graph for fine tuning. Default is False.
    **kwargs
        passed to `matplotlib.pyplot.savefig`.

    Returns
    -------
    a `nxviz.CircosPlot`.
    """
    import networkx as nx
    try:
        import nxviz as nxv
    except:
        raise(ImportError("Unable to import module `nxviz`. Have you done install nxviz? Try pip install git+https://github.com/zktuong/nxviz.git"))

    if min_clone_size is None:
        min_size = 2
    else:
        min_size = int(min_clone_size)

    if clone_key is None:
        clone_ = 'clone_id'
    else:
        clone_ = clone_key

    if self.__class__ == AnnData:
        data = self.obs.copy()
        # get rid of problematic rows that appear because of category conversion?
        data = data[~(data[clone_].isin(
            [np.nan, 'nan', 'NaN', 'No_contig', 'unassigned', None]))]
        if 'clone_overlap' in self.uns:
            overlap = self.uns['clone_overlap'].copy()
        else:
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

            overlap = pd.crosstab(data[clone_], data[groupby])

            if min_size == 0:
                raise ValueError('min_size must be greater than 0.')
            elif min_size > 2:
                overlap[overlap < min_size] = 0
                overlap[overlap >= min_size] = 1
            elif min_size == 2:
                overlap[overlap >= min_size] = 1

            overlap.index.name = None
            overlap.columns.name = None
    elif self.__class__ == Dandelion:
        data = self.metadata.copy()
        # get rid of problematic rows that appear because of category conversion?
        data = data[~(data[clone_].isin(
            [np.nan, 'nan', 'NaN', 'No_contig', 'unassigned', None]))]

        # prepare a summary table
        datc_ = data[clone_].str.split('|', expand=True).stack()
        datc_ = pd.DataFrame(datc_)
        datc_.reset_index(drop=False, inplace=True)
        datc_.columns = ['cell_id', 'tmp', clone_]
        datc_.drop('tmp', inplace=True, axis=1)
        dictg_ = dict(data[groupby])
        datc_[groupby] = [dictg_[l] for l in datc_['cell_id']]

        overlap = pd.crosstab(data[clone_], data[groupby])

        if min_size == 0:
            raise ValueError('min_size must be greater than 0.')
        elif min_size > 2:
            overlap[overlap < min_size] = 0
            overlap[overlap >= min_size] = 1
        elif min_size == 2:
            overlap[overlap >= min_size] = 1

        overlap.index.name = None
        overlap.columns.name = None

    edges = {}
    for x in overlap.index:
        if overlap.loc[x].sum() > 1:
            edges[x] = [y + ({str(clone_): x},) for y in list(combinations(
                [i for i in overlap.loc[x][overlap.loc[x] == 1].index], 2))]

    # create graph
    G = nx.Graph()
    # add in the nodes
    G.add_nodes_from([(p, {str(colorby): d})
                      for p, d in zip(data[groupby], data[colorby])])

    # unpack the edgelist and add to the graph
    for edge in edges:
        G.add_edges_from(edges[edge])
    groupby_dict = dict(zip(data[groupby], data[colorby]))
    if color_mapping is None:
        if self.__class__ == AnnData:
            if pd.api.types.is_categorical_dtype(self.obs[groupby]):
                try:
                    colorby_dict = dict(zip(
                        list(self.obs[str(colorby)].cat.categories), self.uns[str(colorby) + '_colors']))
                except:
                    pass
    else:
        if type(color_mapping) is dict:
            colorby_dict = color_mapping
        else:
            if pd.api.types.is_categorical_dtype(data[groupby]):
                colorby_dict = dict(
                    zip(list(data[str(colorby)].cat.categories), color_mapping))
            else:
                colorby_dict = dict(
                    zip(sorted(list(set(data[str(colorby)]))), color_mapping))
    df = data[[groupby, colorby]]
    if groupby == colorby:
        df = data[[groupby]]
        df = df.sort_values(groupby).drop_duplicates(
            subset=groupby, keep="first").reset_index(drop=True)
    else:
        df = df.sort_values(colorby).drop_duplicates(
            subset=groupby, keep="first").reset_index(drop=True)

    c = nxv.CircosPlot(G,
                       node_color=colorby,
                       node_grouping=colorby,
                       node_labels=node_labels,
                       node_label_layout=node_label_layout,
                       group_label_position=group_label_position,
                       group_label_offset=group_label_offset,
                       figsize=figsize)
    c.nodes = list(df[groupby])
    if 'colorby_dict' in locals():
        c.node_colors = [colorby_dict[groupby_dict[c]] for c in c.nodes]
    c.compute_group_label_positions()
    c.compute_group_colors()
    c.draw()
    if save is not None:
        plt.savefig(save, bbox_inches='tight', **kwargs)
    if return_graph:
        return(c)
