import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import logging
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy as np
from statsmodels.tsa.stattools import acf
from matplotlib.gridspec import GridSpec
import scipy as sp

logger = logging.getLogger(__name__)

def boxplot_with_points(x, y, hue=None, data=None, legend=True, ax=None,
                        dist_plotter=sns.boxplot, dist_kwargs={},
                        points_plotter=sns.stripplot, points_kwargs={},
                        **kwargs):
    if not ax:
        fig, ax = plt.subplots()

    if data is None:
        data = {}
        data['x'] = x
        data['y'] = y
        x = 'x'
        y = 'y'
        if hue:
            data['hue'] = hue
            hue = 'hue'
        data = pd.DataFrame(data)

    data_kwargs = dict(data=data, x=x, y=y, hue=hue, **kwargs)

    distkw = data_kwargs.copy()
    distkw.update(dict(showfliers=False, ax=ax, saturation=0.35))
    distkw.update(dist_kwargs)
    dist_plotter(**distkw)
    handles, labels = ax.get_legend_handles_labels()

    pointskw = data_kwargs.copy()
    pointskw.update(dict(split=True, linewidth=1, ax=ax))
    pointskw.update(points_kwargs)
    points_plotter(**pointskw)
    if legend:
        ax.legend(handles=handles)

    return ax

def residuals_plot(fit, ax=None):
    if not ax:
        fig, ax = plt.subplots()

    x = fit.predict()
    y = fit.resid_pearson

    # Add trend line, but x must be strictly increasing for UnivariateSpline
    # so here's a serious hack. :-/
    delta = lambda x: x[-1:] - x[:-1]
    min_delta = delta(np.unique(x)).min()
    x = x + np.random.uniform(-min_delta/1000, min_delta/1000, len(x))
    x, y = zip(*sorted(zip(x, y)))
    spline = sp.interpolate.UnivariateSpline(x, y)

    ax.scatter(x, y)
    ax.plot(x, spline(x))
    ax.set_ylabel('Standardized Residuals')
    ax.set_xlabel('Fitted Value')

def pcaplot(data, x=0, y=1, z=None, hue=None,
            features=None, norm=True,
            vectors=True,
            pallete=None, legend=True, ax=None,
            vect_scale=4, vect_label_scale=8,
            **kwargs):

    data = data.copy()
    if not ax:
        fig = plt.figure()
        if z:
            ax = fig.add_subplot(1, 1, 1, projection='3d')
        else:
            ax = fig.add_subplot(1, 1, 1)
    if not features:
        features = data.columns
    if norm:
        data[features] = (data[features] - data[features].mean()) / data[features].std()


    pca = PCA().fit(data[features])
    labels = ['PC%d' % (i + 1) for i in range(pca.components_.shape[0])]
    projection = pd.DataFrame(pca.transform(data[features]),
                        index=data.index,
                        columns=labels)
    projection[data.columns] = data  # For grouping and such downstream
    components = pd.DataFrame(pca.components_,
                              index=labels,
                              columns=features)

    if hue is not None:
        groups = projection.groupby(hue)
    else:
        groups = [(None, data)]
    if not pallete:
        pallete = {name: color
                   for (name, _), color
                   in zip(groups, sns.color_palette(n_colors=len(groups)))}

    out = [ax]
    for name, _data in groups:
        if z:
            art = ax.scatter(_data[labels[x]], _data[labels[y]], _data[labels[z]],
                            label=name, c=pallete[name], **kwargs)
        else:
            art = ax.scatter(_data[labels[x]], _data[labels[y]],
                            label=name, color=pallete[name], **kwargs)
        out.append(art)

    if vectors is True:
        vectors = features
    elif vectors is False:
        vectors = []
    else:
        pass  # Allow :vectors: to be a list of names

    for vect_name in vectors:
        _x = components[vect_name][labels[x]]
        _y = components[vect_name][labels[y]]
        if z:
            _z = components[vect_name][labels[z]]
            ax.quiver(0, 0, 0,
                      _x * vect_scale, _y * vect_scale, _z * vect_scale,
                      color='black', pivot='tail')
            ax.text(_x * vect_label_scale,
                    _y * vect_label_scale,
                    _z * vect_label_scale,
                    s=10, zdir=(_x, _y, _z), text=vect_name)
            print(_x, _y, _z)
        else:
            ax.arrow(0, 0, _x * vect_scale, _y * vect_scale, head_width=0.08, color='black')
            ax.annotate(vect_name,
                        xy=(_x * vect_label_scale, _y * vect_label_scale),
                        fontsize=10,
                        horizontalalignment='center', verticalalignment='center')

    ax.set_xlabel('{} ({}%)'.format(labels[x], int(pca.explained_variance_ratio_[x] * 100)))
    ax.set_ylabel('{} ({}%)'.format(labels[y], int(pca.explained_variance_ratio_[y] * 100)))
    if z:
        ax.set_zlabel('{} ({}%)'.format(labels[z], int(pca.explained_variance_ratio_[z] * 100)))
    if legend:
        ax.legend()
    ax.set_aspect('equal', 'datalim')

    return(out)

def traceplot(x, axs=None, alpha=0.05, **kwargs):
    if axs is None:
        gs = GridSpec(2, 2, width_ratios=(4, 1), height_ratios=(2, 1), wspace=0.05, hspace=0.1)
        ax0 = plt.subplot(gs[0,0])
        ax1 = plt.subplot(gs[0,1], sharey=ax0)
        plt.setp(ax1.get_yticklabels(), visible=False)
        kdeplot_kwargs = {'vertical': True}
        ax2 = plt.subplot(gs[1,0], sharex=ax0)
        plt.setp(ax0.get_xticklabels(), visible=False)
        axs = np.array([ax0, ax1, ax2])
    else:
        axs = np.asarray(axs).flatten()
        assert len(axs) == 3

    axs[0].plot(x, **kwargs)
    sns.kdeplot(x, ax=axs[1], **kdeplot_kwargs, **kwargs)
    axs[1].set_xticklabels([])

    lags = np.arange(len(x))
    nlags = len(lags) - 1
    acf_x, confint = acf(x, nlags=nlags, alpha=alpha)
    lags, acf_x, confint = [a[1:] for a in [lags, acf_x, confint]]
    axs[2].vlines(lags, [0], acf_x, **kwargs)
    axs[2].axhline(y=0, color='k', **kwargs)
    axs[2].fill_between(lags, confint[:,0] - acf_x, confint[:,1] - acf_x, alpha=0.25)

    return axs

def load_style(style_name):
    # Defaults:
    color_map = {
                'control': 'steelblue',
                'acarbose': 'goldenrod',
                'male': 'skyblue',
                'female': 'palevioletred',
                'C2013': 'violet',
                'TJL': np.array([0.93834679, 0.00947328, 0.49879277, 1.]),
                'UT': 'forestgreen',
                'UM': np.array([0.22526721, 0.42026913, 0.68868898, 1.]),
                }

    mark_map = {
                'C2013': '<',
                'C2014': 'o',
                'C2012': 'd',
                'Glenn': '8',
                'control': 'o',
                'acarbose': '<',
                'inulin': 'd',
                'male': '^',
                'female': 'o',
                'TJL': 'd',
                'UT': '<',
                'UM': 'o',
                }

    def assign_significance_symbol(pvalue):
        if pvalue >= 0.10:
            return ''
        elif pvalue >= 0.05:
            return '†'
        elif pvalue >= 0.001:
            return '*'
        else:
            return '**'

    def savefig(fig, basename, fmts=['png', 'pdf'], **kwargs):
        final_kwargs = dict(dpi=300, bbox_inches='tight', pad_inches=0.01)
        final_kwargs.update(kwargs)
        for f in fmts:
            fig.savefig('{basename}.{f}'.format(basename=basename, f=f),
                        **final_kwargs)

    _mm_to_inch = lambda x: x / 25.4
    halfwidth = _mm_to_inch(85)
    fullwidth = _mm_to_inch(170)
    fullheight = _mm_to_inch(225)

    # Style name specific:
    if style_name == 'paper':
        sns.set_context('paper')
        sns.set_style('white')
        mpl.rcParams['savefig.dpi'] = 150
        mpl.rcParams['figure.dpi'] = 150
        mpl.rcParams['figure.figsize'] = (4, 3)
        mpl.rcParams['axes.linewidth'] = 0.5
        mpl.rcParams['xtick.direction'] = 'in'
        mpl.rcParams['ytick.direction'] = 'in'
        mpl.rcParams['xtick.major.size'] = 2
        mpl.rcParams['ytick.major.size'] = 2
        mpl.rcParams['xtick.major.width'] = 0.5
        mpl.rcParams['ytick.major.width'] = 0.5
    elif style_name == 'poster':
        sns.set_context('poster')
        sns.set_style('dark')
        mpl.rcParams['savefig.dpi'] = 250
        mpl.rcParams['figure.dpi'] = 250
        mpl.rcParams['figure.figsize'] = (4, 3)
        mpl.rcParams['axes.linewidth'] = 1.5
        mpl.rcParams['xtick.direction'] = 'in'
        mpl.rcParams['ytick.direction'] = 'in'
        mpl.rcParams['xtick.major.size'] = 4
        mpl.rcParams['ytick.major.size'] = 4
        mpl.rcParams['xtick.major.width'] = 1
        mpl.rcParams['ytick.major.width'] = 1
    else:
        raise NotImplementedError(f'Style "{style_name}" has not been implemented.')

    return { 'color_map': color_map
            , 'mark_map': mark_map
            , 'assign_significance_symbol': assign_significance_symbol
            , 'savefig': savefig
            , 'fullwidth': fullwidth
            , 'halfwidth': halfwidth
            , 'fullheight': fullheight
            }
