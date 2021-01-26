#!/usr/bin/env python3

# Running this script from scripts folder
# for d in ../../data/plotting/*; do DS=$(echo $d | rev | cut -d"i ./" -f1 | rev) && echo $DS && ./mkplot.py -d $DS -i ../../data/plotting/$DS -o ./delete; done

import sys
import argparse
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import anndata
import numpy as np
import pandas as  pd
import os
import pickle
import warnings
warnings.filterwarnings("ignore")

markersize = 20
alpha=1
linewidth=5
xmax = 10**6

kallisto_color = "#377eb8"
alevin_color = "#e41a1c"
fsize=20
gridalpha = 0.2

plt.rcParams.update({'font.size': fsize})

from matplotlib.ticker import LogLocator, NullFormatter

def fix_ticks(ax, axis=['x']):
    locmin = LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=25)
    locmaj = LogLocator(base=10,numticks=25)
    if 'x' in axis:
        ax.xaxis.set_minor_locator(locmin)
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.xaxis.set_major_locator(locmaj)
    if 'y' in axis:
        ax.yaxis.set_minor_locator(locmin)
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.set_major_locator(locmaj)
    return ax

def nd(arr):
    return np.asarray(arr).reshape(-1)

def _plt_color(lst):
    cols=[]
    for l in lst:
        if l>0.25 or l<-0.25:
            cols.append("red")
        elif l<=0.25 and l>=-0.25:
            cols.append('black')
    return cols

import matplotlib.patches as mpatches

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

def yex(ax):
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    return ax

from mpl_toolkits.axes_grid1 import make_axes_locatable

def make_hist(A, B, orientation="vertical", ax=None):
    hist, concat_bins = np.histogram(np.concatenate((A,B)), bins='auto')
    hist, A_bins =  np.histogram(A, bins='auto')
    hist, B_bins =  np.histogram(B, bins='auto')

    best_bins = min([A_bins,concat_bins,B_bins], key=len) # may need to change to max
    best_bins = concat_bins
    ax.hist(A, bins=best_bins, orientation=orientation, color=kallisto_color, label="kallisto", alpha=1)
    ax.hist(B, bins=best_bins, orientation=orientation, color=alevin_color, label="alevin", alpha=1)
    return ax

def make_scatter_hist(dist_AA, dist_BB, dist_AB, dist_BA, ax=None):
    x = dist_AA
    y = dist_AB

    xx = dist_BA
    yy = dist_BB

    # the scatter plot:
    ax.scatter(x, y, label="kallisto", color=kallisto_color, s=markersize)
    ax.scatter(xx, yy, label="alevin", color=alevin_color, s=markersize)
    ax.set_aspect(1.)

    # create new axes on the right and on the top of the current axes
    # The first argument of the new_vertical(new_horizontal) method is
    # the height (width) of the axes to be created in inches.
    divider = make_axes_locatable(ax)
    axHistx = divider.append_axes("top", 1.5, pad=0.075, sharex=ax)
    axHisty = divider.append_axes("right", 1.5, pad=0.075, sharey=ax)

    # make some labels invisible
    axHistx.xaxis.set_tick_params(labelbottom=False)
    axHisty.yaxis.set_tick_params(labelleft=False)


    ## Right histogram,  cellranger-cellranger, kallisto-cellranger,
    axHisty = make_hist(y, yy, orientation="horizontal", ax=axHisty)

    # kallisto-kallisto, cellranger-kallisto
    axHistx = make_hist(x, xx, orientation="vertical", ax=axHistx)


    # the xaxis of axHistx and yaxis of axHisty are shared with axScatter,
    # thus there is no need to manually adjust the xlim and ylim of these
    # axis.
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    ax.set_xlabel("$\ell_1$ to nearest kallisto", fontsize=fsize)
    ax.set_ylabel("$\ell_1$ to nearest alevin", fontsize=fsize)
    axHistx.set_ylabel("Barcode counts", fontsize=fsize-8)
    axHisty.set_xlabel("Barcode counts", fontsize=fsize-8)

    axHistx.set_title("E.1", fontweight='bold', loc = 'left' )

    axHistx.legend(fontsize=fsize-5, loc="upper right")
    return axHistx

def QQ_hgmm(path, ax):
    dataset_shortname = os.path.basename(path).split(".")[0]
    df = pd.read_csv(path)
    df.ontology = df.ontology.astype("category")

    from sklearn import preprocessing
    le = preprocessing.LabelEncoder()

    from matplotlib.lines import Line2D
    hg = df[df.mapping.str.contains("Hs")]
    mm = df[df.mapping.str.contains("Mm")]
    legend_elements = [Line2D([0], [0], marker='o', color="w",alpha=0.4, label='Human',markerfacecolor='k', markersize=10),
                      Line2D([0], [0], marker='s', color='w', alpha=0.2, label='Mouse', markerfacecolor='grey', markersize=10)]

    #c = le.fit_transform(df.ontology.values)

    c = le.fit_transform(hg.ontology.values)
    scatter = ax.scatter(hg.uniform_log, hg.p_log, c=c)

    c = le.fit_transform(mm.ontology.values)
    ax.scatter(mm.uniform_log, mm.p_log, c=c, marker='s')

    l1 = ax.legend(handles=legend_elements, loc="lower right", title="Species", fontsize=fsize-5, title_fontsize=fsize-5)

    ax.plot(hg.uniform_log, hg.upper, color='k', alpha=0)
    ax.plot(hg.uniform_log, hg.lower, color='k', alpha=0)

    ax.plot(mm.uniform_log, mm.upper, color='grey', alpha=0)
    ax.plot(mm.uniform_log, mm.lower, color='grey', alpha=0)

    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    ax.fill_between(hg.uniform_log, hg.lower, hg.upper, color='black', alpha='0.4')
    ax.fill_between(mm.uniform_log, mm.lower, mm.upper, color='grey', alpha='0.2')



    # Produce a legend for the ranking (colors). Even though there are 40 different
    # rankings, we only want to show 5 of them in the legend.

    l2 = ax.legend(*(scatter.legend_elements()[0], le.classes_), loc="upper left", title="Ontology", fontsize=fsize-5, title_fontsize=fsize-5)
    ax.add_artist(l1)


    ax.set_xlabel("Expected -log$_{10}$(p)", fontsize=fsize)
    ax.set_ylabel("Observed -log$_{10}$(p)", fontsize=fsize)
    df = df[df.label.astype(str).values != 'nan']
    ax.set_title("QQ", fontweight='bold', fontsize = fsize, loc = 'left' )
    return ax

def QQ_plot(path, ax):
    dataset_shortname = os.path.basename(path).split(".")[0]
    
    if "mm" in dataset_shortname:
        return QQ_hgmm(path, ax)
    
    df = pd.read_csv(path)
    df.ontology = df.ontology.astype("category")

    from sklearn import preprocessing
    le = preprocessing.LabelEncoder()


    c = le.fit_transform(df.ontology.values)

    scatter = ax.scatter(df.uniform_log, df.p_log, c=c)


    ax.plot(df.uniform_log, df.upper, color='k', alpha=0)
    ax.plot(df.uniform_log, df.lower, color='k', alpha=0)

    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    ax.fill_between(df.uniform_log, df.lower, df.upper, color='black', alpha=0.4)

    # Produce a legend for the ranking (colors). Even though there are 40 different
    # rankings, we only want to show 5 of them in the legend.

    l2 = ax.legend(*(scatter.legend_elements()[0], le.classes_), loc="best", title="Ontology", fontsize=fsize-5, title_fontsize=fsize-5)


    ax.set_xlabel("Expected -log$_{10}$(p)", fontsize=fsize)
    ax.set_ylabel("Observed -log$_{10}$(p)", fontsize=fsize)
    df = df[df.label.astype(str).values != 'nan']
    base = os.path.dirname(path)
    with open(os.path.join(base, f"GO_{dataset_shortname}.txt"), 'w') as f:


        for i, txt in enumerate(df["rank"]):
            f.write(str(i+1) + "\t" + df.label.values[i] + "\n")
            print(df.label.values[i])
            ax.annotate(txt, (df.uniform_log[i], df.p_log[i]), fontsize=15)
    return ax

def DE_plot(bar_path, go_terms_path, ax=None):


    fold_change_change_colors = {
    '(4,5]': lighten_color(kallisto_color, 1.4),
    '(3,4]': lighten_color(kallisto_color, 1.1),
    '(2,3]': lighten_color(kallisto_color, 0.8),
    '(1,2]': lighten_color(kallisto_color, 0.5),
    '(0,1]': lighten_color(kallisto_color, 0.2),
    '(-1,0]': lighten_color(alevin_color, 0.2),
    '(-2,-1]':lighten_color(alevin_color, 0.5),
    '(-3,-2]':lighten_color(alevin_color, 0.8),
    '(-4,-3]':lighten_color(alevin_color, 1.1),
    '(-5,-4]':lighten_color(alevin_color, 1.4),
    }
    
    df = pd.read_csv(bar_path)
    go = pd.read_csv(go_terms_path)
    go.index = go['GO.ID']

    df['Term'] = df['GO.ID'].map(go['GO_full_name']).astype(str)

    # For some datasets there are no DE genes, so we need this check to just write a text and make not plot
    df = df[~df.Term.str.contains('nan')]
    genesets = np.sort(df.Term.unique())
    
    if len(genesets)==0:
        ax.text(0.5*(1), 0.5*(1), 'No significant gene sets found ',
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=20, color='black',
            transform=ax.transAxes)
        ax.axis('off')

    # If there are DE genes, then we make a plot
    elif len(genesets)>0:
        gb = (df.groupby(["change", "Term"])[['gene']].nunique())
        new_df = gb.unstack(fill_value=0).sort_index(level=0, key=lambda x: [int(i[1:].split(',')[0]) for i in x])['gene']

        itvs = new_df.index.values
        names = new_df.columns.values
        mtx = new_df.values

        ind = np.arange(names.shape[0])
        bottom = np.zeros_like(names)

        for itv, row in zip(itvs, mtx):
            ax.bar(ind, row, bottom = bottom, label = itv,
                    alpha = 0.5, width=0.8, color = fold_change_change_colors[itv])
            bottom+= row

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], title='Fold change')

        labels = [i.replace(" ", "\n") for i in names]

        ax.set(**{
            "xticks": ind,
            "ylabel": "Number of genes",
            "xlabel": "Gene set"
        })
        ax.set_xticklabels(labels, fontsize=10)

        for label in ax.get_xmajorticklabels():
            label.set_horizontalalignment("center")

        ax.grid(color='dimgrey', linestyle='-', linewidth=0.5, which="both", alpha = 0.5)
    return ax


def load_data(dir):

    kb_common = anndata.read_h5ad(os.path.join(dir, 'kb_common.h5ad'))
    alevin_common = anndata.read_h5ad(os.path.join(dir, 'alevin_common.h5ad'))
    kb_common_cr = anndata.read_h5ad(os.path.join(dir, 'kb_common_cr.h5ad'))
    alevin_common_cr = anndata.read_h5ad(os.path.join(dir, 'alevin_common_cr.h5ad'))

    with open(os.path.join(dir, 'cc_raw.pkl'), 'rb') as fp:
        cc_raw = pickle.load(fp)
    with open(os.path.join(dir, 'A_AB.pkl'), 'rb') as fp:
        A_AB = pickle.load(fp)
    with open(os.path.join(dir, 'M_AB.pkl'), 'rb') as fp:
        M_AB = pickle.load(fp)
    with open(os.path.join(dir, 'dist_AA.pkl'), 'rb') as fp:
        dist_AA = pickle.load(fp)
    with open(os.path.join(dir, 'dist_AB.pkl'), 'rb') as fp:
        dist_AB = pickle.load(fp)
    with open(os.path.join(dir, 'dist_BA.pkl'), 'rb') as fp:
        dist_BA = pickle.load(fp)
    with open(os.path.join(dir, 'dist_BB.pkl'), 'rb') as fp:
        dist_BB = pickle.load(fp)
        
    return (kb_common, alevin_common, kb_common_cr, alevin_common_cr, cc_raw, A_AB, M_AB, dist_AA, dist_AB, dist_BA, dist_BB)

def main(ds, indir, outdir):
    dataset_name = ds
    
    kb_common, alevin_common, kb_common_cr, alevin_common_cr, cc_raw, A_AB, M_AB, dist_AA, dist_AB, dist_BA, dist_BB = load_data(indir)
    
    fig = plt.figure(figsize=(25, 25),constrained_layout=False)
    st = fig.suptitle(dataset_name, fontweight='bold', x=0.4)
    gs = fig.add_gridspec(8, 4)

    ax_a = fig.add_subplot(gs[0:2, 0])
    ax_b = fig.add_subplot(gs[2:4, 0])
    ax_c = fig.add_subplot(gs[4:6, 0])
    ax_d = fig.add_subplot(gs[6:8, 0])

    ax_e_left   = fig.add_subplot(gs[0:2, 1])
    ax_e_right  = fig.add_subplot(gs[0:2, 2])

    ax_f_left   = fig.add_subplot(gs[2:4, 1])
    ax_f_right  = fig.add_subplot(gs[2:4, 2])
    ax_g_left   = fig.add_subplot(gs[4:6, 1])
    ax_g_right  = fig.add_subplot(gs[4:6, 2])

    ax_h        = fig.add_subplot(gs[6:8, 1:3])


    ############################################################################## KNEE plot
    ax = ax_a

    title = "A"

    x = np.sort(nd(kb_common_cr.X.sum(axis=1)))[::-1]
    y = np.arange(x.shape[0])
    
    minkb = min(x)

    ax.plot(x, y, color=kallisto_color, label="kallisto", linewidth=linewidth)

    # ## Alevin
    x = np.sort(nd(alevin_common_cr.X.sum(axis=1)))[::-1]
    y = np.arange(x.shape[0])
    
    minal = min(x)

    cutoff = min(minkb, minal)

    ax.axvline(x=cutoff, color="grey", linestyle="--")

    ax.plot(x, y, color = alevin_color, label="alevin", linewidth=linewidth, zorder=-1)

    ax.set(**{
        'xscale': 'log',
        'xlim': (1, xmax),
        'ylim': 1,
        'yscale': 'log',
        'xlabel': "kallisto UMI counts",
        'ylabel': "Cumulative number of barcodes"
     })

    fix_ticks(ax, ['x', 'y'])

    # Customize the major grid
    ax.grid(which='both', linestyle='-', linewidth='0.5', color='dimgrey', alpha=0.2)
    
    ax.set_title(title, fontweight='bold', loc = 'left' )
    ax.legend()

    ############################################################################## PSEUDOBULK
    ax = ax_b

    title = "B"

    x = kb_common.obs['counts'].values
    y = alevin_common.obs['counts'].values

    cutoff_mask = x>=cutoff

    xx = x[cutoff_mask]
    yy = y[cutoff_mask]
    ax.scatter(xx,yy, alpha = 0.05, s=markersize, color="k")

    xx = x[~cutoff_mask]
    yy = y[~cutoff_mask]
    ax.scatter(xx,yy, alpha = 0.05, s=markersize, color="lightgrey")


    ax.axvline(x=cutoff, color="grey", linestyle="--")


    ax.set(**{
        "xlim": (1, xmax),
         'xscale': 'log',
         'yscale': 'log',
        'xlabel': 'kallisto UMI counts',
        'ylabel': 'alevin UMI counts'
    })
    yex(ax)

    fix_ticks(ax, ['x', 'y'])

    # Customize the major grid
    ax.grid(which='both', linestyle='-', linewidth='0.5', color='dimgrey', alpha=0.2)

    ax.set_title(title, fontweight='bold', loc = 'left' )

    ############################################################################## GENES DETECTED
    ax = ax_c

    title = "C"

    x = kb_common.obs['counts'].values
    y = kb_common.obs['ngenes'].values

    xx = x[cutoff_mask]
    yy = y[cutoff_mask]
    ax.scatter(xx,yy, alpha =0.05, s=markersize, label='kallisto', color=kallisto_color)

    xx = x[~cutoff_mask]
    yy = y[~cutoff_mask]
    ax.scatter(xx,yy, alpha = 0.05, s=markersize, color=lighten_color(kallisto_color))

    ax.axvline(x=cutoff, color="grey", linestyle="--")

    x = alevin_common.obs['counts'].values
    y = alevin_common.obs['ngenes'].values

    xx = x[cutoff_mask]
    yy = y[cutoff_mask]
    ax.scatter(xx,yy, alpha =0.05, s=markersize, label='alevin', color=alevin_color, zorder = -1)

    xx = x[~cutoff_mask]
    yy = y[~cutoff_mask]
    ax.scatter(xx,yy, alpha = 0.05, s=markersize, color=lighten_color(alevin_color), zorder=-1)

    ax.set(**{
        'xscale': 'log',
        'yscale': 'log',
        'xlim': (1, xmax),
        'ylim': (1),
        'xlabel': 'UMI counts',
        'ylabel': 'Genes detected'
    })

    A_patch = mpatches.Patch(color=kallisto_color, label="kallisto")
    B_patch = mpatches.Patch(color=alevin_color, label="alevin")
    ax.legend(handles=[A_patch, B_patch], loc="lower right")

    fix_ticks(ax, ['x', 'y'])

    # Customize the major grid
    ax.grid(which='both', linestyle='-', linewidth='0.5', color='dimgrey', alpha=0.2)

    ax.set_title(title, fontweight='bold', loc = 'left' )

    ############################################################################## CORRELATION
    ax = ax_d

    title = "D"

    x = kb_common.obs['counts'].values
    y = cc_raw

    xx = x[cutoff_mask]
    yy = y[cutoff_mask]
    ax.scatter(xx,yy, s=markersize, alpha=0.05, color="k")

    xx = x[~cutoff_mask]
    yy = y[~cutoff_mask]
    ax.scatter(xx,yy, s=markersize, alpha=0.05, color="lightgray")

    ax.axvline(x=cutoff, color="grey", linestyle="--")

    ax.set(**{
        'xscale': 'log',
        "xlim": (1, xmax),
        'ylim': (0,1),
        'xlabel': 'kallisto UMI counts',
        'ylabel': 'Pearson Correlation'
    })

    fix_ticks(ax, ['x'])

    # Customize the major grid
    ax.grid(which='both', linestyle='-', linewidth='0.5', color='dimgrey', alpha=0.2)
    
    ax.set_title(title, fontweight='bold', loc = 'left' )

    ############################################################################## l1
    ax = ax_e_left

    title = "E.1"

    make_scatter_hist(dist_AA, dist_BB, dist_AB, dist_BA, ax=ax)

    ax.set_title(title, fontweight='bold', loc = 'left' )


    #############
    ax = ax_e_right

    title = "E.2"

    # cnts = nd(kb_common_cr.layers["log1p"].sum(1))
    cnts = nd(kb_common_cr.obs['counts'].values)
    ax.scatter(cnts, dist_AB, color=alevin_color, alpha=alpha, label="alevin", s=markersize)
    
    ax.scatter(cnts, dist_AA, color=kallisto_color, alpha=alpha, label="kallisto", s=markersize)


    ax.set(**{
        "xlabel": "kallisto UMI counts",
        "ylabel": "$\ell_1$ distance",
        'xscale': 'log',
        'yscale': 'log',
#        'xlim': (1, xmax)
    })

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='upper left', markerscale=4)
    
    fix_ticks(ax, ['x'])
    ax.grid(which='both', linestyle='-', linewidth='0.5', color='dimgrey', alpha=0.2)

    ax.set_title(title, fontweight='bold', loc = 'left' )




    ############################################################################## TSNE
    ax = ax_f_left

    title = "F.1"

    label = "kallisto"
    color = kallisto_color
    x = kb_common_cr.obsm['TSNE'][:,0]
    y = kb_common_cr.obsm['TSNE'][:,1]

    ax.scatter(x, y, s =markersize, c = color, alpha = 1, edgecolors = 'none', label = label)

    ax.set(**{
        "xlabel": "t-SNE 1",
        "ylabel": "t-SNE 2",
        "xticklabels": [],
        "yticklabels": [],

    })

    ax.set_title(title, fontweight='bold', fontsize = fsize, loc = 'left' )
    ax.tick_params(axis=u'both', which=u'both',length=0)
    ax.legend(fontsize=fsize, markerscale=4, loc="upper left")

    ################
    ax = ax_f_right

    title = "F.2"

    label = "alevin"
    color = alevin_color

    x = alevin_common_cr.obsm['TSNE'][:,0]
    y = alevin_common_cr.obsm['TSNE'][:,1]

    ax.scatter(x, y, s =markersize, c = color, alpha = 1, edgecolors = 'none', label = label)

    ax.set(**{
        "xlabel": "t-SNE 1",
        "ylabel": "t-SNE 2",
        "xticklabels": [],
        "yticklabels": [],

    })

    ax.set_title(title, fontweight='bold', fontsize = fsize, loc = 'left' )
    ax.tick_params(axis=u'both', which=u'both',length=0)
    ax.legend(fontsize=fsize, markerscale=4, loc="upper left")


    ############################################################################## MA
    ax = ax_g_left

    title = "G.1"

    cols = _plt_color(M_AB)
    ax.scatter(A_AB, M_AB, alpha=0.05, c=cols, s=markersize)

    ax.set(**{
        "ylabel": "log$_2$(Fold change)",
        "xlabel": "log$_2$(Average gene count)",
        "ylim": (-5, 5)
    })


    ax.set_title(title, fontweight='bold', loc = 'left' )

    A_patch = mpatches.Patch(color=kallisto_color, label="kallisto")
    B_patch = mpatches.Patch(color=alevin_color, label="alevin")
    same = mpatches.Patch(color='white', label='log$_2$(FC)$\leq$ 0.25 ({:.3f})'.format(M_AB[M_AB<=0.25].shape[0]/M_AB.shape[0]))
    ax.arrow(0, 1, 0, 1.5, length_includes_head=True, width=.05, color=kallisto_color)
    ax.arrow(0, -1, 0, -1.5, length_includes_head=True, width=.05, color=alevin_color)
    ax.legend(handles=[A_patch, B_patch, same], fontsize=fsize, loc="upper right")

    ################################################################################# QQ plot
    ax = ax_g_right

    title = "G.2"
    QQ_plot(f"../../data/gsea_qq/{dataset_name}.csv", ax)
    ax.set_title(title, fontweight='bold', loc='left')

    ################################################################################ GSEA Plot
    ax = ax_h

    REF_map = {
        "mouse": "Mus musculus",
        "human": "Homo sapiens",
        'rat': "Rattus norvegicus",
        "arabidopsis": "Arabidopsis thaliana",
        "fly": "Drosophila melanogaster",
        "worm": "Caenorhabditis elegans",
        "zebrafish": "Danio rerio",
        'human_mouse': "Homo sapiens Mus musculus"
    }

    species = ds.split('-')[0]
    species = "_".join(REF_map.get(species).split(" "))

    title = "H"
    DE_plot(f"../../data/gsea_bar/{ds}.csv", f"../../reference/GO/go_def{species}.csv", ax)

    ax.set_title(title, fontweight='bold', loc='left')

    plt.tight_layout()
    st.set_y(0.95)
    fig.subplots_adjust(top=0.925)
    
    plt.savefig(os.path.join(outdir, ds+".png"), dpi=300, bbox_inches='tight')
    return None

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="processes all data")
    parser.add_argument("-d", help="dataset name")
    parser.add_argument("-i", help="input dir with all plotting data")
    parser.add_argument("-o", help="output dir to save pic")
    
    
    args = parser.parse_args()

    if not os.path.exists(args.o):
        os.makedirs(args.o)
    
    main(args.d, args.i, args.o)