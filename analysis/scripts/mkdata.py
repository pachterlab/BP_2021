#!/usr/bin/env python3

# Running this script from scripts folder
# for d in ../../data/kallisto_out/*; do DS=$(echo $d | rev | cut/da -d"/" -f1 | rev) && ./mkdata.py -d $DS -o ../../data/plotting/$DS; done

import sys
import argparse

from kb_python.utils import import_matrix_as_anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import pickle
import datetime
import time

# for d in ../../data/kallisto_out/*; do  DS=$(echo $d | cut -d'/' -f5) && echo $DS && ./mkdata.py -o ../../data/plotting/$DS -d $DS; done


from sklearn.decomposition import TruncatedSVD
from sklearn.manifold import TSNE

from sklearn.metrics.pairwise import manhattan_distances


REF = {
    "arabidopsis-SRR8257100_v2": "arabidopsis-tair10",
    "fly-SRR8513910_v2": "fly-dm6",
    "human_mouse-hgmm1k_v2": "human_mouse-hg19_mm10",
    "human-SRR8327928_v2": "human-grch38",
    "human-SRR8524760_v2": "human-grch38",
    "mouse-EMTAB7320_v2": "mouse-mm10",
    "mouse-heart1k_v2": "mouse-mm10",
    "mouse-SRR6998058_v2": "mouse-mm10",
    "mouse-SRR8206317_v2": "mouse-mm10",
    "mouse-SRR8599150_v2": "mouse-mm10",
    "mouse-SRR8639063_v2": "mouse-mm10",
    "rat-SRR7299563_v2": "rat-rnor6",
    "worm-SRR8611943_v2": "worm-ws260",
    "zebrafish-SRR6956073_v2": "zebrafish-dr82",
    "human_mouse-hgmm10k_v3": "human_mouse-hg19_mm10",
    "human_mouse-hgmm1k_v3": "human_mouse-hg19_mm10", # redo
    "human-pbmc10k_v3": "human-grch38",
    "human-pbmc1k_v3": "human-grch38",
    "mouse-heart1k_v3": "mouse-mm10",
    "mouse-neuron10k_v3": "mouse-mm10",
}

def nd(arr):
    return np.asarray(arr).reshape(-1)

def basic_process(A):
    '''
    sum counts per cell
    count ngenes per cell
    filter out cells with zero genes
    '''
    adata = A.copy()
    adata.obs['counts'] = nd(adata.X.sum(axis=1))
    adata.obs['ngenes'] = nd((adata.X > 0).sum(axis=1))
    adata = adata[adata.obs['counts'] > 0]
    adata.layers['log1p'] = np.log1p(adata.X)
    print(adata)

    return adata

def MA(AX, BX):
    '''
        Computes MA for MA plot
        X: A.var["gene_count"]
        Y: B.var["gene_count"]
    '''
    X = nd(AX.mean(axis=0))
    Y = nd(BX.mean(axis=0))

    M_AB = np.log2(X + 1) - np.log2(Y + 1)
    A_AB = 0.5*(np.log2(X + 1) + np.log2(Y + 1))
    return A_AB, M_AB

# Correlations
def _sparse_M_std(X):
    n = X.shape[0]
    return np.sqrt(n * X.multiply(X).sum(0) - np.multiply(X.sum(0), X.sum(0)))

def sparse_M_corr(X, Y):
    '''
        Computes Pearson correlation between X and Y (both in sparse format). Must be same shape.
        X: A_raw[common_obs.index].layers['log1p'] # raw
        Y: B_raw[common_obs.index].layers['log1p']# raw
        X: A.layers['log1p'] # filtered
        Y: B.layers['log1p'] # filtered
        Notes: I changed the axis in sum and shape, need to check if right
    '''
    X_std = _sparse_M_std(X)
    Y_std = _sparse_M_std(Y)
    XY_std = np.multiply(X_std, Y_std)
    n = X.shape[0]
    XY_cov = n*X.multiply(Y).sum(0) - np.multiply(X.sum(0), Y.sum(0))
    R = np.divide(XY_cov, XY_std)
    return np.squeeze(np.asarray(R))

def compute_tsvd(X):
    tsvd = TruncatedSVD(n_components=10)
    Y = tsvd.fit_transform(X)
    return Y

def compute_tsne(X):
    tsne = TSNE(perplexity=30, metric="euclidean", n_jobs=10, random_state=42, n_iter=750 )
    Y = tsne.fit_transform(X)
    return Y

def l1_dist(X, Y):
    '''
    computes manhattan distance between corresponding cells, and nearest cells
        X: A.layers['log1p']
        Y: B.layers['log1p']
    '''
    dist_AA = manhattan_distances(X, X)
    dist_AB = manhattan_distances(X, Y)

    # nkc are the kallisto-alevin distances
    dist_AB = np.diagonal(dist_AB)

    # ncc are the kallisto-kallisto distances
    AA = []
    for row in dist_AA:
        val = np.partition(row, 1)[1]
        AA.append(val)
    dist_AA = AA

    return dist_AA, dist_AB

def main(ds, out_dir):
    start = time.time()
    # load cell ranger barcodes
    t = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{t}] load cr barcodes")
    
    
    cellranger_barcodes_file = f"../../data/cellranger_barcodes/{ds}.txt"
    c = pd.read_csv(cellranger_barcodes_file, header=None, names=['bcs'])
    c = c.bcs.apply(lambda x: x.split('-')[0]).values
    
    # load whitelists
    kb_wl_file = f"../../data/kallisto_out/{ds}/whitelist.txt"    
    al_wl_file = f"../../data/alevin_out/{ds}/whitelist/permit_freq.tsv"
    
    k = nd(pd.read_csv(kb_wl_file, header=None).values)
    a = nd(pd.read_csv(al_wl_file, header=None, sep="\t", names = ['bcs', 'cnt'])['bcs'].values)
    
    cr_barcodes = np.intersect1d(np.intersect1d(k, a), c)
    
    
    # load alevin data
    t = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{t}] load alevin mtx")
    
    alevin_path = f"../../data/alevin_out/{ds}/quant/alevin/"
    
    alevin_raw_decoys = import_matrix_as_anndata(
        barcodes_path=os.path.join(alevin_path, "quants_mat_rows.txt"), 
        genes_path=os.path.join(alevin_path, "quants_mat_cols.txt"), 
        matrix_path=os.path.join(alevin_path, "quants_mat.mtx"))
    
    # load decoys
    t = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{t}] load alevin decoys")

    ref = REF[ds]
    decoys = pd.read_csv(f"../../reference/{ref}/salmon/decoys.txt", header = None, names=['decoys'])
    alevin_raw = alevin_raw_decoys[:, ~alevin_raw_decoys.var.index.isin(decoys['decoys'])].copy()
    
    # load kb
    t = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{t}] load kb mtx")
    
    kb_path = f"../../data/kallisto_out/{ds}/count/"
    kb_raw = import_matrix_as_anndata(
        barcodes_path=os.path.join(kb_path, "output.barcodes.txt"), 
        genes_path=   os.path.join(kb_path, "output.genes.txt"), 
        matrix_path=  os.path.join(kb_path, "output.mtx"))
    
    # basic process
    t = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{t}] basic process")
    
    kb = basic_process(kb_raw)
    alevin = basic_process(alevin_raw)
    
    # common barcodes and genes
    t = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{t}] common barcodes and genes")

    kb_barcodes = kb.obs.index.values
    alevin_barcodes = alevin.obs.index.values
    common_barcodes = np.intersect1d(kb_barcodes, alevin_barcodes)
    
    kb_genes = kb.var.index.values
    alevin_genes = alevin.var.index.values
    common_genes = np.intersect1d(kb_genes, alevin_genes)
    
    # common matrices --> SAVE
    kb_common = kb[common_barcodes][:, common_genes]
    kb_common.write_h5ad(os.path.join(out_dir, "kb_common.h5ad"))
    
    alevin_common = alevin[common_barcodes][:, common_genes]
    alevin_common.write_h5ad(os.path.join(out_dir, "alevin_common.h5ad"))
    
    # filtered by CR --> SAVE AFTER TSNE
    kb_common_cr     = kb_common[cr_barcodes]
    alevin_common_cr = alevin_common[cr_barcodes]
    
    # MA data --> SAVE
    t = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{t}] ma")
    
    A_AB, M_AB = MA(kb_common_cr.X, alevin_common_cr.X)
    
    with open(os.path.join(out_dir, "A_AB.pkl"), 'wb') as fp:
        pickle.dump(A_AB, fp)
    with open(os.path.join(out_dir, "M_AB.pkl"), 'wb') as fp:
        pickle.dump(M_AB, fp)
    
    # correlations --> SAVE
    t = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{t}] correlation")
    
    cc_raw = sparse_M_corr(
        kb_common.layers['log1p'].T, 
        alevin_common.layers['log1p'].T)

    with open(os.path.join(out_dir, "cc_raw.pkl"), 'wb') as fp:
        pickle.dump(cc_raw, fp)
    
    # TSVD
    t = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{t}] tsvd")
    
    kb_common_cr.obsm['TSVD']     = compute_tsvd(kb_common_cr.layers['log1p'])
    alevin_common_cr.obsm['TSVD'] = compute_tsvd(alevin_common_cr.layers['log1p'])
    
    # TSNE
    t = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{t}] tsne")
    
    kb_common_cr.obsm['TSNE']     = compute_tsne(kb_common_cr.obsm['TSVD'])
    kb_common_cr.write_h5ad(os.path.join(out_dir, "kb_common_cr.h5ad"))
    
    alevin_common_cr.obsm['TSNE'] = compute_tsne(alevin_common_cr.obsm['TSVD'])
    alevin_common_cr.write_h5ad(os.path.join(out_dir, "alevin_common_cr.h5ad"))
    
    # L1 DIST --> SAVE
    t = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{t}] l1 distance")
    dist_AA, dist_AB = l1_dist(kb_common_cr.layers['log1p'], alevin_common_cr.layers['log1p'])
    dist_BB, dist_BA = l1_dist(alevin_common_cr.layers['log1p'], kb_common_cr.layers['log1p'])

    with open(os.path.join(out_dir, "dist_AA.pkl"), 'wb') as fp:
        pickle.dump(dist_AA, fp)
    with open(os.path.join(out_dir, "dist_AB.pkl"), 'wb') as fp:
        pickle.dump(dist_AB, fp)
    with open(os.path.join(out_dir, "dist_BA.pkl"), 'wb') as fp:
        pickle.dump(dist_BA, fp)
    with open(os.path.join(out_dir, "dist_BB.pkl"), 'wb') as fp:
        pickle.dump(dist_BB, fp)

    t = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')   
    print(f"[{t}] Done, took {(time.time()-start)/60} min.")
    
    return None

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="processes all data")
    parser.add_argument("-d", help="dataset name")
    parser.add_argument("-o", help="output dir")
    
    args = parser.parse_args()

    if not os.path.exists(args.o):
        os.makedirs(args.o)
    
    main(args.d, args.o)