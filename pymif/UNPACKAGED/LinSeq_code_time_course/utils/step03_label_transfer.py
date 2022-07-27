import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

def load_scrnaseq(file_scrnaseq, hvg=False):
    """
    

    Parameters
    ----------
    file_scrnaseq : TYPE
        DESCRIPTION.

    Returns
    -------
    data : TYPE
        DESCRIPTION.
    gem : TYPE
        DESCRIPTION.
    genes : TYPE
        DESCRIPTION.
    n_scrna_cells : TYPE
        DESCRIPTION.
    n_scrna_genes : TYPE
        DESCRIPTION.

    """
    data = ad.read_h5ad(file_scrnaseq)
    if hvg:
        sc.pp.highly_variable_genes(data, min_mean=0.00125, max_mean=30, min_disp=0.15)
        data = data[:, data.var.highly_variable]
    data.var_names_make_unique()
    n_scrna_cells = len(data.obs)
    n_scrna_genes = len(data.var)
    
    return data, n_scrna_cells, n_scrna_genes

def visualize_marker_distribution_adatas(adata1, adata2, genes, gene_highlight=[], n_cols=5, figsize=(8,3.5), ylog=False):
    """
    

    Parameters
    ----------
    df : TYPE
        DESCRIPTION.
    adata : TYPE
        DESCRIPTION.
    gene_markers : TYPE
        DESCRIPTION.
    n_cols : TYPE, optional
        DESCRIPTION. The default is 5.
    figsize : TYPE, optional
        DESCRIPTION. The default is (8,3.5).

    Returns
    -------
    fig : TYPE
        DESCRIPTION.

    """
    
    assert all([i in adata1.var_names for i in genes]), 'A gene is not present in dataset1'
    assert all([i in adata2.var_names for i in genes]), 'A gene is not present in dataset2'
    assert all([i in genes for i in gene_highlight]), 'A highlight gene is not present in the list of genes shown!'
    

    ### visualize distributions
    nrows = (len(genes)-1)//n_cols + 1
    fig, ax = plt.subplots(nrows=nrows, ncols=n_cols, sharey=True, figsize=figsize)
    if nrows>1:
        ax=ax.flatten()
    for i in range(len(genes), len(ax)):
        fig.delaxes(ax[i])
    fig.subplots_adjust(left=0.1, right=0.95, bottom=0.01, wspace=0, hspace=1.0, top=0.97)
    for i, gene in enumerate(genes):
        
        scrna1_gene = np.array(adata1[:,gene].X)[:,0].copy()
        # loc_gene = (loc_gene-np.min(loc_gene))/(np.max(loc_gene)-np.min(loc_gene))
        
        scrna2_gene = np.array(adata2[:,gene].X)[:,0].copy()
        # scrna_gene = (scrna_gene-np.min(scrna_gene))/(np.max(scrna_gene)-np.min(scrna_gene))
        
        ax[i].hist([scrna1_gene, scrna2_gene], bins=10, label=['scRNAseq1', 'scRNAseq2'], density=False)
        if ylog:
            ax[i].set_yscale('log')
        color = 'black'
        if gene in gene_highlight:
            color = 'red'
        ax[i].set_title(gene, color=color)
        # ax[i].set_xlabel('Expression\n level')
        # ax[i].set_xticks([])
    # ax[0].set_ylabel('Cell count')
    ax[0].legend(loc='upper right', frameon=False, fontsize=6)
    
    return fig
