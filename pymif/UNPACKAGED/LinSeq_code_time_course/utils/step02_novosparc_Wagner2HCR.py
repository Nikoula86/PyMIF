import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist

import _GWadjusted as gwan
import _reconstruction as rec

from importlib import reload
reload(gwan)
reload(rec)

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

def visualize_marker_distribution(df, adata, gene_markers, n_cols=5, figsize=(8,3.5), ylog=False):
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
    
    assert all([i in df.keys() for i in gene_markers])
    assert all([i in adata.var_names for i in gene_markers])
    
    marker_hcr = df[gene_markers].to_numpy()

    ### visualize distributions
    nrows = (len(gene_markers)-1)//n_cols + 1
    fig, ax = plt.subplots(nrows=nrows, ncols=n_cols, sharey=True, figsize=figsize)
    if nrows>1:
        ax=ax.flatten()
    for i in range(len(gene_markers), len(ax)):
        fig.delaxes(ax[i])
    fig.subplots_adjust(left=0.1, right=0.95, bottom=0.2, wspace=0, hspace=0.5)
    for i, gene in enumerate(gene_markers):
        
        loc_gene = marker_hcr[:,i].copy()
        # loc_gene = (loc_gene-np.min(loc_gene))/(np.max(loc_gene)-np.min(loc_gene))
        
        scrna_gene = np.array(adata[:,gene].X)[:,0].copy()
        # scrna_gene = (scrna_gene-np.min(scrna_gene))/(np.max(scrna_gene)-np.min(scrna_gene))
        
        ax[i].hist([loc_gene, scrna_gene], bins=10, label=['locations', 'scRNAseq'], density=False)
        if ylog:
            ax[i].set_yscale('log')
        ax[i].title.set_text(gene)
        ax[i].set_xlabel('Expression\n level')
        # ax[i].set_ylim((0,100000))
    ax[0].set_ylabel('Cell count')
    ax[0].legend(loc='upper right', frameon=False, fontsize=6)
    
    return fig

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
    
    return fig, ax

def check_data_dimensions_ok(locations, data, num_locs_cells, n_scrna_cells, gene_markers):
    genes = data.var_names
    assert n_scrna_cells==len(data.obs), "Something is wrong with the scRNAseq data"
    assert all([i in genes for i in gene_markers]), "Not all markers in scRNAseq dataset"
    assert locations.shape[0]==num_locs_cells, "Something wrong with the locations"

    print('\nNumber of cells in marker dataset:',num_locs_cells)
    print('Number of cells in scRNAseq dataset:',n_scrna_cells)
    print('Number of genes in scRNAseq dataset:',len(genes))
    print("All seems good.\n")

def find_marker_expression_in_scrna_data(data, gene_markers):
    
    try:
        marker_scRNA = np.array(data[:,gene_markers].X.todense())
    except:
        marker_scRNA = np.array(data[:,gene_markers].X)
    
    return marker_scRNA

def setup_reconstruction(gem, locations, 
                         marker_scRNA, marker_gt, 
                         num_marker_cells, num_scRNAcells,
                         num_neighbors_source = 5, num_neighbors_target = 5):
    print('\nSetup reconstruction...\n')
    
    cost_expression, cost_locations = rec.setup_for_OT_reconstruction( 
                                                        gem,
                                                        locations,
                                                        num_neighbors_source = num_neighbors_source,
                                                        num_neighbors_target = num_neighbors_target 
                                                        )
    
    cost_marker_genes = cdist( 
                                marker_scRNA/np.amax(marker_scRNA),
                                marker_gt/np.amax(marker_gt)
                                    )

    # Distributions at target and source spaces
    p_locations, p_expression = rec.create_space_distributions(num_marker_cells, num_scRNAcells)
    
    return cost_expression, cost_locations, cost_marker_genes, p_locations, p_expression

def scrnaseq2locs(locs_df, gene_markers, adata, num_locs_cells, n_scrna_cells, alpha, 
                tol=10e-9, epsilon=1e-3,
                num_neighbors_source=5, num_neighbors_target=5,
                space_keys=["z", "y", "x"]):
    """
    

    Parameters
    ----------
    locs_df : TYPE
        DESCRIPTION.
    gene_markers : TYPE
        DESCRIPTION.
    adata : TYPE
        DESCRIPTION.
    num_locs_cells : TYPE
        DESCRIPTION.
    n_scrna_cells : TYPE
        DESCRIPTION.
    alpha : TYPE
        DESCRIPTION.

    Returns
    -------
    adata_mapped : TYPE
        DESCRIPTION.

    """
    locs_array = locs_df[space_keys].to_numpy()
    marker_gt = locs_df[gene_markers].to_numpy()

    try:
        gem = np.array(adata.X.todense())
    except:
        gem = np.array(adata.X)

    check_data_dimensions_ok(locs_array, adata, num_locs_cells, n_scrna_cells, gene_markers)
    
    ### find marker gene in scrna data and setup reconstruction
    marker_scrna = find_marker_expression_in_scrna_data(adata, gene_markers)
    cost_exp, cost_loc, cost_marker_genes, p_loc, p_exp = setup_reconstruction(gem, locs_array, 
                                                                    marker_scrna, marker_gt, 
                                                                    num_locs_cells, n_scrna_cells,
                                                                    num_neighbors_source = num_neighbors_source, 
                                                                    num_neighbors_target = num_neighbors_target)
    
    print('\nReconstructing with alpha='+str(alpha))
    gw = gwan.gromov_wasserstein_adjusted_norm(cost_marker_genes, 
                                            cost_exp, cost_loc,
                                            alpha, 
                                            p_exp, p_loc,
                                            'square_loss', 
                                            epsilon=epsilon,
                                            tol=tol, 
                                            verbose=True)

    print('Done.')
        
    ### plot
    print('Transporting scRNA values onto locations')
    try:
        gem = adata.X.todense()
    except:
        gem = adata.X
    sgem = np.dot(gem.T, gw)

    adata_mapped = ad.AnnData(X=sgem.T)
    adata_mapped.var_names = adata.var_names
    adata_mapped.var = adata.var
    adata_mapped.obs = locs_df
    
    return adata_mapped

def scrnaseq2locs_adatas(adata1, gene_markers, adata, num_locs_cells, n_scrna_cells, alpha, 
                tol=10e-9, epsilon=1e-3,
                num_neighbors_source=5, num_neighbors_target=5,
                space_keys=["x", "y", "z"]):
    """
    

    Parameters
    ----------
    locs_df : TYPE
        DESCRIPTION.
    gene_markers : TYPE
        DESCRIPTION.
    adata : TYPE
        DESCRIPTION.
    num_locs_cells : TYPE
        DESCRIPTION.
    n_scrna_cells : TYPE
        DESCRIPTION.
    alpha : TYPE
        DESCRIPTION.

    Returns
    -------
    adata_mapped : TYPE
        DESCRIPTION.

    """
    locs_array = adata1.obs[space_keys].to_numpy()
    marker_gt = adata1[:,gene_markers].X

    try:
        gem = np.array(adata.X.todense())
    except:
        gem = np.array(adata.X)

    check_data_dimensions_ok(locs_array, adata, num_locs_cells, n_scrna_cells, gene_markers)
    
    ### find marker gene in scrna data and setup reconstruction
    marker_scrna = find_marker_expression_in_scrna_data(adata, gene_markers)
    cost_exp, cost_loc, cost_marker_genes, p_loc, p_exp = setup_reconstruction(gem, locs_array, 
                                                                    marker_scrna, marker_gt, 
                                                                    num_locs_cells, n_scrna_cells,
                                                                    num_neighbors_source = num_neighbors_source, 
                                                                    num_neighbors_target = num_neighbors_target)
    
    print('\nReconstructing with alpha='+str(alpha))
    gw = gwan.gromov_wasserstein_adjusted_norm(cost_marker_genes, 
                                            cost_exp, cost_loc,
                                            alpha, 
                                            p_exp, p_loc,
                                            'square_loss', 
                                            epsilon=epsilon,
                                            tol=tol, 
                                            verbose=True)

    print('Done.')
        
    ### plot
    print('Transporting scRNA values onto locations')
    try:
        gem = adata.X.todense()
    except:
        gem = adata.X
    sgem = np.dot(gem.T, gw)

    adata_mapped = ad.AnnData(X=sgem.T)
    adata_mapped.var_names = adata.var_names
    adata_mapped.var = adata.var
    adata_mapped.obs = adata1.obs
    
    return adata_mapped

def visualize_3d(adata, col_names, ncols=3, elev=15, azim_init=-60, figsize=(15,15), save=False, folder='', alpha=0.2):
    """
    

    Parameters
    ----------
    df : dataframe
        DESCRIPTION.
    col_names : list, str
        should contain a valid column value.
    ncols : int
    elev : float
    azim : float
    figsize : tuple
    save : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    
    # check that all keys provided are valid
    assert len(col_names)>0
    assert all([i in adata.var_names for i in col_names])
    
    fig = plt.figure(figsize=figsize)
    axs = []
    nrows = (len(col_names)-1)//ncols + 1
    
    gem = adata[:,col_names].X
    
    i = 0
    for col_name in col_names:
        ax = fig.add_subplot(nrows,ncols,i+1,projection='3d')
        ax.scatter(
            adata.obs.x, adata.obs.y, adata.obs.z, 
            c = gem[:,i], 
            cmap = 'hot',
            alpha=alpha,
            linewidths=0.
           )
            
        ax.title.set_text(col_name)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.zaxis.set_ticklabels([])
        ax.set_xlim(-1,1)
        ax.set_ylim(-1,1)
        ax.set_zlim(-1,1)
        ax.view_init(elev=elev, azim=azim_init)
        axs.append(ax)
        
        i+=1
        
    plt.tight_layout()
        
    return fig