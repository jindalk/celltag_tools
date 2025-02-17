import igraph as ig
import numpy as np
import pandas as pd
import scipy

from .celltag_data import CellTagData, celltag_mtx_dict


def jaccard_similarities(mat):
    """Computes the Jaccard similarity for all pairs of columns in a given
    sparse matrix.

    Args:
        mat (scipy.sparse.spmatrix):
            A binary sparse matrix (rows x columns). Each column is a feature vector
            to be compared with every other column.

    Returns:
        scipy.sparse.spmatrix:
            A sparse matrix (same shape as mat.T * mat) where each entry (i, j)
            represents the Jaccard similarity between columns i and j of the input
            matrix. The diagonal is set to 0.
    """

    cols_sum = mat.getnnz(axis=0)
    ab = mat.T * mat

    # for rows
    aa = np.repeat(cols_sum, ab.getnnz(axis=0))
    # for columns
    bb = cols_sum[ab.indices]

    similarities = ab.copy()
    similarities.data /= aa + bb - ab.data

    return similarities


def table_to_spmtx(row_data, col_data, count_data):
    """Converts row, column, and count data into a CSR (Compressed Sparse Row)
    matrix.

    Args:
        row_data (array-like):
            Row labels (e.g., cell barcodes).
        col_data (array-like):
            Column labels (e.g., CellTag identifiers).
        count_data (array-like):
            Counts or other values to populate the sparse matrix.

    Returns:
        tuple:
            A tuple (celltag_mat, cells, celltags) where:
                - celltag_mat (scipy.sparse.csr_matrix): The constructed sparse matrix
                  of shape (len(unique_rows), len(unique_columns)).
                - cells (numpy.ndarray): Sorted unique row labels.
                - celltags (numpy.ndarray): Sorted unique column labels.
    """
    cb_u = list(np.sort(np.unique(row_data)))
    celltag_u = list(np.sort(np.unique(col_data)))

    data = count_data.tolist()
    row = pd.Categorical(row_data, categories=cb_u).codes
    col = pd.Categorical(col_data, categories=celltag_u).codes
    celltag_mat = scipy.sparse.csr_matrix(
        (data, (row, col)), shape=(len(cb_u), len(celltag_u))
    )

    cells = np.array(cb_u)
    celltags = np.array(celltag_u)

    return (celltag_mat, cells, celltags)


def check_mtx_dict(target_mtx_dict):
    """Validates that the provided matrix dictionary conforms to the expected
    structure for CellTagData matrices (e.g., allow_mtx, bin_mtx, metric_mtx).

    Args:
        target_mtx_dict (celltag_mtx_dict):
            A dictionary-like object expected to contain:
                - 'mtx': A scipy.sparse.spmatrix
                - 'cells': A numpy.ndarray of cell identifiers
                - 'celltags': A numpy.ndarray of cell tag identifiers

    Raises:
        ValueError: If `target_mtx_dict` is not a celltag_mtx_dict, if it does not have
            exactly three keys ('mtx', 'cells', 'celltags'), or if the types of those
            values are incorrect.
    """
    if not isinstance(target_mtx_dict, celltag_mtx_dict):
        raise ValueError(
            "Supplied matrix dictionary attribute is not the right data type"
        )

    if not len(target_mtx_dict) == 3:
        raise ValueError("Supplied matrix dictionary is not of length 3")

    # check keys
    if set(target_mtx_dict.keys()) != {"mtx", "cells", "celltags"}:
        raise ValueError(
            "Supplied matrix dictionary should exactly have 'mtx', 'cells', and 'celltags' as keys"
        )

    # check dtypes
    if not isinstance(target_mtx_dict["mtx"], scipy.sparse.spmatrix):
        raise ValueError("cell x celltag matrix should be a scipy sparse array")

    if not isinstance(target_mtx_dict["cells"], np.ndarray):
        raise ValueError("cells should be a numpy array in matrix dictionary")

    if not isinstance(target_mtx_dict["celltags"], np.ndarray):
        raise ValueError("celltags should be a numpy array in matrix dictionary")


def find_homoplasy(n_cells, moi, barcode_abundance, ct_min=2, ct_max=25, n_iters=1000, verbose=False):
    """
    Simulates CellTag signatures in a population of cells to estimate the rate of CellTag signature duplication (homoplasy) across unrelated cells (i.e. false clones).
    
    In each iteration:
    1. A Poisson-distributed random count of CellTags is assigned to each cell (mean = `moi`).
    2. Cells with CellTag counts outside [ct_min, ct_max] are filtered out.
    3. CellTags are sampled from the provided abundance distribution and assigned to each remaining cell.
    4. The duplication rate is computed as the fraction of cell pairs sharing the exact same CellTag signature.
    
    Args:
        n_cells (int):
            The number of cells to simulate in each iteration (prior to filtering).
        moi (float):
            The mean of the Poisson distribution from which the CellTag counts per cell are drawn.
        barcode_abundance (pd.DataFrame | list):
            A DataFrame containing CellTag abundances (first column) with barcodes as the index,
            or a list of barcodes (assumed uniform abundance).
        ct_min (int, optional):
            The minimum allowed number of CellTags in a cell (inclusive). Defaults to 2.
        ct_max (int, optional):
            The maximum allowed number of CellTags in a cell (inclusive). Defaults to 25.
        n_iters (int, optional):
            The number of Monte Carlo simulation iterations to run. Defaults to 1000.
        verbose (bool, optional):
            If True, prints progress messages every 10 iterations. Defaults to False.
    
    Returns:
        list[float]:
            A list of duplication rates (homoplasy) across the simulation iterations.
            Each entry represents the duplication rate in one iteration.
    
    Raises:
        ValueError:
            If `barcode_abundance` is neither a DataFrame nor a list.
    
    Example:
        >>> # Using a uniform abundance of barcodes
        >>> homoplasy_rates = find_homoplasy(
        ...     n_cells=1000,
        ...     moi=5,
        ...     barcode_abundance=["tagA", "tagB", "tagC"],
        ...     ct_min=2,
        ...     ct_max=25,
        ...     n_iters=10,
        ...     verbose=True
        ... )
        >>> print(homoplasy_rates)
    
    Notes:
        - The duplication rate is the proportion of pairs of cells that share the exact
          same set of CellTags. It's computed as:
            net_dup_pairs / comb(len(filtered_cells), 2).
        - `comb(x, 2)` is shorthand for binomial coefficient C(x, 2) = x*(x-1)/2.
    """

    print("Simulating celltag data")
    duplication_rate_list = []

    #Handle barcode abundance table
    if(isinstance(barcode_abundance, pd.DataFrame)):
        print("Using user provided barcode abundances, please make sure CellTag barcodes are indexes and abundance is the first column")
        barcode_abundance_df = barcode_abundance.copy()
        
    elif(isinstance(barcode_abundance, list)):
        print("User provided a list of CellTag barcodes, assuming uniform abundance")
        barcode_abundance_df = pd.DataFrame(np.ones(len(barcode_abundance))/len(barcode_abundance),
                                                    index = barcode_abundance)
        
    else:
        raise ValueError("Barcode abundance should either be a DataFrame or a list")
        
    
    for iteration_curr in range(n_iters):
        if(verbose & iteration_curr%10==0):
            print(f"Iteration: {iteration_curr}")
    
        #simulate n_cells
        cells = np.random.poisson(moi, n_cells)

        #filter out cells outside ct_min and ct_max constraints
        filter_1 = cells >= ct_min
        filter_2 = cells <= ct_max

        filter_final = filter_1 & filter_2
        cells = cells[filter_final]

        #assign tags to each cell
        celltag_sigs = {}
        seen_lens = set()

        #Generating CellTag signatures
        for i in range(len(cells)):
            #simulate a celltag signature based on barcode abundance
            celltag_sig_curr = np.sort(np.random.choice(barcode_abundance_df.index,
                                                        size=cells[i],
                                                        p=barcode_abundance_df.iloc[:,0]))
            celltag_sig_curr = "/".join(celltag_sig_curr)

            #add sorted signature to dict - based on length of signature
            if(cells[i] in celltag_sigs.keys()):
                celltag_sigs[cells[i]].append(celltag_sig_curr)
            else:
                celltag_sigs[cells[i]] = []
                celltag_sigs[cells[i]].append(celltag_sig_curr)

        #check duplication rates
        net_dup_pairs = 0
        for i in celltag_sigs.keys():
            df_curr = pd.DataFrame(celltag_sigs[i])
            df_count = df_curr.value_counts()
            df_count = df_count[df_count > 1].copy()

            if(len(df_count) > 0):
                dup_pairs = sum([comb(x,2) for x in df_count.values])
                net_dup_pairs += dup_pairs
        duplication_rate_list.append(net_dup_pairs/comb(len(cells),2))
    
    print("Finished!")
    return(duplication_rate_list)


def get_clone_cell_embed(adata_obj, ct_obj, clone_weight = 1):

    """
    Creates a combined AnnData object containing both single-cell RNA data and
    clone-level “pseudo-cells” co-embedded in a knowledge graph, based on the connectivities in `adata_obj` and
    the clone assignments in `ct_obj.clone_table`.

    The new connectivity graph is constructed by:
    1. Scaling down the original `adata_obj.obsp['connectivities']` by `1 / clone_weight`
       if `clone_weight >= 1`.
    2. Building a sparse clone-cell connectivity matrix from `ct_obj.clone_table`.
    3. Combining the two connectivity matrices into a larger graph with
       rows/columns for both cells and clones.
    4. Storing the result in `adata_obj_coembed.obsp['connectivities']`.

    Args:
        adata_obj (anndata.AnnData):
            The AnnData object containing single-cell data and a precomputed neighbors
            graph in `adata_obj.obsp['connectivities']`.
        ct_obj (CellTagData):
            A valid CellTagData object containing a `clone_table` with columns for
            clone IDs and cell barcodes.
        clone_weight (float, optional):
            A scaling factor for weighting or penalizing the clone-cell connections
            relative to cell-cell connections. Defaults to 1.

    Returns:
        anndata.AnnData:
            A new AnnData object containing:
            - `.obs_names`: The concatenation of the original cell barcodes and the
              clone IDs.
            - `.obsp['connectivities']`: The merged connectivity matrix for cells and
              clones.
            - `.uns['neighbors']`: Copied parameters from the original `adata_obj`.

    Raises:
        ValueError: If `ct_obj` is invalid or missing `clone_table`, or if
            `adata_obj.obsp['connectivities']` is empty.
    """
    
    import scanpy as sc
    
    if not (isinstance(ct_obj, CellTagData)):
        raise ValueError("Please provide a valid CellTagData object")

    if not isinstance(ct_obj.clone_table, pd.DataFrame):
        raise ValueError("Clone table is either not available or not a DataFrame, please re-run tl.call_clones")

    if not 'connectivities' in adata_obj.obsp:
        raise ValueError("'connectivities' slot is empty in the AnnData obsp slot, was pp.neighbors run?")

    clone_table = CellTagData.clone_table.copy()
    all_clones = clone_table['clone_id'].unique()
    all_obs = pd.concat((pd.Series(adata_obj.obs_names), pd.Series(all_clones)))
    cell_cols = adata_obj.obs_names

    #create coembed object
    adata_obj_coembed = sc.AnnData(np.zeros((len(all_obs),100)))

    print("Imputing connectivities matrix")
    #create new connectivities based on edge list
    adata_obj_connectivities = adata_obj.obsp['connectivities']

    if(clone_weight >= 1):
        adata_obj_connectivities = adata_obj_connectivities/clone_weight

    #create clone-cell connectivities

    if(clone_weight < 1):
        imputed_connectivities, new_rows, new_cols = table_to_spmtx(clone_table['clone_id'],
                                                                             clone_table['cell.bc'],clone_weight*np.ones(len(clone_table)))
    else:
        imputed_connectivities, new_rows, new_cols = table_to_spmtx(clone_table['clone_id'],
                                                                             clone_table['cell.bc'],np.ones(len(clone_table)))

    #add missing cells (columns)
    missing_cells = np.array([*set.difference(set(cell_cols), set(new_cols))])
    new_cols = np.hstack((new_cols,missing_cells))
    imputed_connectivities = scipy.sparse.hstack((imputed_connectivities,
                                                  scipy.sparse.csr_matrix(np.zeros((len(all_clones),len(missing_cells)))))).tocsr()

    # reorder 
    idx = []
    for i in cell_cols:
        new_idx = np.where(i==new_cols)
        idx.append(new_idx[0][0])

    assert((new_cols[idx] == cell_cols).all())

    imputed_connectivities = imputed_connectivities.tocsr()[:,idx]
    new_cols = cell_cols

    conn1 = scipy.sparse.vstack((adata_obj_connectivities, imputed_connectivities))
    conn2 = scipy.sparse.vstack((imputed_connectivities.transpose(), scipy.sparse.csr_matrix(np.zeros((len(new_rows),len(new_rows))))))
    final_connectivities = scipy.sparse.hstack((conn1,conn2)).tocsr()

    #assert connectivities mtx properties
    assert((final_connectivities - final_connectivities.T).sum() == 0)
    assert((final_connectivities.diagonal()==0).all())

    #add imputed connectivities to coembed object
    adata_obj_coembed.uns['neighbors'] = dict()
    adata_obj_coembed.uns['neighbors']['params'] = adata_obj.uns['neighbors']['params']
    adata_obj_coembed.obsp['connectivities'] = final_connectivities

    all_obs = np.hstack((np.array(new_cols),new_rows))
    adata_obj_coembed.obs_names = all_obs

    
    return(adata_obj_coembed)



def merge_nn(nn_graph, all_cells, cell_list):
    """
    Merges a given list of cells with their nearest neighbors as defined by a
    nearest-neighbor graph.

    For each cell in `cell_list`, the function retrieves its neighbors from
    `nn_graph` (row corresponding to that cell in `all_cells`) and unions them
    into a set.

    Args:
        nn_graph (scipy.sparse.spmatrix or numpy.ndarray):
            A nearest-neighbor matrix where row i contains nonzero entries at
            the columns corresponding to the neighbors of cell i.
        all_cells (array-like):
            A list or array of all cell identifiers, matching the rows/columns
            of `nn_graph`.
        cell_list (array-like):
            A list of cell identifiers whose neighbors should be collected
            together.

    Returns:
        set:
            A set of cell identifiers including all `cell_list` cells plus
            any of their nearest neighbors found in `nn_graph`.
    """

    nn_cell_set = set(cell_list)
    
    for i in cell_list:
        if i in all_cells:
            cells_curr = all_cells[nn_graph[np.where(all_cells == i)[0][0],:].nonzero()[1]]
            nn_cell_set = nn_cell_set.union(cells_curr)
    return nn_cell_set


def _call_clones_util(jac_mat_low, cells, jac_th=0.6, return_graph=True):
    """Identifies clonal groups based on a Jaccard similarity matrix and
    returns a clone table, with an option to also return the underlying clone
    graph.

    Args:
        jac_mat_low (scipy.sparse.spmatrix):
            A lower-triangular Jaccard similarity matrix (e.g., from `jaccard_similarities`).
        cells (array-like):
            An array of cell identifiers corresponding to the rows/columns of the Jaccard matrix.
        jac_th (float, optional):
            Threshold above which cells are considered connected/clonally related. Defaults to 0.6.
        return_graph (bool, optional):
            If True, returns the igraph Graph object in addition to the clone table. Defaults to True.

    Returns:
        tuple | pandas.DataFrame:
            - If `return_graph=True`, returns (g, clone_table):
                - g (igraph.Graph): Graph where vertices represent cells and edges
                  represent similarity > `jac_th`.
                - clone_table (pandas.DataFrame): Table mapping cells to their clone IDs
                  and containing the edge density within each subgraph.
            - If `return_graph=False`, returns only `clone_table`.

    Raises:
        ValueError:
            If there is an issue constructing the graph or the input data is invalid.
    """

    jac_mat = jac_mat_low > jac_th

    edge_df = pd.DataFrame(np.array(jac_mat.nonzero()).T)
    vertex_df = pd.DataFrame(cells)
    vertex_df["cb"] = vertex_df[0]
    vertex_df[0] = vertex_df.index
    vertex_df.columns = ['name','cb']

    # create clone graph
    g = ig.Graph.DataFrame(edges=edge_df, directed=False, vertices=vertex_df)
    g.vs.select(_degree=0).delete()

    # call clones
    clones = g.components()

    # evaluate cols for clone table
    edge_den = [[g.induced_subgraph(i).density()] * len(i) for i in clones]
    clones_bc = [g.induced_subgraph(i).vs["cb"] for i in clones]
    clone_id = [[j + 1] * len(i) for j, i in enumerate(clones_bc)]

    # create and return clone table
    clone_table = pd.concat(
        [
            pd.DataFrame((i, j, k), index=["clone.id", "cell.bc", "edge.den"]).T
            for i, j, k in zip(clone_id, clones_bc, edge_den)
        ]
    )

    if return_graph:
        return (g, clone_table)

    return clone_table
    
def _assign_fate_util(i, fate_col = 'day', fate_key='d5', cell_type_key = 'cell_type2'):

    """
    Internal helper for fate assignment of clones. This function
    is used by `assign_fate` to process each group of cells belonging to a single clone.

    Args:
        i (pd.DataFrame):
            A subset of the clone table for a single clone.
        fate_col (str, optional):
            Column name designating the fate time point or condition. Defaults to 'day'.
        fate_key (str, optional):
            Value in `fate_col` that specifies which rows belong to the fate subset.
            Defaults to 'd5'.
        cell_type_key (str, optional):
            Column name for cell types. Defaults to 'cell_type2'.

    Returns:
        pd.DataFrame:
            The input DataFrame with added 'fate' and 'fate_pct' columns, specifying
            the clone’s fate assignment and the percentage of cells matching it.

    Notes:
        - This function references `clone_mtx` inside, which must be accessible
          in the caller's scope. The code snippet suggests usage of `clone_mtx`
          but does not show how it's defined. Adjust accordingly for your implementation.
    """
    
    clone_mtx_curr = clone_mtx[clone_mtx[fate_col] == fate_key].copy()
    if(len(clone_mtx_curr) == 0):
        clone_mtx['fate'] = 'no_fate_cells'
        clone_mtx['fate_pct'] = 0
    else:
        clone_mtx['fate'] = clone_mtx_curr[cell_type_key].value_counts(dropna=False).sort_index().idxmax()
        clone_mtx['fate_pct'] = 100*clone_mtx_curr[cell_type_key].value_counts(dropna=False).sort_index().max()/len(clone_mtx_curr)

    return(clone_mtx.copy())
    
