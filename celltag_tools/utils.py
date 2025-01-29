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


def call_clones_util(jac_mat_low, cells, jac_th=0.6, return_graph=True):
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
