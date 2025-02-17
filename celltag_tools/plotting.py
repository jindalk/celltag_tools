# plotting
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .celltag_data import CellTagData, celltag_mtx_dict
from .utils import check_mtx_dict


def diagnostic_plots(ct_obj, mtx_use="metric"):
    """Generates diagnostic scatter and histogram plots for a specified matrix
    (allow, bin, or metric) within a CellTagData object.

    Args:
        ct_obj (CellTagData):
            The CellTagData object containing the matrix to be visualized.
        mtx_use (str, optional):
            The type of matrix to plot. Must be one of {"allow", "bin", "metric"}.
            Defaults to "metric".

    Raises:
        ValueError:
            If `ct_obj` is not a valid CellTagData object, or if `mtx_use` is not in
            {"metric", "allow", "bin"}.

    Notes:
        The function checks whether the specified matrix dictionary (`mtx_use_mtx`) is valid,
        then creates a 2x2 figure showing:
        - A scatter plot of row sums (CellTags per cell).
        - A scatter plot of column sums (cells per CellTag).
        - A histogram of row sums.
        - A histogram of column sums.
    """
    if not (isinstance(ct_obj, CellTagData)):
        raise ValueError("Please provide a valid CelltagData object")

    if mtx_use not in ["metric", "allow", "bin"]:
        raise ValueError("Please define mtx_use as either 'allow', 'bin', or 'metric'")

    matrix_type = f"{mtx_use}_mtx"
    mtx_dict = getattr(ct_obj, matrix_type)

    check_mtx_dict(mtx_dict)

    r2 = mtx_dict["mtx"]
    fig, ax = plt.subplots(2, 2, figsize=(10, 10))
    ax[0, 0].scatter(
        np.expand_dims(np.arange(0, r2.shape[0]), axis=1),
        np.array(r2.sum(axis=1)),
        rasterized=True,
    )
    ax[0, 1].scatter(
        np.expand_dims(np.arange(0, r2.shape[1]), axis=1),
        np.array(r2.sum(axis=0)),
        rasterized=True,
    )
    ax[1, 0].hist(np.array(r2.sum(axis=1)), bins=100, rasterized=True)
    # ax[1,0].set_ylim((0,50))
    ax[1, 1].hist(np.array(r2.sum(axis=0).transpose()), bins=100, rasterized=True)

    ax[0, 0].set_title("celltags/cell")
    ax[0, 1].set_title("cells/celltag")

    fig.suptitle("metric filtered mtx QC", fontsize=19)
    fig.tight_layout()



def plot_size_by_den(ct_obj, highlight_clones=None, ax=None, **kwargs):
    """
    Creates a scatter plot of clone size versus edge density for visualization,
    with the option to highlight specific clones by coloring or labeling them.

    Args:
        clone_meta (pd.DataFrame):
            A DataFrame containing per-clone metadata, expected to have columns:
            - 'edge.den': The clone's edge density.
            - 'size': The clone's size (number of cells).
            - 'clone.id': The clone ID for labeling (if highlighting).
        highlight_clones (list | np.ndarray | None, optional):
            A list of clone IDs to highlight and label in red. Defaults to None.
        ax (matplotlib.axes.Axes | None, optional):
            A matplotlib Axes object on which to draw the plot. If None, the current
            Axes (`plt.gca()`) is used. Defaults to None.
        **kwargs:
            Additional keyword arguments passed to `sns.scatterplot`.

    Returns:
        matplotlib.axes.Axes:
            The Axes object containing the scatter plot.

    Raises:
        ValueError: If `clone_info` is missing or invalid, or if `highlight_clones`
            is not a valid list/array when provided.
    """

    import adjustText

    if not (isinstance(ct_obj, CellTagData)):
        raise ValueError("Please provide a valid CellTagData object")

    if not (isinstance(ct_obj.clone_info, pd.DataFrame)):
        raise ValueError("Clone table is either not available or not a DataFrame, please re-run tl.call_clones")

    clone_meta = ct_obj.clone_info.copy()
    
    if(ax==None):
        ax=plt.gca()
    ax = sns.scatterplot(x = "edge.den", y="size", data=clone_meta,color='royalblue', **kwargs)

    #highlight clones if any
    if(highlight_clones is not None):

        if not isinstance(highlight_clones, (np.ndarray, list)):
            raise TypeError("highlight_clones should be either a list or a numpy array")
        
        if not highlight_clones:
            print("highlight_clones is set to an emplty list")
            return ax
        
            
        clone_ids_to_label = np.array(highlight_clones)
        clone_meta_sub = clone_meta[clone_meta['clone.id'].isin(clone_ids_to_label)].copy()
        texts = [plt.text(x['edge.den'], x['size'], "clone "+str(int(x['clone.id'])),
                          ha='center', va='center',
                          fontsize=12) for _,x in clone_meta_sub.iterrows()]
        adjustText.adjust_text(texts, expand_points=(1.5, 1.5), arrowprops=dict(arrowstyle="-", color='k', lw=0.5))

    return(ax)
