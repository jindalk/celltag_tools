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
