import pickle


class CellTagData:
    """A container class for storing and managing CellTag data, including read
    data, thresholds, various matrices, and clone information. Provides methods
    for serialization and easy attribute access.

    Attributes:
        ct_reads (pd.DataFrame | None):
            The raw or processed CellTag read data.
        thresholds (dict | None):
            A dictionary of thresholds used at various processing steps
            (e.g., 'starcode', 'triplet', 'binarization', etc.).
        seq_sat (float | None):
            Sequencing saturation value.
        clone_graph (igraph.Graph | None):
            Graph representing clone relationships among cells.
        allow_mtx (celltag_mtx_dict):
            Dictionary-like container for the allow matrix and its axes.
        bin_mtx (celltag_mtx_dict):
            Dictionary-like container for the binarized matrix and its axes.
        metric_mtx (celltag_mtx_dict):
            Dictionary-like container for the metric-filtered matrix and its axes.
        jaccard_mtx (scipy.sparse.spmatrix | None):
            Jaccard similarity matrix.
        clone_table (pd.DataFrame | None):
            Table mapping cells to their clones.
        clone_info (pd.DataFrame | None):
            Table containing clone level metadata.
            Defaults to None.
    """

    def __init__(
        self,
        ct_reads=None,
        thresholds=None,
        seq_sat=None,
        clone_graph=None,
        jaccard_mtx=None,
        clone_table=None,
        clone_info = None
    ):
        """Initializes the CellTagData object with optional parameters for
        immediate data assignment.

        Args:
            ct_reads (pd.DataFrame | None, optional):
                Raw or processed CellTag read data. Defaults to None.
            thresholds (dict | None, optional):
                Dictionary of thresholds used at various processing steps. Defaults to None.
            seq_sat (float | None, optional):
                Sequencing saturation value. Defaults to None.
            clone_graph (igraph.Graph | None, optional):
                Graph representing clone relationships. Defaults to None.
            jaccard_mtx (scipy.sparse.spmatrix | None, optional):
                Jaccard similarity matrix. Defaults to None.
            clone_table (pd.DataFrame | None, optional):
                Table mapping cells to their clones.
                Defaults to None.
            clone_info (pd.DataFrame | None, optional):
                Table containing clone level metadata.
                Defaults to None.
        """

        self.ct_reads = ct_reads
        self.thresholds = thresholds
        self.seq_sat = seq_sat
        self.clone_graph = clone_graph

        # matrices
        self.allow_mtx = celltag_mtx_dict(
            {"mtx": None, "cells": None, "celltags": None}
        )
        self.bin_mtx = celltag_mtx_dict({"mtx": None, "cells": None, "celltags": None})
        self.metric_mtx = celltag_mtx_dict(
            {"mtx": None, "cells": None, "celltags": None}
        )
        self.jaccard_mtx = jaccard_mtx
        self.clone_table = clone_table
        self.clone_info = clone_info

    def save(self, path):
        """Saves the current CellTagData object to a file using pickle
        serialization.

        Args:
            path (str):
                The file path where the serialized object will be stored.
        """

        with open(path, "wb") as f:
            pickle.dump(self, f)

    def __repr__(self):
        return f"{self.__class__.__name__}"


class celltag_mtx_dict:
    """A specialized dictionary-like container for a sparse matrix ('mtx') and
    its associated row ('cells') and column ('celltags') labels, specifically
    tailored for the CellTagData workflow.

    This class strictly maintains three keys:
        - 'mtx': A scipy.sparse matrix representing the cell-tag matrix.
        - 'cells': A numpy.ndarray of cell identifiers (e.g., barcodes).
        - 'celltags': A numpy.ndarray of tag identifiers (e.g., CellTags).

    Any attempt to add keys beyond these three raises a KeyError. Standard dictionary
    operations (e.g., indexing, iteration) are supported, but the structure is
    ensured to remain consistent with these fixed keys.

    Key Features:
        • Only three keys ('mtx', 'cells', 'celltags') are permitted.
        • You can retrieve or set each key via subscript notation (e.g., `obj['mtx']`).
        • The `is_empty()` method returns True if all three keys are set to None.
        • The length (`len(obj)`) is always 3.
        • The `__repr__` method provides a concise, formatted view of the data.

    Example Usage:
        ctdict = celltag_mtx_dict({'mtx': None, 'cells': None, 'celltags': None})
        if ctdict.is_empty():
            print("All entries are currently None.")
        ctdict['mtx'] = my_sparse_matrix
        ctdict['cells'] = my_cell_array
        ctdict['celltags'] = my_celltag_array
    """


    def __init__(self, initial_data):
        """Initializes the celltag_mtx_dict with a predefined structure.

        Args:
            initial_data (dict):
                A dictionary containing exactly the keys 'mtx', 'cells', and 'celltags'.
                The values should be compatible with the CellTagData workflow
                (e.g., sparse matrices, numpy arrays).
        """

        self._data = dict(initial_data)

    def __getitem__(self, key):
        return self._data[key]

    def __setitem__(self, key, value):
        if key not in self._data:
            raise KeyError(
                f"Cannot add new key: '{key}'. Only existing keys can be modified."
            )
        self._data[key] = value

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return iter(self._data)

    def keys(self):
        return self._data.keys()

    def values(self):
        return self._data.values()

    def items(self):
        return self._data.items()

    def __repr__(self):
        return f"{self.__class__.__name__}({self._data})"

    def is_empty(self):
        """Check if all keys are set to None."""
        return all(value is None for value in self._data.values())
