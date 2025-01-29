.. celltag_tools documentation master file, created by
   sphinx-quickstart on Tue Jan 28 15:11:51 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CellTag-tools: Computational framework for Clonal analysis
===========================

.. toctree::
   :maxdepth: 1
   :caption: Tutorials
   :hidden:

   simple_tutorial
   allowlist

.. toctree::
   :maxdepth: 1
   :caption: API Reference
   :hidden:

   celltag_data
   tl
   pl
   utils



CellTag-tools is a highly scalable and efficient software package for analysing multi-modal single cell lineage tracing data. Read more about CellTagging at `Jindal et al. Nat. Biotech. (2023) <https://www.nature.com/articles/s41587-023-01931-4>`_ , `Kong et al. Nat. Protocols (2020) <https://www.nature.com/articles/s41596-019-0247-2>`_, and `Biddy et al. Nature (2018) <https://www.nature.com/articles/s41586-018-0744-4>`_


CellTag-tools identifies clonal relationships within single-cell datasets, starting from CellTag reads enabling insights into cellular lineage and dynamics.

Key Features
------------

- **CellTag Read Extraction**  
  Extracts CellTag sequences from single-cell BAM files.
  
- **Error Correction and Filtering**  
  Refines CellTag data for high accuracy downstream analysis.
  
- **Clone Identification**  
  Identifies clonal populations based on shared CellTag signatures.

Workflow Overview
-----------------

1. **Parsing Single-Cell BAM Files**  
   Extract CellTag sequences from BAM files, supporting pipelines like CellRanger and CellRanger-ATAC.

   * Requirements:
     - Sample configuration file (specifies sample details, paths, cell barcodes, assay types, CellTag versions).
     - A pipeline processes the samples, as outlined on our `GitHub repository <https://github.com/morris-lab/newCloneCalling>`_.

2. **Processing CellTag Reads to Identify Clones**
   After parsing, the tool performs:

   - **Filtering and Error Correction**  
     Filters out low-quality reads and corrects sequencing errors.

   - **Allowlisting**  
     Validates CellTag sequences using predefined lists.

   - **Matrix Construction**  
     Constructs cell-by-CellTag matrices, which are filtered and binarized.

   - **Clone Calling**  
     Identifies clones based on shared CellTag patterns using similarity matrices.

Installation
------------

To install the package, clone the repository and use ``pip``:

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/jindalk/celltag_tools.git

   # Navigate to the repository
   cd celltag_tools

   # Install the package
   pip install .

