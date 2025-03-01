import os
from setuptools import setup, find_packages

install_requires=[
        "numpy>=1.21.0",
        "pandas>=1.3.0",
        "matplotlib>=3.4.0",
        "python-igraph",
        "scipy>=1.11.4"
    ]

setup(
    name="celltag_tools",
    version="0.1.0",
    description="Software for celltag clonal analysis",
    author="Kunal Jindal",
    author_email="kjkjindal@gmail.com",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "celltag_tools": ["data/*.txt"]
    },
    install_requires=install_requires,
    python_requires=">=3.6",
)