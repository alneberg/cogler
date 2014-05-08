#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os

version = '0.1.0'

setup(name='cogler',
      version=version,
      description="Summarise and visualise COG annotations for metagenomic contigs",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='Python Scilifelab Metagenomics Contig',
      author='Johannes Alneberg',
      author_email='johannes.alneberg@scilifelab.se',
      url='https://github.com/alneberg/cogler',
      license='Creative Commons',
      packages=find_packages(exclude=['examples', 'tests']),
      scripts=["scripts/COG_table.py", 
          "scripts/COGPlot.R"],
      include_package_data=True,
      zip_safe=False,
      install_requires= ['biopython>=1.62b',
                        'bcbio-gff>=0.4',
                        'nose==1.3.0'],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )

