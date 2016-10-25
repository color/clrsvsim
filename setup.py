import os

try:
  from setuptools import setup
except:
  from distutils.core import setup


setup(name = "clr-sv-sim",
      version = "0.0.1",
      description = "Color Genomics Structural Variant Simulator",
      author = "Color Genomics",
      author_email = "dev@color.com",
      url = "https://github.com/ColorGenomics/clr-sv-sim",
      packages = ["clr-sv-sim"],
      install_requires=[
        'cigar==0.1.3',
        'mock==2.0.0',
        'nose==1.3.7',
        'numpy==1.10.1',
        'preconditions==0.1',
        'pyfasta==0.5.2',
        'pysam==0.9.0',
      ],
      license = "Apache-2.0",
      )
