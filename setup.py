#!/usr/bin/env python

from setuptools import setup



setup(
    name = 'HOME',
    version = '0.4',
    description="HOME: Histogram Of MEthylation",
    author = 'akanksha srivastava',
    install_requires = [
        'numpy',
        'pandas==0.17.1',
        'scipy==0.16.0',
        'scikit-learn==0.16.1',
        'statsmodels==0.6.1',
    ],
    author_email = 'akanksha.srivastava@research.uwa.edu.au',
    
    scripts = ["scripts/HOME-pairwise","scripts/HOME-timeseries"],
    packages = ['HOME'],
)
