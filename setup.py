#!/usr/bin/env python

from setuptools import setup



setup(
    name = 'HOME',
    version = '0.5',
    description="HOME: Histogram Of MEthylation",
    author = 'akanksha srivastava',
    install_requires = [
        'numpy>=1.9.3',
        'pandas>=0.17.1',
        'scipy>=0.16.0',
        'scikit-learn>=0.16.1',
        'statsmodels>=0.6.1',
    ],
    author_email = 'akanksha.srivastava@research.uwa.edu.au',
    python_requires='>=2.7,!=3.0.*,!=3.1.*,!=3.2.*'
    scripts = ["scripts/HOME-pairwise","scripts/HOME-timeseries"],
    packages = ['HOME'],
)
