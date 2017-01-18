#!/usr/bin/env python

from setuptools import setup



setup(
    name = 'HOME',
    version='0.0',
    description="HOME: Histogram Of MEthylation",
    author = 'akanksha srivastava',
    install_requires = [
        'numpy==1.10.1',
        'pandas==0.17.1',
        'rpy2==2.7.7',
        'scikit-learn==0.16.1',
        'scipy==0.16.0',
        'statsmodels==0.6.1',
    ],
    author_email = 'akanksha.iitkgp.2012@gmail.com',
    
    scripts = ["scripts/HOME-pairwise","scripts/HOME-timeseries"],
    packages = ['HOME'],
)
