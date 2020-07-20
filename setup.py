#!/usr/bin/python

from setuptools import setup, find_packages
from glob import glob
import os

__copyright__ = "Copyright 2018-2019, ProBioPred"
__license__ = "PyPA"
__version__ = "0.0.1"
__maintainer__ = "Dattatray Mongad"

long_description = ("MicFunPred prediction tool")

setup(name='MicFunPred',
      version=__version__,
      description=('MicFunPred prediction'),
      maintainer=__maintainer__,
      url='https://github.com/microDM/MicFunPred',
      packages=['micfunpreDefinitions'],
      scripts=glob('scripts/*py'),
      install_requires=['biopython', 'numpy'],
      package_data={'micfunpreDefinitions':
                    ['data/db/*',
                     'data/other/*',
                     'data/phenotype/*',
                     'data/*',
                     'data/db/blastDB/*']},
      long_description=long_description)
