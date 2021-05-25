#!/usr/bin/python

from setuptools import setup, find_packages
from glob import glob
import os

__copyright__ = "Copyright 2018-2019, ProBioPred"
__license__ = "PyPA"
__version__ = "0.0.1"
__maintainer__ = "Dattatray Mongad"

long_description = ("MicFunPred prediction tool")

files = [os.path.join(r,i).replace('micfunpreDefinitions/','') for r,d,f in os.walk('micfunpreDefinitions/data/') for i in f ]

setup(name='MicFunPred',
      version=__version__,
      description=('MicFunPred prediction'),
      maintainer=__maintainer__,
      url='https://github.com/microDM/MicFunPred',
      packages=['micfunpreDefinitions'],
      scripts=glob('scripts/*py'),
      install_requires=['biopython==1.76', 'numpy==1.18.1','pandas==1.0.5','plotly==4.9.0','pyarrow'],
      package_data={'micfunpreDefinitions':
                    files},
      #data_files=files,
      long_description=long_description)
