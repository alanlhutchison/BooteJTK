from distutils.core import setup
from Cython.Build import cythonize
#import numpy as np


setup(
    ext_modules = cythonize('get_stat_probs.pyx'))

