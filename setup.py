from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize('get_stat_probs1.pyx')
    )
