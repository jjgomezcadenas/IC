from distutils.core import setup
from Cython.Build import cythonize
import os
import numpy as np

# NPPATH path to numpy includes
np_path = '/'.join(np.__file__.split('/')[:-1]) + '/core/include'

setup(ext_modules=cythonize('cBLR.pyx'),
      include_dirs=[np_path])

