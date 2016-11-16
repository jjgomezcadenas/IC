from distutils.core import setup
from Cython.Build import cythonize
import os
print os.environ['NPPATH']
np_path = os.environ['NPPATH']

# NPPATH path to numpy includes
# eg=~/anaconda2/lib/python2.7/site-packages/numpy/core/include/

setup(ext_modules=cythonize('cBLR.pyx'),
      include_dirs=[np_path])

os.system('cp ./build/lib.macosx-10.6-x86_64-2.7/Cython/cBLR.so .')
