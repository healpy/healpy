#!/usr/bin/env python

import platform

TARGET_DICT = {
    'linux': 'healpy',
    'darwin': 'healpy_osx'
}

SYSTEM_STRING = platform.system().lower ()
try:
    HEALPIX_TARGET=TARGET_DICT[SYSTEM_STRING]
    print 'Using Healpix configuration "%s" for system "%s"' % \
            (HEALPIX_TARGET, SYSTEM_STRING)
except KeyError:
    raise AssertionError ('Unsupported platform: %s' % SYSTEM_STRING)

from distutils.core import setup, Extension
from os.path import join,isdir
import sys

try:
  from Cython.Distutils import build_ext
  import Cython.Compiler.Version
  version = [int(v) for v in Cython.Compiler.Version.version.split(".")]
  assert version[1]>=12
  ext = "pyx"
except Exception:
  from distutils.command.build_ext import build_ext
  ext = "c"
  print "No Cython>0.12 found, defaulting to pregenerated c version."
  
from numpy import get_include
numpy_inc = get_include()

def compile_healpix_cxx(target):
    import os
    print "Compiling healpix_cxx (this may take a while)"
    compil_result = os.system('cd hpbeta && '
                              'HEALPIX_TARGET=%s make '%(target))
    if compil_result != 0:
        raise Exception('Error while compiling healpix_cxx')

def get_version():
    try:
        exec(open('healpy/version.py'))
    except Exception, e:
        print e
        raise ValueError('Error getting revision number from '
                         'healpy/version.py')
    return __version__


healpy_pixel_lib_src = ['_healpy_pixel_lib.cc']
healpy_spht_src = ['_healpy_sph_transform_lib.cc']
healpy_fitsio_src = ['_healpy_fitsio_lib.cc']


################################################
#
#    Healpix data (pixel window and ring files
healpix_cxx_dir='hpbeta/%s'%HEALPIX_TARGET
healpix_cxx_inc = healpix_cxx_dir+'/include'
healpix_cxx_lib = healpix_cxx_dir+'/lib'

if sys.argv[1] != 'sdist':
    compile_healpix_cxx(HEALPIX_TARGET)
    if not ( isdir(healpix_cxx_dir+'/include') and
             isdir(healpix_cxx_dir+'/lib') ):
        raise IOError("No include and lib directory : needed for healpy !")

###############################################

healpix_libs =['healpix_cxx','cxxsupport','psht','fftpack','c_utils','cfitsio','gomp']
healpix_args =['-fopenmp']

#start with base extension
pixel_lib = Extension('healpy._healpy_pixel_lib',
                      sources=[join('healpy','src',s)
                               for s in healpy_pixel_lib_src],
                      include_dirs=[numpy_inc,healpix_cxx_inc],
                      library_dirs=[healpix_cxx_lib],
                      libraries=healpix_libs,
                      extra_compile_args=healpix_args
                      )

spht_lib = Extension('healpy._healpy_sph_transform_lib',
                     sources=[join('healpy','src',s) for s in healpy_spht_src],
                     include_dirs=[numpy_inc,healpix_cxx_inc],
                     library_dirs=[healpix_cxx_lib],
                     libraries=healpix_libs,
                     extra_compile_args=healpix_args
                     )

hfits_lib = Extension('healpy._healpy_fitsio_lib',
                      sources=[join('healpy','src',s)
                               for s in healpy_fitsio_src],
                      include_dirs=[numpy_inc,healpix_cxx_inc],
                      library_dirs=[healpix_cxx_lib],
                      libraries=healpix_libs,
                      extra_compile_args=healpix_args
                      )

# 
setup(name='healpy',
      version=get_version(),
      description='Healpix tools package for Python',
      author='C. Rosset',
      author_email='rosset@lal.in2p3.fr',
      url='http://code.google.com/p/healpy',
      packages=['healpy'],
      py_modules=['healpy.pixelfunc','healpy.sphtfunc',
                  'healpy.visufunc','healpy.fitsfunc',
                  'healpy.projector','healpy.rotator',
                  'healpy.projaxes','healpy.version'],
      cmdclass = {'build_ext': build_ext},
      ext_modules=[pixel_lib,spht_lib,hfits_lib,
                   Extension("pshyt", ["pshyt/pshyt."+ext],
                             include_dirs = [numpy_inc,healpix_cxx_inc],
                             libraries = ['psht','gomp','fftpack','c_utils'],
                             library_dirs = [healpix_cxx_lib])
                   ],
      package_data={'healpy': ['data/*.fits']},
      license='GPLv2'
      )


