#!/usr/bin/env python

import platform
import os
import sys
import shutil

TARGET_DICT = {
    'linux': 'healpy',
    'darwin': 'healpy_osx'
}

# Option to use OpenMP in healpix routines
opt_without_openmp = '--without-openmp'
without_openmp = opt_without_openmp in sys.argv
if without_openmp:
    sys.argv.remove(opt_without_openmp)
with_openmp = not without_openmp

# Option to use compiler flag "-march=native"
opt_without_native = '--without-native'
without_native = opt_without_native in sys.argv
if without_native:
    sys.argv.remove(opt_without_native)
with_native = not without_native

SYSTEM_STRING = platform.system().lower ()
try:
    HEALPIX_TARGET=TARGET_DICT[SYSTEM_STRING]
    HEALPIX_EXTRAFLAGS=""
    if with_openmp:
        HEALPIX_EXTRAFLAGS += "-fopenmp "
    if with_native:
        HEALPIX_EXTRAFLAGS += "-march=native "
    print 'Using Healpix configuration "%s" for system "%s"' % \
            (HEALPIX_TARGET, SYSTEM_STRING)
    print 'Extra flags used: "%s"' % (HEALPIX_EXTRAFLAGS)
except KeyError:
    raise AssertionError ('Unsupported platform: %s' % SYSTEM_STRING)

if 'distclean' in sys.argv:
    # Remove build directory of healpy and hpbeta
    build_hpx = os.path.join('hpbeta', 'build.' + HEALPIX_TARGET)
    print 'Removing ', build_hpx, ' directory...'
    shutil.rmtree(build_hpx, True)
    hpx = os.path.join('hpbeta', HEALPIX_TARGET)
    print 'Removing ', hpx, 'directory...'
    shutil.rmtree(hpx, True)
    hpy = 'build'
    print 'Removing ', hpy, ' directory...'
    shutil.rmtree(hpy, True)
    sys.exit(0)

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
  print "No Cython >= 0.12 found, defaulting to pregenerated c version."
  
from numpy import get_include
numpy_inc = get_include()

def compile_healpix_cxx(target):
    import os
    print "Compiling healpix_cxx (this may take a while)"
    compil_result = os.system('cd hpbeta && '
                              'HEALPIX_TARGET=%s HEALPIX_EXTRAFLAGS="%s" make '%(target,HEALPIX_EXTRAFLAGS) )
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

library_dirs = [healpix_cxx_lib]
include_dirs = [numpy_inc, healpix_cxx_inc]
extra_link = []

if 'CFITSIO_EXT_PREFIX' in os.environ:
    cfitsio_inc_dir = os.path.join(os.environ['CFITSIO_EXT_PREFIX'], 'include')
    cfitsio_lib_dir = os.path.join(os.environ['CFITSIO_EXT_PREFIX'], 'lib')
    include_dirs.append(cfitsio_inc_dir)
    #library_dirs.append(cfitsio_lib_dir)
    extra_link.append(os.path.join(cfitsio_lib_dir, 'libcfitsio.a'))
if 'CFITSIO_EXT_INC' in os.environ:
    cfitsio_inc_dir = os.environ['CFITSIO_EXT_INC']
    include_dirs.append(cfitsio_inc_dir)
if 'CFITSIO_EXT_LIB' in os.environ:
    cfitsio_lib_dir = os.environ['CFITSIO_EXT_LIB']
    #library_dirs.append(cfitsio_lib_dir)
    extra_link.append(os.path.join(cfitsio_lib_dir, 'libcfitsio.a'))

healpix_libs =['healpix_cxx','cxxsupport','psht','fftpack','c_utils']
if with_openmp:
    healpix_libs.append('gomp')

if not extra_link:
    healpix_libs.append('cfitsio')

healpix_args =['-fpermissive']
if with_openmp:
    healpix_args.append('-fopenmp')

#start with base extension
pixel_lib = Extension('healpy._healpy_pixel_lib',
                      sources=[join('healpy','src',s)
                               for s in healpy_pixel_lib_src],
                      include_dirs=include_dirs,
                      library_dirs=library_dirs,
                      libraries=healpix_libs,
                      extra_compile_args = healpix_args,
                      extra_link_args = extra_link
                      )

spht_lib = Extension('healpy._healpy_sph_transform_lib',
                     sources=[join('healpy','src',s) for s in healpy_spht_src],
                     include_dirs=include_dirs,
                     library_dirs=library_dirs,
                     libraries=healpix_libs,
                     extra_compile_args=healpix_args,
                     extra_link_args = extra_link
                     )

hfits_lib = Extension('healpy._healpy_fitsio_lib',
                      sources=[join('healpy','src',s)
                               for s in healpy_fitsio_src],
                      include_dirs=include_dirs,
                      library_dirs=library_dirs,
                      libraries=healpix_libs,
                      extra_compile_args=healpix_args,
                      extra_link_args = extra_link
                      )

# 
setup(name='healpy',
      version=get_version(),
      description='Healpix tools package for Python',
      author='C. Rosset',
      author_email='cyrille.rosset@apc.univ-paris-diderot.fr',
      url='http://www.healpy.org',
      packages=['healpy'],
      py_modules=['healpy.pixelfunc','healpy.sphtfunc',
                  'healpy.visufunc','healpy.fitsfunc',
                  'healpy.projector','healpy.rotator',
                  'healpy.projaxes','healpy.version'],
      cmdclass = {'build_ext': build_ext},
      ext_modules=[pixel_lib,spht_lib,hfits_lib,
                 #  Extension("healpy.pshyt", ["pshyt/pshyt."+ext],
                 #            include_dirs = [numpy_inc,healpix_cxx_inc],
                 #            libraries = ['psht','gomp','fftpack','c_utils'],
                 #            library_dirs = library_dirs)
                   ],
      package_data={'healpy': ['data/*.fits']},
      license='GPLv2'
      )


