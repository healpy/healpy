#!/usr/bin/env python

import platform
import os
import sys
import shutil

TARGET_DICT = {
    'linux': 'healpy',
    'darwin': 'healpy_osx'
}

FLAGS_DICT = {
    'openmp' : '-fopenmp',
    'native' : '-march=native'
}

DEFAULT_OPT_DICT = {
    'linux': {'openmp' : True, 'native' : True},
    'darwin' : {'openmp' : True, 'native' : False}
}

SYSTEM_STRING = platform.system().lower ()

# For each option, check is it is set or unset
# check first in defaults, then in environment and finally on command line
try:
    default_options = DEFAULT_OPT_DICT[SYSTEM_STRING]
    HEALPIX_TARGET = TARGET_DICT[SYSTEM_STRING]
except KeyError:
    raise AssertionError ('Unsupported platform: %s' % SYSTEM_STRING)

# Command distclean to remove the build directory (both healpy and in hpbeta)
# 
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


options = []
for option in FLAGS_DICT:
    # Get default value
    opt_val = default_options.get(option, False)
    # Check environment variable
    env_with = 'HEALPY_WITH_' + option.upper()
    env_without = 'HEALPY_WITHOUT_' + option.upper()
    if env_with in os.environ and env_without in os.environ:
        raise ValueError('Both %s and %s environment variable are set !' %
                         (env_with, env_without))
    if env_with in os.environ:
        opt_val = True
    elif env_without in os.environ:
        opt_val = False
    # Check command line arguments
    opt_with = '--with-' + option.lower()
    opt_without = '--without-' + option.lower()
    if opt_with in sys.argv and opt_without in sys.argv:
        raise ValueError('Both %s and %s options are given on command line !' %
                         (opt_with, opt_without))
    if opt_with in sys.argv:
        opt_val = True
        while opt_with in sys.argv:
            sys.argv.remove(opt_with)
    elif opt_without in sys.argv:
        opt_val = False
        while opt_without in sys.argv:
            sys.argv.remove(opt_without)
    if opt_val:
        options.append(option)

HEALPIX_EXTRAFLAGS = ' '.join([FLAGS_DICT[opt] for opt in options])

print 'Using Healpix configuration %s for system "%s"' % (HEALPIX_TARGET,
                                                          SYSTEM_STRING)
print 'Extra flags used: "%s"' % (HEALPIX_EXTRAFLAGS)

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
    # Export necessary environment variables
    os.environ['HEALPIX_TARGET'] = target
    os.environ['HEALPIX_EXTRAFLAGS'] = HEALPIX_EXTRAFLAGS
    # launch compilation
    compil_result = os.system('cd hpbeta && make ')
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
if 'openmp' in options:
    healpix_libs.append('gomp')

if not extra_link:
    healpix_libs.append('cfitsio')

healpix_args =['-fpermissive']
if 'openmp' in options:
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
                   Extension("healpy.pshyt", ["pshyt/pshyt."+ext],
                             include_dirs = [numpy_inc,healpix_cxx_inc],
                             libraries = ['psht','gomp','fftpack','c_utils'],
                             library_dirs = library_dirs)
                   ],
      package_data={'healpy': ['data/*.fits']},
      license='GPLv2'
      )


