#!/usr/bin/env python

import platform
import os
import sys
import shutil

def is_clang_the_default_cc():
    """Check if the cc command runs clang or not. Return true if it does.
    """
    import subprocess
    import re

    try:
        cc_output = subprocess.check_output(['cc', '--version'],
                                            stderr = subprocess.STDOUT)
    except:
        return False

    return re.search('clang', cc_output) is not None

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

# For ReadTheDocs, do not build the extensions, only install .py files
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

# For each option, check is it is set or unset
# check first in defaults, then in environment and finally on command line
try:
    default_options = DEFAULT_OPT_DICT[SYSTEM_STRING]
    HEALPIX_TARGET = TARGET_DICT[SYSTEM_STRING]
except KeyError:
    raise AssertionError ('Unsupported platform: %s' % SYSTEM_STRING)

try:
    if is_clang_the_default_cc():
        print ("Detected clang compiler, disabling openMP, as it is currently unsupported")
        default_options['openmp'] = False
except:
    print ("Cannot check if compiler is clang")

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
    if sys.argv[1] in ['sdist']: #we need to distribute also c and cpp sources
        raise
    if not os.path.exists('healpy/src/_query_disc.pyx'): #pypi source package does contain the pyx files
        raise
    from Cython.Distutils import build_ext
    import Cython.Compiler.Version
    version = [int(v) for v in Cython.Compiler.Version.version.split(".")]
    assert version[1]>=14
    ext = "pyx"
    extcpp = "pyx"
except:
  from distutils.command.build_ext import build_ext
  ext = "c"
  extcpp = "cpp"
  print "No Cython >= 0.14 found, defaulting to pregenerated c version."
  
if on_rtd:
    numpy_inc = ''
else:
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

# Do we need to compile healpix_cxx ?
do_compile = (sys.argv[1] in ['build', 'build_ext', 'build_clib',
                              'bdist', 'bdist_dumb', 'bdist_rpm',
                              'bdist_wininst',
                              'install', 'install_lib', 'upload']
              and 'help' not in [x.strip(' -') for x in sys.argv]
              and not on_rtd)

if do_compile:
    compile_healpix_cxx(HEALPIX_TARGET)
    if not ( isdir(healpix_cxx_dir+'/include') and
             isdir(healpix_cxx_dir+'/lib') ):
        raise IOError("No include and lib directory : needed for healpy !")

###############################################

# Standard system libraries in /usr are included, as in most linux distribution
# the cfitsio package installed via package manager is located there;
# this way install should work often without specifying CFITSIO_EXT_PREFIX
library_dirs = ['/usr/lib', healpix_cxx_lib]
include_dirs = ['/usr/include', numpy_inc, healpix_cxx_inc]
extra_link = []

if 'CFITSIO_EXT_PREFIX' in os.environ:
    cfitsio_inc_dir = os.path.join(os.environ['CFITSIO_EXT_PREFIX'], 'include')
    cfitsio_lib_dir = os.path.join(os.environ['CFITSIO_EXT_PREFIX'], 'lib')
    include_dirs.append(cfitsio_inc_dir)
    extra_link.append(os.path.join(cfitsio_lib_dir, 'libcfitsio.a'))
if 'CFITSIO_EXT_INC' in os.environ:
    cfitsio_inc_dir = os.environ['CFITSIO_EXT_INC']
    include_dirs.append(cfitsio_inc_dir)
if 'CFITSIO_EXT_LIB' in os.environ:
    cfitsio_lib_dir = os.environ['CFITSIO_EXT_LIB']
    extra_link.append(os.path.join(cfitsio_lib_dir, 'libcfitsio.a'))

healpix_libs =['healpix_cxx','cxxsupport','psht','fftpack','c_utils']
healpix_pshyt_libs = ['psht','fftpack','c_utils']

if 'openmp' in options:
    healpix_libs += ['gomp']
    healpix_pshyt_libs += ['gomp']

if not extra_link:
    healpix_libs += ['cfitsio']

healpix_args =['-fpermissive']
if 'openmp' in options:
    healpix_args += ['-fopenmp']

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

if on_rtd:
    extension_list = []
else:
    extension_list = [pixel_lib, spht_lib, hfits_lib,
                      Extension("healpy.pshyt", ["pshyt/pshyt."+ext],
                                include_dirs = [numpy_inc,healpix_cxx_inc],
                                libraries = healpix_pshyt_libs,
                                library_dirs = library_dirs),
                      Extension("healpy._query_disc", 
                                ['healpy/src/_query_disc.'+extcpp],
                                include_dirs = [numpy_inc, healpix_cxx_inc],
                                libraries = healpix_libs,
                                library_dirs = library_dirs,
                                language='c++'),
                      Extension("healpy._sphtools", 
                                ['healpy/src/_sphtools.'+extcpp],
                                include_dirs = [numpy_inc, healpix_cxx_inc],
                                libraries = healpix_libs,
                                library_dirs = library_dirs,
                                extra_compile_args = healpix_args,
                                extra_link_args = extra_link,
                                language='c++')
                      ]
    for e in extension_list[-3:]: #extra setup for Cython extensions
        e.pyrex_directives = {"embedsignature": True}

setup(name='healpy',
      version=get_version(),
      description='Healpix tools package for Python',
      author='C. Rosset, A. Zonca',
      author_email='cyrille.rosset@apc.univ-paris-diderot.fr',
      url='http://github.com/healpy',
      packages=['healpy','healpy.test'],
      py_modules=['healpy.pixelfunc','healpy.sphtfunc',
                  'healpy.visufunc','healpy.fitsfunc',
                  'healpy.projector','healpy.rotator',
                  'healpy.projaxes','healpy.version'],
      cmdclass = {'build_ext': build_ext},
      ext_modules = extension_list,
      package_data = {'healpy': ['data/*.fits', 'data/totcls.dat']},
      license='GPLv2'
      )
