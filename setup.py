#!/usr/bin/env python

import platform
import os
import sys
import subprocess

def is_clang_or_llvm_the_default_cc():
    """Check if the cc command runs clang or not. Return true if it does.
    """
    from distutils import sysconfig
    from distutils import ccompiler
    compiler = ccompiler.new_compiler()
    sysconfig.customize_compiler(compiler)
    cc = compiler.compiler

    try:
        p = subprocess.Popen(cc + ['--version'], stdout=subprocess.PIPE)
    except OSError:
        return False

    cc_output, _ = p.communicate()

    if p.returncode:
        return False

    return ('clang' in cc_output or 'llvm' in cc_output)

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

if is_clang_or_llvm_the_default_cc():
    print ("Detected clang/llvm compiler, disabling openMP, as it is currently unsupported")
    default_options['openmp'] = False


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
from distutils.command.build_clib import build_clib
from os.path import join
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

class build_healpix(build_clib):
    def build_libraries(self, libraries):
        if 'HEALPIX_EXT_PREFIX' in os.environ:
            print "not building HEALPix; will use an external HEALPix build"
        else:
            cc = self.compiler.compiler[0]
            cxx = self.compiler.compiler_cxx[0]
            build_temp = os.path.realpath(self.build_temp)
            build_clib = os.path.realpath(self.build_clib)
            cmdline = ['make', '-w', '-C', 'hpbeta',
                'HEALPIX_TARGET=%s' % HEALPIX_TARGET,
                'HEALPIX_EXTRAFLAGS=%s' % HEALPIX_EXTRAFLAGS,
                'CC=%s' % cc,
                'CXX=%s' % cxx,
                'CL=%s' % cc,
                'CXXL=%s' % cxx,
                'BLDROOT=%s' % build_temp,
                'LIBDIR=%s' % build_clib,
                'INCDIR=%s' % build_clib,
                'BINDIR=%s' % build_clib,
            ]
            print " ".join(cmdline)
            subprocess.check_call(cmdline)

class custom_build_ext(build_ext):
    def run(self):
        # If we were asked to build any C/C++ libraries, add the directory
        # where we built them to the include path. (It's already on the library
        # path.)
        if self.distribution.has_c_libraries():
            build_clib = self.get_finalized_command('build_clib')
            self.include_dirs.append(build_clib.build_clib)
        build_ext.run(self)

def get_version():
    try:
        exec(open('healpy/version.py'))
    except Exception, e:
        print e
        raise ValueError('Error getting revision number from '
                         'healpy/version.py')
    return __version__


healpy_pixel_lib_src = '_healpy_pixel_lib.cc'
healpy_spht_src = '_healpy_sph_transform_lib.cc'
healpy_fitsio_src = '_healpy_fitsio_lib.cc'

library_dirs = []
include_dirs = [numpy_inc]
# cfitsio needs to be last in order to correctly link the healpix libraries
# so needs to go in extra_link instead of libraries
extra_link = ['-lcfitsio']
healpix_libs = []

if 'CFITSIO_EXT_PREFIX' in os.environ:
    cfitsio_inc_dir = os.path.join(os.environ['CFITSIO_EXT_PREFIX'], 'include')
    cfitsio_lib_dir = os.path.join(os.environ['CFITSIO_EXT_PREFIX'], 'lib')
    include_dirs.append(cfitsio_inc_dir)
    library_dirs.append(cfitsio_lib_dir)

# Standard system libraries in /usr are included, as in most linux distribution
# the cfitsio package installed via package manager is located there;
# this way install should work often without specifying CFITSIO_EXT_PREFIX
if not library_dirs:
    print ("""CFITSIO_EXT_PREFIX environment variable not set,
    trying with the default /usr/,
    if healpy fails at runtime, please see INSTALL""")
    include_dirs.append("/usr/include/")
    library_dirs.append("/usr/lib/")

if 'HEALPIX_EXT_PREFIX' in os.environ:
    healpix_inc_dir = os.path.join(os.environ['HEALPIX_EXT_PREFIX'], 'include')
    healpix_lib_dir = os.path.join(os.environ['HEALPIX_EXT_PREFIX'], 'lib')
    include_dirs.append(healpix_inc_dir)
    library_dirs.append(healpix_lib_dir)

if 'openmp' in options:
    extra_link.append('-lgomp')

healpix_args =['-fpermissive']
if 'openmp' in options:
    healpix_args.append('-fopenmp')

#start with base extension
pixel_lib = Extension('healpy._healpy_pixel_lib',
                      sources=[join('healpy','src', healpy_pixel_lib_src)],
                      include_dirs=include_dirs,
                      library_dirs=library_dirs,
                      libraries=healpix_libs,
                      extra_compile_args = healpix_args,
                      extra_link_args = extra_link
                      )

spht_lib = Extension('healpy._healpy_sph_transform_lib',
                     sources=[join('healpy','src', healpy_spht_src)],
                     include_dirs=include_dirs,
                     library_dirs=library_dirs,
                     libraries=healpix_libs,
                     extra_compile_args=healpix_args,
                     extra_link_args = extra_link
                     )

hfits_lib = Extension('healpy._healpy_fitsio_lib',
                      sources=[join('healpy','src', healpy_fitsio_src)],
                      include_dirs=include_dirs,
                      library_dirs=library_dirs,
                      libraries=healpix_libs,
                      extra_compile_args=healpix_args,
                      extra_link_args = extra_link
                      )

# 

if on_rtd:
    libraries = []
    cmdclass = {}
    extension_list = []
else:
    libraries = [('healpix_cxx', {'sources':[]}),
                 ('cxxsupport', {'sources':[]}),
                 ('psht', {'sources':[]}),
                 ('fftpack', {'sources':[]}),
                 ('c_utils', {'sources':[]})]
    cmdclass = {'build_ext': custom_build_ext, 'build_clib': build_healpix}
    extension_list = [pixel_lib, spht_lib, hfits_lib,
                      Extension("healpy.pshyt", ["pshyt/pshyt."+ext],
                                include_dirs = [numpy_inc] + include_dirs,
                                library_dirs = library_dirs,
                                extra_link_args = extra_link),
                      Extension("healpy._query_disc",
                                ['healpy/src/_query_disc.'+extcpp],
                                include_dirs = [numpy_inc] + include_dirs,
                                libraries = healpix_libs,
                                library_dirs = library_dirs,
                                extra_link_args = extra_link,
                                language='c++'),
                      Extension("healpy._sphtools", 
                                ['healpy/src/_sphtools.'+extcpp],
                                include_dirs = [numpy_inc] + include_dirs,
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
      libraries=libraries,
      py_modules=['healpy.pixelfunc','healpy.sphtfunc',
                  'healpy.visufunc','healpy.fitsfunc',
                  'healpy.projector','healpy.rotator',
                  'healpy.projaxes','healpy.version'],
      cmdclass = cmdclass,
      ext_modules = extension_list,
      package_data = {'healpy': ['data/*.fits', 'data/totcls.dat']},
      license='GPLv2'
      )
