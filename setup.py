#!/usr/bin/env python

# Bootstrap setuptools installation.
# Setuptools 0.6.10 is the oldest version with which we have tested Healpy,
# and is also the least common denominator present on Scientific Linux 6.
# If the user has setuptools >= 0.6.10, just take it.
# Otherwise, let use_setuptools() download its default, more recent version. 
try:
    import pkg_resources
    pkg_resources.require("setuptools >= 0.6.10")
except:
    from ez_setup import use_setuptools
    use_setuptools()

import os
from os.path import join
import errno
import sys
import shlex
from distutils.sysconfig import get_config_var
from setuptools import setup, Extension
from distutils.command.build_clib import build_clib
from distutils.errors import DistutilsExecError
from distutils.dir_util import mkpath
from distutils import log


#
# FIXME: Copied from Python 2.7's subprocess.check_output,
# but with the output= argument to CalledProcessError, which also
# dates to Python 2.7, removed.
#
# When Python 2.6 becomes unsupported, replace this with:
#   from subprocess import check_output, CalledProcessError, check_call
#
from subprocess import Popen, PIPE, CalledProcessError, check_call
def check_output(*popenargs, **kwargs):
    r"""Run command with arguments and return its output as a byte string.

    If the exit code was non-zero it raises a CalledProcessError.  The
    CalledProcessError object will have the return code in the returncode
    attribute and output in the output attribute.

    The arguments are the same as for the Popen constructor.  Example:

    >>> check_output(["ls", "-l", "/dev/null"])
    'crw-rw-rw- 1 root root 1, 3 Oct 18  2007 /dev/null\n'

    The stdout argument is not allowed as it is used internally.
    To capture standard error in the result, use stderr=STDOUT.

    >>> check_output(["/bin/sh", "-c",
    ...               "ls -l non_existent_file ; exit 0"],
    ...              stderr=STDOUT)
    'ls: non_existent_file: No such file or directory\n'
    """
    if 'stdout' in kwargs:
        raise ValueError('stdout argument not allowed, it will be overridden.')
    process = Popen(stdout=PIPE, *popenargs, **kwargs)
    output, unused_err = process.communicate()
    retcode = process.poll()
    if retcode:
        cmd = kwargs.get("args")
        if cmd is None:
            cmd = popenargs[0]
        raise CalledProcessError(retcode, cmd)
    return output
#
# FIXME: end section copied from Python 2.7's subprocess.check_output
#


# For ReadTheDocs, do not build the extensions, only install .py files
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

try:
    if sys.argv[1] in ['sdist']: #we need to distribute also c and cpp sources
        raise
    if not os.path.exists('healpy/src/_query_disc.pyx'): #pypi source package does contain the pyx files
        raise
    from Cython.Distutils import build_ext
    import Cython.Compiler.Version as CythonVersion
    from distutils.version import LooseVersion
    assert LooseVersion(CythonVersion.version) >= LooseVersion("0.16")
    ext = "pyx"
    extcpp = "pyx"
except:
  from distutils.command.build_ext import build_ext
  ext = "c"
  extcpp = "cpp"
  print "No Cython >= 0.16 found, defaulting to pregenerated c version."
  
if on_rtd:
    numpy_inc = ''
else:
    from numpy import get_include
    numpy_inc = get_include()


# Test if pkg-config is present. If not, fall back to pykg-config.
try:
    check_output(['pkg-config', '--version'])
    setup_requires = []
except OSError, e:
    if e.errno != errno.ENOENT:
        raise ValueError
    log.warn('pkg-config is not installed, falling back to pykg-config')
    setup_requires = ['pykg-config >= 1.2.0']
    os.environ['PKG_CONFIG'] = sys.executable + ' ' + os.path.abspath('pykg_config.py')


class build_external_clib(build_clib):
    """Subclass of Distutils' standard build_clib subcommand. Adds support for
    libraries that are installed externally and detected with pkg-config, with
    an optional fallback to build from a local configure-make-install style
    distribution."""

    def __init__(self, dist):
        build_clib.__init__(self, dist)
        self.build_args = {}

    @property
    def _environ(self):
        """Construct an environment dictionary suitable for having pkg-config
        pick up .pc files in the build_clib directory."""
        pkg_config_path = os.path.join(os.path.realpath(self.build_clib), 'lib', 'pkgconfig')
        try:
            pkg_config_path += ':' + os.environ['PKG_CONFIG_PATH']
        except KeyError:
            pass
        return dict(os.environ, PKG_CONFIG_PATH=pkg_config_path)

    def pkgconfig(self, *packages):
        PKG_CONFIG = tuple(shlex.split(os.environ.get('PKG_CONFIG', 'pkg-config')))
        kw = {}
        index_key_flag = (
            (2, '--cflags-only-I', ('include_dirs',)),
            (0, '--cflags-only-other', ('extra_compile_args', 'extra_link_args')),
            (2, '--libs-only-L', ('library_dirs', 'runtime_library_dirs')),
            (2, '--libs-only-l', ('libraries',)),
            (0, '--libs-only-other', ('extra_link_args',)))
        for index, flag, keys in index_key_flag:
            cmd = PKG_CONFIG + (flag,) + tuple(packages)
            log.debug('%s', ' '.join(cmd))
            args = [token[index:] for token in check_output(cmd, env=self._environ).split()]
            if args:
                for key in keys:
                    kw.setdefault(key, []).extend(args)
        return kw

    def build_library(self, library, pkg_config_name, url=None):
        log.info("checking if library '%s' is installed", library)
        try:
            build_args = self.pkgconfig(pkg_config_name)
            log.info("found '%s' installed, using it", library)
        except CalledProcessError:

            # If local_source is not specified, then immediately fail.
            if url is None:
                raise DistutilsExecError("library '%s' is not installed", library)

            from setuptools.compat import urlopen, urlparse
            from setuptools.archive_util import unpack_tarfile
            from shutil import copyfileobj
            import tarfile

            # Download tarball if necessary
            local_filename = os.path.basename(urlparse(url).path)
            if not os.path.exists(local_filename):
                log.info("fetching source for library '%s' from '%s'", library, url)
                remote_file = urlopen(url)
                local_file = open(local_filename, 'wb')
                copyfileobj(remote_file, local_file)
                remote_file.close()
                local_file.close()

            # Determine base directory from tarball table of contents
            tf = tarfile.open(local_filename)
            extract_name = tf.next().name
            tf.close()
            if extract_name.startswith('/') or '..' in extract_name.split('/'):
                raise DistutilsExecError("tarball for library '%s' seems to contain files with relative paths", library)
            build_temp = os.path.realpath(os.path.join(self.build_temp, extract_name))

            if not os.path.isdir(build_temp):
                log.info("extracting source for library '%s'", library)
                unpack_tarfile(local_filename, os.path.realpath(self.build_temp))

            log.info("building library '%s' from source", library)

            # Determine which compilers we are to use.
            cc = self.compiler.compiler[0]
            cxx = self.compiler.compiler_cxx[0]

            # Destination for headers and libraries is build_clib.
            build_clib = os.path.realpath(self.build_clib)

            # Create build directories if they do not yet exist.
            mkpath(build_clib)

            # This flag contains the architecture flags, if any.
            # On 32-bit Python on Mac OS, this will map to something like
            # '-arch i386  -DNDEBUG -g -O3  -arch i386'.
            basecflags = get_config_var('BASECFLAGS')
            cflags = get_config_var('CFLAGS')
            cflags.partition(basecflags)[-1]

            # Run configure.
            cmd = ['./configure',
                '--prefix=' + build_clib,
                'CC=' + cc,
                'CXX=' + cxx,
                'CFLAGS=' + cflags,
                'CXXFLAGS=' + cflags,
                '--disable-shared',
                '--with-pic']

            log.info('%s', ' '.join(cmd))
            check_call(cmd, cwd=build_temp, env=self._environ)

            # Run make install.
            cmd = ['make', 'install']
            log.info('%s', ' '.join(cmd))
            check_call(cmd, cwd=build_temp, env=self._environ)

            build_args = self.pkgconfig(pkg_config_name)

        return build_args
        # Done!

    def get_source_files(self):
        """Copied from Distutils' own build_clib, but modified so that it is not
        an error for a build_info dictionary to lack a 'sources' key. If there
        is no 'sources' key, then all files contained within the path given by
        the 'local_sources' value are returned."""
        self.check_library_list(self.libraries)
        filenames = []
        for (lib_name, build_info) in self.libraries:
            sources = build_info.get('sources')
            if sources is None or not isinstance(sources, (list, tuple)):
                sources = []
            filenames.extend(sources)
        return filenames

    def build_libraries(self, libraries):
        # Build libraries that have no 'sources' key, accumulating the output
        # from pkg-config.
        for lib_name, build_info in libraries:
            if 'sources' not in build_info:
                for key, value in self.build_library(lib_name, **build_info).iteritems():
                    if key in self.build_args:
                        self.build_args[key].extend(value)
                    else:
                        self.build_args[key] = value

        # Use parent method to build libraries that have a 'sources' key.
        build_clib.build_libraries(self, ((lib_name, build_info)
            for lib_name, build_info in libraries if 'sources' in build_info))


class custom_build_ext(build_ext):
    def run(self):
        # If we were asked to build any C/C++ libraries, add the directory
        # where we built them to the include path. (It's already on the library
        # path.)
        if self.distribution.has_c_libraries():
            self.run_command('build_clib')
            build_clib = self.get_finalized_command('build_clib')
            for key, value in build_clib.build_args.iteritems():
                for ext in self.extensions:
                    if not hasattr(ext, key) or getattr(ext, key) is None:
                        setattr(ext, key, value)
                    else:
                        getattr(ext, key).extend(value)
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

#start with base extension
pixel_lib = Extension('healpy._healpy_pixel_lib',
                      sources=[join('healpy','src', healpy_pixel_lib_src)],
                      include_dirs=[numpy_inc],
                      language='c++'
                      )

spht_lib = Extension('healpy._healpy_sph_transform_lib',
                     sources=[join('healpy','src', healpy_spht_src)],
                     include_dirs=[numpy_inc],
                     language='c++'
                     )

hfits_lib = Extension('healpy._healpy_fitsio_lib',
                      sources=[join('healpy','src', healpy_fitsio_src)],
                      include_dirs=[numpy_inc],
                      language='c++'
                      )

if on_rtd:
    libraries = []
    cmdclass = {}
    extension_list = []
else:
    cmdclass = {'build_ext': custom_build_ext, 'build_clib': build_external_clib}
    libraries = [
        ('cfitsio', {
        'pkg_config_name': 'cfitsio',
        'url': 'http://localhost:8000/cfitsio3340.tar.gz'}),
        # FIXME: replace with official tarball URL when U.S. government restarts
        ('healpix_cxx', {
        'pkg_config_name': 'healpix_cxx',
        'url': 'http://localhost:8000/healpix_cxx-3.11.1.tar.gz'})
        # FIXME: replace with official Autotools tarball URL
    ]
    extension_list = [pixel_lib, spht_lib, hfits_lib,
                      Extension("healpy._query_disc",
                                ['healpy/src/_query_disc.'+extcpp],
                                include_dirs=[numpy_inc],
                                language='c++'),
                      Extension("healpy._sphtools", 
                                ['healpy/src/_sphtools.'+extcpp],
                                include_dirs=[numpy_inc],
                                language='c++'),
                      Extension("healpy._pixelfunc", 
                                ['healpy/src/_pixelfunc.'+extcpp],
                                include_dirs=[numpy_inc],
                                language='c++'),
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
      package_data = {'healpy': ['data/*.fits', 'data/totcls.dat', 'test/data/*.fits', 'test/data/*.sh']},
      setup_requires=setup_requires,
      install_requires=['pyfits'],
      license='GPLv2'
      )
