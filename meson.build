project(
    'healpy',
    'cython', 'cpp',
    meson_version: '1.4.0',
    default_options: [
        'cpp_std=c++11',
        'cython_language=cpp',
        'cython_args=-Xembedsignature=true'])

py_mod = import('python')
py = py_mod.find_installation(pure: false, modules:'numpy')

add_project_arguments(
    '-DNPY_NO_DEPRECATED_API=NPY_1_19_API_VERSION',
    language: 'cpp')

add_project_dependencies(
    dependency('numpy'),
    dependency('healpix_cxx'),
    declare_dependency(include_directories: 'src'),
    language: 'cpp')

py.install_sources(
    'lib/healpy/__init__.py',
    'lib/healpy/conftest.py',
    'lib/healpy/cookbook.py',
    'lib/healpy/fitsfunc.py',
    'lib/healpy/newvisufunc.py',
    'lib/healpy/pixelfunc.py',
    'lib/healpy/projaxes.py',
    'lib/healpy/projector.py',
    'lib/healpy/rotator.py',
    'lib/healpy/sphtfunc.py',
    'lib/healpy/visufunc.py',
    'lib/healpy/zoomtool.py',
    subdir: 'healpy'
)

py.install_sources(
    'lib/healpy/utils/__init__.py',
    'lib/healpy/utils/deprecation.py',
    subdir: 'healpy/utils'
)

py.install_sources(
    'lib/healpy/data/pixel_window_n0002.fits',
    'lib/healpy/data/pixel_window_n0004.fits',
    'lib/healpy/data/pixel_window_n0008.fits',
    'lib/healpy/data/pixel_window_n0016.fits',
    'lib/healpy/data/pixel_window_n0032.fits',
    'lib/healpy/data/pixel_window_n0064.fits',
    'lib/healpy/data/pixel_window_n0128.fits',
    'lib/healpy/data/pixel_window_n0256.fits',
    'lib/healpy/data/pixel_window_n0512.fits',
    'lib/healpy/data/pixel_window_n1024.fits',
    'lib/healpy/data/pixel_window_n2048.fits',
    'lib/healpy/data/pixel_window_n4096.fits',
    'lib/healpy/data/pixel_window_n8192.fits',
    'lib/healpy/data/planck_cmap.dat',
    'lib/healpy/data/planck_log_cmap.dat',
    'lib/healpy/data/totcls.dat',
    'lib/healpy/data/weight_ring_n00002.fits',
    'lib/healpy/data/weight_ring_n00004.fits',
    'lib/healpy/data/weight_ring_n00008.fits',
    'lib/healpy/data/weight_ring_n00016.fits',
    'lib/healpy/data/weight_ring_n00032.fits',
    'lib/healpy/data/weight_ring_n00064.fits',
    'lib/healpy/data/weight_ring_n00128.fits',
    'lib/healpy/data/weight_ring_n00256.fits',
    'lib/healpy/data/weight_ring_n00512.fits',
    'lib/healpy/data/weight_ring_n01024.fits',
    'lib/healpy/data/weight_ring_n02048.fits',
    'lib/healpy/data/weight_ring_n04096.fits',
    'lib/healpy/data/weight_ring_n08192.fits',
    'lib/healpy/data/wmap_cmap.dat',
    subdir: 'healpy/data'
)

py.extension_module(
    'healpy._healpy_pixel_lib', 
    'src/_healpy_pixel_lib.cc')
py.extension_module(
    'healpy._healpy_sph_transform_lib', 
    'src/_healpy_sph_transform_lib.cc')
py.extension_module(
    'healpy._query_disc',
    'src/_query_disc.pyx')
py.extension_module(
    'healpy._sphtools',
    'src/_sphtools.pyx')
py.extension_module(
    'healpy._pixelfunc',
    'src/_pixelfunc.pyx')
py.extension_module(
    'healpy._masktools',
    'src/_masktools.pyx')
py.extension_module(
    'healpy._hotspots',
    ['src/_hotspots.pyx', 'src/_healpy_hotspots_lib.cc'])
py.extension_module(
    'healpy._line_integral_convolution',
    'src/_line_integral_convolution.pyx')
