project(
    'healpy',
    'cython',
    'cpp',
    version: run_command('python', '-m', 'setuptools_scm', check: true).stdout().strip(),
    meson_version: '>=1.4.0',
    default_options: [
        'cpp_std=c++11',
        'cython_language=cpp',
        'cython_args=-Xembedsignature=true',
    ],
)

py = import('python').find_installation(pure: false)

add_project_arguments(
    '-DNPY_TARGET_VERSION=NPY_2_0_API_VERSION',
    '-DNPY_NO_DEPRECATED_API=NPY_2_0_API_VERSION',
    language: 'cpp',
)

deps = [
    py.dependency(),
    dependency('numpy', version: '>=2.0.0'),
    dependency('healpix_cxx', version: '>=3.83.0'),
    dependency('openmp'),
]
inc = include_directories('src')

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
    'lib/healpy/superpixels.py',
    'lib/healpy/test_pixelfunc_extra.py',
    'lib/healpy/visufunc.py',
    'lib/healpy/zoomtool.py',
    subdir: 'healpy',
)

py.install_sources(
    'lib/healpy/utils/__init__.py',
    'lib/healpy/utils/deprecation.py',
    subdir: 'healpy/utils',
)

py.install_sources(
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
    subdir: 'healpy/data',
)

foreach name, sources : {
    '_healpy_pixel_lib': 'src/_healpy_pixel_lib.cc',
    '_healpy_sph_transform_lib': 'src/_healpy_sph_transform_lib.cc',
    '_query_disc': 'src/_query_disc.pyx',
    '_sphtools': 'src/_sphtools.pyx',
    '_pixelfunc': 'src/_pixelfunc.pyx',
    '_masktools': 'src/_masktools.pyx',
    '_hotspots': ['src/_hotspots.pyx', 'src/_healpy_hotspots_lib.cc'],
    '_line_integral_convolution': 'src/_line_integral_convolution.pyx',
}
    py.extension_module(
        name,
        sources,
        include_directories: inc,
        dependencies: deps,
        install: true,
        subdir: 'healpy',
    )
endforeach

install_data(
    'bin/healpy_get_wmap_maps.sh',
    install_dir: get_option('bindir'),
    install_mode: 'rwxr-xr-x',
)
