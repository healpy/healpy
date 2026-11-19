project('libsharp', 'c', default_options: ['c_std=c99'])
inc = include_directories('.')
meson.override_dependency(
    'libsharp',
    declare_dependency(
        link_with: shared_library(
            'sharp',
            'c_utils/c_utils.c',
            'c_utils/walltime_c.c',
            'pocketfft/pocketfft.c',
            'libsharp/sharp.c',
            'libsharp/sharp_almhelpers.c',
            'libsharp/sharp_core.c',
            'libsharp/sharp_geomhelpers.c',
            'libsharp/sharp_legendre_roots.c',
            'libsharp/sharp_ylmgen_c.c',
            dependencies: [
                meson.get_compiler('c').find_library('m'),
                dependency('openmp'),
            ],
            install: true,
        ),
        include_directories: inc,
    ),
)
