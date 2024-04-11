project('healpix_cxx', 'cpp')
add_project_dependencies(dependency('cfitsio'), language: 'cpp')
meson.override_dependency(
    'healpix_cxx',
    import(
        'unstable-external_project'
    ).add_project(
        'configure'
    )
)
