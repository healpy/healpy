project('cfitsio', 'c')
meson.override_dependency(
    'cfitsio',
    import(
        'unstable-external_project'
    ).add_project(
        'configure'
    ).dependency(
        'cfitsio'
    )
)
