project('LidDrivenCavity', 'cpp',
        version : '1.0.0',
        default_options : ['warning_level=3', 'cpp_std=c++20'])

sfml_dep = dependency('sfml-all')
open_mp_dep = dependency('openmp')

# Source files
sources = files(
    'renderer.cpp',
    'main.cpp'
)

LidDrivenCavity = executable('LidDrivenCavity',
                             sources,
                             dependencies : [sfml_dep, open_mp_dep],
                             install : true)