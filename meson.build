project('LidDrivenCavity', 'cpp',
        version : '1.0.0',
        default_options : ['warning_level=3', 'cpp_std=c++20'])

sfml_dep = dependency('sfml-all')
open_mp_dep = dependency('openmp')

toml_dir = include_directories('subprojects/tomlplusplus')
toml_dep = declare_dependency(include_directories: toml_dir)

tomlplusplus = subproject('tomlplusplus')
# Source files
sources = files(
    'src/renderer.cpp',
    'src/main.cpp',
    'src/LidDrivenCavity.cpp'
)

configure_file(
    input: 'config.toml',
    output: 'config.toml',
    copy: true
)

configure_file(
    input: 'src/plot.py',
    output: 'plot.py',
    copy: true
)

LidDrivenCavity = executable('LidDrivenCavity',
                             sources,
                             dependencies : [sfml_dep, open_mp_dep, toml_dep],
                             install : true)