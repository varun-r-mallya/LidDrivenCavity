project('LidDrivenCavity', 'cpp',
        version : '1.0.0',
        default_options : ['warning_level=3', 'cpp_std=c++20'])



#LidDrivenCavity = executable('LidDrivenCavity', 'main.cpp', install : true)

#test('test', LidDrivenCavity)

# SFML dependency
sfml_dep = dependency('sfml-all')
open_mp_dep = dependency('openmp')

# Include directories
#incdir = include_directories('include')

# Source files
sources = files(
        'renderer.cpp',
        'main.cpp'
)

# Executable
LidDrivenCavity = executable('cavity_sim',
           sources,
#           include_directories : incdir,
           dependencies : [sfml_dep, open_mp_dep],
           install : true)

# Data installation (for assets like fonts)
install_data('assets/arial.ttf', install_dir : 'bin')