from os import getcwd, chdir, path, popen, system
import sys
import f90nml
import ninja_syntax
# see https://gist.github.com/syoyo/778d2294fd5534e0f923 for examples
# on how to use ninja_syntax.py

def build_and_run(
    target, sources, srcdir=None, builddir=None, rundir=None, config_file=None
    ):
    fc = 'gfortran'  # TODO: make configurable in dict or explicitly

    if rundir is None:
        rundir = getcwd()
    if builddir is None:
        builddir = path.abspath(path.join(rundir, '..', 'BUILD'))
    if srcdir is None:
        srcdir = path.abspath(path.join(rundir, '..', 'SRC'))
    if config_file is None:
        config_file = f'{target}.in'

    config = f90nml.read(config_file)['config']
    exe = f'{target}.x'

    sourcepaths = []
    for k in range(len(sources)):
        sourcepaths.append(path.join(srcdir, sources[k].format_map(config)))

    status_build = -1
    status_run = -1
    try:
        chdir(builddir)
        with open('build.ninja', 'w') as ninjafile:
            ninja = ninja_syntax.Writer(ninjafile)
            ninja.rule('fc', f'{fc} $in -o $out')
            ninja.build(exe, 'fc', sourcepaths)
        status_build = system('ninja')

    finally:
        chdir(rundir)

    if status_build == 0:  # build has succeeded
        try:
            status_run = system(path.join(builddir, exe))
        finally:
            if status_run != 0:
                print(f'Error: {exe} exited with code {status_run}.')
    else:
        print(f'Error: Exiting due to build error for {exe}.')
