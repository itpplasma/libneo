from os import getcwd, chdir, path, popen, system
import sys
import f90nml
import ninja_syntax
# see https://gist.github.com/syoyo/778d2294fd5534e0f923 for examples
# on how to use ninja_syntax.py

def read_config(config_file):
    return f90nml.read(config_file)['config']

def build_and_run(
    target, sources, srcdir=None, builddir=None, rundir=None, config=None,
    fc = None, fflags = ''
    ):

    if fc is None:
        fc = 'gfortran'
    if rundir is None:
        rundir = getcwd()
    if builddir is None:
        builddir = path.abspath(path.join(rundir, '..', 'BUILD'))
    if srcdir is None:
        srcdir = path.abspath(path.join(rundir, '..', 'SRC'))
    if config is None:
        config = read_config(f'{target}.in')

    exe = f'{target}.x'

    sourcepaths = []
    objpaths = []
    for k in range(len(sources)):
        sourcefile = sources[k].format_map(config)
        base, ext = path.splitext(path.basename(sourcefile))
        objfile = f'{base}.o'
        sourcepaths.append(path.join(srcdir, sourcefile))
        objpaths.append(path.join(builddir, objfile))

    status_build = -1
    status_run = -1
    try:
        chdir(srcdir)
        system('git show --pretty=oneline -s')
        chdir(builddir)
        with open('build.ninja', 'w') as ninjafile:
            ninja = ninja_syntax.Writer(ninjafile)
            ninja.rule('compile', f'{fc} -o $out -c $in {fflags}')
            ninja.rule('link', f'{fc} -o $out $in {fflags}')
            for k, source in enumerate(sourcepaths):
                ninja.build(objpaths[k], 'compile', source)
            ninja.build(exe, 'link', objpaths)
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
