#!/bin/bash
# The same native test selection runs in a build container or a heavy Slurm tier.
set -euo pipefail
src=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$src"
if [ -n "${SLURM_JOB_ID:-}" ]; then
    root=$(dirname "$src")
    tier=$(cat "$root/tier.txt")
    . "$root/task.env"
    . /home/ert/opt/kin6d-env.sh
    export PATH=/home/ert/.local/bin:$PATH
else
    tier=${1:?usage: ci/run.sh TIER}
    root=$(mktemp -d "${RUNNER_TEMP:-/tmp}/native-ci.XXXXXX")
    export FC=/usr/bin/gfortran CC=/usr/bin/gcc CXX=/usr/bin/g++
    export OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
fi
unset DISPLAY
export CMAKE_TOOLCHAIN_FILE=$src/ci/cluster-toolchain.cmake
if [ -d "$root/mirrors" ]; then
    export GIT_CONFIG_GLOBAL=$root/gitconfig
    : > "$GIT_CONFIG_GLOBAL"
    git config --global protocol.file.allow always
    shopt -s nullglob
    for mirror in "$root"/mirrors/*/*.git; do
        slug=${mirror#"$root/mirrors/"}; slug=${slug%.git}
        git config --global --add "url.file://$mirror.insteadOf" "https://github.com/$slug.git"
        git config --global --add "url.file://$mirror.insteadOf" "https://github.com/$slug"
        git config --global --add "url.file://$mirror.insteadOf" "git@github.com:$slug.git"
    done
fi
export GITHUB_WORKSPACE=$src
export GITHUB_ENV=$root/github-env GITHUB_PATH=$root/github-path GITHUB_OUTPUT=$root/github-output
: > "$GITHUB_ENV"; : > "$GITHUB_PATH"; : > "$GITHUB_OUTPUT"
export SCI_CPUS=${SCI_CPUS:-8}
export CMAKE_BUILD_PARALLEL_LEVEL=$SCI_CPUS CTEST_PARALLEL_LEVEL=$SCI_CPUS
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export RUNNER_TEMP=$root/tmp
mkdir -p "$RUNNER_TEMP" ci-out
# Compiler probes and isolated Python builds must execute outside noexec /tmp.
export TMPDIR=$RUNNER_TEMP
{ hostname; /usr/bin/gfortran --version | head -1; git rev-parse HEAD; } | tee ci-out/toolchain.txt
load_job_environment() {
    while IFS='=' read -r key value; do
        [[ $key =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] && export "$key=$value"
    done < "$GITHUB_ENV"
    while IFS= read -r item; do [ -z "$item" ] || export PATH="$item:$PATH"; done < "$GITHUB_PATH"
}
case "$tier" in
  main-GNU)
    python3 -m venv "$root/venv"
    . "$root/venv/bin/activate"
    export LIBNEO_TESTING="1"
    echo 'Check dependency pins'
    (
python ci/sync_dependency_pins.py --check
python ci/test_dependency_pins.py
python ci/test_downstream_gate.py
python ci/test_nightly_triage.py

    )
    load_job_environment
    echo 'Install Python build dependencies'
    (
set -euo pipefail
python -m pip install --upgrade pip
python -m pip install -c constraints/ci-python.txt numpy

    )
    load_job_environment
    echo Build
    (
if [ "${SCI_MATRIX_BUILD_TYPE}" = "Release" ]; then
  make
else
  build_type="${SCI_MATRIX_BUILD_TYPE}"
  cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE="$build_type"
  cmake --build build --parallel "${SCI_CPUS:-8}"
fi

    )
    load_job_environment
    echo 'Install Python dependencies for tests'
    (
python -m pip install -c constraints/ci-python.txt ".[dev]"

    )
    load_job_environment
    echo Test
    (
# SCI_CTEST_ARGS="-LE slow" leaves the slow tests to the nightly run.
cd build && ctest --output-on-failure ${SCI_CTEST_ARGS:-}

    )
    load_job_environment
    echo 'Editable install of libneo'
    if [ "${SCI_CONDITION_9}" = true ]; then
    (
python -m pip install -e .

    )
    fi
    load_job_environment
    echo 'Run Python tests (pytest)'
    if [ "${SCI_CONDITION_10}" = true ]; then
    (
export PYTHONFAULTHANDLER="1"
export LIBNEO_DEBUG_BOOZER="1"
python -m pytest test/python -q
python -m pytest python/tests -vv -s

    )
    fi
    load_job_environment
    echo 'Build ASCOT5 and install a5py'
    if [ "${SCI_CONDITION_11}" = true ]; then
    (
git clone --depth 1 --branch 5.6.1 https://github.com/ascot4fusion/ascot5.git build/test/ascot5_compare/ascot5
make -C build/test/ascot5_compare/ascot5 libascot
python -m pip install --force-reinstall build/test/ascot5_compare/ascot5

    )
    fi
    load_job_environment
    echo 'Compare with ASCOT5 reference'
    if [ "${SCI_CONDITION_12}" = true ]; then
    (
export LD_LIBRARY_PATH="$GITHUB_WORKSPACE/build/test/ascot5_compare/ascot5/build"
python test/scripts/test_ascot5_compare.py --build-dir build --output-dir build/test/ascot5_compare

    )
    fi
    load_job_environment
    ;;
  *) echo "unknown CI tier: $tier" >&2; exit 2 ;;
esac
