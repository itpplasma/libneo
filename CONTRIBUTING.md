# Contributing to libneo

Coding conventions are described in [CODING_STANDARD.md](CODING_STANDARD.md),
and architectural notes in [DESIGN.md](DESIGN.md).

## Dependency pins

Fortio and Fortnum revisions are defined in `fpm.toml`; both the main CMake
build and the standalone `h5merge` tool read them. To test another revision,
configure with `-DFORTIO_REF=<ref>` or `-DFORTNUM_REF=<ref>`; local source
overrides remain available through CMake's FetchContent options.

After changing a revision in `fpm.toml`, run
`python3 ci/sync_dependency_pins.py` to regenerate the Nix inputs and lockfile.
CI checks their consistency with `python3 ci/sync_dependency_pins.py --check`.

## Merging

- Pull requests, including every layer of a stacked PR, run the fast CI tier.
  Label a PR `full-ci` to run the full tier on it before queueing.
- Review a stack layer by layer. When it is ready, collapse it into one PR into
  `main` (retarget the top layer to `main`); the lower PRs then show as merged.
- Click "Merge when ready" (or `gh pr merge --merge --auto`). The PR enters the
  merge queue, the full CI tier runs on the exact merge result, and on success
  it lands on `main` as a merge commit. `main` itself is not re-tested on push.
- A branch need not be up to date with `main`. Merge `main` into a PR only
  when it conflicts.
