The ODEINT golden test fetches the latest public libneo tag by default. Its
generated Fortran source is written into the build directory.

For offline or captured builds, provide the already adapted public reference
as `golden_odeint_dynamic.f90` in a separate directory and configure with
`-DLIBNEO_ODEINT_GOLDEN_SOURCE_DIR=/path/to/reference` and
`-DLIBNEO_ODEINT_GOLDEN_SHA256=<expected SHA256>`. The hash is checked at
configuration and before copying the input into the build directory. The
reference file is also declared as a native CMake source input for capture.

For example, the adapted reference from public tag v2026.10.16 (commit
c92b33667b299062473855f7ebe0fe9a59b627ec) has SHA256
5df050fa0fe6581b92c4ba6b6c70619e2a9bdb1a3f737c9204abdc6f26b18d8b.
The reference implementation and numerical comparison tolerances are unchanged.
