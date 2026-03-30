{
  description = "libneo - Fortran library for plasma physics codes";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };

        python = pkgs.python3.withPackages (ps: [
          ps.numpy
          ps.scipy
        ]);

        # Nix splits HDF5 into lib, dev, and bin outputs with different store
        # paths. The bundled hdf5-config.cmake hardcodes the lib output hash
        # in PACKAGE_PREFIX_DIR and then calls set_and_check on ${prefix}/bin
        # which does not exist in the lib-only output. Fix: merge all outputs
        # and patch the cmake config to resolve relative to the merged tree.
        hdf5 = pkgs.hdf5-fortran;
        hdf5-merged = pkgs.symlinkJoin {
          name = "hdf5-merged";
          paths = [ hdf5 hdf5.dev hdf5.bin ];
          postBuild = ''
            rm "$out/lib/cmake/hdf5-config.cmake"
            sed \
              -e 's|"''${CMAKE_CURRENT_LIST_DIR}/\.\./\.\./\.\./[^"]*"|"'"$out"'"|' \
              -e 's|''${PACKAGE_PREFIX_DIR}//nix/store/[a-z0-9]*-hdf5-cpp-fortran-[^/]*/|''${PACKAGE_PREFIX_DIR}/|g' \
              -e 's|"//nix/store/[a-z0-9]*-hdf5-cpp-fortran-[^/]*/|"'"$out"'/|g' \
              "${hdf5.dev}/lib/cmake/hdf5-config.cmake" \
              > "$out/lib/cmake/hdf5-config.cmake"

            for f in "$out"/lib/cmake/hdf5-targets*.cmake; do
              rm "$f"
              sed \
                -e 's|${hdf5}|'"$out"'|g' \
                -e 's|${hdf5.dev}|'"$out"'|g' \
                "${hdf5.dev}/lib/cmake/$(basename "$f")" \
                > "$f"
            done
          '';
        };

        buildDeps = [
          pkgs.openmpi
          pkgs.openblas
          pkgs.lapack
          pkgs.fftw
          pkgs.gsl
          hdf5-merged
          pkgs.netcdf
          pkgs.netcdffortran
          pkgs.zlib
          pkgs.curl
        ];

        nativeDeps = [
          pkgs.cmake
          pkgs.ninja
          pkgs.gfortran
          pkgs.pkg-config
        ];
      in
      {
        packages.default = pkgs.stdenv.mkDerivation {
          pname = "libneo";
          version = "0.0.0-dev";
          src = self;

          nativeBuildInputs = nativeDeps;
          buildInputs = buildDeps ++ [ python ];

          cmakeFlags = [
            "-DPREFER_SYSTEM_LIBS=ON"
            "-DLIBNEO_BUILD_TESTING=OFF"
            "-DHDF5_DIR=${hdf5-merged}/lib/cmake"
          ];
        };

        devShells.default = pkgs.mkShell {
          name = "libneo-dev";

          packages = nativeDeps ++ buildDeps ++ [
            pkgs.fftw.dev
            pkgs.gsl.dev
            pkgs.git
            python
          ];

          HDF5_DIR = "${hdf5-merged}/lib/cmake";

          shellHook = ''
            echo "libneo dev shell (all deps from nix)"
            echo "  cmake -S . -B build -G Ninja && cmake --build build"
          '';
        };
      });
}
