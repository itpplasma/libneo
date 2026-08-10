{
  description = "libneo - Fortran library for plasma physics codes";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    fortio = {
      url = "github:lazy-fortran/fortio/b6a0c57e7a7577612637ed1619c3bc50d26c2556";
      flake = false;
    };
    fortnum = {
      url = "github:lazy-fortran/fortnum/62a559c0b8e1b28bd63550164a09ed9644e659b1";
      flake = false;
    };
  };

  outputs = { self, nixpkgs, flake-utils, fortio, fortnum }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };

        python = pkgs.python3.withPackages (ps: [
          ps.numpy
          ps.scipy
        ]);

        buildDeps = [
          pkgs.openmpi
          pkgs.openblas
          pkgs.lapack
          pkgs.fftw
          pkgs.gsl
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
            "-DLIBNEO_BUILD_TESTING=OFF"
            "-DFETCHCONTENT_SOURCE_DIR_FORTIO=${fortio}"
            "-DFETCHCONTENT_SOURCE_DIR_FORTNUM=${fortnum}"
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

          shellHook = ''
            echo "libneo dev shell (all deps from nix)"
            echo "  cmake -S . -B build -G Ninja && cmake --build build"
          '';
        };
      });
}
