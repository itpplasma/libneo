{
  description = "libneo - Fortran library for plasma physics codes";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    fortio = {
      url = "github:lazy-fortran/fortio/6929247e03cc8657cc3f6aabf210525f16237fd7";
      flake = false;
    };
    fortnum = {
      url = "github:lazy-fortran/fortnum/30c14fe2b4674655a3cac9b2525c60aa7bc31738";
      flake = false;
    };
    fortfem = {
      url = "github:lazy-fortran/fortfem/e871b6c5ab14158e5a482cf93c8c17e855621a52";
      flake = false;
    };
    # FortFEM's own CMake pin of its sparse solver dependency.
    fortsparse = {
      url = "github:lazy-fortran/fortsparse/d2d912b95139871b9c8b2fdd695574496bc98810";
      flake = false;
    };
  };

  outputs = { self, nixpkgs, flake-utils, fortio, fortnum, fortfem, fortsparse }:
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
            "-DFETCHCONTENT_SOURCE_DIR_FORTFEM=${fortfem}"
            "-DFETCHCONTENT_SOURCE_DIR_FORTSPARSE=${fortsparse}"
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
