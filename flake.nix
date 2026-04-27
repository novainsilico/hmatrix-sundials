{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";

    treefmt-nix.url = "github:numtide/treefmt-nix";
  };

  outputs = { nixpkgs, treefmt-nix, ... }:
    let
      pkgs = nixpkgs.legacyPackages.x86_64-linux;
      sundials = pkgs.sundials.overrideAttrs (old: rec {
        version = "7.6.0";
        src = pkgs.fetchFromGitHub {
          owner = "LLNL";
          repo = "sundials";
          tag = "v${version}";
          hash = "sha256-DdVZXFfQXpJ9z5ikaK1ZQ/ZkL/vAGdlNsE9MJsIkLdM=";
        };
      });

      treefmt = treefmt-nix.lib.evalModule pkgs {
        projectRootFile = "flake.nix";
        programs.nixpkgs-fmt.enable = true;
        programs.ormolu = {
          enable = true;

          # Setup the default extensions used
          # In editor, this setting is done through HLS
          # ghcOpts = import ./buildlib/default_extensions.nix;

        };
      };
    in
    {
      formatter.x86_64-linux = treefmt.config.build.wrapper;

      packages.x86_64-linux.default = pkgs.callPackage ./default.nix { inherit sundials; };
      devShells.x86_64-linux.default = pkgs.mkShell {
        buildInputs = [
          pkgs.cabal-install
          pkgs.ghc
          pkgs.haskellPackages.haskell-language-server
          sundials
          pkgs.blas
          pkgs.lapack
          pkgs.suitesparse
        ];
      };
    };
}
