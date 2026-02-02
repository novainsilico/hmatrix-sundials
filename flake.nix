{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  };

  outputs = { nixpkgs, ... }:
    let pkgs = nixpkgs.legacyPackages.x86_64-linux;
        sundials = pkgs.sundials.overrideAttrs (old: rec {
          version = "7.6.0";
         src = pkgs.fetchFromGitHub {
           owner = "LLNL";
           repo = "sundials";
           tag = "v${version}";
           hash = "sha256-DdVZXFfQXpJ9z5ikaK1ZQ/ZkL/vAGdlNsE9MJsIkLdM=";
         };
        });
    in {
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
