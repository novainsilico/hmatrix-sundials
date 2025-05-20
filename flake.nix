{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  };

  outputs = { nixpkgs, ... }:
    let pkgs = nixpkgs.legacyPackages.x86_64-linux;
    in {
      packages.x86_64-linux.default = pkgs.callPackage ./default.nix { };
      devShells.x86_64-linux.default = pkgs.mkShell {
        buildInputs = with pkgs; [
          cabal-install
          ghc
          haskellPackages.haskell-language-server
          sundials
          blas
          lapack
          suitesparse
        ];
      };
    };
}
