{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs";
  };

  outputs = {nixpkgs, ...}:
    let pkgs = nixpkgs.legacyPackages.x86_64-linux;
    in {
    devShells.x86_64-linux.default = pkgs.mkShell {
      buildInputs = with pkgs; [ cabal-install ghc haskell-language-server sundials
      blas lapack
      suitesparse

    ];
    };
  };
}
