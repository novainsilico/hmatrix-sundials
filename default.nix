{ suitesparse, sundials, haskellPackages, ... }:
  haskellPackages.callCabal2nix "hmatrix-sundials" ./. {
    klu = suitesparse;
    suitesparseconfig = suitesparse;
    sundials_arkode = sundials;
    sundials_cvode = sundials;
    sundials_core = sundials;
    sundials_sunlinsolklu = sundials;
    sundials_sunmatrixsparse = sundials;
  }
