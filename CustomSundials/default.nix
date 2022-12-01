{ stdenv
, lib
, cmake
, fetchFromGitHub
, llvmPackages
, python
, liblapack
, gfortran
, suitesparse
, lapackSupport ? true }:

let liblapackShared = liblapack.override {
  shared = true;
};

in stdenv.mkDerivation rec {
  pname = "sundials";
  version = "6.4.1.patched";

  buildInputs = lib.optionals (lapackSupport) [ gfortran ];
  nativeBuildInputs =  [ cmake ];

  src = fetchFromGitHub {
    owner = "LLNL";
    repo = pname;
    rev = "cbeb1766c23169a3306f7e82484ead644062b97c";
    sha256 = "sha256-diGykyavV365bekiW9Wdty4Et7LA+D9cx4+j+/0YaA0=";
  };

  cmakeFlags = [
    "-DEXAMPLES_INSTALL_PATH=${placeholder "out"}/share/examples"
  ] ++ lib.optionals (lapackSupport) [

    "-DSUNDIALS_INDEX_SIZE=64"

    "-DBUILD_SHARED_LIBS=OFF"
    "-DBUILD_STATIC_LIBS=ON"

    "-DBUILD_CVODE=ON"
    "-DBUILD_CVODES=OFF"
    "-DBUILD_IDA=OFF"
    "-DBUILD_IDAS=OFF"
    "-DBUILD_ARKODE=ON"
    "-DBUILD_KINSOL=OFF"
    "-DBUILD_TESTING=ON"
    "-DEXAMPLES_ENABLE_C=OFF"
    "-DEXAMPLES_ENABLE_CXX=OFF"
    "-DEXAMPLES_INSTALL=OFF"

    "-DENABLE_KLU=ON"
    "-DKLU_INCLUDE_DIR=${suitesparse.dev}/include"
    "-DKLU_LIBRARY_DIR=${suitesparse}/lib"

    "-DENABLE_LAPACK=ON"
    "-DLAPACK_LIBRARIES=${liblapackShared}/lib/liblapack${stdenv.hostPlatform.extensions.sharedLibrary};${liblapackShared}/lib/libblas${stdenv.hostPlatform.extensions.sharedLibrary}"
  ];

  doCheck = false;

  meta = with lib; {
    description = "Suite of nonlinear differential/algebraic equation solvers";
    homepage    = https://computation.llnl.gov/projects/sundials;
    platforms   = platforms.all;
    maintainers = with maintainers; [ flokli idontgetoutmuch ];
    license     = licenses.bsd3;
  };
}
