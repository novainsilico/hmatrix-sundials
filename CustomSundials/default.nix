{ stdenv
, lib
, cmake
, fetchFromGitHub
, llvmPackages
, python
, liblapack
, gfortran
, suitesparse
  # Whether to compile sundials with lapack support.
, lapackSupport ? true
  # The index type used by sundials (probably "int64_t" or "int32_t").
  # See https://sundials.readthedocs.io/en/latest/Install_link.html#cmakeoption-SUNDIALS_INDEX_TYPE
, indexType ? null
  # The index size used by sundials ("64" or "32").
  # See https://sundials.readthedocs.io/en/latest/Install_link.html#cmakeoption-SUNDIALS_INDEX_SIZE
, indexSize ? null
}:

let
  liblapackShared = liblapack.override {
    shared = true;
  };

in
stdenv.mkDerivation rec {
  pname = "sundials";
  version = "6.4.1.patched";

  buildInputs = lib.optionals (lapackSupport) [ gfortran ];

  nativeBuildInputs = [ cmake ];

  src = fetchFromGitHub {
    owner = "LLNL";
    repo = pname;
    # develop branch of LLNL/sundials corresponding to v6.4.1 + https://github.com/LLNL/sundials/pull/230
    rev = "cbeb1766c23169a3306f7e82484ead644062b97c";
    sha256 = "sha256-diGykyavV365bekiW9Wdty4Et7LA+D9cx4+j+/0YaA0=";
  };

  cmakeFlags =
    [
      "-DBUILD_SHARED_LIBS=ON"
      "-DBUILD_STATIC_LIBS=OFF"

      "-DBUILD_ARKODE=ON"
      "-DBUILD_CVODE=ON"
      "-DBUILD_CVODES=OFF"
      "-DBUILD_IDA=OFF"
      "-DBUILD_IDAS=OFF"
      "-DBUILD_KINSOL=OFF"

      "-DENABLE_KLU=ON"
      "-DKLU_INCLUDE_DIR=${suitesparse.dev}/include"
      "-DKLU_LIBRARY_DIR=${suitesparse}/lib"

      "-DEXAMPLES_ENABLE_C=OFF"
      "-DEXAMPLES_ENABLE_CXX=OFF"
      "-DEXAMPLES_INSTALL=OFF"
    ] ++
    lib.optionals (indexType != null) [ "-DSUNDIALS_INDEX_TYPE=${indexType}" ] ++
    lib.optionals (indexSize != null) [ "-DSUNDIALS_INDEX_SIZE=${indexSize}" ] ++
    lib.optionals (lapackSupport) [
      "-DENABLE_LAPACK=ON"
      "-DLAPACK_LIBRARIES=${liblapackShared}/lib/liblapack${stdenv.hostPlatform.extensions.sharedLibrary};${liblapackShared}/lib/libblas${stdenv.hostPlatform.extensions.sharedLibrary}"
    ];

  doCheck = false;

  meta = with lib; {
    description = "Suite of nonlinear differential/algebraic equation solvers";
    homepage = https://computation.llnl.gov/projects/sundials;
    platforms = platforms.all;
    maintainers = with maintainers; [ flokli idontgetoutmuch ];
    license = licenses.bsd3;
  };
}
