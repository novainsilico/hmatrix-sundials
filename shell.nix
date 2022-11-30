{ compiler ? "default"
, doBenchmark ? false
}:

let

  drv = (import ./default.nix { inherit doBenchmark compiler; });

in

drv.env.overrideAttrs (old: {
  nativeBuildInputs = (old.nativeBuildInputs or [ ]) ++ [ drv.haskellPackages.cabal-install ];
})
