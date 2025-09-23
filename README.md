These are Haskell bindings to the [sundials][1] ODE solver suite
(specifically, the CVODE, ARKODE, IDA solvers) used at
[Novadinsilico][2] in the jinko platform [3].

This repository is a fork of the original `hmatrix-sundials` and includes a few
changes:

- Support for IDA solver (algebraic rules)
- A complete rewrite of the main loop from `inline-c` to pure haskell. The main
  benefit is a better develloper experience (e.g. hacking pure haskell with
  haskell-language-server is "funnier" than editing `inline-c` block) as well
  as a better support for asynchronous exceptions. In the future we also hope
  to unlock better API, such as a `Stream` based timepoint producer.
- Extended (but not complete) sundial bindings.
- Better logging
- A callback usable by the caller in order to get notified of solving progress
- Support for editing with Haskell-language-server

[1]: https://computing.llnl.gov/projects/sundials
[2]: https://www.novainsilico.ai
[3]: https://www.jinko.ai
