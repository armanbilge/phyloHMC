# phyloHMC

An implementation of [Probabilistic Path Hamiltonian Monte Carlo](https://arxiv.org/abs/1702.07814) for Bayesian phylogenetic inference.

# Prerequisites

phyloHMC is written in Scala. To compile and run, you will need:
* Java 8
* [sbt](http://www.scala-sbt.org/)

Optionally (but for significant speed-up):
* [libpll](https://github.com/xflouris/libpll)

To create a standalone jar, compile with the command `sbt assembly`.

# Example

An example dataset is provided.
```
sbt "run examples/surrogate.phmc examples/primates.fst 1000 1E-2 50"
```
