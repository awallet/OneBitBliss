# OneBitBliss

Source code for the timing attack against the
[BLISS](https://eprint.iacr.org/2013/383) lattice-based signature scheme
described in the paper [*One Bit is All It Takes: A Devastating Timing
Attack on BLISS’s Non-Constant Time Sign
Flips*](https://eprint.iacr.org/2019/898). Based on the [original
implementation of BLISS](http://bliss.di.ens.fr/) due to Léo Ducas and
Tancrède Lepoint. We of course take no credit for their code.

Requires the [Eigen3](http://eigen.tuxfamily.org) library. Users of
Debian-like distributions can install it with:

  \# apt-get install libeigen3-dev

The attack demo also relies on the [rlutil](https://tapiov.net/rlutil/)
header-only library (included, as per the terms of the [rlutil
license](https://tapiov.net/rlutil/docs/License.txt)).
