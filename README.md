# **LevelUp** - Find what is in your gamma ray spectrum!

Maintainer: Jonathan Williams


## Goals

* Parse ENSDF files and build up a list of gamma decay chains.
* Read in a list of gamma-ray energies present in a spectrum
* See which of the known decay chains are present in that list.
* Report which species are present in the spectrum.

## How to Install

Use 'make' to compile.  To run the program from anywhere, move the resulting 'levelup' executable to any directory under your $PATH environment variable.

Tested using gcc and GNU make on Ubuntu 14.04.  The code is self-contained and should work on more or less any Linux distro.
