# **LevelUp**

Maintainer: Jonathan Williams

## Description

An offline viewer of ENSDF data.

The program parses the plaintext ENSDF data files into a binary database, which can then be queried for information on various nuclei (levels, cascades, gamma-rays). 

## Features

Currently implemented features include:

* Listing of levels and gamma rays corresponding to a given nucleus.
* Listing of gamma-ray cascades corresponding to a given nucleus.
* Lookup of gamma-ray cascades entered by the user, in order to determine the nucleus from which the cascade originates.

## How to Install

Use `make` to compile.  To run the program from anywhere, move the resulting `levelup` executable to any directory under your `PATH` environment variable.

The environment variable `ENSDF` pointing to a directory containing unzipped ENSDF data files must be defined in order for the program to run.  Running the program without this variable defined will display an error message.

Tested using gcc and GNU make on Ubuntu 14.04 and Arch Linux (as of March 2017).  The code is self-contained and should work on more or less any Linux distro.
