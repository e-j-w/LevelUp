# **LevelUp**

Maintainer: Jonathan Williams

## Description

An offline viewer of ENSDF data, with extensions for analyzing gamma-ray spectra.

The program parses plaintext ENSDF data files (available at: https://www-nds.iaea.org/ensdf_base_files/) into a binary database, which can then be queried for information on various nuclei (levels, cascades, gamma-rays).

## Features

Right now the program's capabilities are fairly limited, but more will be added as I find use for them in my research.  Currently implemented features include:

* Listing of levels and gamma rays corresponding to a given nucleus.
* Listing of gamma-ray cascades corresponding to a given nucleus.
* Checking for overlapping/nearby gamma-ray energies in two nuclei, or in one nucleus and nuclei in the region of another nucleus.
* Lookup of gamma-ray cascades entered by the user, in order to determine the nucleus from which the cascade originates.
* Determination of nuclei present in gamma ray spectra based on cascades present in the spectrum.

## How to Install

### Dependencies

* gcc
* make
* readline

The current version has been tested under Centos 7.

### Instructions

Use `make` to compile.  To run the program from anywhere, move the resulting `levelup` executable to any directory under your `PATH` environment variable.

The environment variable `ENSDF` pointing to a directory containing unzipped ENSDF data files must be defined in order for the program to run.  Running the program without this variable defined will display an error message.

## Notes

### Input Data Types

Some functions in this program allow for analysis of gamma-ray spectra.  The program recognizes the following file formats for gamma-ray spectra:

**.mca** - An .mca file is simply a 2D array of integers, with the first index denoting a spectrum number (up to 100) and the second index denoting a bin number (up to 32768).

**.spe** - An .spe file is the data type written by radware when using the 'ws' command in gf3.


