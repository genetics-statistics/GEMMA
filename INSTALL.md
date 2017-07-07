# INSTALL GEMMA: Genome-wide Efficient Mixed Model Association

## Check version

Simply run gemma once installed

    gemma

and it should give you the version.

## GEMMA dependencies

GEMMA runs on Linux and MAC OSX and the runtime has the following
dependencies:

* C++ tool chain
* GNU Science library (GSL)
* blas
* [Eigen library](http://eigen.tuxfamily.org/dox/)
* zlib

## Install GEMMA

### Bioconda

Recent versions of GEMMA can be installed with
[BioConda](http://ddocent.com/bioconda/) without root permissions using the following
command

    conda install gemma

### GNU Guix

The GNU Guix package manager can install recent versions of GEMMA
using the following command

    guix package -i gemma


### Install from source
