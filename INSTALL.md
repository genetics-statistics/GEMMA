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

(Note Bioconda install is a work in [progress](https://github.com/xiangzhou/GEMMA/issues/52)

Recent versions of GEMMA can be installed with
[BioConda](http://ddocent.com/bioconda/) without root permissions using the following
command

    conda install gemma

### GNU Guix

The GNU Guix package manager can install recent versions of [GEMMA](https://www.gnu.org/software/guix/packages/g.html)
using the following command

    guix package -i gemma

### Install from source

Install listed dependencies and run

    make

if you get an Eigen error you may need to override the include
path. E.g. on GNU Guix with shared libs this may work

    make EIGEN_INCLUDE_PATH=~/.guix-profile/include/eigen3 FORCE_DYNAMIC=1
