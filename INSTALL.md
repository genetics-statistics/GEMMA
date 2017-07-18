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

    make EIGEN_INCLUDE_PATH=~/.guix-profile/include/eigen3 FORCE_DYNAMIC=1 WITH_OPENBLAS=1

to run GEMMA tests

    make check

## Run tests

GEMMA uses the shunit2 test framework (version 2.0) and can be found
[here](https://github.com/genenetwork/shunit2)

In the source tree:

    git clone https://github.com/genenetwork/shunit2 contrib/shunit2

and run

    make check

or

    ./run_tests.sh
