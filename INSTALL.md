# INSTALL GEMMA: Genome-wide Efficient Mixed Model Association

## Check version

Simply run gemma once installed

    gemma

and it should give you the version.

## GEMMA dependencies

GEMMA runs on Linux, MAC OSX and Windows (with Docker). The runtime
has the following dependencies:

* C++ tool chain >= 5.5.0 (see Travis CI and we test with file .guix-dev-gcc-older)
* GNU Science library (GSL) 2.x (GEMMA dropped support for GSL 1.x)
* blas/openblas
* lapack
* zlib

See below for installation on Guix.

## Install GEMMA

### Debian and Ubuntu

Travis-CI uses Ubuntu for testing. Check the test logs for version numbers.

[![Build Status](https://travis-ci.org/genetics-statistics/GEMMA.svg?branch=master)](https://travis-ci.org/genetics-statistics/GEMMA)

Current settings can be found in [travis.yml](.travis.yml).

### Bioconda

(Note Bioconda install is a work in [progress](https://github.com/genetics-statistics/GEMMA/issues/52)

Recent versions of GEMMA can be installed with
[BioConda](http://ddocent.com/bioconda/) without root permissions using the following
command

    conda install gemma

### FreeBSD

Recent editions of FreeBSD ports include [GEMMA](https://www.freebsd.org/cgi/ports.cgi?query=gemma&stype=all)

### GNU Guix

The GNU Guix package manager can install recent versions of [GEMMA](https://www.gnu.org/software/guix/packages/g.html)
using the following command

    guix package -i gemma

A more recent version may be found in the guix-bioinformatics channel
which is maintained by the authors. See the
[README](http://git.genenetwork.org/guix-bioinformatics/guix-bioinformatics), e.g.

    env GUIX_PACKAGE_PATH=./guix-bioinformatics guix package -A gemma

To build GEMMA from source you can opt to install the build tools with
GNU Guix, the current build container is in [guix-dev](./.guix-dev)

    source .guix-dev
    make

Guix allows for easy versioning. To build with an older gcc, for
example:

    guix environment -C guix --ad-hoc gcc-toolchain@9.3.0 gdb gsl openblas zlib bash ld-wrapper perl vim which

### Install with Docker

Recent version of GEMMA come with a 64-bit Docker image that should run
on Linux, Windows and MacOS.

### Install from source

Install listed dependencies (you may want to take hints from
the Travis-CI [tests](./.travis.yml)) and run

	make -j 4

(the -j switch builds on 4 cores).

	time make check

You can run gemma in the debugger with, for example

	gdb --args \
		./bin/gemma -g example/mouse_hs1940.geno.txt.gz \
		-p example/mouse_hs1940.pheno.txt -a example/mouse_hs1940.anno.txt \
		-snps example/snps.txt -nind 400 -loco 1 -gk -debug -o myoutput

Note that if you get <optimized out> warnings on inspecting variables you
should compile with GCC_FLAGS="" to disable optimizations (-O3). E.g.

    make WITH_OPENBLAS=1 GCC_FLAGS=

Other options, such as compiling with warnings, are listed in the
Makefile.

### GNU Guix commands used

Some development examples.  With git bisect build the older versions
of gemma with openblas

    ~/.config/guix/current/bin/guix environment -C guix --ad-hoc gcc gdb gfortran:lib gsl lapack openblas zlib bash ld-wrapper perl ldc
    make clean ; make WITH_OPENBLAS=1 FORCE_DYNAMIC=1 -j 8

or with atlas

    ~/.config/guix/current/bin/guix environment -C guix --ad-hoc gcc gdb gfortran:lib gsl lapack atlas zlib bash ld-wrapper perl ldc
    make clean ; make WITH_OPENBLAS= FORCE_DYNAMIC=1 -j 25

## Run tests

GEMMA includes the shunit2 test framework (version 2.0).

    make check

or

    ./run_tests.sh

## Releases

### Docker release

To distribute GEMMA I made static versions of the binary. A container
can be made instead with, for example

```sh
env GUIX_PACKAGE_PATH=~/guix-bioinformatics ~/.config/guix/current/bin/guix \
  pack -f docker gemma-gn2 -S /bin=bin
```

which created a container in of size 51MB. Tiny! For more information
see
[GUIX-NOTES](http://git.genenetwork.org/guix-bioinformatics/guix-notes/CONTAINERS.org).


### Static release

To create a static release, locate the gfortran lib and use

    source .guix-dev-static
    make WITH_GFORTRAN=1 EXTRA_FLAGS=-L/gnu/store/741057r2x06zwg6zcmqmdyv51spm6n9i-gfortran-7.5.0-lib/lib static

otherwise OpenBlas will complain with

    undefined reference to `_gfortran_concat_string'
