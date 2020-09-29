# INSTALL GEMMA: Genome-wide Efficient Mixed Model Association

## Check version

Simply run gemma once installed

    gemma

and it should give you the version.

## GEMMA dependencies

GEMMA runs on Linux and MAC OSX and the runtime has the following
dependencies:

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

To build GEMMA from source you can opt to install the build tools with GNU Guix

    guix package -i make gcc linux-libre-headers gsl openblas lapack glibc ld-wrapper

The current build container is in [guix-dev](../guix-dev)

    guix environment -C guix --ad-hoc gcc-toolchain gdb gsl openblas zlib bash ld-wrapper perl vim which

To build with an older gcc, for example:

    guix environment -C guix --ad-hoc gcc-toolchain@9.3.0 gdb gsl openblas zlib bash ld-wrapper perl vim which

### Install from source

Note: Eigen is no longer required!

Install listed dependencies and run

	make -j 4

(the -j switch builds on 4 cores).

if you get an Eigen error you may need to override the include
path. E.g. to build GEMMA on GNU Guix with shared libs the following
may work

    make EIGEN_INCLUDE_PATH=~/.guix-profile/include/eigen3 WITH_OPENBLAS=1

another example overriding optimization and LIB flags (so as to link
against gslv1) would be

    make EIGEN_INCLUDE_PATH=~/.guix-profile/include/eigen3 WITH_OPENBLAS=1 GCC_FLAGS="-Wall" LIBS="$HOME/opt/gsl1/lib/libgsl.a $HOME/opt/gsl1/lib/libgslcblas.a -L$HOME/.guix-profile/lib -pthread -llapack -lblas -lz"

to run GEMMA tests

	time make check

You can run gemma in the debugger with, for example

	gdb --args \
		./bin/gemma -g example/mouse_hs1940.geno.txt.gz \
		-p example/mouse_hs1940.pheno.txt -a example/mouse_hs1940.anno.txt \
		-snps example/snps.txt -nind 400 -loco 1 -gk -debug -o myoutput

Note that if you get <optimized out> warnings on inspecting variables you
should compile with GCC_FLAGS="" to disable optimizations (-O3). E.g.

    make EIGEN_INCLUDE_PATH=~/.guix-profile/include/eigen3 WITH_OPENBLAS=1 GCC_FLAGS=

Other options, such as compiling with warnings, are listed in the
Makefile.

### GNU Guix commands used

With git bisect build the older versions of gemma with openblas

    ~/.config/guix/current/bin/guix environment -C guix --ad-hoc gcc gdb gfortran:lib gsl eigen lapack openblas zlib bash ld-wrapper perl ldc
    make clean ; make EIGEN_INCLUDE_PATH=$GUIX_ENVIRONMENT/include/eigen3 WITH_OPENBLAS=1 FORCE_DYNAMIC=1 -j 8

or with atlas

    ~/.config/guix/current/bin/guix environment -C guix --ad-hoc gcc gdb gfortran:lib gsl eigen lapack atlas zlib bash ld-wrapper perl ldc
    make clean ; make EIGEN_INCLUDE_PATH=$GUIX_ENVIRONMENT/include/eigen3/Eigen/ WITH_OPENBLAS= FORCE_DYNAMIC=1 -j 25

You may need to symlink Eigen in some older versions

    ln -s $GUIX_ENVIRONMENT/include/eigen3/Eigen src/Eigen


## Run tests

GEMMA includes the shunit2 test framework (version 2.0).

    make check

or

    ./run_tests.sh
