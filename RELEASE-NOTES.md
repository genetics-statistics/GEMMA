For contributions
see
[contributors](https://github.com/genetics-statistics/GEMMA/graphs/contributors)
and
[commits](https://github.com/genetics-statistics/GEMMA/commits/master).

## ChangeLog v0.98.4 (2021/01/29)

* Fix error on free with randomizer, see #239
* Moved travis-ci.org to travis-ci.com
* GEMMA builds on ARM and other architectures, see #189 and https://buildd.debian.org/status/package.php?p=gemma (thanks @tillea)
* Fixed static build with 00480e8549987b6cae7100b28bcead2a2d501177 - requires gfortran path for OpenBLAS
* Updated README's and Manual
* Added `-lmm 9` switch which shows ~beta/se~ with ~lmle~ and ~plrt~, see #237

## ChangeLog v0.98.3 (2020/11/28)

Maintenance release

* Fix Travis build with gcc 5.5 (OpenBLAS related round-offs)
* Fix Travis build on OSX (brew related)
* GEMMA installs on FreeBSD (thanks @outpaddling)
* Added github issue templates to ascertain the github issue
  tracker is only used for reporting bugs
* Added more debug output creating the GRM
* Remove info on the floating point version (gemmaf).
* Sane randomization handling: GEMMA now honours the -seed option
  (mostly for bslmm). It also allows GSL_RNG_SEED and GSL_RNG_TYPE to
  be used. See the
  [docs](https://www.gnu.org/software/gsl/doc/html/rng.html).
* The tests now use a fixed seed for the randomizer

A docker binary that runs on Linux, MaxOS and Windows can be downloaded from

http://ipfs.genenetwork.org/ipfs/Qmaq1q73ox53ykKdRF6tYDXL9bEKJQfnGCqBxFdo1fcYPb/gemma-0.98.3-AMD64-Guix-docker-release.tgz

After loading the image into Docker, run with something like

    docker run -w /run -v ${PWD}:/run ed5bf7499691 gemma -gk -bfile example/mouse_hs1940


## ChangeLog v0.98.2 (2019/05/28)

GCC 10.1 fix release

* Fix build on gcc 10.1 (mostly BLAS include files)
* Removed Eigenlib dependencies and modifed test results to match openblas output
* Also tested with guix gcc-toolchain@8.4.0 gdb gsl openblas zlib bash ld-wrapper perl vim which

## ChangeLog v0.98.1 (2018/12/10)

Bug fix release

* Fixes regression on Plink analysis with missing data #188 (thank you @voichek)

To install the image, download and

```sh
md5sum gemma-0.98.1.gz
8e7eda8091e4c4e587d91cf9d94ed147  gemma-0.98.1.gz
gzip -d gemma-0.98.1.gz
chmod a+x gemma-0.98.1
./gemma-0.98.1

GEMMA 0.98.1 (2018-12-10) by Xiang Zhou and team (C) 2012-2018

type ./gemma -h [num] for detailed help
```

The binary images were reproducibly built on x86_64 with

```sh
Generation 1    Nov 03 2018 12:30:38    (current)
  guix ff34941
    repository URL: https://git.savannah.gnu.org/git/guix.git
    commit: ff349415b27cc764fd7168ef35ca76c3b8b05889

guix environment -C guix --ad-hoc gcc gdb gfortran:lib gsl eigen openblas zlib bash ld-wrapper perl
make clean && make -j 16 && make fast-check
for x in `ldd bin/gemma|cut -d ' ' -f 3` ; do realpath $x ; done
/gnu/store/kf8pva0pgwg6yrcpa52iri293j8fc56q-gsl-2.5/lib/libgsl.so.23.1.0
/gnu/store/fxiwj2wpp11sif613axdax7gmwzsg6kp-zlib-1.2.11/lib/libz.so.1.2.11
/gnu/store/zpm494j02m6snmvcjxcdqxkgwx43nkmj-openblas-0.3.2/lib/libopenblasp-r0.3.2.so
/gnu/store/4ik1aw34nxs4372xlimvnaq2ilhclwpw-gfortran-8.2.0-lib/lib/libgfortran.so.5.0.0
/gnu/store/bmaxmigwnlbdpls20px2ipq1fll36ncd-gcc-8.2.0-lib/lib/libquadmath.so.0.0.0
/gnu/store/bmaxmigwnlbdpls20px2ipq1fll36ncd-gcc-8.2.0-lib/lib/libstdc++.so.6.0.25
/gnu/store/l4lr0f5cjd0nbsaaf8b5dmcw1a1yypr3-glibc-2.27/lib/libm-2.27.so
/gnu/store/bmaxmigwnlbdpls20px2ipq1fll36ncd-gcc-8.2.0-lib/lib/libgcc_s.so.1
/gnu/store/l4lr0f5cjd0nbsaaf8b5dmcw1a1yypr3-glibc-2.27/lib/libpthread-2.27.so
/gnu/store/l4lr0f5cjd0nbsaaf8b5dmcw1a1yypr3-glibc-2.27/lib/libc-2.27.so
/gnu/store/1yym4xrvnlsvcnbzgxy967cg6dlb19gq-gfortran-5.5.0-lib/lib/libgfortran.so.3.0.0

# build static image
make clean && make WITH_GFORTRAN=1 static -j 16 && make check
```

## ChangeLog v0.98 (2018/09/28)

With the v0.98 release GEMMA has stabilized, is faster than ever, and
contains extensive error checking. This release contains quite a few
bug fixes, hardware-based floating point checking and speedups.

* GEMMA is faster than ever #136
* Fixes log standing 'GSL ERROR: matrix is singular in lu.c at line
  266' mvlmm regression with correlated phenotypes #179 with
  https://github.com/genetics-statistics/GEMMA/commit/99527865c00b74a3a48daa2e1e5eb7c71bd861b5 -
  (thank you @HannahVMeyer and @xiangzhou)
* Adds x86 hardware based floating point error checking #161, see 70f419673d5d3e49a3eada70c70c2d284b502d7b
* Provide a static release of GEMMA again for Linux #162
* Fixes clang++ build on MacOS #160
* Enforce failed for 0.00 in src/io.cpp at line 1431 in BimbamKin #149
* -no-check is default now
* Force SIGINT on error so debuggers can get a stack trace
* Improved CalcPab functions to not break on division by zero
  https://github.com/genenetwork/GEMMA/commit/8010061e8af476d66a0ca6fb6d509b36acdb9b9a
  (thank you @xiangzhou)
* Add Windows MingW compilation support (thanks @DannyArends)
* Fully utilizing GNU Guix containers for a build system (see INSTALL.md)
* Updated manual by @xiangzhou)

Note: This is the last purely C/C++ compilable release because we are
integrating faster-lmm-d code for new functionality in the next
version. Also we are working on a Python and R interface.

To install the image, download and

```sh
md5sum gemma-0.98.gz
875cde6d37fb96014356b15dc77ebf93  gemma-0.98.gz
gzip -d gemma-0.98.gz
chmod a+x gemma-0.98
./gemma-0.98

GEMMA 0.98 (2018-09-28) by Xiang Zhou and team (C) 2012-2018

type ./gemma -h [num] for detailed help
```

The binary images were reproducibly built on x86_64 with

```sh
guix pull -l
Generation 4    Sep 25 2018 10:16:39    (current)
  guix 932839f
    repository URL: https://git.savannah.gnu.org/git/guix.git
    branch: origin/master
    commit: 932839ff124ff3b0dd3070914fb1c5beec69bf32

guix environment -C guix --ad-hoc gcc gdb gfortran:lib gsl eigen openblas zlib bash ld-wrapper perl
make clean && make -j 16 && make fast-check
for x in `ldd bin/gemma|cut -d ' ' -f 3` ; do realpath $x ; done
  /gnu/store/kf8pva0pgwg6yrcpa52iri293j8fc56q-gsl-2.5/lib/libgsl.so.23.1.0
  /gnu/store/fxiwj2wpp11sif613axdax7gmwzsg6kp-zlib-1.2.11/lib/libz.so.1.2.11
  /gnu/store/zpm494j02m6snmvcjxcdqxkgwx43nkmj-openblas-0.3.2/lib/libopenblasp-r0.3.2.so
  /gnu/store/4ik1aw34nxs4372xlimvnaq2ilhclwpw-gfortran-8.2.0-lib/lib/libgfortran.so.5.0.0
  /gnu/store/bmaxmigwnlbdpls20px2ipq1fll36ncd-gcc-8.2.0-lib/lib/libquadmath.so.0.0.0
  /gnu/store/bmaxmigwnlbdpls20px2ipq1fll36ncd-gcc-8.2.0-lib/lib/libstdc++.so.6.0.25
  /gnu/store/l4lr0f5cjd0nbsaaf8b5dmcw1a1yypr3-glibc-2.27/lib/libm-2.27.so
  /gnu/store/bmaxmigwnlbdpls20px2ipq1fll36ncd-gcc-8.2.0-lib/lib/libgcc_s.so.1
  /gnu/store/l4lr0f5cjd0nbsaaf8b5dmcw1a1yypr3-glibc-2.27/lib/libpthread-2.27.so
  /gnu/store/l4lr0f5cjd0nbsaaf8b5dmcw1a1yypr3-glibc-2.27/lib/libc-2.27.so
  /gnu/store/1yym4xrvnlsvcnbzgxy967cg6dlb19gq-gfortran-5.5.0-lib/lib/libgfortran.so.3.0.0

# build static image
make clean && make FORCE_STATIC=1 -j 16 && make check
```


## ChangeLog v0.97 (2017/12/19)

This is a massive bug fix release with many improvements.

### Speedup of GEMMA by using optimized OpenBlas

* Binary release with OpenBlas optimization for generic x86_64 and for Intel Haswell
* Dropped using standard lapack and gslcblas libs
* Fixed NaN bug with GSL2 and made recent libraries the default
* Minimized use of Eigenlib libraries (single threaded and slow compilation)
* -legacy switch provides v0.96 behaviour (incl. eigenlib)

### Added Leave One Chromosome Out (LOCO) support for Bimbam (K and LMM)

* See 449d882a3b33ef81ef4f0127c3932b01fa796dbb
* -snps [filename] option allow selecting a subset of SNPs for analysis
* -loco [chr] option for K and LMM computations
* added [gemma-wrapper](https://github.com/genetics-statistics/gemma-wrapper) to make using LOCO easy
* LOCO examples in https://github.com/genetics-statistics/GEMMA/blob/master/test/dev_test_suite.sh

### Added checks for matrices

* #72 and #45 implements
  1. Fail if K has negative eigen values
  2. Fail if K is not symmetric
  3. Fail if K is not positive definite
  4. Warn in eigen values are very small
  5. Warn if K is ill conditioned
* Check for NaN values

### Added test framework and unit tests

* Added integration and unit tests, as well as
  [Travis-CI](https://travis-ci.org/genenetwork/GEMMA) support
* Improved debug information and testing of input files

### Other

* #81 printing out beta and se(beta) under -lmm 2 as well as logl_H1
* Improved README and INSTALL docs
* Added support info and code of conduct
* Reformatted the full source tree with 3935ba39d30666dd7d4a831155631847c77b70c4
* Merged LMM computation for Plink and Bimbam formats
* Fixed progressbar issues
* #46 removed support for Oxford format
* Got rid of all compiler warnings
* Updated copyright banner, info and license information for included software
* Started a [discussion list](https://groups.google.com/forum/#!forum/gemma-discussion)

See also [commits](https://github.com/genetics-statistics/GEMMA/commits/master).

## GEMMA 0.96.0

+ First stable release.

## GEMMA 0.95.2

+ Resolved Issue #36.

## GEMMA 0.95.1

+ Created first release of GEMMA 0.95a following request in Issue #33.

## GEMMA 0.94.1

+ Fixed a bug (the predict option for multiple phenotype imputation
was not recoginzed with PLINK files).

## GEMMA 0.94.0

+ Implemented the multivariate linear mixed model.

## GEMMA 0.93

+ Implemented the Bayesian sparse linear mixed model.

## GEMMA 0.92

+ Fixed a few typos.

+ Now allows for missing values in the covariates file.

+ Included REMLE estimate for lambda in the output .log file.

+ Added small GWAS example dataset

+ Added detailed user manual.

## GEMMA 0.91

+ Fixed a bug (BIMBAM annotation file not recognized).

## GEMMA 0.90

+ Initial pre-release.

See also https://github.com/genetics-statistics/GEMMA/releases
