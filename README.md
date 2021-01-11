![Genetic associations identified in CFW mice using GEMMA (Parker et al,
Nat. Genet., 2016)](cfw.gif)

# GEMMA: Genome-wide Efficient Mixed Model Association

[![Build Status](https://travis-ci.org/genetics-statistics/GEMMA.svg?branch=master)](https://travis-ci.com/genetics-statistics/GEMMA) [![Anaconda-Server Badge](https://anaconda.org/bioconda/gemma/badges/installer/conda.svg)](https://anaconda.org/bioconda/gemma) [![DL](https://anaconda.org/bioconda/gemma/badges/downloads.svg)](https://anaconda.org/bioconda/gemma) [![BrewBadge](https://img.shields.io/badge/%F0%9F%8D%BAbrew-gemma--0.98-brightgreen.svg)](https://github.com/brewsci/homebrew-bio) [![GuixBadge](https://img.shields.io/badge/gnuguix-gemma-brightgreen.svg)](https://www.gnu.org/software/guix/packages/G/) [![DebianBadge](https://badges.debian.net/badges/debian/testing/gemma/version.svg)](https://packages.debian.org/search?keywords=gemma&searchon=names&suite=all&section=all)

GEMMA is a software toolkit for fast application of linear mixed
models (LMMs) and related models to genome-wide association studies
(GWAS) and other large-scale data sets.

Check out [RELEASE-NOTES.md](RELEASE-NOTES.md) to see what's new in
each GEMMA release.

Please post suspected bugs to
[Github issues](https://github.com/genetics-statistics/GEMMA/issues). For
questions or other discussion, please post to the
[GEMMA Google Group](https://groups.google.com/group/gemma-discussion). We
also encourage contributions, for example, by forking the repository,
making your changes to the code, and issuing a pull request.

Currently, GEMMA provides a runnable Docker container for 64-bit
MacOS, Windows and Linux platforms. GEMMA can be installed with
Debian, Conda, Homebrew and GNU Guix. With Guix you find the latest
version
[here](http://git.genenetwork.org/guix-bioinformatics/guix-bioinformatics)
as it is the version we use every day on http://genenetwork.org. For
installation instructions see also [INSTALL.md](INSTALL.md).  We use
continous integration builds on Travis-CI for Linux (amd64 & arm64)
and MacOS (amd64). GEMMA builds on multiple architectures, see the
[Debian build farm](https://buildd.debian.org/status/package.php?p=gemma).

*(The above image depicts physiological and behavioral trait
loci identified in CFW mice using GEMMA, from [Parker et al, Nature
Genetics, 2016](https://doi.org/10.1038/ng.3609).)

* [Key features](#key-features)
* [Installation](#installation)
  * [Precompiled binaries](#precompiled-binaries)
* [Run GEMMA](#run-gemma)
  * [Debugging and optimization](#debugging-and-optimization)
* [Help](#help)
* [Citing GEMMA](#citing-gemma)
* [License](#license)
* [Optimizing performance](#optimizing-performance)
* [Building from source](#building-from-source)
* [Input data formats](#input-data-formats)
* [Reporting a GEMMA bug or issue](#reporting-a-gemma-bug-or-issue)
  * [Check list:](#check-list)
* [Code of conduct](#code-of-conduct)
* [Credits](#credits)

## Key features

1. Fast assocation tests implemented using the univariate linear mixed
model (LMM). In GWAS, this can correct for population structure and
sample non-exchangeability. It also provides estimates of the
proportion of variance in phenotypes explained by available genotypes
(PVE), often called "chip heritability" or "SNP heritability".

2. Fast association tests for multiple phenotypes implemented using a
multivariate linear mixed model (mvLMM). In GWAS, this can correct for
population structure and sample (non)exchangeability - jointly in
multiple complex phenotypes.

3. Bayesian sparse linear mixed model (BSLMM) for estimating PVE,
phenotype prediction, and multi-marker modeling in GWAS.

4. Estimation of variance components ("chip/SNP heritability") partitioned
by different SNP functional categories from raw (individual-level)
data or summary data. For raw data, HE regression or the REML AI
algorithm can be used to estimate variance components when
individual-level data are available. For summary data, GEMMA uses the
MQS algorithm to estimate variance components.

## Installation

To install GEMMA you can

1. Download the precompiled or Docker binaries
   from [releases](https://github.com/genetics-statistics/GEMMA/releases).

2. Use existing package managers, see [INSTALL.md](INSTALL.md).

3. Compile GEMMA from source, see [INSTALL.md](INSTALL.md).

Compiling from source takes more work, but can potentially boost
performance of GEMMA when using specialized C++ compilers and
numerical libraries.

### Precompiled binaries

1. Fetch the [latest stable release][latest_release] and download the
   file appropriate for your platform.

2. For Docker images, install Docker, load the image into Docker and
   run with something like

        docker run -w /run -v ${PWD}:/run ed5bf7499691 gemma -gk -bfile example/mouse_hs1940

3. For .gz files run `gunzip gemma.linux.gz` or `gunzip
gemma.linux.gz` to unpack the file. And make sure it is executable with

        chmod u+x gemma-linux
        ./gemma-linux

## Run GEMMA

GEMMA is run from the command line. To run gemma

```sh
gemma -h
```

a typical example would be

```sh
# compute Kinship matrix
gemma -g ../example/mouse_hs1940.geno.txt.gz -p ../example/mouse_hs1940.pheno.txt \
    -gk -o mouse_hs1940
# run univariate LMM
gemma -g ../example/mouse_hs1940.geno.txt.gz \
    -p ../example/mouse_hs1940.pheno.txt -n 1 -a ../example/mouse_hs1940.anno.txt \
    -k ./output/mouse_hs1940.cXX.txt -lmm -o mouse_hs1940_CD8_lmm
```

Above example files can be downloaded from
[github](https://github.com/genetics-statistics/GEMMA/tree/master/example).

### Debugging and optimization

GEMMA has a wide range of debugging options which can be viewed with

```
 DEBUG OPTIONS

 -check                   enable checks (slower)
 -no-fpe-check            disable hardware floating point checking
 -strict                  strict mode will stop when there is a problem
 -silence                 silent terminal display
 -debug                   debug output
 -debug-data              debug data output
 -nind       [num]        read up to num individuals
 -issue      [num]        enable tests relevant to issue tracker
 -legacy                  run gemma in legacy mode
```

typically when running gemma you should use -debug which includes
relevant checks. When compiled for debugging the debug version of
GEMMA gives more information.

For performance you may want to use the -no-check option. Also check
the build optimization notes in [INSTALL.md](INSTALL.md).

## Help

+ [The GEMMA manual](doc/manual.pdf).

+ [Detailed example with HS mouse data](example/demo.txt).

+ [Tutorial on GEMMA for genome-wide association
analysis](https://github.com/rcc-uchicago/genetic-data-analysis-2).

## Citing GEMMA

If you use GEMMA for published work, please cite our paper:

+ Xiang Zhou and Matthew Stephens (2012). [Genome-wide efficient
mixed-model analysis for association studies.](http://doi.org/10.1038/ng.2310)
*Nature Genetics* **44**, 821–824.

If you use the multivariate linear mixed model (mvLMM) in your
research, please cite:

+ Xiang Zhou and Matthew Stephens (2014). [Efficient multivariate linear
mixed model algorithms for genome-wide association
studies.](http://doi.org/10.1038/nmeth.2848)
*Nature Methods* **11**, 407–409.

If you use the Bayesian sparse linear mixed model (BSLMM), please cite:

+ Xiang Zhou, Peter Carbonetto and Matthew Stephens (2013). [Polygenic
modeling with bayesian sparse linear mixed
models.](http://doi.org/10.1371/journal.pgen.1003264) *PLoS Genetics*
**9**, e1003264.

And if you use of the variance component estimation using summary
statistics, please cite:

+ Xiang Zhou (2016). [A unified framework for variance component
estimation with summary statistics in genome-wide association
studies.](https://doi.org/10.1101/042846) *Annals of Applied Statistics*, in press.

## License

Copyright (C) 2012–2020, Xiang Zhou and team.

The *GEMMA* source code repository is free software: you can
redistribute it under the terms of the
[GNU General Public License](http://www.gnu.org/licenses/gpl.html). All
the files in this project are part of *GEMMA*. This project is
distributed in the hope that it will be useful, but **without any
warranty**; without even the implied warranty of **merchantability or
fitness for a particular purpose**. See file [LICENSE](LICENSE) for
the full text of the license.

Both the source code for the
[gzstream zlib wrapper](http://www.cs.unc.edu/Research/compgeom/gzstream/)
and [shUnit2](https://github.com/genenetwork/shunit2) unit testing
framework included in GEMMA are distributed under the
[GNU Lesser General Public License](contrib/shunit2-2.0.3/doc/LGPL-2.1),
either version 2.1 of the License, or (at your option) any later
revision.

The source code for the included [Catch](http://catch-lib.net) unit
testing framework is distributed under the
[Boost Software Licence version 1](https://github.com/philsquared/Catch/blob/master/LICENSE.txt).

## Optimizing performance

Precompiled binaries and libraries may not be optimal for your particular
hardware. See [INSTALL.md](INSTALL.md) for speeding up tips.

## Building from source

More information on source code, dependencies and installation can be
found in [INSTALL.md](INSTALL.md).

## Input data formats

Currently GEMMA takes two types of input formats

1. BIMBAM format (preferred)
2. PLINK format

See this [example](./doc/example/data-munging.org) where we convert some
spreadsheets for use in GEMMA.

## Reporting a GEMMA bug or issue

For bugs GEMMA has an
[issue tracker](https://github.com/genetics-statistics/GEMMA/issues)
on github. For general support GEMMA has a mailing list at
[gemma-discussion](https://groups.google.com/forum/#!forum/gemma-discussion)

Before posting an issue search the issue tracker and mailing list
first. It is likely someone may have encountered something
similiar. Also try running the latest version of GEMMA to make sure it
has not been fixed already. Support/installation questions should be
aimed at the mailing list - it is the best resource to get answers.

The issue tracker is specifically meant for development issues around
the software itself. When reporting an issue include the output of the
program and the contents of the .log.txt file in the output directory.

### Check list:

1. [X] I have found an issue with GEMMA
2. [ ] I have searched for it on the [issue tracker](https://github.com/genetics-statistics/GEMMA/issues?q=is%3Aissue) (incl. closed issues)
3. [ ] I have searched for it on the [mailing list](https://groups.google.com/forum/#!forum/gemma-discussion)
4. [ ] I have tried the latest [release](https://github.com/genetics-statistics/GEMMA/releases) of GEMMA
5. [ ] I have read and agreed to below code of conduct
6. [ ] If it is a support/install question I have posted it to the [mailing list](https://groups.google.com/forum/#!forum/gemma-discussion)
7. [ ] If it is software development related I have posted a new issue on the [issue tracker](https://github.com/genetics-statistics/GEMMA/issues) or added to an existing one
8. [ ] In the message I have included the output of my GEMMA run
9. [ ] In the message I have included the relevant .log.txt file in the output directory
10. [ ] I have made available the data to reproduce the problem (optional)

To find bugs the GEMMA software developers may ask to install a
development version of the software. They may also ask you for your
data and will treat it confidentially.  Please always remember that
GEMMA is written and maintained by volunteers with good
intentions. Our time is valuable too. By helping us as much as
possible we can provide this tool for everyone to use.

## Code of conduct

By using GEMMA and communicating with its communtity you implicitely
agree to abide by the
[code of conduct](https://software-carpentry.org/conduct/) as
published by the Software Carpentry initiative.

## Credits

The *GEMMA* software was developed by:

[Xiang Zhou](http://www.xzlab.org)<br>
Dept. of Biostatistics<br>
University of Michigan<br>

Peter Carbonetto, Tim Flutre, Matthew Stephens,
[Pjotr Prins](http://thebird.nl/) and
[others](https://github.com/genetics-statistics/GEMMA/graphs/contributors)
have also contributed to the development of this software.

[latest_release]: https://github.com/genetics-statistics/GEMMA/releases "Most recent stable releases"
