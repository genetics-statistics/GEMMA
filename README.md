![Genetic associations identified in CFW mice using GEMMA (Parker et al,
Nat. Genet., 2016)](cfw.gif)

# GEMMA: Genome-wide Efficient Mixed Model Association

GEMMA is a software toolkit for fast application of linear mixed
models (LMMs) and related models to genome-wide association studies
(GWAS) and other large-scale data sets.

Please post comments, feature requests or suspected bugs to
[Github issues](https://github.com/xiangzhou/GEMMA/issues).

*Note: The above image depicts physiological and behavioral trait
loci identified in CFW mice using GEMMA ([Parker et al, Nature
Genetics, 2006](https://doi.org/10.1038/ng.3609)).*

## Key features

1. Fast assocation tests implemented using the univariate linear mixed
model (LMM). In GWAS, this can correct for population structure and
sample nonexchangeability. It also provides estimates of the
proportion of variance in phenotypes explained by available genotypes
(PVE), often called "chip heritability" or "SNP heritability".

2. Fast association tests for multiple phenotypes implemented using a
multivariate linear mixed model (mvLMM). In GWAS, this can correct for
populations tructure and sample nonexchangeability jointly in multiple
complex phenotypes.

3. Bayesian sparse linear mixed model (BSLMM) for estimating PVE,
phenotype prediction, and multi-marker modeling in GWAS.

4. Estimation of variance components ("chip heritability") partitioned
by different SNP functional categories from raw (individual-level)
data or summary data. For raw data, HE regression or the REML AI
algorithm can be used to estimate variance components when
individual-level data are available. For summary data, GEMMA uses the
MQS algorithm to estimate variance components.

## Quick start

1. Download and install the software. *Give more details here.*

2. Work through the tutorial. *Give more details here.*

3. Read the manual. *Give more details.*

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
studies.](https://doi.org/10.1101/042846) *bioRxiv* 042846.

## License

Copyright (C) 2012–2017, Xiang Zhou.

The *GEMMA* source code repository is free software: you can
redistribute it under the terms of the
[GNU General Public License](http://www.gnu.org/licenses/gpl.html). All
the files in this project are part of *GEMMA*. This project is
distributed in the hope that it will be useful, but **without any
warranty**; without even the implied warranty of **merchantability or
fitness for a particular purpose**. See file [LICENSE](LICENSE) for
the full text of the license.

## Setup 

### Using precompiled executables

### Building from source

## Credits

*Add text here.*
