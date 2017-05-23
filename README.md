# GEMMA: Genome-wide Efficient Mixed Model Association

GEMMA is a software toolkit for fast application of linear mixed
models and related models to genome-wide association studies (GWAS)
and other large-scale data sets.

![Genetic associations discovered in CFW mice using GEMMA (Parker et al,
Nat. Genet., 2016)](cfw.gif)

Features include:

+ Fast assocation tests implemented using the univariate linear mixed
model (LMM). In GWAS, this can correct for account for population
stratification and sample nonexchangeability. It also provides
estimates of the proportion of variance in phenotypes explained (PVE)
by available genotypes (often called "chip heritability" or "SNP
heritability").

+ Fast association tests for multiple phenotypes implemented using a
multivariate linear mixed model (lvLMM).

It fits a multivariate linear mixed model (mvLMM) for testing marker
associations with multiple phenotypes simultaneously while controlling
for population stratification, and for estimating genetic correlations
among complex phenotypes.

+ It fits a Bayesian sparse linear mixed model (BSLMM) using Markov
chain Monte Carlo (MCMC) for estimating PVE by typed genotypes,
predicting phenotypes, and identifying associated markers by jointly
modeling all markers while controlling for population structure.

+ It estimates variance component/chip heritability, and partitions it
by different SNP functional categories. In particular, it uses HE
regression or REML AI algorithm to estimate variance components when
individual-level data are available. It uses MQS to estimate variance
components when only summary statisics are available.

*Add note here about posting questions, comments or bug reports to
Issues.*

### Citing GEMMA

*Add text here.*

### License

Copyright (C) 2012–2017, Xiang Zhou.

### Quick start

*Add text here.*

### Setup

*Add text here.*

