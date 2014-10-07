## Detailed description of the data set is available in the online GEMMA user manual
## Each of the following steps may take over one minute to run






## To calculate a centered relatedness matrix:
../bin/gemma -g mouse_hs1940.geno.txt.gz -p mouse_hs1940.pheno.txt -a mouse_hs1940.anno.txt -gk -o mouse_hs1940

# The estimated relatedness matrix should look like this:
0.3350590  -0.0227226  0.0103535 ...
-0.0227226  0.3035960 -0.0253762 ...
0.0103535  -0.0253762  0.3536100 ...
....................................










## To perform association tests with a univariate linear mixed model:

../bin/gemma -g mouse_hs1940.geno.txt.gz -p mouse_hs1940.pheno.txt -n 1 -a mouse_hs1940.anno.txt -k ./output/mouse_hs1940.cXX.txt -lmm -o mouse_hs1940_CD8_lmm

# The result for top 5 SNPs should look like this:
chr	rs	ps	n_miss	allele1	allele0	af	beta	se	l_remle	p_wald
1	rs3683945	3197400	0	A	G	0.443	-7.788665e-02	6.193502e-02	4.317993e+00	2.087616e-01
1	rs3707673	3407393	0	G	A	0.443	-6.654282e-02	6.210234e-02	4.316144e+00	2.841271e-01
1	rs6269442	3492195	0	A	G	0.365	-5.344241e-02	5.377464e-02	4.323611e+00	3.204804e-01
1	rs6336442	3580634	0	A	G	0.443	-6.770154e-02	6.209267e-02	4.315713e+00	2.757541e-01
1	rs13475700	4098402	0	A	C	0.127	-5.659089e-02	7.175374e-02	4.340145e+00	4.304306e-01

# If you do a manhattan plot in R, you will see a strong signal in chr17

# The log file also contains pve estimates and its standard error
pve estimate in the null model = 0.608801
se(pve) in the null model = 0.032774











## To perform association tests with a multivariate linear mixed model, for two phenotypes CD8 (column 1) and MCH (column 6):
## Notice that the number of individuals in this analysis is different from that above, so the allele frequencies are different between the two analyses

../bin/gemma -g mouse_hs1940.geno.txt.gz -p mouse_hs1940.pheno.txt -n 1 6 -a mouse_hs1940.anno.txt -k ./output/mouse_hs1940.cXX.txt -lmm -o mouse_hs1940_CD8MCH_lmm

# The result for top 5 SNPs should look like this:
chr	rs	ps	n_miss	allele1	allele0	af	beta_1	beta_2	Vbeta_1_1	Vbeta_1_2	Vbeta_2_2	p_wald
1	rs3683945	3197400	0	A	G	0.451	-9.611213e-02	8.165302e-02	3.966873e-03	-2.526118e-04	5.540032e-03	1.862363e-01
1	rs3707673	3407393	0	G	A	0.451	-8.464470e-02	7.130876e-02	3.986286e-03	-2.593467e-04	5.571616e-03	2.757067e-01
1	rs6269442	3492195	0	A	G	0.377	-7.146771e-02	5.179252e-02	3.157023e-03	-7.187157e-05	4.276041e-03	3.317712e-01
1	rs6336442	3580634	0	A	G	0.451	-8.502513e-02	6.813728e-02	3.985054e-03	-2.577585e-04	5.568602e-03	2.835426e-01
1	rs13475700	4098402	0	A	C	0.128	-6.727883e-02	1.685363e-01	5.597160e-03	-1.366799e-04	7.574216e-03	1.060482e-01

# The log file also contains Vg and Ve estimates and their standard errors
## REMLE estimate for Vg in the null model: 
1.39398	
-0.226714	2.08168	
## se(Vg): 
0.156661	
0.136319	0.235858	
## REMLE estimate for Ve in the null model: 
0.348882	
0.0490525	0.414433	
## se(Ve): 
0.0206226	
0.0166233	0.0266869

# Since there are individuals with partially missing phenotypes, one can impute these missing values before association tests
../bin/gemma -g mouse_hs1940.geno.txt.gz -p mouse_hs1940.pheno.txt -n 1 6 -a mouse_hs1940.anno.txt -k ./output/mouse_hs1940.cXX.txt -predict -o mouse_hs1940_CD8MCH_prdt

../bin/gemma -g mouse_hs1940.geno.txt.gz -p ./output/mouse_hs1940_CD8MCH_prdt.prdt.txt -n 1 2 -a mouse_hs1940.anno.txt -k ./output/mouse_hs1940.cXX.txt -lmm -o mouse_hs1940_CD8MCH_prdt_lmm











## To fit BSLMM in the training set:

## To fit a quantitative trait
../bin/gemma -g mouse_hs1940.geno.txt.gz -p mouse_hs1940.pheno.txt -n 2 -a mouse_hs1940.anno.txt -bslmm -o mouse_hs1940_CD8_bslmm -w 1000 -s 10000 -seed 1

# the following three files may be of most importance:
# the *.hyp.txt contains a column for pve and pge
# the *.param.txt contains estimates for betas, gammas and alphas
# the *.bv.txt contains breeding value estimates

## To fit a binary trait using a linear model
../bin/gemma -g mouse_hs1940.geno.txt.gz -p mouse_hs1940.pheno.txt -n 4 -a mouse_hs1940.anno.txt -bslmm -o mouse_hs1940_CD8_bslmm_cc1 -w 1000 -s 10000 -seed 1

## To fit a binary trait using a probit model instead
../bin/gemma -g mouse_hs1940.geno.txt.gz -p mouse_hs1940.pheno.txt -n 4 -a mouse_hs1940.anno.txt -bslmm 3 -o mouse_hs1940_CD8_bslmm_cc3 -w 1000 -s 10000 -seed 1

# The pve estimates in the log file are based on the standard linear model (i.e. on the observed scale), and so you will need to properly transform it back to the liability scale

## To generate relatedness matrix based on the training data
../bin/gemma -g mouse_hs1940.geno.txt.gz -p mouse_hs1940.pheno.txt -n 2 -a mouse_hs1940.anno.txt -gk 1 -o mouse_hs1940_CD8_train

# This matrix will only be required if you want to do prediction based on estimated breeding values
# Prediction can also be done without using the breeding values but instead using the alphas. This later approach does not appear to lose much accuracy in many examples we have encountered, although this may not be the case in your data.










## To obtain predicted values for the test set using estimates from BSLMM
## To do prediction in the test set for quantitative traits
../bin/gemma -g mouse_hs1940.geno.txt.gz -p mouse_hs1940.pheno.txt -n 2 -epm ./output/mouse_hs1940_CD8_bslmm.param.txt -emu ./output/mouse_hs1940_CD8_bslmm.log.txt -ebv ./output/mouse_hs1940_CD8_bslmm.bv.txt -k ./output/mouse_hs1940_CD8_train.cXX.txt -predict -o mouse_hs1940_CD8_prdt_k

## or use the alphas instead of breeding values
../bin/gemma -g mouse_hs1940.geno.txt.gz -p mouse_hs1940.pheno.txt -n 2 -epm ./output/mouse_hs1940_CD8_bslmm.param.txt -emu ./output/mouse_hs1940_CD8_bslmm.log.txt -predict -o mouse_hs1940_CD8_prdt

# The results will be inside ./output/*.prdt.txt
# If you load both results in R and check the mean squared error or correlation, you will find that both ways give very similar results. Both the correlation and the mean squared error should be around 0.65

## Now, do prediction in the test set for the binary traits
## If the traits were fitted using the linear model, then:
../bin/gemma -g mouse_hs1940.geno.txt.gz -p mouse_hs1940.pheno.txt -n 4 -epm ./output/mouse_hs1940_CD8_bslmm_cc1.param.txt -emu ./output/mouse_hs1940_CD8_bslmm_cc1.log.txt -predict -o mouse_hs1940_CD8_prdt_cc1

## If the traits were fitted using the probit model, then use predict option 2:
../bin/gemma -g mouse_hs1940.geno.txt.gz -p mouse_hs1940.pheno.txt -n 4 -epm ./output/mouse_hs1940_CD8_bslmm_cc3.param.txt -emu ./output/mouse_hs1940_CD8_bslmm_cc3.log.txt -predict 2 -o mouse_hs1940_CD8_prdt_cc3

# You will find that fitting the binary traits using either the linear version or the probit version of BSLMM gives similar results. The brier scores should be around 0.19 and the area under the curve (AUC) should be around 0.78.
