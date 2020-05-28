#!/bin/sh

gemma=../bin/gemma
# gemmaopts="-debug -strict"
gemmaopts="-debug -check"

testLinearModel() {
    $gemma $gemmaopts -g ../example/mouse_hs1940.geno.txt.gz \
           -p ../example/mouse_hs1940.pheno.txt \
           -n 1 \
           -a ../example/mouse_hs1940.anno.txt \
           -lm \
           -o mouse_hs1940_CD8_lm
    assertEquals 0 $?
    grep "total computation time" < output/mouse_hs1940_CD8_lm.log.txt
    assertEquals 0 $?
    outfn=output/mouse_hs1940_CD8_lm.assoc.txt
    assertEquals "118459" `wc -w < $outfn`
    assertEquals "4053667109.67" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

# Related to https://github.com/genetics-statistics/GEMMA/issues/78
testBXDStandardRelatednessMatrixKSingularError() {
    outn=BXDerr
    rm -f output/$outn.*
    $gemma $gemmaopts \
           -g ../example/BXD_geno.txt.gz \
           -p ../example/BXD_pheno.txt \
           -c ../example/BXD_covariates.txt \
           -a ../example/BXD_snps.txt \
           -gk \
           -no-check \
           -o $outn
    # assertEquals 130 $? # should show singular error FIXME
}

testBXDStandardRelatednessMatrixK() {
    outn=BXD
    rm -f output/$outn.*
    $gemma $gemmaopts -g ../example/BXD_geno.txt.gz \
           -p ../example/BXD_pheno.txt \
           -c ../example/BXD_covariates2.txt \
           -a ../example/BXD_snps.txt \
           -gk \
           -o $outn
    assertEquals 0 $?
    outfn=output/$outn.cXX.txt
    assertEquals "198" `wc -l < $outfn`
    assertEquals "-116.13" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

testBXDLMLikelihoodRatio() {
    outn=BXD_LM_LR
    $gemma $gemmaopts -g ../example/BXD_geno.txt.gz \
           -p ../example/BXD_pheno.txt \
           -c ../example/BXD_covariates2.txt \
           -a ../example/BXD_snps.txt \
           -k ./output/BXD.cXX.txt \
           -lm 4 -maf 0.1 \
           -o $outn
    assertEquals 0 $?

    outfn=output/$outn.assoc.txt
    assertEquals "95134" `wc -w < $outfn`
    assertEquals "3089042886.08" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

testBXDLMMLikelihoodRatio() {
    outn=BXD_LMM_LR
    $gemma $gemmaopts -g ../example/BXD_geno.txt.gz \
           -p ../example/BXD_pheno.txt \
           -c ../example/BXD_covariates2.txt \
           -a ../example/BXD_snps.txt \
           -k ./output/BXD.cXX.txt \
           -lmm 2 -no-check -maf 0.1 \
           -o $outn
    assertEquals 0 $?

    outfn=output/$outn.assoc.txt
    assertEquals "73180" `wc -w < $outfn`
    assertEquals "3088458212.86" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

testCenteredRelatednessMatrixissue188() {
    outn=issue188
    rm -f output/$outn.*
    $gemma $gemmaopts -b data/issue188/2000 -gk -o $outn
    assertEquals 0 $?
    outfn=output/$outn.cXX.txt
    assertEquals "193.78" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

testLMMissue188() {
    outn=issue188
    $gemma $gemmaopts -b data/issue188/2000 -lmm 2 -k output/$outn.cXX.txt -maf 0.01 -o $outn -n 1
    assertEquals 0 $?
    outfn=output/$outn.assoc.txt
    assertEquals "338154001.76" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

testCenteredRelatednessMatrixKLOCO1() {
    outn=mouse_hs1940_LOCO1
    rm -f output/$outn.*
    $gemma $gemmaopts -g ../example/mouse_hs1940.geno.txt.gz -p ../example/mouse_hs1940.pheno.txt \
           -a ../example/mouse_hs1940.anno.txt -snps ../example/mouse_hs1940_snps.txt -nind 400 -loco 1 -gk -o $outn
    assertEquals 0 $?
    grep "total computation time" < output/$outn.log.txt
    outfn=output/$outn.cXX.txt
    assertEquals 0 $?
    assertEquals "400" `wc -l < $outfn`
    assertEquals "0.312" `head -c 5 $outfn`
    assertEquals "71.03" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

testUnivariateLinearMixedModelLOCO1() {
    outn=mouse_hs1940_CD8_LOCO1_lmm
    rm -f output/$outn.*
    $gemma $gemmaopts -g ../example/mouse_hs1940.geno.txt.gz \
           -p ../example/mouse_hs1940.pheno.txt \
	   -n 1 \
	   -loco 1 \
           -a ../example/mouse_hs1940.anno.txt \
           -k ./output/mouse_hs1940_LOCO1.cXX.txt \
	   -snps ../example/mouse_hs1940_snps.txt -lmm \
	   -nind 400 -no-check \
           -o $outn
    assertEquals 0 $?
    grep "total computation time" < output/$outn.log.txt
    assertEquals 0 $?
    outfn=output/$outn.assoc.txt
    assertEquals "68" `wc -l < $outfn`
    assertEquals "15465346.22" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

testPlinkCenteredRelatednessMatrixKLOCO1() {
    return 0
    outn=mouse_hs1940_Plink_LOCO1
    rm -f output/$outn.*
    $gemma $gemmaopts -bfile ../example/mouse_hs1940 \
           -a ../example/mouse_hs1940.anno.txt \
           -snps ../example/mouse_hs1940_snps.txt \
           -nind 400 \
           -loco 1 \
           -gk \
           -o $outn
    assertEquals 0 $?
    grep "total computation time" < output/$outn.log.txt
    outfn=output/$outn.cXX.txt
    assertEquals 0 $?
    assertEquals "400" `wc -l < $outfn`
    assertEquals "0.312" `head -c 5 $outfn`
    assertEquals "71.03" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}


testPlinkUnivariateLinearMixedModelLOCO1() {
    return 0
    outn=mouse_hs1940_CD8_Plink_LOCO1_lmm
    rm -f output/$outn.*
    $gemma $gemmaopts -bfile ../example/mouse_hs1940 \
	   -n 1 \
	   -loco 1 \
           -k ./output/mouse_hs1940_Plink_LOCO1.cXX.txt \
           -a ../example/mouse_hs1940.anno.txt \
	   -snps ../example/mouse_hs1940_snps.txt -lmm \
	   -nind 400 \
           -o $outn
    assertEquals 0 $?
    grep "total computation time" < output/$outn.log.txt
    assertEquals 0 $?
    outfn=output/$outn.assoc.txt
    assertEquals "68" `wc -l < $outfn`
    assertEquals "15465346.22" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}


testCorrelatedPhenotypesMvLLM() {
    # https://github.com/genetics-statistics/GEMMA/issues/179
    outn=corrpheno
    $gemma $gemmaopts -p data/correlated_phenotypes/Ysim_reg_gemma.txt \
           -g data/correlated_phenotypes/Genotypes_gemma.csv \
           -d data/correlated_phenotypes/Kinship_eigenval_gemma.txt \
           -u data/correlated_phenotypes/Kinship_eigenvec_gemma.txt \
           -lmm 2 -n 1 9 4 6 10 -o $outn
    assertEquals 0 $?
    outfn=output/$outn.assoc.txt
    assertEquals "101" `wc -l < $outfn`
    # assertEquals "777.32" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}


shunit2=`which shunit2`

if [ -x "$shunit2" ]; then
    echo run system shunit2
    . $shunit2
elif [ -e ../contrib/shunit2-2.0.3/src/shell/shunit2 ]; then
    echo run shunit2 provided in gemma repo
    . ../contrib/shunit2-2.0.3/src/shell/shunit2
else
    echo "Can not find shunit2 - see INSTALL.md"
fi
