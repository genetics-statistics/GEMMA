#!/usr/bin/env bash

gemma=../bin/gemma

# Related to https://github.com/genetics-statistics/GEMMA/issues/78
testBXDStandardRelatednessMatrixKSingularError() {
    outn=BXDerr
    rm -f output/$outn.*
    $gemma -g ../example/BXD_geno.txt.gz \
           -p ../example/BXD_pheno.txt \
           -c ../example/BXD_covariates.txt \
           -a ../example/BXD_snps.txt \
           -gk \
           -debug -o $outn
    assertEquals 22 $? # should show singular error
}

testBXDStandardRelatednessMatrixK() {
    outn=BXD
    rm -f output/$outn.*
    $gemma -g ../example/BXD_geno.txt.gz \
           -p ../example/BXD_pheno.txt \
           -c ../example/BXD_covariates2.txt \
           -a ../example/BXD_snps.txt \
           -gk \
           -debug \
           -o $outn
    assertEquals 0 $?
    outfn=output/$outn.cXX.txt
    assertEquals "198" `wc -l < $outfn`
    assertEquals "-116.11" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

testBXDLMMLikelihoodRatio() {
    outn=BXD_LMM_LR
    $gemma -g ../example/BXD_geno.txt.gz \
           -p ../example/BXD_pheno.txt \
           -c ../example/BXD_covariates2.txt \
           -a ../example/BXD_snps.txt \
           -k ./output/BXD.cXX.txt \
           -lmm 2 -maf 0.1 \
           -debug \
           -o $outn
    assertEquals 0 $?

    outfn=output/$outn.assoc.txt
    assertEquals "73180" `wc -w < $outfn`
    assertEquals "3088458212.93" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

testCenteredRelatednessMatrixKLOCO1() {
    outn=mouse_hs1940_LOCO1
    rm -f output/$outn.*
    $gemma -g ../example/mouse_hs1940.geno.txt.gz -p ../example/mouse_hs1940.pheno.txt \
           -a ../example/mouse_hs1940.anno.txt -snps ../example/mouse_hs1940_snps.txt -nind 400 -loco 1 -gk -debug -o $outn
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
    $gemma -g ../example/mouse_hs1940.geno.txt.gz \
           -p ../example/mouse_hs1940.pheno.txt \
	   -n 1 \
	   -loco 1 \
           -a ../example/mouse_hs1940.anno.txt \
           -k ./output/mouse_hs1940_LOCO1.cXX.txt \
	   -snps ../example/mouse_hs1940_snps.txt -lmm \
	   -nind 400 \
	   -debug \
           -o $outn
    assertEquals 0 $?
    grep "total computation time" < output/$outn.log.txt
    assertEquals 0 $?
    outfn=output/$outn.assoc.txt
    assertEquals "68" `wc -l < $outfn`
    assertEquals "15465346.22" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
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
