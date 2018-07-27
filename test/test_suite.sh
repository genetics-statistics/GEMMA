#!/usr/bin/env bash

gemma=../bin/gemma
gemmaopts="-debug"

testBslmm1() {
    outn=mouse_hs1940_CD8_bslmm
    $gemma $gemmaopts -g ../example/mouse_hs1940.geno.txt.gz \
           -p ../example/mouse_hs1940.pheno.txt \
           -n 2 -a ../example/mouse_hs1940.anno.txt \
           -bslmm -no-check \
           -o $outn -w 1000 -s 10000 -seed 1
    assertEquals 0 $?
    outfn1=output/$outn.hyp.txt
    outfn2=output/$outn.param.txt
    # assertEquals "45181" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.0f",(substr($x,,0,6))) } END { printf "%.0f",$sum }' $outfn1`
    # assertEquals "4043967139.42" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn2`
}

testBslmm2() {
    outn=mouse_hs1940_CD8_train
    $gemma $gemmaopts -g ../example/mouse_hs1940.geno.txt.gz \
           -p ../example/mouse_hs1940.pheno.txt \
           -n 2 \
           -a ../example/mouse_hs1940.anno.txt \
           -gk 1 -o $outn
    assertEquals 0 $?
    outfn=output/$outn.cXX.txt
    assertEquals "579.66" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

testBslmm3() {
    ## Fit a binary trait using a linear model
    outn=mouse_hs1940_CD8_bslmm_cc1
    $gemma $gemmaopts \
           -g ../example/mouse_hs1940.geno.txt.gz \
           -p ../example/mouse_hs1940.pheno.txt \
           -n 4 \
           -a ../example/mouse_hs1940.anno.txt \
           -bslmm \
           -o $outn \
           -w 1000 -s 10000 -seed 1 -no-check
    assertEquals 0 $?
    outfn=output/$outn.hyp.txt
    # assertEquals "291" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.0f",(substr($x,,0,6))) } END { printf "%.0f",100*$sum }' $outfn`
}

testBslmm4() {
    outn=mouse_hs1940_CD8_prdt_k
    $gemma $gemmaopts -g ../example/mouse_hs1940.geno.txt.gz \
           -p ../example/mouse_hs1940.pheno.txt \
           -n 2 \
           -epm ./output/mouse_hs1940_CD8_bslmm.param.txt \
           -emu ./output/mouse_hs1940_CD8_bslmm.log.txt \
           -ebv ./output/mouse_hs1940_CD8_bslmm.bv.txt \
           -k ./output/mouse_hs1940_CD8_train.cXX.txt \
           -predict -no-check \
           -o $outn
    assertEquals 0 $?
    outfn=output/$outn.prdt.txt
    # assertEquals "-60.33" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

testBslmm5() {
    ## Now, do prediction in the test set for the binary traits
    ## If the traits were fitted using the linear model, then:
    outn=mouse_hs1940_CD8_prdt_cc1
    $gemma $gemmaopts \
           -g ../example/mouse_hs1940.geno.txt.gz \
           -p ../example/mouse_hs1940.pheno.txt \
           -n 4 \
           -epm ./output/mouse_hs1940_CD8_bslmm_cc1.param.txt \
           -emu ./output/mouse_hs1940_CD8_bslmm_cc1.log.txt \
           -predict \
           -o $outn
    assertEquals 0 $?
    outfn=output/$outn.prdt.txt
    assertEquals "550.67" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

testCenteredRelatednessMatrixKFullLOCO1() {
    outn=mouse_hs1940_full_LOCO1
    $gemma $gemmaopts -g ../example/mouse_hs1940.geno.txt.gz \
           -p ../example/mouse_hs1940.pheno.txt \
           -a ../example/mouse_hs1940.anno.txt \
           -loco 1 -gk -o $outn
    assertEquals 0 $?
    outfn=output/$outn.cXX.txt
    assertEquals "1940" `wc -l < $outfn`
    assertEquals "2246.57" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

testUnivariateLinearMixedModelFullLOCO1() {
    outn=mouse_hs1940_CD8_full_LOCO1_lmm
    $gemma $gemmaopts -g ../example/mouse_hs1940.geno.txt.gz \
           -p ../example/mouse_hs1940.pheno.txt \
	   -n 1 \
	   -loco 1 \
           -a ../example/mouse_hs1940.anno.txt \
           -k ./output/mouse_hs1940_full_LOCO1.cXX.txt \
	   -lmm -no-check \
           -o $outn
    assertEquals 0 $?
    grep "total computation time" < output/$outn.log.txt
    assertEquals 0 $?
    outfn=output/$outn.assoc.txt
    assertEquals "951" `wc -l < $outfn`
    assertEquals "267507851.98" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

testCenteredRelatednessMatrixK() {
    $gemma $gemmaopts -g ../example/mouse_hs1940.geno.txt.gz \
           -p ../example/mouse_hs1940.pheno.txt \
           -gk -o mouse_hs1940
    assertEquals 0 $?
    outfn=output/mouse_hs1940.cXX.txt
    assertEquals "1940" `wc -l < $outfn`
    assertEquals "3763600" `wc -w < $outfn`
    assertEquals "0.335" `head -c 5 $outfn`
    assertEquals "1119.64" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

testUnivariateLinearMixedModel() {
    $gemma $gemmaopts -g ../example/mouse_hs1940.geno.txt.gz \
           -p ../example/mouse_hs1940.pheno.txt \
           -n 1 \
           -a ../example/mouse_hs1940.anno.txt \
           -k ./output/mouse_hs1940.cXX.txt \
           -lmm \
           -o mouse_hs1940_CD8_lmm
    assertEquals 0 $?
    grep "total computation time" < output/mouse_hs1940_CD8_lmm.log.txt
    assertEquals 0 $?
    outfn=output/mouse_hs1940_CD8_lmm.assoc.txt
    assertEquals "129228" `wc -w < $outfn`
    assertEquals "4038540440.86" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

testLinearMixedModelPhenotypes() {
    $gemma $gemmaopts -g ../example/mouse_hs1940.geno.txt.gz \
           -p ../example/mouse_hs1940.pheno.txt \
           -n 1 6 \
           -a ../example/mouse_hs1940.anno.txt \
           -k ./output/mouse_hs1940.cXX.txt \
           -lmm -no-check \
           -o mouse_hs1940_CD8MCH_lmm
    assertEquals 0 $?

    outfn=output/mouse_hs1940_CD8MCH_lmm.assoc.txt
    assertEquals "139867" `wc -w < $outfn`
    assertEquals "4029037056.54" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

testPlinkStandardRelatednessMatrixK() {
    testname=testPlinkStandardRelatednessMatrixK
    datadir=../example
    outfn=output/$testname.sXX.txt
    rm -f $outfn
    $gemma $gemmaopts -bfile $datadir/HLC \
           -gk 2 -o $testname
    assertEquals 0 $?
    assertEquals "427" `wc -l < $outfn`
    assertEquals "-358.07" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

# Test for https://github.com/genetics-statistics/GEMMA/issues/58
# fixed GSLv2 NaN's that appeared with covariates.
testPlinkLinearMixedModelCovariates() {
    testname=testPlinkLinearMixedModelCovariates
    datadir=../example
    $gemma $gemmaopts -bfile $datadir/HLC \
           -k output/testPlinkStandardRelatednessMatrixK.sXX.txt \
           -lmm 1 \
           -maf 0.1 \
           -c $datadir/HLC_covariates.txt \
           -no-check \
           -o $testname
    assertEquals 0 $?
    outfn=output/$testname.assoc.txt
    assertEquals "223243" `wc -l < $outfn`
    assertEquals "89757159113.77" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
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
