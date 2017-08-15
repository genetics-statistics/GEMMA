#!/usr/bin/env bash

gemma=../bin/gemma

# Test for https://github.com/genetics-statistics/GEMMA/issues/26
# Always getting 'pve estimate =0.99xxx se(pve) =-nan'
testIssue26() {
    outn=issue26
    rm -f output/$outn.*
    $gemma -bfile data/issue26/mydata -k data/issue26/mydata_kinship.sXX.txt \
           -miss 1 -maf 0.01 -r2 1 -lmm \
           -debug -issue 26 \
           -o $outn
    assertEquals 0 $?
    outfn=output/$outn.assoc.txt
    grep "total computation time" < output/$outn.log.txt
    assertEquals 0 $?
    assertEquals "2001" `wc -l < $outfn`
    assertEquals "1582899231.18" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
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
    $gemma -g ../example/mouse_hs1940.geno.txt.gz -p ../example/mouse_hs1940.pheno.txt \
	   -n 1 \
	   -loco 1 \
           -a ../example/mouse_hs1940.anno.txt -k ./output/mouse_hs1940_LOCO1.cXX.txt \
	   -snps ../example/mouse_hs1940_snps.txt -lmm \
	   -nind 400 \
	   -debug \
           -o $outn
    assertEquals 0 $?
    grep "total computation time" < output/$outn.log.txt
    assertEquals 0 $?
    outfn=output/$outn.assoc.txt
    assertEquals "68" `wc -l < $outfn`
    assertEquals "15465553.30" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

shunit2=`which shunit2`

if [ -x "$shunit2" ]; then
    echo run system shunit2
    . $shunit2
elif [ -e shunit2-2.0.3/src/shell/shunit2 ]; then
    echo run shunit2 provided in gemma repo
    . shunit2-2.0.3/src/shell/shunit2
else
    echo "Can not find shunit2 - see INSTALL.md"
fi
