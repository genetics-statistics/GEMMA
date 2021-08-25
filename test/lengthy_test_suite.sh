#!/usr/bin/env bash
#
# Long running tests go here

echo "WARNING: THIS TEST SUITE IS NO LONGER USED"

gemma=../bin/gemma
export GSL_RNG_SEED=100

testPlinkStandardRelatednessMatrixK() {
    testname=testPlinkStandardRelatednessMatrixK
    datadir=../example
    outfn=output/$testname.sXX.txt
    rm -f $outfn
    $gemma -bfile $datadir/HLC \
           -gk 2 -o $testname \
           -debug
    assertEquals 0 $?
    assertEquals "427" `wc -l < $outfn`
    assertEquals "-358.07" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
}

testPlinkMultivariateLinearMixedModelMultiplePhenotypes_Issue58() {
    echo "Long running test!"
    # This test passes, but takes over 30 minutes to run!
    # n=2 is original pheno in fam file
    # n=1 is causal1
    # n=3..12 is causal2
    # n=13..22 is causal3
    # -n 1 2 3 15 is independent
    testname=testPlinkMultivariateLinearMixedModelMultiplePhenotypes
    datadir=../example
    $gemma -bfile $datadir/HLC \
           -p $datadir/HLC.simu.pheno.txt \
           -k output/testPlinkStandardRelatednessMatrixK.sXX.txt \
           -lmm 1 \
           -maf 0.1 \
           -n 1 2 3 15 \
           -c $datadir/HLC_covariates.txt \
           -debug \
           -o $testname
    assertEquals 0 $?
    outfn=output/$testname.assoc.txt
    assertEquals "223243" `wc -l < $outfn`
    assertEquals "89754977983.69" `perl -nle 'foreach $x (split(/\s+/,$_)) { $sum += sprintf("%.2f",(substr($x,,0,6))) } END { printf "%.2f",$sum }' $outfn`
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
