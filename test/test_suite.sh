#!/usr/bin/env bash

gemma=../bin/gemma

testCenteredRelatednessMatrixK() {
    $gemma -g ../example/mouse_hs1940.geno.txt.gz -p ../example/mouse_hs1940.pheno.txt \
           -a ../example/mouse_hs1940.anno.txt -gk -o mouse_hs1940
    assertEquals 0 $?
    grep "total computation time" < output/mouse_hs1940.log.txt
    assertEquals 0 $?
    assertEquals "3763600" `wc -w < output/mouse_hs1940.cXX.txt`
    # assertEquals "15f680c" `md5sum < output/mouse_hs1940.cXX.txt | head -c 7`
    assertEquals "0.335" `head -c 5 output/mouse_hs1940.cXX.txt`
    # FIXME: The following test fails in the Guix build system (https://github.com/xiangzhou/GEMMA/issues/55)
    assertEquals "24.9799" `awk '{s+=substr($1,0,6)}END{print s}' output/mouse_hs1940.cXX.txt`
}

testUnivariateLinearMixedModel() {
    $gemma -g ../example/mouse_hs1940.geno.txt.gz -p ../example/mouse_hs1940.pheno.txt -n 1 \
           -a ../example/mouse_hs1940.anno.txt -k ./output/mouse_hs1940.cXX.txt -lmm \
           -o mouse_hs1940_CD8_lmm
    assertEquals 0 $?
    grep "total computation time" < output/mouse_hs1940_CD8_lmm.log.txt
    assertEquals 0 $?
    assertEquals "118459" `wc -w < output/mouse_hs1940_CD8_lmm.assoc.txt`
    assertEquals "92047" `awk '{s+=substr($1,0,6)}END{print s}' output/mouse_hs1940_CD8_lmm.assoc.txt`
}

testMultivariateLinearMixedModel() {
    $gemma -g ../example/mouse_hs1940.geno.txt.gz -p ../example/mouse_hs1940.pheno.txt \
           -n 1 6 -a ../example/mouse_hs1940.anno.txt -k ./output/mouse_hs1940.cXX.txt \
           -lmm -o mouse_hs1940_CD8MCH_lmm
    assertEquals 0 $?
    grep "total computation time" < output/mouse_hs1940_CD8MCH_lmm.log.txt
    assertEquals 0 $?

    outfn=output/mouse_hs1940_CD8MCH_lmm.assoc.txt
    assertEquals "139867" `wc -w < $outfn`
    assertEquals "92079" `awk '{s+=substr($1,0,6)}END{print s}' $outfn`
}

shunit2=`which shunit2`
if [ -e "../contrib/shunit2/source/2.0/src/shell/shunit2" ]; then
    echo try to run the locally installed shunit2
    . ../contrib/shunit2/source/2.0/src/shell/shunit2
elif [ -e "shunit2-2.0.3/src/shell/shunit2" ]; then
    echo try to run the older locally installed shunit2
    . shunit2-2.0.3/src/shell/shunit2
elif [ -x "$shunit2" ]; then
    echo run system shunit2
    . $shunit2
else
    echo "Can not find shunit2 - see INSTALL.md"
fi
