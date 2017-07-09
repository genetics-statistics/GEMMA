#!/usr/bin/env bash

gemma=../bin/gemma

testCenteredRelatednessMatrix() {
    $gemma -g ../example/mouse_hs1940.geno.txt.gz -p ../example/mouse_hs1940.pheno.txt \
                 -a ../example/mouse_hs1940.anno.txt -gk -o mouse_hs1940
    assertEquals "3763600" `wc -w < output/mouse_hs1940.cXX.txt`
}

shunit2=`which shunit2`
if [ -x "$shunit2" ]; then
    . $shunit2
else
    # try to run the locally installed shunit2
    . ../shunit2-2.0.3/src/shell/shunit2
fi
