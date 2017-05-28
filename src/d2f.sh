#!/bin/bash

for prefix in gemma param io lm lmm mvlmm bslmm mathfunc prdt vc
do
for extension in cpp h
do
cp ${prefix}.${extension} ${prefix}_float.${extension}
sed -i.bak 's/_vector_/_vector_float_/g' ${prefix}_float.${extension}
sed -i.bak 's/_vector /_vector_float /g' ${prefix}_float.${extension}
sed -i.bak 's/_matrix_/_matrix_float_/g' ${prefix}_float.${extension}
sed -i.bak 's/_matrix /_matrix_float /g' ${prefix}_float.${extension}
sed -i.bak 's/ddot/dsdot/g' ${prefix}_float.${extension}
sed -i.bak 's/dtrsv/strsv/g' ${prefix}_float.${extension}
sed -i.bak 's/dtrsy/strsy/g' ${prefix}_float.${extension}
sed -i.bak 's/dgemm/sgemm/g' ${prefix}_float.${extension}
sed -i.bak 's/dgemv/sgemv/g' ${prefix}_float.${extension}
sed -i.bak 's/dsyr/ssyr/g' ${prefix}_float.${extension}
sed -i.bak 's/dsyr2/ssyr2/g' ${prefix}_float.${extension}
sed -i.bak 's/ddot/sdot/g' ${prefix}_float.${extension}
sed -i.bak 's/dger/sger/g' ${prefix}_float.${extension}
sed -i.bak 's/dsyrk/ssyrk/g' ${prefix}_float.${extension}
sed -i.bak 's/daxpy/saxpy/g' ${prefix}_float.${extension}
rm ${prefix}_float.${extension}.bak
done
done
