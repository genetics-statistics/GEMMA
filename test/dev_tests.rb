
require 'minitest/autorun'
require 'lib/test-helpers'

class TestQuick < MiniTest::Test

  include TestHelpers

  def setup
    ENV['GSL_RNG_SEED']=100.to_s
  end

  def test_linear_model
    gemma("-g ./example/mouse_hs1940.geno.txt.gz \
           -p ./example/mouse_hs1940.pheno.txt \
           -n 1 \
           -a ./example/mouse_hs1940.anno.txt \
           -lm \
           -o mouse_hs1940_CD8_lm")

    expect('output/mouse_hs1940_CD8_lm.assoc.txt',[[2,1,"rs3707673"],
                                                   [2,10,"5.252187e-05"],
                                                   [3,9,"3.909916e-02"]])
  end

  def test_BXD
    gemma("-g ./example/BXD_geno.txt.gz \
           -p ./example/BXD_pheno.txt \
           -c ./example/BXD_covariates2.txt \
           -a ./example/BXD_snps.txt \
           -gk \
           -o BXD")

    gemma("-g ./example/BXD_geno.txt.gz \
           -p ./example/BXD_pheno.txt \
           -c ./example/BXD_covariates2.txt \
           -a ./example/BXD_snps.txt \
           -k ./output/BXD.cXX.txt \
           -lmm 2 -no-check -maf 0.1 \
           -o BXDLMM")

    expect("output/BXDLMM.assoc.txt",[[2,9,"1.234747e-01"],
                                      [:max,"p_lrt","9.997119e-01"]])

    gemma("-g ./example/BXD_geno.txt.gz \
           -p ./example/BXD_pheno.txt \
           -c ./example/BXD_covariates2.txt \
           -a ./example/BXD_snps.txt \
           -k ./output/BXD.cXX.txt \
           -lmm 9 -no-check -maf 0.1 \
           -o BXDLMM9")

    expect("output/BXDLMM9.assoc.txt",[[:max,"l_mle","0.7531109"],
                                       [:max,"p_lrt","9.997119e-01"]])
  end

  def test_mouse_hs1940_loco
    gemma("-g ./example/mouse_hs1940.geno.txt.gz \
           -p ./example/mouse_hs1940.pheno.txt \
           -a ./example/mouse_hs1940.anno.txt \
           -snps ./example/mouse_hs1940_snps.txt \
           -nind 400 -loco 1 -gk \
           -o mouse_hs1940_loco")

    gemma("-g ./example/mouse_hs1940.geno.txt.gz \
           -p ./example/mouse_hs1940.pheno.txt \
	   -n 1 \
	   -loco 1 \
           -a ./example/mouse_hs1940.anno.txt \
           -k ./output/mouse_hs1940_loco.cXX.txt \
	   -snps ./example/mouse_hs1940_snps.txt \
           -lmm \
	   -nind 400 -no-check \
           -o mouse_hs1940_loco")
    expect("output/mouse_hs1940_loco.assoc.txt",[[2,9,"-3.073643e+02"],
                                                 [:max,"p_wald","9.716047e-01"]])
  end

  # Test for https://github.com/genetics-statistics/GEMMA/issues/58
  # fixed GSLv2 NaN's that appeared with covariates.
  def test_plink_covariates_lmm
    gemma("-bfile example/HLC -gk 2 -o testPlinkStandardRelatednessMatrixK")

    gemma("-bfile example/HLC \
           -k output/testPlinkStandardRelatednessMatrixK.sXX.txt \
           -lmm 1 \
           -maf 0.1 \
           -c example/HLC_covariates.txt \
           -no-check \
           -o plink_lmm1_cov")
    expect("output/plink_lmm1_cov.assoc.txt",[[100,"p_wald","5.189953e-01"],
                                              [:max,"logl_H1","279.2689"],
                                              [:max,"l_remle","1.686062"],
                                              [:max,"p_wald","0.9999996"]])
  end


end
