# Profiling

gperftools (formerly the Google profiler) is included in the .guix-dev
startup script. Compile gemma for profiling:

    make clean
    make profile

Run the profiler

     env CPUPROFILE=/tmp/prof.out ./bin/gemma -g ./example/mouse_hs1940.geno.txt.gz -p ./example/mouse_hs1940.pheno.txt     -gk -o mouse_hs1940
     pprof ./bin/gemma /tmp/prof.out

and `top` shows

```
Welcome to pprof!  For help, type 'help'.
(pprof) top
Total: 720 samples
     103  14.3%  14.3%      103  14.3% dgemm_kernel_ZEN
      39   5.4%  19.7%       79  11.0% ____strtod_l_internal
      37   5.1%  24.9%       53   7.4% __printf_fp_l
      36   5.0%  29.9%       36   5.0% __sched_yield
      34   4.7%  34.6%       34   4.7% __strlen_avx2
      31   4.3%  38.9%       31   4.3% __strspn_sse42
      26   3.6%  42.5%      116  16.1% ReadFile_geno
      25   3.5%  46.0%       26   3.6% _int_malloc
      23   3.2%  49.2%       23   3.2% gsl_vector_set
      18   2.5%  51.7%       18   2.5% __strcspn_sse42
```
