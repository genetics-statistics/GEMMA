;; To use this file to build HEAD of gemma:
;;
;;   guix build -f guix.scm
;;
;; To get a development container (e.g., run in emacs shell).
;;
;;   guix shell -C -D -f guix.scm
;;

(use-modules
  ((guix licenses) #:prefix license:)
  (guix gexp)
  (guix packages)
  (guix git-download)
  (guix build-system gnu)
  (guix utils)
  (gnu packages algebra)
  (gnu packages base)
  (gnu packages compression)
  (gnu packages bioinformatics)
  (gnu packages build-tools)
  (gnu packages check)
  (gnu packages curl)
  (gnu packages gcc)
  (gnu packages gdb)
  (gnu packages llvm)
  (gnu packages maths)
  (gnu packages ninja)
  (gnu packages parallel)
  (gnu packages perl)
  ;; (gnu packages perl6)
  (gnu packages ruby)
  (gnu packages pkg-config)
  ;; (gnu packages shell)  ;; for shunit2
  (srfi srfi-1)
  (ice-9 popen)
  (ice-9 rdelim))

(define %source-dir (dirname (current-filename)))

(define %git-commit
    (read-string (open-pipe "git show HEAD | head -1 | cut -d ' ' -f 2" OPEN_READ)))

(define %gemma-version
    (read-string (open-pipe "cat VERSION" OPEN_READ)))

(define-public openblas-pangemma
  (package
    (name "openblas-pangemma")
    (version "0.3.21")
    (source
     (origin
       (method git-fetch)
       (uri (git-reference
             (url "https://github.com/xianyi/OpenBLAS")
             (commit (string-append "v" version))))
       (file-name (git-file-name name version))
       (sha256
        (base32
         "0yx1axiki12y0xz0d5s76vvl7ds36k0npv1sww08k2qslhz1g9qp"))))
    (build-system gnu-build-system)
    (properties `((tunable? . #t)))
    (arguments
     (list
      #:tests? #f ;; skip tests
      #:test-target "test"
      ;; No default baseline is supplied for powerpc-linux.
      #:substitutable? (not (target-ppc32?))
      #:make-flags
      #~(list (string-append "PREFIX=" #$output)
              (string-append "CFLAGS=-O3 -g -Wno-incompatible-pointer-types -Wno-error=implicit-function-declaration")
              "COPT="
              "COMMON_OPT="
              "DYNAMIC_ARCH="
              "SHELL=bash"
              "MAKE_NB_JOBS=0"          ;use jobserver for submakes

              ;; This is the maximum number of threads OpenBLAS will ever use (that
              ;; is, if $OPENBLAS_NUM_THREADS is greater than that, then NUM_THREADS
              ;; is used.)  If we don't set it, the makefile sets it to the number
              ;; of cores of the build machine, which is obviously wrong.
              "NUM_THREADS=128"

              ;; DYNAMIC_ARCH is only supported on some architectures.
              ;; DYNAMIC_ARCH combined with TARGET=GENERIC provides a library
              ;; which uses the optimizations for the detected CPU.  This can
              ;; be overridden at runtime with the environment variable
              ;; OPENBLAS_CORETYPE=<type>, where "type" is a supported CPU
              ;; type.  On other architectures we target only the baseline CPU
              ;; supported by Guix.
              #$@(cond
                    ((or (target-x86-64?)
                         (target-x86-32?)
                         (target-ppc64le?)
                         (target-aarch64?))
                     ;; Dynamic older enables a few extra CPU architectures
                     ;; on x86_64 that were released before 2010.
                     '("DYNAMIC_ARCH=1" "TARGET=GENERIC"))
                     ;; '("DYNAMIC_ARCH=" "TARGET_CORE=ZEN"))
                    ;; On some of these architectures the CPU type can't be detected.
                    ;; We list the oldest CPU core we want to have support for.
                    ;; On MIPS we force the "SICORTEX" TARGET, as for the other
                    ;; two available MIPS targets special extended instructions
                    ;; for Loongson cores are used.
                    ((target-mips64el?)
                     '("TARGET=SICORTEX"))
                    ((target-arm32?)
                     '("TARGET=ARMV7"))
                    ((target-riscv64?)
                     '("TARGET=RISCV64_GENERIC"))
                    (else '())))
      ;; no configure script
      #:phases
      #~(modify-phases %standard-phases
          (delete 'configure)
          (add-before 'build 'set-extralib
            (lambda* (#:key inputs #:allow-other-keys)
              ;; Get libgfortran found when building in utest.
              (setenv "FEXTRALIB"
                      (string-append
                       "-L"
                       (dirname
                        (search-input-file inputs "/lib/libgfortran.so")))))))))
    (inputs
     (list `(,gfortran "lib")))
    (native-inputs
     (list cunit gfortran perl))
    (home-page "https://www.openblas.net/")
    (synopsis "Optimized BLAS library based on GotoBLAS")
    (description
     "OpenBLAS is a BLAS library forked from the GotoBLAS2-1.13 BSD version.")
    (license license:bsd-3)))


(define-public gemma-git
  (package
    (name "gemma-git")
    (version (git-version "0.98.5" "HEAD" %git-commit))
    (source (local-file %source-dir #:recursive? #t))
    (build-system gnu-build-system)
    (inputs
     `(
       ("catch" ,catch2)
       ("gdb" ,gdb)
       ;; ("gsl-static" ,gsl-static)
       ("gsl" ,gsl)
       ("openblas" ,openblas-pangemma)
       ("ruby" ,ruby) ;; for testing
       ("zlib:static" ,zlib "static")
       ("zlib" ,zlib)))

    (native-inputs ; for running tests
     `(("perl" ,perl)
       ("which" ,which)
       ))
    (arguments
     `(#:phases
       (modify-phases %standard-phases
                      (delete 'configure)
                      (delete 'validate-runpath)
                      (add-before 'build 'bin-mkdir
                                  (lambda _
                                    (mkdir-p "bin")
                                    ))
                      (replace 'install
                               (lambda* (#:key outputs #:allow-other-keys)
                                 (let ((out (assoc-ref outputs "out")))
                                   (install-file "bin/gemma" (string-append out "/bin"))))))
       #:tests? #t
       #:parallel-tests? #f))
    (home-page "https://github.com/genetics-statistics")
    (synopsis "Tool for genome-wide efficient mixed model association")
    (description "Genome-wide Efficient Mixed Model Association (GEMMA)
provides a standard linear mixed model resolver with application in
genome-wide association studies (GWAS).")
    (license license:gpl3)))

gemma-git
