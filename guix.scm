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
  (gnu packages algebra)
  (gnu packages base)
  (gnu packages compression)
  (gnu packages bioinformatics)
  (gnu packages build-tools)
  (gnu packages check)
  (gnu packages curl)
  (gnu packages gdb)
  (gnu packages llvm)
  (gnu packages maths)
  (gnu packages ninja)
  (gnu packages parallel)
  (gnu packages perl)
  ;; (gnu packages perl6)
  (gnu packages ruby)
  (gnu packages pkg-config)
  (pjotr packages openblas) ;; we use this for the static builds
  ;; (gnu packages shell)  ;; for shunit2
  (srfi srfi-1)
  (ice-9 popen)
  (ice-9 rdelim))

(define %source-dir (dirname (current-filename)))

(define %git-commit
    (read-string (open-pipe "git show HEAD | head -1 | cut -d ' ' -f 2" OPEN_READ)))

(define %gemma-version
    (read-string (open-pipe "cat VERSION" OPEN_READ)))

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
       ("gsl-static" ,gsl-static)
       ("gsl" ,gsl)
       ("openblas" ,openblas)
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
