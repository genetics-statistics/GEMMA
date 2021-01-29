;; To use this file to build HEAD of gemma:
;;
;;   guix build -f guix.scm
;;
;; To get a development container (e.g., run in emacs shell).
;;
;;   guix environment -C -l guix.scm

(use-modules
  ((guix licenses) #:prefix license:)
  (guix gexp)
  (guix packages)
  (guix git-download)
  (guix build-system meson)
  (gnu packages algebra)
  (gnu packages base)
  (gnu packages compression)
  (gnu packages bioinformatics)
  (gnu packages build-tools)
  (gnu packages curl)
  (gnu packages llvm)
  (gnu packages maths)
  (gnu packages ninja)
  (gnu packages parallel)
  (gnu packages perl)
  (gnu packages perl6)
  (gnu packages pkg-config)
  ;; (gnu packages shell)  ;; for shunit2
  (srfi srfi-1)
  (ice-9 popen)
  (ice-9 rdelim))

(define %source-dir (dirname (current-filename)))

(define %git-commit
    (read-string (open-pipe "git show HEAD | head -1 | cut -d ' ' -f 2" OPEN_READ)))

(define-public gemma-git
  (package
    (name "gemma-git")
    (version (git-version "0.98.4" "HEAD" %git-commit))
    (source (local-file %source-dir #:recursive? #t))
    (build-system meson-build-system)
    (inputs `(
              ("gsl" ,gsl)
              ;; ("shunit2" ,shunit2) ;; comes with gemma
              ("openblas" ,openblas)
              ("zlib" ,zlib)
              ))
    (native-inputs ; for running tests
     `(("perl" ,perl)
       ("which" ,which)
       ))
    (home-page "https://github.com/genetics-statistics")
    (synopsis "Tool for genome-wide efficient mixed model association")
    (description "Genome-wide Efficient Mixed Model Association (GEMMA)
provides a standard linear mixed model resolver with application in
genome-wide association studies (GWAS).")
    (license license:gpl3)))

gemma-git
