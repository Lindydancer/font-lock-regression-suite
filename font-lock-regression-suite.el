;;; font-lock-regression-suite.el --- Test suite for font-lock.

;; Copyright (C) 2015-2017  Anders Lindgren

;; Author: Anders Lindgren
;; Keywords: faces
;; Version: 0.0.1
;; URL: https://github.com/Lindydancer/font-lock-regression-suite
;; Package-Requires: ((faceup "0.0.4"))

;; This program is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.

;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with this program.  If not, see <http://www.gnu.org/licenses/>.

;;; Commentary:

;; A collection of example source files for a large number of
;; programming languages, with ERT tests to ensure that syntax
;; highlighting does not accidentally change.
;;
;; For each source file, font-lock reference files are provided for
;; various Emacs versions.  The reference files contains a plain-text
;; representation of source file with syntax highlighting, using the
;; format "faceup".
;;
;; Of course, the collection source file can be used for other kinds
;; of testing, not limited to font-lock regression testing.

;; Copyright note:
;;
;; The Copyright at the beginning of this file applies to the files
;; that drive the regression suite. It does not apply to the source
;; examples.  See the individual source files for information
;; regarding copyright and licensing terms.

;; Usage:
;;
;; Run `M-x font-lock-regression-suite-add-testcases RET'. This will
;; add a number of ERT test cases to verify that source files are
;; highlighted according to the reference files.
;;
;; Run, for example, `M-x ert RET t RET' to run all tests.
;;
;; You can bind `font-lock-regression-suite-reference-version' to
;; another Emacs version, to see what the changes are compared to that
;; version.
;;
;; Reference files for several major Emacs versions are provided.
;; You can compare the files to see how syntax highlighting has
;; evolved over the years.  If you find the "faceup" format hard to
;; read, you can run `M-x faceup-render-view-buffer RET' to see how
;; Emacs used to highlight the buffer (given that all relevant faces
;; are defined).

;; See also:
;;
;; - [Comments on syntax highlighting support provided by various
;;   major modes](doc/CommentsOnMajorModes.org)
;;
;; - [The origin of the packages used as test
;;   examples](doc/PackageSources.org)

;; Using the source files in other contexts:
;;
;; The function `font-lock-regression-suite-each-src-ref-file' can be
;; used to traverse all the files in the suite. It will accept one
;; argument, a function that will be called with four arguments: A
;; name, the source file name, the reference file name, and a mode.
;;
;; Today the mode is a single symbol. However, to be future
;; compatible, this can be a list of symbols, which should be called
;; in order. (Think of this as a major mode and a number of minor
;; modes.)
;;
;; Example:
;;
;; The following piece of code will traverse all source file and echo
;; the source names:
;;
;;     (font-lock-regression-suite-each-src-ref-file
;;      (lambda (name src-file ref-file mode)
;;        (message src-file)))
;;
;; Real-world examples:
;;
;; This package is used to test the packages `font-lock-profiler' and
;; `font-lock-studio', to ensure that they behaves like the normal
;; font-lock engine, for non-trivial examples.

;; Other Font Lock Tools:
;;
;; This package is part of a suite of font-lock tools.  The other
;; tools in the suite are:
;;
;;
;; Font Lock Studio:
;;
;; Interactive debugger for font-lock keywords (Emacs syntax
;; highlighting rules).
;;
;; Font Lock Studio lets you *single-step* Font Lock keywords --
;; matchers, highlights, and anchored rules, so that you can see what
;; happens when a buffer is fontified. You can set *breakpoints* on or
;; inside rules and *run* until one has been hit. When inside a rule,
;; matches are *visualized* using a palette of background colors. The
;; *explainer* can describe a rule in plain-text English. Tight
;; integration with *Edebug* allows you to step into Lisp expressions
;; that are part of the Font Lock keywords.
;;
;;
;; Font Lock Profiler:
;;
;; A profiler for font-lock keywords.  This package measures time and
;; counts the number of times each part of a font-lock keyword is
;; used.  For matchers, it counts the total number and the number of
;; successful matches.
;;
;; The result is presented in table that can be sorted by count or
;; time.  The table can be expanded to include each part of the
;; font-lock keyword.
;;
;; In addition, this package can generate a log of all font-lock
;; events.  This can be used to verify font-lock implementations,
;; concretely, this is used for back-to-back tests of the real
;; font-lock engine and Font Lock Studio, an interactive debugger for
;; font-lock keywords.
;;
;;
;; Highlight Refontification:
;;
;; Minor mode that visualizes how font-lock refontifies a buffer.
;; This is useful when developing or debugging font-lock keywords,
;; especially for keywords that span multiple lines.
;;
;; The background of the buffer is painted in a rainbow of colors,
;; where each band in the rainbow represent a region of the buffer
;; that has been refontified.  When the buffer is modified, the
;; rainbow is updated.
;;
;;
;; Faceup:
;;
;; Emacs is capable of highlighting buffers based on language-specific
;; `font-lock' rules. This package makes it possible to perform
;; regression test for packages that provide font-lock rules.
;;
;; The underlying idea is to convert text with highlights ("faces")
;; into a plain text representation using the Faceup markup
;; language. This language is semi-human readable, for example:
;;
;;     «k:this» is a keyword
;;
;; By comparing the current highlight with a highlight performed with
;; stable versions of a package, it's possible to automatically find
;; problems that otherwise would have been hard to spot.
;;
;; This package is designed to be used in conjunction with Ert, the
;; standard Emacs regression test system.
;;
;; The Faceup markup language is a generic markup language, regression
;; testing is merely one way to use it.

;;; Code:

;; TODO:
;;
;; * Rake script for regenerating all faceup files.
;;
;; * List orphaned reference files.

(require 'faceup)

(defvar font-lock-regression-suite-languages
  '(("ada"         ada-mode)
    ("autoconf"    autoconf-mode)
    ("bat"         bat-mode)
    ("C"           c-mode ("CWarn" c-mode cwarn-mode))
    ("elisp"       emacs-lisp-mode)
    ("f90"         f90-mode)
    ("JavaScript"  js-mode)
    ("m4"          m4-mode ("ac" autoconf-mode))
    ("make"        makefile-mode)
    ("Objective-C" objc-mode)
    ("Perl"        perl-mode ("CPerl" cperl-mode))
    ("PostScript"  ps-mode)
    ("prolog"      prolog-mode)
    ("Python"      python-mode)
    ("Ruby"        ruby-mode))
  "List of directories and corresponding modes.

Each entry in the list has the following format:

    (DIR MODE-OR-MODES ...)

Where DIR is a directory in the source tree and MODE-OR-MODES is
a mode or a list of modes that should be used.")

(defvar font-lock-regression-suite-dir (faceup-this-file-directory))


(defvar font-lock-regression-suite-reference-version
  (if (string-match "\\([0-9]+\\.[0-9]+\\.[0-9]+\\)\\.[0-9]+"
                    emacs-version)
      (concat (match-string 1 emacs-version) ".x")
    emacs-version)
  "The version of the revision files to use, defaults to Emacs version.

For released Emacs versions, the same as `emacs-versions', else simplified.

When building an Emacs from source, the fourth version number is
increased for every build. This normalized this by replacing the
fourth version number with an `x'.")


;; -------------------------------------------------------------------
;; Support functions.
;;

(defun font-lock-regression-suite-each-src-ref-file--internal
    (src-dir ref-dir ref-start path func &rest args)
  "Helper function for `font-lock-regression-suite-each-src-ref-file'.

Traverse SRC-DIR recursively, relative to PATH and call FUNC with
name, source file, reference files (rooted in REF-DIR), and
ARGS."
  (dolist (f (directory-files (if path
                                  (concat src-dir path)
                                src-dir)))
    (unless (string-match "^\\." f)     ; ".", "..", ".nosearch" etc.
      (let ((new-path (if path
                          (concat path "/" f)
                        f)))
        (if (file-directory-p (concat src-dir new-path))
            (apply 'font-lock-regression-suite-each-src-ref-file--internal
                   src-dir ref-dir ref-start new-path func args)
          (apply func
                 ;; Note: In theory, the same name could be generated
                 ;; twice. For example, if both "A-B" and "A/B"
                 ;; exists, they will both be mapped to "A-B". One way
                 ;; to work-around this would be to map "-" into "--".
                 (font-lock-regression-suite-dashify
                  (concat ref-start "/" new-path))
                 (concat src-dir new-path)
                 (concat ref-dir ref-start "/" new-path ".faceup")
                 args))))))


(defun font-lock-regression-suite-each-src-ref-file (func &rest args)
  "Call FUNC with each name, source file, reference file, mode, and ARGS.

The name is a unique identifier representing the file.  The
reference file may not exist.

Modes is a function to call or a list of function to call.  You
can use `font-lock-regression-suite-apply-modes' to enable the
modes.

`font-lock-regression-suite-dir' contains the root of the source
files and `font-lock-regression-suite-languages' contains a list
of subdirectories and corresponding modes.

When non-nil `font-lock-regression-suite-reference-version', the
reference files of that version of Emacs is used.  When nil, the
reference files of the current Emacs version is used.

Example:

    (font-lock-regression-suite-each-src-ref-file
     (lambda (name src-file ref-file modes)
       (message \"%s: %s %s in %s\" name src-file ref-file mode)))"
  (let ((seen-ids '()))
    (dolist (entry font-lock-regression-suite-languages)
      (dolist (mode-or-modes (cdr entry))
        (let* ((ref-dir (font-lock-regression-suite-reference-directory))
               (src-start (nth 0 entry))
               (ref-start-base (if (stringp (car-safe mode-or-modes))
                                   (pop mode-or-modes)
                                 (nth 0 entry)))
               (ref-start ref-start-base)
               (count 2))
          ;; Make the reference directory unique.
          (while (member ref-start seen-ids)
            (setq ref-start (format "%s-%d" ref-start-base count))
            (setq count (+ count 1)))
          (push ref-start seen-ids)
          (apply
           #'font-lock-regression-suite-each-src-ref-file--internal
           (concat font-lock-regression-suite-dir
                   "src/"
                   src-start "/")
           ref-dir
           ref-start
           nil
           func
           mode-or-modes
           args))))))


(defun font-lock-regression-suite-dashify (path)
  "Convert PATH to something suitable to be part of an elisp identifier."
  (setq path (file-name-sans-extension path))
  (while (string-match "/" path)
    (setq path (replace-match "-" nil nil path)))
  path)


(defun font-lock-regression-suite-reference-directory ()
  "The root of the reference directory, with a trailing slash."
  (concat font-lock-regression-suite-dir
          "ref/"
          font-lock-regression-suite-reference-version "/"))


(defun font-lock-regression-suite-apply-modes (modes)
  "Apply all modes in MODES.

Modes can be a function to call or a list of functions.

Return nil if any of the function isn't defined, non-nil otherwise."
  ;; Note: Both are needed to recognize lambda expressions and symbols
  ;; referring to undefined functions.
  (when (or (symbolp modes)
            (functionp modes))
    (setq modes (list modes)))
  (let ((res t))
    (while (and res
                modes)
      (let ((m (pop modes)))
        (if (functionp m)
            (funcall m)
          (setq res nil))))
    res))


;; ----------------------------------------------------------------------
;; Add ERT test cases.
;;
;; This generates one ERT test case for each source file. This allows
;; you to use the ERT selection mechanism to test a subset of the files.

(defun font-lock-regression-suite-add-testcases ()
  (interactive)
  (font-lock-regression-suite-each-src-ref-file
   (lambda (name
            src-file
            ref-file
            mode)
     (eval `(ert-deftest
                ,(intern (concat "font-lock-regression-suite--" name))
                ()
              (if (file-exists-p ,ref-file)
                  (should (faceup-test-font-lock-file
                           (quote ,mode)
                           ,src-file
                           ,ref-file))
                (error "The reference file `%s' does not exist."
                       ,ref-file)))))))


;; ----------------------------------------------------------------------
;; Regenerate
;;

(defun font-lock-regression-suite-regenerate (&optional force)
  "Regenerate all reference files.

When C-u prefix, or when FORCE is non-nil, only regenerate missing files."
  (interactive (list current-prefix-arg))
  (font-lock-regression-suite-each-src-ref-file
   (lambda (name
            src-file
            ref-file
            modes)
     (when (or force
               (not (file-exists-p ref-file)))
       (with-temp-buffer
         (insert-file-contents src-file)
         (when (font-lock-regression-suite-apply-modes modes)
           ;; Don't generate a reference file when the font-lock
           ;; keywords have triggered an error. (For example,
           ;; prolog-mode on Emacs 23.3 throws a "No match 3 in
           ;; highlight" error.)
           (when (condition-case nil
                     (progn
                       (font-lock-fontify-region (point-min) (point-max))
                       t)
                   (error nil))
             (make-directory (file-name-directory ref-file) t)
             (faceup-write-file ref-file))))))))


;; ----------------------------------------------------------------------
;; Testing
;;

(defun font-lock-regression-suite-list ()
  "Echo all source files in the regression suite."
  (interactive)
  (with-output-to-temp-buffer "*FontLockRegressionSuite*"
    (font-lock-regression-suite-each-src-ref-file
     (lambda (name src-file ref-file mode)
       (princ (format "%s:\n  %s\n  %s\n  %s\n" name src-file ref-file mode))))
    (display-buffer (current-buffer))))


;; ----------------------------------------------------------------------
;; The End
;;

(provide 'font-lock-regression-suite)

;;; font-lock-regression-suite.el ends here
