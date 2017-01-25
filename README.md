# font-lock-regression-suite - Test suite for font-lock

*Author:* Anders Lindgren<br>
*URL:* [https://github.com/Lindydancer/font-lock-regression-suite](https://github.com/Lindydancer/font-lock-regression-suite)<br>

A collection of example source files for a large number of
programming languages, with ERT tests to ensure that syntax
highlighting does not accidentally change.

For each source file, font-lock reference files are provided for
various Emacs versions.  The reference files contains a plain-text
representation of source file with syntax highlighting, using the
format "faceup".

Of course, the collection source file can be used for other kinds
of testing, not limited to font-lock regression testing.

## Copyright note

The Copyright at the beginning of this file applies to the files
that drive the regression suite. It does not apply to the source
examples.  See the individual source files for information
regarding copyright and licensing terms.

## Usage

Run <kbd>M-x font-lock-regression-suite-add-testcases RET</kbd>. This will
add a number of ERT test cases to verify that source files are
highlighted according to the reference files.

Run, for example, <kbd>M-x ert RET t RET</kbd> to run all tests.

You can bind `font-lock-regression-suite-reference-version` to
another Emacs version, to see what the changes are compared to that
version.

Reference files for several major Emacs versions are provided.
You can compare the files to see how syntax highlighting has
evolved over the years.  If you find the "faceup" format hard to
read, you can run <kbd>M-x faceup-render-view-buffer RET</kbd> to see how
Emacs used to highlight the buffer (given that all relevant faces
are defined).

## See also

- [Comments on syntax highlighting support provided by various
  major modes](doc/CommentsOnMajorModes.org)
- [The origin of the packages used as test
  examples](doc/PackageSources.org)

## Using the source files in other contexts

The function `font-lock-regression-suite-each-src-ref-file` can be
used to traverse all the files in the suite. It will accept one
argument, a function that will be called with four arguments: A
name, the source file name, the reference file name, and a mode.

Today the mode is a single symbol. However, to be future
compatible, this can be a list of symbols, which should be called
in order. (Think of this as a major mode and a number of minor
modes.)

### Example

The following piece of code will traverse all source file and echo
the source names:

        (font-lock-regression-suite-each-src-ref-file
         (lambda (name src-file ref-file mode)
           (message src-file)))

### Real-world examples

This package is used to test the packages `font-lock-profiler` and
`font-lock-studio`, to ensure that they behaves like the normal
font-lock engine, for non-trivial examples.

## Other Font Lock Tools

This package is part of a suite of font-lock tools.  The other
tools in the suite are:

### [Font Lock Studio](https://github.com/Lindydancer/font-lock-studio)

Interactive debugger for font-lock keywords (Emacs syntax
highlighting rules).

Font Lock Studio lets you *single-step* Font Lock keywords --
matchers, highlights, and anchored rules, so that you can see what
happens when a buffer is fontified. You can set *breakpoints* on or
inside rules and *run* until one has been hit. When inside a rule,
matches are *visualized* using a palette of background colors. The
*explainer* can describe a rule in plain-text English. Tight
integration with *Edebug* allows you to step into Lisp expressions
that are part of the Font Lock keywords.

### [Font Lock Profiler](https://github.com/Lindydancer/font-lock-profiler)

A profiler for font-lock keywords.  This package measures time and
counts the number of times each part of a font-lock keyword is
used.  For matchers, it counts the total number and the number of
successful matches.

The result is presented in table that can be sorted by count or
time.  The table can be expanded to include each part of the
font-lock keyword.

In addition, this package can generate a log of all font-lock
events.  This can be used to verify font-lock implementations,
concretely, this is used for back-to-back tests of the real
font-lock engine and Font Lock Studio, an interactive debugger for
font-lock keywords.

### [Highlight Refontification](https://github.com/Lindydancer/highlight-refontification)

Minor mode that visualizes how font-lock refontifies a buffer.
This is useful when developing or debugging font-lock keywords,
especially for keywords that span multiple lines.

The background of the buffer is painted in a rainbow of colors,
where each band in the rainbow represent a region of the buffer
that has been refontified.  When the buffer is modified, the
rainbow is updated.

### [Faceup](https://github.com/Lindydancer/faceup)

Emacs is capable of highlighting buffers based on language-specific
`font-lock` rules. This package makes it possible to perform
regression test for packages that provide font-lock rules.

The underlying idea is to convert text with highlights ("faces")
into a plain text representation using the Faceup markup
language. This language is semi-human readable, for example:

    «k:this» is a keyword

By comparing the current highlight with a highlight performed with
stable versions of a package, it's possible to automatically find
problems that otherwise would have been hard to spot.

This package is designed to be used in conjunction with Ert, the
standard Emacs regression test system.

The Faceup markup language is a generic markup language, regression
testing is merely one way to use it.


---
Converted from `font-lock-regression-suite.el` by [*el2markdown*](https://github.com/Lindydancer/el2markdown).
