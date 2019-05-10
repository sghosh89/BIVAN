# The importance of a complete statistical description of dependence between variables for ecology and related fields : Introduction to the Repository of All Analyses Supporting the Paper

Shyamolina Ghosh, University of Kansas 

Lawrence W. Sheppard, University of Kansas

Mark T. Holder, University of Kansas

Terrance E. Loecke, University of Kansas

Philip C. Reid, University of Plymouth

James D. Bever, University of Kansas

Daniel C. Reuman, University of Kansas

# Introduction

This repository can be used to reproduce the complete analyses behind the paper 
"The importance of a complete statistical description of dependence between variables for ecology and related fields" and 
to recompile the paper itself along with Supporting Information. Some data are also included in the repository.

# How to compile
Knit the makefile.Rmd using R markdown. If all dependencies are in place (see next section) this should re-compute all analyses from data to paper, resulting in three pdfs: **Paper.pdf** (the main text of the paper), **SupportingInformation.pdf** (the supporting information file for the paper), and **makefile.pdf** (notes on the compilation process - can be useful for error mitigation in the event of failure).

The knit may take a several hours or a few days, depending on your computer speed, number of cores used for parallel computing (we used all cores except 2, you can change this number in r chunk named *setup* of SupportingInformation.Rmd), the value of nsurrogs in SupportingInformation.Rmd, and other factors. Subsequent knits, if any, can be faster because packages will be installed (see below) and because intermediate results are cached.

If you try to knit Paper.Rmd or SupportingInformation.Rmd directly, you may have some success, but cross-document references and other features may fail so this is not recommended.

To compile the documents from the command line, use the following: Rscript -e "library(knitr); knit('makefile.Rmd')".

# Dependencies

## Core dependencies
   - R 
   - R markdown
   - R studio
   - latex 
   - bibtex

## Dependencies on the R checkpoint package

Code uses the R *checkpoint* package. This is set up in the master file makefile.Rmd in r chunk named *checkpoint_chunk*, which contains the following line of code specifying a date :

checkpoint("2018-03-29",checkpointLocation = "./")

The checkpoint package then automatically scans through other files looking for other required R packages. It then downloads and installs the newest versions of those packages available on the given date. This helps ensure that re-compiling the document uses exactly the same code that was originally used. This can take some time on first run (you are warned) but it is faster on subsequent runs because the packages are already installed. This also means that R package dependencies should only be the checkpoint package, since that package should scan for other packages and install them locally. Quite a few MB disk space are used (approx. 300).

## Dependencies on pandoc
The open source program pandoc converts documents from one format to another. 
Here, the knitr package uses it to convert the markdown files into latex format so that 
they can then be turned into PDF files. Installers for multiple operating systems are available here: https://pandoc.org/installing.html.

## Dependencies on pdflatex
The makefile makes a system call to pdflatex, so software supporting that needs to be installed:

 - On Windows, you can use Miktex (https://miktex.org/howto/install-miktex),
 - On Linux, install latex (e.g., sudo apt-get install texlive), and
 - On Mac, use the MacTeX installer (http://www.tug.org/mactex/)

## Additional dependencies?
If you find additional dependencies were needed on your system, please let us know: reuman@ku.edu. The compilation process was tested by Reuman on........... using R version ......... and R studio version ......... ; and by Ghosh on Ubuntu 16.04.5 LTS using R version 3.4.4 and R studio version 1.1.463. It has not been tested on Mac. We have endeavored to list all dependencies we can think of above, but we have only compiled on our own machines, so we cannot guarantee that additional dependencies will not also be needed on other machines. This repository is intended to record a workflow, and is not designed or tested for distribution and wide use on multiple machines. It is not guaranteed to work on the first try without any hand-holding on arbitrary computing setups.

# Intermediate files:
Knitting the makefile automatically produces a lot of 'intermediate' files. Files ending in .tex are the converted documents from .Rmd including all the R code output and the rest (files ending .log, .aux, .lof, .lot, .toc and .out ) are intermediate files that pdflatex uses to keep track of various parts of the document. Some of these can be useful for diagnosing problems, if any.

# Acknowlegements :

We thank the many contributors to the large datasets we used; D. Stevens and P. Verrier for data extraction; and Joel E. Cohen, Jonathan Walter, Thomas Anderson, and Lei Zhao for helpful suggestions. We thank James Bell of the Rothamsted Insect Survey (RIS). The Rothamsted Insect Survey, a National Capability, is funded by the Biotechnology and Biological Sciences Research Council under the Core Capability Grant BBS/E/C/000J0200. SG, LWS and DCR were partly funded by US National Science Foundation grants 0075516 and 1442595 and the James S McDonnell Foundation. 













