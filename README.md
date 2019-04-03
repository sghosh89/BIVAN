# The importance of a complete statistical description of dependence between variables for ecology and related fields : 
Introduction to the Repository of All Analyses Supporting the Paper

Shyamolina Ghosh
Lawrence W. Sheppard
Mark T. Holder
Terrance E. Loecke
Philip C. Reid
Daniel C. Reuman 

# Introduction

This repository can be used to reproduce the complete analyses behind the paper 
"The importance of a complete statistical description of dependence between variables for ecology and related fields" and 
to recompile the paper itself along with Supporting Information. Some data are also included in the repository.

# How to compile
Knit the makefile.Rmd using R markdown. If all dependencies are in place (see next section) this should re-compute all analyses from data to paper, resulting in three pdfs: Paper.pdf (the main text of the paper), SupportingInformation.pdf (the supporting information file for the paper), and makefile.pdf (notes on the compilation process - can be useful for error mitigation in the event of failure).

The knit may take a several hours or a few days, depending on your computer speed, number of cores used for parallell computing, the value of nsurrogs in SupportingInformation.Rmd, and other factors. Subsequent knits, if any, can be faster because packages will be installed (see below) and because intermediate results are cached.

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

Code uses the R 'checkpoint' package. This is set up in the master file makefile.Rmd in r chunk named checkpoint_chunk, which contains the following line of code specifying a date :

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



















