# LDhelmet, version 1.10
-----
 Copyright (C) 2012  Andrew H. Chan, Paul A. Jenkins, Yun S. Song

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 <http://www.gnu.org/licenses/>

 email: P.Jenkins@warwick.ac.uk or yss@berkeley.edu

-----

## SUMMARY:

LDhelmet is a software program for statistical inference of fine-scale crossover recombination rates from population genetic data. LDhelmet expects the input data to be phased and aligned DNA sequences. The statistical model assumes the sample of DNA sequences is randomly drawn from a single population, and that the population followed a neutral evolution with constant population size. However, the statistical inference can still be accurate even with mild deviations from these assumptions. LDhelmet can handle sample sizes of up to 50 individuals (haplotypes), and is suitable for whole-genome sequence analysis.

LDhelmet is an implementation of the method described in

  * Andrew H. Chan, Paul A. Jenkins, Yun S. Song (2012).  Genome-Wide Fine-Scale Recombination Rate Variation in Drosophila melanogaster.  PLoS Genet 8(12): e1003090. https://doi.org/10.1371/journal.pgen.1003090

An additional module (convert_table) has been added (along with minor bug fixes) in 2016.  This module should be used in conjunction with lookup tables created by LDpop described in

  * Kamm, J.A., Spence, J.P., Chan, J., and Song, Y.S.  Two-Locus Likelihoods under Variable Population Size and Fine-Scale Recombination Rate Estimation.  Genetics, Vol. 203 No. 3 (2016) 1381-1399. http://dx.doi.org/10.1534/genetics.115.184820

We recommend that you also try using pyrho (https://github.com/popgenmethods/pyrho), a fast demography-aware inference of fine-scale recombination rates based on fused-LASSO, described in 

  * Spence, J.P. and Song, Y.S.  Inference and analysis of population-specific fine-scale recombination maps across 26 diverse human populations.   Science Advances, Vol. 5, No. 10, eaaw9206 (2019).  https://doi.org/10.1126/sciadv.aaw9206

## LICENSES:

The source code is released under the GNU General Public License, version 3. The full text of the license can be found in LICENSE_GPLv3.txt, which should have been included with this README.

The documentation is released under GNU Free Documentation License. The full text of the license can be found in LICENSE_GFDLv1.3.txt, which should have been included with this README.

## INSTALL:

The accompanying manual (ldhelmet_manual.pdf) provides detailed installation instructions for LDhelmet. Please refer to the manual for additional installation instructions.

LDhelmet depends on
a) the Boost C++ Libraries
b) the GNU Scientific Library.

These libraries must first be installed, and if necessary, the Makefile for LDhelmet must be modified to contain path information for the libraries (see comments in the Makefile).

Once the dependencies are installed, simply type

make

in the top-level directory of LDhelmet. The resulting binary is called ldhelmet and will be placed in the top-level directory of LDhelmet.

## USAGE:

Please refer to the accompanying manual (ldhelmet_manual.pdf) for usage instructions.

## CONTACT:

Please contact P.Jenkins@warwick.ac.uk or yss@berkeley.edu with bugs, comments, or questions regarding this version of the software.

## HISTORY:
* 1.10 Fixed a bug which led to a segmentation fault on very long MCMC runs.
* 1.9 Fixed a bug related to likelihood lookup table interpolation.
* 1.8 Added a module (convert_table) to convert LDHat-like tables to LDHelmet tables. Bug fix: Unfortunately, previous versions computed rjMCMC acceptance probabilities incorrectly. These have now been fixed. In previous versions the issue seems to have been ameliorated by using a large block penalty (around 50 in our analyses), so for this version onwards we recommend experimenting with lower block penalties.
Other minor bug fixes. Added likelihood lookup interpolation.
* 1.7: Corrected typos in manual.
* 1.6: Made minor changes to fix problems with some compilers.
* 1.5: Updated copyright notice in header.
* 1.4: Fixed bug in error checking in table_gen.
* 1.3: Fixed bug in post_to_text.
* 1.2: Fixed bug in SNP partitioning.
* 1.1: Fixed bug in max_lk command-line options.
* 1.0: Initial release.
