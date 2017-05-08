---
layout: page
title: NJOY2016 Test Descriptions
---
## Test Problem 1

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/01/input)]

 This one runs on carbon to check `RECONR` for materials with no resonance parameters, to do Doppler broadening to 300K, to calculate heating KERMA and radiation damage, and to calculate the thermal scattering from graphite.  The pointwise (PENDF) results are converted into multigroup form.  The important things to look at are the various counts in `RECONR` and the multigroup values from `GROUPR` on the [output file](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/01/referenceOutput).  To see even more detail, check the reference [PENDF file](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/01/referenceTape25).

## Test Problem 2

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/02/input)]

 This is a demonstration of preparing a library for LMFBR applications using methods active in the 80's.  It is still useful for checking `RECONR` resonance reconstruction, Doppler broadening to more than one temperature, the generation of unresolved resonance data, and multigroup averaging with self shielding.  The `CCCC` output provides an alternative look at the cross sections and group-to-group matrices.  Pay attention to the output from `RECONR` to see if the same number of resonance points are produced.  Note the number of points at each temperature and the thermal quantities printed by `BROADR`.  The unresolved self shielding shows up in the `UNRESR` output and again in `GROUPR` and `CCCCR`.  The multigroup constants printed on the [output file](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/02/referenceOutput) provide plenty of opportunities to check your installation.  For even more detail, check the reference [PENDF file](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/02/referenceTape28).

## Test Problem 3

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/03/input)]

 This one demonstrates processing photoatomic data and the use of MATXS output files.  Newer versions of the photoatomic data are all on one file, rather than the two shown here.  Note that two materials are processed to show the limiting behaviors.  The multigroup constants printed by `GAMINR` on the [output file](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/02/referenceOutput) can be checked.  Two different output formats are generated: the DTF format and the MATXS format.  The `DTFR` numbers are given on the listing.  `DTFR` also generates Postscript plots of its results; the [Postscript](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/03/referenceTape37) file has been converted to [PDF](https://raw.githubusercontent.com/njoy/NJOY2016/master/docs/tests/03/referenceTape37.pdf) format for easy viewing.

## Test Problem 4

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/04/input)]

 This test problem computes cross section covariances using the `ERRORR` module for <sup>235</sup>U from ENDF/B-V.  The performance of `ERRORR` can be tested by looking at the values on the [output file](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/04/referenceOutput).

## Test Problem 5

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/05/input)]

 This one demonstrates formatting covariance data using our "boxer" format, and it generates detailed Postscript plots of the variances and correlations.  The [Postscript](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/05/referenceTape35) file has been converted to [PDF](https://raw.githubusercontent.com/njoy/NJOY2016/master/docs/tests/05/referenceTape35.pdf) form for easy viewing.

## Test Problem 6

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/06/input)]

 This case demonstrates and tests some of the features of the `PLOTR` and `VIEWR` modules.  The [output file](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/06/referenceOutput) is not very interesting.  Look at [PDF](https://raw.githubusercontent.com/njoy/NJOY2016/master/docs/tests/06/referenceTape32.pdf) for the results.

## Test Problem 7

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/07/input)]

 This example demonstrates the production of a library for the MCNP continuous-energy Monte Carlo code for <sup>238</sup>Pu from ENDF/B-IV.  Older versions of MCNP required a set of multigroup photon production tables to construct the 30x20 grid of photon emission bins used in those days.  This problem demonstrates how that was done, but the method is rarely needed these days.  Still, these results provide another useful set of tests for your installation.  The most interesting numbers are those in the `ACER` [output file](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/07/referenceOutput).  We provide both the [PENDF file](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/07/referenceTape26) and the [ACE file](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/07/referenceTape26) to allow for very detailed checks.

## Test Problem 8

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/08/input)]

 This case was added to check the processing of a typical ENDF/B-VI material using Reich-Moore resonances and File 6 for energy-angle distributions.  Pay special attention to the resonance results as shown by the counts on the `RECONR` listing and the values on the [PENDF](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/08/referenceTape28) and [ACE](https://raw.githubusercontent.com/njoy/NJOY2016/master/docs/tests/08/referenceTape25) files.  In addition, check the energy-angle distributions as printed out by `ACER` on the [ACE file](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/08/referenceOutput).

## Test Problem 9

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/09/input)]

 This one demonstrates the use of `LEAPR` to generate a scattering kernel for water.  We have used the ENDF physics model, but we reduced the alpha and beta ranges to make the case run faster with less output.  Note that `RECONR` and `BROADR` are run to prepare a base for the thermal data at 296K.  `THERMR` was run with its long printout option to provide plenty of numbers on the [output file](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/07/referenceOutput) for comparisons.  For additional details, look at the actual `LEAPR` [output](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/09/referenceTape24).

## Test Problem 10

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/10/input)]

 This sample problem demonstrates the production of unresolved resonance probability tables for MCNP.  We run both `UNRESR` and `PURR` using the same sigma0 grid to allow additional checks of both modules and to allow for comparisons between the deterministic and probabilistic approaches to computing Bondarenko cross sections.  We provide the [output](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/10/referenceOutput) listing, the [PENDF](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/10/referenceTape28) file, and the [ACE](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/10/referenceTape26) file.  Lots of results to check.

## Test Problem 11

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/11/input)]

 This one demonstrates the production of a library for the WIMS reactor lattice code using <sup>238</sup>Pu from ENDF/B-IV.  Check both `GROUPR` output and `WIMSR` in the [output file](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/11/referenceOutput).  The [WIMS library file](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/11/referenceTape27) produced is also provided.

## Test Problem 12

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/12/input)]

  This one shows how the gas production capability works with <sup>61</sup>Ni from ENDF/B-VI.  To see the detailed results for gas production, look for MT=203 and 207 in File 3 of the [PENDF tape](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/12/referenceTape22).  Also, check out the [plots](https://raw.githubusercontent.com/njoy/NJOY2016/master/docs/tests/12/referenceTape24.pdf) of resonance cross sections and gas production cross sections.

## Test Problem 13

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/13/input)]

  This case demonstrates the modern MCNP formats using <sup>61</sup>Ni from ENDF/B-VI.  Note that `ACER` is run twice, once to prepare the ACE file, and again to do consistency checking and to prepare detailed [plots](https://raw.githubusercontent.com/njoy/NJOY2016/master/docs/tests/13/referenceTape36.pdf). 

## Test Problem 14

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/14/input)]

  This problem was developed to demonstrate ACE incident proton data and the MCNPX charged-particle format.  Once again, two `ACER` runs are used to do checks and prepare [plots](https://raw.githubusercontent.com/njoy/NJOY2016/master/docs/tests/14/referenceTape36.pdf).

## Test Problem 15

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/15/input)]

  This problem demonstrates advanced covariance processing using <sup>238</sup>U from JENDL-3.3.  Covariances from resonance parameters are included.  Multiple `ERRORR`, `COVR`, `GROUPR`, and `VIEWR` runs are stacked together in order to generate covariances for the multigroup nubar (see [reference plot 45]()), the multigroup cross sections (see [reference plot 46]()), and for the elastic angular P1 angular coefficient (see [reference plot 47]()).  The script is complicated, but it is nice to get all the results in one run.

## Test Problem 16

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/16/input)]

  This problem is similar to 15, except it omits the `GROUPR` module, demonstrating that covariances can be calculated directly from PENDF data.  It also demonstrates some advanced plotting techniques.

## Test Problem 17

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/17/input)]

  This problem uses <sup>235</sup>U, <sup>238</sup>U, and <sup>239</sup>Pu from JENDL-3.3 in order to demonstrate cross-material covariances.  This run is pretty slow.

## Test Problem 18

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/18/input)]

  This problem demonstrates covariance processing for a secondary energy distribution.  An artificial version of <sup>252</sup>Cf was prepared, and the spontaneous fission spectrum covariances were produced.

## Test Problem 19

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/19/input)]

  This case is for a more modern actinide evaluation, <sup>241</sup>Pu from ENDF/B-VI.

## Test Problem 20

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/20/input)]

  This test case uses a preliminary ENDF/B-VII.1 evaluation for <sup>35</sup>Cl, which uses the advanced Reich-Moore-Limited format for the resonances.  This allows more complex covariances to be represented.  In this case, covariances involving the (n,p) reaction are calculated.

## Test Problem 21

[[input](https://raw.githubusercontent.com/njoy/NJOY2016/master/test/21/input)]
 
  This test case was developed to check if NaNs were "calculated" in `PURR`. This problem was discovered by Dave Brown and Paul Romano and was fixed in [Pull Request 18](https://github.com/njoy/NJOY2016/pull/18). This test was created in an attempt to prevent this from happening again. (No guarantees.)

