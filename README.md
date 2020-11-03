# SNL-NJOY-2016
Nuclear data processing with enhanced version of the legacy NJOY-2016 code.

The NJOY-2016 code is an Open Source code, developed by LANL, and hosted on GitHub code repository. 
Sandia has made changes to the baseline NJOY-2016 code in order to support applications by the 
radiation damage community. The application area that our enhancements address are in the development 
of the energy-dependent radiation response functions for describing radiation damage to materials. 

A typical application is iron embrittlement of the critical weld in the pressure vessel for a light 
water reactor. This energy dependent response permits the research community to correlate observed 
damage between exposures to the neutron environments in different reactors (typically comparison 
between commercial reactors and material test reactors). Thus, this enhanced code capability is 
useful to the academic community and the ASTM reactor pressure vessel embrittlement safety community. 

We tried to incorporate our code changes into the current LANL Open Source version of NJOY-2016, but the 
NJOY code system has moved on to the NJOY21 code system and there were some compatibility issues with 
incorporating our required capability enhancements into the NJOY21 code.  LANL has frozen their work on NJOY-2016. 
To distinguish our modified version of NJOY from the current LANL GitHub version, we are calling our code 
SNL-NJOY-2016 - and have also made our version available to the general public via GitHub. 
Our version includes our enhancements to the baseline LANL NJOY-2016 code using the Git software configuration 
control system - so configuration control between our version and LANL's version is rigorously maintained.  
LANL (Jeremy Conlin) will examine these enhancements and try to find ways to incorporate our changes 
into a future LANL version of the NJOY21 code. 

"Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of 
Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software."


                                                 NOTICE:

For five (5) years from 6/11/2020 the United States Government is granted for itself and others acting on its 
behalf a paid-up, nonexclusive, irrevocable worldwide license in this data to reproduce, prepare derivative works, 
and perform publicly and display publicly, by or on behalf of the Government. There is provision for the possible 
extension of the term of this license. Subsequent to that period or any extension granted, the United States Government 
is granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable worldwide license in this 
data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, 
and to permit others to do so. The specific term of the license can be identified by inquiry made to National Technology 
and Engineering Solutions of Sandia, LLC or DOE.
 
NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR NATIONAL TECHNOLOGY AND ENGINEERING 
SOLUTIONS OF SANDIA, LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL 
RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, 
OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
 
Any licensee of this software has the obligation and responsibility to abide by the applicable export control laws, 
regulations, and general prohibitions relating to the export of technical data. Failure to obtain an export control 
license or other authority from the Government may result in criminal liability under U.S. laws.

This software has been assigned SCR# 2499.0. This number is an internal tracking number and is a useful reference 
when contacting the Legal Technology Transfer Center or Licensing and IP Management. 

This software was originally assigned Export Control Classification Number (ECCN) EAR 99. However, since this software 
is to be released as OSS, the software is deemed to be Publicly Available. 

