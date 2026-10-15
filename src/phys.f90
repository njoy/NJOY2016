module physics
   ! Provides pi and physics constants taken from CODATA'2014 as
   ! given on the NIST site: 
   !      https://physics.nist.gov/cuu/Constants
   ! namely, bk (Boltzmann's constant),
   ! amassn (the neutron mass in amu), amu (the amu value itself),
   ! hbar (Planck's constant), ev (the conversion to eV), and
   ! clight (the speed of light).
   use locale ! provides kr
   implicit none
!
! Constants per ENDF-102, Appendix H (February 1, 2018) edition, also
! identified as SVN Commit:  Revision 215).  This is the manual version
! available from http://www.nndc.bnl.gov/csewg/docs/endf-manual.pdf in
! the Summer 2018.
!  - numerical values are given for pi, Boltzmann's constant, eV
!    speed of light, atomic mass unit in eV, Planck's constant and
!    the fine structure constant.
!  - other values are given in terms of the above.
   real(kr),parameter,public::pi=3.141592653589793238e0_kr     !
   real(kr),parameter,public::euler=0.57721566490153286e0_kr !
   real(kr),parameter,public::bk=8.617333262e-5_kr         !eV/degK
   real(kr),parameter,public::ev=1.602176634e-12_kr     !erg/eV
   real(kr),parameter,public::clight=2.99792458e10_kr   !cm/s
   real(kr),parameter,public::amu=931.49410242e6_kr*ev/&
                                        (clight*clight) !g/amu
   real(kr),parameter,public::hbar=6.582119569e-16_kr*ev !Planck/2pi, erg
   real(kr),parameter,public::finstri=1.e16_kr*hbar/(ev*ev*clight) !inv fine str
! ****************************************************************
! * Light particle masses (in amu), per ENDF-102, Appendix H:    *
! * - note, these are particle masses, not atomic masses.        *
! * - we use "a" to start the variable name due to the legacy    *
! *   fortran naming convention that would have considered a     *
! *   name of the form "massn" to be an integer variable.        *
! ****************************************************************
   real(kr),parameter,public::amassn=1.00866491595e0_kr !neutron
   real(kr),parameter,public::amassp=1.007276466621e0_kr !proton
   real(kr),parameter,public::amassd=2.013553212745e0_kr !deuteron
   real(kr),parameter,public::amasst=3.01550071621e0_kr  !triton
   real(kr),parameter,public::amassh=3.014932247175e0_kr !helion (3)
   real(kr),parameter,public::amassa=4.001506179127e0_kr  !alpha
   real(kr),parameter,public::amasse=5.48579909065e-4_kr  !electron
   real(kr),parameter,public::pnratio=amassp/amassn ! proton/neutron mass
   real(kr),parameter,public::dnratio=amassd/amassn ! deuteron/neutron mass
   real(kr),parameter,public::tnratio=amasst/amassn ! triton/neutron mass
   real(kr),parameter,public::hnratio=amassh/amassn ! helion/neutron mass
   real(kr),parameter,public::anratio=amassa/amassn ! alpha/neutron mass

   real(kr),parameter,public::epair=amasse*amu*clight*clight/ev
end module physics

