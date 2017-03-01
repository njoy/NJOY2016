module physics
   ! Provides values of pi (taken from http://mathworld.wolfram.com/Pi.html)
   ! and physics contants taken from CODATA 2014 
   ! (taken from http://physics.nist.gov/cuu/Constants/). Both sites visted on
   ! 3/14/2017 (Pi Day).
   ! The contstants (units) used are:
   !  - bk: Boltzmann's constant (ev/K)
   !  - amassn: Neutron mass (amu)
   !  - amu: atomic mass unit (g)
   !  - hbar: Planck's constant (erg s) (units in cgs)
   !  - ev: Conversion from erg to eV
   !  - clight: Speed of light (cm/s)
   use locale ! provides kr
   implicit none
   real(kr),parameter,public::pi=3.14159265358979323846264338327950e0_kr
   real(kr),parameter,public::bk=8.6173303e-5_kr
   real(kr),parameter,public::amassn=1.00866491588_kr
   real(kr),parameter,public::amu=1.660539040e-24_kr
   real(kr),parameter,public::hbar=1.054571800e-27_kr
   real(kr),parameter,public::ev=1.6021766208e-19_kr
   real(kr),parameter,public::clight=2.99792458e10_kr
end module physics

