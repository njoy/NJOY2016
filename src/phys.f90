module physics
   ! Provides pi and physics constants taken from CODATA'89 as
   ! given on the NIST site; namely, bk (Boltzmann's constant),
   ! amassn (the neutron mass in amu), amu (the amu value itself),
   ! hbar (Planck's constant), ev (the conversion to eV), and
   ! clight (the speed of light).
   use locale ! provides kr
   implicit none
   real(kr),parameter,public::pi=3.14159265358979e0_kr
   real(kr),parameter,public::bk=8.617385e-5_kr
   real(kr),parameter,public::amassn=1.008664904e0_kr
   real(kr),parameter,public::amu=1.6605402e-24_kr
   real(kr),parameter,public::hbar=1.05457266e-27_kr
   real(kr),parameter,public::ev=1.60217733e-12_kr
   real(kr),parameter,public::clight=2.99792458e10_kr
end module physics

