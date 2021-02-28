module snl
   ! Special SNL parameters
   use locale ! provides kr
   implicit none

   ! global variables
   integer,public:: icntrl(40), imode(40), icode
   character(80),public::run_title='Default'
   character*80:: file_damage_name, ldir, optical
   real(kr):: break_new, al_new, zl_new, dummy, al_old, zl_old
   integer:: ipun, len_file, ndamage, nldir, noptical, nfile2
   integer:: idam_fnc
   real(kr):: epoint, ypoint, xthresh
   character(2):: media, izn
   character(1):: dot
   character(100):: file, filename
   real(kr):: displace_th, bgr_energy, bgr_let, damage_threshold
   real(kr):: a_crit, z_crit, value
   real(kr):: energy_let(175), xinc_let(175,25) 
   real(kr):: upper_let(175,25)
   real(kr):: energy_dam(1000), value_dam(1000)
   real(kr):: pe_archive, pe_store, pe_trial
end module snl
