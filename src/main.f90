program njoy
!-----------------------------------------------------------------------
!
!    NJOY Nuclear Data Processing System
!    Version 2016.20
!
!-----------------------------------------------------------------------
!
! NJOY is a system of processing modules intended to convert evaluated
! nuclear data in the ENDF format into forms useful for practical
! applications.  It consists of a number of independent processing
! modules that perform various processing tasks and this main program
! (njoy) that runs the modules in the desired order.  The modules
! currently available are the following:
!
! reconr...reconstruct pointwise cross sections from ENDF/B resonance
!          parameters and interpolation schemes.
!
! broadr...doppler broaden and thin pointwise cross sections.
!
! unresr...compute effective pointwise self-shielded cross sections in
!          the unresolved energy range.
!
! heatr....compute heat production cross sections (kerma) and damage
!          energy production.
!
! thermr...generate neutron scattering cross sections and
!          point-to-point scattering kernels in the thermal range for
!          free or bound atoms.
!
! groupr...generate self-shielded multigroup cross sections and
!          group-to-group scattering and photon production matrices.
!
! gaminr...compute multigroup photon interaction cross sections,
!          scattering matrices, and heat production.
!
! errorr...construct multigroup covariance matrices.
!
! covr.....process covariance data from errorr.
!
! moder....convert between ENDF/B standard coded mode and the njoy
!          blocked binary mode.
!
! dtfr.....output and plot multigroup data for discrete ordinates
!          transport codes.
!
! ccccr....format multigroup data into the CCCC standard interface
!          files ISOTXS, BRKOXS, and DLAYXS.
!
! matxsr...convert multigroup data into the comprehensive MATXS
!          cross section interface format.
!
! resxsr...prepare a CCCC-like file of pointwise resonance
!          cross sections for thermal flux calculations.
!
! acer.....prepare library for the Los Alamos continuous energy
!          Monte Carlo code MCNP.
!
! powr.....convert multigroup data into libraries for the thermal
!          power reactor codes EPRI-CELL and EPRI-CPM.
!
! wimsr....convert multigroup data into libraries for the reactor
!          assembly codes WIMS-D or WIMS-E.
!
! plotr....plot ENDF, PENDF, GENDF, or experimental cross sections,
!          distributions, or matrices.
!
! viewr....view plots from plotr, dtfr, covr, etc. in Postscript.
!
! mixr.....mix file 3 cross sections (for example, to make elemental
!          cross sections) for plotting, etc.
!
! purr.....generate unresolved-resonance probability tables for the
!          MCNP Monte Carlo code.
!
! leapr....generate S(alpha,beta) for thermal moderators.
!
! gaspr....add gas production (MT203-207) to PENDF tape.
!
! Each processing module is implemented as a Fortran-90 module with
! only one public subroutine named as above.  The modules have the
! names as given above with "m" instead of "r".  Processing modules
! communicate through I/O units, and their internal data and
! subroutines are kept private to reduce the possibility of side
! effects between modules.  There are also some modules providing
! common capabilities to all the processing modules:
!
!   module      description
!   ---------   --------------------------------------------
!   version     contains the current version information
!   locale      contains localization parameters
!   mainio      contains the main I/O unit numbers
!   util        contains utility subroutines used by modules
!   endf        contains routines for handling ENDF files
!   math        contains math routines used by other modules
!   physics     contains common physics constants
!   snl         contains parameters that support SNL-specific enhanced capabilities
!
!---input specifications (free format)----------------------------------
!
! card 1       module option
!
!    module    six character module name, e.g., reconr.
!              it is not necessary to use quotes.
!
!          repeat card 1 for each module desired, and
!          use the name "stop" to terminate the program.
!
! See the comments at the start of each module for its specific input
! instructions.
!
!-----------------------------------------------------------------------
!
! Add SNL-specific modifications
!
! Item 0: Title card
!         special flag for SNL new control integers flagged by the use of the "Enh:" in the first four characters
!
!         Description of optional control imode flags
!
!         ************************************************
!
!         SNL_enhanced_input_format = 0     Enhanced format off
!         SNL_enhanced_input_format = 1     Enhanced format o  
!
! Item 1: implement optional read of flow control flags based on title card format (imode)
!
!         Description of optional control imode flags
!
!         ************************************************
!
!         imode(1) =   0    not used
!         imode(2) =   0    not used
!         imode(3) =   0    default
!                     <0    enhanced levels of debug print-out
!
! Item 2: implement optional read of control flags based on title card format (icntrl)
!
!         Description of optional control icntrl flags
!
!         ************************************************
!
!         icntrl(1) =  0    output displacement energy (with lower imtegration bound of Ed) - NJOY default
!                               (redundant with icntrl(1)=8 option)
!                      1    output recoil energy instead of damage energy in df
!                      2    BGR recoil energy threshold option
!                      3    BGR LET threshold option
!                      4    output LET instead of damage energy 
!                      5    output NRT damage energy (break = 2*Ed/beta)
!                      6    output original 3-level Kinchin-Pease damage energy (break = 2*Ed, no beta=0.8 factor))
!                      7    output displacement kerma (i.e. Ed = 0 eV, no break)
!                      8    output sharp transition Kinchin-Pease damage energy - no transistion region (break = Ed)
!
!         icntrl(2) =  1    punch pka spectra at epoint (na)
!         icntrl(3) =  9    TENDL-2010 bypass File 32 format issue
!         icntrl(4) =  1    apply damage efficiency factor (to be read-in) to displacement partition
!         icntrl(5) =  1    modify displacement threshold and target mass/charge 
!                      2    read-in BGR recoil threshold energy
!                      3    read-in BGR LET threshold value
!         icntrl(6) =  0    for icntrl(1)=4, score LET for particles > d (not implemented)
!                      1                                   particles > a (not implemented)
!                      2                                   all particles (not implemented)
!         icntrl(7) =  1    list reaction source for damage increments
!         icntrl(8) =  0    use build-in Robinson function
!                      1    use tabular function (to be read-in) with variable threshold energy
!         icntrl(9)         Recoil particle-dependent damage partition function
!                   =  0    full Robinson treatment   
!                   =  i    ignore damage energy from charged particles
!                             with atomic mass "i" less than the lattice atom
!         icntrl(10)=  0    report damage energy - default
!                   =  1    report dpa - converted from damage energy
!
!         ***************************************************
!
!
!
   use version ! provides vers,vday
   use locale  ! provides lab,mx
   use mainio  ! provides nsysi,nsyso,nsyse
   use util    ! provides dater,wclock,timer,error
   use endf    ! provides npage,iverf

   use reconm  ! provides reconr
   use broadm  ! provides broadr
   use unresm  ! provides unresr
   use heatm   ! provides heatr
   use thermm  ! provides thermr
   use groupm  ! provides groupr
   use gaminm  ! provides gaminr
   use errorm  ! provides errorr
   use covm    ! provides covr
   use modem   ! provides moder
   use dtfm    ! provides dtfr
   use ccccm   ! provides ccccr
   use matxsm  ! provides matxsr
   use resxsm  ! provides resxsr
   use acem    ! provides acer
   use powm    ! provides powr
   use wimsm   ! provides wimsr
   use plotm   ! provides plotr
   use viewm   ! provides viewr
   use mixm    ! provides mixr
   use purm    ! provides purr
   use leapm   ! provides leapr
   use gaspm   ! provides gaspr
   
   use snl     ! provides icntrl

   implicit none
   character(6)::module
   character(8)::date
   character(8)::time
   character(28)::strng
   real(kr)::secs
   integer:: jk, module_loop_number

   !--set page size for blocked binary mode in module endf
   npage=(npage/102)*102

   !--open the NJOY listing file
   !--using nsyso from module mainio
   open(nsyso,file='output')

   !--write njoy banner on main output listing
   !--using information from modules locale and version.
   call dater(date)
   call wclock(time)
   write(nsyso,'(''1'',77(''*'')/ &
    &'' *'',41x,''*'',14x,''*'',18x,''*''/ &
    &'' *'',3x,''$$'',3x,''$$'',7x,''$$'',3x,''$$$$$'',3x, &
    &''$$'',4x,''$$'',3x,''*'',14x,''*'',2x, &
    &''vers: '',a8,2x,''*''/ &
    &'' *'',3x,''$$$'',2x,''$$'',7x,''$$'',2x,''$$$$$$$'',3x, &
    &''$$'',2x,''$$'',4x,''*'',2x,''nuclear'',5x,''*'',2x, &
    &''vday: '',a8,2x,''*''/ &
    &'' *'',3x,''$$$$ $$'',7x,''$$'',2x,''$$'',3x,''$$'',4x, &
    &''$$$$'',5x,''*'',2x,''data'',8x,''*'',2x, &
    &''site: '',a8,2x,''*''/ &
    &'' *'',3x,''$$ $$$$'',2x,''$$'',3x,''$$'',2x,''$$'',3x, &
    &''$$'',5x,''$$'',6x,''*'',2x,''processing'',2x,''*'',2x, &
    &''mach: '',a8,2x,''*''/ &
    &'' *'',3x,''$$'',2x,''$$$'',2x,''$$$$$$$'',2x,''$$$$$$$'', &
    &5x,''$$'',6x,''*'',2x,''system'',6x,''*'',2x, &
    &''date: '',a8,2x,''*''/ &
    &'' *'',3x,''$$'',3x,''$$'',3x,''$$$$$'',4x,''$$$$$'',6x, &
    &''$$'',6x,''*'',14x,''*'',2x, &
    &''time: '',a8,2x,''*''/ &
    &'' *'',41x,''*'',14x,''*'',18x,''*''/ &
    &1x,77(''*''))') vers,vday,lab,mx,date,time

   !--write tty banner on terminal output
   !--using nsyse from module mainio.
   write(nsyse,'(/'' njoy '',a8,1x,a8,38x,a8,1x,a8/ &
    &1x,77(''*''))') vers,vday,date,time

   !--loop over the requested modules
   !--using nsysi from module mainio.
   
   run_title = 'Def: Placeholder for the default NJOY-2016 title'
   do jk = 1, 40
     icntrl(jk) = 0
     imode(jk) = 0
   enddo
   
   read (nsysi,'(a)') run_title
   write (nsyse,'(/,''Run Title:: '', a80,/)') run_title
   write (nsyso,'(/,''Run Title:: '', a80,/)') run_title
   SNL_enhanced_input_format = 0
   if ( run_title(1:4) .eq. 'Enh:') then
   SNL_enhanced_input_format = 1
     read (nsysi,'(3i2)') (imode(jk), jk=1,3)
     read (nsysi,'(10i2)') (icntrl(jk), jk=1,10)
     write (nsyso,'(''Enhanced control logic inputs:'')')
     write (nsyso,'(/,''    imode:  '', 3i2)') (imode(jk), jk=1,3)
     write (nsyso,'(  ''    icntrl: '', 10i2,/)') (icntrl(jk), jk=1,10)
   endif

   if ( run_title(1:5) .eq. 'moder') then
     SNL_enhanced_input_format = -1
     write (nsyso,'(''Ols style control logic inputs:'')')
   endif

   do jk = 1, 2 
     if ( imode(jk) .ne. 0) then 
        write (nsyso,'(''imode control logic input error for flag: '', i3)') jk
     endif
   enddo
   if ( imode(3) .ne. 0) then 
        write (nsyso,'(/,''imode(3) debug printer control logic flag set: '', i3,/)') imode(3)
   endif

   if ( icntrl(1) .eq. 0 .and. SNL_enhanced_input_format .eq. 1) then 
        write (nsyso,'(''icntrl(1) control logic flag set for NJOY default (spKP-DE): '', i3)') icntrl(1)
   elseif (icntrl(1) .eq. 1) then
        write (nsyso,'(''icntrl(1) control logic flag not implemented 1: '', i3)') icntrl(1)
   elseif (icntrl(1) .eq. 2) then
        write (nsyso,'(''icntrl(1) control logic flag not implemented 1: '', i3)') icntrl(1)
   elseif (icntrl(1) .eq. 3) then
        write (nsyso,'(''icntrl(1) control logic flag not implemented 1: '', i3)') icntrl(1)
   elseif (icntrl(1) .eq. 4) then
        write (nsyso,'(''icntrl(1) control logic flag not implemented 1: '', i3)') icntrl(1)
   elseif (icntrl(1) .eq. 5) then
        write (nsyso,'(''icntrl(1) control logic flag set for NRT-DE: '', i3)') icntrl(1)
   elseif (icntrl(1) .eq. 6) then
        write (nsyso,'(''icntrl(1) control logic flag set for origKP-DE: '', i3)') icntrl(1)
   elseif (icntrl(1) .eq. 7) then
        write (nsyso,'(''icntrl(1) control logic flag set for displacement kerma: '', i3)') icntrl(1)
   elseif (icntrl(1) .eq. 8) then
        write (nsyso,'(''icntrl(1) control logic flag set for spKP-DE: '', i3)') icntrl(1)
   else
        write (nsyso,'(''icntrl(1) control logic flag not implemented:    1 '', i5)') icntrl(1)
   endif
   do jk = 2, 3
     if ( icntrl(jk) .ne. 0) then 
        write (nsyso,'(''icntrl control logic flag not implemented yet: '', 2i5)') jk, icntrl(jk)
     endif
   enddo
   if ( icntrl(4) .eq. 1) then 
!     read-in recoil atom energy-dependent efficiency factor
        write (nsyso,'(''icntrl(4) control logic flag set to damage efficiency factor read-in: '', i5)') icntrl(4)
   elseif ( icntrl(4) .gt. 1) then
        write (nsyso,'(''icntrl(4) control logic flag not implemented:    4 '', i5)') icntrl(5)
   elseif ( icntrl(4) .lt. 0) then
        write (nsyso,'(''icntrl(4) control logic flag not implemented:    4 '', i5)') icntrl(5)
   endif
   if ( icntrl(5) .eq. 1) then 
!     read-in replacement lattice ion mass/charge
        write (nsyso,'(''icntrl(5) control logic flag set to lattice atom read-in: '', i5)') icntrl(5)
   elseif ( icntrl(5) .gt. 1) then
        write (nsyso,'(''icntrl(5) control logic flag (a) not implemented:    5 '', i5)') icntrl(5)
   elseif ( icntrl(5) .lt. 0) then
        write (nsyso,'(''icntrl(5) control logic flag (b) not implemented:    5 '', i5)') icntrl(5)
   endif

   do jk = 6, 7
     if ( icntrl(jk) .ne. 0) then 
        write (nsyso,'(''icntrl(6/7) control logic flag (c) not implemented: '', 2i5)') jk, icntrl(jk)
     endif
   enddo

   if ( icntrl(8) .eq. 1) then
      write (nsyso,'(''icntrl(8) control logic flag set for user-input of damage parition function: '', i3)') icntrl(8)
   elseif ( icntrl(8) .ne. 0 .and. icntrl(8) .ne. 1) then 
      write (nsyso,'(''icntrl(8) control logic flag (d) not implemented: '', 2i5)') jk, icntrl(jk)
   endif

   do jk = 9, 9
     if ( icntrl(jk) .ne. 0) then 
        write (nsyso,'(''icntrl(9) control logic flag for low mass residual particles: '', i5)') icntrl(jk)
     endif
   enddo
   if ( icntrl(10) .eq. 0 .and. SNL_enhanced_input_format .eq. 1) then 
        write (nsyso,'(''icntrl(10) control logic flag set for damage energy: '', i3)') icntrl(10)
   elseif (icntrl(10) .eq. 1) then
        write (nsyso,'(''icntrl(10) control logic flag set for dpa: '', i3)') icntrl(10)
   else
        write (nsyso,'(''icntrl(10) control logic flag not implemented:    10 '', i5)') icntrl(10)
   endif
   
   module_loop_number = 0
   do
      module_loop_number = module_loop_number + 1
      if ( module_loop_number .eq. 1 .and. SNL_enhanced_input_format .eq. -1) then 
         module = run_title(1:6)
      else
         read(nsysi,*) module
      endif
      if (module.eq.'stop') exit
      select case(module)

      case ('#')     ! pre-process instruction; nothing to do

      case ('@')     ! post-process instruction; nothing to do

      case('reconr') ! reconstruct pointwise cross-sections
         call reconr

      case('broadr') ! doppler broaden point cross sections
         call broadr

      case('unresr') ! compute unresolved resonance cross-sections
         call unresr

      case('heatr')  ! compute heating and damage
         call heatr

      case('thermr') ! compute thermal scattering cross sections
         call thermr

      case('groupr')  ! compute self-shielded multigroup cross sections
         call groupr

      case('gaminr') ! compute gamma interaction cross sections.
         call gaminr

      case('errorr') ! compute cross section covariances
         call errorr

      case('covr')   ! process covariance data from errorr
         call covr

      case('moder')  ! change mode of an endf tape
         call moder

      case('dtfr')   ! produce output in dtf format
         call dtfr

      case('ccccr')  ! produce output in cccc-iii interface format
         call ccccr

      case('matxsr') ! produce output in matxs interface format
         call matxsr

      case('resxsr') ! prepare point cross sections for flux calcs
         call resxsr

      case('acer')   ! prepare data for continuous energy monte carlo
         call acer

       case('powr')  ! produce input tapes for epri codes
         call powr

      case('wimsr')  ! produce libraries for wims
         call wimsr

      case('plotr')  ! plot cross sections
         call plotr

      case('viewr')  ! view plots in postscript
         call viewr

      case('mixr')   ! combine cross sections into mixtures or elements
         call mixr

      case('purr')   ! generate unresolved-region probability tables
         call purr

      case('leapr')  ! generate s(alpha,beta) for thermal moderators
         call leapr

      case('gaspr')  ! add gas production (mt203-207) to pendf
         call gaspr

      case ('--')    ! comment card; nothing to do

      case default
         write(strng,'(''illegal module name---'',a6)') module
         call error('njoy',strng,' ')

      end select
   end do

   !--njoy complete
   call timer(secs)
   write(nsyse,'(69x,f8.1,''s''/1x,77(''*''))') secs

end program njoy

