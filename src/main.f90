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

   implicit none
   character(6)::module
   character(8)::date
   character(8)::time
   character(28)::strng
   real(kr)::secs

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
   do
      read(nsysi,*) module
      if (module.eq.'stop') exit
      select case(module)

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

