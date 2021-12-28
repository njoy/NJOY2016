module heatm
   ! provides heatr for NJOY2016
   use locale
   implicit none
   private
   public heatr

   ! global variables

   ! units for heatr
   integer::nendf,nin,nout,nplot

   ! heatr global parameters
   integer,parameter::nbuf=1000
   integer,parameter::nqamax=30
   integer,parameter::ilmax=200
   integer,parameter::maxmf6=320

   ! global variables
   integer::matd,mgam,npk,mtp(28),nqa,mta(nqamax),mt303,mt19
   integer::ne,kchk,iprint,lqs(nqamax)
   real(kr),dimension(:),allocatable::qbar
   real(kr)::qa(nqamax),efirst,elast,za,awr,elist(ilmax)
   real(kr)::break
   real(kr)::qdel,etabmax
   integer::mt103,mt104,mt105,mt106,mt107,mt16
   integer::miss4(250),nmiss4
   integer::i6,mt6(maxmf6),i6g,mt6no(maxmf6)
   integer::i6p,mt6yp(maxmf6)
   integer::jp,jpn,jpp
   real(kr)::q,zat,awrt,zap,awp
   integer::lct
   integer::idame
   real(kr)::ebot,etop
   integer::mt458,nply,lfc,ifc1,ifc2,ifc3,ifc4,ifc5,ifc6
   real(kr),dimension(:),allocatable::c458,cpoly,hpoly,afr,anp,agp
   real(kr)::emc2,tm,rtm

   ! target
   integer::izat

   ! projectile
   integer::izap

contains

   subroutine heatr
   !-------------------------------------------------------------------
   !
   ! Compute heating kerma (kinetic energy release in material)
   ! and radiation damage energy production.
   !
   ! The prompt kerma is computed pointwise on the grid of the
   ! total cross section from the input pendf tape and written
   ! onto the output PENDF tape at infinite dilution using the
   ! 300 series of MT numbers.  All temperatures on the input PENDF
   ! tape for the desired material are processed.  The dictionary
   ! is revised.  Reaction Q values are obtained from the ENDF
   ! tape unless the user enters his own value.  Partial kermas
   ! can be requested for self-shielding calculations or other
   ! purposes.  The code uses the energy balance method where
   ! photon files are available and deposits all photon energy
   ! locally when files are not available.  This assures
   ! consistency between neutron heating and energy deposition by
   ! subsequent photon interactions.  An exception is made for
   ! capture where recoil is computed by momentum conservation.
   ! Photon files are used to estimate the average photon momentum
   ! when available.  A diagnostic message is printed if the
   ! momentum calculation leads to a significant error in
   ! energy conservation.
   !
   ! If desired, the energy-balance kerma factors can be compared
   ! with conservative kinematic limits (set iprint=2).
   ! A plot file for viewr can be automatically prepared.
   !
   ! Damage energy is computed using the Lindhard electronic
   ! screening damage function with a displacement threshold
   ! from a table of default values for important elements
   ! or a value provided by the user.
   !
   !---input specifications (free format)---------------------------
   !
   ! card 1
   !    nendf    unit for endf tape
   !    nin      unit for input pendf tape
   !    nout     unit for output pendf tape
   !    nplot    unit for graphical check output
   ! card 2
   !    matd     material to be processed
   !    npk      number of partial kermas desired (default=0)
   !    nqa      number of user q values (default=0)
   !    ntemp    number of temperatures to process
   !             (default=0, meaning all on pendf)
   !    local    0/1=gamma rays transported/deposited locally
   !             (default=0)
   !    iprint   print (0 min, 1 max, 2 check) (default=0)
   !    ed       displacement energy for damage
   !             (if negative, use default from built-in table)
   !             Note, unlike in the original NJOY-2016 formulation, zero is now a valid entry 
   !                   since it is used to generate the displacement kerma
   !    idam_fnc flag for alternate treatment of damage partition function
   !             (default is 0 for built-in Robinson partition function)
   !                 0 = Robinson partition function /
   !                 1 = read-in function / 
   !                 2 = read-in efficiency factor and apply to Robinson partition function /
   !                 3 = read-in partition function and efficiency factor and use product /
   !                Note, icntrl(4) entry:
   !                  icntrl(4) = 0 & icntrl(8) = 0   corresponds to idam_fnc = 0 - NJOY default
   !                  icntrl(4) = 0 & icntrl(8) = 1   corresponds to idam_fnc = 1 - read-in damage partition function
   !                  icntrl(4) = 1 & icntrl(8) = 0   corresponds to idam_fnc = 2 - read-in efficiency function and apply on top of Robinson partition function
   !                  icntrl(4) = 1 & icntrl(8) = 1   corresponds to idam_fnc = 3 - read-in efficiency function and apply on top of read-in partition function
   !    al_new   lattice atom atomic mass
   !             (default from ENDF file header cards)
   !                if positive, units of amu
   !                if negative, units of neutron mass
   !                Note, NJOY HEATR module uses logic based on ENDF-6 format/definitions,
   !                i.e. neutron mass is used.
   !    zl_new   lattice atom atomic number
   !             (default from ENDF file header cards)
   !    filename_damage_efficiency  
   !             filename for user-defined damage efficiency function, max lenght 60 characaters
   !    filename_partition_func
   !             filename for user-defined damage partition function if not Robinson partition function, max lenght 60 characaters
   ! card 3      for npk gt 0 only
   !    mtk      mt numbers for partial kermas desired
   !             total (mt301) will be provided automatically.
   !             partial kerma for reaction mt is mt+300
   !             and may not be properly defined unless
   !             a gamma file for mt is on endf tape.
   !             special values allowed--
   !               303   non-elastic (all but mt2)
   !               304   inelastic (mt51 thru 91)
   !               318   fission (mt18 or mt19, 20, 21, 38)
   !               401   disappearance (mt102 thru 120)
   !               442   total photon ev-barns
   !               443   total kinematic kerma (high limit)
   !             damage energy production values--
   !               444   total
   !               445   elastic (mt2)
   !               446   inelastic (mt51 thru 91)
   !               447   disappearance (mt102 thru 120)
   !             Note, the damage energy formulation is determined by icntrl(1) entry:
   !               icntrl(1) =  0    output displacement energy (with lower imtegration bound of Ed) - NJOY default
   !                                     (redundant with icntrl(1)=8 option)
   !                            5    output NRT damage energy (breakpoint = 2*Ed/beta)
   !                            6    output original 3-level Kinchin-Pease damage energy (breakpoint = 2*Ed, no beta=0.8 factor)
   !                            7    output displacement kerma (i.e. Ed = 0 eV, no breakpoint)
   !                            8    output sharp transition Kinchin-Pease damage energy - no transistion region (breakpoint = Ed)
   !          cards 4 and 5 for nqa gt 0 only
   ! card 4
   !    mta      mt numbers for users q values
   ! card 5
   !    qa       user specified q values (ev)
   !               (if qa.ge.99.e6, read in variable qbar
   !                  for this reaction)
   ! card 5a     variable qbar (for reactions with qa flag only)
   !    qbar      tab1 record giving qbar versus e (1000 words max)
   !
   !-------------------------------------------------------------------
   use mainio  ! provides nsysi,nsyso
   use util    ! provides timer,openz,error,mess,repoz,closz
   use endf    ! provides endf routines and variables
   use physics ! provides amassn,amu,ev,clight
   use snl     ! provides SNL
   ! internals
   integer::ntemp,local,npkk,i,nz,j,isave,iz,loc,nzq,nz0
   integer::nscr,nend4,nend6,iold,inew,itemp,idone,nb,nw
   real(kr)::time,flag
   integer,parameter::npkmax=28
   integer::mtk(2+npkmax)
   real(kr)::z(17)
   real(kr),dimension(:),allocatable::tmp
   real(kr),parameter::qtest=99.e6_kr
   real(kr),parameter::qflag=-1.e-9_kr
   real(kr),parameter::zero=0
   integer,parameter::maxqbar=10000
   real(kr) ::   al_new_temp, amu_over_nmass, al_new_amu, al_new_nmass
   character(60) :: filename_damage_efficiency, filename_partition_func
   integer:: lk

   !--start
   call timer(time)
   write(nsyso,&
     '(/'' heatr...prompt kerma'',48x,f8.1,''s'')') time
   write(nsyse,'(/'' heatr...'',60x,f8.1,''s'')') time

   if (imode(3) .lt. 0) then
     write (nsyso,2301) imode(3)
2301 format (/,1x, 'HEATR debug print invoked ', i3, /)
   endif

   !--read user input.
   nplot=0
   read(nsysi,*) nendf,nin,nout,nplot
   call openz(nendf,0)
   call openz(nin,0)
   call openz(nout,1)
   if (nplot.ne.0) call openz(nplot,1)
   npk=0
   nqa=0
   ntemp=0
   local=0
   iprint=0
   break=-1.0
   idam_fnc = 0
   filename_damage_efficiency = 'not_defined'
   filename_partition_func = 'not_defined'
   al_new_temp = 0.0
   zl_new = 0.0
   read(nsysi,*) matd,npk,nqa,ntemp,local,iprint,break, idam_fnc, al_new_temp, zl_new,  &
 &               filename_damage_efficiency, filename_partition_func
   break_new = break
   kchk=0
   if (iprint.eq.2) kchk=1
   if (iprint.eq.2) iprint=1
   npkk=npk+3
   if (kchk.eq.1) npkk=3*npk+7
   if (npkk.gt.npkmax) then 
      write (nsyso, 3981) iprint, kchk, npkk, npkmax, npk 
 3981 format (1x, 'HEATR requests too many kermas: ', 5i8)
      call error('heatr',&
        'requested too many kerma mt-s (6+mt301 allowed).',' ')
   endif
   if (nqa.gt.nqamax)&
     call error('heatr','requested too many q values.',' ')
   if (npk.gt.0) then
      read(nsysi,*) (mtk(i),i=1,npk)
   endif
   allocate(tmp(maxqbar))
   loc=1
   if (nqa.ne.0) then
      read(nsysi,*) (mta(i),i=1,nqa)
      read(nsysi,*) (qa(i),i=1,nqa)
      do i=1,nqa
         lqs(i)=0
         if (qa(i).ge.qtest) then
            lqs(i)=loc
            nz=1000
            flag=qflag
            do j=1,nz
               tmp(loc+j-1)=flag
            enddo
            read(nsysi,*) (tmp(loc+j-1),j=1,nz)
            nz0=0
            do j=1,nz
               if (tmp(loc+j-1).gt.flag) nz0=j
            enddo
            loc=loc+nz0
            if (loc.gt.maxqbar) call error('heatr',&
              'too much energy-dependent q data',' ')
         endif
      enddo
   endif
   nzq=loc-1
   if (nzq.gt.0) then
      allocate(qbar(nzq))
      do i=1,nzq
         qbar(i)=tmp(i)
      enddo
   endif
   deallocate(tmp)
   if ((nin.lt.0.and.nout.gt.0).or.(nin.gt.0.and.nout.lt.0))&
     call error('heatr',&
     'mode conversion not allowed between nin and nout.',' ')

   !--check the partial kermas
   if (npk.gt.0) then
      isave=0
      do i=1,npk
         if (mtk(i).eq.443.and.kchk.eq.0) kchk=2
         if (mtk(i).eq.301) isave=i
      enddo
      if (isave.gt.0) then
         call mess('heatr','mt301 always calculated',&
           '--you do not need to ask for it.')
         npk=npk-1
         if (isave.le.npk) then
            do i=isave,npk
               mtk(i)=mtk(i+1)
            enddo
         endif
      endif
      do i=1,npk
         mtp(npk+3-i)=mtk(i)
      enddo
   endif
   mtp(1)=0
   mtp(2)=301
   npk=npk+2
   ! there are always at least 2 elements in the npk array,
   ! (1)=energy, (2)=mt301 (total)
   if (npk.gt.3) call horder(mtp,npk)

   !--list all input parameters
   write(nsyso,'(/&
     &'' input endf unit ...................... '',i10/&
     &'' input pendf unit ..................... '',i10/&
     &'' output pendf unit .................... '',i10)')&
     nendf,nin,nout
   write(nsyso,'(&
     &'' mat to be processed .................. '',i10)') matd
   write(nsyso,'(&
     &'' no. temperatures (0=all) ............. '',i10/&
     &'' gamma heat (0 nonlocal, 1 local) ..... '',i10/&
     &'' print option (0 min, 1 more, 2 chk) .. '',i10)')&
     ntemp,local,iprint

!  over-ride previous convevntion and permit a zero displacement threshold energy
   if (break.lt.zero) then
      write(nsyso,'(&
        &'' damage displacement energy ...........    default'')')
   else
      write(nsyso,'(&
        &'' damage displacement energy ........... '',f11.2,&
        &'' ev'')') break
   endif
   displace_th = break_new

   if ( idam_fnc .eq. 0) then 
      if ( icntrl(8) .ne. 0 .or. icntrl(4) .ne. 0) then 
         write(nsyso,'(&
           &'' ERROR: idam_fnc flag = 0 input conflict'', 3i5)') icntrl(4), icntrl(8), idam_fnc
         call error('heatr','idam_fnc in conflict with icntrl(4); icntrl(8)',' ')
      endif
   elseif ( idam_fnc .eq. 1) then 
      if ( icntrl(8) .ne. 1 .or. icntrl(4) .ne. 0) then 
         write(nsyso,'(&
           &'' ERROR: idam_fnc = 1 flag input conflict'', 3i5)') icntrl(4), icntrl(8), idam_fnc
         call error('heatr','idam_fnc in conflict with icntrl(4); icntrl(8)',' ')
      endif
   elseif ( idam_fnc .eq. 2) then
      if ( icntrl(8) .ne. 0 .or. icntrl(4) .ne. 1) then 
         write(nsyso,'(&
           &'' ERROR: idam_fnc=2 flag input conflict'', 3i5)') icntrl(4), icntrl(8), idam_fnc
         call error('heatr','idam_fnc in conflict with icntrl(4); icntrl(8)',' ')
      endif
   elseif ( idam_fnc .eq. 3) then
      if ( icntrl(8) .ne. 1 .or. icntrl(4) .ne. 1) then 
         write(nsyso,'(&
           &'' ERROR: idam_fnc =3 flag input conflict'', 3i5)') icntrl(4), icntrl(8), idam_fnc
         call error('heatr','idam_fnc in conflict with icntrl(4); icntrl(8)',' ')
      endif
   else
         write(nsyso,'(&
           &'' ERROR: illegal idam_fnc input flag'', 3i5)') icntrl(4), icntrl(8), idam_fnc
         call error('heatr','idam_fnc /= icntrl(4), icntrl(8), idam_fnc',' ')
   endif
!
! Damage efficiency function input/modification options
!
   if ( icntrl(4) .eq. 1) then 
         write(nsyso,'(&
          &'' Damage efficiency filename ='', a60)') filename_damage_efficiency
         open(unit=34,status='old',file=filename_damage_efficiency)
         read (34, *) npoints_eff
         write(nsyso,'(&
        &'' Number of points in damage efficiency ='', 5x, i5)') npoints_eff
        open(unit=34,status='old',file=filename_damage_efficiency)
        if ( npoints_eff .gt. 250) then 
            write(nsyso,'(&
           &'' ERROR: illegal npoints_eff, >250 '', i5)') npoints_eff
            call error('heatr','efficiency curve has too many points ',' ')
        endif
        write(nsyso,'(&
        &''                '')')
        write(nsyso,'(&
        &'' Tabulated damage efficiency points: '')')
        write(nsyso,'(28x, &
        &'' PKA Energy (eV)      Correction Factor '')')
        write(nsyso,'(&
        &''                '')')
        do lk = 1, npoints_eff
          read(34,*) eff_eng(lk), eff_value(lk)
         write(nsyso,'(&
        &''           damage efficiency ='', g14.7, 6x, g14.7)') eff_eng(lk), eff_value(lk)
        open(unit=34,status='old',file=filename_damage_efficiency)
        enddo
        close(unit=34)
        write(nsyso,'(&
        &''                '')')
   elseif (icntrl(4) .gt. 1) then 
         write(nsyso,'(&
        &'' ERROR: illegal icntrl(4) input flag'', 3i5)') icntrl(4), icntrl(8), idam_fnc
      call error('heatr','icntrl(4) efficiency function option not implemented',' ')
   elseif (icntrl(4) .lt. 0) then
         write(nsyso,'(&
        &'' ERROR: illegal icntrl(4) input flag'', 3i5)') icntrl(4), icntrl(8), idam_fnc
      call error('heatr','icntrl(4) efficiency function option not implemented',' ')
   endif
!
! Damage partition function modifications
!
   if ( icntrl(8) .ne. 0 .and. icntrl(8) .ne. 1) then 
         write(nsyso,'(&
        &'' ERROR: illegal icntrl(8) input flag'', 3i5)') icntrl(4), icntrl(8), idam_fnc
          write(nsyso,'(&
        &'' Damage partition function filename = '', a60)') filename_partition_func
      call error('heatr','icntrl(8) partition function read not implemented',' ')
   endif
   if ( icntrl(8) .eq. 1) then 
         write(nsyso,'(&
          &'' Partition function filename = '', a60)') filename_partition_func
         open(unit=34,status='old',file=filename_partition_func)
         read (34, *) ndamage
         write(nsyso,'(&
        &'' Number of points in partition function ='', 5x, i5)') ndamage
        open(unit=34,status='old',file=filename_partition_func)
        if ( ndamage .gt. 1000) then 
            write(nsyso,'(&
           &'' ERROR: illegal ndamage, >1000 '', i5)') ndamage
            call error('heatr','partition function curve has too many points ',' ')
        endif
        write(nsyso,'(&
        &''                '')')
        write(nsyso,'(&
        &'' Tabulated damage partition function points: '')')
        write(nsyso,'(28x, &
        &'' PKA Energy (eV)      Partition Fraction '')')
        write(nsyso,'(&
        &''                '')')
        do lk = 1, ndamage
          read(34,*) energy_dam(lk), value_dam(lk)
         write(nsyso,'(&
        &''           partition function = '', g14.7, 6x, g14.7)') energy_dam(lk), value_dam(lk)
        open(unit=34,status='old',file=filename_partition_func)
        enddo
        close(unit=34)
        write(nsyso,'(&
        &''                '')')
   elseif (icntrl(8) .gt. 1) then 
         write(nsyso,'(&
        &'' ERROR: illegal icntrl(8) input flag'', 3i5)') icntrl(4), icntrl(8), idam_fnc
      call error('heatr','icntrl(8) partition function option not implemented',' ')
   elseif (icntrl(8) .lt. 0) then
         write(nsyso,'(&
        &'' ERROR: illegal icntrl(8) input flag'', 3i5)') icntrl(4), icntrl(8), idam_fnc
      call error('heatr','icntrl(8) partition function option not implemented',' ')
   endif

   if (zl_new.eq.zero) then
      write(nsyso,'(&
        &'' lattice atomic number ................    default'')')
   else
      write(nsyso,'(&
        &'' lattice atom atomic number ........... '',f11.2&
        & )') zl_new
   endif

   if (al_new_temp.eq.zero) then
      write(nsyso,'(&
        &'' lattice atomic mass ..................    default'')')
   else
      write(nsyso,'(&
        &'' lattice atom atomic mass ............. '',g14.7&
        &)') al_new_temp
   endif
   amu_over_nmass =  931.49410242 / 939.56542052
   if ( al_new_temp .lt. 0.0) then
!      input nmass
        al_new_amu = abs(al_new_temp) * amu_over_nmass
        al_new_nmass = abs(al_new_temp)
   else
!      input amu
       al_new_amu = al_new_temp
       al_new_nmass = al_new_temp / amu_over_nmass
   endif
   al_new = al_new_nmass
   write (nsyso, 1903) break_new, al_new_amu,  &
 &         al_new_nmass, zl_new
 1903  format (/,12x, 'dpa parameters changed = ', g14.7,&
     &          '  ev ',/, &
     &       12x, 'al parameter changed   = ', g14.7, ' amu',/,&
     &       12x, '                       = ', g14.7, ' nmass',/, &
     &       12x, 'zl parameter changed   = ', g14.7, &
     &       '  Z  ',/)
   if ( zl_new.ne.zero .or. al_new_temp.ne.zero) then
     if ( icntrl(5) .ne. 1) then 
        write (nsyso, 7903) icntrl(5)
 7903   format (12x, 'WARNING: icntrl(5) not set to 1, but lattice atom information changed ', i5)
     endif 
   endif

   if (ntemp.eq.0) ntemp=100
   if (npk.ge.3) write(nsyso,'(&
     &'' partial kerma mt-s desired ........... '',i10)') mtp(3)
   if (npk.ge.4) write(nsyso,'(40x,i10)')(mtp(i),i=4,npk)
   if (nqa.ne.0) then
      write(nsyso,'(&
        &'' auxiliary reactions .................. '',i10)') mta(1)
      if (nqa.gt.1) write(nsyso,'(40x,i10)') (mta(i),i=2,nqa)
      write(nsyso,'(&
        &'' auxiliary q values ................... '',1p,e12.4)') qa(1)
      if (nqa.gt.1) write(nsyso,'(40x,1p,e12.4)') (qa(i),i=2,nqa)
      if (nzq.ne.0) then
         do i=1,nqa
            if (qa(i).ge.qtest) then
               nsh=0
               math=1
               mfh=1
               mth=mta(i)
               write(nsyso,'('' input q ----'')')
               call tab1io(0,nsyso,0,qbar,nb,nw)
            endif
         enddo
      endif
   endif

   !--assign scratch units
   nscr=10
   if (nendf.lt.0) nscr=-nscr
   nend4=13
   if (nendf.lt.0) nend4=-nend4
   nend6=14
   if (nendf.lt.0) nend6=-nend6
   iold=11
   inew=12
   call openz(-iold,1)
   call openz(-inew,1)
   call openz(nscr,1)
   call openz(nend4,1)
   call openz(nend6,1)

   !--process input pendf tape
   call repoz(nin)
   call repoz(nout)
   nsh=0
   call tpidio(nin,nout,0,z,nb,nw)
   call findf(matd,1,451,nin)
   itemp=0

   !--loop over desired temperatures for this mat
   idone=0
   do while (idone.eq.0)
      call contio(nin,0,0,z,nb,nw)
      if (math.ne.matd) then
         idone=1
      else
         itemp=itemp+1
         if (itemp.gt.ntemp) then
            idone=1
         else
            za=c1h
            awr=c2h
            emc2=amassn*amu*clight*clight/ev
            tm=emc2*(awr+1)
            rtm=1/tm

            !--check on default damage displacement energy
            if (break.lt.zero) then
               iz=nint(za/1000)
               if (iz.eq.4) then
                  break=31
               else if (iz.eq.6) then
                  break=31
               else if (iz.eq.12) then
                  break=25
               else if (iz.eq.13) then
                  break=27
               else if (iz.eq.14) then
                  break=25
               else if (iz.eq.20) then
                  break=40
               else if (iz.ge.22.and.iz.le.29) then
                  break=40
               else if (iz.eq.40) then
                  break=40
               else if (iz.eq.41) then
                  break=40
               else if (iz.eq.42) then
                  break=60
               else if (iz.eq.47) then
                  break=60
               else if (iz.eq.73) then
                  break=90
               else if (iz.eq.74) then
                  break=90
               else if (iz.eq.79) then
                  break=30
               else if (iz.eq.82) then
                  break=25
               else
                  break=25
               endif
               write(nsyso,'(/'' default damage energy ='',f5.1,&
                 &'' ev'')') break
            endif

            !--analyze dictionary and read energy grids.
            call hinit(iold,nend4,nend6,nscr)

            !--calculate heating from neutron files.
            call timer(time)
            write(nsyso,'(&
              &61x,''temp'',i2,2x,f8.1,''s'')') itemp,time
            call nheat(iold,inew,nscr,nend4,nend6,local)

            !--correct for gamma production.
            if (mgam.gt.0.and.mgam.ne.10.and.local.eq.0) then
               call timer(time)
               write(nsyso,'(&
                 &61x,''temp'',i2,2x,f8.1,''s'')') itemp,time
               call gheat(iold,inew,nscr)
            endif

            !--write kermas on pendf tape.
            call hout(iold)
         endif
      endif
   enddo

   !--heatr is finished.
   if (allocated(c458)) deallocate(c458)
   if (allocated(cpoly)) deallocate(cpoly)
   if (allocated(hpoly)) deallocate(hpoly)
   if (allocated(afr)) deallocate(afr)
   if (allocated(anp)) deallocate(anp)
   if (allocated(agp)) deallocate(agp)
   call atend(nout,0)
   call repoz(nout)
   call repoz(nin)
   call repoz(nendf)
   call closz(nin)
   call closz(nout)
   call closz(nendf)
   call closz(-iold)
   call closz(-inew)
   call closz(nscr)
   call closz(nend4)
   call closz(nend6)
   call timer(time)
   write(nsyso,'(/69x,f8.1,''s''/&
     &1x,7(''**********''),''*******'')') time
   return
   end subroutine heatr

   subroutine horder(iarray,n)
   !-------------------------------------------------------------------
   ! Arrange the n elements of iarray in ascending order.
   !-------------------------------------------------------------------
   ! externals
   integer::n,iarray(n)
   ! internals
   integer::i,j,isave

   do i=1,n
      do j=i,n
         if (iarray(i).gt.iarray(j)) then
            isave=iarray(i)
            iarray(i)=iarray(j)
            iarray(j)=isave
         endif
      enddo
   enddo
   return
   end subroutine horder

   subroutine hinit(iold,nend4,nend6,nscr)
   !-------------------------------------------------------------------
   ! Analyze dictionary and set necessary flags.
   ! Read energy grid of total cross section.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util   ! provides repoz,error,mess,loada
   use endf   ! provides endf routines and variables
   ! externals
   integer::iold,nend4,nend6,nscr
   ! internals
   integer::ifiss,i519,nb,nw,n6,nx,i,ifini,ii6,lf,n,nk,nfc,ifc
   integer::mfd,mtd,j,idone,mf4,mt4,ie,nmu,imu,npkk,idis
   integer::ielem,ik,nr,idnx,nen,law
   integer::lrel,nmod
   integer::i6t,nl,l
   real(kr)::zar,awrr,e,enext,y,yld,test,efix
   real(kr)::c(30)
   character(60)::strng
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::bufo
   real(kr),parameter::fact=.99999e0_kr
   real(kr),parameter::rup=1.000000001e0_kr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::zero=0
   real(kr),parameter::mevev=1.e-6_kr
   real(kr),parameter::tol=1.e-5_kr

   !--analyze dictionary
   lrel=0
   nmod=0
   i6=0
   mt103=0
   mt104=0
   mt105=0
   mt106=0
   mt107=0
   mt16=0
   nmiss4=0
   mgam=0
   ifiss=0
   mt19=0
   i519=0
   mt458=0
   qdel=0
   jp=0
   jpn=0
   jpp=0
   etop=20000000
   ebot=1
   ebot=ebot/100000
   etabmax=-1.0
   call repoz(nendf)
   call tpidio(nendf,0,0,c(1),nb,nw)
   call findf(matd,1,451,nendf)
   call contio(nendf,0,0,c(1),nb,nw)
   n6=n2h
   call contio(nendf,0,0,c(1),nb,nw)
   if (n1h.ne.0) then
      iverf=4
      nx=n6
   else if (n2h.eq.0) then
      iverf=5
   else
      iverf=6
   endif
   call skiprz(nendf,-1)
   if (iverf.ge.5) call contio(nendf,0,0,c(1),nb,nw)
   if (iverf.eq.6) call contio(nendf,0,0,c(1),nb,nw)
   if (iverf.eq.6) then
      if (c2h.gt.etop) etop=c2h
      lrel=l1h
      nmod=n2h
   endif
   allocate(scr(npage+50))
   call hdatio(nendf,0,0,scr,nb,nw)
   do while (nb.ne.0)
      call moreio(nendf,0,0,scr,nb,nw)
   enddo
   if (iverf.ge.5) nx=n2h
   deallocate(scr)
   nw=6*nx
   if (nw.lt.npage+50) nw=npage+50
   allocate(scr(nw))
   nw=nx
   if (nw.ne.0) then
      call dictio(nendf,0,0,scr,nb,nw)
      i=1
      ifini=0
      do while (i.le.nw.and.ifini.eq.0)
         mfd=nint(scr(i+2))
         mtd=nint(scr(i+3))
         if (mfd.ne.2) then
            if (mfd.eq.12.and.mtd.ne.460) mgam=1
            if (mfd.eq.13) mgam=1
            if (mfd.ge.13) then
               ifini=1
            else
               if (mfd.le.3) then
                  if (mfd.eq.1.and.mtd.eq.458) mt458=1
                  if (mtd.eq.19) mt19=1
                  if (mtd.eq.18.or.mtd.eq.19) ifiss=1
                  if (iverf.ge.6) then
                     if (mtd.ge.600.and.mtd.lt.650) mt103=1
                     if (mtd.ge.650.and.mtd.lt.700) mt104=1
                     if (mtd.ge.700.and.mtd.lt.750) mt105=1
                     if (mtd.ge.750.and.mtd.lt.800) mt106=1
                     if (mtd.ge.800.and.mtd.lt.850) mt107=1
                     if (mtd.ge.875.and.mtd.lt.891) mt16=1
                     if (mtd.ge.600.and.mtd.le.850) then
                        j=i
                        idone=0
                        do while (j.le.nw.and.idone.eq.0)
                           mf4=nint(scr(j+2))
                           mt4=nint(scr(j+3))
                           if (mf4.eq.4.and.mt4.eq.mtd) then
                              idone=1
                           else
                              if (mf4.eq.6.and.mt4.eq.mtd) then
                                 idone=1
                              else
                                 j=j+6
                              endif
                           endif
                        enddo
                        if (idone.eq.0) then
                           write(strng,'(&
                             &''mf4 and 6 missing, '',&
                             &''isotropy assumed for mt '',i3)') mtd
                           call mess('hinit',strng,' ')
                           nmiss4=nmiss4+1
                           miss4(nmiss4)=mtd
                        endif
                     endif
                  else
                     if (mtd.ge.700.and.mtd.lt.720) mt103=1
                     if (mtd.ge.720.and.mtd.lt.740) mt104=1
                     if (mtd.ge.740.and.mtd.lt.760) mt105=1
                     if (mtd.ge.760.and.mtd.lt.780) mt106=1
                     if (mtd.ge.780.and.mtd.lt.800) mt107=1
                     if (mtd.ge.700.and.mtd.le.800) then
                        j=i
                        idone=0
                        do while (j.le.nw.and.idone.eq.0)
                           mf4=nint(scr(j+2))
                           mt4=nint(scr(j+3))
                           if (mf4.eq.4.and.mt4.eq.mtd) then
                              idone=1
                           else
                              if (mf4.eq.6.and.mt4.eq.mtd) then
                                 idone=1
                              else
                                 j=j+6
                              endif
                           endif
                        enddo
                        if (idone.eq.0) then
                           write(strng,'(&
                             &''mf4 and 6 missing, '',&
                             &''isotropy assumed for mt '',i3)') mtd
                           call mess('hinit',strng,' ')
                           nmiss4=nmiss4+1
                           miss4(nmiss4)=mtd
                        endif
                     endif
                  endif
               else
                  if (mfd.eq.5) then
                     if (mtd.eq.19) i519=1
                  else if (mfd.eq.6) then
                     i6=i6+1
                     if (i6.gt.maxmf6) call error('hinit',&
                       'too many mf6 reactions',' ')
                     mt6(i6)=mtd
                  endif
               endif
            endif
         endif
         i=i+6
      enddo
   endif
   if (mt19.eq.1.and.i519.eq.0) mt19=2
   if (mt19.eq.1) call mess('hinit','mt18 is redundant',&
     '19, 20, 21, ... will be used.')
   if (mt19.eq.2) call mess('hinit','mt19 has no spectrum',&
     'mt18 spectrum will be used.')

   !--calculate q change for fissionables in mt 458
   if (ifiss.ne.0) then
      deallocate(scr)
      nw=10000
      allocate(scr(nw))
      if (allocated(c458)) deallocate(c458)
      if (allocated(cpoly)) deallocate(cpoly)
      if (allocated(hpoly)) deallocate(hpoly)
      if (allocated(afr)) deallocate(afr)
      if (allocated(anp)) deallocate(anp)
      if (allocated(agp)) deallocate(agp)
      if (mt458.eq.1) then
         call findf(matd,1,458,nendf)
         call contio(nendf,0,0,scr,nb,nw)
         lfc=nint(scr(4))
         nfc=nint(scr(6))
         ifc1=0
         ifc2=0
         ifc3=0
         ifc4=0
         ifc5=0
         ifc6=0
         call listio(nendf,0,0,scr,nb,nw)
         nply=nint(scr(4))
         if (lfc.eq.0) then
            allocate(cpoly(0:nply))
            allocate(hpoly(0:nply))
            cpoly=0
            hpoly=0
         endif
         allocate(c458(nint(scr(5))))
         ! save polynomial coefficients ... c0 & c1 terms are ok,
         ! but for endf/b-vii.1, c2 and higher are mistakenly given
         ! as if the energy will be given in MeV and therefore need
         ! correction.
         write(nsyso,'(/,&
           &'' fission energy components''/)')
         c458(1:18)=scr(7:24)
         if (lfc.eq.0) then      ! thermal point or polynomial
            qdel=c458(5)+c458(9)+c458(11) ! delayed fission Q at 0 eV
            if (nply.ge.1) then
               c458(19:36)=scr(25:42)
               if (nply.gt.1) then
                  do n=2,nply
                     if (nmod.ne.7.or.lrel.ne.1) then
                        efix=1
                     else
                        efix=mevev**(n-1)
                     endif
                     c458(1+18*n:18+18*n)=scr(7+18*n:24+18*n)*efix
                  enddo
               endif
               write(nsyso,'(&
                 &''   fission products: polynomial of order '',i2/&
                 &''   prompt neutrons : polynomial of order '',i2/&
                 &''   prompt gammas   : polynomial of order '',i2)')&
                 nply,nply,nply
            else
               write(nsyso,'(&
                 &''   fission products: thermal point''/&
                 &''   prompt neutrons : thermal point''/&
                 &''   prompt gammas   : thermal point'')')
            endif
            do  i=0,nply
               cpoly(i)=c458(1+i*18)+c458(3+i*18)+c458(7+i*18)
               hpoly(i)=c458(1+i*18)+c458(7+i*18)
            enddo
         else if (lfc.eq.1) then ! tabulated components
            qdel=0
            do i=1,nfc
               l=1
               call tab1io(nendf,0,0,scr,nb,nw)
               do while (nb.ne.0)
                  l=l+nw
                  call moreio(nendf,0,0,scr(l),nb,nw)
               enddo
               ifc=l2h
               nl=6+2*n1h+2*n2h
               if (etabmax.lt.0) then
                  etabmax=scr(nl-1)
               else
                  if ((abs(etabmax-scr(nl-1))/etabmax).gt.tol) then
                     write(strng,'(''upper energy mismatch for ifc='',i2,&
                        &'' in mt=458.'')')ifc
                     call error('hinit',strng,' ')
                  endif
               endif
               if (ifc.eq.1) then      ! EFR is tabulated
                  ifc1=1
                  allocate(afr(1:nl))
                  afr(1:nl)=scr(1:nl)
               else if (ifc.eq.2) then ! ENP is tabulated
                  ifc2=1
                  allocate(anp(1:nl))
                  anp(1:nl)=scr(1:nl)
               else if (ifc.eq.3) then ! END is tabulated
                  ifc3=1
                  qdel=qdel+scr(7)
               else if (ifc.eq.4) then ! EGP is tabulated
                  ifc4=1
                  allocate(agp(1:nl))
                  agp(1:nl)=scr(1:nl)
               else if (ifc.eq.5) then ! EGD is tabulated
                  ifc5=1
                  qdel=qdel+scr(7)
               else if (ifc.eq.6) then ! EB is tabulated
                  ifc6=1
                  qdel=qdel+scr(7)
               endif
            enddo
            if (ifc3.eq.0) qdel=qdel+c458(5)
            if (ifc5.eq.0) qdel=qdel+c458(9)
            if (ifc6.eq.0) qdel=qdel+c458(11)
            if (ifc1.eq.1) then
               write(nsyso,'(''   fission products: tabulated values'')')
            else
               write(nsyso,'(''   fission products: thermal point'')')
            endif
            if (ifc2.eq.1) then
               write(nsyso,'(''   prompt neutrons : tabulated values'')')
            else
               write(nsyso,'(''   prompt neutrons : thermal point'')')
            endif
            if (ifc4.eq.1) then
               write(nsyso,'(''   prompt gammas   : tabulated values'')')
            else
               write(nsyso,'(''   prompt gammas   : thermal point'')')
            endif
            if (etabmax.lt.0) then
               write(strng,'(''no tabulated fission q components found'',&
                  &'' in mt=458.'')')
               call error('nheat',strng,' ')
            endif
         else
            call error('hinit','bad LFC in mt=458.',' ')
         endif
      else
         lfc=0
         nply=0
         allocate(cpoly(0:nply))
         cpoly=0
         call mess('hinit','mt458 is missing for this mat',' ')
      endif
   endif

   !--prepare a scratch tape for file 6
   !--check for missing recoil subsections
   !--and check for photon data
   i6p=0
   mt6yp=0
   if (i6.gt.0) then
      call repoz(nend6)
      nsp=1
      math=1
      call afend(nend6,0)
      call findf(matd,6,0,nendf)
      ii6=0
   endif
   mfh=6
   do while (mfh.eq.6)
      call contio(nendf,0,0,scr,nb,nw)
      if (mth.eq.18) then
         jp=l1h
         jpn=mod(jp,10)
         jpp=jp-10*jpn
      endif
      if (mth.ne.18 .or. (mth.eq.18 .and. jp.eq.0)) then
         call contio(0,nend6,0,scr,nb,nw)
      else
         !--delete mt=18 from mt6 array and skip this section (for now)
         i6t=i6
         do i=1,i6t
            if (mt6(i).lt.18) cycle
            if (mt6(i).eq.18) i6=i6-1
            if (mt6(i).gt.18 .and. i6.ne.i6t) mt6(i-1)=mt6(i)
         enddo
         call tosend(nendf,0,0,scr)
         cycle
      endif
      if (mfh.ne.0) then
         if (mth.ne.0) then
            nk=n1h
            ii6=ii6+1
            mt6no(ii6)=0
            zar=c1h+1
            awrr=c2h+1
            ielem=mod(nint(c1h),1000)
            do ik=1,nk
               call tab1io(nendf,nend6,0,scr,nb,nw)
               nr=n1h
               zap=c1h
               if (nint(zap).eq.0) then
                  i6p=i6p+1
                  mt6yp(i6p)=mth
               endif
               lf=l1h
               if (mth.eq.5) zar=0
               if (mth.ne.5) then
                  yld=scr(6+2*nr+2)
                  zar=zar-zap*yld
                  awrr=awrr-c2h*yld
                  if (zap.eq.zero.and.(ik.gt.1.or.nk.eq.1)) then
                     if (nk.ne.1.or.zap.ne.zero) then
                        if (zar.ne.zero) then
                           test=3000
                           if (ielem.ne.0.or.abs(zar).ge.test) then
                              write(strng,'(''mf6, mt'',i3,&
                                &'' does not give recoil za='',i6)')&
                                mth,nint(zar)
                              if (mth.eq.102) then
                                 call mess('hinit',strng,&
                                   'photon momentum recoil used.')
                              else if (mth.ne.10) then
                                 call mess('hinit',strng,&
                                   'one-particle recoil approx. used.')
                              endif
                              mt6no(ii6)=ik-1
                              zar=0
                           endif
                        endif
                     endif
                  endif
               endif
               if (mth.eq.102.and.zap.gt.0.and.lf.eq.0) then
                  write(strng,'(''mf6, mt'',i3,&
                   &'' has recoil with no spectrum'')') mth
                  call mess('hinit',strng,&
                                  'photon momentum recoil used.')
                  mt6no(ii6)=nk
               endif
               if (zap.eq.zero) mgam=10+mod(mgam,10)
               law=l2h
               do while (nb.ne.0)
                  call moreio(nendf,nend6,0,scr,nb,nw)
               enddo
               if (law.eq.6) then
                  call contio(nendf,nend6,0,scr,nb,nw)
               else if (law.eq.1.or.law.eq.2.or.law.eq.5) then
                  call tab2io(nendf,nend6,0,scr,nb,nw)
                  ne=n2h
                  do ie=1,ne
                     call listio(nendf,nend6,0,scr,nb,nw)
                     do while (nb.ne.0)
                        call moreio(nendf,nend6,0,scr,nb,nw)
                     enddo
                  enddo
               else if (law.eq.7) then
                  call tab2io(nendf,nend6,0,scr,nb,nw)
                  ne=n2h
                  do ie=1,ne
                     call tab2io(nendf,nend6,0,scr,nb,nw)
                     nmu=n2h
                     do imu=1,nmu
                        call tab1io(nendf,nend6,0,scr,nb,nw)
                        do while (nb.ne.0)
                           call moreio(nendf,nend6,0,scr,nb,nw)
                        enddo
                     enddo
                  enddo
               else if (law.lt.0) then
                  call skip6(nendf,nend6,0,scr,law)
               endif
            enddo
            if (zar.ne.zero) then
               test=3000
               if (ielem.ne.0.or.abs(zar).ge.test) then
                  write(strng,'(''mf6, mt'',i3,&
                    &'' does not give recoil za='',i6)') mth,nint(zar)
                  if (mth.eq.102) then
                     call mess('hinit',strng,&
                       'photon momentum recoil used.')
                  else if (mth.ne.10) then
                     call mess('hinit',strng,&
                       'one-particle recoil approx. used.')
                  endif
                  mt6no(ii6)=nk
               endif
            endif
         endif
      endif
   enddo
   call amend(nend6,0)
   call atend(nend6,0)
   call repoz(nend6)
   if (mgam.eq.0) write(nsyso,'(/&
     &'' no photon production files...'',&
     &''all photon energy will be deposited locally.'')')

   !--copy endf tape to nscr for use in conbar and gheat
   !--and to nend4 for use in disbar
   !--convert photon transition probability arrays if present
   call repoz(nendf)
   call repoz(nend4)
   call repoz(nscr)
   nsp=1
   nsc=1
   math=1
   call afend(nend4,nscr)
   call findf(matd,1,451,nendf)
   call hconvr(nendf,nend4,nscr)
   call amend(nend4,nscr)
   call atend(nend4,nscr)
   call repoz(nend4)
   call repoz(nscr)
   call repoz(nendf)

   !--store energy grid of the total cross section.
   call findf(matd,3,1,nin)
   npkk=npk
   if (kchk.eq.1) npkk=3*npk-2
   if (mgam.gt.0) npkk=npkk+3
   if (mgam.eq.0.and.i6.gt.0) npkk=npkk+1
   do i=1,npkk
      c(i)=0
   enddo
   call contio(nin,0,0,scr,nb,nw)
   e=0
   call gety1(e,enext,idis,y,nin,scr)
   enext=enext*rup
   efirst=enext
   idnx=0
   allocate(bufo(nbuf))
   ne=0
   nen=1
   do while (nen.gt.0)
      ne=ne+1
      test=enext*fact
      if (idnx.gt.0.and.test.gt.e) enext=test
      e=enext
      call gety1(e,enext,idis,y,nin,scr)
      nen=ne
      if (enext.gt.big-big/100) nen=-nen
      c(1)=e
      call loada(nen,c,npkk,iold,bufo,nbuf)
      if (idis.gt.0.and.idnx.gt.0) idnx=0
      if (idis.gt.0.and.idnx.eq.0) idnx=1
      if (idis.eq.0) idnx=0
   enddo
   elast=e
   call tosend(nin,0,0,scr)
   deallocate(bufo)
   deallocate(scr)
   mt303=0
   if (npk.lt.3) return
   do i=3,npk
      if (mtp(i).eq.303) mt303=i
   enddo
   return
   end subroutine hinit

   subroutine nheat(iold,inew,nscr,nend4,nend6,local)
   !-------------------------------------------------------------------
   ! Calculate contributions to heating from neutron files.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util   ! provides error,mess,loada,finda,skiprz,sigfig
   use endf   ! provides endf routines and variables
   use snl     ! provides SNL
   ! externals
   integer::iold,inew,nscr,nend4,nend6,local
   ! internals
   integer::nb,nw,npkk,ipk,n6,j6,irec,jrec,last6,new6
   integer::lr,icon,mtd,i,iimt,nmt,na,nwa,nwm,idis
   integer::ipx,irx,ilist,nn,iii,lq0,idx,ii,indxx
   integer::ipfr,irfr,ipnp,irnp,ipgp,irgp
   integer::npkt,intt,isave,npkkk,nwmax,iflag
   real(kr)::awfac,aw1fac,e,thresh,y,t,qendf,pnue
   real(kr)::qs,q0,h0,ebar,dame,yld,qsave
   real(kr)::elst,yld0,xxx,test,enext,enx,ebal6,h
   real(kr)::ebarp,yldp,yp,hp,damep,ebal6p,hmin,hmax,tt
   real(kr)::qfr,qnp,qgp
   real(kr)::c(30)
   integer::imt(30)
   character(60)::strng,strng1,strng2
   real(kr),dimension(:),allocatable::b
   real(kr),dimension(:),allocatable::a
   real(kr),dimension(:),allocatable::d
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::bufo,bufn
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::efis=15.e6_kr
   real(kr),parameter::fq1=8.07e6_kr
   real(kr),parameter::fq2=1.307e0_kr
   real(kr),parameter::one=1
   real(kr),parameter::ten=10
   real(kr),parameter::qtest=99.e6_kr
   real(kr),parameter::smin=1.e-9_kr
   real(kr),parameter::break1=1.e6_kr
   real(kr),parameter::break2=2.e6_kr
   real(kr),parameter::break3=5.e6_kr
   real(kr),parameter::break4=20.e6_kr
   real(kr),parameter::step1=.2e6_kr
   real(kr),parameter::step2=.5e6_kr
   real(kr),parameter::step3=1.e6_kr
   real(kr),parameter::step4=5.e6_kr
   real(kr),parameter::up=1.1e0_kr
   real(kr),parameter::zero=0
   real(kr),parameter::qsmall=1e-6_kr
   real(kr),parameter::tol=1.e-5_kr

   !--allocate storage
   allocate(b(npage+50))
   nwmax=250000
   allocate(d(nwmax))
   na=500000
   allocate(a(na))
   allocate(scr(npage+50))
   allocate(bufo(nbuf))
   allocate(bufn(nbuf))

   if (imode(3) .lt. 0) then
     write (nsyso,2301) 
2301 format (/,1x, 'NHEAT entry ', /)
   endif

   awfac=1/awr
   aw1fac=1/(awr+1)
   npkk=npk
   if (kchk.eq.1) npkk=3*npk-2
   if (mgam.gt.0) npkk=npkk+3
   if (mgam.eq.0.and.i6.gt.0) npkk=npkk+1
   npkkk=0
   do ipk=3,npk
      if (mtp(ipk).eq.442) npkkk=ipk
   enddo
   ! check for requested dame energy calculation
   idame=0
   ipk=2
   do while (idame.eq.0.and.ipk.lt.npk)
      ipk=ipk+1
      if (mtp(ipk).ge.444) idame=1
   enddo
   n6=0
   j6=0
   i6g=0
   irec=0
   jrec=0
   last6=0
   new6=0

   !--loop over non-redundant mt-s in file 3.
  105 continue

   if (imode(3) .lt. 0) then
     write (nsyso,3402) irec
3402 format (1x, ' nheat-105 ', i6)
   endif

   call contio(nin,0,0,scr,nb,nw)
   if (mfh.eq.0) go to 300
   if (mth.eq.0) go to 105
   if (mth.eq.3.or.mth.eq.4) go to 110
   if (mth.eq.10) go to 110
   if (mth.eq.26.or.mth.eq.27) go to 110
   if (mth.eq.18.and.mt19.eq.1) go to 110
   if (mth.ge.19.and.mth.le.21.and.mt19.eq.2) go to 110
   if (mth.eq.38.and.mt19.eq.2) go to 110
   if (mth.eq.101) go to 110
   if (mth.eq.103.and.mt103.gt.0) go to 110
   if (mth.eq.104.and.mt104.gt.0) go to 110
   if (mth.eq.105.and.mt105.gt.0) go to 110
   if (mth.eq.106.and.mt106.gt.0) go to 110
   if (mth.eq.107.and.mt107.gt.0) go to 110
   if (mth.eq.16.and.mt16.gt.0) go to 110
   if (iverf.le.5.and.mth.eq.719) go to 110
   if (iverf.le.5.and.mth.eq.739) go to 110
   if (iverf.le.5.and.mth.eq.759) go to 110
   if (iverf.le.5.and.mth.eq.779) go to 110
   if (iverf.le.5.and.mth.eq.799) go to 110
   if (mth.gt.120.and.mth.lt.152) go to 110
   if (mth.gt.200.and.mth.lt.600) go to 110
   zat=c1h
   awrt=c2h
   izat=nint(zat)
   go to 120
  110 continue
   call tosend(nin,0,0,scr)
   go to 105

   !--compute heating from this section.
  120 continue
   e=0
   call gety1(e,thresh,idis,y,nin,b)
   thresh=sigfig(thresh,7,0)
   t=c1h
   q=c2h
   qendf=q
   lr=l2h
   qs=0
   q0=0
   h0=0
   icon=0
   ebar=0
   dame=0
   mtd=mth

   !--check for an endf-6 file 6
   if (j6.gt.0) go to 123
   do i=1,i6
      if (mtd.eq.mt6(i)) go to 121
   enddo
   go to 124
  121 continue
   call findf(matd,6,mtd,nend6)
   call contio(nend6,0,0,scr,nb,nw)
   n6=nint(scr(5))
   lct=nint(scr(4))
   jp=nint(scr(3))
   jpn=mod(jp,10)
   jpp=jp-10*jpn
   if (mth.eq.18 .and. jpn.ne.0) then
      call tosend(nend6,0,0,scr)
      go to 105
   endif
   i6g=i6g+1
  123 continue
   j6=j6+1
   if (irec.gt.0.and.mt6no(i6g).gt.0) j6=j6-1
   last6=0
   if (j6.eq.n6.and.irec.gt.0) last6=1
   if (j6.eq.n6.and.mt6no(i6g).eq.0) last6=1
   new6=0
   if (j6.eq.1.and.irec.eq.0) new6=1
  124 continue

   !--adjust fission q to e(in)=0 prompt value
   if ((mth.ge.18.and.mth.le.21).or.mth.eq.38) then
      if (mt458.eq.1) then
         qendf=c458(15)
         q=qendf-qdel
         write(strng1,&
              &'(''changed q from '',1p,e14.6,'' to '',1p,e14.6)')&
              qendf,q
         write(strng2,'(''for mt '',i3,'' by taking out delayed components'')') mth
         call mess('nheat',strng1,strng2)
      else
         q=qendf
         cpoly(0)=q
      endif
   endif

   !--choose appropriate ground state q value.
   !--initialize heating routines.
   !--define neutron yield for this reaction.
   iimt=0
   if (nqa.eq.0) go to 140
   do 130 i=1,nqa
      if (mta(i).ne.mth) go to 130
      qs=qa(i)
      iimt=i
      go to 140
  130 continue
  140 continue
   if (j6.gt.0) go to 175
   if (mth.gt.15.and.mth.lt.51) go to 150
   if (mth.gt.100) go to 170
   icon=1
   q0=0
   if (lr.ne.0.and.lr.ne.31) q0=t
   if (iimt.ne.0) q0=qs
   if (mth.eq.91) go to 160
   yld=1
   ! accomodate neutron producing pseudo-levels
   if (lr.eq.16.or.lr.eq.21.or.lr.eq.24.or.lr.eq.26.or.lr.eq.30) yld=2
   if (lr.eq.17.or.lr.eq.25.or.lr.eq.38) yld=3
   if (lr.eq.37) yld=4
   nwa=na
   call disbar(e,ebar,dame,q,za,awr,nend4,matd,mtd,a,nwa)
   go to 180
  150 continue
   q0=q
   if (iimt.ne.0) q0=qs
  160 continue
   icon=2
   nwm=nwmax
   nwa=na
   call conbar(e,ebar,dame,nendf,nscr,matd,mtd,a,nwa,d,nwm,yld,q)
   go to 180
  170 continue
   if (iverf.le.5.and.mth.ge.700.and.mth.le.799) go to 171
   if (iverf.ge.6.and.mth.ge.600.and.mth.le.849) go to 171
   q0=q
   if (idame.gt.0) call capdam(e,dame,q,za,awr,mtd)
   go to 180
  171 continue
   if (iverf.le.5.and.mth.eq.718) go to 174
   if (iverf.le.5.and.mth.eq.738) go to 174
   if (iverf.le.5.and.mth.eq.758) go to 174
   if (iverf.le.5.and.mth.eq.778) go to 174
   if (iverf.le.5.and.mth.eq.798) go to 174
   if (iverf.ge.6.and.mth.eq.649) go to 174
   if (iverf.ge.6.and.mth.eq.699) go to 174
   if (iverf.ge.6.and.mth.eq.749) go to 174
   if (iverf.ge.6.and.mth.eq.799) go to 174
   if (iverf.ge.6.and.mth.eq.849) go to 174
   ! discrete charged-particle levels
   icon=1
   if (iverf.le.5.and.mth.eq.700) qsave=q
   if (iverf.le.5.and.mth.eq.720) qsave=q
   if (iverf.le.5.and.mth.eq.740) qsave=q
   if (iverf.le.5.and.mth.eq.760) qsave=q
   if (iverf.le.5.and.mth.eq.780) qsave=q
   if (iverf.le.5) q0=qsave
   if (iverf.ge.6) q0=t
   yld=1
   nwa=na
   call disbar(e,ebar,dame,q,za,awr,nend4,matd,mtd,a,nwa)
   go to 180
   ! continuum charged-particle emission
  174 continue
   yld=1
   if (iverf.le.5) q0=qsave
   if (iverf.ge.6) q0=t
   if (idame.gt.0) call capdam(e,dame,q,za,awr,mtd)
   go to 180
   ! file 6 reactions
  175 continue
   nwm=nwmax
   nwa=na
   call sixbar(e,ebar,yld,dame,nend6,a,nwa,nscr,d,nwm,n6,j6,irec,jrec,iflag)
   if (ebar.lt.zero) go to 291
   icon=-1
   q0=t
   if (iimt.ne.0) q0=qs
   if ((mtd.ge.18.and.mtd.le.21).or.mtd.eq.38) q0=q
  180 continue
   nmt=0
   if (npk.gt.2) call indx(npk,mtp,mt19,mtd,imt,nmt)
   if (icon.lt.0) go to 179
   if (iprint.ne.1) go to 182
   if (qs.ge.qtest) go to 181
   if ((mtd.ge.18.and.mtd.le.21).or.mtd.eq.38) go to 183
   if (icon.lt.0) go to 179
   if (idame.eq.0) write(nsyso,'(/&
     &'' neutron heating for mt'',i3,''   q0 ='',1p,e12.4,5x,&
     &''q ='',e12.4/14x,''e'',10x,''ebar'',9x,''yield'',10x,&
     &''xsec'',7x,''heating'')') mtd,q0,qendf
   if (idame.gt.0) write(nsyso,'(/&
     &'' neutron heating for mt'',i3,''   q0 ='',1p,e12.4,5x,&
     &''q ='',e12.4/14x,''e'',10x,''ebar'',9x,''yield'',10x,&
     &''xsec'',7x,''heating'',8x,''damage'')') mtd,q0,qendf
   go to 182
  183 continue
   if (idame.eq.0) write(nsyso,'(/&
     &'' neutron heating for mt'',i3,&
     &''   q0= variable fission q''/&
     &14x,''e'',12x,''q0'',10x,''ebar'',9x,''yield'',10x,&
     &''xsec'',7x,''heating'')') mtd
   if (idame.gt.0) write(nsyso,'(/&
     &'' neutron heating for mt'',i3,&
     &''   q0= variable fission q''/&
     &14x,''e'',12x,''q0'',10x,''ebar'',9x,''yield'',10x,&
     &''xsec'',7x,''heating'',8x,''damage'')') mtd
   if (irec.gt.0.and.icon.lt.0) izap=100
   go to 182
  181 continue
   if (idame.eq.0) write(nsyso,'(/&
     &'' neutron heating for mt'',i3,''   q0= variable'',5x,&
     &''q ='',e12.4/14x,''e'',10x,''qbar'',10x,''ebar'',9x,&
     &''yield'',10x,''xsec'',7x,''heating'')') mtd,qendf
   if (idame.gt.0) write(nsyso,'(/&
     &'' neutron heating for mt'',i3,''   q0= variable'',5x,&
     &''q ='',e12.4/14x,''e'',10x,''qbar'',10x,''ebar'',9x,&
     &''yield'',10x,''xsec'',7x,''heating'',8x,''damage'')')&
     &mtd,qendf
   go to 182
  179 continue
   izap=nint(zap)
   if (irec.gt.0) then
      if (mod(izat,1000).eq.0) then
         izap=1000*(izat/1000-izap/1000)
      else
         izap=nint(zat+1-zap)
      endif
   endif
   if (iprint.ne.1) go to 182
   if (qs.ge.qtest) go to 178
   if (idame.eq.0) write(nsyso,'(/&
     &'' file six heating for mt'',i3,'', particle ='',i6,&
     &5x,''q = '',1p,e12.4/14x,''e'',10x,''ebar'',9x,&
     &''yield'',10x,''xsec'',7x,''heating'')') mtd,izap,q0
   if (idame.gt.0) write(nsyso,'(/&
     &'' file six heating for mt'',i3,'', particle ='',i6,&
     &5x,''q = '',1p,e12.4/14x,''e'',10x,''ebar'',9x,&
     &''yield'',10x,''xsec'',7x,''heating'',8x,''damage'')')&
     mtd,izap,q0
   go to 182
  178 continue
   if (idame.eq.0) write(nsyso,'(/&
     &'' file six heating for mt'',i3,'', particle ='',i6,&
     &5x,''q0 = variable''/&
     &14x,''e'',10x,''qbar'',10x,''ebar'',9x,''yield'',10x,&
     &''xsec'',7x,''heating'')') mtd,izap
   if (idame.gt.0) write(nsyso,'(/&
     &'' file six heating for mt'',i3,'', particle ='',i6,&
     &5x,''q0 = variable''/&
     &14x,''e'',10x,''qbar'',10x,''ebar'',9x,''yield'',10x,&
     &''xsec'',7x,''heating'',8x,''damage'')') mtd,izap
  182 continue

   if (imode(3) .lt. 0 .and. mtd.eq. 102 ) then
     write (nsyso,3309) icon, irec
3309 format (/,1x, 'file six nheat mt=102, title write ', 2i6,/)
   endif

!
!  Temporary code fix to reset irec and correct multiplicity for MF6 MT102 recoils [PJG 9/14/2020]
!  No - this puts it into an infinite loop
!
!   irec = 0

   ipx=2
   irx=1
   ipfr=2
   irfr=1
   ipnp=2
   irnp=1
   ipgp=2
   irgp=1
   elst=0
   ilist=1
   if (icon.eq.0) yld=0
   yld0=0

   !--compute heating for energies in union grid.
   nn=0
  190 continue
   nn=nn+1
   call finda(nn,c,npkk,iold,bufo,nbuf)
   e=c(1)
   e=sigfig(e,9,0)
   if (new6.eq.1) c(npkk)=0

   if (imode(3) .lt. 0 .and. e .le. 1.01E-5 ) then
     write (nsyso,7309) e, mth, mtd, mfh, matd
7309 format (/,1x, 'file six nheat heating loop ', g14.6, 6i6,/)
   endif

   if (mtd.ne.2) go to 193
   if (e.lt.elst.and.nn.lt.ne) go to 193
   if (ilist.gt.ilmax) call error('nheat','storage exceeded.',' ')
   elist(ilist)=e
   elist(ilist+1)=big
   if (e.lt.break1) then
      iii=int(log10(e*(1+small))+small)
      if (log10(e*(1+small)).lt.zero) iii=iii-1
      xxx=ten**iii
      if (2*xxx.gt.e*(1+small).and.e.ge.one*(1-small)) then
         elst=2*xxx
      else if (5*xxx.gt.e*(1+small).and.e.ge.one*(1-small)) then
         elst=5*xxx
      else
         elst=10*xxx
      endif
   else
      test=break4
      test=sigfig(test,9,1)
      if (e.ge.test) elst=step4*int(e/step4+up)
      if (e.lt.test) elst=step3*int(e/step3+up)
      test=break3
      test=sigfig(test,9,1)
      if (e.lt.test) elst=step2*int(e/step2+up)
      test=break2
      test=sigfig(test,9,1)
      if (e.lt.test) elst=step1*int(e/step1+up)
   endif
   elst=elst-elst/10000000
  193 continue

   if (imode(3) .lt. 0 .and. mtd .eq. 102 ) then
     write (nsyso,3307) e, thresh, icon
3307 format (/,1x, 'nheat heating increment-B thresh', 2g14.7, i6)
   endif

   if (e.lt.thresh) go to 290
   call gety1(e,enext,idis,y,nin,b)
   ! check for energy-dependent q
   if (iimt.gt.0) then
      if (qa(iimt).ge.qtest) then
         lq0=lqs(iimt)
         call terpa(q0,e,enx,idx,qbar(lq0),ipx,irx)
         q=q0
      endif
   endif

   if (imode(3) .lt. 0) then
     write (nsyso,3306) e, ebar, dame, za, awr, icon
3306 format (1x, 'NHEAT location check ', 5g14.7, i6)
   endif

   if (icon.eq.1) call disbar(e,ebar,dame,q,za,awr,&
     nend4,matd,mtd,a,nwa)
   if (icon.eq.2) call conbar(e,ebar,dame,nendf,nscr,matd,mtd,&
     a,nwa,d,nwm,yld,q)
   if (icon.eq.0.and.idame.gt.0) call capdam(e,dame,q,za,awr,mtd)
   if (icon.lt.0)&
     call sixbar(e,ebar,yld,dame,nend6,a,nwa,nscr,d,nwm,n6,j6,irec,jrec,iflag)
   if (yld0.eq.zero) yld0=yld
   ! correct fission q for current incident energy
   if ((mtd.ge.18.and.mtd.le.21).or.mtd.eq.38) then
      if (lfc.eq.0) then     ! thermal point or polynomial
         if (nply.eq.0) then
            q0=cpoly(0)-(yld-yld0)*fq1+fq2*e
            if (mt458.eq.1) then
               h0=hpoly(0)
            else
               h0=q0-ebar*yld
            endif
         else
            q0=cpoly(nply)
            h0=hpoly(nply)
            do i=nply-1,0,-1
               q0=q0*e+cpoly(i)
               h0=h0*e+hpoly(i)
            enddo
         endif
         pnue=q0-h0
      else                   ! tabulated
         if (etabmax.lt.e.and.((abs(etabmax-e)/etabmax).gt.tol)) then
            write(strng,'(''upper energy tabulated fission q components'',&
               &'' is too low.'')')
            call error('nheat',strng,' ')
         endif
         if (ifc1.eq.1) then
            call terpa(qfr,e,enx,idx,afr,ipfr,irfr)
         else
            qfr=c458(1)
         endif
         if (ifc2.eq.1) then
            call terpa(qnp,e,enx,idx,anp,ipnp,irnp)
         else
            qnp=c458(3)-(yld-yld0)*fq1+fq2*e
         endif
         if (ifc4.eq.1) then
            call terpa(qgp,e,enx,idx,agp,ipgp,irgp)
         else
            qgp=c458(7)
         endif
         pnue=qnp
         h0=qfr+qgp
         q0=h0+pnue
      endif
      q0=q0-e
   endif
   if (icon.ge.0) then
      ebal6=0
      dame=dame*y
      if ((mtd.lt.18.or.mtd.gt.21).and.mtd.ne.38) then
         h=(e+q0-ebar*yld)*y
      else
         h=h0*y
      endif

   if (imode(3) .lt. 0 .and. mtd .eq. 102 .and. e .lt. 1.e-4) then
     write (nsyso,4303) e, h, ebar, yld, y, q0, h0, dame, icon
4303 format (/,1x, 'nheat heating increment-C ', 8g14.7, i6)
   endif

   else
      if ((mtd.lt.18.or.mtd.gt.21).and.mtd.ne.38) then
         h=ebar*yld*y
      else
         h=h0*y
      endif

   if (imode(3) .lt. 0 .and. mtd .eq. 102 .and. e .lt. 1.e-4) then
     write (nsyso,2303) e, h, ebar, yld, y, dame, icon, npkk
2303 format (/,1x, 'nheat heating increment-A ', 6g14.7, 2i6)
   elseif (imode(3) .lt. 0 .and. mtd .eq. 102 .and. e .gt. 9.6E6 .and. e .lt. 10.2E6) then
     write (nsyso,2303) e, h, ebar, yld, y, dame, icon, npkk
   endif

      dame=dame*y
      c(npkk)=c(npkk)+h
      if (izap.eq.0) c(npkk-1)=c(npkk-1)+h
      if (izap.eq.0.and.npkkk.gt.0) c(npkkk)=c(npkkk)+h
      if (izap.eq.1) h=0
      if (izap.eq.0.and.local.eq.0) h=0
      if (izap.le.1) dame=0
      ebal6=0
      if (last6.eq.1.and.mth.ne.5.and.mth.ne.102)&
        ebal6=(e+q0)*y-c(npkk)
      if ((mtd.ge.18.and.mtd.le.21).or.mtd.eq.38) ebal6=0
   endif
   if (mtd.ge.46.and.mtd.le.49) then
      ! treat second neutron from sequential (n,2n), mt=46-49.
      if (qs.ge.zero) call error('nheat',&
        'binding energy for sequential n2n needed.',' ')
      h=(qs-ebar)*y
   endif

   !--accumulate total and partial heating factors.
   c(2)=c(2)+h+ebal6
   if ((mtd.ge.18.and.mtd.le.21).or.mtd.eq.38) ebar=pnue/yld
   ebarp=sigfig(ebar,9,0)
   yldp=sigfig(yld,9,0)
   yp=sigfig(y,9,0)
   hp=sigfig(h,9,0)
   damep=sigfig(dame,9,0)

   if (imode(3) .lt. 0 .and. mtd .eq. 102 .and. e .lt. 1.e-4) then
     write (nsyso,2302) e, yp, hp, damep, c(2), h, ebal6, c(npkk), c(npkk-1), c(npkkk), npkk, npkkk
2302 format (1x, 'file six nheat mt=102, e <= 1.e-5 accumulated heating ', 10g14.7, 2i6)
   elseif (imode(3) .lt. 0 .and. mtd .eq. 102 .and. e .gt. 9.6e6 .and. e .lt. 10.2E6) then
     write (nsyso,2302) e, yp, hp, damep, c(2), h, ebal6, c(npkk), c(npkk-1), c(npkkk), npkk, npkkk
   endif

   ebal6p=sigfig(ebal6,9,0)
   if (iprint.ne.1.or.abs(e-elist(ilist)).gt.small*e&
     .or.y.le.smin) go to 199
   if (qs.ge.qtest) go to 198
   if ((mtd.ge.18.and.mtd.le.21).or.mtd.eq.38) go to 198
   if (idame.eq.0) then
      write(nsyso,'(1x,1p,7e14.4)') e,ebarp,yldp,yp,hp
   else
      write(nsyso,'(1x,1p,7e14.4)') e,ebarp,yldp,yp,hp,damep
   endif

   if (imode(3) .lt. 0 .and. mtd.eq. 102 ) then
     write (nsyso,3308)  e, ebarp, yldp, yp, hp
3308 format (1x, 'file six nheat mt=102, before 198 ', 5g14.7)
   endif

   go to 199
  198 continue

   if (imode(3) .lt. 0 .and. mtd .eq. 102 .and. e .le. 1.e-5) then
     write (nsyso,3302) e, yp, hp, damep
3302 format (1x, 'file six nheat mt=102-B, e <= 1.e-5 accumulated heating ', 4g14.7)
   elseif (imode(3) .lt. 0 .and. mtd .eq. 102 .and. e .gt. 9.6e6 .and. e .lt. 10.2e6) then
     write (nsyso,3302) e, yp, hp, damep
   endif

   if (idame.eq.0) then
      write(nsyso,'(1x,1p,7e14.4)') e,q0,ebarp,yldp,yp,hp
   else
      write(nsyso,'(1x,1p,7e14.4)') e,q0,ebarp,yldp,yp,hp,damep
   endif
  199 continue

   if (abs(ebal6).ge.100*y) then
      if (iprint.eq.1.and.abs(e-elist(ilist)).le.small*e&
        .and.y.gt.smin) then
         if (qs.lt.qtest) then
            write(nsyso,'(53x,''ebal'',1p,e14.4)') ebal6p
         else
            write(nsyso,'(67x,''ebal'',1p,e14.4)') ebal6p
         endif
      endif
   endif
   if (nmt.ne.0) then
      ! partial kermas
      do ii=1,nmt
         indxx=imt(ii)

   if (imode(3) .lt. 0 .and. e .le. 1.e-5) then
     write (nsyso,8302) e, c(indxx), dame, indxx
8302 format (1x, 'file six nheat mt=447-A:  ', 3g14.7, 2x, i6)
   elseif (imode(3) .lt. 0 .and. e .gt. 9.6e6 .and. e .lt. 10.2e6) then
     write (nsyso,8302) e, c(indxx), dame, indxx
   endif

         if (mtp(indxx).ge.444) c(indxx)=c(indxx)+dame
         if (mtp(indxx).lt.442) c(indxx)=c(indxx)+h+ebal6

   if (imode(3) .lt. 0 .and. e .le. 1.e-5) then
     write (nsyso,8572) e, c(indxx), dame, h, ebal6, indxx
8572 format (1x, 'file six nheat look ', 5g14.7, 2x, i6)
   elseif (imode(3) .lt. 0 .and. e .gt. 9.6e6 .and. e .lt. 10.2e6) then
     write (nsyso,8572) e, c(indxx), dame, h, ebal6, indxx
   endif

      enddo
   endif

   !--compute kinematic limits on heating, if desired.
   if (kchk.eq.0) go to 290
   if (icon.lt.0) go to 275
   if (mtd.ne.2) go to 215
   hmin=h
   hmax=hmin
   go to 285
  215 continue
   if (mtd.ne.102) go to 220
   hmin=y*e*aw1fac
   tt=q+awr*e*aw1fac
   hmax=hmin+y*tt*tt*rtm/2

   if (imode(3) .lt. 0 .and. e .le. 1.e-5) then
     write (nsyso,8372) e, hmax, hmin, y, aw1fac, awr, q, tt, rtm
8372 format (1x, 'nheat capture heating limits: ', 9g14.7, 2x, i6)
   elseif (imode(3) .lt. 0 .and. e .gt. 9.5E6 .and. e .lt. 10.1E6) then
     write (nsyso,8372) e, hmax, hmin, y, aw1fac, awr, q, tt, rtm
   endif

   go to 285
  220 continue
   if (mtd.ge.51.and.mtd.lt.91) go to 221
   if (iverf.le.5.and.mtd.ge.700.and.mtd.lt.718) go to 221
   if (iverf.le.5.and.mtd.ge.720.and.mtd.lt.738) go to 221
   if (iverf.le.5.and.mtd.ge.740.and.mtd.lt.758) go to 221
   if (iverf.le.5.and.mtd.ge.760.and.mtd.lt.778) go to 221
   if (iverf.le.5.and.mtd.ge.780.and.mtd.lt.798) go to 221
   if (iverf.ge.6.and.mtd.ge.600.and.mtd.lt.649) go to 221
   if (iverf.ge.6.and.mtd.ge.650.and.mtd.lt.699) go to 221
   if (iverf.ge.6.and.mtd.ge.700.and.mtd.lt.749) go to 221
   if (iverf.ge.6.and.mtd.ge.750.and.mtd.lt.799) go to 221
   if (iverf.ge.6.and.mtd.ge.800.and.mtd.lt.849) go to 221
   go to 230
  221 continue
   hmin=(e+q-ebar*yld)*y
   hmin=(e+q-ebar*yld)*y
   hmax=hmin
   if (lr.ne.0.and.lr.ne.31) hmax=(e+q0-ebar*yld)*y
   go to 285
  230 continue
   if (mtd.ne.91) go to 240
   hmin=y*(e+ebar)*awfac
   hmax=hmin
   if (lr.ne.0.and.lr.ne.31) hmax=(e+q0-ebar)*y
   go to 285
  240 continue
   if (mtd.lt.103) go to 250
   hmin=y*e*awfac
   hmax=h
   go to 285
  250 continue
   if (mtd.ne.18.and.mtd.ne.20.and.mtd.ne.21.and.mtd.ne.38) go to 260
   hmin=y*(e+q0-ebar*yld/2-efis)
   hmax=h
   go to 285
  260 continue
   test=10
   if (mtd.ne.16.or.awr.lt.test) go to 265
   hmin=0
   hmax=y*(e+ebar)/(awr-1)
   if (hmax.gt.h) hmax=h
   go to 285
  265 continue
   test=10
   if (mtd.ne.17.or.awr.lt.test) go to 270
   hmin=0
   hmax=y*(e+2*ebar)/(awr-2)
   if (hmax.gt.h) hmax=h
   go to 285
  270 continue
   hmin=0
   hmax=h
   go to 285
  275 continue
   hmin=h
   hmax=h
  285 continue
   npkt=npk-1
   if (kchk.ne.1) go to 283
   c(2+npkt)=c(2+npkt)+hmin
   c(2+2*npkt)=c(2+2*npkt)+hmax
  283 continue
   if (nmt.eq.0) go to 290
   do 286 ii=1,nmt
   indxx=imt(ii)
   if (mtp(indxx).eq.442) go to 286
   if (mtp(indxx).eq.443) go to 284
   if (mtp(indxx).ge.444) go to 286
   if (kchk.ne.1) go to 286
   c(indxx+npkt)=c(indxx+npkt)+hmin
   c(indxx+2*npkt)=c(indxx+2*npkt)+hmax
   go to 286
  284 continue
   c(indxx)=c(indxx)+hmax

   if (imode(3) .lt. 0) then
     write (nsyso,3319) c(indxx), c(indxx+2*npkt), c(indxx+npkt), hmax, hmin, indxx , npkt
3319 format (1x, ' c-look-1 ', 5g14.7, 2i9)
   endif

  286 continue

   !--store results on scratch and go to next energy.
  290 continue
   intt=nn
   if (intt.eq.ne) intt=-intt
   call loada(intt,c,npkk,inew,bufn,nbuf)
   if (e.eq.elist(ilist)) ilist=ilist+1
   if (nn.lt.ne) go to 190

   !--for file 6, loop over subsections
   isave=iold
   iold=inew
   inew=isave
  291 continue
   if (j6.ge.n6) then
      ! make sure we're at the end of section (will not be so if
      ! there was no particle distribution for the last subsection
      ! in this section).
      if (ebar.lt.zero) then
         call tosend(nin,0,0,b)
         call skiprz(nin,-1)
      endif
      go to 296
   endif
   call skiprz(nin,-2)
   call findf(matd,3,mtd,nin)
   go to 297
  296 continue
   if (j6.eq.0) go to 297
   if (irec.gt.0.and.mt6no(i6g).gt.0) go to 293
   if (iprint.eq.1.and.i6g.gt.0.and.izap.ne.0.and.iflag.eq.0)&
     write(nsyso,'(/''   no explicit file 6 photon production for mt'',i3)')&
     mtd
   if (mt6no(i6g).ne.0) then
      if (mth.ne.102 ) write(nsyso,'(/&
        &''   generating recoil with one-particle approx.'')')
      if (mth.eq.102 ) write(nsyso,'(/&
        &''   generating recoil from photon momentum.'')')

   if (imode(3) .lt. 0) then
     write (nsyso,3301) irec, mt6no(i6g)
3301 format (1x, ' nheat enerating recoil ', 2i6)
   endif

      irec=mt6no(i6g)
      call skiprz(nin,-2)
      call findf(matd,3,mtd,nin)
      go to 297
   endif
  293 continue
   irec=0
   jrec=0
   j6=0
   n6=0
  297 continue
   go to 105

   !--nheat is finished.
  300 continue
   deallocate(b)
   deallocate(d)
   deallocate(a)
   deallocate(scr)
   deallocate(bufn)
   deallocate(bufo)
   call skiprz(nin,-4)
   return
   end subroutine nheat

   subroutine indx(npk,mtp,mt19,mt,imt,nmt)
   !-------------------------------------------------------------------
   ! Find the index of a given MT in the stored desired MT array.
   ! Heating identifiers
   !     300 + MT for reaction
   ! Special values allowed are
   !     303=nonelastic (all but MT2)
   !     304=inelastic (MTs1-91)
   !     318=fission (MT18 or MT19-21 and MT38)
   !     401=disappearance (MT102-120)
   !     442=total photon ev-barns in kerma
   !     443=kinematic kerma
   ! Damage energy identifiers
   !     444=total
   !     445=elastic
   !     446=inelastic (MT51-MT91)
   !     447=total capture (MT102-MT120)
   !-------------------------------------------------------------------
   use endf ! provides iverf
   ! externals
   integer::npk,mtp(*),mt19,mt,imt(*),nmt
   ! internals
   integer::i,mtpi,iflag

   nmt=0
   if (npk.lt.3) return
   do i=3,npk
      mtpi=mtp(i)
      iflag=0
      if (mtpi.eq.442) iflag=1
      if (mtpi.eq.443) iflag=1
      if (mtpi.eq.444) iflag=1
      if (mtpi.eq.445.and.mt.eq.2) iflag=1
      if (mtpi.eq.446.and.mt.ge.51.and.mt.le.91) iflag=1
      if (mtpi.eq.447.and.mt.ge.102.and.mt.le.120) iflag=1
      if (iverf.lt.6) then
         if (mtpi.eq.447.and.mt.ge.700.and.mt.le.799) iflag=1
      else
         if (mtpi.eq.447.and.mt.ge.600.and.mt.le.849) iflag=1
      endif
      if (mtpi.eq.303.and.mt.ne.2.and.mt.ne.3) iflag=1
      if (mtpi.eq.304.and.mt.ge.51.and.mt.le.91) iflag=1
      if (mtpi.eq.318) then
         if (mt19.eq.0.and.mt.eq.18) iflag=1
         if (mt19.eq.1.and.(mt.eq.19.or.&
           mt.eq.20.or.mt.eq.21.or.mt.eq.38)) iflag=1
         if (mt19.eq.2.and.mt.eq.18) iflag=1
      endif
      if (mtpi.eq.401) then
         if (mt.ge.102.and.mt.le.120) iflag=1
         if (iverf.lt.6.and.mt.ge.700.and.mt.lt.800) iflag=1
         if (iverf.ge.6.and.mt.ge.600.and.mt.lt.850) iflag=1
      endif
      if (iverf.lt.6) then
         if (mtpi.eq.403.and.mt.ge.700.and.mt.le.719) iflag=1
         if (mtpi.eq.404.and.mt.ge.720.and.mt.le.739) iflag=1
         if (mtpi.eq.405.and.mt.ge.740.and.mt.le.759) iflag=1
         if (mtpi.eq.406.and.mt.ge.760.and.mt.le.779) iflag=1
         if (mtpi.eq.407.and.mt.ge.780.and.mt.le.799) iflag=1
      else
         if (mtpi.eq.403.and.mt.ge.600.and.mt.le.649) iflag=1
         if (mtpi.eq.404.and.mt.ge.650.and.mt.le.699) iflag=1
         if (mtpi.eq.405.and.mt.ge.700.and.mt.le.749) iflag=1
         if (mtpi.eq.406.and.mt.ge.750.and.mt.le.799) iflag=1
         if (mtpi.eq.407.and.mt.ge.800.and.mt.le.849) iflag=1
      endif
      if (mtpi.eq.mt+300) iflag=1
      if (iflag.eq.1) then
         nmt=nmt+1
         imt(nmt)=i
      endif
   enddo
   return
   end subroutine indx

   subroutine capdam(ee,dame,q,za,awr,mtd)
   !-------------------------------------------------------------------
   ! Compute damage energy for capture reactions.
   !-------------------------------------------------------------------
   use endf    ! provides iverf,terp1
   use mainio  ! provides nsysi,nsyso
   use snl     ! provides SNL
   ! externals
   integer::mtd
   real(kr)::ee,dame,q,za,awr
   ! internals
   integer::iq,iz
   real(kr)::e,zx,ax,denom,z,en,aw1fac
   real(kr)::damn,daml,el,er,ea,ec,et
   real(kr):: er_low, er_high
   real(kr),dimension(4),parameter::qp=(/&
      -.86114e0_kr,-.33998e0_kr,.33998e0_kr,.86114e0_kr/)
   real(kr),dimension(4),parameter::qw=(/&
     .34785e0_kr,.65215e0_kr,.65215e0_kr,.34785e0_kr/)
   real(kr),parameter::third=.333333333e0_kr
   real(kr),parameter::econ=1.029e6_kr
   real(kr),parameter::eps=1.e-10_kr
   real(kr),parameter::step=1.1e0_kr
   integer,parameter::nq=4
   real(kr),parameter::zero=0
   save en,damn,el,daml
   save zx,ax,denom,z,aw1fac
   save er_low, er_high

   if (imode(3) .lt. 0) then
     write (nsyso,2301) ee, en, q, za, awr, mtd
2301 format (1x, 'CAPDAM entry-mod ', 5g14.7, i6)
   endif
   !--initialize
   e=ee
   if (e.eq.zero) then
      zx=0
      if (mtd.gt.102) zx=1
      if (mtd.eq.106.or.mtd.eq.107) zx=2
      if (iverf.le.5.and.(mtd.eq.778.or.mtd.eq.798)) zx=2
      ax=0
      if (mtd.gt.102) ax=1
      if (mtd.eq.104) ax=2
      if (iverf.le.5.and.mtd.eq.738) ax=2
      if (mtd.eq.105.or.mtd.eq.106) ax=3
      if (iverf.le.5.and.(mtd.eq.758.or.mtd.eq.778)) ax=3
      if (mtd.eq.107) ax=4
      if (iverf.le.5.and.mtd.eq.798) ax=4
      denom=1/awr**third
      if (ax.ne.zero) denom=1/(ax**third+awr**third)
      iz=nint(za/1000)
      z=iz
      en=0
      damn=df(en,z-zx,awr+1-ax,z,awr)
      dame=0
      aw1fac=1/(awr+1)

      er_low = 0
      er_high = 0

      return
   endif

   !--normal entry
   !--if energy not in current range,
   !--construct new interpolation range.
   if (ee.ge.en*(1-eps)) then
      el=en
      daml=damn

      er_low = er_high

      e=step*el
      if (e.lt.ee*(1-eps)) e=(1-eps)*ee
      dame=0
      if (mtd.eq.102) then
         ! capture followed by gamma emission
         er=e*aw1fac
         ea=awr*er+q

   if (imode(3) .lt. 0) then
     write (nsyso,9302) er, e, aw1fac, awr, q, ea
9302 format (1x, 'CAPDAM damage partial ', 8g14.7)
   endif

         er=er+ea*ea*rtm/2
         dame=df(er,z,awr+1,z,awr)

         er_high = er

   if (imode(3) .lt. 0) then
     write (nsyso,4302) e, er, dame, aw1fac, awr, rtm, ea, q
4302 format (1x, 'CAPDAM damage default increment ', 8g14.7)
   endif

      else if (ax.ne.zero) then
         ! capture followed by particle emission
         ec=econ*zx*z*denom
         ea=q+awr*e*aw1fac
         if (ea.ge.zero) then
            et=(awr+1-ax)*e*aw1fac
            if (ea.gt.ec*(1+eps)) ea=ec
            do iq=1,nq
               er=(et-2*sqrt(et*ax*ea)*qp(iq)+ax*ea)*aw1fac
               dame=dame+qw(iq)*df(er,z-zx,awr+1-ax,z,awr)/2

               if (imode(3) .lt. 0) then
                   write (nsyso,4301) iq, er, dame
4301               format (1x, 'CAPDAM damage increment ', i6, 2g14.7)
               endif

            enddo
         endif
      endif
      en=e
      damn=dame

      er_high = er

   endif

   !--interpolate for damage energy
   call terp1(el,daml,en,damn,ee,dame,2)

   !-- if new grid point not defined, e is not defined, so, interpolate for recoil energy
   if ( er .lt. 1.0E-10) then
       call terp1(el,er_low,en,er_high,ee,er,2)
   endif

   if (imode(3) .lt. 0) then
     write (nsyso,5301) ee, er, dame, nq
5301 format (1x, 'CAPDAM exit ', 3g14.7, i6)
   endif

   return
   end subroutine capdam

   subroutine disbar(ee,ebar,dame,q,za,awr,ntape,matd,mtd,a,nwa)
   !-------------------------------------------------------------------
   ! Compute ebar and damage energy for discrete scattering.
   !-------------------------------------------------------------------
   use endf ! provides iverf,terp1
   use snl     ! provides SNL
   use mainio  ! provide nsyso
   use mathm ! provides legndr
   use physics !get global physics and light particle mass constants
   ! externals
   integer::ntape,matd,mtd,nwa
   real(kr)::ee,ebar,dame,q,za,awr,a(*)
   real(kr) :: edis_replace
   ! internals
   integer::imiss,i,mfd,nld,lcd,idis,iz,iq,il
   real(kr)::e,thresh,awp,afact,arat,en,cn,el,cl,z
   real(kr)::daml,damn,enx,enext,u,f,e2,ce,r,b,g,test,wbar
   real(kr)::fl(65),p(65)
   integer,parameter::nq=64
   real(kr),dimension(64),parameter::qp=(/&
     -9.99305042E-01_kr,-9.96340117E-01_kr,-9.91013371E-01_kr,&
     -9.83336254E-01_kr,-9.73326828E-01_kr,-9.61008800E-01_kr,&
     -9.46411375E-01_kr,-9.29569172E-01_kr,-9.10522137E-01_kr,&
     -8.89315446E-01_kr,-8.65999398E-01_kr,-8.40629296E-01_kr,&
     -8.13265315E-01_kr,-7.83972359E-01_kr,-7.52819907E-01_kr,&
     -7.19881850E-01_kr,-6.85236313E-01_kr,-6.48965471E-01_kr,&
     -6.11155355E-01_kr,-5.71895646E-01_kr,-5.31279464E-01_kr,&
     -4.89403146E-01_kr,-4.46366017E-01_kr,-4.02270158E-01_kr,&
     -3.57220158E-01_kr,-3.11322872E-01_kr,-2.64687162E-01_kr,&
     -2.17423644E-01_kr,-1.69644420E-01_kr,-1.21462819E-01_kr,&
     -7.29931218E-02_kr,-2.43502927E-02_kr, 2.43502927E-02_kr,&
      7.29931218E-02_kr, 1.21462819E-01_kr, 1.69644420E-01_kr,&
      2.17423644E-01_kr, 2.64687162E-01_kr, 3.11322872E-01_kr,&
      3.57220158E-01_kr, 4.02270158E-01_kr, 4.46366017E-01_kr,&
      4.89403146E-01_kr, 5.31279464E-01_kr, 5.71895646E-01_kr,&
      6.11155355E-01_kr, 6.48965471E-01_kr, 6.85236313E-01_kr,&
      7.19881850E-01_kr, 7.52819907E-01_kr, 7.83972359E-01_kr,&
      8.13265315E-01_kr, 8.40629296E-01_kr, 8.65999398E-01_kr,&
      8.89315446E-01_kr, 9.10522137E-01_kr, 9.29569172E-01_kr,&
      9.46411375E-01_kr, 9.61008800E-01_kr, 9.73326828E-01_kr,&
      9.83336254E-01_kr, 9.91013371E-01_kr, 9.96340117E-01_kr,&
      9.99305042E-01_kr/)
   real(kr),dimension(64),parameter::qw=(/&
      1.78328072E-03_kr, 4.14703326E-03_kr, 6.50445797E-03_kr,&
      8.84675983E-03_kr, 1.11681395E-02_kr, 1.34630479E-02_kr,&
      1.57260305E-02_kr, 1.79517158E-02_kr, 2.01348232E-02_kr,&
      2.22701738E-02_kr, 2.43527026E-02_kr, 2.63774697E-02_kr,&
      2.83396726E-02_kr, 3.02346571E-02_kr, 3.20579284E-02_kr,&
      3.38051618E-02_kr, 3.54722133E-02_kr, 3.70551285E-02_kr,&
      3.85501532E-02_kr, 3.99537411E-02_kr, 4.12625632E-02_kr,&
      4.24735151E-02_kr, 4.35837245E-02_kr, 4.45905582E-02_kr,&
      4.54916279E-02_kr, 4.62847966E-02_kr, 4.69681828E-02_kr,&
      4.75401657E-02_kr, 4.79993886E-02_kr, 4.83447622E-02_kr,&
      4.85754674E-02_kr, 4.86909570E-02_kr, 4.86909570E-02_kr,&
      4.85754674E-02_kr, 4.83447622E-02_kr, 4.79993886E-02_kr,&
      4.75401657E-02_kr, 4.69681828E-02_kr, 4.62847966E-02_kr,&
      4.54916279E-02_kr, 4.45905582E-02_kr, 4.35837245E-02_kr,&
      4.24735151E-02_kr, 4.12625632E-02_kr, 3.99537411E-02_kr,&
      3.85501532E-02_kr, 3.70551285E-02_kr, 3.54722133E-02_kr,&
      3.38051618E-02_kr, 3.20579284E-02_kr, 3.02346571E-02_kr,&
      2.83396726E-02_kr, 2.63774697E-02_kr, 2.43527026E-02_kr,&
      2.22701738E-02_kr, 2.01348232E-02_kr, 1.79517158E-02_kr,&
      1.57260305E-02_kr, 1.34630479E-02_kr, 1.11681395E-02_kr,&
      8.84675983E-03_kr, 6.50445797E-03_kr, 4.14703326E-03_kr,&
      1.78328072E-03_kr/)
   real(kr),parameter::edis=25.e0_kr
   real(kr),parameter::step=1.1e0_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0
   save imiss,thresh,awp,afact,arat
   save en,cn,el,cl,z,daml,damn,enx,enext

   if (imode(3) .lt. -1) then
     write (nsyso,2301) ee,q,za,awr,matd,mtd
2301 format (/,1x, 'DISBAR entry ', 4g14.7, 2i6,/)
   endif

   !--initialize if e=0.
   e=ee
   if (e.eq.zero) then
      edis_replace = edis
      imiss=0
      do i=1,nmiss4
         if (miss4(i).eq.mtd) imiss=1
      enddo
      if (imiss.ne.1) then
         mfd=4
         nld=60
         lcd=2
         call hgtfle(e,enext,idis,fl,nld,lcd,matd,mfd,mtd,&
           ntape,a,nwa)
      else
         nld=1
         enext=etop
      endif
      thresh=((awr+1)/awr)*(-q)
      awp=1
      if (iverf.le.5.and.mtd.ge.700.and.mtd.lt.718) awp=pnratio
      if (iverf.le.5.and.mtd.ge.720.and.mtd.lt.738) awp=dnratio
      if (iverf.le.5.and.mtd.ge.740.and.mtd.lt.758) awp=tnratio
      if (iverf.le.5.and.mtd.ge.760.and.mtd.lt.778) awp=hnratio
      if (iverf.le.5.and.mtd.ge.780.and.mtd.lt.798) awp=anratio
      if (iverf.ge.6.and.mtd.ge.600.and.mtd.lt.649) awp=pnratio
      if (iverf.ge.6.and.mtd.ge.650.and.mtd.lt.699) awp=dnratio
      if (iverf.ge.6.and.mtd.ge.700.and.mtd.lt.749) awp=tnratio
      if (iverf.ge.6.and.mtd.ge.750.and.mtd.lt.799) awp=hnratio
      if (iverf.ge.6.and.mtd.ge.800.and.mtd.lt.849) awp=anratio
      afact=awp/((awr+1)**2)
      arat=awp/(awr+1-awp)
      en=0
      cn=0
      el=0
      cl=0
      enx=0
      if ( edis_replace .gt. break_new .and. SNL_enhanced_input_format .eq. 1) then 
         write (nsyso, 4581) edis_replace, break_new
 4581    format (/,1x, 'edis parameter in disbar over-ridden by displacement', &
  &          ' threshold energy: ', 2g14.7)
         edis_replace = break_new
      endif
      if (idame.ne.0) then
         iz=nint(za/1000)
         z=iz
         damn=df(e,z,awr,z,awr)
         daml=damn
         if (mtd.eq.2) then
            daml=0
!            enx=edis*arat/(4*afact)
            enx=edis_replace*arat/(4*afact)
!           write (nsyso, 2891) enx, edis_replace, arat, afact
 2891       format (1x, 'disbar enx set: ', 4g14.7)
            if ( icntrl(1) .eq. 7) then
               enx = 0.0
!              write (nsyso, 6714) enx
 6714          format (1x, 'disbar enx over-ride ', g14.7)
            endif
            damn=0
         endif
      endif
      return
   endif

   !--normal entry.
   if (ee.le.en*(1+small)) go to 130
   el=en
   cl=cn
   if (idame.ne.0) daml=damn
   e=step*el
   if (enext.lt.e*(1-small)) e=enext
   if (ee.gt.e*(1+small)) e=ee
   r=1
   if (mtd.ne.2) then
      if (thresh.ge.e*(1-small)) then
         r=0
      else
         r=sqrt(1-thresh/e)
      endif
   endif
   if (imiss.ne.1) then
      mfd=4
      nld=60
      lcd=2
      call hgtfle(e,enext,idis,fl,nld,lcd,matd,mfd,mtd,ntape,a,nwa)
   else
      nld=1
      fl(1)=1
      fl(2)=0
   endif
   en=e
   b=r*sqrt(awr/arat)
   g=b*arat
   cn=0
   test=1
   test=test/1000000
   if (abs(awp-1).gt.test) go to 120

   !--compute ebar
   wbar=fl(2)
   if (wbar.gt.qp(64)) wbar=qp(64)
   cn=(1+2*b*wbar+b*b)*afact
   if (idame.eq.0) go to 130
!  write (nsyso, 7814) e, enx 
 7814   format (1x, 'disbar angle intg start', 4g14.7, 2i6, g14.7 )
   if (e.lt.enx*(1-small)) go to 130

   !--angle integration by gauss-legendre method.
  120 continue
   dame=0
   do iq=1,nq
      u=qp(iq)
      call legndr(u,p,nld)
      f=0
      do il=1,nld
         f=f+(2*il-1)*fl(il)*p(il)/2
      enddo
      e2=e*(1-2*g*u+g*g)*afact/arat
!     if ( iq .eq. 1) write (nsyso, 89) e, iq, nq, e2, awp, awr, &
!     &   afact, arat,    &
!     &   4.0*afact/arat, (1-2*g*u+g*g), g, u
 89   format (1x, 'Discrete scatter recoil energy: ', g14.7, 1x, &
      & 2i5, 9g14.7)
      dame=dame+qw(iq)*f*df(e2,z,awr,z,awr)
   enddo
   if (dame.lt.zero) dame=zero
   damn=dame

   !--interpolate for results.
  130 continue
   call terp1(el,cl,en,cn,ee,ce,2)
   ebar=ee*ce
   if (idame.gt.0) then
      call terp1(el,daml,en,damn,ee,dame,2)
   endif
!   write (nsyso, 7813) damn, dame, ebar, ee, e, ce
 7813   format (1x, 'disbar exit', 6g14.7 )
   return
   end subroutine disbar

   real(kr) function df(e,zr_in,ar_in,zl_in,al_in)
   !-------------------------------------------------------------------
   ! Damage function using the Lindhard partition of
   ! energy between atomic and electronic motion.
   ! Call with e=0 for each reaction to precompute the constants.
   !-------------------------------------------------------------------
   use snl     ! provides SNL
   use mainio  ! provides nsysi,nsyso,nsyse
   ! externals
   real(kr):: zr_in, ar_in, zl_in, al_in
   real(kr)::e,zr,ar,zl,al
   real(kr) :: threshold_factor
   integer:: jloop_al, jloop_zl, jloop_de
   data jloop_al, jloop_zl, jloop_de /0, 0, 0/
   real(kr):: over_ride, critical_mass
   ! internals
!
!  break_i is the intermediate break point for the dpa model
!  e       is the energy of the recoil ion
! 
   real(kr):: break_i, dam_tabular, damage_energy, dam_effic
   real(kr)::el,rel,denom,fl,ep,dam
   real(kr),parameter::twothd=.666666667e0_kr
   real(kr),parameter::threeq=.75e0_kr
   real(kr),parameter::sixth=.166666667e0_kr
   real(kr),parameter::onep5=1.5e0_kr
   real(kr),parameter::c1=30.724e0_kr
   real(kr),parameter::c2=.0793e0_kr
   real(kr),parameter::c3=3.4008e0_kr
   real(kr),parameter::c4=.40244e0_kr
   real(kr),parameter::zero=0
   save rel,fl

   if (imode(3) .lt. -1) then
     write (nsyso,2301) e,zr_in,ar_in,zl_in,al_in
2301 format (/,1x, 'DF entry ', 5g14.7,/)
   endif

   zr = zr_in
   ar = ar_in
   zl = zl_in
   al = al_in
!  write (nsyso, 6712) zl, al, zr, ar, zl_new, al_new
 6712 format (1x,'DF entry: ', 6g14.7)
!
!  added logic for expaned damage functions
! 
   if ( icntrl(5) .ne. 0 ) then
!        over-ride displacement threshold if icntrl(5) set
      break = displace_th
!        notify user first time over-ride takes place
      if (jloop_de .le. 0) write (nsyso,2034) break
      if (jloop_de .le. 0) write (6,2034) break
 2034  format (1x, 'Threshold displacement energy replaced with ', g14.7)
       jloop_de = jloop_de + 1
   endif
   if (break .lt. 0.0) then
       write (nsyso,9034) break, displace_th
       write (6,9034) break, displace_th
 9034  format (1x, 'Invalid displace_th replaced ', 2g14.7)
       displace_th = 25.
       break = 25.0
   endif
   if ( icntrl(1) .eq. 7) then
!    Output displacement kerma in MT=444, i.e. over-ride Ed with 0.0 eV
       break = 0.0
   endif
   if ( al_new .gt. 0.0) then
       al = al_new
       jloop_al = jloop_al + 1
       if ( jloop_al .le. 1) write (6,3034) al
 3034  format (1x, 'Lattice atomic mass replaced with ', g14.7 ' nmass')
   endif
   if ( zl_new .gt. 0.0) then
       zl = zl_new
       jloop_zl = jloop_zl + 1
       if ( jloop_zl .le. 1) write (6,8034) zl
 8034  format (1x, 'Lattice atomic number replaced with ', g14.7)
   endif
   break_i = break
   if ( icntrl(1) .eq. 6) then
      break_i = 2.0*break
   else if (icntrl(1) .eq. 5) then
      break_i = 2.0*break/0.8
   else if ( icntrl(1) .eq. 8) then
      break_i = break
   else if ( icntrl(1) .eq.0) then
     break_i = 2.0*break/0.8
   endif
!
! end of added damage logic section
!
   if (zr.eq.zero) then
      df=0
      dam = 0.0
   else if (e.le.zero) then
      el=c1*zr*zl*sqrt(zr**twothd+zl**twothd)*(ar+al)/al
      rel=1/el
      denom=(zr**twothd+zl**twothd)**threeq*ar**onep5*sqrt(al)
      fl=c2*zr**twothd*sqrt(zl)*(ar+al)**onep5/denom
      df=0
      dam = 0.0
   else if (e.lt.break) then
      df=0
      dam = 0.0
   else
      el=c1*zr*zl*sqrt(zr**twothd+zl**twothd)*(ar+al)/al
      rel=1/el
      denom=(zr**twothd+zl**twothd)**threeq*ar**onep5*sqrt(al)
      fl=c2*zr**twothd*sqrt(zl)*(ar+al)**onep5/denom
      ep=e*rel
      dam=e/(1+fl*(c3*ep**sixth+c4*ep**threeq+ep))
      df=dam
!     
!     Over-ride damage partition function with user-supplied damage partition function
!     
      if ( icntrl(8) .eq.1) then
!         dam_tabular = 1.0
!        damage partition function interpolation - log (fraction) with log (PKA energy)
         icode = 5
         dam_tabular = fitmd(e, ndamage, energy_dam, value_dam, icode)
         if ( imode(3) .lt. -1) then
            write (nsyso, 7823) e, dam/e, dam_tabular, zr, ar, zl, al
 7823       format (1x, 'icntrl(8) damage partition ', 7g14.7)
         endif
         dam = dam_tabular*e
         df = dam
      endif
!  
!     Save the user-specified damage efficiency function on top of damage partition function
!       set default to unity; compute efficiency if user-flag is set
!
      dam_effic = 1.0
      if ( icntrl(4) .eq. 1) then 
!        damage efficiency interpolation - log (efficiency) with log (PKA energy)
         icode = 5
         if ( e .lt. eff_eng(1)) then 
            dam_effic = eff_value(1)
         elseif ( e .gt. eff_eng(npoints_eff)) then 
            dam_effic = eff_value(npoints_eff)
         else
            dam_effic = fitmd(e, npoints_eff, eff_eng, eff_value, icode)
         endif
         if ( imode(3) .lt. -1) then 
            write (nsyso, 3823) e, dam, dam_effic
 3823       format (1x, 'icntrl(4) damage efficiency ', 7g14.7)
         endif
      endif
!
!     Check on option to inhibit low recoil mass contributions to damage
!
      if ( icntrl(9) .gt. 0) then
         critical_mass = al - icntrl(9)*1.0
         if (ar .le. critical_mass) then
 !          turn-off damage energy from secondary charged particles
 !          with mass less "i" less than the lattice atom
            dam = 0.0
            df = 0.0
         endif
      endif
!
   endif
   damage_energy = dam
!
!  Add logic to address various damage energy interpretations
!
      if ( icntrl(1) .eq. 0) then
!
!          return the damage energy 
!  
        threshold_factor = 1.0
        df = dam
!
      elseif ( icntrl(1) .eq. 5) then
!
!          return the NRT damage energy model
!               three interval threshold model
!               transitions at Ed and 2*Ed/beta, where beta = 0.8
!
        if ( damage_energy .lt. break) then
           dam = 0.0
           df  = 0.0
           threshold_factor = 1.0
        elseif ( damage_energy .lt. break_i) then
           if ( dam .le. 0.0) then
             threshold_factor = 1.0
             df = 0.0
             dam = 0.0
           else
!
!            for DE corresponding to 1 dpa, 
!                remove damage energy and set to energy for one Frenkel pair, i.e. 2*Ed/beta
!
             threshold_factor = 2.0*break/0.8/dam
           endif
           df = threshold_factor*dam
           dam = df
        else 
           threshold_factor = 1.0
           df = threshold_factor*dam
           dam = df
        endif
!
      elseif ( icntrl(1) .eq. 6) then
!
!       Kinchin-Pease damage energy model
!           three interval threshold model
!           transitions at Ed and 2*Ed
!
        if ( damage_energy .lt. break) then
           dam = 0.0
           df  = 0.0
           threshold_factor = 1.0
        elseif ( damage_energy .lt. 2.0*break) then
           if ( dam .le. 0.0) then
             threshold_factor = 0.0
           else
             threshold_factor = 2.0*break/dam
           endif
           df = threshold_factor*dam
        else
           threshold_factor = 1.0
           df = threshold_factor*dam
        endif
!
      elseif (icntrl(1) .eq. 7) then
!
!       Displacement Kerma as output metric
!           no threshold model
!           full non-ionizing damage is scored for all enegies
!
        threshold_factor = 1.0
        if ( break .gt. 0.0) then
          write (nsyso, 3891) break, dam, icntrl(1)
 3891     format (1x, 'ERROR: displacement kerma requires zero ', &
  &          'displacement threshold energy ', 2g14.7, i5)
        endif
        df = dam
!        write (nsyso, 3441) e, dam
 3441   format (1x, 'df icntrl(1)=7 set: ', g14.7)
!
      elseif (icntrl(1) .eq. 8) then
!
!       Sharp transition Kinchin-Pease damage energy model
!           two interval threshold model
!           transitions at Ed
!
        if ( damage_energy .lt. break) then
           dam = 0.0
           df  = 0.0
           threshold_factor = 1.0
        else
           threshold_factor = 1.0
           df = threshold_factor*dam
        endif
      endif
!
!    Add over-ride logic to produce a dpa rather than a damage energy
!
      if ( icntrl(1) .eq. 6 .and. icntrl(10) .eq. 1 .and. break .gt. 0.0) then
!        conversion from damage energy to Kinchin-Pease dpa
         df = dam/(2.0*break)
      elseif ( icntrl(1) .eq. 5 .and. icntrl(10) .eq. 1 .and. break .gt. 0.0) then
!        conversion from damage energy to NRT dpa
         df = 0.8*dam/(2.0*break)
      elseif ( icntrl(1) .eq. 8 .and. icntrl(10) .eq. 1 .and. break .gt.0.0) then
         df = dam/(2.0*break)
      endif
      if ( icntrl(8) .eq. 1) then
         over_ride = 1.0
!         if( jloop_de .lt. 20) write (nsyso, 8001) jloop_de, &
!     &           e, dam, df, over_ride, break
 8001    format (1x, 'Function df icntrl(8) process: ', i5, 1x, 5g14.7)
      endif
      if ( imode(3) .le. -3) then
        write (nsyso, 6711) e, df, dam, dam/e, break, break_i, zl, &
  &                  al, zr, ar, threshold_factor, icntrl(1)
 6711   format (1x, 'DF exit debug: ', 11g14.7, 2x, i5)
      endif
!
!     Now, apply the damage efficiency, if it was selected, on top of the treshold and dpa treatment
!
      df = df*dam_effic
!
   return
   end function df

   real(kr) function fitmd( x, nx, xarray, yarray,icode)
   use endf    ! provides iverf,terp1
   real(kr):: x
   integer:: nx, icode,idirection,ipick,ipick1,i
   real(kr):: xarray(nx), yarray(nx)
   real(kr):: x1, x2, y1, y2, y,place
   real(kr):: log, exp
!
!  determine is xarray is increasing or deceasing
!
   ipick = 1
 101  continue
   if ( xarray(ipick+1) .gt. xarray(ipick)) then
         idirection = 1
   elseif ( xarray(ipick+1) .lt. xarray(ipick)) then
         idirection = -1
   else
         ipick = ipick + 1
         if ( ipick .ge. nx) then
!           all values are the same - return an answer
            if ( x .eq. xarray(1)) then
                fitmd = yarray(1)
                return
            endif
            write (6,901)
901         format (1x, 'error in fitmd - constant x ')
            stop 'fitmd-cnst'
         endif
         go to 101
   endif
!
!  check array bounds
!
   if ( (x .lt. xarray(1) .and. idirection .eq.  1) .or. &
  &      (x .gt. xarray(1) .and. idirection .eq. -1)) then
         fitmd = 0.0
         return
   endif
   if ( (x .lt. xarray(nx) .and. idirection .eq. -1) .or. &
  &      (x .gt. xarray(nx) .and. idirection .eq.  1)) then
         fitmd = 0.0
         return
   endif
!  bracket the array points
!
   ipick = 0
   if ( idirection .eq. 1) then
      do i=1,nx
         place = x - xarray(i)
         if ( place .lt. 0.0) go to 909
         ipick = i
      enddo
   elseif ( idirection .eq. -1) then
      do i=1,nx
         place = xarray(i) - x
         if ( place .lt. 0.0) go to 909
         ipick = i
      enddo
   endif
   stop 'fitmd-err'
909 continue
!
!  perform interpolation
!
   ipick1 = ipick + 1
   call terp1(xarray(ipick), yarray(ipick), xarray(ipick1), &
  &           yarray(ipick1), x, y, icode)
   fitmd = y
   return
   end function fitmd

   subroutine conbar(e,ebar,dame,nin,nscr,matd,mtd,&
     c,ncmax,b,nbmax,yld,q)
   !-------------------------------------------------------------------
   ! Calculate neutron ebar and damage energy
   ! for continuum scattering reactions.
   ! Damage for (n,np), (n,na), etc., is currently treated like
   ! damage for (n,p), (n,a), etc.  neutron recoil is ignored.
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   use mainio  ! provides nsysi,nsyso
   use snl     ! provides SNL
   ! externals
   integer::nin,nscr,matd,mtd,ncmax,nbmax
   real(kr)::e,ebar,dame,c(ncmax),b(nbmax),yld,q
   ! internals
   integer::l,mf1,mt1,mfd,idis,nb,nw,iz,il,ik,nktot
   integer::ip,ir,idisc,nr,np,ll,nnr,inn,ie,idone,ii,mtt
   integer::lf,nne,lnow,iraw,ne
   real(kr)::enext,t1,t2,ehi,fhi,dhi,pe,eihi,etp,f,d,s
   real(kr)::awr,awfac,aw1fac,z,za
   integer,parameter::nkmax=12
   integer::loc(nkmax)
   real(kr)::d1(nkmax),d2(nkmax),e1(nkmax),e2(nkmax)
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::step=1.5e0_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0
   save nktot,awr,awfac,aw1fac,z,za,mtt
   save loc,d1,d2,e1,e2,mf1,mt1

   if (imode(3) .lt. 0) then
     write (nsyso,2301) e,matd,mtd
2301 format (/,1x, 'CONBAR entry ', g14.7, 2i6,/)
   endif

   !--initialize if e=0.
   if (e.eq.zero) then
      l=1
      if ((mtd.lt.18.or.mtd.gt.21).and.mtd.ne.38) then
         yld=1
         if (mtd.eq.16.or.mtd.eq.30) yld=2
         if (mtd.eq.24.or.mtd.eq.26) yld=2
         if (mtd.eq.11.or.mtd.eq.41) yld=2
         if (mtd.eq.17.or.mtd.eq.25) yld=3
         if (mtd.eq.38.or.mtd.eq.42) yld=3
         if (mtd.eq.37) yld=4
      else
         mf1=1
         mt1=452
         if (nply.gt.0) mt1=456
         ! delayed neutrons are treated the same as fast neutrons.
         ! this will cause a slight overestimate of en for fission.
         ! ... but when nply.ne.0, use prompt neutron data only
         call hgtyld(e,enext,idis,yld,matd,mf1,mt1,nscr,b,nbmax)
      endif
      mfd=5
      call findf(matd,mfd,mtd,nin)
      call contio(nin,0,0,c(l),nb,nw)
      nktot=nint(c(5))
      if (nktot.gt.nkmax) call error('conbar','nktot gt nkmax.',' ')
      awr=c(2)
      awfac=1/awr
      aw1fac=1/(awr+1)
      iz=nint(c(1)/1000)
      z=iz
      ! treat damage for (n,np) as (n,p),etc.
      za=c(1)
      mtt=0
      if (mtd.eq.22) mtt=107
      if (mtd.eq.28) mtt=103
      if (mtd.eq.32) mtt=104
      if (mtd.eq.33) mtt=105
      if (mtd.eq.34) mtt=106
      if (mtt.gt.0) call capdam(e,dame,q,za,awr,mtt)
      ne=0
      il=1
      ik=0
      do while (ik.lt.nktot)
         ik=ik+1
         d2(ik)=0
         e2(ik)=0
         loc(ik)=il
         l=il
         call tab1io(nin,0,0,c(l),nb,nw)
         l=l+nw
         do while (nb.ne.0)
            call moreio(nin,0,0,c(l),nb,nw)
            l=l+nw
         enddo
         lf=nint(c(il+3))
         il=l

         !--analytic subsection.
         if (lf.gt.1) then
            if (lf.eq.3) then
               t1=(awr*awr+1)*aw1fac*aw1fac
               t2=awr*aw1fac
            else
               call tab1io(nin,0,0,c(l),nb,nw)
               l=l+nw
               do while (nb.ne.0)
                  if (l.gt.ncmax) call error('conbar',&
                    'insufficient storage for raw endf data.',' ')
                  call moreio(nin,0,0,c(l),nb,nw)
                  l=l+nw
               enddo
               il=l
               if (lf.eq.5.or.lf.eq.11) then
                  call tab1io(nin,0,0,c(l),nb,nw)
                  l=l+nw
                  do while (nb.ne.0)
                     if (l.gt.ncmax) call error('conbar',&
                       'insufficient storage for raw endf data.',' ')
                     call moreio(nin,0,0,c(l),nb,nw)
                     l=l+nw
                  enddo
                  il=l
               endif
            endif

         !--tabulated subsection
         else
            call tab2io(nin,0,0,c(l),nb,nw)
            ne=nint(c(l+5))
            l=l+nw
            nne=0
            lnow=l
            iraw=lnow+3*ne

            !--read spectrum for each incident energy
            !--compute heating and damage contributions
            do while (nne.lt.ne)
               l=iraw
               call tab1io(nin,0,0,c(l),nb,nw)
               l=l+nw
               do while (nb.ne.0)
                  if (l.gt.ncmax) call error('conbar',&
                    'insufficient storage for raw endf data.',' ')
                  call moreio(nin,0,0,c(l),nb,nw)
                  l=l+nw
               enddo
               nne=nne+1
               ehi=c(iraw+1)
               call tabbar(fhi,c(iraw),lf)
               dhi=0
               if (idame.ne.0) then
                  if (mtt.ne.0) then
                     call capdam(ehi,dhi,q,za,awr,mtt)
                  else
                     call tabdam(ehi,dhi,z,awr,awfac,c(iraw))
                  endif
               endif
               c(lnow+3*(nne-1))=ehi
               c(lnow+3*(nne-1)+1)=fhi
               c(lnow+3*(nne-1)+2)=dhi
            enddo
            il=lnow+3*ne
         endif
      enddo
      return
   endif

   !--normal entry.
   enext=emax
   idis=0
   ebar=0
   dame=0

   !--get energy dependent yld.
   if (mtd.eq.18.or.mtd.eq.19.or.mtd.eq.20.or.mtd.eq.21.or.mtd.eq.38) then
      mf1=1
      call hgtyld(e,enext,idis,yld,matd,mf1,mt1,nscr,b,nbmax)
   endif

   !--loop over subsections.
   do ik=1,nktot
      lnow=loc(ik)
      lf=nint(c(lnow+3))

      !--interpolate for fractional probability.
      ip=2
      ir=1
      call terpa(pe,e,eihi,idisc,c(lnow),ip,ir)
      if (abs(eihi-enext).lt.small*enext.and.idisc.gt.idis) idis=idisc
      if (eihi.lt.enext*(1-small)) idis=idisc
      if (eihi.lt.enext*(1-small)) enext=eihi
      nr=nint(c(lnow+4))
      np=nint(c(lnow+5))
      ll=lnow+6+2*nr+2*(np-1)
      etp=c(ll)
      if (pe.gt.zero) then

         !--tabulated subsection.
         !--interpolate for results
         if (lf.eq.1) then
            lnow=lnow+6+2*nr+2*np
            ne=nint(c(lnow+5))
            nnr=nint(c(lnow+4))
            inn=2
            lnow=lnow+6+2*nnr
            ie=0
            idone=0
            do while (ie.lt.ne.and.idone.eq.0)
               ie=ie+1
               if (c(lnow+3*ie).ge.e*(1-small)) then
                  ii=lnow+3*(ie-1)
                  idone=1
               endif
            enddo
            if (idone.eq.0) ii=lnow+3*(ne-1)
            call terp1(c(ii),c(ii+1),c(ii+3),c(ii+4),e,f,inn)
            ebar=ebar+f*pe
            if (idame.ne.0) then
               call terp1(c(ii),c(ii+2),c(ii+3),c(ii+5),e,d,inn)
               dame=dame+d*pe
            endif

         !--analytic subsection.
         else
            call anabar(s,e,c(lnow),t1,t2)
            ebar=ebar+s*pe
            if (idame.ne.0) then
               if (mtt.ne.0) then
                  call capdam(e,dame,q,za,awr,mtt)
               else
                  if (e.ge.e2(ik)*(1-small)) then
                     idone=0
                     do while (idone.eq.0)
                        d1(ik)=d2(ik)
                        e1(ik)=e2(ik)
                        e2(ik)=step*e1(ik)
                        if (e2(ik).gt.etp*(1+small)) e2(ik)=etp
                        if (e1(ik).eq.zero) e2(ik)=e
                        call anadam(e2(ik),d2(ik),z,awr,c(lnow))
                        if (abs(e2(ik)-etp).lt.small*etp) then
                           idone=1
                        else
                           if (e1(ik).ne.zero.and.&
                             e.le.e2(ik)*(1+small)) idone=1
                        endif
                     enddo
                  endif
                  call terp1(e1(ik),d1(ik),e2(ik),d2(ik),e,d,2)
                  dame=dame+d*pe
               endif
            endif
         endif
      endif
   enddo
   end subroutine conbar

   subroutine hgtyld(e,enext,idis,yld,mat,mf,mt,itape,a,na)
   !-------------------------------------------------------------------
   ! Retrieve or compute yield for requested mat,mf,mt on itape.
   ! If mt=455, retrieve delayed neutron time constants.
   ! Initialize if e=0.
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   use mainio  ! provides nsysi,nsyso
   use snl     ! provides SNL
   ! externals
   integer::idis,mat,mf,mt,itape,na
   real(kr)::e,enext,yld,a(*)
   ! internals
   integer::nb,nw,lnu,lnd,loc,ip,ir,nr,i,nc
   real(kr)::term
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::zero=0
   save lnu,ip,ir

   if (imode(3) .lt. 0) then
     write (nsyso,2301) e, enext, yld, mt
2301 format (1x, 'HGTYLD entry ', 3g14.7, i6)
   endif

   !--initialize.
   if (e.eq.zero) then
      call findf(mat,mf,mt,itape)
      call contio(itape,0,0,a,nb,nw)
      lnu=nint(a(4))
      if (mth.eq.455) then
         call listio(itape,0,0,a,nb,nw)
         lnd=nint(a(5))
         if (lnd.ne.6.and.lnd.ne.8)&
           call error('hgtyld','illegal lnd, must be 6 or 8',' ')
         loc=1
         call tab1io(itape,0,0,a(loc),nb,nw)
         loc=loc+nw
         if (loc.gt.na) call error('hgtyld','storage exceeded.',' ')
         do while (nb.ne.0)
            call moreio(itape,0,0,a(loc),nb,nw)
            loc=loc+nw
            if (loc.gt.na) call error('hgtyld','storage exceeded.',' ')
         enddo
         na=loc-1
         ip=2
         ir=1
         nr=nint(a(5))
         enext=a(7+2*nr)
      else
         lnd=0
         loc=1
         if (lnu.ne.1) then
            call tab1io(itape,0,0,a(loc),nb,nw)
            loc=loc+nw
            if (loc.gt.na) call error('hgtyld','storage exceeded.',' ')
            do while (nb.ne.0)
               call moreio(itape,0,0,a(loc),nb,nw)
               loc=loc+nw
               if (loc.gt.na) call error('hgtyld','storage exceeded.',' ')
            enddo
            na=loc-1
            ip=2
            ir=1
            nr=nint(a(5))
            enext=a(7+2*nr)
         else
            call listio(itape,0,0,a,nb,nw)
            na=nw
            enext=emax
         endif
      endif
      idis=0
      return
   endif

   !--tabulated data.
   if (lnu.ne.1) then
      call terpa(yld,e,enext,idis,a,ip,ir)

   !--polynomial data.
   else
      yld=a(7)
      term=1
      nc=nint(a(5))
      do i=2,nc
         term=term*e
         yld=yld+a(i+6)*term
      enddo
      enext=emax
      idis=0
   endif
   return
   end subroutine hgtyld

   subroutine anabar(s,e,a,t1,t2)
   !-------------------------------------------------------------------
   ! Compute mean energy for an analytic section of File 5.
   !-------------------------------------------------------------------
   use physics ! provides pi
   use endf    ! provides terpa
   ! externals
   real(kr)::s,e,a(*),t1,t2
   ! internals
   integer::lf,nr,np,lnext,ip,ir,idis,nrn,npn,lnxt
   real(kr)::u,theta,enext,b1,b2,b3,expa,test,b,rpi2
   real(kr),parameter::thrhaf=1.5e0_kr
   real(kr),parameter::zero=0
   rpi2=sqrt(pi)/2

   !--retrieve common factors.
   s=0
   u=a(1)
   if (e.le.u) return
   lf=nint(a(4))
   nr=nint(a(5))
   np=nint(a(6))
   lnext=7+2*nr+2*np
   ip=2
   ir=1
   call terpa(theta,e,enext,idis,a(lnext),ip,ir)

   !--law 3--discrete energy spectrum.
   if (lf.eq.3) then
      s=e*t1-t2*theta

   !--law 5--general evaporation spectrum.
   else if (lf.eq.5) then
      nrn=nint(a(lnext+4))
      npn=nint(a(lnext+5))
      lnxt=lnext+2*nrn+2*npn+6
      call tabbar(s,a(lnxt),lf)
      s=s*theta

   !--law 7--simple fission spectrum.
   else if (lf.eq.7) then
      if (theta.ne.zero) then
         b1=(e-u)/theta
         b2=sqrt(b1)
         b3=exp(-b1)
         expa=exp(-b2*b2)
         s=theta*(thrhaf-b1*b2*b3/(rpi2*(1-expa*rerfc(b2))-b2*b3))
      endif

   !--law 9--evaporation spectrum.
   else if (lf.eq.9) then
      if (theta.ne.zero) then
         b1=(e-u)/theta
         b2=sqrt(b1)
         b3=exp(-b1)
         test=1
         test=test/1000
         if (b1.ge.test) then
            s=theta*(2-b1*b1*b3/(1-(b1+1)*b3))
         else
            s=4*(e-u)/3
         endif
      endif

   !--law 11--energy-dependent watt spectrum.
   else if (lf.eq.11) then
      nrn=nint(a(lnext+4))
      npn=nint(a(lnext+5))
      lnxt=lnext+6+2*nrn+2*npn
      ip=2
      ir=1
      call terpa(b,e,enext,idis,a(lnxt),ip,ir)
      s=theta*(thrhaf+theta*b/4)

   !--law 12--energy dependent fission neutron spectrum (madland-nix)
   else if (lf.eq.12) then
          s=((a(lnext)+a(lnext+1))/2)+(4*theta/3)
   endif
   return

   contains
      real(kr) function r(z)
      real(kr)::z
      real(kr),parameter::a0=.3275911e0_kr
      r=1/(1+a0*abs(z))
      return
      end function r

      real(kr) function rerfc(z)
      real(kr)::z
      real(kr),parameter::a1=.254829592e0_kr
      real(kr),parameter::a2=-.284496736e0_kr
      real(kr),parameter::a3=1.421413741e0_kr
      real(kr),parameter::a4=-1.453152027e0_kr
      real(kr),parameter::a5=1.061405429e0_kr
      rerfc=r(z)*(a1+r(z)*(a2+r(z)*(a3+r(z)*(a4+a5*r(z)))))
      return
      end function rerfc

   end subroutine anabar

   subroutine anadam(e,d,z,awr,a)
   !-------------------------------------------------------------------
   ! Compute damage energy for an analytic section of File 5
   !-------------------------------------------------------------------
   use endf ! provides terpa
   ! externals
   real(kr)::e,d,z,awr,a(*)
   ! internals
   integer::nr,np,lf,lnext,ip,ir,idis,i,iq,j,iflag
   real(kr)::u,awfac,test,theta,enext,ep,er,xlast,ylast,xm,ym,yt
   integer,parameter::imax=10
   real(kr)::x(imax),y(imax)
   integer,parameter::nq=4
   real(kr),dimension(nq),parameter::qp=(/&
     -.86114e0_kr,-.33998e0_kr,.33998e0_kr,.86114e0_kr/)
   real(kr),dimension(nq),parameter::qw=(/&
     .34785e0_kr,.65215e0_kr,.65215e0_kr,.34785e0_kr/)
   real(kr),parameter::eps=.05e0_kr
   real(kr),parameter::delta=1.e-7_kr

   !--retrieve spectrum temperature
   d=0
   u=a(1)
   nr=nint(a(5))
   np=nint(a(6))
   lf=nint(a(4))
   awfac=1/awr
   if (lf.ne.9) return
   if (e.le.u) return
   if ((e-u).lt.delta*e) return
   test=10
   if ((e-u).lt.test) return
   lnext=7+2*nr+2*np
   ip=2
   ir=1
   call terpa(theta,e,enext,idis,a(lnext),ip,ir)

   !--adaptive energy integration
   !--gauss-legendre angle integration
   i=4
   do while (i.gt.1)
      i=i-1
      if (i.eq.3) ep=1
      if (i.eq.2) ep=(e-u)/2
      if (i.eq.2.and.theta.lt.e-u) ep=theta
      if (i.eq.1) ep=e-u
      x(i)=ep
      y(i)=0
      do iq=1,nq
         er=(e-2*sqrt(e*ep)*qp(iq)+ep)*awfac
         y(i)=y(i)+qw(iq)*df(er,z,awr,z,awr)/2
      enddo
      y(i)=y(i)*sed(e,ep,theta,u)
   enddo
   i=3
   j=0
   xlast=0
   ylast=0
   do while (i.gt.1)
      ! test for convergence
      iflag=0
      if (i.gt.1.and.i.lt.imax) then
         xm=(x(i-1)+x(i))/2
         ym=(y(i-1)+y(i))/2
         yt=0
         do iq=1,nq
            er=(e-2*sqrt(e*xm)*qp(iq)+xm)*awfac
            yt=yt+qw(iq)*df(er,z,awr,z,awr)/2
         enddo
         yt=yt*sed(e,xm,theta,u)
         if (abs(yt-ym).gt.eps*yt) iflag=1
      endif
      ! fails test.
      ! add midpoint to stack and continue.
      if (iflag.eq.1) then
         i=i+1
         x(i)=x(i-1)
         x(i-1)=xm
         y(i)=y(i-1)
         y(i-1)=yt
      ! passes test.
      ! use top point off the stack.
      else
         j=j+1
         if (j.gt.1) d=d+(x(i)-xlast)*(y(i)+ylast)/2
         xlast=x(i)
         ylast=y(i)
         i=i-1
      endif
   enddo
   return
   end subroutine anadam

   real(kr) function sed(e,ep,theta,u)
   !-------------------------------------------------------------------
   ! Secondary energy distribution function for anadam.
   ! A Taylor expansion is used for small values to avoid overflow.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::e,ep,theta,u
   ! internals
   real(kr)::xeu,aeu,test

   xeu=(e-u)/theta
   aeu=abs(xeu)
   test=1
   test=test/10000
   if (aeu.ge.test) then
      sed=ep*exp(-ep/theta)/(theta*theta*(1-exp(-xeu)*(1+xeu)))
   else
      sed=ep*exp(-ep/theta)/(theta*theta*xeu*xeu*(1-xeu)/2)
   endif
   return
   end function sed

  subroutine tabbar(f,a,law)
   !-------------------------------------------------------------------
   ! Compute mean energy for a tabulated section of File 5.
   ! (If law<0, the mean energy of a subsection of File 6
   ! can be computed using abs(law) as the interpolation law).
   !-------------------------------------------------------------------
   use util ! provides error
   use mainio  ! provides nsysi,nsyso
   use snl     ! provides SNL
   ! externals
   integer::law
   real(kr)::f,a(*)
   ! internals
   integer::nr,np,ibase,ir,nbt,inn,ncyc,i
   real(kr)::xh,xhsq,yh,xl,xlsq,delxsq,yl,b
   real(kr)::alxl,alxh,alyl,alyh
   real(kr)::photon_count

   !--initialize.
   if (law.gt.0.and.law.ne.1.and.law.ne.5)&
     call error('tabbar','coded for lf=1 and lf=5 only.',' ')
   if (law.ge.0) then
      f=0
      photon_count = 0.0
      nr=nint(a(5))
      np=nint(a(6))
      ibase=6+2*nr
      ir=1
      nbt=nint(a(7))
      inn=nint(a(8))
      ncyc=2
   else
      f=0
      photon_count = 0.0
      np=nint(a(6))
      inn=iabs(law)
      ncyc=nint(a(4))+2
      ibase=6
   endif
   xh=a(ibase+1)
   xhsq=xh*xh
   yh=a(ibase+2)

   !--accumulate contributions to integral.
   do i=2,np
      xl=xh
      xh=a(ibase+ncyc*(i-1)+1)
      xlsq=xhsq
      xhsq=xh*xh
      delxsq=xhsq-xlsq
      yl=yh
      yh=a(ibase+ncyc*(i-1)+2)

      if (law.gt.0.and.i.gt.nbt) then
         ir=ir+1
         nbt=nint(a(ibase+2*ir-1))
         inn=nint(a(ibase+2*ir))
      endif
      if (xl.ne.xh) then

         !--y is a constant.
         if (inn.eq.1) then
            f=f+yl*delxsq/2
            photon_count = photon_count + (xh-xl)*(yh+yl)/2.

         !--y is linear in x.
         else if (inn.eq.2) then
            b=(yh-yl)/(xh-xl)
            f=f+((yl-xl*b)*delxsq/2+b*(xh*xhsq-xl*xlsq)/3)
            photon_count = photon_count + (xh-xl)*(yh+yl)/2.

         !--y is linear in ln(x)
         else if (inn.eq.3) then
            alxl=log(xl)
            alxh=log(xh)
            b=(yh-yl)/(alxh-alxl)
            f=f+(yl-b*alxl+b/2)*delxsq/2+b*(xhsq*alxh-xlsq*alxl)/2
            photon_count = photon_count + (xh-xl)*(yh+yl)/2.

         !--ln(y) is linear in x
         else if (inn.eq.4) then
            alyl=log(yl)
            alyh=log(yh)
            b=(alyh-alyl)/(xh-xl)
            f=f+(yh*(b*xh-1)-yl*(b*xl-1))/b**2
            photon_count = photon_count + (xh-xl)*(yh+yl)/2.

         !--ln(y) is linear in ln(x)
         else if (inn.eq.5) then
            alxl=log(xl)
            alxh=log(xh)
            alyl=log(yl)
            alyh=log(yh)
            b=(alyh-alyl)/(alxh-alxl)
            f=f+yl*xlsq*((xh/xl)**(b+2)-1)/(b+2)
            photon_count = photon_count + (xh-xl)*(yh+yl)/2.
         endif
      endif

      if (imode(3) .lt. 0) then
        write (nsyso,2306) i, inn, xl, xh, yl, yh, f, photon_count
2306    format (1x, 'TABBAR accumulation ', 2i5, 6g14.7)
      endif

   enddo
   return
   end subroutine tabbar

   subroutine tabdam(e,d,z,awr,awfac,a)
   !-------------------------------------------------------------------
   ! Compute damage energy for a tabulated section of File 5.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::e,d,z,awr,awfac,a(*)
   ! internals
   integer::nr,np,ibase,i,iq
   real(kr)::xh,yh,fl,xl,f,er
   integer,parameter::nq=4
   real(kr),dimension(nq),parameter::qp=(/&
     -.86114e0_kr,-.33998e0_kr,.33998e0_kr,.86114e0_kr/)
   real(kr),dimension(nq),parameter::qw=(/&
     .34785e0_kr,.65215e0_kr,.65215e0_kr,.34785e0_kr/)
   real(kr),parameter::zero=0

   !--initialize
   if (e.eq.zero) then
      d=df(e,z,awr,z,awr)
      return
   endif

   !--trapazoidal energy integration
   !--gauss-legendre angle integration
   d=0
   nr=nint(a(5))
   np=nint(a(6))
   ibase=6+2*nr
   xh=a(ibase+1)
   yh=a(ibase+2)
   fl=0
   do i=1,np
      xl=xh
      xh=a(ibase+2*i-1)
      yh=a(ibase+2*i)
      f=0
      do iq=1,nq
         er=(e-2*sqrt(e*xh)*qp(iq)+xh)*awfac
         f=f+qw(iq)*df(er,z,awr,z,awr)/2
      enddo
      f=f*yh
      if (i.gt.1) d=d+(xh-xl)*(f+fl)/2
      fl=f
   enddo
   return
   end subroutine tabdam

   subroutine sixbar(e,ebar,yld,dame,nin,c,ncmax,nscr,b,nbmax,&
     n2,j6,irec,jrec,iflag)
   !-------------------------------------------------------------------
   ! Calculate charged-particle ebar and damage energy
   ! for reactions represented using ENDF-6 File 6.
   ! If recoil subsection is not given, attempt to compute
   ! from the distribution for the heavy emitted particle.
   ! For normal particles, irec is zero.  For recoil cases,
   ! irec is the index of the section for the driving particle.
   !-------------------------------------------------------------------
   use util ! provides skiprz,mess
   use endf ! provides endf routines and variables
   use mainio  ! provides nsysi,nsyso
   use snl     ! provides SNL
   ! externals
   integer::nin,n2,j6,irec,jrec,ncmax,nbmax,nscr,iflag
   real(kr)::e,ebar,yld,dame
   real(kr)::c(ncmax),b(nbmax)
   ! internals
   integer::matd,mfd,mtd,l,nb,nw,ik,nnt,nmu,imu,mf1,mt1,idis
   integer::ip,ir,idisc,iraw,nne,ne,law,lang,lep,intl
   real(kr)::disc102,zp,zt,ap,at,ztt,pe,eihi,f,d,s,enext
   real(kr)::elo,ehi,flo,fhi,dlo,dhi
   character(60)::strng
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::zero=0
   save iraw,nne,ne,law,lang,lep,intl
   save elo,ehi,flo,fhi,dlo,dhi
   save disc102,zp,ap,zt,at

   if (imode(3) .lt. 0) then
     write (nsyso,2301) e, disc102, elo,ehi,flo,fhi,dlo,dhi, pe, yld, disc102, mth, irec, jrec
2301 format (1x, 'SIXBAR entry ', 11g14.7, 3i6)
   endif

   !--initialize for this subsection when e=0.
   if (e.gt.zero) go to 300
   if (irec.eq.0) go to 110

   if (imode(3) .lt. 0) then
     write (nsyso,5865) e, irec, jrec
5865 format (1x, 'SIXBAR irec step-D ', g14.7, 2i6)
   endif

   if (jrec.gt.0) irec=jrec

   !--for recoils, back up to the driving particle,
   !--or restore position after doing recoil.
  100 continue
   matd=math
   mfd=6
   mtd=mth
   call skiprz(nin,-2)
   call findf(matd,mfd,mtd,nin)
   l=1
   call contio(nin,0,0,c(l),nb,nw)
   ik=0
   do while (ik.lt.irec-1)
      ik=ik+1
      call tab1io(nin,0,0,c(l),nb,nw)
      law=l2h
      call skip6(nin,0,0,c(l),law)
   enddo
   if (jrec.eq.irec) then
      irec=0
      jrec=0
   endif

   if (imode(3) .lt. 0) then
     write (nsyso,5861) e, c(1), c(2), irec, jrec
5861 format (1x, 'SIXBAR irec step-A ', 3g14.7, 2i6)
   endif

   !--start reading in the desired subsection
  110 continue
   l=1
   call tab1io(nin,0,0,c(l),nb,nw)
   zap=c1h
   awp=c2h
   law=l2h
   l=l+nw
   do while (nb.ne.0)
      call moreio(nin,0,0,c(l),nb,nw)
      l=l+nw
   enddo
   iflag=0
   disc102=0
   if (zap.eq.zero) then
      iflag=1
      disc102=awp
      awp=0
   endif
   ebar=0
   yld=0
   if (law.gt.0) go to 150

   !--warning for law=0.
   write(strng,'(''no distribution for mt'',i3,'' particle '',i5)')&
     mth,nint(zap)
   call mess('sixbar',strng,' ')
   ebar=-2
   return

   !--watch for 2-body recoils
  150 continue
   if (law.ne.4) go to 210
   irec=j6-1
   jrec=j6+1

   if (imode(3) .lt. 0) then
     write (nsyso,5862) e, irec, jrec
5862 format (1x, 'SIXBAR irec step-B ', g14.7, 2i6)
   endif

   go to 100

   !--continue processing this subsection
  210 continue
   zp=int(zap/1000)
   zt=int(zat/1000)
   if (irec.gt.0) zp=zt-zp
   ap=awp
   if (irec.gt.0) ap=awr+1-awp
   at=awr
   if (zap.eq.zero) then
       ap=awr+1
       zp=zt
   endif
   dame=df(e,zp,ap,zt,at)
   if (disc102.gt.zero) go to 295
   if (law.ne.3.and.law.ne.6) then
      call tab2io(nin,0,0,c(l),nb,nw)
      lang=nint(c(l+2))
      lep=nint(c(l+3))
      ne=n2h
      nnt=l+6
      intl=nint(c(nnt+1))
      intl=2
      l=l+nw
   endif
   nne=0
   iraw=l
   if (law.eq.3) go to 230
   if (law.eq.6) go to 240
   if (law.eq.7) go to 250

   !--laws 1, 2, and 5.
   call listio(nin,0,0,c(l),nb,nw)
   l=l+nw
   do while (nb.ne.0)
      call moreio(nin,0,0,c(l),nb,nw)
      l=l+nw
   enddo
   nne=nne+1
   elo=c(iraw+1)
   ehi=0
   dhi=0
   fhi=0
   go to 290

   !--law 3 -- isotropic distribution
  230 continue
   ehi=emax
   lang=0
   return

   !--law 6 -- phase-space distribution
  240 continue
   call contio(nin,0,0,c(l),nb,nw)
   lang=0
   lep=2
   l=l+nw
   ehi=emax
   dhi=0
   fhi=0
   return

   !--law 7 -- laboratory angle-energy law
  250 continue
   call tab2io(nin,0,0,c(l),nb,nw)
   nmu=n2h
   l=l+nw
   do imu=1,nmu
      call tab1io(nin,0,0,c(l),nb,nw)
      l=l+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,c(l),nb,nw)
         l=l+nw
      enddo
   enddo
   nne=nne+1
   elo=c(iraw+1)
   ehi=0
   dhi=0
   fhi=0

   !--compute low values of ebar and dame
  290 continue

   if (imode(3) .lt. 0 .and. mth .eq. 102) then
     write (nsyso,4342) elo, flo, dlo, ebar, yld, zap, c(1), c(2), c(iraw), irec, iraw
4342 format (1x, 'Before GETSIX call from SIXBAR at low return: ', 9g14.7, 2i6)
   endif

   if (mth.ne.102.or.(zap.eq.0.and.irec.eq.0)) then
      call getsix(elo,flo,dlo,c(iraw),law,lang,lep,irec)

      if (imode(3) .lt. 0) then
        write (nsyso,4341) elo, flo, dlo, ebar, yld
4341    format (1x, 'SIXBAR low return after GETSIX ', 5g14.7)
      endif

   else
      ztt=int(zat/1000)

      if (mth .eq. 102 .and. irec .gt. 0) then 
          if (imode(3) .lt. 0) then
             write (nsyso,5341) pe_archive, pe_store
5341         format (1x, 'SIXBAR low return pe reset ', 2g14.7)
          endif

!          ip=2
!          ir=1
!          call terpa(pe_trial,e,eihi,idisc,c(1),ip,ir)

          if (imode(3) .lt. 0 .and. e .le. 1.E-4) then
             write (nsyso,9306) e, irec, pe_trial, disc102, mth, eihi, c(1), idisc, ip, ir
9306         format (1x, 'SIXBAR low pe_trial terpa reset ', e14.7, i6, 2g14.7, i6, 2g14.7, 3i6)
          endif
         
      endif


   if (imode(3) .lt. 0) then
     write (nsyso,2331) e, flo, dlo, ebar, yld, pe_archive, pe_store
2331 format (1x, 'Before TABSQ6 low call in SIXBAR ', 7g14.7)
   endif

      call tabsq6(e, flo,dlo,c(iraw),law,ztt,awrt)

   if (imode(3) .lt. 0) then
     write (nsyso,2311) e, flo, dlo, ebar, yld, c(iraw), iraw
2311 format (1x, 'TABSQ6 low return in SIXBAR ', 6g14.7, i6)
   endif

   endif
   if ((mth.ge.18.and.mth.le.21).or.mth.eq.38) then
      matd=math
      mf1=1
      mt1=452
      if (nply.gt.0) mt1=456
      ! delayed neutrons are treated the same as fast neutrons.
      ! this will cause a slight overestimate of en for fission.
      ! ... but when nply.ne.0, use prompt neutron data only.
      call hgtyld(e,enext,idis,yld,matd,mf1,mt1,nscr,b,nbmax)
   endif

   if (imode(3) .lt. 0) then
     write (nsyso,2371) flo, dlo, ebar, yld, irec
2371 format (1x, 'SIXBAR low value exit ', 4g14.7, i6)
   endif

   return

  295 continue
   yld=1
   call skip6(nin,0,0,c(l),law)
   return

   !--normal entry
  300 continue
      if (disc102.gt.zero) go to 430

   !--interpolate for particle yield
   ip=2
   ir=1
   call terpa(pe,e,eihi,idisc,c(1),ip,ir)

   if (imode(3) .lt. 0 .and. e .le. 1.E-3) then
     write (nsyso,8306) e, irec, pe, disc102, mth, eihi, c(1), c(2), idisc, ip, ir
8306 format (1x, 'SIXBAR pe over-ride after terpa ', e14.7, i6, 2g14.7, i6, 3g14.7, 3i6)
   endif

!   Temporary change to fix MF12 capture gamma total kerma and displacement kerma
!   inhibit normal coding to stop multiplicity change in MF12 capture gamma case 
!   but - permit this for MF6 cases.
!   Inhibited!
!   [PJG, 9/14/2020]
!
!     Legacy code:
       pe_archive = pe
       if (irec.gt.0) pe=1
       pe_store = pe
!     Replacement code:
!        if ( irec.gt.0 .and. mth .ne. 102) pe = 1
!        if ( irec.gt. 0 .and. mth .eq. 102 .and. mfd .eq. 0) pe = 1
        if ( irec.gt.0 .and. mth .eq. 102 .and. e .le. 0.101E-4) then
            write (nsyso,5356) e, pe, pe_archive, irec, mth, mfd
5356       format (1x, '*** Warning ***: temporary pe over-ride xonsidered ', 3e14.7, 3i6)
        endif 
!   End of temporary change

   if (imode(3) .lt. 0 .and. e .le. 1.E-4) then
     write (nsyso,8305) e, irec, pe, pe_archive
8305 format (1x, 'SIXBAR loc-8305 ', e14.7, i6, 2g14.7)
   endif

   !--is desired energy in current panel
  305 continue

   if (imode(3) .lt. 0 .and. e .le. 1.e-4) then
     write (nsyso,2305) e, disc102, ehi, nne, ne, law, irec
2305 format (1x, 'SIXBAR loc-305 ', 3g14.7, 4i6)
   endif

   if (e.lt.ehi*(1-small)) go to 400
   if (nne.eq.1) go to 310
   if (nne.eq.ne.and.e.le.ehi*(1+small)) go to 400
   if (nne.eq.ne) go to 450

   !--no. slide high-energy data into low-energy positions
   !--and read in new high energy data.
   elo=ehi
   flo=fhi
   dlo=dhi

   if (imode(3) .lt. 0) then
     write (nsyso,4311) elo, flo, dlo
4311 format (1x, 'SIXBAR low fill ', 2g14.7)
   endif

   !--read in new data for high energy
   !--laws 1, 2, and 5.
  310 continue
   if (law.eq.1.or.law.eq.2.or.law.eq.5) then
      l=iraw
      call listio(nin,0,0,c(l),nb,nw)
      l=l+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,c(l),nb,nw)
         l=l+nw
      enddo
      nne=nne+1
      ehi=c(iraw+1)

   !--law 6.
   else if (law.eq.6) then
      l=iraw
         call contio(nin,0,0,c(l),nb,nw)

   !--law 7.
   else if (law.eq.7) then
      l=iraw
      call tab2io(nin,0,0,c(l),nb,nw)
      nmu=n2h
      l=l+nw
      do imu=1,nmu
         call tab1io(nin,0,0,c(l),nb,nw)
         l=l+nw
         do while (nb.ne.0)
            call moreio(nin,0,0,c(l),nb,nw)
            l=l+nw
         enddo
      enddo
      ehi=c(iraw+1)
      nne=nne+1
   endif

   !--compute high values of ebar and dame
   if (mth.ne.102.or.(zap.eq.0.and.irec.eq.0)) then
      call getsix(ehi,fhi,dhi,c(iraw),law,lang,lep,irec)
   else
      ztt=int(zat/1000)
      call tabsq6(e, fhi,dhi,c(iraw),law,ztt,awrt)

      if (imode(3) .lt. 0) then
        write (nsyso,2401) e, fhi, dhi, c(iraw), law
2401    format (1x, 'TABSQ6 high return in SIXBAR ', 4g14.7, i6)
      endif

!     Code change: PJG 9/21/2020
!     flo/dlo was not set correctly in above logic because pe was not properly defined at call
!     so, we reset flo/dlo to the fhi/dhi values here
      call tabsq6(e, flo,dlo,c(iraw),law,ztt,awrt)
!     end code change

   endif
   go to 305

   !--yes.  compute result.
  400 continue
   if (law.ne.3.and.law.ne.6) go to 420
   c(iraw+1)=e
   call getsix(e,f,d,c(iraw),law,lang,lep,irec)

   if (imode(3) .lt. 0) then
     write (nsyso,5863) e, irec, jrec
5863 format (1x, 'SIXBAR irec step-C ', g14.7, 2i6)
   endif

   ebar=f
   yld=pe
   dame=pe*d
   return

   !--yes. interpolate for results
  420 continue
   if (nne.eq.1) go to 450
   call terp1(elo,flo,ehi,fhi,e,s,intl)

   if (imode(3) .lt. 0) then
     write (nsyso,7311) elo, flo, ehi, fhi, e, s, pe, pe_archive
7311 format (1x, 'SIXBAR terp for s ', 8g14.7)
   endif

   ebar=s
   yld=pe
   if ((mth.ge.18.and.mth.le.21).or.mth.eq.38) then
      matd=math
      mf1=1
      mt1=452
      call hgtyld(e,enext,idis,yld,matd,mf1,mt1,nscr,b,nbmax)
   endif
   call terp1(elo,dlo,ehi,dhi,e,d,intl)

   if (imode(3) .lt. 0) then
     write (nsyso,7312) elo, dlo, ehi, dhi, e, d, pe, pe_archive
7312 format (1x, 'SIXBAR terp for d ', 8g14.7)
   endif

   dame=pe*d

   if (imode(3) .lt. 0) then
     write (nsyso,2304) e, ebar, yld, dame, pe, d, s, irec
2304 format (1x, 'SIXBAR exit ', 7g14.7, i6)
   endif

   return

   !--discrete relativistic capture gamma
  430 call hgam102(e,ebar,dame,disc102,c,irec,zp,ap,zt,at)
   yld=1

   if (imode(3) .lt. 0) then
     write (nsyso,2302) e, ebar, dame, yld
2302 format (1x, 'SIXBAR discrete capture gamma exit ', 4g14.7)
   endif

   return

   !--return zeros outside range of table
  450 continue
   ebar=0
   dame=0
   return
   end subroutine sixbar

   subroutine getsix(e,ebar,dame,c,law,lang,lep,irec)
   !-------------------------------------------------------------------
   ! Compute ebar and dame for one particular incident energy
   ! using data in ENDF-6 File 6 format.  For irec.eq.0, return the
   ! numbers for the emitted particle.  For irec.gt.0, return the
   ! numbers for the recoil due to a single emission of the particle.
   !-------------------------------------------------------------------
   use mathm ! provides legndr
   use mainio  ! provides nsysi,nsyso
   use snl     ! provides SNL
   ! externals
   integer::law,lang,lep,irec
   real(kr)::e,ebar,dame,c(*)
   ! internals
   integer::nl,l,i,nd,na,nep,ncyc,iq,ia,nld,il,nmu,imu
   integer::nr,np,next,ibase,iint
   real(kr)::zp,zt,ap,at,summ,epnext,el,fl,dy,da,xm
   real(kr):: delta_d, delta_h
   real(kr)::epn,ym,test,en,f,f1,h,d,xl,yl,xx,yy
   real(kr)::x2,zpp,ztt,b,er,t1,t2,thresh,beta
   real(kr)::afact,arec,u,e2,ul,fn,dl,hl,f2
   real(kr)::x(10),y(10,2)
   real(kr)::term(2)
   real(kr)::p(65)
   integer,parameter::nq=64
   real(kr),dimension(64),parameter::qp=(/&
     -9.99305042E-01_kr,-9.96340117E-01_kr,-9.91013371E-01_kr,&
     -9.83336254E-01_kr,-9.73326828E-01_kr,-9.61008800E-01_kr,&
     -9.46411375E-01_kr,-9.29569172E-01_kr,-9.10522137E-01_kr,&
     -8.89315446E-01_kr,-8.65999398E-01_kr,-8.40629296E-01_kr,&
     -8.13265315E-01_kr,-7.83972359E-01_kr,-7.52819907E-01_kr,&
     -7.19881850E-01_kr,-6.85236313E-01_kr,-6.48965471E-01_kr,&
     -6.11155355E-01_kr,-5.71895646E-01_kr,-5.31279464E-01_kr,&
     -4.89403146E-01_kr,-4.46366017E-01_kr,-4.02270158E-01_kr,&
     -3.57220158E-01_kr,-3.11322872E-01_kr,-2.64687162E-01_kr,&
     -2.17423644E-01_kr,-1.69644420E-01_kr,-1.21462819E-01_kr,&
     -7.29931218E-02_kr,-2.43502927E-02_kr, 2.43502927E-02_kr,&
      7.29931218E-02_kr, 1.21462819E-01_kr, 1.69644420E-01_kr,&
      2.17423644E-01_kr, 2.64687162E-01_kr, 3.11322872E-01_kr,&
      3.57220158E-01_kr, 4.02270158E-01_kr, 4.46366017E-01_kr,&
      4.89403146E-01_kr, 5.31279464E-01_kr, 5.71895646E-01_kr,&
      6.11155355E-01_kr, 6.48965471E-01_kr, 6.85236313E-01_kr,&
      7.19881850E-01_kr, 7.52819907E-01_kr, 7.83972359E-01_kr,&
      8.13265315E-01_kr, 8.40629296E-01_kr, 8.65999398E-01_kr,&
      8.89315446E-01_kr, 9.10522137E-01_kr, 9.29569172E-01_kr,&
      9.46411375E-01_kr, 9.61008800E-01_kr, 9.73326828E-01_kr,&
      9.83336254E-01_kr, 9.91013371E-01_kr, 9.96340117E-01_kr,&
      9.99305042E-01_kr/)
   real(kr),dimension(64),parameter::qw=(/&
      1.78328072E-03_kr, 4.14703326E-03_kr, 6.50445797E-03_kr,&
      8.84675983E-03_kr, 1.11681395E-02_kr, 1.34630479E-02_kr,&
      1.57260305E-02_kr, 1.79517158E-02_kr, 2.01348232E-02_kr,&
      2.22701738E-02_kr, 2.43527026E-02_kr, 2.63774697E-02_kr,&
      2.83396726E-02_kr, 3.02346571E-02_kr, 3.20579284E-02_kr,&
      3.38051618E-02_kr, 3.54722133E-02_kr, 3.70551285E-02_kr,&
      3.85501532E-02_kr, 3.99537411E-02_kr, 4.12625632E-02_kr,&
      4.24735151E-02_kr, 4.35837245E-02_kr, 4.45905582E-02_kr,&
      4.54916279E-02_kr, 4.62847966E-02_kr, 4.69681828E-02_kr,&
      4.75401657E-02_kr, 4.79993886E-02_kr, 4.83447622E-02_kr,&
      4.85754674E-02_kr, 4.86909570E-02_kr, 4.86909570E-02_kr,&
      4.85754674E-02_kr, 4.83447622E-02_kr, 4.79993886E-02_kr,&
      4.75401657E-02_kr, 4.69681828E-02_kr, 4.62847966E-02_kr,&
      4.54916279E-02_kr, 4.45905582E-02_kr, 4.35837245E-02_kr,&
      4.24735151E-02_kr, 4.12625632E-02_kr, 3.99537411E-02_kr,&
      3.85501532E-02_kr, 3.70551285E-02_kr, 3.54722133E-02_kr,&
      3.38051618E-02_kr, 3.20579284E-02_kr, 3.02346571E-02_kr,&
      2.83396726E-02_kr, 2.63774697E-02_kr, 2.43527026E-02_kr,&
      2.22701738E-02_kr, 2.01348232E-02_kr, 1.79517158E-02_kr,&
      1.57260305E-02_kr, 1.34630479E-02_kr, 1.11681395E-02_kr,&
      8.84675983E-03_kr, 6.50445797E-03_kr, 4.14703326E-03_kr,&
      1.78328072E-03_kr/)
   real(kr),parameter::tol=0.02e0_kr
   real(kr),parameter::small=1.e-8_kr
   real(kr),parameter::emax=1.e10_kr
   integer,parameter::imax=10
   real(kr),parameter::zero=0
   real(kr),parameter::one=1

   if (imode(3) .lt. 0) then
     write (nsyso,2301) e, ebar, zap, law, lct
2301 format (1x, 'GETSIX entry ', 3g14.7, 2i6)
   endif

   !--branch to appropriate law
   zp=int(zap/1000)
   zt=int(zat/1000)
   if (irec.gt.0) zp=zt-zp
   ap=awp
   if (irec.gt.0) ap=awrt+1-awp
   at=awrt
   if (law.eq.7) go to 510
   if (law.eq.1.and.lct.eq.1) go to 400
   if (law.eq.1.and.lct.eq.3.and.zap.gt.2004) go to 400
   if (law.eq.1.and.zap.eq.0.) go to 400
   if (law.eq.2.or.law.eq.3) go to 450
   if (law.eq.6) lang=0

   !--law 1 in the cm system.
   !--adaptive integration over lab energy.
   nl=1
   if (irec.gt.0) nl=2
   summ=0
   ebar=0
   dame=0
   x(2)=0
   call h6cm(x(2),epnext,term,nl,lang,lep,c)
   do l=1,nl
      y(2,l)=term(l)
   enddo
   el=x(2)*y(2,1)
   fl=0
  120 continue
   x(1)=epnext
   call h6cm(x(1),epnext,term,nl,lang,lep,c)
   do l=1,nl
      y(1,l)=term(l)
   enddo
   i=2
   dy=(y(1,1)+y(2,1))/10000
  140 continue
   if (i.eq.imax) go to 190
   da=(x(i-1)-x(i))*(y(i-1,1)+y(i,1))/2
   if (da.lt.small) go to 190
   xm=(x(i-1)+x(i))/2
   if (xm.ge.x(i-1).or.xm.le.x(i)) go to 190
   epn=x(i)
   call h6cm(xm,epn,term,nl,lang,lep,c)
   ym=(y(i-1,1)+y(i,1))/2
   test=2*tol*abs(ym)+2*dy
   if (abs(term(1)-ym).gt.test) go to 160
   go to 190
   ! fails test
  160 continue
   i=i+1
   x(i)=x(i-1)
   do l=1,nl
      y(i,l)=y(i-1,l)
   enddo
   x(i-1)=xm
   do l=1,nl
      y(i-1,l)=term(l)
   enddo
   go to 140
   ! passes test
  190 continue
   if (i.eq.1) go to 220
   summ=summ+(x(i-1)-x(i))*(y(i-1,1)+y(i,1))/2
   if (irec.le.0) then
      en=y(i-1,1)*x(i-1)
      ebar=ebar+(x(i-1)-x(i))*(en+el)/2
      f=y(i-1,1)*df(x(i-1),zp,ap,zt,at)
      dame=dame+(x(i-1)-x(i))*(f+fl)/2
   else
      f1=0
      if (y(i-1,1).ne.zero) f1=y(i-1,2)/y(i-1,1)
      en=(e-2*sqrt(e*awp*x(i-1))*f1+awp*x(i-1))/(awrt+1-awp)
      f=y(i-1,1)*df(en,zp,ap,zt,at)
      en=en*y(i-1,1)
      ebar=ebar+(x(i-1)-x(i))*(en+el)/2
      dame=dame+(x(i-1)-x(i))*(f+fl)/2
   endif
   el=en
   fl=f
   i=i-1
   if (i.gt.1) go to 140
   test=emax-emax/100
   if (epnext.ge.test) go to 190
   x(2)=x(1)
   do l=1,nl
      y(2,l)=y(1,l)
   enddo
   go to 120
  220 continue
   return

   !--law 1 for lab system.
   !--trapazoidal energy integration, and
   !--gauss-legendre angle integration, if needed.
  400 continue
   ebar=0
   dame=0
   nd=nint(c(3))
   na=nint(c(4))
   nep=nint(c(6))
   ncyc=na+2
   h=0
   d=0
   xl=0
   yl=0
   el=0
   fl=0
   do i=1,nep
      l=7+ncyc*(i-1)
      xx=c(l)
      yy=c(l+1)
      x2=xx
      if (irec.eq.0) then

   if (imode(3) .lt. 0) then
     write (nsyso,3329) i, xx, yy
3329 format (1x, 'GETSIX energy partial increment-A1 ', i6, 2g14.7)
   endif

         en=yy*xx
         zpp=int(zap/1000)
         ztt=int(zat/1000)
         f2=df(xx,zpp,awp,ztt,awrt)
         f=yy*f2
      else
         en=0
         f=0
         do iq=1,nq
            b=0
            call legndr(qp(iq),p,na)
            do ia=1,na+1
               b=b+(2*ia-1)*c(l+ia)*p(ia)/2
            enddo
            if (yy.ne.zero) b=b/yy
            er=(e-2*sqrt(e*awp*xx)*qp(iq)+awp*xx)/(awrt+1-awp)
            en=en+qw(iq)*b*er
            f=f+qw(iq)*b*df(er,zp,ap,zt,at)
         enddo
         x2=en

   if (imode(3) .lt. 0) then
     write (nsyso,2329) i, en, yy
2329 format (1x, 'GETSIX energy partial increment-A2 ', i6, 2g14.7)
   endif

         en=en*yy
         f2=f
         f=f*yy
      endif
      if (i.le.nd) then
         h=h+en
         d=d+f

   if (imode(3) .lt. 0) then
     write (nsyso,2309) i, nd, h, en, d, f
2309 format (1x, 'GETSIX energy increment-A ', 2i6, 4g14.7)
   endif

      else if (i.eq.nd+1) then
         xl=xx
         yl=yy
         el=en
         fl=f
      else
         if (lep.eq.1) t1=x2*yl+el
         if (lep.gt.1) t1=en+el

   if (imode(3) .lt. 0) then
     write (nsyso,2379) i, lep, t1, x2, yl, el, xx, xl, t1/2, (xx-xl)
2379 format (1x, 'GETSIX energy partial increment-B ', 2i6, 8g14.7)
   endif


         delta_h = (xx-xl)*t1/2
         h=h+(xx-xl)*t1/2
         if (lep.eq.1) t2=f2*yl+fl
         if (lep.gt.1) t2=f+fl
         delta_d = (xx-xl)*t2/2
         d=d+(xx-xl)*t2/2

   if (imode(3) .lt. 0) then
     write (nsyso,2319) i, nd, lep, h, delta_h, d, delta_d, xx, xl, t1
2319 format (1x, 'GETSIX energy increment-B ', 3i6, 7g14.7)
   endif

         xl=xx
         yl=yy
         el=en
         fl=f
      endif
   enddo
   ebar=h
   dame=d

   if (imode(3) .lt. 0) then
     write (nsyso,2305) e, ebar, dame, nd, na, nep, irec
2305 format (1x, 'GETSIX exit-400 ', 3g14.7, 4i6)
   endif

   return

   !--laws 2 and 3.
   !--discrete-level scattering
  450 continue
   ebar=0
   dame=0
   lang=0
   if (law.eq.2) lang=nint(c(3))
   nld=0
   if (law.eq.2) nld=nint(c(6))
   thresh=((awrt+1)/awrt)*(-q)
   if (e.lt.thresh) return
   beta=sqrt(awrt*(awrt+1-awp)*(1-thresh/e)/awp)
   afact=awp/(1+awrt)**2
   arec=awp/(awrt+1-awp)
   if (lang.eq.0) then
      ! legendre coefficients
      ! uses gauss-legendre integration for damage
      f1=0
      if (law.ne.3) f1=c(7)
      if (irec.gt.0) then
         f1=-f1
         beta=arec*beta
         afact=afact/arec
      endif
      ebar=afact*e*(1+2*beta*f1+beta*beta)
      if (idame.eq.0) return
      do iq=1,nq
         u=qp(iq)
         if (irec.gt.0) u=-u
         f=one/2
         if (law.ne.3) then
            call legndr(u,p,nld)
            do il=1,nld
               f=f+(2*il+1)*c(6+il)*p(il+1)/2
            enddo
         endif
         e2=afact*e*(1+2*beta*u+beta*beta)
         dame=dame+qw(iq)*f*df(e2,zp,ap,zt,at)
      enddo
   else
      ! tabulated distribution
      ! uses trapazoidal integration for heat and damage
      nmu=nint(c(6))
      if (irec.gt.0) then
         beta=arec*beta
         afact=afact/arec
      endif
      el=0
      fl=0
      ul=0
      do imu=1,nmu
         l=6+2*imu-1
         if (irec.gt.0) l=6+2*nmu-2*imu+1
         u=c(l)
         if (irec.gt.0) u=-u
         yy=c(l+1)
         e2=afact*e*(1+2*beta*u+beta*beta)
         en=yy*e2
         fn=yy*df(e2,zp,ap,zt,at)
         if (imu.gt.1) then
            ebar=ebar+(u-ul)*(en+el)/2
            dame=dame+(u-ul)*(fn+fl)/2
         endif
         ul=u
         el=en
         fl=fn
      enddo
   endif
   return

   !--law 7.
   !--trapazoidal integration over energy for each angle
   !--then trapazoidal integration over angle
  510 continue
   ebar=0
   dame=0
   dl=0
   hl=0
   ul=0
   nmu=nint(c(6))
   l=7+2*nint(c(5))
   do imu=1,nmu
      u=c(l+1)
      nr=nint(c(l+4))
      np=nint(c(l+5))
      iint=nint(c(l+7))
      next=l+6+2*nr+2*np
      ibase=l+5+2*nr
      h=0
      d=0
      do i=1,np
         xx=c(ibase+2*i-1)
         yy=c(ibase+2*i)
         en=yy*xx
         if (irec.gt.0) then
            xx=(e-2*sqrt(e*awp*xx)*u+awp*xx)/(awrt+1-awp)
         endif
         if (i.gt.1) then
            if (iint.eq.1) then
               h=h+(xx-xl)*el
            else
               h=h+(xx-xl)*(en+el)/2
            endif
         endif
         f=yy*df(xx,zp,ap,zt,at)
         if (i.gt.1) then
            if (iint.eq.1) then
               d=d+(xx-xl)*fl
            else
               d=d+(xx-xl)*(f+fl)/2
            endif
         endif
         xl=xx
         el=en
         fl=f
      enddo
      l=next
      if (imu.gt.1) ebar=ebar+(u-ul)*(h+hl)/2
      if (imu.gt.1) dame=dame+(u-ul)*(d+dl)/2
      ul=u
      hl=h
      dl=d
   enddo

   if (imode(3) .lt. 0) then
     write (nsyso,2302) e, ebar, dame, law
2302 format (1x, 'GETSIX exit ', 3g14.7, i6)
   endif

   return
   end subroutine getsix

   subroutine h6cm(ep,epnext,term,nl,lang,lep,cnow)
   !-------------------------------------------------------------------
   ! Compute the Legendre coefficients for the double-differential
   ! cross section for secondary energy ep (in the lab system)
   ! using the ENDF-6 File 6 data in the cm system located
   ! in the array cnow.  Call with ep=0 to initialize.
   ! Thereafter, ep can be requested in any order.
   !-------------------------------------------------------------------
   use mathm ! provides legndr
   ! externals
   integer::nl,lang,lep
   real(kr)::ep,epnext,term(*),cnow(*)
   ! internals
   integer::i,na,l,j,ll,ndnow,npnow,ncnow
   real(kr)::aprime,epnn,epm,epn,w,t,epx,xx,cc,c,umin,dm
   real(kr)::umax,un,u,yy,test,epp,s,dx,ym,epnxt
   real(kr)::wlo,whi,f,eps,xc,elmax,e,epmax
   real(kr)::da
   integer,parameter::imax=10
   integer,parameter::mxlg=65
   real(kr)::p(mxlg)
   real(kr)::x(imax),y(imax,mxlg),yt(mxlg)
   real(kr),parameter::tol=0.02e0_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::step=0.05e0_kr
   real(kr),parameter::check=.99999e0_kr
   real(kr),parameter::rdn=.999995e0_kr
   real(kr),parameter::rup=1.000005e0_kr
   real(kr),parameter::tiny=1.e-9_kr
   real(kr),parameter::dn=.999999e0_kr
   real(kr),parameter::up=1.000001e0_kr
   real(kr),parameter::thou=1.e-3_kr
   real(kr),parameter::amil=1.e-6_kr
   real(kr),parameter::one=1.e0_kr
   real(kr),parameter::onep5=1.5e0_kr
   real(kr),parameter::du=.0001e0_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0
   save na,eps,xc,ndnow,npnow,ncnow,elmax,e,epmax

   !--initialize.
   if (ep.gt.zero) go to 200
   e=cnow(2)
   ndnow=nint(cnow(3))
   npnow=nint(cnow(6))
   ncnow=nint(cnow(5))/npnow
   aprime=awp
   xc=aprime/(awrt+1)**2
   epnext=emax
   elmax=0
   epnn=emax
   epmax=0
   epm=0
   epn=0
   w=0
   if (npnow.ne.ndnow) then
      if (lang.gt.0) then
         t=h6ddx(ep,epn,epm,w,cnow,lang,lep)
      else
         t=h6psp(ep,epn,epm,w,e,cnow)
      endif
      if (epn.lt.epnn*(1-small)) epnn=epn
      if (epm.gt.epmax*(1+small)) epmax=epm
      elmax=e*(sqrt(epmax/e)+sqrt(xc))**2
      elmax=up*elmax
      epn=step*elmax
      if(epmax*up.lt.xc*e) epn=e*(sqrt(epmax/e)-sqrt(xc))**2
      if (epn.lt.epnext) epnext=epn
      eps=tiny/elmax
      epnext=dn*epnext
   endif
   if (ndnow.ne.0) then
      i=0
      t=h6dis(i,epn,epm,w,cnow,lang)
      do i=1,ndnow
         t=h6dis(i,epn,epm,w,cnow,lang)
         if (epn.lt.epnn*(1-small)) epnn=epn
         if (epn.gt.epmax*(1+small)) epmax=epn
         epx=dn*e*(sqrt(epn/e)-sqrt(xc))**2
         if (epx.lt.epnext*(1-small)) epnext=epx
         epx=up*e*(sqrt(epn/e)+sqrt(xc))**2
         if (epx.gt.elmax*(1+small)) elmax=epx
      enddo
   endif
   na=nl-1
   do l=1,nl
      term(l)=0
   enddo
   return

   !--normal entry.
  200 continue
   do l=1,nl
      term(l)=0
   enddo
   if (ndnow.eq.npnow) go to 330
   xx=ep/e
   cc=xc/xx
   c=sqrt(cc)
   umin=(1+cc-epmax/ep)/(2*c)
   if (umin.lt.check) go to 215
   if (c.le.one) go to 420
   if (c.gt.one) go to 320
  215 continue
   if (umin.lt.-one) umin=-one
   dm=(1-umin)/1000+amil
   if (dm.lt.thou/5) dm=thou/5
   epnn=emax

   !--for the continuum part of the distribution,
   !--integrate over the trajectory ep=constant
   !--using adaptive integration on lab cosine.
   umax=1
   un=umax
   x(2)=un
  220 continue
   if (x(2).le.umin) go to 310
   j=2
   do while (j.eq.2)
      u=un
      yy=1+cc-2*c*u
      test=amil
      if (yy.lt.test) then
         yy=amil
         c=u-thou
         cc=c*c
      endif
      epp=yy*ep
      w=(u-c)/sqrt(yy)
      if (lang.gt.0) then
         s=h6ddx(epp,epn,epm,w,cnow,lang,lep)
      else
         s=h6psp(epp,epn,epm,w,e,cnow)
      endif
      if (u.eq.umax.and.epn.lt.epnn) epnn=epn
      call legndr(u,p,na)
      j=1
      if (u.eq.umax) j=2
      x(j)=u
      do l=1,nl
         y(j,l)=p(l)*s/sqrt(yy)
      enddo
      un=u-(epn-epp)/(2*c*ep)
      if (un.gt.u-du) un=u-du
      if (un.lt.umin+du/10) un=umin
   enddo
   i=2
  240 continue
   if (i.eq.imax) go to 290
   dx=x(i)-x(i-1)
   if (dx.lt.dm) go to 290
   da=dx*(y(i,1)+y(i-1,1))/2
   if (da.lt.eps) go to 290
   u=(x(i-1)+x(i))/2
   yy=1+cc-2*c*u
   epp=yy*ep
   w=(u-c)/sqrt(yy)
   if (lang.gt.0) s=h6ddx(epp,epn,epm,w,cnow,lang,lep)
   if (lang.eq.0) s=h6psp(epp,epn,epm,w,e,cnow)
   call legndr(u,p,na)
   do l=1,nl
      yt(l)=p(l)*s/sqrt(yy)
   enddo
   do 250 l=1,nl
   ym=(y(i-1,l)+y(i,l))/2
   if (abs(yt(l)-ym).gt.l*tol*abs(ym)+eps) go to 260
  250 continue
   go to 290
   ! fails test
  260 continue
   i=i+1
   x(i)=x(i-1)
   do l=1,nl
      y(i,l)=y(i-1,l)
   enddo
   x(i-1)=u
   do l=1,nl
      y(i-1,l)=yt(l)
   enddo
   go to 240
   ! passes test
  290 continue
   do l=1,nl
      term(l)=term(l)+(x(i)-x(i-1))*(y(i,l)+y(i-1,l))/2
   enddo
   i=i-1
   if (i.gt.1) go to 240
   x(2)=x(1)
   do l=1,nl
      y(2,l)=y(1,l)
   enddo
   if (i.eq.1) go to 220
  310 continue
   epnxt=e*(sqrt(epnn/e)-sqrt(xc))**2
   if (epnxt.le.ep*(1+small)) epnxt=e*(sqrt(epnn/e)+sqrt(xc))**2
   if (epnxt.le.ep*(1+small)) epnxt=onep5*ep
   if (epnxt.gt.elmax*(1+small).and.elmax.gt.ep*(1+small))epnxt=elmax
   go to 400
  320 continue
   epnxt=(epnext+elmax)/2
   go to 400

   !--check for contributions from delta functions
  330 continue
   if (ndnow.eq.0) go to 400
   epnxt=emax
   do i=1,ndnow
      ll=7+ncnow*(i-1)
      epp=cnow(ll)
      w=(ep-epp-xc*e)/(2*sqrt(xc*e)*sqrt(epp))
      u=(xc*e+ep-epp)/(2*sqrt(xc*e)*sqrt(ep))
      wlo=-1
      whi=1
      if (w.ge.wlo.and.w.le.whi) then
         call legndr(u,p,na)
         f=h6dis(i,epn,epm,w,cnow,lang)
         f=f/(2*sqrt(xc*e*epp))
         epn=e*(sqrt(epp/e)+sqrt(xc))**2
         if (ep.lt.check*epn) then
            epn=rdn*epn
         else
            epn=rup*epn
         endif
         if (epn.lt.epnxt) epnxt=epn
         do l=1,nl
            term(l)=term(l)+f*p(l)
         enddo
      else if (w.lt.wlo) then
         epn=e*(sqrt(epp/e)-sqrt(xc))**2
         if (ep.lt.check*epn) then
            epn=rdn*epn
         else
            epn=rup*epn
         endif
         if (epn.lt.epnxt) epnxt=epn
      endif
   enddo

   !--select epnext
  400 continue
   epnext=epnxt
   return
  420 continue
   epnext=emax
   return
   end subroutine h6cm

   real(kr) function h6ddx(ep,epnext,epmax,w,cnow,lang,lep)
   !-------------------------------------------------------------------
   ! Retrieve the double differential cross section for
   ! secondary energy ep and cosine w (in the cm system) from a
   ! subsection in ENDF-6 File 6 format stored in cnow.  Call with
   ! ep=0 to initialize for each new e. Thereafter, ep can be
   ! requested in any order.  Only the continuum part is
   ! returned.  Delta functions, if given, must be handled
   ! separately.
   !-------------------------------------------------------------------
   use mathm ! provides legndr
   use util  ! provides error
   use endf  ! provides terp1
   ! externals
   integer::lang,lep
   real(kr)::ep,epnext,epmax,w,cnow(*)
   ! internals
   integer::ndnow,npnow,ncnow,mnow,inow,lnow,nl,na,l,lll
   integer::iza2,inn,ii,jj,ia,illdef
   real(kr)::enow,efirst,eplast,t,tt,s,r,aa,tii,tjj
   real(kr)::x1,x2,y1,y2
   integer,parameter::nlmax=65
   real(kr)::p(nlmax)
   real(kr),parameter::up=1.00001e0_kr
   real(kr),parameter::dn=.99999e0_kr
   real(kr),parameter::off=.999995e0_kr
   real(kr),parameter::step=0.05e0_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0
   save efirst,enow,nl,inow,mnow,lnow,ncnow,na,illdef

   !--initialize
   if (ep.eq.zero) then
      enow=cnow(2)
      ndnow=nint(cnow(3))
      npnow=nint(cnow(6))
      ncnow=nint(cnow(5))/npnow
      mnow=7+(npnow-1)*ncnow
      inow=7+ncnow*ndnow
      lnow=inow
      if (cnow(lnow).le.zero) lnow=lnow+ncnow
      epnext=cnow(lnow)
      epmax=cnow(mnow)
      nl=ncnow-1
      if (nl.gt.nlmax) call error('h6ddx',&
        'too many legendre terms',' ')
      na=nl-1
      efirst=0
      h6ddx=0
      illdef=0
      return
   endif

   !--check for desired interpolation range
   if (ep.lt.efirst*(1-small)) go to 200
   if (ep.gt.cnow(mnow)) go to 190
  130 continue
   epnext=cnow(lnow)
   if (ep.lt.off*epnext) go to 140
   if (lnow.ge.mnow) go to 135
   lnow=lnow+ncnow
   go to 130
  135 continue
   go to 150
  140 continue
   if (lnow.le.ncnow) go to 200
   eplast=cnow(lnow-ncnow)
   if (ep.ge.off*eplast) go to 150
   if (lnow.le.inow+ncnow) go to 150
   lnow=lnow-ncnow
   go to 130
  150 continue
   go to 210

   !--return zero above table
  190 continue
   t=0
   epnext=emax
   go to 350

   !--use analytic form below efirst.
  200 continue
   t=0
   go to 350

   !--legendre coefficients
  210 continue
   if (lang.eq.1) then
      call legndr(w,p,na)
      t=0
      do l=1,nl
         tt=0
         if (l.le.ncnow-1&
             .and.ep.ge.cnow(inow)*(1-small)&
             .and.ep.le.cnow(mnow)*(1+small)) then
            lll=lnow+l
            x1=cnow(lnow-ncnow)
            x2=cnow(lnow)
            if (x1.eq.x2.and.lep.gt.1) then
               x2=sigfig(x2,6,1)
               if (illdef.eq.0) then
                  illdef=1
               endif
            endif
            y1=cnow(lll-ncnow)
            y2=cnow(lll)
            call terp1(x1,y1,x2,y2,ep,tt,lep)
         endif
         t=t+p(l)*(2*l-1)*tt/2
      enddo
      if (t.lt.zero) t=0

   !--kalbach systematics
   else if (lang.eq.2) then
      s=0
      r=0
      if (ep.ge.cnow(inow)*(1-small).and.&
          ep.le.cnow(mnow)*(1+small)) then
         x1=cnow(lnow-ncnow)
         x2=cnow(lnow)
         if (x1.eq.x2.and.lep.gt.1) then
            x2=sigfig(x2,6,1)
            if (illdef.eq.0) then
               illdef=1
            endif
         endif
         y1=cnow(lnow-ncnow+1)
         y2=cnow(lnow+1)
         call terp1(x1,y1,x2,y2,ep,s,lep)
         y1=cnow(lnow-ncnow+2)
         y2=cnow(lnow+2)
         call terp1(x1,y1,x2,y2,ep,r,lep)
      endif
      iza2=nint(zap)
      aa=bacha(izap,iza2,izat,enow,ep)
      t=aa*(cosh(aa*w)+r*sinh(aa*w))/(2*sinh(aa))
      t=t*s
      if (t.lt.zero) t=0

   !--tabulated angular distribution
   else if (lang.ge.11.and.lang.le.15) then
      inn=lang-10
      ii=lnow-ncnow
      jj=lnow
      x1=cnow(ii)
      x2=cnow(jj)
      if (x1.eq.x2.and.lep.gt.1) then
         x2=sigfig(x2,6,1)
         if (illdef.eq.0) then
            illdef=1
         endif
      endif
      y1=cnow(ii+1)
      y2=cnow(jj+1)
      call terp1(x1,y1,x2,y2,ep,s,lep)
      tii=0
      do ia=1,na,2
         if (w.ge.cnow(ii+1+ia).and.w.le.cnow(ii+3+ia)) then
            call terp1(cnow(ii+1+ia),cnow(ii+2+ia),cnow(ii+3+ia),&
              cnow(ii+4+ia),w,tii,inn)
         endif
      enddo
      tjj=0
      do ia=1,na,2
         if (w.ge.cnow(jj+1+ia).and.w.le.cnow(jj+3+ia)) then
            call terp1(cnow(jj+1+ia),cnow(jj+2+ia),cnow(jj+3+ia),&
              cnow(jj+4+ia),w,tjj,inn)
         endif
      enddo
      call terp1(x1,tii,x2,tjj,ep,t,lep)
      t=t*s

   !--bad value for lang
   else
      call error('h6ddx','illegal lang.',' ')
   endif

   !--return final value
  350 continue
   h6ddx=t
   if (lep.gt.1) return
   if (epnext.eq.emax) return
   if (ep.lt.dn*dn*epnext) then
      epnext=dn*epnext
   else
      epnext=up*epnext
   endif
   return
   end function h6ddx

   real(kr) function h6dis(i,epnext,epmax,w,cnow,lang)
   !-------------------------------------------------------------------
   ! retrieve the double differential cross section for
   ! discrete-energy i and cosine w (in the cm system) from a
   ! subsection in endf-6 file 6 format stored in cnow.  call with
   ! i=0 to initialize for each new e. thereafter, i can be
   ! requested in any order.
   !-------------------------------------------------------------------
   use util  ! provides error
   use mathm ! provides legndr
   use endf  ! provides terp1
   ! externals
   integer::i,lang
   real(kr)::epnext,epmax,w,cnow(*)
   ! internals
   integer::ndnow,npnow,ncnow,inow,lnow,mnow,nl,na,l,iza2,inn,ia
   integer,parameter::mxlg=65
   real(kr)::enow,t,s,r,aa
   real(kr)::p(mxlg)
   real(kr),parameter::zero=0
   save nl,inow,lnow,mnow,ncnow,na,enow

   !--initialize
   if (i.eq.0) then
      enow=cnow(2)
      ndnow=nint(cnow(3))
      npnow=nint(cnow(6))
      ncnow=nint(cnow(5))/npnow
      inow=7
      lnow=inow
      mnow=lnow+(ndnow-1)*ncnow
      epnext=cnow(lnow)
      epmax=cnow(mnow)
      nl=ncnow-1
      na=nl-1
      h6dis=0
      return
   endif

   !--select desired discrete line
   inow=7+(i-1)*ncnow
   epnext=cnow(inow)

   !--legendre coefficients
   if (lang.eq.1) then
      call legndr(w,p,na)
      t=0
      do l=1,nl
         if (l.le.ncnow-1) then
            t=t+(2*l-1)*p(l)*cnow(inow+1+l)/2
         endif
      enddo
      t=t*cnow(inow+1)
      if (t.lt.zero) t=0

   !--kalbach systematics
   else if (lang.eq.2) then
      s=cnow(inow+1)
      r=cnow(inow+2)
      iza2=nint(zap)
      aa=bacha(izap,iza2,izat,enow,epnext)
      t=aa*(cosh(aa*w)+r*sinh(aa*w))/(2*sinh(aa))
      t=t*s
      if (t.lt.zero) t=0

   !--tabulated angular distribution
   else if (lang.ge.11.and.lang.le.15) then
      inn=lang-10
      t=0
      do ia=1,na-2,2
         if (w.ge.cnow(inow+1+ia).and.w.le.cnow(inow+3+ia)) then
            call terp1(cnow(inow+1+ia),cnow(inow+2+ia),&
              cnow(inow+3+ia),cnow(inow+4+ia),w,t,inn)
         endif
      enddo
      t=t*cnow(inow+1)

   !--illegal lang
   else
      call error('h6dis','illegal lang.',' ')
   endif

   !--return final value
   h6dis=t
   return
   end function h6dis

   real(kr) function bacha(iza1i,iza2,izat,e,ep)
   !-------------------------------------------------------------------
   ! Compute the Kalbach a parameter.
   !-------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::iza1i,iza2,izat
   real(kr)::e,ep
   ! internals
   integer::iza1,iza
   real(kr)::aa,za,ac,zc,ab,zb,nc,nb,na,sa,sb
   real(kr)::ecm,ea,eb,x1,x3,fa,fb
   character(60)::strng
   real(kr),parameter::third=.333333333e0_kr
   real(kr),parameter::twoth=.666666667e0_kr
   real(kr),parameter::fourth=1.33333333e0_kr
   real(kr),parameter::c1=15.68e0_kr
   real(kr),parameter::c2=-28.07e0_kr
   real(kr),parameter::c3=-18.56e0_kr
   real(kr),parameter::c4=33.22e0_kr
   real(kr),parameter::c5=-0.717e0_kr
   real(kr),parameter::c6=1.211e0_kr
   real(kr),parameter::s2=2.22e0_kr
   real(kr),parameter::s3=8.48e0_kr
   real(kr),parameter::s4=7.72e0_kr
   real(kr),parameter::s5=28.3e0_kr
   real(kr),parameter::brk1=130.e0_kr
   real(kr),parameter::brk2=41.e0_kr
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::b1=0.04e0_kr
   real(kr),parameter::b2=1.8e-6_kr
   real(kr),parameter::b3=6.7e-7_kr
   real(kr),parameter::tomev=1.e-6_kr
   real(kr),parameter::zero=0

   iza1=iza1i
   if (iza1i.eq.0) iza1=1
   iza=izat
   if (iza.eq.6000) iza=6012
   if (iza.eq.12000) iza=12024
   if (iza.eq.14000) iza=14028
   if (iza.eq.16000) iza=16032
   if (iza.eq.17000) iza=17035
   if (iza.eq.19000) iza=19039
   if (iza.eq.20000) iza=20040
   if (iza.eq.22000) iza=22048
   if (iza.eq.23000) iza=23051
   if (iza.eq.24000) iza=24052
   if (iza.eq.26000) iza=26056
   if (iza.eq.28000) iza=28058
   if (iza.eq.29000) iza=29063
   if (iza.eq.31000) iza=31069
   if (iza.eq.40000) iza=40090
   if (iza.eq.42000) iza=42096
   if (iza.eq.48000) iza=48112
   if (iza.eq.49000) iza=49115
   if (iza.eq.50000) iza=50120
   if (iza.eq.63000) iza=63151
   if (iza.eq.72000) iza=72178
   if (iza.eq.74000) iza=74184
   if (iza.eq.82000) iza=82208
   aa=mod(iza,1000)
   if (aa.eq.zero) then
      write(strng,'(''dominant isotope not known for '',i8)') iza
      call error('bacha',strng,' ')
   endif
   za=int(iza/1000)
   ac=aa+mod(iza1,1000)
   zc=za+int(iza1/1000)
   ab=ac-mod(iza2,1000)
   zb=zc-int(iza2/1000)
   na=nint(aa-za)
   nb=nint(ab-zb)
   nc=nint(ac-zc)
   sa=c1*(ac-aa)&
     +c2*((nc-zc)**2/ac-(na-za)**2/aa)&
     +c3*(ac**twoth-aa**twoth)&
     +c4*((nc-zc)**2/ac**fourth-(na-za)**2/aa**fourth)&
     +c5*(zc**2/ac**third-za**2/aa**third)&
     +c6*(zc**2/ac-za**2/aa)
   if (iza1.eq.1002) sa=sa-s2
   if (iza1.eq.1003) sa=sa-s3
   if (iza1.eq.2003) sa=sa-s4
   if (iza1.eq.2004) sa=sa-s5
   sb=c1*(ac-ab)&
     +c2*((nc-zc)**2/ac-(nb-zb)**2/ab)&
     +c3*(ac**twoth-ab**twoth)&
     +c4*((nc-zc)**2/ac**fourth-(nb-zb)**2/ab**fourth)&
     +c5*(zc**2/ac**third-zb**2/ab**third)&
     +c6*(zc**2/ac-zb**2/ab)
   if (iza2.eq.1002) sb=sb-s2
   if (iza2.eq.1003) sb=sb-s3
   if (iza2.eq.2003) sb=sb-s4
   if (iza2.eq.2004) sb=sb-s5
   ecm=aa*e/ac
   ea=ecm*tomev+sa
   eb=tomev*ep*ac/ab+sb
   x1=eb
   if (ea.gt.brk1) x1=brk1*eb/ea
   x3=eb
   if (ea.gt.brk2) x3=brk2*eb/ea
   fa=1
   if (iza1.eq.2004) fa=0
   fb=1
   if (iza2.eq.1) fb=half
   if (iza2.eq.2004) fb=2
   bacha=b1*x1+b2*x1**3+b3*fa*fb*x3**4
   return
   end function bacha

   real(kr) function h6psp(ep,epnext,epmax,w,e,c)
   !-------------------------------------------------------------------
   ! Compute the double differential cross section for this incident
   ! energy (in the lab system), and secondary energy ep and cosine
   ! w (in the CM system) using the n-body phase-space formula.
   ! Call with ep=0 to initialize the retrieval routine for each
   ! new e.  Thereafter, ep can be requested in any order.
   !-------------------------------------------------------------------
   use endf ! provides contio
   use util ! provides error
   ! externals
   real(kr)::ep,epnext,epmax,w,e,c(*)
   ! internals
   integer::npsx
   real(kr)::apsx,f1,f2,test,s,cn,ex,eimax
   real(kr),parameter::c3=1.2732e0_kr
   real(kr),parameter::c4=3.2813e0_kr
   real(kr),parameter::c5=5.8205e0_kr
   real(kr),parameter::thrhaf=1.5e0_kr
   real(kr),parameter::sevhaf=3.5e0_kr
   real(kr),parameter::step=0.05e0_kr
   real(kr),parameter::rndup=1.000001e0_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::zero=0
   save cn,ex,eimax

   !--initialize
   if (ep.eq.zero) then
      apsx=c(1)
      npsx=nint(c(6))
      f1=(apsx-awp)/apsx
      f2=awrt/(awrt+1)
      ex=thrhaf*npsx-4
      eimax=f1*(f2*e+q)
      if (eimax.le.zero) eimax=1
      if (npsx.eq.3) then
         cn=c3/eimax**2
      else if (npsx.eq.4) then
         cn=c4/eimax**sevhaf
      else if (npsx.eq.5) then
         cn=c5/eimax**5
      else
         call error('h6psp','3, 4, or 5 particles only.',' ')
      endif
      h6psp=0
      epnext=step*eimax
      test=1
      if (eimax.eq.test) epnext=emax
      epmax=eimax
      return
   endif

   !--compute cross section
   s=0
   if (ep.lt.eimax*(1-small)) s=cn*sqrt(ep)*(eimax-ep)**ex
   h6psp=s
   epnext=ep+step*eimax
   if (epnext.gt.eimax*(1+small)) epnext=rndup*eimax
   if (ep.ge.eimax*(1-small)) epnext=emax
   return
   end function h6psp

   subroutine tabsq6(e_incident, g,h,a,law,z,awr)
   !-------------------------------------------------------------------
   ! Compute average of photon recoil energy from capture and
   ! corresponding damage energy for File 6 capture photons
   !-------------------------------------------------------------------
   use endf ! provides terp1
   use mainio  ! provides nsysi,nsyso
   use snl     ! provides SNL
   ! externals
   integer::law
   real(kr)::g,h,a(*),z,awr
   real(kr):: e_incident
   ! internals
   integer::nd,np,ncyc,ibase,inn,nc,i,j
   real(kr)::ein,rein,x,y,xr,awc,xh,yh,xl,yl,dx,s
   real(kr):: g_prime, h_prime, max_intensity, xr_prime
   integer,parameter::nq=4
   integer:: ichange_flag
   real(kr):: total_discrete_gamma_energy, total_continuum_gamma_energy, total_gamma_energy
   real(kr):: partial_sum, discrete_sum, s_discrete, s_continuum
   real(kr):: g_prime_discrete, h_prime_discrete, g_prime_continuum, h_prime_continuum
   real(kr),dimension(nq),parameter::qp=(/&
     -.86114e0_kr,-.33998e0_kr,.33998e0_kr,.86114e0_kr/)
   real(kr),dimension(nq),parameter::qw=(/&
     .34785e0_kr,.65215e0_kr,.65215e0_kr,.34785e0_kr/)
   data ichange_flag /0/

   if (imode(3) .lt. 0) then
     write (nsyso,2301) g, h, a(1), law, z, awr, pe_archive, pe_store
2301 format (1x, 'TABSQ6 entry ', 3g14.7, i6, 4g14.7)
   endif

   !--initialize
   ein=2*tm
   rein=1/ein
   g=0
   g_prime = 0
   h_prime = 0
   h=0
   s=0
   nd=nint(a(3))
   np=nint(a(6))
   ncyc=nint(a(5))/np
   ibase=6
   inn=2

   !--accumulate contributions from discrete levels
   total_gamma_energy = 0.0
   total_discrete_gamma_energy = 0.0
   total_continuum_gamma_energy = 0.0
   if (nd.ne.0) then

!    precompute discrete intensity sum and maximum intensity
     discrete_sum = 0.0
     max_intensity = 0.0
     do i = 1,nd
       discrete_sum = discrete_sum + a(ibase+ncyc*(i-1)+2)
       max_intensity = max(max_intensity,a(ibase+ncyc*(i-1)+2))
     enddo

      do i=1,nd
         x=a(ibase+ncyc*(i-1)+1)
         y=a(ibase+ncyc*(i-1)+2)
         xr=x*x*rein
         xr_prime = xr
!       
!        over-ride xr to include base capture recoil energy
!        xr_prime = xr + e_incident/(awr+1.)
!
         g=g+xr*y
         g_prime = g_prime + xr_prime*y
         awc=awr+1
         h=h+df(xr,z,awc,z,awr)*y
         h_prime = h_prime + df(xr_prime,z,awc,z,awr)*y
         s=s+y

         if (imode(3) .lt. 0) then
           write (nsyso,2306) i, x, y, xr, rein, g, h, s, g_prime, h_prime, discrete_sum,  &
              max_intensity, pe_archive, pe_store
2306       format (1x, 'TABSQ6 discrete level entry ', i5, 13g14.7)
         endif
         total_discrete_gamma_energy = total_discrete_gamma_energy + x*y

      enddo
   endif

   if (imode(3) .lt. 0) then
     write (nsyso,2303) g, h, s, total_discrete_gamma_energy, nd, np
2303 format (1x, 'TABSQ6 discrete level sum ', 4g14.7, 2i6)
   endif
   g_prime_discrete = g_prime
   h_prime_discrete = h_prime
   s_discrete = s
   g_prime_continuum = 0.0
   h_prime_continuum = 0.0
   s_continuum = 0.0

   !--accumulate contributions from continuum
   if (np.gt.nd) then
      nc=np-nd
      ibase=ibase+ncyc*nd
      xh=a(ibase+1)
      yh=a(ibase+2)
      do i=2,nc
         xl=xh
         xh=a(ibase+ncyc*(i-1)+1)
         yl=yh
         yh=a(ibase+ncyc*(i-1)+2)
         partial_sum = 0.0
         if (xl.ne.xh) then
            dx=xh-xl
            do j=1,nq
               x=xl+(1+qp(j))*dx/2
               call terp1(xl,yl,xh,yh,x,y,inn)
               xr=x*x*rein
               g=g+qw(j)*xr*dx*y
               h=h+qw(j)*df(xr,z,awr+1,z,awr)*dx*y

!              Note: proper g_prime may not be computed by baseline appraoch.
!                    Tally discrete and continuum g/h contributions separately
!                      since they scale slightly differently.

               xr_prime = xr
!       
!              over-ride xr to include base capture recoil energy
!              xr_prime = xr + e_incident/(awr+1.)
!
               g_prime = g_prime + qw(j)*xr_prime*dx*y
               h_prime = h_prime + qw(j)*df(xr_prime,z,awr+1,z,awr)*dx*y

               g_prime_continuum = g_prime_continuum + qw(j)*xr_prime*dx*y
               h_prime_continuum = h_prime_continuum + qw(j)*df(xr_prime,z,awr+1,z,awr)*dx*y

               s_continuum=s_continuum+qw(j)*dx*y

               s=s+qw(j)*dx*y
               partial_sum = partial_sum + qw(j)*dx*y

               if (imode(3) .lt. 0 ) then
!                   write (nsyso,8336) j, g, g_prime, g_prime_continuum, h, h_prime, & 
!                        h_prime_continuum, qw(j), xr, xr_prime, &
!                        dx, y, s, s_continuum
8336               format (1x, 'TABSQ6 continuum level g/h_prime/continuum loop ', i5, 13g14.7)
               endif

            enddo
            total_continuum_gamma_energy = total_continuum_gamma_energy + partial_sum
         endif

         if (imode(3) .lt. 0 ) then
            write (nsyso,2336) i, xl, xh, yl, yh, g, g_prime, g_prime_continuum, h, h_prime, & 
                  h_prime_continuum, s, &
                  total_continuum_gamma_energy, partial_sum
2336        format (1x, 'TABSQ6 continuum level entry ', i5, 13g14.7)
         endif

      enddo
   endif
   
   total_gamma_energy = total_discrete_gamma_energy + total_continuum_gamma_energy
   if (imode(3) .lt. 0) then
     write (nsyso,2304) g, g_prime, h, h_prime, s, total_gamma_energy, nc
2304 format (1x, 'TABSQ6 discrete+continuum level sum ', 6g14.7, i6)
     write (nsyso, 2314) g_prime_discrete, h_prime_discrete, g_prime_continuum, h_prime_continuum, s, &
        s_discrete, s_continuum, pe_archive
2314 format (1x, 'TABSQ6 altered sum: ', 8g14.7)
   endif

   !--finished
   g=g/s
   g_prime = g_prime_discrete*pe_archive + g_prime_continuum*pe_archive/s
!   g_prime = g_prime_discrete*pe_archive + g_prime_continuum
   h=h/s
   h_prime = h_prime_discrete*pe_archive + h_prime_continuum*pe_archive/s
!   h_prime = h_prime_discrete*pe_archive + h_prime_continuum

   if (imode(3) .lt. 0 .and. ichange_flag .lt. 10) then
     write (nsyso,8702) ichange_flag, g_prime_discrete*pe_archive, g_prime_continuum*pe_archive/s, &
         h_prime_discrete*pe_archive, &
         h_prime_continuum*pe_archive/s, s, pe_archive
8702 format (1x, 'Note: TABSQ6 disc/cont partition of recoil and damage energies: ', i6, 9g14.7)
   endif

!  code modifications to correct multiphoton MF6 kerma calculation
!  change by PJG 9/19/2020

   if (imode(3) .lt. 0 .and. ichange_flag .lt. 10) then
     ichange_flag = ichange_flag + 1
     write (nsyso,2702) ichange_flag, g, g_prime, h, h_prime, s, pe_archive, pe_store, xr, xr_prime
2702 format (1x, 'Warning: code modification taking affect for neutron kerma ', i6, 9g14.7)
   endif
   g = g_prime
   h = h_prime

!  end code modification section

   if (imode(3) .lt. 0) then
     write (nsyso,2302) g, g_prime, h, h_prime, s, xr, xr_prime
2302 format (1x, 'TABSQ6 exit ', 9g14.7)
   endif

   return
   end subroutine tabsq6

   subroutine hgtfle(e,enext,idis,fle,nle,lcd,matd,mfd,mtd,nin,c,nc)
   !-------------------------------------------------------------------
   ! Retrieve or compute Legendre coefficients at e.
   ! Return isotropic distribution below first point.
   !-------------------------------------------------------------------
   use util ! provides error,mess
   use endf ! provides endf routines
   ! externals
   integer::idis,nle,lcd,matd,mfd,mtd,nin,nc
   real(kr)::e,enext,fle(*),c(nc)
   ! internals
   integer::nb,iso,lvt,ltt,ltt3,lttn,nwc,i,nlmax,il,l
   integer::iraw,ir,nne,ne,inn,innt,lct,nlo,nhi
   real(kr)::test,awr,elo,ehi
   character(60)::strng
   real(kr)::flo(65),fhi(65)
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::zero=0
   save iso,iraw,ir,nne,ne,inn
   save elo,ehi,nlo,nhi,flo,fhi,ltt,ltt3,lttn
   save lct,awr

   !--initialize
   idis=0
   if (e.gt.zero) go to 110
   iso=0
   call findf(matd,mfd,mtd,nin)
   call contio(nin,0,0,c,nb,nwc)
   awr=c(2)
   lvt=nint(c(3))
   ltt=nint(c(4))
   ltt3=ltt
   if (ltt.eq.3) then
      ltt=1
      lttn=1
   endif
   if (lvt.eq.0) then
      call contio(nin,0,0,c,nb,nwc)
   else
      call listio(nin,0,0,c,nb,nwc)
      do while (nb.ne.0)
         call moreio(nin,0,0,c(7),nb,nwc)
      enddo
   endif
   iso=nint(c(3))
   lct=nint(c(4))
   ! check for lab distribution for two-body scattering.
   test=10
   if (lct.eq.1.and.awr.ge.test) then
      if (mtd.ge.51.and.mtd.le.90) then
         write(strng,&
           '(''lab distribution changed to cm for mt '',i3)') mtd
         call mess('hgtfle',strng,' ')
         lct=2
      endif
   endif
   if (iso.ne.1) then
      call tab2io(nin,0,0,c(1),nb,nwc)
      ne=nint(c(6))
      ! read in raw data at first incident energy point
      ! retrieve or compute legendre coefficients.
      iraw=1+nwc
      if (ltt.eq.1) call listio(nin,0,0,c(iraw),nb,nwc)
      if (ltt.eq.2) call tab1io(nin,0,0,c(iraw),nb,nwc)
      if (nb.ne.0) then
         l=iraw
         do while (nb.ne.0)
           l=l+nwc
           call moreio(nin,0,0,c(l),nb,nwc)
         enddo
      endif
      elo=c(iraw+1)
      nlo=nle
      call hgetco(flo,nlo,lcd,c(iraw),lct,ltt,awr,idis)
      if (ltt.eq.1) call listio(nin,0,0,c(iraw),nb,nwc)
      if (ltt.eq.2) call tab1io(nin,0,0,c(iraw),nb,nwc)
      if (nb.ne.0) then
         l=iraw
         do while (nb.ne.0)
           l=l+nwc
           call moreio(nin,0,0,c(l),nb,nwc)
         enddo
      endif
      ehi=c(iraw+1)
      nhi=nle
      call hgetco(fhi,nhi,lcd,c(iraw),lct,ltt,awr,idis)
      ! check for total isotropy.
      if (ne.eq.2.and.nlo.eq.1.and.nhi.eq.1) iso=1
   endif
   if (iso.eq.1) then
      nle=1
      enext=etop
   else
      enext=elo
   endif
   nne=2
   ir=1
   inn=nint(c(6+2*ir))
   return

   !--normal entry
   !--is desired energy in current panel
  110 continue
   if (nle.eq.1) go to 150
   if (iso.eq.1) go to 150
   if (e.lt.ehi*(1-small)) go to 130
   if (nne.eq.ne.and.e.lt.ehi+ehi/100) go to 130

   !--no.  slide raw high energy data into low energy positions.
  120 continue
   do i=1,nhi
      flo(i)=fhi(i)
   enddo
   nlo=nhi
   elo=ehi

   !--read in new raw data at high energy.
   if (nne.eq.ne.and.ltt3.eq.3.and.lttn.eq.1) then
      call tab2io(nin,0,0,c(1),nb,nwc)
      ne=nint(c(6))
      nne=0
      ir=1
      ltt=2
      lttn=2
   else if (nne.eq.ne) then
      call error('hgtfle','desired energy above highest given.',' ')
   endif
   if (ltt.eq.1) call listio(nin,0,0,c(iraw),nb,nwc)
   if (ltt.eq.2) call tab1io(nin,0,0,c(iraw),nb,nwc)
   if (nb.ne.0) then
      l=iraw
      do while (nb.ne.0)
        l=l+nwc
        call moreio(nin,0,0,c(l),nb,nwc)
      enddo
   endif
   ehi=c(iraw+1)
   nhi=nle
   call hgetco(fhi,nhi,lcd,c(iraw),lct,ltt,awr,idis)
   nne=nne+1
   if (nne.gt.nint(c(5+2*ir))) ir=ir+1
   inn=nint(c(6+2*ir))
   if (ehi.le.e*(1+small).and.nne.lt.ne) go to 120

   !--yes.  interpolate for coefficients at desired energy.
  130 continue
   if (e.lt.elo*(1-small)) go to 140
   nlmax=nlo
   if (nlmax.lt.nhi) nlmax=nhi
   do i=1,nle
      if (i.le.nlmax) then
         innt=inn
         if (innt.eq.4.or.innt.eq.5) then
            if (flo(i)*fhi(i).le.zero) innt=innt-2
         endif
         call terp1(elo,flo(i),ehi,fhi(i),e,fle(i),innt)
      else
         fle(i)=0
      endif
   enddo
   nle=nlmax
   enext=ehi
   if (inn.eq.1) idis=1
   return

   !--return isotropic distribution below first point.
  140 continue
   fle(1)=1
   do il=2,nle
      fle(il)=0
   enddo
   nle=1
   enext=ehi
   idis=1
   return

   !--isotropic distribution or request.
  150 continue
   fle(1)=1
   if (nle.gt.1) then
      do il=2,nle
         fle(il)=0
      enddo
      nle=1
   endif
   enext=etop
   idis=0
   return
   end subroutine hgtfle

   subroutine hgetco(fl,nl,lcd,c,lct,ltt,awr,idis)
   !-------------------------------------------------------------------
   ! Retrieve or compute Legendre coefficients using raw data in
   ! ENDF format.  Conversion between reference systems and
   ! conversion of tabulated data to coefficients is performed by
   ! direct Gaussian quadrature of order 64. Significance of the
   ! results is reduced to reflect precision of calculation.
   !      nl   on entry, maximum number of coefficients desired
   !      nl   on return, number of coefficients required
   !      lcd  system desired (1 for lab, 2 for cm)
   !      lct  system of raw data
   !      ltt  format of raw data (1 for coefficients, 2 for tabulated)
   !      awr  atomic weight ratio
   !-------------------------------------------------------------------
   use util  ! provides error
   use mathm ! provides legndr
   use endf  ! provides terpa
   ! externals
   integer::nl,lcd,lct,ltt,idis
   real(kr)::fl(*),c(*),awr
   ! internals
   integer::l,np,il,ipc,irc,iq,ip,nlz,j
   real(kr)::x,y,xnext,xn
   real(kr)::p(65)
   real(kr),dimension(64),parameter::qp=(/&
     -9.99305042E-01_kr,-9.96340117E-01_kr,-9.91013371E-01_kr,&
     -9.83336254E-01_kr,-9.73326828E-01_kr,-9.61008800E-01_kr,&
     -9.46411375E-01_kr,-9.29569172E-01_kr,-9.10522137E-01_kr,&
     -8.89315446E-01_kr,-8.65999398E-01_kr,-8.40629296E-01_kr,&
     -8.13265315E-01_kr,-7.83972359E-01_kr,-7.52819907E-01_kr,&
     -7.19881850E-01_kr,-6.85236313E-01_kr,-6.48965471E-01_kr,&
     -6.11155355E-01_kr,-5.71895646E-01_kr,-5.31279464E-01_kr,&
     -4.89403146E-01_kr,-4.46366017E-01_kr,-4.02270158E-01_kr,&
     -3.57220158E-01_kr,-3.11322872E-01_kr,-2.64687162E-01_kr,&
     -2.17423644E-01_kr,-1.69644420E-01_kr,-1.21462819E-01_kr,&
     -7.29931218E-02_kr,-2.43502927E-02_kr, 2.43502927E-02_kr,&
      7.29931218E-02_kr, 1.21462819E-01_kr, 1.69644420E-01_kr,&
      2.17423644E-01_kr, 2.64687162E-01_kr, 3.11322872E-01_kr,&
      3.57220158E-01_kr, 4.02270158E-01_kr, 4.46366017E-01_kr,&
      4.89403146E-01_kr, 5.31279464E-01_kr, 5.71895646E-01_kr,&
      6.11155355E-01_kr, 6.48965471E-01_kr, 6.85236313E-01_kr,&
      7.19881850E-01_kr, 7.52819907E-01_kr, 7.83972359E-01_kr,&
      8.13265315E-01_kr, 8.40629296E-01_kr, 8.65999398E-01_kr,&
      8.89315446E-01_kr, 9.10522137E-01_kr, 9.29569172E-01_kr,&
      9.46411375E-01_kr, 9.61008800E-01_kr, 9.73326828E-01_kr,&
      9.83336254E-01_kr, 9.91013371E-01_kr, 9.96340117E-01_kr,&
      9.99305042E-01_kr/)
   real(kr),dimension(64),parameter::qw=(/&
      1.78328072E-03_kr, 4.14703326E-03_kr, 6.50445797E-03_kr,&
      8.84675983E-03_kr, 1.11681395E-02_kr, 1.34630479E-02_kr,&
      1.57260305E-02_kr, 1.79517158E-02_kr, 2.01348232E-02_kr,&
      2.22701738E-02_kr, 2.43527026E-02_kr, 2.63774697E-02_kr,&
      2.83396726E-02_kr, 3.02346571E-02_kr, 3.20579284E-02_kr,&
      3.38051618E-02_kr, 3.54722133E-02_kr, 3.70551285E-02_kr,&
      3.85501532E-02_kr, 3.99537411E-02_kr, 4.12625632E-02_kr,&
      4.24735151E-02_kr, 4.35837245E-02_kr, 4.45905582E-02_kr,&
      4.54916279E-02_kr, 4.62847966E-02_kr, 4.69681828E-02_kr,&
      4.75401657E-02_kr, 4.79993886E-02_kr, 4.83447622E-02_kr,&
      4.85754674E-02_kr, 4.86909570E-02_kr, 4.86909570E-02_kr,&
      4.85754674E-02_kr, 4.83447622E-02_kr, 4.79993886E-02_kr,&
      4.75401657E-02_kr, 4.69681828E-02_kr, 4.62847966E-02_kr,&
      4.54916279E-02_kr, 4.45905582E-02_kr, 4.35837245E-02_kr,&
      4.24735151E-02_kr, 4.12625632E-02_kr, 3.99537411E-02_kr,&
      3.85501532E-02_kr, 3.70551285E-02_kr, 3.54722133E-02_kr,&
      3.38051618E-02_kr, 3.20579284E-02_kr, 3.02346571E-02_kr,&
      2.83396726E-02_kr, 2.63774697E-02_kr, 2.43527026E-02_kr,&
      2.22701738E-02_kr, 2.01348232E-02_kr, 1.79517158E-02_kr,&
      1.57260305E-02_kr, 1.34630479E-02_kr, 1.11681395E-02_kr,&
      8.84675983E-03_kr, 6.50445797E-03_kr, 4.14703326E-03_kr,&
      1.78328072E-03_kr/)
   real(kr),parameter::toler=1.e-6_kr
   real(kr),parameter::rtoler=1.e6_kr
   real(kr),parameter::half=0.5e0_kr
   integer,parameter::nqp=64
   integer,parameter::nlmax=65
   real(kr),parameter::zero=0

   !--determine representation
   l=nl-1
   if (nl.gt.nlmax)&
     call error('hgetco','limited to 64 legendre coefficients.',' ')
   if (ltt.eq.2) go to 110
   if (lcd.eq.1.and.lct.ne.1) go to 110
   if (lcd.eq.2.and.lct.lt.2) go to 110

   !--retrieve coefficients from raw data
   np=nint(c(5))+1
   fl(1)=1
   do il=2,nl
      if (il.le.np) then
         fl(il)=c(il+5)
      else
         fl(il)=0
      endif
   enddo
   if (np.lt.nl) nl=np
   return

   !--integrate for fl in desired system
  110 continue
   fl(1)=1
   do il=2,nl
      fl(il)=0
   enddo
   ipc=2
   irc=1
   do iq=1,nqp
      x=qp(iq)

      !--compute scattering probability
      if (ltt.ne.2) then
         ! from coefficients
         np=nint(c(5))
         call legndr(x,p,np)
         y=half
         do ip=1,np
            y=(2*ip+1)*p(ip+1)*c(ip+6)/2+y
         enddo
      else
         ! from tabular data
         call terpa(y,x,xnext,idis,c,ipc,irc)
      endif

      !--multiply scattering probability by legendre polynomials
      !--in desired system and accumulate
      xn=x
      if (lcd.eq.2.and.lct.eq.1)&
        call error('hgetco','lab to cm conversion not coded.',' ')
      if (lcd.eq.1.and.lct.ge.2)&
        xn=(1+awr*x)/sqrt(1+awr*awr+2*awr*x)
      call legndr(xn,p,l)
      y=y*qw(iq)
      do il=2,nl
         fl(il)=fl(il)+y*p(il)
      enddo
   enddo

   !--reduce number of significant figures
   !--and count number of coefficients required
   nlz=1
   do il=2,nl
      j=nint(fl(il)*rtoler)
      fl(il)=j*toler
      if (abs(fl(il)).gt.zero) nlz=il
   enddo
   nl=nlz
   return
   end subroutine hgetco

   subroutine hconvr(nin,nout,nscr)
   !-------------------------------------------------------------------
   ! Convert photon transition probability arrays (LO=2), if any,
   ! to photon yields (LO=1).  Add MT456 if necessary.
   ! copy all other sections in this material.
   !-------------------------------------------------------------------
   use util ! provides error,mess
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr
   ! internals
   integer::l2flg,imax,lmax,i,istor,nb,nw,j,lg,mt0,mt0old,n,jm1
   integer::kk,k,ii,kp1,l,ja,jb,lm1,ip1,mtl,no455,lnu
   integer::mf,mt,l1,l2,n1,n2,nnu,idone,imax2
   integer::i10,maths,mtnow,mtmess,mttst,nn
   integer::m1,m2,nmf3,nmf12
   integer,parameter::nqmx=450
   real(kr)::g,ei,p,ysum,eja,ejb,yy,tsave
   character(60)::strng
   integer::ngam(300)
   integer::mtq(nqmx)
   real(kr)::eeq(nqmx)
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::e,eg,es,y,aa,r
   real(kr),dimension(:),allocatable::nu,tmpnu
   real(kr),parameter::zero=0

   !--initialize
   mt0=0
   mt0old=0
   l2flg=0
!   imax=49
   imax=50
   lmax=500
   do i=1,imax
      ngam(i)=0
   enddo
   istor=0
   allocate(scr(npage+50))
   nsh=0
   nsc=0

   !--check for consistancy between mf3 mt's and, if present, mf12
   !  mt's.  if any discrete inelastic levels are defined, we expect
   !  to a sequential list of mt values from 51 to a maximum of 90.
   !  Furthermore, for version 6 formatted endf files there may be
   !  sequential mt values from 600 to a maximum of 648, 650 to a
   !  maximum of 698, 700 to a maximum of 748, 750 to a maximum of
   !  798 and 800 to a maximum of 848 and 876 to 891 for outgoing
   !  protons, deuterons, tritons, 3He, alpha particles and the (n,2n)
   !  reaction, respectively.  The allowed mt interval for version 5
   !  formatted files differs, and is accounted for in the coding that
   !  follows.  The absence of an expected mf3 mt section is flagged
   !  as an error condition; the absence of an expected mf12 mt section
   !  is noted in a message since, while unusual, may be normal.
   mtq=0
   eeq=0
   mtnow=0
   mtmess=0
   i10=0
   nmf3=0
   nmf12=0
  10 continue
   call contio(nin,0,0,scr,nb,nw)
   if (i10.eq.0) then
      i10=1
      maths=math
   endif
   if (math.eq.0) go to 20
   if (mfh.lt.3 .or. (mfh.gt.3.and.mfh.lt.12)) then
       call tofend(nin,0,0,scr)
       go to 10
   endif
   if (mfh.gt.12) go to 20
   if (mfh.eq.3) then
      if (mth.eq.51.or.                                                &
          (iverf.ge.6.and.(mth.eq.601.or.mth.eq.651.or.mth.eq.701.or.  &
                           mth.eq.751.or.mth.eq.801.or.mth.eq.876)).or.&
          (iverf.le.5.and.(mth.eq.701.or.mth.eq.721.or.mth.eq.741.or.  &
                           mth.eq.761.or.mth.eq.781))) then
         mtnow=mth
      elseif ((mth.gt.51.and.mth.lt.91).or.                 &
              (iverf.ge.6.and.mth.gt.601.and.mth.lt.648).or.&
              (iverf.ge.6.and.mth.gt.651.and.mth.lt.698).or.&
              (iverf.ge.6.and.mth.gt.701.and.mth.lt.748).or.&
              (iverf.ge.6.and.mth.gt.751.and.mth.lt.798).or.&
              (iverf.ge.6.and.mth.gt.801.and.mth.lt.848).or.&
              (iverf.ge.6.and.mth.gt.876.and.mth.lt.890).or.&
              (iverf.le.5.and.mth.gt.701.and.mth.lt.718).or.&
              (iverf.le.5.and.mth.gt.721.and.mth.lt.738).or.&
              (iverf.le.5.and.mth.gt.741.and.mth.lt.758).or.&
              (iverf.le.5.and.mth.gt.761.and.mth.lt.778).or.&
              (iverf.le.5.and.mth.gt.781.and.mth.lt.798)) then
         if (mtnow+1.ne.mth) then
            mtmess=1
            write(strng,'(''mf3, mt'',i2,'' is missing'')')mth-1
            call mess('hconvr',strng,'')
         endif
         mtnow=mth
      endif
      call contio(nin,0,0,scr,nb,nw)
      nmf3=nmf3+1
      mtq(nmf3)=mth
      eeq(nmf3)=scr(2)
      call tosend(nin,0,0,scr)
      go to 10
   elseif (mfh.eq.12) then
      if (mth.eq.51.or.                                                &
          (iverf.ge.6.and.(mth.eq.601.or.mth.eq.651.or.mth.eq.701.or.  &
                           mth.eq.751.or.mth.eq.801.or.mth.eq.876)).or.&
          (iverf.le.5.and.(mth.eq.701.or.mth.eq.721.or.mth.eq.741.or.  &
                           mth.eq.761.or.mth.eq.781))) then
         mtnow=mth
      elseif ((mth.gt.51.and.mth.lt.91).or.                 &
              (iverf.ge.6.and.mth.gt.601.and.mth.lt.648).or.&
              (iverf.ge.6.and.mth.gt.651.and.mth.lt.698).or.&
              (iverf.ge.6.and.mth.gt.701.and.mth.lt.748).or.&
              (iverf.ge.6.and.mth.gt.751.and.mth.lt.798).or.&
              (iverf.ge.6.and.mth.gt.801.and.mth.lt.848).or.&
              (iverf.ge.6.and.mth.gt.876.and.mth.lt.890).or.&
              (iverf.le.5.and.mth.gt.701.and.mth.lt.718).or.&
              (iverf.le.5.and.mth.gt.721.and.mth.lt.738).or.&
              (iverf.le.5.and.mth.gt.741.and.mth.lt.758).or.&
              (iverf.le.5.and.mth.gt.761.and.mth.lt.778).or.&
              (iverf.le.5.and.mth.gt.781.and.mth.lt.798)) then
         if (mtnow+1.ne.mth) then
            write(strng,'(''mf12, mt'',i2,'' may be missing'')')mth-1
            call mess('hconvr',strng,&
                               'discrete photon data may be incomplete')
         endif
         mtnow=mth
      endif
      nmf12=nmf12+1
      call tosend(nin,0,0,scr)
      go to 10
   endif
   if (mtmess.ne.0) then
      call error('hconvr',&
                 'missing mf3 mt''s, probable endf error','')
   endif

   !--mt checks are done, rewind the endf tape and position
   !  at the first record for this material
  20 continue
   call repoz(nin)
   call findf(maths,1,451,nin)

   !--loop over all sections of this material
  110 continue
   call contio(nin,0,0,scr,nb,nw)
   if (mfh.eq.0.and.math.ne.0) go to 130
   if (mfh.eq.1) go to 550
   if (mfh.eq.12.and.mth.ne.460) go to 120
   if (math.eq.1149.and.mfh.eq.13.and.mth.ge.51.and.iverf.lt.5)&
     go to 140
   if (math.eq.1150.and.mfh.eq.13.and.mth.le.54.and.iverf.lt.5)&
     go to 140
   if (mfh.eq.14.and.l2flg.eq.1) go to 500
   call contio(0,nout,nscr,scr,nb,nw)
   if (math.eq.0) go to 600
   call tosend(nin,nout,nscr,scr)
   go to 110
  120 continue
   if (l1h.eq.2) go to 150
   call contio(0,nout,nscr,scr,nb,nw)
   call tosend(nin,nout,nscr,scr)
   go to 110
  130 continue
   call contio(0,nout,nscr,scr,nb,nw)
   go to 110
  140 continue
   call tosend(nin,0,0,scr)
   write(strng,'(''gamma prod patch made for mt '',i3)') mth
   call mess('hconvr',strng,' ')
   go to 110

   !--allocate storage
  150 continue
   if (istor.eq.0) then
      istor=1
      imax2=imax*imax
      allocate(e(imax))
      allocate(eg(lmax))
      allocate(es(lmax))
      allocate(y(lmax))
      allocate(aa(imax2))
      allocate(r(imax2))
      e(1)=0
      do i=1,imax2
         r(i)=0
         aa(i)=0
      enddo
      do i=1,imax
         r(imax*(i-1)+i)=1
      enddo
   endif

   !--convert transition probability array to yields
   za=c1h
   awr=c2h
   lg=l2h
   g=1
   l2flg=1
   call listio(nin,0,0,scr,nb,nw)

   !--set base value for mt0
   if (mth.ge.51.and.mth.le.91.and.mt0.ne.49) mt0=49
   if (iverf.ge.6) then
      if (mth.ge.601.and.mth.le.649.and.mt0.ne.599) mt0=599
      if (mth.ge.651.and.mth.le.699.and.mt0.ne.649) mt0=649
      if (mth.ge.701.and.mth.le.749.and.mt0.ne.699) mt0=699
      if (mth.ge.751.and.mth.le.799.and.mt0.ne.749) mt0=749
      if (mth.ge.801.and.mth.le.849.and.mt0.ne.799) mt0=799
      if (mth.ge.876.and.mth.le.891.and.mt0.ne.874) mt0=874
   elseif (iverf.le.5) then
      if (mth.ge.701.and.mth.le.719.and.mt0.ne.699) mt0=699
      if (mth.ge.721.and.mth.le.739.and.mt0.ne.719) mt0=719
      if (mth.ge.741.and.mth.le.759.and.mt0.ne.739) mt0=739
      if (mth.ge.761.and.mth.le.779.and.mt0.ne.759) mt0=759
      if (mth.ge.781.and.mth.le.799.and.mt0.ne.779) mt0=779
   endif

   !--load the e(_) array with mf3 -q.  Will overwrite with mf12
   !  mt data when possible (which should be always but if there
   !  is a missing mf12 mt, we're covered).
   if (mt0.ne.mt0old) then
      e=0
      mt0old=mt0
      m1=mt0+2
      if (mt0.eq.49) then
         m2=91
      elseif (iverf.ge.6.and.mt0.ne.874) then
         m2=m1+48
      elseif (iverf.ge.6.and.mt0.eq.874) then
         m2=m1+15
      elseif (iverf.le.5) then
         m2=m1+17
      endif
      nn=1
      mttst=-1
      do while (mttst.lt.m2.and.nn.le.nmf3)
         mttst=mtq(nn)
         if (mttst.ge.m1.and.mttst.le.m2) e(mttst-mt0)=-eeq(nn)
         nn=nn+1
      enddo
   endif
   j=mth-mt0
   e(j)=scr(1)
   n=n2h
   jm1=j-1
   do kk=1,jm1
      k=jm1-kk+1
      ii=0
      idone=0
      do while (ii.lt.n.and.idone.eq.0)
         ii=ii+1
         i=ii
         ei=scr(6+(lg+1)*i-lg)
         if (ei.eq.zero.and.e(k).eq.zero) then
            idone=1
         elseif (ei.ne.zero) then
            if (abs(ei-e(k))/ei.le.0.0001_kr) idone=1
         endif
      enddo
      if (idone.eq.0) then
         aa((k-1)*imax+j)=0
         r((k-1)*imax+j)=0
      else
         p=scr(7+(lg+1)*i-lg)
         g=1
         if (lg.eq.2) g=scr(6+(lg+1)*i)
         aa((k-1)*imax+j)=p*g
         r((k-1)*imax+j)=p
      endif
      if (k.ne.jm1) then
         kp1=k+1
         do ii=kp1,jm1
            r((k-1)*imax+j)=r((k-1)*imax+j)&
              +r((ii-1)*imax+j)*aa((k-1)*imax+ii)
         enddo
      endif
   enddo
   l=0
   ysum=0
   do i=2,j
      ja=j+2-i
      eja=e(ja)
      do ii=1,jm1
         jb=j-ii
         ejb=e(jb)
         yy=aa((jb-1)*imax+ja)*r((ja-1)*imax+j)
         if (yy.gt.zero) then
            l=l+1
            if (l.gt.lmax) then
               call error('hconvr','too many lo=2 gammas',' ')
            endif
            eg(l)=eja-ejb
            es(l)=eja
            y(l)=yy
            ysum=ysum+y(l)
         endif
      enddo
   enddo
   e(1)=0

   !--arrange gamma energies in descending order
   if (l.gt.1) then
      lm1=l-1
      do i=1,lm1
         ip1=i+1
         do ii=ip1,l
            if (eg(i).lt.eg(ii)) then
               tsave=eg(i)
               eg(i)=eg(ii)
               eg(ii)=tsave
               tsave=y(i)
               y(i)=y(ii)
               y(ii)=tsave
               tsave=es(i)
               es(i)=es(ii)
               es(ii)=tsave
            endif
         enddo
      enddo
   endif

   !--output computed yields in endf lo=1 format
   scr(1)=za
   scr(2)=awr
   scr(3)=1
   scr(4)=0
   scr(5)=l
   scr(6)=0
   mtl=49
   if (iverf.le.5.and.mth.ge.700) mtl=649
   if (iverf.ge.6.and.mth.ge.600) mtl=549
   ngam(mth-mtl)=l
   call contio(0,nout,nscr,scr,nb,nw)
   ! output tab1 sum record
   if (l.gt.1) then
      scr(1)=0
      scr(2)=0
      scr(3)=0
      scr(5)=1
      scr(6)=2
      scr(7)=2
      scr(8)=2
      scr(9)=ebot
      scr(10)=ysum
      scr(11)=etop
      scr(12)=ysum
      nw=12
      call tab1io(0,nout,nscr,scr,nb,nw)
   endif
   ! output rest of tab1 records
   do i=1,l
      scr(1)=eg(i)
      scr(2)=es(i)
      scr(3)=0
      scr(4)=2
      scr(5)=1
      scr(6)=2
      scr(7)=2
      scr(8)=2
      scr(9)=ebot
      scr(10)=y(i)
      scr(11)=etop
      scr(12)=y(i)
      nw=12
      call tab1io(0,nout,nscr,scr,nb,nw)
   enddo

   !--copy send record and loop over remaining sections.
   call contio(nin,nout,nscr,scr,nb,nw)
   go to 110

   !--angular distributions
  500 continue
   if (mth.ge.51.and.mth.le.90) then
      scr(3)=1
      scr(4)=0
      scr(6)=0
      scr(5)=ngam(mth-50)
      call contio(0,nout,nscr,scr,nb,nw)
      call tosend(nin,0,0,scr)
      call asend(nout,nscr)
   else
      call contio(0,nout,nscr,scr,nb,nw)
      call tosend(nin,nout,nscr,scr)
   endif
   go to 110

   !--file1.  add 456=452 if necessary.
  550 continue
   call contio(0,nout,nscr,scr,nb,nw)
   call tosend(nin,nout,nscr,scr)
   no455=0
   call contio(nin,nout,nscr,scr,nb,nw)
   if (mfh.eq.0) go to 110
   if (mth.ne.452) go to 595
!   nnu=9999
!   allocate(nu(nnu))
   if(allocated(tmpnu))deallocate(tmpnu)
   allocate(tmpnu(npage+46))
   l=1
   lnu=l2h
!   if (lnu.eq.1) call listio(nin,nout,nscr,nu(l),nb,nw)
!   if (lnu.eq.2) call tab1io(nin,nout,nscr,nu(l),nb,nw)
   if (lnu.eq.1) then
      call listio(nin,nout,nscr,tmpnu(l),nb,nw)
   elseif (lnu.eq.2) then
      call tab1io(nin,nout,nscr,tmpnu(l),nb,nw)
   endif
   nnu=nb+nw+1
   if (allocated(nu))deallocate(nu)
   allocate(nu(nnu))
   nu(1:nw)=tmpnu(1:nw)
   do while (nb.ne.0)
      l=l+nw
      if (l.gt.nnu) call error('hconvr',&
        'exceeded storage for nubar',' ')
      call moreio(nin,nout,nscr,nu(l),nb,nw)
   enddo
   call tosend(nin,nout,nscr,scr)
  565 continue
   call contio(nin,0,0,scr,nb,nw)
   if (mth.eq.456) go to 590
   if (mth.eq.0.and.no455.eq.0) go to 570
   if (mfh.eq.0) go to 590
   if (mth.eq.455) no455=1
   if (mth.gt.456.and.no455.eq.0) go to 570
   call contio(0,nout,nscr,scr,nb,nw)
   call tosend(nin,nout,nscr,scr)
   go to 565
  570 continue
   if (no455.eq.1) go to 590
   mf=mfh
   mt=mth
   l1=l1h
   l2=l2h
   n1=n1h
   n2=n2h
   mfh=1
   mth=456
   scr(3)=0
   scr(4)=lnu
   scr(5)=0
   scr(6)=0
   l=1
   call contio(0,nout,nscr,scr,nb,nw)
   if (lnu.eq.1) call listio(0,nout,nscr,nu(l),nb,nw)
   if (lnu.eq.2) call tab1io(0,nout,nscr,nu(l),nb,nw)
   do while (nb.ne.0)
      l=l+nw
      if (l.gt.nnu) call error('hconvr',&
        'exceeded storage for nubar',' ')
      call moreio(0,nout,nscr,nu(l),nb,nw)
   enddo
   call asend(nout,nscr)
   mfh=mf
   mth=mt
   scr(3)=l1
   scr(4)=l2
   scr(5)=n1
   scr(6)=n2
  590 continue
   call contio(0,nout,nscr,scr,nb,nw)
   deallocate(nu)
   if (mfh.eq.0) go to 110
  595 continue
   call tofend(nin,nout,nscr,scr)
   go to 110

   !--conver is finished.
  600 continue
   if (istor.gt.0) then
      deallocate(e)
      deallocate(eg)
      deallocate(es)
      deallocate(y)
      deallocate(aa)
      deallocate(r)
   endif
   deallocate(scr)
   return
   end subroutine hconvr

   subroutine hgam102(e,ebar,dame,disc102,c,irec,zp,ap,zt,at)
   !-------------------------------------------------------------------
   ! Process the relativistic discrete gamma or its recoil as
   ! given in mf6/mt102 for ENDF/B-VII neutron + H-1.
   !-------------------------------------------------------------------
   ! externals
   use mainio  ! provides nsysi,nsyso
   use snl     ! provides SNL
   integer::irec
   real(kr)::e,ebar,dame,disc102,zp,ap,zt,at
   real(kr)::c(*)
   ! internals
   real(kr)::er,eg2

   if (imode(3) .lt. 0) then
     write (nsyso,2301) e
2301 format (/,1x, 'HGAM102 entry ', 4g14.7, /)
   endif

   if (irec.eq.0) then
      ebar=disc102+e*awr/(1+awr)
      dame=0
   else
      !--include recoil energy plus photon "kick" energy
      er=e/(awr+1)
      eg2=disc102*disc102*rtm/2
      ebar=er+eg2
      dame=df(ebar,zp,ap,zt,at)
   endif
   return
   end subroutine hgam102

   subroutine gheat(iold,inew,nscr)
   !-------------------------------------------------------------------
   ! Subtract energy carried away by gammas.
   !-------------------------------------------------------------------
   use util ! provides error,mess,loada,finda,sigfig
   use mainio ! provides nsyso
   use snl     ! provides SNL
   use endf ! provides endf routines and variables
   ! externals
   integer::iold,inew,nscr
   ! internals
   integer::nb,nw,mfd,mgam,nk,ik,mfc,mtc,npkk,npkt
   integer::lo,mfx,mtx,mtd,jdis,ilist,nwd,lp,lqx,idis
   integer::i,ipx,irx,nmt,mty,idx,ii,indxx,it,isave,jerr,nd
   real(kr)::afact,aw1fac,z,e,dame
   real(kr):: h_prior
   real(kr)::cerr,enxt,el,elow,ehigh,test,thresh
   real(kr)::egkr,ebar,egam,edam,damn,enext,enx,egk
   real(kr)::h,hk,cfix,eava,ebarp,xp,yp,hp,egamp,edamp
   real(kr)::damep,hx,hxp,subtot,elo,ehi,x,y
   real(kr)::c(30)
   integer::imt(30)
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::b
   real(kr),dimension(:),allocatable::d
   real(kr),dimension(:),allocatable::bufo,bufn
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::smin=1.e-9_kr
   real(kr),parameter::qtest=99.e6_kr
   real(kr),parameter::zero=0
   integer::mf6flg
   character(len=70)::strng1,strng2

   !--allocate buffers for loada/finda
   allocate(bufo(nbuf))
   allocate(bufn(nbuf))

   if (imode(3) .lt. 0) then
     write (nsyso,2301)
2301 format (/,1x, 'GHEAT entry ',/)
   endif

   !--find first section of mf12 or mf13
   allocate(scr(npage+50))
   mfd=3
   afact=awr/(awr+1)
   aw1fac=1/(awr+1)
   z=nint(za/1000)
   e=0
   dame=df(e,z,awr,z,awr)
   mgam=2
  100 continue
   call contio(nscr,0,0,scr,nb,nw)
   if (mfh.eq.12) mgam=1
   if (mfh.ge.12) go to 102
   call tofend(nscr,0,0,scr)
   go to 100
  102 continue
   mfc=mfh
   mtc=mth
   call findf(matd,mfc,mtc,nscr)
   npkk=npk
   if (kchk.eq.1) npkk=3*npk-2
   if (mgam.gt.0) npkk=npkk+3
   npkt=npk-1

   !--sum all reactions found in mf12 and mf13.
   !--reconstruct redundant cross sections if needed.
   allocate(b(npage+50))
   nd=10000
   allocate(d(nd))
   if (mgam.eq.2) go to 210
  105 continue
   call contio(nscr,0,0,scr,nb,nw)
   if (mfh.gt.13) go to 250
   if (mfh.eq.0) go to 220
   if (mth.eq.0) go to 105
   if (mfh.eq.12.and.mth.eq.460) then
      call tosend(nscr,0,0,scr)
      go to 105
   endif
   ! skip over this mf/mt if photon data were already found in mf6
   if (i6p.gt.0) then
      mf6flg=0
      do i=1,i6p
         if (mt6yp(i).ne.mth) cycle
         mf6flg=1
         exit
      enddo
      if (mf6flg.ne.0) then
         write(strng1,'(''skipping mf'',i2,''/mt = '',i3)')mfh,mth
         write(strng2,'(''photons were already processed in mf6'')')
         call mess('gheat ',strng1,strng2)
         call tosend(nscr,0,0,scr)
         go to 105
      endif
   endif
   lo=l1h
   if (lo.eq.2) call error('gheat','lo=2 not coded.',' ')
   mfx=mfh
   mtd=mth
   mtx=mth
   if (mth.eq.3) mtd=1
   if (iprint.eq.1.and.mfh.eq.12.and.idame.gt.0.and.mtd.eq.102)&
     write(nsyso,'(/&
       &'' photon energy (from yields) mf'',i2,'', mt'',i3/&
       &14x,''e'',6x,''ebar/err'',10x,''egam'',10x,''edam'',10x,&
       &''xsec'',9x,''yield'',7x,''heating'',8x,''damage'')')&
       mfh,mth
   if (iprint.eq.1.and.mfh.eq.12.and.(idame.eq.0.or.mtd.ne.102))&
     write(nsyso,'(/&
       &'' photon energy (from yields) mf'',i2,'', mt'',i3/&
       &14x,''e'',10x,''ebar'',10x,''xsec'',9x,''yield'',&
       &8x,''energy'',7x,''heating'')') mfh,mth
   if (iprint.eq.1.and.mfh.eq.13) write(nsyso,'(/&
     &'' photon energy (from xsecs) mf'',i2,'', mt'',i3/&
     &14x,''e'',10x,''ebar'',10x,''xsec'',8x,''energy'',7x,&
     &''heating'')') mfh,mth
  110 continue
   nk=n1h
   ik=0
   cerr=0
   ! check for leading or trailing zeroes to locate ranges
   ! where this partial is not defined
   e=0
   call gety1(e,enxt,jdis,y,nscr,scr)
   e=enxt
   el=sigfig(e,7,0)
   elow=big
   ehigh=0
   test=big-big/100
   do while (e.lt.test)
      call gety1(e,enxt,jdis,y,nscr,scr)
      if (y.gt.zero.and.e.lt.elow) elow=el
      if (y.gt.zero) ehigh=0
      if (ehigh.eq.zero.and.y.eq.zero) ehigh=sigfig(e,7,0)
      el=sigfig(e,7,0)
      e=enxt
   enddo
   if (ehigh.eq.zero) ehigh=sigfig(e,7,0)
   if (nk.eq.1) then
      call findf(matd,mfx,mtx,nscr)
      call contio(nscr,0,0,scr,nb,nw)
   endif
  130 continue
   e=0
   call gety1(e,thresh,jdis,y,nscr,scr)
   ilist=1
   egkr=c1h
   nwd=nd
   if (egkr.eq.zero)&
     call gambar(e,ebar,egam,edam,nendf,matd,mtx,d,nwd)
   lp=l1h
   ik=ik+1
   if (iprint.eq.1.and.egkr.eq.zero) write(nsyso,&
     '(i6,''  continuum gammas'')') ik
   if (iprint.eq.1.and.egkr.gt.zero) write(nsyso,&
     '(i6,1p,e13.4,'' ev gamma'')') ik,egkr
   lqx=0
   if (mfd.eq.13) go to 150
   call findf(matd,mfd,mtd,nin)
   call contio(nin,0,0,b,nb,nw)
   call gety2(e,enxt,idis,x,nin,b)

   if (imode(3) .lt. 0 .and. e .lt. 1.e-4) then
     write (nsyso,4231) e, x, b
4231 format (1x, 'GHEAT x from gety2 ', 3g14.7)
   endif

   if (mth.eq.102) then
      q=c2h
      if (nqa.gt.0) then
         do i=1,nqa
            if (mta(i).eq.mth) then
               q=qa(i)
               if (q.ge.qtest) lqx=lqs(i)
            endif
         enddo
      endif
      call capdam(e,damn,q,za,awr,mth)
      call disgam(e,egam,edam,z,awr)
      ipx=2
      irx=1
   endif
   if (enxt.gt.thresh) thresh=enxt

   !--loop over all energies in grid
  150 continue
   if (ik.gt.1) go to 160
   nmt=0
   if (npk.lt.3) go to 160
   mty=mtx
   if (mty.eq.3) mty=1
   if (mty.eq.4) mty=51
   if (mty.eq.18.and.mt19.eq.1) mty=19
   call indx(npk,mtp,mt19,mty,imt,nmt)
  160 continue
   i=0
  190 continue
   i=i+1
   call finda(i,c,npkk,iold,bufo,nbuf)
   e=sigfig(c(1),9,0)
   if (e.lt.elow*(1-small).or.e.gt.ehigh*(1+small)) go to 180
   call gety1(e,enext,jdis,y,nscr,scr)
   x=zero
   if (mfd.eq.3) then
      call gety2(e,enext,idis,x,nin,b)
      if (x.eq.zero) go to 195
   endif
   if (lqx.ne.0) call terpa(q,e,enx,idx,qbar(lqx),ipx,irx)
   egk=egkr
   if (lp.eq.2) egk=egkr+e*afact
   ebar=egk
   dame=0
   if (ik.eq.1) c(npkk)=0
   if (egkr.eq.zero)&
     call gambar(e,ebar,egam,edam,nendf,matd,mtx,d,nwd)
   if (imode(3) .lt. 0 .and. e .lt. 1.e-4) then
     write (nsyso,7331) e, ebar, egam, edam
7331 format (1x, 'GHEAT gambar return ', 4g14.7)
   endif
   if (mfd.eq.13) go to 170
   if (mth.ne.102) go to 164
   ! photon recoil correction
   if (egkr.ne.zero) call disgam(egkr,egam,edam,z,awr)
   if (imode(3) .lt. 0 .and. e .lt. 1.e-4) then
     write (nsyso,7332) e, egkr, egam, edam
7332 format (1x, 'GHEAT disgam return ', 4g14.7)
   elseif (imode(3) .lt. 0 .and. e .gt. 9.51e6 .and. e .lt. 10.1E6) then
     write (nsyso,7332) e, egkr, egam, edam
   endif
   h=egam*x*y
   hk=h
   c(npkk-1)=c(npkk-1)+x*y*ebar

   if (imode(3) .lt. 0 .and. e .le. 1.e-4) then
     write (nsyso,8332) e, c(npkk-1), x, y, ebar, egam, edam, h, y*egam, y*edam, npkk
8332 format (1x, 'GHEAT c(npkk-1)a: ', 10g14.7, 2x, i6)
   elseif (imode(3) .lt. 0 .and. e .ge. 9.52e6 .and. e .le. 10.1E6) then
     write (nsyso,8332) e, c(npkk-1), x, y, ebar, egam, edam, h, y*egam, y*edam, npkk
   endif

   c(npkk)=c(npkk)+y*ebar
   dame=edam*x*y

   if (imode(3) .lt. 0 .and. e .lt. 1.e-4) then
     write (nsyso,6331) dame, edam, x, y, idame, ik, nk
6331 format (1x, 'GHEAT dame calc ', 4g14.7, 3i6)
   elseif (imode(3) .lt. 0 .and. e .gt. 9.53E6 .and. e .lt. 10.1E6) then
     write (nsyso,6331) dame, edam, x, y, idame, ik, nk
   endif

   if (ik.eq.nk) then
      if (idame.gt.0) then
         call capdam(e,damn,q,za,awr,mth)
         dame=dame-damn*x
      endif
      hk=hk-x*rtm*(q+e*awr*aw1fac)**2/2
      cfix=(e+q-e*aw1fac)*x
      h_prior = h
      h=h-cfix
      c(npkk-2)=cfix-x*c(npkk)
      eava=afact*e+q
      cerr=100*(c(npkk)-eava)/eava
      if (imode(3) .lt. 0 .and. e .lt. 1.e-4) then
        write (nsyso,6332) e, cfix, q, aw1fac, x, h, h_prior, c(npkk-2), x*c(npkk), hk, eava, npkk
6332    format (1x, 'GHEAT hk calc ', 11g14.7, 2i6)
      elseif (imode(3) .lt. 0 .and. e .gt. 9.53E6 .and. e .lt. 10.1E6) then
        write (nsyso,6332) e, cfix, q, aw1fac, x, h, h_prior, c(npkk-2), x*c(npkk), hk, eava, npkk
      endif
   endif
   ebarp=sigfig(ebar,9,0)
   xp=sigfig(x,9,0)
   yp=sigfig(y,9,0)
   hp=sigfig(h,9,0)
   egamp=sigfig(egam,9,0)
   edamp=sigfig(edam,9,0)
   damep=sigfig(dame,9,0)
   cerr=sigfig(cerr,4,0)
   if (iprint.ne.1.or.abs(e-elist(ilist)).gt.small*e&
     .or.y.le.smin) go to 166
   if (idame.eq.0) then
      write(nsyso,'(1x,1p,8e14.4)') e,ebarp,xp,yp,hp
   else
      write(nsyso,'(1x,1p,8e14.4)')&
        e,ebarp,egamp,edamp,xp,yp,hp,damep

   if (imode(3) .lt. 0 .and. e .lt. 1.e-4) then
     write (nsyso,2331) ik, egkr, e, yp, hp, damep
2331 format (1x, 'GHEAT component echo ', i6, 5g14.7)
   endif

   endif
   if (ik.eq.nk) write(nsyso,&
     '(1x,1p,e14.4,0p,f10.1,'' pc'')') e,cerr
   go to 166
   ! other mf12 reactions
  164 continue
   if (mtd.ne.2.and.mtd.ne.460) h=-y*x*ebar
   if (mtd.eq.2) h=y*x*ebar
   c(npkk-1)=c(npkk-1)-h
   c(npkk)=c(npkk)+h

   if (imode(3) .lt. 0 .and. e .lt. 1.e-4) then
     write (nsyso,8331) e, c(npkk), h, egkr
8331 format (1x, 'GHEAT h increment ', 4g14.7)
   elseif (imode(3) .lt. 0 .and. e .ge. 9.52e6 .and. e .le. 10.1E6) then
     write (nsyso,8331) e, c(npkk), h, egkr
   endif

  166 continue
   hx=h
   c(2)=c(2)+h
   if (nmt.gt.0) then
      do ii=1,nmt
         indxx=imt(ii)
         if (mtp(indxx).ge.444) c(indxx)=c(indxx)+dame

         if ( imode(3) .lt. 0 .and. mtp(indxx) .eq. 401 .and. c(1) .lt. 1.03e-5) then 
              write (nsyso,6301) c(1), c(2), c(indxx), h, indxx, dame, y, x, ebar
6301          format (1x, 'GHEAT partial sum-a-401: ', 4g14.7, i6, 4g14.6)
         elseif (imode(3) .lt. 0 .and. mtp(indxx) .eq. 401 .and. c(1) .ge. 9.52e6 .and. c(1) .le. 10.1E6) then
              write (nsyso,6301) c(1), c(2), c(indxx), h, indxx, dame, y, x, ebar
         endif
         if ( imode(3) .lt. 0 .and. mtp(indxx) .eq. 442 .and. c(1) .lt. 1.03e-5) then 
              write (nsyso,6302) c(1), c(2), c(indxx), h, indxx, dame, y, x, ebar
6302          format (1x, 'GHEAT partial sum-a-442: ', 4g14.7, i6, 4g14.6)
         elseif (imode(3) .lt. 0 .and. mtp(indxx) .eq. 442 .and. c(1) .ge. 9.52e6 .and. c(1) .le. 10.1E6) then
              write (nsyso,6302) c(1), c(2), c(indxx), h, indxx, dame, y, x, ebar
         endif

         if (mtp(indxx).lt.442) c(indxx)=c(indxx)+h
         if (mtp(indxx).eq.442) c(indxx)=c(indxx)-h
         if (mtp(indxx).eq.443.and.mth.eq.102) c(indxx)=c(indxx)+hk

   if (imode(3) .lt. 0 .and. e .le. 1.e-5) then
     write (nsyso,8302) e, c(indxx), h, hk, indxx, mth, mtp(indxx)
8302 format (1x, 'c-look-2 ', 4g14.7, 2x, 3i6)
   endif

      enddo
   endif

         if ( imode(3) .lt. 0 .and. mtp(indxx) .eq. 401 .and. c(1) .lt. 1.03e-5) then 
              write (nsyso,6391) c(1), c(indxx), h
6391          format (1x, 'GHEAT MF401 calc: ', 4g14.7)
         elseif (imode(3) .lt. 0 .and. mtp(indxx) .eq. 401 .and. c(1) .ge. 9.52e6 .and. c(1) .le. 10.1E6) then
              write (nsyso,6391) c(1), c(indxx), h
         endif
         if ( imode(3) .lt. 0 .and. mtp(indxx) .eq. 442 .and. c(1) .lt. 1.03e-5) then 
              write (nsyso,6392) c(1), c(indxx), h
6392          format (1x, 'GHEAT MF442 calc: ', 4g14.7)
        elseif (imode(3) .lt. 0 .and. mtp(indxx) .eq. 442 .and. c(1) .ge. 9.52e6 .and. c(1) .le. 10.1E6) then
              write (nsyso,6392) c(1), c(indxx), h
         endif
         if ( imode(3) .lt. 0 .and. mtp(indxx) .eq. 447 .and. c(1) .lt. 1.03e-5) then 
              write (nsyso,6393) c(1), c(indxx), h
6393          format (1x, 'GHEAT MF447 calc: ', 4g14.7)
        elseif (imode(3) .lt. 0 .and. mtp(indxx) .eq. 447 .and. c(1) .ge. 9.52e6 .and. c(1) .le. 10.1E6) then
              write (nsyso,6393) c(1), c(indxx), h
         endif

   if (iprint.eq.0.or.e.ne.elist(ilist)) go to 180
   if (mth.eq.102) go to 180
   if (ik.lt.nk.and.y.lt.smin) go to 180
   if (nk.eq.1.and.y.lt.smin) go to 180
   test=-1
   test=test/10000
   if (nk.gt.1.and.ik.eq.nk.and.c(npkk).gt.test) go to 180
   ebarp=sigfig(ebar,9,0)
   xp=sigfig(x,9,0)
   yp=sigfig(y,9,0)
   hxp=sigfig(hx,9,0)
   hp=sigfig(h,9,0)
   write(nsyso,'(1x,1p,8e14.4)') e,ebarp,xp,yp,-hxp,hp
   if (ik.eq.nk.and.nk.gt.1) then
      subtot=-c(npkk)
      subtot=sigfig(subtot,9,0)
      write(nsyso,'(48x,''subtot   '',1p,e14.4)') subtot
   endif
   go to 180
   ! mf13 reactions
  170 continue
   h=-y*ebar
   c(npkk-1)=c(npkk-1)-h
   c(npkk)=c(npkk)+h
   hx=h
   c(2)=c(2)+h
   if (nmt.gt.0) then
      do ii=1,nmt
         indxx=imt(ii)

         if ( imode(3) .lt. 0 .and. mtp(indxx) .eq. 401 .and. c(1) .lt. 1.03e-5) then 
              write (nsyso,6102) c(1), c(2), c(indxx), h, indxx
6102          format (1x, 'GHEAT partial sum-b: ', 4g14.7, i6)
        elseif (imode(3) .lt. 0 .and. mtp(indxx) .eq. 401 .and. c(1) .ge. 9.52e6 .and. c(1) .le. 10.1E6) then
              write (nsyso,6102) c(1), c(2), c(indxx), h, indxx
         endif


         if (mtp(indxx).lt.442) c(indxx)=c(indxx)+h
         if (mtp(indxx).eq.442) c(indxx)=c(indxx)-h

         if ( imode(3) .lt. 0 .and. mtp(indxx) .eq. 401 .and. c(1) .lt. 1.03e-5) then 
              write (nsyso,6502) c(1), c(2), c(indxx), h, indxx
6502          format (1x, 'GHEAT partial sum-f: ', 4g14.7, i6)
        elseif (imode(3) .lt. 0 .and. mtp(indxx) .eq. 401 .and. c(1) .ge. 9.52e6 .and. c(1) .le. 10.1E6) then
              write (nsyso,6502) c(1), c(2), c(indxx), h, indxx
         endif

      enddo
   endif
   if (iprint.eq.0.or.abs(e-elist(ilist)).gt.small*e) go to 180
   if (ik.lt.nk.and.y.lt.smin) go to 180
   if (nk.eq.1.and.y.lt.smin) go to 180
   test=-1
   test=test/10000
   if (nk.gt.1.and.ik.eq.nk.and.c(npkk).gt.test) go to 180
   ebarp=sigfig(ebar,9,0)
   yp=sigfig(y,9,0)
   hxp=sigfig(hx,9,0)
   hp=sigfig(h,9,0)
   write(nsyso,'(1x,1p,8e14.4)') e,ebarp,yp,-hxp,hp
   if (ik.eq.nk.and.nk.gt.1) then
      subtot=-c(npkk)
      subtot=sigfig(subtot,9,0)
      write(nsyso,'(48x,''subtot   '',1p,e14.4)') subtot
   endif
  180 continue
  195 continue
   it=i
   if (i.eq.ne) it=-it
   call loada(it,c,npkk,inew,bufn,nbuf)
   if (abs(e-elist(ilist)).lt.small*e) ilist=ilist+1
   if (i.lt.ne) go to 190

   !--loop over all subsections and reactions.
   isave=iold
   iold=inew
   inew=isave
   if (ik.lt.nk) go to 130
   if (mfd.eq.13) go to 105
   if (mtx.eq.3.and.mtd.eq.1) go to 200
   go to 105
  200 continue
   mtd=2
   call findf(matd,mfx,mtx,nscr)
   call contio(nscr,0,0,scr,nb,nw)
   go to 110
  210 continue
   call mess('gheat','no file 12 for this material.',' ')
   call findf(matd,13,0,nscr)

   !--finished with present file.
  220 continue
   if (mfd.eq.13) go to 250
   mfd=13
   mtd=0
   go to 105
  250 continue
   if (kchk.ne.1) go to 300
   if (mt303.eq.0) go to 300

   !--print kinematic test for total photon production
   ilist=1
   write(nsyso,'(/'' photon energy production check''/&
     &14x,''e'',6x,''ev-barns'',11x,''min'',11x,''max'')')
   do i=1,ne
      call finda(i,c,npkk,iold,bufo,nbuf)
      e=sigfig(c(1),9,0)
      if (abs(e-elist(ilist)).le.small*e) then
         ilist=ilist+1
         elo=c(npkk-1)+c(mt303)-c(mt303+2*npkt)+c(npkk-2)
         ehi=c(npkk-1)+c(mt303)-c(mt303+npkt)+c(npkk-2)
         jerr=0
         if (c(npkk-1).gt.ehi+ehi/10) jerr=+1
         if (c(npkk-1).lt.elo-elo/10) jerr=-1
         if (jerr.eq.0) write(nsyso,'(1x,1p,8e14.4)')&
           e,c(npkk-1),elo,ehi
         if (jerr.gt.0) write(nsyso,'(1x,1p,4e14.4,''   ++++'')')&
           e,c(npkk-1),elo,ehi
         if (jerr.lt.0) write(nsyso,'(1x,1p,4e14.4,''   ----'')')&
           e,c(npkk-1),elo,ehi
      endif
   enddo
  300 continue
   deallocate(d)
   deallocate(b)
   deallocate(scr)
   deallocate(bufn)
   deallocate(bufo)
   return
   end subroutine gheat

   subroutine gambar(e,ebar,esqb,esqd,nin,matd,mtd,a,nwamax)
   !-------------------------------------------------------------------
   ! Calculate ebar for continuous spectra and the photon
   ! recoil correction for capture.
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   use mainio  ! provides nsysi,nsyso
   use snl     ! provides SNL
   ! externals
   integer::nin,matd,mtd,nwamax
   real(kr)::e,ebar,esqb,esqd,a(*)
   ! internals
   integer::nb,nw,mfd,l,iz,lf,iend,j,il
   integer::jstart,jend,istart,ilo,ihi,nnt,ne,nne,inn,nbt
   real(kr)::z,awr,elo,ehi,flo,fhi,glo,ghi,hlo,hhi,sen
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::up=1.000001e0_kr
   real(kr),parameter::zero=0
   save lf,z,awr
   save jstart,jend,istart,ilo,ihi,nnt,ne,nne,inn,nbt
   save elo,ehi,flo,fhi,glo,ghi,hlo,hhi

   if (imode(3) .lt. 0) then
     write (nsyso,2301) e, matd, mtd
2301 format (/,1x, 'GAMBAR entry ', g14.7, 2i6,/)
   endif

   !--initialize if e=0.
   if (e.gt.zero) go to 110
   mfd=15
   call findf(matd,mfd,mtd,nin)
   l=1
   call contio(nin,0,0,a(l),nb,nw)
   iz=nint(c1h/1000)
   z=iz
   awr=c2h
   call tab1io(nin,0,0,a(l),nb,nw)
   l=l+nw
   do while (nb.ne.0)
      if (l.gt.nwamax) call error('gambar',&
        'storage exceeded in a.',' ')
      call moreio(nin,0,0,a(l),nb,nw)
      l=l+nw
   enddo
   lf=1

   !--read in first incident energy point.
   l=1
   call tab2io(nin,0,0,a(l),nb,nw)
   ne=nint(a(l+5))
   nnt=l+6
   nbt=nint(a(nnt))
   inn=nint(a(nnt+1))
   inn=2
   istart=l+nw
   l=istart
   call tab1io(nin,0,0,a(l),nb,nw)
   l=l+nw
   do while (nb.ne.0)
      if (l.gt.nwamax) call error('gambar',&
        'storage exceeded in a.',' ')
      call moreio(nin,0,0,a(l),nb,nw)
      l=l+nw
   enddo
   iend=l-1
   jstart=iend+1
   jend=0
   ehi=a(istart+1)
   elo=ehi
   nne=1
   ilo=0
   ihi=0
   fhi=0
   ghi=0
   hhi=0
   return

   !--normal entry.  is desired energy in current panel.
  110 continue
   if (e.lt.ehi*(1-small)) go to 140
   if (nne.eq.1) go to 130
   if (nne.eq.ne) go to 160

   !--no.  slide high energy data into low energy positions.
  120 continue
   j=-1
   do il=jstart,jend
      j=j+1
      a(j+istart)=a(il)
   enddo
   iend=istart+j
   jstart=iend+1
   elo=ehi
   flo=fhi
   if (mtd.eq.102) glo=ghi
   if (mtd.eq.102) hlo=hhi
   ilo=ihi

   !--read new high energy data record.
  130 continue
   if (nne.eq.ne)&
     call error('gambar','requested energy gt highest given.',' ')
   l=jstart
   call tab1io(nin,0,0,a(l),nb,nw)
   l=l+nw
   do while (nb.ne.0)
      if (l.gt.nwamax) call error('gambar',&
        'storage exceeded in a.',' ')
      call moreio(nin,0,0,a(l),nb,nw)
      l=l+nw
   enddo
   jend=l-1
   nne=nne+1
   ihi=0
   ehi=a(jstart+1)
   if (ehi.lt.e*(1-small).and.nne.lt.ne) go to 120

   !--integrate at high and low energies.
   if (ilo.le.0) then
      call tabbar(flo,a(istart),lf)
      if (mtd.eq.102) call tabsqr(glo,hlo,a(istart),lf,z,awr)
   endif
   call tabbar(fhi,a(jstart),lf)
   if (mtd.eq.102) call tabsqr(ghi,hhi,a(jstart),lf,z,awr)
   ihi=1

   !--yes.  interpolate for mean energy.
  140 continue
   if (nne.eq.1) go to 150
   if (nne.gt.nbt) then
      nnt=nnt+2
      nbt=nint(a(nnt))
      inn=nint(a(nnt+1))
      inn=2
   endif
   if (elo.eq.ehi) then
       elo=sigfig(elo,7,-1)
       ehi=sigfig(ehi,7,+1)
   endif
   call terp1(elo,flo,ehi,fhi,e,sen,inn)
   if (mtd.eq.102) call terp1(elo,glo,ehi,ghi,e,esqb,inn)
   if (mtd.eq.102) call terp1(elo,hlo,ehi,hhi,e,esqd,inn)
   ebar=sen
   return

   !--return zero below threshold.
  150 continue
   ebar=0
   esqb=0
   esqd=0
   return

   !--return zero or fhi above table.
  160 continue
   if (inn.ne.1.and.e.ge.up*ehi) then
      ebar=0
      esqb=0
      esqd=0
   else
      ebar=fhi
      esqb=ghi
      esqd=hhi
   endif
   return
   end subroutine gambar

   subroutine tabsqr(g,h,a,law,z,awr)
   !-------------------------------------------------------------------
   ! Compute average of photon recoil energy from capture and
   ! corresponding damage energy for a tabulated section of File 15.
   !-------------------------------------------------------------------
   use endf ! provides terp1
   use mainio  ! provides nsysi,nsyso
   use snl     ! provides SNL
   ! externals
   integer::law
   real(kr)::g,h,a(*),z,awr
   ! internals
   integer::nr,np,ibase,ir,nbt,inn,i,j
   real(kr)::ein,rein,xl,yl,xh,yh,dx,x,y,xr,s
   integer,parameter::nq=4
   real(kr),dimension(nq),parameter::qp=(/&
     -.86114e0_kr,-.33998e0_kr,.33998e0_kr,.86114e0_kr/)
   real(kr),dimension(nq),parameter::qw=(/&
     .34785e0_kr,.65215e0_kr,.65215e0_kr,.34785e0_kr/)

   if (imode(3) .lt. 0) then
     write (nsyso,2301) awr
2301 format (/,1x, 'TABSQR entry ', g14.7)
   endif

   !--initialize
   ein=2*tm
   rein=1/ein
   g=0
   h=0
   s=0
   nr=nint(a(5))
   np=nint(a(6))
   ibase=6+2*nr
   ir=1
   nbt=nint(a(7))
   inn=nint(a(8))
   xh=a(ibase+1)
   yh=a(ibase+2)

   !--accumulate contributions to integrals
   do i=2,np
      xl=xh
      xh=a(ibase+2*i-1)
      yl=yh
      yh=a(ibase+2*i)
      if (i.gt.nbt) then
         ir=ir+1
         nbt=nint(a(ibase+2*ir-1))
         inn=nint(a(ibase+2*ir))
      endif
      if (xl.ne.xh) then
         dx=xh-xl
         do j=1,nq
            x=xl+(1+qp(j))*dx/2
            call terp1(xl,yl,xh,yh,x,y,inn)
            xr=x*x*rein
            g=g+qw(j)*xr*dx*y
            h=h+qw(j)*df(xr,z,awr+1,z,awr)*dx*y
            s=s+qw(j)*dx*y

            if (imode(3) .lt. 0) then
              write (nsyso,2392) xl, xh, yl, yh, g, h, s, xr
2392          format (1x, 'TABSQR sum: ', 8g14.7)
           endif

         enddo
      endif
   enddo
   g=g/s
   h=h/s

   if (imode(3) .lt. 0) then
     write (nsyso,2391) g, h, s
2391 format (/,1x, 'TABSQR exit ', 7g14.7, /)
   endif

   return
   end subroutine tabsqr

   subroutine disgam(e,egam,edam,z,awr)
   !-------------------------------------------------------------------
   ! Compute recoil and damage energy of discrete capture photons.
   !-------------------------------------------------------------------
   ! externals
   use mainio  ! provides nsysi,nsyso
   use snl     ! provides SNL
   real(kr)::e,egam,edam,z,awr
   ! internals
   real(kr),parameter::zero=0

   if (imode(3) .lt. 0 .and. e .lt. 1.e-4) then
     write (nsyso,2301) e, z, awr, egam, edam
2301 format (1x, 'DISGAM entry ', 5g14.7)
   endif

   if (e.eq.zero) then
      return
   else
      egam=e*e*rtm/2
      edam=df(egam,z,awr+1,z,awr)
   endif

   if (imode(3) .lt. 0 .and. e .lt. 1.e-4) then
     write (nsyso,2302) e, egam, edam, rtm
2302 format (1x, 'DISGAM exit ', 4g14.7)
   endif

   return
   end subroutine disgam

   subroutine hout(iold)
   !-------------------------------------------------------------------
   ! Write the output tape with corrected dictionary and
   ! desired MT-s added.
   !-------------------------------------------------------------------
   use util ! provides openz,repoz,closz,error,sigfig
   use mainio ! provides nsyso
   use endf ! provides endf routines and variables
   use snl     ! provides SNL
   ! externals
   integer::iold
   ! internals
   integer::nscr,npkk,npkt,npkd,npktd,i,ilist,inpk,nb,nw
   integer::nxn,k,mfi,mti,mtn,ia,nn
   integer::idict,idict1,nx,inow,nwd,n,ilo,ihi,j,ibase
   integer::iowr,mfnin,mtnin,mfnscr,mtnscr,inow6
   real(kr)::x,y,xlast,ylo,yhi,thin,rat,elo,test,e
   real(kr)::c(30)
   integer::ncds(30)
   character(4)::klo(9),khi(9)
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::b
   real(kr),dimension(:),allocatable::buf
   character(60)::strng
   character(4)::iblank='    '
   character(4)::idown='low '
   character(4)::iup='high'
   character(1),parameter::qu=''''
   real(kr),parameter::small=1.e-10_kr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::emax=20.e6_kr
   real(kr),parameter::ecut=.01e6_kr
   real(kr),parameter::zero=0
   real(kr),parameter::ten=10
   integer,parameter::nwscr=10000

   ! allocate scratch storage for endf records and buffer
   allocate(scr(nwscr))
   allocate(buf(nbuf))

   if (imode(3) .lt. 0 .and. c(1) .lt. 1.e-4) then
     write (nsyso,1302) iold
1302 format (1x, 'hout entry ', 3i6)
   endif

   !--construct tab1 records for the kerma factors
   nscr=15
   if (nout.lt.0) nscr=-nscr
   call openz(nscr,1)
   call repoz(nscr)
   if (iprint.eq.1) write(nsyso,&
     '(/'' final kerma factors''/14x,''e'',8(11x,i3))')&
     (mtp(i),i=2,npk)
   npkk=npk
   if (kchk.eq.1) npkk=3*npk-2
   if (mgam.gt.0) npkk=npkk+3
   if (mgam.eq.0.and.i6.gt.0) npkk=npkk+1
   npkt=npk-1
   npkd=npk
   npktd=npkt
   do i=2,npk
      if (mtp(i).ge.442) then
         npkd=npkd-1
         npktd=npktd-1
      endif
   enddo
   ilist=1
   nsc=0
   math=matd
   mfh=3
   inpk=1
  105 continue
   inpk=inpk+1
   mth=mtp(inpk)
   scr(1)=za
   scr(2)=awr
   scr(3)=0
   scr(4)=0
   scr(5)=0
   scr(6)=0
   call contio(0,0,nscr,scr,nb,nw)
   nn=0
   scr(1)=0
   scr(2)=0
   scr(5)=1
   scr(6)=0
   scr(7)=0
   scr(8)=2
   ibase=8
   i=0
   nw=0
  110 continue
   nn=nn+1
   if (nn.gt.ne) go to 125
   call finda(nn,c,npkk,iold,buf,nbuf)
   do j=1,npkk
      if (c(j).gt.big) then
         c(j)=sigfig(c(j),6,0)
      else
         c(j)=sigfig(c(j),7,0)
      endif
   enddo
   e=sigfig(c(1),9,0)
   if (iprint.eq.1.and.abs(e-elist(ilist)).le.small*e) then
      ihi=0
      ilo=0
      if (kchk.eq.1) then
         ihi=0
         ilo=0
         do j=2,npkd
            read(iblank,'(a4)') klo(j)
            test=c(j+npkt)-c(j+npkt)/10
            if (c(j).lt.test) then
               read(idown,'(a4)') klo(j)
               ilo=1
            endif
            read(iblank,'(a4)') khi(j)
            test=c(j+2*npkt)+c(j+2*npkt)/10
            if (c(j).gt.test) then
               read(iup,'(a4)') khi(j)
               ihi=1
            endif
         enddo
      endif
      if (ilo.eq.1) write(nsyso,'('' '')')
      if (kchk.eq.1) write(nsyso,'(15x,8(6x,a4,4x))')&
        (klo(j),j=2,npkd)
      if (kchk.eq.1) write(nsyso,'(12x,''min'',1p,8e14.4)')&
        (c(npk+j),j=1,npktd)
      write(nsyso,'(1x,1p,9e14.4)') (c(j),j=1,npk)

   if (imode(3) .lt. 0 .and. c(1) .lt. 1.e-4) then
     write (nsyso,9302) c(1), c(2), c(3), c(4), c(5), c(6), c(npk+2), c(npk+npkt+2), npk, npktd
9302 format (1x, 'final kerma 401 look ', 8g14.7, 2x, 2i6)
   endif

      if (kchk.eq.1) write(nsyso,'(12x,''max'',1p,8e14.4)')&
        (c(npk+npkt+j),j=1,npktd)
      if (ihi.eq.1) write(nsyso,'(15x,8(6x,a4,4x))')&
        (khi(j),j=2,npkd)
   endif
   if (e.ge.elist(ilist)) ilist=ilist+1
   if (i.gt.0) go to 118
   if (c(inpk).eq.zero.and.i.eq.0.and.ibase.eq.8) go to 110
   if (nw.gt.0) go to 118
   n=ne-nn+1
   scr(6)=n
   scr(7)=n
  118 continue
   i=i+2
   scr(ibase+i-1)=c(1)
   scr(ibase+i)=c(inpk)
   if (i.lt.npage.and.nn.ne.ne) go to 110
   if (ibase.eq.0) go to 120
   call tab1io(0,0,nscr,scr,nb,nw)
   if (nn.eq.ne) go to 130
   ibase=0
   i=0
   go to 110
  120 continue
   call moreio(0,0,nscr,scr,nb,nw)
   if (nn.eq.ne) go to 130
   i=0
   go to 110
   ! all zeroes for this mt
  125 continue
   scr(6)=2
   scr(7)=2
   scr(9)=efirst
   scr(10)=0
   scr(11)=elast
   scr(12)=0
   nw=12
   call tab1io(0,0,nscr,scr,nb,nw)
   n=2
  130 continue
   call asend(0,nscr)
   ncds(inpk)=3+(n+2)/3
   if (inpk.lt.npk) go to 105
   call afend(0,nscr)
   call repoz(nscr)

   !--update dictionary in file 1
   call findf(matd,1,451,nin)
   call contio(nin,0,0,scr,nb,nw)
   nx=n2h
   inow=7
   if (iverf.ge.5) then
      call contio(nin,0,0,scr(inow),nb,nw)
      inow=inow+6
      if (iverf.eq.6) then
         call contio(nin,0,0,scr(inow),nb,nw)
         inow=inow+6
      endif
   endif
   call hdatio(nin,0,0,scr(inow),nb,nw)
   nwd=nint(scr(inow+4))
   if (iverf.ge.5) nx=nint(scr(inow+5))
   inow=inow+nw
   do while (nb.ne.0)
      call moreio(nin,0,0,scr(inow),nb,nw)
      inow=inow+nw
   enddo
   nw=nx
   idict=inow
   idict1=idict-1
   call dictio(nin,0,0,scr(idict),nb,nw)
   nw=inow+nw
   nw=6*(nx+20)
   if (nw.lt.npage+50) nw=npage+50
   allocate(b(nw))
   inpk=2
   mtn=mtp(inpk)
   j=0

   !--loop over existing dictionary records
   do i=1,nx
      iowr=0
      ia=6*(i-1)
      mfi=nint(scr(idict1+ia+3))
      mti=nint(scr(idict1+ia+4))
      !--insert heating mt dictionary data (or overwrite
      !  original dictionary data with revised data).
      do while ((mfi.eq.3 .and. mtn.le.mti) .or.&
                (mfi.gt.3 .and. inpk.le.npk))
         if (mfi.eq.3 .and. mtn.eq.mti) iowr=1
         b(j+1)=0
         b(j+2)=0
         b(j+3)=3
         b(j+4)=mtp(inpk)
         b(j+5)=ncds(inpk)
         b(j+6)=0
         j=j+6
         inpk=inpk+1
         if (inpk.le.npk) then
            mtn=mtp(inpk)
         else
            mtn=10000
         endif
      enddo
      !--insert original dictionary data (unless the heating
      !  mt entry above was an overwrite of existing data).
      if (mfi.ne.3 .or. (mfi.eq.3 .and. iowr.eq.0)) then
         b(j+1)=0
         b(j+2)=0
         b(j+3)=mfi
         b(j+4)=mti
         b(j+5)=scr(idict1+ia+5)
         b(j+6)=scr(idict1+ia+6)
         j=j+6
      endif
   enddo
   !--insert rest of heating mt dictionary data, if any.
   do i=inpk,npk
      b(j+1)=0
      b(j+2)=0
      b(j+3)=3
      b(j+4)=mtp(i)
      b(j+5)=ncds(i)
      b(j+6)=0
      j=j+6
   enddo
   mth=451
   nxn=j/6
   inow=1
   if (iverf.le.4) scr(inow+5)=nxn
   call contio(0,nout,0,scr(inow),nb,nw)
   inow=inow+6
   if (iverf.eq.6) call contio(0,nout,0,scr(inow),nb,nw)
   if (iverf.eq.6) inow=inow+6
   if (iverf.gt.4) then
      call contio(0,nout,0,scr(inow),nb,nw)
      inow=inow+6
      scr(inow+5)=nxn
   endif
   nw=6+nwd
   call hdatio(0,nout,0,scr(inow),nb,nw)
   inow=inow+nw
   do while (nb.ne.0)
      call moreio(0,nout,0,scr(inow),nb,nw)
      inow=inow+nw
   enddo
   nw=nxn
   call dictio(0,nout,0,b,nb,nw)
   call tofend(nin,nout,0,scr)
   call tofend(nin,nout,0,scr)

   !--add new reactions to file 3.
   call contio(nin,0,0,scr,nb,nw)
   if (mfh.ne.3 .or. mth.eq.0) then
      write(strng,'("nin out of order.  read mfh,mth = ",2i5)')mfh,mth
      call error('hout',strng,' ')
   endif
   mfnin=mfh
   mtnin=mth
   inow=7
   call contio(nscr,0,0,scr(inow),nb,nw)
   if (mfh.ne.3 .or. mth.eq.0) then
      write(strng,'("nscr out of order.  read mfh,mth = ",2i5)')mfh,mth
      call error('hout',strng,' ')
   endif
   mfnscr=mfh
   mtnscr=mth
   inow6=inow+6
   do while (mfnin.eq.3 .and. mfnscr.eq.3)
      if (mtnin.lt.mtnscr) then
         mfh=mfnin
         mth=mtnin
         call contio(0,nout,0,scr,nb,nw)
         call tosend(nin,nout,0,scr(inow6))
         call contio(nin,0,0,scr,nb,nw)
         mfnin=mfh
         mtnin=mth
      elseif (mtnin.eq.mtnscr) then
         mfh=mfnscr
         mth=mtnscr
         call contio(0,nout,0,scr(inow),nb,nw)
         call tosend(nscr,nout,0,scr(inow6))
         call tosend(nin,0,0,scr(inow6))
         call contio(nin,0,0,scr,nb,nw)
         mfnin=mfh
         mtnin=mth
         call contio(nscr,0,0,scr(inow),nb,nw)
         mfnscr=mfh
         mtnscr=mth
      else
         do while (mtnin.ge.mtnscr .and. mfnscr.eq.3)
            mfh=mfnscr
            mth=mtnscr
            call contio(0,nout,0,scr(inow),nb,nw)
            call tosend(nscr,nout,0,scr(inow6))
            call contio(nscr,0,0,scr(inow),nb,nw)
            mfnscr=mfh
            mtnscr=mth
         enddo
      endif
   enddo
   if (mfnin.eq.3 .and. mfnscr.eq.0) then
      mfh=mfnin
      mth=mtnin
      call contio(0,nout,0,scr,nb,nw)
   elseif (mfnin.eq.0 .and. mfnscr.eq.3) then
      mfh=mfnscr
      mth=mtnscr
      call contio(0,nout,0,scr(inow),nb,nw)
      call tofend(nscr,nout,0,scr(inow6))
   elseif (mfnin.eq.0 .and. mfnscr.eq.0) then
      mfh=mfnin
      mth=mtnin
      call contio(0,nout,0,scr,nb,nw)
   endif
   call tomend(nin,nout,0,scr)
   deallocate(b)
   deallocate(scr)
   call closz(nscr)

   !--plot heating, photon energy production, and kinematic limits
   if (nplot.eq.0) go to 500
   write(nplot,'(''/'')')

   !--total heating
   !--log-log
   elo=1
   elo=elo/1000
   rat=emax/elo
   thin=ten**(log10(rat)/2500)
   ylo=big
   yhi=-big
   xlast=elo
   do i=1,ne
      call finda(i,c,npkk,iold,buf,nbuf)
      x=c(1)
      if (x.ge.thin*xlast) then
         if (c(2).gt.yhi) yhi=c(2)
         if (c(2).lt.ylo) ylo=c(2)
         if (c(2+npkt).gt.yhi) yhi=c(2+npkt)
         if (c(2+npkt).lt.ylo) ylo=c(2+npkt)
         if (c(2+2*npkt).gt.yhi) yhi=c(2+2*npkt)
         if (c(2+2*npkt).lt.ylo) ylo=c(2+2*npkt)
         xlast=x
      endif
   enddo
   if (ylo.le.zero) ylo=small
   if (yhi.le.zero) yhi=small
   do k=1,3
      if (k.eq.1) then
         write(nplot,'(''1/'')')
         write(nplot,'(a,''Energy-Balance Check'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(''4 0 2 1/'')')
         write(nplot,'(''/'')')
         write(nplot,'(a,''<e>nergy (e<v>)'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(a,''<h>eating (e<v>-barns)'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(''0 0 4/'')')
         write(nplot,'(a,''<u>pper limit'',a,''/'')') qu,qu
      else if (k.eq.2) then
         write(nplot,'(''2/'')')
         write(nplot,'(''/'')')
         write(nplot,'(''/'')')
         write(nplot,'(a,''<e>nergy-balance heating'',a,''/'')') qu,qu
      else
         write(nplot,'(''3/'')')
         write(nplot,'(''/'')')
         write(nplot,'(''0 0 4/'')')
         write(nplot,'(a,''<l>ower limit'',a,''/'')') qu,qu
      endif
      write(nplot,'(''0/'')')
      xlast=elo
      j=0
      do i=1,ne
         call finda(i,c,npkk,iold,buf,nbuf)
         x=c(1)
         if (k.eq.1) y=c(2+2*npkt)
         if (k.eq.2) y=c(2)
         if (k.eq.3) y=c(2+npkt)
         if (x.ge.thin*xlast.and.j.lt.2500) then
            j=j+1
            if (y.le.zero) y=small
            write(nplot,'(1p,2e14.6,''/'')') x,y
            xlast=x
         endif
      enddo
      if (k.eq.1) write(nplot,'(1p,2e14.6,''/'')') x,ylo
      if (k.eq.1) write(nplot,'(1p,2e14.6,''/'')') x,yhi
      write(nplot,'(''/ end'')')
   enddo

   !--total heating
   !--lin-lin
   thin=emax/2500
   ylo=big
   yhi=-big
   xlast=ecut
   do i=1,ne
      call finda(i,c,npkk,iold,buf,nbuf)
      x=c(1)
      if (x.ge.xlast+thin) then
         if (c(2).gt.yhi) yhi=c(2)
         if (c(2).lt.ylo) ylo=c(2)
         if (c(2+npkt).gt.yhi) yhi=c(2+npkt)
         if (c(2+npkt).lt.ylo) ylo=c(2+npkt)
         if (c(2+2*npkt).gt.yhi) yhi=c(2+2*npkt)
         if (c(2+2*npkt).lt.ylo) ylo=c(2+2*npkt)
         xlast=x
      endif
   enddo
   do k=1,3
      if (k.eq.1) then
         write(nplot,'(''1/'')')
         write(nplot,'(a,''Energy-Balance Check'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(''1 0 2 1/'')')
         write(nplot,'(''/'')')
         write(nplot,'(a,''<e>nergy (e<v>)'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(a,''<h>eating (e<v>-barns)'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(''0 0 4/'')')
         write(nplot,'(a,''<u>pper limit>'',a,''/'')') qu,qu
      else if (k.eq.2) then
         write(nplot,'(''2/'')')
         write(nplot,'(''/'')')
         write(nplot,'(''/'')')
         write(nplot,'(a,''<e>nergy-balance heating'',a,''/'')') qu,qu
      else
         write(nplot,'(''3/'')')
         write(nplot,'(''/'')')
         write(nplot,'(''0 0 4/'')')
         write(nplot,'(a,''<l>ower limit'',a,''/'')')qu,qu
      endif
      write(nplot,'(''0/'')')
      xlast=ecut
      j=0
      do i=1,ne
         call finda(i,c,npkk,iold,buf,nbuf)
         x=c(1)
         if (k.eq.1) y=c(2+2*npkt)
         if (k.eq.2) y=c(2)
         if (k.eq.3) y=c(2+npkt)
         if (x.ge.xlast+thin.and.j.lt.2500) then
            j=j+1
            write(nplot,'(1p,2e14.6,''/'')') x,y
            xlast=x
         endif
      enddo
      if (k.eq.1) write(nplot,'(1p,2e14.6,''/'')') x,ylo
      if (k.eq.1) write(nplot,'(1p,2e14.6,''/'')') x,yhi
      write(nplot,'(''/ end'')')
   enddo

   !--photon energy production
   !--log-log
   if (mt303.eq.0) go to 490
   elo=1
   elo=elo/1000
   rat=emax/elo
   thin=ten**(log10(rat)/2500)
   ylo=big
   yhi=-big
   xlast=elo
   do i=1,ne
      call finda(i,c,npkk,iold,buf,nbuf)
      x=c(1)
      if (x.ge.thin*xlast) then
         y=c(npkk-1)+c(mt303)-c(mt303+npkt)+c(npkk-2)
         if (y.gt.yhi) yhi=y
         if (y.lt.ylo) ylo=y
         y=c(npkk-1)
         if (y.gt.yhi) yhi=y
         if (y.lt.ylo) ylo=y
         y=c(npkk-1)+c(mt303)-c(mt303+2*npkt)+c(npkk-2)
         if (y.gt.yhi) yhi=y
         if (y.lt.ylo) ylo=y
         xlast=x
      endif
   enddo
   if (ylo.le.zero) ylo=small
   if (yhi.le.zero) yhi=small
   do k=1,3
      if (k.eq.1) then
         write(nplot,'(''1/'')')
         write(nplot,'(a,''Energy-Balance Check'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(''4 0 2 1/'')')
         write(nplot,'(''/'')')
         write(nplot,'(a,''<e>nergy (e<v>)'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(a,''<p>hoton energy prod (e<v>-barns)'',&
           &a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(''0 0 4/'')')
         write(nplot,'(a,''<u>pper limit'',a,''/'')') qu,qu
      else if (k.eq.2) then
         write(nplot,'(''2/'')')
         write(nplot,'(''/'')')
         write(nplot,'(''/'')')
         write(nplot,'(a,''<p>hoton energy production'',a,''/'')')&
           qu,qu
      else
         write(nplot,'(''3/'')')
         write(nplot,'(''/'')')
         write(nplot,'(''0 0 4/'')')
         write(nplot,'(a,''<l>ower limit'',a,''/'')') qu,qu
      endif
      write(nplot,'(''0/'')')
      xlast=elo
      j=0
      do i=1,ne
         call finda(i,c,npkk,iold,buf,nbuf)
         x=c(1)
         if (k.eq.1) y=c(npkk-1)+c(mt303)-c(mt303+npkt)+c(npkk-2)
         if (k.eq.2) y=c(npkk-1)
         if (k.eq.3) y=c(npkk-1)+c(mt303)-c(mt303+2*npkt)+c(npkk-2)
         if (x.ge.thin*xlast.and.j.lt.2500) then
            j=j+1
            if (y.le.zero) y=small
            write(nplot,'(1p,2e14.6,''/'')') x,y
            xlast=x
         endif
      enddo
      if (k.eq.1) write(nplot,'(1p,2e14.6,''/'')') x,ylo
      if (k.eq.1) write(nplot,'(1p,2e14.6,''/'')') x,yhi
      write(nplot,'(''/ end'')')
   enddo

   !--photon energy production
   !--lin-lin
   thin=emax/2500
   ylo=big
   yhi=-big
   xlast=ecut
   do i=1,ne
      call finda(i,c,npkk,iold,buf,nbuf)
      x=c(1)
      if (x.ge.xlast+thin) then
         y=c(npkk-1)+c(mt303)-c(mt303+npkt)+c(npkk-2)
         if (y.gt.yhi) yhi=y
         if (y.lt.ylo) ylo=y
         y=c(npkk-1)
         if (y.gt.yhi) yhi=y
         if (y.lt.ylo) ylo=y
         y=c(npkk-1)+c(mt303)-c(mt303+2*npkt)+c(npkk-2)
         if (y.gt.yhi) yhi=y
         if (y.lt.ylo) ylo=y
         xlast=x
      endif
   enddo
   do k=1,3
      if (k.eq.1) then
         write(nplot,'(''1/'')')
         write(nplot,'(a,''Energy-Balance Check'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(''1 0 2 1/'')')
         write(nplot,'(''/'')')
         write(nplot,'(a,''<e>nergy (e<v>)'',a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(a,''<p>hoton energy prod (e<v>-barns)'',&
           &a,''/'')') qu,qu
         write(nplot,'(''/'')')
         write(nplot,'(''0 0 4/'')')
         write(nplot,'(a,''<u>pper limit'',a,''/'')') qu,qu
      else if (k.eq.2) then
         write(nplot,'(''2/'')')
         write(nplot,'(''/'')')
         write(nplot,'(''/'')')
         write(nplot,'(a,''<p>hoton energy production'',a,''/'')')&
           qu,qu
      else
         write(nplot,'(''3/'')')
         write(nplot,'(''/'')')
         write(nplot,'(''0 0 4/'')')
         write(nplot,'(a,''<l>ower limit'',a,''/'')') qu,qu
      endif
      write(nplot,'(''0/'')')
      xlast=ecut
      j=0
      do i=1,ne
         call finda(i,c,npkk,iold,buf,nbuf)
         x=c(1)
         if (k.eq.1) y=c(npkk-1)+c(mt303)-c(mt303+npkt)+c(npkk-2)
         if (k.eq.2) y=c(npkk-1)
         if (k.eq.3) y=c(npkk-1)+c(mt303)-c(mt303+2*npkt)+c(npkk-2)
         if (x.ge.xlast+thin.and.j.lt.2500) then
            j=j+1
            write(nplot,'(1p,2e14.6,''/'')') x,y
            xlast=x
         endif
      enddo
      if (k.eq.1) write(nplot,'(1p,2e14.6,''/'')') x,ylo
      if (k.eq.1) write(nplot,'(1p,2e14.6,''/'')') x,yhi
      write(nplot,'(''/ end'')')
   enddo

   !--plotting finished
  490 continue
   write(nplot,'(''99/'')')
   call closz(nplot)
  500 continue
   deallocate(buf)
   return
   end subroutine hout

end module heatm
