module dragm
   ! provides subroutine dragr for NJOY2016
   use locale
   use xsmm
   implicit none
   private
   type(xsm_file),pointer :: draglib
   public dragr

   ! random number generator for unresolved energy domain
   integer :: kk
contains

   subroutine dragr
   !-----------------------------------------------------------------
   !
   !     Produce a draglib interface file from Njoy intermediate
   !   cross-section library.
   !
   !     The draglib format provide an efficient way to store multi-
   !   group isotopic nuclear data to be used in a lattice code.
   !   For example the draglib can then be used by the Dragon code.
   !
   !     Working from thermr and groupr output tape,this module
   !   produce the standard interface file on the xsm file. The xsm
   !   data base is used to organize the data in a hierarchical
   !   format. Therefore, it will be easy to convert back and forth
   !   between the binary direct access format (efficient during a
   !   calculation) and the ascii or binary format (useful for
   !   backup and exchange purposes).
   !
   ! Authors: Hasan Saygin (original version -- 1992)
   !          Alain Hebert (reprogrammed in 2003)
   !          Alain Hebert (adapted to Njoy2012 in 2014)
   !
   ! Copyright:
   !  Copyright (C) 2003 Ecole Polytechnique de Montreal
   !  This library is free software; you can redistribute it and/or
   !  modify it under the terms of the BSD license.
   !
   !---input specifications (free format)----
   !
   ! card 1 file units
   !   nendf     input endf unit
   !   npendf    input pendf unit
   !   ngendf    input gendf unit
   !   nfp       input endf unit for fission yield data
   !   ndcy      input endf unit for radioactive decay data
   !   nimpo     input draglib unit. A negative value corresponds to
   !             a binary direct access file named drag**. Set to zero
   !             if the file is to be created.
   !   nexpo     output draglib unit. A negative value corresponds to
   !             a binary direct access file named drag**. Set to nimpo
   !             to modify its contents in place (both nimpo and nexpo
   !             must be negative in this case).
   !   pfflag    scattering storage flag: = 0/1/2. pfflag=1 to avoid
   !             storing scattering matrices; pfflag=2 to avoid
   !             storing scattering matrices above thermal energies
   !             (default=0).
   !   prflag    probability table flag: off/on = 0/1. prflag=1 if
   !             purr information is available. (default=0)
   !   idecay    fission product decay heat option (0 include, 1 don't)
   !   iprint    long print option (0/1=minimum/maximum) (default=0)
   ! card 2 hollerith identification for the library (nimpo=0 only)
   !   labell    72 character identification for the library
   ! card 3 material data (one card per material)
   !   matno     integer material identifier (endf/b mat number)
   !   hmat      hollerith material identifier (up to 8 characters
   !             each). By default, an ascii identifier is
   !             constructed.
   !   nbesp     zero (no energy-dependent fission spectra) or
   !             number of energy-dependent fission spectra)
   ! card 4 material readme comment (up to 72 characters)
   ! card 5 energy mesh for energy-dependent fission spectra
   ! present if and only if nbesp.ne.0
   !   e1, e2, e3, etc. nbesp+1 increasing values of energies (eV)
   ! card 6 energy limits for dilution-dependent xs data. This data
   !        is important to avoid self-shielding in very low and
   !        very high energy group xs and to reduce the size of the
   !        draglib. (one card per material)
   !   eres0     lower limit of resolved resonance xs data (ev)
   !   eres1     upper limit of unresolved resonance xs data (ev)
   ! card 7 thermal cutoff used by thermr, reskr, dragr and acer.
   !   ether     thermal energy cutoff (ev)
   ! card 8 autolib energy limits. This data is used with advanced
   !        self-shielding models such as the Sanchez-Coste, Ribon
   !        extended and RSE models. It is highly recommended to use
   !        the same data for all resonant isotopes. (one card per
   !        material)
   ! present if and only if npendf.ne.0
   !   eaut0     lower limit of autolib resonance xs data (ev)
   !   eaut1     upper limit of autolib resonance xs data (ev)
   !   deli      elementary lethargy width
   !
   !            repeat cards 3 to 8 for all materials desired
   !            matno=0/ terminates dragr run.
   !
   ! card 9 burnup chain data (one or two cards per isotope)
   ! present if and only if nfp.ne.0 and ndcy.ne.0. Repeat card 9 for
   ! all isotopes present in the burnup chain.
   !   hich      hollerith isotope identifier (same as hmat below)
   !             hich must be constructed in a way compatible with
   !             subroutine dranam (ex: Am242m). If an isotope is
   !             missing in this burnup chain, an error message of
   !             the type 'isotope xxx should not be lumped' may be
   !             issued by subroutine dralum.
   !   hrch(1)   hollerith neutron-induced reaction identifier for
   !             first daughter
   !   br(1)     branching ratio for an isomeric daughter (=0.000 if
   !             there is no isomeric daughter)
   !   hrch(2)   hollerith neutron-induced reaction identifier for
   !             second daughter
   !   br(2)     branching ratio for an isomeric daughter (=0.000 if
   !             there is no isomeric daughter)
   !  ...
   !
   ! card 10 optional incident neutron energies (one card per incident
   ! neutron energy). By default, the most thermal neutron energy is
   ! used.
   !
   !  'interpol' keyword to set the incident neutron energy
   !  'at'       mandatory keyword to make dragr happy
   !   energy    incident neutron energy (eV) with MT454
   !
   ! hich='end'/ terminates the burnup data.
   !
   !-----------------------------------------------------------------
   use mainio ! provides nsysi,contio,nsyso,nsyse
   use endf   ! provides endf routines and variables
   use util   ! provides timer,openz,repoz,error
   real(kr)::time,ether,eres0,eres1,eaut0,eaut1,deli,urlimit
   integer :: maxa,maxgr,maxesp
   parameter (maxa=2000,maxgr=2000,maxesp=4)
   type(xsm_file),pointer :: draglib_in
   character text72*72,text6*6,text8*8,text12*12
   integer ngen,matno,ng
   integer iesp(maxesp+1)
   character(len=4) labell(18)
   character(len=150) hsmg
   integer i,ig,igaut0,igaut1,igrest,igres0,igres1,igecco,iig,ipflag,iza, &
   & nb,nbesp,ndcy,nendf,nexpo,nfp,nimpo,npen,nw,ilong,ilong1,ityxsm,irflag, &
   & idecay,iprint,ner,igar(1),impy
   logical lsame,lurr,lkerma,lexist,lres
   real(kr) ener(maxgr+1),eesp(maxesp+1),delecco,delig,eh
   real eespi(maxesp+1),gar(1)
   real,allocatable,dimension(:) :: gar1,gar2
   real(kr),allocatable,dimension(:) :: scr
   !
   call timer(time)
   write(nsyso,'(/'' dragr...produce a draglib format output file'',24x,f8.1, &
   & ''s'')') time
   write(nsyse,'(/'' dragr...'',60x,f8.1,''s'')') time
   !
   !**allocate xsm file containing the output draglib
   allocate(draglib)
   !
   !**read users input/output unit numbers.
   nendf=0
   npen=0
   ngen=0
   nfp=0
   ndcy=0
   nimpo=0
   nexpo=0
   ipflag=0
   irflag=0
   idecay=0
   iprint=0
   read(nsysi,*) nendf,npen,ngen,nfp,ndcy,nimpo,nexpo,ipflag,irflag,idecay, &
   & iprint
   impy=max(0,iprint-1)
   allocate(scr(maxa))
   !
   !**open the input/output files
   write(nsyso,'(/ &
   &  '' input endf unit ............................. '',i10/ &
   &  '' input pendf unit ............................ '',i10/ &
   &  '' input gendf unit ............................ '',i10/ &
   &  '' input endf unit for fission yield data ...... '',i10/ &
   &  '' input endf unit for radioactive decay data .. '',i10/ &
   &  '' input draglib unit .......................... '',i10/ &
   &  '' output draglib unit ......................... '',i10/ &
   &  '' scattering storage flag ..................... '',i10/ &
   &  '' purr information flag ....................... '',i10/ &
   &  '' FP decay heat option (0 include, 1 dont) .... '',i10/ &
   &  '' print option (0 min, 1 max) ................. '',i10)') &
   &  nendf,npen,ngen,nfp,ndcy,nimpo,nexpo,ipflag,irflag,idecay,iprint
   if(nendf.ne.0) call openz(nendf,0)
   if(npen.ne.0) call openz(npen,0)
   if(ngen.ne.0) call openz(ngen,0)
   if((nfp.ne.0).and.(ndcy.ne.0)) then
     call openz(nfp,0)
     call openz(ndcy,0)
   endif
   !
   !**recover the input draglib on unit nimpo and copy it on unit nexpo
   lexist=.false.
   if(nexpo.eq.0) then
     call error('dragr','output draglib unit cannot be zero',' ')
   else if(nimpo.eq.0) then
     !**open a new xsm file and create the signature
     write(text6,'(4hdrag,i2.2)') abs(nexpo)
     inquire(file=text6,exist=lres)
     if(lres) then
       write(hsmg,'(13hdraglib file ,a6,15h already exists)') text6
       call error('dragr',hsmg,' ')
     endif
     call xsmop(draglib,text6,0,1)
     text12='L_DRAGLIB'
     read(text12,'(3a4)') (labell(i),i=1,3)
     labell(1)=text12(1:4)
     labell(2)=text12(5:8)
     labell(3)=text12(9:12)
     call xsmput(draglib,'SIGNATURE',labell(1:3))
     text12='RELEASE_2025'
     labell(1)=text12(1:4)
     labell(2)=text12(5:8)
     labell(3)=text12(9:12)
     call xsmput(draglib,'VERSION',labell(1:3))
   else if((nimpo.lt.0).and.(nexpo.eq.nimpo)) then
     !**open an existing direct access draglib in modification mode
     write(text6,'(4hdrag,i2.2)') -nimpo
     call xsmop(draglib,text6,1,1)
     call xsmlen(draglib,'ENERGY',ilong1,ityxsm)
     if(ilong1.eq.0) then
       call xsmlib(draglib)
       call error('dragr','no ENERGY record on existing draglib(1)',' ')
     endif
     allocate(gar1(ilong1))
     call xsmget(draglib,'ENERGY',gar1)
     lexist=.true.
   else if((nimpo.lt.0).and.(nexpo.lt.0).and.(nexpo.ne.nimpo)) then
     !**open an existing direct access draglib in read-only mode
     write(text6,'(4hdrag,i2.2)') -nimpo
     allocate(draglib_in)
     call xsmop(draglib_in,text6,2,1)
     call xsmlen(draglib_in,'ENERGY',ilong1,ityxsm)
     if(ilong1.eq.0) then
       call xsmlib(draglib_in)
       call error('dragr','no ENERGY record on existing draglib(2)',' ')
     endif
     allocate(gar1(ilong1))
     call xsmget(draglib_in,'ENERGY',gar1)
     !
     write(text6,'(4hdrag,i2.2)') -nexpo
     inquire(file=text6,exist=lres)
     if(lres) then
       call xsmop(draglib,text6,1,1) ! open draglib in modification mode
     else
       call xsmop(draglib,text6,0,1) ! open draglib in creation mode
     endif
     call xsmequ(draglib_in,draglib)
     call xsmcl(draglib_in,1)
     deallocate(draglib_in)
     lexist=.true.
   else
     !**open and import an existing draglib
     if(nexpo.lt.0) then
       write(text6,'(4hdrag,i2.2)') -nexpo
       inquire(file=text6,exist=lres)
       if(lres) then
         call xsmop(draglib,text6,1,1) ! open draglib in modification mode
       else
         call xsmop(draglib,text6,0,1) ! open draglib in creation mode
       endif
       call openz(nimpo,0)
       call xsmexp(draglib,nimpo,2,2,impy)
       call closz(nimpo)
     else if(nexpo.gt.0) then
       if(nimpo.lt.0) then
         write(text6,'(4hdrag,i2.2)') -nimpo
         call xsmop(draglib,text6,1,1)
       else
         write(text6,'(4hdrag,i2.2)') nexpo
         call xsmop(draglib,text6,0,1)
         call openz(nimpo,0)
         call xsmexp(draglib,nimpo,2,2,impy) ! import
         call closz(nimpo)
       endif
       write(text6,'(4htape,i2.2)') nexpo
       inquire(file=text6,exist=lres)
       if(lres) then
         !** concatenate nimpo to existing nexpo
         call openz(nexpo,0)
         allocate(draglib_in)
         call xsmop(draglib_in,'dummy',0,1)
         call xsmexp(draglib_in,nexpo,2,2,impy) ! import
         call xsmequ(draglib,draglib_in)
         call xsmcl(draglib,2)
         deallocate(draglib)
         draglib=>draglib_in
         call closz(nexpo)
       endif
     endif
     call xsmlen(draglib,'ENERGY',ilong1,ityxsm)
     if(ilong1.eq.0) then
       call xsmlib(draglib)
       call error('dragr','no ENERGY record on existing draglib(3)',' ')
     endif
     allocate(gar1(ilong1))
     call xsmget(draglib,'ENERGY',gar1)
     lexist=.true.
   endif
   !**skip draglib feeding in trivial case
   if((nendf.eq.0).and.(npen.eq.0).and.(ngen.eq.0).and. &
     &(nfp.eq.0).and.(ndcy.eq.0)) go to 30
   !
   !**read identification for the library
   if(nimpo.eq.0) then
     text72=' '
     read(nsysi,*) text72
     if(text72.ne.' ') then
       do i=1,18
         labell(i)=text72((i-1)*4+1:i*4)
       enddo
       write(nsyso,'('' identification for the library .............. ''/1x, &
       & a)') text72
       call xsmput(draglib,'README',labell(1:18))
     endif
   endif
   !
   !**material loop
   if(ngen.eq.0) go to 20
   do
     matno=0
     text8=' '
     read(nsysi,*) matno,text8,nbesp
     if(nbesp.gt.maxesp) call error('dragr','maxesp overflow',' ')
     if(matno.eq.0) then
       go to 10
     else if(text8.eq.' ') then
       ! automatic composition of the ascii material name.
       if(nendf.eq.0) call error('dragr','endf unit not provided',' ')
       call findf(matno,1,451,nendf)
       call contio(nendf,0,0,scr,nb,nw)
       iza=nint(c1h+0.1)
       call contio(nendf,0,0,scr,nb,nw)
       call dranam(10*iza+l2h,text8)
     endif
     write(nsyso,'('' --> processing material '',a, &
     &  ''  (endf identification ='',i8,'')'')') text8,matno
     !
     urlimit=0.0
     if(nendf.ne.0) then
       ! ***recover the upper limit of the resolved energy domain.
       call repoz(nendf)
       call findf(matno,2,151,nendf)
       call contio(nendf,0,0,scr,nb,nw)
       call contio(nendf,0,0,scr,nb,nw)
       ner=n1h
       lurr=.false.
       eh=-1.0e10
       do
         if((c1h.eq.0.0).and.(c2h.eq.0.0).and.(n1h.eq.0).and.(n2h.eq.0)) exit
         call contio(nendf,0,0,scr,nb,nw)
         if((l1h.eq.0).and.(l2h.eq.0)) then
           urlimit=c2h
           exit
         else if((l1h.eq.1).and.(l2h.gt.0)) then
           eh=c2h
           urlimit=c2h
         else if(c1h.eq.eh) then
           lurr=.true.
           exit
         endif
         call listio(nendf,0,0,scr,nb,nw)
         if(nw.gt.maxa) call error('dragr','endf input size exceeded(1)',' ')
         do while (nb.ne.0)
           call moreio(nendf,0,0,scr,nb,nw)
           if(nw.gt.maxa) call error('dragr','endf input size exceeded(2)',' ')
         enddo
       enddo
       write(nsyso,'('' number of resonance ranges .................. '', &
       & i10)') ner
       write(nsyso,'('' resolved energy domain upper limit .......'',1p, &
       & e11.4,'' eV'')') urlimit
       write(nsyso,'('' unresolved energy domain exists ............. '',1p, &
       & l10)') lurr
       if(.not.lurr) irflag=0
     endif
     !
     ! ***energy group information for the library
     call drahd(ngen,matno,ng,maxgr,ener)
     text72=' '
     read(nsysi,*) text72
     do i=1,18
       labell(i)=text72((i-1)*4+1:i*4)
     enddo
     write(nsyso,'('' identification for the material ............. ''/1x, &
     & a)') text72
     !
     ! ***recover energy limits for energy-dependent fission spectra
     if(nbesp.ne.0) then
       read(nsysi,*) (eesp(i),i=1,nbesp+1)
       iig=nbesp+1
       do ig=1,ng+1
         if(iig.eq.0) call error('dragr','bad limits for energy-dep' &
         &  //'endent fission spectra(1)',' ')
         if(eesp(iig).ge.0.999d0*ener(ig)) then
           iesp(nbesp+2-iig)=ig-1
           iig=iig-1
         endif
       enddo
       if(iig.ne.0) call error('dragr','bad limits for energy-depen' &
       & //'dent fission spectra(2)',' ')
       do iig=1,nbesp+1
         eespi(iig)=real(eesp(nbesp+2-iig))
       enddo
       call xsmput(draglib,'CHI-ENERGY',eespi(1:nbesp+1))
       call xsmput(draglib,'CHI-LIMITS',iesp(1:nbesp+1))
     else
       nbesp=1
       iesp(1)=0
       iesp(2)=ng
     endif
     !
     ! ***step on material
     call xsmsix(draglib,text8,1)
     call xsmput(draglib,'README',labell(1:18))
     !
     ! ***recover energy limits for dilution-dependent xs data
     read(nsysi,*) eres0,eres1
     if(eres0.eq.0.0d0) eres0=ener(ng+1)
     if(eres1.eq.0.0d0) eres1=ener(1)
     !
     ! ***recover thermal cutoff in eV
     read(nsysi,*) ether
     !
     ! ***compute igrest, igres0 and igres1
     igrest=0
     igres0=1
     igres1=ng
     do ig=ng,2,-1
       if(ether.ge.0.999d0*ener(ig)) igrest=igrest+1
     enddo
     do ig=1,ng
       if(eres1.le.1.001d0*ener(ig+1)) igres0=ig+1
     enddo
     do ig=ng+1,2,-1
       if(eres0.ge.0.999d0*ener(ig)) igres1=ig-1
     enddo
     write(nsyso,'(/ &
     &  '' number of thermal groups .................... '',i10/ &
     &  '' first group of unresolved resonance xs data . '',i10/ &
     &  '' last group of resolved resonance xs data .... '',i10)') &
     &  igrest,igres0,igres1
     !
     ! ***validate energy mesh for option ipflag=2 (eccolib mesh)
     igecco=0
     if(ipflag.eq.2) then
       delecco=log(ener(1)/ener(2))
       do ig=1,ng
         if(ig.gt.ng-igrest+1) exit
         delig=log(ener(ig)/ener(ig+1))
         if(abs(delig-delecco).le.1.0e-5) then
           igecco=ig
         else if(ener(ig+1).gt.10.0) then
           write(hsmg,'(44hunable to use ipflag=2 (not an eccolib mesh))')
           call mess('dragr',hsmg,' ')
           write(hsmg,'(23hlethargy width of group,i5,2h (,1p,e10.2,7h) is no, &
           & 29ht equal to width of group 1 (,e10.2,1h))') ig,delig,delecco
           call error('dragr',hsmg,' ')
         else
           exit
         endif
       enddo
       write(nsyso,'( &
       &    '' last eccolib group with equal lethargy width  '',i10/ &
       &    '' eccolib lethargy width ...................... '',1p,d10.4)') &
       &    igecco,delecco
       if((igecco.lt.igres1).and.(igres1.lt.ng)) then
           write(hsmg,'(38hlast group with equal lethargy width (,i5,4h) is, &
           & 33h lower than last resolved group (,i5,1h))') igecco,igres1
           call error('dragr',hsmg,' ')
       endif
       igar(1)=igecco
       call xsmput(draglib,'ONFLIGHTIGR',igar)
       gar(1)=real(delecco)
       call xsmput(draglib,'ONFLIGHTDEL',gar)
     endif
     !
     ! ***generate draglib file for this material
     call dramat(nendf,ngen,matno,ng,igrest,igres0,igres1,igecco,ipflag,nbesp, &
     & iesp,ener,lkerma)
     write(nsyso,'(/ &
     &    '' heatr kerma factor availability ............. '',l10)') lkerma
     if(nendf.ne.0) then
       ! ***construct the "isotopic specification line" for the burnup data.
       call dradat(nendf,matno,text8)
     endif
     !
     ! ***put autolib data in draglib file for this material
     if(npen.ne.0) then
       read(nsysi,*) eaut0,eaut1,deli
       igaut0=1
       igaut1=ng
       do ig=1,ng
         if(eaut1.le.1.001d0*ener(ig+1)) igaut0=ig+1
       enddo
       do ig=ng+1,2,-1
         if(eaut0.ge.0.999d0*ener(ig)) igaut1=ig-1
       enddo
       write(nsyso,'(/ &
       &    '' first group of autolib xs data.............. '',i10/ &
       &    '' last group of autolib xs data............... '',i10/ &
       &    '' elementary lethargy width................... '',1p,d10.4)') &
       &    igaut0,igaut1,deli
       if((eaut1.gt.urlimit).and.(irflag.eq.0).and.(ner.gt.1)) then
         write(hsmg,'(48h The upper limit of the resolved energy domain (, &
         & 1p,e10.3,51h eV) is lower that the upper limit of the autolib (, &
         & 1p,e10.3,5h eV).)') urlimit,eaut1
         call mess('dragr',hsmg,' ')
       endif
       call drauto(npen,matno,ng,ener,eaut0,eaut1,igaut0,igaut1,deli,irflag, &
       & urlimit,iprint)
     endif
     call xsmsix(draglib,' ',2)
   enddo
   !
   ! **finished
   10 call closz(ngen)
   call closz(npen)
   call closz(nendf)
   !
   ! **compute depletion-related data
   20 if((nfp.ne.0).and.(ndcy.ne.0)) then
     call xsmsix(draglib,'DEPL-CHAIN',1)
     call dradep(nfp,ndcy,idecay,iprint)
     call xsmsix(draglib,' ',2)
     call closz(ndcy)
     call closz(nfp)
   endif
   !
   ! **check the energy mesh consistency
   30 if(lexist) then
     call xsmlen(draglib,'ENERGY',ilong,ityxsm)
     if(ilong.eq.0) then
       call xsmlib(draglib)
       call error('dragr','no ENERGY record on new draglib',' ')
     endif
     allocate(gar2(ilong))
     call xsmget(draglib,'ENERGY',gar2)
     if(ilong.ne.ilong1)  call error('dragr','number of energy' &
     & //' group do not match for concatenate',' ')
     lsame=.true.
     do ig=1,ilong-1
       lsame=lsame.and.(abs(gar1(ig)-gar2(ig)).le.1.0e-3*gar1(ig))
     enddo
     if(.not.lsame) then
       do ig=1,ilong-1
         write(nsyso,*) ig,gar1(ig),gar2(ig), &
         & abs(gar1(ig)-gar2(ig))/gar1(ig)*100,'%'
       enddo
       call error('dragr','inconsistent energy mesh',' ')
     endif
     deallocate(gar1,gar2)
   endif
   !
   !** export/close the output draglib
   if(nexpo.lt.0) then
     call xsmcl(draglib,1) ! keep the file
   else if(nexpo.gt.0) then
     ! **open and export the new ascii draglib
     call openz(nexpo,1)
     call xsmexp(draglib,nexpo,2,1,impy)
     call closz(nexpo)
     if(nimpo.lt.0) then
       call xsmcl(draglib,1) ! close the direct access file
     else
       call xsmcl(draglib,2) ! erase the direct access file
     endif
   endif
   deallocate(draglib)
   !
   call closz(nexpo)
   call timer(time)
   write(nsyso,'(69x,f8.1,''s''/1x,7(''**********''),''*******'')') time
   deallocate(scr)
   return
   end subroutine dragr
   !
   subroutine dradat(nendf,matno,text8)
   !-----------------------------------------------------------------
   ! Construct the "isotopic specification line" for the burnup data.
   !-----------------------------------------------------------------
   use mainio ! provides nsysi,contio,nsyso,nsyse
   use endf   ! provides endf routines and variables
   use util   ! provides timer,openz,repoz,error
   integer :: maxa,maxedi
   parameter(maxa=2000,maxedi=12)
   integer nendf,matno
   character hline*80,text8*8
   integer ied,iof,iof2,iza,nb,nw
   logical lfind
   integer,save,dimension(maxedi) :: malist1= &
   & (/ 16, 17, 18, 28, 37,102,103,104,105,106,107,108 /)
   integer,save,dimension(maxedi) :: malist2= &
   & (/-99,-99,-99,-99,-99,-99,600,650,700,750,800,-99 /)
   character(len=8),save,dimension(maxedi) :: namedi= &
   & (/'n2n     ','n3n     ','nftot   ','nnp     ','n4n     ','ng      ', &
   &   'np      ','nd      ','nt      ','nhe3    ','na      ','n2a     '/)
   real(kr),allocatable,dimension(:) :: scr
   !
   allocate(scr(maxa))
   write(nsyso,'(/'' isotopic specification line for material'',i8, &
   & '':''/1x,80(''-''))') matno
   hline=' '
   iof=1
   call repoz(nendf)
   call findf(matno,1,451,nendf)
   call contio(nendf,0,0,scr,nb,nw)
   iza=nint(c1h+0.1)
   call contio(nendf,0,0,scr,nb,nw)
   write(hline(iof:),'(a8,1x)') text8
   iof=iof+9
   do ied=1,maxedi
     call repoz(nendf)
     lfind=.false.
     do while (.not.lfind)
       if(nendf.lt.0) then
         read(-nendf,end=50) math,mfh,mth,nb,nw
       else if(nendf.gt.0) then
         read(nendf,'(66x,i4,i2,i3,i5)',end=50) math,mfh,mth,nsp
       endif
       lfind=(math.eq.matno).and.(mfh.eq.3).and.((mth.eq.malist1(ied)).or. &
       & (mth.eq.malist2(ied)))
     enddo
     if(.not.lfind) go to 50
     if(iof+15.gt.80) then
       write(nsyso,'(1x,a80)') hline
       hline=' '
       iof=10
     endif
     iof2=index(namedi(ied),' ')
     write(hline(iof:),'(a,1x)') namedi(ied)(:iof2-1)
     iof=iof+iof2
     ! the isomeric branching ratio is set to 0.0. It may be corrected before
     ! calling dradep()
     write(hline(iof:),'('' 0.000 '')')
     iof=iof+7
     50 continue
   enddo
   if(iof+1.gt.80) then
     write(nsyso,'(1x,a80)') hline
     hline=' '
     iof=10
   endif
   write(hline(iof:),'(''/'')')
   write(nsyso,'(1x,a80)') hline
   write(nsyso,'(1x,80(''-''))')
   deallocate(scr)
   return
   end subroutine dradat
   !
   subroutine drahd(ngen,matno,ng,maxgr,ener)
   !-----------------------------------------------------------------
   ! write energy group data to the xsm file
   !-----------------------------------------------------------------
   use endf   ! provides endf routines and variables
   use util   ! provides error
   integer :: maxa
   parameter(maxa=3000)
   integer ngen,matno,ng,maxgr
   real(kr) ener(maxgr+1)
   integer ig,ilong,ityxsm,loc,nb,nw,nz
   logical lsame
   real,allocatable,dimension(:) :: gar
   real(kr),allocatable,dimension(:) :: scr
   !
   allocate(scr(maxa))
   call findf(matno,1,451,ngen)
   call contio(ngen,0,0,scr,nb,nw)
   nz=nint(scr(4))
   call listio(ngen,0,0,scr(1),nb,nw)
   loc=1+nw
   do while (nb.ne.0)
     if(loc+302.gt.maxa) call error('drahd','endf input size exceeded',' ')
     call moreio(ngen,0,0,scr(loc),nb,nw)
     loc=loc+nw
   enddo
   ng=nint(scr(3))
   if(ng.gt.maxgr) call error('drahd','maxgr overflow',' ')
   call xsmlen(draglib,'ENERGY',ilong,ityxsm)
   allocate(gar(ng+1))
   if(ilong.eq.0) then
     do ig=1,ng+1
       ener(ng-ig+2)=scr(7+nz+ig)
     enddo
     gar(1:ng+1)=real(ener(1:ng+1))
     call xsmput(draglib,'ENERGY',gar(1:ng+1))
   else
     if(ilong.ne.ng+1) call error('drahd','invalid number of groups',' ')
     call xsmget(draglib,'ENERGY',gar)
     ener(1:ng+1)=gar(1:ng+1)
     lsame=.true.
     do ig=1,ng
       lsame=lsame.and.(abs(ener(ng-ig+2)-scr(7+nz+ig)).le.1.0e-3*ener(ng-ig+2))
     enddo
     if(.not.lsame) call error('drahd','inconsistent energy mesh',' ')
   endif
   deallocate(gar)
   deallocate(scr)
   return
   end subroutine drahd
   !
   subroutine dramat(nendf,ngen,matno,ng,igrest,igres0,igres1,igecco,ipflag, &
   & nbesp,iesp,ener,lkerma)
   !-----------------------------------------------------------------
   !   write draglib data to the xsm (direct access) file for endf
   !   material matno
   !   the dilution-dependent xs deltas are recovered between groups
   !   igres0 and igres1.
   !-----------------------------------------------------------------
   use endf   ! provides endf routines and variables
   use util   ! provides error
   integer :: maxa,maxgr,maxnl,maxnz,maxtmp,maxedi,maxedi2,lz
   parameter (maxa=3000,maxgr=2000,maxnl=8,maxnz=30,maxtmp=100,maxedi=17, &
   & maxedi2=57,lz=6)
   integer nendf,ngen,matno,ng,igrest,igres0,igres1,igecco,ipflag,nbesp, &
   & iesp(nbesp+1)
   real(kr) aa(6),ener(ng+1)
   logical lkerma,lfind,lover,exist,exist2,exist3,lfiss
   character cd*4,hsmg*131,hmt*8
   real awr(1),dilut(maxnz),olddil(maxnz),oldtmp(maxtmp),vector(maxgr)
   integer igfirs(maxnl),iglast(maxnl)
   integer ied,ig,ig1,ig2,igmax,il,ilong,ityxsm,imt,isbmat,itm,itm0,iz,iz0, &
   & loc,nb,ng0,nl,nl2,nlgar,ntmp,nw,nz,nz0,nz0bis,nzgar,ned,alloc_ok
   real(kr) temps(maxtmp),ytemp,kap,kap2
   integer, save, dimension(maxedi) :: malist= &
   & (/ 1,2,4,5,16,17,18,28,37,102,103,104,105,106,107,108,301 /)
   integer, save, dimension(maxedi2) :: malist2= &
   & (/ 1,2,5,16,17,18,28,37,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66, &
   & 67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90, &
   & 91,102,103,104,105,106,107,108,301 /)
   character(len=8), save, dimension(maxedi) :: namedi= &
   & (/ 'NTOT0   ','NELAS   ','NINEL   ','NX      ','N2N     ','N3N     ', &
   &    'NFTOT   ','NNP     ','N4N     ','NG      ','NP      ','ND      ', &
   &    'NT      ','NHE3    ','NA      ','N2A     ','H-FACTOR'/)
   character(len=8), save, dimension(maxedi2) :: namedi2= &
   & (/ 'NTOT0   ','NELAS   ','NX      ','N2N     ','N3N     ','NFTOT   ', &
   &    'NNP     ','N4N     ','NINEL001','NINEL002','NINEL003','NINEL004', &
   &    'NINEL005','NINEL006','NINEL007','NINEL008','NINEL009','NINEL010', &
   &    'NINEL011','NINEL012','NINEL013','NINEL014','NINEL015','NINEL016', &
   &    'NINEL017','NINEL018','NINEL019','NINEL020','NINEL021','NINEL022', &
   &    'NINEL023','NINEL024','NINEL025','NINEL026','NINEL027','NINEL028', &
   &    'NINEL029','NINEL030','NINEL031','NINEL032','NINEL033','NINEL034', &
   &    'NINEL035','NINEL036','NINEL037','NINEL038','NINEL039','NINEL040', &
   &    'NINEL041','NG      ','NP      ','ND      ','NT      ','NHE3    ', &
   &    'NA      ','N2A     ','H-FACTOR'/)
   real(kr),allocatable,dimension(:) :: scr,deltau
   real(kr),allocatable,dimension(:,:,:) :: rv2,rv3,rv4,flux,rv
   real(kr),allocatable,dimension(:,:,:,:) :: rm2
   !
   ! scratch storage allocation
   allocate(scr(maxa))
   allocate(deltau(maxgr),rv2(maxgr,maxnl,1),flux(maxgr,maxnl,maxnz), &
   & rv(maxgr,maxnl,maxnz))
   allocate(rm2(maxgr,maxgr,maxnl,1),stat=alloc_ok)
   if(alloc_ok /= 0) call error('dramat','unable to allocate rm2',' ')
   !
   ! recover the new dilution values.
   call findf(matno,1,451,ngen)
   call contio(ngen,0,0,scr,nb,nw)
   awr(1)=real(scr(2))
   call xsmput(draglib,'AWR',awr)
   nz=nint(scr(4))
   if(nz.gt.maxnz) call error('dramat','maxnz overflow',' ')
   call listio(ngen,0,0,scr(1),nb,nw)
   loc=1+nw
   do while(nb.ne.0)
     if(loc+302.gt.maxa) call error('dramat','endf input size exceeded',' ')
     call moreio(ngen,0,0,scr(loc),nb,nw)
     loc=loc+nw
   enddo
   do iz=1,nz
     dilut(nz-iz+1)=real(scr(7+iz))
   enddo
   if(dilut(nz).lt.1.0e10) then
     call error('dramat','missing infinite dilution value',' ')
   endif
   do ig=1,ng
     deltau(ig)=log(scr(7+nz+(ng-ig+2))/scr(7+nz+(ng-ig+1)))
   enddo
   !
   ! recover the temperature values.
   ntmp=0
   10 lfind=.false.
   do while(.not.lfind)
     if(ngen.lt.0) then
       read(-ngen,end=20) math,mfh,mth,nb,nw
     else if(ngen.gt.0) then
       read(ngen,'(6e11.0,i4,i2,i3,i5)',end=20) aa,math,mfh,mth,nsp
     endif
     lfind=(math.eq.matno).and.(mfh.eq.3).and.(mth.eq.1)
   enddo
   if(.not.lfind) go to 20
   call listio(ngen,0,0,scr(1),nb,nw)
   ntmp=ntmp+1
   if(ntmp.gt.maxtmp) call error('dramat','maxtmp overflow',' ')
   temps(ntmp)=scr(1)
   call tomend(ngen,0,0,scr)
   go to 10
   20 if(ntmp.eq.0) call error('dramat','no temperatures on gendf',' ')
   !
   ! set the 'TEMPERATURE' record on xsm.
   call xsmlen(draglib,'TEMPERATURE',ilong,ityxsm)
   if(ilong.gt.0) then
     if(ilong.gt.maxtmp) call error('dramat','maxtmp overflow',' ')
     call xsmget(draglib,'TEMPERATURE',oldtmp)
     lover=.true.
     do itm=1,ntmp
       lover=lover.and.(temps(itm).gt.oldtmp(ilong))
     enddo
     if(lover) then
       if(ntmp+ilong.gt.maxtmp) call error('dramat','maxtmp overflow',' ')
       do itm=1,ntmp
          ilong=ilong+1
          oldtmp(ilong)=real(temps(itm))
       enddo
     endif
   else
     ilong=ntmp
     do itm=1,ntmp
        oldtmp(itm)=real(temps(itm))
     enddo
   endif
   call xsmput(draglib,'TEMPERATURE',oldtmp(1:ilong))
   !
   ! loop over the temperatures.
   do itm=1,ntmp
     do itm0=1,ilong
       if(abs(temps(itm)-oldtmp(itm0)).le.1.0e-5*temps(itm)) then
         isbmat=itm0
         go to 30
       endif
     enddo
     call error('dramat','missing temperature',' ')
     30   ytemp=temps(itm)
     write (cd,'(i4.4)') isbmat
     call xsmsix(draglib,'SUBTMP'//cd,1)
     !
     ! ***store the new dilutions
     call xsmlen(draglib,'DILUTION',nz0,ityxsm)
     if(nz0.gt.0) then
       if(nz0.gt.maxnz) call error('dramat','maxnz overflow',' ')
       call xsmget(draglib,'DILUTION',olddil)
       lover=.true.
       do iz=1,nz-1
         lover=lover.and.(dilut(iz).gt.olddil(nz0))
       enddo
       if(lover) then
         if(nz+nz0-1.gt.maxnz) call error('dramat','maxnz overflow',' ')
         do iz=1,nz-1
           olddil(nz0+iz)=dilut(iz)
         enddo
         call xsmput(draglib,'DILUTION',olddil(1:nz+nz0-1))
       else
         nz0bis=0
         do iz0=1,nz0
           if(abs(dilut(1)-olddil(iz0)).le.1.0e-5*dilut(1)) then
             nz0bis=iz0-1
             go to 35
           endif
         enddo
         call error('dramat','inconsistent dilutions(1).',' ')
         35 if(nz0bis+nz-1.gt.nz0) then
           call error('dramat','inconsistent dilutions(2).',' ')
         endif
         do iz=1,nz-1
           if(abs(dilut(iz)-olddil(nz0bis+iz)).gt.1.0e-5*dilut(iz)) &
           & call error('dramat','inconsistent dilutions(3).',' ')
         enddo
         nz0=nz0bis
       endif
     else if(nz.gt.1) then
       call xsmput(draglib,'DILUTION',dilut(1:nz-1))
     endif
     !
     !   ***recover thermal correction x-sections and compute the
     !      correction term.
     !   NOTE: Apparently, there is a bug in THERMR so that the matrix
     !         221 to 250 reactions does not sum to the corresponding
     !         vector reactions. We will therefore avoid using the
     !         vector reactions for thermal scattering.
     !    call fvect(maxgr,maxnl,1,ngen,matno,imt,ytemp,ng,nl2,nzgar,rv2,
     ! &  flux,exist2)
     rv2(:maxgr,:maxnl,1)=0.0
     exist3=.false.
     do imt=250,221,-1
       call fmatr(maxgr,maxnl,1,ngen,matno,imt,ytemp,ng0,nl2,nzgar,rm2,flux, &
       & exist2)
       if(exist2) then
         exist3=.true.
         nlgar=nl2
         if(ng0.ne.ng) call error('dramat','inconsistent ng(1).',' ')
         do il=1,nl2
           do ig1=1,ng
             do ig2=1,ng
               rv2(ig1,il,1)=rv2(ig1,il,1)+rm2(ig1,ig2,il,1)
             enddo
           enddo
         enddo
       endif
     enddo
     if(exist3) then
       call fvect(maxgr,maxnl,maxnz,ngen,matno,2,ytemp,ng0,nl,nzgar, &
       & rv,flux,exist)
       if(nl.gt.nlgar) call error('dramat','inconsistent number o'// &
       & 'f Legendre orders(1)',' ')
       do il=1,nl
         do ig=1,ng
           if(rv2(ig,il,1).ne.0.0) then
             rv2(ig,il,1)=rv2(ig,il,1)-rv(ig,il,1)
           endif
         enddo
       enddo
     endif
     !
     ! ***loop over the vector cross section sets.
     lkerma=.false.
     lfiss=.false.
     igfirs(:maxnl)=0
     iglast(:maxnl)=0
     if(igecco.eq.0) then
       ned=maxedi
     else
       ned=maxedi2
     endif
     do ied=1,ned
       !  ***process x-sections for reaction malist(ied).
       if(igecco.eq.0) then
         ! store the sum of inelastic scattering channels
         imt=malist(ied)
         hmt=namedi(ied)
       else
         ! store all channels of inelastic scattering (used for eccolibs)
         imt=malist2(ied)
         hmt=namedi2(ied)
       endif
       call fvect(maxgr,maxnl,maxnz,ngen,matno,imt,ytemp,ng0,nl,nzgar,rv,flux, &
       & exist)
       if(.not.exist) cycle
       if(ng0.ne.ng) call error('dramat','inconsistent ng(2).',' ')
       if(imt.eq.18) lfiss=.true.
       !
       if((imt.le.2).and.exist3) then
         ! ***perform thermal correction on 'NTOT0' and 'NELAS' x-sections.
         do iz=1,nzgar
           do ig=1,ng
             rv(ig,1,iz)=rv(ig,1,iz)+rv2(ig,1,1)
           enddo
         enddo
       else if(imt.eq.301) then
         ! ***add H-FACTOR information on Draglib based on heatr kermas.
         lkerma=.true.
         allocate(rv3(maxgr,maxnl,maxnz),rv4(maxgr,maxnl,maxnz))
         if(lfiss) then
           call fvect(maxgr,maxnl,maxnz,ngen,matno,318,ytemp,ng0,nl,nzgar, &
           & rv3,flux,exist)
           !
           !   ***recompute the fission kerma according to the Serpent2 recipe
           !   ***fission -- recover kappa info in mf=1, mt=458
           call repoz(nendf)
           lfind=.false.
           do while (.not.lfind)
             if(nendf.lt.0) then
               read(-nendf,end=50) math,mfh,mth,nb,nw
             else if(nendf.gt.0) then
               read(nendf,'(66x,i4,i2,i3,i5)',end=50) math,mfh,mth,nsp
             endif
             if(math.eq.-1) then
               call repoz(nendf)
               go to 50
             endif
             lfind=(math.eq.matno).and.(mfh.eq.1).and.(mth.eq.458)
           enddo
           if(.not.lfind) go to 50
           call listio(nendf,0,0,scr,nb,nw)
           ! ***remove promp and delayed neutron energy from pseudo Q
           kap=scr(lz+15)-scr(lz+3)-scr(lz+5)
           call fvect(maxgr,maxnl,maxnz,ngen,matno,18,ytemp,ng0,nl,nzgar, &
           & rv4,flux,exist)
           do iz=1,nzgar
             do ig=1,ng
               kap2=kap+0.5*(ener(ig)+ener(ig+1)) ! add incident neutron energy
               rv(ig,1,iz)=rv(ig,1,iz)-rv3(ig,1,iz)+kap2*rv4(ig,1,iz)
             enddo
           enddo
         endif
         50 deallocate(rv4,rv3)
       endif
       !
       igmax=ng
       do ig=ng,1,-1
         igmax=ig
         if(rv(ig,1,1).ne.0.0) go to 60
       enddo
       cycle
       60 do ig=1,igmax
         vector(ig)=real(rv(ig,1,1))
       enddo
       call xsmput(draglib,hmt,vector(1:igmax))
       if(nzgar.eq.1) cycle
       !
       ! ***find the last self-shielded group.
       do il=1,maxnl
         igfirs(il)=igres0
         iglast(il)=igres1+1
       enddo
       do iz=1,nz-1
         write (cd,'(i4.4)') nz0+iz
         call xsmsix(draglib,'SUBMAT'//cd,1)
         !
         ! ***process finite dilution x-sections
         vector(1:iglast(1)-1)=0.0
         igmax=ng
         do ig=iglast(1)-1,igfirs(1),-1
           igmax=ig
           if((flux(ig,1,1).ne.0.0).or.(rv(ig,1,1).ne.0.0)) go to 70
         enddo
         go to 80
         70 do ig=igfirs(1),igmax
           flux(ig,1,nz-iz+1)=flux(ig,1,nz-iz+1)/flux(ig,1,1)
           if(flux(ig,1,nz-iz+1).gt.10.0) then
             write(hsmg,'(19hinconsistent flux (,1p,e10.3,5h) in , &
             & 5hgroup,i5,13h at dilution=,e10.2,5h barn)') &
             & flux(ig,1,nz-iz+1),ig,dilut(iz)
             call error('dramat',hsmg,' ')
           endif
           vector(ig)=real(rv(ig,1,nz-iz+1)*flux(ig,1,nz-iz+1)-rv(ig,1,1))
         enddo
         call xsmput(draglib,hmt,vector(1:iglast(1)-1))
         !
         80 if(ied.eq.1) then
           ! ***process 'NWT0' finite dilution fluxes
           do ig=igfirs(1),igmax
             vector(ig)=real(flux(ig,1,nz-iz+1)-1.0d0)
           enddo
           call xsmput(draglib,'NWT0',vector(1:iglast(1)-1))
         endif
         !
         call xsmsix(draglib,' ',2)
       enddo
     enddo
     !
     ! ***create H-FACTOR if lkerma=.false.
     if(.not.lkerma) then
       ! use the legacy Dragon formula based on Q and pseudo-Q values
       call draq(nendf,ngen,matno,ytemp,nz0,ng,igfirs(1),iglast(1))
     endif
     !
     ! ***process fission x-sections.
     call drafis(ngen,matno,ytemp,nz0,igfirs(1),iglast(1),nbesp,iesp)
     !
     ! ***process scattering x-sections.
     call drasc(ngen,matno,ytemp,nz0,igrest,igecco,igfirs,iglast,ipflag)
     !
     call xsmsix(draglib,' ',2)
   enddo
   !
   ! scratch storage deallocation
   deallocate(rm2,rv,flux,rv2,deltau,scr)
   return
   end subroutine dramat
   !
   subroutine drauto(npen,matno,ng,ener,eaut0,eaut1,igaut0,igaut1,deli,irflag, &
   & urlimit,iprint)
   !-----------------------------------------------------------------
   !  put autolib data in draglib file for material matno
   !-----------------------------------------------------------------
   use mainio ! provides nsysi,contio,nsyso,nsyse
   use endf   ! provides endf routines and variables
   use util   ! provides timer,openz,repoz,error
   ! externals
   integer :: npen,matno,ng,igaut0,igaut1,irflag,iprint
   real(kr) ener(ng+1),eaut0,eaut1,deli,urlimit
   ! internals
   integer :: maxa,maxgr,maxtmp
   parameter (maxa=2000,maxgr=2000,maxtmp=100)
   integer :: igar,i,ibin,ibin0,idis,ig,ig1,ig2,ilfiss,isbmat,itm,ityxsm,nb, &
   & nbfine,nbin,lssf,nbdil,ntmp,nw,loc,ndata,nunr,nbinpt,irflag_2
   real oldtmp(maxtmp),gar1t(maxgr),gar1s(maxgr),gar1f(maxgr),deli3(1)
   integer nfs(maxgr),ijj(maxgr),njj(maxgr)
   real(kr) deltau,gar,em,ep,emlog,aa(6),enext,eplog,err,err1,err2,err3,err4, &
   & errmax,fact1,fact2,fact3,sum1,sum2,sum3,sum4,temps,tm,tp,uuu,zbin, &
   & deli2
   real(kr), allocatable, dimension(:) :: bener,bsig1,bsig2,bsig3,scr1,scr2a, &
   & scr2b,scr,eneurr,erandom
   real, allocatable, dimension(:) :: bgar,gar2
   logical lfind,lfiss
   character cd*4
   !
   allocate(scr(maxa))
   !**find the autolib energy limits.
   nbfine=0
   errmax=0.0
   do ig=igaut0,igaut1
     deli2=deli
     if(ener(ig+1).ge.2.0d2) then
       deli2=deli/4.0d0 ! use thinner bins above 0.1 keV
     endif
     zbin=log(ener(ig)/ener(ig+1))
     nbin=int(zbin/deli2+0.01)
     err1=abs(real(nbin)*deli2-zbin)/abs(zbin)
     err=err1
     if(ener(ig+1).lt.2.0d2) then
       err2=abs(real(nbin)*deli2+deli2/4.0d0-zbin)/abs(zbin)
       err3=abs(real(nbin)*deli2+deli2/2.0d0-zbin)/abs(zbin)
       err4=abs(real(nbin)*deli2+3.0d0*deli2/4.0d0-zbin)/abs(zbin)
       err=min(err1,err2,err3,err4)
     endif
     if(err.ne.err1) nbin=nbin+1
     nbfine=nbfine+nbin
     errmax=max(errmax,err)
   enddo
   write(nsyso,'( '' number of bin energy groups in autolib...... '',i10/ &
   & '' percent accuracy of the bin mesh............ '',1p,d9.2, &
   & ''%'')') nbfine,errmax*100.0
   allocate(bener(nbfine+1))
   bener(1)=ener(igaut0)
   nfs(:ng)=0
   uuu=log(ener(1)/ener(igaut0))
   deli3(1)=real(deli)
   ibin=1
   do ig=igaut0,igaut1
     deli2=deli
     if(ener(ig+1).ge.2.0d2) then
       deli2=deli/4.0d0 ! use thinner bins above 0.1 keV
       deli3(1)=min(real(deli/4.0d0),deli3(1))
     endif
     uuu=log(ener(1)/ener(ig))
     zbin=log(ener(ig)/ener(ig+1))
     nbin=int(zbin/deli2+0.01)
     do i=1,nbin
       ibin=ibin+1
       uuu=uuu+deli2
       bener(ibin)=ener(1)*exp(-uuu)
     enddo
     err1=abs(real(nbin)*deli2-zbin)/abs(zbin)
     err2=1.0e10
     err3=1.0e10
     err4=1.0e10
     err=err1
     if(ener(ig+1).lt.2.0d2) then
       err2=abs(real(nbin)*deli2+deli2/4.0d0-zbin)/abs(zbin)
       err3=abs(real(nbin)*deli2+deli2/2.0d0-zbin)/abs(zbin)
       err4=abs(real(nbin)*deli2+3.0d0*deli2/4.0d0-zbin)/abs(zbin)
       err=min(err1,err2,err3,err4)
     endif
     if(err.eq.err2) then
       nbin=nbin+1
       ibin=ibin+1
       uuu=uuu+deli2/4.0d0
       bener(ibin)=ener(ig+1)
       deli3(1)=min(real(deli2/4.0d0),deli3(1))
     else if(err.eq.err3) then
       nbin=nbin+1
       ibin=ibin+1
       uuu=uuu+deli2/2.0d0
       bener(ibin)=ener(ig+1)
       deli3(1)=min(real(deli2/4.0d0),deli3(1))
     else if(err.eq.err4) then
       nbin=nbin+1
       ibin=ibin+1
       uuu=uuu+3.0d0*deli2/4.0d0
       bener(ibin)=ener(ig+1)
       deli3(1)=min(real(deli2/4.0d0),deli3(1))
     endif
     if(ibin.gt.nbfine+1) call error('drauto','bin overflow',' ')
     nfs(ig)=nbin
   enddo
   allocate(bsig1(nbfine),bsig2(nbfine),bsig3(nbfine),bgar(nbfine+1))
   do ibin=1,nbfine+1
     bgar(ibin)=real(bener(ibin))
   enddo
   write(nsyso,'( &
   & '' adjusted elementary lethargy width .........'',1p,d11.4)') &
   & deli3(1)
   call xsmput(draglib,'BIN-DELI',deli3)
   call xsmput(draglib,'BIN-NFS',nfs(1:ng))
   call xsmput(draglib,'BIN-ENERGY',bgar(1:nbfine+1))
   !
   !**recover the temperature values.
   call repoz(npen)
   10 lfind=.false.
   do while(.not.lfind)
     if(npen.lt.0) then
       read(-npen,end=30) math,mfh,mth,nb,nw,aa
     else if(npen.gt.0) then
       read(npen,'(6e11.0,i4,i2,i3,i5)',end=30) aa,math,mfh,mth,nsp
     endif
     lfind=(math.eq.matno).and.(mfh.eq.1).and.(mth.eq.451)
   enddo
   if(.not.lfind) go to 30
   call contio(npen,0,0,scr,nb,nw)
   call contio(npen,0,0,scr,nb,nw)
   call contio(npen,0,0,scr,nb,nw)
   temps=c1h
   if(iprint.gt.0) then
     write(nsyso,'(/29h drauto: process temperature=,1p,e12.4)') temps
   endif
   call xsmlen(draglib,'TEMPERATURE',ntmp,ityxsm)
   if(ntmp.eq.0) call error('drauto','missing temperatures on xsm',' ')
   if(ntmp.gt.maxtmp) call error('drauto','oldtmp overflow',' ')
   call xsmget(draglib,'TEMPERATURE',oldtmp)
   isbmat=0
   do itm=1,ntmp
     if(abs(temps-oldtmp(itm)).le.1.0e-3*temps) then
        isbmat=itm
        go to 20
     endif
   enddo
   20 if(isbmat.eq.0) call error('drauto','missing temperature',' ')
   write (cd,'(i4.4)') isbmat
   call xsmsix(draglib,'SUBTMP'//cd,1)
   allocate(scr1(maxa))
   !
   !**recover the purr data (probability tables in the UR domain).
   lssf=0
   irflag_2=irflag
   if(irflag.eq.1) then
     math=matno
     call findf(math,2,152,npen)
     if(math.eq.-1) then
       irflag_2=0
       go to 40
     endif
     allocate(erandom(nbfine))
     kk=-matno ! use a different random sequence for each isotope
     do i=1,nbfine
       erandom(i)=draran(kk)
     enddo
     call contio(npen,0,0,scr1,nb,nw)
     lssf=l1h
     call listio(npen,0,0,scr1,nb,nw)
     if(nw.gt.maxa) call error('drauto','endf input size exceeded(1)',' ')
     nbdil=l2h
     ndata=n1h
     nunr=n2h
     allocate(scr2a(ndata))
     scr2a(:nw-6)=scr1(7:nw)
     loc=1+nw-6
     do while (nb.ne.0)
       call moreio(npen,0,0,scr2a(loc),nb,nw)
       loc=loc+nw
       if(loc-1.gt.n1h) call error('drauto','endf input size exceeded(2)',' ')
     enddo
     !
     call findf(matno,2,153,npen)
     call contio(npen,0,0,scr1,nb,nw)
     nbinpt=n2h
     call listio(npen,0,0,scr1,nb,nw)
     if(nw.gt.maxa) call error('drauto','endf input size exceeded(3)',' ')
     ndata=n1h
     nunr=n2h
     allocate(scr2b(ndata))
     scr2b(:nw-6)=scr1(7:nw)
     loc=1+nw-6
     do while (nb.ne.0)
       call moreio(npen,0,0,scr2b(loc),nb,nw)
       loc=loc+nw
       if(loc-1.gt.n1h) call error('drauto','endf input size exceeded(4)',' ')
     enddo
     allocate(eneurr(nunr))
     loc=1
     do i=1,nunr
       eneurr(i)=scr2b(loc)
       loc=loc+1+6*nbinpt
     enddo
     write(nsyso,'('' unresolved energy limits (eV) ..............'')')
     write(nsyso,'(5x,1p,10e12.4)') eneurr(:)
   endif
   !
   !**recover the pointwise total cross sections.
   40 call findf(matno,3,1,npen)
   call contio(npen,0,0,scr1,nb,nw)
   ep=0.0d0
   call gety1(ep,enext,idis,tp,npen,scr1)
   em=eaut0
   call gety1(em,enext,idis,tm,npen,scr1)
   ibin=nbfine
   gar=0.0d0
   emlog=log(eaut1/em)
   do while(ep*(1.0d0+1.0d-6).lt.eaut1)
     if(ibin.le.0) call error('drauto','invalid index',' ')
     ep=min(enext,bener(ibin))
     eplog=log(eaut1/ep)
     call gety1(ep,enext,idis,tp,npen,scr1)
     if((irflag_2.eq.0).or.(ep.le.urlimit)) then
       gar=gar+0.5d0*(tm+tp)*(emlog-eplog)
     else
       ! select a random value in the probability table
       gar=gar+draurr(1,ibin,nbinpt,nunr,nbdil,lssf,ep,erandom(ibin),eneurr, &
       & scr2a,scr2b)*(emlog-eplog)
     endif
     if(ep.eq.bener(ibin)) then
       deltau=log(bener(ibin)/bener(ibin+1))
       bsig1(ibin)=gar/deltau
       ibin=ibin-1
       gar=0.0d0
     endif
     em=ep
     tm=tp
     emlog=eplog
   enddo
   !
   !**recover the pointwise scattering cross sections.
   call findf(matno,3,2,npen)
   call contio(npen,0,0,scr1,nb,nw)
   ep=0.0d0
   call gety1(ep,enext,idis,tp,npen,scr1)
   em=eaut0
   call gety1(em,enext,idis,tm,npen,scr1)
   ibin=nbfine
   gar=0.0d0
   emlog=log(eaut1/em)
   do while(ep*(1.0d0+1.0d-6).lt.eaut1)
     if(ibin.le.0) call error('drauto','invalid index',' ')
     ep=min(enext,bener(ibin))
     eplog=log(eaut1/ep)
     call gety1(ep,enext,idis,tp,npen,scr1)
     if((irflag_2.eq.0).or.(ep.le.urlimit)) then
       gar=gar+0.5d0*(tm+tp)*(emlog-eplog)
     else
       ! select a random value in the probability table
       gar=gar+draurr(2,ibin,nbinpt,nunr,nbdil,lssf,ep,erandom(ibin),eneurr, &
       & scr2a,scr2b)*(emlog-eplog)
     endif
     if(ep.eq.bener(ibin)) then
       deltau=log(bener(ibin)/bener(ibin+1))
       bsig2(ibin)=gar/deltau
       ibin=ibin-1
       gar=0.0d0
     endif
     em=ep
     tm=tp
     emlog=eplog
   enddo
   !
   !**recover the pointwise fission cross sections.
   call xsmlen(draglib,'NUSIGF',ilfiss,ityxsm)
   if(ilfiss.gt.0) then
     call findf(matno,3,18,npen)
     call contio(npen,0,0,scr1,nb,nw)
     ep=0.0d0
     call gety1(ep,enext,idis,tp,npen,scr1)
     em=eaut0
     call gety1(em,enext,idis,tm,npen,scr1)
     ibin=nbfine
     gar=0.0d0
     emlog=log(eaut1/em)
     do while(ep*(1.0d0+1.0d-6).lt.eaut1)
       if(ibin.le.0) call error('drauto','invalid index',' ')
       ep=min(enext,bener(ibin))
       eplog=log(eaut1/ep)
       call gety1(ep,enext,idis,tp,npen,scr1)
       if((irflag_2.eq.0).or.(ep.le.urlimit)) then
         gar=gar+0.5d0*(tm+tp)*(emlog-eplog)
       else
         ! select a random value in the probability table
         gar=gar+draurr(3,ibin,nbinpt,nunr,nbdil,lssf,ep,erandom(ibin), &
         & eneurr,scr2a,scr2b)*(emlog-eplog)
       endif
       if(ep.eq.bener(ibin)) then
         deltau=log(bener(ibin)/bener(ibin+1))
         bsig3(ibin)=gar/deltau
         ibin=ibin-1
         gar=0.0d0
       endif
       em=ep
       tm=tp
       emlog=eplog
     enddo
   endif
   if(irflag_2.eq.1) deallocate(scr2b,scr2a,eneurr,erandom)
   !
   !**normalization of BIN cross sections.
   allocate(gar2(ng**2))
   gar1t(:ng)=0.0
   gar1s(:ng)=0.0
   gar1f(:ng)=0.0
   gar2(:ng*ng)=0.0
   call xsmget(draglib,'NTOT0',gar1t)
   call xsmget(draglib,'IJJS00',ijj)
   call xsmget(draglib,'NJJS00',njj)
   call xsmget(draglib,'SCAT00',gar2)
   if(ilfiss.gt.0) call xsmget(draglib,'NUSIGF',gar1f)
   igar=0
   do ig2=1,ng
     do ig1=ijj(ig2),ijj(ig2)-njj(ig2)+1,-1
       igar=igar+1
       gar1s(ig1)=gar1s(ig1)+gar2(igar)
     enddo
   enddo
   deallocate(gar2)
   ibin0=0
   do ig=igaut0,igaut1
     sum1=0.0
     sum2=0.0
     sum3=0.0
     sum4=0.0
     do ibin=ibin0+1,ibin0+nfs(ig)
       deltau=log(bener(ibin)/bener(ibin+1))
       sum1=sum1+deltau
       sum2=sum2+bsig1(ibin)*deltau
       sum3=sum3+bsig2(ibin)*deltau
       if(ilfiss.gt.0) sum4=sum4+bsig3(ibin)*deltau
     enddo
     fact1=gar1t(ig)*(sum1/sum2)
     fact2=gar1s(ig)*(sum1/sum3)
     fact3=0.0
     if((ilfiss.gt.0).and.(sum4.gt.0.0)) fact3=gar1f(ig)*(sum1/sum4)
     if(iprint.gt.0) then
       write(nsyso,'('' Autolib normalization in group'',i5,'':'',3f7.4)') &
       & ig,fact1,fact2,fact3
     endif
     do ibin=ibin0+1,ibin0+nfs(ig)
       bsig1(ibin)=fact1*bsig1(ibin)
       bsig2(ibin)=fact2*bsig2(ibin)
       if(ilfiss.gt.0) bsig3(ibin)=fact3*bsig3(ibin)
     enddo
     ibin0=ibin0+nfs(ig)
   enddo
   !
   !**store the BIN cross sections on xsm.
   do ibin=1,nbfine
     bgar(ibin)=real(bsig1(ibin))
   enddo
   call xsmput(draglib,'BIN-NTOT0',bgar(1:nbfine))
   if((isbmat.eq.1).and.(iprint.gt.1)) then
     write(nsyso,'('' -----------------------------------'')')
     write(nsyso,100) (ener(i),i=igaut0,igaut1+1)
     write(nsyso,'('' --BIN-ENERGY-----------------------'')')
     write(nsyso,100) (bener(i),i=1,nbfine+1)
     write(nsyso,'('' --BIN-NTOT0------------------------'')')
     write(nsyso,100) (bsig1(i),i=1,nbfine),bsig1(nbfine)
   endif
   do ibin=1,nbfine
     bgar(ibin)=real(bsig2(ibin))
   enddo
   call xsmput(draglib,'BIN-SIGS00',bgar(1:nbfine))
   if(ilfiss.gt.0) then
     lfiss=.false.
     do ibin=1,nbfine
       lfiss=lfiss.or.(bsig3(ibin).ne.0.0)
       bgar(ibin)=real(bsig3(ibin))
     enddo
     if(lfiss) call xsmput(draglib,'BIN-NUSIGF',bgar(1:nbfine))
   endif
   call xsmsix(draglib,' ',2)
   !
   call tomend(npen,0,0,scr)
   deallocate(scr1)
   go to 10
   30 deallocate(bgar,bsig3,bsig2,bsig1,bener)
   deallocate(scr)
   return
   100 format(10(1p,e12.5,1h,))
   end subroutine drauto
   !
   subroutine draq(nendf,ngen,matno,ytemp,nz0,ng,igfirs,iglast)
   !-----------------------------------------------------------------
   !   recompute kerma factors using the legacy Dragon method based
   !   on Q and pseudo-Q information.
   !-----------------------------------------------------------------
   use endf   ! provides endf routines and variables
   use util   ! provides error
   integer, intent(in) :: nendf,ngen,matno,nz0,ng,igfirs,iglast
   real(kr), intent(in) :: ytemp
   integer :: maxa,maxgr,maxnl,maxnz,maxedi,nb,nw,ngtmp,nl,nz,lz
   parameter (maxa=3000,maxgr=2000,maxnl=8,maxnz=30,maxedi=12,lz=6)
   integer ied,ig,igmax,iz,nztmp
   logical lfind,exist,lfiss
   real(kr) kap
   character cd*4
   integer,save,dimension(maxedi) :: malist1= &
   & (/ 16, 17, 18, 28, 37,102,103,104,105,106,107,108 /)
   integer,save,dimension(maxedi) :: malist2= &
   & (/-99,-99,-99,-99,-99,-99,600,650,700,750,800,-99 /)
   real,allocatable,dimension(:) :: vector
   real(kr),allocatable,dimension(:) :: scr
   real(kr),allocatable,dimension(:,:) :: sigq
   real(kr),allocatable,dimension(:,:,:) :: flux,rv
   !
   ! scratch storage allocation
   allocate(scr(maxa))
   allocate(sigq(maxgr,maxnz),flux(maxgr,maxnl,maxnz),rv(maxgr,maxnl,maxnz))
   !
   ! initialize the sigq array.
   sigq(:maxgr,:maxnz)=0.0
   !
   call repoz(nendf)
   call findf(matno,1,451,nendf)
   call contio(nendf,0,0,scr,nb,nw)
   call contio(nendf,0,0,scr,nb,nw)
   nz=0
   do ied=1,maxedi
     !  ***process x-sections for reaction malist1(ied).
     call fvect(maxgr,maxnl,maxnz,ngen,matno,malist1(ied),ytemp,ngtmp, &
     & nl,nztmp,rv,flux,exist)
     if(.not.exist) cycle
     if(ngtmp.ne.ng) call error('draq','inconsistent ng.',' ')
     if(nztmp.gt.maxnz) call error('draq','maxnz overflow.',' ')
     if(malist1(ied).eq.18) lfiss=.true.
     nz=max(nz,nztmp)
     do iz=1,nztmp
       do ig=igfirs,iglast-1
         if(flux(ig,1,1).eq.0.0) cycle
         flux(ig,1,nztmp-iz+1)=flux(ig,1,nztmp-iz+1)/flux(ig,1,1)
       enddo
     enddo
     !
     call repoz(nendf)
     lfind=.false.
     do while (.not.lfind)
       if(nendf.lt.0) then
         read(-nendf,end=10) math,mfh,mth,nb,nw
       else if(nendf.gt.0) then
         read(nendf,'(66x,i4,i2,i3,i5)',end=10) math,mfh,mth,nsp
       endif
       lfind=(math.eq.matno).and.(mfh.eq.3).and.((mth.eq.malist1(ied)).or. &
       & (mth.eq.malist2(ied)))
     enddo
     if(.not.lfind) cycle
     ! recover the Q or pseudo-Q of the reaction (may be negative)
     call contio(nendf,0,0,scr,nb,nw)
     kap=c1h
     if(mth.eq.18) then
       ! ***fission -- recover pseudo-Q kappa info in mf=1, mt=458 if exists
       call repoz(nendf)
       lfind=.false.
       do while (.not.lfind)
         if(nendf.lt.0) then
           read(-nendf,end=30) math,mfh,mth,nb,nw
         else if(nendf.gt.0) then
           read(nendf,'(66x,i4,i2,i3,i5)',end=30) math,mfh,mth,nsp
         endif
         lfind=(math.eq.matno).and.(mfh.eq.1).and.(mth.eq.458)
       enddo
       if(lfind) then
         call listio(nendf,0,0,scr,nb,nw)
         kap=scr(lz+15)
       endif
     endif
     if(kap.eq.0.0) cycle
     do ig=1,ng
       sigq(ig,1)=sigq(ig,1)+kap*rv(ig,1,1)
     enddo
     do iz=2,nztmp
       ! process finite-dilution h-factors
       do ig=igfirs,iglast-1
         if(rv(ig,1,1).eq.0.0) cycle
         sigq(ig,iz)=sigq(ig,iz)+kap*(rv(ig,1,iz)*flux(ig,1,iz)-rv(ig,1,1))
       enddo
     enddo
     10 continue
   enddo ! ied
   !
   igmax=ng
   do ig=ng,1,-1
     igmax=ig
     if(sigq(ig,1).ne.0.0) go to 20
   enddo
   go to 30
   20 allocate(vector(ng))
   do ig=1,igmax
     vector(ig)=real(sigq(ig,1))
   enddo
   call xsmput(draglib,'H-FACTOR',vector(1:igmax))
   do iz=1,nz-1
     write (cd,'(i4.4)') nz0+iz
     call xsmsix(draglib,'SUBMAT'//cd,1)
     do ig=1,iglast-1
       vector(ig)=real(sigq(ig,nz-iz+1))
     enddo
     call xsmput(draglib,'H-FACTOR',vector(1:iglast-1))
     call xsmsix(draglib,' ',2)
   enddo ! iz
   deallocate(vector)
   30 deallocate(rv,flux,sigq)
   deallocate(scr)
   return
   end subroutine draq
   !
   subroutine drafis(nin,matd,ytemp,nz0,igfirs,iglast,nbesp,iesp)
   !-----------------------------------------------------------------
   !   utility routine for recovering fission-related information from a
   !   gendf file.
   !   input parameters:
   !   nin     file unit number of the gendf file.
   !   matd    material number.
   !   ytemp   absolute temperature (Kelvin).
   !   nz0     previous dilution index.
   !   igfirs  fastest group index with self-shielding effect.
   !   iglast  fastest thermal group index with no self-shielding effect.
   !   nbesp   number of energy-dependent fission spectra.
   !   iesp    energy limit indices of the energy-dependent fission spectra.
   !-----------------------------------------------------------------
   use mainio ! provides nsysi,contio,nsyso,nsyse
   use endf   ! provides endf routines and variables
   use util   ! provides timer,openz,repoz,error
   integer :: maxa,maxgr,maxnl,maxnz,maxdel,lz
   real(kr) :: ytemp
   parameter (maxa=20000,maxgr=2000,maxnl=8,maxnz=30,maxdel=10,lz=6)
   integer nin,matd,nz0,igfirs,iglast,nbesp,iesp(nbesp+1)
   logical lfind,exist,exist2
   character text2*2,cd*4
   real(kr) aa(6),ssum(maxdel),vsum(maxgr)
   integer idel,ig,ig1,ig2,igmax,il,isp,iz,jg,loc,nb,ndel,nf,ng,ng2, &
   & nl,nw,nz,nz2,nzm,nzv,alloc_ok
   real(kr) all
   real gar1(maxdel),gar2(maxgr)
   integer, save, dimension(4) :: mflist= (/ 19,20,21,38 /)
   real(kr),allocatable,dimension(:) :: scr
   real(kr),allocatable,dimension(:,:) :: chidel,sigf2
   real(kr),allocatable,dimension(:,:,:) :: rv,sigf,flux
   real(kr),allocatable,dimension(:,:,:,:) :: sigfm,sigfm2
   allocate(scr(maxa))
   allocate(rv(maxgr,maxnl,maxnz),sigf(maxgr,maxnl,maxnz), &
   & flux(maxgr,maxnl,maxnz),chidel(maxgr,maxdel),sigf2(maxgr,maxnz))
   allocate(sigfm(maxgr,maxgr,maxnl,maxnz),stat=alloc_ok)
   if(alloc_ok /= 0) call error('drafis','unable to allocate sigfm',' ')
   allocate(sigfm2(maxgr,maxgr,maxnl,maxnz),stat=alloc_ok)
   if(alloc_ok /= 0) call error('drafis','unable to allocate sigfm2',' ')
   !
   ! initialize the sigf and sigfm arrays.
   sigf(:maxgr,:maxnl,:maxnz)=0.0
   sigfm(:maxgr,:maxgr,:maxnl,:maxnz)=0.0
   sigfm2(:maxgr,:maxgr,:maxnl,:maxnz)=0.0
   !
   ! recover the 'NFTOT' vector fission cross sections.
   call fvect(maxgr,maxnl,maxnz,nin,matd,18,ytemp,ng,nl,nzv,rv,flux,exist)
   if(.not.exist) go to 60
   !
   ! recover the 'NFTOT' matrix fission cross sections.
   call fmatr(maxgr,maxnl,maxnz,nin,matd,18,ytemp,ng,nl,nzm,sigfm,flux,exist)
   if(.not.exist) then
     ! try to recover a fission matrix by summing partial fission matrices
     ! (if they exist).
     nzm=0
     do nf=1,4
       call fmatr(maxgr,maxnl,maxnz,nin,matd,mflist(nf),ytemp,ng,nl, &
       & nz2,sigfm2,flux,exist2)
       exist=exist.or.exist2
       if(exist2) then
         nzm=max(nzm,nz2)
         do iz=1,nz2
           do il=1,nl
             do ig=1,ng
               do jg=1,ng
                 sigfm(jg,ig,il,iz)=sigfm(jg,ig,il,iz)+sigfm2(jg,ig,il,iz)
               enddo
             enddo
           enddo
         enddo
       endif
     enddo
   endif
   if(.not.exist) then
     write(nsyso,'('' missing fission matrix  '', &
     & ''(endf identification ='',i8,'')'')') matd
     return
   endif
   !
   ! recover the delayed spectra.
   chidel(:maxgr,:maxdel)=0.0
   call repoz(nin)
   10 lfind=.false.
   do while (.not.lfind)
     if(nin.lt.0) then
       read(-nin,end=20) math,mfh,mth,nb,nw
     else if(nin.gt.0) then
       read(nin,'(6e11.0,i4,i2,i3,i5)',end=20) aa,math,mfh,mth,nsp
     endif
     lfind=(math.eq.matd).and.(mfh.eq.5).and.(mth.eq.455)
   enddo
   if(.not.lfind) go to 20
   if(nin.lt.0) then
     call skiprz(nin,-1)
     read(-nin) math,mfh,mth,nb,nw,aa
   endif
   ndel=nint(aa(3))
   nz=nint(aa(4))
   ng=nint(aa(6))
   call listio(nin,0,0,scr(1),nb,nw)
   if(abs(scr(1)-ytemp).gt.1.0e-3) then
     call tomend(nin,0,0,scr)
     go to 10
   endif
   if(ng.gt.maxgr) call error('drafis','maxgr overflow',' ')
   if(nz.ne.1) call error('drafis','no dilution dependence expected',' ')
   if(ndel.gt.maxdel) call error('drafis','maxdel overflow',' ')
   ng2=l1h
   loc=1+nw
   do while (nb.ne.0)
     if(loc+302.gt.maxa) call error('drafis','endf input size exceeded',' ')
     call moreio(nin,0,0,scr(loc),nb,nw)
     loc=loc+nw
   enddo
   do idel=1,ndel
     gar1(idel)=real(scr(lz+idel))
   enddo
   if(ndel.gt.0) call xsmput(draglib,'LAMBDA-D',gar1(1:ndel))
   do idel=1,ndel
     ssum(idel)=0.0
     do ig=1,ng2-1
       ssum(idel)=ssum(idel)+scr(lz+ndel+(ig-1)*ndel+idel)
     enddo
     do ig=1,ng2-1
       chidel(ng-ig+1,idel)=scr(lz+ndel+(ig-1)*ndel+idel)/ssum(idel)
     enddo
     do ig=ng2,ng
       chidel(ng-ig+1,idel)=0.0
     enddo
     do ig=1,ng
       gar2(ig)=real(chidel(ig,idel))
     enddo
     write(text2,'(i2.2)') idel
     call xsmput(draglib,'CHI'//text2,gar2(1:ng))
   enddo
   call fvect(maxgr,maxnl,maxnz,nin,matd,455,ytemp,ng,nl,nz,rv,flux,exist)
   call fvect(maxgr,maxnl,maxnz,nin,matd,18,ytemp,ng,nl,nz,sigf,flux,exist)
   do iz=1,nz
     do ig=igfirs,iglast-1
       if(flux(ig,1,1).eq.0.0) cycle
       flux(ig,1,nz-iz+1)=flux(ig,1,nz-iz+1)/flux(ig,1,1)
     enddo
   enddo
   do idel=1,ndel
     do ig=1,ng
       gar2(ig)=real(rv(ig,1,1)*sigf(ig,1,1)*ssum(idel))
     enddo
     igmax=ng
     do ig=ng,1,-1
       igmax=ig
       if(gar2(ig).ne.0.0) go to 15
     enddo
     15 write(text2,'(i2.2)') idel
     call xsmput(draglib,'NUSIGF'//text2,gar2(1:igmax))
     do iz=1,nz-1
       !
       ! ***process finite dilution x-sections
       gar2(1:iglast-1)=0.0
       do ig=igfirs,iglast-1
         if(sigf(ig,1,1).eq.0.0) cycle
         gar2(ig)=real(rv(ig,1,1)*ssum(idel)*(sigf(ig,1,nz-iz+1)* &
         & flux(ig,1,nz-iz+1)-sigf(ig,1,1)))
       enddo
       igmax=iglast-1
       do ig=iglast-1,1,-1
         igmax=ig
         if(gar2(ig).ne.0.0) go to 16
       enddo
       go to 17
       16 write (cd,'(i4.4)') nz0+iz
       call xsmsix(draglib,'SUBMAT'//cd,1)
       call xsmput(draglib,'NUSIGF'//text2,gar2(1:igmax))
       call xsmsix(draglib,' ',2)
       17 continue
     enddo
     do iz=1,nz
       do ig2=1,ng
         do ig1=1,ng
           sigfm(ig1,ig2,1,iz)=sigfm(ig1,ig2,1,iz)+rv(ig1,1,1)* &
           & sigf(ig1,1,iz)*ssum(idel)*chidel(ig2,idel)
         enddo
       enddo
     enddo
   enddo
   !
   ! process the prompt 'NUSIGF' and 'CHI' vectors.
   20 call fvect(maxgr,maxnl,maxnz,nin,matd,18,ytemp,ng,nl,nz,sigf,flux,exist)
   if((nzm.ne.nz).and.(nzm.ne.1)) then
     call error('drafis','invalid number of dilutions with mf=6, mt='//'18',' ')
   endif
   do 25 isp=1,nbesp
   do iz=1,nzm
     if(iz.eq.1) then
       all=0.0
       do ig2=1,ng
         vsum(ig2)=0.0
         do ig1=iesp(isp)+1,iesp(isp+1)
           vsum(ig2)=vsum(ig2)+sigfm(ig1,ig2,1,iz)*flux(ig1,1,iz)
         enddo
         do ig1=iesp(isp)+1,iesp(isp+1)
           all=all+sigfm(ig1,ig2,1,iz)*flux(ig1,1,iz)
         enddo
       enddo
       do ig1=1,ng
         gar2(ig1)=real(vsum(ig1)/all)
       enddo
       if(nbesp.eq.1) then
         call xsmput(draglib,'CHI',gar2(1:ng))
       else if(all.ne.0.0) then
         write(text2,'(i2.2)') isp
         call xsmput(draglib,'CHI--'//text2,gar2(1:ng))
       endif
     endif
   enddo
   25 continue
   do iz=1,nzm
     do ig1=1,ng
       vsum(ig1)=0.0
       do ig2=1,ng
         vsum(ig1)=vsum(ig1)+sigfm(ig1,ig2,1,iz)
       enddo
     enddo
     if(nzm.gt.1) then
       do ig=1,ng
         sigf2(ig,iz)=vsum(ig)
       enddo
     endif
   enddo
   if(nzm.eq.1) then
     do iz=nz,1,-1
       do ig=1,ng
         if(sigf(ig,1,1).ne.0.0) then
           sigf2(ig,iz)=sigf(ig,1,iz)*vsum(ig)/sigf(ig,1,1)
         endif
       enddo
     enddo
   endif
   igmax=ng
   do ig=ng,1,-1
     igmax=ig
     if(sigf2(ig,1).ne.0.0) go to 30
   enddo
   30 do ig=1,igmax
     gar2(ig)=real(sigf2(ig,1))
   enddo
   call xsmput(draglib,'NUSIGF',gar2(1:igmax))
   do iz=1,nz-1
     !
     ! ***process 'NUSIGF' finite dilution x-sections
     do ig=1,iglast-1
       gar2(ig)=0.0d0
     enddo
     do ig=igfirs,iglast-1
       if(flux(ig,1,1).eq.0.0) cycle
       if(sigf2(ig,1).eq.0.0) cycle
       flux(ig,1,nz-iz+1)=flux(ig,1,nz-iz+1)/flux(ig,1,1)
       gar2(ig)=real(sigf2(ig,nz-iz+1)*flux(ig,1,nz-iz+1)-sigf2(ig,1))
     enddo
     igmax=iglast-1
     do ig=iglast-1,1,-1
       igmax=ig
       if(gar2(ig).ne.0.0) go to 40
     enddo
     go to 50
     40 write (cd,'(i4.4)') nz0+iz
     call xsmsix(draglib,'SUBMAT'//cd,1)
     call xsmput(draglib,'NUSIGF',gar2(1:igmax))
     call xsmsix(draglib,' ',2)
     50 continue
   enddo
   !
   ! scratch storage deallocation
   60 deallocate(sigfm2,sigfm)
   deallocate(sigf2,chidel,flux,sigf,rv)
   deallocate(scr)
   return
   end subroutine drafis
   !
   subroutine drasc(nin,matd,ytemp,nz0,igrest,igecco,igfirs,iglast,ipflag)
   !-----------------------------------------------------------------
   !   utility routine for recovering scattering-related information from
   !   a gendf file.
   !   input parameters:
   !   nin     file unit number of the gendf file.
   !   matd    material number.
   !   ytemp   absolute temperature (Kelvin).
   !   nz0     previous dilution index.
   !   igrest  number of thermal groups.
   !   igecco  last eccolib group with equal lethargy width.
   !   igfirs  fastest group index with self-shielding effect.
   !   iglast  fastest thermal group index with no self-shielding effect.
   !   ipflag  fission product flag.
   !-----------------------------------------------------------------
   use mainio ! provides nsysi,contio,nsyso,nsyse
   use endf   ! provides endf routines and variables
   use util   ! provides timer,openz,repoz,error
   integer :: maxa,maxgr,maxnl,maxnz,maxmat,maxgar
   parameter (maxa=2000)
   parameter(maxgr=2000,maxnl=8,maxnz=30,maxmat=47,maxgar=maxgr*maxgr)
   integer nin,matd,nz0,igrest,igecco,ijj(maxgr),njj(maxgr),igfirs(maxnl), &
   & iglast(maxnl),ipflag
   logical exist,exist2,exist3,lfind
   character cd*4,cdl*2
   integer i,ig,ig1,ig2,igar,igmax,igmin,il,imat,imax,imt,iz,loc,nb,ng, &
   & ngther,nl,nl2,nlgar,nw,nz,nze,alloc_ok
   real(kr) ytemp,aa(6),delta,dgar,garmax
   real,allocatable,dimension(:) :: gar
   real(kr),allocatable,dimension(:) :: scr
   real(kr),allocatable,dimension(:,:,:) :: flux
   real(kr),allocatable,dimension(:,:,:,:) :: rm,sigsm
   integer,save,dimension(maxmat) :: mtlist= &
   & (/ 2,5,16,17,28,37,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67, &
   & 68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91 /)
   !
   ! scratch storage allocation
   allocate(scr(maxa),gar(maxgar))
   allocate(rm(maxgr,maxgr,maxnl,maxnz),stat=alloc_ok)
   if(alloc_ok /= 0) call error('drasc','unable to allocate rm',' ')
   allocate(flux(maxgr,maxnl,maxnz),stat=alloc_ok)
   if(alloc_ok /= 0) call error('drasc','unable to allocate flux',' ')
   allocate(sigsm(maxgr,maxgr,maxnl,maxnz),stat=alloc_ok)
   if(alloc_ok /= 0) call error('drasc','unable to allocate sigsm',' ')
   !
   ! recover the number of thermal-corrected groups (ngther).
   ngther=0
   call repoz(nin)
   10 lfind=.false.
   call skiprz(nin,1)
   do while (.not.lfind)
     if(nin.lt.0) then
       read(-nin,end=20) math,mfh,mth,nb,nw
     else if(nin.gt.0) then
       read(nin,'(6e11.0,i4,i2,i3,i5)',end=20) aa,math,mfh,mth,nsp
     endif
     lfind=(math.eq.matd).and.(mfh.eq.6).and.(mth.ge.221).and.(mth.le.250)
   enddo
   if(.not.lfind) go to 20
   if(nin.lt.0) then
     call skiprz(nin,-1)
     read(-nin) math,mfh,mth,nb,nw,aa
   endif
   ng=nint(aa(6))
   call listio(nin,0,0,scr(1),nb,nw)
   if(abs(scr(1)-ytemp).gt.1.0e-3) then
     call tomend(nin,0,0,scr)
     go to 10
   endif
   lfind=.true.
   do while (mth.ne.0) ! loop over energy groups
     if(.not.lfind) call listio(nin,0,0,scr(1),nb,nw)
     if(n2h.eq.ng) go to 20
     ngther=max(ngther,n2h)
     lfind=.false.
     loc=1+nw
     do while (nb.ne.0)
       if(loc+302.gt.maxa) call error('drasc','endf input size exceeded',' ')
       call moreio(nin,0,0,scr(loc),nb,nw)
       loc=loc+nw
     enddo
   enddo
   20 ngther=min(ngther,igrest)
   write(nsyso,'(/ &
   & '' number of energy groups affected by thermalization for '', &
   & ''material'',i6,'' ='',i6)') matd,ngther
   ngther=ng-ngther
   !
   ! **recover thermal correction x-sections and compute the
   ! correction term for groups indices higher than ngther.
   do iz=1,maxnz
     do il=1,maxnl
       do ig2=1,maxgr
         do ig1=1,maxgr
          sigsm(ig1,ig2,il,iz)=0.0
         enddo
       enddo
     enddo
   enddo
   exist3=.false.
   do imt=250,221,-1
     call fmatr(maxgr,maxnl,1,nin,matd,imt,ytemp,ng,nl2,nz,rm,flux,exist2)
     if(exist2) then
       exist3=.true.
       nlgar=nl2
       do iz=1,maxnz
         do il=1,nl2
           do ig2=1,ng
             do ig1=ngther+1,ng
               sigsm(ig1,ig2,il,iz)=sigsm(ig1,ig2,il,iz)+rm(ig1,ig2,il,1)
             enddo
           enddo
         enddo
       enddo
     endif
   enddo
   if(exist3) then
     call fmatr(maxgr,maxnl,maxnz,nin,matd,2,ytemp,ng,nl,nz,rm,flux,exist)
     if(nl.ne.nlgar) call error('drasc','inconsistent number of L'// &
     & 'egendre orders(2)',' ')
     do iz=1,maxnz
       do il=1,nl
         do ig2=1,ng
           do ig1=ngther+1,ng
             sigsm(ig1,ig2,il,iz)=sigsm(ig1,ig2,il,iz)-rm(ig1,ig2,il,1)
           enddo
         enddo
       enddo
     enddo
   endif
   !
   nze=0
   do imat=1,maxmat
     call fmatr(maxgr,maxnl,maxnz,nin,matd,mtlist(imat),ytemp,ng,nl, &
     & nz,rm,flux,exist)
     if(mtlist(imat).eq.2) nze=nz
     if(exist) then
       if(nz.eq.1) then
         do iz=1,nze
           do il=1,nl
             do ig2=1,ng
               do ig1=1,ng
                 sigsm(ig1,ig2,il,iz)=sigsm(ig1,ig2,il,iz)+rm(ig1,ig2,il,1)
               enddo
             enddo
           enddo
         enddo
       else if(nz.eq.nze) then
         do iz=1,nz
           do il=1,nl
             do ig2=1,ng
               do ig1=1,ng
                 sigsm(ig1,ig2,il,iz)=sigsm(ig1,ig2,il,iz)+rm(ig1,ig2,il,iz)
               enddo
             enddo
           enddo
         enddo
       else
         call error('drasc','invalid number of dilutions.',' ')
       endif
     endif
   enddo
   !
   ! store the scattering matrices on the xsm file
   call fmatr(maxgr,maxnl,maxnz,nin,matd,2,ytemp,ng,nl,nz,rm,flux,exist)
   do il=1,nl
     !
     ! ***lump the scattering matrices
     if(ipflag.eq.1) then
       do iz=1,nz
         do ig1=1,ng
           dgar=0.0
           do ig2=1,ng
             dgar=dgar+sigsm(ig1,ig2,il,iz)
           enddo
           sigsm(ig1,:ng,il,iz)=0.0
           sigsm(ig1,ig1,il,iz)=dgar
         enddo
       enddo
     else if(ipflag.eq.2) then
       do iz=1,nz
         do ig1=1,igecco
           if(il.eq.1) then
             dgar=0.0
             do ig2=1,ng
               dgar=dgar+sigsm(ig1,ig2,1,iz)
             enddo
           else
             dgar=sigsm(ig1,ig1,1,iz)
           endif
           sigsm(ig1,:igecco,il,iz)=0.0
           sigsm(ig1,ig1,il,iz)=dgar
         enddo
       enddo
     endif
     !
     write(cdl,'(i2.2)') il-1
     igar=0
     garmax=0.0
     ! ***ig2 is the secondary group
     do ig2=1,ng
       igmin=ig2
       igmax=ig2
       !
       ! ***ig1 is the primary group
       do ig1=ng,1,-1
         if(sigsm(ig1,ig2,il,1).ne.0.0) then
           igmin=min(igmin,ig1)
           igmax=max(igmax,ig1)
         endif
       enddo
       ijj(ig2)=igmax
       njj(ig2)=igmax-igmin+1
       do ig1=igmax,igmin,-1
         igar=igar+1
         if(igar.gt.maxgar) call error('drasc','gar overflow.',' ')
         gar(igar)=real(sigsm(ig1,ig2,il,1))
         if(abs(gar(igar)).gt.garmax) garmax=abs(gar(igar))
       enddo
     enddo
     !
     ! ***write on the xsm file
     call xsmput(draglib,'SCAT'//cdl,gar(1:igar))
     call xsmput(draglib,'NJJS'//cdl,njj(1:ng))
     call xsmput(draglib,'IJJS'//cdl,ijj(1:ng))
     do iz=1,nz-1
       write (cd,'(i4.4)') nz0+iz
       call xsmsix(draglib,'SUBMAT'//cd,1)
       !
       ! ***process finite dilution scattering information
       do ig=1,ng
         if(flux(ig,il,1).eq.0.0) cycle
         flux(ig,il,nz-iz+1)=flux(ig,il,nz-iz+1)/flux(ig,il,1)
       enddo
       igar=0
       garmax=0.0
       ! ***ig2 is the secondary group
       do ig2=1,igfirs(il)-1
         ijj(ig2)=ig2
         njj(ig2)=1
         igar=igar+1
         if(igar.gt.maxgar) call error('drasc','gar overflow.',' ')
         gar(igar)=0.0d0
       enddo
       do ig2=igfirs(il),ng
         igmin=ig2
         igmax=ig2
         !
         ! ***ig1 is the primary group
         do ig1=iglast(il)-1,igfirs(il),-1
           delta=sigsm(ig1,ig2,il,nz-iz+1)*flux(ig1,il,nz-iz+1)- &
           & sigsm(ig1,ig2,il,1)
           if(abs(delta).gt.1.0e-6*abs(sigsm(ig1,ig1,il,1))) then
             igmin=min(igmin,ig1)
             igmax=max(igmax,ig1)
           endif
         enddo
         ijj(ig2)=igmax
         njj(ig2)=igmax-igmin+1
         do ig1=igmax,igmin,-1
           igar=igar+1
           if(igar.gt.maxgar) call error('drasc','gar overflow.',' ')
           if((ig1.ge.igfirs(il)).and.(ig1.lt.iglast(il))) then
             gar(igar)=real(sigsm(ig1,ig2,il,nz-iz+1)*flux(ig1,il,nz-iz+1)- &
             & sigsm(ig1,ig2,il,1))
             if(abs(gar(igar)).gt.garmax) garmax=abs(gar(igar))
           else
             gar(igar)=0.0d0
           endif
         enddo
       enddo
       !
       ! ***write on the xsm file
       imax=igar
       do i=imax,1,-1
         if(abs(gar(i)).ge.1.0e-8*garmax) go to 30
         igar=igar-1
       enddo
       30 if(igar.gt.0) then
         call xsmput(draglib,'SCAT'//cdl,gar(1:igar))
         call xsmput(draglib,'NJJS'//cdl,njj(1:ng))
         call xsmput(draglib,'IJJS'//cdl,ijj(1:ng))
       endif
       call xsmsix(draglib,' ',2)
     enddo
   enddo
   !**end of Legendre order loop
   !
   ! scratch storage deallocation
   deallocate(sigsm,flux,rm)
   deallocate(gar,scr)
   return
   end subroutine drasc
   !
   subroutine fmatr(maxgr,maxnl,maxnz,nin,matd,mt,ytemp,ng,nl,nz,rm,flux,exist)
   !-----------------------------------------------------------------
   !   utility routine for recovering a matrix reaction from a gendf file.
   !   input parameters:
   !   maxgr   first and second dimensions of matrix rm. Maximum number
   !           of energy groups.
   !   maxnl   third dimension of matrix rm. Maximum number of Legendre orders.
   !   maxnz   4-th dimensions of matrix rm. Maximum number of dilutions.
   !   nin     file unit number of the gendf file.
   !   matd    material number.
   !   mt      reaction index.
   !   ytemp   absolute temperature (Kelvin).
   !   output parameters:
   !   ng      number of energy groups.
   !   nl      number of Legendre orders.
   !   nz      number of dilutions.
   !   rm      matrix reaction (primary->secondary,Legendre,dilution).
   !   flux    flux (primary,Legendre,dilution).
   !   exist   matrix reaction existence flag.
   !-----------------------------------------------------------------
   use endf   ! provides endf routines and variables
   use util   ! provides error
   use util   ! provides timer,openz,repoz,error
   integer :: maxa,lz
   parameter (maxa=200000,lz=6)
   integer maxgr,maxnl,maxnz,nin,matd,mt,ng,nl,nz,ig,ig2,ig2lo,il,iz,jg, &
   & k,loc,loca,locf,nb,ng2,nw
   logical exist,lfind
   real(kr) ytemp,aa(6),rm(maxgr,maxgr,maxnl,maxnz),flux(maxgr,maxnl,maxnz)
   real(kr),allocatable,dimension(:) :: scr,cspc
   !
   ! scratch storage allocation
   allocate(scr(maxa),cspc(maxgr))
   !
   rm(:maxgr,:maxgr,:maxnl,:maxnz)=0.0
   flux(:maxgr,:maxnl,:maxnz)=0.0
   cspc(:maxgr)=0.0
   call repoz(nin)
   10 lfind=.false.
   call skiprz(nin,1)
   do while (.not.lfind)
     if(nin.lt.0) then
       read(-nin,end=900) math,mfh,mth,nb,nw
     else if(nin.gt.0) then
       read(nin,'(6e11.0,i4,i2,i3,i5)',end=900) aa,math,mfh,mth,nsp
     endif
     lfind=(math.eq.matd).and.(mfh.eq.6).and.(mth.eq.mt)
   enddo
   if(.not.lfind) go to 900
   if(nin.lt.0) then
     call skiprz(nin,-1)
     read(-nin) math,mfh,mth,nb,nw,aa
   endif
   nl=nint(aa(3))
   nz=nint(aa(4))
   ng=nint(aa(6))
   call listio(nin,0,0,scr(1),nb,nw)
   if(abs(scr(1)-ytemp).gt.1.0e-3) then
     call tomend(nin,0,0,scr)
     go to 10
   endif
   if(ng.gt.maxgr) call error('fmatr','maxgr overflow',' ')
   if(nl.gt.maxnl) call error('fmatr','maxnl overflow',' ')
   if(nz.gt.maxnz) call error('fmatr','maxnz overflow',' ')
   lfind=.true.
   do while (mth.ne.0) ! loop over energy groups
     if(.not.lfind) call listio(nin,0,0,scr(1),nb,nw)
     ng2=l1h
     ig2lo=l2h
     ig=n2h
     lfind=.false.
     loc=1+nw
     do while (nb.ne.0)
       if(loc+302.gt.maxa) call error('fmatr','endf input size exceeded',' ')
       call moreio(nin,0,0,scr(loc),nb,nw)
       loc=loc+nw
     enddo
     if(ig.ne.0) then
       jg=ng-ig+1
       do il=1,nl
         do iz=1,nz
           locf=1+lz+nl*(iz-1)+(il-1)
           flux(jg,il,iz)=scr(locf)
           do k=2,ng2
             loca=locf+nl*nz*(k-1)
             if(ig2lo.ne.0) then
               ! matrix part
               ig2=ig2lo+k-2
               rm(jg,ng-ig2+1,il,iz)=rm(jg,ng-ig2+1,il,iz)+scr(loca)
             else
               ! spectrum part
               do ig2=1,ng
                 rm(jg,ng-ig2+1,il,iz)=rm(jg,ng-ig2+1,il,iz)+cspc(ig2)*scr(loca)
               enddo
             endif
           enddo
         enddo
       enddo
     else
       ! save constant spectrum
       do k=1,ng2
         cspc(ig2lo+k-1)=scr(1+lz+nl*nz*(k-1))
       enddo
     endif
   enddo
   exist=.true.
   deallocate(cspc,scr)
   return
   !
   900 exist=.false.
   ng=0
   nl=0
   nz=0
   !
   ! scratch storage deallocation
   deallocate(cspc,scr)
   return
   end subroutine fmatr
   !
   subroutine fvect(maxgr,maxnl,maxnz,nin,matd,mt,ytemp,ng,nl,nz,rv,flux,exist)
   !-----------------------------------------------------------------
   !   utility routine for recovering a vector reaction from a gendf file.
   !   input parameters:
   !   maxgr   first dimension of vector rv. Maximum number of energy groups.
   !   maxnl   second dimension of vector rv. Maximum number of Legendre orders.
   !   maxnz   third dimensions of vector rv. Maximum number of dilutions.
   !   nin     file unit number of the gendf file.
   !   matd    material number.
   !   mt      reaction index.
   !   ytemp   absolute temperature (Kelvin).
   !   output parameters:
   !   ng      number of energy groups.
   !   nl      number of Legendre orders.
   !   nz      number of dilutions.
   !   rv      vector reaction (group,Legendre,dilution).
   !   flux    flux (group,Legendre,dilution).
   !   exist   vector reaction existence flag.
   !-----------------------------------------------------------------
   use endf   ! provides endf routines and variables
   use util   ! provides error
   use util   ! provides timer,openz,repoz,error
   integer :: maxa,lz
   parameter (maxa=2000,lz=6)
   integer maxgr,maxnl,maxnz,nin,matd,mt,ng,nl,nz,ig,il,iz,jg,loc,locf,nb, &
   & ng2, nw
   logical exist,lfind
   real(kr) ytemp,rv(maxgr,maxnl,maxnz),flux(maxgr,maxnl,maxnz),aa(6)
   real(kr),allocatable,dimension(:) :: scr
   !
   allocate(scr(maxa))
   rv(:maxgr,:maxnl,:maxnz)=0.0
   flux(:maxgr,:maxnl,:maxnz)=0.0
   call repoz(nin)
   10 lfind=.false.
   call skiprz(nin,1)
   do while (.not.lfind)
     if(nin.lt.0) then
       read(-nin,end=900) math,mfh,mth,nb,nw
     else if(nin.gt.0) then
       read(nin,'(6e11.0,i4,i2,i3,i5)',end=900) aa,math,mfh,mth,nsp
     endif
     lfind=(math.eq.matd).and.(mfh.eq.3).and.(mth.eq.mt)
   enddo
   if(.not.lfind) go to 900
   if(nin.lt.0) then
     call skiprz(nin,-1)
     read(-nin) math,mfh,mth,nb,nw,aa
   endif
   nl=nint(aa(3))
   nz=nint(aa(4))
   ng=nint(aa(6))
   call listio(nin,0,0,scr(1),nb,nw)
   if(abs(scr(1)-ytemp).gt.1.0e-3) then
     call tomend(nin,0,0,scr)
     go to 10
   endif
   if(ng.gt.maxgr) call error('fvect','maxgr overflow',' ')
   if(nl.gt.maxnl) call error('fvect','maxnl overflow',' ')
   if(nz.gt.maxnz) call error('fvect','maxnz overflow',' ')
   lfind=.true.
   do while (mth.ne.0) ! loop over energy groups
     if(.not.lfind) call listio(nin,0,0,scr(1),nb,nw)
     ng2=l1h
     ig=n2h
     lfind=.false.
     loc=1+nw
     do while (nb.ne.0)
       if(loc+302.gt.maxa) call error('fvect','endf input size exceeded', &
       & ' ')
       call moreio(nin,0,0,scr(loc),nb,nw)
       loc=loc+nw
     enddo
     jg=ng-ig+1
     do il=1,nl
       do iz=1,nz
         locf=1+lz+nl*(iz-1)+(il-1)
         flux(jg,il,iz)=scr(locf)
         rv(jg,il,iz)=scr(locf+nl*nz)
       enddo
     enddo
   enddo
   exist=.true.
   deallocate(scr)
   return
   !
   900 exist=.false.
   ng=0
   nl=0
   nz=0
   deallocate(scr)
   return
   end subroutine fvect
   !
   subroutine dradep(nfp,ndcy,idecay,iprint)
   !-----------------------------------------------------------------
   !   compute depletion-related information.
   !-----------------------------------------------------------------
   use mainio ! provides nsysi,contio,nsyso,nsyse
   use endf   ! provides endf routines and variables
   use util   ! provides timer,openz,repoz,error
   integer :: maxa,maxiso,nreac,nfath,maxch,maxen,maxfis,maxter
   parameter(maxa=10000,maxiso=4000,nreac=14,nfath=50,maxch=800,maxen=4, &
   & maxfis=100,maxter=4)
   integer nfp,ndcy,idecay,iprint,izae,nbch,nbfiss,nbfp,nbfpch,nbiso,nw
   integer mylist(maxiso,3)
   character(len=4) hiso(3,maxiso)
   integer i,i1,ia,ifp,ifps,igar,ii,ile,ind,iof,ipos,iso,itext4,iz,j, &
   & jpos,k,lep1,loc,maxfp,nb,nbdpf,ninter
   real(kr) awr,energy(maxter),za,eeee,t1,t2,t3
   character hname*8,hname2*8,hsmg*131,hich(maxch)*8,text4*4
   integer,allocatable,dimension(:,:) :: idreac,ipreac
   real(kr),allocatable,dimension(:) :: ddeca,scr
   real(kr),allocatable,dimension(:,:) :: br,dener,prate
   real(kr),allocatable,dimension(:,:,:) :: yield,terp
   character(len=8),allocatable,dimension(:,:) :: hrch
   real(kr), dimension(maxen) :: enrgs
   !
   !**read the specification lines for the isotopes present in the
   !  burnup chain. The specification is:
   !   [[
   !   hich [[ hrch br ]] /
   !   ]]
   !   [[ 'interpol' 'at' energy / ]]
   !   end /
   !   where hich  : character*8 name of the Draglib isotope
   !         hrch  : character*8 name of a neutron induced reaction (not
   !                 a scattering type reaction)
   !         br    : branching ratio to an isomeric daughter
   !         energy: incident neutron energy with MT454 (default = 0.0).
   allocate(scr(maxa),br(nfath,maxch),hrch(nfath,maxch))
   mylist(:maxiso,:3)=0
   nbch=0
   nbiso=0
   ninter=0 ! number of incident neutron energies for fission yields
   do
     nbch=nbch+1
     if(nbch.gt.maxch) call error('dradep','maxch overflow',' ')
     hrch(:nfath,nbch)=' '
     br(:nfath,nbch)=0.0
     read(nsysi,*) hich(nbch),(hrch(i,nbch),br(i,nbch),i=1,nfath)
     if(hich(nbch).eq.'end') then
       go to 90
     else if(hich(nbch).eq.'interpol') then
       ! special card to set the incident neutron energy with MT454
       if(hrch(1,nbch).ne.'at') call error('dradep','at keyword expected',' ')
       ninter=ninter+1
       if(ninter.gt.maxter) call error('dradep','maxter overflow',' ')
       energy(ninter)=br(1,nbch)
       nbch=nbch-1
       cycle
     else
       do i=1,nfath
         if(br(i,nbch).ne.0.0) then
           write(nsyso,'('' isomeric '',a,'' branching ratio='',f7.3, &
           & '' set for '',a)') hrch(i,nbch),br(i,nbch),hich(nbch)
         endif
       enddo
     endif
     izae=0
     i1=index(hich(nbch),'_')
     if(i1.eq.0) then
       hname=hich(nbch)
     else
       hname=hich(nbch)(:i1-1)
     endif
     call dranam(izae,hname)
     nbiso=nbiso+1
     if(nbiso.gt.maxiso) call error('dradep','maxiso overflow',' ')
     mylist(nbiso,1)=izae
     mylist(nbiso,3)=nbch
   enddo
   90 nbch=nbch-1
   if(ninter.eq.0) then
     ! assume thermal fission yields by default
     ninter=1
     energy(1)=0.0
   endif
   if(nbch.eq.0) return
   !
   !**set the list of depleting isotopes
   allocate(terp(maxen,maxfis,ninter))
   nbfpch=0
   call repoz(nfp)
   100 if(nfp.gt.0) then
     read(nfp,'(66x,i4,i2,i3)') math,mfh,mth
   else
     read(-nfp) math,mfh,mth
   endif
   if(math.eq.0) go to 100
   nbfiss=0
   if((mfh.ne.0).and.(mth.ne.0)) call skiprz(nfp,-1)
   do
     150 call contio(nfp,0,0,scr,nb,nw)
     if(math.eq.-1) then
       go to 300
     else if(mfh.ne.1) then
       call tofend(nfp,0,0,scr)
       go to 150
     else if(mth.ne.451) then
       call tosend(nfp,0,0,scr)
       go to 150
     endif
     nbfiss=nbfiss+1
     if(nbfiss.gt.maxfis) call error('dradep','maxfis overflow',' ')
     za=c1h
     awr=c2h
     call contio(nfp,0,0,scr,nb,nw)
     izae=10*nint(za+0.1)+l2h
     if(izae.eq.0) call error('dradep','izae error',' ')
     call draind(nbiso,izae,mylist(1,1),ind)
     if(ind.le.0) then
       nbiso=nbiso+1
       if(nbiso.gt.maxiso) call error('dradep','maxiso overflow',' ')
       mylist(nbiso,1)=izae
     endif
     call dranam(izae,hname)
     call findf(math,8,454,nfp)
     call contio(nfp,0,0,scr,nb,nw)
     lep1=l1h
     if(lep1.gt.maxen) then
       write(hsmg,'(44hinvalid number of incident neutron energies., &
       & 18h Increase maxen to,i3,1h.)') lep1
       call error('dradep',hsmg,' ')
     endif
     do ile=1,lep1
       call listio(nfp,0,0,scr(1),nb,nw)
       nbfp=n2h
       enrgs(ile)=c1h
       loc=1+nw
       do while (nb.ne.0)
         if(loc+302.gt.maxa) call error('dradep','endf input size exceeded',' ')
         call moreio(nfp,0,0,scr(loc),nb,nw)
         loc=loc+nw
       enddo
       if(enrgs(ile).gt.2.0e6) cycle
       do ifp=1,nbfp
         iof=6+(ifp-1)*4+1
         if(scr(iof+2).gt.1.0e-10) then
           iz=nint(scr(iof)/1000+0.1)
           ifps=nint(scr(iof+1)+0.1)
           ia=mod(nint(scr(iof)+0.1),1000)
           izae=10000*iz+10*ia+ifps
           call draind(nbiso,izae,mylist(1,1),ind)
           if(ind.le.0) then
             nbiso=nbiso+1
             if(nbiso.gt.maxiso) call error('dradep','maxiso overflow',' ')
             mylist(nbiso,1)=izae
             call dranam(izae,hname2)
             do i=1,nbch
               if(hich(i).eq.hname2) then
                 nbfpch=nbfpch+1
                 go to 170
               endif
             enddo
           endif
           170 continue
         endif
       enddo
     enddo
     ! compute terp factors for MT454
     do i=1,ninter
       eeee=energy(i)
       terp(:maxen,nbfiss,i)=0.0
       if((lep1.ne.1).and.(energy(i).ne.0.0)) then
         if((energy(i).lt.enrgs(1)).or.(energy(i).gt.enrgs(lep1))) then
           write(hsmg,'(a,30h: invalid interpolation energy,1p,e11.3, &
           & 28h for MT454 not in interval (,2e11.3,2h).)') trim(hname), &
           & energy(i),enrgs(1),enrgs(lep1)
           call mess('dradep',hsmg,' ')
         endif
       endif
       if((lep1.eq.1).or.(energy(i).lt.enrgs(1))) then
         eeee=enrgs(1)
         terp(1,nbfiss,i)=1.0
       else if(energy(i).gt.enrgs(lep1)) then
         eeee=enrgs(lep1)
         terp(lep1,nbfiss,i)=1.0
       else
         t1=log(1.0d7/energy(i))
         do ile=2,lep1
           if(energy(i).le.enrgs(ile)) then
             t2=log(1.0d7/enrgs(ile-1))
             t3=log(1.0d7/enrgs(ile))
             terp(ile-1,nbfiss,i)=(t3-t1)/(t3-t2)
             terp(ile,nbfiss,i)=1.0-terp(ile-1,nbfiss,i)
             exit
           endif
         enddo
       endif
       if(iprint.gt.0) then
         write(nsyso,'( &
         &  '' incident neutron energy with MT454 ..........  '',i4,1x,a8,1h:, &
         &  1p,e11.3,3h eV)') nbfiss,hname,eeee
       endif
     enddo
     call tomend(nfp,0,0,scr)
   enddo
   300 nbdpf=nbiso-nbfiss
   nbfpch=nbch
   !
   call repoz(ndcy)
   400 if(ndcy.gt.0) then
     read(ndcy,'(66x,i4,i2,i3)') math,mfh,mth
   else
     read(-ndcy) math,mfh,mth
   endif
   if(math.eq.0) go to 400
   if((mfh.ne.0).and.(mth.ne.0)) call skiprz(ndcy,-1)
   do
     500 call contio(ndcy,0,0,scr,nb,nw)
     if(math.eq.-1) then
       go to 600
     else if(mfh.ne.1) then
       call tofend(ndcy,0,0,scr)
       go to 500
      else if(mth.ne.451) then
       call tosend(ndcy,0,0,scr)
       go to 500
     endif
     za=c1h
     call contio(ndcy,0,0,scr,nb,nw)
     izae=10*nint(za+0.1)+l2h
     call draind(nbiso,izae,mylist(1,1),ind)
     if(ind.le.0) then
       nbiso=nbiso+1
       if(nbiso.gt.maxiso) call error('dradep','maxiso overflow',' ')
       mylist(nbiso,1)=izae
       mylist(nbiso,2)=math
     else
       if((mylist(ind,2).ne.math).and.(mylist(ind,2).ne.0)) then
         write(hsmg,'(''invalid mat index (='',i6,''; expected='', &
         & i6,'' izae='',i8,'')'')') math,mylist(ind,2),izae
         call error('dradep',hsmg,' ')
       endif
       mylist(ind,2)=math
     endif
     call tomend(ndcy,0,0,scr)
   enddo
   !
   !*sort the izae entries.
   600 do k=2,nbiso
     i=k-1
     ipos=0
     if(mylist(i,3).gt.0) ipos=index(hich(mylist(i,3)),'_')+1
     do j=k,nbiso
       jpos=0
       if(mylist(j,3).gt.0) jpos=index(hich(mylist(j,3)),'_')+1
       if((mylist(i,1).gt.mylist(j,1)).or. &
       & ((mylist(i,1).eq.mylist(j,1)).and.(ipos.gt.jpos))) then
         do ii=1,3
           igar=mylist(i,ii)
           mylist(i,ii)=mylist(j,ii)
           mylist(j,ii)=igar
         enddo
       endif
     enddo
   enddo
   !
   text4=' '
   read(text4,'(a4)') itext4
   do iso=1,nbiso
     if(mylist(iso,3).ne.0) then
       hname=hich(mylist(iso,3))
     else
       call dranam(mylist(iso,1),hname)
     endif
     hiso(1,iso)=hname(1:4)
     hiso(2,iso)=hname(5:8)
     hiso(3,iso)=' '
   enddo
   !
   maxfp=nbdpf+50 ! reserve 50 location for lumped fp daughters
   allocate(idreac(nreac,nbiso),ipreac(nfath,nbiso))
   allocate(dener(nreac,nbiso),ddeca(nbiso),prate(nfath,nbiso), &
   & yield(ninter,nbfiss,maxfp))
   call draevo(maxfp,nbiso,nbfiss,nbdpf,ninter,nreac,nfath,mylist(1,1),nfp, &
   & ndcy,idreac,dener,ddeca,ipreac,prate,yield,terp,iprint)
   deallocate(terp)
   !
   !**lump the burnup chain from nbiso to nbch isotopes
   call dralum(maxfp,nbiso,nbfiss,nbdpf,ninter,nreac,nfath,mylist(1,1),hiso, &
   & nbch,nbfpch,hich,hrch,br,idreac,dener,ddeca,ipreac,prate,yield,idecay, &
   & energy,iprint)
   deallocate(yield,prate,ddeca,dener)
   deallocate(ipreac,idreac)
   deallocate(hrch,br,scr)
   return
   end subroutine dradep
   !
   subroutine draevo(maxfp,nbiso,nbfiss,nbdpf,ninter,nreac,nfath,mylist,nfp, &
   & ndcy,idreac,dener,ddeca,ipreac,prate,yield,terp,iprint)
   !-----------------------------------------------------------------
   !   compute blocks 'DEPLETE-IDEN', 'DEPLETE-REAC', 'DEPLETE-ENER',
   !   'DEPLETE-DECA', etc. for the non-lumped chain
   !-----------------------------------------------------------------
   use mainio ! provides nsysi,contio,nsyso,nsyse
   use endf   ! provides endf routines and variables
   use util   ! provides timer,openz,repoz,error
   integer :: maxa,maxen,maxfis,lz
   parameter(maxa=10000,maxen=4,maxfis=100,lz=6)
   integer maxfp,nbiso,nbfiss,nbdpf,ninter,nreac,nfath,nfp,ndcy,i,ia,idy, &
   & ifath,ifis,ifp,ifpp,ifps,ifpss,ile,ind,iof,iz,iza,izae,ja,jnd,jz,jzae, &
   & lep1,loc,nb,nbdy,nbfp,nw,iprint
   real(kr) za,rtyp
   integer mylist(nbiso),idreac(nreac,nbiso),ipreac(nfath,nbiso)
   real(kr) awr,energy,dener(nreac,nbiso),ddeca(nbiso),prate(nfath,nbiso), &
   & yield(ninter,nbfiss,maxfp),terp(maxen,maxfis,ninter)
   integer, allocatable, dimension(:) :: indpf
   real(kr), allocatable, dimension(:) :: summ,scr
   character text6*6,hname*8
   !
   allocate(indpf(nbdpf),scr(maxa),summ(ninter))
   ddeca(:nbiso)=0.0
   idreac(:nreac,:nbiso)=0
   dener(:nreac,:nbiso)=0.0
   ipreac(:nfath,:nbiso)=0
   prate(:nfath,:nbiso)=0.0
   indpf(:nbdpf)=0
   yield(:ninter,:nbfiss,:nbdpf)=0.0
   !
   call repoz(nfp)
   ifpss=0
   100 if(nfp.gt.0) then
     read(nfp,'(66x,i4,i2,i3)') math,mfh,mth
   else
     read(-nfp) math,mfh,mth
   endif
   if(math.eq.0) go to 100
   if((mfh.ne.0).and.(mth.ne.0)) call skiprz(nfp,-1)
   ifis=0
   if(iprint.gt.0) then
     write(nsyso,'(/50h draevo: check initial fission yield normalization)')
   endif
   do
     150 call contio(nfp,0,0,scr,nb,nw)
     if(math.eq.-1) then
       go to 300
     else if(mfh.ne.1) then
       call tofend(nfp,0,0,scr)
       go to 150
     else if(mth.ne.451) then
       call tosend(nfp,0,0,scr)
       go to 150
     endif
     ifis=ifis+1
     if(ifis.gt.nbfiss) call error('draevo','nbfiss overflow',' ')
     za=c1h
     awr=c2h
     call contio(nfp,0,0,scr,nb,nw)
     izae=10*nint(za+0.1)+l2h
     if(izae.eq.0) call error('draevo','izae error(1)',' ')
     call draind(nbiso,izae,mylist,ind)
     if(ind.le.0) then
       call dranam(izae,hname)
       call error('draevo','missing isotope(1):'//hname,' ')
     endif
     idreac(2,ind)=ifis*100+4 ! ind is producing fp
     call findf(math,8,454,nfp)
     call contio(nfp,0,0,scr,nb,nw)
     lep1=l1h
     summ(:ninter)=0.0
     do ile=1,lep1
       call listio(nfp,0,0,scr(1),nb,nw)
       nbfp=n2h
       energy=c1h
       loc=1+nw
       do while (nb.ne.0)
         if(loc+302.gt.maxa) call error('draevo','endf input size exceeded',' ')
         call moreio(nfp,0,0,scr(loc),nb,nw)
         loc=loc+nw
       enddo
       if(energy.le.2.0e6) then
         do ifp=1,nbfp
           iof=6+(ifp-1)*4+1
           if(scr(iof+2).gt.1.0e-10) then
             iz=nint(scr(iof)/1000+0.1)
             ifps=nint(scr(iof+1)+0.1)
             ia=mod(nint(scr(iof)+0.1),1000)
             jzae=10000*iz+10*ia+ifps
             if(jzae.eq.0) call error('draevo','jzae error',' ')
             call draind(nbiso,jzae,mylist,jnd)
             if(jnd.le.0) then
               call dranam(jzae,hname)
               call error('draevo','missing isotope(2):'//hname,' ')
             endif
             do i=1,ifpss
               if(indpf(i).eq.jzae) then
                 ifpp=i
                 go to 200
               endif
             enddo
             ifpss=ifpss+1
             if(ifpss.gt.nbdpf) call error('draevo','nbdpf overflow',' ')
             indpf(ifpss)=jzae
             ifpp=ifpss
             200 idreac(2,jnd)=ifpp*100+5 ! jnd is a fission fragment
             do i=1,ninter
               yield(i,ifis,ifpp)=yield(i,ifis,ifpp)+terp(ile,ifis,i)*scr(iof+2)
               summ(i)=summ(i)+terp(ile,ifis,i)*scr(iof+2)
             enddo
           endif
         enddo
       endif
     enddo
     call tomend(nfp,0,0,scr)
     if(iprint.gt.0) then
       call dranam(izae,hname)
       write(nsyso,'(1x,i5,1x,a9,10f9.4)') ifis,hname,summ(:ninter)
     endif
   enddo
   !
   300 call repoz(ndcy)
   400 if(ndcy.gt.0) then
     read(ndcy,'(66x,i4,i2,i3)') math,mfh,mth
   else
     read(-ndcy) math,mfh,mth
   endif
   if(math.eq.0) go to 400
   if((mfh.ne.0).and.(mth.ne.0)) call skiprz(ndcy,-1)
   do
     500 call contio(ndcy,0,0,scr,nb,nw)
     if(math.eq.-1) then
       go to 700
     else if(mfh.ne.1) then
       call tofend(ndcy,0,0,scr)
       go to 500
     else if(mth.ne.451) then
       call tosend(ndcy,0,0,scr)
       go to 500
     endif
     za=c1h
     awr=c2h
     call contio(ndcy,0,0,scr,nb,nw)
     izae=10*nint(za+0.1)+l2h
     if(izae.eq.0) call error('draevo','izae error(2)',' ')
     call draind(nbiso,izae,mylist,ind)
     if(ind.le.0) then
       call dranam(izae,hname)
       call error('draevo','missing isotope(3):'//hname,' ')
     endif
     idreac(1,ind)=1 ! ind can decay
     call findf(math,8,457,ndcy)
     call contio(ndcy,0,0,scr,nb,nw)
     iza=nint(c1h+0.1)
     iz=nint(c1h/1000+0.1)
     ia=mod(iza,1000)
     call listio(ndcy,0,0,scr,nb,nw)
     if(c1h.eq.0.0) then
       ddeca(ind)=0.0
     else
       ddeca(ind)=1.0e8*log(2.0)/c1h
     endif
     ! ***recover energy for n1h/2 decay reactions
     dener(1,ind)=0.0
     do i=1,n1h,2
       dener(1,ind)=dener(1,ind)+scr(lz+i)*1.0e-6
     enddo
     call listio(ndcy,0,0,scr(1),nb,nw)
     nbdy=n2h
     loc=1+nw
     do while (nb.ne.0)
       if(loc+302.gt.maxa) call error('draevo','endf input size exceeded',' ')
       call moreio(ndcy,0,0,scr(loc),nb,nw)
       loc=loc+nw
     enddo
     do idy=1,nbdy
       rtyp=scr(6*idy+1)
       jz=iz
       ja=ia
       if(abs(rtyp-1.0).le.0.00001) then
         ! **beta- decay
         jz=iz+1
       else if(abs(rtyp-1.1).le.0.00001) then
         ! **two beta- decays
         jz=iz+2
       else if(abs(rtyp-1.4).le.0.00001) then
         ! **beta- decay followed by alpha emission
         jz=iz-1
         ja=ia-4
       else if(abs(rtyp-1.5).le.0.00001) then
         ! **beta- decay followed by neutron emission (delayed neutron)
         jz=iz+1
         ja=ia-1
       else if(abs(rtyp-1.55).le.0.00001) then
         ! **beta- decay followed by two neutrons emission (delayed neutrons)
         jz=iz+1
         ja=ia-2
       else if(abs(rtyp-1.555).le.0.00001) then
         ! **beta- decay followed by three neutrons emission (delayed neutrons)
         jz=iz+1
         ja=ia-3
       else if(abs(rtyp-2.0).le.0.00001) then
         ! **beta+ decay or electron capture
         jz=iz-1
       else if(abs(rtyp-2.4).le.0.00001) then
         ! **beta+ decay followed by alpha emission
         jz=iz-3
         ja=ia-4
       else if(abs(rtyp-2.7).le.0.00001) then
         ! **beta+ decay followed by proton emission
         jz=iz-2
         ja=ia-1
       else if(abs(rtyp-2.77).le.0.00001) then
         ! **beta+ decay followed by two protons emission
         jz=iz-3
         ja=ia-2
       else if(abs(rtyp-3.0).le.0.00001) then
         ! **isomeric transition
       else if(abs(rtyp-3.4).le.0.00001) then
         ! **isomeric transition followed by alpha decay
         jz=iz-2
         ja=ia-4
       else if(abs(rtyp-3.5).le.0.00001) then
         ! **isomeric transition followed by neutron emission (delayed neutron)
         ja=ia-1
       else if(abs(rtyp-4.0).le.0.00001) then
         ! **alpha decay
         jz=iz-2
         ja=ia-4
       else if(abs(rtyp-5.0).le.0.00001) then
         ! **neutron emission (not delayed neutron)
         ja=ia-1
       else if(abs(rtyp-5.5).le.0.00001) then
         ! **two neutrons emission (not delayed neutron)
         ja=ia-2
       else if(abs(rtyp-5.55).le.0.00001) then
         ! **three neutrons emission (not delayed neutron)
         ja=ia-3
       else if(abs(rtyp-6.0).le.0.00001) then
         ! **spontaneous fission
         jz=0
       else if(abs(rtyp-7.0).le.0.00001) then
         ! **proton emission
         jz=iz-1
         ja=ia-1
       else if(abs(rtyp-7.7).le.0.00001) then
         ! **two protons emission
         jz=iz-2
         ja=ia-2
       else if(abs(rtyp-10.0).le.0.00001) then
         ! **unknown origin
         jz=0
       else
         call dranam(izae,hname)
         write(text6,'(f6.4)') rtyp
         write(nsyso,'(A)')'unknown type of decay: '//text6//' for '//hname
         jz=0
       endif
       if(jz.ne.0) then
         izae=10000*jz+10*ja+nint(scr(6*idy+2)+0.1)
         call draind(nbiso,izae,mylist,jnd)
         if(jnd.gt.0) then
           ifath=0
           do i=1,nfath
             if(ipreac(i,jnd).eq.0) then
               ifath=i
               go to 600
             endif
           enddo
           call error('draevo','nfath overflow',' ')
           600 ipreac(ifath,jnd)=ind*100+1 ! decay father identification
           prate(ifath,jnd)=scr(6*idy+5)
         endif
       endif
     enddo
     call tomend(ndcy,0,0,scr)
   enddo
   700 deallocate(summ,scr,indpf)
   return
   end subroutine draevo
   !
   subroutine draind(nb,izae,mylist,ind)
   !-----------------------------------------------------------------
   !   find the position ind of material izae in array mylist.
   !-----------------------------------------------------------------
   integer nb,izae,mylist(nb),ind,i
   do i=1,nb
     ind=i
     if(mylist(i).eq.izae) return
   enddo
   ind=-1
   return
   end subroutine draind
   !
   subroutine dralum(maxfp,nbiso,nbfiss,nbdpf,ninter,nreac,nfath,mylist,hiso, &
   & nbch,nbfpch,hich,hrch,br,idreac,dener,ddeca,ipreac,prate,yield,idecay, &
   & energy,iprint)
   !-----------------------------------------------------------------
   !   complete and lump the burnup chain from nbiso to nbch isotopes.
   !   write the lumped chain on the xsm file.
   !-----------------------------------------------------------------
   use mainio ! provides nsysi,contio,nsyso,nsyse
   use util   ! provides error
   integer :: maxrea,nstate,maxit,nfath
   parameter (maxrea=14,nstate=40,maxit=20)
   integer maxfp,nbiso,nbfiss,nbdpf,ninter,nreac,nbch,nbfpch,idecay,iprint, &
   & ireac,iter,iso,isoo,i,j,ia,ja,ibfp,ida,ifa,ifath,ifi,ifp,ifps,iii,im, &
   & ind,ipgar,iz,j0,jfath,jfp,jnd,jnd1,jnd2,jso,jz,k,knd,kreac,kt,nbheav, &
   & nn,nlump,ja2,jz2,indifi
   real(kr) prgar,ymax
   real(kr) br(nfath,nbch)
   character hrch(nfath,nbch)*8,hich(nbch)*8,hname*8,text4*4,hsmg*131
   integer mylist(nbiso),idreac(nreac,nbiso),ipreac(nfath,nbiso),istate(nstate)
   character(len=4) hiso(3,nbiso),hreac(2,maxrea),in(2)
   real(kr) over_half,dener(nreac,nbiso),ddeca(nbiso),prate(nfath,nbiso), &
   & yield(ninter,nbfiss,maxfp),energy(ninter)
   integer, allocatable, dimension(:,:) :: jpreac,jdreac,ipos
   character(len=4), allocatable, dimension(:,:) :: hhhh
   real, allocatable, dimension(:) :: dddd,eeee
   real, allocatable, dimension(:,:) :: rrate,eener
   real, allocatable, dimension(:,:,:) :: eyiel
   real(kr), allocatable, dimension(:) :: delayed,cyield
   character(len=8), save, dimension(maxrea) :: reac= &
   & (/ 'DECAY   ','NFTOT   ','NG      ','N2N     ','N3N     ','N4N     ', &
   &    'NA      ','NP      ','N2A     ','NNP     ','ND      ','NT      ', &
   &    'NHE3    ','NODECAY ' /)
   !
   if(nreac.ne.maxrea) call error('dralum','maxrea overflow',' ')
   allocate(jpreac(nfath,nbch),jdreac(nreac-1,nbch),ipos(nbch,2),hhhh(3,nbch))
   allocate(rrate(nfath,nbch),eener(nreac-1,nbch),eyiel(ninter,nbfiss,nbfpch), &
   & dddd(nbch))
   !
   ! **find the position of the lumped isotopes in the complete chain
   do iso=1,nbch
     j0=0
     in(1)=hich(iso)(1:4)
     in(2)=hich(iso)(5:8)
     do jso=1,nbiso
       j0=jso
       if((in(1).eq.hiso(1,jso)).and.(in(2).eq.hiso(2,jso))) go to 10
     enddo
     call error('dralum','unable to find '//hich(iso),' ')
     10 ipos(iso,1)=j0
   enddo
   !
   if(iprint.gt.0) then
     write(nsyso,'(/50h dralum: check initial fission yield normalization)')
     do iso=1,nbch
       ida=ipos(iso,1)
       if(mod(idreac(2,ida),100).ne.4) cycle
       ifi=idreac(2,ida)/100
       if(ifi.eq.0) cycle
       do i=1,ninter
         write(nsyso,'(i3,1x,2a4,f9.4)') i,hiso(:2,ida),sum(yield(i,ifi,:maxfp))
       enddo
     enddo
   endif
   !
   ! **complete idreac, dener, ipreac and prate.
   do iso=1,nbch
     ind=ipos(iso,1)
     iz=mylist(ind)/10000
     ia=mod(mylist(ind),10000)/10
     ifps=mod(mylist(ind),10)
     do i=1,nfath
       if(hrch(i,iso).eq.' ') go to 30
       ireac=0
       jz=iz
       ja=0
       jz2=0
       ja2=0
       if(hrch(i,iso).eq.'nftot') then
         jz=0
         if(idreac(2,ind).eq.0) idreac(2,ind)=3
         ireac=2
       else if(hrch(i,iso).eq.'NODECAY') then
         jz=0
         if(idreac(2,ind).eq.0) idreac(2,ind)=3
         ireac=14
       else if(hrch(i,iso).eq.'ng') then
         ja=ia+1
         idreac(3,ind)=1
         ireac=3
       else if(hrch(i,iso).eq.'n2n') then
         ja=ia-1
         idreac(4,ind)=1
         ireac=4
       else if(hrch(i,iso).eq.'n3n') then
         ja=ia-2
         idreac(5,ind)=1
         ireac=5
       else if(hrch(i,iso).eq.'n4n') then
         ja=ia-3
         idreac(6,ind)=1
         ireac=6
       else if(hrch(i,iso).eq.'na') then
         jz=iz-2
         ja=ia-3
         idreac(7,ind)=1
         ireac=7
         jz2=2 ; ja2=4 ! tertiary particle=He4
       else if(hrch(i,iso).eq.'np') then
         jz=iz-1
         ja=ia
         idreac(8,ind)=1
         ireac=8
         jz2=1 ; ja2=1 ! tertiary particle=H1
       else if(hrch(i,iso).eq.'n2a') then
         jz=iz-4
         ja=ia-7
         idreac(9,ind)=1
         ireac=9
         jz2=2 ; ja2=4 ! tertiary particle=He4
       else if(hrch(i,iso).eq.'nnp') then
         jz=iz-1
         ja=ia-1
         idreac(10,ind)=1
         ireac=10
         jz2=1 ; ja2=1 ! tertiary particle=H1
       else if(hrch(i,iso).eq.'nd') then
         jz=iz-1
         ja=ia-1
         idreac(11,ind)=1
         ireac=11
         jz2=1 ; ja2=2 ! tertiary particle=H2
       else if(hrch(i,iso).eq.'nt') then
         jz=iz-1
         ja=ia-2
         idreac(12,ind)=1
         ireac=12
         jz2=1 ; ja2=3 ! tertiary particle=H3
       else if(hrch(i,iso).eq.'nhe3') then
         jz=iz-2
         ja=ia-2
         idreac(13,ind)=1
         ireac=13
         jz2=2 ; ja2=3 ! tertiary particle=He3
       endif
       if(ireac.eq.0) then
         write(hsmg,'(8hisotope ,2A4,16h has no reaction)') hiso(:2,ind)
         call error('dralum',hsmg,' ')
       endif
       !
       ! tertiary production of a light nucleus
       if((jz2.ne.0).and.(ja2.ne.0).and.(iz.le.10)) then
         ! only consider tertiary production for lighter father isotopes
         call draind(nbiso,10000*jz2+10*ja2,mylist,jnd)
         if(jnd.eq.-1) then
           write(hsmg,'(21hno isotope with izae=,i10)') 10000*jz2+10*ja2
           call error('dralum',hsmg,' ')
         else if((ja.eq.ja2).and.(jz.eq.jz2)) then
           ! two identical secondary nuclei
           do ifath=1,nfath
             if(ipreac(ifath,jnd).eq.ind*100+ireac) then
               prate(ifath,jnd)=2.0*(1.0-br(i,iso))
               exit
             endif
           enddo
         else
           do ifath=1,nfath
             if(ipreac(ifath,jnd).eq.0) then
               ipreac(ifath,jnd)=ind*100+ireac
               prate(ifath,jnd)=1.0-br(i,iso)
               exit
             endif
           enddo
         endif
       endif
       !
       ! secondary nucleus
       if(jz.ne.0) then
         call draind(nbiso,10000*jz+10*ja,mylist,jnd)
         if(jnd.gt.0) then
           do ifath=1,nfath
             if(ipreac(ifath,jnd).eq.0) then
               ipreac(ifath,jnd)=ind*100+ireac
               prate(ifath,jnd)=1.0-br(i,iso)
               go to 15
             endif
           enddo
           call error('dralum','nfath overflow-1',' ')
         endif
       endif
       15 if((jz.ne.0).and.(br(i,iso).ne.0.0)) then
         ! *** neutron induced production of an isomer
         call draind(nbiso,10000*jz+10*ja+1,mylist,jnd)
         if(jnd.gt.0) then
           do ifath=1,nfath
             if(ipreac(ifath,jnd).eq.0) then
               ipreac(ifath,jnd)=ind*100+ireac
               prate(ifath,jnd)=br(i,iso)
               go to 20
             endif
           enddo
           call error('dralum','nfath overflow-2',' ')
         endif
       endif
       20 continue
     enddo
     30 continue
   enddo
   !
   ! compute cumulative decay heat from lumped isotopes. Assume thermal
   ! incident neutrons.
   allocate(cyield(nbiso),delayed(nbfiss))
   delayed(:nbfiss)=0.0d0
   do ifi=1,nbfiss
     cyield(:nbiso)=0.0d0 ! cumulative yields
     indifi=0
     do ida=1,nbiso
       if(mod(idreac(2,ida),100).eq.4) then
         if(idreac(2,ida)/100.eq.ifi) indifi=ida
       else if(mod(idreac(2,ida),100).eq.5) then
         jfp=idreac(2,ida)/100
         cyield(ida)=yield(1,ifi,jfp)
       endif
     enddo ! ida
     if(indifi.eq.0) call error('dralum','no fissile index',' ')
     loop1: do ida=1,nbiso
       do ifath=1,nfath
         if(ipreac(ifath,ida).eq.0) exit
         if(mod(ipreac(ifath,ida),100).ne.1) cycle
         ind=ipreac(ifath,ida)/100
         cyield(ida)=cyield(ida)+prate(ifath,ida)*cyield(ind)
       enddo
       if(cyield(ida)*dener(1,ida).eq.0.0) cycle
       delayed(ifi)=delayed(ifi)+cyield(ida)*dener(1,ida)
       if(idecay.eq.0) then
         ! decay heat correction to EGD+EB+END for Draglib fission products
         do iso=1,nbch
           ind=ipos(iso,1)
           if(ind.eq.ida) then
             ! ida is a delayed energy precursor which is not lumped. Remove
             ! this energy deposition from DEPLETE-ENER
             dener(2,indifi)=dener(2,indifi)-cyield(ida)*dener(1,ida)
             cycle loop1
           endif
         enddo
       endif
     enddo loop1
   enddo ! ifi
   if(iprint.gt.0) then
     write(nsyso,'(/52h dralum: delayed fission energy from cumulative yiel, &
     & 2hds/12x,10hEND+EGD+EB,5x,26hcorrection in DEPLETE-ENER)')
     do ida=1,nbiso
       if(mod(idreac(2,ida),100).ne.4) cycle
       ifi=idreac(2,ida)/100
       if(ifi.eq.0) cycle
       write(nsyso,'(1x,2a4,f9.4,4h MeV,3x,f9.4,4h MeV)') hiso(:2,ida), &
       & delayed(ifi),dener(2,ida)
     enddo
   endif
   deallocate(delayed,cyield)
   !
   ! lump the decay chains, including independent yields
   iter=0
   40 iter=iter+1
   if(iter.gt.maxit) call error('dralum','too many iterations',' ')
   nlump=0
   do iso=1,nbch
     ind=ipos(iso,1)
     do ifath=1,nfath
       if(ipreac(ifath,ind).eq.0) exit
       if(mod(ipreac(ifath,ind),100).ne.1) go to 50
       jnd=ipreac(ifath,ind)/100
       if(mylist(jnd).eq.0) go to 50
       do j=1,nbch
         if(ipos(j,1).eq.jnd) go to 50
       enddo
       nlump=nlump+1 ! isotope jnd is lumped
       do ida=1,nbiso
         do ifa=1,nfath
           ipgar=ipreac(ifa,ida)
           if((ipgar/100.eq.jnd).and.(mod(ipgar,100).eq.1)) then
             if(mylist(ida).eq.0) go to 50
           endif
         enddo
       enddo
       do ida=1,nbiso
         do ifa=1,nfath
           ipgar=ipreac(ifa,ida)
           if((ipgar/100.eq.jnd).and.(mod(ipgar,100).eq.1)) then
             if(ida.eq.jnd) call error('dralum','bug',' ')
             prgar=prate(ifa,ida)
             do im=ifa,nfath-1
               ipreac(im,ida)=ipreac(im+1,ida)
               prate(im,ida)=prate(im+1,ida)
             enddo
             ipreac(nfath,ida)=0
             prate(nfath,ida)=0.0
             do jfath=1,nfath
               if(ipreac(jfath,jnd).eq.0) exit
               im=nfath+1
               do k=nfath,1,-1
                 if(ipreac(k,ida).eq.ipreac(jfath,jnd)) then
                   prate(k,ida)=prate(k,ida)+prgar*prate(jfath,jnd)
                   go to 45
                 endif
                 if(ipreac(k,ida).eq.0) im=k
               enddo
               if(im.gt.nfath) then
                 write(hsmg,'(25hnfath overflow-3 isotope=,2a4)') hiso(:2,ida)
                 call error('dralum',hsmg,' ')
               endif
               ipreac(im,ida)=ipreac(jfath,jnd)
               prate(im,ida)=prgar*prate(jfath,jnd)
               45 continue
             enddo
             if(mod(idreac(2,jnd),100).eq.5) then
               jfp=idreac(2,jnd)/100
               if(mod(idreac(2,ida),100).eq.5) then
                 ifp=idreac(2,ida)/100
               else
                 nbdpf=nbdpf+1
                 if(nbdpf.gt.maxfp) call error('dralum','maxfp overflow',' ')
                 ifp=nbdpf
                 yield(:ninter,:nbfiss,ifp)=0.0
               endif
               do ifi=1,nbfiss
                 do i=1,ninter
                   yield(i,ifi,ifp)=yield(i,ifi,ifp)+yield(i,ifi,jfp)*prgar
                 enddo
               enddo
               idreac(2,ida)=ifp*100+5
             endif
           endif
         enddo
       enddo
       ! add decay heat to Draglib isotopes
       do jfath=1,nfath
         if(ipreac(jfath,jnd).eq.0) exit
         kt=mod(ipreac(jfath,jnd),100) ! reaction
         if(kt.eq.2) call error('dralum','fission forbidden',' ')
         knd=ipreac(jfath,jnd)/100 ! new father
         dener(kt,knd)=dener(kt,knd)+prate(jfath,jnd)*dener(1,jnd)
         ipreac(jfath,jnd)=0
         prate(jfath,jnd)=0.0
       enddo
       ymax=0.0
       if(mod(idreac(2,jnd),100).eq.5) then
         jfp=idreac(2,jnd)/100
         do ifi=1,nbfiss
           do i=1,ninter
             ymax=max(ymax,abs(yield(i,ifi,jfp)))
           enddo
         enddo
       endif
       dener(1,jnd)=0.0
       idreac(1,jnd)=0
       if(ymax.gt.0.0) then
         over_half=86400.0*ddeca(jnd)/(1.0e8*log(2.0))
         if(over_half.eq.0.0) then
           call dranam(mylist(jnd),hname)
           write(nsyso,'('' warning: isotope '',a8,'' is lumped and is'', &
           & '' stable. Max fission yield='',1p,e8.1,''%'')') hname,ymax*100.0
           if(ymax.gt.1.0E-2) call error('dralum','isotope '//hname// &
           & ' should not be lumped(1)',' ')
         else if((1.0/over_half.gt.40.0).and.(1.0/over_half.lt.999999.99)) then
           call dranam(mylist(jnd),hname)
           write(nsyso,'('' warning: isotope '',a8,'' is lumped and'', &
           & '' has a half-life of'',f10.2,'' days. Max fission yield='', &
           & 1p,e8.1,''%'')') hname,1.0/over_half,ymax*100.0
           if(ymax.gt.1.0E-2) call error('dralum','isotope '//hname// &
           & ' should not be lumped(2)',' ')
         else if(1.0/over_half.gt.40.0) then
           call dranam(mylist(jnd),hname)
           write(nsyso,'('' warning: isotope '',a8,'' is lumped and'', &
           & '' has a half-life of'',1p,e10.3,'' days. Max fission yi'', &
           & ''eld='',e8.1,''%'')') hname,1.0/over_half,ymax*100.0
           if(ymax.gt.1.0E-2) call error('dralum','isotope '//hname// &
           & ' should not be lumped(3)',' ')
         endif
       endif
       ddeca(jnd)=0.0
       mylist(jnd)=0
       50 continue
     enddo
   enddo
   write(nsyse,'('' ......... nlump='',i5)') nlump
   write(nsyso,'('' ......... nlump='',i5)') nlump
   if(nlump.gt.0) go to 40
   !
   ! if idecay=1, set decay heat of fission products to zero
   if(idecay.eq.1) then
     do iso=1,nbch
       ind=ipos(iso,1)
       ia=mod(mylist(ind),10000)/10
       if(ia.le.210) dener(1,ind)=0.0
     enddo
   endif
   !
   ! *write vectors 'PRODUCE-REAC' and 'PRODUCE-RATE' to the xsm file
   do iso=1,nbch
     do ifath=1,nfath
       jpreac(ifath,iso)=0
       rrate(ifath,iso)=0.0
     enddo
     ind=ipos(iso,1)
     nn=0
     do ifath=1,nfath
       if(ipreac(ifath,ind).ne.0) then
         do j=1,ifath-1
           if(ipreac(ifath,ind).eq.ipreac(j,ind)) then
             jnd1=ipreac(ifath,ind)/100
             jnd2=ipreac(j,ind)/100
             write(nsyso,'(/27h dralum: duplicate fathers:,2a4, &
             & 1x,2a4)') hiso(1,jnd1),hiso(2,jnd1),hiso(1,jnd2),hiso(2,jnd2)
             write(hname,'(2a4)') hiso(1,ind),hiso(2,ind)
             call error('dralum','duplicate fathers for '//hname,' ')
           endif
         enddo
         call draind(nbch,ipreac(ifath,ind)/100,ipos(1,1),jso)
         if(jso.le.0) then
           jnd=ipreac(ifath,ind)/100
           write(nsyso,'(/24h dralum: unknown father ,2a4,5h for , &
           & 2a4)') hiso(1,jnd),hiso(2,jnd),hiso(1,ind),hiso(2,ind)
         else
           nn=nn+1
           if(nn.gt.nfath) then
             write(text4,'(i4)') nn
             call error('dralum','nfath overflow nn='//text4,' ')
           endif
           jpreac(nn,iso)=100*jso+mod(ipreac(ifath,ind),100)
           rrate(nn,iso)=real(prate(ifath,ind))
         endif
       endif
     enddo
   enddo
   call xsmput(draglib,'PRODUCE-REAC',reshape(jpreac,(/ nfath*nbch /)))
   call xsmput(draglib,'PRODUCE-RATE',reshape(rrate,(/ nfath*nbch /)))
   !
   ! *write the isotope ascii names and header information on xsm file
   do i=1,nreac-1
     hreac(1,i)=reac(i)(1:4)
     hreac(2,i)=reac(i)(5:8)
   enddo
   call xsmput(draglib,'DEPLETE-IDEN',reshape(hreac,(/ 2*(nreac-1) /)))
   !
   ! *write the lumped fission yield matrix to the xsm file
   ibfp=0
   do iso=1,nbch
     ind=ipos(iso,1)
     ipos(iso,2)=0
     if(mod(idreac(2,ind),100).eq.5) then
       ibfp=ibfp+1
       if(ibfp.gt.nbfpch) call error('dralum','nbfpch overflow',' ')
       ipos(iso,2)=ibfp
       do ifi=1,nbfiss
         do i=1,ninter
           eyiel(i,ifi,ibfp)=real(yield(ninter-i+1,ifi,idreac(2,ind)/100))
         enddo
       enddo
     endif
   enddo
   if(ibfp.gt.0) then
     ! check lumped fission yield normalization
     if(iprint.gt.0) then
       write(nsyso,'(/49h dralum: check lumped fission yield normalization)')
       do iso=1,nbch
         ida=ipos(iso,1)
         if(mod(idreac(2,ida),100).ne.4) cycle
         ifi=idreac(2,ida)/100
         if(ifi.eq.0) cycle
         do i=1,ninter
          write(nsyso,'(i3,1x,2a4,f9.4)') i,hiso(:2,ida),sum(eyiel(i,ifi,:ibfp))
         enddo
       enddo
     endif
     call xsmput(draglib,'FISSIONYIELD',reshape(eyiel,(/ ninter*nbfiss*ibfp /)))
     allocate(eeee(ninter))
     do i=1,ninter
       eeee(i)=real(energy(ninter-i+1))
     enddo
     call xsmput(draglib,'PENERG-YIELD',eeee)
     deallocate(eeee)
   endif
   !
   ! *write vectors 'DEPLETE-REAC' and 'DEPLETE-ENER' to the xsm file
   do iso=1,nbch
     ind=ipos(iso,1)
     do i=1,nreac-1
       if(idreac(i,ind)/100.gt.0) then
         kreac=mod(idreac(i,ind),100)
         if((kreac.le.0).or.(kreac.gt.5)) then
           call error('dralum','invalid reaction',' ')
         endif
       endif
       if((i.eq.2).and.(mod(idreac(i,ind),100).eq.5)) then
         jdreac(i,iso)=ipos(iso,2)*100+5
       else
         jdreac(i,iso)=idreac(i,ind)
       endif
       eener(i,iso)=real(dener(i,ind))
     enddo
     if(dener(14,ind).ne.0.0) then
        do isoo=1,nbch
           eener(1,isoo)=0.0
           do iii=3,nreac-1
              eener(iii,isoo)=0.0
           enddo
          enddo
        eener(2,iso)=real(dener(14,ind))
     endif
   enddo
   call xsmput(draglib,'DEPLETE-REAC',reshape(jdreac,(/ (nreac-1)*nbch /)))
   call xsmput(draglib,'DEPLETE-ENER',reshape(eener,(/ (nreac-1)*nbch /)))
   !
   ! *write vectors 'CHARGEWEIGHT', 'DEPLETE-DECA', 'ISOTOPESDEPL'
   ! and 'STATE-VECTOR' to the xsm file
   nbheav=0
   do iso=1,nbch
     hhhh(1,iso)=hiso(1,ipos(iso,1))
     hhhh(2,iso)=hiso(2,ipos(iso,1))
     hhhh(3,iso)=hiso(3,ipos(iso,1))
     dddd(iso)=real(ddeca(ipos(iso,1)))
     ipos(iso,1)=mylist(ipos(iso,1))
     if(ipos(iso,1).ge.900000) nbheav=nbheav+1
   enddo
   call xsmput(draglib,'ISOTOPESDEPL',reshape(hhhh,(/ 3*nbch /)))
   call xsmput(draglib,'CHARGEWEIGHT',ipos(1:nbch,1))
   call xsmput(draglib,'DEPLETE-DECA',dddd)
   istate(:nstate)=0
   istate(1)=nbch
   istate(2)=nbfiss
   istate(3)=ibfp
   istate(4)=nbheav
   istate(5)=nbch-nbheav
   istate(8)=nreac-1
   istate(9)=nfath
   istate(10)=ninter
   call xsmput(draglib,'STATE-VECTOR',istate)
   deallocate(dddd,eyiel,eener,rrate)
   deallocate(hhhh,ipos,jdreac,jpreac)
   return
   end subroutine dralum
   !
   subroutine dranam(izae,hname)
   !-----------------------------------------------------------------
   !   utility to compose a material ascii (character*8) name or to
   !   recover the izae for an existing name.
   !-----------------------------------------------------------------
   use util   ! provides error
   integer izae,i,ia,iz,ifps,iof1,iof2,nia,niz
   character hname*8,hformat*12
   character(len=2), save, dimension(0:111) :: cs= &
    (/    'n ','H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
   & 'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti', &
   & 'V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se', &
   & 'Br','Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd', &
   & 'Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce', &
   & 'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
   & 'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb', &
   & 'Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U ','Np','Pu', &
   & 'Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg', &
   & 'Bh','Hs','Mt','Ds','Rg' /)
   character(len=2), save, dimension(0:3) :: cm= (/ '  ','m ','m2','m3' /)
   !
   iof1=0
   if(izae.eq.0) then
     ! **compute the izae index from name
     iz=0
     ifps=0
     do i=1,111
       if(hname(:2).eq.cs(i)) then
         iz=i
         iof1=3
         go to 10
       endif
     enddo
     do i=1,111
       if((hname(:1).eq.cs(i)(:1)).and.(cs(i)(2:).eq.' ')) then
         iz=i
         iof1=2
         go to 10
       endif
     enddo
     call error('dranam','invalid name:'//hname,' ')
     10 iof2=index(hname(iof1:),'m')
     if(iof2.gt.0) then
       do i=1,3
         if(hname(iof1+iof2-1:).eq.cm(i)) then
           ifps=i
           go to 20
         endif
       enddo
       call error('dranam','invalid suffix:'//hname,' ')
     else
       iof2=index(hname(iof1:),' ')
     endif
     20   read(hname(iof1:iof1+iof2-2),'(i3)') ia
     izae=iz*10000+ia*10+ifps
   else
     ! **compute the name from izae index
     iz=izae/10000
     ia=mod(izae,10000)/10
     ifps=mod(izae,10)
     if(iz.gt.111) call error('dranam','cs overflow',' ')
     if(ifps.gt.3) call error('dranam','cm overflow',' ')
     niz=2
     if(cs(iz)(2:).eq.' ') niz=1
     nia=3
     if(ia.lt.100) nia=2
     if(ia.lt.10) nia=1
     write(hformat,'(''(a'',i1,'',i'',i1,'',a2)'')') niz,nia
     write(hname,hformat) cs(iz)(:niz),ia,cm(ifps)
   endif
   return
   end subroutine dranam
   !
   real(kr) function draran(idum)
   !-------------------------------------------------------------------
   ! Random number generator.
   ! Set idum negative to reset the seed.
   !-------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::idum
   ! internals
   integer,parameter::m=714025
   integer,parameter::ia=1366
   integer,parameter::ic=150889
   real(kr),parameter::zero=0
   real(kr),parameter::one=1
   real(kr),parameter::rm=one/m
   integer::j,iy
   integer::ir(97)
   integer::iff=0
   save iff,iy,ir

   ! Do not allow zero
   100 continue
   if(idum.lt.0.or.iff.eq.0) then
      iff=1
      idum=mod(ic-idum,m)
      do j=1,97
         idum=mod(ia*idum+ic,m)
         ir(j)=idum
      enddo
      idum=mod(ia*idum+ic,m)
      iy=idum
   endif
   j=1+(97*iy)/m
   if(j.gt.97.or.j.lt.1) call error('draran','failed',' ')
   iy=ir(j)
   draran=iy*rm
   idum=mod(ia*idum+ic,m)
   ir(j)=idum
   if(draran.eq.zero) go to 100
   return
   end function draran
   !
   real(kr) function draurr(ityp,ibin,nbinpt,nunr,nbdil,lssf,ep,rand, &
   & eneurr,scr2a,scr2b)
   !-------------------------------------------------------------------
   ! Draw a random cross section in the unresolved energy range
   !-------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer :: ityp,ibin,nbinpt,nunr,nbdil,lssf
   real(kr) :: ep,rand,eneurr(nunr),scr2a(nbdil+(1+5*nbdil)*nunr), &
   & scr2b((1+6*nbinpt)*nunr)
   ! internals
   integer :: i,iunr,iof,inor
   real(kr) :: ww,siginf
   real(kr), allocatable, dimension(:) :: www,sig
   character(len=131) hsmg
   !
   draurr=0.0d0
   allocate(www(nbinpt),sig(nbinpt))
   iunr=0
   do i=nunr,1,-1
     if(ep.ge.eneurr(i)) then
       iunr=i
       go to 10
     else if((i.eq.1).and.(ep.ge.eneurr(1)-1.0e-2)) then
       iunr=i
       go to 10
     endif
   enddo
   write(hsmg,'(7henergy=,1p,e12.4,12h outside urr)') ep
   call error('draurr',hsmg,' ')
   10 siginf=scr2a(nbdil+(iunr-1)*(1+5*nbdil)+(ityp-1)*nbdil+2)
   iof=1+(iunr-1)*(1+6*nbinpt)
   www(:nbinpt)=scr2b(iof+1:iof+nbinpt)
   iof=iof+ityp*nbinpt
   sig(:nbinpt)=scr2b(iof+1:iof+nbinpt)
   ww=0.0d0
   do inor=1,nbinpt
     ww=ww+www(inor)
     if(rand.le.ww+1.0e-6) then
       if(lssf.eq.1) then
         draurr=sig(inor)*siginf
       else
         draurr=sig(inor)
       endif
       go to 20
     endif
   enddo
   write(hsmg,'(43hweight normalization issue in autolib group,i5)') ibin
   call error('draurr',hsmg,' ')
   20 deallocate(sig,www)
   return
   end function draurr
end module dragm
