module calenm
   ! provides subroutine calenr for NJOY2016
   use locale
   implicit none
   private
   public calenr

   ! units
   integer::nendf,npendf,nin,nout
   integer::nb
contains

   subroutine calenr
   !-----------------------------------------------------------------
   !
   ! Include CALENDF probability tables in the gendf tape.
   !
   ! Author: Alain Hebert
   !
   ! Copyright:
   !  Copyright (C) 2025 Polytechnique Montreal
   !  This library is free software; you can redistribute it and/or
   !  modify it under the terms of the BSD license.
   !
   !---input specifications (free format)----
   !
   ! card 1 tape units
   !   nendf   input endf unit
   !   npendf  input pendf unit
   !   nin     unit for input gendf tape
   !   nout    unit for output gendf tape
   !   ipreci  accuracy criterion for tables (=1/=2/=3/=4)
   !   iprint  long print option (0/1=minimum/maximum) (default=0)
   ! card 2
   !   matd    material to be processed
   !           matd=0 terminates calenr
   !   ntemp   no. of temperatures (default=1)
   ! card 3
   !   temp    temperatures in kelvin (including zero)
   ! card 4 energy limits for dilution-dependent xs data. This data
   !        is important to avoid self-shielding in very low and
   !        very high energy group xs and to reduce the size of the
   !        probability tables. (one card per material)
   !   eres0     lower limit of resolved resonance xs data (ev)
   !   eres1     upper limit of unresolved resonance xs data (ev)
   !
   !            repeat cards 2 to 4 for all materials desired
   !            matd=0/ terminates calenr run.
   !-----------------------------------------------------------------
   use mainio ! provides nsysi,contio,nsyso,nsyse
   use endf   ! provides endf routines and variables
   use util   ! provides timer,openz,repoz,error
   ! internals
   integer :: maxa,maxgr,maxnl,maxnz,maxnor,npar
   ! maxnor: maximum order of a CALENDF probability table
   ! npar: number of partial cross section types
   parameter (maxa=2000,maxgr=2000,maxnl=8,maxnz=30,maxnor=12,npar=6)
   integer::matd,itemp,ntemp,irflag,iprint,ner,nw,ng,ng0,ipreci,mtin,ipar, &
   & nl,nz
   logical lurr,lfind,exist
   real(kr)::time,eres0,eres1,urlimit,eh,aa(6)
   real(kr) ener(maxgr+1)
   real(kr),dimension(:),allocatable :: scr,temp
   real(kr),dimension(:,:,:),allocatable :: proref,flux,rv
   ! list of mt corresponding to partial (correlated) reactions
   integer, save, dimension(npar) :: malist= &
   & (/ 2,4,16,17,18,102 /) ! must be in increasing order
!
   call timer(time)
   write(nsyso,'(/'' calenr...produce a gendf format output tape'',24x,f8.1, &
   & ''s'')') time
   write(nsyse,'(/'' calenr...'',59x,f8.1,''s'')') time
   !
   !**read users input/output unit numbers and input data
   nendf=0
   npendf=0
   nin=0
   nout=0
   ipreci=4
   iprint=0
   read(nsysi,*) nendf,npendf,nin,nout,ipreci,iprint
   if(npendf.eq.0) call error('calenr','no pendf tape defined',' ')
   allocate(scr(maxa))
   !
   !**open the input/output tapes
   write(nsyso,'(/ &
   &  '' input endf unit ............................. '',i10/ &
   &  '' input pendf unit ............................ '',i10/ &
   &  '' input gendf unit ............................ '',i10/ &
   &  '' output gendf unit ........................... '',i10/ &
   &  '' accuracy criterion for tables (=1/=2/=3/=4).. '',i10/ &
   &  '' print option (0 min, 1 max) ................. '',i10)') &
   &  nendf,npendf,nin,nout,ipreci,iprint
   if((nin.lt.0.and.nout.gt.0).or.(nin.gt.0.and.nout.lt.0)) &
     & call error('calenr', &
     & 'mode conversion not allowed between nin and nout',' ')
   if(nendf.ne.0) call openz(nendf,0)
   if(npendf.ne.0) call openz(npendf,0)
   if(nin.ne.0) call openz(nin,0)
   if(nout.ne.0) call openz(nout,1)
   call repoz(nendf)
   call repoz(npendf)
   call repoz(nin)
   call repoz(nout)
   !
   !--loop over requested materials
   do
     matd=0
     iprint=1
     ntemp=1
     read(nsysi,*) matd,ntemp
     if(matd.eq.0) exit
     !
     !--process this material
     call timer(time)
     write(nsyso,'(/'' mat = '',i4,58x,f8.1,''s'')') matd,time
     write(nsyse,'(/'' mat = '',i4,58x,f8.1,''s'')') matd,time
     !
     !--copy through to desired material
     call tpidio(nin,nout,0,scr,nb,nw)
     10 call contio(nin,0,0,scr,nb,nw)
     if(math == -1) go to 120
     if(math == matd) go to 120
     if(math < 0) call error('calenr','desired material not on gendf tape',' ')
     call contio(0,nout,0,scr,nb,nw)
     call tomend(nin,nout,0,scr)
     go to 10
     !
     ! ***energy group information for the library
     120 call repoz(nin)
     call calhd(nin,matd,ng,maxgr,ener)
     !
     !--recover temperature information
     if(allocated(temp)) deallocate(temp)
     allocate(temp(ntemp))
     read(nsysi,*) (temp(itemp),itemp=1,ntemp)
     write(nsyso,'(&
       &'' temperatures ................................ '',1p,e10.3)') temp(1)
     if(ntemp.gt.1) write(nsyso,'(47x,1p,e10.3)') (temp(itemp),itemp=2,ntemp)
     urlimit=0.0
     if(nendf.ne.0) then
       ! ***recover the upper limit of the resolved energy domain.
       call repoz(nendf)
       call findf(matd,2,151,nendf)
       call contio(nendf,0,0,scr,nb,nw)
       call contio(nendf,0,0,scr,nb,nw)
       ner=n1h
       lurr=.false.
       eh=-1.0E10
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
         if(nw.gt.maxa) call error('calenr','endf input size exceeded(1)',' ')
         do while (nb.ne.0)
           call moreio(nendf,0,0,scr,nb,nw)
           if(nw.gt.maxa) call error('calenr','endf input size exceeded(2)',' ')
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
     ! ***recover energy limits for the complete self-shielding rang
     read(nsysi,*) eres0,eres1
     if(eres0.eq.0.0d0) eres0=ener(ng+1)
     if(eres1.eq.0.0d0) eres1=ener(1)
     !
     ! check if purr data is there
     irflag=0
     if(lurr) then
       math=matd
       lfind=.false.
       do while (.not.lfind)
         if(npendf.lt.0) then
           read(-npendf,end=20) math,mfh,mth,nb,nw
         else if(npendf.gt.0) then
           read(npendf,'(6e11.0,i4,i2,i3,i5)',end=20) aa,math,mfh,mth,nsp
         endif
         if(math.eq.0) go to 20
         lfind=(math.eq.matd).and.(mfh.eq.2).and.(mth.eq.153)
       enddo
       if(lfind) irflag=1
       20 call repoz(npendf)
     endif
     write(nsyso,'( &
     &  '' purr information flag ....................... '',i10)') irflag
     !
     ! recover infinite-dilution GENDF cross sections
     allocate(proref(npar+1,ng,ntemp))
     proref(:npar+1,:ng,:ntemp)=0.d0
     allocate(flux(maxgr,maxnl,maxnz),rv(maxgr,maxnl,maxnz))
     do itemp=1,ntemp
       do ipar=1,npar+1
         if(ipar.eq.1) then
           mtin=1 ! total
         else
           mtin=malist(max(1,ipar-1)) ! partial reaction
         endif
         call cfvect(maxgr,maxnl,maxnz,nin,matd,mtin,temp(itemp),ng0,nl,nz, &
         & rv,flux,exist)
         if(.not.exist) cycle
         if(ng0.ne.ng) call error('calenr','inconsistent ng',' ')
         proref(ipar,:ng,itemp)=rv(:ng,1,1)
       enddo
     enddo
     deallocate(rv,flux)
     !
     ! ***put calenf probability tables in output gendf tape for this material
     call repoz(nin)
     call findf(matd,1,451,nin)
     call calpt(npendf,matd,ng,ener,eres0,eres1,irflag,urlimit,ntemp,temp, &
     & maxnor,ipreci,npar,malist,proref,iprint)
     deallocate(proref)
   enddo ! material loop
   deallocate(temp)
   !
   ! **finished
   call atend(nout,0)
   call closz(nout)
   call closz(nin)
   call closz(npendf)
   call timer(time)
   write(nsyso,'(69x,f8.1,''s''/1x,7(''**********''),''*******'')') time
   deallocate(scr)
   return
   end subroutine calenr
   !
   subroutine calhd(ngn,matd,ng,maxgr,ener)
   !-----------------------------------------------------------------
   ! recover energy group data from gendf tape
   ! Group 1 is the most thermal group.
   !-----------------------------------------------------------------
   use endf   ! provides endf routines and variables
   use util   ! provides error
   integer :: maxa
   parameter(maxa=3000)
   integer ngn,matd,ng,maxgr
   real(kr) ener(maxgr+1)
   integer ig,loc,nw,nz
   real(kr),allocatable,dimension(:) :: scr
   !
   allocate(scr(maxa))
   call findf(matd,1,451,ngn)
   call contio(ngn,0,0,scr,nb,nw)
   nz=nint(scr(4))
   call listio(ngn,0,0,scr(1),nb,nw)
   loc=1+nw
   do while (nb.ne.0)
     if(loc+302.gt.maxa) call error('calhd','endf input size exceeded',' ')
     call moreio(ngn,0,0,scr(loc),nb,nw)
     loc=loc+nw
   enddo
   ng=nint(scr(3))
   if(ng.gt.maxgr) call error('calhd','maxgr overflow',' ')
   do ig=1,ng+1
     ener(ig)=scr(7+nz+ig)
   enddo
   deallocate(scr)
   return
   end subroutine calhd
   !
   subroutine calurr(maxnor,nbinpt,nunr,nbdil,npar,lssf,ep,eneurr,scr2a, &
   & scr2b,factor,momt,momp,ndil,dilut,proref,sefr)
   !-------------------------------------------------------------------
   ! Compute cross section moments in the unresolved resonance rang
   !-------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer :: maxnor,nbinpt,nunr,nbdil,npar,lssf,ndil
   real(kr) :: ep,eneurr(nunr),scr2a(nbdil+(1+5*nbdil)*nunr), &
   & scr2b((1+6*nbinpt)*nunr),factor,momt(2*maxnor),momp(maxnor,npar), &
   & dilut(ndil),proref(npar+1),sefr(npar+2,ndil)
   ! internals
   integer :: i,iunr,iof,inor,jnor,idil,ipar,ityp
   real(kr) :: tt,sigt,sigx,bond,prodif
   integer, save, dimension(6) :: matype= &
   & (/ 1,0,0,0,2,3 /) ! purr/malist correspondance
   real(kr), allocatable, dimension(:) :: siginf,www
   real(kr), allocatable, dimension(:,:) :: sig
   character(len=131) hsmg
   !
   allocate(siginf(npar+1),www(nbinpt),sig(nbinpt,npar+1))
   iunr=0
   if(ep.ge.eneurr(nunr)) go to 20 ! continuum
   do i=nunr,1,-1
     if(ep.ge.eneurr(i)) then
       iunr=i
       go to 10
     else if((i.eq.1).and.(ep.ge.eneurr(1)-1.0e-2)) then
       iunr=i
       go to 10
     endif
   enddo
   call error('calurr',hsmg,' ')
   !
   ! process URR
   10 do ityp=1,4 ! assume 4 types of PURR basepoints
     siginf(ityp)=scr2a(nbdil+(iunr-1)*(1+5*nbdil)+(ityp-1)*nbdil+2)
   enddo
   iof=1+(iunr-1)*(1+6*nbinpt)
   www(:nbinpt)=scr2b(iof+1:iof+nbinpt)
   do ityp=1,4
     iof=iof+nbinpt
     sig(:nbinpt,ityp)=scr2b(iof+1:iof+nbinpt)
   enddo
   do inor=1,nbinpt
     if(lssf.eq.1) then
       sigt=sig(inor,1)*siginf(1)
     else
       sigt=sig(inor,1)
     endif
     tt=www(inor)
     do jnor=maxnor,2*maxnor
       momt(jnor)=momt(jnor)+factor*tt
       tt=tt*sigt
     enddo
     tt=www(inor)/sigt
     do jnor=maxnor-1,1,-1
       momt(jnor)=momt(jnor)+factor*tt
       tt=tt/sigt
     enddo
     do idil=1,ndil
       bond=factor*dilut(idil)/(sigt+dilut(idil))
       sefr(1,idil)=sefr(1,idil)+www(inor)*bond
       sefr(2,idil)=sefr(2,idil)+www(inor)*sigt*bond
     enddo
     do ipar=1,npar
       ! sigx is the base point of a partial reaction
       ityp=matype(ipar)
       if(ityp.gt.0) then
         if(lssf.eq.1) then
           sigx=sig(inor,ityp+1)*siginf(ityp+1)
         else
           sigx=sig(inor,ityp+1)
         endif
       else
         ! purr doesn't process malist(ipar); use an approximation in the URR
         if(lssf.eq.1) then
           sigx=sig(inor,1)*siginf(1)-sig(inor,2)*siginf(2)-sig(inor,3)* &
           & siginf(3)-sig(inor,4)*siginf(4)
         else
           sigx=sig(inor,1)-sig(inor,2)-sig(inor,3)-sig(inor,4)
         endif
         ! 1: total; 2: elastic; 6: fission; 7: capture
         prodif=proref(1)-proref(2)-proref(6)-proref(7)
         if(prodif.gt.1.0E-4) then
           sigx=sigx*proref(ipar+1)/prodif
         else
           sigx=proref(ipar+1)
         endif
       endif
       tt=www(inor)*sigx
       do jnor=maxnor/2+1,maxnor
         momp(jnor,ipar)=momp(jnor,ipar)+factor*tt
         tt=tt*sigt
       enddo ! inor
       tt=www(inor)*sigx/sigt
       do jnor=maxnor/2,1,-1
         momp(jnor,ipar)=momp(jnor,ipar)+factor*tt
         tt=tt/sigt
       enddo ! jnor
       do idil=1,ndil
         bond=factor*dilut(idil)/(sigt+dilut(idil))
         sefr(ipar+2,idil)=sefr(ipar+2,idil)+www(inor)*sigx*bond
       enddo
     enddo ! ipar
   enddo ! inor
   go to 30
   !
   ! process continuum
   20 sigt=proref(1) ! sigt is the total cross section
   tt=1.0
   do jnor=maxnor,2*maxnor
     momt(jnor)=momt(jnor)+factor*tt
     tt=tt*sigt
   enddo
   tt=1.0/sigt
   do jnor=maxnor-1,1,-1
     momt(jnor)=momt(jnor)+factor*tt
     tt=tt/sigt
   enddo
   do idil=1,ndil
     bond=factor*dilut(idil)/(sigt+dilut(idil))
     sefr(1,idil)=bond
     sefr(2,idil)=sigt*bond
   enddo
   do ipar=1,npar
     sigx=proref(ipar+1) ! sigx is the base point cross section
     tt=sigx
     do jnor=maxnor/2+1,maxnor
       momp(jnor,ipar)=momp(jnor,ipar)+factor*tt
       tt=tt*sigt
     enddo ! inor
     tt=sigx/sigt
     do jnor=maxnor/2,1,-1
       momp(jnor,ipar)=momp(jnor,ipar)+factor*tt
       tt=tt/sigt
     enddo ! jnor
     do idil=1,ndil
       bond=factor*dilut(idil)/(sigt+dilut(idil))
       sefr(ipar+2,idil)=sigx*bond
     enddo
   enddo
   ! deallocate scratch memory
   30 deallocate(sig,www,siginf)
   return
   end subroutine calurr
   !
   subroutine calpt(npendf,matno,ng,ener,eres0,eres1,irflag,urlimit, &
   & ntemp,temp,maxnor,ipreci,npar,malist,proref,iprint)
   !-----------------------------------------------------------------
   !  put CALENDF probability tables on GENDF tape for material matno
   !-----------------------------------------------------------------
   use mainio ! provides nsysi,contio,nsyso,nsyse
   use endf   ! provides endf routines and variables
   use util   ! provides openz,repoz,error
   ! externals
   integer :: npendf,matno,ng,irflag,ntemp,maxnor,ipreci,npar,malist(npar), &
   & iprint
   real(kr) ener(ng+1),eres0,eres1,urlimit,temp(ntemp),proref(npar+1,ng,ntemp)
   ! internals
   integer :: maxa,maxtmp,maxbin,ndil
   parameter (maxa=3000,maxtmp=10,maxbin=1000000,ndil=5)
   integer i,j,ig,itm,mt,irflag_2,lssf,nbdil,ndata,nunr,loc,nw,nbinpt,igres0, &
   & igres1,idis,iener,nener,inor,idil,ipar,ibase,lim,nl,nz,icount,igmax
   character(len=131) :: hsmg
   real(kr) temps,aa(6),deltau,em,emlog,ep,eplog,enext,tm,tp,sm,sp,tt,sigt, &
   & sigx,errbst,factor,za,awr
   logical lnoraj,lfind
   integer, allocatable, dimension(:) :: npanel,nor
   real(kr), allocatable, dimension(:) :: scr,scr1,scr2a,scr2b,eneurr
   real(kr), allocatable, dimension(:,:) :: momt,ekep
   real(kr), allocatable, dimension(:,:,:) :: momp,sefr,prosig
   logical, allocatable, dimension(:) :: lfindv
   real(kr),dimension(ndil),parameter :: dilut = &
     & (/ 1.e10, 1.e5, 1.e4, 1.e3, 1.e2 /)
   allocate(scr(maxa),nor(ng),lfindv(npar))
   nor(:ng)=0
   !
   ! ***compute igres0 and igres1
   igres0=1
   igres1=ng
   do ig=ng,2,-1
     if(eres1.le.0.999d0*ener(ig+1)) igres1=ig-1
   enddo
   do ig=1,ng
     if(eres0.ge.1.001d0*ener(ig)) igres0=ig+1
   enddo
   write(nsyso,'(/ &
   &  '' first group of resolved resonance xs data ... '',i10/ &
   &  '' last group of unresolved resonance xs data .. '',i10)') &
   &  igres0,igres1
   !
   !**recover the temperature values.
   lnoraj=.true.
   if(ntemp.gt.maxtmp) call error('calpt','ntemp overflow',' ')
   do itm=1,ntemp
     call repoz(npendf)
     call skiprz(npendf,1)
     10 lfind=.false.
     do while(.not.lfind)
       if(npendf.lt.0) then
         read(-npendf,end=100) math,mfh,mth,nb,nw,aa
       else if(npendf.gt.0) then
         read(npendf,'(6e11.0,i4,i2,i3,i5)',end=100) aa,math,mfh,mth,nsp
       endif
       if(math.eq.-1) go to 100
       lfind=(math.eq.matno).and.(mfh.eq.1).and.(mth.eq.451)
     enddo
     if(.not.lfind) go to 100
     za=aa(1)
     awr=aa(2)
     call contio(npendf,0,0,scr,nb,nw)
     call contio(npendf,0,0,scr,nb,nw)
     call contio(npendf,0,0,scr,nb,nw)
     call skiprz(npendf,1) ! skip the comment line
     if(abs(c1h-temp(itm)).gt.1.0e-3*c1h) then
       call tomend(npendf,0,0,scr)
       go to 10
     endif
     temps=c1h
     if(iprint.gt.0) then
       write(nsyso,'(/28h calpt: process temperature=,1p,e12.4)') temps
     endif
     allocate(scr1(maxa))
     call tosend(nin,nout,0,scr)
     !
     !**recover the purr data (probability tables in the UR domain).
     lssf=0
     irflag_2=irflag
     nunr=0
     nbdil=0
     if(irflag_2.eq.1) then
       math=matno
       call findf(math,2,152,npendf)
       if(math.eq.-1) then
         irflag_2=0
         go to 40
       endif
       call contio(npendf,0,0,scr1,nb,nw)
       lssf=l1h
       call listio(npendf,0,0,scr1,nb,nw)
       if(nw.gt.maxa) call error('calpt','endf input size exceeded(1)',' ')
       nbdil=l2h
       ndata=n1h
       nunr=n2h
       allocate(scr2a(ndata))
       scr2a(:nw-6)=scr1(7:nw)
       loc=1+nw-6
       do while (nb.ne.0)
         call moreio(npendf,0,0,scr2a(loc),nb,nw)
         loc=loc+nw
         if(loc-1.gt.n1h) call error('calpt','endf input size exceeded(2)',' ')
       enddo
       !
       call findf(matno,2,153,npendf)
       call contio(npendf,0,0,scr1,nb,nw)
       nbinpt=n2h
       call listio(npendf,0,0,scr1,nb,nw)
       if(nw.gt.maxa) call error('calpt','endf input size exceeded(3)',' ')
       ndata=n1h
       nunr=n2h
       allocate(scr2b(ndata))
       scr2b(:nw-6)=scr1(7:nw)
       loc=1+nw-6
       do while (nb.ne.0)
         call moreio(npendf,0,0,scr2b(loc),nb,nw)
         loc=loc+nw
         if(loc-1.gt.n1h) call error('calpt','endf input size exceeded(4)',' ')
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
     40 allocate(prosig(maxnor,2+npar,ng))
     allocate(momt(2*maxnor,ng),momp(maxnor,npar,ng),sefr(npar+2,ndil,ng), &
      & npanel(ng))
     call findf(matno,3,1,npendf)
     call contio(npendf,0,0,scr1,nb,nw)
     ep=0.0d0
     call gety1(ep,enext,idis,tp,npendf,scr1)
     em=eres0
     call gety1(em,enext,idis,tm,npendf,scr1)
     allocate(ekep(2,maxbin)) ! allocate the big array containing sigt values
     nener=1
     ekep(1,1)=em ; ekep(2,1)=tm
     ig=1
     momt(:2*maxnor,:ng)=0.0d0
     momp(:maxnor,:npar,:ng)=0.0d0
     sefr(:npar+2,:ndil,:ng)=0.0d0
     npanel(:ng)=0
     emlog=log(eres1/em)
     do while(ep*(1.0d0+1.0d-6).lt.eres1)
       if(ig.gt.ng) call error('calpt','invalid index(1)',' ')
       ep=min(enext,ener(ig+1))
       eplog=log(eres1/ep)
       call gety1(ep,enext,idis,tp,npendf,scr1)
       if((irflag_2.eq.0).or.(ep.le.urlimit)) then
         npanel(ig)=npanel(ig)+1
         nener=nener+1
         if(nener.gt.maxbin) then
           write(hsmg,'(24hmaxbin overflow in group,i5)') ig
           call error('calpt',hsmg,' ')
         endif
         ekep(1,nener)=ep ; ekep(2,nener)=tp
         sigt=max(0.001,0.5d0*(tm+tp))
         tt=1.0d0
         do inor=maxnor,2*maxnor
           momt(inor,ig)=momt(inor,ig)+(emlog-eplog)*tt
           tt=tt*sigt
         enddo
         tt=1.0d0/sigt
         do inor=maxnor-1,1,-1
           momt(inor,ig)=momt(inor,ig)+(emlog-eplog)*tt
           tt=tt/sigt
         enddo
         do idil=1,ndil
           factor=(emlog-eplog)*dilut(idil)/(sigt+dilut(idil))
           sefr(1,idil,ig)=sefr(1,idil,ig)+factor
           sefr(2,idil,ig)=sefr(2,idil,ig)+factor*sigt
         enddo
       else
         npanel(ig)=-1
         factor=emlog-eplog
         call calurr(maxnor,nbinpt,nunr,nbdil,npar,lssf,ep,eneurr,scr2a, &
         & scr2b,factor,momt(1,ig),momp(1,1,ig),ndil,dilut,proref(1,ig,itm), &
         & sefr(1,1,ig))
       endif
       if(ep.eq.ener(ig+1)) ig=ig+1
       em=ep ; tm=tp ; sm=sp
       emlog=eplog
     enddo
     !
     lfindv(:npar)=.false.
     do ipar=1,npar
       !**recover the pointwise cross sections fot malist(ipar) (if exists).
       mt=malist(ipar)
       icount=0
       do while (.not.lfindv(ipar))
         if(npendf.lt.0) then
           read(-npendf,end=50) math,mfh,mth,nb,nw
         else if(npendf.gt.0) then
           read(npendf,'(6e11.0,i4,i2,i3,i5)',end=50) aa,math,mfh,mth,nsp
         endif
         icount=icount+1
         if(math.eq.0) then
           call skiprz(npendf,-icount)
           go to 50
         endif
         lfindv(ipar)=(math.eq.matno).and.(mfh.eq.3).and.(mth.eq.mt)
       enddo
       if(lfindv(ipar)) then
         call findf(matno,3,mt,npendf)
         call contio(npendf,0,0,scr1,nb,nw)
         ep=0.0d0
         call gety1(ep,enext,idis,tp,npendf,scr1)
         em=eres0
         call gety1(em,enext,idis,tm,npendf,scr1)
         iener=1
         sm=ekep(2,1)
         ig=1
         emlog=log(eres1/em)
         do while(ep*(1.0d0+1.0d-6).lt.eres1)
           if(ig.gt.ng) then
             write(hsmg,'(21hinvalid index for mt=,i3)') mt
             call error('calpt',hsmg,' ')
           endif
           if((irflag_2.eq.0).or.(ep.le.urlimit)) then
             iener=iener+1
             if(iener.gt.nener) exit
             ep=ekep(1,iener)
             eplog=log(eres1/ep)
             call gety1(ep,enext,idis,tp,npendf,scr1)
             sp=ekep(2,iener)
             sigt=max(0.001,0.5d0*(sm+sp))
             sigx=0.5d0*(tm+tp)
             tt=sigx
             do inor=maxnor/2+1,maxnor
               momp(inor,ipar,ig)=momp(inor,ipar,ig)+(emlog-eplog)*tt
               tt=tt*sigt
             enddo
             tt=sigx/sigt
             do inor=maxnor/2,1,-1
               momp(inor,ipar,ig)=momp(inor,ipar,ig)+(emlog-eplog)*tt
               tt=tt/sigt
             enddo
             do idil=1,ndil
               factor=(emlog-eplog)*dilut(idil)/(sigt+dilut(idil))
               sefr(ipar+2,idil,ig)=sefr(ipar+2,idil,ig)+factor*sigx
             enddo
           endif
           if(ep.eq.ener(ig+1)) ig=ig+1
           em=ep ; tm=tp ; sm=sp
           emlog=eplog
         enddo
       endif
       50 cycle
     enddo ! ipar
     deallocate(ekep) ! deallocate the big array containing sigt values
     if(irflag_2.eq.1) deallocate(scr2b,scr2a,eneurr)
     !
     ! compute calendf probability tables
     if(iprint.gt.0) then
       write(nsyso,'(/41h calpt: calendf information for material=,i8, &
       & 13h temperature=,1p,e12.4/2x,5hgroup,6x,6henergy,3x,6hpanels,3x, &
       & 5horder,7x,5herror)') matno,temps
     endif
     do ig=igres0,igres1
       deltau=log(ener(ig+1)/ener(ig))
       do inor=1,2*maxnor
         momt(inor,ig)=momt(inor,ig)/deltau
       enddo
       if(npar.gt.0) then
         do inor=1,maxnor
           do ipar=1,npar
             momp(inor,ipar,ig)=momp(inor,ipar,ig)/deltau
           enddo
         enddo
       endif
       do idil=1,ndil
         do ipar=2,npar+2
           sefr(ipar,idil,ig)=sefr(ipar,idil,ig)/sefr(1,idil,ig)
         enddo
         sefr(1,idil,ig)=sefr(1,idil,ig)/deltau
       enddo
       call calcat(maxnor,npar,ndil,momt(1,ig),momp(1,1,ig),ipreci,lnoraj, &
       & dilut,sefr(1,1,ig),nor(ig),prosig(1,1,ig),errbst)
       if(iprint.gt.0) then
         write(nsyso,'(1x,i6,1p,e12.4,i9,i8,e12.4)') ig,ener(ig),npanel(ig), &
         & nor(ig),errbst
       endif
     enddo
     lnoraj=.false.
     deallocate(npanel,sefr,momp,momt)
     !
     ! normalize probability tables with respect to unit nin information
     do ipar=2,npar+2
       if(ipar.eq.2) then
         mth=1 ! total
       else
         mth=malist(max(1,ipar-2)) ! partial reaction
         if(.not.lfindv(max(1,ipar-2))) cycle
       endif
       do ig=igres0,igres1
         if(proref(ipar-1,ig,itm).ne.0.0) go to 60
       enddo
       cycle
       60 if(iprint.gt.0) then
         write(nsyso,'(/48h calpt: normalization of probability tables for , &
         & 3hmt=,i4/13x,5hgendf,5x,7hcalendf,6x,6hfactor,3x,5horder)') mth
       endif
       do ig=igres0,igres1
         sigt=0.d0
         do inor=1,nor(ig)
           sigt=sigt+prosig(inor,1,ig)*prosig(inor,ipar,ig)
         enddo
         if(sigt.ne.0.d0) then
           factor=proref(ipar-1,ig,itm)/sigt
         else
           factor=1.d0
         endif
         prosig(:nor(ig),ipar,ig)=prosig(:nor(ig),ipar,ig)*factor
         if((iprint.gt.0).and.(proref(ipar-1,ig,itm).ne.0.0)) then
           write(nsyso,'(1x,i5,1p,3e12.4,i8,e12.4)') ig,proref(ipar-1,ig,itm), &
           & sigt,factor,nor(ig)
         endif
       enddo
     enddo
     !
     ! store probability tables on unit nout
     if(nout.ne.0) then
       call tomend(nin,nout,0,scr)
       call skiprz(nout,-1)
       mfh=50 ; nz=1 ; nl=1
       do ipar=2,npar+2
         if(ipar.eq.2) then
           mth=1 ! total
         else
           mth=malist(max(1,ipar-2))
         endif
         igmax=0
         do ig=igres0,igres1
           if(proref(ipar-1,ig,itm).ne.0.0) igmax=ig
         enddo
         if(igmax.eq.0) cycle
         scr(1)=za
         scr(2)=awr
         scr(3)=1
         scr(4)=1
         scr(5)=0
         scr(6)=ng
         call contio(0,nout,0,scr,nb,nw)
         do ig=igres0,igres1
           if(ipar.gt.2) then
             if(proref(ipar-1,ig,itm).eq.0.0) cycle
           endif
           nw=6
           lim=2*nl*nz*nor(ig)
           scr(1)=temps
           scr(2)=0
           scr(3)=2
           scr(4)=nor(ig)
           scr(5)=lim
           scr(6)=ig
           j=6
           ibase=6
           do inor=1,2*nor(ig)
             j=j+1
             if (inor.le.nor(ig)) then
                scr(j)=prosig(inor,1,ig) ! weight
             else
                scr(j)=prosig(inor-nor(ig),ipar,ig) ! base point
             endif
             if((j.ge.npage+ibase).or.(j-6.eq.lim)) then
               if(ibase.ne.0) then
                 call listio(0,nout,0,scr,nb,j)
                 ibase=0
                 j=0
               else
                 call moreio(0,nout,0,scr,nb,nw)
                 j=0
               endif
             endif
           enddo ! inor
         enddo ! ig
         if(igmax.lt.ng) then
           ! add one group with nor=1 to make ecco happy
           scr(1)=temps
           scr(2)=0
           scr(3)=2
           scr(4)=1
           scr(5)=2
           scr(6)=ng
           scr(7)=1.0
           scr(8)=proref(ipar-1,ng,itm)
           j=8
           call listio(0,nout,0,scr,nb,j)
         endif
         call asend(nout,0)
       enddo ! ipar
       call amend(nout,0)
     endif
     deallocate(prosig)
     call tomend(npendf,0,0,scr)
     deallocate(scr1)
   enddo ! itm
   100 deallocate(lfindv,nor,scr)
   return
   end subroutine calpt
   !
   subroutine alprtb(maxnor,nor,iini,demt,ier,weight,basept)
   !
   !-----------------------------------------------------------------------
   !
   !purpose:
   ! compute a probability table preserving 2*nor moments of a function
   ! using the modified Ribon approach.
   !
   !author(s): A. Hebert
   !
   !parameters: input
   ! maxnor  maximum order of a probability table.
   ! nor     the number of moments to preserve is 2*nor.
   ! iini    minimum order of the moment we want to preserve. we must have
   !         2-2*nor <= iini <= 0 (order 0 and 1 moments are always preserved).
   ! demt    moments.
   !
   !parameters: output
   ! ier     error flag (=0/=1 success/failure of the algorithm).
   ! weight  weights of the probability table.
   ! basept  base points of the probability table.
   !
   !-----------------------------------------------------------------------
   !
   use util   ! provides error
   implicit real(kr)(a-h,o-z)
   !----
   ! subroutine arguments
   !----
   integer maxnor,nor,iini,ier
   real(kr) demt(iini:2*nor+iini-1)
   real(kr) weight(nor),basept(nor)
   !----
   ! local variables
   !----
   integer i,j,k,ior,jor,j0
   real(kr) ds(maxnor+1,maxnor+1),dda(0:maxnor),dd,dsigx
   complex*16 roots(maxnor),ccc,dcc,xcc
   complex cgar
   logical lfail
   !
   if(nor.gt.maxnor) call error('alprtb','storage overflow',' ')
   if(nor.le.0) call error('alprtb','negative or zero value of nor',' ')
   if((2-2*nor.gt.iini).or.(iini.gt.0)) call error('alprtb','inconsistent va' &
     & //'lue of iini.',' ')
   !
   ! build the matrix.
   do 15 ior=1,nor
   ds(ior,nor+1)=-demt(nor+ior+iini-1)
   do 10 jor=1,ior
   ds(ior,jor)=demt(ior+jor+iini-2)
   ds(jor,ior)=demt(ior+jor+iini-2)
   10 continue
   15 continue
   !
   ! l-d-l(t) factorization of the matrix.
   do 40 i=1,nor
   do 30 j=1,i-1
   ds(j,i)=ds(i,j)
   do 20 k=1,j-1
   ds(j,i)=ds(j,i)-ds(k,i)*ds(j,k)
   20 continue
   ds(i,j)=ds(j,i)*ds(j,j)
   ds(i,i)=ds(i,i)-ds(j,i)*ds(i,j)
   30 continue
   if(ds(i,i).eq.0.d0) then
      ier=1
      return
   endif
   ds(i,i)=1.d0/ds(i,i)
   40 continue
   !
   ! solution of the factorized system to obtain the denominator of the
   ! Pade approximation.
   do 55 i=1,nor
   do 50 k=1,i-1
   ds(i,nor+1)=ds(i,nor+1)-ds(i,k)*ds(k,nor+1)
   50 continue
   55 continue
   do 60 i=1,nor
   ds(i,nor+1)=ds(i,nor+1)*ds(i,i)
   60 continue
   do 71 i=nor,1,-1
   do 70 k=i+1,nor
   ds(i,nor+1)=ds(i,nor+1)-ds(k,i)*ds(k,nor+1)
   70 continue
   71 continue
   ds(nor+1,nor+1)=1.0d0
   !
   ! compute the base points as the roots of the denominator.
   call alroot(ds(1,nor+1),nor,roots,lfail)
   if(lfail) call error('alprtb','polynomial root finding failure',' ')
   do 80 i=1,nor
   !
   ! Newton improvement of the roots.
   ccc=0.0d0
   xcc=1.0d0
   do 74 j=0,nor
   ccc=ccc+ds(j+1,nor+1)*xcc
   xcc=xcc*roots(i)
   74 continue
   dcc=0.0d0
   xcc=1.0d0
   do 75 j=1,nor
   dcc=dcc+ds(j+1,nor+1)*xcc*real(j)
   xcc=xcc*roots(i)
   75 continue
   roots(i)=roots(i)-ccc/dcc
   !
   cgar=cmplx(roots(i))
   if(abs(aimag(cgar)).gt.1.0e-4*abs(real(cgar))) then
      ier=1
      return
   else
      basept(i)=real(cmplx(roots(i)))
   endif
   80 continue
   !
   ! compute the weights.
   do 130 i=1,nor
   dsigx=dble(roots(i))
   dda(0)=1.0d0
   j0=0
   do 100 j=1,nor
   if(j.eq.i) go to 100
   j0=j0+1
   dda(j0)=dda(j0-1)
   do 90 k=1,j0-1
   dda(j0-k)=dda(j0-k-1)-dda(j0-k)*dble(roots(j))
   90 continue
   dda(0)=-dda(0)*dble(roots(j))
   100 continue
   dd=0.0d0
   do 110 j=0,nor-1
   dd=dd+dda(j)*demt((iini-1)/2+j)
   110 continue
   do 120 j=1,nor
   if(j.ne.i) dd=dd/(dble(roots(j))-dsigx)
   120 continue
   weight(i)=real(((-1.0d0)**(nor-1))*dd*dsigx**((1-iini)/2))
   130 continue
   !
   ! test the consistency of the solution.
   do 140 i=1,nor
   if((weight(i).le.0.0).or.(basept(i).le.0.0)) then
      ier=1
      return
   endif
   140 continue
   ier=0
   return
   end subroutine alprtb
   !
   subroutine alroot(a,m,roots,lfail)
   !
   !-----------------------------------------------------------------------
   !
   !purpose:
   ! find the roots of a polynomial.
   !
   !author(s): A. Hebert
   !
   !parameters: input
   ! a    polynomial coefficients.
   ! m    polynomial order.
   !
   !parameters: output
   ! roots   complex roots.
   ! lfail   flag set to .true. in case of failure.
   !
   !-----------------------------------------------------------------------
   !
   use util   ! provides error
   implicit real(kr)(a-h,o-z)
   !----
   !  subroutine arguments
   !----
   integer m
   real(kr) a(m+1)
   complex*16 roots(m)
   logical lfail
   !----
   !  local variables
   !----
   integer i,j,jj,maxm,its
   complex*16 aaa,bbb,sq1,test,sqrtm3
   parameter (sqrtm3=(0.0,1.73205080756888))
   real(kr) eps
   parameter (eps=1.d-6,maxm=101)
   complex*16 ad(maxm),x,b,c
   !
   lfail=.false.
   if(m+1.gt.maxm) call error('alroot','insufficient storage',' ')
   if(a(m+1).eq.0.0d0) call error('alroot','invalid coefficient',' ')
   if(m.eq.1) then
      roots(1)=-a(1)/a(2)
   else if(m.eq.2) then
      cq=a(2)/a(3)
      dq=a(1)/a(3)
      aaa=cq*cq-4.0d0*dq
      aaa=sqrt(aaa)
      roots(1)=-0.5d0*(cq+aaa)
      roots(2)=-0.5d0*(cq-aaa)
   else if(m.eq.3) then
      if(a(1).eq.0.0) then
      cq=a(3)/a(4)
      dq=a(2)/a(4)
      aaa=cq*cq-4.0d0*dq
      aaa=sqrt(aaa)
      roots(1)=0.0
      roots(2)=-0.5d0*(cq+aaa)
      roots(3)=-0.5d0*(cq-aaa)
      else
      bq=a(3)/a(4)
      cq=a(2)/a(4)
      dq=a(1)/a(4)
      aa=(3.0d0*cq-bq**2)/3.0d0
      bb=(2.0d0*bq**3-9.0d0*bq*cq+27.0d0*dq)/27.0d0
      sq1=bb**2/4.0d0+aa**3/27.0d0
      test=bb/2.0d0-sqrt(sq1)
      if(dble(test).eq.0.0) then
         aaa=0.0d0
      else if(dble(test).gt.0.0) then
         aaa=-(test)**(1.0d0/3.0d0)
      else
         aaa=(-test)**(1.0d0/3.0d0)
      endif
      test=bb/2.0d0+sqrt(sq1)
      if(dble(test).eq.0.0) then
         bbb=0.0d0
      else if(dble(test).gt.0.0) then
         bbb=-(test)**(1.0d0/3.0d0)
      else
         bbb=(-test)**(1.0d0/3.0d0)
      endif
      roots(1)=aaa+bbb-bq/3.0d0
      roots(2)=-(aaa+bbb)/2.0d0+(aaa-bbb)*sqrtm3/2.0d0-bq/3.0d0
      roots(3)=-(aaa+bbb)/2.0d0-(aaa-bbb)*sqrtm3/2.0d0-bq/3.0d0
      endif
   else if(m.eq.4) then
      call alquar(a,roots)
   else
      do j=1,m+1
        ad(j)=cmplx(a(j),0.d0,kind=kind(ad))
      enddo
      do j=m,1,-1
        x=cmplx(0.d0,0.d0,kind=kind(x))
        call alguer(ad,j,x,its,lfail)
        if(lfail) return
        if(abs(dimag(x)).le.2.d0*eps**2*abs(dble(x))) &
         & x=cmplx(dble(x),0.d0,kind=kind(x))
        roots(j)=x
        b=ad(j+1)
        do jj=j,1,-1
          c=ad(jj)
          ad(jj)=b
          b=x*b+c
        enddo
      enddo ! j
      do j=1,m+1
        ad(j)=cmplx(a(j),0.d0,kind=kind(ad))
      enddo
      do j=1,m
        call alguer(ad,m,roots(j),its,lfail)
        if(lfail) return
      enddo
   endif
   !
   do 70 j=2,m
   x=roots(j)
   do 50 i=j-1,1,-1
   if(dble(roots(i)).le.dble(x)) goto 60
   roots(i+1)=roots(i)
   50 continue
   i=0
   60 roots(i+1)=x
   70 continue
   return
   end subroutine alroot
   !
   subroutine alquar(a,roots)
   !
   !-----------------------------------------------------------------------
   !
   !purpose:
   ! compute the roots of the real quartic polynomial defined as
   ! a(1)+a(2)*z + ... + a(5)*z**4.
   ! note: it is assumed that a(5) is non-zero. no test is made here.
   !
   !author(s): A. H. Morris, W. L. Davis, A. Miller, and R. L. Carmichael
   !
   !parameters: input
   ! a   polynomial coefficients
   !
   !parameters: output
   ! roots  complex roots
   !
   !-----------------------------------------------------------------------
   !
   implicit real(kr)(a-h,o-z)
   !----
   !  subroutine arguments
   !----
   real(kr), intent(in) :: a(5)
   complex*16, intent(out) :: roots(4)
   !----
   !  local variables
   !----
   integer :: i, j
   integer,dimension(1) :: k
   real(kr) :: b, b2, c, d, e, h, p, q, r, t, work, u, v, v1, v2, x, xx(3),  &
     & y, bq, cq, dq, aa, bb
   complex*16 :: w, aaa, bbb, sq1, test, sqrtm3
   real(kr),dimension(4) :: temp
   parameter (sqrtm3=(0.0,1.73205080756888))
   !
   if(a(1)==0.0) then
     if(a(2).eq.0.0) then
        cq=a(4)/a(5)
        dq=a(3)/a(5)
        aaa=cq*cq-4.0d0*dq
        aaa=sqrt(aaa)
        roots(1)=0.0
        roots(2)=0.0
        roots(3)=-0.5d0*(cq+aaa)
        roots(4)=-0.5d0*(cq-aaa)
     else
        bq=a(4)/a(5)
        cq=a(3)/a(5)
        dq=a(2)/a(5)
        aa=(3.0d0*cq-bq**2)/3.0d0
        bb=(2.0d0*bq**3-9.0d0*bq*cq+27.0d0*dq)/27.0d0
        sq1=bb**2/4.0d0+aa**3/27.0d0
        test=bb/2.0d0-sqrt(sq1)
        if(dble(test).eq.0.0) then
        aaa=0.0d0
        else if(dble(test).gt.0.0) then
       aaa=-(test)**(1.0d0/3.0d0)
        else
       aaa=(-test)**(1.0d0/3.0d0)
        endif
        test=bb/2.0d0+sqrt(sq1)
        if(dble(test).eq.0.0) then
        bbb=0.0d0
        else if(dble(test).gt.0.0) then
        bbb=-(test)**(1.0d0/3.0d0)
        else
        bbb=(-test)**(1.0d0/3.0d0)
        endif
        roots(1)=0.0
        roots(2)=aaa+bbb-bq/3.0d0
        roots(3)=-(aaa+bbb)/2.0d0+(aaa-bbb)*sqrtm3/2.0d0-bq/3.0d0
        roots(4)=-(aaa+bbb)/2.0d0-(aaa-bbb)*sqrtm3/2.0d0-bq/3.0d0
     endif
     return
   endif
   !----
   !  solve a quartic equation
   !----
   b = a(4)/(4.0d0*a(5))
   c = a(3)/a(5)
   d = a(2)/a(5)
   e = a(1)/a(5)
   b2 = b*b

   p = 0.5d0*(c - 6.0d0*b2)
   q = d - 2.0d0*b*(c - 4.0d0*b2)
   r = b2*(c - 3.0d0*b2) - b*d + e
   !----
   !  solve the resolvent cubic equation. the cubic has at least one
   !  nonnegative real root.  if w1, w2, w3 are the roots of the cubic
   !  then the roots of the original equation are
   !     roots = -b + csqrt(w1) + csqrt(w2) + csqrt(w3)
   !  where the signs of the square roots are chosen so
   !  that csqrt(w1) * csqrt(w2) * csqrt(w3) = -q/8.
   !----
   temp(1) = -q*q/64.0d0
   temp(2) = 0.25d0*(p*p - r)
   temp(3) =  p
   temp(4) = 1.0d0
   bq=temp(3)
   cq=temp(2)
   dq=temp(1)
   aa=(3.0d0*cq-bq**2)/3.0d0
   bb=(2.0d0*bq**3-9.0d0*bq*cq+27.0d0*dq)/27.0d0
   sq1=bb**2/4.0d0+aa**3/27.0d0
   test=bb/2.0d0-sqrt(sq1)
   if(dble(test).eq.0.0) then
      aaa=0.0d0
   else if(dble(test).gt.0.0) then
      aaa=-(test)**(1.0d0/3.0d0)
   else
      aaa=(-test)**(1.0d0/3.0d0)
   endif
   test=bb/2.0d0+sqrt(sq1)
   if(dble(test).eq.0.0) then
      bbb=0.0d0
   else if(dble(test).gt.0.0) then
      bbb=-(test)**(1.0d0/3.0d0)
   else
      bbb=(-test)**(1.0d0/3.0d0)
   endif
   roots(1)=aaa+bbb-bq/3.0d0
   roots(2)=-(aaa+bbb)/2.0d0+(aaa-bbb)*sqrtm3/2.0d0-bq/3.0d0
   roots(3)=-(aaa+bbb)/2.0d0-(aaa-bbb)*sqrtm3/2.0d0-bq/3.0d0
   if(aimag(roots(2)).ne.0.0d0) go to 60
   !----
   !  the resolvent cubic has only real roots.
   !  reorder the roots in increasing order.
   !----
   xx(1) = dble(roots(1))
   xx(2) = dble(roots(2))
   xx(3) = dble(roots(3))
   do 25 j=2,3
   x=xx(j)
   do 10 i=j-1,1,-1
   if(xx(i).le.x) goto 20
   xx(i+1)=xx(i)
   10 continue
   i=0
   20 xx(i+1)=x
   25 continue

   u = 0.0d0
   if(xx(3).gt.0.0d0) u = sqrt(xx(3))
   if(xx(2).le.0.0d0) go to 41
   if(xx(1).ge.0.0d0) go to 30
   if(abs(xx(1)).gt.xx(2)) go to 40
   xx(1) = 0.0d0

   30 xx(1) = sqrt(xx(1))
   xx(2) = sqrt(xx(2))
   if(q.gt.0.0d0) xx(1) = -xx(1)
   temp(1) = (( xx(1) + xx(2)) + u) - b
   temp(2) = ((-xx(1) - xx(2)) + u) - b
   temp(3) = (( xx(1) - xx(2)) - u) - b
   temp(4) = ((-xx(1) + xx(2)) - u) - b

   do j=1,3
     k=minloc(temp(j:))
     if(j.ne.k(1)) then
       work = temp(j)
       temp(j) = temp(k(1))
       temp(k(1)) = work
     endif
   enddo

   if(abs(temp(1)).ge.0.1d0*abs(temp(4))) go to 31
   t = temp(2)*temp(3)*temp(4)
   if(t.ne.0.0d0) temp(1) = e/t
   31 roots(1) = cmplx(temp(1), 0.0d0, kind=kind(roots))
   roots(2) = cmplx(temp(2), 0.0d0, kind=kind(roots))
   roots(3) = cmplx(temp(3), 0.0d0, kind=kind(roots))
   roots(4) = cmplx(temp(4), 0.0d0, kind=kind(roots))
   return

   40 v1 = sqrt(abs(xx(1)))
   v2 = 0.0d0
   go to 50
   41 v1 = sqrt(abs(xx(1)))
   v2 = sqrt(abs(xx(2)))
   if(q < 0.0d0) u = -u

   50 x = -u - b
   y = v1 - v2
   roots(1) = cmplx(x, y, kind=kind(roots))
   roots(2) = cmplx(x,-y, kind=kind(roots))
   x =  u - b
   y = v1 + v2
   roots(3) = cmplx(x, y, kind=kind(roots))
   roots(4) = cmplx(x,-y, kind=kind(roots))
   return
   !----
   !  the resolvent cubic has complex roots.
   !----
   60 t = dble(roots(1))
   x = 0.0d0
   if(t < 0.0d0) then
     go to 61
   else if(t.eq.0.0d0) then
     go to 70
   else
     go to 62
   endif
   61 h = abs(dble(roots(2))) + abs(aimag(roots(2)))
   if(abs(t).le.h) go to 70
   go to 80
   62 x = sqrt(t)
   if(q.gt.0.0d0) x = -x

   70 w = sqrt(roots(2))
   u = 2.0d0*dble(w)
   v = 2.0d0*abs(aimag(w))
   t =  x - b
   xx(1) = t + u
   xx(2) = t - u
   if(abs(xx(1)).le.abs(xx(2))) go to 71
   t = xx(1)
   xx(1) = xx(2)
   xx(2) = t
   71 u = -x - b
   h = u*u + v*v
   if(xx(1)*xx(1) < 0.01d0*min(xx(2)*xx(2),h)) xx(1) = e/(xx(2)*h)
   roots(1) = cmplx(xx(1), 0.0d0, kind=kind(roots))
   roots(2) = cmplx(xx(2), 0.0d0, kind=kind(roots))
   roots(3) = cmplx(u, v, kind=kind(roots))
   roots(4) = cmplx(u,-v, kind=kind(roots))
   return

   80 v = sqrt(abs(t))
   roots(1) = cmplx(-b, v, kind=kind(roots))
   roots(2) = cmplx(-b,-v, kind=kind(roots))
   roots(3) = roots(1)
   roots(4) = roots(2)
   return
   end subroutine alquar
   !
   subroutine alguer(a,m,x,its,lfail)
   !
   !-----------------------------------------------------------------------
   !
   !purpose:
   ! find one root of a polynomial.
   !
   !author(s): A. Hebert
   !
   !parameters: input
   ! a    polynomial coefficients. dimension a(m+1)
   ! m    polynomial order.
   !
   !parameters: output
   ! x       complex single root.
   ! its     number of iterations.
   ! lfail   set to .true. in case of failure.
   !
   !-----------------------------------------------------------------------
   !
   implicit real(kr)(a-h,o-z)
   !----
   !  subroutine arguments
   !----
   integer m,its
   complex*16 a(m+1),x
   logical lfail
   !----
   !  local variables
   !----
   integer maxit,mr,mt
   real(kr) epss
   parameter (epss=2.d-7,mr=8,mt=10,maxit=mt*mr)
   integer iter,j
   real(kr) abx,abp,abm,err,frac(mr)
   complex*16 dx,x1,b,d,f,g,h,sq,gp,gm,g2,tmp
   save frac
   data frac /.5d0,.25d0,.75d0,.13d0,.38d0,.62d0,.88d0,1.d0/
   !
   lfail=.false.
   do 12 iter=1,maxit
     its=iter
     b=a(m+1)
     err=abs(b)
     d=cmplx(0.d0,0.d0,kind=kind(d))
     f=cmplx(0.d0,0.d0,kind=kind(f))
     abx=abs(x)
     do j=m,1,-1
       f=x*f+d
       d=x*d+b
       b=x*b+a(j)
       err=abs(b)+abx*err
     enddo
     err=epss*err
     if(abs(b).le.err) then
       return
     else
       g=d/b
       g2=g*g
       h=g2-2.d0*f/b
       sq=sqrt((m-1)*(m*h-g2))
       gp=g+sq
       gm=g-sq
       abp=abs(gp)
       abm=abs(gm)
       if(abp.lt.abm) gp=gm
       if(max(abp,abm).gt.0.d0) then
      dx=m/gp
       else
      tmp=cmplx(log(1.d0+abx),dble(iter),kind=kind(tmp))
      dx=exp(tmp)
       endif
     endif
     x1=x-dx
     if(x.eq.x1)return
     if(mod(iter,mt).ne.0) then
       x=x1
     else
       x=x-dx*frac(iter/mt)
     endif
   12 continue
   lfail=.true.
   return
   end subroutine alguer
   !
   subroutine calmpa(maxnor,nor,jini,weight,basept,demp,sp)
   !
   !-----------------------------------------------------------------------
   !
   !purpose:
   ! compute a set of partial base points preserving nor partial moments
   ! of a function using the modified Ribon approach.
   !
   !author(s): A. Hebert
   !
   !parameters: input
   ! maxnor  maximum order of a probability table.
   ! nor     number of partial moments to preserve.
   ! jini    minimum order of the partial moment we want to preserve. We must
   !         have 1-nor <= jini <= 0 (order 0 moment is always preserved).
   ! weight  weights of the probability table.
   ! basept  base points of the probability table.
   ! demp    partial moments.
   !
   !parameters: output
   ! sp      base points for the partial cross section.
   !
   !-----------------------------------------------------------------------
   !
   use util   ! provides error
   implicit real(kr)(a-h,o-z)
   !----
   ! subroutine arguments
   !----
   integer maxnor,nor,jini
   real(kr) demp(jini:nor+jini-1)
   real(kr) weight(nor),basept(nor),sp(nor)
   !----
   ! local variables
   !----
   integer i,j,k,j0
   real(kr) dda(0:maxnor),dd,dsigx
   !
   if(nor.gt.maxnor) call error('calmpa','storage overflow',' ')
   if(nor.le.0) call error('calmpa','negative or zero value of nor',' ')
   if((1-nor.gt.jini).or.(jini.gt.0)) call error('calmpa','inconsistent va' &
     & //'lue of jini',' ')
   !
   do 50 i=1,nor
   dsigx=dble(basept(i))
   dda(0)=1.0d0
   j0=0
   do 20 j=1,nor
   if(j.eq.i) go to 20
   j0=j0+1
   dda(j0)=dda(j0-1)
   do 10 k=1,j0-1
   dda(j0-k)=dda(j0-k-1)-dda(j0-k)*dble(basept(j))
   10 continue
   dda(0)=-dda(0)*dble(basept(j))
   20 continue
   dd=0.0d0
   do 30 j=0,nor-1
   dd=dd+dda(j)*demp(j+jini)
   30 continue
   do 40 j=1,nor
   if(j.ne.i) dd=dd/(dble(basept(j))-dsigx)
   40 continue
   sp(i)=real(((-1.0d0)**(nor-1))*dd*dsigx**(-jini))/weight(i)
   50 continue
   return
   end subroutine calmpa
   !
   subroutine calcat(maxnor,npar,ndil,demt,demp,ipreci,lnoraj,sigerd,sefr, &
   & nor,prosig,errbst)
   !
   !-----------------------------------------------------------------------
   !
   !purpose:
   ! Compute the weights and base points for a principal cross-section type
   ! and the partial base points for npar partial reactions.
   !
   !author(s): A. Hebert
   !
   !parameters: input
   ! maxnor  maximum order of a probability table.
   ! npar    number of partial cross sections types.
   ! ndil    number of dilutions used to test the accuracy of the table.
   ! demt    moments of the principal cross section.
   ! demp    moments for each partial cross section.
   ! ipreci  accuracy criterion for the table (=1/=2/=3).
   ! lnoraj  algorithm flag (=.true.: find an order nor.le.maxnor
   !         corresponding to accuracy ipreci; =.false.: compute the table at
   !         order nor. if this is impossible, try an order smaller than nor).
   ! sigerd  list of dilutions used to test the accuracy of the table.
   ! sefr    list of reference self-shielded cross sections corresponding
   !         to each cross-section type and each dilution.
   !
   !parameters: input/output
   ! nor     input order of the table if lnoraj=.false. and
   !         output order of the table if lnoraj=.true.
   !
   !parameters: output
   ! prosig  probability table.
   !         prosig(inor,1): weights;
   !         prosig(inor,2): base points for the principal x-s;
   !         prosig(inor,3): base points for a partial x-s;
   !         etc.
   ! errbst  probability table error.
   !
   !-----------------------------------------------------------------------
   !
   use util   ! provides error
   implicit real(kr)(a-h,o-z)
   !----
   !  subroutine arguments
   !----
   integer maxnor,npar,ndil,ipreci,nor
   real(kr) sigerd(ndil),sefr(npar+2,ndil),prosig(maxnor,2+npar),errbst
   real(kr) demt(1-maxnor:maxnor),demp(-maxnor/2:(maxnor-1)/2,npar)
   logical lnoraj
   !----
   ! local variables
   !----
   integer i,idil,inor,ior,ipar,jini,ier
   !----
   !  allocatable arrays
   !----
   real(kr), allocatable, dimension(:,:) :: seff
   real(kr), allocatable, dimension(:,:,:) :: prosic
   !
   eps=0.2**(1+ipreci)
   if(.not.lnoraj) then
     ! compute a table for an imposed value of nor.
     if(nor.gt.maxnor) call error('calcat','invalid input order',' ')
     10 call alprtb(maxnor,nor,1-nor,demt(1-nor),ier,prosig(1,1),prosig(1,2))
     if(ier.ne.0) then
       nor=nor-1
       go to 10
     endif
     jini=-nor/2
     do ipar=1,npar
       call calmpa(maxnor,nor,jini,prosig(1,1),prosig(1,2),demp(jini,ipar), &
       & prosig(1,ipar+2))
     enddo
     !----
     ! compute the self-shielded cross sections from the table.
     !----
     allocate(seff(npar+2,ndil))
     do idil=1,ndil
       seff(:npar+2,idil)=0.0
       do ior=1,nor
         astpd=sigerd(idil)*prosig(ior,1)/(prosig(ior,2)+sigerd(idil))
         seff(1,idil)=seff(1,idil)+astpd
         do ipar=2,npar+2
           seff(ipar,idil)=seff(ipar,idil)+astpd*prosig(ior,ipar)
         enddo
       enddo
       do ipar=2,npar+2
         seff(ipar,idil)=seff(ipar,idil)/seff(1,idil)
       enddo
     enddo ! idil
     !----
     ! compute the table accuracy.
     !----
     errbst=0.0
     do idil=1,ndil
       do i=2,npar+2
         errbst=max(errbst,abs(seff(i,idil)-sefr(i,idil))/abs(sefr(i,ndil)))
       enddo
     enddo
     deallocate(seff)
     return
   endif
   !----
   ! compute the table for each available order.
   !----
   allocate(prosic(maxnor,2+npar,maxnor),seff(npar+2,ndil))
   errbst=1.0e10
   do inor=1,maxnor
     call alprtb(maxnor,inor,1-inor,demt(1-inor),ier,prosic(1,1,inor), &
     & prosic(1,2,inor))
     if(ier.ne.0) cycle
     jini=-inor/2
     do ipar=1,npar
       call calmpa(maxnor,inor,jini,prosic(1,1,inor),prosic(1,2,inor), &
       & demp(jini,ipar),prosic(1,ipar+2,inor))
     enddo ! ipar
     !----
     ! compute the self-shielded cross sections from the table.
     !----
     do idil=1,ndil
       seff(:npar+2,idil)=0.0
       do ior=1,inor
         astpd=sigerd(idil)*prosic(ior,1,inor)/(prosic(ior,2,inor)+sigerd(idil))
         seff(1,idil)=seff(1,idil)+astpd
         do ipar=2,npar+2
           seff(ipar,idil)=seff(ipar,idil)+astpd*prosic(ior,ipar,inor)
         enddo
       enddo
       do ipar=2,npar+2
         seff(ipar,idil)=seff(ipar,idil)/seff(1,idil)
       enddo
     enddo ! idil
     !----
     ! compute the table accuracy.
     !----
     ermax=0.0
     do idil=1,ndil
       do i=2,npar+2
         ermax=max(ermax,abs(seff(i,idil)-sefr(i,idil))/abs(sefr(i,ndil)))
       enddo
     enddo
     if(1.2*ermax.lt.errbst) then
        nor=inor
        errbst=ermax
        if(ermax.lt.eps) go to 100
     endif
   enddo ! inor
   !----
   ! select the order nor table.
   !----
   100 do ipar=1,2+npar
     do i=1,nor
       prosig(i,ipar)=prosic(i,ipar,nor)
     enddo
   enddo
   !----
   !  scratch storage deallocation
   !----
   deallocate(seff,prosic)
   return
   end subroutine calcat
   !
   subroutine cfvect(maxgr,maxnl,maxnz,nin,matd,mt,ytemp,ng,nl,nz,rv,flux,exist)
   !-----------------------------------------------------------------
   !   utility routine for recovering a vector reaction from a gendf tape.
   !   Group 1 is the most thermal group.
   !   input parameters:
   !   maxgr   first dimension of vector rv. Maximum number of energy groups.
   !   maxnl   second dimension of vector rv. Maximum number of Legendre orders.
   !   maxnz   third dimensions of vector rv. Maximum number of dilutions.
   !   nin     unit number of the gendf tape.
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
   integer maxgr,maxnl,maxnz,nin,matd,mt,ng,nl,nz,ig,il,iz,loc,locf,nb,ng2,nw
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
   if(ng.gt.maxgr) call error('cfvect','maxgr overflow',' ')
   if(nl.gt.maxnl) call error('cfvect','maxnl overflow',' ')
   if(nz.gt.maxnz) call error('cfvect','maxnz overflow',' ')
   lfind=.true.
   do while (mth.ne.0) ! loop over energy groups
     if(.not.lfind) call listio(nin,0,0,scr(1),nb,nw)
     ng2=l1h
     ig=n2h
     if(ig.eq.0) exit
     lfind=.false.
     loc=1+nw
     do while (nb.ne.0)
       if(loc+302.gt.maxa) call error('cfvect','endf input size exceeded', &
       & ' ')
       call moreio(nin,0,0,scr(loc),nb,nw)
       loc=loc+nw
     enddo
     do il=1,nl
       do iz=1,nz
         locf=1+lz+nl*(iz-1)+(il-1)
         flux(ig,il,iz)=scr(locf)
         rv(ig,il,iz)=scr(locf+nl*nz)
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
   end subroutine cfvect
end module calenm
