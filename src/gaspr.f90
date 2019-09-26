module gaspm
   ! provides subroutine gaspr for NJOY2016
   use locale
   implicit none
   private
   public gaspr
contains

   subroutine gaspr
   !-------------------------------------------------------------------
   !
   ! Add gas production reactions (mt203-207) to the pendf tape.
   ! Any old gas sections on the input pendf tape are deleted.
   ! The directory is updated to show the new reactions.
   ! This module can be run anywhere in the pendf preparation
   ! sequence as long as it is somewhere after broadr.
   !
   ! If the input pendf tape omits mt103 to mt107, but it does
   ! have the corresponding charged-particle cross sections, they
   ! are processed and will appear in the appropriate mt20x section.
   !
   !---input specifications (free format)--------------------------
   !
   ! card 1
   !    nendf    unit for endf tape
   !    nin      unit for input pendf tape
   !    nout     unit for output pendf tape
   !
   !-------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso,nsyse
   use endf   ! provides endf routines and variables
   use util   ! provides timer,error,openz,sigfig
   ! internals
   integer::nscr1,iprint,nendf,npend,noutp,nb,nw
   integer::n6,nx,nsub,matd,mf6mt5,i,idone,mfi,mti,mfd,mtd
   integer::nk,lsix,l203,l204,l205,l206,l207
   integer::ik,izap,law,ll,ne,ie,nmu,imu,itemp,mfb
   integer::izg,izr,idis,lr,ngas,j,ip,ir,nsec
   integer::n203,n204,n205,n206,n207,jg,ii,k,nold
   integer::ni,nm,mtb,ifini,np,istart,iend,ib
   integer::mpmin,mpmax,mdmin,mdmax,mtmin,mtmax,m3min,m3max,m4min,m4max
   integer::maxg
   real(kr)::time,za,awr,zain,thrg,en,enext,y
   real(kr)::y203,y204,y205,y206,y207,yyy,xnext
   real(kr)::a(10000)
   real(kr)::b(6)
   real(kr),dimension(:),allocatable::egas
   real(kr),dimension(:,:),allocatable::sgas
   real(kr)::six(5000)
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::zero=0
   real(kr),parameter::hund=100

   nscr1=10
   ! set iprint nonzero for gas production printout on listing
   iprint=0

   !--read user input
   call timer(time)
   write(nsyso,'(/'' gaspr...'',&
     &''add gas production cross sections'',&
     &27x,f8.1,''s'')') time
   write(nsyse,'(/'' gaspr...'',60x,f8.1,''s'')') time
   read(nsysi,*) nendf,npend,noutp
   if (npend.lt.0) nscr1=-nscr1
   if (npend*noutp.lt.0) call error('gaspr',&
     'npend and noutp must both be binary',&
     'or both be coded')
   write(nsyso,'(/'' units:'',3i6)') nendf,npend,noutp
   call openz(nendf,0)
   call openz(npend,0)
   call openz(noutp,1)
   call openz(nscr1,1)

   !--check endf tape for mf6,mt5.
   !--also set flags for the absence or presence of charged-particle
   !--reactions and the endf version dependent partial cross section
   !--mt ranges.
   call repoz(nendf)
   call tpidio(nendf,0,0,a(1),nb,nw)
   call contio(nendf,0,0,a(1),nb,nw)
   za=c1h
   awr=c2h
   n6=n2h
   call contio(nendf,0,0,a(1),nb,nw)
   if (n1h.ne.0) then
      iverf=4
      nx=n6
   else if (n2h.eq.0) then
      iverf=5
   else
      iverf=6
   endif
   call skiprz(nendf,-1)
   if (iverf.ge.5) call contio(nendf,0,0,a(1),nb,nw)
   nsub=10
   zain=1
   if (iverf.ge.6) then
      call contio(nendf,0,0,a(1),nb,nw)
      nsub=n1h
      zain=int(nsub/10)
   endif
   call hdatio(nendf,0,0,a(1),nb,nw)
   matd=math
   do while (nb.ne.0)
      call moreio(nendf,0,0,a(1),nb,nw)
   enddo
   if (iverf.ne.4) nx=n2h
   nw=nx
   call dictio(nendf,0,0,a(1),nb,nw)
   mf6mt5=0
   i=1
   idone=0
   if (iverf .ge. 6) then
       mpmin=600
       mpmax=649
       mdmin=650
       mdmax=699
       mtmin=700
       mtmax=749
       m3min=750
       m3max=799
       m4min=800
       m4max=849
   else
       mpmin=700
       mpmax=718
       mdmin=720
       mdmax=738
       mtmin=740
       mtmax=758
       m3min=760
       m3max=768
       m4min=780
       m4max=798
   endif
   do while (i.le.nw.and.idone.eq.0)
      mfi=nint(a(1+i+1))
      mti=nint(a(1+i+2))
      if (mfi.gt.6) then
         idone=1
      else
         if (mfi.eq.6.and.mti.eq.5) mf6mt5=1
         i=i+6
      endif
   enddo
   call repoz(nendf)
   if (mf6mt5.ne.0) then
      write(nsyso,'(/'' mf6,mt5 found'')')
      mfd=6
      mtd=5
      call findf(matd,mfd,mtd,nendf)
      call contio(nendf,0,0,a(1),nb,nw)
      nk=n1h
      lsix=1
      l203=0
      l204=0
      l205=0
      l206=0
      l207=0
      do ik=1,nk
         call tab1io(nendf,0,0,a(1),nb,nw)
         izap=nint(c1h)
         law=l2h
         ll=1+nw
         do while (nb.ne.0)
            call moreio(nendf,0,0,a(ll),nb,nw)
            ll=ll+nw
         enddo
         if (izap.eq.1001) then
            l203=lsix
            do i=1,nw
               six(lsix+i-1)=a(i)
            enddo
            lsix=lsix+nw
         else if (izap.eq.1002) then
            l204=lsix
            do i=1,nw
               six(lsix+i-1)=a(i)
            enddo
            lsix=lsix+nw
         else if (izap.eq.1003) then
            l205=lsix
            do i=1,nw
               six(lsix+i-1)=a(i)
            enddo
            lsix=lsix+nw
         else if (izap.eq.2003) then
            l206=lsix
            do i=1,nw
               six(lsix+i-1)=a(i)
            enddo
            lsix=lsix+nw
         else if (izap.eq.2004) then
            l207=lsix
            do i=1,nw
               six(lsix+i-1)=a(i)
            enddo
            lsix=lsix+nw
         endif
         if (law.eq.6) then
            call contio(nendf,0,0,a,nb,nw)
         else if (law.eq.1.or.law.eq.2.or.law.eq.5) then
            call tab2io(nendf,0,0,a,nb,nw)
            ne=n2h
            do ie=1,ne
               call listio(nendf,0,0,a,nb,nw)
               do while (nb.ne.0)
                  call moreio(nendf,0,0,a,nb,nw)
               enddo
            enddo
         else if (law.eq.7) then
            call tab2io(nendf,0,0,a,nb,nw)
            ne=n2h
            do ie=1,ne
               call tab2io(nendf,0,0,a,nb,nw)
               nmu=n2h
               do imu=1,nmu
                  call tab1io(nendf,0,0,a,nb,nw)
                  do while (nb.ne.0)
                     call moreio(nendf,0,0,a,nb,nw)
                  enddo
               enddo
            enddo
         endif
      enddo
   endif

   !--loop over all temperatures on the pendf tape
   call repoz(npend)
   call repoz(noutp)
   call tpidio(npend,noutp,0,a(1),nb,nw)
   itemp=0
  225 continue
   itemp=itemp+1
   call contio(npend,0,0,b(1),nb,nw)
   if (mfh.eq.0) go to 550

   !--copy data through file3, include any partial
   !-charged-particle cross sections.
   call repoz(nscr1)
   call contio(0,nscr1,0,b(1),nb,nw)
   call tofend(npend,nscr1,0,a(1))
   call tofend(npend,nscr1,0,a(1))
   idone=0
   do while (idone.eq.0)
      call contio(npend,0,0,b(1),nb,nw)
      mfb=mfh
      mtb=mth
      if (mth.gt.m4max.or.mth.eq.0) then
         idone=1
      else
         call contio(0,nscr1,0,b(1),nb,nw)
         call tosend(npend,nscr1,0,a(1))
      endif
   enddo
   call contio(0,nscr1,0,b(1),nb,nw)

   !--reposition npend to the location where gas production
   !--cross sections would go.
   call repoz(npend)
   call tpidio(npend,0,0,a(1),nb,nw)
   if (itemp.gt.1) then
      do i=2,itemp
         call tomend(npend,0,0,a(1))
      enddo
   endif
   call tofend(npend,0,0,a(1))
   call tofend(npend,0,0,a(1))
   idone=0
   do while (idone.eq.0)
      call contio(npend,0,0,b(1),nb,nw)
      if (mth.gt.200.or.mth.eq.0) then
         idone=1
      else
         call tosend(npend,0,0,a(1))
      endif
   enddo

   !--check the reactions in file 3 for gas threshold
   call repoz(nscr1)
   call tofend(nscr1,0,0,a(1))
   call tofend(nscr1,0,0,a(1))
   thrg=emax
  240 continue
   call contio(nscr1,0,0,a(1),nb,nw)
   if (mth.gt.200.and.mth.lt.mpmin) go to 245
   if (mth.gt.m4max.or.mth.eq.0) go to 250
   izg=0
   izr=nint(za+zain)
   en=0
   call gety1(en,enext,idis,y,nscr1,a(1))
   lr=l2h
   if (mth.ge.1.and.mth.le.4) go to 245
   if (mth.eq.5.and.mf6mt5.eq.1) izg=1
   if (mth.ge.6.and.mth.le.10) go to 245
   if (mth.eq.11) izg=1
   if (mth.eq.11) izr=izr-1004
   if (mth.ge.12.and.mth.le.15) go to 245
   if (mth.eq.16) izr=izr-2
   if (mth.eq.17) izr=izr-3
   if (mth.ge.18.and.mth.le.21) go to 245
   if (mth.ge.22.and.mth.le.25) izg=1
   if (mth.eq.22) izr=izr-2005
   if (mth.eq.23) izr=izr-6013
   if (mth.eq.24) izr=izr-2006
   if (mth.eq.25) izr=izr-2007
   if (mth.ge.26.and.mth.le.27) go to 245
   if (mth.ge.28.and.mth.le.32) izg=1
   if (mth.eq.28) izr=izr-1002
   if (mth.eq.29) izr=izr-4009
   if (mth.eq.30) izr=izr-4010
   if (mth.eq.31) go to 245
   if (mth.ge.32.and.mth.le.36) izg=1
   if (mth.eq.32) izr=izr-1003
   if (mth.eq.33) izr=izr-1004
   if (mth.eq.34) izr=izr-2004
   if (mth.eq.35) izr=izr-5011
   if (mth.eq.36) izr=izr-5012
   if (mth.eq.37) izr=izr-4
   if (mth.ge.38.and.mth.le.40) go to 245
   if (mth.ge.41.and.mth.le.42) izg=1
   if (mth.eq.41) izr=izr-1003
   if (mth.eq.42) izr=izr-1004
   if (mth.eq.43) go to 245
   if (mth.ge.44.and.mth.le.45) izg=1
   if (mth.eq.44) izr=izr-2003
   if (mth.eq.45) izr=izr-3006
   if (mth.ge.46.and.mth.le.50) go to 245
   if (mth.ge.51.and.mth.le.91) then
      izr=izr-1
      if (lr.ge.22.and.lr.le.25) izg=1
      if (lr.eq.22) izr=izr-2004
      if (lr.eq.23) izr=izr-6012
      if (lr.eq.24) izr=izr-2005
      if (lr.eq.25) izr=izr-2006
      if (lr.ge.28.and.lr.le.30) izg=1
      if (lr.eq.28) izr=izr-1001
      if (lr.eq.29) izr=izr-4008
      if (lr.eq.30) izr=izr-4009
      if (lr.ge.32.and.lr.le.36) izg=1
      if (lr.eq.32) izr=izr-1002
      if (lr.eq.33) izr=izr-1003
      if (lr.eq.34) izr=izr-2003
      if (lr.eq.35) izr=izr-5010
      if (lr.eq.36) izr=izr-5011
   endif
   if (mth.ge.92.and.mth.le.101) go to 245
   if (mth.ge.103.and.mth.le.117) izg=1
   if (mth.eq.103) izr=izr-1001
   if (mth.eq.104) izr=izr-1002
   if (mth.eq.105) izr=izr-1003
   if (mth.eq.106) izr=izr-2003
   if (mth.eq.107) izr=izr-2004
   !--always skip levels for mt103-107
   if (mth.ge.mpmin.and.mth.le.mpmax) go to 245
   if (mth.ge.mdmin.and.mth.le.mdmax) go to 245
   if (mth.ge.mtmin.and.mth.le.mtmax) go to 245
   if (mth.ge.m3min.and.mth.le.m3max) go to 245
   if (mth.ge.m4min.and.mth.le.m4max) go to 245
   if (mth.eq.108) izr=izr-4008
   if (mth.eq.109) izr=izr-6012
   if (mth.eq.111) izr=izr-2002
   if (mth.eq.112) izr=izr-3005
   if (mth.eq.113) izr=izr-5011
   if (mth.eq.114) izr=izr-5010
   if (mth.eq.115) izr=izr-2003
   if (mth.eq.116) izr=izr-2004
   if (mth.eq.117) izr=izr-3006
   if (mth.ge.154.and.mth.le.159) izg=1
   if (mth.ge.162.and.mth.le.200) izg=1
   if (mth.eq.152) izr=izr-5
   if (mth.eq.153) izr=izr-6
   if (mth.eq.154) izr=izr-1005
   if (mth.eq.155) izr=izr-3007
   if (mth.eq.156) izr=izr-1005
   if (mth.eq.157) izr=izr-1005
   if (mth.eq.158) izr=izr-3007
   if (mth.eq.159) izr=izr-3007
   if (mth.eq.160) izr=izr-7
   if (mth.eq.161) izr=izr-8
   if (mth.eq.162) izr=izr-1006
   if (mth.eq.163) izr=izr-1007
   if (mth.eq.164) izr=izr-1008
   if (mth.eq.165) izr=izr-2008
   if (mth.eq.166) izr=izr-2009
   if (mth.eq.167) izr=izr-2010
   if (mth.eq.168) izr=izr-2011
   if (mth.eq.169) izr=izr-1006
   if (mth.eq.170) izr=izr-1007
   if (mth.eq.171) izr=izr-1008
   if (mth.eq.172) izr=izr-1006
   if (mth.eq.173) izr=izr-1007
   if (mth.eq.174) izr=izr-1008
   if (mth.eq.175) izr=izr-1009
   if (mth.eq.176) izr=izr-2005
   if (mth.eq.177) izr=izr-2006
   if (mth.eq.178) izr=izr-2007
   if (mth.eq.179) izr=izr-2005
   if (mth.eq.180) izr=izr-4011
   if (mth.eq.181) izr=izr-3008
   if (mth.eq.182) izr=izr-2005
   if (mth.eq.183) izr=izr-2004
   if (mth.eq.184) izr=izr-2005
   if (mth.eq.185) izr=izr-2006
   if (mth.eq.186) izr=izr-3005
   if (mth.eq.187) izr=izr-3006
   if (mth.eq.188) izr=izr-3007
   if (mth.eq.189) izr=izr-3008
   if (mth.eq.190) izr=izr-2004
   if (mth.eq.191) izr=izr-3004
   if (mth.eq.192) izr=izr-3005
   if (mth.eq.193) izr=izr-4007
   if (mth.eq.194) izr=izr-2006
   if (mth.eq.195) izr=izr-4012
   if (mth.eq.196) izr=izr-3009
   if (mth.eq.197) izr=izr-3003
   if (mth.eq.198) izr=izr-3004
   if (mth.eq.199) izr=izr-4009
   if (mth.eq.200) izr=izr-2007
   if (izr.eq.4008) izg=1
   if (izg.eq.0.and.&
     (izr.gt.2004.or.izr.le.0)) go to 245
   if (enext.lt.thrg) thrg=enext
  245 continue
   call tosend(nscr1,0,0,a(1))
   go to 240
  250 continue
   if (itemp.eq.1)&
      write(nsyso,'(/'' the gas production threshold is'',&
      &1pe12.4,'' ev'')') thrg

   !--read through the scratch tape
   !--collecting gas production in the process
   !--use the energy grid of mt1 starting at thrg
   call repoz(nscr1)
   call tofend(nscr1,0,0,a(1))
   call tofend(nscr1,0,0,a(1))
   call contio(nscr1,0,0,a(1),nb,nw)
   en=0
   call gety1(en,enext,idis,y,nscr1,a(1))
   maxg=nint(a(6))
    if (allocated(egas)) then
       deallocate(egas)
       deallocate(sgas)
    endif
   allocate (egas(maxg))
   allocate (sgas(5,maxg))
   enext=thrg
   i=0
   do while (enext.lt.emax)
      i=i+1
      if (i.gt.maxg) call error('gaspr',&
        'too many gas production energy points',' ')
      en=enext
      call gety1(en,enext,idis,y,nscr1,a(1))
      egas(i)=en
   enddo
   ngas=i
   if (itemp.eq.1) write(nsyso,'(/'' found'',i6,'' points'')') ngas
   call tosend(nscr1,0,0,a(1))
   if (itemp.eq.1) write(nsyso,&
     '(/'' pendf mt  mt203  mt204  mt205  mt206  mt207''/&
     &'' ________  _____  _____  _____  _____  _____'')')
   do i=1,ngas
      do j=1,5
         sgas(j,i)=0
      enddo
   enddo

   !--loop over other reactions and
   !--sum up gas production values
  270 continue
   call contio(nscr1,0,0,a(1),nb,nw)
   if (mth.gt.200.and.mth.lt.mpmin) go to 310
   if (mth.gt.m4max.or.mth.eq.0) go to 330
   if (mth.le.4) go to 310
   if (mth.ge.6.and.mth.le.10) go to 310
   if (mth.ge.12.and.mth.le.15) go to 310
   if (mth.ge.18.and.mth.le.21) go to 310
   if (mth.ge.38.and.mth.le.40) go to 310
   if (mth.eq.43) go to 310
   if (mth.ge.46.and.mth.le.50) go to 310
   if (mth.ge.92.and.mth.le.101) go to 310
   !--always skip levels for mt103-107
   if (mth.ge.mpmin.and.mth.le.mpmax) go to 310
   if (mth.ge.mdmin.and.mth.le.mdmax) go to 310
   if (mth.ge.mtmin.and.mth.le.mtmax) go to 310
   if (mth.ge.m3min.and.mth.le.m3max) go to 310
   if (mth.ge.m4min.and.mth.le.m4max) go to 310
   if (mth.eq.152.or.mth.eq.153.or.mth.eq.160.or.mth.eq.161)go to 310
   en=0
   call gety1(en,enext,idis,y,nscr1,a(1))
   lr=l2h
   izr=nint(za+zain)
   y203=0
   y204=0
   y205=0
   y206=0
   y207=0
   if (mth.eq.5.and.mf6mt5.eq.1) then
      if (l203.gt.0) y203=111
      if (l204.gt.0) y204=111
      if (l205.gt.0) y205=111
      if (l206.gt.0) y206=111
      if (l207.gt.0) y207=111
   else if (mth.eq.11) then
      izr=izr-1004
      y204=1
   else if (mth.eq.16) then
      izr=izr-2
   else if (mth.eq.17) then
      izr=izr-3
   else if (mth.eq.22) then
      izr=izr-2005
      y207=1
   else if (mth.eq.23) then
      izr=izr-6013
      y207=3
   else if (mth.eq.24) then
      izr=izr-2006
      y207=1
   else if (mth.eq.25) then
      izr=izr-2007
      y207=1
   else if (mth.eq.28) then
      izr=izr-1002
      y203=1
   else if (mth.eq.29) then
      izr=izr-4009
      y207=2
   else if (mth.eq.30) then
      izr=izr-4010
      y207=2
   else if (mth.eq.32) then
      izr=izr-1003
      y204=1
   else if (mth.eq.33) then
      izr=izr-1004
      y205=1
   else if (mth.eq.34) then
      izr=izr-2004
      y206=1
   else if (mth.eq.35) then
      izr=izr-5011
      y204=1
      y207=2
   else if (mth.eq.36) then
      izr=izr-5012
      y205=1
      y207=2
   else if (mth.eq.37) then
      izr=izr-4
   else if (mth.eq.41) then
      izr=izr-1003
      y203=1
   else if (mth.eq.42) then
      izr=izr-1004
      y203=1
   else if (mth.eq.44) then
      izr=izr-2003
      y203=2
   else if (mth.eq.45) then
      izr=izr-3006
      y203=1
      y207=1
   else if (mth.ge.51.and.mth.le.91) then
      izr=izr-1
      if (lr.eq.22) then
         izr=izr-2004
         y207=1
      else if (lr.eq.23) then
         izr=izr-6012
         y207=3
      else if (lr.eq.24) then
         izr=izr-2005
         y207=1
      else if (lr.eq.25) then
         izr=izr-2006
         y207=1
      else if (lr.eq.28) then
         izr=izr-1001
         y203=1
      else if (lr.eq.29) then
         izr=izr-4008
         y207=2
      else if (lr.eq.30) then
         izr=izr-4009
         y207=2
      else if (lr.eq.32) then
         izr=izr-1002
         y204=1
      else if (lr.eq.33) then
         izr=izr-1003
         y205=1
      else if (lr.eq.34) then
         izr=izr-2003
         y206=1
      else if (lr.eq.35) then
         izr=izr-5010
         y204=1
         y207=2
      else if (lr.eq.36) then
         izr=izr-5011
         y205=1
         y207=2
      else if (lr.eq.39) then
         izr=izr
      else if (lr.eq.40) then
         izr=izr
      endif
   else if (mth.eq.103) then
      izr=izr-1001
      y203=1
   else if (mth.eq.104) then
      izr=izr-1002
      y204=1
   else if (mth.eq.105) then
      izr=izr-1003
      y205=1
   else if (mth.eq.106) then
      izr=izr-2003
      y206=1
   else if (mth.eq.107) then
      izr=izr-2004
      y207=1
   else if (mth.eq.108) then
      izr=izr-4008
      y207=2
   else if (mth.eq.109) then
      izr=izr-6012
      y207=3
   else if (mth.eq.111) then
      izr=izr-2002
      y203=2
   else if (mth.eq.112) then
      izr=izr-3005
      y203=1
      y207=1
   else if (mth.eq.113) then
      izr=izr-5011
      y205=1
      y207=2
   else if (mth.eq.114) then
      izr=izr-5010
      y204=1
      y207=2
   else if (mth.eq.115) then
      izr=izr-2003
      y203=1
      y204=1
   else if (mth.eq.116) then
      izr=izr-2004
      y203=1
      y205=1
   else if (mth.eq.117) then
      izr=izr-3006
      y204=1
      y207=1
   else if (mth.eq.152) then
      izr=izr-5
   else if (mth.eq.153) then
      izr=izr-6
   else if (mth.eq.154) then
      izr=izr-1005
      y205=1
   else if (mth.eq.155) then
      izr=izr-3007
      y205=1
      y207=1
   else if (mth.eq.156) then
      izr=izr-1005
      y203=1
   else if (mth.eq.157) then
      izr=izr-1005
      y204=1
   else if (mth.eq.158) then
      izr=izr-3007
      y204=1
      y207=1
   else if (mth.eq.159) then
      izr=izr-3007
      y203=1
      y207=1
   else if (mth.eq.160) then
      izr=izr-7
   else if (mth.eq.161) then
      izr=izr-8
   else if (mth.eq.162) then
      izr=izr-1006
      y203=1
   else if (mth.eq.163) then
      izr=izr-1007
      y203=1
   else if (mth.eq.164) then
      izr=izr-1008
      y203=1
   else if (mth.eq.165) then
      izr=izr-2008
      y207=1
   else if (mth.eq.166) then
      izr=izr-2009
      y207=1
   else if (mth.eq.167) then
      izr=izr-2010
      y207=1
   else if (mth.eq.168) then
      izr=izr-2011
      y207=1
   else if (mth.eq.169) then
      izr=izr-1006
      y204=1
   else if (mth.eq.170) then
      izr=izr-1007
      y204=1
   else if (mth.eq.171) then
      izr=izr-1008
      y204=1
   else if (mth.eq.172) then
      izr=izr-1006
      y205=1
   else if (mth.eq.173) then
      izr=izr-1007
      y205=1
   else if (mth.eq.174) then
      izr=izr-1008
      y205=1
   else if (mth.eq.175) then
      izr=izr-1009
      y205=1
   else if (mth.eq.176) then
      izr=izr-2005
      y206=1
   else if (mth.eq.177) then
      izr=izr-2006
      y206=1
   else if (mth.eq.178) then
      izr=izr-2007
      y206=1
   else if (mth.eq.179) then
      izr=izr-2005
      y203=2
   else if (mth.eq.180) then
      izr=izr-4011
      y207=2
   else if (mth.eq.181) then
      izr=izr-3008
      y203=1
      y207=1
   else if (mth.eq.182) then
      izr=izr-2005
      y204=1
      y205=1
   else if (mth.eq.183) then
      izr=izr-2004
      y203=1
      y204=1
   else if (mth.eq.184) then
      izr=izr-2005
      y203=1
      y205=1
   else if (mth.eq.185) then
      izr=izr-2006
      y204=1
      y205=1
   else if (mth.eq.186) then
      izr=izr-3005
      y203=1
      y206=1
   else if (mth.eq.187) then
      izr=izr-3006
      y204=1
      y206=1
   else if (mth.eq.188) then
      izr=izr-3007
      y205=1
      y206=1
   else if (mth.eq.189) then
      izr=izr-3008
      y205=1
      y207=1
   else if (mth.eq.190) then
      izr=izr-2004
      y203=2
   else if (mth.eq.191) then
      izr=izr-3004
      y203=1
      y206=1
   else if (mth.eq.192) then
      izr=izr-3005
      y204=1
      y206=1
   else if (mth.eq.193) then
      izr=izr-4007
      y206=1
      y207=1
   else if (mth.eq.194) then
      izr=izr-2006
      y203=2
   else if (mth.eq.195) then
      izr=izr-4012
      y207=2
   else if (mth.eq.196) then
      izr=izr-3009
      y203=1
      y207=1
   else if (mth.eq.197) then
      izr=izr-3003
      y203=3
   else if (mth.eq.198) then
      izr=izr-3004
      y203=3
   else if (mth.eq.199) then
      izr=izr-4009
      y203=2
      y207=1
   else if (mth.eq.200) then
      izr=izr-2007
      y203=2
   endif
   if (izr.eq.1001) y203=y203+1
   if (izr.eq.1002) y204=y204+1
   if (izr.eq.1003) y205=y205+1
   if (izr.eq.2003) y206=y206+1
   if (izr.eq.2004) y207=y207+1
   if (izr.eq.4008) y207=y207+2
   if (y203.eq.zero.and.y204.eq.zero.and.y205.eq.zero&
     .and.y206.eq.zero.and.y207.eq.zero) go to 310
   enext=thrg
   if (itemp.eq.1)&
     write(nsyso,'(i8,5(4x,f3.1))') mth,y203,y204,y205,y206,y207
   i=0
   do while (i.lt.ngas)
      i=i+1
      en=egas(i)
      call gety1(en,enext,idis,y,nscr1,a(1))
      ip=2
      ir=1
      if (y203.gt.hund) then
         call terpa(yyy,egas(i),xnext,idis,six(l203),ip,ir)
         sgas(1,i)=sgas(1,i)+yyy*y
      else
         sgas(1,i)=sgas(1,i)+y203*y
      endif
      if (y204.gt.hund) then
         call terpa(yyy,egas(i),xnext,idis,six(l204),ip,ir)
         sgas(2,i)=sgas(2,i)+yyy*y
      else
         sgas(2,i)=sgas(2,i)+y204*y
      endif
      if (y205.gt.hund) then
         call terpa(yyy,egas(i),xnext,idis,six(l205),ip,ir)
         sgas(3,i)=sgas(3,i)+yyy*y
      else
         sgas(3,i)=sgas(3,i)+y205*y
      endif
      if (y206.gt.hund) then
         call terpa(yyy,egas(i),xnext,idis,six(l206),ip,ir)
         sgas(4,i)=sgas(4,i)+yyy*y
      else
         sgas(4,i)=sgas(4,i)+y206*y
      endif
      if (y207.gt.hund) then
         call terpa(yyy,egas(i),xnext,idis,six(l207),ip,ir)
         sgas(5,i)=sgas(5,i)+yyy*y
      else
         sgas(5,i)=sgas(5,i)+y207*y
      endif
   enddo
  310 continue
   call tosend(nscr1,0,0,a(1))
   go to 270
  330 continue
   if (itemp.eq.1)&
     write(nsyso,'(/'' *** means that the yield is '',&
     &''energy dependent'')')

   !--print out the gas-production results
   if (iprint.ne.0) then
      if (itemp.eq.1) then
         write(nsyso,'(/'' gas production versus energy'')')
         do i=1,ngas
           write(nsyso,'(1p,6e12.4)') egas(i),sgas(1,i),sgas(2,i),&
             sgas(3,i),sgas(4,i),sgas(5,i)
         enddo
      endif
   endif

   !--update the directory in file 1
   call repoz(nscr1)
   nsec=0
   n203=0
   n204=0
   n205=0
   n206=0
   n207=0
   do jg=1,5
      idone=0
      i=0
      do while (i.lt.ngas.and.idone.eq.0)
         i=i+1
         ii=i
         if (sgas(jg,i).ne.zero) idone=1
      enddo
      if (idone.eq.1) then
         ii=ii-1
         nsec=nsec+1
         if (jg.eq.1) n203=ngas-ii+1
         if (jg.eq.2) n204=ngas-ii+1
         if (jg.eq.3) n205=ngas-ii+1
         if (jg.eq.4) n206=ngas-ii+1
         if (jg.eq.5) n207=ngas-ii+1
      endif
   enddo
   call contio(nscr1,0,0,a(1),nb,nw)
   k=1+nw
   if (iverf.eq.4) nx=n2h
   if (iverf.ge.5) then
      call contio(nscr1,0,0,a(k),nb,nw)
      k=k+nw
   endif
   if (iverf.ge.6) then
      call contio(nscr1,0,0,a(k),nb,nw)
      k=k+nw
   endif
   call hdatio(nscr1,0,0,a(k),nb,nw)
   k=k+nw
   if (iverf.ne.4) nx=n2h
   do while (nb.ne.0)
      call moreio(nscr1,0,0,a(k),nb,nw)
      k=k+nw
   enddo
   nw=nx
   call dictio(nscr1,0,0,a(k),nb,nw)
   nold=0
   do i=1,nx
      j=k-1+6*(i-1)
      mfi=nint(a(j+3))
      mti=nint(a(j+4))
      if (mfi.eq.3.and.mti.ge.203.and.mti.le.207) nold=nold+1
   enddo
   if (nold.gt.0)&
     write(nsyso,'(/'' gas data on input pendf tape deleted'')')
   i=0
   idone=0
   do while (i.lt.nx.and.idone.eq.0)
      i=i+1
      j=k-1+6*(i-1)
      if (nint(a(j+3)).gt.3) idone=1
      if (nint(a(j+3)).eq.3.and.nint(a(j+4)).gt.200) idone=1
   enddo
   if (idone.eq.0) then
      j=k-1+6*nx
   else
      ni=6*(nsec-nold)
      nm=6*nx-j+k+5
      do i=1,nm
         a(6*nx+ni+k-i)=a(6*nx+k-i)
      enddo
   endif
   if (n203.ne.0) then
      a(j+1)=0
      a(j+2)=0
      a(j+3)=3
      a(j+4)=203
      a(j+5)=int((n203+2)/3)
      a(j+6)=1
      j=j+6
   endif
   if (n204.ne.0) then
      a(j+1)=0
      a(j+2)=0
      a(j+3)=3
      a(j+4)=204
      a(j+5)=int((n204+2)/3)
      a(j+6)=1
      j=j+6
   endif
   if (n205.ne.0) then
      a(j+1)=0
      a(j+2)=0
      a(j+3)=3
      a(j+4)=205
      a(j+5)=int((n205+2)/3)
      a(j+6)=1
      j=j+6
   endif
   if (n206.ne.0) then
      a(j+1)=0
      a(j+2)=0
      a(j+3)=3
      a(j+4)=206
      a(j+5)=int((n206+2)/3)
      a(j+6)=1
      j=j+6
   endif
   if (n207.ne.0) then
      a(j+1)=0
      a(j+2)=0
      a(j+3)=3
      a(j+4)=207
      a(j+5)=int((n207+2)/3)
      a(j+6)=1
      j=j+6
   endif
   nw=6
   k=1
   if (iverf.eq.4) a(6)=nx+nsec-nold
   call contio(0,noutp,0,a(k),nb,nw)
   k=k+nw
   if (iverf.ge.5) then
      call contio(0,noutp,0,a(k),nb,nw)
      k=k+nw
   endif
   if (iverf.ge.6) then
      call contio(0,noutp,0,a(k),nb,nw)
      k=k+nw
   endif
   if (iverf.ne.4) a(k+5)=nx+nsec-nold
   call hdatio(0,noutp,0,a(k),nb,nw)
   k=k+nw
   do while (nb.ne.0)
      call moreio(0,noutp,0,a(k),nb,nw)
      k=k+nw
   enddo
   nw=nx+nsec-nold
   nw=nx+nsec-nold
   call dictio(0,noutp,0,a(k),nb,nw)
   call tofend(nscr1,noutp,0,a(1))
   call tofend(nscr1,noutp,0,a(1))

   !--copy file 3 down to the gas production area
   idone=0
   do while (idone.eq.0)
      call contio(nscr1,0,0,b(1),nb,nw)
      mfb=mfh
      mtb=mth
      if (mth.gt.200.or.mth.eq.0) then
         idone=1
      else
         call contio(0,noutp,0,b(1),nb,nw)
         call tosend(nscr1,noutp,0,a(1))
      endif
   enddo

   !--write the gas production sections in file 3
   do jg=1,5
      i=0
      ifini=0
      do while (i.lt.ngas.and.ifini.eq.0)
         i=i+1
         ii=i
         if (sgas(jg,i).ne.zero) ifini=1
      enddo
      if (ifini.eq.1) then
         i=ii-1
         if (i.eq.0) i=1
         if (jg.eq.1) mth=203
         if (jg.eq.2) mth=204
         if (jg.eq.3) mth=205
         if (jg.eq.4) mth=206
         if (jg.eq.5) mth=207
         math=matd
         mfh=3
         a(1)=za
         a(2)=awr
         a(3)=0
         a(4)=0
         a(5)=0
         a(6)=0
         nw=6
         call contio(0,noutp,0,a(1),nb,nw)
         np=ngas-i+1
         a(1)=0
         a(2)=0
         a(3)=0
         a(4)=0
         a(5)=1
         a(6)=np
         a(7)=np
         a(8)=2
         k=8
         istart=i
         idone=0
         do while (idone.eq.0)
            iend=ngas
            if ((iend-istart).ge.npage/2) iend=istart+npage/2-1
            j=k-1
            ib=istart-1
            do while (ib.lt.iend)
               j=j+2
               ib=ib+1
               a(j)=egas(ib)
               a(j+1)=sigfig(sgas(jg,ib),7,0)
            enddo
            nw=j+1
            if (k.ne.0) then
               k=0
               call tab1io(0,noutp,0,a(1),nb,nw)
            else
               call moreio(0,noutp,0,a(1),nb,nw)
            endif
            if (nb.eq.0) then
               idone=1
            else
               istart=iend+1
            endif
         enddo
         call asend(noutp,0)
      endif
   enddo

   !--copy rest of this temperature to output file
   !--delete any gas sections on the input pend file
   mfh=mfb
   mth=mtb
   nw=6
  520 continue
   if (mfh.eq.0) go to 530
   if (mth.ge.203.and.mth.le.207) go to 525
   call contio(0,noutp,0,b(1),nb,nw)
   call tosend(npend,noutp,0,a(1))
   call contio(npend,0,0,b(1),nb,nw)
   go to 520
  525 continue
   call tosend(npend,0,0,a(1))
   call contio(npend,0,0,b(1),nb,nw)
   go to 520
  530 continue
   call contio(0,noutp,0,b(1),nb,nw)
   call tomend(npend,noutp,0,a(1))
   go to 225
  550 continue
   call contio(0,noutp,0,b(1),nb,nw)

   !--finished
   call closz(nendf)
   call closz(npend)
   call closz(noutp)
   call closz(nscr1)
   write(nsyso,'(/'' found'',i2,'' temperatures'')') itemp-1
   call timer(time)
   write(nsyso,'(69x,f8.1,''s''/1x,7(''**********''),&
     &''*******'')') time

   deallocate(egas)
   deallocate(sgas)

   return
   end subroutine gaspr

end module gaspm

