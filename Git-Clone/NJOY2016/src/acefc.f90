module acefc
   ! provides fast continuous options for acer
   use locale
   implicit none
   private

   !--Public routines
   public acetop,acefix

   !--Private global variables

   ! header parameters for ace continuous format
   character(13)::hz
   character(10)::hd
   character(10)::hm
   real(kr)::aw0,tz

   ! parameters for continuous nxs block
   integer::len2,izaid,nes,ntr,nr,ntrp,ntype,ndnf,nxsd(8)

   ! parameters for continuous jxs block
   integer::esz,nu,mtr,lqr,tyr,lsig,sig,land,and,ldlw,dlw,&
     gpd,mtrp,lsigp,sigp,landp,andp,ldlwp,dlwp,yp,fis,end,&
     iurpt,nud,dndat,ldnd,dnd,jxsd(2),ptype,ntro,ploct

   ! index of sections
   integer,parameter::nxcmax=500
   integer::nxc,mfs(nxcmax),mts(nxcmax),ncs(nxcmax)

   ! scratch units
   integer::mscr,iold,inew,nscr

   ! various bin sizes
   integer::nbina,nbinp,negn

   ! parameters obtained from the input files
   real(kr)::za,awr
   integer::mt19,mf1x(3)
   real(kr)::elast
   real(kr)::elim
   integer::ngmt,nned
   integer::mt16,mt455
   integer::mt5n,mt5p,mt5d,mt5t,mt5he3,mt5a
   integer::mt103,mt104,mt105,mt106,mt107
   integer::mpmin,mpmax,mdmin,mdmax,mtmin,mtmax,m3min,m3max,m4min,m4max

   ! record parameters for Type-2 binary files
   integer::ner,nbw

   ! sizes for equally probable bins
   integer::mcoars,npt

   ! File 6 parameters
   integer::ipnu,jpnu,nkk
   real(kr),allocatable,dimension(:)::e1456,p1456nu
   real(kr),allocatable,dimension(:)::p16456nu,p6456nu
   integer::nsix,n16

   ! particle production information (from common/ace8/)
   real(kr)::t201,t203,t204,t205,t206,t207
   integer,parameter::maxpr=300
   integer::nprod,kprod(maxpr),mprod(maxpr),iprod(maxpr),lprod(maxpr)

   ! incident particle
   real(kr)::awi
   integer::izai

   ! photon production reactions
   integer,parameter::maxpp=250
   integer::ntrpp,nf12s,mf12s(maxpp),nf16s,mf16s(maxpp)

   ! storage for reaction thresholds and gamma discontinuities
   real(kr),dimension(:),allocatable::ethr
   real(kr),dimension(:),allocatable::disc

   ! storage for group structure
   real(kr),dimension(:),allocatable::egn

   ! storage for gamma mt values
   real(kr),dimension(:),allocatable::gmt

   ! storage for ptleg data
   real(kr),dimension(:),allocatable::xat

   ! main container array for fast continuous data
   integer,parameter::nxss=20000000
   real(kr)::xss(nxss)

contains

   subroutine acetop(nendf,npend,ngend,nace,ndir,iprint,itype,mcnpx,&
     suff,hk,izn,awn,matd,tempd,newfor,iopp,ismooth,thin)
   !--------------------------------------------------------------------
   ! Prepare an ACE fast continuous file.
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use util   ! provides openz,mess,closz
   use endf   ! provides endf routines and variables
   ! externals
   integer::nendf,npend,ngend,nace,ndir,iprint,itype,matd,newfor,iopp,ismooth
   integer::mcnpx
   real(kr)::suff
   character(70)::hk
   integer::izn(16)
   real(kr)::awn(16)
   real(kr)::tempd
   real(kr)::thin(4)
   ! internals
   integer::nb,nw,nethr,nedis,itab,nscr2
   real(kr)::b(17)

   nxsd=0
   jxsd=0
   xss=0

   !--assign input units
   call openz(nendf,0)
   call openz(npend,0)
   call openz(ngend,0)

   !--determine what endf version is being used
   call repoz(nendf)
   call tpidio(nendf,0,0,b,nb,nw)
   call contio(nendf,0,0,b,nb,nw)
   call contio(nendf,0,0,b,nb,nw)
   if (n1h.ne.0) then
      iverf=4
   else if (n2h.eq.0) then
      iverf=5
   else
      iverf=6
   endif
   write(nsyso,'(/'' using endf-'',i1,'' format'')') iverf

   !--assign scratch files
   mscr=10
   mscr=iabs(mscr)
   if (npend.lt.0) mscr=-mscr
   call openz(mscr,1)
   iold=11
   inew=12
   nscr=13
   nscr2=14

   ! assign parameters for binary files
   ner=512
   nbw=8

   ! assign bin sizes for equally probable representations
   nbina=32
   nbinp=20

   !--prepare new files 1 and 2.
   nethr=0
   call first(nendf,npend,mscr,nethr,itab,matd,tempd,newfor,iopp)

   !--convert photon production representations.
   nedis=0
   call convr(nendf,npend,nscr2,0,nedis,nethr,matd)

   !--prepare new file 3 on unionized grid.
   call unionx(nendf,npend,mscr,matd,nedis,nethr,thin)

   !--prepare new files 4, 5, and 6.
   call topfil(nendf,mscr,matd,newfor)

   !--prepare new file 13 on unionized grid.
   if (ngmt.ne.0) call gamsum(npend,mscr,nscr2,matd,iopp)

   !--prepare new files 14 and 15 in equal probability bins.
   if (ngmt.ne.0) call gamout(ngend,nendf,mscr,nscr2,matd)

   !--terminate the scratch file.
   call amend(mscr,0)
   call atend(mscr,0)

   !--load ace data into memory.
   call acelod(mscr,nedis,suff,matd,tempd,newfor,mcnpx,ismooth)

   !--print ace file.
   if (iprint.gt.0) call aceprt(hk)

   !--write output file in desired mode.
   if (nace.ne.0) call aceout(itype,nace,ndir,hk,izn,awn,mcnpx)

   !--acetop is finished.
   deallocate(disc)
   deallocate(ethr)
   deallocate(gmt)
   if (allocated(xat)) deallocate(xat)
   if (allocated(egn)) deallocate(egn)
   call closz(nendf)
   call closz(npend)
   call closz(nscr2)
   call closz(ngend)
   call closz(mscr)
   return
   end subroutine acetop

   subroutine first(nendf,npend,nout,nethr,itab,matd,tempd,newfor,iopp)
   !--------------------------------------------------------------------
   ! Use tape id and descriptive data from PENDF tape.
   ! Convert MT452, 455, and 456 to tabular form if necessary.
   ! Remove other sections from MF1.
   ! Write new File 2 appropriate for no resonance parameters.
   ! Copy probability tables from MF2/Mt153 if found.
   !--------------------------------------------------------------------
   use mainio ! provides nsyso
   use util   ! provides repoz,closz,error
   use endf   ! provides endf routines and variables
   ! externals
   integer::nendf,npend,nout,nethr,matd,newfor,iopp
   real(kr)::tempd
   ! internals
   integer::nwscr,nb,nw,nxp,jscr,itab,nlib,nbt,nwt,nwd,nx
   integer::ncds,i,mfd,mtd,iflag,nprod3,nprodt,n
   integer::mf,mt,izr,ln,mtt,nc2,jethr,isix,nk,ik
   integer::izap,law,ll,nn,jj,ii,j,isort1,isort2,isave1
   integer::isave2,isave3,isave4,idis
   integer::nwtst,newnw
   real(kr)::delta,temp,y201,y203,y204,y205,y206,y207
   real(kr)::e,enext,y,xsthr,e201,e203,e204,e205,e206,e207
   real(kr)::b(350)
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::dict,dictp
   integer,parameter::ngmtmx=500
   integer,parameter::itape=11
   real(kr),parameter::etop=1.e10_kr
   real(kr),parameter::elo=1.e-5_kr
   real(kr),parameter::ehi=2.e7_kr
   real(kr),parameter::zero=0
   real(kr),parameter::ttol=1
   real(kr),parameter::err=.01e0_kr

   !--copy tape id and descriptive data from pendf tape.
   nethr=300
   allocate(ethr(nethr))
   allocate(gmt(ngmtmx))
   ngmt=0
   ntrpp=0
   nf12s=0
   nwscr=2*npage+50
   allocate(scr(nwscr))
   nsh=0
   nsc=0
   itab=0
   call repoz(nendf)
   call repoz(npend)
   call repoz(nout)
   call tpidio(npend,nout,0,scr,nb,nw)
   call tpidio(nendf,0,0,scr,nb,nw)
   call findf(matd,1,451,npend)
   call findf(matd,1,451,nendf)
   delta=ttol
   do while (delta.ge.ttol)
      call contio(npend,0,0,scr,nb,nw)
      if (math.ne.matd)&
        call error('first','desired temperature not found.',' ')
      za=scr(1)
      za=nint(za)
      awr=scr(2)
      nxp=n2h
      jscr=7
      if (iverf.ge.5) call contio(npend,0,0,scr(jscr),nb,nw)
      if (iverf.ge.5) jscr=jscr+6
      if (iverf.eq.6) call contio(npend,0,0,scr(jscr),nb,nw)
      if (iverf.eq.6) jscr=jscr+6
      awi=1
      izai=1
      if (iverf.eq.6) then
         awi=scr(13)
         nlib=nint(scr(17))
         izai=nlib/10
         if (izai.ne.1) write(nsyso,'(/&
           &'' charged incident particle:'',i6)') izai
      endif
      call hdatio(npend,0,0,scr(jscr),nb,nw)
      if (iverf.ge.5) nxp=n2h
      temp=scr(jscr)
      delta=abs(temp-tempd)
      if (delta.ge.ttol) call tomend(npend,0,0,scr)
   enddo
   call contio(0,nout,0,scr,nbt,nwt)
   if (iverf.ge.5) call contio(0,nout,0,scr(7),nbt,nwt)
   if (iverf.eq.6) call contio(0,nout,0,scr(13),nbt,nwt)
   call hdatio(0,nout,0,scr(jscr),nb,nw)
   nwd=nint(scr(jscr+4))
   do while (nb.ne.0)
      call moreio(npend,nout,0,scr,nb,nw)
   enddo
   call asend(nout,0)

   !--process dictionary from endf tape.
   call contio(nendf,0,0,scr,nb,nw)
   nx=n2h
   if (iverf.ge.5) call contio(nendf,0,0,scr,nb,nw)
   if (iverf.eq.6) call contio(nendf,0,0,scr,nb,nw)
   elim=ehi
   if (iverf.eq.6.and.scr(2).gt.ehi) elim=scr(2)
   call hdatio(nendf,0,0,scr,nb,nw)
   if (iverf.ge.5) nx=n2h
   do while (nb.ne.0)
      call moreio(nendf,0,0,scr,nb,nw)
   enddo
   ! read dictionary.
   nw=nx*6
   allocate(dict(nw))
   mt5n=1   !default is no neutron production in mt5
   mt5p=1   !  "           proton
   mt5d=1   !  "           deuteron
   mt5t=1   !  "           triton
   mt5he3=1 !  "           3He
   mt5a=1   !  "           alpha
   mt16=0
   mt455=0
   mt19=0
   mf1x(1)=0
   mf1x(2)=0
   mf1x(3)=0
   nsix=0
   nxc=1
   mfs(nxc)=1
   mts(nxc)=451
   ! save number of card images for descriptive data only for mt=451
   if (npend.gt.0) ncds=nwd+2
   if (npend.lt.0) ncds=nwd/17+2
   ncs(nxc)=ncds
   nw=nx
   if (nw.ne.0) then
      call dictio(nendf,0,0,dict,nb,nw)
      nprod=0
      do i=1,nw,6
         mfd=nint(dict(i+2))
         mtd=nint(dict(i+3))
         if (mfd.ge.3) then
            if (mtd.ge.875.and.mtd.le.891) mt16=1
            if (mfd.eq.5.and.mtd.eq.455) mt455=1
            if (mfd.ge.4.and.mfd.le.6.and.mtd.eq.19) mt19=1
            if (mfd.eq.6) nsix=nsix+1
            if (mfd.eq.12.and.iopp.ne.0) mf1x(1)=mf1x(1)+1
            if (mfd.eq.12.and.(mtd.lt.5.or.mtd.gt.600)) then
               nf12s=nf12s+1
               if (nf12s.gt.maxpp)&
                 call error('first','storage exceeded.',' ')
               mf12s(nf12s)=mtd
            endif
            if (mfd.eq.13.and.iopp.ne.0) mf1x(2)=mf1x(2)+1
            if (mfd.eq.15.and.iopp.ne.0) mf1x(3)=mf1x(3)+1
            if (iopp.ne.0.and.&
                ((mfd.eq.12.and.mtd.ne.460).or.(mfd.eq.13))) then
               ngmt=ngmt+1
               gmt(ngmt)=mfd*1000+mtd
            endif
            iflag=0
            if (mfd.eq.3.and.mtd.eq.102.and.izai.eq.1&
              .and.za+izai.lt.2004.1) iflag=1
            if (mfd.eq.4.and.mtd.eq.2) iflag=1
            if (mfd.eq.4.and.mtd.ge.600.and.mtd.le.648) iflag=1
            if (mfd.eq.4.and.mtd.ge.650.and.mtd.le.698) iflag=1
            if (mfd.eq.4.and.mtd.ge.700.and.mtd.le.748) iflag=1
            if (mfd.eq.4.and.mtd.ge.750.and.mtd.le.798) iflag=1
            if (mfd.eq.4.and.mtd.ge.800.and.mtd.le.848) iflag=1

            ! for the new format, check for particle production.
            ! note that acer will only include production from mts
            ! that appear in file 4 two-body reactions, file 6, and
            ! neutron mt102 with light targets.  thus, contributions
            ! from mf3 lr flags and any other mf3 sections with no
            ! distribution given will be omitted.  as a result, the
            ! acer particle production cross section may be less
            ! than the gas production cross sections in mt=203-207.

            if (iflag.ne.0.and.newfor.ne.0) then
               if (nprod+6.gt.maxpr) call error('first',&
                 'too many production items',' ')
               y201=0
               y203=0
               y204=0
               y205=0
               y206=0
               y207=0
               mt=mtd
               mf=mfd
               izr=0
               if (mt.eq.2) then
                  izr=nint(za+izai)-izai
               else if (mt.ge.50.and.mt.le.91) then
                  izr=nint(za+izai)-1
                  y201=1
               else if (mt.eq.102) then
                  izr=nint(za+izai)
               else if (mt.ge.600.and.mt.le.649) then
                  izr=nint(za+izai)-1001
                  y203=1
               else if (mt.ge.650.and.mt.le.699) then
                  izr=nint(za+izai)-1002
                  y204=1
               else if (mt.ge.700.and.mt.le.749) then
                  izr=nint(za+izai)-1003
                  y205=1
               else if (mt.ge.750.and.mt.le.799) then
                  izr=nint(za+izai)-2003
                  y206=1
               else if (mt.ge.800.and.mt.le.849) then
                  izr=nint(za+izai)-2004
                  y207=1
               endif
               if (izr.eq.1001) y203=y203+1
               if (izr.eq.1002) y204=y204+1
               if (izr.eq.1003) y205=y205+1
               if (izr.eq.2003) y206=y206+1
               if (izr.eq.2004) y207=y207+1
               if (izr.eq.4008) y207=y207+2
               if (y201.gt.zero.and.izai.ne.1) then
                  nprod=nprod+1
                  mprod(nprod)=mt
                  iprod(nprod)=1
                  kprod(nprod)=mf
                  lprod(nprod)=0
               else if (y203.gt.zero.and.izai.ne.1001) then
                  nprod=nprod+1
                  mprod(nprod)=mt
                  iprod(nprod)=1001
                  kprod(nprod)=mf
                  lprod(nprod)=0
               else if (y204.gt.zero.and.izai.ne.1002) then
                  nprod=nprod+1
                  mprod(nprod)=mt
                  iprod(nprod)=1002
                  kprod(nprod)=mf
                  lprod(nprod)=0
               else if (y205.gt.zero.and.izai.ne.1003) then
                  nprod=nprod+1
                  mprod(nprod)=mt
                  iprod(nprod)=1003
                  kprod(nprod)=mf
                  lprod(nprod)=0
               else if (y206.gt.zero.and.izai.ne.2003) then
                  nprod=nprod+1
                  mprod(nprod)=mt
                  iprod(nprod)=2003
                  kprod(nprod)=mf
                  lprod(nprod)=0
               else if (y207.gt.zero.and.izai.ne.2004) then
                  nprod=nprod+1
                  mprod(nprod)=mt
                  iprod(nprod)=2004
                  kprod(nprod)=mf
                  lprod(nprod)=0
               endif
            endif
         endif
      enddo
   endif
   nw=nxp*6
   allocate(dictp(nw))
   call dictio(npend,0,0,dictp,nb,nxp)
   if (ngmt.gt.ngmtmx) call error('first','storage exceeded.',' ')
   call tosend(nendf,0,0,scr)

   !--process fission yield sections from endf tape.
   !--remove all other sections.
   mfh=1
   do while (mfh.ne.0)
      call contio(nendf,0,0,b,nb,nw)
      if (mfh.eq.1) then
         if (mth.ne.452.and.mth.ne.455.and.mth.ne.456) then
            call tosend(nendf,0,0,scr)
         else
            ln=l2h
            mt=mth
            call tabize(mt,ln,nendf,nout,itape,err,b)
            itab=1
         endif
      endif
   enddo
   call afend(nout,0)

   !--write a file 2 with no resonance parameters.
   mfh=2
   mth=151
   mtt=151
   scr(1)=za
   scr(2)=awr
   scr(3)=0
   scr(4)=0
   scr(5)=1
   scr(6)=0
   call contio(0,nout,0,scr,nb,nw)
   scr(2)=1
   call contio(0,nout,0,scr,nb,nw)
   scr(1)=elo
   scr(2)=ehi
   scr(5)=0
   call contio(0,nout,0,scr,nb,nw)
   scr(1)=0
   scr(2)=0
   call contio(0,nout,0,scr,nb,nw)
   nc2=4
   call asend(nout,0)

   !--copy probability tables if found
   mfh=2
   do while (mfh.ne.0)
      call contio(npend,0,0,scr,nb,nw)
      if (mfh.gt.0) then
         if (mth.eq.153) then
            call contio(0,nout,0,scr,nb,nw)
            call tosend(npend,nout,0,scr)
         else
            call tosend(npend,0,0,scr)
         endif
      endif
   enddo

   !--finished with file 2.
   call afend(nout,0)
   nxc=nxc+1
   mfs(nxc)=2
   mts(nxc)=mtt
   ncs(nxc)=nc2
   nxc=nxc+1
   mfs(nxc)=3
   mts(nxc)=1
   ncs(nxc)=0

   !--make ordered list of all thresholds in file 3.
   !--find particle production thresholds.
   call findf(matd,3,2,npend)
   call tosend(npend,0,0,scr)
   t201=etop
   t203=etop
   t204=etop
   t205=etop
   t206=etop
   t207=etop
   jethr=0
   mfh=1
   do while (mfh.ne.0)
      call contio(npend,0,0,scr,nb,nw)
      if (mfh.ne.0) then
         e=0
         call gety1(e,enext,idis,y,npend,scr)
         ethr(1+jethr)=enext
         jethr=jethr+1
         if (nprod.gt.0) then
            do i=1,nprod
               if (mprod(i).eq.mth) then
                  if (iprod(i).eq.1.and.enext.lt.t201) t201=enext
                  if (iprod(i).eq.1001.and.enext.lt.t203) t203=enext
                  if (iprod(i).eq.1002.and.enext.lt.t204) t204=enext
                  if (iprod(i).eq.1003.and.enext.lt.t205) t205=enext
                  if (iprod(i).eq.2003.and.enext.lt.t206) t206=enext
                  if (iprod(i).eq.2004.and.enext.lt.t207) t207=enext
               endif
            enddo
         endif
         call tosend(npend,0,0,scr)
         if (jethr.ge.nethr) then
            call error('first','too many thresholds',' ')
         endif
      endif
   enddo
   if (jethr.gt.1) call aordr(jethr,nethr,ethr)
   nethr=jethr
   call skiprz(npend,-2)

   !--for new format...
   if (newfor.ne.0) then

      !--check particle productions in mf6.
      nprod3=nprod
      if (nsix.gt.0) then
         call findf(matd,6,0,nendf)
         isix=0
         newnw=0
         do while (isix.lt.nsix)
            isix=isix+1
            call contio(nendf,0,0,scr,nb,nw)
            nk=n1h
            do ik=1,nk
               call tab1io(nendf,0,0,scr,nb,nw)
               if (mth.eq.5) then
                  !--test for light particle production; 0=yes
                  if (abs(c1h-1).lt.0.0001) mt5n=0
                  if (abs(c1h-1001).lt.0.0001) mt5p=0
                  if (abs(c1h-1002).lt.0.0001) mt5d=0
                  if (abs(c1h-1003).lt.0.0001) mt5t=0
                  if (abs(c1h-2003).lt.0.0001) mt5he3=0
                  if (abs(c1h-2004).lt.0.0001) mt5a=0
               endif
               law=l2h
               nwtst=6+2*(nint(scr(5))+nint(scr(6)))
               newnw=max(nwscr,nwtst,newnw)
               do while (nb.ne.0)
                  call moreio(nendf,0,0,scr,nb,nw)
               enddo
               call skip6(nendf,0,0,scr,law)
            enddo
            call contio(nendf,0,0,scr,nb,nw) !read eos record
         enddo
         if (newnw.gt.nwscr) then
            deallocate(scr)
            nwscr=newnw
            allocate(scr(nwscr))
         endif
         call repoz(nendf)
         call findf(matd,6,0,nendf)
      endif
      isix=0
      do while (isix.lt.nsix)
         isix=isix+1
         call contio(nendf,0,0,scr,nb,nw)
         nk=n1h
         mtd=mth
         call findf(matd,3,mtd,npend)
         call contio(npend,0,0,b,nb,nw)
         e=0
         call gety1(e,enext,idis,y,npend,b)
         xsthr=enext
         do ik=1,nk
            call tab1io(nendf,0,0,scr,nb,nw)
            izap=nint(c1h)
            law=l2h
            ll=1+nw
            do while (nb.ne.0)
               call moreio(nendf,0,0,scr(ll),nb,nw)
               ll=ll+nw
            enddo
            nn=nint(scr(6))
            jj=0
            do ii=1,nn
               if (jj.eq.0.and.scr(8+2*ii).gt.zero) jj=ii
            enddo
            if (jj.gt.1) jj=jj-1
            if (nprod+6.gt.maxpr) call error('first',&
              'too many production items',' ')
            if (izap.eq.1.and.izai.ne.1) then
               e201=scr(7+2*jj)
               if (xsthr.gt.e201) e201=xsthr
               if (e201.lt.t201) t201=e201
               nprod=nprod+1
               mprod(nprod)=mth
               iprod(nprod)=1
               kprod(nprod)=6
               lprod(nprod)=ik
            else if (izap.eq.1001.and.izai.ne.1001) then
               e203=scr(7+2*jj)
               if (xsthr.gt.e203) e203=xsthr
               if (e203.lt.t203) t203=e203
               nprod=nprod+1
               mprod(nprod)=mth
               iprod(nprod)=1001
               kprod(nprod)=6
               lprod(nprod)=ik
            else if (izap.eq.1002.and.izai.ne.1002) then
               e204=scr(7+2*jj)
               if (xsthr.gt.e204) e204=xsthr
               if (e204.lt.t204) t204=e204
               nprod=nprod+1
               mprod(nprod)=mth
               iprod(nprod)=1002
               kprod(nprod)=6
               lprod(nprod)=ik
            else if (izap.eq.1003.and.izai.ne.1003) then
               e205=scr(7+2*jj)
               if (xsthr.gt.e205) e205=xsthr
               if (e205.lt.t205) t205=e205
               nprod=nprod+1
               mprod(nprod)=mth
               iprod(nprod)=1003
               kprod(nprod)=6
               lprod(nprod)=ik
            else if (izap.eq.2003.and.izai.ne.2003) then
               e206=scr(7+2*jj)
               if (xsthr.gt.e206) e206=xsthr
               if (e206.lt.t206) t206=e206
               nprod=nprod+1
               mprod(nprod)=mth
               iprod(nprod)=2003
               kprod(nprod)=6
               lprod(nprod)=ik
            else if (izap.eq.2004.and.izai.ne.2004) then
               e207=scr(7+2*jj)
               if (xsthr.gt.e207) e207=xsthr
               if (e207.lt.t207) t207=e207
               nprod=nprod+1
               mprod(nprod)=mth
               iprod(nprod)=2004
               kprod(nprod)=6
               lprod(nprod)=ik
            endif
            call skip6(nendf,0,0,scr,law)
         enddo
         call tosend(nendf,0,0,scr)
      enddo

      !--sort the production contributions
      if (nprod.gt.0) then
         do i=1,nprod-1
            do j=i,nprod
               isort1=100000*mprod(i)+10*iprod(i)+lprod(i)
               isort2=100000*mprod(j)+10*iprod(j)+lprod(j)
               if (isort1.gt.isort2) then
                  isave1=iprod(j)
                  isave2=mprod(j)
                  isave3=kprod(j)
                  isave4=lprod(j)
                  iprod(j)=iprod(i)
                  mprod(j)=mprod(i)
                  kprod(j)=kprod(i)
                  lprod(j)=lprod(i)
                  iprod(i)=isave1
                  mprod(i)=isave2
                  kprod(i)=isave3
                  lprod(i)=isave4
               endif
            enddo
         enddo

         ! If light recoil particle production (p,d,t,3He,a) was
         ! found in file 6 it may be redundant with light particle
         ! production that can be inferred from mf=3, mt=102 to
         ! 107.  Therefore, check for redundant mt values in
         ! mprod(i) & mprod(i+1), redundant zap values in iprod(i)
         ! & iprod(i+1), but different mf values in kprod(i) &
         ! kprod(i+1). If found, delete the ith entry, move i+1
         ! and later entries down one array location and decrement
         ! nprod.
         if (nprod.ne.nprod3.and.nprod.gt.1) then
            nprodt=2
            do while (nprodt.le.nprod)
               if (mprod(nprodt).eq.mprod(nprodt-1) .and.&
                 iprod(nprodt).eq.iprod(nprodt-1) .and.&
                 kprod(nprodt).ne.kprod(nprodt-1)) then
                  nprod=nprod-1
                  if (nprodt.le.nprod) then
                     nprodt=nprodt-1
                     do n=nprodt,nprod
                        iprod(n)=iprod(n+1)
                        mprod(n)=mprod(n+1)
                        kprod(n)=kprod(n+1)
                        lprod(n)=lprod(n+1)
                     enddo
                  endif
               endif
               nprodt=nprodt+1
            enddo
         endif

      endif
   endif

   !--finished with first.
   deallocate(dict)
   deallocate(dictp)
   deallocate(scr)
   return
   end subroutine first

   subroutine tabize(mti,ln,nin,nout,itape,err,b)
   !--------------------------------------------------------------------
   ! Convert polynomial expansions to tabulated representations,
   ! or convert other interpolation formulas to linear-linear,
   ! if necessary, and write TAB1 records on nout.
   !--------------------------------------------------------------------
   use util ! provides openz,loada,finda,sigfig
   use endf ! provides endf routines and variables
   ! externals
   integer::mti,ln,nin,nout,itape
   real(kr)::err,b(*)
   ! internals
   integer::nb,nw,i,nonlin,int,ia,ir,in,nbt,ibase,id,ne,idone,nen
   integer::nwm,nwscr
   integer::j
   real(kr)::dy,xm,yy,ym,test
   real(kr),dimension(:),allocatable::scr,scr2
   real(kr),dimension(:),allocatable::buf
   real(kr),dimension(:),allocatable::x,y
   integer,parameter::nbuf=1000
   integer,parameter::ndim=100
   integer,parameter::maxscr=2000
   real(kr),parameter::xmax=2.e7_kr
   real(kr),parameter::errlim=1.e-5_kr
   real(kr),parameter::stpmax=1.25e0_kr
   integer,parameter::nc=4

   !--initialize
   if (allocated(scr)) deallocate(scr)
   allocate(scr(maxscr))
   nxc=nxc+1
   mfs(nxc)=1
   mts(nxc)=mti
   ncs(nxc)=0
   nb=0
   nw=6

   !--copy polynomials (ln=1), don't convert them.
   if (ln.eq.1) then
      call contio(0,nout,0,b,nb,nw)
      call listio(nin,nout,0,scr,nb,nw)
      ncs(nxc)=ncs(nxc)+1+nint((scr(5)+5)/6)
      if (mti.eq.455) then
         do while (nb.ne.0)
            call moreio(nin,nout,0,scr,nb,nw)
         enddo
         call listio(nin,nout,0,scr,nb,nw)
         ncs(nxc)=ncs(nxc)+1+nint((scr(5)+5)/6)
      endif
      call tosend(nin,nout,0,scr)
      deallocate(scr)
      return
   endif

   !--this section contains tabulated data (ln=2).
   b(4)=2
   ncs(nxc)=ncs(nxc)+1
   call contio(0,nout,0,b,nb,nw)
   ! copy precursor record for mt 455.
   if (mti.eq.455) then
      call listio(nin,nout,0,scr,nb,nw)
      ncs(nxc)=ncs(nxc)+1+nint((scr(5)+5)/6)
      do while (nb.ne.0)
         call moreio(nin,nout,0,scr,nb,nw)
      enddo
   endif
   call tab1io(nin,0,0,scr,nb,nw)
   nr=nint(scr(5))
   i=0
   nonlin=0
   do while (i.lt.nr.and.nonlin.eq.0)
      i=i+1
      int=nint(scr(6+2*i))
      if (int.gt.2) nonlin=1
   enddo

   !--the section is already linear.  copy and return.
   if (nonlin.eq.0) then
      call tab1io(0,nout,0,scr,nb,nw)
      ncs(nxc)=ncs(nxc)+1+nint((scr(5)+2)/3)+nint((scr(6)+2)/3)
      if (mti.eq.456) then
         ipnu=nint(scr(6))
         if (allocated(e1456)) then
            deallocate(e1456)
            deallocate(p1456nu)
         endif
         allocate(e1456(ipnu))
         allocate(p1456nu(ipnu))
         ibase=7+2*nr
         i=0
         do j=ibase,nw,2
            i=i+1
            e1456(i)=scr(j)
            p1456nu(i)=scr(j+1)
         enddo
      endif
      do while (nb.ne.0)
         call moreio(nin,nout,0,scr,nb,nw)
         if (mti.eq.456) then
            do j=1,nw,2
               i=i+1
               e1456(i)=scr(j)
               p1456nu(i)=scr(j+1)
            enddo
         endif
      enddo
      if (allocated(e1456)) then
         deallocate(e1456)
         deallocate(p1456nu)
      endif
      call contio(nin,nout,0,scr,nb,nw)
      deallocate(scr)
      return
   endif

   !--linearize the tabulated data.
   call openz(-itape,1)
   nwscr=npage+50
   allocate(scr2(nwscr))
   allocate(buf(nbuf))
   allocate(x(ndim))
   allocate(y(ndim))
   ia=1+6+2*nr
   ir=1
   in=1
   nbt=nint(scr(5+2*ir))
   int=nint(scr(6+2*ir))
   ibase=6+2*nr
   x(2)=scr(ia)
   y(2)=scr(ia+1)
   scr2(1)=0
   scr2(2)=0
   scr2(3)=0
   scr2(4)=0
   scr2(5)=1
   scr2(6)=0
   scr2(7)=0
   scr2(8)=2
   id=8
   ne=0
   idone=0
   do while (idone.eq.0)
      i=2
      ia=ia+2
      if (ia.gt.nw.and.nb.eq.0) then
         idone=1
      else
         in=in+1
         if (in.gt.nbt) then
            ir=ir+1
            nbt=nint(scr(5+2*ir))
            int=nint(scr(6+2*ir))
         endif
         x(1)=scr(ia)
         y(1)=scr(ia+1)
         if (ia.ge.nw.and.nb.ne.0) then
            call moreio(nin,0,0,scr(1+ibase),nb,nwm)
            nw=ibase+nwm
         endif
         do while (i.gt.1)
            dy=0
            xm=(x(i)+x(i-1))/2
            xm=sigfig(xm,7,0)
            if (xm.gt.x(i).and.xm.lt.x(i-1)) then
               call terp1(x(i),y(i),x(i-1),y(i-1),xm,yy,int)
               call terp1(x(i),y(i),x(i-1),y(i-1),xm,ym,2)
               dy=abs(yy-ym)
            endif
            if (i.eq.ndim) dy=0
            test=err*ym
            if (test.lt.errlim) test=errlim
            if (dy.gt.test) then
            ! not converged.
            ! add the midpoint to the stack and continue.
               i=i+1
               x(i)=x(i-1)
               y(i)=y(i-1)
               x(i-1)=xm
               y(i-1)=ym
            else
               ! converged.
               ! use the top point off the stack.
               ne=ne+1
               b(1)=x(i)
               b(2)=y(i)
               call loada(ne,b,2,itape,buf,nbuf)
               i=i-1
            endif
         enddo
         if (x(1).eq.xmax) idone=1
         x(2)=x(1)
         y(2)=y(1)
      endif
   enddo
   ne=ne+1
   nen=-ne
   b(1)=x(1)
   b(2)=y(1)
   call loada(nen,b,2,itape,buf,nbuf)

   !--finished with linearization.
   !--write new tab1-formatted mt.
   scr2(6)=ne
   scr2(7)=ne
   if (mti.eq.456) then
      ipnu=ne
      if (allocated(e1456)) then
         deallocate(e1456)
         deallocate(p1456nu)
      endif
      allocate(e1456(ipnu))
      allocate(p1456nu(ipnu))
   endif
   do i=1,ne
      call finda(i,b,2,itape,buf,nbuf)
      id=id+1
      scr2(id)=b(1)
      id=id+1
      scr2(id)=b(2)
      if (mti.eq.456) then
         e1456(i)=b(1)
         p1456nu(i)=b(2)
      endif
      if (id.ge.npage.or.i.eq.ne) then
         if (ibase.ne.0) then
            call tab1io(0,nout,0,scr2,nb,id)
            id=0
            ibase=0
         else
            call moreio(0,nout,0,scr2,nb,id)
            id=0
         endif
      endif
   enddo
   if (allocated(e1456)) then
      deallocate(e1456)
      deallocate(p1456nu)
   endif
   call contio(nin,nout,0,scr,nb,nw)
   call closz(-itape)
   deallocate(y)
   deallocate(x)
   deallocate(buf)
   deallocate(scr2)
   deallocate(scr)

   !--finished
   return
   end subroutine tabize

   subroutine unionx(nendf,nin,nout,matd,nedis,nethr,thin)
   !--------------------------------------------------------------------
   ! Write out all File 3 reactions on a unionized energy grid.
   ! The grid of MT1 is used, therefore the routine will work on
   ! non-unionized input tapes at reduced accuracy.
   ! Redundant reactions and ratio quantities are deleted.
   ! Discontinuities are removed by moving first point down slightly.
   !--------------------------------------------------------------------
   use util ! provides openz,repoz,loada,finda,sigfig,error
   use endf ! provides endf routines and variables
   use mainio ! provides nsyso
   ! externals
   integer::nendf,nin,nout,matd,nedis,nethr
   real(kr)::thin(4)
   ! internals
   integer::nws,nwscr,ithopt,iwtt,npts,iskp,j,mtcomp
   integer::nb,nw,iinel,iedis,jethr,lt,lr
   integer::i,jt,ne,isave,iter,idone,nsave,ilast,inext
   integer::limit,ii,k,nc,nen,ll,nee,iee,mtn,itest
   integer::kbase,np,ifrst,it,ibase,idis,nold,nxcs
   integer::nins
   real(kr)::rsigz,e1,e2,e,enext,y,test,egl,egh,egd
   real(kr)::eet,edl,edh,s,q,el,rsum,r,rl,csum,rmin,ee
   real(kr)::eg,es,rs,cl,cs,rtot,ctot,sum,cap
   real(kr)::en,rn,rinc,rlin,cn,cinc,clin,er
   real(kr)::egrid,thresh,thrx,errr,ethrr
   character(52)::messs
   real(kr)::b(6),c(4)
   integer::jsave(10)
   real(kr)::esave(10),rsave(10,10),csave(10,10)
   real(kr),dimension(:),allocatable::scr,scr2
   real(kr),dimension(:),allocatable::buf,bufn
   real(kr),parameter::fac=2.e0_kr
   real(kr),parameter::rfact=.0025e0_kr
   real(kr),parameter::errm=.95e0_kr
   real(kr),parameter::eps=1.e-10_kr
   real(kr),parameter::etop=1.e10_kr
   real(kr),parameter::elow=1.e-5_kr
   real(kr),parameter::step=1.2e0_kr
   integer,parameter::miter=9
   integer,parameter::nbuf=1000
   real(kr),parameter::zero=0
   real(kr),parameter::one=1

   !--initialize
   inew=iabs(inew)
   nscr=iabs(nscr)
   if (nin.lt.0) nscr=-nscr
   call openz(-inew,1)
   call openz(-iold,1)
   call openz(nscr,1)
   call repoz(nscr)
   nxcs=nxc
   nws=2*npage+50
   allocate(scr(nws))
   nwscr=npage+50
   allocate(scr2(nwscr))
   allocate(buf(nbuf))
   allocate(bufn(nbuf))
   ithopt=nint(thin(4))
   if (ithopt.eq.2) iwtt=nint(thin(1))
   if (ithopt.eq.2) npts=nint(thin(2))
   if (ithopt.eq.2) rsigz=thin(3)
   iskp=0
   if (ithopt.eq.1) iskp=nint(thin(3))
   if (ithopt.eq.1) e1=thin(1)
   if (ithopt.eq.1) e2=thin(2)
   j=0
   c(2)=0
   c(3)=0
   ethrr=0
   mt103=0
   mt104=0
   mt105=0
   mt106=0
   mt107=0
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

   !--check for mt4 overlap in probability tables
   mtcomp=0
   if (izai.eq.1) then
      call findf(matd,2,151,nin)
      mfh=2
      do while (mfh.eq.2.and.mth.ne.153)
         call tosend(nin,0,0,scr)
         call contio(nin,0,0,scr,nb,nw)
      enddo
      if (mfh.ne.0) then
         !--new competition flag (can only be -1, 51, 91 or 4)
         iinel=l1h
         !--continue but check for file with old competition flag
         call listio(nin,0,0,scr,nb,nw)
         if (iinel.eq.zero.and.l2h.gt.zero) iinel=mod(l2h,1000)
         if (iinel.eq.4) mtcomp=4
      endif
   endif

   !--watch for discontinuities and thresholds
   iedis=1
   egl=etop
   egh=etop
   if (iedis.le.nedis) then
      egd=disc(iedis)
      egl=sigfig(egd,7,-1)
      egh=sigfig(egd,7,+1)
   endif
   jethr=1
   eet=0
   if (jethr.le.nethr) eet=ethr(jethr)
   edl=0
   edh=0

   !--for incident neutrons:
   !--store energy grid of the total cross section
   !--with thinning, if desired
   if (izai.eq.1) then
      call findf(matd,3,1,nin)
      call contio(nin,0,0,b,nb,nw)
      e=0
      call gety1(e,enext,idis,y,nin,scr)
      test=elow+elow/1000
      if (enext.gt.test) enext=elow
      s=scr(1)
      q=scr(2)
      lt=nint(scr(3))
      lr=nint(scr(4))
      el=0
      rsum=0
      r=0
      rl=0
      i=0
      jt=1
      do while (jt.gt.0)
         i=i+1
         e=enext
         call gety1(e,enext,idis,y,nin,scr)
         if (e.ge.(1-eps)*eet) then
            i=iskp
            jethr=jethr+1
            eet=0
            if (jethr.le.nethr) eet=ethr(jethr)
         endif
         if (abs(e-egl).lt.eps*egl) then
            enext=egh
            i=iskp
         else if (abs(e-egh).lt.eps*egh) then
            i=iskp
            egl=etop
            egh=etop
            iedis=iedis+1
            if (iedis.le.nedis) then
               egd=disc(iedis)
               egl=sigfig(egd,7,-1)
               egh=sigfig(egd,7,+1)
            endif
            if (idis.ne.0) then
               edl=sigfig(enext,7,-1)
               edh=sigfig(enext,7,+1)
               enext=edl
            endif
         else if (abs(e-edl).lt.eps*edl) then
            enext=edh
            i=iskp
         else if (abs(e-edh).lt.eps*edh) then
            i=iskp
            edl=0
            edh=0
            if (idis.ne.0) then
               edl=sigfig(enext,7,-1)
               edh=sigfig(enext,7,+1)
               enext=edl
            endif
         else if (enext.ge.(1-eps)*egl) then
            enext=egl
         endif
         c(1)=e
         c(2)=y
         if (ithopt.eq.2) then
            r=y*awt(iwtt,e)/(1+y/rsigz)
            if (j.gt.1) rsum=rsum+(e-el)*(r+rl)/2
            c(3)=r
            rl=r
         endif
         el=e
         itest=0
         if (iskp.le.0.or.i.ge.iskp) then
             itest=1
         endif
         if (itest.eq.0.and.ithopt.eq.1) then
             if (e.lt.(1-eps)*e1.or.e.gt.(1+eps)*e2) then
                itest=1
             endif
         endif
         if (itest.eq.1) then
            i=0
            j=j+1
            jt=j
            test=etop-etop/100
            if (enext.gt.test) jt=-jt
            call loada(jt,c,4,iold,buf,nbuf)
         endif
      enddo
      ne=j

      !--for integral thinning, get the
      !--corresponding capture cross sections
      if (ithopt.eq.2) then
         call findf(matd,3,102,nin)
         call contio(nin,0,0,b,nb,nw)
         e=0
         call gety1(e,enext,idis,y,nin,scr)
         j=0
         csum=0
         el=elow
         rl=0
         jt=1
         do while (jt.gt.0)
            j=j+1
            call finda(j,c,4,iold,buf,nbuf)
            e=c(1)
            call gety1(e,enext,idis,y,nin,scr)
            r=y*awt(iwtt,e)/(1+c(2)/rsigz)
            if (j.gt.1) csum=csum+(e-el)*(r+rl)/2
            c(4)=r
            rl=r
            el=e
            jt=j
            if (j.eq.ne) jt=-j
            call loada(jt,c,4,inew,bufn,nbuf)
         enddo
         isave=inew
         inew=iold
         iold=isave
         rmin=rfact*rsum/npts
         write(nsyso,'(/&
           &'' original grid='',i6,'' with integrals '',1p,2e12.4)')&
           ne,rsum,csum

         !--generate the thinned grid using the integral method
         !--on the total and capture cross sections
         iter=1
         jethr=1
         ee=0
         if (jethr.le.nethr) ee=ethr(jethr)
         iedis=1
         egl=etop
         egh=etop
         if (iedis.le.nedis) then
            eg=disc(iedis)
            egl=sigfig(eg,7,-1)
            egh=sigfig(eg,7,+1)
         endif
         idone=0
         do while (idone.eq.0)
            i=1
            j=1
            nsave=1
            call finda(i,c,4,iold,buf,nbuf)
            call loada(j,c,4,inew,bufn,nbuf)
            el=c(1)
            es=c(1)
            rl=c(3)
            rs=c(3)
            cl=c(4)
            cs=c(4)
            rsum=0
            csum=0
            rtot=0
            ctot=0
            sum=0
            cap=0
            errr=0
            jt=1
            do while (jt.gt.0)
               i=i+1
               call finda(i,c,4,iold,buf,nbuf)
               en=c(1)
               rn=c(3)
               rinc=(rn+rl)*(en-el)/2
               rsum=rsum+rinc
               rlin=(rn+rs)*(en-es)/2
               cn=c(4)
               cinc=(cn+cl)*(en-el)/2
               csum=csum+cinc
               clin=(cn+cs)*(en-es)/2
               if (i.eq.nsave*ne/10) then
                  rsave(nsave,iter+1)=sum+rlin
                  csave(nsave,iter+1)=cap+clin
                  jsave(nsave)=j
                  if (iter.eq.1) then
                     esave(nsave)=c(1)
                     rsave(nsave,1)=rtot+rsum
                     csave(nsave,1)=ctot+csum
                  endif
                  test=100*(rsave(nsave,iter+1)-rsave(nsave,1))&
                    /rsave(nsave,1)
                  if (test.gt.errr) errr=test
                  test=100*(csave(nsave,iter+1)-csave(nsave,1))&
                    /csave(nsave,1)
                  if (test.gt.errr) errr=test
                  nsave=nsave+1
               endif
               ! check for gamma discontinuities
               ilast=0
               inext=1
               if (abs(ee-egl).ge.eps*egl) then
                  if (abs(ee-egh).lt.eps*egh) then
                     iedis=jethr+1
                     egl=etop
                     egh=etop
                     if (iedis.le.nedis) then
                        eg=disc(iedis)
                        egl=sigfig(eg,7,-1)
                        egh=sigfig(eg,7,+1)
                     endif
                  else
                     ! check for thresholds
                     if (ee.ne.zero.and.en.ge.(1-eps)*ee) then
                        jethr=jethr+1
                        ee=0
                        if (jethr.le.nethr) ee=ethr(jethr)
                     else
                        ! check for convergence
                        if (rsum.ge.rmin.or.rlin.ge.rmin) then
                           ilast=1
                        else
                           if (abs(rsum-rlin).ge.rlin/10&
                             .or.abs(csum-clin).ge.clin/10) then
                              ilast=1
                           else
                              if (i.lt.ne.and.en.lt.es+es/20) then
                                 cl=cn
                                 rl=rn
                                 el=en
                                 inext=0
                              endif
                           endif
                        endif
                     endif
                  endif
               endif
               ! accept last point (el)
               if (ilast.eq.1) then
                  if (abs(el-es).ge.eps*es) then
                     sum=sum+(rl+rs)*(el-es)/2
                     rtot=rtot+rsum-rinc
                     rsum=rinc
                     cap=cap+(cl+cs)*(el-es)/2
                     ctot=ctot+csum-cinc
                     csum=cinc
                     c(1)=el
                     c(3)=rl
                     c(4)=cl
                  else
                     sum=sum+rlin
                     rtot=rtot+rsum
                     rsum=0
                     cap=cap+clin
                     ctot=ctot+csum
                     csum=0
                     c(1)=en
                     c(3)=rn
                     c(4)=cn
                  endif
               endif
               ! accept next point (en)
               if (inext.eq.1) then
                  sum=sum+rlin
                  rtot=rtot+rsum
                  rsum=0
                  cap=cap+clin
                  ctot=ctot+csum
                  csum=0
                  c(1)=en
                  c(3)=rn
                  c(4)=cn
               endif
               ! save this energy point
               if (ilast.eq.1.or.inext.eq.1) then
                  j=j+1
                  jt=j
                  if (i.eq.ne) jt=-j
                  call loada(jt,c,4,inew,bufn,nbuf)
                  es=c(1)
                  rs=c(3)
                  cs=c(4)
                  cl=cn
                  rl=rn
                  el=en
               endif
            enddo

            !--print iteration and check for end
            write(nsyso,'(/&
              &'' new grid='',i6,'' with integrals '',1p,2e12.4)')&
              j,sum,cap
            if (iter.eq.1.and.j.lt.3*npts/4) rmin=rmin/(fac**6)
            if (iter.eq.1.and.j.lt.3*npts/4) j=npts+1
            if (j.le.npts) then
               idone=1
            else
               if (j.eq.ne) then
                  idone=1
               else
                  if (iter.eq.miter) then
                     idone=1
                  else
                     if (errr.gt.errm) then
                        idone=1
                     else
                        jethr=1
                        ee=0
                        if (jethr.le.nethr) ee=ethr(jethr)
                        rmin=fac*rmin
                        iter=iter+1
                     endif
                  endif
               endif
            endif
         enddo
         isave=iold
         iold=inew
         inew=isave
         e=c(1)

         !--print detailed table of integrals
         limit=iter+1
         nsave=nsave-1
         do i=2,nsave
            ii=nsave-i+2
            jsave(ii)=jsave(ii)-jsave(ii-1)
            do k=1,limit
               rsave(ii,k)=rsave(ii,k)-rsave(ii-1,k)
               csave(ii,k)=csave(ii,k)-csave(ii-1,k)
            enddo
         enddo
         write(nsyso,'(/'' total'')')
         do i=1,nsave
            do k=2,limit
               if (rsave(i,1).ne.zero) then
                  rsave(i,k)=100*(rsave(i,k)-rsave(i,1))/rsave(i,1)
               endif
            enddo
            write(nsyso,'(1x,1p,2e12.4,0p,9f6.1)')&
              esave(i),(rsave(i,k),k=1,limit)
         enddo
         write(nsyso,'(/'' capture'')')
         do i=1,nsave
            do k=2,limit
               if (csave(i,1).ne.zero) then
                  csave(i,k)=100*(csave(i,k)-csave(i,1))/csave(i,1)
               endif
            enddo
            write(nsyso,'(1x,1p,2e12.4,0p,9f6.1)')&
              esave(i),(csave(i,k),k=1,limit)
         enddo
         write(nsyso,'(/1x,10i6)') (jsave(i),i=1,nsave)
      endif
      nc=4
      nen=j
      elast=e
      call findf(matd,3,2,nin)

   !--for incident charged particles:
   !--generate a union grid for incident charged particles
   !--first get the union grid of file 3
   else
      call findf(matd,3,0,nin)
      lt=0
      lr=0
      nold=0
      mfh=3
      do while (mfh.ne.0)
         if (mth.eq.2) then
            call contio(nin,0,0,b,nb,nw)
         else
            call contio(nin,0,0,scr,nb,nw)
         endif
         if (mfh.ne.0) then
            k=0
            j=0
            if (nold.gt.0) then
               k=k+1
               call finda(k,c,2,iold,bufn,nbuf)
               eg=c(1)
            else
               eg=etop
            endif
            er=0
            call gety1(er,enext,idis,y,nin,scr)
            er=enext
            jt=1
            do while (jt.gt.0)
               if (er.gt.eg) then
                  j=j+1
                  jt=j
                  test=etop-etop/100
                  if (enext.gt.test.and.k.eq.nold) jt=-jt
                  c(1)=eg
                  c(2)=0
                  call loada(jt,c,2,inew,buf,nbuf)
                  if (k.lt.nold) then
                     k=k+1
                     call finda(k,c,2,iold,bufn,nbuf)
                     eg=c(1)
                  else
                     eg=etop
                  endif
               else
                  call gety1(er,enext,idis,y,nin,scr)
                  j=j+1
                  jt=j
                  test=etop-etop/100
                  if (enext.gt.test.and.k.eq.nold) jt=-jt
                  c(1)=er
                  c(2)=0
                  call loada(jt,c,2,inew,buf,nbuf)
                  if (er.eq.eg.and.k.lt.nold) then
                     k=k+1
                     call finda(k,c,2,iold,bufn,nbuf)
                     eg=c(1)
                  endif
                  er=enext
               endif
            enddo
            call tosend(nin,0,0,scr)
            isave=iold
            iold=inew
            inew=isave
            nold=j
         endif
      enddo

      !--then add in the incident energy points from mf6/mt2
      !  save current value of nin and switch to nendf
      nins=nin
      nin=nendf
      call findf(matd,6,2,nin)
      call contio(nin,0,0,scr,nb,nw)
      call tab1io(nin,0,0,scr,nb,nw)
      ll=1+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,scr(ll),nb,nw)
         ll=ll+nw
      enddo
      call tab2io(nin,0,0,scr,nb,nw)
      nee=nint(scr(6))
      j=0
      k=1
      call finda(k,c,2,iold,bufn,nbuf)
      eg=c(1)
      do iee=1,nee
         call listio(nin,0,0,scr,nb,nw)
         ll=1+nb
         do while (nb.ne.0)
            call moreio(nin,0,0,scr(ll),nb,nw)
            ll=ll+nw
         enddo
         ee=scr(2)
         do while (eg.lt.ee)
            j=j+1
            call loada(j,c,2,inew,buf,nbuf)
            k=k+1
            call finda(k,c,2,iold,bufn,nbuf)
            eg=c(1)
         enddo
         j=j+1
         jt=j
         if (iee.eq.nee) jt=-jt
         c(1)=ee
         c(2)=0
         call loada(jt,c,2,inew,buf,nbuf)
         if (eg.eq.ee.and.jt.gt.0) then
            k=k+1
            call finda(k,c,2,iold,bufn,nbuf)
            eg=c(1)
         endif
      enddo
      isave=iold
      iold=inew
      inew=isave
      nold=j

      !--now add extra points for panels with large steps
      j=0
      k=1
      call finda(k,c,2,iold,bufn,nbuf)
      eg=c(1)
      do while (k.lt.nold)
         k=k+1
         call finda(k,c,2,iold,bufn,nbuf)
         if (k.lt.nold) then
            er=c(1)
            egrid=0
            do while (egrid.lt.er)
               j=j+1
               c(1)=eg
               c(2)=0
               call loada(j,c,2,inew,buf,nbuf)
               egrid=sigfig(step*eg,2,0)
               if (egrid.lt.er) eg=egrid
            enddo
            eg=er
         endif
      enddo
      j=j+1
      jt=-j
      call loada(jt,c,2,inew,buf,nbuf)
      isave=iold
      iold=inew
      inew=isave
      nold=j
      nc=2
      nen=nold
      elast=c(1)
      !--switch back to original nin tape
      nin=nins
      call findf(matd,3,0,nin)
   endif

   !--write other desired reactions on this grid.
   mfh=3
   do while (mfh.ne.0)
      call contio(nin,0,0,scr,nb,nw)
      if (mfh.ne.0) then

      !--eliminate redundant reactions
      if (mth.eq.3.or.(mth.eq.4.and.izai.eq.1)) then
         idone=0
         if (nf12s.ne.0) then
            i=0
            do while (i.lt.nf12s.and.idone.eq.0)
                  i=i+1
                  if (mth.eq.mf12s(i)) idone=1
               enddo
            endif
            if (idone.eq.0.and.nf16s.ne.0) then
               idone=0
               i=0
               do while (i.lt.nf16s.and.idone.eq.0)
                  i=i+1
                  if (mth.eq.mf16s(i)) idone=1
               enddo
            endif
            ! keep mt=4 if it is needed for ur competition
            if (mtcomp.eq.4.and.mth.eq.4) idone=1
         else
            idone=1
            if ((mt19.eq.1.and.mth.eq.18).or.&
              (mt19.eq.0.and.&
              (mth.eq.19.or.mth.eq.20.or.mth.eq.21.or.mth.eq.38)).or.&
              (mth.eq.26.or.mth.eq.27).or.&
              (mth.eq.101.or.mth.eq.120).or.&
              (mth.eq.151).or.&
              (mth.gt.207.and.mth.lt.221).or.&
              (mth.ge.221.and.mth.le.260).or.&
              (mth.gt.250.and.mth.le.300).or.&
              (mth.gt.301.and.mth.lt.444).or.&
              (iverf.lt.6.and.(mth.gt.444.and.mth.lt.700)).or.&
              (iverf.lt.6.and.mth.gt.800).or.&
              (iverf.ge.6.and.(mth.gt.444.and.mth.lt.600)).or.&
              (iverf.ge.6.and.(mth.ge.851.and.mth.le.870)).or.&
              (iverf.ge.6.and.mth.gt.900)) then
               idone=0
            endif
            if (mth.eq.10) then
               idone=0
               call mess('unionx','redundant mt=10 found',&
                 'cross section and distribution excluded')
            endif
         endif

         !--skip this section
         if (idone.eq.0) then
            call tosend(nin,0,0,scr)

         !--include this section
         else
            mtn=mth
            if (mth.eq.103) mt103=1
            if (mth.eq.104) mt104=1
            if (mth.eq.105) mt105=1
            if (mth.eq.106) mt106=1
            if (mth.eq.107) mt107=1
            call contio(0,0,nscr,scr,nb,nw)
            e=0
            call gety1(e,thresh,idis,y,nin,scr)
            thrx=(awr+awi)*(-c2h)/awr
            if (thrx.lt.zero) thrx=0
            test=0
            if (thrx.ne.zero) test=thresh/thrx
            write(messs,'(i6,1p,3e15.7)') mth,thresh,thrx,test
            if (test.lt.one.and.test.ne.zero)&
              call mess('unionx','threshold error',messs)
            if (mth.eq.2.or.mth.eq.3.or.mth.eq.102.or.mth.eq.301)&
               thresh=elow
            if (mth.eq.444) thresh=elow
            k=8
            kbase=k
            scr2(1)=scr(1)
            scr2(2)=scr(2)
            scr2(3)=scr(3)
            scr2(4)=scr(4)
            scr2(5)=1
            scr2(6)=0
            scr2(7)=0
            scr2(8)=2
            np=0
            do i=1,nen
               if (i.eq.1) ifrst=0
               call finda(i,c,nc,iold,buf,nbuf)
               if (c(1).ge.(1-eps)*thresh) then
                  e=c(1)
                  if (i.eq.1) e=e*(1+eps)
                  call gety1(e,enext,idis,y,nin,scr)
                  if (mth.eq.2.and.e.lt.ethrr) y=0
                  if (i.gt.1.and.ifrst.eq.0) y=0
                  ifrst=1
                  if (mth.le.45.or.mth.ge.50) then
                     if ((mth.le.120)&
                       .or.(mt103.eq.0.and.mth.ge.mpmin.and.mth.le.mpmax)&
                       .or.(mt104.eq.0.and.mth.ge.mdmin.and.mth.le.mdmax)&
                       .or.(mt105.eq.0.and.mth.ge.mtmin.and.mth.le.mtmax)&
                       .or.(mt106.eq.0.and.mth.ge.m3min.and.mth.le.m3max)&
                       .or.(mt107.eq.0.and.mth.ge.m4min.and.mth.le.m4max)) then
                        if (mth.eq.2) c(2)=0
                        c(2)=c(2)+y
                     endif
                  endif
               endif
               it=i
               if (i.eq.nen) it=-it
               call loada(it,c,2,inew,bufn,nbuf)
               if (c(1).ge.(1-eps)*thresh) then
                  if (np.le.0) then
                     np=nen-i+1
                     scr2(6)=np
                     scr2(7)=np
                  endif
                  k=k+2
                  scr2(k-1)=c(1)
                  scr2(k)=y
                  if (k.ge.(npage+kbase).or.i.eq.nen) then
                     if (kbase.ne.0) then
                        call tab1io(0,0,nscr,scr2,nb,k)
                        k=0
                        kbase=0
                     else
                        call moreio(0,0,nscr,scr2,nb,k)
                        k=0
                     endif
                  endif
               endif
            enddo
            call tosend(nin,0,0,scr)
            call asend(nscr,0)
            isave=iold
            nxc=nxc+1
            mfs(nxc)=mfh
            mts(nxc)=mtn
            ncs(nxc)=3+(np+2)/3
            iold=inew
            inew=isave
            nc=2
         endif
      endif
   enddo

   !--all reactions added.  write out file end card
   call contio(0,0,nscr,scr,nb,nw)

   !--write out new total cross section.
   mfh=3
   mth=1
   call contio(0,nout,0,b,nb,nw)
   call repoz(nscr)
   scr(1)=0
   scr(2)=0
   scr(3)=lt
   scr(4)=lr
   scr(5)=1
   scr(6)=nen
   scr(7)=nen
   scr(8)=2
   ibase=8
   j=ibase
   do i=1,nen
      call finda(i,c,2,iold,buf,nbuf)
      j=j+2
      scr(j-1)=c(1)
      scr(j)=sigfig(c(2),7,0)
      if (j.ge.(npage+ibase).or.i.eq.nen) then
         if (ibase.ne.0) then
            call tab1io(0,nout,0,scr,nb,j)
            j=0
            ibase=0
         else
            call moreio(0,nout,0,scr,nb,j)
            j=0
         endif
      endif
   enddo
   ncs(nxcs)=3+(nen+2)/3
   call asend(nout,0)

   !--copy unionized data for other reactions
   !--from scratch tape to output tape
   mfh=3
   do while (mfh.ne.0)
      call contio(nscr,nout,0,scr,nb,nw)
      if (mfh.ne.0) then
         if (mth.ne.0) then
            call tab1io(nscr,nout,0,scr,nb,nw)
            do while (nb.ne.0)
               call moreio(nscr,nout,0,scr,nb,nw)
            enddo
         endif
      endif
   enddo
   deallocate(scr)
   deallocate(buf)
   deallocate(bufn)
   deallocate(scr2)
   call closz(nscr)
   if (ngmt.eq.0) then  !will close these in gamsum when ngmt != 0
      call closz(iold)
      call closz(inew)
   endif
   return
   end subroutine unionx

   real(kr) function awt(iwtt,e)
   !--------------------------------------------------------------------
   ! Compute weight function for integral thinning.
   !--------------------------------------------------------------------
   ! externals
   integer::iwtt
   real(kr)::e
   ! internals
   real(kr)::test

   if (iwtt.le.1) then
      awt=1
   else
      test=1
      test=test/10
      if (e.lt.test) awt=10
      if (e.ge.test) awt=1/e
   endif
   return
   end function awt

   subroutine topfil(nin,nout,matd,newfor)
   !-------------------------------------------------------------------
   ! Prepare Files 4, 5, and 6 for further processing.  For the old
   ! format, angular distributions are converted to equally probable
   ! bin form.  Sections of File 6 using lab or cm LF=1 Legendre or
   ! tabulated data are converted to LF=7 angle-energy form.  Other
   ! sections are just copied.  For the new formats, just copy
   ! MF4, 5, and 6.  The conversion of distributions will be done
   ! in acelod.  There is an option triggered by no7=1 that will
   ! cause any law=7 sections to be converted to law=1 with tabulated
   ! angular distributions.  These sections will then end up using
   ! ACE LAW=61 in the final library.  It is set by default.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util ! provides repoz,error,openz,closz
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,matd,newfor
   ! internals
   integer::iso,ltt,lf,lvt,ltt3,lttn,nwb,no,ne,idone,ie
   integer::now,ne1,nk,new6,ik,mtd,l,l2,nmu,l3,imu,l1,ir,ip
   integer::nep,ll,iep,nb,nw,iend,nwmax,nin0,lang,i,np
   integer::intmu,loc,intep,igrd,ngrd,m,n,ians,ipp,irr,idis
   integer::nt1w,loct1,loct2,ne2,nn,nr,nf,locmx,npe
   integer::lis,jp,jpn,jpp,jpnut
   character(60)::string
   real(kr)::dzap,test,zap,e1,e2,f,ei,ep,epn,ss,ff,dmu
   real(kr)::b(50)
   real(kr)::y,enext
   real(kr),dimension(:),allocatable::tab1
   real(kr),dimension(:),allocatable::scr
   integer,parameter::nwmaxn=1000000
   real(kr),parameter::big=1.e9_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1

   !--this flag says convert law=7 to law=1
   integer::no7=1

   !--initialize and compute coefficients.
   jp=0
   jpn=0
   jpp=0
   jpnu=0
   nkk=0
   nt1w=20000
   allocate(tab1(nt1w))
   nsix=19
   nin0=0
   mcoars=nbina
   nsc=1
   npt=mcoars+1
   nwmax=nwmaxn
   call ptinit
   allocate(scr(nwmax))
   write(nsyso,'(/)')

   !--read through sections of endf tape.
   call repoz(nin)
   call tpidio(nin,0,0,scr,nb,nw)
   if (izai.ne.1) then
      call findf(matd,6,0,nin)
   else
      call findf(matd,4,0,nin)
   endif
   iend=0
   do while (iend.eq.0)
      iso=0
      call contio(nin,0,0,scr,nb,nw)
      ltt=0
      lf=0
      if (mfh.gt.6.or.math.eq.0) then
         iend=1

      !--just copy send or fend records
      else if (mth.eq.0.or.mfh.eq.0) then
         call contio(0,nout,0,scr,nb,nw)

      !--skip over some particular reaction sections
      else if ((mt19.eq.1.and.mth.eq.18).or.(mt19.eq.0.and.&
        (mth.eq.19.or.mth.eq.20.or.mth.eq.21.or.mth.eq.38)).or.&
        (mfh.eq.6.and.mth.eq.10).or.&
        (mfh.eq.5.and.mth.gt.900)) then
         call tosend(nin,0,0,scr)

      !--for the new format, just copy mf4 sections.
      else if (mfh.eq.4.and.newfor.eq.1) then
         call contio(0,nout,0,scr,nb,nw)
         nxc=nxc+1
         if (nxc.gt.nxcmax) call error('topfil','nxc.gt.nxcmax.',' ')
         mfs(nxc)=mfh
         mts(nxc)=mth
         ncs(nxc)=1
         call tosend(nin,nout,0,scr)

      !--for the old format, process file 4.
      !--skip transformation matrix in mf4.
      else if (mfh.eq.4) then
         lvt=l1h
         ltt=l2h
         scr(3)=0
         scr(4)=2
         nxc=nxc+1
         if (nxc.gt.nxcmax) call error('topfil','nxc.gt.nxcmax.',' ')
         mfs(nxc)=mfh
         mts(nxc)=mth
         ncs(nxc)=0
         call contio(0,nout,0,scr,nb,nw)
         ltt3=ltt
         lttn=0
         if (ltt3.eq.3) then
            ltt=1
            lttn=1
         endif
         ncs(nxc)=ncs(nxc)+1
         if (lvt.ne.0) then
            call listio(nin,0,0,scr,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,0,0,scr(7),nb,nw)
            enddo
         else
            call contio(nin,0,0,scr,nb,nw)
            iso=nint(scr(3))
         endif
         scr(5)=0
         scr(6)=0
         call contio(0,nout,0,scr,nb,nw)
         ncs(nxc)=ncs(nxc)+1

         !--convert distributions to equal probability form.
         if (iso.ne.1) then
            call tab2io(nin,0,0,b,nb,nwb)
            no=nout
            if (ltt3.eq.3) then
               no=nscr
               call repoz(nscr)
            else
               call tab2io(0,nout,0,b,nb,nwb)
            endif
            ne=nint(b(6))
            ncs(nxc)=ncs(nxc)+2
            idone=0
            do while (idone.eq.0)
               if (ltt3.eq.3.and.lttn.eq.2) ltt=2
               lttn=lttn+1
               do ie=1,ne
                  if (ltt.ne.2) then
                     call listio(nin,0,0,scr,nb,nw)
                     now=1+nw
                     do while (nb.ne.0)
                        call moreio(nin,0,0,scr(now),nb,nw)
                        now=now+nw
                     enddo
                     now=now-1
                     call ptleg(no,scr)
                  else
                     call tab1io(nin,0,0,scr,nb,nw)
                     now=1+nw
                     do while (nb.ne.0)
                        call moreio(nin,0,0,scr(now),nb,nw)
                        now=now+nw
                     enddo
                     call pttab(ltt,scr,no)
                  endif
               enddo
               idone=1
               if (ltt3.eq.3.and.lttn.eq.2) then
                  ne1=ne
                  call tab2io(nin,0,0,scr,nb,nw)
                  ne=nint(scr(6))
                  idone=0
               else if (ltt3.eq.3.and.lttn.eq.3) then
                  ne=ne+ne1
                  b(6)=ne
                  b(7)=ne
                  call tab2io(0,nout,0,b,nw,nwb)
                  call asend(nscr,0)
                  call repoz(nscr)
                  call tosend(nscr,nout,0,scr)
               endif
            enddo
         endif
         call tosend(nin,nout,0,scr)

      !--for file 5, get mf, mt, and tab1 lf values.
      !--if lf.ne.1, copy as is to nout
      !--if lf.eq.1, check secondary spectrum tab1 functions for
      !--multiple [e,f(e)=0.] data pairs.  For non-histogram
      !--interpolation, eliminate the lower energy pairs before
      !--writing the function to nout.  For histogram interpolation,
      !--eliminate all low-energy f(e)=0 data pairs.
      else if (mfh.eq.5) then
         call contio(0,nout,0,scr,nb,nw)
         nxc=nxc+1
         if (nxc.gt.nxcmax) call error('topfil','nxc.gt.nxcmax.',' ')
         mfs(nxc)=mfh
         mts(nxc)=mth
         ncs(nxc)=1
         nk=n1h
         ik=0
     111 ik=ik+1
         loct1=1
         ! read the initial tab1 and get lf.
         call tab1io(nin,0,0,tab1(loct1),nb,nw)
         do while (nb.ne.0)
            loct1=loct1+nw
            if (loct1.gt.nt1w) then
               call error('topfil','tab1 allocation is too small',' ')
            endif
            call moreio(nin,0,0,tab1(loct1),nb,nw)
         enddo
         lf=nint(tab1(4))
         ! move this tab1 to nout (all the time).
         loct1=1
         call tab1io(0,nout,0,tab1(loct1),nb,nw)
         loct1=loct1+nw
         do while (nb.ne.0)
            call moreio(0,nout,0,tab1(loct1),nb,nw)
            loct1=loct1+nw
         enddo
         ! if not lf=1, write the rest of this section to nout
         if (lf.ne.1) then
            call tosend(nin,nout,0,scr)
         else
            ! lf=1, read the tab2 record
            loct2=1
            call tab2io(nin,0,0,tab1(loct2),nb,nw)
            ne2=nint(tab1(6))
            ! copy this tab2 to nout (all the time).
            call tab2io(0,nout,0,tab1(loct2),nb,nw)
            ! check secondary tab1 records
            do nn=1,ne2
               loct1=1
               call tab1io(nin,0,0,tab1(loct1),nb,nw)
               do while (nb.ne.0)
                  loct1=loct1+nw
                  if (loct1.gt.nt1w) then
                     call error('topfil','tab1 allocation is too small2',' ')
                  endif
                  call moreio(nin,0,0,tab1(loct1),nb,nw)
               enddo
               nr=nint(tab1(5))
               nf=nint(tab1(6))
               ! check tab1 for multiple f(e)=0 data.
               ! if the first f(e) is non-zero, nothing else to do.
               ! if histogram interpolation, check from the first
               ! f(e) value, if not check from the second f(e) value
               loc=8+2*nr
               if (tab1(loc).eq.zero) then
                  if (nint(tab1(8)).ne.1) loc=loc+2
                  locmx=7+2*nr+2*nf
                  npe=0
                  do while (tab1(loc).eq.zero)
                     if (loc.gt.locmx)then
                        write(string,'(a,i3,a)')'mf=5,mt=',mth,&
                          ', entire tab1 function is zero.'
                        call error('topfil',string,' ')
                     endif
                     loc=loc+2
                     npe=npe+1
                  enddo
                  if (npe.ne.0) then
                     ! yes multiple zero data found, eliminate them.
                     ! fix nf;
                     ! fix the pointer array (interpolation code
                     !  array remains the same);
                     ! shift tab1 array.
                     tab1(6)=tab1(6)-float(npe)
                     loc=7
                     do n=loc,loc+2*nr,2
                        tab1(n)=tab1(n)-float(npe)
                     enddo
                     loc=6+2*nr
                     do n=1,nint(2*tab1(6))
                        tab1(loc+n)=tab1(loc+n+2*npe)
                     enddo
                  endif
               endif
               ! write the modified tab1 record to nout.
               loct1=1
               call tab1io(0,nout,0,tab1(loct1),nb,nw)
               loct1=loct1+nw
               do while (nb.ne.0)
                  call moreio(0,nout,0,tab1(loct1),nb,nw)
                  loct1=loct1+nw
               enddo
            enddo
            if (ik.lt.nk) go to 111
            ! copy section end record to nout.
            call contio(nin,nout,0,scr,nb,nw)
         endif

      !--check for conversion to law7
      else if (mfh.eq.6) then
         if (newfor.eq.0) then
            nk=n1h
            new6=0
            dzap=1
            test=1
            test=test/1000
            ik=0
            do while (dzap.gt.test.and.ik.lt.nk)
               ik=ik+1
               call tab1io(nin,0,0,scr,nb,nw)
               zap=c1h
               lf=l2h
               dzap=abs(zap-1)
               do while(nb.ne.0)
                  call moreio(nin,0,0,scr,nb,nw)
               enddo
               if (lf.eq.6) then
                  call contio(nin,0,0,scr,nb,nw)
               else if (lf.eq.1.or.lf.eq.2.or.lf.eq.5) then
                  call tab2io(nin,0,0,scr,nb,nw)
                  lang=l1h
                  if (dzap.le.test.and.lf.eq.1.and.lang.ne.2) new6=1
                  ne=n2h
                  do ie=1,ne
                     call listio(nin,0,0,scr,nb,nw)
                     do while (nb.ne.0)
                        call moreio(nin,0,0,scr,nb,nw)
                     enddo
                  enddo
               else if (lf.eq.7) then
                  call tab2io(nin,0,0,scr,nb,nw)
                  ne=n2h
                  do ie=1,ne
                     call tab2io(nin,0,0,scr,nb,nw)
                     nmu=n2h
                     do imu=1,nmu
                        call tab1io(nin,0,0,scr,nb,nw)
                        do while (nb.ne.0)
                           call moreio(nin,0,0,scr,nb,nw)
                        enddo
                     enddo
                  enddo
               endif
            enddo
            mtd=mth
            call findf(matd,6,mtd,nin)

            !--convert mf6,law1 lab or cm distributions to law7
            if (new6.ne.0) then
               if (nin.lt.0) nsix=-iabs(nsix)
               call openz(nsix,1)
               call fix6(nin,nsix,ik)
               call repoz(nsix)
               nin0=nin
               nin=nsix
            endif
            call contio(nin,0,0,scr,nb,nw)
         endif

         !--work on file 6
         ltt=0
         jp=l1h
         jpn=mod(jp,10)
         jpp=(jp-jpn)/10
         nk=n1h
         ik=0
         scr(5)=nk
         nxc=nxc+1
         if (nxc.gt.nxcmax) call error('topfil','nxc.gt.nxcmax.',' ')
         mfs(nxc)=mfh
         mts(nxc)=mth
         ncs(nxc)=0
         call contio(0,nout,0,scr,nb,nw)
         ncs(nxc)=ncs(nxc)+1
         do while (ik.lt.nk)
            ik=ik+1
            call tab1io(nin,0,0,scr,nb,nw)
            ncs(nxc)=ncs(nxc)+1+nint((scr(5)+2)/3)&
              +nint((scr(6)+2)/3)
            lis=l1h
            lf=l2h
            if (nint(scr(1)).eq.1.and.&
                (jpn.eq.2 .or. (jpn.eq.1.and.ik.ge.2))) then
               nkk=nkk+1
               jpnut=6+2*(n1h+n2h)
               if (jpnut.gt.jpnu) jpnu=jpnut
            endif
            if (newfor.eq.1.and.lf.eq.7.and.no7.eq.1) scr(4)=1
            call tab1io(0,nout,0,scr,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,nout,0,scr,nb,nw)
            enddo

            ! negative laws, law 0, 3 and 4 have no law dependent
            ! structure so no need to do anything
            if (lf.eq.6) then
               call contio(nin,nout,0,scr,nb,nw)
            else if (lf.eq.7.and.newfor.eq.1.and.no7.eq.1) then
               ! law=7 for newfor=1 -- convert the law7
               ! data into law1 format.
               call tab2io(nin,0,0,b,nb,nw)
               ne=nint(b(6))
               do ie=1,ne
                  ! read in the data
                  call tab2io(nin,0,0,scr,nb,nw)
                  ei=scr(2)
                  intmu=nint(scr(8))
                  nmu=n2h
                  loc=1+nmu
                  do imu=1,nmu
                     scr(imu)=loc
                     call tab1io(nin,0,0,scr(loc),nb,nw)
                     intep=nint(scr(loc+7))
                     loc=loc+nw
                     do while (nb.ne.0)
                        call moreio(nin,0,0,scr(loc),nb,nw)
                        loc=loc+nw
                     enddo
                  enddo
                  ! fix up the tab2 for law1
                  if (ie.eq.1) then
                     b(3)=10+intmu
                     b(4)=intep
                     call tab2io(0,nout,0,b,nb,nw)
                     ncs(nxc)=ncs(nxc)+2
                  endif
                  ! construct a union grid for eprime
                  igrd=loc
                  ngrd=0
                  do imu=1,nmu
                     loc=nint(scr(imu))
                     m=nint(scr(loc+4))
                     n=nint(scr(loc+5))
                     do iep=1,n
                        ngrd=ngrd+1
                        scr(igrd+ngrd-1)=scr(loc+4+2*m+2*iep)
                     enddo
                  enddo
                  call ordr(scr(igrd),ngrd)
                  ! interpolate for angular distributions
                  ! on the union eprime grid to construct
                  ! the law1 distribution.
                  ians=igrd+ngrd
                  scr(ians)=0
                  scr(ians+1)=ei
                  scr(ians+2)=0
                  scr(ians+3)=2*nmu
                  scr(ians+4)=ngrd*(2+2*nmu)
                  scr(ians+5)=ngrd
                  ll=ians+6
                  do iep=1,ngrd
                     ep=scr(igrd+iep-1)
                     scr(ll)=ep
                     ss=0
                     do imu=1,nmu
                        loc=nint(scr(imu))
                        ipp=2
                        irr=1
                        call terpa(ff,ep,epn,idis,scr(loc),ipp,irr)
                        scr(ll+2*imu)=scr(loc+1)
                        scr(ll+1+2*imu)=ff
                        if (imu.gt.1) then
                           dmu=scr(ll+2*imu)-scr(ll+2*imu-2)
                           if (intmu.eq.1) then
                              ss=ss+dmu*scr(ll+1+2*imu-2)
                           else
                              ss=ss+dmu*&
                                (scr(ll+1+2*imu)+scr(ll+1+2*imu-2))/2
                           endif
                        endif
                     enddo
                     scr(ll+1)=ss
                     do imu=1,nmu
                        if (ss.ne.zero) then
                           scr(ll+1+2*imu)=scr(ll+1+2*imu)/ss
                        else
                           scr(ll+1+2*imu)=one/2
                        endif
                     enddo
                     ll=ll+2+2*nmu
                  enddo
                  call listio(0,nout,0,scr(ians),nb,nw)
                  ll=ians+nw
                  do while (nb.ne.0)
                     call moreio(0,nout,0,scr(ll),nb,nw)
                     ll=ll+nw
                  enddo
                  nw=ngrd*(2+2*nmu)
                  nw=(nw+5)/6
                  ncs(nxc)=ncs(nxc)+1+nw
               enddo
            else if (lf.eq.1.or.lf.eq.2.or.lf.eq.5.or.lf.eq.7) then
               call tab2io(nin,nout,0,b,nb,nw)
               ne=nint(b(6))
               ncs(nxc)=ncs(nxc)+2
               do ie=1,ne
                  ! law=1 -- copy the subsection,
                  ! but also check interpolation
                  if (lf.eq.1) then
                     call listio(nin,0,0,scr,nb,nw)
                     now=1+nw
                     do while (nb.ne.0)
                        call moreio(nin,0,0,scr(now),nb,nw)
                        now=now+nw
                     enddo
                     call cptab(nout,scr)
                  ! law=2 for newfor=1 - copy the subsection
                  else if (lf.eq.2.and.newfor.eq.1) then
                     call listio(nin,nout,0,scr,nb,nw)
                     now=1+nw
                     do while (nb.ne.0)
                        call moreio(nin,nout,0,scr(now),nb,nw)
                        now=now+nw
                     enddo
                  ! law=2 for newfor=0 - convert to probability bins
                  else if (lf.eq.2.and.newfor.eq.0) then
                     call listio(nin,0,0,scr,nb,nw)
                     now=1+nw
                     do while (nb.ne.0)
                        call moreio(nin,0,0,scr(now),nb,nw)
                        now=now+nw
                     enddo
                     now=now-1
                     lang=nint(scr(3))
                     if (lang.eq.0) then
                        ! legendre coefficients
                        call ptleg(nout,scr)
                     else
                        ! tabulated angular distribution
                        do i=1,now
                           scr(now+2-i+1)=scr(now-i+1)
                        enddo
                        np=nint(scr(8))
                        scr(1)=scr(3)
                        scr(2)=scr(4)
                        scr(3)=0
                        scr(4)=0
                        scr(5)=1
                        scr(6)=np
                        scr(7)=np
                        scr(8)=lang-10
                        call pttab(ltt,scr,nout)
                     endif
                  ! law=5 -- copy the subsection
                  else if (lf.eq.5) then
                     call listio(nin,nout,0,scr,nb,nw)
                     now=1+nw
                     do while (nb.ne.0)
                        call moreio(nin,nout,0,scr(now),nb,nw)
                        now=now+nw
                     enddo
                  ! law=7 -- the tab2 is converted to a special tab1
                  ! containing the overall angular distribution and
                  ! the angle-energy data are copied
                  else if (lf.eq.7) then
                     l=1
                     call tab2io(nin,0,0,scr(l),nb,nw)
                     l=l+nw
                     l2=l
                     nmu=n2h
                     l=l+66
                     l3=l
                     do imu=1,nmu
                        l1=l
                        call tab1io(nin,0,0,scr(l),nb,nw)
                        l=l+nw
                        do while (nb.ne.0)
                           call moreio(nin,0,0,scr(l),nb,nw)
                           l=l+nw
                        enddo
                        e1=0
                        e2=big
                        ir=1
                        ip=2
                        call intega(f,e1,e2,scr(l1),ip,ir)
                        scr(l2+2*imu-2)=scr(l1+1)
                        scr(l2+2*imu-1)=f
                        nr=nint(scr(l1+4))
                        nep=nint(scr(l1+5))
                        ll=l1+5+2*nr
                        do iep=1,nep
                           if (f.ne.0) scr(ll+2*iep)=scr(ll+2*iep)/f
                        enddo
                     enddo
                     scr(4)=nmu
                     if (newfor.eq.0) then
                        call pttab(3,scr,nout)
                     else
                        call tab1io(0,nout,0,scr,nb,nw)
                     endif
                     l=l3
                     do imu=1,nmu
                        call tab1io(0,nout,0,scr(l),nb,nw)
                        l=l+nw
                        do while (nb.ne.0)
                           call moreio(0,nout,0,scr(l),nb,nw)
                           l=l+nw
                        enddo
                     enddo
                  endif
               enddo
            endif
         enddo
         call tosend(nin,nout,0,scr)
         if (nin.eq.nsix) then
            call closz(nin)
            nin=nin0
         endif
      endif
   enddo

   !--topfil is finished.
   deallocate(scr)
   deallocate(tab1)
   return
   end subroutine topfil

   subroutine ordr(x,n)
   !-------------------------------------------------------------------
   ! sort the n elements of x into ascending order
   ! removing any duplicate elements
   !-------------------------------------------------------------------
   ! externals
   real(kr)::x(*)
   integer::n
   ! internals
   integer::i,j,k,m
   real(kr)::tsave
   real(kr),parameter::small=1.e-10_kr

   if (n.le.2) return
   ! sort
   i=0
  110 i=i+1
   j=i
  120 j=j+1
   if (x(j).lt.x(i)) then
      tsave=x(j)
      x(j)=x(i)
      x(i)=tsave
   endif
   if (j.lt.n) go to 120
   if (i.lt.n-1) go to 110
   ! remove duplicates
   m=n
   i=1
   do while (i.lt.m)
      i=i+1
      if (abs(x(i)-x(i-1)).le.small*x(i)) then
         m=m-1
         do k=i,m
            x(k)=x(k+1)
         enddo
         i=i-1
      endif
   enddo
   n=m
   return
   end subroutine ordr

   subroutine ptinit
   !-------------------------------------------------------------------
   ! Initialize the calculation of equal probability bins
   ! from Legendre coefficients.  See ptleg.
   !-------------------------------------------------------------------
   ! internals
   integer::nwords,n,no,ist,ifi,jfi,ivar,kxat,k
   real(kr)::afac,asign,cf,aa,b,c,d
   integer,parameter::ni=64

   !--assign storage
   nwords=ni*(1+ni/2)
   allocate(xat(nwords))

   !--calculate xaterm(n,k) - the coefficients in
   !--the series for the integral of p-sub-n up to order ni.
   !--k is the number of terms in the series.
   afac=1
   asign=-1
   do n=1,ni
      if (n.gt.1) afac=afac*(2*n-1)/n
      cf=1
      no=1+n/2
      ist=0
      ifi=1
      jfi=1
      ivar=2
      aa=2*n+1
      aa=aa/2
      xat(n)=cf*afac*aa/(n+jfi)
      if (no.ge.2) then
         kxat=n-ni
         do k=2,no
            jfi=jfi-2
            b=(n-ist)*(n-ifi)
            c=ivar*(2*n-ifi)
            d=n+jfi
            cf=cf*asign*b/c
            xat(ni*k+kxat)=cf*afac*aa/d
            ist=ist+2
            ifi=ifi+2
            ivar=ivar+2
         enddo
      endif
   enddo
   return
   end subroutine ptinit

   subroutine ptleg(nout,scr)
   !-------------------------------------------------------------------
   ! This subroutine translates ENDF Legendre ang dist data into
   ! tabulated form with equal probability mu intervals.
   ! Borrowed from ETOPL.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use endf ! provides endf routines and variables
   use util ! provides sigfig
   ! externals
   integer::nout
   real(kr)::scr(*)
   ! internals
   integer::mmfne,mfine,nord,j,iso,nz,nb,nw,int,ibase,mm,nexp,iloop
   integer::l,n,no,index,k,idone,mfmid,mvar,nwds,jscr
   real(kr)::e,delmuf,varmu,afunc,temp,sum,amumid,delmu1,delmu2
   real(kr)::avar,test,dprob,ap,b,c,d,g,dum
   real(kr)::afnc(1001),acum(1001),amu(1001)
   real(kr)::prob(1000),pro(100)
   real(kr)::acmu(64),fl(64),xxmu(65),pmu(65),expm(65)
   integer,parameter::ni=64
   real(kr),parameter::half=0.5e0_kr
   real(kr),parameter::oh2=0.02e0_kr
   real(kr),parameter::tiny=1.e-50_kr
   real(kr),parameter::zero=0

   ! set mmfne and mfine - to be defined later.
   mmfne=100
   mfine=1000

   !--work with list record read in topfil
   nord=n1h
   do j=1,nord
      fl(j)=scr(6+j)
   enddo

   !--test for isotropy.
   !--write isotropic section and return.
   iso=1
   if (fl(1).ne.zero) iso=0
   do nz=1,nord
      if (fl(nz).ne.zero) iso=0
   enddo
   if (iso.eq.1) then
      scr(5)=1
      scr(6)=2
      scr(7)=2
      scr(8)=2
      scr(9)=-1
      scr(10)=half
      scr(11)=1
      scr(12)=half
      nw=12
      call tab1io(0,nout,0,scr,nb,nw)
      ncs(nxc)=ncs(nxc)+3
      return
   endif

   !--process anisotropic distribution.
   int=1
   e=scr(2)
   scr(5)=1
   scr(6)=npt
   scr(7)=npt
   scr(8)=1
   ibase=8

   !--to compute the equal prob distributions, we first divide the
   !--mu interval from -1 to +1 into 100 (mmfne) equal subintervals.
   !--we then integrate and find the mu value at which the
   !--cumulative area first reaches a value above 0.5.  call ]
   !--this amumid.  then we divide the interval from -1 to amumid
   !--into 500 equal subintervals and do the same from amumid
   !--to +1.  so we have a total of 1000 (mfine) subintervals.
   delmuf=oh2
   mm=mmfne+1
   amu(1)=-1
   do j=2,mm
      amu(j)=amu(j-1)+delmuf
   enddo
   nexp=nord+1
   iloop=1
   do while (iloop.le.2)
      do j=1,mm
         varmu=amu(j)
         expm(1)=varmu
         do l=2,nexp
            expm(l)=expm(l-1)*varmu
            if (abs(expm(l)).lt.tiny) expm(l)=0
         enddo
         afunc=varmu/2
         do n=1,nord
            no=1+n/2
            index=n+1
            temp=fl(n)
            do k=1,no
               afunc=afunc+xat((k-1)*ni+n)*expm(index)*temp
               index=index-2
            enddo
         enddo
         afnc(j)=afunc
      enddo
      sum=0
      idone=0
      j=0
      do while (j.lt.mmfne.and.idone.eq.0)
         j=j+1
         prob(j)=afnc(j+1)-afnc(j)
         if (prob(j).le.zero) write(nsyso,'(&
           &'' ---message from ptleg---negative area between mu='',&
           &1p,e12.5,'' and '',e12.5,'',  '',e12.5,'',  e='',e12.5/&
           &''    mat='',i4,'', mf='',i2,'', mt='',i3)')&
           amu(j),amu(j+1),prob(j),e,math,mfh,mth
         sum=sum+prob(j)
         if (iloop.le.1.and.sum.gt.half) idone=1
         acum(j)=sum
      enddo
      if (iloop.eq.2) then
         iloop=3
      else
         amumid=amu(j)
         mfmid=mfine/2
         delmu1=(amumid+1)/mfmid
         amu(1)=-1
         mvar=mfmid+1
         do j=2,mvar
            amu(j)=amu(j-1)+delmu1
         enddo
         delmu2=(1-amumid)/mfmid
         mvar=mvar+1
         mm=mfine+1
         do j=mvar,mm
            amu(j)=amu(j-1)+delmu2
         enddo
         iloop=2
         ! note that after finding amumid we set mmfne=mfine -
         ! i.e. change from 100 to 1000.
         mmfne=mfine
      endif
   enddo

   !--sum is the integrated area computed using the 1000
   !--subintervals described above.  if it differs from 1.0
   !--by more than 1.0e-5, we print out a non-fatal diagnostic.
   avar=sum
   test=1
   test=test/100000
   if (abs(1-sum).gt.test) write(nsyso,'(&
     &'' ---message from ptleg---integrated area of legendre '',&
     &''dist'',/,20x,'' using 1000 subintervals is '',1p,e12.5/&
     &21x,''mat='',i4,'' mf='',i2,'' mt='',i3,'' e='',e12.6)')&
     avar,math,mfh,mth,e

   !--compute equal probability figure as sum (see comment above)
   !--divided by the number of coarse bins specified.
   dprob=sum/mcoars

   !--compute the mu values for the equal probability coarse bins.
   j=1
   k=1
   varmu=-1
   pro(k)=0
   idone=0
   do while (idone.eq.0)
      pro(k)=acum(j)
      varmu=amu(j+1)
      if (pro(k).eq.sum) then
         idone=1
      else
         ap=pro(k)-dprob*k
         if (ap.eq.zero) then
            acmu(k)=varmu
            j=j+1
            k=k+1
         else if (ap.gt.zero) then
            b=dprob*k
            c=amu(j+1)-amu(j)
            d=b-acum(j-1)
            g=acum(j)-acum(j-1)
            varmu=amu(j)+(d*c/g)
            acmu(k)=varmu
            k=k+1
         else
            j=j+1
         endif
      endif
   enddo
   acmu(k)=varmu
   xxmu(1)=-1
   do l=1,mcoars
      xxmu(l+1)=acmu(l)
   enddo

   !--knowing the mu values for the equal probability bins, we
   !--next compute the corresponding p(mu) values appropriate
   !--to a histogram representations of the data and write out
   !--the table.
   do j=1,mcoars
      pmu(j)=dprob/(xxmu(j+1)-xxmu(j))
   enddo
   pmu(npt)=0
   do l=1,npt
      scr(2*l-1+ibase)=sigfig(xxmu(l),7,0)
      scr(2*l+ibase)=sigfig(pmu(l),7,0)
   enddo
   nw=ibase+2*npt-1
   nwds=nw
   if (j.gt.npage) nw=ibase+npage
   call tab1io(0,nout,0,scr,nb,nw)
   do while (nb.ne.0)
      nwds=nwds-nw
      jscr=1+nw
      nw=nwds
      if (nw.gt.npage) nw=npage
      call moreio(0,nout,0,scr(jscr),nb,nw)
   enddo
   ncs(nxc)=ncs(nxc)+2+(npt+2)/3

   !--now compute area for mu-pmu histogram.  if area differs
   !--from 1.0 by more than 0.0001 print non-fatal diagnostic
   call summer(xxmu,pmu,npt,int,2,dum,e)

   !--finished
   return
   end subroutine ptleg

   subroutine cptab(nout,a)
   !-------------------------------------------------------------------
   ! Copies tabulated energy distributions as is, except sections
   ! with nr>1 and int>2 are linearized to nr=1.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   ! externals
   integer::nout
   real(kr)::a(*)
   ! internals
   integer::nr,np,int,nw,now,ll,k,i,ir,jnt,j,jn,nb
   real(kr)::dy,xm,ym,xn,xl,yt
   integer,parameter::imax=10
   real(kr)::x(imax),y(imax)
   real(kr),parameter::err=0.01e0_kr

   !--work with tab1 or list record read in topfil
   if (mfh.eq.6) then
      nr=1
      int=2
      nw=6+nint(a(5))
   else
      nr=nint(a(5))
      np=nint(a(6))
      int=nint(a(8))
      nw=6+2*nr+2*np
   endif
   now=nw+1
   ll=1
   if (nr.ne.1.or.int.gt.2) then

      !--linearize the secondary energy distribution
      ll=now
      a(ll)=a(1)
      a(ll+1)=a(2)
      a(ll+2)=a(3)
      a(ll+3)=a(4)
      k=0
      ! prime the adaptive stack
      i=2
      x(2)=a(6+2*nr+1)
      y(2)=a(6+2*nr+2)
      x(1)=a(6+2*nr+2*np-1)
      y(1)=a(6+2*nr+2*np)
      ! do adaptive reconstruction
      do while (i.gt.0)
         dy=0
         if (i.gt.1.and.i.lt.imax) then
            xm=(x(i-1)+x(i))/2
            ym=(y(i-1)+y(i))/2
            ir=1
            jnt=nint(a(6+2*(ir-1)+1))
            int=nint(a(6+2*(ir-1)+2))
            xn=xm
            j=0
            do while (j.lt.np-1.and.xm.ge.xn)
               j=j+1
               jn=j
               if (j.ge.jnt) then
                  ir=ir+1
                  jnt=nint(a(6+2*(ir-1)+1))
                  int=nint(a(6+2*(ir-1)+2))
               endif
               xn=a(6+2*nr+2*j+1)
            enddo
            xl=a(6+2*nr+2*(jn-1)+1)
            call terp1(xl,a(6+2*nr+2*(jn-1)+2),xn,a(6+2*nr+2*jn+2),&
              xm,yt,int)
            dy=abs(yt-ym)
         endif
         if (dy.gt.err*yt) then
            ! not converged.
            ! add the midpoint to the stack and continue.
            i=i+1
            x(i)=x(i-1)
            y(i)=y(i-1)
            x(i-1)=xm
            y(i-1)=yt
         else
            ! converged.
            ! use the top point off the stack.
            k=k+1
            a(ll+5+2+2*(k-1)+1)=x(i)
            a(ll+5+2+2*(k-1)+2)=y(i)
            i=i-1
         endif
      enddo
      nw=6+2+2*k
      a(ll+4)=1
      a(ll+5)=k
      a(ll+6)=k
      a(ll+7)=2
   endif

   !--copy subsection to output unit.
   if (mfh.eq.5) call tab1io(0,nout,0,a(ll),nb,nw)
   if (mfh.eq.6) call listio(0,nout,0,a(ll),nb,nw)
   now=nw+ll
   do while (nb.ne.0)
      call moreio(0,nout,0,a(now),nb,nw)
      now=now+nw
   enddo

   return
   end subroutine cptab

   subroutine pttab(ltt,a,nout)
   !-------------------------------------------------------------------
   ! Translates ENDF tabulated ang dist data into equal
   ! probability intervals.  Borrowed from ETOPL.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   integer::ltt,nout
   real(kr)::a(*)
   ! internals
   integer::lt,nr,np,nmu,intmu,int,ibase,i,nbt,nb,nw,l,j,k
   integer::idone
   real(kr)::t,e,dum,sum,dprob,area1,aneed,area,test,slpe
   real(kr)::aa,b,c,w,d,x1,x2
   integer,parameter::npmax=300
   real(kr)::amu(npmax),p(npmax),tbmu(npmax),pmu(npmax)
   ! set area of back angles to be omitted from distribution
   real(kr),parameter::aback=1.e-4_kr
   real(kr),parameter::zero=0

   !--work on tab1 record read (or constructed) by topfil
   t=a(1)
   e=a(2)
   lt=nint(a(3))
   nr=nint(a(5))
   np=nint(a(6))
   if (ltt.eq.3) nmu=np
   if (ltt.eq.3) intmu=nint(a(8))
   if (np.gt.npmax) call error('pttab','storage exceeded.',' ')
   ! if number of interpolation ranges gt. 1, print diagnostic and
   ! exit.
   if (nr.gt.2) call error('pttab',&
     'tab ang dist has more than one terp range.',' ')
   if (nr.gt.1) call chekit(a,nw)
   nr=nint(a(5))
   np=nint(a(6))
   nbt=nint(a(7))
   int=nint(a(8))
   ibase=8
   if (int.eq.3.or.int.eq.5)&
     call error('pttab','tab ang dist not allowed for.',' ')
   do i=1,np
      amu(i)=a(2*i-1+ibase)
      p(i)=a(2*i+ibase)
   enddo

   !--integrate distribution as given in endf files.
   !--if area differs from 1.0 by more than 0.0001,
   !--print non-fatal diagnostic.
   call summer(amu,p,np,int,1,sum,e)

   !--test for isotropy.
   !--write out mu-pmu pairs for an isotropic distribution
   !--and return.
   if (np.eq.2.and.p(1).eq.p(2)) then
      a(1)=t
      a(2)=e
      a(3)=lt
      a(4)=0
      if (ltt.eq.3) a(3)=intmu
      if (ltt.eq.3) a(4)=nmu
      a(5)=nr
      a(6)=np
      a(7)=nbt
      a(8)=int
      ibase=8
      do i=1,np
         a(2*i-1+ibase)=amu(i)
         a(2*i+ibase)=p(i)
      enddo
      call tab1io(0,nout,0,a,nb,nw)
      do while (nb.ne.0)
         call moreio(0,nout,0,a,nb,nw)
      enddo
      ncs(nxc)=ncs(nxc)+2+(np+2)/3
      return
   endif

   !--compute equal prob figure - delprob - as integrated
   !--area divided by number of coarse bins specified.
   dprob=sum/mcoars

   !--process according to interpolation scheme.
   !--int=1 -- histogram interpolation
   if (int.eq.1) then
      l=2
      j=2
      k=1
      area1=0
      aneed=dprob
      idone=0
      tbmu(1)=amu(1)
      do while (idone.eq.0)
         area=p(k)*(amu(j)-amu(k))
         if (area.gt.aneed.and.l.lt.npt) then
            if (p(k).le.zero) then
               p(k)=p(k-1)/1000
            endif
            tbmu(l)=(aneed/p(k))+amu(k)
            amu(k)=tbmu(l)
            area1=0
            aneed=dprob
            l=l+1
         else
            area1=area1+area
            aneed=aneed-area
            if (l.eq.2.and.area1.lt.aback) tbmu(1)=amu(j)
            if (j.ge.np) then
               idone=1
            else
               k=k+1
               j=j+1
            endif
         endif
      enddo
      tbmu(l)=amu(np)
      area1=area1-dprob
      test=1
      test=test/100000
      if (abs(area1).gt.test) write(nsyso,'(/&
         &'' ---message from pttab---'',&
         &''integrated area of distribution''/&
         &20x,''is off by '',1p,e12.5,'' at e='',e12.4/&
         &20x,''mat='',i4,'' mf='',i2,'' mt='',i3)')&
         area1,e,math,mfh,mth

   !--int=2 -- linear-linear interpolation
   else if (int.eq.2) then
      l=2
      j=2
      k=1
      area1=0
      aneed=dprob
      idone=0
      tbmu(1)=amu(1)
      do while (idone.eq.0)
         area=(p(j)+p(k))*(amu(j)-amu(k))/2
         if (area.gt.aneed.and.l.lt.npt) then
            test=1
            test=test/10000
            if (aneed.lt.test) then
               tbmu(l)=amu(k)+aneed/p(k)
            else
               slpe=(p(j)-p(k))/(amu(j)-amu(k))
               if (abs(slpe).lt.test) then
                  if (p(k).le.zero) p(k)=p(k-1)/1000
                  tbmu(l)=aneed/p(k)+amu(k)
               else
                  aa=slpe/2
                  b=p(k)-slpe*amu(k)
                  c=(slpe*amu(k)/2-p(k))*amu(k)-aneed
                  w=b*b-4*aa*c
                  if (w.le.zero) then
                     write(nsyso,'(/&
                       &'' ---message from pttab---'',&
                       &''neg arg in sqrt for int=2''/&
                       &''  mat='',i4,''  mf='',i2,''  mt='',i3,&
                       &''  e='',1p,e12.5)') math,mfh,mth,e
                     write(nsyso,'(/&
                       &''    values are as follows '',1p,4e12.5)')&
                        p(j),p(k),amu(j),amu(k)
                  endif
                  d=sqrt(w)
                  x1=(-b+d)/(2*aa)
                  x2=(-b-d)/(2*aa)
                  if (amu(k).lt.x1.and.x1.le.amu(j)) tbmu(l)=x1
                  if (amu(k).lt.x2.and.x2.le.amu(j)) tbmu(l)=x2
               endif
            endif
            p(k)=p(k)+slpe*(tbmu(l)-amu(k))
            amu(k)=tbmu(l)
            area1=0
            aneed=dprob
            l=l+1
         else
            area1=area1+area
            aneed=aneed-area
            if (l.eq.2.and.area1.lt.aback) tbmu(1)=amu(j)
            if (j.ge.np) then
               idone=1
            else
               k=k+1
               j=j+1
            endif
         endif
      enddo
      tbmu(l)=amu(np)
      area1=area1-dprob
      test=1
      test=test/100000
      if (abs(area1).gt.test) write(nsyso,'(/&
        &'' ---message from pttab---'',&
        &''integrated area of distribution''/&
        &20x,''is off by '',1p,e12.5,'' at e='',e12.4/&
        &20x,''mat='',i4,'' mf='',i2,'' mt='',i3)')&
        area1,e,math,mfh,mth

   !--int=4 -- lin-log interpolation
   else if (int.eq.4) then
      l=2
      j=2
      k=1
      area1=0
      aneed=dprob
      idone=0
      tbmu(1)=amu(1)
      do while (idone.eq.0)
         b=log(p(j)/p(k))/(amu(j)-amu(k))
         test=1
         test=test/10000
         if (abs(b).lt.test) then
            area=p(k)*(amu(j)-amu(k))
         else
            area=p(k)*(exp(b*(amu(j)-amu(k)))-1)/b
         endif
         if (area.gt.aneed.and.l.lt.npt) then
            test=1
            test=test/10000
            if (abs(b).lt.test) then
               tbmu(l)=aneed/p(k)+amu(k)
            else
               tbmu(l)=log(1+b*aneed/p(k))/b+amu(k)
            endif
            p(k)=p(k)*exp(b*(tbmu(l)-amu(k)))
            amu(k)=tbmu(l)
            area1=0
            aneed=dprob
            l=l+1
         else
            area1=area1+area
            aneed=aneed-area
            if (l.eq.2.and.area1.lt.aback) tbmu(1)=amu(j)
            if (j.ge.np) then
               idone=1
            else
               k=k+1
               j=j+1
            endif
         endif
      enddo
      tbmu(l)=amu(np)
      area1=area1-dprob
      test=1
      test=test/100000
      if (abs(area1).gt.test) write(nsyso,'(/&
        &'' ---message from pttab---'',&
        &''integrated area of distribution''/&
        &20x,''is off by '',1p,e12.5,'' at e='',e12.4/&
        &20x,''mat='',i4,'' mf='',i2,'' mt='',i3)')&
        area1,e,math,mfh,mth
   endif

   !--knowing the mu values for the equal probability bins,
   !--we next compute the corresponding p(mu) values appropriate
   !--to a histogram representation of the data for input
   !--interpolation schemes = 1, 2, and 4.
   do j=1,mcoars
      pmu(j)=dprob/(tbmu(j+1)-tbmu(j))
   enddo
   test=1
   if (tbmu(npt).ne.test.and.mfh.eq.4) write(nsyso,'(/&
     &'' ---message from pttab---last value of mu ne 1 at e='',&
     &1p,e12.4,'' for mt='',i3)') e,mth
   pmu(npt)=pmu(npt-1)

   !--write out mu-pmu pairs for equal probability bins.
   !--used for interpolation schemes 1, 2, and 4.
   a(1)=t
   a(2)=e
   a(3)=lt
   a(4)=0
   if (ltt.eq.3) a(3)=intmu
   if (ltt.eq.3) a(4)=nmu
   a(5)=nr
   a(6)=npt
   a(7)=npt
   a(8)=1
   ibase=8
   do i=1,npt
      a(2*i-1+ibase)=sigfig(tbmu(i),7,0)
      a(2*i+ibase)=sigfig(pmu(i),7,0)
   enddo
   nw=2*npt+ibase
   call tab1io(0,nout,0,a,nb,nw)
   do while (nb.ne.0)
      call moreio(0,nout,0,a,nb,nw)
   enddo
   ncs(nxc)=ncs(nxc)+2+(npt+2)/3

   !--integrate distribution as computed for equal prob bins.
   !--if area differs from 1.0 by more than 0.0001,
   !--print non-fatal diagnostic.
   call summer(tbmu,pmu,npt,1,2,dum,e)

   !--finished
   return
   end subroutine pttab

   subroutine chekit(a,nw)
   !-------------------------------------------------------------------
   ! Convert certain types of File 5 MT's with nr=2 to nr=1.
   !-------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::nw
   real(kr)::a(*)
   ! internals
   integer::nbt1,int1,nbt2,nwm,i
   real(kr),parameter::zero=0

   nbt1=nint(a(7))
   int1=nint(a(8))
   nbt2=nint(a(9))
   if (nbt2.ne.nbt1+1)&
     call error('chekit','wrong type of nr=2 file 5 mt.',' ')
   if (int1.ne.1)&
     call error('chekit','wrong type of nr=2 file 5 mt.',' ')
   if (a(nw).ne.zero.or.a(nw-1).ne.a(nw-3).or.a(nw-2).ne.a(nw-4))&
     call error('chekit','wrong type of nr=2 file 5 mt.',' ')
   a(6)=a(6)-1
   a(5)=1
   a(7)=a(6)
   nw=nw-4
   nwm=nw-1
   do i=9,nwm
      a(i)=a(i+2)
   enddo
   a(nw)=0
   return
   end subroutine chekit

   subroutine summer(x,y,nopts,jnt,itype,area,e)
   !-------------------------------------------------------------------
   ! Integrates lin-lin or log-lin tabular ang dist.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use endf ! provides math,mfh,mth,gral
   ! externals
   integer::nopts,jnt,itype
   real(kr)::x(*),y(*),area,e
   ! internals
   integer::j,jj
   real(kr),parameter::test=1.e-4_kr

   area=0
   do j=2,nopts
      jj=j-1
      area=area+gral(x(jj),y(jj),x(j),y(j),x(jj),x(j),jnt)
   enddo
   if (itype.ne.2) then
      if (abs(1-area).gt.test) write(nsyso,'(/&
         &'' ---message from summer---for distribution as per '',&
         &''endf'',/&
         &6x,''area='',1p,e12.5,''  mat='',i4,''  mf='',i2,&
         &''  mt='',i3,''  e='',e12.5)')&
         area,math,mfh,mth,e
   else
      if (abs(1-area).gt.test) write(nsyso,'(/&
         &'' ---message from summer---for distr with equal prob '',&
         &''bins, '',/,6x,''area='',1p,e12.5,''  mat='',i4,&
         &''  mf='',i2,''  mt='',i3,''  e='',e12.5)')&
         area,math,mfh,mth,e
   endif
   return
   end subroutine summer

   subroutine fix6(nin,nout,nk)
   !-------------------------------------------------------------------
   ! Convert a section of File 6 using Legendre or tabulated
   ! angular distributions into Law 7 format with 33 angles.
   ! Convert just the first nk subsections of this MF6 section.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use endf ! provides endf routines and variables
   use mathm ! provides legndr
   use util ! provides error
   ! externals
   integer::nin,nout,nk
   ! internals
   integer::ndebug,ncos,i,nb,nw,lct,ik,lang,lep,ne,l,ll
   integer::nep,ncyc,nl,nmu,imu,j,idone,l1,l2,iep,kmu
   integer::nnep,j2,mmu
   real(kr)::dmu,awr,ein,acos,ep,csn,elb,drv,clb,cmn,qq,aw1
   real(kr)::fmu
   integer,parameter::namax=9000
   real(kr)::a(namax)
   real(kr)::amu(50)
   real(kr)::p(65)
   real(kr),parameter::zero=0

   !--start the conversion process
   ndebug=nsyso
   ndebug=0
   ncos=33
   dmu=2
   dmu=dmu/(ncos-1)
   do i=1,ncos
      amu(i)=-1+(i-1)*dmu
   enddo
   call contio(nin,0,0,a,nb,nw)
   awr=c2h
   lct=l2h
   ! force laboratory coordinate system
   a(4)=1
   a(5)=nk
   write(nsyso,'(/'' converting mf=6, mt='',i3,&
     &'' to law 7 format'')') mth
   if (nk.gt.1) write(nsyso,'(i6,'' neutron subsections found'')') nk
   call contio(0,nout,ndebug,a,nb,nw)

   !--loop over the desired subsections
   ik=0
   do while (ik.lt.nk)
      ik=ik+1
      call tab1io(nin,0,0,a,nb,nw)
      a(4)=7
      call tab1io(0,nout,ndebug,a,nb,nw)
      do while (nb.ne.0)
         call moreio(nin,nout,ndebug,a,nb,nw)
      enddo
      call tab2io(nin,0,0,a,nb,nw)
      lang=l1h
      lep=l2h
      ne=n2h
      a(3)=0
      a(4)=0
      call tab2io(0,nout,ndebug,a,nb,nw)
      l1=1

      !--loop over incident energies
      do i=1,ne
         l=l1
         call listio(nin,0,0,a(l),nb,nw)
         l=l+nw
         do while (nb.ne.0)
            if (l.gt.namax) call error('fix6','storage exceeded',' ')
            call moreio(nin,0,0,a(l),nb,nw)
            l=l+nw
         enddo
         l2=l
         ein=a(l1+1)
         nep=nint(a(l1+5))
         nw=nint(a(l1+4))
         ncyc=nw/nep
         nl=ncyc-1
         nmu=ncos
         if (nl.eq.1) nmu=2

         !--write the tab2 for the cosine grid
         a(l2)=0
         a(l2+1)=a(l1+1)
         a(l2+2)=0
         a(l2+3)=0
         a(l2+4)=1
         a(l2+5)=nmu
         a(l2+6)=nmu
         a(l2+7)=2
         nw=8
         call tab2io(0,nout,ndebug,a(l2),nb,nw)

         !--loop over the cosines for law 7
         do imu=1,nmu
            if (nl.eq.1) then
               if (imu.eq.1) acos=-1
               if (imu.eq.2) acos=+1
            else
               acos=amu(imu)
            endif

            !--reconstruct the energy distribution for this cosine
            a(l2)=0
            a(l2+1)=acos
            a(l2+2)=0
            a(l2+3)=0
            a(l2+4)=1
            a(l2+5)=nep
            a(l2+6)=nep
            a(l2+7)=lep
            j=l2+8
            idone=0
            iep=0
            do while (idone.eq.0)
               iep=iep+1
               ep=a(l1+6+ncyc*(iep-1))
               csn=acos
               elb=ep
               drv=1
               if (lct.ne.1) then
                  aw1=awr+1
                  clb=csn
                  ! minimum lab cosine (= zero particle energy in cm)
                  cmn=-1
                  qq=1-aw1*aw1*ep/ein
                  if (qq.gt.zero) cmn=sqrt(qq)
                  if (clb.lt.cmn) then
                     clb=cmn
                     ! zero distribution when far outside
                     ! valid cosine range
                     if (imu.lt.nmu.and.amu(imu+1).le.cmn) drv=0
                  endif
                  ! outgoing particle energy in the laboratory system
                  qq=ep-ein*(1-clb*clb)/(aw1*aw1)
                  if (qq.lt.zero) qq=0
                  elb=clb*sqrt(ein)/aw1+sqrt(qq)
                  elb=elb*elb
                  ! calculate corresponding cosine in the cm system
                  if (ep.gt.zero) then
                     csn=clb*sqrt(elb/ep)-sqrt(ein/ep)/aw1
                  endif
                  ! check the limits
                  qq=-1
                  if (csn.lt.qq) csn=qq
                  qq=1
                  if (csn.gt.qq) csn=qq
                  ! include jacobian for cm-to-lab transformation
                  if (ep.ne.zero) drv=drv*sqrt(elb/ep)
               endif
               ! calculate the probability from Legendre polynomials
               if (lang.le.10) then
                  call legndr(csn,p,nl)
                  fmu=0
                  do ll=1,nl
                     fmu=fmu+(2*ll-1)*p(ll)*a(l1+6+ncyc*(iep-1)+ll)/2
                  enddo
               ! calculate the probability from
               ! pointwise representation
               else
                  mmu=(ncyc-1)/2
                  fmu=0
                  do kmu=1,mmu-1
                     ll=l1+5+2*kmu+ncyc*(iep-1)
                     if (csn.ge.a(ll).and.csn.le.a(ll+2))&
                       call terp1(a(ll),a(ll+1),a(ll+2),a(ll+3),&
                       csn,fmu,lang-10)
                  enddo
               endif
               if (j.le.l2+8) then
                  a(j)=elb
                  a(j+1)=fmu*drv
                  j=j+2
               else if (elb.gt.a(j-2)) then
                  a(j)=elb
                  a(j+1)=fmu*drv
                  j=j+2
               endif
               if (j.ge.namax-1) call error('fix6',&
                 'storage in a exceeded',' ')
               if (iep.eq.nep) idone=1
            enddo
            nnep=(j-(l2+8))/2
            if (nnep.eq.1) then
               a(l2+10)=2*a(l2+8)
               a(l2+11)=0
               nnep=2
            endif
            a(l2+5)=nnep
            a(l2+6)=nnep
            j2=l2
            call tab1io(0,nout,ndebug,a(j2),nb,nw)
            do while (nb.ne.0)
               j2=j2+nw
               call moreio(0,nout,ndebug,a(j2),nb,nw)
            enddo

         !--continue the loop over cosines
         enddo

      !--continue the loop over incident energies
      enddo

   !--continue the loop over subsections
   enddo

   !--finished with this section
   call tosend(nin,nout,ndebug,a)
   return
   end subroutine fix6

   subroutine gamsum(npend,nout,nf12c,matd,iopp)
   !-------------------------------------------------------------------
   ! Sum all gamma production cross sections computed from
   ! MF12*MF3 or MF16*MF3, or read from MF13, on unionized grid
   ! used for neutron cross sections.  Write sum as MF13/MT1.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util ! provides repoz,finda,loada,error
   use endf ! provides endf routines and variables
   ! externals
   integer::npend,nout,nf12c,matd,iopp
   ! internals
   integer::jscr,nin,isave,jt,j,nen,kw,mfd,nb,nw,i16,ntape
   integer::l,nk,mtd,mtx,ik,jdis,idone,idis,jw,k,kk,jscr2,i
   integer::ii,it,nw2,nsave,nwscr,nwscr2,ntemp
   real(kr)::e,enxt,x,thresh,y,enext
   real(kr)::c(2)
   real(kr),dimension(:),allocatable::scr,scr2
   real(kr),dimension(:),allocatable::bufn,bufo
   integer,parameter::nbuf=1000
   real(kr),parameter::eps=1.e-10_kr

   !--set up units and storage
   call repoz(nf12c)
   allocate(bufo(nbuf))
   allocate(bufn(nbuf))
   nwscr=2*(npage+50)
   allocate(scr(nwscr))
   nwscr2=npage+50
   allocate(scr2(nwscr2))
   jscr=1+npage+50
   nin=npend
   nscr=13

   !--use the energy grid from unionx
   isave=iold
   iold=inew
   inew=isave
   jt=1
   j=0
   do while (jt.gt.0)
      j=j+1
      call finda(j,c,2,inew,bufn,nbuf)
      c(2)=0
      jt=j
      if (c(1).ge.(1-eps)*elast) jt=-j
      call loada(jt,c,2,iold,bufo,nbuf)
   enddo
   nen=j
   kw=0

   !--sum all reactions found in mf12, mf16, or mf13.
   !--reconstruct redundant cross sections if needed.
   mfd=3
   call repoz(nf12c)
   call tpidio(nf12c,0,0,scr,nb,nw)
   i16=0
   ntape=0
   l=1
   if (mf1x(1).eq.0.and.iopp.ne.0) write(nsyso,&
     '(/'' message from gamsum---file 12 not found.'')')
   if (mf1x(1).eq.0) go to 240
   call findf(matd,12,0,nf12c)
  150 continue
   call contio(nf12c,0,0,scr,nb,nw)
   if (mfh.eq.16) i16=1
   if (mfh.eq.0) go to 240
   if (mth.eq.0) go to 150
   if (mfh.eq.12.and.mth.eq.460) then
      call tosend(nf12c,0,0,scr)
      go to 150
   endif
   nk=n1h
   if (iopp.gt.0) ntrpp=ntrpp+nk
   mtd=mth
   mtx=mth
   if (mth.eq.3) mtd=1
  160 continue
   call findf(matd,mfd,mtd,nin)
  170 continue
   call contio(nin,0,ntape,scr,nb,nw)
   if (mfh.eq.0.or.mfh.gt.13) go to 245
   ik=0
   e=0
   if (mfd.eq.13) then
      nk=n1h
      if (iopp.gt.0) ntrpp=ntrpp+nk
      mtd=mth
   else
      call gety1(e,enxt,jdis,x,nf12c,scr)
   endif
   idone=0
   do while (idone.eq.0)
      call gety2(e,thresh,idis,y,nin,scr(jscr))
      if (iopp.ne.0) then
         jw=0
         k=0
         kk=0
         scr2(1)=scr(l)
         scr2(2)=scr(l+1)
         scr2(3)=scr(l+2)
         scr2(4)=scr(l+3)
         scr2(5)=1
         scr2(6)=0
         scr2(7)=0
         scr2(8)=2
         jscr2=7
      endif
      do i=1,nen
         call finda(i,c,2,iold,bufo,nbuf)
         if (c(1).ge.(1-eps)*thresh) then
            e=c(1)
            if (i.eq.1) e=e*(1+eps)
            if (i.eq.nen) e=e*(1-eps)
            if (mfd.eq.3) call gety1(e,enxt,jdis,x,nf12c,scr)
            call gety2(e,enext,idis,y,nin,scr(jscr))
            if (i.gt.1.and.e.lt.thresh*(1+eps)) y=0
            if (mfd.ne.13) then
               if (mtd.eq.2) c(2)=c(2)-x*y
               if (mtd.ne.2) c(2)=c(2)+x*y
            else
               if (nk.eq.1.or.ik.ne.0) c(2)=c(2)+y
               if (iopp.ne.0) then
                  y=sigfig(y,7,0)
                  k=k+1
                  kk=kk+1
                  scr2(jscr2+2*kk)=c(1)
                  scr2(jscr2+1+2*kk)=y
                  if (2*kk.ge.npage.or.(i.ge.nen.and.jw.ne.0)) then
                     jw=jw+1
                     kw=kw+1
                     nw=2*kk
                     if (jw.eq.1) nw=nw+8
                     write(ntemp) nw,(scr2(ii),ii=1,nw)
                     kk=0
                     jscr2=-1
                  endif
               endif
            endif
         endif
         it=i
         if (i.eq.nen) it=-it
         call loada(it,c,2,inew,bufn,nbuf)
      enddo
      if (mfd.ne.13.or.iopp.ne.1) then
         call tosend(nf12c,0,0,scr)
         idone=1
      else
         if (jw.ne.0) then
            call repoz(-ntemp)
            read(ntemp) nw2,(scr2(ii),ii=1,nw2)
         endif
         scr2(6)=k
         scr2(7)=k
         math=matd
         mfh=13
         mth=mtd
         call tab1io(0,0,nscr,scr2,nb,nw)
         do while (nb.ne.0)
            read(ntemp) nw2,(scr2(ii),ii=1,nw2)
            call moreio(0,0,nscr,scr2,nb,nw)
         enddo
         call repoz(-ntemp)
         if (nk.eq.1) then
            idone=1
         else
            ik=ik+1
            if (ik.eq.nk+1) then
               idone=1
            else
               isave=iold
               iold=inew
               inew=isave
               e=0
            endif
         endif
      endif
   enddo
   call tosend(nin,0,ntape,scr)
   isave=iold
   iold=inew
   inew=isave
   if (mfd.eq.13) go to 170
   if (mtx.eq.3.and.mtd.eq.1) go to 235
   go to 150
  235 continue
   mtd=2
   call findf(matd,12,mtx,nf12c)
   call contio(nf12c,0,0,scr,nb,nw)
   go to 160

   !--finished with file 12.  check for mf16.
  240 continue
  if (i16.eq.1) go to 242
   if (n16.eq.0) go to 242
   call findf(matd,16,0,nf12c)
   go to 150

   !--read and add file 13 using same loop
  242 continue
   if (mf1x(2).eq.0) go to 245
   mfd=13
   mtd=0
   if (iopp.eq.0) go to 160
   nscr=13
   if (nout.lt.0) nscr=-nscr
   call openz(nscr,1)
   call repoz(nscr)
   ntemp=18
   call openz(-ntemp,1)
   call repoz(-ntemp)
   nsc=0
   ntape=nscr
   nsave=nin
   nin=nf12c
   l=jscr
   go to 160
   ! copy file 12 from nf12c to nout
  245 continue
   if (mf1x(1).ne.0) then
      call repoz(nf12c)
      call findf(matd,12,0,nf12c)
      call contio(nf12c,nout,0,scr,nb,nw)
      call tofend(nf12c,nout,0,scr)
   endif
   if (mf1x(2).ne.0) nin=nsave
   if (kw.gt.0) call closz(-ntemp)

   !--write new mf13, mt1 on nout
   math=matd
   mfh=13
   mth=1
   nxc=nxc+1
   if (nxc.gt.nxcmax)&
     call error('gamsum','exceeded storage in dictionary.',' ')
   mfs(nxc)=mfh
   mts(nxc)=mth
   scr(1)=za
   scr(2)=awr
   scr(3)=0
   scr(4)=0
   scr(5)=1
   scr(6)=0
   call contio(0,nout,0,scr,nb,nw)
   scr(1)=0
   scr(2)=0
   scr(3)=0
   scr(4)=0
   scr(5)=1
   scr(6)=nen
   scr(7)=nen
   scr(8)=2
   j=8
   do i=1,nen
      call finda(i,c,2,iold,bufo,nbuf)
      j=j+2
      scr(j-1)=c(1)
      scr(j)=sigfig(c(2),7,0)
      if (i.gt.(npage/2)) then
         if (j.eq.npage.or.i.eq.nen) then
            call moreio(0,nout,0,scr,nb,j)
            j=0
         endif
      else
         if (i.ge.(npage/2).or.i.eq.nen) then
            call tab1io(0,nout,0,scr,nb,j)
            j=0
         endif
      endif
   enddo
   call asend(nout,0)
   if (mf1x(2).eq.0) call afend(nout,0)
   if (mf1x(2).ne.0) then
      ! copy file 13 from nscr to nout
      call afend(0,nscr)
      call repoz(nscr)
      call contio(nscr,nout,0,scr,nb,nw)
      call tofend(nscr,nout,0,scr)
      call closz(nscr)
   endif
   call closz(-iold)
   call closz(-inew)
   call repoz(nf12c)
   deallocate(scr2)
   deallocate(scr)
   deallocate(bufn)
   deallocate(bufo)
   return
   end subroutine gamsum

   subroutine convr(nin,npend,nout,nscr,nedis,nethr,matd)
   !-------------------------------------------------------------------
   ! Convert MF12 photon transition probability arrays (LO=2)
   ! to photon yields (LO=1).  Copy MF13.  Fix up MF14 for changes
   ! made to MF12.  copy MF15.  Move any MF6 photon production.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util ! provides openz,repoz,closz,sigfig
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,npend,nout,nscr,nedis,nethr,matd
   ! internals
   integer::kt,mf12,nw,i,l2flg,j,nk,mfd,mtd,ik,np,nwt,jtot
   integer::jp,jpn,jpp
   integer::nb,l,iedis,nemax,imu,nmu,ie,ne,nbt
   integer::imax,lmax,imaxsq,nwscr,kscr,idis,kprnt
   integer::lg,mt0,mt0old,n,jm1,kk,k,idone,ii,kp1,ja,jb
   integer::lm1,ip1,mtl,igm,law,mf,nnth,jdis
   integer::m1,m2,maths,mtnow,mttst,nn
   integer,parameter::nqmx=450
   real(kr)::e,e1,y,enext,test,diff,g,ei,p,ysum,eja,ejb,save,zap,egamma
   real(kr)::enxt,x,elow
   integer::ngam(300)
   integer::mtth(nqmx)
   real(kr)::eeth(nqmx)
   real(kr)::gamid(17)
   character(4)::t(17)
   character(60)::strng
   equivalence(t(1),gamid(1))
   real(kr),dimension(:),allocatable::ee,eg,es,yy,aa,rr
   real(kr),dimension(:),allocatable::tot
   integer,parameter::nwtot=5000
   real(kr),dimension(:),allocatable::yold
   integer,parameter::nyold=1000
   real(kr),dimension(:),allocatable::scr
   character(4)::blank='    '
   real(kr),parameter::etop=1.e10_kr
   real(kr),parameter::emax=2.e7_kr
   real(kr),parameter::zero=0

   !--set up size and storage, and assign i/o units
   nout=iabs(nout)
   if (nin.lt.0) nout=-nout
   call openz(nout,1)
   mt0=0
   mt0old=0
!  imax=49
   imax=50
   imaxsq=imax*imax
   lmax=100
   nned=50
   allocate(disc(nned))
   allocate(ee(imax))
   allocate(eg(lmax))
   allocate(es(lmax))
   allocate(yy(lmax))
   allocate(aa(imaxsq))
   allocate(rr(imaxsq))
   nwscr=npage+50
   allocate(scr(nwscr))
   allocate(yold(nyold))
   allocate(tot(nwtot))
   n16=0
   nf16s=0
   kscr=inew
   if (nin.lt.0) kscr=-kscr
   kt=0
   mf12=mf1x(1)

   !--write a tapeid on nout
   math=1
   mfh=0
   mth=0
   nw=17
   do i=1,nw
      read(blank,'(a4)') t(i)
   enddo
   call tpidio(0,nout,0,gamid,nb,nw)

   !--initialize
   l2flg=0
   do i=1,300
      ngam(i)=0
   enddo
   ee(1)=0
   do i=1,imax
      do j=1,imax
         rr((j-1)*imax+i)=0
         if (i.eq.j) rr((j-1)*imax+i)=1
         aa((j-1)*imax+i)=0
      enddo
   enddo
   nsh=0
   nsc=0
   if (mf1x(1).eq.0.and.mf1x(2).eq.0) go to 515

   !--get thresholds vs mt number
   call findf(matd,3,0,npend)
   nnth=0
  101 call contio(npend,0,0,scr,nb,nw)
   if (mfh.eq.0) go to 102
   e=0
   call gety1(e,enxt,jdis,x,npend,scr)
   nnth=nnth+1
   eeth(nnth)=enxt
   mtth(nnth)=mth
   call tosend(npend,0,0,scr)
   go to 101
  102 continue
   if (mf12.eq.0) go to 92

   !--check for missing mf12 mt's ... advise user if found
   call findf(matd,12,0,nin)
   mtnow=-1
   maths=matd
  90 continue
   call contio(nin,0,0,scr,nb,nw)
   if (mfh.eq.0.or.math.eq.0.or.math.eq.-1) go to 91
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
         write(strng,'('' mf12, mt'',i2,'' may be missing'')')mth-1
         call mess('convr',strng,&
                           'discrete photon data may be incomplete')
      endif
      mtnow=mth
   endif
   call tosend(nin,0,0,scr)
   go to 90
  91 continue
   call repoz(nin)
   matd=maths

   !--loop over sections of mf12, mf13, and mf14
  92 continue
   mf=12
   if (mf12.eq.0) mf=13
   call findf(matd,mf,0,nin)
  110 continue
   call contio(nin,0,0,scr,nb,nw)
   if (mfh.eq.0.and.math.ne.0) go to 190
   if (mfh.eq.13) go to 125
   if ((mfh.eq.12.and.mth.eq.460).or.(mfh.eq.14.and.mth.eq.460)) then
       call tosend(nin,0,0,scr)
       go to 110
   endif
   if (mfh.eq.12) go to 120
   if (mfh.eq.14.and.l2flg.eq.1) go to 500
   if (mfh.gt.15) go to 515
   call contio(0,nout,0,scr,nb,nw)
   if (math.eq.0.or.mfh.gt.15) go to 515
   call tofend(nin,nout,nscr,scr)
   go to 110
  120 continue
   if (l1h.eq.2) go to 200
  125 continue
   call contio(0,nout,nscr,scr,nb,nw)
   nk=n1h
   mfd=mfh
   mtd=mth
   if (nk.eq.1) go to 180
   kt=kt+1
   if (kt.eq.1) call openz(kscr,1)
   ik=0
   ! copy total yield or xs to a scratch file
   call repoz(kscr)
   nsc=0
   call tab1io(nin,0,kscr,scr,nb,nw)
   np=n2h
   nr=n1h
   nwt=6+2*nr+2*np
   if (nwt.gt.nwtot)&
     call error('convr','storage exceeded for photon data',' ')
   call tosend(nin,0,kscr,scr)
   call repoz(kscr)
   jtot=1
   call tab1io(kscr,0,0,tot(jtot),nb,nw)
   do while (nb.ne.0)
      jtot=jtot+nw
      call moreio(kscr,0,0,tot(jtot),nb,nw)
      if (jtot.gt.nwtot) call error('convr',&
        'exceeded tot array',' ')
   enddo
   ! supplement grid with discontinuities
   l=6+2*nr
   do i=1,np
      if (i.eq.1.and.mfh.eq.13) then
         ethr(1+nethr)=tot(l+2*i-1)
         nethr=nethr+1
      else if (i.ne.1.and.tot(l+2*i-1).eq.tot(l+2*i-3)) then
         nedis=nedis+1
         if (nedis.gt.nned)&
           &call error('convr','storage exceeded for edis',' ')
         disc(nedis)=tot(l+2*i-3)
         tot(l+2*i-3)=sigfig(tot(l+2*i-3),7,-1)
         tot(l+2*i-1)=sigfig(tot(l+2*i-1),7,1)
      endif
      yold(i)=tot(l+2*i)
      tot(l+2*i)=0
   enddo
   ! compute corrected total on the new grid
   do while (ik.lt.nk)
      e=0
      call gety1(e,e1,idis,y,kscr,scr)
      enext=e1
      test=etop-etop/100
      do i=1,np
         e=tot(l+2*i-1)
         if (e.ge.e1.and.enext.le.test) then
            call gety1(e,enext,idis,y,kscr,scr)
            tot(l+2*i)=tot(l+2*i)+y
         endif
      enddo
      ik=ik+1
   enddo
   ! write corrected total to nout and nscr
   mfh=mfd
   mth=mtd
   jtot=1
   call tab1io(0,nout,nscr,tot(jtot),nb,nw)
   do while (nb.ne.0)
      jtot=jtot+nw
      call moreio(0,nout,nscr,tot(jtot),nb,nw)
   enddo
   call repoz(kscr)
   call tab1io(kscr,0,0,scr,nb,nw)
   do while (nb.ne.0)
      call moreio(kscr,0,0,scr,nb,nw)
   enddo
   ! copy rest of mt to nout and nscr
   call tab1io(kscr,nout,nscr,scr,nb,nw)
   call tosend(kscr,nout,nscr,scr)
   ! check for magnitude of differences found
   kprnt=0
   do i=1,np
      if (yold(i).ne.zero) then
         diff=(abs(tot(l+2*i)-yold(i))/yold(i))*100
         test=1
         test=test/10
         if (diff.gt.test) then
            kprnt=kprnt+1
            if (kprnt.eq.1) write(nsyso,'(/)')
            write(nsyso,'('' ('',i2,''/'',i3,'') at energy'',&
              &1p,e13.5,'', endf yield='',e13.5,'', correct yield='',&
              &e13.5,'', p.c. diff='',e13.5)')&
              mfd,mtd,tot(l+2*i-1),yold(i),tot(l+2*i),diff
         endif
      endif
   enddo
   go to 110
   ! handle the case nk=1 for mf12 or mf13
  180 continue
   jtot=1
   call tab1io(nin,nout,nscr,tot(jtot),nb,nw)
   np=n2h
   nr=n1h
   nwt=6+2*nr+2*np
   if (nwt.gt.nwtot)&
     call error('convr','storage exceeded for photon data',' ')
   do while (nb.ne.0)
      jtot=jtot+nw
      call moreio(nin,nout,nscr,tot(jtot),nb,nw)
   enddo
   nr=n1h
   if (mf.eq.13) then
      ethr(1+nethr)=tot(7+2*nr)
      nethr=nethr+1
   endif
   ! check for discontinuities
   l=6+2*nr
   do i=1,np
      if (i.ne.1.and.tot(l+2*i-1).eq.tot(l+2*i-3)) then
         nedis=nedis+1
         if (nedis.gt.nned)&
           call error('convr','storage exceeded for edis',' ')
         disc(nedis)=tot(l+2*i-3)
      endif
   enddo
   call tosend(nin,nout,nscr,scr)
   go to 110
  190 continue
   call contio(0,nout,nscr,scr,nb,nw)
   go to 110

   !--convert transition probability array to yields
  200 continue
   za=scr(1)
   awr=scr(2)
   lg=l2h
   g=1
   l2flg=1
   call listio(nin,0,0,scr,nb,nw)

   !--make sure mt0 is correct for this range of mt's.
   if (mth.ge.51.and.mth.le.90.and.mt0.ne.49) mt0=49
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

   !--load the ee(_) array with mf3 -q.  Will overwrite with mf12
   !  mt data when possible (which should be always but if there is
   !  a missing mf12 mt, we're covered).
   if (mt0.ne.mt0old) then
      ee=0
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
      do while (mttst.lt.m2.and.nn.le.nnth)
         mttst=mtth(nn)
         if (mttst.ge.m1.and.mttst.le.m2)&
                                 ee(mttst-mt0)=eeth(nn)*awr/(awr+1)
         nn=nn+1
      enddo
   endif

   j=mth-mt0
   ee(j)=scr(1)
   n=nint(scr(6))
   jm1=j-1
   do kk=1,jm1
      k=jm1-kk+1
      idone=0
      ii=0
      do while (ii.lt.n.and.idone.eq.0)
         ii=ii+1
         i=ii
         ei=scr(6+(lg+1)*i-lg)
         if (ei.eq.zero.and.ee(k).eq.zero) idone=1
         if (ei.ne.zero) then
            if (abs(ei-ee(k))/ei.lt.0.0001) idone=1
         endif
      enddo
      if (idone.eq.0) then
         aa((k-1)*imax+j)=0
         rr((k-1)*imax+j)=0
      else
         p=scr(7+(lg+1)*i-lg)
         g=1
         if (lg.eq.2) g=scr(6+(lg+1)*i)
         aa((k-1)*imax+j)=p*g
         rr((k-1)*imax+j)=p
      endif
      if (k.ne.jm1) then
         kp1=k+1
         do ii=kp1,jm1
            rr((k-1)*imax+j)=rr((k-1)*imax+j)&
              +rr((ii-1)*imax+j)*aa((k-1)*imax+ii)
         enddo
      endif
   enddo
   l=0
   ysum=0
   do i=2,j
      ja=j+2-i
      eja=ee(ja)
      do ii=1,jm1
         jb=j-ii
         ejb=ee(jb)
         y=aa((jb-1)*imax+ja)*rr((ja-1)*imax+j)
         if (y.ne.zero) then
            l=l+1
            if (l.gt.lmax) call error('convr',&
              'too many lo=2 photons',' ')
            eg(l)=eja-ejb
            es(l)=eja
            yy(l)=y
            ysum=ysum+yy(l)
         endif
      enddo
   enddo
   ee(1)=0

   !--arrange gamma energies in descending order
   if (l.gt.1) then
      lm1=l-1
      do i=1,lm1
         ip1=i+1
         do ii=ip1,l
            if (eg(i).lt.eg(ii)) then
               save=eg(i)
               eg(i)=eg(ii)
               eg(ii)=save
               save=yy(i)
               yy(i)=yy(ii)
               yy(ii)=save
               save=es(i)
               es(i)=es(ii)
               es(ii)=save
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
   elow=0
   do i=1,nnth
      if (mtth(i).eq.mth) elow=eeth(i)
   enddo
   ! output tab1 sum record
   if (l.gt.1) then
      scr(1)=0
      scr(2)=0
      scr(3)=0
      scr(5)=1
      scr(6)=2
      scr(7)=2
      scr(8)=2
      scr(9)=elow
      scr(10)=ysum
      scr(11)=elim
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
      scr(9)=elow
      scr(10)=yy(i)
      scr(11)=elim
      scr(12)=yy(i)
      nw=12
      call tab1io(0,nout,nscr,scr,nb,nw)
   enddo

   !--copy send record and loop over remaining sections.
   call contio(nin,nout,nscr,scr,nb,nw)
   go to 110

   !--angular distributions
  500 continue
   mt0=0
   if (mth.ge.50.and.mth.le.90) mt0=50
   if (iverf.le.5.and.mth.ge.700.and.mth.le.717) mt0=700
   if (iverf.le.5.and.mth.ge.720.and.mth.le.737) mt0=720
   if (iverf.le.5.and.mth.ge.740.and.mth.le.757) mt0=740
   if (iverf.le.5.and.mth.ge.760.and.mth.le.777) mt0=760
   if (iverf.le.5.and.mth.ge.780.and.mth.le.797) mt0=780
   if (iverf.ge.6.and.mth.ge.600.and.mth.le.648) mt0=600
   if (iverf.ge.6.and.mth.ge.650.and.mth.le.698) mt0=650
   if (iverf.ge.6.and.mth.ge.700.and.mth.le.748) mt0=700
   if (iverf.ge.6.and.mth.ge.750.and.mth.le.798) mt0=750
   if (iverf.ge.6.and.mth.ge.800.and.mth.le.848) mt0=800
   if (mt0.ne.0) then
      mtl=49
      if (iverf.le.5.and.mth.ge.700) mtl=649
      if (iverf.ge.6.and.mth.ge.600) mtl=549
      scr(3)=1
      scr(4)=0
      scr(6)=0
      scr(5)=ngam(mth-mtl)
      call contio(0,nout,nscr,scr,nb,nw)
      call tosend(nin,0,0,scr)
      call asend(nout,nscr)
   else
      call contio(0,nout,nscr,scr,nb,nw)
      call tosend(nin,nout,nscr,scr)
   endif
   go to 110

   !--move mf6 photons to mf16
  515 continue
   if (nsix.ne.0) then
      call findf(matd,6,0,nin)
      mfh=6
      do while (mfh.ne.0)
         call contio(nin,0,0,scr,nb,nw)
         jp=l1h
         jpn=mod(jp,10)
         jpp=jp-10*jpn
         if (mfh.ne.0 .and. mth.eq.18 .and. jpp.ne.0) then
            call tosend(nin,0,0,scr) !skip past mf6/mt18 (for now)
         else if (mfh.ne.0) then
            igm=0
            nk=n1h
            ik=0
            do while (ik.lt.nk)
               call tab1io(nin,0,0,scr(7),nb,nw)
               zap=c1h
               law=l2h
               egamma=0
               if (zap.eq.zero .and. law.ne.zero) then
                  if (law.eq.2) then
                     call mess('convr',&
                       'discrete anisotropic photon',&
                       'treated as simple primary photon')
                     egamma=c2h
                  endif
                  scr(5)=1
                  mfh=16
                  n16=n16+1
                  if (mth.ge.600) nf16s=nf16s+1
                  if (mth.ge.600) mf16s(nf16s)=mth
                  ngmt=ngmt+1
                  gmt(ngmt)=1000*mfh+mth
                  igm=1+(nw+5)/6
                  call contio(0,nout,0,scr,nbt,nw)
                  call tab1io(0,nout,0,scr(7),nbt,nw)
               endif
               do while (nb.ne.0)
                  call moreio(nin,0,0,scr(7),nb,nw)
                  if (zap.eq.zero .and. law.ne.zero) then
                     mfh=16
                     call moreio(0,nout,0,scr(7),nbt,nw)
                  endif
               enddo
               if (law.eq.1.or.law.eq.2.or.law.eq.5) then
                  call tab2io(nin,0,0,scr(7),nb,nw)
                  ne=n2h
                  mfh=16
                  if (zap.eq.zero)&
                    call tab2io(0,nout,0,scr(7),nbt,nw)
                  if (zap.eq.zero) igm=igm+(nw+5)/6
                  do ie=1,ne
                     call listio(nin,0,0,scr(7),nb,nw)
                     mfh=16
                     if (egamma.gt.0) then
                        e=c2h
                        scr(9)=1
                        scr(10)=0
                        scr(11)=2
                        scr(12)=1
                        scr(13)=egamma+awr*e/(awr+1)
                        scr(14)=1
                     endif
                     if (zap.eq.0.)&
                       call listio(0,nout,0,scr(7),nbt,nw)
                     if (zap.eq.0.) igm=igm+(nw+5)/6
                     do while (nb.ne.0)
                        call moreio(nin,0,0,scr(7),nb,nw)
                        if (zap.eq.zero) then
                           mfh=16
                           call moreio(0,nout,0,scr(7),nbt,nw)
                        endif
                     enddo
                  enddo
               else if (law.eq.6) then
                  call contio(nin,0,0,scr(7),nb,nw)
               else if (law.eq.7) then
                  call tab2io(nin,0,0,scr(7),nb,nw)
                  ne=n2h
                  do ie=1,ne
                     call tab2io(nin,0,0,scr(7),nb,nw)
                     nmu=n2h
                     do imu=1,nmu
                        call tab1io(nin,0,0,scr(7),nb,nw)
                        do while (nb.ne.0)
                           call moreio(nin,0,0,scr(7),nb,nw)
                        enddo
                     enddo
                  enddo
               endif
               ik=ik+1
            enddo
            call tosend(nin,0,0,scr)
            if (igm.ne.0) then
               mfh=16
               call asend(nout,0)
            endif
         endif
      enddo
      call contio(0,nout,0,scr,nb,nw)
   endif

   !--convr is finished.
   call amend(nout,0)
   call atend(nout,0)
   ! put energies in order and eliminate duplicate points
   if (nedis.gt.1) call aordr(nedis,nned,disc)
   if (nethr.gt.1) call aordr(nethr,300,ethr)
   nemax=nedis
   if (nemax.gt.nned) call error('convr','storage exceeded.',' ')
   if (kt.gt.0) call closz(kscr)
   if (nedis.gt.0) write(nsyso,'(/&
     &'' energy discontinuities found in gamma files''/&
     &(10x,1p,e13.5))')&
     (disc(iedis),iedis=1,nedis)
   deallocate(tot)
   deallocate(yold)
   deallocate(scr)
   deallocate(rr)
   deallocate(aa)
   deallocate(yy)
   deallocate(es)
   deallocate(eg)
   deallocate(ee)
   return
   end subroutine convr

   subroutine aordr(n,max,e)
   !-------------------------------------------------------------------
   ! Put the array e in ascending order, and eliminate duplicate
   ! elements.  The value of n may be changed on return.
   !-------------------------------------------------------------------
   ! externals
   integer::n,max
   real(kr)::e(max)
   ! internals
   integer::nm1,k,i,i1,j
   real(kr)::save
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::eps=1.e-10_kr

   if (n.lt.2) return
   nm1=n-1
   k=0
   do i=1,nm1
      i1=i+1
      do j=i1,n
         if (e(j).le.(1+eps)*e(i)) then
            if (e(j).lt.(1-eps)*e(i)) then
               save=e(i)
               e(i)=e(j)
               e(j)=save
            else if (e(j).lt.emax) then
               e(j)=emax
               k=k+1
            endif
         endif
      enddo
   enddo
   n=n-k
   return
   end subroutine aordr

   subroutine gamout(ngend,nendf,nout,nf12c,matd)
   !-------------------------------------------------------------------
   ! Write Files 14 and 15.  If multigroup photon data are available,
   ! create an isotropic MF14/MT1.  Read in the multigroup energy
   ! distributions, convert them into equally probable photons, and
   ! write them out as MF15/MT1 using a specially defined law LF3.
   ! Convert the normal ENDF photon angular distributions from MF14
   ! into equal-probability bins in a new version of MF14.  copy
   ! the normal ENDF MF15 and the specially converted MF16 photons.
   !-------------------------------------------------------------------
   use util ! provides error,repoz,sigfig
   use endf ! provides endf routines and variables
   ! externals
   integer::ngend,nendf,nout,nf12c,matd
   ! internals
   integer::nin,nb,nw
   integer::li,ltt,nk,ik,ne,ie,now,n15,ip,i,ig,ngn,nl,nz,ng,ng2
   integer::ig2lo,nwa,lz,ibase,ngg,ngnp1,nggp1,jg,jpos,loca,locs
   integer::iabase,ia,ign,jgghi,jgg,jb,ncds,n,nbinp1,nbinp6
   integer::nsmax,nwamax,nwords,nbin,ni
   real(kr)::sum,sumn,rsumn,elast,r,s,dist
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::egn,egg
   real(kr),dimension(:),allocatable::ee,eb
   real(kr),dimension(:),allocatable::sig
   integer,parameter::maxsig=10000
   real(kr),parameter::elo=1.e-5_kr
   real(kr),parameter::ehi=2.e7_kr
   real(kr),parameter::zero=0

   nin=ngend
   !--set up storage
   nbin=nbinp
   nbinp1=nbin+1
   nbinp6=nbin+6
   nwords=241
   allocate(egn(nwords))
   nwords=151
   allocate(egg(nwords))
   nwamax=3*npage+50
   allocate(scr(nwamax))
   nsmax=5000
   allocate(sig(nsmax))

   !--write the isotropic mf14/mt1 for the multigroup photons
   if (nin.ne.0) then
      scr(1)=za
      scr(2)=awr
      scr(3)=1
      scr(4)=0
      scr(5)=1
      scr(6)=0
      math=matd
      mfh=14
      mth=1
      call contio(0,nout,0,scr,nb,nw)
      nxc=nxc+1
      mfs(nxc)=mfh
      mts(nxc)=mth
      ncs(nxc)=1
      call asend(nout,0)
   endif

   !--convert the endf mf14 distributions to equal-prob bins
   if (mf1x(1).eq.0.and.mf1x(2).eq.0) go to 550
   call findf(matd,14,0,nf12c)
   mfh=14
   do while (mfh.ne.0)
      call contio(nf12c,nout,0,scr,nb,nw)
      if (mfh.ne.0) then
         li=nint(scr(3))
         if (li.ne.0) then
            ! this entire reaction is isotropic
            call tosend(nf12c,nout,0,scr)
         else
            ! this reaction contains anisotropic photons
            ltt=nint(scr(4))
            nk=nint(scr(5))
            ni=nint(scr(6))
            do ik=1,nk
               if (ik.le.ni) then
                  ! this subsection is isotropic
                  call contio(nf12c,nout,0,scr,nb,nw)
               else
                  ! this subsection is anisotropic.  convert it.
                  call tab2io(nf12c,nout,0,scr,nb,nw)
                  ne=nint(scr(6))
                  do  ie=1,ne
                     if (ltt.ne.2) then
                        call listio(nf12c,0,0,scr,nb,nw)
                        call ptleg(nout,scr)
                     else
                        call tab1io(nf12c,0,0,scr,nb,nw)
                        now=1+nw
                        do while (nb.ne.0)
                           call moreio(nf12c,0,0,scr(now),nb,nw)
                           now=now+nw
                        enddo
                        call pttab(ltt,scr,nout)
                     endif
                  enddo
               endif
            enddo
            call contio(nf12c,nout,0,scr,nb,nw)
            if (mth.ne.0) call error('gamout',&
              'expected send card while reading mf14',' ')
         endif
      endif
   enddo

   !--check for multigroup mf15 data
   n15=mf1x(3)
   if (nin.eq.0) go to 500
   ip=1
   do i=1,nsmax
      sig(i)=0
   enddo

   !--find desired material on input tape.
   call repoz(nin)
   call tpidio(nin,0,0,scr,nb,nw)
   ig=0
   ngn=1
  110 continue
   call contio(nin,0,0,scr,nb,nw)
   if (math.eq.1) go to 110
   if (math.eq.-1) call error('gamout','mat not found.',' ')
   if (math.eq.matd) go to 116
   call tomend(nin,0,0,scr)
   go to 110

   !--process all photon production sections.
  115 continue
   if (ig.lt.ngn) go to 120
   call contio(nin,0,0,scr,nb,nw)
   if (math.eq.0) go to 300
   if (mfh.eq.0) go to 115
   if (mth.eq.0) go to 115
  116 continue
   nl=l1h
   nz=l2h
   ng=n2h
  120 continue
   call listio(nin,0,0,scr,nb,nw)
   ng2=l1h
   ig2lo=l2h
   nw=n1h
   ig=n2h
   nwa=nw+1
   lz=6
   do while (nb.ne.0)
      call moreio(nin,0,0,scr(nwa),nb,nw)
      nwa=nwa+nw
      if (nwa.gt.nwamax)&
         call error('gamout','storage in a exceeded',' ')
      enddo

   !--store information from header record.
   if (mfh.eq.1) then
      ibase=lz+ng+nz
      ngn=ng2
      ngg=ig2lo
      if (ngg.eq.0)&
        call error('gamout','no gamma groups on ngend.',' ')
      if ((ngn*ngg).gt.nsmax)&
        call error('gamout','storage in sig exceeded.',' ')
      ngnp1=ngn+1
      do i=1,ngnp1
         egn(i)=scr(i+ibase)
      enddo
      ibase=ibase+ngnp1
      nggp1=ngg+1
      do i=1,nggp1
         egg(i)=scr(i+ibase)
      enddo
      ig=ngn
      allocate(ee(nbinp1))
      allocate(eb(nbinp1))

   !--sum all gamma production matrices.
   else if (mfh.eq.16.or.mfh.eq.17) then
      jg=ig-1
      if (jg.gt.0) then
         do i=2,ng2
            jpos=ig2lo+i-2
            if (jpos.ge.1.and.jpos.le.ngg) then
               loca=lz+ip+nl*nz*(i-1)
               locs=jpos+ngg*jg
               sig(locs)=sig(locs)+scr(loca)
               if (locs.gt.nsmax) call error('gamout',&
                 'exceeded size of sig array',' ')
            endif
         enddo
      endif
   endif
   go to 115

   !--compute break energies and equal probability bins.
  300 continue
   iabase=13+nbinp6
   ia=iabase
   ee(1)=egg(1)
   do 310 ign=1,ngn
   sum=0
   jpos=ngg*(ign-1)
   jgghi=1
   do jgg=1,ngg
      locs=jgg+jpos
      sum=sum+sig(locs)
      if (sig(locs).ne.zero) jgghi=jgg
   enddo
   if (sum.eq.zero) go to 344
   sumn=sum/nbin
   rsumn=1/sumn
   jgg=1
   elast=ee(1)
   do 330 jb=2,nbinp1
   r=sumn
   s=0
  325 continue
   locs=jgg+jpos
   dist=sig(locs)/(egg(jgg+1)-egg(jgg))
   if (jgg.eq.jgghi) go to 335
   if ((dist*(egg(jgg+1)-elast)).lt.r) go to 340
  335 continue
   ee(jb)=elast+r/dist
   eb(jb-1)=rsumn*(s+dist*(ee(jb)*ee(jb)-elast*elast)/2)
   elast=ee(jb)
   go to 330
  340 continue
   r=r-dist*(egg(jgg+1)-elast)
   s=s+dist*(egg(jgg+1)*egg(jgg+1)-elast*elast)/2
   elast=egg(jgg+1)
   jgg=jgg+1
   go to 325
  330 continue
   ee(nbinp1)=egg(nggp1)
   do jb=1,nbin
      eb(jb)=sigfig(eb(jb),7,0)
   enddo
   go to 360
   ! if no photons for this group, set all eb to zero.
  344 continue
   do jb=1,nbin
      eb(jb)=0
   enddo
   ! save new array in list record format.
  360 continue
   scr(ia+1)=0
   scr(ia+2)=egn(ign)
   scr(ia+3)=0
   scr(ia+4)=0
   scr(ia+5)=nbin
   scr(ia+6)=0
   ia=ia+6
   do jb=1,nbin
      ia=ia+1
      if (ia.gt.nwamax)&
        &call error('gamout','storage in a exceeded.',' ')
      scr(ia)=eb(jb)
   enddo
  310 continue

   !--write file 15, mt 1
   mfh=15
   mth=1
   math=matd
   scr(1)=za
   scr(2)=awr
   scr(3)=0
   scr(4)=0
   scr(5)=1
   scr(6)=0
   nw=6
   call contio(0,nout,0,scr,nb,nw)
   ncds=1
   scr(1)=0
   scr(2)=0
   scr(3)=0
   scr(4)=0
   scr(5)=1
   scr(6)=2
   scr(7)=2
   scr(8)=2
   scr(9)=elo
   scr(10)=1
   scr(11)=ehi
   scr(12)=1
   nw=12
   call tab1io(0,nout,0,scr,nb,nw)
   ncds=ncds+3
   scr(1)=0
   scr(2)=0
   scr(3)=0
   scr(4)=0
   scr(5)=1
   scr(6)=ngn
   scr(7)=ngn
   scr(8)=1
   nw=8
   call tab2io(0,nout,0,scr,nb,nw)
   ncds=ncds+2
   ia=iabase
   do n=1,ngn
      do jb=1,nbinp6
         ia=ia+1
         scr(jb)=scr(ia)
      enddo
      nw=nbinp6
      call listio(0,nout,0,scr,nb,nw)
      ncds=ncds+1+(nbin+5)/6
   enddo
   nxc=nxc+1
   mfs(nxc)=mfh
   mts(nxc)=mth
   ncs(nxc)=ncds
   mf1x(3)=mf1x(3)+1
   call asend(nout,0)

   !--copy the endf mf15 data
  500 continue
   if (n15.eq.0) go to 550
   call findf(matd,15,0,nendf)
   call contio(nendf,nout,0,scr,nb,nw)
   call tofend(nendf,nout,0,scr)

   !--copy the mf16 data derived from the endf mf6, if any
  550 continue
   if (n16.ne.0) then
      call findf(matd,16,0,nf12c)
      call contio(nf12c,nout,0,scr,nb,nw)
      call tofend(nf12c,nout,0,scr)
   endif

   !--the gamma angle and energy distributions are ready
   call amend(nout,0)
   call atend(nout,0)
   if (allocated(ee)) deallocate(ee)
   if (allocated(eb)) deallocate(eb)
   deallocate(scr)
   deallocate(sig)
   deallocate(egn)
   deallocate(egg)
   return
   end subroutine gamout

   subroutine acelod(nin,nedis,suff,matd,tempd,newfor,mcnpx,ismooth)
   !-------------------------------------------------------------------
   ! Load data in ace format from the input file.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use physics ! provides pi,bk,amassn,amu,hbar,ev,clight
   use util ! repoz,dater,error,skiprz,sigfig
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nedis,matd,newfor,mcnpx,ismooth
   real(kr)::suff,tempd
   ! internals
   integer::nwscr,nnu,nnup,kfis,mtnr,mtntr,i,nnud,nnf
   integer::nurd,idone,mta,nb,nw,lnu,n,m,jnt,j
   integer::lssf,iinel,iabso,nunr,ncyc,i1,idis
   integer::k,it,ic,ie,ih,next,keep3,keep4,keep,ir,iskip
   integer::nnex,keep1,keep2,l,ij,lct,lvt,ltt,ltt3,lttn
   integer::jscr,iso,ne,law,lidp,last,il,ja,jb,nure,nurb
   integer::jj,ll,ib,iza,mf,mt,lend,lendp,inow
   integer::lff,lxx,nn,mm,iint,loc,ix
   integer::mt418,mt518
   real(kr)::urlo,urhi,e,enext,s,test,awp,spi,q,x,teste,zaid
   real(kr)::xxmin,xxmax,sumup,ex,fx,val
   character(8)::hdt
   real(kr),dimension(8)::dntc
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::nut,nup,nudn
   real(kr),dimension(:),allocatable::urd
   real(kr),parameter::eps=1.e-10_kr
   real(kr),parameter::emev=1.e6_kr
   real(kr),parameter::shake=1.e8_kr
   real(kr),parameter::etop=1.e10_kr
   real(kr),parameter::zero=0

   !--initialize
   inow=0
   nnex=0
   nwscr=1000000
   allocate(scr(nwscr))
   do i=1,8
      nxsd(i)=0
   enddo
   do i=1,2
      jxsd(i)=0
   enddo
   call repoz(nin)
   iza=nint(za)
   aw0=awr
   write(hm,'(''   mat'',i4)') matd
   tz=tempd*bk/emev
   zaid=iza+suff
   if (mcnpx.eq.0) then
      if (izai.eq.1) then
         write(hz,'(f9.2,''c'')') zaid
      else if (izai.eq.1001) then
         write(hz,'(f9.2,''h'')') zaid
      else if (izai.eq.1002) then
         write(hz,'(f9.2,''o'')') zaid
      else if (izai.eq.1003) then
         write(hz,'(f9.2,''r'')') zaid
      else if (izai.eq.2003) then
         write(hz,'(f9.2,''s'')') zaid
      else if (izai.eq.2004) then
         write(hz,'(f9.2,''a'')') zaid
      else
        call error('acelod',&
          'not coded for this incident particle',' ')
      endif
   else
      if (izai.eq.0) then
         write(hz,'(f10.3,''pc '')') zaid
      else if (izai.eq.1) then
         write(hz,'(f10.3,''nc '')') zaid
      else if (izai.eq.1001) then
         write(hz,'(f10.3,''hc '')') zaid
      else if (izai.eq.1002) then
         write(hz,'(f10.3,''dc '')') zaid
      else if (izai.eq.1003) then
         write(hz,'(f10.3,''tc '')') zaid
      else if (izai.eq.2003) then
         write(hz,'(f10.3,''sc '')') zaid
      else if (izai.eq.2004) then
         write(hz,'(f10.3,''ac '')') zaid
      else
        call error('acelod',&
          'not coded for this incident particle',' ')
      endif
   endif
   ntrp=ntrpp
   izaid=nint(za)
   call dater(hdt)
   hd='  '//hdt

   !--count reactions and set flags
   !--keep all reactions that survived unionx
   nnu=0
   nnud=0
   nnup=0
   ntr=0
   nr=0
   kfis=0
   mtnr=0
   mtntr=0
   do i=1,nxc
      mf=mfs(i)
      mt=mts(i)
      if (mf.eq.1.and.mt.eq.452) kfis=1
      if (mf.eq.1.and.mt.eq.455) kfis=2
      if (izai.eq.1) then
         if (mf.eq.3) then
            if (mt.ne.1.and.mt.ne.2.and.mt.ne.301.and.&
                (mt.ne.5.or.(mt.eq.5.and.mt5n.eq.0))) then
               ntr=ntr+1
               if ((mt.ge.5.and.mt.le.91).or.&
                 (mt.ge.152.and.mt.le.154).or.&
                 (mt.ge.156.and.mt.le.181).or.&
                 (mt.ge.183.and.mt.le.190).or.&
                 (mt.ge.194.and.mt.le.196).or.&
                 (mt.ge.198.and.mt.le.200).or.&
                 (mt.ge.875.and.mt.le.899)) then
                  nr=nr+1
                  if (mt16.gt.0.and.mt.eq.16) nr=nr-1
                  mtnr=mt
               endif
            endif
            mtntr=mt
         endif
      else if (izai.eq.1001) then
         if (mf.eq.3) then
            if (mt.ne.1.and.mt.ne.2.and.&
                (mt.ne.5.or.(mt.eq.5.and.mt5p.eq.0))) then
               ntr=ntr+1
               if (mt.eq.5.or.mt.eq.28.or.mt.eq.41.or.&
                 mt.eq.42.or.mt.eq.44.or.mt.eq.45.or.&
                 mt.eq.103.or.mt.eq.111.or.&
                 mt.eq.112.or.mt.eq.115.or.mt.eq.116.or.&
                 mt.eq.156.or.mt.eq.159.or.mt.eq.162.or.&
                 mt.eq.163.or.mt.eq.164.or.mt.eq.179.or.&
                 mt.eq.181.or.mt.eq.183.or.mt.eq.184.or.&
                 mt.eq.186.or.mt.eq.190.or.mt.eq.191.or.&
                 mt.eq.194.or.mt.eq.196.or.mt.eq.197.or.&
                 mt.eq.198.or.mt.eq.199.or.mt.eq.200.or.&
                 (mt.ge.600.and.mt.le.649)) then
                  nr=nr+1
                  mtnr=mt
               endif
            endif
            mtntr=mt
         endif
      else if (izai.eq.1002) then
         if (mf.eq.3) then
            if (mt.ne.1.and.mt.ne.2.and.&
                (mt.ne.5.or.(mt.eq.5.and.mt5d.eq.0))) then
               ntr=ntr+1
               if (mt.eq.5.or.mt.eq.11.or.mt.eq.32.or.mt.eq.35.or.&
                 mt.eq.104.or.mt.eq.114.or.mt.eq.115.or.mt.eq.117.or.&
                 mt.eq.157.or.mt.eq.158.or.mt.eq.169.or.&
                 mt.eq.170.or.mt.eq.171.or.mt.eq.182.or.&
                 mt.eq.183.or.mt.eq.185.or.mt.eq.187.or.&
                 mt.eq.192.or.&
                 (mt.ge.650.and.mt.le.699)) then
                  nr=nr+1
                  mtnr=mt
               endif
            endif
            mtntr=mt
         endif
      else if (izai.eq.1003) then
         if (mf.eq.3) then
            if (mt.ne.1.and.mt.ne.2.and.&
                (mt.ne.5.or.(mt.eq.5.and.mt5t.eq.0))) then
               ntr=ntr+1
               if (mt.eq.5.or.mt.eq.33.or.mt.eq.36.or.&
                 mt.eq.105.or.mt.eq.113.or.mt.eq.116.or.&
                 mt.eq.154.or.mt.eq.155.or.mt.eq.172.or.&
                 mt.eq.173.or.mt.eq.174.or.mt.eq.175.or.&
                 mt.eq.182.or.mt.eq.184.or.mt.eq.185.or.&
                 mt.eq.188.or.mt.eq.189.or.&
                 (mt.ge.700.and.mt.le.749)) then
                  nr=nr+1
                  mtnr=mt
               endif
            endif
            mtntr=mt
         endif
      else if (izai.eq.2003) then
         if (mf.eq.3) then
            if (mt.ne.1.and.mt.ne.2.and.&
                (mt.ne.5.or.(mt.eq.5.and.mt5he3.eq.0))) then
               ntr=ntr+1
               if (mt.eq.5.or.mt.eq.34.or.mt.eq.106.or.&
                  mt.eq.176.or.mt.eq.177.or.mt.eq.178.or.&
                  mt.eq.186.or.mt.eq.187.or.mt.eq.188.or.&
                  mt.eq.191.or.mt.eq.192.or.mt.eq.193.or.&
                  (mt.ge.750.and.mt.le.799)) then
                  nr=nr+1
                  mtnr=mt
               endif
            endif
            mtntr=mt
         endif
      else if (izai.eq.2004) then
         if (mf.eq.3) then
            if (mt.ne.1.and.mt.ne.2.and.&
                (mt.ne.5.or.(mt.eq.5.and.mt5a.eq.0))) then
               ntr=ntr+1
               if (mt.eq.5.or.(mt.ge.22.and.mt.le.25).or.&
                 mt.eq.29.or.mt.eq.30.or.&
                 mt.eq.35.or.mt.eq.36.and.mt.eq.45.or.&
                 mt.eq.107.or.mt.eq.108.or.mt.eq.109.or.&
                 mt.eq.112.or.mt.eq.113.or.&
                 mt.eq.114.or.mt.eq.117.or.&
                 mt.eq.155.or.mt.eq.158.or.mt.eq.159.or.&
                 mt.eq.165.or.mt.eq.166.or.mt.eq.167.or.&
                 mt.eq.168.or.mt.eq.180.or.mt.eq.181.or.&
                 mt.eq.189.or.mt.eq.193.or.mt.eq.195.or.&
                  mt.eq.196.or.mt.eq.199.or.&
                 (mt.ge.800.and.mt.le.849)) then
                  nr=nr+1
                  mtnr=mt
               endif
            endif
            mtntr=mt
         endif
      endif
   enddo
   nurd=0
   iurpt=0
   urhi=0

   !--process file 1 fission nu data
   if (kfis.ne.0) then
      idone=0
      mta=452
      do while (idone.eq.0)
         call findf(matd,1,mta,nin)
         call contio(nin,0,0,scr,nb,nw)
         lnu=nint(scr(4))
         if (lnu.eq.1.and.(mta.eq.452.or.mta.eq.456)) then
            ! polynomial representation
            call listio(nin,0,0,scr,nb,nw)
            n=nint(scr(5))
            nw=n+2
            if (mta.eq.452) then
               allocate(nut(nw))
               nut(1)=lnu
               nut(2)=n
               do i=1,n
                  nut(i+2)=scr(i+6)
                  if (i.gt.1) nut(i+2)=nut(i+2)*emev*(i-1)
               enddo
               nnu=nw
            endif
            if (mta.eq.456) then
               allocate(nup(nw))
               nup(1)=lnu
               nup(2)=n
               do i=1,n
                  nup(i+2)=scr(i+6)
                  if (i.gt.1) nup(i+2)=nup(i+2)*emev*(i-1)
               enddo
               nnup=nw
            endif
         elseif (lnu.eq.1.and.mta.eq.455) then
            ! polynomial expansion for spontaneous fission
            call listio(nin,0,0,scr,nb,nw)
            nnf=n1h
            do i=1,nnf
               dntc(i)=scr(6+i)
            enddo
            call listio(nin,0,0,scr,nb,nw)
            n=n1h
            nw=3+2*(n+1)
            allocate(nudn(nw))
            nudn(1)=2
            nudn(2)=0
            nudn(3)=n
            nudn(4)=1.e-11_kr
            nudn(5)=scr(7)
            nudn(6)=etop
            nudn(7)=scr(7)
            nnud=nw
         else
            ! tabular representation
            if (mta.eq.455) then
               call listio(nin,0,0,scr,nb,nw)
               nnf=n1h
               do i=1,nnf
                  dntc(i)=scr(6+i)
               enddo
            endif
            call tab1io(nin,0,0,scr,nb,nw)
            inow=1+nw
            do while (nb.ne.0)
               call moreio(nin,0,0,scr(inow),nb,nw)
               inow=inow+nw
            enddo
            m=nint(scr(5))
            n=nint(scr(6))
            jnt=nint(scr(8))
            nw=3+2*n
            if (m.ne.1.or.jnt.ne.2) then
               nw=nw+2*m
            endif
            if (mta.eq.452) then
               allocate(nut(nw))
               nut(1)=lnu
               j=2
               nut(j)=0
               if (m.ne.1.or.jnt.ne.2) then
                  nut(j)=m
                  do i=1,m
                     nut(i+j)=scr(2*i+5)
                     nut(i+m+j)=scr(2*i+6)
                  enddo
                  j=j+2*m
               endif
               j=j+1
               nut(j)=n
               do i=1,n
                  nut(i+j)=sigfig(scr(2*i+2*m+5)/emev,7,0)
                  nut(i+j+n)=scr(2*i+2*m+6)
               enddo
               nnu=nw
            endif
            if (mta.eq.455) then
               allocate(nudn(nw))
               nudn(1)=lnu
               j=2
               nudn(j)=0
               if (m.ne.1.or.jnt.ne.2) then
                  nudn(j)=m
                  do i=1,m
                     nudn(i+j)=scr(2*i+5)
                     nudn(i+m+j)=scr(2*i+6)
                  enddo
                  j=j+2*m
               endif
               j=j+1
               nudn(j)=n
               do i=1,n
                  nudn(i+j)=sigfig(scr(2*i+2*m+5)/emev,7,0)
                  nudn(i+j+n)=scr(2*i+2*m+6)
               enddo
               nnud=nw
            endif
            if (mta.eq.456) then
               allocate(nup(nw))
               nup(1)=lnu
               j=2
               nup(j)=0
               if (m.ne.1.or.jnt.ne.2) then
                  nup(j)=m
                  do i=1,m
                     nup(i+j)=scr(2*i+5)
                     nup(i+m+j)=scr(2*i+6)
                  enddo
                  j=j+2*m
               endif
               j=j+1
               nup(j)=n
               do i=1,n
                  nup(i+j)=sigfig(scr(2*i+2*m+5)/emev,7,0)
                  nup(i+j+n)=scr(2*i+2*m+6)
               enddo
               nnup=nw
            endif
         endif
         if (mta.eq.452.and.kfis.eq.2) then
            mta=455
         else if (mta.eq.455) then
            mta=456
         else
            idone=1
         endif
      enddo
      call tofend(nin,0,0,scr)
   endif

   !--process unresolved probability tables if found
   call findf(matd,2,0,nin)
   mfh=2
   do while (mfh.ne.0)
      call contio(nin,0,0,scr,nb,nw)
      if (mth.eq.153) then
         write(nsyso,'(/'' found mt=153 with unresolved-range'',&
           &'' probability tables'')')
         iinel=l1h !--new competition flag (can only be -1, 51, 91 or 4)
         iabso=l2h !--new competition flag (can be -1, 0 or positive)
         call listio(nin,0,0,scr,nb,nw)
         lssf=l1h
         if (iinel.eq.zero) then !--old competition flags are used
            iinel=-1
            iabso=-1
            if (mod(l2h,1000).ne.zero) iinel=mod(l2h,1000)
            if (l2h/1000.ne.zero) iabso=l2h/1000
         endif
         nunr=n2h
         ncyc=n1h/nunr
         nurd=n1h+6
         allocate(urd(nurd))
         i1=1
         do i=1,nw
            urd(i1-1+i)=scr(i)
         enddo
         do while (nb.ne.0)
            i1=i1+nw
            call moreio(nin,0,0,scr,nb,nw)
            do i=1,nw
               urd(i1-1+i)=scr(i)
            enddo
         enddo
         urlo=urd(7)
         urhi=urd(7+ncyc*(nunr-1))
         write(nsyso,'(''   energy range: '',1p,e11.4,&
           &'' - '',e11.4,'' ev'')') urlo,urhi
         if (lssf.eq.0)&
           write(nsyso,'(''   tables are cross sections'')')
         if (lssf.eq.1)&
           write(nsyso,'(''   tables are factors'')')
         if (iinel.lt.0.and.iabso.lt.0) then
            write(nsyso,'(''   no competition'')')
         else
            write(nsyso,'(''   inelastic competition ='',i3/&
              &''   absorption competition ='',i3)') iinel,iabso
         endif
      endif
      if (mfh.ne.0) then
         call tosend(nin,0,0,scr)
      endif
   enddo

   !--locate and store energy grid of total xsec
   esz=0
   call findf(matd,3,1,nin)
   call contio(nin,0,0,scr,nb,nw)
   e=0
   call gety1(e,enext,idis,s,nin,scr)
   k=0
   test=etop-etop/100
   do while (enext.lt.test)
      e=enext
      call gety1(e,enext,idis,s,nin,scr)
      k=k+1
      xss(esz+k)=e
   enddo
   nes=k
   if ((esz+5*nes).gt.nxss)&
     call error('acelod','insufficient storage for esz block.',' ')
   call tosend(nin,0,0,scr)

   !--assign locators for esz block
   it=esz+nes
   ic=it+nes
   ie=ic+nes
   ih=ie+nes
   next=ih+nes
   next=next+1

   !--assign locators for nu thru sig blocks
   nu=0
   if (kfis.ne.0) then
      nu=next
      next=nu+nnu
      if (kfis.ne.1) then
         next=next+1+nnup
      endif
   endif
   mtr=next
   lqr=mtr+ntr
   tyr=lqr+ntr
   lsig=tyr+ntr
   sig=lsig+ntr
   keep3=0
   keep4=0
   keep=0
   next=sig
   ir=0

   !--read and store cross sections producing incident particle
   mt=-1
   do while (mt.lt.mtnr)
      call contio(nin,0,0,scr,nb,nw)
      mt=mth
      if (mth.eq.3) keep3=1
      if (mth.eq.4.and.izai.eq.1) keep4=1
      if (izai.eq.1) then
         iskip=0
         if (mt.eq.3.or.mt.eq.4) iskip=1
         if (mt.eq.5.and.mt5n.eq.1) iskip=1
         if (mt.gt.91.and.mt.le.151) iskip=1
         if (mt.eq.155.or.mt.eq.182.or.mt.eq.191) iskip=1
         if (mt.eq.192.or.mt.eq.193.or.mt.eq.197) iskip=1
         if (mt.gt.200.and.mt.le.849) iskip=1
         if (mt16.gt.0.and.mt.eq.16) iskip=1
      else if (izai.eq.1001) then
         iskip=1
         if (mt.eq.2.or.mt.eq.28.or.mt.eq.41.or.&
           mt.eq.42.or.mt.eq.44.or.mt.eq.45.or.&
           mt.eq.103.or.mt.eq.111.or.&
           mt.eq.112.or.mt.eq.115.or.mt.eq.116.or.&
           mt.eq.156.or.mt.eq.159.or.mt.eq.162.or.&
           mt.eq.163.or.mt.eq.164.or.mt.eq.179.or.&
           mt.eq.181.or.mt.eq.183.or.mt.eq.184.or.&
           mt.eq.186.or.mt.eq.190.or.mt.eq.191.or.&
           mt.eq.194.or.mt.eq.196.or.mt.eq.197.or.&
           mt.eq.198.or.mt.eq.199.or.mt.eq.200.or.&
           (mt.ge.600.and.mt.le.649)) iskip=0
         if (mt.eq.5.and.mt5p.eq.0) iskip=0
      else if (izai.eq.1002) then
         iskip=1
         if (mt.eq.2.or.mt.eq.11.or.mt.eq.32.or.mt.eq.35.or.&
           mt.eq.104.or.mt.eq.114.or.mt.eq.115.or.mt.eq.117.or.&
           mt.eq.157.or.mt.eq.158.or.mt.eq.169.or.&
           mt.eq.170.or.mt.eq.171.or.mt.eq.182.or.&
           mt.eq.183.or.mt.eq.185.or.mt.eq.187.or.&
           mt.eq.192.or.&
           (mt.ge.650.and.mt.le.699)) iskip=0
         if (mt.eq.5.and.mt5d.eq.0) iskip=0
      else if (izai.eq.1003) then
         iskip=1
         if (mt.eq.2.or.mt.eq.33.or.mt.eq.36.or.&
           mt.eq.105.or.mt.eq.113.or.mt.eq.116.or.&
           mt.eq.154.or.mt.eq.155.or.mt.eq.172.or.&
           mt.eq.173.or.mt.eq.174.or.mt.eq.175.or.&
           mt.eq.182.or.mt.eq.184.or.mt.eq.185.or.&
           mt.eq.188.or.mt.eq.189.or.&
           (mt.ge.700.and.mt.le.749)) iskip=0
         if (mt.eq.5.and.mt5t.eq.0) iskip=0
      else if (izai.eq.2003) then
         iskip=1
         if (mt.eq.2.or.mt.eq.34.or.mt.eq.106.or.&
           mt.eq.176.or.mt.eq.177.or.mt.eq.178.or.&
           mt.eq.186.or.mt.eq.187.or.mt.eq.188.or.&
           mt.eq.191.or.mt.eq.192.or.mt.eq.193.or.&
           (mt.ge.750.and.mt.le.799)) iskip=0
         if (mt.eq.5.and.mt5he3.eq.0) iskip=0
      else if (izai.eq.2004) then
         iskip=1
         if (mt.eq.2.or.(mt.ge.22.and.mt.le.25).or.&
           mt.eq.29.or.mt.eq.30.or.&
           mt.eq.35.or.mt.eq.36.and.mt.eq.45.or.&
           mt.eq.107.or.mt.eq.108.or.mt.eq.109.or.&
           mt.eq.112.or.mt.eq.113.or.&
           mt.eq.114.or.mt.eq.117.or.&
           mt.eq.155.or.mt.eq.158.or.mt.eq.159.or.&
           mt.eq.165.or.mt.eq.166.or.mt.eq.167.or.&
           mt.eq.168.or.mt.eq.180.or.mt.eq.181.or.&
           mt.eq.189.or.mt.eq.193.or.mt.eq.195.or.&
           mt.eq.196.or.mt.eq.199.or.&
           (mt.ge.800.and.mt.le.849)) iskip=0
         if (mt.eq.5.and.mt5a.eq.0) iskip=0
      endif
      if (iskip.eq.0) then
         e=0
         call gety1(e,enext,idis,s,nin,scr)
         j=1
         n=-1

         !--check for reaction cross sections
         if (mth.ne.2) then

            !--locate energy index for threshold
            j=nes
            idone=0
            do while (idone.eq.0)
               if (j.ge.1) then
                  if (enext.le.(1+eps)*xss(esz+j)) then
                     j=j-1
                  else
                     idone=1
                  endif
               else
                  idone=1
               endif
            enddo

            !--store reaction parameters
            j=j+1
            xss(mtr+ir)=mth
            xss(lqr+ir)=sigfig(scr(2)/emev,7,0)
            xss(tyr+ir)=0
            xss(lsig+ir)=next-sig+1
            ir=ir+1
            if (mth.eq.18) fis=next
            xss(next)=j
            nnex=next+1
            next=next+2
            n=0
         endif

         !--store cross sections
         do while (j.le.nes)
            e=xss(esz+j)
            call gety1(e,enext,idis,s,nin,scr)
            test=100
            if (j.gt.1.and.n.eq.0.and.enext.gt.test) s=0
            s=sigfig(s,7,0)
            if (mth.eq.2) then
               xss(ie+j)=s
               xss(it+j)=xss(it+j)+s
            else
               if (mth.ge.6.and.mth.le.9) s=sigfig(s/2,7,0)
               if (mth.ge.46.and.mth.le.49) s=sigfig(s/2,7,0)
               xss(next)=s
               if (mth.ge.5) xss(it+j)=xss(it+j)+s
               n=n+1
               next=next+1
               if (next.gt.nxss) call error('acelod',&
                 'insufficient space for cross sections',' ')
            endif
            j=j+1
         enddo
         if (mth.ne.2) xss(nnex)=n
      endif
      call tosend(nin,0,0,scr)
   enddo

   !--read and store cross sections not producing incident particle
   call findf(matd,3,2,nin)
   mt=2
   do while (mt.lt.mtntr)
      call contio(nin,0,0,scr,nb,nw)
      mt=mth
      if (izai.eq.1) then
         iskip=1
         if (mt.gt.91.and.mt.lt.152) iskip=0
         if (mt.eq.155.or.mt.eq.182.or.mt.eq.191) iskip=0
         if (mt.eq.192.or.mt.eq.193.or.mt.eq.197) iskip=0
         if (mt.gt.200.and.mt.le.849) iskip=0
         if (mt16.gt.0.and.mt.eq.16) iskip=0
      else if (izai.eq.1001) then
         iskip=0
         if (mt.eq.2.or.mt.eq.5.or.mt.eq.28.or.mt.eq.41.or.&
           mt.eq.42.or.mt.eq.44.or.mt.eq.45.or.&
           mt.eq.103.or.mt.eq.111.or.&
           mt.eq.112.or.mt.eq.115.or.mt.eq.116.or.&
           mt.eq.156.or.mt.eq.159.or.mt.eq.162.or.&
           mt.eq.163.or.mt.eq.164.or.mt.eq.179.or.&
           mt.eq.181.or.mt.eq.183.or.mt.eq.184.or.&
           mt.eq.186.or.mt.eq.190.or.mt.eq.191.or.&
           mt.eq.194.or.mt.eq.196.or.mt.eq.197.or.&
           mt.eq.198.or.mt.eq.199.or.mt.eq.200.or.&
           (mt.ge.600.and.mt.le.649)) iskip=1
      else if (izai.eq.1002) then
         iskip=0
         if (mt.eq.2.or.mt.eq.5.or.mt.eq.11.or.&
           mt.eq.32.or.mt.eq.35.or.&
           mt.eq.104.or.mt.eq.114.or.mt.eq.115.or.mt.eq.117.or.&
           mt.eq.157.or.mt.eq.158.or.mt.eq.169.or.&
           mt.eq.170.or.mt.eq.171.or.mt.eq.182.or.&
           mt.eq.183.or.mt.eq.185.or.mt.eq.187.or.&
           mt.eq.192.or.&
           (mt.ge.650.and.mt.le.699)) iskip=1
      else if (izai.eq.1003) then
         iskip=0
         if (mt.eq.2.or.mt.eq.5.or.mt.eq.33.or.mt.eq.36.or.&
           mt.eq.105.or.mt.eq.113.or.mt.eq.116.or.&
           mt.eq.154.or.mt.eq.155.or.mt.eq.172.or.&
           mt.eq.173.or.mt.eq.174.or.mt.eq.175.or.&
           mt.eq.182.or.mt.eq.184.or.mt.eq.185.or.&
           mt.eq.188.or.mt.eq.189.or.&
           (mt.ge.700.and.mt.le.749)) iskip=1
      else if (izai.eq.2003) then
         iskip=0
         if (mt.eq.2.or.mt.eq.5.or.mt.eq.34.or.mt.eq.106.or.&
           mt.eq.176.or.mt.eq.177.or.mt.eq.178.or.&
           mt.eq.186.or.mt.eq.187.or.mt.eq.188.or.&
           mt.eq.191.or.mt.eq.192.or.mt.eq.193.or.&
           (mt.ge.750.and.mt.le.799)) iskip=1
      else if (izai.eq.2004) then
         iskip=0
         if (mt.eq.2.or.mt.eq.5.or.(mt.ge.22.and.mt.le.25).or.&
           mt.eq.29.or.mt.eq.30.or.&
           mt.eq.35.or.mt.eq.36.and.mt.eq.45.or.&
           mt.eq.107.or.mt.eq.108.or.mt.eq.109.or.&
           mt.eq.112.or.mt.eq.113.or.&
           mt.eq.114.or.mt.eq.117.or.&
           mt.eq.155.or.mt.eq.158.or.mt.eq.159.or.&
           mt.eq.165.or.mt.eq.166.or.mt.eq.167.or.&
           mt.eq.168.or.mt.eq.180.or.mt.eq.181.or.&
           mt.eq.189.or.mt.eq.193.or.mt.eq.195.or.&
           mt.eq.196.or.mt.eq.199.or.&
           (mt.ge.800.and.mt.le.849)) iskip=1
      endif
      if (iskip.eq.0) then
         e=0
         call gety1(e,enext,idis,s,nin,scr)
         j=1
         n=-1

         !--check for reaction cross sections
         if (mth.ne.301) then

            !--locate energy index for threshold
            j=nes
            idone=0
            do while (idone.eq.0)
               if (j.ge.1) then
                  if (enext.le.(1+eps)*xss(esz+j)) then
                     j=j-1
                  else
                     idone=1
                  endif
               else
                  idone=1
               endif
            enddo

            !--store reaction parameters
            j=j+1
            xss(mtr+ir)=mth
            xss(lqr+ir)=sigfig(scr(2)/emev,7,0)
            xss(tyr+ir)=0
            xss(lsig+ir)=next-sig+1
            ir=ir+1
            if (mth.eq.18) fis=next
            xss(next)=j
            nnex=next+1
            next=next+2
            n=0
         endif

         !--store cross sections
         do while (j.le.nes)
            e=xss(esz+j)
            call gety1(e,enext,idis,s,nin,scr)
            test=100
            if (j.gt.1.and.n.eq.0.and.enext.gt.test) s=0
            s=sigfig(s,7,0)
            if (mth.eq.301) then
               if (xss(it+j).eq.zero) then
                  xss(ih+j)=0
               else
                  xss(ih+j)=sigfig(s/emev/xss(it+j),7,0)
               endif
            else
               if (izai.eq.1) then
                  if ((mth.ge.102.and.mth.le.150)&
                    .or.mth.eq.155.or.mth.eq.182.or.mth.eq.191&
                    .or.mth.eq.192.or.mth.eq.193.or.mth.eq.197&
                    .or.(mt103.eq.0.and.mth.ge.mpmin.and.mth.le.mpmax)&
                    .or.(mt104.eq.0.and.mth.ge.mdmin.and.mth.le.mdmax)&
                    .or.(mt105.eq.0.and.mth.ge.mtmin.and.mth.le.mtmax)&
                    .or.(mt106.eq.0.and.mth.ge.m3min.and.mth.le.m3max)&
                    .or.(mt107.eq.0.and.mth.ge.m4min.and.mth.le.m4max)) then
                     xss(ic+j)=xss(ic+j)+s
                  endif
               else
                  if (mth.le.200.or.mth.ge.600) then
                     xss(ic+j)=xss(ic+j)+s
                  endif
               endif
               if (mth.eq.444) s=sigfig(s/emev,7,0)
               xss(next)=s
               if (izai.eq.1) then
                  if ((mth.ge.5.and.mth.le.150)&
                    .or.(mth.ge.152.and.mth.le.200)&
                    .or.(mt103.eq.0.and.mth.ge.mpmin.and.mth.le.mpmax)&
                    .or.(mt104.eq.0.and.mth.ge.mdmin.and.mth.le.mdmax)&
                    .or.(mt105.eq.0.and.mth.ge.mtmin.and.mth.le.mtmax)&
                    .or.(mt106.eq.0.and.mth.ge.m3min.and.mth.le.m3max)&
                    .or.(mt107.eq.0.and.mth.ge.m4min.and.mth.le.m4max)) then
                     xss(it+j)=xss(it+j)+s
                  endif
               else
                  if (mth.le.200.or.mth.ge.600) then
                     xss(it+j)=xss(it+j)+s
                  endif
               endif
               n=n+1
               next=next+1
               if (next.gt.nxss) call error('acelod',&
                 'insufficient space for cross sections',' ')
            endif
            j=j+1
         enddo
         if (mth.ne.2.and.mth.ne.301) xss(nnex)=n
      endif
      call tosend(nin,0,0,scr)
   enddo

   !--go back and add mt3 and/or 4, if needed
   if (keep3.eq.1.or.keep4.eq.1) then
      keep1=3
      if (keep3.eq.0) keep1=4
      keep2=4
      if (keep4.eq.0) keep2=3
      keep=keep1
      do while (keep.le.keep2)
         call findf(matd,3,keep,nin)
         call contio(nin,0,0,scr,nb,nw)
         e=0
         call gety1(e,enext,idis,s,nin,scr)
         j=1
         n=-1

         !--locate energy index for threshold
         j=nes
         idone=0
         do while (idone.eq.0)
            if (j.ge.1) then
               if (enext.le.(1+eps)*xss(esz+j)) then
                  j=j-1
               else
                  idone=1
               endif
            else
               idone=1
            endif
         enddo

         !--store reaction parameters
         j=j+1
         xss(mtr+ir)=mth
         xss(lqr+ir)=sigfig(scr(2)/emev,7,0)
         xss(tyr+ir)=0
         xss(lsig+ir)=next-sig+1
         ir=ir+1
         if (mth.eq.18) fis=next
         xss(next)=j
         nnex=next+1
         next=next+2
         n=0

         !--store cross sections
         do while (j.le.nes)
            e=xss(esz+j)
            call gety1(e,enext,idis,s,nin,scr)
            test=100
            if (j.gt.1.and.n.eq.0.and.enext.gt.test) s=0
            s=sigfig(s,7,0)
            xss(next)=s
            n=n+1
            next=next+1
            if (next.gt.nxss) call error('acelod',&
              'insufficient space for cross sections',' ')
            j=j+1
         enddo
         xss(nnex)=n
         keep=keep+1
      enddo
   endif

   !--store fission nu data
   if (nu.ne.0) then
      l=nu-1
      if (nnup.ne.0) then
         l=l+1
         xss(l)=-nnup
         do i=1,nnup
            l=l+1
            xss(l)=nup(i)
         enddo
      endif
      do i=1,nnu
         l=l+1
         xss(l)=nut(i)
      enddo
   endif

   !--loop over angular distributions
   land=next
   and=land+1+nr
   next=and
   do ij=1,nr+1
      mt418=0
      ir=ij-1
      if (ij.eq.1) then
         mt=2
      else
         mt=nint(xss(mtr+ij-2))
      endif
      mf=0
      iso=0
      do k=1,nxc
         if (mts(k).eq.mt) then
            if (mt.ne.18.and.(mfs(k).eq.4.or.mfs(k).eq.6)) then
               mf=mfs(k)
            else
               if (mt.eq.18.and.mfs(k).eq.4.and.mt418.eq.0) then
                  mt418=1
                  mf=mfs(k)
                  exit
               else if (mt.eq.18.and.mfs(k).eq.6) then
                  mf=mfs(k)
                  exit
               endif
            endif
         endif
      enddo
      if (mf.eq.4.or.mf.eq.6) then
         call findf(matd,mf,mt,nin)
         call contio(nin,0,0,scr,nb,nw)
         awr=c2h
         iso=0
         if (mfh.ne.6) then
            law=-999
            lvt=nint(scr(3))
            ltt=nint(scr(4))
            ltt3=ltt
            lttn=1
            if (lvt.eq.0) then
               call contio(nin,0,0,scr,nb,nw)
            else
               jscr=1
               call listio(nin,0,0,scr(jscr),nb,nw)
               jscr=jscr+nw
               do while (nb.ne.0)
                  call moreio(nin,0,0,scr(jscr),nb,nw)
                  jscr=jscr+nw
               enddo
            endif
            iso=nint(scr(3))
            lct=nint(scr(4))
            if (iso.ne.1) then
               call tab2io(nin,0,0,scr,nb,nw)
               ne=nint(scr(6))
            endif
         else
            lct=l2h
            call tab1io(nin,0,0,scr,nb,nw)
            awp=c2h
            ltt=1
            law=l2h
            call tab2io(nin,0,0,scr,nb,nw)
            ne=nint(scr(6))
            spi=scr(1)
            lidp=nint(scr(3))
         endif

         !--coupled energy-angle distributions
         if (mfh.eq.6.and.(law.eq.1.or.law.eq.6)) then
            xss(land+ir)=-1

         !--regular angular distributions
         else
            if (mfh.eq.6.and.(law.eq.2.or.law.eq.3)) then
               do k=1,nxc
                  if (mts(k).eq.mth.and.mfs(k).eq.6) mfs(k)=64
               enddo
               if (law.eq.3) iso=1
            endif
            last=next
            xss(land+ir)=next-and+1
            n=1
            if (mth.ge.6.and.mth.le.9) n=2
            if (mth.eq.16.or.mth.eq.24.or.mth.eq.26.or.mth.eq.0) n=2
            if (mth.ge.46.and.mth.le.49) n=2
            if (mth.eq.11.or.mth.eq.41) n=2
            if (mth.eq.17.or.mth.eq.25) n=3
            if (mth.eq.42) n=3
            if (mth.eq.37) n=4
            if (mth.eq.18) n=19
            if ((mth.ge.19.and.mth.le.21).or.mth.eq.38) n=19
            if (lct.ge.2.and.mth.ne.18) n=-n
            if (mth.ne.2.and.mth.le.91) xss(tyr+ir-1)=n

            !--shorten table for all-isotropic reactions.
            if (iso.eq.1) then
               xss(land+ir)=0
               next=last

            !--treat anisotropic reactions.
            else
               xss(next)=ne
               ie=next
               il=ie+ne
               next=il+ne+1
               iso=1

               !--store neutron angular distributions
               if (izai.eq.1) then
                  call acensd(ir,next,scr,nin,awp,ltt3,lttn,&
                    ltt,last,law,ne,ie,il,iso,newfor)

               !--treat charged-particle elastic
               else
                  call acecpe(next,scr,nin,awr,awp,&
                    spi,ne,lidp,ie,il,nes)
               endif
            endif
         endif
      endif
   enddo
   call skiprz(nin,-1)

   !--store energy distribution data
   ldlw=next
   next=ldlw+nr
   dlw=next
   if (nr.ne.0) then
      do i=1,nr
         mt518=0
         mt=nint(xss(mtr+i-1))
         q=xss(lqr+i-1)
         do k=1,nxc
            if (mts(k).eq.mt) then
               if (mt.eq.18.and.mfs(k).eq.5.and.mt518.eq.0) then
                  mt518=1
                  mf=mfs(k)
                  exit
               else if (mt.eq.18.and.mfs(k).eq.6) then
                  mf=mfs(k)
                  exit
               else
                  mf=mfs(k)
               endif
            endif
         enddo
         if (mf.eq.5) then
            call acelf5(next,i,matd,mt,q,nin,ismooth)
         else if (mf.eq.6) then
            if (mt518.eq.0) then
               call acelf6(next,i,matd,mt,q,iza,izai,nin,newfor,ismooth)
            endif
         else
            if ((next+11).gt.nxss) call error('acelod',&
              'insufficient space for energy distributions',' ')
            xss(ldlw+i-1)=next-dlw+1
            xss(next)=0
            xss(next+1)=3
            xss(next+2)=next-dlw+10
            xss(next+3)=0
            xss(next+4)=2
            l=nint(xss(lsig+i-1))
            ja=nint(xss(sig+l-1))
            jb=nint(xss(sig+l)+ja-1)
            xss(next+5)=sigfig(xss(esz+ja)/emev,7,0)
            xss(next+6)=sigfig(xss(esz+jb)/emev,7,0)
            xss(next+7)=1
            xss(next+8)=1
            x=(awr+1)/awr
            xss(next+9)=sigfig(x*(-q),7,0)
            xss(next+10)=sigfig(1/x**2,7,0)
            if (xss(next+9).gt.xss(next+5)) then
               xss(next+9)=sigfig(xss(next+5)*0.999998_kr,7,0)
            endif
            next=next+11
         endif
      enddo
   endif
   call tofend(nin,0,0,scr)

   !--store unresolved-range probability tables
   !--after energy distributions
   if (nurd.ne.0) then
      iurpt=next
      nure=nint(urd(6))
      xss(next)=nure
      nurb=nint(urd(5))
      nurb=(nurb/nure-1)/6
      xss(next+1)=nurb
      xss(next+2)=2
      xss(next+3)=iinel
      xss(next+4)=iabso
      xss(next+5)=lssf
      next=next+6
      do ie=1,nure
         jj=7+(ie-1)*(1+6*nurb)
         xss(next-1+ie)=sigfig(urd(jj)/emev,7,0)
         ll=next-1+nure+(ie-1)*6*nurb
         do ib=1,nurb
            if (ib.eq.1) xss(ll+ib)=urd(jj+ib)
            if (ib.gt.1) xss(ll+ib)=xss(ll+ib-1)+urd(jj+ib)
            xss(ll+nurb+ib)=urd(jj+nurb+ib)
            xss(ll+2*nurb+ib)=urd(jj+2*nurb+ib)
            xss(ll+3*nurb+ib)=urd(jj+3*nurb+ib)
            xss(ll+4*nurb+ib)=urd(jj+4*nurb+ib)
            if (lssf.eq.0) then
               xss(ll+5*nurb+ib)=urd(jj+5*nurb+ib)/emev
            else
               xss(ll+5*nurb+ib)=urd(jj+5*nurb+ib)
            endif
         enddo
         do ib=1,nurb
            xss(ll+ib)=sigfig(xss(ll+ib),7,0)
         enddo
         xss(ll+nurb)=1
      enddo
      next=next+nure*(1+6*nurb)
   endif

   !--store delayed neutron data
   ndnf=0
   if (nnud.gt.0.and.mt455.eq.0) write(nsyso,'(/&
     &''  a delayed nubar section was found, but''/&
     &''   no delayed neutron spectra were found:''/&
     &''   delayed neutron data supressed'')')
   if (nnud.gt.0.and.mt455.eq.1) then
      write(nsyso,'(/''  adding delayed neutron data'')')
      nud=next
      l=next-1

      !--fission delayed nubar data
      do i=1,nnud
         l=l+1
         xss(l)=nudn(i)
      enddo
      next=l+1

      !--locate the delayed neutron data in file 5
      call findf(matd,5,455,nin)
      call contio(nin,0,0,scr,nb,nw)
      ndnf=n1h

      !--dndat block
      dndat=next

      !--read through the section to load the dndat block
      lff=dndat
      do i=1,ndnf
         call tab1io(nin,0,0,scr,nb,nw)
         law=l2h
         nn=n1h
         n=n2h
         !--dndat entry
         xss(lff)=dntc(i)/shake
         lff=lff+1
         if (nn.eq.1.and.nint(scr(8)).eq.2) then
            xss(lff)=0
            lff=lff+1
         else
            xss(lff)=nn
            do j=1,nn
               xss(lff+j)=scr(5+2*j)
               xss(lff+nn+j)=scr(6+2*j)
            enddo
            lff=lff+1+2*nn
         endif
         xss(lff)=n
         do j=1,n
            xss(lff+j)=sigfig(scr(5+2*nn+2*j)/emev,7,0)
            xss(lff+n+j)=sigfig(scr(6+2*nn+2*j),7,0)
         enddo
         lff=lff+1+2*n
         if (law.eq.1) then
            !--law=1
            call tab2io(nin,0,0,scr,nb,nw)
            ne=n2h
            do ie=1,ne
               call tab1io(nin,0,0,scr,nb,nw)
               do while (nb.ne.0)
                  call moreio(nin,0,0,scr,nb,nw)
               enddo
            enddo
         else if (law.eq.5) then
            !--law=5
            call tab1io(nin,0,0,scr,nb,nw)
            call tab1io(nin,0,0,scr,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,0,0,scr,nb,nw)
            enddo
         endif
      enddo
      next=lff

      !--ldnd block
      ldnd=next
      next=ldnd+ndnf

      !--dnd block
      dnd=next

      !--go back to the start of the sections
      call skiprz(nin,-2)
      call findf(matd,5,455,nin)
      call contio(nin,0,0,scr,nb,nw)

      !--store the data
      do i=1,ndnf
         call tab1io(nin,0,0,scr,nb,nw)
         law=l2h
         nn=n1h
         n=n2h
         xxmin=scr(7+2*nn)
         xxmax=scr(5+2*nn+2*n)
         !--ldnd entry
         xss(ldnd-1+i)=next-dnd+1
         !--dnd data
         !--there is only one law per family
         xss(next)=0
         xss(next+1)=4
         xss(next+2)=next-dnd+10
         xss(next+3)=0
         xss(next+4)=2
         xss(next+5)=sigfig(xxmin/emev,7,0)
         xss(next+6)=sigfig(xxmax/emev,7,0)
         xss(next+7)=1
         xss(next+8)=1
         if (law.eq.1) then
            !--law=1
            call tab2io(nin,0,0,scr,nb,nw)
            ne=n2h
            xss(next+9)=0
            xss(next+10)=ne
            next=next+11
            lxx=next
            next=next+2*ne
            do ie=1,ne
               call tab1io(nin,0,0,scr,nb,nw)
               nn=n1h
               mm=n2h
               iint=nint(scr(8))
               xss(lxx+ie-1)=sigfig(c2h/emev,7,0)
               xss(lxx+ne+ie-1)=next-dnd+1
               xss(next)=iint
               xss(next+1)=mm
               loc=1+nw
               do while (nb.ne.0)
                  call moreio(nin,0,0,scr(loc),nb,nw)
                  loc=loc+nw
               enddo
               l=next+1
               sumup=0
               do j=1,mm
                  xss(l+j)=sigfig(scr(5+2*nn+2*j)/emev,7,0)
                  xss(l+mm+j)=sigfig(scr(5+2*nn+2*j+1)*emev,7,0)
                  xss(l+2*mm+j)=sigfig(sumup,9,0)
                  ll=5+2*nn+2*j
                  if (j.lt.mm.and.iint.eq.1) then
                     sumup=sumup+(scr(ll+2)-scr(ll))*scr(ll+1)
                  else if (j.lt.mm.and.iint.eq.2) then
                     sumup=sumup+(scr(ll+2)-scr(ll))*(scr(ll+3)+scr(ll+1))/2
                  endif
               enddo
               if (100000000*abs(sumup-1).gt.1) then
                  write(nsyso,'(&
                    &'' renormalizing delayed spectrum:'',&
                    &'' precursor'',i2,&
                    &'' norm='',f11.8)') i,sumup
                  do j=1,mm
                     xss(l+mm+j)=sigfig(xss(l+mm+j)/sumup,7,0)
                     xss(l+2*mm+j)=sigfig(xss(l+2*mm+j)/sumup,9,0)
                  enddo
               endif
               next=next+2+3*mm
            enddo
         else if (law.eq.5) then
            !--law=5
            call tab1io(nin,0,0,scr,nb,nw)
            call tab1io(nin,0,0,scr,nb,nw)
            loc=1+nw
            do while (nb.ne.0)
               call moreio(nin,0,0,scr(loc),nb,nw)
               loc=loc+nw
            enddo
            nn=n1h
            mm=n2h

            ! extend lowest delayed bin using sqrt(e) shape
            if (ismooth.gt.0.and.nint(scr(8)).eq.1) then
               ex=40
               fx=.8409
               write(nsyso,&
                 '('' extending lowest delayed bin using sqrt(E)'')')
               do while (scr(11).gt.ex)
                  do ix=2*mm,1,-1
                     scr(10+ix)=scr(8+ix)
                  enddo
                  scr(11)=fx*scr(13)
                  val=scr(10)
                  scr(10)=sqrt(fx)*val
                  scr(12)=(1-fx*sqrt(fx))*val/(1-fx)
                  mm=mm+1
               enddo
            endif

            !--there is no incident energy dependence, we represent
            !--this by two energies with duplicated distributions
            xss(next+9)=0
            xss(next+10)=2
            xss(next+11)=sigfig(xxmin/emev,7,0)
            xss(next+12)=sigfig(xxmax/emev,7,0)
            xss(next+13)=next+15-dnd+1
            xss(next+14)=next+15+2+3*mm-dnd+1
            iint=nint(scr(8))
            xss(next+15)=iint
            xss(next+15+2+3*mm)=iint
            xss(next+16)=mm
            xss(next+16+2+3*mm)=mm
            l=next+16
            sumup=0
            do j=1,mm
               xss(l+j)=sigfig(scr(5+2*nn+2*j)/emev,7,0)
               xss(l+j+2+3*mm)=sigfig(scr(5+2*nn+2*j)/emev,7,0)
               xss(l+mm+j)=sigfig(scr(5+2*nn+2*j+1)*emev,7,0)
               xss(l+mm+j+2+3*mm)=sigfig(scr(5+2*nn+2*j+1)*emev,7,0)
               xss(l+2*mm+j)=sigfig(sumup,9,0)
               xss(l+2*mm+j+2+3*mm)=sigfig(sumup,9,0)
               ll=5+2*nn+2*j
               if (j.lt.mm.and.iint.eq.1) then
                  sumup=sumup+(scr(ll+2)-scr(ll))*scr(ll+1)
               else if (j.lt.mm.and.iint.eq.2) then
                  sumup=sumup+(scr(ll+2)-scr(ll))*(scr(ll+3)+scr(ll+1))/2
               endif
            enddo
            if (100000000*abs(sumup-1).gt.1) then
               write(nsyso,'(&
                 &'' renormalizing single delayed spectrum:'',&
                 &'' precursor'',i2,&
                 &'' norm='',f11.8)') i,sumup
               do j=1,mm
                  xss(l+mm+j)=sigfig(xss(l+mm+j)/sumup,7,0)
                  xss(l+mm+j+2+3*mm)=&
                    sigfig(xss(l+mm+j+2+3*mm)/sumup,7,0)
                  xss(l+2*mm+j)=sigfig(xss(l+2*mm+j)/sumup,9,0)
                  xss(l+2*mm+j+2+3*mm)=&
                   sigfig(xss(l+2*mm+j+2+3*mm)/sumup,9,0)
                  if (xss(l+2*mm+j).ge.0.999999997_kr.and.&
                         xss(l+2*mm+j).le.1.000100000_kr) then
                     xss(l+2*mm+j)=1
                     xss(l+2*mm+j+2+3*mm)=1
                  endif
               enddo
               if (xss(l+3*mm).lt.1) then
                  xss(l+3*mm)=1
                  xss(l+6*mm+2)=1
               endif
            endif
            next=next+15+2*(2+3*mm)
         endif
      enddo
   endif

   !--store photon production cross section
   if (mf1x(1).ne.0.or.mf1x(2).ne.0.or.n16.ne.0) then
      if (next+nes.gt.nxss) call error('acelod',&
        'insufficient storage for photon production',' ')
      call findf(matd,13,1,nin)
      call contio(nin,0,0,scr,nb,nw)
      e=0
      call gety1(e,enext,idis,s,nin,scr)
      gpd=next
      do i=1,nes
        e=xss(esz+i)
        teste=e/10000+eps
        if (i.eq.1.and.abs(e-enext).lt.teste) e=enext
        call gety1(e,enext,idis,s,nin,scr)
        xss(gpd+i-1)=s
      enddo
      next=gpd+nes

      !--store photon spectra
      !--omit obsolete 30x20 photon spectrum
      negn=0
      if (mf1x(3).ne.0) then
         if (next+600.gt.nxss) call error('acelod',&
           'insufficient storage for photon spectra.',' ')
         call findf(matd,15,0,nin)
         call contio(nin,0,0,scr,nb,nw)
         if (mth.eq.1) then
            call tab1io(nin,0,0,scr,nb,nw)
            call tab2io(nin,0,0,scr,nb,nw)
            negn=nint(scr(6))
            if (negn.ne.30) call error('acelod',&
              '30 groups are required for 30x20 photon matrix',' ')
            allocate(egn(negn))
            do i=1,negn
               call listio(nin,0,0,scr,nb,nw)
               egn(i)=sigfig(scr(2)/emev,7,0)
               do j=1,nbinp
                  xss(next)=sigfig(scr(6+j)/emev,7,0)
                  next=next+1
               enddo
            enddo
         endif
      endif
   endif
   lend=next-1
   lendp=lend

   !--store detailed photon production
   if (ntrp.ne.0) then
      call acelpp(next,matd,ngmt,nin)
      lendp=next-1
   endif

   !--change energies to mev, and
   !--fix up sig figs in summation cross sections.
   do i=1,nes
      xss(esz+i)=sigfig(xss(esz+i)/emev,9,0)
      xss(it+i)=sigfig(xss(it+i),9,0)
      xss(ic+i)=sigfig(xss(ic+i),9,0)
   enddo
   esz=esz+1
   end=lendp
   len2=lendp

   !--load particle production information
   call acelcp(next,matd,nin,za,awr)
   if (allocated(urd)) deallocate(urd)
   if (allocated(nudn)) deallocate(nudn)
   if (allocated(nup)) deallocate(nup)
   if (allocated(nut)) deallocate(nut)
   deallocate(scr)

   return
   end subroutine acelod

   subroutine acensd(ir,next,scr,nin,awp,ltt3,lttn,ltt,last,law,&
      ne,ie,il,iso,newfor)
   !-------------------------------------------------------------------
   ! Process this neutron scattering distribution.
   !-------------------------------------------------------------------
   use util ! provides error
   use mainio ! provides nsyso
   use endf ! provides endf routines and variables
   use acecm ! provides ptleg2,pttab2
   ! externals
   integer::ir,next,nin,ltt3,lttn,ltt,last,law,ne,ie,il,iso,newfor
   real(kr)::scr(*)
   real(kr)::awp
   ! internals
   integer::idone,j,nb,nw,lang,iint,nn,kk,nmu,m,n,i,ne1,ii,ll
   real(kr)::sum,renorm
   real(kr),parameter::emev=1.e6_kr
   real(kr),parameter::rmin=1.e-30_kr

   idone=0
   ne1=0
   do while (idone.eq.0)
      do j=1,ne
         if (newfor.eq.0) then
            call tab1io(nin,0,0,scr,nb,nw)
         else
            if (law.eq.7) then
               call tab1io(nin,0,0,scr,nb,nw)
            else if (ltt.eq.2) then
               call tab1io(nin,0,0,scr,nb,nw)
               ll=1+nw
               do while (nb.ne.0)
                  call moreio(nin,0,0,scr(ll),nb,nw)
                  ll=ll+nw
               enddo
               call pttab2(scr)
            else
               call listio(nin,0,0,scr,nb,nw)
               if (mfh.eq.6) then
                  lang=nint(scr(3))
                  if (lang.eq.0) then
                     call ptleg2(scr)
                  else
                     iint=lang-10
                     nn=nint(scr(6))
                     do kk=1,nw-4
                        scr(3+nw-kk)=scr(1+nw-kk)
                     enddo
                     scr(5)=1
                     scr(6)=nn
                     scr(7)=nn
                     scr(8)=iint
                     call pttab2(scr)
                  endif
               else
                  call ptleg2(scr)
               endif
            endif
         endif
         xss(ie+j)=sigfig(scr(2)/emev,7,0)
         if (law.eq.7) nmu=nint(scr(4))
         m=nint(scr(5))
         n=nint(scr(6))
         if (n.le.2.and.newfor.ne.1) then
            xss(il+j)=0
         else
            xss(il+j)=next-and+1
            if (newfor.ne.0) then
               xss(il+j)=-xss(il+j)
               iso=0
               iint=nint(scr(8))
               xss(next)=iint
               xss(next+1)=n
               if (next+2+3*n.gt.nxss) call error('acensd',&
                 'insufficient storage for angular distributions.',&
                 ' ')
               do i=1,n
                  xss(next+1+i)=sigfig(scr(5+2*m+2*i),7,0)
                  xss(next+1+n+i)=sigfig(scr(6+2*m+2*i),7,0)
                  if (xss(next+1+n+i).lt.rmin) xss(next+1+n+i)=0
                  if (i.eq.1) xss(next+1+2*n+i)=0
                  if (i.gt.1.and.iint.eq.1) then
                     sum=xss(next+1+2*n+i-1)+xss(next+1+n+i-1)&
                       *(xss(next+1+i)-xss(next+1+i-1))
                     xss(next+1+2*n+i)=sigfig(sum,7,0)
                  endif
                  if (i.gt.1.and.iint.eq.2) then
                     sum=xss(next+1+2*n+i-1)&
                       +(xss(next+1+n+i)+xss(next+1+n+i-1))&
                       *(xss(next+1+i)-xss(next+1+i-1))/2
                     xss(next+1+2*n+i)=sigfig(sum,7,0)
                  endif
               enddo
               renorm=1/xss(next+1+3*n)
               do i=1,n
                  xss(next+1+2*n+i)=renorm*xss(next+1+2*n+i)
                  xss(next+1+2*n+i)=sigfig(xss(next+1+2*n+i),9,0)
               enddo
               next=next+2+3*n
            else
               iso=0
               do i=1,n
                  xss(next)=sigfig(scr(2*i+2*m+5),6,0)
                  next=next+1
               enddo
               if (next.gt.nxss) call error('acensd',&
                 'insufficient storage for angular distributions.',&
                 ' ')
            endif
         endif
         if (law.eq.7) then
            do i=1,nmu
               call tab1io(nin,0,0,scr,nb,nw)
               do while (nb.ne.0)
                  call moreio(nin,0,0,scr,nb,nw)
               enddo
            enddo
         endif
      enddo

      !--for ltt=3, there are two series, one for the legendre
      !--data, and one for the tabulated distributions.  after
      !--the first series, read the the tab2 for the second,
      !--move the already loaded data up to make room, and
      !--start reading the tabulated distributions.
      idone=1
      if (ltt3.eq.3.and.lttn.eq.1) then
         ne1=ne
         call tab2io(nin,0,0,scr,nb,nw)
         ne=nint(scr(6))
         ltt=2
         lttn=2
         do i=il+ne1+1,next
            ii=next+il+ne1+1-i
            xss(ii+2*ne)=xss(ii)
         enddo
         do i=il+1,il+ne1
            ii=il+ne1+il+1-i
            xss(ii+ne)=xss(ii)-2*ne
         enddo
         next=next+2*ne
         il=il+ne1+ne
         ie=ie+ne1
         xss(last)=ne1+ne
         idone=0
      else if (ltt3.eq.3.and.lttn.eq.2) then
         ne=ne+ne1
      endif
   enddo

   !--shorten table for all-isotropic reactions.
   if (iso.ne.0) then
      xss(land+ir)=0
      next=last
   endif
   return
   end subroutine acensd

   subroutine acecpe(next,scr,nin,awr,awp,spi,ne,lidp,ie,il,nes)
   !-------------------------------------------------------------------
   ! Process this charged-particle elastic distribution.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use endf ! provides endf routines and variables
   use physics ! provides amassn,amu,hbar,pi,ev,clight
   use util ! provides sigfig
   ! externals
   integer::next,nin,ne,lidp,ie,il,nes
   real(kr)::scr(*),awr,awp,spi
   ! internals
   integer::llht,lld,j,ll,nb,nw,ltp,nl,ien,kk,iterp,jl,ione,itwo
   integer::i2s,ii,nrr,npp,jj,idis
   real(kr)::amass,e,f,xelas,cumm,amul,smul,sigcl,ratr,ratrl
   real(kr)::eht,amuu,pmu,ai,at,zt,zi,cc1,ee,cc2,wn,eta
   real(kr)::sigc,signi,signow,h,en
   real(kr)::xxs(200),yys(200)
   real(kr),parameter::emev=1.e6_kr
   real(kr),parameter::fm=1.e-12_kr

   write(nsyso,'(/'' working on charged-particle elastic'')')
   amass=awr/awp
   llht=1
   scr(llht)=0
   scr(llht+1)=0
   scr(llht+2)=0
   scr(llht+3)=0
   scr(llht+4)=1
   scr(llht+5)=ne
   scr(llht+6)=ne
   scr(llht+7)=5
   lld=llht+8+2*ne
   do j=1,ne
      ll=lld
      call listio(nin,0,0,scr(ll),nb,nw)
      ll=ll+nw
      do while (nb.ne.0)
         call moreio(nin,0,0,scr(ll),nb,nw)
         ll=ll+nw
      enddo
      e=scr(lld+1)
      ltp=nint(scr(lld+2))
      scr(lld+3)=lidp
      nl=nint(scr(lld+5))
      if (ltp.lt.12) then
         call ptlegc(scr(lld),awi,izai,awr,nint(za),spi)
         nl=nint(scr(lld+5))
         ll=lld+6+2*nl
      endif
      xss(ie+j)=sigfig(e/emev,7,0)
      xss(il+j)=-(next-and+1)
      ien=1
      do while (xss(esz+ien+1).lt.e.and.ien.lt.nes)
         ien=ien+1
      enddo
      f=(xss(esz+ien+1)-e)/(xss(esz+ien+1)-xss(esz+ien))
      xelas=xss(esz+3*nes+ien)*f+xss(esz+3*nes+ien+1)*(1-f)
      write(nsyso,'('' e,elas='',1p,2e12.4)') e,xelas
      write(nsyso,'(15x,''mu'',7x,''signi'',8x,''sigc'',8x,&
       &''sige'',8x,''ratr'',8x,''cumm'')')
      cumm=0
      amul=0
      smul=0
      sigcl=0
      ratr=0
      ratrl=0
      kk=0
      iterp=0
      eht=0
      do jl=1,nl
         ione=0
         do while (ione.eq.0)
            amuu=scr(lld+6+2*(jl-1))
            pmu=scr(lld+7+2*(jl-1))
            itwo=0
            do while (itwo.eq.0)
               i2s=nint(2*spi)
               ai=awi*amassn
               at=awr*amassn
               zt=nint(za/1000)
               zi=int(izai/1000)
               cc1=2*amu*ev*fm**2/hbar**2
               ee=(ev/10000000)*(clight/10)
               cc2=ee**4*amu/(2*hbar**2*ev)
               wn=at*sqrt(cc1*e*ai)/(ai+at)
               eta=zt*zi*sqrt(cc2*ai/e)
               sigc=0
               if (lidp.eq.0) sigc=(eta**2/wn**2)/(1-amuu)**2
               if (lidp.eq.1) sigc=((2*eta**2/wn**2)&
                 /(1-amuu**2))*((1+amuu**2)/(1-amuu**2)&
                 +(-1**i2s)*cos(eta*log((1+amuu)/(1-amuu)))&
                 /(2*spi+1))
               if (ltp.lt.12) pmu=pmu-sigc
               if (iterp.eq.1) then
                  signi=(ratr-1)*sigc
               else
                  signi=pmu*xelas
                  if (signi.lt.-sigc) signi=-sigc
               endif
               ratr=(sigc+signi)/sigc
               itwo=1
               if (jl.gt.1.and.iterp.eq.0.and.sigc.gt.abs(signi)&
                 .and.sigc.gt.(2*sigcl)) then
                  iterp=1
                  amuu=(amul+amuu)/2
                  ratr=(ratrl+ratr)/2
                  itwo=0
               endif
            enddo
            if (jl.gt.1) cumm=cumm+(amuu-amul)*(signi+sigc+smul)/2
            write(nsyso,'(5x,1p,6e12.4)')&
              amuu,signi,sigc,sigc+signi,ratr,cumm
            kk=kk+1
            scr(ll+3*(kk-1))=amuu
            scr(ll+1+3*(kk-1))=signi+sigc
            scr(ll+2+3*(kk-1))=cumm
            if (jl.gt.1) eht=eht&
              +(amuu-amul)*((1-amuu)*(signi+sigc)+(1-amul)*smul)/2
            amul=amuu
            sigcl=sigc
            smul=signi+sigc
            ratrl=ratr
            ione=1
            if (iterp.eq.1) then
               iterp=0
               ione=0
            endif
         enddo
      enddo
      scr(llht+6+2*j)=e
      scr(llht+7+2*j)=2*amass*(e/emev)*eht*2*pi/(1+amass)**2
      if (izai.eq.nint(za)) scr(llht+7+2*j)=0
      xxs(j)=e
      yys(j)=cumm*2*pi
      xss(next)=2
      xss(next+1)=kk
      do ii=1,kk
         xss(next+1+ii)=sigfig(scr(ll+3*(ii-1)),7,0)
         xss(next+1+kk+ii)=sigfig(scr(ll+1+3*(ii-1))/cumm,7,0)
         xss(next+1+2*kk+ii)=sigfig(scr(ll+2+3*(ii-1))/cumm,9,0)
      enddo
      next=next+2+3*kk
   enddo

   !--adjust the elastic and total cross sections
   !--use log-log interpolation for the elastic and heat
   !--note that the heating for the elastic recoil
   !--particles up through alphas comes from acelcp
   !--and is reported as zero here.
   write(nsyso,'(/'' new total, elastic, and heat'')')
   nrr=1
   npp=2
   do j=1,nes
      e=xss(esz+j)
      signi=xss(esz+3*nes+j)
      xss(esz+nes+j)=xss(esz+nes+j)-signi
      jj=1
      do while (xxs(jj+1).lt.e.and.jj.lt.ne)
         jj=jj+1
      enddo
      call terp1(xxs(jj),yys(jj),xxs(jj+1),yys(jj+1),e,signow,5)
      signow=sigfig(signow,7,0)
      xss(esz+nes+j)=sigfig(xss(esz+nes+j)+signow,9,0)
      xss(esz+3*nes+j)=signow
      call terpa(h,e,en,idis,scr(llht),npp,nrr)
      h=h/xss(esz+nes+j)
      if (izai.le.2004) h=0
      write(nsyso,'(5x,1p,4e12.4)') xss(esz+j),&
        xss(esz+nes+j),xss(esz+3*nes+j),h
      xss(esz+4*nes+j)=sigfig(h,7,0)
   enddo
   return
   end subroutine acecpe

   subroutine acelf5(next,i,matd,mt,q,nin,ismooth)
   !-------------------------------------------------------------------
   ! Process this reaction from File 5.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util ! provides sigfig
   use endf ! provides endf routines and variables
   ! externals
   integer::next,i,matd,mt,nin,ismooth
   real(kr)::q
   ! internals
   integer::nb,nw,nk,k,lf,m,n,jnt,ja,jb,j,l,nextn,nexd,ne,jscr,ki
   integer::ie,js,is,last,ix,ixx
   real(kr)::u,e,ep,renorm,efl,efh,emin,emax,tme,tmt
   real(kr)::dy,test,xm,ym,yt,dele,ta11
   character(60)::strng
   integer,parameter::ismax=20
   integer,parameter::jsmax=200
   real(kr)::xs(ismax),ys(ismax),xxs(jsmax),yys(jsmax)
   real(kr),dimension(:),allocatable::scr
   integer,parameter::nwscr=5000
   real(kr),parameter::emev=1.e6_kr
   real(kr),parameter::tol=.02e0_kr
   real(kr),parameter::tmin=1.e-12_kr
   real(kr),parameter::rmin=1.e-30_kr
   real(kr),parameter::zero=0

   !--allocate scratch storage area
   allocate(scr(nwscr))

   !--find desired reaction
   call findf(matd,5,mt,nin)
   call contio(nin,0,0,scr,nb,nw)
   nk=nint(scr(5))
   xss(ldlw+i-1)=next-dlw+1

   !--store each subsection
   do k=1,nk
      call tab1io(nin,0,0,scr,nb,nw)
      u=scr(1)
      lf=nint(scr(4))
      m=nint(scr(5))
      n=nint(scr(6))
      if (next+2*m+2*n+3.gt.nxss) call error('acelf5',&
        'insufficient space for energy distributions',' ')
      jnt=nint(scr(8))
      if (k.gt.1) xss(last)=next-dlw+1
      last=next
      xss(next)=0
      xss(next+1)=lf
      if (lf.eq.1) xss(next+1)=4
      if (lf.eq.12) xss(next+1)=4
      if (m.ne.1.or.jnt.ne.2) then
         ja=next+3
         xss(ja)=m
         jb=ja+m
         do j=1,m
            xss(ja+j)=scr(5+2*j)
            xss(jb+j)=scr(6+2*j)
         enddo
         next=jb+m+1
      else
         xss(next+3)=0
         next=next+4
      endif
      xss(next)=n
      l=5+2*m
      do j=1,n
         xss(j+next)=sigfig(scr(2*j+l)/emev,7,0)
         xss(j+n+next)=scr(2*j+1+l)
      enddo
      next=next+2*n+1
      xss(last+2)=next-dlw+1

      !--store law data according to type
      !--laws 1, 5, 7, 9, 11 only

      !--law 1...tabulated distributions
      if (lf.eq.1) then
         call tab2io(nin,0,0,scr,nb,nw)
         m=nint(scr(5))
         n=nint(scr(6))
         if (next+2*m+1.gt.nxss) call error('acelf5',&
           'insufficient space for energy distributions',' ')
         jnt=nint(scr(8))
         jnt=mod(jnt,10)
         if (jnt.gt.2) jnt=2
         if (m.ne.1.or.jnt.ne.2) then
            xss(next)=m
            do j=1,m
               xss(j+next)=scr(5+2*j)
               jnt=nint(scr(6+2*j))
               jnt=mod(jnt,10)
               if (jnt.gt.2) jnt=2
               xss(j+m+next)=jnt
            enddo
            next=next+1+2*m
         else
            xss(next)=0
            next=next+1
         endif
         xss(next)=n
         nextn=next+n
         nexd=nextn+n+1
         ne=n
         do j=1,ne
            call tab1io(nin,0,0,scr,nb,nw)
            jscr=1
            do while (nb.ne.0)
               jscr=jscr+nw
               if (jscr.gt.nwscr) call error('acelf5',&
                 'scratch storage exceeded reading lf=1',' ')
               call moreio(nin,0,0,scr(jscr),nb,nw)
            enddo
            e=c2h
            xss(next+j)=sigfig(e/emev,7,0)
            xss(nextn+j)=nexd-dlw+1
            m=n1h
            n=n2h
            jnt=nint(scr(6+2*m))
            if (ismooth.gt.0.and.mt.eq.18.and.jnt.eq.2) then
               write(nsyso,'('' may supplement the fission '',&
                 &''grid above 10 MeV using exponential shape '',&
                 &''if delta-E exceeds 200 keV.'')')
               ix=1
               do while (ix.lt.n)
                  jscr=5+2*m+2*ix
                  if (scr(jscr).lt.9.99e6_kr) then
                     ix=ix+1
                  else
                     dele=scr(jscr+2)-scr(jscr)
                     if (dele.gt.2.e5_kr) then
                        do ixx=n,ix+1,-1
                           scr(6+2*m+2*ixx+8)=scr(6+2*m+2*ixx)
                           scr(5+2*m+2*ixx+8)=scr(5+2*m+2*ixx)
                        enddo
                        ta11=zero
                        if (scr(jscr+11).eq.zero) then
                           ta11=scr(jscr+11)
                           scr(jscr+11)=1.e-6_kr*scr(jscr+1)
                        endif
                        do ixx=1,4
                           scr(jscr+2*ixx)=scr(jscr)+ixx*dele/5
                           call terp1(scr(jscr),scr(jscr+1),&
                                      scr(jscr+10),scr(jscr+11),&
                                      scr(jscr+2*ixx),scr(jscr+2*ixx+1),4)
                        enddo
                        if (ta11.ne.zero) scr(jscr+11)=ta11
                        n=n+4
                        ix=ix+5
                     else
                        ix=ix+1
                     endif
                  endif
               enddo
            endif
            xss(nexd)=jnt
            xss(nexd+1)=n
            nexd=nexd+1
            xss(nexd+1+2*n)=0
            do ki=1,n
               ep=scr(5+2*m+2*ki)
               if (ep.gt.e.and.q.lt.zero) then
                   write(nsyso,'(/'' ---warning from acelod---'',&
                     &6x,''mf5 ep.gt.e with negative q''/&
                     &6x,''mt='',i3,'' e='',1p,e12.4,'' ep='',e12.4/&
                     &6x,''patching...'')') mt,e/emev,ep/emev
                   ep=e-(n-ki)*1000
                   scr(5+2*m+2*ki)=ep
               endif
               xss(ki+nexd)=sigfig(scr(5+2*m+2*ki)/emev,7,0)
               xss(ki+n+nexd)=sigfig(scr(6+2*m+2*ki)*emev,7,0)
               if (xss(ki+n+nexd).lt.rmin) xss(ki+n+nexd)=0
               if (ki.gt.1.and.jnt.eq.1) xss(ki+2*n+nexd)=&
                 xss(ki+2*n-1+nexd)+scr(6+2*m+2*(ki-1))&
                 *(scr(5+2*m+2*ki)-scr(5+2*m+2*(ki-1)))
               if (ki.gt.1.and.jnt.eq.2) xss(ki+2*n+nexd)=&
                 xss(ki+2*n-1+nexd)+((scr(6+2*m+2*(ki-1))&
                 +scr(6+2*m+2*ki))/2)*(scr(5+2*m+2*ki)&
                 -scr(5+2*m+2*(ki-1)))
            enddo
            ! renormalize
            renorm=1
            if (xss(3*n+nexd).ne.zero) renorm=1/xss(3*n+nexd)
            do ki=1,n
               xss(ki+n+nexd)=sigfig(xss(ki+n+nexd)*renorm,7,0)
               xss(ki+2*n+nexd)=sigfig(xss(ki+2*n+nexd)*renorm,9,0)
            enddo
            nexd=nexd+3*n+1
         enddo
         next=nexd

      !--law 5...generalized evaporation spectrum
      else if (lf.eq.5) then
         call error('acelf5','sorry. acer cannot handle lf=5.',&
           'you will have to patch the evaluation to use lf=1.')
         call tab1io(nin,0,0,scr,nb,nw)
         m=nint(scr(5))
         n=nint(scr(5))
         if (next+2*m+2*n+2.gt.nxss) call error('acelf5',&
           'insufficient space for energy distributions',' ')
         jnt=nint(scr(8))
         if (m.ne.1.or.jnt.ne.2) then
            xss(next)=m
            do j=1,m
               xss(j+next)=scr(5+2*j)
               xss(j+m+next)=scr(6+2*j)
            enddo
            next=next+1+2*m
         else
            xss(next)=0
            next=next+1
         endif
         xss(next)=n
         do j=1,n
            xss(j+next)=sigfig(scr(5+2*j+2*m)/emev,7,0)
            xss(j+n+next)=sigfig(scr(6+2*j+2*m)/emev,7,0)
         enddo
         next=next+1+2*n
         call tab1io(nin,0,0,scr,nb,nw)
         m=nint(scr(5))
         n=nint(scr(6))
         if (next+n+1.gt.nxss) call error('acelf5',&
           'insufficient space for energy distributions',' ')
         xss(next)=n
         do j=1,n
            xss(j+next)=sigfig(scr(5+2*j+2*m)/emev,7,0)
         enddo
         next=next+1+n

      !--law 7...simple fission spectrum
      !--law 9...evaporation spectrum
      else if (lf.eq.7.or.lf.eq.9) then
         call tab1io(nin,0,0,scr,nb,nw)
         m=nint(scr(5))
         n=nint(scr(6))
         if (next+2*m+2*n+2.gt.nxss) call error('acelf5',&
           'insufficient space for energy distributions',' ')
         jnt=nint(scr(8))
         if (m.ne.1.or.jnt.ne.2) then
            xss(next)=m
            do j=1,m
               xss(j+next)=scr(5+2*j)
               xss(j+m+next)=scr(6+2*j)
            enddo
            next=next+2*m+1
         else
            xss(next)=0
            next=next+1
         endif
         xss(next)=n
         do j=1,n
            xss(j+next)=sigfig(scr(5+2*j+2*m)/emev,7,0)
            xss(j+n+next)=sigfig(scr(6+2*j+2*m)/emev,7,0)
         enddo
         next=next+1+2*n
         xss(next)=u/emeV
         next=next+1

      !--law 11...watt fission spectrum
      else if (lf.eq.11) then
         call tab1io(nin,0,0,scr,nb,nw)
         m=nint(scr(5))
         n=nint(scr(6))
         jnt=nint(scr(8))
         if (next+2*m+2*n+2.gt.nxss) call error('acelf5',&
           'insufficient space for energy distributions',' ')
         if (m.ne.1.or.jnt.ne.2) then
            xss(next)=m
            do j=1,m
               xss(j+next)=sigfig(scr(5+2*j)/emev,7,0)
               xss(j+m+next)=sigfig(scr(6+2*j)/emev,7,0)
            enddo
            next=next+2*m+1
         else
            xss(next)=0
            next=next+1
         endif
         xss(next)=n
         do j=1,n
            xss(j+next)=sigfig(scr(5+2*j+2*m)/emev,7,0)
            xss(j+n+next)=sigfig(scr(6+2*j+2*m)/emev,7,0)
         enddo
         next=next+1+2*n
         call tab1io(nin,0,0,scr,nb,nw)
         m=nint(scr(5))
         n=nint(scr(6))
         if (next+2*m+2*n+2.gt.nxss) call error('acelf5',&
           'insufficient space for energy distributions',' ')
         jnt=nint(scr(8))
         if (m.ne.1.or.jnt.ne.2) then
            xss(next)=m
            do j=1,m
               xss(j+next)=sigfig(scr(5+2*j+2*m)/emev,7,0)
               xss(j+n+next)=sigfig(scr(6+2*j+2*m)/emev,7,0)
            enddo
            next=next+1+2*m
         else
            xss(next)=0
            next=next+1
         endif
         xss(next)=n
         do j=1,n
            xss(j+next)=sigfig(scr(5+2*j+2*m)/emev,7,0)
            xss(j+n+next)=sigfig(scr(6+2*j+2*m)*emev,7,0)
         enddo
         next=next+1+2*n
         xss(next)=sigfig(u/emev,7,0)
         next=next+1

      !--law 12...madland-nix fission spectrum
      else if (lf.eq.12) then
         call tab1io(nin,0,0,scr,nb,nw)
         efl=scr(1)
         efh=scr(2)
         m=nint(scr(5))
         n=nint(scr(6))
         jnt=nint(scr(8))
         emin=scr(7+2*m)
         emax=scr(7+2*m+2*n-2)
         xss(next)=0
         next=next+1
         xss(next)=n
         nextn=next+n
         nexd=nextn+n+1
         ! convert madland-nix to ace law=4 using the given e grid
         ! and developing the e' grid adaptively.
         do ie=1,n
            tme=scr(5+2*m+2*ie)
            tmt=scr(6+2*m+2*ie)
            ! prime the adaptive stack
            renorm=0
            js=0
            is=3
            xs(3)=emin
            ys(3)=fmn(emin,efl,efh,tmt)
            xs(2)=emev
            ys(2)=fmn(xs(2),efl,efh,tmt)
            xs(1)=emax
            ys(1)=fmn(emax,efl,efh,tmt)
            ! carry out the adaptive linearization
            do while (is.gt.0)
               dy=0
               test=1
               if (is.gt.1.and.is.lt.ismax) then
                  xm=(xs(is-1)+xs(is))/2
                  ym=(ys(is-1)+ys(is))/2
                  yt=fmn(xm,efl,efh,tmt)
                  test=tol*abs(yt)+tmin
                  dy=abs(yt-ym)
               endif
               if (dy.gt.test) then
                  ! not converged.
                  ! add midpoint to stack and continue.
                  is=is+1
                  xs(is)=xs(is-1)
                  ys(is)=ys(is-1)
                  xs(is-1)=xm
                  ys(is-1)=yt
               else
                  ! converged.
                  ! move top point in stack to temporary spectrum.
                  js=js+1
                  if (js.gt.jsmax) js=jsmax
                  xxs(js)=xs(is)
                  yys(js)=ys(is)
                  if (js.gt.1) renorm=renorm&
                    +(xxs(js)-xxs(js-1))*(yys(js)+yys(js-1))/2
                  is=is-1
               endif
            enddo
            ! transfer spectrum to law 4
            renorm=1/renorm
            jnt=2
            xss(next+ie)=sigfig(tme/emev,7,0)
            xss(nextn+ie)=nexd-dlw+1
            xss(nexd)=jnt
            nexd=nexd+1
            xss(nexd)=js
            xss(nexd+1+2*js)=0
            do ki=1,js
               xss(ki+nexd)=sigfig(xxs(ki)/emev,7,0)
               xss(ki+js+nexd)=sigfig(yys(ki)*renorm*emev,7,0)
               if (ki.gt.1.and.jnt.eq.1) xss(ki+2*js+nexd)=&
                 xss(ki+2*js-1+nexd)&
                 +renorm*yys(ki)*(xxs(ki)-xxs(ki-1))
               if (ki.gt.1.and.jnt.eq.2) xss(ki+2*js+nexd)=&
                 xss(ki+2*js-1+nexd)+renorm*(yys(ki)+yys(ki-1))&
                 *(xxs(ki)-xxs(ki-1))/2
            enddo
            nexd=nexd+3*js+1
         enddo
         next=nexd

      !--illegal value of lf
      else
         write(strng,'(''illegal lf='',i2)') lf
         call error('acelf5',strng,' ')
      endif

   !--end of loop over laws
   enddo
   call tosend(nin,0,0,scr)
   deallocate(scr)
   return
   end subroutine acelf5

   subroutine acelf6(next,i,matd,mt,q,iza,izai,nin,newfor,ismooth)
   !-------------------------------------------------------------------
   ! Prepare generalized yields and energy-angle distributions for
   ! this reaction from File 6.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util ! provides error,mess
   use endf ! provides endf routines and variables
   use acecm ! provides bachaa,ptleg2,pttab2
   ! externals
   integer::next,i,matd,mt,iza,izai,nin,newfor,ismooth
   real(kr)::q
   ! internals
   integer::nb,nw,lct,nk,jscr,ivar,ik,idone,ikk,law,m,n,jnt
   integer::nnn,j,ii,ntyr,igyl,k,ishift,ip,ir,idis,ngyl,lgyl
   integer::npsx,nn,lang,lep,nextn,nexd,ne,nd,na,ncyc,nexcd
   integer::ki,iso,ik3,ii1,ia,ll,intmu,nmu,imu,mus,npep,intep
   integer::last,nx,ix
   integer::jp,jpn,jpp
   integer::nxyz1,nxyz2,nxyz3,nxyz4
   real(kr)::test,eemx,yield,xnext,xx,yy,y,xn,eyl,gyl,en
   real(kr)::apsx,step1,step2,xl,pl,yn,pn,rn,sum,ee
   real(kr)::ep,e,bzro,sfe,sfo,bbi,fbarcm,delfcm,akal,rkal
   real(kr)::emu1,emu2,fbl,ffl,fbcm,ffcm,akak,del,av,renorm
   real(kr)::zap,aa,test1,test2,test3,ex,fx,cx,cxx,val,dx
   real(kr)::e1,p1,e2,p2
   integer::loc(5)
   character(60)::strng
   real(kr),dimension(:),allocatable::scr
   integer,parameter::nwscr=1000000
   real(kr),parameter::emev=1.e6_kr
   real(kr),parameter::small=1.e-30_kr
   real(kr),parameter::eps=.001e0_kr
   real(kr),parameter::tmin=1.e-6_kr
   real(kr),parameter::etop=1.e10_kr
   real(kr),parameter::up=1.00001e0_kr
   real(kr),parameter::elow=1.e-5_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1
   real(kr),parameter::ten=10

   !--allocate scratch storage
   allocate(scr(nwscr))

   !--generalized yield from file 6
   call findf(matd,6,mt,nin)
   call contio(nin,0,0,scr,nb,nw)
   jp=nint(scr(3))
   jpn=mod(jp,10)
   jpp=(jp-jpn)/10
   lct=nint(scr(4))
   if (lct.gt.2) lct=2
   nk=nint(scr(5))
   jscr=1
   ivar=0
   ik=0
   idone=0
   ikk=0
   do while (idone.eq.0)
      ikk=ikk+1
      call tab1io(nin,0,0,scr(jscr),nb,nw)
      zap=scr(jscr)
      test=1
      test=test/1000
      law=nint(scr(jscr+3))
      if (zap.eq.zero.or.zap-izai.gt.test.or.&
          jpn.eq.2.or.(jpn.eq.1.and.ikk.gt.1)) then
         idone=1
      else

         !--wrong zap.  skip this subsection
         if (abs(zap-izai).gt.test) then
            do while (nb.ne.0)
               call moreio(nin,0,0,scr(jscr),nb,nw)
            enddo
            call skip6a(nin,0,0,scr(jscr),law)

         !--got desired zap
         else
            ik=ik+1
            loc(ik)=jscr
            if (law.ne.1.and.law.ne.6.and.law.ne.7)&
              call error('acelf6',&
              'illegal law for endf6 file6 neutrons',' ')
            m=nint(scr(jscr+4))
            n=nint(scr(jscr+5))
            jnt=nint(scr(jscr+5+2*m))
            if (next+2*m+2*n+3.gt.nxss) call error('acelf6',&
              'insufficient space for mf6 neutron yield',' ')
            jscr=jscr+nw
            do while (nb.ne.0)
               call moreio(nin,0,0,scr(jscr),nb,nw)
               jscr=jscr+nw
               if (jscr.gt.nwscr) call error('acelf6',&
                 'exceeded scratch storage',' ')
            enddo
            nnn=n
            if (jnt.eq.1) nnn=nnn-1
            if (nnn.ge.2) then
               do j=2,nnn
                  ii=loc(ik)+5+2*m+2*j
                  if (scr(ii)-scr(ii-2).ne.zero) ivar=1
               enddo
               eemx=scr(ii-1)
            endif
            ii=loc(ik)+5+2*m+2
            test=1
            test=test/1000
            if (ivar.eq.0.and.abs(scr(ii)-nint(scr(ii))).gt.test) ivar=2
            call skip6a(nin,0,0,scr(jscr),law)
         endif
      endif
      if (ikk.eq.nk) idone=1
   enddo
   if (ik.gt.1) write(nsyso,&
     '(/'' multiple mf6 subsections found for mt='',i3)') mth
   if (ivar.eq.1) write(nsyso,&
     '(/'' energy-dependent mf6 yields found for mt='',i3)') mth
   if (ivar.eq.2) write(nsyso,&
     '(/'' non-integer mf6 yields found for mt='',i3)') mth

   !--constant integer yield or yields
   if (ivar.le.0) then
      yield=0
      do j=1,ik
         ii=loc(j)+7+2*m
         yield=yield+scr(ii)
      enddo
      ntyr=nint(yield)
      if (mth.ge.6.and.mth.le.9) ntyr=2*ntyr
      if (mth.ge.46.and.mth.le.49) ntyr=2*ntyr
      if (mth.eq.18) ntyr=19
      if (mth.eq.18) lct=1 ! forces lab system for fission
      xss(tyr+i-1)=(3-2*lct)*ntyr
      if (law.eq.6) xss(tyr+i-1)=-ntyr

   !--generalized yield
   else
      ntyr=100+next-dlw+1
      xss(tyr+i-1)=(3-2*lct)*ntyr
      igyl=next
      xnext=0
      k=-1
      do while (k.le.0.or.xnext.lt.etop)
         k=k+1
         xx=xnext
         yy=0
         xnext=etop
         ishift=500
         do j=1,ik
            ii=loc(j)
            ip=2
            ir=1
            call terpa(y,xx,xn,idis,scr(ii),ip,ir)
            if (k.eq.0.and.xn.lt.xnext) xnext=xn
            if (k.ne.0) then
               yy=yy+y
               if (idis.eq.1.and.xx.lt.xn-xn/10000) xn=sigfig(xn,7,-1)
               if (xx.lt.eemx) xnext=xn
            endif
         enddo
         if (k.gt.0) then
            xss(next+2*k+ishift)=xx
            xss(next+2*k+1+ishift)=sigfig(yy,7,0)
            if (2*k+2.ge.ishift) call error('acelf6',&
              'storage exceeded for generalized yield',' ')
         endif
      enddo
      xss(next)=0
      xss(next+1)=k
      do j=1,k
         eyl=xss(next+2*j+ishift)
         xss(next+1+j)=sigfig(eyl/emev,7,0)
         xss(next+1+k+j)=sigfig(xss(next+2*j+1+ishift),7,0)
      enddo
      next=next+2*k+2
   endif

   !--return to the beginning of the section
   call findf(matd,6,mt,nin)
   call contio(nin,0,0,scr,nb,nw)
   xss(ldlw+i-1)=next-dlw+1
   jp=nint(scr(3))
   jpn=mod(jp,10)
   jpp=(jp-jpn)/10
   ik=0
   ikk=0
  110 continue
   ikk=ikk+1
   if (ikk.gt.nk) go to 130
   jscr=1
   call tab1io(nin,0,0,scr(jscr),nb,nw)
   zap=scr(1)
   test=1
   test=test/1000
   law=nint(scr(4))
   if (zap.eq.zero.or.zap-izai.gt.test) go to 130
   if (jpn.eq.2.or.(jpn.eq.1.and.ikk.gt.1)) go to 130
   if (abs(zap-izai).lt.test) go to 120
   ! skip this subsection
   do while (nb.ne.0)
      call moreio(nin,0,0,scr(jscr),nb,nw)
   enddo
   call skip6a(nin,0,0,scr(jscr),law)
   go to 110
  120 continue
   ik=ik+1
   m=nint(scr(5))
   n=nint(scr(6))
   jnt=nint(scr(8))
   jscr=jscr+nw
   do while (nb.ne.0)
      call moreio(nin,0,0,scr(jscr),nb,nw)
      jscr=jscr+nw
      if (jscr.gt.nwscr) call error('acelf6',&
        'exceeded scratch storage',' ')
   enddo

   !--store yield
   if (ik.gt.1) xss(last)=next-dlw+1
   last=next
   xss(next)=0
   if (law.eq.1) then
      xss(next+1)=44
   else if (law.eq.7) then
      xss(next+1)=67
   else
      xss(next+1)=66
   endif
   next=next+3
   if (ivar.le.0) then
      if (m.ne.1.or.jnt.ne.2) then
         xss(next)=m
         do k=1,m
            xss(next+k)=scr(5+2*k)
            xss(next+m+k)=scr(6+2*k)
         enddo
         next=next+1+2*m
      else
         xss(next)=0
         next=next+1
      endif
      xss(next)=n
      do k=1,n
         xss(next+k)=sigfig(scr(5+2*m+2*k)/emev,7,0)
         xss(next+n+k)=sigfig(scr(6+2*m+2*k)/yield,7,0)
      enddo
      next=next+1+2*n
   else
      ngyl=nint(xss(igyl+1))
      xss(next)=0
      next=next+1
      xss(next)=ngyl
      do j=1,ngyl
         eyl=xss(igyl+1+j)*emev
         if (j.eq.1) eyl=up*eyl
         if (j.eq.ngyl) eyl=eyl/up
         gyl=xss(igyl+1+ngyl+j)
         lgyl=igyl+2
         if (j.eq.1) call terp1(xss(lgyl),xss(lgyl+ngyl),&
           xss(lgyl+1),xss(lgyl+1+ngyl),eyl/emev,gyl,2)
         lgyl=igyl+1+ngyl
         if (j.eq.ngyl) call terp1(xss(lgyl-1),xss(lgyl-1+ngyl),&
           xss(lgyl),xss(lgyl+ngyl),eyl/emev,gyl,2)
         ir=1
         ip=2
         call terpa(y,eyl,en,idis,scr,ip,ir)
         if (j.eq.1) eyl=eyl/up
         if (j.eq.ngyl) eyl=up*eyl
         xss(next+j)=sigfig(eyl/emev,7,0)
         yy=y
         if (gyl.gt.0) yy=yy/gyl
         xss(next+ngyl+j)=sigfig(yy,7,0)
      enddo
      next=next+1+2*ngyl
   endif

   !--phase-space distribution
   if (law.eq.6) then
      call contio(nin,0,0,scr,nb,nw)
      apsx=scr(1)
      npsx=nint(scr(6))
      xss(last+2)=next-dlw+1
      xss(next)=npsx
      xss(next+1)=apsx
      xss(next+2)=2
      next=next+3
      step1=ten**(one/5)
      step2=one/50
      xx=elow
      n=1
      test1=one+one/100000
      test2=one/10-one/1000000
      test3=one-one/10000
      do while (xx.lt.test1)
         n=n+1
         if (xx.lt.test2) then
             xx=xx*step1
         else
            xx=xx+step2
         endif
      enddo
      nn=n
      xss(next)=nn
      xl=0
      pl=0
      yn=0
      n=1
      xss(next+n)=xl
      xss(next+nn+n)=pl
      xss(next+2*nn+n)=yn
      xx=elow
      do while (xx.lt.test1)
         n=n+1
         if (xx.gt.test3) then
            xx=1
            pn=0
         else
            rn=3
            pn=sqrt(xx)*(1-xx)**(rn*npsx/2-4)
         endif
         yn=yn+(xx-xl)*(pn+pl)/2
         xss(next+n)=sigfig(xx,7,0)
         xss(next+nn+n)=pn
         xss(next+2*nn+n)=yn
         xl=xx
         pl=pn
         test=one/10
         test=test-test/10000
         if (xx.lt.test) then
            xx=xx*step1
         else
            xx=xx+step2
         endif
      enddo
      sum=yn
      do k=1,nn
         xss(next+nn+k)=sigfig(xss(next+nn+k)/sum,7,0)
         xss(next+2*nn+k)=sigfig(xss(next+2*nn+k)/sum,9,0)
      enddo
      next=next+1+3*nn

   !--other energy-angle distributions
   else
      call tab2io(nin,0,0,scr,nb,nw)
      lang=nint(scr(3))
      if (law.eq.1.and.(lang.lt.1.or.(lang.gt.2.and.lang.lt.11)&
        .or.lang.gt.13)) call error('acelf6',&
        'only lang=1,2,11-13 allowed for endf-6 file 6 neutrons',&
        ' ')
      if (newfor.eq.0.and.law.eq.1.and.lct.eq.2.and.lang.ne.2) then
         write(nsyso,'(/'' converting to kalbach:'',&
           &'' mt ='',i4)') mth
         xss(tyr+i-1)=-abs(xss(tyr+i-1))
      endif
      if (newfor.eq.1.and.law.eq.1.and.lang.ne.2) xss(last+1)=61
      lep=nint(scr(4))
      m=nint(scr(5))
      n=nint(scr(6))
      if (next+2*m+1.gt.nxss) call error('acelf6',&
        'insufficient space for mf6 tab2',' ')
      jnt=nint(scr(8))
      jnt=mod(jnt,10)
      if (jnt.gt.2) jnt=2
      xss(last+2)=next-dlw+1
      if (m.ne.1.or.jnt.ne.2) then
         xss(next)=m
         do k=1,m
            xss(next+k)=scr(5+2*k)
            jnt=nint(scr(6+2*k))
            jnt=mod(jnt,10)
            if (jnt.gt.2) jnt=2
            xss(next+m+k)=jnt
         enddo
         next=next+1+2*m
      else
         xss(next)=0
         next=next+1
      endif
      xss(next)=n
      nextn=next+n
      nexd=nextn+n+1
      ne=n

      !--loop over incident energies
      do j=1,ne
         if (law.ne.7) then
            jscr=1
            call listio(nin,0,0,scr(jscr),nb,nw)
            jscr=jscr+nw
            do while (nb.ne.0)
               call moreio(nin,0,0,scr(jscr),nb,nw)
               jscr=jscr+nw
            enddo
            xss(next+j)=sigfig(c2h/emev,7,0)
            ee=c2h/emev
            xss(nextn+j)=nexd-dlw+1
            nd=nint(scr(3))
            na=nint(scr(4))
            ncyc=na+2
            nx=nint(scr(5))
            n=nint(scr(6))
            if (next+5*n.gt.nxss) call error('acelf6',&
              'insufficient storage for energy distributions.',' ')

            ! extend low histogram bins as sqrt(e) using log energy scale
            ! only do this for outgoing neutrons with law=1, lang=2
            ! only do this if there are no discrete data
            if (ismooth.gt.0.and.law.eq.1.and.lang.eq.2.and.&
                lep.eq.1.and.zap.eq.1.and.nd.eq.0) then
               fx=.8409
               ex=40
               cx=scr(7+ncyc)*scr(8)
               do while (n.gt.2)
                  cxx=cx+scr(8+ncyc)*(scr(7+2*ncyc)-scr(7+ncyc))
                  if (abs(cxx/scr(7+2*ncyc)**1.5-cx/scr(7+ncyc)**1.5)&
                    .gt.cx/scr(7+ncyc)**1.5/50) exit
                  scr(8)=(scr(8)*scr(7+ncyc)&
                    +scr(8+ncyc)*(scr(7+2*ncyc)&
                    -scr(7+ncyc)))/scr(7+2*ncyc)
                  do ix=1,nx-2*ncyc
                     scr(6+ix+ncyc)=scr(6+ix+2*ncyc)
                  enddo
                  cx=cxx
                  nx=nx-ncyc
                  n=n-1
               enddo
               write(nsyso,'('' extending histograms as sqrt(E) below'',&
                 &1p,e10.2,'' MeV for E='',e10.2,'' MeV mt='',i3)')&
                 scr(7+ncyc)/emev,ee,mt
               do while (scr(7+ncyc).gt.ex)
                  do ix=nx,1,-1
                     scr(6+ncyc+ix)=scr(6+ix)
                  enddo
                  scr(7+ncyc)=fx*scr(7+2*ncyc)
                  scr(7+ncyc)=sigfig(scr(7+ncyc),6,0)
                  val=scr(8)
                  scr(8)=sqrt(fx)*val
                  scr(8)=sigfig(scr(8),6,0)
                  scr(8+ncyc)=(1-fx*sqrt(fx))*val/(1-fx)
                  scr(8+ncyc)=sigfig(scr(8+ncyc),6,0)
                  nx=nx+ncyc
                  jscr=jscr+ncyc
                  n=n+1
               enddo
            ! extend to lower energy as sqrt(e) using linear interpol
            ! only do this for outgoing neutrons with law=1, lang=2
            ! only do this if there are no discrete data
            else if (ismooth.gt.0.and.law.eq.1.and.lang.eq.2.and.&
                     lep.eq.2.and.zap.eq.1.and.n.gt.3.and.nd.eq.0) then
               ex=40
               fx=0.50
               nn=0
               ! make room for initial data at zero energy, then
               ! insert those zero energy data
               if (scr(7).gt.ex) then
                  write(nsyso,'('' extending lin-lin as sqrt(E) '',&
                   &''below'',1p,e10.2,'' MeV for E='',e10.2,'' MeV mt='',i3)&
                   &')scr(7)/emev,ee,mt
                  do ix=nx,1,-1
                     scr(6+ncyc+ix)=scr(6+ix)
                  enddo
                  do ii=1,ncyc
                     scr(6+ii)=0
                  enddo
                  nx=nx+ncyc
                  jscr=jscr+ncyc
                  n=n+1
               endif
               ! insert new data between the original E1 and zero.
               ! continue until reach an energy below ex.
               do while (scr(7+ncyc).gt.ex)
                  nn=nn+1
                  do ix=nx,ncyc+1,-1
                     scr(6+ncyc+ix)=scr(6+ix)
                  enddo
                  scr(7+ncyc)=fx*scr(7+2*ncyc)
                  scr(7+ncyc)=sigfig(scr(7+ncyc),6,0)
                  do ii=2,ncyc
                     scr(6+ncyc+ii)=scr(6+2*ncyc+ii)*&
                       sqrt(scr(7+ncyc)/scr(7+2*ncyc))
                  enddo
                  nx=nx+ncyc
                  jscr=jscr+ncyc
                  n=n+1
                  scr(5)=nx
                  scr(6)=n
               enddo
               ! if new data were inserted, need to renormalize
               ! the emission distribution back to unity.
               if (nn.gt.0) then
                  cxx=0
                  e1=scr(7)
                  p1=scr(8)
                  do ix=1+ncyc,nx,ncyc
                     e2=scr(6+ix)
                     p2=scr(7+ix)
                     cxx=cxx+(e2-e1)*(p1+p2)/2
                     e1=e2
                     p1=p2
                  enddo
                  if (cxx.ne.0) then
                     cxx=1/cxx
                     do ix=2+ncyc,nx,ncyc
                        scr(6+ix)=cxx*scr(6+ix)
                     enddo
                  endif
               endif
            endif

            xss(nexd)=lep+10*nd
            xss(nexd+1)=n
            nexd=nexd+1
            xss(nexd+1+2*n)=0
            nexcd=nexd+4*n+1

            !--loop over secondary energies
            do ki=1,n
               ep=scr(7+ncyc*(ki-1))
               e=ee*emev
               if (ep.gt.e-e/1000.and.ki.lt.n.and.mth.ne.5.and.q.lt.zero) then
                  write(nsyso,'(/'' ---warning from acelf6---'',&
                    &6x,''mf6 ep.gt.e with negative q''/&
                    &6x,''mt='',i3,'' e='',1p,e12.4,'' ep='',e12.4/&
                    &6x,''patching...'')') mt,e/emev,ep/emev
                  ep=e-(n-ki)*1000
                  scr(7+ncyc*(ki-1))=ep
               else if (ep.gt.e.and.ki.eq.n.and.mth.ne.5.and.q.lt.zero) then
                  write(nsyso,'(/'' ---warning from acelf6---'',&
                    &6x,''mf6 ep.gt.e with negative q''/&
                    &6x,''mt='',i3,'' e='',1p,e12.4,'' ep='',e12.4/&
                    &6x,''patching...'')') mt,e/emev,ep/emev
                  ep=e-(n-ki)*1000
                  scr(7+ncyc*(ki-1))=ep
               else if (ep.gt.e.and.mth.eq.5) then
                  write(nsyso,'(/'' ---warning from acelf6---'',&
                    &6x,''mf6/mt5 ep.gt.e''/&
                    &6x,''mt='',i3,'' e='',1p,e12.4,'' ep='',e12.4/&
                    &6x,''leaving it as is...'')') mt,e/emev,ep/emev
               endif
               xss(ki+nexd)=sigfig(scr(7+ncyc*(ki-1))/emev,7,0)
               if (ki.le.nd) xss(ki+n+nexd)=scr(8+ncyc*(ki-1))
               if (ki.gt.nd) xss(ki+n+nexd)&
                 =sigfig(scr(8+ncyc*(ki-1))*emev,7,0)
               test=xss(ki+n+nexd)
               if (test.gt.0.and.test.lt.small) xss(ki+n+nexd)=small
               bzro=xss(nexd+ki+n)/emev
               if (ki.le.nd) xss(ki+2*n+nexd)=xss(ki+2*n-1+nexd)&
                 +scr(8+ncyc*(ki-1))
               if (nd.gt.0.and.ki.eq.nd+1)&
                  xss(ki+2*n+nexd)=xss(ki+2*n-1+nexd)
               if (ki.gt.nd+1.and.lep.eq.1)&
                 xss(ki+2*n+nexd)=xss(ki+2*n-1+nexd)&
                 +scr(8+ncyc*(ki-2))*(scr(7+ncyc*(ki-1))&
                 -scr(7+ncyc*(ki-2)))
               if (ki.gt.nd+1.and.lep.eq.2)&
                 xss(ki+2*n+nexd)=xss(ki+2*n-1+nexd)&
                 +((scr(8+ncyc*(ki-2))+scr(8+ncyc*(ki-1)))/2)&
                 *(scr(7+ncyc*(ki-1))-scr(7+ncyc*(ki-2)))
               ep=xss(ki+nexd)

               !--process according to lang

               !--distribution given in kalbach format
               if (lang.eq.2) then
                  xss(ki+3*n+nexd)=scr(9+ncyc*(ki-1))
                  if (na.eq.2) then
                     aa=scr(10+ncyc*(ki-1))
                  else
                     aa=bachaa(1,1,iza,ee,ep)
                  endif
                  xss(ki+4*n+nexd)=sigfig(aa,7,0)

               !--convert legendre distribution to kalbach form
               else if (lang.eq.1.and.newfor.ne.1) then
                  iso=1
                  do ik3=1,na
                     if (scr(8+ik3+ncyc*(ki-1)).ne.zero) iso=0
                  enddo
                  if (bzro.eq.0.) iso=1

                  !--isotropic case
                  if (iso.ne.0) then
                     xss(ki+3*n+nexd)=0
                     xss(ki+4*n+nexd)=elow

                  !--nonisotropic case
                  else
                     sfe=1
                     sfe=sfe/2
                     sfo=0
                     do  ii1=1,na
                        bbi=(2*ii1+1)*scr(8+ii1+ncyc*(ki-1))&
                          /(2*bzro)
                        if (ii1.eq.2*(ii1/2)) then
                           sfe=sfe+bbi
                        else
                           sfo=sfo+2*bbi
                        endif
                     enddo
                     fbarcm=sfe
                     delfcm=sfo
                     ! don't accept unreasonable distributions
                     if (delfcm.lt.zero) then
                        write(strng,'('' for mt='',i3,''  e='',1p,&
                          &e12.4,''  eprime='',e12.4)') mth,ee,ep
                        call mess('acelf6',&
                          'converted cm distribution unreasonable',&
                          strng)
                        delfcm=0
                     endif
                     call fndar1(akal,rkal,fbarcm,delfcm,ee,ep)
                     xss(ki+3*n+nexd)=rkal
                     xss(ki+4*n+nexd)=akal
                  endif

               !--convert tabulated distribution to kalbach form
               else if (lang.gt.2.and.newfor.eq.0) then
                  emu1=scr(9+ncyc*(ki-1))
                  emu2=scr(7+na+ncyc*(ki-1))
                  fbl=scr(10+ncyc*(ki-1))
                  ffl=scr(8+na+ncyc*(ki-1))
                  test=1
                  test=test/1000
                  if (abs(emu1+1).ge.test.or.abs(emu2-1).ge.test) then
                     fbcm=fbl
                     ffcm=ffl
                     call fndar2(akak,rkal,fbcm,ffcm,emu1,emu2,ee,ep)
                     call mess('acelf6',&
                       'tabulated angular distribution',&
                       'does not extend over entire cosine range.')
                  else
                     fbarcm=(ffl+fbl)/2
                     delfcm=ffl-fbl
                     if (delfcm.lt.zero) then
                        write(strng,'('' for mt='',i3,''  e='',1p,&
                          &e12.4,''  eprime='',e12.4)') mth,ee,ep
                        call mess('acelf6',&
                          'converted cm distribution unreasonable',&
                          strng)
                        delfcm=0
                     endif
                     call fndar1(akal,rkal,fbarcm,delfcm,ee,ep)
                     xss(ki+3*n+nexd)=sigfig(rkal,7,0)
                     xss(ki+4*n+nexd)=sigfig(akal,7,0)
                  endif

               !--convert legendre distribution to law 61
               else if (lang.eq.1.and.newfor.eq.1) then
                  scr(jscr)=0
                  scr(jscr+1)=ep
                  scr(jscr+2)=0
                  scr(jscr+3)=0
                  scr(jscr+4)=na
                  scr(jscr+5)=0
                  do ia=1,na
                     ll=8+ncyc*(ki-1)
                     scr(jscr+5+ia)=0
                     if (scr(ll).ne.zero) then
                        scr(jscr+5+ia)=scr(ll+ia)/scr(ll)
                     endif
                  enddo
                  call ptleg2(scr(jscr))
                  xss(ki+3*n+nexd)=nexcd-dlw+1
                  intmu=2
                  xss(nexcd)=intmu
                  nmu=nint(scr(jscr+5))
                  xss(nexcd+1)=nmu
                  do imu=1,nmu
                     xss(nexcd+1+imu)=sigfig(scr(jscr+6+2*imu),7,0)
                     xss(nexcd+1+nmu+imu)=sigfig(scr(jscr+7+2*imu),7,0)
                     if (imu.eq.1) then
                        xss(nexcd+1+2*nmu+imu)=0
                     else if (imu.eq.nmu) then
                        xss(nexcd+1+2*nmu+imu)=1
                     else
                        del=scr(jscr+6+2*imu)-scr(jscr+4+2*imu)
                        av=(scr(jscr+7+2*imu)+scr(jscr+5+2*imu))/2
                        xss(nexcd+1+2*nmu+imu)=&
                          xss(nexcd+1+2*nmu+imu-1)+del*av
                        xss(nexcd+1+2*nmu+imu)=&
                          sigfig(xss(nexcd+1+2*nmu+imu),7,0)
                     endif
                  enddo
                  nexcd=nexcd+2+3*nmu

               !--convert tabulated distribution to law 61
               else if (lang.gt.2.and.newfor.eq.1) then
                  xss(ki+3*n+nexd)=nexcd-dlw+1
                  ll=7+ncyc*(ki-1)
                  intmu=lang-10
                  xss(nexcd)=intmu
                  nmu=na/2
                  xss(nexcd+1)=nmu
                  do imu=1,nmu
                     xss(nexcd+1+imu)=scr(ll+2*imu)
                     xss(nexcd+1+nmu+imu)=scr(ll+2*imu+1)
                     if (imu.eq.1) then
                        sum=0
                        xss(nexcd+1+2*nmu+imu)=0
                     else
                        del=scr(ll+2*imu)-scr(ll+2*imu-2)
                        if (intmu.eq.1) then
                           sum=sum+del*scr(ll+1+2*imu-2)
                           xss(nexcd+1+2*nmu+imu)=sum
                        else
                           av=(scr(ll+1+2*imu)+scr(ll+1+2*imu-2))/2
                           sum=sum+del*av
                           xss(nexcd+1+2*nmu+imu)=sum
                        endif
                     endif
                  enddo
                  do imu=1,nmu
                     xss(nexcd+1+imu)=&
                       sigfig(xss(nexcd+1+imu)/sum,7,0)
                     xss(nexcd+1+nmu+imu)=&
                       sigfig(xss(nexcd+1+nmu+imu)/sum,7,0)
                     xss(nexcd+1+2*nmu+imu)=&
                       sigfig(xss(nexcd+1+2*nmu+imu)/sum,7,0)
                  enddo
                  nexcd=nexcd+2+3*nmu
               endif
            enddo
            if (xss(3*n+nexd).ne.zero) then
               renorm=1/xss(3*n+nexd)
            else
               renorm=1
            endif
            do ki=1,n
               xss(ki+n+nexd)=sigfig(xss(ki+n+nexd)*renorm,7,0)
               xss(ki+2*n+nexd)=sigfig(xss(ki+2*n+nexd)*renorm,9,0)
            enddo
            nexd=nexd+5*n+1
            if (nint(xss(last+1)).eq.61) nexd=nexcd

         !--law 7, angle-energy format
         else
            jscr=1
            call tab1io(nin,0,0,scr(jscr),nb,nw)
            intmu=l1h
            nmu=l2h
            ee=c2h/emev
            xss(next+j)=sigfig(ee,7,0)
            xss(nextn+j)=nexd-dlw+1
            xss(nexd)=intmu
            xss(nexd+1)=nmu
            nexd=nexd+2
            mus=nexd
            nexd=nexd+2*nmu
            do imu=1,nmu
               jscr=1
               call tab1io(nin,0,0,scr(jscr),nb,nw)
               npep=n2h
               intep=nint(scr(jscr+7))
               jscr=jscr+nw
               do while (nb.ne.0)
                  call moreio(nin,0,0,scr(jscr),nb,nw)
                  jscr=jscr+nw
               enddo
               xss(mus+imu-1)=c2h
               xss(mus+nmu+imu-1)=nexd-dlw+1
               xss(nexd)=intep
               xss(nexd+1)=npep
               nexd=nexd+1
               xss(nexd+1+2*npep)=0
               do ki=1,npep
                  xss(nexd+ki)=sigfig(scr(9+2*(ki-1))/emev,7,0)
                  xss(nexd+npep+ki)=&
                    sigfig(scr(10+2*(ki-1))*emev,7,0)
                  if (ki.ne.1) then
                     if (intep.eq.1) xss(nexd+2*npep+ki)=&
                       xss(nexd+2*npep+ki-1)+scr(10+2*(ki-2))&
                       *(scr(9+2*(ki-1))-scr(9+2*(ki-2)))
                     if (intep.eq.2) xss(nexd+2*npep+ki)=&
                       xss(nexd+2*npep+ki-1)&
                       +(scr(10+2*(ki-2))+scr(10+2*(ki-1)))&
                       *(scr(9+2*(ki-1))-scr(9+2*(ki-2)))/2
                  endif
               enddo
               renorm=1/xss(nexd+3*npep)
               do ki=1,npep
                  xss(nexd+npep+ki)=&
                    sigfig(renorm*xss(nexd+npep+ki),7,0)
                  xss(nexd+2*npep+ki)=&
                    sigfig(renorm*xss(nexd+2*npep+ki),9,0)
               enddo
               nexd=nexd+3*npep+1
            enddo
         endif
      enddo
      next=nexd
   endif
   go to 110
  130 continue
   call tosend(nin,0,0,scr)

   !--finished
   deallocate(scr)
   return
   end subroutine acelf6

   subroutine skip6a(nin,nout,nscr,a,law)
   !-------------------------------------------------------------------
   ! Special version of skip6 for special version of File 6 used
   ! in ACER.  Law=7 has a TAB1 containing the angular distribution
   ! instead of the normal TAB2 for each incident energy.
   ! Skip the next subsection in the current section (MT).
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr,law
   real(kr)::a(*)
   ! internals
   integer::nb,nw,ne,ie,nmu,imu

   if (law.eq.6) then
      call contio(nin,nout,nscr,a(1),nb,nw)
   else if (law.eq.1.or.law.eq.2.or.law.eq.5) then
      call tab2io(nin,nout,nscr,a(1),nb,nw)
      ne=n2h
      do ie=1,ne
         call listio(nin,nout,nscr,a(1),nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a(1),nb,nw)
         enddo
      enddo
   else if (law.eq.7) then
      call tab2io(nin,nout,nscr,a(1),nb,nw)
      ne=n2h
      do ie=1,ne
         call tab1io(nin,nout,nscr,a(1),nb,nw)
         nmu=nint(a(4))
         do imu=1,nmu
            call tab1io(nin,nout,nscr,a(1),nb,nw)
            do while (nb.ne.0)
               call moreio(nin,nout,nscr,a(1),nb,nw)
            enddo
         enddo
      enddo
   endif
   return
   end subroutine skip6a

   subroutine ptlegc(c,awp,izai,awt,iza,spi)
   !-------------------------------------------------------------------
   ! For charged particle Legendre coefficient representations,
   ! reconstruct the angular distribution adaptively.  The
   ! distribution given is the actual elastic cross section.
   !-------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::izai,iza
   real(kr)::c(*),awp,awt,spi
   ! internals
   integer::ii,i,nn,idone,j,jj,k
   real(kr)::dy,dm,xm,yt,test,check,f,diff,ym,dco
   integer,parameter::imax=20
   real(kr)::x(imax),y(imax)
   integer,parameter::maxang=2000
   real(kr)::aco(maxang),cprob(maxang)
   real(kr),parameter::tol1=.001e0_kr
   real(kr),parameter::tol2=.01e0_kr
   real(kr),parameter::one=1.e0_kr
   real(kr),parameter::half=.5e0_kr
   real(kr),parameter::hund=.01e0_kr
   real(kr),parameter::umin=.96e0_kr
   real(kr),parameter::zero=0

   !--adaptive reconstruction of angular distribution
   ii=0
   ! prime the adaptive stack
   i=2
   x(2)=-1
   if (iza.eq.izai) x(2)=-umin
   call coul(x(2),y(2),c,awp,izai,awt,iza,spi)
   x(1)=umin
   call coul(x(1),y(1),c,awp,izai,awt,iza,spi)
   ! carry out the adaptive reconstruction
   do while (i.gt.0)
      dy=0
      if (i.gt.1.and.i.lt.imax) then
         dm=x(i-1)-x(i)
         xm=half*(x(i-1)+x(i))
         xm=sigfig(xm,3,0)
         if (xm.gt.x(i).and.xm.lt.x(i-1)) then
            ym=half*(y(i-1)+y(i))
            call coul(xm,yt,c,awp,izai,awt,iza,spi)
            test=tol1*abs(yt)
            dy=abs(yt-ym)
            if (dm.gt.hund) dy=2*test
            if (ym.ne.zero.and.yt/ym.gt.one*5) dy=2*test
            if (ym.ne.zero.and.yt/ym.lt.one/5) dy=2*test
         endif
      endif
      if (dy.gt.test) then
         ! not converged.
         ! add midpoint to stack and continue.
         i=i+1
         x(i)=x(i-1)
         y(i)=y(i-1)
         x(i-1)=xm
         y(i-1)=yt
      else
         ! converged.
         ! use top point in stack.
         ii=ii+1
         if (ii.gt.maxang) call error('ptlegc',&
           'too many coulomb angles',' ')
         aco(ii)=x(i)
         cprob(ii)=y(i)
         i=i-1
      endif
   enddo
   nn=ii-1

   !--now thin the distribution to a coarser tolerance
   i=1
   ii=1
   idone=0
   dco=0
   do while (i.lt.nn-1.and.idone.eq.0)
      j=i+1
      check=0
      do while (j.lt.nn+1.and.check.le.0.and.dco.le.one/2)
         j=j+1
         jj=j-1
         dco=aco(j)-aco(i)
         if (dco.le.one/2) then
            k=i
            do while (k.lt.j-1.and.check.le.0)
               k=k+1
               f=(aco(j)-aco(k))/dco
               test=f*cprob(i)+(1-f)*cprob(j)
               diff=tol2*abs(cprob(k))
               check=abs(test-cprob(k))-diff
            enddo
         endif
      enddo
      if (dco.gt.one/2.or.check.gt.zero) then
         i=jj
         ii=ii+1
         aco(ii)=aco(i)
         cprob(ii)=cprob(i)
         dco=0
      else
         idone=1
      endif
   enddo
   i=nn+1
   ii=ii+1
   aco(ii)=aco(i)
   cprob(ii)=cprob(i)

   !--load the new distribution back into the input array
   c(5)=2*ii
   c(6)=ii
   do i=1,ii
      c(5+2*i)=aco(i)
      c(6+2*i)=cprob(i)
   enddo
   return
   end subroutine ptlegc

   subroutine coul(x,y,c,awp,izap,awt,izat,spi)
   !-------------------------------------------------------------------
   ! Compute the charged-particle elastic scattering cross section
   ! and the normal Coulomb cross section from the LTP=1 or LTP=2
   ! representations of File 6.
   !-------------------------------------------------------------------
   use physics ! provides amu,hbar,ev,clight,amassn
   use mathm ! provides legndr
   ! externals
   integer::izap,izat
   real(kr)::x,y,c(*),awp,awt,spi
   ! internals
   integer::ltp,lidp,i2s,nt,np,ip,it
   real(kr)::e,ai,at,zt,zi,cc1,ee,cc2,eta,wn,sigc
   real(kr)::sigr,sigi,sgn
   complex(kr)::carg1,carg2,cs1,cs2
   real(kr)::p(65)
   real(kr),parameter::fm=1.e-12_kr
   real(kr),parameter::zero=0
   real(kr),parameter::uno=1

   !--initialize constants
   e=c(2)
   ltp=nint(c(3))
   lidp=nint(c(4))
   ai=awp*amassn
   at=awt*amassn
   zt=int(izat/1000)
   zi=int(izap/1000)
   i2s=nint(2*spi)
   cc1=2*amu*ev*fm**2/hbar**2
   ee=(ev/10000000)*(clight/10)
   cc2=ee**4*amu/(2*hbar**2*ev)
   eta=zt*zi*sqrt(cc2*ai/e)
   wn=at*sqrt(cc1*ai*e)/(ai+at)
   sigc=0
   if (lidp.eq.0) sigc=(eta**2/wn**2)/(1-x)**2
   if (lidp.eq.1) sigc=((2*eta**2/wn**2)&
     /(1-x**2))*((1+x**2)/(1-x**2)&
     +(-1**i2s)*cos(eta*log((1+x)/(1-x)))/(2*spi+1))
   nt=nint(c(6))
   np=2*nt
   call legndr(x,p,np)

   !--ltp=1: nuclear amplitude expansion
   if (ltp.eq.1) then
      if (lidp.ne.1) then
         sigr=c(7)/2
         do ip=1,np
            sigr=sigr+(2*ip+1)*p(ip+1)*c(ip+7)/2
         enddo
         cs1=dcmplx(c(8+np),c(9+np))/2
         do it=1,nt
            cs1=cs1+(2*it+1)*p(it+1)&
              *dcmplx(c(8+np+2*it),c(9+np+2*it))/2
         enddo
         carg1=dcmplx(zero,uno)*eta*log((1-x)/2)
         sigi=(-2*eta/(1-x))*dble(cs1*exp(carg1))
         y=sigc+sigr+sigi
      else
         sigr=c(7)/2
         do it=1,nt
            sigr=sigr+(4*it+1)*p(2*it+1)*c(it+7)/2
         enddo
         cs1=dcmplx(c(8+nt),c(9+nt))/2
         cs2=cs1
         sgn=-1
         do it=1,nt
            cs1=cs1+(2*it+1)*p(it+1)&
              *dcmplx(c(8+nt+2*it),c(9+nt+2*it))/2
            cs2=cs2+sgn*(2*it+1)*p(it+1)&
              *dcmplx(c(8+nt+2*it),c(9+nt+2*it))/2
            sgn=-sgn
         enddo
         carg1=dcmplx(zero,uno)*eta*log((1-x)/2)
         carg2=dcmplx(zero,uno)*eta*log((1+x)/2)
         sigi=(-2*eta/(1-x*x))*dble(cs1*(1+x)*exp(carg1)&
           +cs2*(1-x)*exp(carg2))
         y=sigc+sigr+sigi
      endif

   !--ltp=2: residual cross section expansion
   else
      if (lidp.ne.1) then
         sigr=c(7)/2
         do ip=1,np
            sigr=sigr+(2*ip+1)*p(ip+1)*c(ip+7)/2
         enddo
         y=sigc+sigr
      else
         sigr=c(7)/2
         do it=1,nt
            sigr=sigr+(4*it+1)*p(2*it+1)*c(it+7)/2
         enddo
         y=sigc+sigr
      endif
   endif
   return
   end subroutine coul

   subroutine acelpp(next,matd,ngmt,nin)
   !-------------------------------------------------------------------
   ! Store detailed photon production data starting at location next.
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   integer::next,matd,ngmt,nin
   ! internals
   integer::nesp,nex,j,nwords,kgmt,mfd,mtd,mtdnc,nb,nw
   integer::nk,mto,ik,ifini,jscr,idone,lf,lp,m,n,nn,nnn,jnt,i
   integer::ie,je,jn,jfirst,jlast,nlast,law,lff,li,ni,ii,mmm
   integer::ne,lc,imu,nexl,nc,ic,nexd,k,lep,nd,na,ncyc,ki
   integer::nyp,mtl,loct,nd0,mtdold
   integer::jp,jpn,jpp
   real(kr)::awr,eg,egamma,ei,en,ep,epu,ef,el,e1,teste,renorm,r
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::dise
   real(kr),dimension(:),allocatable::tdise
   character(66)::strng
   integer,parameter::nwscr=1000000
   integer,parameter::ndise=5000
   real(kr),dimension(:),allocatable::phot
   real(kr),parameter::emev=1.e6_kr
   real(kr),parameter::rmin=1.e-30_kr
   real(kr),parameter::eps=4.e-6_kr
   real(kr),parameter::zero=0

   !--initialize
   mtrp=next
   nesp=nes
   lsigp=mtrp+ntrp
   sigp=lsigp+ntrp
   nex=sigp
   j=0
   nwords=5*ntrp
   allocate(phot(nwords))
   allocate(scr(nwscr))
   allocate(dise(ndise))
   dise=zero

   !--loop over photon production reactions
   kgmt=0
   do while (kgmt.lt.ngmt)
      kgmt=kgmt+1
      mfd=int(gmt(kgmt)/1000)
      mtd=nint(gmt(kgmt)-mfd*1000)
      call findf(matd,mfd,mtd,nin)
      call contio(nin,0,0,scr,nb,nw)
      jp=nint(scr(3))
      jpn=mod(jp,10)
      jpp=(jp-jpn)/10
      nk=n1h
      mto=mtd*10000
      ik=0
      ifini=0
      do while (ifini.eq.0)
         call tab1io(nin,0,0,scr,nb,nw)
         law=l2h
         jscr=1
         idone=0
         do while (idone.eq.0)
            if (nb.eq.0) idone=1
            if (mfd.eq.13.and.(ik.gt.0.or.nk.eq.1)) idone=1
            if (idone.eq.0) then
               if (mfd.ne.13) jscr=jscr+nw
               if (jscr.gt.nwscr)&
                 call error('acelpp','unit error.',' ')
               call moreio(nin,0,0,scr(jscr),nb,nw)
            endif
         enddo
         if (ik.ne.0.or.nk.eq.1) then
            lf=l2h
            j=j+1
            xss(lsigp-1+j)=nex-sigp+1
            xss(nex)=mfd
            nex=nex+1
            if (mfd.eq.12) xss(nex)=mtd
            if (mfd.eq.12) nex=nex+1
            eg=c1h
            en=c2h
            lp=l1h
            m=n1h
            n=n2h
            if (mfd.eq.16) then
               xss(nex)=mtd
               nex=nex+1
               eg=0
               en=0
               lp=0
            endif
            jnt=nint(scr(6+2*m))

            !--file 12
            if (mfd.ne.13) then
               ef=scr(7+2*m)
               el=scr(5+2*m+2*n)
               xss(nex)=m
               if (m.ne.1.or.jnt.ne.2) then
                  do i=1,m
                     xss(nex+i)=scr(5+2*i)
                     xss(i+nex+m)=scr(6+2*i)
                  enddo
                  nex=nex+2*m+1
               else
                  xss(nex)=0
                  nex=nex+1
               endif
               xss(nex)=n
               jscr=5+2*m
               do i=1,n
                  xss(nex+i)=sigfig(scr(jscr+2*i)/emev,7,0)
                  xss(i+n+nex)=scr(2*i+jscr+1)
               enddo
               nex=nex+2*n+1
               if (nex-1.gt.nxss) call error('acelpp',&
                 'insufficient storage for photon production.',' ')

            !--file 13
            else
               e1=scr(7+2*m)
               teste=e1/100000000
               idone=0
               ie=0
               do while (ie.lt.nes.and.idone.eq.0)
                  ie=ie+1
                  if (abs(xss(ie+esz)-e1).lt.teste) idone=1
               enddo
               xss(nex)=ie
               nex=nex+1
               je=0
               jn=0
               jfirst=0
               jlast=0
               jscr=6+2*m
               nlast=(nw-(6+2*m))/2
               idone=0
               i=0
               do while (i.lt.n.and.idone.eq.0)
                  i=i+1
                  if (i.gt.nlast.and.nb.ne.0) then
                     call moreio(nin,0,0,scr,nb,nw)
                     jscr=1-2*nlast-1
                     nlast=nlast+nw/2
                  endif
                  teste=xss(esz+ie+je)/100000000
                  if (abs(scr(jscr+2*i-1)-xss(esz+ie+je)).le.teste) then
                     je=je+1
                     if (scr(jscr+2*i).eq.zero.and.jlast.eq.0.and.&
                       jn.eq.1) jn=0
                     if (jn.eq.0) jfirst=ie+je-1
                     jn=jn+1
                     if (scr(jscr+2*i).ne.zero) jlast=jn
                     xss(jn+nex)=scr(jscr+2*i)
                     if ((ie+je).gt.nesp) idone=1
                  endif
               enddo
               if (jlast.lt.jn) jlast=jlast+1
               n=jlast
               xss(nex-1)=jfirst
               if (ie.ne.jfirst) ie=jfirst
               ef=xss(esz+ie)
               el=xss(esz+ie+n-1)
               xss(nex)=n
               nex=nex+n+1
               if (nex-1.gt.nxss) call error('acelpp',&
                 'insufficient storage for photon production',' ')
            endif
            mto=mto+10
            xss(mtrp-1+j)=mto
            if (lf.eq.1) xss(mtrp-1+j)=mto+1
            if (mfh.eq.16) xss(mtrp-1+j)=mto+2
            phot(5*(j-1)+1)=eg
            phot(5*(j-1)+2)=en
            phot(5*(j-1)+3)=lp
            phot(5*(j-1)+4)=ef
            phot(5*(j-1)+5)=el
         endif
         ik=ik+1
         if (nk.eq.1) ifini=1
         if (ik.gt.nk) ifini=1
      enddo
   enddo
   if (j.ne.ntrp)&
     call error('acelpp','no. of gamma energies not complete.',' ')

   !--store the photon angular distributions
   landp=nex
   nex=nex+ntrp
   andp=nex
   i=0
   mtdold=-1
   do while (i.lt.ntrp)
      i=i+1
      lff=mod(nint(xss(i-1+mtrp)),10)
      if (lff.le.1) then
         mtd=int(xss(i-1+mtrp)/10000)
         if (mtd.ne.mtdold) then
            mtdold=mtd
         else
            call error('acelpp','mf14/mt infinite loop',&
                       'probable endf error')
         endif
         call findf(matd,14,mtd,nin)
         call contio(nin,0,0,scr,nb,nw)
         li=l1h
         nk=n1h
         if (li.ne.0) then
            ! all gammas isotropic for this reaction
            do ik=1,nk
               xss(i+landp-1)=0
               i=i+1
            enddo
            i=i-1
         else
            ! some of the gammas are anisotropic
            ni=n2h
            ik=0
            i=i-1
            if (ni.ne.0) then
               do ii=1,ni
                  ik=ik+1
                  i=i+1
                  call contio(nin,0,0,scr,nb,nw)
                  eg=c1h
                  do j=1,ntrp
                     mmm=nint(xss(mtrp+j-1))
                     mmm=mmm/10000
                     if (mmm.eq.mtd.and.&
                       abs(eg-phot(5*(j-1)+1)).le.eg/10000000)&
                       xss(j+landp-1)=0
                  enddo
               enddo
            endif
            do while (ik.lt.nk)
               ik=ik+1
               i=i+1
               call tab2io(nin,0,0,scr,nb,nw)
               eg=c1h
               ne=n2h
               do j=1,ntrp
                  mmm=nint(xss(mtrp+j-1))
                  mmm=mmm/10000
                  if (mmm.eq.mtd.and.&
                    abs(eg-phot(5*(j-1)+1)).le.eg/10000000)&
                    xss(j+landp-1)=nex-andp+1
               enddo
               xss(nex)=ne
               lc=nex+2*ne
               do ie=1,ne
                  call tab1io(nin,0,0,scr,nb,nw)
                  xss(nex+ie)=sigfig(c2h/emev,7,0)
                  if (n2h.ne.2) then
                     xss(nex+ne+ie)=lc-andp+2
                     do imu=1,33
                        xss(lc+imu)=scr(7+2*imu)
                     enddo
                     lc=lc+33
                  else
                     xss(nex+ne+ie)=0
                  endif
               enddo
               nex=lc+1
            enddo
         endif
      endif
   enddo

   !--store photon energy distributions
   ldlwp=nex
   nex=ldlwp+ntrp
   dlwp=nex
   do i=1,ntrp
      xss(i-1+ldlwp)=nex-dlwp+1
      lff=mod(nint(xss(mtrp-1+i)),10)
      xss(mtrp-1+i)=xss(mtrp-1+i)/10

      !--discrete energies
      !--law=2
      if (lff.eq.0) then
         eg=phot(5*(i-1)+1)
         en=phot(5*(i-1)+2)
         lp=nint(phot(5*(i-1)+3))
         ef=phot(5*(i-1)+4)
         el=phot(5*(i-1)+5)
         xss(nex+1)=2
         xss(nex+2)=nex+9-dlwp+1
         xss(nex+3)=0
         xss(nex+4)=2
         xss(nex+5)=sigfig(ef/emev,7,0)
         xss(nex+6)=sigfig(el/emev,7,0)
         xss(nex+7)=1
         xss(nex+8)=1
         xss(nex+9)=lp
         xss(nex+10)=sigfig(eg/emev,7,0)
         nex=nex+11

      !--data from file 15
      !--law=4
      else if (lff.eq.1) then
         mtd=int(xss(i-1+mtrp)/1000)
         mtdnc=int(xss(i-1+mtrp))-mtd*1000
         xss(nex)=0
         xss(nex+1)=4
         nexl=nex+2
         call findf(matd,15,mtd,nin)
         call contio(nin,0,0,scr,nb,nw)
         nc=n1h
         if (mtdnc.gt.nc) mtdnc=1
         do ic=1,nc
            if (ic.ne.mtdnc) then
            ! dummy read of the unmatched probability record
               call tab1io(nin,0,0,scr,nb,nw)
               do while (nb.ne.0)
                  call moreio(nin,0,0,scr,nb,nw)
               enddo
               call tab2io(nin,0,0,scr,nb,nw)
               ne=n2h
               do ie=1,ne
                  call tab1io(nin,0,0,scr,nb,nw)
                  do while (nb.ne.0)
                     call moreio(nin,0,0,scr,nb,nw)
                  enddo
               enddo
            else
               ! read probability record
               call tab1io(nin,0,0,scr,nb,nw)
               jscr=1
               do while (nb.ne.0)
                  jscr=jscr+nw
                  call moreio(nin,0,0,scr(jscr),nb,nw)
               enddo
               m=n1h
               n=n2h
               jnt=nint(scr(6+2*m))
               if (m.ne.1.or.jnt.ne.2) then
                  xss(nex+3)=m
                  do j=1,m
                     xss(j+3+nex)=scr(5+2*j)
                     xss(j+3+m+nex)=scr(6+2*j)
                  enddo
                  nex=nex+4+2*m
               else
                  xss(nex+3)=0
                  nex=nex+4
               endif
               xss(nex)=n
               do j=1,n
                  xss(j+nex)=sigfig(scr(5+2*j+2*m)/emev,7,0)
                  xss(j+n+nex)=scr(6+2*j+2*m)
               enddo
               nex=nex+2*n+1
               xss(nexl)=nex-dlwp+1
               call tab2io(nin,0,0,scr,nb,nw)
               m=n1h
               jnt=nint(scr(6+2*m))
               if (m.ne.1.or.jnt.ne.2) then
                  xss(nex)=m
                  do j=1,m
                     xss(j+nex)=scr(5+2*j)
                     xss(j+m+nex)=scr(6+2*j)
                  enddo
                  nex=nex+2*m+1
               else
                  xss(nex)=0
                  nex=nex+1
               endif
               ne=n2h
               xss(nex)=ne
               nex=nex+1
               nexd=nex+2*ne
               do ie=1,ne
                  call tab1io(nin,0,0,scr,nb,nw)
                  jscr=1
                  do while (nb.ne.0)
                     jscr=jscr+nw
                     if (jscr.gt.nwscr) call error('acelpp',&
                       'insufficient storage for input photon data.',&
                       ' ')
                     call moreio(nin,0,0,scr(jscr),nb,nw)
                  enddo
                  xss(nex-1+ie)=sigfig(c2h/emev,7,0)
                  xss(nex-1+ne+ie)=nexd-dlwp+1
                  m=n1h
                  n=n2h
                  jnt=nint(scr(6+2*m))
                  xss(nexd)=jnt
                  xss(nexd+1)=n
                  nexd=nexd+1
                  xss(nexd+1+2*n)=0
                  do k=1,n
                     xss(k+nexd)=sigfig(scr(7+2*k)/emev,7,0)
                     xss(k+n+nexd)=sigfig(scr(8+2*k)*emev,7,0)
                     if (xss(k+n+nexd).lt.rmin) xss(k+n+nexd)=0
                     if (k.gt.1.and.jnt.eq.1)&
                        xss(k+2*n+nexd)=xss(k+2*n-1+nexd)&
                        +scr(8+2*(k-1))&
                        *(scr(7+2*k)-scr(7+2*(k-1)))
                     if (k.gt.1.and.jnt.eq.2)&
                        xss(k+2*n+nexd)=xss(k+2*n-1+nexd)&
                        +((scr(8+2*(k-1))+scr(8+2*k))/2)&
                        *(scr(7+2*k)-scr(7+2*(k-1)))
                  enddo
                  ! renormalize cdf to sum to 1, and pdf correspondingly
                  renorm=1/xss(3*n+nexd)
                  do k=1,n
                     xss(k+n+nexd)=sigfig(xss(k+n+nexd)*renorm,9,0)
                     xss(k+2*n+nexd)=sigfig(xss(k+2*n+nexd)*renorm,9,0)
                  enddo
                  nexd=nexd+3*n+1
               enddo
               nex=nex+2*ne+1
            endif
         enddo
         nex=nexd
         if (nex.gt.nxss) call error('acelpp',&
           'insufficient space for cross sections',' ')

      !--data from file 16
      !--law=4
      else
         mtd=int(xss(mtrp-1+i)/1000)
         xss(nex)=0
         xss(nex+1)=4
         xss(nex+2)=nex-dlwp+10
         call findf(matd,16,mtd,nin)

         !--create a union list of all discrete photons for this mtd
         call contio(nin,0,0,scr,nb,nw)
         awr=c2h
         call tab1io(nin,0,0,scr,nb,nw)
         law=nint(scr(4))
         egamma=zero
         if (law.eq.2) egamma=scr(2)
         do while (nb.ne.0)
            call moreio(nin,0,0,scr,nb,nw)
         enddo
         call tab2io(nin,0,0,scr,nb,nw)
         ne=n2h
         nd0=0
         !--loop over incident energies
         do ie=1,ne
            call listio(nin,0,0,scr,nb,nw)
            ei=scr(2)
            nd=nint(scr(3))
            jscr=1
            do while (nb.ne.0)
               jscr=jscr+nw
               call moreio(nin,0,0,scr(jscr),nb,nw)
            enddo
            if (nd.gt.0) then
               !--create an initial list of discrete photons.
               !--use nd0 as a list counter.  if a discrete photon
               !--energy is zero, reset it to a small non-zero
               !--value (for now), but this should not happen
               !--and indicates an error in the evaluated file.
               !--Also do not allow degenerate energy (shouldn't
               !--happen, but has in the past and will cause an
               !--array overflow error later on).
               if (nd0.eq.0) then
                  do ki=1,nd
                     ep=scr(5+2*ki)
                     if (ep.eq.zero) then
                        call mess('acelpp',&
                                 '1discrete photon energy must .ne. 0',&
                                 'reset to 1.e-5 eV')
                        ep=1.e-5_kr
                     endif
                     if (law.eq.2) ep=ep-awr*ei/(awr+1)
                     dise(ki)=ep
                     if (ki.gt.1) then
                        if (dise(ki-1).eq.dise(ki))&
                         dise(ki)=dise(ki)*(1._kr+2*eps)
                     endif
                  enddo
                  nd0=nd
               else
                  !--compare discrete photon list at higher incident
                  !--neutron energies with a union list from lower
                  !--incident neutron energies.  endf formats allow
                  !--these to differ but mcnp doesn't.  let the
                  !--dise array accumulate a union list.
                  do ki=1,nd
                     ep=scr(5+2*ki)
                     if (ep.eq.zero) then
                        call mess('acelpp',&
                                 '2discrete photon energy must .ne. 0',&
                                 'reset to 1.e-5 eV')
                        ep=1.e-5_kr
                     endif
                     if (law.eq.2) ep=ep-awr*ei/(awr+1)
                     if (ki.gt.1.) then
                        if (ep.eq.scr(5+2*(ki-1)))&
                           ep=ep*(1._kr+2*eps)
                     endif
                     do m=1,nd0
                        r=abs(ep/dise(m)-1)
                        if (r.le.eps) go to 111
                     enddo
                     !--found a new discrete energy.  insert it into
                     !--the existing dise array, making sure to
                     !--maintain a highest to lowest energy order.
                     if (abs(scr(5+2*ki)).gt.abs(dise(1))) then
                        do m=nd0,1,-1
                           dise(m+1)=dise(m)
                        enddo
                        dise(1)=ep
                        nd0=nd0+1
                     elseif (abs(ep).lt.abs(dise(nd0))) then
                        dise(nd0+1)=ep
                        nd0=nd0+1
                     else
                        do m=1,nd0-1
                           if (abs(ep).lt.abs(dise(m)).and.&
                               abs(ep).gt.abs(dise(m+1))) then
                              do j=nd0,m+1,-1
                                 dise(j+1)=dise(j)
                              enddo
                              dise(m+1)=ep
                              go to 110
                           endif
                        enddo
                        write(strng,'(''mtd='',i3,''  mt='',i6,&
                          &''  ie='',i4,i5,''  nd='',3i4,'' ed='',&
                          &1p,e12.5)')&
                          mtd,int(xss(mtrp-1+i)),ie,ne,nd,nd0,ki,&
                          scr(5+2*ki)
                        call error('acelpp',&
                          'a discrete energy was not consistent.',&
                          strng)
  110                   continue
                        nd0=nd0+1
                     endif
  111                continue
                  enddo
               endif
               if (nd0.gt.ndise) then
                  write(strng,'(''nwords is'',i6,'' but need '',i6)')&
                        nwords,nd0
                  call error('acelpp',&
                             'too many discrete photons found',strng)
               endif
            endif
         enddo
         !--union list of discrete photons is complete.

         call repoz(nin)
         call findf(matd,16,mtd,nin)
         call contio(nin,0,0,scr,nb,nw)
         jscr=1
         call tab1io(nin,0,0,scr(jscr),nb,nw)
         do while (nb.ne.0)
            jscr=jscr+nw
            call moreio(nin,0,0,scr(jscr),nb,nw)
         enddo
         m=n1h
         n=n2h
         xss(nex+3)=0
         xss(nex+4)=2
         xss(nex+5)=sigfig(scr(7+2*m)/emev,7,0)
         xss(nex+6)=sigfig(scr(5+2*m+2*n)/emev,7,0)
         xss(nex+7)=1
         xss(nex+8)=1
         nex=nex+9
         call tab2io(nin,0,0,scr,nb,nw)
         lep=l2h
         m=n1h
         jnt=nint(scr(6+2*m))
         if (m.ne.1.or.jnt.ne.2) then
            xss(nex)=m
            do j=1,m
               xss(nex+j)=scr(5+2*j)
               xss(nex+m+j)=scr(6+2*j)
            enddo
            nex=nex+2*m+1
         else
            xss(nex)=0
            nex=nex+1
         endif
         ne=n2h
         xss(nex)=ne
         nex=nex+1
         nexd=nex+2*ne
         !--loop over incident energies
         if (nd0.gt.0) then
            allocate(tdise(2*ndise))
            tdise=0
         endif
         do ie=1,ne
            call listio(nin,0,0,scr,nb,nw)
            ei=scr(2)
            nd=nint(scr(3))
            na=nint(scr(4))
            ncyc=na+2
            jscr=1
            do while (nb.ne.0)
               jscr=jscr+nw
               call moreio(nin,0,0,scr(jscr),nb,nw)
            enddo
            xss(nex-1+ie)=sigfig(c2h/emev,7,0)
            xss(nex-1+ne+ie)=nexd-dlwp+1
            n=n2h

            !--if discrete photons are present make sure all energies
            !--are non-zero, then save a copy of these energies and
            !--their probabilities.  Increment any degenerate energy
            !--as was done before.
            if (nd.ne.0) then
               do nn=1,nd
                  if (scr(5+2*nn).eq.zero) scr(5+2*nn)=1.e-5_kr
                  tdise(2*(nn-1)+1)=scr(5+2*nn)
                  if (nn.gt.1) then
                     if (tdise(2*(nn-1)+1).eq.tdise(2*(nn-2)+1))&
                        tdise(2*(nn-1)+1)=tdise(2*(nn-1)+1)*(1._kr+2*eps)
                  endif
                  tdise(2*nn)=scr(6+2*nn)
               enddo
            endif

            !--make sure all (nd0) discrete photons are included for
            !--all incident energies.  nd is the number of discrete
            !--photons for the current incident energy.  when nd=nd0,
            !--only need to check law and whether this is a primary
            !--photon.
            if (nd0.ne.0.and.nd.eq.nd0) then
               do nn=1,nd
                  if (law.eq.1.and.scr(5+2*nn).lt.zero)&
                               scr(5+2*nn)=-scr(5+2*nn)+ei*awr/(awr+1)
               enddo
            elseif (nd0.ne.0.and.nd.ne.nd0) then
               !--if nd=0 then must insert all discrete photons
               !--with zero probability
               if (nd.eq.0) then
                  !--move continuous data, then insert discrete data
                  do ki=n*2,1,-1
                     scr(6+2*nd0+ki)=scr(6+ki)
                  enddo
                  do ki=1,nd0
                     epu=dise(ki)
                     if (law.eq.1.and.epu.lt.zero)&
                                      epu=-epu+ei*awr/(awr+1)
                     if (law.eq.2) epu=epu+ei*awr/(awr+1)
                     scr(5+2*ki)=epu
                     scr(6+2*ki)=0
                  enddo
                  !--update length of list record
                  n=n+nd0

               !--already have some, but not all discrete photons
               !--tabulated.  compare current list versus union
               !--and insert missing data
               else
                  !--move continuous data, then insert discrete data
                  do ki=n*2,nd*2+1,-1
                     scr(6+2*nd0+ki-2*nd)=scr(6+ki)
                  enddo
                  !--loop over union list of photons, inserting
                  !--missing energies with zero probability.  also
                  !--check law and/or sign of photon energy to know
                  !--if this is a primary or secondary photon.  if
                  !--a primary photon, its energy must be increased
                  !--to account for the incident neutron energy.
                  do m=nd0,1,-1
                     ep=dise(m)
                     if (law.eq.1.and.ep.lt.zero)ep=-ep+ei*awr/(awr+1)
                     if (law.eq.2)ep=ep+ei*awr/(awr+1)
                     scr(5+2*m)=ep
                     scr(6+2*m)=zero
                     nn=nd+1
                     do while (nn.gt.1)
                        nn=nn-1
                        r=(abs(tdise(2*(nn-1)+1))/abs(dise(m)))-1
                        if (abs(r).lt.0.001*eps) then
                           scr(6+2*m)=tdise(2*nn)
                           nn=-1
                        elseif (r.gt.zero) then
                           nn=-1
                        endif
                     enddo
                  enddo
                  !--update length of list record
                  n=n-nd+nd0
               endif
               !--update number of discrete photons
               nd=nd0
            endif
            !--done with discrete photon corrections

            xss(nexd)=lep+10*nd
            xss(nexd+1)=n
            nexd=nexd+1
            !--cdf starts at zero, unless there are discrete
            !  photons in which case the initial cdf equals
            !  that first photon's pdf.
            xss(nexd+1+2*n)=0
            do ki=1,n
               xss(ki+nexd)=sigfig(scr(7+ncyc*(ki-1))/emev,7,0)
               if (ki.le.nd) xss(ki+n+nexd)=scr(8+ncyc*(ki-1))
               if (ki.gt.nd) xss(ki+n+nexd)=&
                 scr(8+ncyc*(ki-1))*emev
               if (xss(ki+n+nexd).lt.rmin) xss(ki+n+nexd)=0
               if (nd.gt.0.and.ki.eq.1) then
                  xss(ki+2*n+nexd)=xss(ki+n+nexd)
               elseif (ki.le.nd) then
                  xss(ki+2*n+nexd)=xss(ki+2*n-1+nexd)+scr(8+ncyc*(ki-1))
               endif
               if (nd.gt.0.and.ki.eq.nd+1) xss(ki+2*n+nexd)=&
                 xss(ki+2*n-1+nexd)
               if (ki.gt.nd+1.and.lep.eq.1) xss(ki+2*n+nexd)=&
                 xss(ki+2*n-1+nexd)+scr(8+ncyc*(ki-2))&
                 *(scr(7+ncyc*(ki-1))-scr(7+ncyc*(ki-2)))
               if (ki.gt.nd+1.and.lep.eq.2) xss(ki+2*n+nexd)=&
                 xss(ki+2*n-1+nexd)+((scr(8+ncyc*(ki-2))&
                 +scr(8+ncyc*(ki-1)))/2)&
                 *(scr(7+ncyc*(ki-1))-scr(7+ncyc*(ki-2)))
            enddo
            renorm=1
            if (xss(3*n+nexd).ne.0.) renorm=1/xss(3*n+nexd)
            do ki=1,n
               xss(ki+n+nexd)=sigfig(xss(ki+n+nexd)*renorm,9,0)
               xss(ki+2*n+nexd)=sigfig(xss(ki+2*n+nexd)*renorm,9,0)
            enddo
            nexd=nexd+3*n+1
         enddo
         if (allocated(tdise)) deallocate(tdise)
         nex=nexd
         if (nex.gt.nxss) call error('acelpp',&
           'insufficient space for cross sections',' ')
      endif
   enddo
   deallocate(phot)
   deallocate(dise)
   yp=nex
   nex=nex+1
   nyp=0
   mtl=0
   do i=1,ntrp
      mtd=int(xss(i-1+mtrp)/1000)
      loct=nint(xss(i-1+lsigp)+sigp-1)
      mfd=nint(xss(loct))
      if (mfd.ne.13.and.mtd.ne.mtl) then
         mtl=mtd
         nyp=nyp+1
         xss(nex)=mtd
         nex=nex+1
      endif
   enddo
   xss(yp)=nyp
   next=nex
   if (next.gt.nxss)call error('acelpp',&
     'insufficient space for cross sections.',' ')
   deallocate(scr)
   return
   end subroutine acelpp

   subroutine fndar1(akal,rkal,fbarc,delfc,ee,ep)
   !-------------------------------------------------------------------
   ! This subroutine finds the Kalbach angular distribution parameters
   ! akal and rkal from the average of the forward and backward
   ! amplitudes fbarc, and the difference of the forward and backward
   ! amplitudes delfc, both defined in the cm frame.
   ! Writen by A. J. Sierk, LANL, 1 March 1990.
   !-------------------------------------------------------------------
   use util ! provides mess
   ! externals
   real(kr)::akal,rkal,fbarc,delfc,ee,ep
   ! internals
   integer::iq1,idone
   real(kr)::fofa,s2a,c2a,dela
   character(60)::strng
   real(kr),parameter::test=1.e-10_kr
   real(kr),parameter::small=1.e-5_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1

   iq1=0
   idone=0
   akal=2*fbarc*(1-4*fbarc**2/3)

   !--use newton's method to converge to a solution of fofa=0,
   !--where fofa=akal/tanh(akal)-2*fbarc
   do while (iq1.lt.30.and.idone.eq.0)
      fofa=akal/tanh(akal)-2*fbarc
      if (abs(fofa).lt.test) then
         idone=1
      else
         s2a=sinh(2*akal)
         c2a=cosh(2*akal)
         dela=(2*fbarc*(c2a-1)-akal*s2a)/(s2a-2*akal)
         akal=akal+dela
         if (akal.lt.zero) akal=-akal
         iq1=iq1+1
      endif
   enddo
   if (idone.eq.0) then
      write(strng,'('' e='',1p,e12.4,''  eprime='',e12.4)') ee,ep
      call mess('fndar1','loop to find kalbach a not converged',&
        strng)
      akal=small
      rkal=0
   else
      if (akal.lt.zero) then
         write(strng,'('' e='',1p,e12.4,''  eprime='',e12.4)') ee,ep
         call mess('fndar1','kalbach a is negative',strng)
         akal=small
         rkal=0
      else
         rkal=0
         if (akal.ne.zero) rkal=delfc/akal
         if (rkal.lt.zero.or.rkal.gt.one) then
            write(strng,'('' e='',1p,e12.4,''  eprime='',e12.4)') ee,ep
            call mess('fndar1','kalbach r is unreasonable',strng)
            if (rkal.lt.zero) rkal=0
            if (rkal.gt.one) rkal=1
         endif
         if (akal.lt.small) akal=small
      endif
   endif
   return
   end subroutine fndar1

   subroutine fndar2(akal,rkal,fbcm,ffcm,emu1,emu2,ee,ep)
   !-------------------------------------------------------------------
   ! This subroutine finds the Kalbach angular distribution parameters
   ! akal and rkal from tabulated angular distributions using the lang
   ! =11--13 options.  This subroutine is only used when the table
   ! does not cover the entire range (-1,1) in cos(theta-cm).
   !-------------------------------------------------------------------
   use util ! provides mess
   ! externals
   real(kr)::akal,rkal,fbcm,ffcm,emu1,emu2,ee,ep
   ! internals
   integer::iq1,idone
   real(kr)::fb,dmu,cam,cap,sam,sap,sa,ca,s2a,cad,sad,fofa
   real(kr)::denom,dela
   character(60)::strng
   real(kr),parameter::test=1.e-10_kr
   real(kr),parameter::small=1.e-5_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1

   iq1=0
   idone=0
   fb=fbcm+ffcm
   dmu=emu2-emu1
   akal=fb*(1-fb**2/3)
   do while (iq1.lt.30.and.idone.eq.0)
      cam=cosh(akal*emu1)
      cap=cosh(akal*emu2)
      sam=sinh(akal*emu1)
      sap=sinh(akal*emu2)
      sa=sinh(akal)
      ca=cosh(akal)
      s2a=sinh(2*akal)
      cad=cosh(akal*dmu)
      sad=sinh(akal*dmu)
      fofa=ffcm*sam-fbcm*sap+akal*sad/(2*sa)
      if (abs(fofa).le.test) then
         idone=1
      else
         denom=sad*(sa-akal*ca)+akal*dmu*cad*sa+&
           2*sa*sa*(ffcm*emu1*cam-fbcm*emu2*cap)
         dela=(2*(fbcm*sap-ffcm*sam)*sa-akal*sad)*sa/denom
         akal=akal+dela
      endif
   enddo
   if (idone.eq.0) then
      write(strng,'('' e='',1p,e12.4,''  eprime='',e12.4)') ee,ep
      call mess('fndar2','loop to find kalbach a not converged',&
        strng)
      akal=small
      rkal=0
   else
      if (akal.lt.zero) then
         write(strng,'('' e='',1p,e12.4,''  eprime='',e12.4)') ee,ep
         call mess('fndar2','kalbach a is negative',strng)
         akal=small
         rkal=0
      else
         rkal=0
         if (akal.ne.zero) rkal=2*sa*(ffcm*cam-fbcm*cap)/(akal*s2a)
         if (rkal.lt.zero.or.rkal.gt.one) then
            write(strng,'('' e='',1p,e12.4,''  eprime='',e12.4)')&
              ee,ep
            call mess('fndar2','kalbach r is unreasonable',strng)
            if (rkal.lt.zero) rkal=0
            if (rkal.gt.one) rkal=1
         endif
         if (akal.lt.small) akal=small
      endif
   endif
   return
   end subroutine fndar2

   real(kr) function fmn(e2,efl,efh,tm)
   !-------------------------------------------------------------------
   ! Compute value of Madland-Nix fission spectrum.
   !-------------------------------------------------------------------
   use mathm ! provides e1
   ! externals
   real(kr)::e2,efl,efh,tm
   ! internals
   real(kr)::u1,u2,g1,g2
   real(kr),parameter::thrhaf=1.5e0_kr

   u1=(sqrt(e2)-sqrt(efl))**2/tm
   u2=(sqrt(e2)+sqrt(efl))**2/tm
   g1=(u2**thrhaf*e1(u2)-u1**thrhaf*e1(u1)+gami(thrhaf,u2)&
     -gami(thrhaf,u1))/(3*sqrt(efl*tm))
   u1=(sqrt(e2)-sqrt(efh))**2/tm
   u2=(sqrt(e2)+sqrt(efh))**2/tm
   g2=(u2**thrhaf*e1(u2)-u1**thrhaf*e1(u1)+gami(thrhaf,u2)&
     -gami(thrhaf,u1))/(3*sqrt(efh*tm))
   fmn=(g1+g2)/2
   return
   end function fmn

   subroutine acelcp(next,matd,nin,za,awr)
   !-------------------------------------------------------------------
   ! Store particle production data.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util ! provides error, mess
   use endf ! provides endf routines and variables
   use physics ! provides amassn,amu,ev,clight
   use acecm ! provides ptleg2,pttab2,eavl
   integer hpd,tyrh,lsigh,sigh,landh,andh,ldlwh,dlwh,yh
   ! externals
   integer::next,matd,nin
   real(kr)::za,awr
   ! internals
   character(60)::strng
   integer::iza,i,ip,j,itype,ipp,it,ie,jp,mt,mf,mtt,ir,k,n
   integer::iaa,nb,nw,lct,nk,ik,izap,law,ll,nrr,npp,idis
   integer::nrint,ne,ipj,loce,leee,llht,int,naa,ltt,il
   integer::iie,m,lly,izarec,lang,nl,iil,iint,nn,kk
   integer::na,nc,nmu,nx,intx,ix,imu,iskip,last,idone
   integer::lawnow,lep,lee,lle,llh,nd,ncyc,ng,nexcd,ig
   integer::ia,lll,intmu,iep,llad,mus,nra,npa,npep,intep
   integer::ki,np,npsx,lld,isocp,llx
   integer::ipt,mtrh,ntrh,ipn
   real(kr)::emc2,thresh,tt,e,y,en,ss,q,amass,ubar,sum,renorm
   real(kr)::h,awp,aprime,th,r1,r2,betasq,awprec,ee
   real(kr)::avadd,avlab,avll,test,rkal,ep,akal,del,av
   real(kr)::eavi,avl,avcm,sign,dele,avav,tt1,tt2
   real(kr)::apsx,step1,step2,xx,test1,test2,test3
   real(kr)::xl,pl,yn,pn,rn,ecm,ea,eim,chek,dmu,e1l,e2l
   real(kr)::chk,pp,disc,v1,v2,e1,e2,pp1,pp2,de,pp1l,pp2l
   real(kr)::eav,suml,chkl,ebar,ad,amuu,amulst,chklst,heat,g
   real(kr)::epl,gl,gammsq,amun,eavlst
   real(kr),dimension(:),allocatable::scr
   integer,parameter::nwscr=10000000
   real(kr),parameter::awr1=.99862e0_kr
   real(kr),parameter::awr2=1.99626e0_kr
   real(kr),parameter::awr3=2.98960e0_kr
   real(kr),parameter::awr4=2.98903e0_kr
   real(kr),parameter::awr5=3.96713e0_kr
   real(kr),parameter::small=1.e-30_kr
   real(kr),parameter::emev=1.e6_kr
   real(kr),parameter::etop=1.e10_kr
   real(kr),parameter::elow=1.e-5_kr
   real(kr),parameter::delt=1.e-10_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1
   real(kr),parameter::ten=10

   !--allocate scratch storage
   allocate(scr(nwscr))

   iza=nint(za)
   emc2=amassn*amu*clight*clight/ev/emev

   !--store particle production data
   ptype=end+1
   ntype=0
   if (t201.lt.etop.and.izai.ne.1) then
      ntype=ntype+1
      xss(ptype+ntype-1)=1
   endif
   if (t203.lt.etop.and.izai.ne.1001) then
      ntype=ntype+1
      xss(ptype+ntype-1)=9
   endif
   if (t204.lt.etop.and.izai.ne.1002) then
      ntype=ntype+1
      xss(ptype+ntype-1)=31
   endif
   if (t205.lt.etop.and.izai.ne.1003) then
      ntype=ntype+1
      xss(ptype+ntype-1)=32
   endif
   if (t206.lt.etop.and.izai.ne.2003) then
      ntype=ntype+1
      xss(ptype+ntype-1)=33
   endif
   if (t207.lt.etop.and.izai.ne.2004) then
      ntype=ntype+1
      xss(ptype+ntype-1)=34
   endif

   !--exit if none are found
   if (ntype.eq.0) then
      ptype=0
      ntro=0
      ploct=0
      deallocate(scr)
      return
   endif

   !--count up productions
   ntro=ptype+ntype
   ploct=ntro+ntype
   do i=1,ntype
      if (nint(xss(ptype+i-1)).eq.1) ip=1
      if (nint(xss(ptype+i-1)).eq.9) ip=1001
      if (nint(xss(ptype+i-1)).eq.31) ip=1002
      if (nint(xss(ptype+i-1)).eq.32) ip=1003
      if (nint(xss(ptype+i-1)).eq.33) ip=2003
      if (nint(xss(ptype+i-1)).eq.34) ip=2004
      ntrh=0
      do j=1,nprod
         if (iprod(j).eq.ip.and.kprod(j).gt.0) ntrh=ntrh+1
      enddo
      xss(ntro+i-1)=ntrh
   enddo
   next=ploct+10*ntype

   !--loop over each of the ntype productions to build data
   do itype=1,ntype
      ipt=nint(xss(ptype+itype-1))
      ntrh=nint(xss(ntro+itype-1))
      if (ipt.eq.1) ip=1
      if (ipt.eq.9) ip=1001
      if (ipt.eq.31) ip=1002
      if (ipt.eq.32) ip=1003
      if (ipt.eq.33) ip=2003
      if (ipt.eq.34) ip=2004

      !--determine the threshold
      !--and assign some locators
      if (ip.eq.1) thresh=t201/emev
      if (ip.eq.1001) thresh=t203/emev
      if (ip.eq.1002) thresh=t204/emev
      if (ip.eq.1003) thresh=t205/emev
      if (ip.eq.2003) thresh=t206/emev
      if (ip.eq.2004) thresh=t207/emev
      ipp=ip
      it=1
      do while (xss(esz+it-1).lt.thresh*(1-delt))
         it=it+1
      enddo
      hpd=next
      xss(ploct+10*(itype-1))=hpd
      xss(hpd)=it
      xss(hpd+1)=nes-it+1
      mtrh=hpd+2+2*(nes-it+1)
      xss(ploct+10*(itype-1)+1)=mtrh
      tyrh=mtrh+ntrh
      xss(ploct+10*(itype-1)+2)=tyrh
      lsigh=tyrh+ntrh
      xss(ploct+10*(itype-1)+3)=lsigh
      sigh=lsigh+ntrh
      xss(ploct+10*(itype-1)+4)=sigh
      next=sigh
      do ie=it,nes
         xss(hpd+2+ie-it)=0
         xss(hpd+2+(nes-it+1)+ie-it)=0
      enddo

      !--find each mt and subsection that contributes
      !--to this production
      jp=0
      do j=1,nprod
         ipn=iprod(j)
         mt=mprod(j)
         mf=kprod(j)
         if (ipn.eq.ip.and.mf.ne.0) then
            jp=jp+1
            xss(mtrh+jp-1)=mt
            xss(lsigh+jp-1)=next-sigh+1
            mtt=0
            ir=0
            do while (mtt.ne.mt)
               ir=ir+1
               mtt=nint(xss(mtr+ir-1))
               k=nint(xss(lsig+ir-1))+sig-1
               n=nint(xss(k+1))
               iaa=nint(xss(k))
            enddo
            if (mt.eq.2) iaa=1

            !--special branch for mf4, neutron elastic recoil,
            !--or neutron capture recoil
            if ((mt.eq.102.and.izai.eq.1).or.mf.eq.4) then
               xss(tyrh+jp-1)=-1
               if (mt.eq.102) xss(tyrh+jp-1)=1
               do ie=iaa,nes
                  if (mt.eq.2) then
                     xss(hpd+2+ie-it)=&
                       sigfig(xss(esz+3*nes+ie-it),7,0)
                  else
                     tt=xss(hpd+2+ie-it)+xss(2+k+ie-iaa)
                     xss(hpd+2+ie-it)=sigfig(tt,7,0)
                  endif
               enddo
               xss(next)=12
               xss(next+1)=mt
               xss(next+2)=0
               xss(next+3)=2
               xss(next+4)=sigfig(xss(esz+it-1),7,0)
               xss(next+6)=1
               xss(next+5)=sigfig(xss(esz+nes-1),7,0)
               xss(next+7)=1
               next=next+8

            !--normal branch for other reactions (mf=6)
            else

               !--first check the subsection to see whether
               !--the distribution is isotropic or not.
               isocp=1
               call findf(matd,mf,mt,nin)
               call contio(nin,0,0,scr,nb,nw)
               nk=n1h
               ik=0
               idone=0
               do while (ik.lt.nk.and.idone.eq.0)
                  ik=ik+1
                  ll=1
                  lly=ll
                  call tab1io(nin,0,0,scr(ll),nb,nw)
                  izap=nint(c1h)
                  awp=c2h
                  law=l2h
                  ll=ll+nw
                  do while (nb.ne.0)
                     call moreio(nin,0,0,scr(ll),nb,nw)
                     ll=ll+nw
                  enddo

                  !--if not the desired particle,
                  !--or not a law=1 subsection,
                  !--skip the subsection
                  if (izap.ne.ip.or.law.ne.1) then
                     call skip6a(nin,0,0,scr,law)
                  else
                     call tab2io(nin,0,0,scr(ll),nb,nw)
                     lang=nint(scr(ll+2))
                     lep=nint(scr(ll+3))
                     ne=nint(scr(ll+5))
                     do ie=1,ne
                        ll=1
                        call listio(nin,0,0,scr(ll),nb,nw)
                        ll=ll+nw
                        do while (nb.ne.0)
                           call moreio(nin,0,0,scr(ll),nb,nw)
                           ll=ll+nw
                        enddo
                        na=nint(scr(4))
                        if (na.gt.0) isocp=0
                     enddo
                  endif
               enddo

               !--go back and process the subsection

               call findf(matd,mf,mt,nin)
               call contio(nin,0,0,scr,nb,nw)
               lct=l2h
               if (lct.eq.1) then
                  xss(tyrh+jp-1)=1
               else
                  xss(tyrh+jp-1)=-1
               endif
               nk=n1h

               !--loop over subsections
               do ik=1,nk
                  call tab1io(nin,0,0,scr,nb,nw)
                  izap=nint(c1h)
                  law=l2h
                  ll=1+nw
                  do while (nb.ne.0)
                     call moreio(nin,0,0,scr(ll),nb,nw)
                     ll=ll+nw
                     if (ll.gt.nwscr) call error('acelcp',&
                        'exceeded scratch storage',' ')
                  enddo
                  if (mt.eq.2.and.nk.eq.1) izap=nint(za)

                  !--accumulate contribution to production and
                  !--store the yield from this subsection
                  if (ik.eq.lprod(j)) then
                     nrr=1
                     npp=2
                     do ie=iaa,nes
                        e=xss(esz+ie-1)*emev
                        call terpa(y,e,en,idis,scr,npp,nrr)
                        if (y.lt.delt) y=0
                        e=e/emev
                        if (mth.eq.2) then
                           ss=xss(esz+3*nes+ie-1)
                        else
                           ss=xss(2+k+ie-iaa)
                        endif
                        tt=xss(hpd+2+ie-it)+y*ss
                        xss(hpd+2+ie-it)=sigfig(tt,7,0)
                     enddo

                     !--store the yield
                     xss(next)=12
                     xss(next+1)=mt
                     next=next+2
                     nrint=nint(scr(5))
                     if (nrint.eq.1.and.nint(scr(8)).eq.2) then
                        xss(next)=0
                     else
                        xss(next)=nrint
                        do i=1,nrint
                           xss(next+i)=nint(scr(5+2*i))
                           xss(next+nrint+i)=nint(scr(6+2*i))
                        enddo
                        next=next+2*nrint
                     endif
                     next=next+1
                     ne=nint(scr(6))
                     xss(next)=ne
                     do i=1,ne
                        xss(next+i)=&
                          sigfig(scr(5+2*nrint+2*i)/emev,7,0)
                        xss(next+i+ne)=&
                          sigfig(scr(6+2*nrint+2*i),7,0)
                     enddo
                     next=next+1+2*ne
                  endif

                  !--skip to the next subsection
                  call skip6a(nin,0,0,scr,law)
               enddo
            endif
         endif

      !--continue searching for the mts that contribute
      enddo

      !--now go back and look for angular distribution data
      landh=next
      xss(ploct+10*(itype-1)+5)=landh
      andh=landh+ntrh
      xss(ploct+10*(itype-1)+6)=andh
      next=andh
      jp=0

      !--loop over production types
      do j=1,nprod
         ipj=iprod(j)
         mt=mprod(j)
         mf=kprod(j)
         if (ipj.eq.ip.and.mf.ne.0) then
            mtt=0
            ir=0
            do while (mtt.ne.mt)
               ir=ir+1
               mtt=nint(xss(mtr+ir-1))
               k=nint(xss(lsig+ir-1))+sig-1
               n=nint(xss(k+1))
               iaa=nint(xss(k))
               q=xss(lqr+ir-1)
            enddo
            jp=jp+1
            if (mt.eq.2) iaa=1

            !--angular representation for elastic recoil particle
            !--just reverse the distribution of scattered particle
            if (mt.eq.2) then
               xss(landh+jp-1)=next-andh+1
               loce=nint(xss(land))+and-1
               ne=nint(xss(loce))
               xss(next)=ne
               leee=next
               do ie=1,ne
                  xss(next+ie)=xss(loce+ie)
               enddo
               nb=next+ne
               next=next+1+2*ne
               loce=loce+1+2*ne
               llht=1
               scr(llht)=0
               scr(llht+1)=0
               scr(llht+2)=0
               scr(llht+3)=0
               scr(llht+4)=1
               scr(llht+5)=ne
               scr(llht+6)=ne
               scr(llht+7)=5
               amass=awr/awi
               do ie=1,ne
                  int=nint(xss(loce))
                  n=nint(xss(loce+1))
                  xss(nb+ie)=-(next-andh+1)
                  xss(next)=int
                  xss(next+1)=n
                  do i=1,n
                     xss(next+1+i)=-xss(loce+2+n-i)
                     xss(next+1+n+i)=xss(loce+2+2*n-i)
                     if (xss(next+1+n+i).lt.small) xss(next+1+n+i)=0
                     if (i.eq.1) then
                        xss(next+1+2*n+1)=0
                        ubar=0
                     endif
                     if (i.gt.1.and.int.eq.1) then
                        sum=xss(next+1+2*n+i-1)&
                          +xss(next+1+n+i-1)*(xss(next+1+i)&
                          -xss(next+1+i-1))
                        xss(next+1+2*n+i)=sigfig(sum,7,0)
                        ubar=ubar&
                          +xss(next+1+n+i-1)*(xss(next+1+i)&
                          -xss(next+1+i-1))&
                          +(xss(next+1+i)+xss(next+1+i-1))/2
                     endif
                     if (i.gt.1.and.int.eq.2) then
                        sum=xss(next+1+2*n+i-1)&
                          +(xss(next+1+n+i)+xss(next+1+n+i-1))&
                          *(xss(next+1+i)-xss(next+1+i-1))/2
                        xss(next+1+2*n+i)=sigfig(sum,7,0)
                        ubar=ubar&
                          +(xss(next+1+n+i)+xss(next+1+n+i-1))&
                         *(xss(next+1+i)-xss(next+1+i-1))&
                         *(xss(next+1+i)+xss(next+1+i-1))/4
                     endif
                  enddo
                  renorm=1/xss(next+1+3*n)
                  do i=1,n
                     xss(next+1+n+i)= &
                       sigfig(renorm*xss(next+1+n+i),7,0)
                     xss(next+1+2*n+i)=&
                        sigfig(renorm*xss(next+1+2*n+i),9,0)
                  enddo
                  next=next+2+3*n
                  loce=loce+2+3*n
                  e=xss(leee+ie)
                  scr(llht+6+2*ie)=e
                  scr(llht+7+2*ie)=2*amass*e*(1+ubar)/(1+amass)**2
               enddo
               ! add in contribution to heating
               naa=nint(xss(hpd+1))
               nrr=1
               npp=2
               do ie=it,nes
                  e=xss(esz+ie-1)
                  call terpa(h,e,en,idis,scr(llht),npp,nrr)
                  xss(hpd+2+naa+ie-it)=xss(hpd+2+naa+ie-it)&
                    +h*xss(hpd+2+ie-it)
               enddo

            !--angular representation for neutron capture recoil
            !--just construct a straight-ahead distribution
            else if (mt.eq.102.and.izai.eq.1) then
               xss(landh+jp-1)=next-andh+1
               xss(next)=2
               xss(next+1)=xss(esz+iaa-1)
               xss(next+2)=xss(esz+nes-1)
               xss(next+3)=-(next-andh+6)
               xss(next+4)=-(next-andh+17)
               next=next+5
               do ie=1,2
                  xss(next)=2
                  xss(next+1)=3
                  xss(next+2)=-1
                     xss(next+5)=0
                  xss(next+8)=0
                  xss(next+3)=one-one/100
                  xss(next+6)=0
                  xss(next+9)=0
                  xss(next+4)=1
                  xss(next+7)=200
                  xss(next+10)=1
                  next=next+11
               enddo

            !--particle distributions using mf=4
            else if (mf.eq.4) then
               call findf(matd,mf,mt,nin)
               call contio(nin,0,0,scr,nb,nw)
               ltt=l2h
               lld=1
               awp=1
               if (mt.ge.600.and.mt.le.649) awp=awr1
               if (mt.ge.650.and.mt.le.699) awp=awr2
               if (mt.ge.700.and.mt.le.749) awp=awr3
               if (mt.ge.750.and.mt.le.799) awp=awr4
               if (mt.ge.800.and.mt.le.849) awp=awr5

               !--legendre distributions
               if (ltt.eq.1) then
                  call contio(nin,0,0,scr,nb,nw)
                  lct=l2h
                  call tab2io(nin,0,0,scr,nb,nw)
                  ne=n2h
                  xss(landh+jp-1)=next-andh+1
                  xss(next)=ne
                  ie=next
                  il=ie+ne
                  next=il+ne+1
                  llht=1
                  lld=llht+8+2*ne
                  scr(llht)=0
                  scr(llht+1)=0
                  scr(llht+2)=0
                  scr(llht+3)=0
                  scr(llht+4)=1
                  scr(llht+5)=ne
                  scr(llht+6)=ne
                  scr(llht+7)=2
                  amass=awr/awi
                  aprime=awp/awi
                  do iie=1,ne
                     ll=lld
                     call listio(nin,0,0,scr(ll),nb,nw)
                     call ptleg2(scr(lld))
                     xss(ie+iie)=sigfig(scr(lld+1)/emev,7,0)
                     m=nint(scr(lld+4))
                     n=nint(scr(lld+5))
                     xss(il+iie)=next-andh+1
                     xss(il+iie)=-xss(il+iie)
                     int=nint(scr(lld+7))
                     xss(next)=int
                     xss(next+1)=n
                     if (next+2+3*n.gt.nxss) call error('acelcp',&
                  'insufficient storage for angular distributions.',&
                     ' ')
                     do i=1,n
                        xss(next+1+i)=&
                          sigfig(scr(lld+4+2*m+2*i),7,0)
                        xss(next+1+n+i)=&
                          sigfig(scr(lld+5+2*m+2*i),7,0)
                        if (i.eq.1) then
                           xss(next+1+2*n+i)=0
                           ubar=0
                        endif
                        if (i.gt.1.and.int.eq.1) then
                           sum=xss(next+1+2*n+i-1)&
                             +xss(next+1+n+i-1)&
                             *(xss(next+1+i)-xss(next+1+i-1))
                           xss(next+1+2*n+i)=sigfig(sum,7,0)
                           ubar=ubar&
                             +xss(next+1+n+i-1)&
                             *(xss(next+1+i)-xss(next+1+i-1))&
                             +(xss(next+1+i)+xss(next+1+i-1))/2
                        endif
                        if (i.gt.1.and.int.eq.2) then
                           sum=xss(next+1+2*n+i-1)&
                            +(xss(next+1+n+i)+xss(next+1+n+i-1))&
                            *(xss(next+1+i)-xss(next+1+i-1))/2
                           xss(next+1+2*n+i)=sigfig(sum,7,0)
                           ubar=ubar&
                             +(xss(next+1+n+i)+xss(next+1+n+i-1))&
                             *(xss(next+1+i)-xss(next+1+i-1))&
                             *(xss(next+1+i)+xss(next+1+i-1))/4
                        endif
                     enddo
                     renorm=1/xss(next+1+3*n)
                     do i=1,n
                        xss(next+1+n+i)=&
                          sigfig(renorm*xss(next+1+n+i),7,0)
                        xss(next+1+2*n+i)=&
                          sigfig(renorm*xss(next+1+2*n+i),9,0)
                     enddo
                     next=next+2+3*n
                     e=xss(ie+iie)
                     th=(1+amass)*q/amass
                     r1=amass*(amass+1-aprime)/aprime
                     r2=aprime/(1+amass)**2
                     betasq=r1*(1+th/e)
                     if (betasq.lt.zero) betasq=0
                     scr(llht+6+2*iie)=e
                     scr(llht+7+2*iie)=e*r2*(betasq+1&
                       +2*ubar*sqrt(betasq))
                  enddo
                  ! add in contribution to heating
                  naa=nint(xss(hpd+1))
                  nrr=1
                  npp=2
                  do ie=it,nes
                     e=xss(esz+ie-1)
                     call terpa(h,e,en,idis,scr(llht),npp,nrr)
                     ss=0
                     if (ie.ge.iaa) ss=xss(2+k+ie-iaa)
                     xss(hpd+2+naa+ie-it)=xss(hpd+2+naa+ie-it)+h*ss
                  enddo

               !--tabulated distributions
               else if (ltt.eq.2) then
                  call contio(nin,0,0,scr,nb,nw)
                  lct=l2h
                  call tab2io(nin,0,0,scr,nb,nw)
                  ne=n2h
                  xss(landh+jp-1)=next-andh+1
                  xss(next)=ne
                  ie=next
                  il=ie+ne
                  next=il+ne+1
                  llht=1
                  lld=llht+8+2*ne
                  scr(llht)=0
                  scr(llht+1)=0
                  scr(llht+2)=0
                  scr(llht+3)=0
                  scr(llht+4)=1
                  scr(llht+5)=ne
                  scr(llht+6)=ne
                  scr(llht+7)=2
                  amass=awr/awi
                  aprime=awp/awi
                  do iie=1,ne
                     ll=lld
                     call tab1io(nin,0,0,scr(ll),nb,nw)
                     call pttab2(scr(lld))
                     xss(ie+iie)=sigfig(scr(lld+1)/emev,7,0)
                     m=nint(scr(lld+4))
                     n=nint(scr(lld+5))
                     xss(il+iie)=next-andh+1
                     xss(il+iie)=-xss(il+iie)
                     int=nint(scr(lld+7))
                     xss(next)=int
                     xss(next+1)=n
                     if (next+2+3*n.gt.nxss) call error('acelcp',&
                   'insufficient storage for angular distributions.',&
                     ' ')
                     do i=1,n
                        xss(next+1+i)=&
                          sigfig(scr(lld+4+2*m+2*i),7,0)
                        xss(next+1+n+i)=&
                          sigfig(scr(lld+5+2*m+2*i),7,0)
                        if (xss(next+1+n+i).lt.small)&
                          xss(next+1+n+i)=0
                        if (i.eq.1) then
                           xss(next+1+2*n+i)=0
                           ubar=0
                        endif
                        if (i.gt.1.and.int.eq.1) then
                           sum=xss(next+1+2*n+i-1)&
                             +xss(next+1+n+i-1)&
                             *(xss(next+1+i)-xss(next+1+i-1))
                           xss(next+1+2*n+i)=sigfig(sum,7,0)
                           ubar=ubar&
                             +xss(next+1+n+i-1)&
                             *(xss(next+1+i)-xss(next+1+i-1))&
                             +(xss(next+1+i)+xss(next+1+i-1))/2
                        endif
                        if (i.gt.1.and.int.eq.2) then
                           sum=xss(next+1+2*n+i-1)&
                            +(xss(next+1+n+i)+xss(next+1+n+i-1))&
                            *(xss(next+1+i)-xss(next+1+i-1))/2
                           xss(next+1+2*n+i)=sigfig(sum,7,0)
                           ubar=ubar&
                             +(xss(next+1+n+i)+xss(next+1+n+i-1))&
                             *(xss(next+1+i)-xss(next+1+i-1))&
                             *(xss(next+1+i)+xss(next+1+i-1))/4
                        endif
                     enddo
                     renorm=1/xss(next+1+3*n)
                     do i=1,n
                        xss(next+1+n+i)=&
                          sigfig(renorm*xss(next+1+n+i),7,0)
                        xss(next+1+2*n+i)=&
                          sigfig(renorm*xss(next+1+2*n+i),9,0)
                     enddo
                     next=next+2+3*n
                     e=xss(ie+iie)
                     th=(1+amass)*q/amass
                     r1=amass*(amass+1-aprime)/aprime
                     r2=aprime/(1+amass)**2
                     betasq=r1*(1+th/e)
                     if (betasq.lt.zero) betasq=0
                     scr(llht+6+2*iie)=e
                     scr(llht+7+2*iie)=e*r2*(betasq+1&
                       +2*ubar*sqrt(betasq))
                  enddo
                  ! add in contribution to heating
                  naa=nint(xss(hpd+1))
                  nrr=1
                  npp=2
                  do ie=it,nes
                     e=xss(esz+ie-1)
                     call terpa(h,e,en,idis,scr(llht),npp,nrr)
                     ss=0
                     if (ie.ge.iaa) ss=xss(2+k+ie-iaa)
                     xss(hpd+2+naa+ie-it)=xss(hpd+2+naa+ie-it)+h*ss
                  enddo
               endif

            !--normal branch for other laws (mf=6)
            else
               call findf(matd,mf,mt,nin)
               call contio(nin,0,0,scr,nb,nw)
               nk=n1h
               do ik=1,nk
                  ll=1
                  lly=ll
                  call tab1io(nin,0,0,scr(ll),nb,nw)
                  izap=nint(c1h)
                  awp=c2h
                  law=l2h
                  ll=ll+nw
                  do while (nb.ne.0)
                     call moreio(nin,0,0,scr(ll),nb,nw)
                     ll=ll+nw
                  enddo
                  if (ik.eq.lprod(j)) then
                     lld=ll
                     xss(landh+jp-1)=-1
                  endif

                  !--special steps for two-body recoil distributions
                  !--back up to corresponding law=2 distribution
                  izarec=0
                  awprec=0
                  if (ik.eq.lprod(j).and.law.eq.4) then
                     izarec=izap
                     awprec=awp
                     call findf(matd,mf,mt,nin)
                     call contio(nin,0,0,scr,nb,nw)
                     ll=1
                     lly=ll
                     call tab1io(nin,0,0,scr(ll),nb,nw)
                     izap=nint(c1h)
                     awp=c2h
                     law=l2h
                     ll=ll+nw
                     do while (nb.ne.0)
                        call moreio(nin,0,0,scr(ll),nb,nw)
                        ll=ll+nw
                     enddo
                  endif

                  !--law 2 angular distribution
                  !--also used for law 4 two-body recoils
                  if ((ik.eq.lprod(j).and.law.eq.2).or.&
                    (izarec.eq.ip.and.law.eq.2)) then
                     ll=lld
                     call tab2io(nin,0,0,scr(ll),nb,nw)
                     xss(landh+jp-1)=next-andh+1
                     ne=nint(scr(ll+5))
                     xss(next)=ne
                     ie=next
                     il=ie+ne
                     next=il+ne+1
                     llht=lld
                     lld=llht+8+2*ne
                     scr(llht)=0
                     scr(llht+1)=0
                     scr(llht+2)=0
                     scr(llht+3)=0
                     scr(llht+4)=1
                     scr(llht+5)=ne
                     scr(llht+6)=ne
                     scr(llht+7)=2
                     amass=awr/awi
                     if (izarec.eq.0) then
                        aprime=awp/awi
                     else
                        aprime=awprec/awi
                     endif
                     do iie=1,ne
                        ll=lld
                        call listio(nin,0,0,scr(ll),nb,nw)
                        lang=nint(scr(lld+2))
                        if (lang.eq.0) then
                           if (izarec.ne.0) then
                              nl=nint(scr(lld+5))
                              do iil=1,nl
                                 if (mod(iil,2).eq.1) then
                                    scr(lld+5+iil)=-scr(lld+5+iil)
                                 endif
                              enddo
                           endif
                           call ptleg2(scr(lld))
                        else
                           iint=lang-10
                           nn=nint(scr(lld+5))
                           do kk=1,nw-4
                              scr(lld+2+nw-kk)=scr(lld+nw-kk)
                           enddo
                           scr(lld+4)=1
                           scr(lld+5)=nn
                           scr(lld+6)=iint
                           scr(lld+7)=iint
                           call pttab2(scr(lld))
                        endif
                        xss(ie+iie)=sigfig(scr(lld+1)/emev,7,0)
                        m=nint(scr(lld+4))
                        n=nint(scr(lld+5))
                        xss(il+iie)=next-andh+1
                        xss(il+iie)=-xss(il+iie)
                        int=nint(scr(lld+7))
                        xss(next)=int
                        xss(next+1)=n
                        if (next+2+3*n.gt.nxss) call error('acelcp',&
                  'insufficient storage for angular distributions.',&
                  ' ')
                        do i=1,n
                           xss(next+1+i)=&
                             sigfig(scr(lld+4+2*m+2*i),7,0)
                              xss(next+1+n+i)=&
                             sigfig(scr(lld+5+2*m+2*i),7,0)
                           if (xss(next+1+n+i).lt.small)&
                             xss(next+1+n+i)=0
                           if (i.eq.1) then
                              xss(next+1+2*n+i)=0
                              ubar=0
                              endif
                           if (i.gt.1.and.int.eq.1) then
                              sum=xss(next+1+2*n+i-1)&
                                +xss(next+1+n+i-1)&
                                *(xss(next+1+i)-xss(next+1+i-1))
                              xss(next+1+2*n+i)=sigfig(sum,7,0)
                              ubar=ubar&
                                +xss(next+1+n+i-1)&
                                *(xss(next+1+i)-xss(next+1+i-1))&
                                +(xss(next+1+i)+xss(next+1+i-1))/2
                           endif
                           if (i.gt.1.and.int.eq.2) then
                              sum=xss(next+1+2*n+i-1)&
                               +(xss(next+1+n+i)+xss(next+1+n+i-1))&
                               *(xss(next+1+i)-xss(next+1+i-1))/2
                              xss(next+1+2*n+i)=sigfig(sum,7,0)
                              ubar=ubar&
                                +(xss(next+1+n+i)+xss(next+1+n+i-1))&
                                *(xss(next+1+i)-xss(next+1+i-1))&
                                *(xss(next+1+i)+xss(next+1+i-1))/4
                           endif
                        enddo
                        renorm=1/xss(next+1+3*n)
                        do i=1,n
                           xss(next+1+n+i)=&
                             sigfig(renorm*xss(next+1+n+i),7,0)
                           xss(next+1+2*n+i)=&
                             sigfig(renorm*xss(next+1+2*n+i),9,0)
                        enddo
                        next=next+2+3*n
                        e=xss(ie+iie)
                        th=(1+amass)*q/amass
                        r1=amass*(amass+1-aprime)/aprime
                        r2=aprime/(1+amass)**2
                        betasq=r1*(1+th/e)
                        if (betasq.lt.zero) betasq=0
                        scr(llht+6+2*iie)=e
                        scr(llht+7+2*iie)=e*r2*(betasq+1&
                          +2*ubar*sqrt(betasq))
                     enddo
                     ! add in contribution to heating
                     naa=nint(xss(hpd+1))
                     nrr=1
                     npp=2
                     do ie=it,nes
                        e=xss(esz+ie-1)
                        call terpa(h,e,en,idis,scr(llht),npp,nrr)
                        ss=0
                        if (ie.ge.iaa) ss=xss(2+k+ie-iaa)
                        xss(hpd+2+naa+ie-it)=xss(hpd+2+naa+ie-it)+h*ss
                     enddo

                  !--law 7 angle-energy distribution
                  else if (ik.eq.lprod(j).and.law.eq.7) then
                     ll=lld
                     call tab2io(nin,0,0,scr(ll),nb,nw)
                     xss(landh+jp-1)=next-andh+1
                     ne=nint(scr(ll+5))
                     xss(next)=ne
                     na=next
                     nc=next+ne
                     next=next+1+2*ne
                     do ie=1,ne
                        ll=lld
                        call tab1io(nin,0,0,scr(ll),nb,nw)
                        ll=ll+nw
                        do while (nb.ne.0)
                           call moreio(nin,0,0,scr(ll),nb,nw)
                           ll=ll+nw
                        enddo
                        nmu=nint(scr(lld+3))
                        nx=nint(scr(lld+5))
                        intx=nint(scr(lld+7))
                        xss(na+ie)=scr(lld+1)
                        xss(nc+ie)=-(next-andh+1)
                        xss(next)=intx
                        xss(next+1)=nx
                        if (next+2+3*nx.gt.nxss) call error('acelcp',&
                  'insufficient storage for angular distributions.',&
                  ' ')
                        do ix=1,nx
                           xss(next+1+ix)=sigfig(scr(lld+6+2*ix),7,0)
                           xss(next+1+nx+ix)=&
                             sigfig(scr(lld+7+2*ix),7,0)
                           if (xss(next+1+nx+ix).lt.small)&
                             xss(next+1+nx+ix)=0
                           if (ix.eq.1) xss(next+1+2*nx+ix)=0
                           if (ix.gt.1) then
                              sum=xss(next+1+2*nx+ix-1)&
                                +(xss(next+1+nx+ix)&
                                +xss(next+1+nx+ix-1))&
                                *(xss(next+1+ix)-xss(next+1+ix-1))/2
                              xss(next+1+2*nx+ix)=sigfig(sum,7,0)
                           endif
                        enddo
                        renorm=1/xss(next+1+3*nx)
                        do ix=1,nx
                           xss(next+1+nx+ix)=&
                             sigfig(renorm*xss(next+1+nx+ix),7,0)
                           xss(next+1+2*nx+ix)=&
                             sigfig(renorm*xss(next+1+2*nx+ix),9,0)
                        enddo
                        next=next+2+3*nx
                        ll=lld
                        do imu=1,nmu
                           call tab1io(nin,0,0,scr(ll),nb,nw)
                           do while (nb.ne.0)
                              call moreio(nin,0,0,scr(ll),nb,nw)
                           enddo
                        enddo
                     enddo

                  !--skip to the next subsection
                  else
                     call skip6a(nin,0,0,scr,law)
                  endif
               enddo
            endif
         endif
      enddo

      !--now go back and get the energy distribution data
      ldlwh=next
      xss(ploct+10*(itype-1)+7)=ldlwh
      dlwh=ldlwh+ntrh
      xss(ploct+10*(itype-1)+8)=dlwh
      next=dlwh
      jp=0
      do j=1,nprod
         ipj=iprod(j)
         mt=mprod(j)
         mf=kprod(j)
         iskip=0
         if (ipj.ne.ip) iskip=1
         if (mf.eq.0) iskip=1
         if (iskip.eq.0) then
            jp=jp+1
            if (mt.eq.2) then
               xss(ldlwh+jp-1)=0
               iskip=1
            endif
         endif
         if (iskip.eq.0) then
            xss(ldlwh+jp-1)=next-dlwh+1
            last=next
            xss(next)=0
            xss(next+1)=0
            next=next+3
            mtt=0
            ir=0
            do while (mtt.ne.mt)
               ir=ir+1
               mtt=nint(xss(mtr+ir-1))
               k=nint(xss(lsig+ir-1))+sig-1
               n=nint(xss(k+1))
               iaa=nint(xss(k))
               q=xss(lqr+ir-1)
            enddo

            !--special branch for neutron mt102
            if (mt.eq.102.and.izai.eq.1) then
               xss(last+1)=33
               xss(next)=0
               xss(next+1)=2
               xss(next+2)=sigfig(xss(esz+it-1),7,0)
               xss(next+3)=1
               xss(next+4)=sigfig(xss(esz+nes-1),7,0)
               xss(next+5)=1
               next=next+2+2*2
               xss(last+2)=next-dlwh+1
               xss(next)=0
               xss(next+1)=awi/(awr+awi)
               next=next+2
               ! add in contribution to heating
               naa=nint(xss(hpd+1))
               do ie=it,nes
                  tt=xss(esz+ie-1)*xss(hpd+2+ie-it)*awi/(awr+awi)
                  xss(hpd+2+naa+ie-it)=xss(hpd+2+naa+ie-it)+tt
               enddo

            !--branch for mf=4 reactions
            else if (mf.eq.4) then
               call findf(matd,mf,mt,nin)
               call contio(nin,0,0,scr,nb,nw)
               ltt=l2h
               awp=1
               if (mt.ge.600.and.mt.le.649) awp=awr1
               if (mt.ge.650.and.mt.le.699) awp=awr2
               if (mt.ge.700.and.mt.le.749) awp=awr3
               if (mt.ge.750.and.mt.le.799) awp=awr4
               if (mt.ge.800.and.mt.le.849) awp=awr5
               xss(last+1)=33
               xss(next)=0
               xss(next+1)=2
               next=next+2+2*2
               xss(last+2)=next-dlwh+1
               amass=awr/awi
               aprime=awp/awi
               xss(next)=sigfig((1+amass)*(-q)/amass,7,0)
               xss(next+1)=&
                 sigfig(amass*(amass+1-aprime)/(1+amass)**2,7,0)
               ! add in contribution to heating
               if (ltt.eq.0) then
                  naa=nint(xss(hpd+1))
                  do ie=it,nes
                     e=xss(esz+ie-1)
                     ss=0
                     if (ie.ge.iaa) ss=xss(2+k+ie-iaa)
                     tt=xss(next+1)*(e-xss(next))*ss
                     xss(hpd+2+naa+ie-it)=xss(hpd+2+naa+ie-it)+tt
                  enddo
               endif
               next=next+2

            !--normal branch for other reactions (mf=6)
            else

               !--first check the subsection to see whether
               !--the distribution is isotropic or not.
               isocp=1
               call findf(matd,mf,mt,nin)
               call contio(nin,0,0,scr,nb,nw)
               nk=n1h
               ik=0
               idone=0
               do while (ik.lt.nk.and.idone.eq.0)
                  ik=ik+1
                  ll=1
                  lly=ll
                  call tab1io(nin,0,0,scr(ll),nb,nw)
                  izap=nint(c1h)
                  awp=c2h
                  law=l2h
                  ll=ll+nw
                  do while (nb.ne.0)
                     call moreio(nin,0,0,scr(ll),nb,nw)
                     ll=ll+nw
                  enddo

                  !--if not the desired particle, skip the subsection
                  if (ik.ne.lprod(j).or.law.ne.1) then
                     call skip6a(nin,0,0,scr,law)

                  !--we only need to check law 1 subsections
                  else
                     call tab2io(nin,0,0,scr(ll),nb,nw)
                     lang=nint(scr(ll+2))
                     lep=nint(scr(ll+3))
                     ne=nint(scr(ll+5))
                     do ie=1,ne
                        ll=lld
                        call listio(nin,0,0,scr(ll),nb,nw)
                        ll=ll+nw
                        do while (nb.ne.0)
                           call moreio(nin,0,0,scr(ll),nb,nw)
                           ll=ll+nw
                        enddo
                        na=nint(scr(lld+3))
                        if (na.gt.0) isocp=0
                     enddo
                  endif
               enddo

               !--go back and process the subsection
               call findf(matd,mf,mt,nin)
               call contio(nin,0,0,scr,nb,nw)
               nk=n1h
               ik=0
               idone=0
               do while (ik.lt.nk.and.idone.eq.0)
                  ik=ik+1
                  ll=1
                  lly=ll
                  call tab1io(nin,0,0,scr(ll),nb,nw)
                  izap=nint(c1h)
                  awp=c2h
                  law=l2h
                  ll=ll+nw
                  do while (nb.ne.0)
                     call moreio(nin,0,0,scr(ll),nb,nw)
                     ll=ll+nw
                  enddo

                  !--if not the desired particle, skip the subsection
                  if (izap.ne.ip) then
                     call skip6a(nin,0,0,scr,law)

                  !--mt102 subsection
                  else if (mt.eq.102) then
                     xss(last+1)=33
                     xss(next)=0
                     xss(next+1)=2
                     xss(next+2)=sigfig(xss(esz+it-1),7,0)
                     xss(next+3)=1
                     xss(next+4)=sigfig(xss(esz+nes-1),7,0)
                     xss(next+5)=1
                     next=next+2+2*2
                     xss(last+2)=next-dlwh+1
                     xss(next)=0
                     xss(next+1)=awi/(awr+awi)
                     next=next+2
                     ! add in contribution to heating
                     naa=nint(xss(hpd+1))
                     do ie=it,nes
                        tt=xss(esz+ie-1)*xss(hpd+2+ie-it)*awi/(awr+awi)
                        xss(hpd+2+naa+ie-it)=xss(hpd+2+naa+ie-it)+tt
                     enddo
                     call skip6a(nin,0,0,scr,law)

                  !--skip if not correct subsection
                  else if (ik.ne.lprod(j)) then
                    call skip6a(nin,0,0,scr,law)

                  !--law 1 or 2 subsections
                  else if (law.eq.1.or.law.eq.2) then
                     call tab2io(nin,0,0,scr(ll),nb,nw)
                     lang=nint(scr(ll+2))
                     lawnow=0
                     if (law.eq.1.and.lang.eq.1.and.isocp.eq.0)&
                       lawnow=61
                     if (law.eq.1.and.lang.eq.1.and.isocp.eq.1)&
                       lawnow=4
                     if (law.eq.1.and.lang.eq.2) lawnow=44
                     if (law.eq.2) lawnow=33
                     if (law.eq.1.and.lang.ge.11) lawnow=61
                     if (lawnow.eq.0) call error('acelcp',&
                       'unsupported law and lang',' ')
                     xss(last+1)=lawnow
                     if (law.eq.1.and.isocp.eq.1) xss(landh+jp-1)=0
                     lep=nint(scr(ll+3))
                     ne=nint(scr(ll+5))
                     xss(next)=0
                     lee=next
                     xss(next+1)=2
                     next=next+2+2*2
                     xss(last+2)=next-dlwh+1
                     if (law.eq.1) then
                        xss(next)=0
                        xss(next+1)=ne
                        lle=next+2
                        next=lle+2*ne
                     endif
                     llh=ll
                     scr(llh)=0
                     scr(llh+1)=0
                     scr(llh+2)=0
                     scr(llh+3)=0
                     scr(llh+4)=1
                     scr(llh+5)=ne
                     scr(llh+6)=ne
                     scr(llh+7)=2
                     lld=llh+8+2*ne
                     do ie=1,ne
                        ll=lld
                        call listio(nin,0,0,scr(ll),nb,nw)
                        ll=ll+nw
                        do while (nb.ne.0)
                           call moreio(nin,0,0,scr(ll),nb,nw)
                           ll=ll+nw
                        enddo
                        if (ie.eq.1) then
                           xss(lee+2)=sigfig(scr(lld+1)/emev,7,0)
                           xss(lee+4)=1
                        else if (ie.eq.ne) then
                           xss(lee+3)=sigfig(scr(lld+1)/emev,7,0)
                           xss(lee+5)=1
                        endif
                        if (law.ne.2) then
                           ! law=1
                           xss(lle+ie-1)=sigfig(scr(lld+1)/emev,7,0)
                           ee=xss(lle+ie-1)
                           xss(lle+ne+ie-1)=next-dlwh+1
                           nd=nint(scr(lld+2))
                           na=nint(scr(lld+3))
                           ncyc=na+2
                           ng=nint(scr(lld+5))
                           if (lawnow.eq.4) then
                              xss(next)=lep
                           else
                              xss(next)=lep+10*nd
                           endif
                           xss(next+1)=ng
                           nexcd=next+4*ng+2
                           amass=awp*emc2
                           avadd=awi*sqrt(2*ee/(emc2*awi))/(awi+awr)
                           avlab=0
                           avll=0
                           do ig=1,ng
                              ! distribution
                              xss(next+1+ig)=&
                                sigfig(scr(lld+6+ncyc*(ig-1))/emev,7,0)
                              xss(next+1+ig+ng)=&
                                sigfig(scr(lld+7+ncyc*(ig-1))*emev,7,0)
                              test=xss(next+1+ig+ng)
                              if (test.gt.0.and.test.lt.small)&
                                xss(next+1+ig+ng)=small
                              xss(next+1+ig+2*ng)=0
                              if (ig.le.nd) xss(next+1+ig+2*ng)=&
                                xss(next+1+ig+2*ng)&
                                +scr(lld+7+ncyc*(ig-1))
                              if (nd.gt.0.and.ig.eq.nd+1)&
                                xss(next+1+ig+2*ng)=xss(next+ig+2*ng)
                              if (ig.gt.nd+1.and.lep.eq.1)&
                                xss(next+1+ig+2*ng)=xss(next+ig+2*ng)&
                                +scr(lld+7+ncyc*(ig-2))&
                                *(scr(lld+6+ncyc*(ig-1))&
                                -scr(lld+6+ncyc*(ig-2)))
                              if (ig.gt.nd+1.and.lep.eq.2)&
                                xss(next+1+ig+2*ng)=xss(next+ig+2*ng)&
                                +((scr(lld+7+ncyc*(ig-2))&
                                +scr(lld+7+ncyc*(ig-1)))/2)&
                                *(scr(lld+6+ncyc*(ig-1))&
                                -scr(lld+6+ncyc*(ig-2)))
                              rkal=0
                              akal=0
                              ! kalbach distribution
                              if (lang.eq.2) then
                                 rkal=scr(lld+8+ncyc*(ig-1))
                                 xss(next+1+ig+3*ng)=sigfig(rkal,7,0)
                                 ep=xss(next+1+ig)
                                 if (na.eq.2) then
                                    akal=scr(lld+9+ncyc*(ig-1))
                                 else
                                    akal=bachaa(izai,izap,iza,ee,ep)
                                 endif
                                 xss(next+1+ig+4*ng)=sigfig(akal,7,0)
                              ! legendre or tabulated distribution
                              else if (lawnow.eq.61) then
                                 ep=xss(next+1+ig)
                                 if (lang.eq.1) then
                                    scr(ll)=0
                                    scr(ll+1)=ep
                                    scr(ll+2)=0
                                    scr(ll+3)=0
                                    scr(ll+4)=na
                                    scr(ll+5)=0
                                    do ia=1,na
                                       lll=lld+7+ncyc*(ig-1)
                                       scr(ll+5+ia)=0
                                       if (scr(lll).ne.zero) then
                                          scr(ll+5+ia)=scr(lll+ia)/scr(lll)
                                       endif
                                    enddo
                                    call ptleg2(scr(ll))
                                    intmu=2
                                    nmu=nint(scr(ll+5))
                                    llx=ll+6
                                 else
                                    intmu=lang-10
                                    nmu=na/2
                                    llx=lld+6
                                 endif
                                 xss(next+1+3*ng+ig)=nexcd-dlwh+1
                                 xss(nexcd)=intmu
                                 xss(nexcd+1)=nmu
                                 do imu=1,nmu
                                    xss(nexcd+1+imu)=scr(llx+2*imu)
                                    xss(nexcd+1+nmu+imu)=scr(llx+1+2*imu)
                                    if (imu.eq.1) then
                                       sum=0
                                       xss(nexcd+1+2*nmu+imu)=0
                                    else
                                       del=&
                                         scr(llx+2*imu)-scr(llx+2*imu-2)
                                       if (intmu.eq.1) then
                                          sum=sum+del*scr(llx+1+2*imu-2)
                                          xss(nexcd+1+2*nmu+imu)=sum
                                       else
                                          av=(scr(llx+1+2*imu)&
                                            +scr(llx+1+2*imu-2))/2
                                          sum=sum+del*av
                                          xss(nexcd+1+2*nmu+imu)=sum
                                       endif
                                    endif
                                 enddo
                                 do imu=1,nmu
                                    xss(nexcd+1+imu)=sigfig(&
                                      xss(nexcd+1+imu)/sum,7,0)
                                    xss(nexcd+1+nmu+imu)=sigfig(&
                                      xss(nexcd+1+nmu+imu)/sum,7,0)
                                    xss(nexcd+1+2*nmu+imu)=sigfig(&
                                      xss(nexcd+1+2*nmu+imu)/sum,7,0)
                                 enddo
                                 nexcd=nexcd+2+3*nmu
                              endif
                              ! average lab energy
                              if (ig.ne.1) then
                                 eavi=xss(next+1+ig)
                                 if (na.eq.0) then
                                    avl=eavi
                                 else
                                    avcm=sqrt(2*eavi/amass)
                                    sign=1
                                    avl=eavl(akal,amass,avcm,&
                                      avadd,rkal,sign)
                                 endif
                                 dele=xss(next+1+ig)-xss(next+ig)
                                 if (lep.eq.1) then
                                    avav=xss(next+ig+ng)*(avll+avl)/2
                                 else
                                    avav=(xss(next+ig+ng)*avll&
                                      +xss(next+ig+1+ng)*avl)/2
                                 endif
                                 avlab=avlab+avav*dele
                                 avll=avl
                              endif
                           enddo
                           renorm=1/xss(next+1+3*ng)
                           do ig=1,ng
                              xss(next+1+ng+ig)=&
                                sigfig(renorm*xss(next+1+ng+ig),7,0)
                              xss(next+1+2*ng+ig)=&
                                sigfig(renorm*xss(next+1+2*ng+ig),9,0)
                           enddo
                           scr(llh+6+2*ie)=ee
                           scr(llh+7+2*ie)=avlab
                           if (lawnow.eq.61) then
                              next=nexcd
                           else
                              next=next+2+(2*na+3)*ng
                           endif
                        endif
                     enddo
                     ! law2
                     if (law.eq.2) then
                        amass=awr/awi
                        aprime=awp/awi
                        xss(next)=sigfig((1+amass)*(-q)/amass,7,0)
                        xss(next+1)=&
                          sigfig(amass*(amass+1-aprime)/(1+amass)**2,&
                          7,0)
                        next=next+2
                     else
                     ! add in contribution to heating
                     ! for this subsection
                        if (lld+nes-it.gt.nwscr)&
                          call error('acelcp',&
                          'scratch array overflowing.',&
                          'reduce the number of energy points. ')
                        nrr=1
                        npp=2
                        do ie=it,nes
                           e=xss(esz+ie-1)*emev
                           call terpa(y,e,en,idis,scr(lly),npp,nrr)
                           ss=0
                           if (ie.ge.iaa) ss=xss(2+k+ie-iaa)
                           scr(lld+ie-it)=y*ss
                        enddo
                        naa=nint(xss(hpd+1))
                        nrr=1
                        npp=2
                        do ie=it,nes
                           e=xss(esz+ie-1)
                           if (law.eq.1) then
                              call terpa(h,e,en,idis,scr(llh),npp,nrr)
                           else
                              tt1=amass*(amass+1-aprime)/aprime
                              tt2=1+((1+amass)/amass)*(q/e)
                              h=e*(1+tt1*tt2)*aprime/(1+amass)**2
                           endif
                           xss(hpd+2+naa+ie-it)=xss(hpd+2+naa+ie-it)&
                             +h*scr(lld+ie-it)
                        enddo
                     endif
                     idone=1

                  !--law 3 isotropic two-body distributions
                  !--and law 4 two-body recoil distributions
                  else if (law.eq.3.or.law.eq.4) then
                     xss(last+1)=33
                     if (law.eq.3) xss(landh+jp-1)=0
                     xss(next)=0
                     xss(next+1)=2
                     next=next+2+2*2
                     xss(last+2)=next-dlwh+1
                     amass=awr/awi
                     aprime=awp/awi
                     xss(next)=sigfig((1+amass)*(-q)/amass,7,0)
                     xss(next+1)=&
                       sigfig(amass*(amass+1-aprime)/(1+amass)**2,7,0)
                     if (law.ne.4) then
                        ! add in contribution to heating
                        naa=nint(xss(hpd+1))
                        do ie=it,nes
                           e=xss(esz+ie-1)
                           ss=0
                           if (ie.ge.iaa) ss=xss(2+k+ie-iaa)
                           tt=xss(next+1)*(e-xss(next))*ss
                           xss(hpd+2+naa+ie-it)=xss(hpd+2+naa+ie-it)+tt
                        enddo
                     endif
                     next=next+2
                     idone=1

                  !--law 6, phase space
                  else if (law.eq.6) then
                     xss(last+1)=66
                     call contio(nin,0,0,scr(ll),nb,nw)
                     apsx=scr(ll)
                     npsx=nint(scr(ll+5))
                     xss(next)=0
                     lee=next
                     xss(next+1)=2
                     next=next+2+2*2
                     xss(last+2)=next-dlwh+1
                     xss(next)=npsx
                     xss(next+1)=apsx
                     xss(next+2)=2
                     next=next+3
                     step1=ten**(one/5)
                     step2=one/50
                     xx=elow
                     n=1
                     test1=one+one/100000
                     test2=one/10-one/1000000
                     test3=one-one/100000
                     do while (xx.lt.test1)
                        n=n+1
                        if (xx.lt.test2) then
                            xx=xx*step1
                        else
                           xx=xx+step2
                        endif
                     enddo
                     nn=n
                     xss(next)=nn
                     xl=0
                     pl=0
                     yn=0
                     n=1
                     xss(next+n)=xl
                     xss(next+nn+n)=pl
                     xss(next+2*nn+n)=yn
                     xx=elow
                     do while (xx.lt.test1)
                        n=n+1
                        if (xx.gt.test3) then
                           xx=1
                           pn=0
                        else
                           rn=3
                           pn=sqrt(xx)*(1-xx)**(rn*npsx/2-4)
                        endif
                        yn=yn+(xx-xl)*(pn+pl)/2
                        xss(next+n)=sigfig(xx,7,0)
                        xss(next+nn+n)=pn
                        xss(next+2*nn+n)=yn
                        xl=xx
                        pl=pn
                        if (xx.lt.test2) then
                           xx=xx*step1
                        else
                           xx=xx+step2
                        endif
                     enddo
                     sum=yn
                     do kk=1,nn
                        xss(next+nn+kk)=&
                          sigfig(xss(next+nn+kk)/sum,7,0)
                        xss(next+2*nn+kk)=&
                          sigfig(xss(next+2*nn+kk)/sum,9,0)
                     enddo
                     naa=nint(xss(hpd+1))
                     nrr=1
                     npp=2
                     ! add in this part of the heating
                     do ie=it,nes
                        e=xss(esz+ie-1)
                        ee=e*emev
                        ecm=e*awi/(awi+awr)
                        ecm=ecm*awp/(awi+awr)
                        ea=e*awr/(awi+awr)+q
                        eim=ea*(apsx-awp)/apsx
                        h=0
                        chek=0
                        dmu=1
                        dmu=dmu/5
                        do imu=1,11
                           amuu=-1+dmu*(imu-1)
                           e1l=0
                           e2l=0
                           sum=0
                           chk=0
                           do iep=1,nn
                              ep=xss(next+iep)*eim
                              pp=xss(next+nn+iep)/eim
                              disc=ep-ecm*(1-amuu**2)
                              if (disc.ge.zero) then
                                 v1=amuu*sqrt(ecm)+sqrt(disc)
                                 v2=amuu*sqrt(ecm)-sqrt(disc)
                                 e1=v1**2
                                 e2=v2**2
                                 if (v2.lt.zero) e2=-e2
                                 pp1=0
                                 if (ep.gt.zero) pp1=pp*sqrt(e1)&
                                   /sqrt(ep)
                                 pp2=0
                                 if (v2.gt.zero.and.ep.gt.zero) then
                                    pp2=pp*sqrt(e2)/sqrt(ep)
                                 endif
                                 if (e1l.gt.zero) then
                                    de=e1-e1l
                                    if (de.lt.zero) de=-de
                                    sum=sum+de*(e1*pp1+e1l*pp1l)/2
                                    chk=chk+de*(pp1+pp1l)/2
                                    if (e2.gt.zero.and.&
                                      e2l.gt.zero) then
                                       de=e2l-e2
                                       if (de.lt.zero) de=-de
                                       sum=sum+de*(e2*pp2+e2l*pp2l)/2
                                       chk=chk+de*(pp2+pp2l)/2
                                    endif
                                 endif
                                 e1l=e1
                                 e2l=e2
                                 pp1l=pp1
                                 pp2l=pp2
                              endif
                           enddo
                           eav=sum
                           if (chk.ne.zero) eav=eav/chk
                           if (imu.gt.1) h=h+(sum+suml)/10
                           if (imu.gt.1) chek=chek+(chk+chkl)/10
                           suml=sum
                           chkl=chk
                        enddo
                        if (chek.gt.zero) h=h/chek
                        call terpa(y,ee,en,idis,scr(lly),npp,nrr)
                        ss=0
                        if (ie.ge.iaa) ss=xss(2+k+ie-iaa)
                        xss(hpd+2+naa+ie-it)=xss(hpd+2+naa+ie-it)&
                          +h*y*ss
                     enddo
                     next=next+1+3*nn
                     idone=1

                  !--law 7, angle-energy
                  else if (law.eq.7) then
                     xss(last+1)=67
                     call tab2io(nin,0,0,scr(ll),nb,nw)
                     lep=nint(scr(ll+3))
                     ne=nint(scr(ll+5))
                     xss(last+2)=next-dlwh+1
                     xss(next)=0
                     xss(next+1)=ne
                     lee=next+2
                     next=lee+2*ne
                     llht=ll
                     llad=llht+8+2*ne
                     scr(llht)=0
                     scr(llht+1)=0
                     scr(llht+2)=0
                     scr(llht+3)=0
                     scr(llht+4)=1
                     scr(llht+5)=ne
                     scr(llht+6)=ne
                     scr(llht+7)=2
                     nrr=1
                     npp=2
                     do ie=1,ne
                        ll=llad
                           call tab1io(nin,0,0,scr(ll),nb,nw)
                        intmu=l1h
                        nmu=l2h
                        e=c2h
                        ee=c2h/emev
                        xss(lee+ie-1)=sigfig(ee,7,0)
                        xss(lee+ne+ie-1)=next-dlwh+1
                        xss(next)=intmu
                        xss(next+1)=nmu
                        next=next+2
                        mus=next
                        next=next+2*nmu
                        lld=llad+nw
                        nra=1
                        npa=2
                        ebar=0
                        chek=0
                        do imu=1,nmu
                           ll=lld
                           call tab1io(nin,0,0,scr(ll),nb,nw)
                           npep=n2h
                           intep=nint(scr(ll+7))
                           ll=ll+nw
                           do while (nb.ne.0)
                              call moreio(nin,0,0,scr(ll),nb,nw)
                              ll=ll+nw
                           enddo
                           xss(mus+imu-1)=c2h
                           amuu=c2h
                           xss(mus+nmu+imu-1)=next-dlwh+1
                           xss(next)=intep
                           xss(next+1)=npep
                           next=next+1
                           xss(next+1+2*npep)=0
                           eav=0
                           chk=0
                           do ki=1,npep
                              xss(next+ki)=&
                                sigfig(scr(lld+8+2*(ki-1))/emev,7,0)
                              xss(next+npep+ki)=&
                                sigfig(scr(lld+9+2*(ki-1))*emev,7,0)
                              if (ki.ne.1) then
                                 if (intep.eq.1) then
                                    xss(next+2*npep+ki)=&
                                      xss(next+2*npep+ki-1)&
                                      +scr(lld+9+2*(ki-2))&
                                      *(scr(lld+8+2*(ki-1))&
                                      -scr(lld+8+2*(ki-2)))
                                    eav=eav&
                                      +scr(lld+9+2*(ki-2))&
                                      *(scr(lld+8+2*(ki-1))&
                                      -scr(lld+8+2*(ki-2)))&
                                      *(scr(lld+8+2*(ki-1))&
                                      +scr(lld+8+2*(ki-2)))/2
                                    chk=chk&
                                      +scr(lld+9+2*(ki-2))&
                                      *(scr(lld+8+2*(ki-1))&
                                      -scr(lld+8+2*(ki-2)))
                                 else if (intep.eq.2) then
                                    xss(next+2*npep+ki)=&
                                      xss(next+2*npep+ki-1)&
                                      +(scr(lld+9+2*(ki-2))&
                                      +scr(lld+9+2*(ki-1)))&
                                      *(scr(lld+8+2*(ki-1))&
                                      -scr(lld+8+2*(ki-2)))/2
                                    eav=eav&
                                      +(scr(lld+9+2*(ki-2))&
                                      *scr(lld+8+2*(ki-2))&
                                      +scr(lld+9+2*(ki-1))&
                                      *scr(lld+8+2*(ki-1)))&
                                      *(scr(lld+8+2*(ki-1))&
                                      -scr(lld+8+2*(ki-2)))/2
                                    chk=chk&
                                      +(scr(lld+9+2*(ki-2))&
                                      +scr(lld+9+2*(ki-1)))&
                                      *(scr(lld+8+2*(ki-1))&
                                      -scr(lld+8+2*(ki-2)))/2
                                 endif
                              endif
                           enddo
                           renorm=1
                           if (xss(next+3*npep).ne.zero)&
                             renorm=1/xss(next+3*npep)
                           do ki=1,npep
                              xss(next+npep+ki)=&
                                sigfig(renorm*xss(next+npep+ki),7,0)
                              xss(next+2*npep+ki)=&
                                sigfig(renorm*xss(next+2*npep+ki),9,0)
                           enddo
                           next=next+3*npep+1
                           call terpa(ad,amuu,amun,idis,&
                             scr(llad),npa,nra)
                           eav=ad*eav
                           chk=ad*chk
                           if (imu.gt.1) then
                              ebar=ebar+(amuu-amulst)*(eav+eavlst)/2
                              chek=chek+(amuu-amulst)*(chk+chklst)/2
                           endif
                           eavlst=eav
                           amulst=amuu
                           chklst=chk
                        enddo
                        call terpa(y,e,en,idis,scr(lly),npp,nrr)
                        scr(llht+6+2*ie)=ee
                        scr(llht+7+2*ie)=y*ebar/emev
                     enddo
                     ! add in this contribution to the heating
                     naa=nint(xss(hpd+1))
                     nrr=1
                     npp=2
                     do ie=it,nes
                        e=xss(esz+ie-1)
                        call terpa(h,e,en,idis,scr(llht),npp,nrr)
                        ss=0
                        if (ie.ge.iaa) ss=xss(2+k+ie-iaa)
                        xss(hpd+2+naa+ie-it)=xss(hpd+2+naa+ie-it)+h*ss
                     enddo
                     idone=1
                  endif
               enddo
            endif
         endif

      !--continue the loop over contributions to this production
      enddo

      !--divide the heating contribution by the total xsec
      naa=nint(xss(hpd+1))
      do ie=it,nes
         if (xss(esz+nes+ie-1).ne.zero) then
            xss(hpd+2+naa+ie-it)=xss(hpd+2+naa+ie-it)&
              /xss(esz+nes+ie-1)
         endif
         if (xss(hpd+2+naa+ie-it).lt.delt) xss(hpd+2+naa+ie-it)=0
         if (ip.eq.1) xss(hpd+2+naa+ie-it)=0
         if (izai.gt.1) then
            xss(esz+4*nes+ie-1)=xss(esz+4*nes+ie-1)&
              +xss(hpd+2+naa+ie-it)
         endif
         xss(hpd+2+naa+ie-it)=sigfig(xss(hpd+2+naa+ie-it),7,0)
         xss(esz+4*nes+ie-1)=sigfig(xss(esz+4*nes+ie-1),7,0)
      enddo

      !--fill in the yh block
      yh=next
      xss(ploct+10*(itype-1)+9)=yh
      xss(yh)=0
      next=next+1
      do j=1,nprod
         if (iprod(j).eq.ipp) then
            xss(next)=mprod(j)
            next=next+1
         endif
      enddo
      xss(yh)=next-yh-1

   !--continue loop over production types
   enddo
   len2=next-1
   if (izai.le.1) then
     deallocate(scr)
     return
   endif

   !--for incident charged particles,
   !--go back through file 6 and get heating from recoils
   call findf(matd,6,0,nin)
   do while (mfh.eq.6)
      call contio(nin,0,0,scr,nb,nw)
      if (mfh.eq.6.and.mth.gt.2) then
         mt=mth
         mtt=0
         ir=0
         do while (mtt.ne.mt)
            ir=ir+1
            mtt=nint(xss(mtr+ir-1))
            k=nint(xss(lsig+ir-1))+sig-1
            n=nint(xss(k+1))
            iaa=nint(xss(k))
            q=xss(lqr+ir-1)
         enddo
         nk=n1h
         lly=1
         do ik=1,nk
            ll=lly
            call tab1io(nin,0,0,scr(ll),nb,nw)
            izap=nint(c1h)
            awp=c2h
            law=l2h
            ll=ll+nw
            do while (nb.ne.0)
               call moreio(nin,0,0,scr(ll),nb,nw)
               ll=ll+nw
            enddo

            !--compute the heating from this recoil nuclide
            if (izap.gt.2004) then

               !--law 1
               if (law.eq.1) then
                  call tab2io(nin,0,0,scr(ll),nb,nw)
                  lep=nint(scr(ll+3))
                  ne=nint(scr(ll+5))
                  llh=ll
                  scr(llh)=0
                  scr(llh+1)=0
                  scr(llh+2)=0
                  scr(llh+3)=0
                  scr(llh+4)=1
                  scr(llh+5)=ne
                  scr(llh+6)=ne
                  scr(llh+7)=2
                  lld=llh+8+2*ne
                  nrr=1
                  npp=2
                  do ie=1,ne
                     ll=lld
                     call listio(nin,0,0,scr(ll),nb,nw)
                     ll=ll+nw
                     do while (nb.ne.0)
                        call moreio(nin,0,0,scr(ll),nb,nw)
                        ll=ll+nw
                     enddo
                     e=c2h
                     if (law.ne.2) then
                        heat=0
                        np=nint(scr(lld+5))
                        call terpa(y,e,en,idis,scr(lly),npp,nrr)
                        do ip=1,np
                           ep=scr(lld+4+2*ip)
                           g=scr(lld+5+2*ip)
                           if (ip.gt.1) then
                              heat=heat+(ep-epl)*gl*(ep+epl)/2
                           endif
                              epl=ep
                           gl=g
                        enddo
                        scr(llh+6+2*ie)=e
                        scr(llh+7+2*ie)=y*heat
                     endif
                  enddo
                  mtt=0
                  ir=0
                  do while (mtt.ne.mth)
                     ir=ir+1
                     mtt=nint(xss(mtr+ir-1))
                     k=nint(xss(lsig+ir-1))+sig-1
                     n=nint(xss(k+1))
                     iaa=nint(xss(k))
                  enddo
                  nrr=1
                  npp=2
                  do ie=iaa,nes
                  e=xss(esz+ie-1)*emev
                     if (law.eq.1) then
                        call terpa(h,e,en,idis,scr(llh),npp,nrr)
                     else
                        amass=awr/awp
                        h=2*amass*e/(1+amass)**2
                     endif
                     h=(h/emev)*xss(2+k+ie-iaa)/xss(esz+nes+ie-1)
                     xss(esz+4*nes+ie-1)=&
                       sigfig(xss(esz+4*nes+ie-1)+h,7,0)
                  enddo
               else if (law.eq.4) then
                  izarec=izap
                  awprec=awp
                  call findf(matd,6,mt,nin)
                  call contio(nin,0,0,scr,nb,nw)
                  nk=n1h
                  lly=1
                  llht=lly
                  llad=llht+8+2*ne
                  scr(llht)=0
                  scr(llht+1)=0
                  scr(llht+2)=0
                  scr(llht+3)=0
                  scr(llht+4)=1
                  scr(llht+5)=ne
                  scr(llht+6)=ne
                  scr(llht+7)=2
                  lly=lly+8+2*ne
                  ll=lly
                  call tab1io(nin,0,0,scr(ll),nb,nw)
                  izap=nint(c1h)
                  awp=c2h
                  law=l2h
                  aprime=awp/awr
                  ll=ll+nw
                  do while (nb.ne.0)
                     call moreio(nin,0,0,scr(ll),nb,nw)
                     ll=ll+nw
                  enddo
                  call tab2io(nin,0,0,scr(ll),nb,nw)
                  ne=nint(scr(ll+5))
                  th=(1+amass)*q/amass
                  r1=amass*(amass+1-aprime)/aprime
                  r2=(amass+1-aprime)/(1+amass)**2
                  betasq=r1*(1+th/e)
                  if (betasq.lt.zero) betasq=0
                  gammsq=(aprime/(amass+1-aprime))**2*betasq
                  do ie=1,ne
                     call listio(nin,0,0,scr(ll),nb,nw)
                     e=scr(ll+1)
                     ubar=scr(ll+6)
                     scr(llht+6+2*ie)=e
                     scr(llht+7+2*ie)=e*r2*(gammsq+1&
                       -2*ubar*sqrt(gammsq))
                  enddo
                  nrr=1
                  npp=2
                  do ie=iaa,nes
                     e=xss(esz+ie-1)*emev
                     call terpa(h,e,en,idis,scr(llht),npp,nrr)
                     h=(h/emev)*xss(2+k+ie-iaa)/xss(esz+nes+ie-1)
                     xss(esz+4*nes+ie-1)=&
                       sigfig(xss(esz+4*nes+ie-1)+h,7,0)
                  enddo
               else if (law.eq.6) then
                  write(nsyso,'('' warning: law=6 heating for '',&
                    &''mt='',i3,'' recoil neglected'')') mth
               else if (law.eq.0) then
                  write(strng,'(''no heating info for recoil '',&
                               &''particle '',i5)')izap
                  call mess('acelcp',strng,' ')
               else
                  call mess('acelcp','unexpected law for recoil',&
                    'subsection skipped')
               endif

            !--skip to the next subsection
            else
               call skip6a(nin,0,0,scr,law)
            endif
         enddo
      endif
      call tosend(nin,0,0,scr)

   !--continue loop over mts
   enddo

   !--finished
   deallocate(scr)
   return
   end subroutine acelcp

   subroutine aceprt(hk)
   !-------------------------------------------------------------------
   ! Print the ACE file from the information in memory.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use acecm ! provides mtname
   ! externals
   character(70)::hk
   ! internals
   integer::lll,k,iaa,ib,ic,i1,i2,i3,i4,i5,id,ig,list,l
   integer::na,nb,nc,j,n,nc2,nkk,m,ii,nr1,nn,itl,ne,nbin1
   integer::int,np,jnu,kf,idone,nure,nurb,lurt,luri,lura,lurf
   integer::ie,ll,i,nlaw,loct,law,loci,intt
   real(kr)::e,eg
   integer::imn(8),imx(8),loc(8)
   character(10)::title(16)
   character(10)::name
   character(15)::kk(40)
   character(6),parameter::blank='      '
   character(12)::dashes='------------'
   character(6),parameter::ek='energy'
   character(10),dimension(3),parameter::hlabl1=(/&
     ' elastic  ','   center ','of mass   '/)
   character(10),dimension(2),parameter::hlabl2=(/&
     '     labor','atory     '/)

   !--print information block
   write(nsyso,'(////////&
     &38x,''zaid'',1x,a13/39x,''awr'',f10.3/&
     &38x,''temp'',1p,e10.2/38x,''date'',a10/39x,''mat'',a10/&
     &6x,''***********************''/&
     &6x,''*                     *'',9x,''len2'',i10/&
     &6x,''*        fast         *'',10x,''nes'',i10/&
     &6x,''*                     *'',10x,''ntr'',i10/&
     &6x,''*   ace format file   *'',11x,''nr'',i10/&
     &6x,''*                     *'',9x,''ntrp'',i10/&
     &6x,''*     processed by    *'',8x,''ntype'',i10/&
     &6x,''*                     *'',9x,''ndnf'',i10/&
     &6x,''*        njoy         *'',10x,''esz'',i10/&
     &6x,''*                     *'',11x,''nu'',i10/&
     &6x,''***********************'',10x,''mtr'',i10/&
     &39x,''lqr'',i10/39x,''tyr'',i10/38x,''lsig'',i10/&
     &39x,''sig'',i10/38x,''land'',i10/39x,''and'',i10/&
     &38x,''ldlw'',i10//39x,''dlw'',i10/39x,''gpd'',i10/&
     &38x,''mtrp'',i10/37x,''lsigp'',i10/38x,''sigp'',i10/&
     &37x,''landp'',i10/38x,''andp'',i10/37x,''ldlwp'',i10/&
     &38x,''dlwp'',i10/40x,''yp'',i10/39x,''fis'',i10/&
     &39x,''end'',i10/37x,''iurpt'',i10/39x,''nud'',i10/&
     &37x,''dndat'',i10/38x,''ldnd'',i10/39x,''dnd'',i10/&
     &37x,''ptype'',i10/38x,''ntro'',i10/37x,''ploct'',i10///&
     &6x,''hk---'',a70)')&
     hz,aw0,tz,hd,hm,len2,nes,ntr,nr,ntrp,ntype,ndnf,esz,&
     nu,mtr,lqr,tyr,lsig,sig,land,and,ldlw,dlw,gpd,mtrp,lsigp,&
     sigp,landp,andp,ldlwp,dlwp,yp,fis,end,iurpt,nud,dndat,ldnd,&
     dnd,ptype,ntro,ploct,hk

   !--temporarily change pointers
   esz=esz-1
   mtr=mtr-1
   tyr=tyr-1
   lsig=lsig-1
   sig=sig-1
   and=and-1
   ldlw=ldlw-1
   dlw=dlw-1
   lqr=lqr-1

   !--print reaction descriptions
   write(nsyso,'(//&
     &7x,''reaction descriptors''/7x,''--------------------'')')
   write(nsyso,'(/&
     &7x,''reaction'',8x,''mt'',5x,''tyr'',6x,''lsig'',6x,&
     &''land'',6x,''ldlw'',3x,&
     &''           emin           emax              q''/&
     &7x,''--------'',8x,''--'',5x,''---'',6x,''----'',6x,&
     &''----'',6x,''----'',3x,&
     &''   ------------   ------------   ------------'')')
   lll=nint(xss(land))
   write(nsyso,'(7x,''elastic'',10x,''2'',16x,i12,13x,&
     &1p,2e15.6)') lll,xss(esz+1),xss(esz+nes)
   if (ntr.ne.0) then
      do i=1,ntr
         k=nint(xss(lsig+i)+sig)
         iaa=nint(xss(k))
         ib=iaa+nint(xss(k+1))-1
         ic=nint(xss(mtr+i))
         call mtname(ic,name,izai)
         if (name(1:1).eq.'(') then
            if (izai.eq.0) name(2:2)='g'
            if (izai.eq.1001) name(2:2)='p'
            if (izai.eq.1002) name(2:2)='d'
            if (izai.eq.1003) name(2:2)='t'
            if (izai.eq.2003) name(2:2)='s'
            if (izai.eq.2004) name(2:2)='a'
         endif
         i1=nint(xss(mtr+i))
         i2=nint(xss(tyr+i))
         i3=nint(xss(lsig+i))
         if (i.le.nr) then
            i4=nint(xss(land+i))
            i5=nint(xss(ldlw+i))
            write(nsyso,'(7x,a10,2i8,3i10,3x,1p,3e15.6)')&
              name,i1,i2,i3,i4,i5,xss(esz+iaa),xss(esz+ib),&
              xss(lqr+i)
         else
            write(nsyso,'(7x,a10,2i8,i10,23x,1p,3e15.6)')&
              name,i1,i2,i3,xss(esz+iaa),xss(esz+ib),xss(lqr+i)
         endif
      enddo
   endif

   !--print esz block plus gpd cross sections
   iaa=esz+nes
   ib=iaa+nes
   ic=ib+nes
   id=ic+nes
   ig=gpd-1
   do i=1,nes
      if (mod(i,57).eq.1) then
         if (gpd.gt.0) write(nsyso,&
           '(''1''/6x,''i'',5x,''energy'',11x,'' total  '',7x,&
           &''absorption'',5x,''elastic'',8x,''heating'',7x,&
           &''gamma prod''/1x,''------'',3x,''--------------'',&
           &5(3x,''------------''))')
         if (gpd.eq.0) write(nsyso,&
           &'(''1''/6x,''i'',5x,''energy'',11x,'' total  '',7x,&
           &''absorption'',5x,''elastic'',8x,''heating''/&
           &1x,''------'',3x,''--------------'',&
           &4(3x,''------------''))')
      endif
      if (gpd.ne.0) then
         write(nsyso,'(1x,i6,1p,e17.8,7e15.6)') i,xss(esz+i),&
           xss(iaa+i),xss(ib+i),xss(ic+i),xss(id+i),xss(ig+i)
      else
         write(nsyso,'(1x,i6,1p,e17.8,7e15.6)') i,xss(esz+i),&
           xss(iaa+i),xss(ib+i),xss(ic+i),xss(id+i)
      endif
   enddo

   !--print nonelastic cross sections
   list=(ntr+5)/6
   if (list.ne.0) then
      do l=1,list
         na=(l-1)*6+1
         nb=min0(na+5,ntr)
         nc=nb-na+1
         iaa=nes
         ib=1
         j=1
         do n=na,nb
            k=nint(xss(lsig+n)+sig)
            imn(j)=nint(xss(k))
            iaa=min0(iaa,imn(j))
            k=k+1
            loc(j)=k
            imx(j)=imn(j)+nint(xss(k))-1
            ib=max0(ib,imx(j))
            k=iabs(nint(xss(mtr+n)))
            call mtname(k,title(j),izai)
            if (title(j)(1:1).eq.'(') then
               if (izai.eq.0) title(j)(2:2)='g'
               if (izai.eq.1001) title(j)(2:2)='p'
               if (izai.eq.1002) title(j)(2:2)='d'
               if (izai.eq.1003) title(j)(2:2)='t'
               if (izai.eq.2003) title(j)(2:2)='s'
               if (izai.eq.2004) title(j)(2:2)='a'
            endif
            j=j+1
         enddo
         nc2=nc
         nkk=nc
         do m=iaa,ib
            if (mod(m+1-iaa,57).eq.1) then
               write(nsyso,'(''1''/6x,''i'',5x,''energy'',11x,a10,&
                 &6(5x,a10))') (title(ii),ii=1,nc2)
               write(nsyso,'(1x,''------'',3x,''--------------'',&
                 &7(3x,a12))') (dashes,ii=1,nc2)
            endif
            do i=1,nc
               if (m.ge.imn(i).and.m.le.imx(i)) then
                  j=loc(i)+m-imn(i)+1
                  write(kk(i),'(1p,e15.6)') xss(j)
               else
                  write(kk(i),'(15x)')
               endif
            enddo
            write(nsyso,'(1x,i6,1p,e17.8,7a15)')&
              m,xss(esz+m),(kk(i),i=1,nkk)
         enddo
      enddo
   endif

   !--print angular distributions
   nr1=nr+1
   do nn=1,nr1
      n=nn-1
      na=nint(xss(land+n))
      if (na.gt.0) then
         na=na+and
         if (n.eq.0) then
            do itl=1,3
               title(itl)=hlabl1(itl)
            enddo
         else
            j=iabs(nint(xss(mtr+n)))
            call mtname(j,title(1),izai)
            if (title(1)(1:1).eq.'(') then
               if (izai.eq.0) title(1)(2:2)='g'
               if (izai.eq.1001) title(1)(2:2)='p'
               if (izai.eq.1002) title(1)(2:2)='d'
               if (izai.eq.1003) title(1)(2:2)='t'
               if (izai.eq.2003) title(1)(2:2)='s'
               if (izai.eq.2004) title(1)(2:2)='a'
            endif
            if (xss(tyr+n).le.0) then
               do itl=2,3
                  title(itl)=hlabl1(itl)
               enddo
            else
               do itl=2,3
                   title(itl)=hlabl2(itl-1)
               enddo
            endif
         endif
         ne=nint(xss(na))
         nb=na+ne

         !--check on elastic format
         k=nint(xss(nb+1))
         if (k.ge.0) then

            !--equally-probable bins format
            list=(ne+7)/8
            nbin1=nbina+1
            do l=1,list
               iaa=(l-1)*8+1
               ib=min0(ne,iaa+7)
               ic=ib-iaa+1
               j=1
               do m=iaa,ib
                  k=nint(xss(m+nb))
                  if (k.gt.0) k=k+and
                  loc(j)=k
                  j=j+1
               enddo
               write(nsyso,'(''1''///&
                 &22x,''angular distributions for '',a10,&
                 &'' reaction in the '',a10,a10,'' system''//)')&
                 (title(i),i=1,3)
               if (l.eq.1) then
                  write(nsyso,'(6x,''ne ='',i4)') ne
                  write(nsyso,'(1x)')
               endif
               write(nsyso,'(6x,8(4x,a6,a4))') (ek,blank,i=1,ic)
               write(nsyso,'(5x,1p,8e14.5)') (xss(i+na),i=iaa,ib)
               write(nsyso,'(/)')
               nkk=ic
               do j=1,nbin1
                  do m=1,ic
                        if (loc(m).ne.0) then
                        i=loc(m)+j-1
                        write(kk(m),'(1p,e14.5)') xss(i)
                     else
                        write(kk(m),'(14x)')
                     endif
                  enddo
                  write(nsyso,'(1x,i4,8a14)') j,(kk(i),i=1,nkk)
               enddo
            enddo
         else

            !--cummulative-distribution format
            write(nsyso,&
              '(''1''///16x,''angular distributions for '',a10,&
              &'' reaction in the '',a10,a10,'' system''//)')&
             (title(i),i=1,3)
            write(nsyso,'(6x,''ne ='',i4)') ne
            do i=1,ne
               e=xss(na+i)
               k=nint(abs(xss(nb+i)))+and
               int=nint(xss(k))
               np=nint(xss(k+1))
               k=k+1
               write(nsyso,'(/5x,'' incident particle energy ='',&
                 &1p,e14.6,''    int ='',i2,''    np ='',i3)')&
                 e,int,np
               write(nsyso,'(/12x,''cosine'',13x,''pdf'',13x,&
                 &''cdf'',10x,''cosine'',13x,''pdf'',13x,''cdf''/&
                 &2x,6(4x,''------------''))')
               do j=1,np,2
                  if (j.lt.np) then
                     write(nsyso,'(1p,2x,6e16.6)')&
                       xss(k+j),xss(k+np+j),xss(k+2*np+j),&
                       xss(k+j+1),xss(k+np+j+1),xss(k+2*np+j+1)
                  else
                     write(nsyso,'(1p,2x,3e16.6)')&
                       xss(k+j),xss(k+np+j),xss(k+2*np+j)
                  endif
               enddo
            enddo
         endif
      endif
   enddo

   !--print fission nu data, if present
   jnu=0
   if (nu.ne.0) then
      write(nsyso,'(''1''/&
        &4x,''fission nu data for any of the reactions'',&
        &'' mt = 18 or 19,20,21,38''/)')
      l=nu
      j=nint(xss(l))
      kf=j
      idone=0
      do while (idone.eq.0)
         if (kf.lt.0) then
            l=l+1
            jnu=jnu+1
            j=nint(xss(l))
            if (jnu.eq.1) write(nsyso,&
              '('' prompt nu''/'' ---------''/)')
            if (jnu.eq.2) write(nsyso,&
              '('' total nu''/'' --------''/)')
         endif

         !--print fission nu, polynomial form
         if (j.ne.2) then
            write(nsyso,'(12x,''lnu = '',i2,14x,&
              &''polynomial nu, nu = sum(c,i)*(e**(i-1))'')') j
            l=l+1
            j=nint(xss(l))
            write(nsyso,'(12x,''nc ='',i4)') j
            write(nsyso,'(12x,''c(i=1,nc) =   '',&
              &1p,6e14.6/(12x,7e14.6))') (xss(l+i),i=1,j)
            l=l+j

         !--print fission nu, tabular form
         else
            write(nsyso,'(12x,''lnu = '',i2,24x,''tabular nu'')') j
            l=l+1
            j=nint(xss(l))
            write(nsyso,'(12x,''nr ='',i4)') j
            if (j.ne.0) then
               write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
               write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
            endif
            l=l+1
            j=nint(xss(l))
            write(nsyso,'(12x,''ne ='',i4)') j
            write(nsyso,'(12x,''e(i=1,ne) =   '',1p,6e14.6/&
              &(12x,7e14.6))') (xss(l+i),i=1,j)
            l=l+j
            write(nsyso,'(12x,''nu(i=1,ne) =  '',1p,6e14.6/&
              &(12x,7e14.6))') (xss(l+i),i=1,j)
            l=l+j
         endif
         idone=1
         if (jnu.eq.1.and.kf.lt.0) idone=0
      enddo
   endif

   !--print energy distributions for secondary neutrons
   if (nr.ne.0) call acepdd(izai)

   !--print unresolved-range probability tables
   if (iurpt.ne.0) then
      write(nsyso,'(''1''/'' unresolved-range probability tables''/&
        &'' -----------------------------------'')')
      nure=nint(xss(iurpt))
      nurb=nint(xss(iurpt+1))
      lurt=nint(xss(iurpt+2))
      luri=nint(xss(iurpt+3))
      lura=nint(xss(iurpt+4))
      lurf=nint(xss(iurpt+5))
      write(nsyso,'(/''     number of energies: '',i6/&
                    &''     number of bins:     '',i6/&
                    &''     interpolation law:  '',i6/&
                    &''     inelastic reaction: '',i6/&
                    &''     absorption reaction:'',i6)')&
         nure,nurb,lurt,luri,lura
      if (lurf.eq.0) write(nsyso,'(&
                    &''     tables are cross sections'')')
      if (lurf.eq.1) write(nsyso,'(&
                    &''     tables are factors'')')
      do ie=1,nure
         write(nsyso,'(/'' energy='',1p,e14.6)') xss(iurpt+5+ie)
         write(nsyso,&
           '(''   bin     prob           tot          elas'',&
           &''          fiss          capt          heat''/&
           &''   ---   ------   -----------   -----------'',&
           &''  ------------   -----------   -----------'')')
         do ib=1,nurb
            ll=iurpt+5+nure+(ie-1)*6*nurb
            write(nsyso,'(i6,f9.4,1p,5e14.6)') ib,xss(ll+ib),&
              xss(ll+nurb+ib),xss(ll+2*nurb+ib),xss(ll+3*nurb+ib),&
              xss(ll+4*nurb+ib),xss(ll+5*nurb+ib)
         enddo
      enddo
   endif

   !--print delayed neutron data
   if (nud.gt.0) then

      !--delayed nubar
      write(nsyso,'(''1''/'' delayed nubar data''/&
                         &'' ------------------''/)')
      l=nud
      j=nint(xss(l))
      write(nsyso,'(12x,''lnu = '',i2,24x,''tabular nu'')') j
      l=l+1
      j=nint(xss(l))
      write(nsyso,'(12x,''nr ='',i4)') j
      if (j.ne.0) then
         write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
           (nint(xss(l+i)),i=1,j)
         l=l+j
         write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
           (nint(xss(l+i)),i=1,j)
         l=l+j
      endif
      l=l+1
      j=nint(xss(l))
      write(nsyso,'(12x,''ne ='',i4)') j
      write(nsyso,'(12x,''e(i=1,ne) =   '',1p,6e14.6/&
        &(12x,7e14.6))')&
        (xss(l+i),i=1,j)
      l=l+j
      write(nsyso,'(12x,''nu(i=1,ne) =  '',1p,6e14.6/&
        &(12x,7e14.6))')&
        (xss(l+i),i=1,j)
      l=l+j

      !--precursor information
      write(nsyso,'(/'' precursor information''/&
                    &'' ---------------------'')')
      l=dndat
      do i=1,ndnf
         write(nsyso,'(/6x,''decay constant'',i3,'' of'',i3,&
           &'' (per shake) ='',1p,e13.5)') i,ndnf,xss(l)
         write(nsyso,'(/6x,''delayed fraction'')')
         l=l+1
         nn=nint(xss(l))
         write(nsyso,'(12x,''nr ='',i4)') nn
         if (nn.ne.0) then
            write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
              (nint(xss(j+l)),j=1,nn)
            write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
              (nint(xss(j+nn+l)),j=1,nn)
            l=l+2*nn
         endif
         l=l+1
         j=nint(xss(l))
         write(nsyso,'(12x,''ne ='',i4)') j
         write(nsyso,'(12x,&
           &''e(i=1,ne) =   '',1p,6e14.6/(12x,7e14.6))')&
           (xss(l+ii),ii=1,j)
         l=l+j
         write(nsyso,'(12x,&
           &''p(i=1,ne) =   '',1p,6e14.6/(12x,7e14.6))')&
           (xss(l+ii),ii=1,j)
         l=l+j+1
      enddo

      !--precursor energy distributions
      write(nsyso,'(/&
        &'' delayed neutron energy distributions by precursor''/&
        &'' -------------------------------------------------'')')
      l=0
      k=3
      do i=1,ndnf
         nlaw=1
         loct=nint(xss(i-1+ldnd)+dnd-1)
         law=nint(xss(loct+1))
         if (law.eq.4) then
            l=l+1
            if (l.gt.1) write(nsyso,'(/)')
            if (l.gt.1) k=1
            write(nsyso,'(//&
              &''   energy distribution for delayed neutrons from '',&
              &''precursor '',i3,'' of'',i3)') i,ndnf
            write(nsyso,'(/&
              &''    law ='',i2,i5,''st of'',i2,'' laws''/)')&
              law,nlaw,nlaw
            k=k+3
            m=nint(xss(loct+3))
            loct=loct+3
            write(nsyso,'(8x,''probability of law'')')
            write(nsyso,'(12x,''nr ='',i4)') m
            if (m.ne.0) then
               write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
                 (nint(xss(j+loct)),j=1,m)
               write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
                 (nint(xss(j+m+loct)),j=1,m)
               k=k+4
               loct=loct+2*m
            endif
            loct=loct+1
            n=nint(xss(loct))
            write(nsyso,'(12x,''ne ='',i4)') n
            write(nsyso,'(12x,''e(i=1,ne) =   '',1p,6e14.6&
              &/(12x,7e14.6))') (xss(j+loct),j=1,n)
            write(nsyso,'(12x,''p(i=1,ne) =   '',1p,6e14.6&
              &/(12x,7e14.6))') (xss(j+n+loct),j=1,n)
            k=k+3
            loct=loct+1+2*n
            write(nsyso,'(/)')
            write(nsyso,'(8x,''data for law'')')
            k=k+2
            m=nint(xss(loct))
            write(nsyso,'(12x,''nr ='',i4)') m
            if (m.ne.0) then
               k=k+4
               write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
                 (nint(xss(j+loct)),j=1,m)
               write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
                 (nint(xss(j+m+loct)),j=1,m)
               loct=loct+2*m
            endif
            loct=loct+1
            ne=nint(xss(loct))
            write(nsyso,'(12x,''ne ='',i4)') ne
            if (m.eq.0) k=k+1
            do ie=1,ne
               eg=xss(ie+loct)
               loci=nint(xss(ie+ne+loct))+dnd-1
               intt=nint(xss(loci))
               n=nint(xss(loci+1))
               loci=loci+1
               if (ie.ne.1.and.k+6+n.ge.57) then
                  write(nsyso,'(/)')
                  k=1
               endif
               write(nsyso,'(/6x,'' incident energy = '',1p,e14.6,&
                 &''   intt ='',i2,''    np = '',i4//&
                 &1x,&
                 &2(''        energy           pdf           cdf'')/&
                 &1x,&
                 &2(''  ------------  ------------  ------------'')/&
                 &(1x,1p,6e14.6))')&
                 eg,intt,n,(xss(j+loci),xss(j+n+loci),&
                 xss(j+2*n+loci),j=1,n)
               k=k+n+6
            enddo
         endif
      enddo
   endif

   !--old-style gamma production matrix
   if (gpd.ne.0.and.negn.ne.0.and.mtrp.gt.gpd+nes) then
      l=gpd+nes-1
      write(nsyso,'(''1''//&
        &'' gamma production energy distribution''/&
        &'' ------------------------------------''/)')
      write(nsyso,'(/&
        &7x,''neutron energy'',35x,&
        &''equally-probable gamma energies''/&
        &7x,''--------------'',35x,&
        &''-------------------------------''/)')
      do i=1,negn
         write(nsyso,&
           '(/i6,1x,1p,e12.4,3x,0p,10f10.4/(22x,10f10.4))')&
           i,egn(i),(xss(l+k),k=1,nbinp)
         l=l+nbinp
      enddo
   endif

   !--print detailed photon production
   if (ntrp.ne.0) call aceppp(izai,nbina)

   !--restore pointers
   esz=esz+1
   mtr=mtr+1
   tyr=tyr+1
   lsig=lsig+1
   sig=sig+1
   and=and+1
   ldlw=ldlw+1
   dlw=dlw+1
   lqr=lqr+1

   !--print out particle production data, if present.
   if (ntype.ne.0) call acepcp(nbina)

   !--finished
   return
   end subroutine aceprt

   subroutine acepdd(izai)
   !--------------------------------------------------------------------
   ! Print energy and angle-energy distribution data.
   !--------------------------------------------------------------------
   use mainio ! provides nsyso
   use acecm ! provides mtname
   ! externals
   integer::izai
   ! internals
   integer::n,mtrn,nlaw,l,j,i,m,law,iaa,ib,ic,ne,ie
   integer::loci,intt,nd,nn,ip,locj,intmu,nmu,imu,npsx
   integer::lein,locmu,k,ll,intep,npep
   real(kr)::e2,apsx,ein,amu
   integer::loc(8)
   character(3)::ordnl(9)=(/'1st','2nd','3rd','4th','5th','6th',&
     '7th','8th','9th'/)
   character(10)::title(2)

   !--loop over all distribution reactions
   write(nsyso,'(/)')
   do n=1,nr
      mtrn=iabs(nint(xss(mtr+n)))
      call mtname(mtrn,title(1),izai)
      if (title(1)(1:1).eq.'(') then
         if (izai.eq.0) title(1)(2:2)='g'
         if (izai.eq.1001) title(1)(2:2)='p'
         if (izai.eq.1002) title(1)(2:2)='d'
         if (izai.eq.1003) title(1)(2:2)='t'
         if (izai.eq.2003) title(1)(2:2)='s'
         if (izai.eq.2004) title(1)(2:2)='a'
      endif

      !--count number of laws and store locators for each
      nlaw=1
      l=nint(xss(ldlw+n)+dlw)
      loc(1)=l
      do while (nint(xss(l)).ne.0)
         l=nint(xss(l))+dlw
         nlaw=nlaw+1
         loc(nlaw)=l
      enddo
      if (izai.eq.1) title(2)='neutrons'
      if (izai.eq.1001) title(2)='protons'
      if (izai.eq.1002) title(2)='deuterons'
      if (izai.eq.1003) title(2)='tritons'
      if (izai.eq.2003) title(2)='he3s'
      if (izai.eq.2004) title(2)='alphas'
      write(nsyso,'(''1''//&
        &'' energy distribution for secondary '',a9,'' from '',&
        &''reaction '',a10,'' (mt ='',i3,'') with '',i2,&
        &'' law(s)'')') title(2),title(1),mtrn,nlaw

      !--print generalized particle yield, if present.
      l=iabs(nint(xss(tyr+n)))
      if (l.ge.100) then
         l=l-100+dlw
         write(nsyso,'(/8x,''particle yield'')')
         j=nint(xss(l))
         write(nsyso,'(12x,''nr ='',i4)') j
         if (j.ne.0) then
            write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
              (nint(xss(l+i)),i=1,j)
            l=l+j
            write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
              (nint(xss(l+i)),i=1,j)
            l=l+j
         endif
         l=l+1
         j=nint(xss(l))
         write(nsyso,'(12x,''ne ='',i4)') j
         write(nsyso,'(12x,&
           &''e(i=1,ne) =   '',1p,6e14.6/(12x,7e14.6))')&
           (xss(l+i),i=1,j)
         l=l+j
         write(nsyso,'(12x,&
           &''nu(i=1,ne) =  '',1p,6e14.6/(12x,7e14.6))')&
           (xss(l+i),i=1,j)
         l=l+j
      endif

      !--print energy distributions by law
      do m=1,nlaw
         l=loc(m)+1
         law=nint(xss(l))
         l=l+2
         write(nsyso,'(/4x,''law = '',i3,5x,a3,'' of '',i2,&
           &'' laws for reaction '',a10)')&
           law,ordnl(m),nlaw,title(1)
         write(nsyso,'(8x,''probability of law'')')
         j=nint(xss(l))
         write(nsyso,'(12x,''nr ='',i4)') j
         if (j.ne.0) then
            write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
              (nint(xss(l+i)),i=1,j)
            l=l+j
            write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
              (nint(xss(l+i)),i=1,j)
            l=l+j
         endif
         l=l+1
         j=nint(xss(l))
         write(nsyso,'(12x,''ne ='',i4)') j
         write(nsyso,'(12x,&
           &''e(i=1,ne) =   '',1p,6e14.6/(12x,7e14.6))')&
           (xss(l+i),i=1,j)
         l=l+j
         write(nsyso,'(12x,&
           &''p(i=1,ne) =   '',1p,6e14.6/(12x,7e14.6))')&
           (xss(l+i),i=1,j)
         l=l+j+1
         write(nsyso,'(8x,''data for law'')')

         !--branch on laws

         !--fission law
         if (law.eq.18) then
            j=nint(xss(l))

            !--print fission nu, polynomial form
            if (j.ne.2) then
               write(nsyso,'(12x,''lnu =  '',15x,&
                 &''polynomial nu, nu = sum(c,i)*(e**(i-1)))'')')
               l=l+1
               j=nint(xss(l))
               write(nsyso,'(12x,''nc ='',i4)')
               write(nsyso,'(12x,''c(i=1,nc) =   '',1p,6e14.6/&
                 &(12x,7e14.6))') (xss(l+i),i=1,j)

            !--print fission nu, tabular form
            else
               write(nsyso,'(12x,''lnu =  '',25x,''tabular nu'')')
               l=l+1
               j=nint(xss(l))
               write(nsyso,'(12x,''nr ='',i4)') j
               if (j.ne.0) then
                  write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
                    (nint(xss(l+i)),i=1,j)
                  l=l+j
                  write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
                    (nint(xss(l+i)),i=1,j)
                  l=l+j
               endif
               l=l+1
               j=nint(xss(l))
               write(nsyso,'(12x,''ne ='',i4)') j
               write(nsyso,'(12x,''e(i=1,ne) =   '',1p,6e14.6/&
                 &(12x,7e14.6))') (xss(l+i),i=1,j)
               l=l+j
                  write(nsyso,'(12x,''nu(i=1,ne) = '',1p,6e14.6/&
                    &(12x,7e14.6))') (xss(l+i),i=1,j)
            endif

         !--law 1
         else if (law.eq.1) then
            j=nint(xss(l))
            write(nsyso,'(12x,''nr ='',i4)') j
            if (j.ne.0) then
               write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
               write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
            endif
            l=l+1
            iaa=nint(xss(l))
            ib=l+iaa+1
            ic=nint(xss(ib))
            write(nsyso,'(12x,''ne ='',i4)') iaa
            do j=1,iaa
               write(nsyso,'(12x,''e('',i3,'' ) = '',1p,e14.6)')&
                 j,xss(l+j)
               write(nsyso,'(12x,''eout(i=1,'',i3,'') ='',13x,&
                 &1p,5e14.6/(12x,7e14.6))') ic,(xss(ib+i),i=1,ic)
               ib=ib+ic
            enddo

         !--law 3
         else if (law.eq.3) then
            write(nsyso,'(12x,''eout = c*(e-ec)   ec ='',1p,e14.6,&
              &5x,''c ='',e14.6)') xss(l),xss(l+1)

         !--law 4
         else if (law.eq.4) then
            j=nint(xss(l))
            write(nsyso,'(12x,''nr ='',i4)') j
            if (j.ne.0) then
               write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
               write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
            endif
            l=l+1
            ne=nint(xss(l))
            write(nsyso,'(12x,''ne ='',i4)') ne
            do ie=1,ne
               e2=xss(ie+l)
               loci=nint(xss(ie+ne+l)+dlw)
               intt=nint(xss(loci))
               nd=nint(xss(loci)/10)
               nn=nint(xss(loci+1))
               loci=loci+1
               write(nsyso,'(/6x,'' incident energy = '',1p,e14.6,&
                 &''   intt ='',i2,''    nd = '',i4,''    np = '',i4//&
                 &1x,&
                 &2(''        energy           pdf           cdf'')/&
                 &1x,&
                 &2(''  ------------  ------------  ------------'')/&
                 &(1x,1p,6e14.6))')&
                 e2,intt,nd,nn,(xss(j+loci),xss(j+nn+loci),&
                 xss(j+2*nn+loci),j=1,nn)
            enddo

         !--law 5
         else if (law.eq.5) then
            j=nint(xss(l))
            write(nsyso,'(12x,''nr ='',i4)') j
            if (j.ne.0) then
               write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
               write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
            endif
            l=l+1
            j=nint(xss(l))
            write(nsyso,'(12x,''ne ='',i4)') j
            write(nsyso,'(12x,''e(i=1,ne) = '',16x,1p,5e14.6/&
              &(12x,7e14.6))') (xss(l+i),i=1,j)
            l=l+j
            write(nsyso,'(12x,''theta(i=1,ne) ='',13x,1p,5e13.5/&
              &(12x,7e14.6))') (xss(l+i),i=1,j)
            l=l+j+1
            j=nint(xss(l))
            write(nsyso,'(12x,''nx ='',i4)') j
            write(nsyso,'(12x,''x(i=1,nx) =   '',1p,6e14.6/&
              &(12x,7e14.6))') (xss(l+i),i=1,j)

         !--laws 7 and 9
         else if (law.eq.7.or.law.eq.9) then
            j=nint(xss(l))
            write(nsyso,'(12x,''nr ='',i4)') j
            if (j.ne.0) then
               write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
               write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
            endif
            l=l+1
            j=nint(xss(l))
            write(nsyso,'(12x,''ne ='',i4)') j
            write(nsyso,'(12x,''e(i=1,ne) = '',16x,1p,5e14.6/&
              &(12x,7e14.6))') (xss(l+i),i=1,j)
            l=l+j
            write(nsyso,'(12x,''theta(i=1,ne) = '',12x,1p,5e14.6/&
              &(12x,7e14.6))') (xss(l+i),i=1,j)
            l=l+j+1
            write(nsyso,'(12x,''u =  '',9x,1p,e14.6)') xss(l)

         !--law 10
         else if (law.eq.10) then
            write(nsyso,'(12x,''rejection scheme constants'',6x,&
              &''b = '',1p,e14.6,5x,''m = '',e14.6,5x,''u = '',&
              &e14.6)') xss(l),xss(l+1),xss(l+2)

         !--law 11
         else if (law.eq.11) then
            j=nint(xss(l))
            write(nsyso,'(12x,''nr ='',i4)') j
            if (j.ne.0) then
               write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
               write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
            endif
            l=l+1
            j=nint(xss(l))
            write(nsyso,'(12x,''ne ='',i4)') j
            write(nsyso,'(12x,''e(i=1,ne) = '',16x,1p,5e14.6/&
              &(12x,7e14.6))') (xss(l+i),i=1,j)
            l=l+j
            write(nsyso,'(12x,''    a(i=1,ne) = '',12x,1p,5e14.6/&
              &(12x,7e14.6))') (xss(l+i),i=1,j)
            l=l+j+1
            j=nint(xss(l))
            write(nsyso,'(12x,''nr ='',i4)') j
            if (j.ne.0) then
               write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
               write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
            endif
            l=l+1
            j=nint(xss(l))
            write(nsyso,'(12x,''ne ='',i4)') j
            write(nsyso,'(12x,''e(i=1,ne) ='',17x,1p,5e14.6/&
              &(12x,7e14.6))') (xss(l+i),i=1,j)
            l=l+j
            write(nsyso,'(12x,''    b(i=1,ne) = '',12x,1p,5e14.6/&
              &(12x,7e14.6))') (xss(l+i),i=1,j)
            l=l+j+1
            write(nsyso,'(12x,''u =  '',9x,1p,e14.6)') xss(l)

         !--law 44 (coupled energy-angle distribution)
         else if (law.eq.44) then
            j=nint(xss(l))
            write(nsyso,'(12x,''nr ='',i4)') j
            if (j.ne.0) then
               write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
               write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
            endif
            l=l+1
            ne=nint(xss(l))
            write(nsyso,'(12x,''ne ='',i4)') ne
            do ie=1,ne
               e2=xss(ie+l)
               loci=nint(xss(ie+ne+l)+dlw)
               intt=mod(nint(xss(loci)),10)
               nd=nint(xss(loci)/10)
               nn=nint(xss(loci+1))
               loci=loci+1
               write(nsyso,'(/6x,'' incident energy = '',1p,e14.6,&
                 &''   intt ='',i2,''    nd = '',i4,''    np = '',i4//&
                 &''         energy           pdf           cdf'',&
                 &''             r             a''/&
                 &''   ------------  ------------  ------------'',&
                 &''  ------------  ------------''/&
                 &(1x,1p,5e14.6))')&
                 e2,intt,nd,nn,(xss(j+loci),xss(j+nn+loci),&
                 xss(j+2*nn+loci),xss(j+3*nn+loci),xss(j+4*nn+loci),&
                 j=1,nn)
            enddo

         !--law 61 (energy-angle cummulative distributions)
         else if (law.eq.61) then
            j=nint(xss(l))
            write(nsyso,'(12x,''nr ='',i4)') j
            if (j.ne.0) then
               write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
               write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
            endif
            l=l+1
            ne=nint(xss(l))
            write(nsyso,'(12x,''ne ='',i4)') ne
            do ie=1,ne
               e2=xss(ie+l)
               loci=nint(xss(ie+ne+l)+dlw)
               intt=mod(nint(xss(loci)),10)
               nd=nint(xss(loci)/10)
               nn=nint(xss(loci+1))
               loci=loci+1
               write(nsyso,'(/6x,'' incident energy = '',1p,e14.6,&
                 &''   intt ='',i2,''    nd = '',i4,''    np = '',&
                 &i4)') e2,intt,nd,nn
               do ip=1,nn
                  locj=nint(xss(ip+3*nn+loci)+dlw)
                  intmu=nint(xss(locj))
                  nmu=nint(xss(locj+1))
                  write(nsyso,'(/&
                    &6x,'' secondary energy = '',1p,e14.6/&
                    &6x,''              pdf = '',e14.6/&
                    &6x,''              cdf = '',e14.6/&
                    &6x,''            intmu = '',i8/&
                    &6x,''              nmu = '',i8/&
                    &''         cosine           pdf           cdf'',&
                    &''        cosine           pdf           cdf''/&
                    &''   ------------  ------------  ------------'',&
                    &''  ------------  ------------  ------------'')')&
                    xss(ip+loci),xss(ip+nn+loci),xss(ip+2*nn+loci),&
                    intmu,nmu
                  do imu=1,nmu,2
                     if (imu.eq.nmu) then
                        write(nsyso,'(1x,1p,3e14.6)')&
                          xss(locj+1+imu),xss(locj+1+nmu+imu),&
                          xss(locj+1+2*nmu+imu)
                     else
                        write(nsyso,'(1x,1p,6e14.6)')&
                          xss(locj+1+imu),xss(locj+1+nmu+imu),&
                          xss(locj+1+2*nmu+imu),xss(locj+1+imu+1),&
                          xss(locj+1+nmu+imu+1),&
                          xss(locj+1+2*nmu+imu+1)
                     endif
                  enddo
               enddo
            enddo

         !--law 66 (phase-space)
         else if (law.eq.66) then
            npsx=nint(xss(l))
            apsx=xss(l+1)
            intt=nint(xss(l+2))
            nn=nint(xss(l+3))
            loci=l+3
            write(nsyso,'(12x,''npsx ='',i2/12x,''apsx ='',f10.4/&
              &12x,''intt ='',i2)') npsx,apsx,intt
            write(nsyso,'(&
              &1x,2(''        energy           pdf           cdf'')/&
              &1x,2(''  ------------  ------------  ------------'')/&
              &(1x,1p,6e14.6))')&
              (xss(j+loci),xss(j+nn+loci),xss(j+2*nn+loci),j=1,nn)

         !--law 67 (angle-energy)
         else if (law.eq.67) then
            j=nint(xss(l))
            write(nsyso,'(12x,''nr ='',i4)') j
            if (j.ne.0) then
               write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
               write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
                 (nint(xss(l+i)),i=1,j)
               l=l+j
            endif
            l=l+1
            ne=nint(xss(l))
            write(nsyso,'(12x,''ne ='',i4)') ne
            do ie=1,ne
               ein=xss(l+ie)
               lein=nint(xss(ie+ne+l)+dlw)
               intmu=nint(xss(lein))
               nmu=nint(xss(lein+1))
               write(nsyso,'(/4x,''ein ='',f10.5,4x,''intmu ='',i2,&
                 &4x,''nmu ='',i3)') ein,intmu,nmu
               locmu=lein+1
               do k=1,nmu
                  amu=xss(locmu+k)
                  ll=nint(xss(locmu+nmu+k)+dlw)
                  intep=nint(xss(ll))
                  ll=ll+1
                  npep=nint(xss(ll))
                  write(nsyso,'(/5x,''mu ='',f10.5,4x,''intep ='',&
                    &i2,4x,''npep ='',i3)') amu,intep,npep
                  write(nsyso,'(&
                    &1x,&
                    &2(''        energy           pdf           cdf'')&
                    &/1x,&
                    &2(''  ------------  ------------  ------------'')&
                    &/(1x,1p,6e14.6))')&
                    (xss(ll+j),xss(ll+npep+j),xss(ll+2*npep+j),&
                    j=1,npep)
               enddo
            enddo
         endif

      !--close loop over number of laws
      enddo
      write(nsyso,'(/)')
      write(nsyso,'(/)')

   !--close loop over reactions
   enddo

   return
   end subroutine acepdd

   subroutine aceppp(izai,nbina)
   !--------------------------------------------------------------------
   ! Print detailed photon production.
   !--------------------------------------------------------------------
   use mainio ! provides nsyso
   use acecm ! provides mtname
   ! externals
   integer::izai,nbina
   ! internals
   integer::ntyp,k,nindx,ipy,nesp,i,loct,mftype,loc1,je,nje
   integer::law,mtmult,m,n,itrp,kl,mtnow,j,list,lst,na,nb,nc
   integer::nkk,ielo,iehi,mtn,ke,ie,l,iflag,ngg,ne,nbin1
   integer::iaa,ib,ic,ii,mti,mtrn,nlaw,loci,intt,nd,nyp
   real(kr)::egamma,eg
   integer::loc(8)
   integer,dimension(:),allocatable::indx
   character(10)::title(16)
   character(15)::kk(40)
   character(6)::ek='energy'
   character(6)::blank='      '
   character(12)::dashes='-----------'

   ntyp=1
   if (ntyp.eq.0) then
      k=26
   else
      write(nsyso,'(''1''//'' ********** photons **********'')')
      k=4
   endif
   write(nsyso,'(/7x,''mt'',4x,''lsigp'',3x,''landp'',3x,&
     &''ldlwp'',2x,''mftype'',2x,''mtmult'',7x,''emin'',10x,&
     &''emax'',9x,''egamma'',8x,''yield'',8x,''lp''/4x,&
     &''------'',3x,''-----'',3x,''-----'',2x,''------'',2x,&
     &''------'',2x,''------'',2x,4(2x,''------------''),3x,&
     &''--'')')
   allocate(indx(ntrp))
   nindx=0
   ipy=0
   nesp=nes

   !--print photon reaction descriptions
   do i=1,ntrp
      loct=nint(xss(i-1+lsigp)+sigp-1)
      mftype=nint(xss(loct))
      loct=loct+1
      loc1=nint(xss(i-1+ldlwp)+dlwp-1)
      k=k+1
      if (mod(k,57).eq.1) then
         write(nsyso,'(''1''/&
     &        7x,''mt'',4x,''lsigp'',3x,''landp'',3x,''ldlwp'',2x,&
     &        ''mftype'',2x,''mtmult'',7x,''emin'',10x,''emax'',9x,&
     &        ''egamma'',8x,''yield'',8x,''lp''/4x,''------'',3x,&
     &        ''-----'',3x,''-----'',3x,''-----'',2x,''------'',2x,&
     &        ''------'',2x,4(2x,''------------''),3x,''--'')')
         k=3
      endif
      if (mftype.eq.13) then
         nindx=nindx+1
         indx(nindx)=i
         je=nint(xss(loct))
         nje=nint(xss(loct+1))
         law=nint(xss(loc1+1))
         if (law.eq.2) then
            write(nsyso,&
              '(4x,i6,2x,i6,2x,i6,1x,i7,2x,i6,10x,1p,3e14.6,16x,i3)')&
              nint(xss(i-1+mtrp)),nint(xss(i-1+lsigp)),&
              nint(xss(i-1+landp)),nint(xss(i-1+ldlwp)),mftype,&
              xss(je+esz),xss(je+nje-1+esz),xss(loc1+10),&
              nint(xss(loc1+9))
         else
            write(nsyso,&
              '(4x,i6,2x,i6,2x,i6,1x,i7,2x,i6,10x,1p,2e14.6,30x,i3)')&
              nint(xss(i-1+mtrp)),nint(xss(i-1+lsigp)),&
              nint(xss(i-1+landp)),nint(xss(i-1+ldlwp)),mftype,&
              xss(je+esz),xss(je+nje-1+esz)
         endif
      else
         mtmult=nint(xss(loct))
         loct=loct+1
         m=nint(xss(loct))
         loct=loct+1
         n=nint(xss(loct))
         law=nint(xss(loc1+1))
         if (law.ne.2.or.(m.ne.0.or.n.gt.2).or.&
           xss(loct+3).ne.xss(loct+4)) then
            nindx=nindx+1
            indx(nindx)=i
            ipy=ipy+1
            loct=loct+2*m
            if (law.eq.2) then
               write(nsyso,'(4x,i6,2x,i6,2x,i6,1x,i7,2x,i6,&
                 &2x,i6,2x,1p,3e14.6,16x,i3)')&
                 nint(xss(i-1+mtrp)),nint(xss(i-1+lsigp)),&
                 nint(xss(i-1+landp)),nint(xss(i-1+ldlwp)),&
                 mftype,mtmult,xss(loct+1),xss(loct+n),&
                 xss(loc1+10),nint(xss(loc1+9))
            else
               write(nsyso,'(4x,i6,2x,i6,2x,i6,1x,i7,2x,i6,&
                 &2x,i6,2x,1p,4e14.6,2x,i3)')&
                 nint(xss(i-1+mtrp)),nint(xss(i-1+lsigp)),&
                 nint(xss(i-1+landp)),nint(xss(i-1+ldlwp)),&
                 mftype,mtmult,xss(loct+1),xss(loct+n)
            endif
         else
            write(nsyso,'(4x,i6,2x,i6,2x,i6,1x,i7,2x,i6,2x,i6,2x,&
              &1p,4e14.6,2x,i3)') nint(xss(i-1+mtrp)),&
              nint(xss(i-1+lsigp)),nint(xss(i-1+landp)),&
              nint(xss(i-1+ldlwp)),mftype,mtmult,&
              xss(loct+1),xss(loct+2),xss(loc1+10),xss(loct+3),&
              nint(xss(loc1+9))
         endif
      endif
   enddo

   !--print energy-dependent photon production yields
   if (nindx.ne.0) then
      if (ipy.gt.0) write(nsyso,'(''1''/&
        &'' photon production yields''/&
        &'' ------------------------'')')
      k=0
      do i=1,nindx
         itrp=indx(i)
         loct=nint(xss(itrp-1+lsigp)+sigp-1)
         mftype=nint(xss(loct))
         loct=loct+1
         if (mftype.eq.13.and.k.eq.0) k=i
         if (mftype.eq.13) kl=i
         if (mftype.ne.13) then
            mtnow=nint(xss(itrp-1+mtrp))
            loc1=nint(xss(itrp-1+ldlwp)+dlwp-1)
            law=nint(xss(loc1+1))
            egamma=0
            if (law.eq.2) egamma=xss(loc1+10)
            mtmult=nint(xss(loct))
            if (law.eq.2) then
               write(nsyso,'(//'' file type = '',i2,6x,''   mt = '',&
                 &i6,''   egamma = '', 1p,e12.5,''   mtmult = '',&
                 &i3)') mftype,mtnow,egamma,mtmult
            else
               write(nsyso,'(//'' file type = '',i2,6x,''   mt = '',&
                 &i6,''   continuum gammas     mtmult = '',i3)')&
                 mftype,mtnow,mtmult
            endif
            loct=loct+1
            m=nint(xss(loct))
            write(nsyso,'(12x,''nr ='',i4)') m
            if (m.ne.0) then
               write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
                 (nint(xss(j+loct)),j=1,m)
               write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
                 (nint(xss(j+m+loct)),j=1,m)
               loct=loct+2*m
            endif
            loct=loct+1
            n=nint(xss(loct))
            write(nsyso,'(12x,''ne ='',i4)') n
            write(nsyso,&
               '(12x,''e(i=1,ne) ='',17x,1p,5e14.6/(12x,7e14.6))')&
               (xss(j+loct),j=1,n)
            write(nsyso,'(12x,''yield(i=1,ne) ='',13x,1p,5e14.6/&
               &(12x,7e14.6))') (xss(j+n+loct),j=1,n)
         endif
      enddo

      !--print photon production cross sections
      if (mftype.ne.12.and.k.ne.0) then
         list=(kl-k+7)/7
         write(nsyso,'(''1''/&
           &'' photon production cross sections''/&
           &'' --------------------------------''/)')
         do lst=1,list
            na=(lst-1)*6+1
            nb=min0(na+5,kl-k+1)
            nc=nb-na+1
            nkk=nc
            ielo=2*nesp
            iehi=-ielo
            do n=na,nb
               j=indx(k+n-1)
               mtn=nint(xss(j-1+mtrp))
               loct=nint(xss(j-1+lsigp)+sigp-1)
               je=nint(xss(loct+1))
               if (je.lt.ielo) ielo=je
               nje=nint(xss(loct+2))
               if ((je+nje-1).gt.iehi) iehi=je+nje-1
               write(title(n-na+1),'(''  mt'',i6)') mtn
            enddo
            ke=1
            do ie=ielo,iehi
               if (mod(ke,57).eq.1) then
                  if (ke.gt.1.or.lst.gt.1) then
                     write(nsyso,'(''1'')')
                  endif
                  write(nsyso,&
                     '(6x,''i'',5x,''energy'',9x,a10,5(5x,a10))')&
                     (title(l),l=1,nc)
                  write(nsyso,'(1x,''------'',3x,''------------'',&
                     &7(3x,a12))') (dashes,l=1,nc)
               endif
               do m=na,nb
                  j=indx(k+m-1)
                  loct=nint(xss(j-1+lsigp)+sigp-1)
                  je=nint(xss(loct+1))
                  nje=nint(xss(loct+2))
                  l=m-na+1
                  if (je.le.ie.and.je+nje-1.ge.ie) then
                     loct=loct+3+ie-je
                     write(kk(l),'(1p,e15.6)') xss(loct)
                  else
                     write(kk(l),'(15x)')
                  endif
               enddo
               ke=ke+1
               write(nsyso,'(1x,i6,1p,e15.6,6a15)')&
                 ie,xss(ie+esz),(kk(l),l=1,nkk)
            enddo
         enddo
         deallocate(indx)
      endif
   endif

   !--print photon angular distributions
   iflag=0
   do i=1,ntrp
      if (nint(xss(i+landp-1)).ne.0) iflag=1
   enddo
   if (iflag.eq.0) then
      write(nsyso,'(''1''/&
        &'' photon angular distributions''/&
        &'' ----------------------------''//&
        &'' all are isotropic.'')')
   else
      do i=1,ntrp
         ngg=nint(xss(i-1+mtrp))
         na=nint(xss(i+landp-1))
         if (na.gt.0) then
            na=na+andp-1
            ne=nint(xss(na))
            list=(ne+7)/8
            nb=na+ne
            nbin1=nbina+1
            do l=1,list
               iaa=(l-1)*8+1
               ib=min0(ne,iaa+7)
               ic=ib-iaa+1
               j=1
               do m=iaa,ib
                  k=nint(xss(m+nb))
                  if (k.gt.0) k=k+andp-1
                  loc(j)=k
                  j=j+1
               enddo
               write(nsyso,'(''1''///22x,&
                 &''angular distributions for photon'',i10//)') ngg
               write(nsyso,'(6x,8(4x,a6,a4))') (ek,blank,ii=1,ic)
               write(nsyso,'(5x,1p,8e14.5)') (xss(ii+na),ii=iaa,ib)
               write(nsyso,'(/)')
               nkk=ic
               do j=1,nbin1
                  do m=1,ic
                     if (loc(m).ne.0) then
                        ii=loc(m)+j-1
                        write(kk(m),'(1p,e14.5)') xss(ii)
                     else
                        write(kk(m),'(''               '')')
                     endif
                  enddo
                  write(nsyso,'(1x,i4,8a14)') j,(kk(ii),ii=1,nkk)
               enddo
            enddo
         endif
      enddo
   endif

   !--print tabulated photon energy distributions
   if (ntrp.ne.0) then
      write(nsyso,'(''1''/&
        &'' photon energy distributions''/&
        &'' ---------------------------'')')
      l=0
      k=3
      do i=1,ntrp
         mti=nint(xss(i-1+mtrp)/1000)
         mtrn=nint(xss(i-1+mtrp))
         call mtname(mti,title(1),izai)
         if (title(1)(1:1).eq.'(') then
            if (izai.eq.0) title(1)(2:2)='g'
            if (izai.eq.1001) title(1)(2:2)='p'
            if (izai.eq.1002) title(1)(2:2)='d'
            if (izai.eq.1003) title(1)(2:2)='t'
            if (izai.eq.2003) title(1)(2:2)='s'
            if (izai.eq.2004) title(1)(2:2)='a'
         endif
         nlaw=1
         loct=nint(xss(i-1+ldlwp)+dlwp-1)
         law=nint(xss(loct+1))
         if (law.eq.4.or.law.eq.44) then
            l=l+1
            if (l.gt.1) write(nsyso,'(/)')
            if (l.gt.1) k=1
            write(nsyso,'(//&
              &'' energy distribution for secondary photons from '',&
              &''reaction '',a10,'' (mt='',i6,'') with '',i1,&
              &'' law(s)'')') title(1),mtrn,nlaw
            write(nsyso,'(/&
              &''    law ='',i2,i5,''st of'',i2,'' laws''/)')&
              law,nlaw,nlaw
            k=k+3
            m=nint(xss(loct+3))
            loct=loct+3
            write(nsyso,'(8x,''probability of law'')')
            write(nsyso,'(/)')
            write(nsyso,'(12x,''nr ='',i4)') m
            if (m.ne.0) then
               write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
                 (nint(xss(j+loct)),j=1,m)
               write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
                 (nint(xss(j+m+loct)),j=1,m)
               k=k+4
               loct=loct+2*m
            endif
            loct=loct+1
            n=nint(xss(loct))
            write(nsyso,'(12x,''ne ='',i4)') n
            write(nsyso,'(12x,''e(i=1,ne) =   '',1p,6e14.6&
              &/(12x,7e14.6))') (xss(j+loct),j=1,n)
            write(nsyso,'(12x,''p(i=1,ne) =   '',1p,6e14.6&
              &/(12x,7e14.6))') (xss(j+n+loct),j=1,n)
            k=k+3
            loct=loct+1+2*n
            write(nsyso,'(/)')
            write(nsyso,'(8x,''data for law'')')
            write(nsyso,'(/)')
            k=k+2
            m=nint(xss(loct))
            write(nsyso,'(12x,''nr ='',i4)') m
            if (m.ne.0) then
               k=k+4
               write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
                 (nint(xss(j+loct)),j=1,m)
               write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
                 (nint(xss(j+m+loct)),j=1,m)
               loct=loct+2*m
            endif
            loct=loct+1
            ne=nint(xss(loct))
            if (m.ne.0) then
               write(nsyso,'(12x,''ne ='',i4)') ne
            endif
            if (m.eq.0) write(nsyso,&
               '(8x,''number of incident energies = '',i3)') ne
            if (m.eq.0) k=k+1
            do ie=1,ne
               eg=xss(ie+loct)
               loci=nint(xss(ie+ne+loct))+dlwp-1
               intt=mod(nint(xss(loci)),10)
               nd=nint(xss(loci)/10)
               n=nint(xss(loci+1))
               loci=loci+1
               if (ie.ne.1.and.k+6+n.ge.57) then
                  write(nsyso,'(/)')
                  k=1
               endif
               write(nsyso,'(/6x,'' incident energy = '',1p,e14.6,&
                 &''   intt ='',i2,''    nd = '',i4,''    np = '',i4//&
                 &1x,&
                 &2(''        energy           pdf           cdf'')/&
                 &1x,&
                 &2(''  ------------  ------------  ------------'')/&
                 &(1x,1p,6e14.6))')&
                 eg,intt,nd,n,(xss(j+loci),xss(j+n+loci),&
                 xss(j+2*n+loci),j=1,n)
               k=k+n+6
            enddo
         endif
      enddo
   endif

   !--print yp block
   nyp=nint(xss(yp))
   if (nyp.ne.0) then
      write(nsyso,'(//&
        &'' neutron mts needed as multipliers for photon yields''/&
        &'' ---------------------------------------------------''//&
        &(6x,12i6))') (nint(xss(i+yp)),i=1,nyp)
   endif
   if (allocated(indx)) deallocate(indx)
   return
   end subroutine aceppp

   subroutine acepcp(nbina)
   !--------------------------------------------------------------------
   ! Print particle production information.
   !--------------------------------------------------------------------
   use mainio ! provides nsyso
   ! externals
   integer::nbina
   ! internals
   integer::i,iaa,naa,imt,ie,mt,l1,l2,nrint,ii,ne,ll,na,nb
   integer::m,k,int,np,j,list,nbin1,l,ib,ic,nkk,law,l3
   integer::loci,intt,nn,nd,ip,locj,intmu,nmu,imu,npsx
   integer::intep,npep,nyh,locmu
   integer::ipt,mtrh,nmtr
   integer hpd,tyrh,lsigh,sigh,landh,andh,ldlwh,dlwh,yh
   real(kr)::e,xs,heat,e2,apsx,amu
   integer::loc(8)
   character(15)::kk(40)
   character(6)::blank='      '
   character(6)::ek='energy'

   !--loop over production types
   do i=1,ntype
      ipt=nint(xss(ptype+i-1))
      hpd=nint(xss(ploct+10*(i-1)))
      iaa=nint(xss(hpd))
      naa=nint(xss(hpd+1))
      if (ipt.eq.1) write(nsyso,'(''1''//&
        &'' ********** neutron production **********'')')
      if (ipt.eq.9) write(nsyso,'(''1''//&
        &'' ********** proton production **********'')')
      if (ipt.eq.31) write(nsyso,'(''1''//&
        &'' ********** deuteron production **********'')')
      if (ipt.eq.32) write(nsyso,'(''1''//&
        &'' ********** triton production **********'')')
      if (ipt.eq.33) write(nsyso,'(''1''//&
        &'' ********** he-3 production **********'')')
      if (ipt.eq.34) write(nsyso,'(''1''//&
        &'' ********** alpha production **********'')')

      !--print reaction information for this type
      mtrh=nint(xss(ploct+10*(i-1)+1))
      nmtr=nint(xss(ntro+i-1))
      tyrh=nint(xss(ploct+10*(i-1)+2))
      lsigh=nint(xss(ploct+10*(i-1)+3))
      sigh=nint(xss(ploct+10*(i-1)+4))
      landh=nint(xss(ploct+10*(i-1)+5))
      andh=nint(xss(ploct+10*(i-1)+6))
      ldlwh=nint(xss(ploct+10*(i-1)+7))
      dlwh=nint(xss(ploct+10*(i-1)+8))
      yh=nint(xss(ploct+10*(i-1)+9))
      write(nsyso,'(/''              nmtr ='',i7/&
                   & ''              mtrh ='',i7/&
                   & ''              tyrh ='',i7/&
                   & ''             lsigh ='',i7/&
                   & ''              sigh ='',i7/&
                   & ''             landh ='',i7/&
                   & ''              andh ='',i7/&
                   & ''             ldlwh ='',i7/&
                   & ''              dlwh ='',i7/&
                   & ''                yh ='',i7)')&
        nmtr,mtrh,tyrh,lsigh,sigh,landh,andh,ldlwh,dlwh,yh
      write(nsyso,'(/&
        &''   mtrh   tyrh   lsigh   landh   ldlwh   mtmult''/&
        &''   ----   ----   -----   -----   -----   ------'')')
      do imt=1,nmtr
         write(nsyso,'(i7,i7,i8,i8,i8,i9)') nint(xss(mtrh+imt-1)),&
           nint(xss(tyrh+imt-1)),nint(xss(lsigh+imt-1)),&
           nint(xss(landh+imt-1)),nint(xss(ldlwh+imt-1)),&
           nint(xss(yh+imt))
      enddo

      !--print production cross section and heating
      write(nsyso,'(/'' production cross section''/&
        &4x,''ie'',9x,''energy'',13x,''xs'',8x,''heating''/&
        &2x,''----'',3(2x,''-------------''))')
      do ie=1,naa
         e=xss(esz+iaa+ie-2)
         xs=xss(hpd+1+ie)
         heat=xss(hpd+1+naa+ie)
         write(nsyso,'(i6,1p,3e15.6)') iaa+ie-1,e,xs,heat
      enddo

      !--print production distributions
      if (ipt.eq.1)&
        write(nsyso,'(''1''/'' neutron production distributions'')')
      if (ipt.eq.9)&
        write(nsyso,'(''1''/'' proton production distributions'')')
      if (ipt.eq.31)&
        write(nsyso,'(''1''/'' deuteron production distributions'')')
      if (ipt.eq.32)&
        write(nsyso,'(''1''/'' triton production distributions'')')
      if (ipt.eq.33)&
        write(nsyso,'(''1''/'' he-3 production distributions'')')
      if (ipt.eq.34)&
        write(nsyso,'(''1''/'' alpha production distributions'')')

      !--loop over mts contributing to this production
      do imt=1,nmtr
         mt=nint(xss(mtrh+imt-1))
         write(nsyso,'(/'' contributions from mt='',i3,&
           &'' ------->'')') mt
         l1=nint(xss(lsigh+imt-1))
         l2=sigh+l1+1
         nrint=nint(xss(l2))
         write(nsyso,'(4x,''nr ='',i4)') nrint
         if (nrint.ne.0) then
            write(nsyso,'(4x,''nbt(i=1,nr) = '',20i5)')&
              (nint(xss(l2+ii)),ii=1,nrint)
            l2=l2+nrint
            write(nsyso,'(4x,''int(i=1,nr) = '',20i5)')&
              (nint(xss(l2+ii)),ii=1,nrint)
            l2=l2+nrint
         endif
         l2=l2+1
         ne=nint(xss(l2))
         write(nsyso,'(4x,''ne ='',i4)') ne
         write(nsyso,'(/&
           &''              energy         yield''/&
           &6x,2(2x,''------------''))')
         do ii=1,ne
            write(nsyso,'(6x,1p,2e14.6)') xss(l2+ii),xss(l2+ne+ii)
         enddo
         ll=nint(xss(landh+imt-1))

         !--angular distributions
         if (ll.gt.0) then
            write(nsyso,'(/''   angular distribution:'')')
            ll=ll+andh-1
            ne=nint(xss(ll))
            write(nsyso,'(/''      ne ='',i4)') ne
            na=ll
            nb=na+ne
            m=nint(xss(nb+1))

            !--cummulative form
            if (m.lt.0) then
               do ie=1,ne
                  e=xss(na+ie)
                  k=nint(abs(xss(nb+ie)))+andh-1
                  int=nint(xss(k))
                  np=nint(xss(k+1))
                  k=k+1
                  write(nsyso,'(/5x,'' incident particle energy ='',&
                    &1p,e14.6,''    int ='',i2,''    np ='',i4)')&
                    e,int,np
                  write(nsyso,'(/12x,''cosine'',13x,''pdf'',13x,&
                    &''cdf'',10x,''cosine'',13x,''pdf'',13x,''cdf''/&
                    &2x,6(4x,''------------''))')
                  do j=1,np,2
                     if (j.lt.np) then
                        write(nsyso,'(1p,2x,6e16.6)')&
                          xss(k+j),xss(k+np+j),xss(k+2*np+j),&
                          xss(k+j+1),xss(k+np+j+1),xss(k+2*np+j+1)
                     else
                        write(nsyso,'(1p,2x,3e16.6)')&
                          xss(k+j),xss(k+np+j),xss(k+2*np+j)
                     endif
                  enddo
               enddo

            !--equally probable bins
            else
               list=(ne+7)/8
               nbin1=nbina+1
               do l=1,list
                  iaa=(l-1)*8+1
                  ib=min0(ne,iaa+7)
                  ic=ib-iaa+1
                  j=1
                  do m=iaa,ib
                     k=nint(xss(m+nb))
                     if (k.gt.0) k=k+and
                     loc(j)=k
                     j=j+1
                  enddo
                  write(nsyso,'(6x,8(4x,a6,a4))') (ek,blank,ii=1,ic)
                  write(nsyso,'(5x,1p,8e14.5)') (xss(i+na),ii=iaa,ib)
                  write(nsyso,'(/)')
                  nkk=ic
                  do j=1,nbin1
                     do m=1,ic
                        if (loc(m).ne.0) then
                           ii=loc(m)+j-1
                           write(kk(m),'(1p,e14.5)') xss(ii)
                        else
                           write(kk(m),'(14x)')
                        endif
                     enddo
                     write(nsyso,'(1x,i4,8a14)') j,(kk(ii),ii=1,nkk)
                  enddo
               enddo
            endif
         endif

         !--energy-angle distributions
         l1=nint(xss(ldlwh+imt-1))
         if (l1.gt.0) then
            l2=dlwh+l1-1
            law=nint(xss(l2+1))
            write(nsyso,&
              '(/''   distributions (law'',i2,''):'')') law
            l3=dlwh+nint(xss(l2+2))-1

            !--law=4
            if (law.eq.4) then
               nrint=nint(xss(l3))
               write(nsyso,'(4x,''nr ='',i4)') nrint
               if (nrint.ne.0) then
                  write(nsyso,'(4x,''nbt(i=1,nr) = '',20i5)')&
                    (nint(xss(l3+ii)),ii=1,nrint)
                  l3=l3+nrint
                  write(nsyso,'(4x,''int(i=1,nr) = '',20i5)')&
                    (nint(xss(l3+ii)),ii=1,nrint)
                  l3=l3+nrint
               endif
               l3=l3+1
               ne=nint(xss(l3))
               write(nsyso,'(4x,''ne ='',i4)') ne
               do ie=1,ne
                  e2=xss(l3+ie)
                  loci=nint(xss(l3+ne+ie))+dlwh-1
                  intt=nint(xss(loci))
                  nn=nint(xss(loci+1))
                  loci=loci+1
                  write(nsyso,'(/6x,&
                    &'' incident energy = '',1p,e14.6,''   intt ='',&
                    &i2,''    np = '',i4//&
                    &1x,8x,''energy'',11x,''pdf'',11x,''cdf'',&
                    &8x,''energy'',11x,''pdf'',11x,''cdf''/&
                    &1x,6(2x,''------------''))') e2,intt,nn
                  do j=1,nn,2
                     if (j.lt.nn) then
                        write(nsyso,'(1x,1p,6e14.6)')&
                          xss(j+loci),xss(j+nn+loci),&
                          xss(j+2*nn+loci),xss(j+1+loci),&
                          xss(j+1+nn+loci),xss(j+1+2*nn+loci)
                     else
                        write(nsyso,'(1x,1p,3e14.6)')&
                       xss(j+loci),xss(j+nn+loci),xss(j+2*nn+loci)
                     endif
                  enddo
               enddo

            !--law=44
            else if (law.eq.44) then
               nrint=nint(xss(l3))
               write(nsyso,'(4x,''nr ='',i4)') nrint
               if (nrint.ne.0) then
                  write(nsyso,'(4x,''nbt(i=1,nr) = '',20i5)')&
                    (nint(xss(l3+ii)),ii=1,nrint)
                  l3=l3+nrint
                  write(nsyso,'(4x,''int(i=1,nr) = '',20i5)')&
                    (nint(xss(l3+ii)),ii=1,nrint)
                  l3=l3+nrint
               endif
               l3=l3+1
               ne=nint(xss(l3))
               write(nsyso,'(4x,''ne ='',i4)') ne
               do ie=1,ne
                  e2=xss(l3+ie)
                  loci=nint(xss(l3+ne+ie))+dlwh-1
                  intt=mod(nint(xss(loci)),10)
                  nd=nint(xss(loci)/10)
                  nd=0
                  nn=nint(xss(loci+1))
                  loci=loci+1
                  write(nsyso,'(/6x,&
                    &'' incident energy = '',1p,e14.6,''   intt ='',&
                    &i2,''   nd = '',i4,''    np = '',i4//&
                    &1x,8x,''energy'',11x,''pdf'',11x,''cdf'',&
                    &13x,''r'',13x,''a''/&
                    &1x,5(2x,''------------'')/(1x,1p,5e14.6))')&
                    e2,intt,nd,nn,(xss(j+loci),xss(j+nn+loci),&
                    xss(j+2*nn+loci),xss(j+3*nn+loci),&
                    xss(j+4*nn+loci),j=1,nn)
               enddo

            !--law=33
            else if (law.eq.33) then
               write(nsyso,'(12x,''eout = c*(e-ec)   ec ='',&
                    &1p,e14.6,5x,''c ='',e14.6)') xss(l3),xss(l3+1)

            !--law=61
            else if (law.eq.61) then
               nrint=nint(xss(l3))
               write(nsyso,'(4x,''nr ='',i4)') nrint
               if (nrint.ne.0) then
                  write(nsyso,'(4x,''nbt(i=1,nr) = '',20i5)')&
                    (nint(xss(l3+ii)),ii=1,nrint)
                  l3=l3+nrint
                  write(nsyso,'(4x,''int(i=1,nr) = '',20i5)')&
                    (nint(xss(l3+ii)),ii=1,nrint)
                  l3=l3+nrint
               endif
               l3=l3+1
               ne=nint(xss(l3))
               write(nsyso,'(4x,''ne ='',i4)') ne
               do ie=1,ne
                  e2=xss(ie+l3)
                  loci=nint(xss(ie+ne+l3)+dlwh-1)
                  intt=mod(nint(xss(loci)),10)
                  nd=nint(xss(loci)/10)
                  nn=nint(xss(loci+1))
                  loci=loci+1
                  write(nsyso,'(/6x,'' incident energy = '',&
                    &1p,e14.6,''   intt ='',i2,''    nd = '',i4,&
                    &''    np = '',i4)') e2,intt,nd,nn
                  do ip=1,nn
                     locj=nint(xss(ip+3*nn+loci)+dlwh-1)
                     intmu=nint(xss(locj))
                     nmu=nint(xss(locj+1))
                     write(nsyso,'(/&
                       &6x,'' secondary energy = '',1p,e14.6/&
                       &6x,''              pdf = '',e14.6/&
                       &6x,''              cdf = '',e14.6/&
                       &6x,''            intmu = '',i8/&
                       &6x,''              nmu = '',i8/&
                  &''         cosine           pdf           cdf'',&
                  &''        cosine           pdf           cdf''/&
                  &''   ------------  ------------  ------------'',&
                  &''  ------------  ------------  ------------'')')&
                       xss(ip+loci),xss(ip+nn+loci),&
                       xss(ip+2*nn+loci),intmu,nmu
                     do imu=1,nmu,2
                        if (imu.eq.nmu) then
                           write(nsyso,'(1x,1p,3e14.6)')&
                             xss(locj+1+imu),xss(locj+1+nmu+imu),&
                             xss(locj+1+2*nmu+imu)
                        else
                           write(nsyso,'(1x,1p,6e14.6)')&
                             xss(locj+1+imu),xss(locj+1+nmu+imu),&
                             xss(locj+1+2*nmu+imu),xss(locj+1+imu+1),&
                             xss(locj+1+nmu+imu+1),&
                             xss(locj+1+2*nmu+imu+1)
                        endif
                     enddo
                  enddo
               enddo

            !--law=66
            else if (law.eq.66) then
               npsx=nint(xss(l3))
               apsx=xss(l3+1)
               intt=nint(xss(l3+2))
               nn=nint(xss(l3+3))
               loci=l3+3
               write(nsyso,'(12x,''npsx ='',i2/12x,''apsx ='',f10.4/&
                 &12x,''intt ='',i2)') npsx,apsx,intt
               write(nsyso,'(1x,2(''        energy'',&
                 &''           pdf           cdf'')/1x,&
                 &2(''  ------------'',&
                 &''  ------------  ------------'')/&
                 &(1x,1p,6e14.6))')&
                 (xss(j+loci),xss(j+nn+loci),xss(j+2*nn+loci),&
                 j=1,nn)

            !--law=67
            else if (law.eq.67) then
               nrint=nint(xss(l3))
               write(nsyso,'(4x,''nr ='',i4)') nrint
               if (nrint.ne.0) then
                  write(nsyso,'(4x,''nbt(i=1,nr) = '',20i5)')&
                    (nint(xss(l3+ii)),ii=1,nrint)
                  l3=l3+nrint
                  write(nsyso,'(4x,''int(i=1,nr) = '',20i5)')&
                    (nint(xss(l3+ii)),ii=1,nrint)
                  l3=l3+nrint
               endif
               l3=l3+1
               ne=nint(xss(l3))
               write(nsyso,'(4x,''ne ='',i4)') ne
               do ie=1,ne
                  e2=xss(l3+ie)
                  loci=nint(xss(l3+ne+ie))+dlwh-1
                  intmu=nint(xss(loci))
                  nmu=nint(xss(loci+1))
                  write(nsyso,'(/4x,''ein ='',f10.5,4x,&
                    &''intmu ='',i2,4x,''nmu ='',i3)') e2,intmu,nmu
                  locmu=loci+1
                  do k=1,nmu
                     amu=xss(locmu+k)
                     ll=nint(xss(locmu+nmu+k)+dlwh-1)
                     intep=nint(xss(ll))
                     ll=ll+1
                     npep=nint(xss(ll))
                     write(nsyso,'(/5x,''mu ='',f10.5,4x,&
                       &''intep ='',i2,4x,''npep ='',i3)')&
                       amu,intep,npep
                     write(nsyso,'(1x,2(''        energy'',&
                       &''           pdf           cdf'')/1x,&
                       &2(''  ------------'',&
                       &''  ------------  ------------'')/&
                       &(1x,1p,6e14.6))')&
                       (xss(ll+j),xss(ll+npep+j),xss(ll+2*npep+j),&
                       j=1,npep)
                  enddo
               enddo
            endif
         endif
      enddo

      !--print out yh block
      nyh=nint(xss(yh))
      if (nyh.ne.0) then
         write(nsyso,&
           '(//'' mts needed as multipliers for particle yields''/&
           &1x,5(''----------''),''-''//(6x,12i6))')&
           (nint(xss(ii+yh)),ii=1,nyh)
      endif

   !--continue loop over production types
   enddo
   return
   end subroutine acepcp

   subroutine aceout(itype,nace,ndir,hk,izn,awn,mcnpx)
   !-------------------------------------------------------------------
   ! Write ACE data out in desired format.
   !-------------------------------------------------------------------
   use util ! provides openz,closz,error
   ! externals
   integer::itype,nace,ndir,mcnpx
   character(70)::hk
   integer::izn(16)
   real(kr)::awn(16)
   ! internals
   integer::nout,ll,nn,n,i,lrec,nern,irec1

   !--branch according to ace file type
   nout=nace

   !--write type-1 header block
   if (itype.eq.1) then
      call openz(nout,1)
      if (mcnpx.eq.0) then
         write(nout,&
           '(a10,f12.6,1x,1p,e11.4,1x,a10/a70,a10)')&
           hz(1:10),aw0,tz,hd,hk,hm
      else
         write(nout,&
           '(a13,f12.6,1x,1p,e11.4,1x,a10/a70,a10)')&
           hz(1:13),aw0,tz,hd,hk,hm
      endif
      write(nout,'(4(i7,f11.0))') (izn(i),awn(i),i=1,16)
      write(nout,'(8i9)')&
        len2,izaid,nes,ntr,nr,ntrp,ntype,ndnf,(nxsd(i),i=1,8),&
        esz,nu,mtr,lqr,tyr,lsig,sig,land,and,ldlw,dlw,&
        gpd,mtrp,lsigp,sigp,landp,andp,ldlwp,dlwp,yp,fis,end,&
        iurpt,nud,dndat,ldnd,dnd,(jxsd(i),i=1,2),ptype,ntro,ploct

      !--write data in memory
      call change (nout)
      lrec=0
      nern=0
      irec1=1
      call closz(nout)

   !--write ace file in type 2 format
   else if (itype.eq.2) then
      call openz(-nout,1)
      if (mcnpx.eq.0) then
         write(nout) hz(1:10),aw0,tz,hd,hk,hm,&
           (izn(i),awn(i),i=1,16),&
           len2,izaid,nes,ntr,nr,ntrp,ntype,ndnf,(nxsd(i),i=1,8),&
           esz,nu,mtr,lqr,tyr,lsig,sig,land,and,ldlw,dlw,&
           gpd,mtrp,lsigp,sigp,landp,andp,ldlwp,dlwp,yp,fis,end,&
           iurpt,nud,dndat,ldnd,dnd,(jxsd(i),i=1,2),ptype,ntro,ploct
      else
         write(nout) hz(1:13),aw0,tz,hd,hk,hm,&
           (izn(i),awn(i),i=1,16),&
           len2,izaid,nes,ntr,nr,ntrp,ntype,ndnf,(nxsd(i),i=1,8),&
           esz,nu,mtr,lqr,tyr,lsig,sig,land,and,ldlw,dlw,&
           gpd,mtrp,lsigp,sigp,landp,andp,ldlwp,dlwp,yp,fis,end,&
           iurpt,nud,dndat,ldnd,dnd,(jxsd(i),i=1,2),ptype,ntro,ploct
      endif
      ll=0
      nn=len2
      do while (nn.gt.0)
         n=nn
         if (n.gt.ner) n=ner
         write(nout) (xss(ll+i),i=1,n)
         ll=ll+n
         nn=nn-n
      enddo
      lrec=ner*nbw
      nern=ner
      irec1=1
      call closz(-nout)
   endif

   !--output directory for mcnp
   call openz(ndir,1)
   if (mcnpx.eq.0) then
      if (iurpt.gt.0) then
         write(ndir,&
           '(a10,f12.6,'' filename route'',i2,i4,1x,i8,2i6,1p,e10.3,&
           &'' ptable'')') hz(1:10),aw0,itype,irec1,len2,lrec,nern,tz
      else
         write(ndir,&
           '(a10,f12.6,'' filename route'',i2,i4,1x,i8,2i6,1p,e10.3)')&
           hz(1:10),aw0,itype,irec1,len2,lrec,nern,tz
      endif
   else
      if (iurpt.gt.0) then
         write(ndir,&
           '(a13,f12.6,'' file route'',i2,i4,1x,i8,2i6,1p,e10.3,&
           &'' ptable'')') hz(1:13),aw0,itype,irec1,len2,lrec,nern,tz
      else
         write(ndir,&
           '(a13,f12.6,'' file route'',i2,i4,1x,i8,2i6,1p,e10.3)')&
          hz(1:13),aw0,itype,irec1,len2,lrec,nern,tz
      endif
   endif
   call closz(ndir)
   return
   end subroutine aceout

   subroutine change(nout)
   !-------------------------------------------------------------------
   ! Change ACE data fields from integer to real or vice versa.
   ! If nout.gt.1, the results are written in Type 1 format
   !    (all fields are assumed to contain real numbers).
   ! If nout.eq.0, real fields are changed to integers in memory
   !    (all fields are assumed to contain real numbers).
   ! If nout.eq.1, integer fields are changed to real in memory
   !    (fields are assumed to contain mixed reals and integers).
   !-------------------------------------------------------------------
   ! externals
   integer::nout
   ! internals
   integer::n,i,l,lnu,m,nc,j,nrr,ne,nn,ll,k,np,nw
   integer::ly,lnw,law,net,nmu,kk,nep,nure,nurb,mftype
   integer::nyp,ntro,jj,ir,nyh,li,ii,ntrh

   !--write or convert esz block
   n=5*nes
   do i=1,n
      call typen(i,nout,2)
   enddo

   !--nu block
   if (nu.ne.0) then
      l=nu
      lnu=nint(xss(l))
      if (lnu.gt.0) then
        m=1
      else
        m=2
        call typen(l,nout,1)
        l=l+1
      endif
      do i=1,m
         lnu=nint(xss(l))
         call typen(l,nout,1)
         l=l+1
         if (lnu.ne.2) then
            nc=nint(xss(l))
            call typen(l,nout,1)
            l=l+1
            do j=1,nc
               call typen(l,nout,2)
               l=l+1
            enddo
         else
            nrr=nint(xss(l))
            call typen(l,nout,1)
            l=l+1
            if (nrr.ne.0) then
               n=2*nrr
               do j=1,n
                  call typen(l,nout,1)
                  l=l+1
               enddo
            endif
            ne=nint(xss(l))
            call typen(l,nout,1)
            l=l+1
            n=2*ne
            do j=1,n
               call typen(l,nout,2)
               l=l+1
            enddo
         endif
      enddo
   endif

   !--cross section data
   if (ntr.ne.0) then

      !--mtr block
      l=mtr
      do i=1,ntr
         call typen(l,nout,1)
         l=l+1
      enddo

      !--lqr block
      l=lqr
      do i=1,ntr
         call typen(l,nout,2)
         l=l+1
      enddo

      !--tyr block
      l=tyr
      do i=1,ntr
         call typen(l,nout,1)
         l=l+1
      enddo

      !--lsig block
      l=lsig
      do i=1,ntr
         call typen(l,nout,1)
         l=l+1
      enddo

      !--sig block
      l=sig
      do i=1,ntr
         call typen(l,nout,1)
         l=l+1
         ne=nint(xss(l))
         call typen(l,nout,1)
         l=l+1
         do j=1,ne
            call typen(l,nout,2)
            l=l+1
         enddo
      enddo
   endif

   !--land block
   n=nr+1
   l=land
   li=l-1
   do i=1,n
      call typen(l,nout,1)
      l=l+1
   enddo

   !--and block
   l=and
   do i=1,n
      nn=nint(xss(li+i))
      if (nn.gt.0) then
         ne=nint(xss(l))
         call typen(l,nout,1)
         l=l+1
         do j=1,ne
            call typen(l,nout,2)
            l=l+1
         enddo
         ll=l-1
         do j=1,ne
            call typen(l,nout,1)
            l=l+1
         enddo
         do j=1,ne
            nn=nint(xss(ll+j))
            if (nn.ne.0) then
               if (nn.ge.0) then
                  do k=1,33
                     call typen(l,nout,2)
                     l=l+1
                  enddo
               else
                  np=nint(xss(iabs(nn)+and))
                  call typen(l,nout,1)
                  l=l+1
                  call typen(l,nout,1)
                  l=l+1
                  nw=3*np
                  do k=1,nw
                     call typen(l,nout,2)
                     l=l+1
                  enddo
               endif
            endif
         enddo
      endif
   enddo

   !--distributions
    if (nr.ne.0) then

      !--ldlw block
      l=ldlw
      do i=1,nr
         call typen(l,nout,1)
         l=l+1
      enddo

      !--dlw block
      l=dlw
      do i=1,nr
         ly=nint(xss(tyr+i-1))
         ly=iabs(ly)
         if (ly.gt.100) then
            l=ly-100+dlw-1
            nrr=nint(xss(l))
            call typen(l,nout,1)
            l=l+1
            if (nrr.gt.0) then
               n=2*nrr
               do j=1,n
                  call typen(l,nout,1)
                  l=l+1
               enddo
            endif
            ne=nint(xss(l))
            call typen(l,nout,1)
            l=l+1
            n=2*ne
            do j=1,n
               call typen(l,nout,2)
               l=l+1
            enddo
         endif

         !--loop over laws
         lnw=1
         do while (lnw.gt.0)
            lnw=nint(xss(l))
            call typen(l,nout,1)
            l=l+1
            law=nint(xss(l))
            call typen(l,nout,1)
            l=l+1
            call typen(l,nout,1)
            l=l+1
            nrr=nint(xss(l))
            call typen(l,nout,1)
            l=l+1
            if (nrr.gt.0) then
               n=2*nrr
               do j=1,n
                  call typen(l,nout,1)
                  l=l+1
               enddo
            endif
            ne=nint(xss(l))
            call typen(l,nout,1)
            l=l+1
            n=2*ne
            do j=1,n
               call typen(l,nout,2)
               l=l+1
            enddo

            !--law 1
            if (law.eq.1) then
               nrr=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               if (nrr.gt.0) then
                  n=2*nrr
                  do j=1,n
                     call typen(l,nout,1)
                     l=l+1
                  enddo
               endif
               ne=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               do j=1,ne
                  call typen(l,nout,2)
                  l=l+1
               enddo
               net=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               do j=1,ne
                  do k=1,net
                     call typen(l,nout,2)
                     l=l+1
                  enddo
               enddo

            !--law 3
            else if (law.eq.3) then
               call typen(l,nout,2)
               l=l+1
               call typen(l,nout,2)
               l=l+1

            !--law 4
            else if (law.eq.4) then
               nrr=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               if (nrr.gt.0) then
                  n=2*nrr
                  do j=1,n
                     call typen(l,nout,1)
                     l=l+1
                  enddo
               endif
               ne=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               do j=1,ne
                  call typen(l,nout,2)
                  l=l+1
               enddo
               do j=1,ne
                  call typen(l,nout,1)
                  l=l+1
               enddo
               do j=1,ne
                  call typen(l,nout,1)
                  l=l+1
                  np=nint(xss(l))
                  call typen(l,nout,1)
                  l=l+1
                  n=3*np
                  do k=1,n
                     call typen(l,nout,2)
                     l=l+1
                  enddo
               enddo

            !--law 5
            else if (law.eq.5) then
               nrr=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               if (nrr.gt.0) then
                  n=2*nrr
                  do j=1,n
                     call typen(l,nout,1)
                     l=l+1
                  enddo
               endif
               ne=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               n=2*ne
               do j=1,n
                  call typen(l,nout,2)
                  l=l+1
               enddo
               net=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               do j=1,net
                     call typen(l,nout,2)
                  l=l+1
               enddo

            !--law 7 or law 9
            else if (law.eq.7.or.law.eq.9) then
               nrr=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               if (nrr.gt.0) then
                  n=2*nrr
                  do j=1,n
                     call typen(l,nout,1)
                     l=l+1
                  enddo
               endif
               ne=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               n=2*ne
               do j=1,n
                  call typen(l,nout,2)
                  l=l+1
               enddo
               call typen(l,nout,2)
               l=l+1

            !--law 11
            else if (law.eq.11) then
               nrr=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               if (nrr.gt.0) then
                  n=2*nrr
                  do j=1,n
                     call typen(l,nout,1)
                     l=l+1
                  enddo
               endif
               ne=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               n=2*ne
               do j=1,n
                  call typen(l,nout,2)
                  l=l+1
               enddo
               nrr=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               if (nrr.gt.0) then
                  n=2*nrr
                  do j=1,n
                     call typen(l,nout,1)
                     l=l+1
                  enddo
               endif
               ne=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               n=2*ne
               do j=1,n
                  call typen(l,nout,2)
                  l=l+1
               enddo
               call typen(l,nout,2)
               l=l+1

            !--law 44
            else if (law.eq.44) then
               nrr=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               if (nrr.gt.0) then
                  n=2*nrr
                  do j=1,n
                     call typen(l,nout,1)
                     l=l+1
                  enddo
               endif
               ne=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               do j=1,ne
                  call typen(l,nout,2)
                  l=l+1
               enddo
               do j=1,ne
                  call typen(l,nout,1)
                  l=l+1
               enddo
               do j=1,ne
                  call typen(l,nout,1)
                  l=l+1
                  np=nint(xss(l))
                  call typen(l,nout,1)
                  l=l+1
                  n=5*np
                  do k=1,n
                     call typen(l,nout,2)
                     l=l+1
                  enddo
               enddo

            !--law 61
            else if (law.eq.61) then
               nrr=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               if (nrr.ne.0) then
                  n=2*nrr
                  do j=1,n
                     call typen(l,nout,1)
                     l=l+1
                  enddo
               endif
               ne=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
                  do j=1,ne
                  call typen(l,nout,2)
                  l=l+1
               enddo
               do j=1,ne
                  call typen(l,nout,1)
                  l=l+1
               enddo
               do j=1,ne
                  call typen(l,nout,1)
                  l=l+1
                  np=nint(xss(l))
                  call typen(l,nout,1)
                  l=l+1
                  n=3*np
                  do k=1,n
                     call typen(l,nout,2)
                     l=l+1
                  enddo
                  do k=1,np
                     call typen(l,nout,1)
                     l=l+1
                  enddo
                  do k=1,np
                     call typen(l,nout,1)
                     l=l+1
                     nmu=nint(xss(l))
                     call typen(l,nout,1)
                     l=l+1
                     nw=3*nmu
                     do kk=1,nw
                        call typen(l,nout,2)
                        l=l+1
                     enddo
                  enddo
               enddo

            !--law 66
            else if (law.eq.66) then
               call typen(l,nout,1)
               l=l+1
               call typen(l,nout,2)
               l=l+1
               call typen(l,nout,1)
               l=l+1
               nn=nint(xss(l))
               n=3*nn
               call typen(l,nout,1)
               l=l+1
               do k=1,n
                  call typen(l,nout,2)
                  l=l+1
               enddo

            !--law 67
            else if (law.eq.67) then
               nrr=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               if (nrr.gt.0) then
                  n=2*nrr
                  do j=1,n
                     call typen(l,nout,1)
                     l=l+1
                  enddo
               endif
               ne=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               do j=1,ne
                  call typen(l,nout,2)
                  l=l+1
               enddo
               do j=1,ne
                  call typen(l,nout,1)
                  l=l+1
               enddo
               do j=1,ne
                  call typen(l,nout,1)
                  l=l+1
                  nmu=nint(xss(l))
                  call typen(l,nout,1)
                  l=l+1
                  do k=1,nmu
                     call typen(l,nout,2)
                     l=l+1
                     enddo
                  do k=1,nmu
                     call typen(l,nout,1)
                     l=l+1
                  enddo
                  do k=1,nmu
                     call typen(l,nout,1)
                     l=l+1
                     nep=nint(xss(l))
                     call typen(l,nout,1)
                     l=l+1
                     nn=3*nep
                     do n=1,nn
                        call typen(l,nout,2)
                        l=l+1
                     enddo
                  enddo
               enddo
            endif
         enddo
      enddo
   endif

   !--unresolved-range probability-table block
   if (iurpt.gt.0) then
      l=iurpt
      nure=nint(xss(l))
      call typen(l,nout,1)
      l=l+1
      nurb=nint(xss(l))
      call typen(l,nout,1)
      l=l+1
      call typen(l,nout,1)
      l=l+1
      call typen(l,nout,1)
      l=l+1
      call typen(l,nout,1)
      l=l+1
      call typen(l,nout,1)
      l=l+1
      n=nure*(1+6*nurb)
      do i=1,n
         call typen(l,nout,2)
         l=l+1
      enddo
   endif

   !--delayed neutron block
   if (ndnf.ne.0) then
      !--delayed nubar
      l=nud
      lnu=nint(xss(l))
      call typen(l,nout,1)
      l=l+1
      nrr=nint(xss(l))
      call typen(l,nout,1)
      l=l+1
      if (nrr.ne.0) then
         n=2*nrr
         do j=1,n
            call typen(l,nout,1)
            l=l+1
         enddo
      endif
      ne=nint(xss(l))
      call typen(l,nout,1)
      l=l+1
      n=2*ne
      do j=1,n
         call typen(l,nout,2)
         l=l+1
      enddo
      !--precursor data
      l=dndat
      do i=1,ndnf
         call typen(l,nout,2)
         l=l+1
      nrr=nint(xss(l))
         call typen(l,nout,1)
         l=l+1
      if (nrr.ne.0) then
         n=2*nrr
            do j=1,n
               call typen(l,nout,1)
               l=l+1
            enddo
         endif
         ne=nint(xss(l))
         call typen(l,nout,1)
         l=l+1
         n=2*ne
         do j=1,n
            call typen(l,nout,2)
            l=l+1
         enddo
      enddo
      !--precursor energy distribution locators
      do i=1,ndnf
         call typen(l,nout,1)
         l=l+1
      enddo
      !--precursor energy distributions
      do i=1,ndnf
         call typen(l,nout,1)
         l=l+1
         call typen(l,nout,1)
         l=l+1
         call typen(l,nout,1)
         l=l+1
      nrr=nint(xss(l))
         call typen(l,nout,1)
         l=l+1
      if (nrr.ne.0) then
         n=2*nrr
            do j=1,n
               call typen(l,nout,1)
               l=l+1
            enddo
         endif
         ne=nint(xss(l))
         call typen(l,nout,1)
         l=l+1
         n=2*ne
         do j=1,n
            call typen(l,nout,2)
            l=l+1
         enddo
         !--law=4 data
      nrr=nint(xss(l))
         call typen(l,nout,1)
         l=l+1
      if (nrr.ne.0) then
         n=2*nrr
            do j=1,n
               call typen(l,nout,1)
               l=l+1
            enddo
         endif
         ne=nint(xss(l))
         call typen(l,nout,1)
         l=l+1
         do j=1,ne
            call typen(l,nout,2)
            l=l+1
         enddo
         do j=1,ne
            call typen(l,nout,1)
            l=l+1
         enddo
         do j=1,ne
            call typen(l,nout,1)
            l=l+1
            np=nint(xss(l))
            call typen(l,nout,1)
            l=l+1
            n=3*np
            do k=1,n
               call typen(l,nout,2)
               l=l+1
            enddo
         enddo
      enddo
   endif

   !--gpd block
   if (gpd.ne.0) then
      l=gpd
      do i=1,nes
         call typen(l,nout,2)
         l=l+1
      enddo
      if (mtrp.gt.gpd+nes) then
         n=20*30
         do i=1,n
            call typen(l,nout,2)
            l=l+1
         enddo
      endif
   endif

   !--detailed photon production
   if (ntrp.ne.0) then

      !--mtrp block
      l=mtrp
      do i=1,ntrp
         call typen(l,nout,1)
         l=l+1
      enddo

      !--lsigp block
      l=lsigp
      do i=1,ntrp
         call typen(l,nout,1)
         l=l+1
      enddo

      !--sigp block
      l=sigp
      do i=1,ntrp
         mftype=nint(xss(l))
         call typen(l,nout,1)
         l=l+1
         if (mftype.ne.12.and.mftype.ne.16) then
            call typen(l,nout,1)
            l=l+1
            ne=nint(xss(l))
            call typen(l,nout,1)
            l=l+1
            do j=1,ne
               call typen(l,nout,2)
               l=l+1
            enddo
         else
            call typen(l,nout,1)
            l=l+1
            nrr=nint(xss(l))
            call typen(l,nout,1)
            l=l+1
            if (nrr.gt.0) then
               n=2*nrr
               do j=1,n
                  call typen(l,nout,1)
                  l=l+1
               enddo
            endif
            ne=nint(xss(l))
            call typen(l,nout,1)
            l=l+1
            n=2*ne
            do j=1,n
               call typen(l,nout,2)
               l=l+1
            enddo
         endif
      enddo

      !--landp block
      l=landp
      li=l-1
      do i=1,ntrp
         call typen(l,nout,1)
         l=l+1
      enddo

      !--andp block
      l=andp
      do i=1,ntrp
         nn=nint(xss(li+i))
         if (nn.gt.0) then
            ne=nint(xss(l))
            call typen(l,nout,1)
            l=l+1
            do j=1,ne
               call typen(l,nout,2)
               l=l+1
            enddo
            ll=l-1
            do j=1,ne
               call typen(l,nout,1)
               l=l+1
            enddo
            do j=1,ne
               nn=nint(xss(ll+j))
               if (nn.gt.0) then
                  do k=1,33
                     call typen(l,nout,2)
                     l=l+1
                  enddo
               endif
            enddo
         endif
      enddo

      !--ldlwp block
      l=ldlwp
      do i=1,ntrp
         call typen(l,nout,1)
         l=l+1
      enddo

      !--dlwp block
      l=dlwp
      do i=1,ntrp
         lnw=1
         do while (lnw.ne.0)
            lnw=nint(xss(l))
            call typen(l,nout,1)
            l=l+1
            law=nint(xss(l))
            call typen(l,nout,1)
            l=l+1
            call typen(l,nout,1)
            l=l+1
            nrr=nint(xss(l))
            call typen(l,nout,1)
            l=l+1
            if (nrr.gt.0) then
               n=2*nrr
               do j=1,n
                  call typen(l,nout,1)
                  l=l+1
               enddo
            endif
            ne=nint(xss(l))
            call typen(l,nout,1)
            l=l+1
            n=2*ne
            do j=1,n
               call typen(l,nout,2)
               l=l+1
            enddo

            !--law 1
            if (law.eq.1) then
               nrr=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               if (nrr.gt.0) then
                  n=2*nrr
                  do j=1,n
                     call typen(l,nout,1)
                     l=l+1
                  enddo
               endif
               ne=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               do j=1,ne
                  call typen(l,nout,2)
                  l=l+1
               enddo
               net=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               do j=1,ne
                  do k=1,net
                     call typen(l,nout,2)
                     l=l+1
                  enddo
               enddo

            !--law 2
            else if (law.eq.2) then
               call typen(l,nout,1)
               l=l+1
               call typen(l,nout,2)
               l=l+1

            !--law 4 and law 44
            else if (law.eq.4.or.law.eq.44) then
               nrr=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               if (nrr.gt.0) then
                  n=2*nrr
                  do j=1,n
                     call typen(l,nout,1)
                     l=l+1
                  enddo
               endif
               ne=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               do j=1,ne
                  call typen(l,nout,2)
                  l=l+1
               enddo
               do j=1,ne
                  call typen(l,nout,1)
                  l=l+1
               enddo
               do j=1,ne
                  call typen(l,nout,1)
                  l=l+1
                  np=nint(xss(l))
                  call typen(l,nout,1)
                  l=l+1
                  n=3*np
                  do k=1,n
                     call typen(l,nout,2)
                     l=l+1
                  enddo
               enddo
            endif
         enddo
      enddo

      !--yp block
        l=yp
        nyp=nint(xss(l))
        call typen(l,nout,1)
        l=l+1
        do i=1,nyp
           call typen(l,nout,1)
           l=l+1
        enddo
     endif

   !--particle production blocks
     if (ntype.gt.0) then

      !--ptype, ntro, and ixs arrays
        do i=1,ntype
           call typen(l,nout,1)
           l=l+1
        enddo
        ntro=l
        do i=1,ntype
           call typen(l,nout,1)
           l=l+1
        enddo
        nw=10*ntype
        do i=1,nw
           call typen(l,nout,1)
           l=l+1
        enddo

      !--loop over particle types
      do i=1,ntype

         !--hpd block
         call typen(l,nout,1)
         l=l+1
         ne=nint(xss(l))
         call typen(l,nout,1)
         l=l+1
         if (ne.ne.0) then
            nw=2*ne
            do j=1,nw
               call typen(l,nout,2)
               l=l+1
            enddo
            ntrh=nint(xss(ntro-1+i))

            !--mtrh block
            do k=1,ntrh
               call typen(l,nout,1)
               l=l+1
            enddo

            !--tyrh block
            do k=1,ntrh
               call typen(l,nout,1)
               l=l+1
            enddo

            !--lsigh block
            do k=1,ntrh
               call typen(l,nout,1)
               l=l+1
            enddo

            !--sigh block
            do j=1,ntrh
               call typen(l,nout,1)
               l=l+1
               call typen(l,nout,1)
               l=l+1
               nrr=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               if (nrr.gt.0) then
                  n=2*nrr
                  do jj=1,n
                     call typen(l,nout,1)
                     l=l+1
                  enddo
               endif
               ne=nint(xss(l))
               call typen(l,nout,1)
               l=l+1
               nw=2*ne
               do k=1,nw
                  call typen(l,nout,2)
                  l=l+1
               enddo
            enddo

            !--landh block
            li=l-1
            do k=1,ntrh
               call typen(l,nout,1)
               l=l+1
            enddo

            !--andh block
            do ir=1,ntrh
               nn=nint(xss(li+ir))
               if (nn.gt.0) then
                  ne=nint(xss(l))
                  call typen(l,nout,1)
                  l=l+1
                  do j=1,ne
                     call typen(l,nout,2)
                     l=l+1
                  enddo
                  ll=l-1
                  do j=1,ne
                     call typen(l,nout,1)
                     l=l+1
                  enddo
                  do j=1,ne
                     nn=nint(xss(ll+j))
                     if (nn.gt.0) then
                        do k=1,33
                           call typen(l,nout,2)
                           l=l+1
                        enddo
                     else if (nn.lt.0) then
                        call typen(l,nout,1)
                        l=l+1
                        np=nint(xss(l))
                        call typen(l,nout,1)
                        l=l+1
                        nw=3*np
                        do k=1,nw
                           call typen(l,nout,2)
                           l=l+1
                        enddo
                     endif
                  enddo
               endif
            enddo

            !--ldlwh block
            li=l-1
            do k=1,ntrh
               call typen(l,nout,1)
               l=l+1
            enddo

            !--dlwh block
            do ii=1,ntrh
               nn=nint(xss(li+ii))
               if (nn.gt.0) then
                  call typen(l,nout,1)
                  l=l+1
                  law=nint(xss(l))
                  call typen(l,nout,1)
                  l=l+1
                  call typen(l,nout,1)
                  l=l+1
                  nrr=nint(xss(l))
                  call typen(l,nout,1)
                  l=l+1
                  if (nrr.ne.0) then
                     nw=2*nrr
                     do k=1,nw
                        call typen(l,nout,1)
                        l=l+1
                     enddo
                  endif
                  ne=nint(xss(l))
                  call typen(l,nout,1)
                  l=l+1
                  nw=2*ne
                  do k=1,nw
                     call typen(l,nout,2)
                     l=l+1
                  enddo
                  if (law.eq.4) then
                     nrr=nint(xss(l))
                     call typen(l,nout,1)
                     l=l+1
                     ne=nint(xss(l))
                     call typen(l,nout,1)
                     l=l+1
                     do k=1,ne
                        call typen(l,nout,2)
                        l=l+1
                     enddo
                     do k=1,ne
                        call typen(l,nout,1)
                        l=l+1
                     enddo
                     do k=1,ne
                        call typen(l,nout,1)
                        l=l+1
                        np=nint(xss(l))
                        call typen(l,nout,1)
                        l=l+1
                        nw=3*np
                        do kk=1,nw
                           call typen(l,nout,2)
                           l=l+1
                        enddo
                     enddo
                  else if (law.eq.44) then
                     nrr=nint(xss(l))
                     call typen(l,nout,1)
                     l=l+1
                     ne=nint(xss(l))
                     call typen(l,nout,1)
                     l=l+1
                     do k=1,ne
                        call typen(l,nout,2)
                        l=l+1
                     enddo
                     do k=1,ne
                        call typen(l,nout,1)
                        l=l+1
                     enddo
                     do j=1,ne
                        call typen(l,nout,1)
                        l=l+1
                        np=nint(xss(l))
                        call typen(l,nout,1)
                        l=l+1
                        nw=5*np
                        do k=1,nw
                           call typen(l,nout,2)
                           l=l+1
                        enddo
                     enddo
                  else if (law.eq.33) then
                     call typen(l,nout,2)
                     l=l+1
                     call typen(l,nout,2)
                     l=l+1
                  else if (law.eq.66) then
                     call typen(l,nout,1)
                     l=l+1
                     call typen(l,nout,2)
                     l=l+1
                     call typen(l,nout,1)
                     l=l+1
                     nn=nint(xss(l))
                     n=3*nn
                     call typen(l,nout,1)
                     l=l+1
                     do k=1,n
                        call typen(l,nout,2)
                        l=l+1
                     enddo
                  else if (law.eq.61) then
                     nrr=nint(xss(l))
                     call typen(l,nout,1)
                     l=l+1
                     if (nrr.ne.0) then
                        n=2*nrr
                        do j=1,n
                           call typen(l,nout,1)
                           l=l+1
                        enddo
                     endif
                     ne=nint(xss(l))
                     call typen(l,nout,1)
                     l=l+1
                     do j=1,ne
                        call typen(l,nout,2)
                        l=l+1
                     enddo
                     do j=1,ne
                        call typen(l,nout,1)
                        l=l+1
                     enddo
                     do j=1,ne
                        call typen(l,nout,1)
                        l=l+1
                        np=nint(xss(l))
                        call typen(l,nout,1)
                        l=l+1
                        n=3*np
                        do k=1,n
                           call typen(l,nout,2)
                           l=l+1
                        enddo
                        do k=1,np
                           call typen(l,nout,1)
                           l=l+1
                        enddo
                        do k=1,np
                           call typen(l,nout,1)
                           l=l+1
                           nmu=nint(xss(l))
                           call typen(l,nout,1)
                           l=l+1
                           nw=3*nmu
                           do kk=1,nw
                              call typen(l,nout,2)
                              l=l+1
                           enddo
                        enddo
                     enddo
                  else if (law.eq.67) then
                     do j=1,ne
                        call typen(l,nout,1)
                        l=l+1
                        nmu=nint(xss(l))
                        call typen(l,nout,1)
                        l=l+1
                        do k=1,nmu
                           call typen(l,nout,2)
                           l=l+1
                        enddo
                        do k=1,nmu
                           call typen(l,nout,1)
                           l=l+1
                        enddo
                        do k=1,nmu
                           call typen(l,nout,1)
                           l=l+1
                           nep=nint(xss(l))
                           call typen(l,nout,1)
                           l=l+1
                           nn=3*nep
                           do n=1,nn
                              call typen(l,nout,2)
                              l=l+1
                           enddo
                        enddo
                     enddo
                  endif
               endif
            enddo

            !--yh block
            nyh=nint(xss(l))
            call typen(l,nout,1)
            l=l+1
            do k=1,nyh
               call typen(l,nout,1)
               l=l+1
            enddo
         endif
      enddo
   endif
   call typen(0,nout,3)

   return
   end subroutine change

   subroutine typen(l,nout,iflag)
   !-------------------------------------------------------------------
   ! Write an integer or a real number to a Type-1 ACE file,
   ! or (if nout=0) convert real to integer for type-3 output,
   ! or (if nout=1) convert integer to real for type-3 input.
   ! Use iflag.eq.1 to write an integer (i20).
   ! Use iflag.eq.2 to write a real number (1pe20.11).
   ! Use iflag.eq.3 to write partial line at end of file.
   !-------------------------------------------------------------------
   ! externals
   integer::l,nout,iflag
   ! internals
   integer::i,j
   character(20)::hl(4)
   save hl,i

   if (iflag.eq.3.and.nout.gt.1.and.i.lt.4) then
      write(nout,'(4a20)') (hl(j),j=1,i)
   else
      i=mod(l-1,4)+1
      if (iflag.eq.1) write(hl(i),'(i20)') nint(xss(l))
      if (iflag.eq.2) write(hl(i),'(1p,e20.11)') xss(l)
      if (i.eq.4) write(nout,'(4a20)') (hl(j),j=1,i)
   endif
   return
   end subroutine typen

   subroutine acefix(nin,itype,nout,ndir,iprint,nplot,suff,&
     nxtra,hk,izn,awn,mcnpx)
   !-------------------------------------------------------------------
   ! Print or edit ACE files.
   !-------------------------------------------------------------------
   use util ! provides openz,closz,error
   ! externals
   integer::nin,itype,nout,ndir,iprint,nplot,nxtra,mcnpx
   real(kr)::suff
   character(70)::hk
   integer::izn(16)
   real(kr)::awn(16)
   ! internals
   integer::n,l,j,iza,max,i
   integer::izo(16)
   real(kr)::zaid
   real(kr)::awo(16)
   character(70)::hko
   character(10)::str
   character(3)::ht
   real(kr),dimension(6),parameter::awit=(/1.0e0_kr,0.99862e0_kr,&
     1.99626e0_kr,2.98960e0_kr,2.98903e0_kr,3.96713e0_kr/)
   real(kr),parameter::zero=0

   !--branch on input file type

   !--read type 1 ace format file
   call openz(nin,0)
   if (itype.eq.1) then
      if (mcnpx.eq.0) then
         read(nin,&
           '(a10,f12.6,1x,1p,e11.4,1x,a10/a70,a10)')&
           hz(1:10),aw0,tz,hd,hko,hm
      else
         read(nin,&
           '(a13,f12.6,1x,1p,e11.4,1x,a10/a70,a10)')&
           hz(1:13),aw0,tz,hd,hko,hm
      endif
      read (nin,'(4(i7,f11.0))') (izo(i),awo(i),i=1,16)
      read (nin,'(8i9)')&
        len2,izaid,nes,ntr,nr,ntrp,ntype,ndnf,(nxsd(i),i=1,8),&
        esz,nu,mtr,lqr,tyr,lsig,sig,land,and,ldlw,dlw,&
        gpd,mtrp,lsigp,sigp,landp,andp,ldlwp,dlwp,yp,fis,end,&
        iurpt,nud,dndat,ldnd,dnd,(jxsd(i),i=1,2),ptype,ntro,ploct
      n=(len2+3)/4
      l=0
      do i=1,n
         read (nin,'(4e20.0)') (xss(l+j),j=1,4)
         l=l+4
      enddo
      call closz(nin)
      if (mcnpx.gt.0) then
         if (hz(11:11).eq.'p') then
            izai=0
            awi=0
         else if (hz(11:11).eq.'n') then
            izai=1
            awi=awit(1)
         else if (hz(11:11).eq.'h') then
            izai=1001
            awi=awit(2)
         else if (hz(11:11).eq.'d') then
            izai=1002
            awi=awit(3)
         else if (hz(11:11).eq.'t') then
            izai=1003
            awi=awit(4)
         else if (hz(11:11).eq.'s') then
            izai=2003
            awi=awit(5)
         else if (hz(11:11).eq.'a') then
            izai=2004
            awi=awit(6)
         else
           call error('acefix',&
             'problem with particle id in zaid',' ')
         endif
      else
         if (hz(10:10).eq.'p') then
            izai=0
            awi=0
         else if (hz(10:10).eq.'u') then
            izai=0
            awi=0
         else if (hz(10:10).eq.'c') then
            izai=1
            awi=awit(1)
         else if (hz(10:10).eq.'t') then
            izai=1
            awi=awit(1)
         else if (hz(10:10).eq.'y') then
            izai=1
            awi=awit(1)
         else if (hz(10:10).eq.'h') then
            izai=1001
            awi=awit(2)
         else if (hz(10:10).eq.'o') then
            izai=1002
            awi=awit(3)
         else if (hz(10:10).eq.'r') then
            izai=1003
            awi=awit(4)
         else if (hz(10:10).eq.'s') then
            izai=2003
            awi=awit(5)
         else if (hz(10:10).eq.'a') then
            izai=2004
            awi=awit(6)
         else
           call error('acefix',&
             'problem with particle id in zaid',' ')
         endif
      endif

   !--read type 2 ace format file
   else if (itype.eq.2) then
      if (mcnpx.eq.0) then
        read(nin) hz(1:10),aw0,tz,hd,hko,hm,(izo(i),awo(i),i=1,16),&
          len2,izaid,nes,ntr,nr,ntrp,ntype,ndnf,(nxsd(i),i=1,8),&
          esz,nu,mtr,lqr,tyr,lsig,sig,land,and,ldlw,dlw,&
          gpd,mtrp,lsigp,sigp,landp,andp,ldlwp,dlwp,yp,fis,end,&
          iurpt,nud,dndat,ldnd,dnd,(jxsd(i),i=1,2),ptype,ntro,ploct
      else
        read(nin) hz(1:13),aw0,tz,hd,hko,hm,(izo(i),awo(i),i=1,16),&
          len2,izaid,nes,ntr,nr,ntrp,ntype,ndnf,(nxsd(i),i=1,8),&
          esz,nu,mtr,lqr,tyr,lsig,sig,land,and,ldlw,dlw,&
          gpd,mtrp,lsigp,sigp,landp,andp,ldlwp,dlwp,yp,fis,end,&
          iurpt,nud,dndat,ldnd,dnd,(jxsd(i),i=1,2),ptype,ntro,ploct
      endif
      n=(len2+ner-1)/ner
      l=0
      do i=1,n
         max=len2-l
         if (max.gt.ner) max=ner
         read (nin) (xss(l+j),j=1,max)
         l=l+max
      enddo
      call closz(-nin)
      if (mcnpx.gt.0) then
         if (hz(8:8).eq.'p') then
            izai=0
            awi=0
         else if (hz(8:8).eq.'n') then
            izai=1
            awi=awit(1)
         else if (hz(8:8).eq.'h') then
            izai=1001
            awi=awit(2)
         else if (hz(8:8).eq.'d') then
            izai=1002
            awi=awit(3)
         else if (hz(8:8).eq.'t') then
            izai=1003
            awi=awit(4)
         else if (hz(8:8).eq.'s') then
            izai=2003
            awi=awit(5)
         else if (hz(8:8).eq.'a') then
            izai=2004
            awi=awit(6)
         else
            call error('acefix',&
             'problem with particle id in zaid',' ')
         endif
      else
         if (hz(10:10).eq.'p') then
            izai=0
            awi=0
         else if (hz(10:10).eq.'u') then
            izai=0
            awi=0
         else if (hz(10:10).eq.'c') then
            izai=1
            awi=awit(1)
         else if (hz(10:10).eq.'t') then
            izai=1
            awi=awit(1)
         else if (hz(10:10).eq.'y') then
            izai=1
            awi=awit(1)
         else if (hz(10:10).eq.'h') then
            izai=1001
            awi=awit(2)
         else if (hz(10:10).eq.'o') then
            izai=1002
            awi=awit(3)
         else if (hz(10:10).eq.'r') then
            izai=1003
            awi=awit(4)
         else if (hz(10:10).eq.'s') then
            izai=2003
            awi=awit(5)
         else if (hz(10:10).eq.'a') then
            izai=2004
            awi=awit(6)
         else
           call error('acefix',&
             'problem with particle id in zaid',' ')
         endif
      endif
   endif

   !--adjust zaid
   if (mcnpx.gt.0) then
       read(hz,'(f10.0,a3)') zaid,ht
   else
       read(hz,'(a9,a1)') str,ht
       if (ht(1:1).ne.'t') then
          read(hz,'(f9.0,a1)') zaid,ht
       endif
   endif
   if (suff.ge.zero.and.ht(1:1).ne.'t') then
      iza=int(zaid+0.001)
      zaid=iza+suff
      if (mcnpx.gt.0) then
          write(hz,'(f10.3,a3)') zaid,ht
      else
          write(hz,'(f9.2,a1)') zaid,ht
      endif
   endif

   !--adjust comment and (iz,aw) list
   if (len_trim(hk).eq.0) then
      hk=hko
   endif
   if (nxtra.ne.0) then
      do i=1,nxtra
         izn(i)=izo(i)
         awn(i)=awo(i)
      enddo
   endif

   !--print, plot, and/or write the file.
   !--check consistency for c files.
   if ((mcnpx.eq.0.and.&
     (ht.eq.'c'.or.ht.eq.'h'.or.ht.eq.'o'.or.ht.eq.'r'&
     .or.ht.eq.'s'.or.ht.eq.'a'))&
     .or.(mcnpx.gt.0.and.ht(2:2).eq.'c')) then
      call consis
      if (iprint.gt.0) call aceprt(hk)
      if (nout.gt.0) call aceout(itype,nout,ndir,hk,izn,awn,mcnpx)
      if (nplot.ne.0) then
         call aplots(nplot,hk)
         call closz(nplot)
      endif
   else
      call error('acefix','illegal file type.',' ')
   endif
   return
   end subroutine acefix

   subroutine consis
   !-------------------------------------------------------------------
   ! Do basic consistency checks on the ACE file in memory.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util ! provides mess
   use acecm ! provides mtname
   ! internals
   integer::nerr,i,k,iaa,iin,nrl,nn,n,na,ic,id,ne,nr1
   integer::nb,nk,ie,k1,im,ll,nlaw,icm,j,m,law,loci,intt
   integer::n2big,ishift,locj,nmu,loct,loc1,mftype,mtmult
   integer::ii,naa,locv,locc,iflag,imt,l1,l2,l3,mt,l
   integer::ipt,mtrh,nmtr
   integer::hpd,lsigh,sigh,ldlwh,dlwh,nj
   integer::locl(20)
   integer::nure,nurb,lurt,luri,lura,lurf,ib
   integer::nerrtb,nerrxt,nerrxe,nerrxf,nerrxg,nerrxh
   real(kr)::thresh,elast,e,aprime,q,epmax,clast,ep,c,p,r
   real(kr)::cclast,co,cc,gsum,ss,y,gg,sum,frac,x,xlast
   character(10)::name
   real(kr),parameter::elow=1.e-11_kr
   real(kr),parameter::oneup=1.0001e0_kr
   real(kr),parameter::oplus=1.000001e0_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1

   write(nsyso,'(/'' ace consistency checks''/&
                 &'' ----------------------'')')
   nerr=0

   !--check threshold against q values
   write(nsyso,'(/'' check reaction thresholds against q values'')')
   do i=1,ntr
      k=nint(xss(lsig-1+i)+sig-1)
      iaa=nint(xss(k))
      iin=nint(xss(mtr-1+i))
      call mtname(iin,name,izai)
      thresh=(aw0+awi)*(-xss(lqr-1+i))/aw0
      if (thresh.lt.elow) thresh=elow
      if (xss(esz-1+iaa).lt.thresh.and.thresh.gt.elow) then
         write(nsyso,'(''   consis: threshold'',1p,e16.8,&
           &'' less than the expected'',e16.8,'' for '',a)')&
           xss(esz-1+iaa),thresh,name
         nerr=nerr+1
      endif
   enddo

   !--check for monotonic energy grid
   write(nsyso,'(/'' check that main energy grid is monotonic'')')
   elast=0
   do i=1,nes
      e=xss(esz-1+i)
      if (e.le.elast) then
         write(nsyso,'(''   consis: energy'',1p,e16.8,&
           &'' less than'',e16.8,'' (see point no.'',i7,'')'')')&
           e,elast,i
         nerr=nerr+1
      endif
      elast=e
   enddo

   !--check angular distributions for correct reference frame
   write(nsyso,'(/'' check angular distributions '',&
     &''for correct reference frame'')')
   nrl=nr+1
   do nn=1,nrl
      n=nn-1
      na=nint(xss(land+n))
      if (na.ge.0) then
         if (na.lt.0) na=-na
         if (n.eq.0) ic=2
         if (n.gt.0) ic=nint(xss(mtr-1+n))
         if (n.eq.0) id=-1
         if (n.gt.0) id=nint(xss(tyr-1+n))
         na=and-1+na
         ne=nint(xss(na))
         call mtname(ic,name,izai)
         if (((ic.ge.51.and.ic.le.90).or.&
           (ic.ge.6.and.ic.le.9)).and.id.gt.0) then
            write(nsyso,'(''   consis: should be cm: '',a10)') name
            nerr=nerr+1
         endif
         if (((ic.ge.16.and.ic.lt.50).or.ic.eq.91).and.id.lt.0) then
            write(nsyso,'(''   consis: should be lab: '',a10)') name
            nerr=nerr+1
         endif
      endif
   enddo

   !--check angular distributions for reasonable cosines
   write(nsyso,'(/'' check angular distributions '',&
     &''for unreasonable cosine values'')')
   nr1=nr+1
   do nn=1,nr1
      n=nn-1
      na=nint(xss(land+n))
      if (n.eq.0) then
         mt=2
         name='elastic'
      else
         mt=iabs(nint(xss(mtr+n-1)))
         call mtname(mt,name,izai)
      endif
      l=len_trim(name)
      if (na.gt.0) then
         na=na+and-1
         ne=nint(xss(na))
         nb=na+ne
         k=nint(xss(nb+1))

         !--equiprobable bins
         if (k.ge.0) then
            nk=0
            do ie=1,ne
               k=nint(xss(ie+nb))
               if (k.ne.0) then
                  k1=k+and-1
                  if (ie.lt.ne) nk=nint(xss(ie+1+nb))-k
                  do im=1,nk
                     if (xss(k1+im-1).lt.-oneup&
                       .or.xss(k1+im-1).gt.oneup) then
                        write(nsyso,'(''   consis:'',&
                          &'' mu out of range for '',a,&
                          &'' at e='',1p,2e14.6)') name(1:l),&
                          xss(na+ie),xss(k1+im-1)
                        nerr=nerr+1
                        endif
                     if (im.gt.1.and.&
                       xss(k1+im-1).le.xss(k1+im-2)) then
                        write(nsyso,&
                          '(''   consis:'',&
                          &'' mu values not monotonic for '',a,&
                          &''at e='',1p,3e14.6)') name(1:l),&
                          xss(na+ie),xss(k1+im-2),xss(k1+im-1)
                        nerr=nerr+1
                     endif
                  enddo
               endif
            enddo

         !--cummulative distributions
         else
            nk=0
            do ie=1,ne
               k=nint(xss(ie+nb))
               if (k.ne.0) then
                  k1=k+and-1
                  if (ie.lt.ne) nk=nint(xss(ie+1+nb))-k
                  do im=1,nk
                     ll=k1+im-1
                     if (xss(ll).lt.-oneup.or.xss(ll).gt.oneup) then
                        write(nsyso,'(''   consis:'',&
                           &'' mu out of range for '',&
                           &a,'' at e='',1p,2e14.6)') name(1:l),&
                           xss(na+ie),xss(ll)
                        nerr=nerr+1
                     endif
                     if (im.gt.1.and.xss(ll).le.xss(ll-1)) then
                        write(nsyso,&
                          '(''   consis:'',&
                          &'' mu values not monotonic for '',a,&
                          &'' at e='',1p,3e14.6)') name(1:l),&
                          xss(na+ie),xss(ll-1),xss(ll)
                        nerr=nerr+1
                     endif
                     ll=k1+2*nk+im-1
                     if (xss(ll).lt.zero.or.xss(ll).gt.one) then
                        write(nsyso,&
                          '(''   consis:'',&
                          &'' cumm prob out of range for '',a,&
                          &'' at e='',1p,3e14.6)') name(1:l),&
                          xss(na+ie),xss(ll)
                        nerr=nerr+1
                     endif
                     if (im.gt.1.and.xss(ll).lt.xss(ll-1)) then
                        write(nsyso,&
                           '(''   consis:'',&
                           &'' cumm probs not monotonic for '',a,&
                           &'' at e='',1p,3e14.6)') name(1:l),&
                           xss(na+ie),xss(ll-1),xss(ll)
                        nerr=nerr+1
                     endif
                  enddo
               endif
            enddo
         endif
      endif
   enddo

   !--check energy distributions
   write(nsyso,'(/'' check energy distributions'')')
   if (nr.ne.0) then
      do n=1,nr
         mt=iabs(nint(xss(mtr+n-1)))
         call mtname(mt,name,izai)
         ll=len_trim(name)
         aprime=1
         q=xss(lqr+n-1)
         nlaw=1
         l=nint(xss(ldlw+n-1)+dlw-1)
         locl(1)=l
         do while (nint(xss(l)).ne.0)
            l=nint(xss(l))
            l=l+dlw-1
            nlaw=nlaw+1
            locl(nlaw)=l
         enddo
         l=iabs(nint(xss(tyr+n-1)))
         icm=0
         if (nint(xss(tyr+n-1)).lt.0) icm=1
         if (l.ge.100) then
            l=l-100+dlw-1
            j=nint(xss(l))
            if (j.ne.0) then
               l=l+2*j
            endif
            l=l+1
            j=nint(xss(l))
            l=l+2*j
         endif
         do m=1,nlaw
            l=locl(m)+1
            law=nint(xss(l))
            l=l+2
            j=nint(xss(l))
            if (j.ne.0) then
               l=l+2*j
            endif
            l=l+1
            j=nint(xss(l))
            l=l+2*j+1

            !--law18
            if (law.eq.18) then
               j=nint(xss(l))
               if (j.eq.2) then
                  l=l+1
                  j=nint(xss(l))
               endif

            !--law4
            else if (law.eq.4) then
               j=nint(xss(l))
               if (j.ne.0) then
                  l=l+2*j
               endif
               l=l+1
               ne=nint(xss(l))
               do ie=1,ne
                  e=xss(ie+l)
                  loci=nint(xss(ie+ne+l)+dlw-1)
                  intt=nint(xss(loci))
                  if (intt.ne.1.and.intt.ne.2) then
                     write(nsyso,'(''   consis:'',&
                       &'' illegal interpolation--'',&
                       &''only int=1 and 2 are allowed'')')
                     nerr=nerr+1
                  endif
                  nn=nint(xss(loci+1))
                  loci=loci+1
                  epmax=e
                  if (icm.eq.1) then
                     epmax=(sqrt(e)-sqrt(aprime*e/(aw0+awi)**2))**2
                     epmax=sigfig(epmax,7,-1)
                  endif
                  n2big=0
                  clast=0
                  do j=1,nn
                     ep=xss(j+loci)
                     c=xss(j+2*nn+loci)
                     if (ep.gt.epmax.and.q.lt.zero) then
                        write(nsyso,&
                          '(''   consis:'',&
                          &'' ep.gt.epmax'',1p,e13.6,&
                          &'' with q.lt.0 for '',a,&
                          &'' at e '',e14.6,'' ->'',e13.6)')&
                          epmax,name(1:ll),e,ep
                        n2big=n2big+1
                        nerr=nerr+1
                     endif
                     if (c.lt.zero.or.c.gt.oplus) then
                        write(nsyso,'(''   consis:'',&
                          &'' bad cumm. prob for '',a,&
                          &'' at '',1p,e14.6,'' ->'',e13.6)')&
                        name(1:ll),e,ep
                        nerr=nerr+1
                     endif
                     if (c.lt.clast) then
                        write(nsyso,'(''   consis:'',&
                           &'' decreasing cumm. prob for '',&
                           &a,'' at '',1p,e14.6,'' ->'',e13.6)')&
                           name(1:ll),e,ep
                        nerr=nerr+1
                     endif
                     clast=c
                  enddo
                  if (n2big.gt.0) then
                     write(nsyso,'(''   consis:'',&
                       &'' shifting eprimes greater than epmax'',&
                       &'' and renorming the distribution'')')
                     do j=nn-n2big+1,nn
                        ishift=j-nn-1
                        ep=xss(j+loci)
                        xss(j+loci)=sigfig(epmax,7,ishift)
                        if (intt.eq.1) then
                           p=(xss(j+2*nn+loci)-xss(j-1+2*nn+loci))&
                             /(xss(j+loci)-xss(j-1+loci))
                           xss(j-1+nn+loci)=p
                           xss(j+nn+loci)=p
                        else
                           p=2*(xss(j+2*nn+loci)-xss(j-1+2*nn+loci))&
                             /(xss(j+loci)-xss(j-1+loci))&
                             -xss(j-1+nn+loci)
                           xss(j+nn+loci)=p
                        endif
                     enddo
                  endif
               enddo

            !--law44
            else if (law.eq.44) then
               j=nint(xss(l))
               if (j.ne.0) then
                  l=l+2*j
               endif
               l=l+1
               ne=nint(xss(l))
               do ie=1,ne
                  e=xss(ie+l)
                  loci=nint(xss(ie+ne+l)+dlw-1)
                  intt=nint(xss(loci))
                  if (intt.ne.1.and.intt.ne.2) then
                     write(nsyso,'(''   consis:'',&
                       &'' illegal interpolation--'',&
                       &''only int=1 and 2 are allowed'')')
                     nerr=nerr+1
                  endif
                  nn=nint(xss(loci+1))
                  loci=loci+1
                  epmax=e
                  if (icm.eq.1) then
                     epmax=(sqrt(e)-sqrt(aprime*e/(aw0+awi)**2))**2
                     epmax=sigfig(epmax,7,-1)
                  endif
                  n2big=0
                  clast=0
                  do j=1,nn
                     ep=xss(j+loci)
                     c=xss(j+2*nn+loci)
                     r=xss(j+3*nn+loci)
                     if (ep.gt.epmax) then
                        if (mt.ne.5.and.q.lt.0) then
                           write(nsyso,'(''   consis:'',&
                             &'' ep.gt.epmax'',1p,e13.6,&
                             &'' with q.lt.0 for '',a,&
                             &'' at e'',e14.6,'' ->'',e13.6)')&
                             epmax,name(1:ll),e,ep
                           n2big=n2big+1
                           nerr=nerr+1
                        else if (mt.eq.5.and.aw0.lt.180.) then
                           write(nsyso,'(''   consis:'',&
                             &'' ep.gt.epmax'',1p,e13.6,&
                             &'' with q.lt.0 for '',a,&
                             &'' at e'',e14.6,'' ->'',e13.6)')&
                             epmax,name(1:ll),e,ep
                           write(nsyso,'(''   consis:'',&
                             &''   awr.lt.180'',&
                             &''---this is probably an error.'')')
                           n2big=n2big+1
                           nerr=nerr+1
                        else if (mt.eq.5.and.aw0.ge.180.) then
                           write(nsyso,'(''   consis:'',&
                             &'' ep.gt.epmax'',1p,e13.6,&
                             &'' with q.lt.0 for '',a,&
                             &'' at e'',e14.6,'' ->'',e13.6)')&
                             epmax,name(1:ll),e,ep
                           write(nsyso,&
                             &'(''   consis: awr.ge.180---'',&
                             &''there could be a legitimate'',&
                             &'' positive-q channel'',&
                             &'' or admixed fission.'')')
                           nerr=nerr+1
                        endif
                     endif
                     if (c.lt.zero.or.c.gt.oplus) then
                        write(nsyso,'(''   consis:'',&
                          &'' bad cumm. prob. for '',a,&
                          &''at'',1p,e14.6,'' ->'',e13.6)')&
                          name(1:ll),e,ep
                        nerr=nerr+1
                     endif
                     if (c.lt.clast) then
                        write(nsyso,'(''   consis:'',&
                          &'' decreasing cumm. prob for '',a,&
                          &'' at '',1p,e14.6,'' ->'',e13.6)')&
                          name(1:ll),e,ep
                        nerr=nerr+1
                     endif
                     clast=c
                     if (r.lt.zero.or.r.gt.oneup) then
                        write(nsyso,'(''   consis:'',&
                          &'' bad kalbach r for '',a,&
                          &''at'',1p,e14.6,'' ->'',e13.6)')&
                          name(1:ll),e,ep
                        nerr=nerr+1
                     endif
                  enddo
                  if (n2big.gt.0) then
                     write(nsyso,'(''   consis:'',&
                       &'' shifting eprimes greater than epmax'',&
                       &'' and renorming the distribution'')')
                     do j=nn-n2big+1,nn
                        ishift=j-nn-1
                        ep=xss(j+loci)
                        xss(j+loci)=sigfig(epmax,7,ishift)
                        if (intt.eq.1) then
                           p=(xss(j+2*nn+loci)-xss(j-1+2*nn+loci))&
                            /(xss(j+loci)-xss(j-1+loci))
                           xss(j-1+nn+loci)=p
                           xss(j+nn+loci)=p
                        else
                           p=2*(xss(j+2*nn+loci)-xss(j-1+2*nn+loci))&
                            /(xss(j+loci)-xss(j-1+loci))&
                             -xss(j-1+nn+loci)
                           xss(j+nn+loci)=p
                        endif
                     enddo
                  endif
               enddo

            !--law61
            else if (law.eq.61) then
               j=nint(xss(l))
               if (j.ne.0) then
                  l=l+2*j
               endif
               l=l+1
               ne=nint(xss(l))
               do ie=1,ne
                  e=xss(ie+l)
                  loci=nint(xss(ie+ne+l))+dlw-1
                  nn=nint(xss(loci+1))
                  loci=loci+1
                  clast=0
                  do j=1,nn
                     ep=xss(j+loci)
                     c=xss(j+2*nn+loci)
                     if (c.lt.zero.or.c.gt.oplus) then
                        write(nsyso,'(''   consis:'',&
                          &'' bad cumm. prob. for '',a,&
                          &'' at'',1p,e14.6,'' ->'',e13.6)')&
                          name(1:ll),e,ep
                        nerr=nerr+1
                     endif
                     if (c.lt.clast) then
                        write(nsyso,'(''   consis:'',&
                          &'' decreasing cumm. prob for '',a,&
                          &'' at '',1p,e14.6,'' ->'',e13.6)')&
                          name(1:ll),e,ep
                        nerr=nerr+1
                     endif
                     clast=c
                     locj=nint(xss(j+3*nn+loci)+dlw-1)
                     nmu=nint(xss(locj+1))
                     cclast=0
                     do k=1,nmu
                        co=xss(locj+1+k)
                        cc=xss(locj+1+2*nmu+k)
                        if (cc.lt.zero.or.cc.gt.oplus) then
                           write(nsyso,'('' consis:'',&
                             &'' bad angular cumm. prob. for '',a,&
                             &''at'',1p,e14.6,'' ->'',e13.6,e14.6)')&
                             name(1:ll),e,ep,co
                           nerr=nerr+1
                        endif
                        if (cc.lt.cclast) then
                           write(nsyso,'(''   consis:'',&
                             &'' decreasing angular'',&
                             &'' cumm. prob for '',a,&
                             &'' at '',1p,e14.6,'' ->'',e13.6,e14.6)')&
                             name(1:ll),e,ep,co
                           nerr=nerr+1
                        endif
                        cclast=cc
                     enddo
                  enddo
               enddo
            endif
         enddo
      enddo
   endif

   !--check unresolved-range probability tables
   if (iurpt.ne.0) then
      write(nsyso,'(/'' check probability tables''/&
      &'' ------------------------'')')
      nure=nint(xss(iurpt))
      nurb=nint(xss(iurpt+1))
      lurt=nint(xss(iurpt+2))
      luri=nint(xss(iurpt+3))
      lura=nint(xss(iurpt+4))
      lurf=nint(xss(iurpt+5))
      write(nsyso,'(/''     number of energies: '',i6/&
                    &''     number of bins:     '',i6/&
                    &''     interpolation law:  '',i6/&
                    &''     inelastic reaction: '',i6/&
                    &''     absorption reaction:'',i6)')&
                    & nure,nurb,lurt,luri,lura
      if (lurf.eq.0) write(nsyso,'(&
                    &''     tables are cross sections (lssf=0)'')')
      if (lurf.eq.1) write(nsyso,'(&
                    &''     tables are factors (lssf=1)'')')
      do ie=1,nure
         write(nsyso,'(/'' energy='',1p,e14.6)') xss(iurpt+5+ie)
         nerrtb=0
         nerrxt=0
         nerrxe=0
         nerrxf=0
         nerrxg=0
         nerrxh=0
         ll=iurpt+5+nure+(ie-1)*6*nurb
         do ib=1,nurb
           if (xss(ll+ib).lt.0.0d0) then
             write(nsyso,'('' consis: negative probability value='',&
             & 1pe14.6,'' for bin '',i4)')xss(ll+ib),ib
             nerrtb=nerrtb+1
           endif
           if (xss(ll+nurb+ib).lt.0.0d0) then
             write(nsyso,'('' consis: negative total cross section/factor='',&
             & 1pe14.6,'' for bin '',i4)')xss(ll+nurb+ib),ib
             nerrxt=nerrxt+1
           endif
           if (xss(ll+2*nurb+ib).lt.0.0d0) then
             write(nsyso,'('' consis: negative elastic cross section/factor='',&
             & 1pe14.6,'' for bin '',i4)')xss(ll+2*nurb+ib),ib
             nerrxe=nerrxe+1
           endif
           if (xss(ll+3*nurb+ib).lt.0.0d0) then
             write(nsyso,'('' consis: negative fission cross section/factor='',&
             & 1pe14.6,'' for bin '',i4)')xss(ll+3*nurb+ib),ib
             nerrxf=nerrxf+1
           endif
           if (xss(ll+4*nurb+ib).lt.0.0d0) then
             write(nsyso,'('' consis: negative capture cross section/factor='',&
             & 1pe14.6,'' for bin '',i4)')xss(ll+4*nurb+ib),ib
             nerrxg=nerrxg+1
           endif
           if (xss(ll+5*nurb+ib).lt.0.0d0) then
             write(nsyso,'('' consis: negative heating cross section/factor='',&
             & 1pe14.6,'' for bin '',i4)')xss(ll+5*nurb+ib),ib
             nerrxh=nerrxh+1
           endif
         enddo
         nerr=nerr+nerrtb+nerrxt+nerrxe+nerrxf+nerrxg+nerrxh
      enddo
   endif

   !--check delayed neutron data
   if (ndnf.gt.0) then
      write(nsyso,'(/'' check delayed neutron fractions'')')
      l=dndat
      sum=0
      do j=1,ndnf
         l=l+1
         nj=nint(xss(l))
         l=l+1+2*nj
         nn=nint(xss(l))
         l=l+nn
         frac=xss(l+1)
         sum=sum+frac
         l=l+nn+1
      enddo
      if (abs(sum-1)*1000.gt.one) then
         write(nsyso,'(''   consis: delayed fractions do not'',&
           &'' sum to one'')')
         nerr=nerr+1
      endif
      write(nsyso,'(/'' check delayed neutron distributions'')')
      do i=1,ndnf
         nlaw=1
         loct=nint(xss(i-1+ldnd)+dnd-1)
         law=nint(xss(loct+1))
         m=nint(xss(loct+3))
         loct=loct+3+2*m
         loct=loct+1
         n=nint(xss(loct))
         loct=loct+1+2*n
         m=nint(xss(loct))
         loct=loct+2*m
         loct=loct+1
         ne=nint(xss(loct))
         loci=nint(xss(1+ne+loct))+dnd-1
         intt=nint(xss(loci))
         n=nint(xss(loci+1))
         loci=loci+1
         do j=1,n
            x=xss(j+loci)
            y=xss(j+loci+n)
            c=xss(j+loci+2*n)
            if (j.gt.1) then
               if (x.lt.xlast) then
                  write(nsyso,'(''   consis: delayed spectrum'',&
                    &'' energies not monotonic'')')
                  nerr=nerr+1
               endif
               if (c.lt.clast) then
                  write(nsyso,'(''   consis: delayed spectrum'',&
                    &'' cummulative probs not monotonic'')')
                  nerr=nerr+1
               endif
            endif
            xlast=x
            clast=c
         enddo
      enddo
   endif

   !--check detailed photon data
   if (ntrp.ne.0) then

      !--total photon production cross section
      write(nsyso,'(/'' check photon production sum'')')
      do ie=1,nes
         e=xss(esz+ie-1)
         gsum=0
         do i=1,ntrp
            loct=nint(xss(i-1+lsigp)+sigp-1)
            loc1=nint(xss(i-1+ldlwp)+dlwp-1)
            mftype=nint(xss(loct))
            if (mftype.eq.13) then
               loct=loct+1
               j=ie-nint(xss(loct))+1
               if (j.ge.1.and.j.le.nint(xss(loct+1))) then
                  ss=xss(loct+1+j)
                  gsum=gsum+ss
               endif
            else
               loct=loct+1
               mtmult=nint(xss(loct))
               do ii=1,ntr
                  if (mtmult.eq.nint(xss(mtr+ii-1))&
                    .or.(mtmult.eq.18.and.nint(xss(mtr+ii-1)).eq.19)&
                    .or.(mtmult.eq.18.and.nint(xss(mtr+ii-1)).eq.20)&
                    .or.(mtmult.eq.18.and.nint(xss(mtr+ii-1)).eq.21)&
                    .or.(mtmult.eq.18.and.nint(xss(mtr+ii-1)).eq.38)&
                    ) then
                     k=nint(xss(lsig+ii-1)+sig-1)
                     iaa=nint(xss(k))
                     naa=nint(xss(k+1))
                     if (ie.ge.iaa.and.ie.lt.iaa+naa) then
                        locv=loct+1
                        locc=locv
                        m=nint(xss(locv))
                        locv=locv+1+2*m
                        n=nint(xss(locv))
                        law=nint(xss(loc1+1))
                        iflag=0
                        if (law.ne.2) iflag=1
                        if (m.ne.0.or.n.gt.2) iflag=1
                        if (xss(locv+3).ne.xss(locv+4)) iflag=1
                        if (iflag.eq.0) then
                           gsum=gsum+xss(locv+3)*xss(2+k+ie-iaa)
                        else
                           call terpc(e,y,xss(locc))
                           gsum=gsum+y*xss(2+k+ie-iaa)
                        endif
                     endif
                  endif
               enddo
            endif
         enddo
         gg=xss(gpd+ie-1)
         if (abs(gg-gsum).gt.gg/10000) then
            nerr=nerr+1
            write(nsyso,'(''   consis: mismatch at'',&
              &1p,e14.6,''  gpd='',e14.6,&
              &''  sum='',e14.6)') e,gg,gsum
         endif
      enddo

      !--photon production distributions
      write(nsyso,'(/'' check photon distributions'')')
      do i=1,ntrp
         loct=nint(xss(i-1+ldlwp)+dlwp-1)
         law=nint(xss(loct+1))
         write(name,'(i6)') nint(xss(i-1+mtrp))
         l=len_trim(name)
         if (law.eq.4) then
            loct=loct+3
            m=nint(xss(loct))
            if (m.ne.0) then
               loct=loct+2*m
            endif
            loct=loct+1
            n=nint(xss(loct))
            loct=loct+1+2*n
            m=nint(xss(loct))
            if (m.ne.0) then
               loct=loct+2*m
            endif
            loct=loct+1
            ne=nint(xss(loct))
            do ie=1,ne
               e=xss(ie+loct)
               loci=nint(xss(ie+ne+loct)+dlwp-1)
               n=nint(xss(loci+1))
               loci=loci+1
               clast=0
               do j=1,n
                  ep=xss(j+loci)
                  c=xss(j+2*n+loci)
                  if (c.lt.zero.or.c.gt.oplus) then
                     write(nsyso,'(''   consis:'',&
                       &'' bad cumm. prob. for '',a,&
                       &''at'',1p,e14.6,'' ->'',e13.6)')&
                       name(1:l),e,ep
                     nerr=nerr+1
                  endif
                  if (c.lt.clast) then
                     write(nsyso,'(''   consis:'',&
                       &''decreasing cumm. prob for '',&
                       &a,'' at '',1p,e14.6,'' ->'',e13.6)')&
                       name(1:l),e,ep
                     nerr=nerr+1
                  endif
                  clast=c
               enddo
            enddo
         endif
      enddo
   endif

   !--do consistency checks for particle production
   if (ntype.ne.0) then
      write(nsyso,'(/'' checking particle production sections'')')
      do i=1,ntype
         ipt=nint(xss(ptype+i-1))
         hpd=nint(xss(ploct+10*(i-1)))
         iaa=nint(xss(hpd))
         if (ipt.eq.1) write(nsyso,'(/''   neutron production:'')')
         if (ipt.eq.9) write(nsyso,'(/''   proton production:'')')
         if (ipt.eq.31) write(nsyso,'(/''   deuteron production:'')')
         if (ipt.eq.32) write(nsyso,'(/''   triton production:'')')
         if (ipt.eq.33) write(nsyso,'(/''   he-3 production:'')')
         if (ipt.eq.34) write(nsyso,'(/''   alpha production:'')')
         write(nsyso,'(/''   checking energy distributions'')')
         mtrh=nint(xss(ploct+10*(i-1)+1))
         nmtr=nint(xss(ntro+i-1))
         lsigh=nint(xss(ploct+10*(i-1)+3))
         sigh=nint(xss(ploct+10*(i-1)+4))
         ldlwh=nint(xss(ploct+10*(i-1)+7))
         dlwh=nint(xss(ploct+10*(i-1)+8))
         do imt=1,nmtr
            mt=nint(xss(mtrh+imt-1))
            call mtname(mt,name,izai)
            l=len_trim(name)
            l1=nint(xss(lsigh+imt-1))
            l2=sigh+l1-1
            ne=nint(xss(l2+3))
            l1=nint(xss(ldlwh+imt-1))
            l2=dlwh+l1-1
            law=nint(xss(l2+1))
            l3=dlwh+nint(xss(l2+2))-1

            !--law=4
            if (law.eq.4) then
               j=nint(xss(l3))
               if (j.ne.0) then
                  l3=l3+2*j
               endif
               l3=l3+1
               ne=nint(xss(l3))
               do ie=1,ne
                  e=xss(l3+ie)
                  loci=nint(xss(l3+ne+ie))+dlwh-1
                  nn=nint(xss(loci+1))
                  loci=loci+1
                  clast=0
                  do j=1,nn
                     ep=xss(j+loci)
                     c=xss(j+2*nn+loci)
                     if (c.lt.zero.or.c.gt.oplus) then
                        write(nsyso,'(''   consis:'',&
                          &'' bad law4 cumm. prob for '',a,&
                          &'' at '',1p,e14.6,'' ->'',e13.6)')&
                          name(1:l),e,ep
                        nerr=nerr+1
                     endif
                     if (c.lt.clast) then
                        write(nsyso,'(''   consis:'',&
                          &'' decreasing law4'',&
                          &'' cumm. prob for '',a,&
                          &'' at '',1p,e14.6,'' ->'',e13.6)')&
                          name(1:l),e,ep
                        nerr=nerr+1
                     endif
                     clast=c
                  enddo
               enddo

            !--law=44
            else if (law.eq.44) then
               j=nint(xss(l3))
               if (j.ne.0) then
                  l3=l3+2*j
               endif
               l3=l3+1
               ne=nint(xss(l3))
               do ie=1,ne
                  e=xss(l3+ie)
                  loci=nint(xss(l3+ne+ie))+dlwh-1
                  nn=nint(xss(loci+1))
                  loci=loci+1
                  clast=0
                  do j=1,nn
                     ep=xss(j+loci)
                     c=xss(j+2*nn+loci)
                     r=xss(j+3*nn+loci)
                     if (c.lt.zero.or.c.gt.oplus) then
                        write(nsyso,'(''   consis:'',&
                           &''   bad law44 cumm. prob for '',a,&
                           &'' at '',1p,e14.6,'' ->'',e13.6)')&
                           name(1:l),e,ep
                        nerr=nerr+1
                     endif
                     if (c.lt.clast) then
                        write(nsyso,'(''   consis:'',&
                          &'' decreasing law44'',&
                          &'' cumm. prob for '',a,&
                          &'' at '',1p,e14.6,'' ->'',e13.6)')&
                          name(1:l),e,ep
                        nerr=nerr+1
                     endif
                     if (r.lt.zero.or.r.gt.oplus) then
                        write(nsyso,&
                          '(''   consis:'',&
                          &'' bad law44 kalbach r for '',a,&
                          &'' at '',1p,e14.6,'' ->'',e13.6)')&
                          name(1:l),e,ep
                        nerr=nerr+1
                     endif
                     clast=c
                  enddo
               enddo

            !--law=61
            else if (law.eq.61) then
               j=nint(xss(l3))
               if (j.ne.0) then
                  l3=l3+2*j
               endif
               l3=l3+1
               ne=nint(xss(l3))
               do ie=1,ne
                  e=xss(l3+ie)
                  loci=nint(xss(l3+ne+ie))+dlwh-1
                  nn=nint(xss(loci+1))
                  loci=loci+1
                  clast=0
                  do j=1,nn
                     ep=xss(j+loci)
                     c=xss(j+2*nn+loci)
                     if (c.lt.zero.or.c.gt.oplus) then
                        write(nsyso,'('' consis:'',&
                          &'' bad law61 cumm. prob for '',a,&
                          &'' at '',1p,e14.6,'' ->'',e13.6)')&
                          name(1:l),e,ep
                        nerr=nerr+1
                     endif
                     if (c.lt.clast) then
                        write(nsyso,'(''   consis:'',&
                          &'' decreasing law61 cumm. prob for '',&
                          &a,'' at '',1p,e14.6,'' ->'',e13.6)')&
                          name(1:l),e,ep
                        nerr=nerr+1
                     endif
                     clast=c
                     locj=nint(xss(j+3*nn+loci)+dlwh-1)
                     nmu=nint(xss(locj+1))
                     cclast=0
                     do k=1,nmu
                        co=xss(locj+1+k)
                        cc=xss(locj+1+2*nmu+k)
                        if (cc.lt.zero.or.cc.gt.oplus) then
                           write(nsyso,'(''   consis:'',&
                             &'' bad law61 angular'',&
                             &'' cumm. prob. for '',a,&
                             &'' at'',1p,e14.6,'' ->'',e13.6,e14.6)')&
                             name(1:l),e,ep,co
                           nerr=nerr+1
                        endif
                        if (cc.lt.cclast) then
                           write(nsyso,'(''   consis:'',&
                             &'' decreasing angular'',&
                             &'' cumm. prob for '',a,&
                             &'' at '',1p,e14.6,'' ->'',e13.6,e14.6)')&
                             name(1:l),e,ep,co
                           nerr=nerr+1
                        endif
                        cclast=cc
                     enddo
                  enddo
               enddo
            endif
         enddo
      enddo
   endif

   if (nerr.eq.0) then
      write(nsyso,'(/'' no problems found'')')
   else
      write(nsyso,'(/i5,'' problems found'')') nerr
      call mess('consis','consistency problems found',' ')
   endif

   return
   end subroutine consis

   subroutine terpc(x,y,a)
   !-------------------------------------------------------------------
   ! Interpolate in an mftype=12 or 16 ACE table
   !-------------------------------------------------------------------
   use endf ! provides terp1
   ! externals
   real(kr)::x,y,a(*)
   ! internals
   integer::i,m,mm,jj,n,j

   y=0
   i=1
   m=nint(a(i))
   i=i+1
   if (m.ne.0) i=i+2*m
   mm=2
   jj=0
   if (m.gt.0) then
      jj=1
      mm=nint(a(1+jj+m))
   endif
   n=nint(a(i))
   do j=1,n-1
      if (m.gt.0.and.j.ge.nint(a(1+jj))) then
         jj=jj+1
         mm=nint(a(1+jj+m))
      endif
      if (a(i+j).le.x.and.a(i+j+1).ge.x) then
         call terp1(a(i+j),a(i+n+j),a(i+j+1),a(i+j+1+n),x,y,mm)
      endif
   enddo
   return
   end subroutine terpc

   subroutine aplots(nout,hk)
   !-------------------------------------------------------------------
   ! Do standard plots of an ACE file
   !-------------------------------------------------------------------
   use util ! provides openz
   use endf ! provides terp1
   use acecm ! provides mtname
   ! externals
   integer::nout
   character(70)::hk
   ! internals
   integer::ipcol,iwcol,i,it,j,maxii,idone,ii1,ii2,iil,nn
   integer::nnf,mt,kf,iif,kc,iic,nofiss,n,k,iaa
   integer::major,minor,mtl,icurv,mtlast,nlev,iflag
   integer::nure,intunr,nurb,lssf,ie,ib,ii,ll,kk,nunu
   integer,parameter::pltumx=10000
   real(kr)::e,tot,abso,elas,gprod,xtag,ytag,thin,abss
   real(kr)::e1,e2,fiss,cap,heat,dam,x,y,xlast
   real(kr)::xmin,xmax,ymin,ymax,xstep,ystep,test
   real(kr)::ee(pltumx),s0(pltumx),s1(pltumx),s2(pltumx)
   real(kr)::f0,f1,f2,c1,c2,cl,dp,pp,pe,capt
   character(1)::qu=''''
   character(10)::name
   character(70)::strng
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::ten=10.e0_kr
   real(kr),parameter::small=1.e-12_kr
   real(kr),parameter::scale=1.e6_kr
   real(kr),parameter::hmin=1.e-10_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1
   integer,parameter::nden=4000

   !--black and white pages
   ipcol=0
   iwcol=0
   !--colored pages
   ipcol=2
   iwcol=3

   !--start the viewr input text
   call openz(nout,1)
   write(nout,'(''1 2 .30'',i3,''/'')') ipcol

   !--plot log-log total, absorption, elastic, and gamma production
   xmin=big
   xmax=0
   ymin=big
   ymax=-big
   do i=1,nes
      e=xss(esz-1+i)
      tot=xss(esz+nes-1+i)
      abso=xss(esz+2*nes-1+i)
      elas=xss(esz+3*nes-1+i)
      if (e.lt.xmin) xmin=e
      if (e.gt.xmax) xmax=e
      if (tot.lt.ymin) ymin=tot
      if (tot.gt.ymax) ymax=tot
      if (abso.lt.ymin) ymin=abso
      if (abso.gt.ymax) ymax=abso
      if (elas.lt.ymin) ymin=elas
      if (elas.gt.ymax) ymax=elas
      if (gpd.ne.0) then
         gprod=xss(gpd-1+i)
         if (gprod.lt.ymin) ymin=gprod
         if (gprod.gt.ymax) ymax=gprod
      endif
   enddo
   call ascll(xmin,xmax)
   if (ymin.lt.ymax/scale) ymin=ymax/scale
   call ascll(ymin,ymax)
   write(nout,'(''1'',i3,''/'')') iwcol
   it=1
   do i=1,70
      if (hk(i:i).ne.' ') it=i
   enddo
   write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
   write(nout,'(a,''<p>rincipal cross sections'',a,''/'')') qu,qu
   xtag=5*xmin
   ytag=7*log10(ymin)/10+3*log10(ymax)/10
   ytag=ten**ytag
   write(nout,'(''4 0 2 1'',2e12.4,''/'')') xtag,ytag
   write(nout,'(1p,3e12.3,''/'')') xmin,xmax,one
   write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
   write(nout,'(1p,3e12.3,''/'')') ymin,ymax,one
   write(nout,'(a,''<c>ross section (barns)'',a,''/'')') qu,qu
   write(nout,'(''/'')')
   write(nout,'(''/'')')
   write(nout,'(a,''total'',a,''/'')') qu,qu
   write(nout,'(''0/'')')
   thin=ten**(log10(xmax/xmin)/nden)
   xlast=small
   j=0
   do i=1,nes
      e=xss(esz-1+i)
      tot=xss(esz+nes-1+i)
      if (tot.lt.ymin) tot=ymin
      if (nes.le.nden.or.e.ge.thin*xlast) then
         j=j+1
         write(nout,'(1p,2e14.6,''/'')') e,tot
         xlast=e
      endif
   enddo
   write(nout,'(''/'')')
   write(nout,'(''2/'')')
   write(nout,'(''/'')')
   if (iwcol.eq.0) then
      write(nout,'(''0 0 1/'')')
   else
      write(nout,'(''0 0 0 3/'')')
   endif
   write(nout,'(a,''absorption'',a,''/'')') qu,qu
   write(nout,'(''0/'')')
   xlast=small
   j=0
   do i=1,nes
      e=xss(esz-1+i)
      abss=xss(esz+2*nes-1+i)
      if (abss.lt.ymin) abss=ymin
      if (nes.le.nden.or.e.ge.thin*xlast) then
         j=j+1
         write(nout,'(1p,2e14.6,''/'')') e,abss
         xlast=e
      endif
   enddo
   write(nout,'(''/'')')
   write(nout,'(''3/'')')
   write(nout,'(''/'')')
   if (iwcol.eq.0) then
      write(nout,'(''0 0 2/'')')
   else
      write(nout,'(''0 0 0 2/'')')
   endif
   write(nout,'(a,''elastic'',a,''/'')') qu,qu
   write(nout,'(''0/'')')
   xlast=small
   j=0
   do i=1,nes
      e=xss(esz-1+i)
      elas=xss(esz+3*nes-1+i)
      if (elas.lt.ymin) elas=ymin
      if (nes.le.nden.or.e.ge.thin*xlast) then
         j=j+1
         write(nout,'(1p,2e14.6,''/'')') e,elas
         xlast=e
      endif
   enddo
   write(nout,'(''/'')')
   if (gpd.ne.0) then
      write(nout,'(''4/'')')
      write(nout,'(''/'')')
      if (iwcol.eq.0) then
         write(nout,'(''0 0 3/'')')
      else
         write(nout,'(''0 0 0 1/'')')
      endif
      write(nout,'(a,''gamma production'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      xlast=small
      j=0
      do i=1,nes
         e=xss(esz-1+i)
         gprod=xss(gpd-1+i)
         if (gprod.lt.ymin) gprod=ymin
         if (nes.le.nden.or.e.ge.thin*xlast) then
            j=j+1
            write(nout,'(1p,2e14.6,''/'')') e,gprod
            xlast=e
         endif
      enddo
      write(nout,'(''/'')')
   endif

   !--plot expanded resonance data for total
   if (nes.ge.1500) then
      maxii=120
      idone=0
      e1=xss(esz)
      ii1=1
      do while (idone.eq.0)
         do while (idone.eq.0)
            ii2=ii1
            e2=10*e1
            iil=nint(log10(e2))
            e2=ten**iil
            do while (ii2.lt.nes.and.xss(esz+ii2-1).lt.e2+small*e2)
               ii2=ii2+1
            enddo
            if (ii2-ii1.gt.maxii) idone=1
            if (ii2.ge.nes) idone=2
            ii2=ii2-1
            if (idone.eq.0) then
               e1=e2
               ii1=ii2
            endif
         enddo
         if (idone.eq.1) then
            e2=xss(esz-1+ii2)
            nn=ii2-ii1+1
            xmin=big
            xmax=small
            ymin=big
            ymax=small
            do i=ii1,ii2
               e=xss(esz-1+i)
               tot=xss(esz+nes-1+i)
               if (e.lt.xmin) xmin=e
               if (e.gt.xmax) xmax=e
               if (tot.gt.zero.and.tot.lt.ymin) ymin=tot
               if (tot.gt.zero.and.tot.gt.ymax) ymax=tot
            enddo
            call ascll(xmin,xmax)
            if (ymin.lt.ymax/scale) ymin=ymax/scale
            call ascll(ymin,ymax)
            write(nout,'(''1'',i3,''/'')') iwcol
            it=1
            do i=1,70
               if (hk(i:i).ne.' ') it=i
            enddo
            write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
            write(nout,&
              '(a,''resonance total cross section'',a,''/'')') qu,qu
            write(nout,'(''4 0 2 1/'')')
            write(nout,'(1p,3e12.3,''/'')') xmin,xmax,one
            write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
            write(nout,'(1p,3e12.3,''/'')') ymin,ymax,one
            write(nout,&
              '(a,''<c>ross section (barns)'',a,''/'')') qu,qu
            write(nout,'(''/'')')
            write(nout,'(''/'')')
            write(nout,'(a,''total'',a,''/'')') qu,qu
            write(nout,'(''0/'')')
            thin=ten**(log10(xmax/xmin)/nden)
            xlast=small
            j=0
            do i=ii1,ii2
               e=xss(esz-1+i)
               if (nn.le.nden.or.e.ge.thin*xlast) then
                  tot=xss(esz+nes-1+i)
                  if (tot.lt.ymin) tot=ymin
                  j=j+1
                  write(nout,'(1p,2e14.6,''/'')') e,tot
                  xlast=e
               endif
            enddo
            write(nout,'(''/'')')
            if (ii2.lt.nes-200.and.e2.lt.ten) then
               e1=e2
               ii1=ii2
               idone=0
            endif
         endif
      enddo
   endif

   !--plot expanded resonance data for fission and capture
   if (nes.ge.1500) then
      maxii=100
      idone=0
      e1=xss(esz)
      ii1=1
      do while (idone.eq.0)
         do while (idone.eq.0)
            ii2=ii1
            e2=10*e1
            iil=nint(log10(e2))
            e2=ten**iil
            do while (ii2.lt.nes.and.xss(esz+ii2-1).lt.e2+small*e2)
               ii2=ii2+1
            enddo
            if (ii2-ii1.gt.maxii) idone=1
            if (ii2.ge.nes) idone=2
            ii2=ii2-1
            if (idone.eq.0) then
               e1=e2
               ii1=ii2
            endif
         enddo
         if (idone.eq.1) then
            e2=xss(esz-1+ii2)
            nn=ii2-ii1+1
            nnf=0
            do i=1,ntr
               mt=nint(xss(mtr+i-1))
               if (mt.eq.18.or.mt.eq.19) then
                  kf=nint(xss(lsig+i-1)+sig-1)
                  nnf=nint(xss(kf+1))
                  iif=nint(xss(kf))
               else if (mt.eq.102) then
                  kc=nint(xss(lsig+i-1)+sig-1)
                  iic=nint(xss(kc))
               endif
            enddo
            xmin=big
            xmax=small
            ymin=big
            ymax=small
            nofiss=1
            do i=ii1,ii2
               e=xss(esz-1+i)
               fiss=0
               if (nnf.gt.0) then
                  if (i.ge.iif) fiss=xss(kf+2+i-iif)
               endif
               if (fiss.gt.zero) nofiss=0
               cap=xss(kc+2+i-iic)
               if (e.lt.xmin) xmin=e
               if (e.gt.xmax) xmax=e
               if (fiss.gt.zero.and.fiss.lt.ymin) ymin=fiss
               if (fiss.gt.zero.and.fiss.gt.ymax) ymax=fiss
               if (cap.gt.zero.and.cap.lt.ymin) ymin=cap
               if (cap.gt.zero.and.cap.gt.ymax) ymax=cap
            enddo
            call ascll(xmin,xmax)
            if (ymin.lt.ymax/scale) ymin=ymax/scale
            call ascll(ymin,ymax)
            write(nout,'(''1'',i3,''/'')') iwcol
            it=1
            do i=1,70
               if (hk(i:i).ne.' ') it=i
            enddo
            write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
            write(nout,&
              '(a,''resonance absorption cross sections'',a,''/'')')&
              qu,qu
            write(nout,'(''4 0 2 1/'')')
            write(nout,'(1p,3e12.3,''/'')') xmin,xmax,one
            write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
            write(nout,'(1p,3e12.3,''/'')') ymin,ymax,one
            write(nout,&
              '(a,''<c>ross section (barns)'',a,''/'')') qu,qu
            write(nout,'(''/'')')
            write(nout,'(''/'')')
            write(nout,'(a,''capture'',a,''/'')') qu,qu
            write(nout,'(''0/'')')
            thin=ten**(log10(xmax/xmin)/nden)
            xlast=small
            j=0
            do i=ii1,ii2
               e=xss(esz-1+i)
               if (nn.le.nden.or.e.ge.thin*xlast) then
                  cap=xss(kc+2+i-iic)
                  if (cap.lt.ymin) cap=ymin
                  j=j+1
                  write(nout,'(1p,2e14.6,''/'')') e,cap
                  xlast=e
               endif
            enddo
            write(nout,'(''/'')')
            if (nofiss.ne.1) then
               write(nout,'(''2/'')')
               write(nout,'(''/'')')
               write(nout,'(''0 0 1/'')')
               write(nout,'(a,''fission'',a,''/'')') qu,qu
               write(nout,'(''0/'')')
               xlast=small
               j=0
               do i=ii1,ii2
                  e=xss(esz-1+i)
                  if (nn.le.nden.or.e.ge.thin*xlast) then
                     fiss=0
                     if (i.ge.iif) fiss=xss(kf+2+i-iif)
                     if (fiss.lt.ymin) fiss=ymin
                     j=j+1
                     write(nout,'(1p,2e14.6,''/'')') e,fiss
                     xlast=e
                  endif
               enddo
               write(nout,'(''/'')')
            endif
            if (ii2.lt.nes-200.and.e2.lt.ten) then
               e1=e2
               ii1=ii2
               idone=0
            endif
         endif
      enddo
   endif

   !--plot ur cross sections
   if (iurpt.ne.0) then
      nure=nint(xss(iurpt))
      if (nure.gt.pltumx) then
         write(strng,'(''1need to redefine pltumx to'',i5)')nure
         call error('aplots',strng,'')
      endif
      intunr=nint(xss(iurpt+2))
      nurb=nint(xss(iurpt+1))
      lssf=nint(xss(iurpt+5))

      !--total
      if (lssf.eq.0) then
         xmin=1e10
         xmax=0
         ymin=1e10
         ymax=0
         do ie=1,nure
            ee(ie)=abs(xss(iurpt+5+ie))
            if (ee(ie).gt.xmax) xmax=ee(ie)
            if (ee(ie).lt.xmin) xmin=ee(ie)
            s0(ie)=0
            s1(ie)=0
            s2(ie)=0
            f0=0
            f1=0
            f2=0
            cl=0
            ll=iurpt+5+nure+(ie-1)*6*nurb
            do ib=1,nurb
               dp=xss(ll+ib)-cl
               s0(ie)=s0(ie)+dp*xss(ll+nurb+ib)
               s1(ie)=s1(ie)+dp*xss(ll+nurb+ib)/(100+xss(ll+nurb+ib))
               s2(ie)=s2(ie)+dp*xss(ll+nurb+ib)/(1+xss(ll+nurb+ib))
               f0=f0+dp
               f1=f1+dp/(100+xss(ll+nurb+ib))
               f2=f2+dp/(1+xss(ll+nurb+ib))
               cl=xss(ll+ib)
            enddo
            s1(ie)=s1(ie)/f1
            s2(ie)=s2(ie)/f2
            if (s0(ie).gt.ymax) ymax=s0(ie)
            if (s1(ie).gt.ymax) ymax=s1(ie)
            if (s2(ie).gt.ymax) ymax=s1(ie)
            if (s0(ie).lt.ymin) ymin=s0(ie)
            if (s1(ie).lt.ymin) ymin=s1(ie)
            if (s2(ie).lt.ymin) ymin=s2(ie)
         enddo
         nunu=nure
      else
         xmin=abs(xss(iurpt+5+1))
         xmax=abs(xss(iurpt+5+nure))
         ymin=1e10
         ymax=0
         ie=0
         do i=1,nes
           e=xss(esz-1+i)
           if (e.lt.xmin.or.e.ge.xmax) cycle
           ie=ie+1
           if (ie.gt.pltumx) cycle
           ee(ie)=e
           tot=xss(esz+nes-1+i)
           s0(ie)=tot
           do ii=1,nure-1
              e1=abs(xss(iurpt+5+ii))
              e2=abs(xss(iurpt+5+ii+1))
              if (e.lt.e1.or.e.ge.e2) cycle
              kk=iurpt+5+nure+(ii-1)*6*nurb
              ll=iurpt+5+nure+(ii+1-1)*6*nurb
              s1(ie)=0
              s2(ie)=0
              f1=0
              f2=0
              c1=0
              c2=0
              do j=1,nurb
                 call terp1(e1,xss(kk+j)-c1,e2,xss(ll+j)-c2,&
                    e,dp,intunr)
                 call terp1(e1,xss(kk+nurb+j),e2,xss(ll+nurb+j),&
                   e,pp,intunr)
                 s1(ie)=s1(ie)+dp*pp*tot/(100+pp*tot)
                 s2(ie)=s2(ie)+dp*pp*tot/(1+pp*tot)
                 f1=f1+dp/(100+pp*tot)
                 f2=f2+dp/(1+pp*tot)
                 c1=xss(kk+j)
                 c2=xss(ll+j)
              enddo
              s1(ie)=s1(ie)/f1
              s2(ie)=s2(ie)/f2
              if (s0(ie).gt.ymax) ymax=s0(ie)
              if (s0(ie).lt.ymin) ymin=s0(ie)
              if (s1(ie).gt.ymax) ymax=s1(ie)
              if (s1(ie).lt.ymin) ymin=s1(ie)
              if (s2(ie).gt.ymax) ymax=s2(ie)
              if (s2(ie).lt.ymin) ymin=s2(ie)
           enddo
         enddo
         if (ie.gt.pltumx) then
            write(strng,'(''2need to redefine pltumx to'',i5)')ie
            call error('aplots',strng,'')
         endif
         nunu=ie
      endif
      ymax=ymax+ymax/10
      if (ymin.lt.ymax/1000) ymin=ymax/1000
      call ascll(xmin,xmax)
      call ascll(ymin,ymax)
      write(nout,'(''1'',i3,''/'')') iwcol
      it=1
      do i=1,70
         if (hk(i:i).ne.' ') it=i
      enddo
      write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
      write(nout,'(a,''UR total cross section'',a,''/'')') qu,qu
      write(nout,'(''4 0 2 1/'')')
      write(nout,'(1p,3e12.3,''/'')') xmin,xmax,one
      write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
      write(nout,'(1p,3e12.3,''/'')') ymin,ymax,one
      write(nout,'(a,''<c>ross section (barns)'',a,''/'')') qu,qu
      write(nout,'(''/'')')
      write(nout,'(''/'')')
      write(nout,'(a,''Inf. Dil.'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      do ie=1,nunu
         write(nout,'(1p,2e14.6,''/'')') ee(ie),s0(ie)
      enddo
      write(nout,'(''/'')')
      write(nout,'(''2/'')')
      write(nout,'(''/'')')
      write(nout,'(''0 0 0 1/'')')
      write(nout,'(a,''100 b'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      do ie=1,nunu
         write(nout,'(1p,2e14.6,''/'')') ee(ie),s1(ie)
      enddo
      write(nout,'(''/'')')
      write(nout,'(''3/'')')
      write(nout,'(''/'')')
      write(nout,'(''0 0 0 2/'')')
      write(nout,'(a,''1 b'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      do ie=1,nunu
         write(nout,'(1p,2e14.6,''/'')') ee(ie),s2(ie)
      enddo
      write(nout,'(''/'')')

      !--elastic
      if (lssf.eq.0) then
         xmin=1e10
         xmax=0
         ymin=1e10
         ymax=0
         do ie=1,nure
            ee(ie)=abs(xss(iurpt+5+ie))
            if (ee(ie).gt.xmax) xmax=ee(ie)
            if (ee(ie).lt.xmin) xmin=ee(ie)
            s0(ie)=0
            s1(ie)=0
            s2(ie)=0
            f0=0
            f1=0
            f2=0
            cl=0
            do ib=1,nurb
               ll=iurpt+5+nure+(ie-1)*6*nurb
               dp=xss(ll+ib)-cl
               s0(ie)=s0(ie)+dp*xss(ll+2*nurb+ib)
               s1(ie)=s1(ie)+dp*xss(ll+2*nurb+ib)/(100+xss(ll+nurb+ib))
               s2(ie)=s2(ie)+dp*xss(ll+2*nurb+ib)/(1+xss(ll+nurb+ib))
               f0=f0+dp
               f1=f1+dp/(100+xss(ll+nurb+ib))
               f2=f2+dp/(1+xss(ll+nurb+ib))
               cl=xss(ll+ib)
            enddo
            s1(ie)=s1(ie)/f1
            s2(ie)=s2(ie)/f2
            if (s0(ie).gt.ymax) ymax=s0(ie)
            if (s1(ie).gt.ymax) ymax=s1(ie)
            if (s2(ie).gt.ymax) ymax=s2(ie)
            if (s0(ie).lt.ymin) ymin=s0(ie)
            if (s1(ie).lt.ymin) ymin=s1(ie)
            if (s2(ie).lt.ymin) ymin=s2(ie)
         enddo
         nunu=nure
      else
         xmin=abs(xss(iurpt+5+1))
         xmax=abs(xss(iurpt+5+nure))
         ymin=1e10
         ymax=0
         ie=0
         do i=1,nes
           e=xss(esz-1+i)
           if (e.lt.xmin.or.e.ge.xmax) cycle
           ie=ie+1
           ee(ie)=e
           tot=xss(esz+nes-1+i)
           elas=xss(esz+3*nes-1+i)
           s0(ie)=elas
           do ii=1,nure-1
              e1=abs(xss(iurpt+5+ii))
              e2=abs(xss(iurpt+5+ii+1))
              if (e.lt.e1.or.e.ge.e2) cycle
              kk=iurpt+5+nure+(ii-1)*6*nurb
              ll=iurpt+5+nure+(ii+1-1)*6*nurb
              s1(ie)=0
              s2(ie)=0
              f1=0
              f2=0
              c1=0
              c2=0
              do j=1,nurb
                 call terp1(e1,xss(kk+j)-c1,e2,xss(ll+j)-c2,&
                    e,dp,intunr)
                 call terp1(e1,xss(kk+nurb+j),e2,xss(ll+nurb+j),&
                   e,pp,intunr)
                 call terp1(e1,xss(kk+2*nurb+j),e2,xss(ll+2*nurb+j),&
                   e,pe,intunr)
                 s1(ie)=s1(ie)+dp*pe*elas/(100+pp*tot)
                 s2(ie)=s2(ie)+dp*pe*elas/(1+pp*tot)
                 f1=f1+dp/(100+pp*tot)
                 f2=f2+dp/(1+pp*tot)
                 c1=xss(kk+j)
                 c2=xss(ll+j)
              enddo
              s1(ie)=s1(ie)/f1
              s2(ie)=s2(ie)/f2
              if (s0(ie).gt.ymax) ymax=s0(ie)
              if (s0(ie).lt.ymin) ymin=s0(ie)
              if (s1(ie).gt.ymax) ymax=s1(ie)
              if (s1(ie).lt.ymin) ymin=s1(ie)
              if (s2(ie).gt.ymax) ymax=s2(ie)
              if (s2(ie).lt.ymin) ymin=s2(ie)
           enddo
         enddo
         nunu=ie
      endif
      ymax=ymax+ymax/10
      if (ymin.lt.ymax/1000) ymin=ymax/1000
      call ascll(xmin,xmax)
      call ascll(ymin,ymax)
      write(nout,'(''1'',i3,''/'')') iwcol
      it=1
      do i=1,70
         if (hk(i:i).ne.' ') it=i
      enddo
      write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
      write(nout,'(a,''UR elastic cross section'',a,''/'')') qu,qu
      write(nout,'(''4 0 2 1/'')')
      write(nout,'(1p,3e12.3,''/'')') xmin,xmax,one
      write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
      write(nout,'(1p,3e12.3,''/'')') ymin,ymax,one
      write(nout,'(a,''<c>ross section (barns)'',a,''/'')') qu,qu
      write(nout,'(''/'')')
      write(nout,'(''/'')')
      write(nout,'(a,''Inf. Dil.'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      do ie=1,nunu
         write(nout,'(1p,2e14.6,''/'')') ee(ie),s0(ie)
      enddo
      write(nout,'(''/'')')
      write(nout,'(''2/'')')
      write(nout,'(''/'')')
      write(nout,'(''0 0 0 1/'')')
      write(nout,'(a,''100 b'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      do ie=1,nunu
         write(nout,'(1p,2e14.6,''/'')') ee(ie),s1(ie)
      enddo
      write(nout,'(''/'')')
      write(nout,'(''3/'')')
      write(nout,'(''/'')')
      write(nout,'(''0 0 0 2/'')')
      write(nout,'(a,''1 b'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      do ie=1,nunu
         write(nout,'(1p,2e14.6,''/'')') ee(ie),s2(ie)
      enddo
      write(nout,'(''/'')')

      !--fission
      if (lssf.eq.0) then
         xmin=1e10
         xmax=0
         ymin=1e10
         ymax=0
         do ie=1,nure
            ee(ie)=abs(xss(iurpt+5+ie))
            if (ee(ie).gt.xmax) xmax=ee(ie)
            if (ee(ie).lt.xmin) xmin=ee(ie)
            s0(ie)=0
            s1(ie)=0
            s2(ie)=0
            f0=0
            f1=0
            f2=0
            cl=0
            do ib=1,nurb
               ll=iurpt+5+nure+(ie-1)*6*nurb
               dp=xss(ll+ib)-cl
               s0(ie)=s0(ie)+dp*xss(ll+3*nurb+ib)
               s1(ie)=s1(ie)+dp*xss(ll+3*nurb+ib)/(100+xss(ll+nurb+ib))
               s2(ie)=s2(ie)+dp*xss(ll+3*nurb+ib)/(1+xss(ll+nurb+ib))
               f0=f0+dp
               f1=f1+dp/(100+xss(ll+nurb+ib))
               f2=f2+dp/(1+xss(ll+nurb+ib))
               cl=xss(ll+ib)
            enddo
            s1(ie)=s1(ie)/f1
            s2(ie)=s2(ie)/f2
            if (s0(ie).gt.ymax) ymax=s0(ie)
            if (s1(ie).gt.ymax) ymax=s1(ie)
            if (s2(ie).gt.ymax) ymax=s2(ie)
            if (s0(ie).lt.ymin) ymin=s0(ie)
            if (s1(ie).lt.ymin) ymin=s1(ie)
            if (s2(ie).lt.ymin) ymin=s2(ie)
         enddo
         nunu=nure
      else
         xmin=abs(xss(iurpt+5+1))
         xmax=abs(xss(iurpt+5+nure))
         ymin=1e10
         ymax=0
         nnf=0
         do i=1,ntr
            mt=nint(xss(mtr-1+i))
            if (mt.eq.18.or.mt.eq.19) then
               kf=nint(xss(lsig-1+i)+sig-1)
               nnf=nint(xss(kf+1))
               iif=nint(xss(kf))
            endif
         enddo
         ie=0
         do i=1,nes
            e=xss(esz-1+i)
            if (e.lt.xmin.or.e.ge.xmax) cycle
            ie=ie+1
            ee(ie)=e
            tot=xss(esz+nes-1+i)
            fiss=0
            if (nnf.gt.0) then
               if (i.ge.iif) fiss=xss(kf+2+i-iif)
            endif
            s0(ie)=fiss
            do ii=1,nure-1
               e1=abs(xss(iurpt+5+ii))
               e2=abs(xss(iurpt+5+ii+1))
               if (e.lt.e1.or.e.ge.e2) cycle
               kk=iurpt+5+nure+(ii-1)*6*nurb
               ll=iurpt+5+nure+(ii+1-1)*6*nurb
               s1(ie)=0
               s2(ie)=0
               f1=0
               f2=0
               c1=0
               c2=0
               do j=1,nurb
                  call terp1(e1,xss(kk+j)-c1,e2,xss(ll+j)-c2,&
                     e,dp,intunr)
                  call terp1(e1,xss(kk+nurb+j),e2,xss(ll+nurb+j),&
                    e,pp,intunr)
                  call terp1(e1,xss(kk+3*nurb+j),e2,xss(ll+3*nurb+j),&
                    e,pe,intunr)
                  s1(ie)=s1(ie)+dp*pe*fiss/(100+pp*tot)
                  s2(ie)=s2(ie)+dp*pe*fiss/(1+pp*tot)
                  f1=f1+dp/(100+pp*tot)
                  f2=f2+dp/(1+pp*tot)
                  c1=xss(kk+j)
                  c2=xss(ll+j)
               enddo
               s1(ie)=s1(ie)/f1
               s2(ie)=s2(ie)/f2
               if (s0(ie).gt.ymax) ymax=s0(ie)
               if (s0(ie).lt.ymin) ymin=s0(ie)
               if (s1(ie).gt.ymax) ymax=s1(ie)
               if (s1(ie).lt.ymin) ymin=s1(ie)
               if (s2(ie).gt.ymax) ymax=s2(ie)
               if (s2(ie).lt.ymin) ymin=s2(ie)
            enddo
         enddo
         nunu=ie
      endif
      if (ymax.gt.zero) then
        ymax=ymax+ymax/10
        if (ymin.lt.ymax/1000) ymin=ymax/1000
        call ascll(xmin,xmax)
        call ascll(ymin,ymax)
        write(nout,'(''1'',i3,''/'')') iwcol
         it=1
         do i=1,70
            if (hk(i:i).ne.' ') it=i
         enddo
         write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
         write(nout,'(a,''UR fission cross section'',a,''/'')') qu,qu
         write(nout,'(''4 0 2 1/'')')
         write(nout,'(1p,3e12.3,''/'')') xmin,xmax,one
         write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
         write(nout,'(1p,3e12.3,''/'')') ymin,ymax,one
         write(nout,'(a,''<c>ross section (barns)'',a,''/'')') qu,qu
         write(nout,'(''/'')')
         write(nout,'(''/'')')
         write(nout,'(a,''Inf. Dil.'',a,''/'')') qu,qu
         write(nout,'(''0/'')')
         do ie=1,nunu
            write(nout,'(1p,2e14.6,''/'')') ee(ie),s0(ie)
         enddo
         write(nout,'(''/'')')
         write(nout,'(''2/'')')
         write(nout,'(''/'')')
         write(nout,'(''0 0 0 1/'')')
         write(nout,'(a,''100 b'',a,''/'')') qu,qu
         write(nout,'(''0/'')')
         do ie=1,nunu
            write(nout,'(1p,2e14.6,''/'')') ee(ie),s1(ie)
         enddo
         write(nout,'(''/'')')
         write(nout,'(''3/'')')
         write(nout,'(''/'')')
         write(nout,'(''0 0 0 2/'')')
         write(nout,'(a,''1 b'',a,''/'')') qu,qu
         write(nout,'(''0/'')')
         do ie=1,nunu
            write(nout,'(1p,2e14.6,''/'')') ee(ie),s2(ie)
         enddo
         write(nout,'(''/'')')
      endif

      !--capture
      if (lssf.eq.0) then
         xmin=1e10
         xmax=0
         ymin=1e10
         ymax=0
         do ie=1,nure
            ee(ie)=abs(xss(iurpt+5+ie))
            if (ee(ie).gt.xmax) xmax=ee(ie)
            if (ee(ie).lt.xmin) xmin=ee(ie)
            s0(ie)=0
            s1(ie)=0
            s2(ie)=0
            f0=0
            f1=0
            f2=0
            cl=0
            do ib=1,nurb
               ll=iurpt+5+nure+(ie-1)*6*nurb
               dp=xss(ll+ib)-cl
               s0(ie)=s0(ie)+dp*xss(ll+4*nurb+ib)
               s1(ie)=s1(ie)+dp*xss(ll+4*nurb+ib)/(100+xss(ll+nurb+ib))
               s2(ie)=s2(ie)+dp*xss(ll+4*nurb+ib)/(1+xss(ll+nurb+ib))
               f0=f0+dp
               f1=f1+dp/(100+xss(ll+nurb+ib))
               f2=f2+dp/(1+xss(ll+nurb+ib))
               cl=xss(ll+ib)
            enddo
            s1(ie)=s1(ie)/f1
            s2(ie)=s2(ie)/f2
            if (s0(ie).gt.ymax) ymax=s0(ie)
            if (s1(ie).gt.ymax) ymax=s1(ie)
            if (s2(ie).gt.ymax) ymax=s2(ie)
            if (s0(ie).lt.ymin) ymin=s0(ie)
            if (s1(ie).lt.ymin) ymin=s1(ie)
            if (s2(ie).lt.ymin) ymin=s2(ie)
         enddo
      else
         xmin=abs(xss(iurpt+5+1))
         xmax=abs(xss(iurpt+5+nure))
         ymin=1e10
         ymax=0
         nnf=0
         do i=1,ntr
            mt=nint(xss(mtr-1+i))
            if (mt.eq.102) then
               kc=nint(xss(lsig-1+i)+sig-1)
               iic=nint(xss(kc))
            endif
         enddo
         ie=0
         do i=1,nes
            e=xss(esz-1+i)
            if (e.lt.xmin.or.e.ge.xmax) cycle
            ie=ie+1
            ee(ie)=e
            tot=xss(esz+nes-1+i)
            capt=xss(kc+2+i-iic)
            s0(ie)=capt
            do ii=1,nure-1
               e1=abs(xss(iurpt+5+ii))
               e2=abs(xss(iurpt+5+ii+1))
               if (e.lt.e1.or.e.ge.e2) cycle
               kk=iurpt+5+nure+(ii-1)*6*nurb
               ll=iurpt+5+nure+(ii+1-1)*6*nurb
               s1(ie)=0
               s2(ie)=0
               f1=0
               f2=0
               c1=0
               c2=0
               do j=1,nurb
                  call terp1(e1,xss(kk+j)-c1,e2,xss(ll+j)-c2,&
                     e,dp,intunr)
                  call terp1(e1,xss(kk+nurb+j),e2,xss(ll+nurb+j),&
                    e,pp,intunr)
                  call terp1(e1,xss(kk+4*nurb+j),e2,xss(ll+4*nurb+j),&
                    e,pe,intunr)
                  s1(ie)=s1(ie)+dp*pe*capt/(100+pp*tot)
                  s2(ie)=s2(ie)+dp*pe*capt/(1+pp*tot)
                  f1=f1+dp/(100+pp*tot)
                  f2=f2+dp/(1+pp*tot)
                  c1=xss(kk+j)
                  c2=xss(ll+j)
               enddo
               s1(ie)=s1(ie)/f1
               s2(ie)=s2(ie)/f2
               if (s0(ie).gt.ymax) ymax=s0(ie)
               if (s0(ie).lt.ymin) ymin=s0(ie)
               if (s1(ie).gt.ymax) ymax=s1(ie)
               if (s1(ie).lt.ymin) ymin=s1(ie)
               if (s2(ie).gt.ymax) ymax=s2(ie)
               if (s2(ie).lt.ymin) ymin=s2(ie)
            enddo
         enddo
         nunu=ie
      endif
      ymax=ymax+ymax/10
      if (ymin.lt.ymax/1000) ymin=ymax/1000
      call ascll(xmin,xmax)
      call ascll(ymin,ymax)
      write(nout,'(''1'',i3,''/'')') iwcol
      it=1
      do i=1,70
         if (hk(i:i).ne.' ') it=i
      enddo
      write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
      write(nout,'(a,''UR capture cross section'',a,''/'')') qu,qu
      write(nout,'(''4 0 2 1/'')')
      write(nout,'(1p,3e12.3,''/'')') xmin,xmax,one
      write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
      write(nout,'(1p,3e12.3,''/'')') ymin,ymax,one
      write(nout,'(a,''<c>ross section (barns)'',a,''/'')') qu,qu
      write(nout,'(''/'')')
      write(nout,'(''/'')')
      write(nout,'(a,''Inf. Dil.'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      do ie=1,nunu
         write(nout,'(1p,2e14.6,''/'')') ee(ie),s0(ie)
      enddo
      write(nout,'(''/'')')
      write(nout,'(''2/'')')
      write(nout,'(''/'')')
      write(nout,'(''0 0 0 1/'')')
      write(nout,'(a,''100 b'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      do ie=1,nunu
         write(nout,'(1p,2e14.6,''/'')') ee(ie),s1(ie)
      enddo
      write(nout,'(''/'')')
      write(nout,'(''3/'')')
      write(nout,'(''/'')')
      write(nout,'(''0 0 0 2/'')')
      write(nout,'(a,''1 b'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      do ie=1,nunu
         write(nout,'(1p,2e14.6,''/'')') ee(ie),s2(ie)
      enddo
      write(nout,'(''/'')')
   endif

   !--plot log-log heating per reaction
   xmin=big
   xmax=0
   ymin=big
   ymax=-big
   do i=1,nes
      e=xss(esz-1+i)
      heat=xss(esz+4*nes-1+i)
      if (heat.lt.hmin) heat=hmin
      if (e.lt.xmin) xmin=e
      if (e.gt.xmax) xmax=e
      if (heat.lt.ymin) ymin=heat
      if (heat.gt.ymax) ymax=heat
   enddo
   if (ymax.gt.ymin) then
      call ascll(xmin,xmax)
      if (ymin.lt.ymax/scale) ymin=ymax/scale
      call ascll(ymin,ymax)
      write(nout,'(''1'',i3,''/'')') iwcol
      it=1
      do i=1,70
         if (hk(i:i).ne.' ') it=i
      enddo
      write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
      write(nout,'(a,''<h>eating'',a,''/'')') qu,qu
      write(nout,'(''4 0 2 1/'')')
      write(nout,'(1p,3e12.3,''/'')') xmin,xmax,one
      write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
      write(nout,'(1p,3e12.3,''/'')') ymin,ymax,one
      write(nout,'(a,''<h>eating (<m>e<v>/reaction)'',a,''/'')') qu,qu
      write(nout,'(''/'')')
      write(nout,'(''/'')')
      write(nout,'(a,''heating'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      thin=ten**(log10(xmax/xmin)/nden)
      xlast=small
      j=0
      do i=1,nes
         e=xss(esz-1+i)
         if (nes.le.nden.or.e.ge.thin*xlast) then
            heat=xss(esz+4*nes-1+i)
            if (heat.lt.hmin) heat=hmin
            if (heat.lt.ymin) heat=ymin
            j=j+1
            write(nout,'(1p,2e14.6,''/'')') e,heat
            xlast=e
         endif
      enddo
      write(nout,'(''/'')')
   endif

   !--plot log-log damage
   xmin=big
   xmax=0
   ymin=big
   ymax=-big
   n=0
   do i=1,ntr
      mt=nint(xss(mtr+i-1))
      if (mt.eq.444) then
         k=nint(xss(lsig+i-1)+sig-1)
         n=nint(xss(k+1))
         iaa=nint(xss(k))
      endif
   enddo
   if (n.ne.0) then
      do i=1,n
         e=xss(iaa+i-1)
         dam=xss(k+2+i-1)
         if (dam.lt.hmin) dam=hmin
         if (e.lt.xmin) xmin=e
         if (e.gt.xmax) xmax=e
         if (dam.lt.ymin) ymin=dam
         if (dam.gt.ymax) ymax=dam
      enddo
      if (ymax.gt.ymin) then
         call ascll(xmin,xmax)
         if (ymin.lt.ymax/scale) ymin=ymax/scale
         call ascll(ymin,ymax)
         write(nout,'(''1'',i3,''/'')') iwcol
         it=1
         do i=1,70
            if (hk(i:i).ne.' ') it=i
         enddo
         write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
         write(nout,'(a,''<d>amage'',a,''/'')') qu,qu
         write(nout,'(''4 0 2 1/'')')
         write(nout,'(1p,3e12.3,''/'')') xmin,xmax,one
         write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
         write(nout,'(1p,3e12.3,''/'')') ymin,ymax,one
         write(nout,'(a,''<d>amage (<m>e<v>-barns)'',a,''/'')') qu,qu
         write(nout,'(''/'')')
         write(nout,'(''/'')')
         write(nout,'(a,''damage'',a,''/'')') qu,qu
         write(nout,'(''0/'')')
         thin=ten**(log10(xmax/xmin)/nden)
         xlast=small
         j=0
         do i=1,n
            e=xss(iaa+i-1)
            if (nes.le.nden.or.e.ge.thin*xlast) then
               dam=xss(k+2+i-1)
               if (dam.lt.hmin) dam=hmin
               if (dam.lt.ymin) dam=ymin
               j=j+1
               write(nout,'(1p,2e14.6,''/'')') e,dam
               xlast=e
            endif
         enddo
         write(nout,'(''/'')')
      endif
   endif

   !--make pages showing the nonthreshold reactions
   mtlast=0

   !--plot log-log non-threshold reactions
   nlev=1
   do while (nlev.gt.0)
      xmin=1000
      xmax=0
      ymin=1000
      ymax=0
      nlev=0
      do i=1,ntr
         mt=nint(xss(mtr+i-1))
         k=nint(xss(lsig+i-1)+sig-1)
         n=nint(xss(k+1))
         iaa=nint(xss(k))
         iflag=0
         if (nlev.eq.5) iflag=1
         if (mt.le.mtlast) iflag=1
         if (mt.gt.207) iflag=1
         if (xss(iaa).gt.1.e-6) iflag=1
         if (iflag.eq.0) then
            nlev=nlev+1
            do j=1,n
               x=xss(iaa+j-1)
               y=xss(k+2+j-1)
               if (y.ne.zero.or.j.le.1) then
                  if (x.lt.xmin) xmin=x
                  if (x.gt.xmax) xmax=x
                  if (y.lt.ymin) ymin=y
                  if (y.gt.ymax) ymax=y
               endif
            enddo
         endif
      enddo
      if (nlev.ne.0) then
         call ascll(xmin,xmax)
         if (ymin.lt.ymax/scale) ymin=ymax/scale
         call ascll(ymin,ymax)
         thin=ten**(log10(xmax/xmin)/nden)
         write(nout,'(''1'',i3,''/'')') iwcol
         it=1
         do i=1,70
            if (hk(i:i).ne.' ') it=i
         enddo
         write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
         write(nout,'(a,''<n>on-threshold reactions'',a,''/'')') qu,qu
         write(nout,'(''4 0 2 1/'')')
         write(nout,'(1p,3e12.3,''/'')') xmin,xmax,one
         write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
         write(nout,'(1p,3e12.3,''/'')') ymin,ymax,one
         write(nout,'(a,''<c>ross section (barns)'',a,''/'')') qu,qu
         write(nout,'(''/'')')
         nlev=0
         do i=1,ntr
            mt=nint(xss(mtr+i-1))
            k=nint(xss(lsig+i-1)+sig-1)
            n=nint(xss(k+1))
            iaa=nint(xss(k))
            iflag=0
            if (nlev.eq.5) iflag=1
            if (mt.le.mtlast) iflag=1
            if (mt.gt.207) iflag=1
            if (xss(iaa).gt.1.e-6) iflag=1
            if (iflag.eq.0) then
               mtl=mt
               nlev=nlev+1
               if (nlev.gt.1) write(nout,'(i2,''/'')') nlev
               if (nlev.gt.1) write(nout,'(''/'')')
               icurv=mod(nlev-1,5)
               if (iwcol.eq.0) then
                  write(nout,'(''0 0'',i2,''/'')') icurv
               else
                  write(nout,'(''0 0 0'',i2,''/'')') icurv
               endif
               call mtname(mt,name,izai)
               if (name(1:1).eq.'(') then
                  if (izai.eq.1) then
                     name(2:2)='n'
                  else if (izai.eq.1001) then
                     name(2:2)='p'
                  else if (izai.eq.1002) then
                     name(2:2)='d'
                  else if (izai.eq.1003) then
                     name(2:2)='t'
                  else if (izai.eq.2003) then
                     name(2:2)='s'
                  else if (izai.eq.2004) then
                     name(2:2)='a'
                  endif
               endif
               write(nout,'(a,a,a,''/'')') qu,name,qu
               write(nout,'(''0/'')')
               xlast=small
               do j=1,n
                  x=xss(iaa+j-1)
                  y=xss(k+2+j-1)
                  if (y.lt.ymin) y=ymin
                  if (n.le.nden.or.x.ge.thin*xlast) then
                     write(nout,'(1p,2e14.6,''/'')') x,y
                     xlast=x
                  endif
               enddo
               write(nout,'(''/'')')
            endif
         enddo
         mtlast=mtl
      endif
   enddo

   !--plot lin-lin total, absorption, elastic, and gamma production
   xmin=big
   xmax=0
   ymin=big
   ymax=-big
   do i=1,nes
      e=xss(esz-1+i)
      if (e.gt.2) then
         tot=xss(esz+nes-1+i)
         abso=xss(esz+2*nes-1+i)
         elas=xss(esz+3*nes-1+i)
         if (e.lt.xmin) xmin=e
         if (e.gt.xmax) xmax=e
         if (tot.lt.ymin) ymin=tot
         if (tot.gt.ymax) ymax=tot
         if (abso.lt.ymin) ymin=abso
         if (abso.gt.ymax) ymax=abso
         if (elas.lt.ymin) ymin=elas
         if (elas.gt.ymax) ymax=elas
         if (gpd.ne.0) then
            gprod=xss(gpd-1+i)
            if (gprod.lt.ymin) ymin=gprod
            if (gprod.gt.ymax) ymax=gprod
         endif
      endif
   enddo
   ymin=0
   call ascle(4,xmin,xmax,major,minor)
   xstep=(xmax-xmin)/major
   call ascle(4,ymin,ymax,major,minor)
   ystep=(ymax-ymin)/major
   write(nout,'(''1'',i3,''/'')') iwcol
   it=1
   do i=1,70
      if (hk(i:i).ne.' ') it=i
   enddo
   write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
   write(nout,'(a,''<p>rincipal cross sections'',a,''/'')') qu,qu
   xtag=35*xmin/100+65*xmax/100
   ytag=ymin/10+9*ymax/10
   write(nout,'(''1 0 2 1'',2e12.4,''/'')') xtag,ytag
   write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
   write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
   write(nout,'(1p,3e12.3,''/'')') ymin,ymax,ystep
   write(nout,'(a,''<c>ross section (barns)'',a,''/'')') qu,qu
   write(nout,'(''/'')')
   write(nout,'(''/'')')
   write(nout,'(a,''total'',a,''/'')') qu,qu
   write(nout,'(''0/'')')
   thin=(xmax-xmin)/nden
   xlast=small
   j=0
   do i=1,nes
      e=xss(esz-1+i)
      test=1
      test=test/5
      if (e.ge.test) then
         if (nes.le.nden.or.i.eq.nes.or.e.ge.xlast+thin) then
            tot=xss(esz+nes-1+i)
            j=j+1
            write(nout,'(1p,2e14.6,''/'')') e,tot
            xlast=e
         endif
      endif
   enddo
   write(nout,'(''/'')')
   write(nout,'(''2/'')')
   write(nout,'(''/'')')
   if (iwcol.eq.0) then
      write(nout,'(''0 0 1/'')')
   else
      write(nout,'(''0 0 0 3/'')')
   endif
   write(nout,'(a,''absorption'',a,''/'')') qu,qu
   write(nout,'(''0/'')')
   xlast=small
   j=0
   do i=1,nes
      e=xss(esz-1+i)
      test=1
      test=test/10
      if (e.ge.test) then
         if (nes.le.nden.or.e.ge.xlast+thin.or.i.eq.nes) then
            abss=xss(esz+2*nes-1+i)
            j=j+1
            write(nout,'(1p,2e14.6,''/'')') e,abss
            xlast=e
         endif
      endif
   enddo
   write(nout,'(''/'')')
   write(nout,'(''3/'')')
   write(nout,'(''/'')')
   if (iwcol.eq.0) then
      write(nout,'(''0 0 2/'')')
   else
      write(nout,'(''0 0 0 2/'')')
   endif
   write(nout,'(a,''elastic'',a,''/'')') qu,qu
   write(nout,'(''0/'')')
   xlast=small
   j=0
   do i=1,nes
      e=xss(esz-1+i)
      test=1
      test=test/10
      if (e.ge.test) then
         if (nes.le.nden.or.e.ge.xlast+thin.or.i.eq.nes) then
            elas=xss(esz+3*nes-1+i)
            j=j+1
            write(nout,'(1p,2e14.6,''/'')') e,elas
            xlast=e
         endif
      endif
   enddo
   write(nout,'(''/'')')
   if (gpd.ne.0) then
      write(nout,'(''4/'')')
      write(nout,'(''/'')')
      if (iwcol.eq.0) then
         write(nout,'(''0 0 3/'')')
      else
         write(nout,'(''0 0 0 1/'')')
      endif
      write(nout,'(a,''gamma production'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      xlast=small
      j=0
      do i=1,nes
         e=xss(esz-1+i)
         test=1
         test=test/10
         if (e.ge.test) then
            if (nes.le.nden.or.e.ge.xlast+thin) then
               gprod=xss(gpd-1+i)
               j=j+1
               write(nout,'(1p,2e14.6,''/'')') e,gprod
               xlast=e
            endif
         endif
      enddo
      write(nout,'(''/'')')
   endif

   !--plot lin-lin heating per reaction
   xmin=big
   xmax=0
   ymin=big
   ymax=-big
   do i=1,nes
      e=xss(esz-1+i)
      heat=xss(esz+4*nes-1+i)
      if (e.lt.xmin) xmin=e
      if (e.gt.xmax) xmax=e
      if (heat.lt.ymin) ymin=heat
      if (heat.gt.ymax.and.e.gt.one) ymax=heat
   enddo
   if (ymin.ne.zero.or.ymax.ne.zero) then
      call ascle(4,xmin,xmax,major,minor)
      xstep=(xmax-xmin)/major
      call ascle(4,ymin,ymax,major,minor)
      ystep=(ymax-ymin)/major
      write(nout,'(''1'',i3,''/'')') iwcol
      it=1
      do i=1,70
         if (hk(i:i).ne.' ') it=i
      enddo
      write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
      write(nout,'(a,''<h>eating'',a,''/'')') qu,qu
      write(nout,'(''1 0 2 1/'')')
      write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
      write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
      write(nout,'(1p,3e12.3,''/'')') ymin,ymax,ystep
      write(nout,'(a,''<h>eating (<m>e<v>/reaction)'',a,''/'')') qu,qu
      write(nout,'(''/'')')
      write(nout,'(''/'')')
      write(nout,'(a,''heating'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      thin=(xmax-xmin)/nden
      xlast=small
      j=0
      do i=1,nes
         e=xss(esz-1+i)
         test=1
         test=test/5
         if (e.ge.test) then
            if (nes.le.nden.or.e.ge.xlast+thin.or.i.eq.nes) then
               heat=xss(esz+4*nes-1+i)
               j=j+1
               write(nout,'(1p,2e14.6,''/'')') e,heat
               xlast=e
            endif
         endif
      enddo
      write(nout,'(''/'')')
   endif

   !--plot lin-lin damage
   xmin=big
   xmax=0
   ymin=big
   ymax=-big
   n=0
   do i=1,ntr
      mt=nint(xss(mtr+i-1))
      if (mt.eq.444) then
         k=nint(xss(lsig+i-1)+sig-1)
         n=nint(xss(k+1))
         iaa=nint(xss(k))
      endif
   enddo
   if (n.ne.0) then
      do i=1,n
         e=xss(iaa+i-1)
         test=1
         test=test/5
         if (e.ge.test) then
            dam=xss(k+2+i-1)
            if (dam.lt.hmin) dam=hmin
            if (e.lt.xmin) xmin=e
            if (e.gt.xmax) xmax=e
            if (dam.lt.ymin) ymin=dam
            if (dam.gt.ymax.and.e.gt.1.) ymax=dam
         endif
      enddo
      if (ymin.ne.zero.or.ymax.ne.zero) then
         call ascle(4,xmin,xmax,major,minor)
         xstep=(xmax-xmin)/major
         call ascle(4,ymin,ymax,major,minor)
         ystep=(ymax-ymin)/major
         write(nout,'(''1'',i3,''/'')') iwcol
         it=1
         do i=1,70
            if (hk(i:i).ne.' ') it=i
         enddo
         write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
         write(nout,'(a,''<d>amage'',a,''/'')') qu,qu
         write(nout,'(''1 0 2 1/'')')
         write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
         write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
         write(nout,'(1p,3e12.3,''/'')') ymin,ymax,ystep
         write(nout,'(a,''<d>amage (<m>e<v>-barns)'',a,''/'')') qu,qu
         write(nout,'(''/'')')
         write(nout,'(''/'')')
         write(nout,'(a,''damage'',a,''/'')') qu,qu
         write(nout,'(''0/'')')
         thin=(xmax-xmin)/nden
         xlast=small
         j=0
         do i=1,n
            e=xss(iaa+i-1)
            test=1
            test=test/5
            if (e.ge.test) then
               if (nes.le.nden.or.e.ge.xlast+thin.or.i.eq.nes) then
                  dam=xss(k+2+i-1)
                  if (dam.lt.hmin) dam=hmin
                  j=j+1
                  write(nout,'(1p,2e14.6,''/'')') e,dam
                  xlast=e
               endif
            endif
         enddo
         write(nout,'(''/'')')
      endif
   endif

   !--make pages showing the nonthreshold reactions
   mtlast=0

   !--plot lin-log non-threshold reactions
   idone=0
   do while (idone.eq.0)
      xmin=1000
      xmax=0
      ymin=1000
      ymax=0
      nlev=0
      do i=1,ntr
         mt=nint(xss(mtr+i-1))
         k=nint(xss(lsig+i-1)+sig-1)
         n=nint(xss(k+1))
         iaa=nint(xss(k))
         iflag=0
         if (nlev.eq.5) iflag=1
         if (mt.le.mtlast) iflag=1
         if (mt.gt.207) iflag=1
         if (xss(iaa).gt.1.e-6) iflag=1
         if (iflag.eq.0) then
            nlev=nlev+1
            do j=1,n
               x=xss(iaa+j-1)
               y=xss(k+2+j-1)
               test=1
               test=test/5
               iflag=0
               if (x.lt.test) iflag=1
               if (y.eq.zero.and.j.gt.1) iflag=1
               if (iflag.eq.0) then
                  if (x.lt.xmin) xmin=x
                  if (x.gt.xmax) xmax=x
                  if (y.lt.ymin) ymin=y
                  if (y.gt.ymax.and.x.gt.one) ymax=y
               endif
            enddo
         endif
      enddo
      if (nlev.eq.0.or.ymax.eq.zero) then
         idone=1
      else
         call ascle(4,xmin,xmax,major,minor)
         xstep=(xmax-xmin)/major
         call ascll(ymin,ymax)
         if (ymax.eq.ymin) then
            idone=1
         else
            thin=(xmax-xmin)/nden
            write(nout,'(''1'',i3,''/'')') iwcol
            it=1
            do i=1,70
               if (hk(i:i).ne.' ') it=i
            enddo
            write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
            write(nout,&
              '(a,''<n>on-threshold reactions'',a,''/'')') qu,qu
            write(nout,'(''2 0 2 1/'')')
            write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
            write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
            write(nout,'(1p,3e12.3,''/'')') ymin,ymax,one
            write(nout,&
              '(a,''<c>ross section (barns)'',a,''/'')') qu,qu
            write(nout,'(''/'')')
            nlev=0
            do i=1,ntr
               mt=nint(xss(mtr+i-1))
               k=nint(xss(lsig+i-1)+sig-1)
               n=nint(xss(k+1))
               iaa=nint(xss(k))
               iflag=0
               if (nlev.eq.5) iflag=1
               if (mt.le.mtlast) iflag=1
               if (mt.gt.207) iflag=1
               if (xss(iaa).gt.1.e-6) iflag=1
               if (iflag.eq.0) then
                  mtl=mt
                  nlev=nlev+1
                  if (nlev.gt.1) write(nout,'(i2,''/'')') nlev
                  if (nlev.gt.1) write(nout,'(''/'')')
                  icurv=mod(nlev-1,5)
                  if (iwcol.eq.0) then
                     write(nout,'(''0 0'',i2,''/'')') icurv
                  else
                     write(nout,'(''0 0 0'',i2,''/'')') icurv
                  endif
                  call mtname(mt,name,izai)
                  if (name(1:1).eq.'(') then
                     if (izai.eq.1) then
                        name(2:2)='n'
                     else if (izai.eq.1001) then
                        name(2:2)='p'
                     else if (izai.eq.1002) then
                        name(2:2)='d'
                     else if (izai.eq.1003) then
                        name(2:2)='t'
                     else if (izai.eq.2003) then
                        name(2:2)='s'
                     else if (izai.eq.2004) then
                        name(2:2)='a'
                     endif
                  endif
                  write(nout,'(a,a,a,''/'')') qu,name,qu
                  write(nout,'(''0/'')')
                  xlast=small
                  do j=1,n
                     x=xss(iaa+j-1)
                     test=1
                     test=test/5
                     if (x.ge.test) then
                        if (n.le.nden.or.x.ge.xlast+thin.or.j.eq.n)&
                          then
                           y=xss(k+2+j-1)
                           if (y.lt.ymin) y=ymin
                           write(nout,'(1p,2e14.6,''/'')') x,y
                           xlast=x
                        endif
                     endif
                  enddo
                  write(nout,'(''/'')')
               endif
            enddo
            mtlast=mtl
         endif
      endif
   enddo

   !--make the pages showing the inelastic level cross sections
   !--and threshold reactions with up to 5 levels per page
   call aplotr(nout,iwcol,hk)

   !--make the page showing the higher fission reactions
   call aplopf(nout,iwcol,hk)

   !--make contour plots for angular distributions
   call aplof4(nout,iwcol,hk)

   !--plot nubar
   if (nu.gt.0) call aplonu(nout,iwcol,hk)

   !--make 3d pictures for energy distributions
   if (nr.ne.0) call aplodd(nout,iwcol,hk)

   !--plot delayed-neutron data
   if (ndnf.gt.0) call aplodn(nout,iwcol,hk)

   !--plot detailed photon data
   call aplopp(nout,iwcol,hk)

   !--plot particle production sections
   if (ntype.gt.0) call aploxp(nout,iwcol,hk)

   !--end the plotr input text
   write(nout,'(''99/'')')
   write(nout,'(''stop'')')
   return
   end subroutine aplots

   subroutine aplotr(nout,iwcol,hk)
   !-------------------------------------------------------------------
   ! Plot the inelastic levels and threshold reactions
   ! showing up to 5 curves per page.
   !-------------------------------------------------------------------
   use acecm ! provides mtname
   ! externals
   integer::nout,iwcol
   character(70)::hk
   ! internals
   integer::mtlast,ilev,ilast,i,mt,k,n,iaa,j,major,minor,mtl
   integer::iflag,icurv,idone,it,nlev
   real(kr)::ylast,x,y,xstep,ystep,thin,xtag,ytag
   real(kr)::xmin,xmax,ymin,ymax,xlast
   character(10)::name
   character(1)::qu=''''
   real(kr),parameter::small=1.e-12_kr
   real(kr),parameter::xsmin=1.e-6_kr
   real(kr),parameter::zero=0
   integer,parameter::nden=4000
   mtlast=0
   ilev=0

   !--loop over the inelastic levels
   !--in pages with 5 curves per page
   ilast=0
   do while (ilast.eq.0)
      xmin=1000
      xmax=0
      ymin=1000
      ymax=0
      nlev=0
      do i=1,ntr
         mt=nint(xss(mtr+i-1))
         if (mt.ge.51+ilev.and.mt.le.90.and.izai.eq.1) then
            if (nlev.lt.5) then
               nlev=nlev+1
               k=nint(xss(lsig+i-1)+sig-1)
               n=nint(xss(k+1))
               iaa=nint(xss(k))
               ylast=-1
               do j=1,n
                  x=xss(iaa+j-1)
                  y=xss(k+2+j-1)
                  if (y.ne.ylast) then
                     if (x.lt.xmin) xmin=x
                     if (x.gt.xmax) xmax=x
                     if (y.lt.ymin) ymin=y
                     if (y.gt.ymax) ymax=y
                  endif
                  ylast=y
               enddo
            endif
         endif
      enddo
      if (nlev.eq.0) then
         ilast=1
      else
         call ascle(4,xmin,xmax,major,minor)
         xstep=(xmax-xmin)/major
         call ascle(4,ymin,ymax,major,minor)
         ystep=(ymax-ymin)/major
         thin=(xmax-xmin)/nden
         write(nout,'(''1'',i3,''/'')') iwcol
         it=1
         do i=1,70
            if (hk(i:i).ne.' ') it=i
         enddo
         write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
         write(nout,'(a,''<i>nelastic levels'',a,''/'')') qu,qu
         xtag=3*xmin/10+7*xmax/10
         ytag=ymin/10+9*ymax/10
         write(nout,'(''1 0 2 1'',2e12.4,''/'')') xtag,ytag
         write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
         write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
         write(nout,'(1p,3e12.3,''/'')') ymin,ymax,ystep
         write(nout,'(a,''<c>ross section (barns)'',a,''/'')') qu,qu
         write(nout,'(''/'')')
         nlev=0
         idone=0
         i=0
         do while (i.lt.ntr.and.idone.eq.0)
            i=i+1
            mt=nint(xss(mtr+i-1))
            if (mt.ge.51+ilev) then
               if (mt.gt.90.or.nlev.ge.5) then
                  idone=1
               else
                  nlev=nlev+1
                  k=nint(xss(lsig+i-1)+sig-1)
                  n=nint(xss(k+1))
                  iaa=nint(xss(k))
                  if (nlev.gt.1) write(nout,'(i2,''/'')') nlev
                  if (nlev.gt.1) write(nout,'(''/'')')
                  icurv=mod(nlev-1,5)
                  if (iwcol.eq.0) then
                     write(nout,'(''0 0'',i2,''/'')') icurv
                  else
                     write(nout,'(''0 0 0'',i2,''/'')') icurv
                  endif
                  call mtname(mt,name,izai)
                  if (name(1:1).eq.'(') then
                     if (izai.eq.1) then
                        name(2:2)='n'
                     else if (izai.eq.1001) then
                        name(2:2)='p'
                     else if (izai.eq.1002) then
                        name(2:2)='d'
                     else if (izai.eq.1003) then
                        name(2:2)='t'
                     else if (izai.eq.2003) then
                        name(2:2)='s'
                     else if (izai.eq.2004) then
                        name(2:2)='a'
                     endif
                  endif
                  write(nout,'(a,a,a,''/'')') qu,name,qu
                  write(nout,'(''0/'')')
                  xlast=small
                  ylast=-1
                  do j=1,n
                     x=xss(iaa+j-1)
                     if (n.le.nden.or.x.ge.xlast+thin.or.j.eq.n) then
                        y=xss(k+2+j-1)
                        write(nout,'(1p,2e14.6,''/'')') x,y
                        xlast=x
                        ylast=y
                     endif
                  enddo
                  write(nout,'(''/'')')
               endif
            endif
         enddo
      endif
      ilev=ilev+nlev
   enddo

   !--loop over the threshold reactions
   idone=0
   mtlast=0
   do while (idone.eq.0)
      xmin=1000
      xmax=0
      ymin=1000
      ymax=0
      nlev=0
      do i=1,ntr
         mt=nint(xss(mtr+i-1))
         k=nint(xss(lsig+i-1)+sig-1)
         n=nint(xss(k+1))
         iaa=nint(xss(k))
         iflag=0
         if (nlev.eq.5) iflag=1
         if (mt.lt.5) iflag=1
         if (mt.ge.18.and.mt.le.21) iflag=1
         if (mt.eq.38) iflag=1
         if (mt.ge.50.and.mt.lt.91.and.izai.eq.1) iflag=1
         if (mt.gt.207.and.mt.lt.600) iflag=1
         if (mt.le.mtlast) iflag=1
         if (xss(iaa).lt.xsmin) iflag=1
         if (iflag.eq.0) then
            nlev=nlev+1
            do j=1,n
               x=xss(iaa+j-1)
               y=xss(k+2+j-1)
               if (y.ne.zero.or.j.le.1) then
                  if (x.lt.xmin) xmin=x
                  if (x.gt.xmax) xmax=x
                  if (y.lt.ymin) ymin=y
                  if (y.gt.ymax) ymax=y
               endif
            enddo
            if (mt.lt.203.and.nint(xss(mtr+i)).ge.203.and.&
              nint(xss(mtr+i)).le.207) nlev=5
         endif
      enddo
      if (nlev.eq.0.or.ymax.eq.zero) then
         idone=1
      else
         call ascle(4,xmin,xmax,major,minor)
         xstep=(xmax-xmin)/major
         call ascle(4,ymin,ymax,major,minor)
         ystep=(ymax-ymin)/major
         thin=(xmax-xmin)/nden
         write(nout,'(''1'',i3,''/'')') iwcol
         it=1
         do i=1,70
            if (hk(i:i).ne.' ') it=i
         enddo
         write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
         write(nout,'(a,''<t>hreshold reactions'',a,''/'')') qu,qu
         write(nout,'(''1 0 2 1/'')')
         write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
         write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
         write(nout,'(1p,3e12.3,''/'')') ymin,ymax,ystep
         write(nout,'(a,''<c>ross section (barns)'',a,''/'')') qu,qu
         write(nout,'(''/'')')
         nlev=0
         do i=1,ntr
            mt=nint(xss(mtr+i-1))
            k=nint(xss(lsig+i-1)+sig-1)
            n=nint(xss(k+1))
            iaa=nint(xss(k))
            iflag=0
            if (nlev.eq.5) iflag=1
            if (mt.lt.5) iflag=1
            if (mt.ge.18.and.mt.le.21) iflag=1
            if (mt.eq.38) iflag=1
            if (mt.ge.50.and.mt.lt.91.and.izai.eq.1) iflag=1
            if (mt.gt.207.and.mt.lt.600) iflag=1
            if (mt.le.mtlast) iflag=1
            if (xss(iaa).lt.xsmin) iflag=1
            if (iflag.eq.0) then
               mtl=mt
               nlev=nlev+1
               if (nlev.gt.1) write(nout,'(i2,''/'')') nlev
               if (nlev.gt.1) write(nout,'(''/'')')
               icurv=mod(nlev-1,5)
               if (iwcol.eq.0) then
                  write(nout,'(''0 0'',i2,''/'')') icurv
               else
                  write(nout,'(''0 0 0'',i2,''/'')') icurv
               endif
               call mtname(mt,name,izai)
               if (name(1:1).eq.'(') then
                  if (izai.eq.1) then
                     name(2:2)='n'
                  else if (izai.eq.1001) then
                     name(2:2)='p'
                  else if (izai.eq.1002) then
                     name(2:2)='d'
                  else if (izai.eq.1003) then
                     name(2:2)='t'
                  else if (izai.eq.2003) then
                     name(2:2)='s'
                  else if (izai.eq.2004) then
                     name(2:2)='a'
                  endif
               endif
               write(nout,'(a,a,a,''/'')') qu,name,qu
               write(nout,'(''0/'')')
               xlast=small
               do j=1,n
                  x=xss(iaa+j-1)
                  iflag=0
                  if (n.le.nden.or.&
                    x.ge.xlast+thin.or.j.eq.1.or.j.eq.n) then
                     y=xss(k+2+j-1)
                     write(nout,'(1p,2e14.6,''/'')') x,y
                     xlast=x
                  endif
               enddo
               write(nout,'(''/'')')
               if (mt.lt.203.and.nint(xss(mtr+i)).ge.203.and.&
                 nint(xss(mtr+i)).le.207) nlev=5
            endif
         enddo
         mtlast=mtl
      endif
   enddo
   return
   end subroutine aplotr

   subroutine aplopf(nout,iwcol,hk)
   !-------------------------------------------------------------------
   ! Plot the higher fission reactions.
   !-------------------------------------------------------------------
   use acecm ! provides mtname
   ! externals
   integer::nout,iwcol
   character(70)::hk
   ! internals
   integer::nlev,i,mt,iflag,k,n,iaa,j,major,minor,it,icurv
   real(kr)::x,y,xstep,ystep,thin,xlast,ylast
   real(kr)::xmin,xmax,ymin,ymax
   character(10)::name
   character(1)::qu=''''
   real(kr),parameter::small=1.e-12_kr
   real(kr),parameter::zero=0
   integer,parameter::nden=4000

   !--search for the higher fission reactions
   xmin=1000
   xmax=0
   ymin=1000
   ymax=0
   nlev=0
   do i=1,ntr
      mt=nint(xss(mtr+i-1))
      iflag=0
      if (mt.lt.20) iflag=1
      if (mt.gt.21.and.mt.lt.38) iflag=1
      if (mt.gt.38) iflag=1
      if (iflag.eq.0) then
         nlev=nlev+1
         k=nint(xss(lsig+i-1)+sig-1)
         n=nint(xss(k+1))
         iaa=nint(xss(k))
         do j=1,n
            x=xss(iaa+j-1)
            y=xss(k+2+j-1)
            if (y.ne.zero.or.j.le.1) then
               if (x.lt.xmin) xmin=x
               if (x.gt.xmax) xmax=x
               if (y.lt.ymin) ymin=y
               if (y.gt.ymax) ymax=y
            endif
         enddo
      endif
   enddo
   if (nlev.ne.0.and.ymax.ne.zero) then
      call ascle(4,xmin,xmax,major,minor)
      xstep=(xmax-xmin)/major
      call ascle(4,ymin,ymax,major,minor)
      ystep=(ymax-ymin)/major
      thin=(xmax-xmin)/nden
      write(nout,'(''1'',i3,''/'')') iwcol
      it=1
      do i=1,70
         if (hk(i:i).ne.' ') it=i
      enddo
      write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
      write(nout,'(a,''<h>igh-order fission reactions'',a,''/'')')&
        qu,qu
      write(nout,'(''1 0 2 1/'')')
      write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
      write(nout,'(a,''<e>nergy (<m>e>v>)'',a,''/'')') qu,qu
      write(nout,'(1p,3e12.3,''/'')') ymin,ymax,ystep
      write(nout,'(a,''<c>ross section (barns)'',a,''/'')') qu,qu
      write(nout,'(''/'')')
      nlev=0
      do i=1,ntr
         mt=nint(xss(mtr+i-1))
         iflag=0
         if (mt.lt.20) iflag=1
         if (mt.gt.21.and.mt.lt.38) iflag=1
         if (mt.gt.38) iflag=1
         if (iflag.eq.0) then
            nlev=nlev+1
            k=nint(xss(lsig+i-1)+sig-1)
            n=nint(xss(k+1))
            iaa=nint(xss(k))
            if (nlev.gt.1) write(nout,'(i2,''/'')') nlev
            if (nlev.gt.1) write(nout,'(''/'')')
            icurv=mod(nlev-1,5)
            if (iwcol.eq.0) then
               write(nout,'(''0 0'',i2,''/'')') icurv
            else
               write(nout,'(''0 0 0'',i2,''/'')') icurv
            endif
            call mtname(mt,name,izai)
            write(nout,'(a,a,a,''/'')') qu,name,qu
            write(nout,'(''0/'')')
            xlast=small
            ylast=-1
            do j=1,n
               x=xss(iaa+j-1)
               if (n.le.nden.or.&
                 x.ge.xlast+thin.or.j.eq.1.or.j.eq.n) then
                  y=xss(k+2+j-1)
                  write(nout,'(1p,2e14.6,''/'')') x,y
                  xlast=x
                  ylast=y
               endif
            enddo
            write(nout,'(''/'')')
         endif
      enddo
   endif
   return
   end subroutine aplopf

   subroutine aplof4(nout,iwcol,hk)
   !-------------------------------------------------------------------
   ! Plot the distribution data in 3D form.
   !-------------------------------------------------------------------
   use util ! provides error
   use acecm ! provides mtname
   ! externals
   integer::nout,iwcol
   character(70)::hk
   ! internals
   integer::nr1,nn,n,na,mt,ne,nb,k,imin,ie,nk,i,j,itwo
   integer::major,minor,np,iflag,intt,it
   real(kr)::x,xstep,ystep,test,break,e,cc,pp,rat,stepm
   real(kr)::elast,ylast
   real(kr)::xmin,xmax,ymin,ymax,zmin,zmax
   integer,parameter::maxe=1200
   integer::loce(maxe)
   character(10)::name
   character(1)::qu=''''
   real(kr),parameter::eps=1.e-5_kr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::dn=.99e0_kr
   real(kr),parameter::up=1.01e0_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1

   !--loop over angular distributions
   nr1=nr+1
   do nn=1,nr1
      n=nn-1
      na=nint(xss(land+n))
      if (na.gt.0) then
         if (n.eq.0) then
            mt=2
            name='elastic'
         else
            mt=iabs(nint(xss(mtr+n-1)))
            call mtname(mt,name,izai)
            if (name(1:1).eq.'(') then
               if (izai.eq.1) then
                  name(2:2)='n'
               else if (izai.eq.1001) then
                  name(2:2)='p'
               else if (izai.eq.1002) then
                  name(2:2)='d'
               else if (izai.eq.1003) then
                  name(2:2)='t'
               else if (izai.eq.2003) then
                  name(2:2)='s'
               else if (izai.eq.2004) then
                  name(2:2)='a'
               endif
            endif
         endif
         na=na+and-1
         ne=nint(xss(na))
         nb=na+ne
         k=nint(xss(nb+1))

         !--plot equiprobable contours
         if (k.ge.0) then
            if (ne.gt.maxe) call error('aplof4',&
              'too many e values in angular distribution',' ')
            xmin=1000
            imin=ne
            do ie=1,ne
               k=nint(xss(ie+nb))
               if (k.gt.0) then
                  x=xss(na+ie)
                  if (x.lt.xmin) then
                     xmin=x
                     imin=ie
                     endif
               endif
               loce(ie)=k+and-1
            enddo
            nk=loce(ne)-loce(ne-1)
            if (imin.ne.ne.and.nk.ge.1) then
               do i=1,70
               if (hk(i:i).ne.' ') it=i
               enddo
               xmax=xss(na+ne)
               call ascle(4,xmin,xmax,major,minor)
               xstep=(xmax-xmin)/major
               write(nout,'(''1'',i3,''/'')') iwcol
               write(nout,&
                  '(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
               write(nout,&
                  '(a,''<s>cattering contours for '',a,a,''/'')')&
                qu,name,qu
               write(nout,'(''1/'')')
               write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
               write(nout,&
                  '(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
               write(nout,'(1p,3e12.3,''/'')') -1.,1.,.5
               write(nout,'(a,''<c>osine>'',a,''/'')') qu,qu
               write(nout,'(''/'')')
               do j=1,nk
                  if (j.gt.1) write(nout,'(i2,''/'')') j
                  if (j.gt.1) write(nout,'(''/'')')
                  write(nout,'(''/'')')
                  write(nout,'(''0'')')
                  do ie=imin,ne
                     k=loce(ie)+j-1
                     write(nout,'(1p,2e14.6,''/'')')&
                        xss(na+ie),xss(k)
                  enddo
                  write(nout,'(''/'')')
               enddo
            endif

         !--plot perspective view of distribution
         !--high-energy distributions are split onto two pages
         else
            xmin=-1
            xmax=1
            xstep=1
            xstep=xstep/2
            ymin=xss(na+1)
            ymax=xss(na+ne)
            itwo=0
            if (ne.gt.4) then
               test=3
               if (ymax.gt.test*xss(na+ne-1)) ymax=xss(na+ne-1)
               test=50
               if (ymax.ge.test) then
                  break=20
                  if (ymin.lt.dn*break.and.ymax.gt.up*break) itwo=1
               endif
            endif
            do while (itwo.ge.0.and.itwo.le.2)
               if (itwo.eq.1) then
                  ymax=break
               else if (itwo.eq.2) then
                  ymin=break
                  ymax=xss(na+ne)
               endif
               call ascle(4,ymin,ymax,major,minor)
               ystep=(ymax-ymin)/major
               zmin=big
               zmax=0
               do i=1,ne
                  e=xss(na+i)
                  if (e.le.(1+eps)*ymax.and.e.ge.ymin) then
                     k=nint(abs(xss(nb+i)))+and-1
                     np=nint(xss(k+1))
                     k=k+1
                     do j=1,np
                        cc=xss(k+j)
                        pp=xss(k+np+j)
                        if (pp.lt.zmin) zmin=pp
                        if (pp.gt.zmax) zmax=pp
                     enddo
                  endif
               enddo
               rat=100000
               if (zmin.le.zero) zmin=zmax/rat
               if (zmax/zmin.gt.rat)  zmin=zmax/rat
               call ascll(zmin,zmax)
               if (zmin.eq.zmax) then
                   zmax=2*zmax
                   zmin=zmax/10
               endif
               do i=1,70
                  if (hk(i:i).ne.' ') it=i
               enddo
               write(nout,'(''1'',i3,''/'')') iwcol
               write(nout,&
                  '(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
               write(nout,&
                 '(a,''angular distribution for '',a,a,''/'')')&
                 qu,name,qu
               write(nout,'(''-1 2/'')')
               write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
               write(nout,'(a,''<c>osine'',a,''/'')') qu,qu
               write(nout,'(1p,3e12.3,''/'')') ymin,ymax,ystep
               write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
               write(nout,'(1p,3e12.3,''/'')') zmin,zmax,one
               write(nout,'(a,''<p>rob/<C>os'',a,''/'')') qu,qu
               write(nout,'(''/'')')
               write(nout,'('' 15. -15. 15. -2.5 6.5 2.5/'')')
               write(nout,'(''1/'')')
               stepm=(ymax-ymin)/150
               elast=0
               do i=1,ne
                  e=xss(na+i)
                  iflag=0
                  if (e.gt.(1+eps)*ymax) iflag=1
                  if (e.lt.ymin) iflag=1
                  if (i.gt.1.and.i.lt.ne.and.ne.gt.150.and.&
                    e-elast.lt.stepm) iflag=1
                  if (iflag.eq.0) then
                     write(nout,'(1p,e14.6,''/'')') e
                     k=nint(abs(xss(nb+i)))+and-1
                     intt=nint(xss(k))
                     np=nint(xss(k+1))
                     k=k+1
                     ylast=zmin
                     do j=1,np
                        cc=xss(k+j)
                        pp=xss(k+np+j)
                        if (pp.lt.zmin) pp=zmin
                        if (intt.eq.1) write(nout,&
                           '(1p,2e14.6,''/'')') cc,ylast
                        write(nout,'(1p,2e14.6,''/'')') cc,pp
                        ylast=pp
                     enddo
                     write(nout,'(''/'')')
                     elast=e
                  endif
               enddo
               write(nout,'(''/'')')
               if (itwo.eq.0) itwo=3
               if (itwo.eq.2) itwo=3
               if (itwo.eq.1) itwo=2
            enddo
         endif
      endif
   enddo
   return
   end subroutine aplof4

   subroutine aplonu(nout,iwcol,hk)
   !-------------------------------------------------------------------
   ! Plot the total fission nubar curve.
   !-------------------------------------------------------------------
   !--externals
   integer::nout,iwcol
   character(70)::hk
   ! internals
   real(kr)::xmin,xmax,ymin,ymax,xstep,ystep,e,emax,x,y,sum
   integer::l,j,kf,n,i,ne,major,minor,it,nr1
   character(1)::qu=''''
   real(kr),parameter::big=1.0e10_kr
   real(kr),parameter::small=1.e-12_kr
   real(kr),parameter::ten=10.e0_kr
   real(kr),parameter::step=1.2e0_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1.e0_kr

   !--set up the page for the total nubar curve
   xmin=big
   xmax=0
   ymin=big
   ymax=-big
   l=nu
   j=nint(xss(l))
   kf=j
   if (kf.lt.0) then
      l=l+iabs(kf)+1
      j=nint(xss(l))
   endif
   if (j.ne.2) then
      e=xss(esz)
      emax=xss(esz+nes-1)
      l=l+1
      n=nint(xss(l))
      xmin=e
      xmax=emax
      ymin=xss(l+1)
      ymax=ymin
      do i=2,n
         ymax=ymax+xss(l+i)*emax**(i-1)
      enddo
   else
      l=l+1
      nr1=nint(xss(l))
      if (nr1.gt.0) l=l+2*nr1
      l=l+1
      ne=nint(xss(l))
      do i=1,ne
         x=xss(l+i)
         y=xss(l+i+ne)
         if (x.lt.xmin) xmin=x
         if (x.gt.xmax) xmax=x
         if (y.lt.ymin) ymin=y
         if (y.gt.ymax) ymax=y
      enddo
   endif
   call ascle(4,ymin,ymax,major,minor)
   ystep=(ymax-ymin)/major
   xstep=(xmax-xmin)/4
   write(nout,'(''1'',i3,''/'')') iwcol
   it=1
   do i=1,70
      if (hk(i:i).ne.' ') it=i
   enddo
   write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
   write(nout,'(a,''<t>otal fission nubar'',a,''/'')') qu,qu
   write(nout,'(''1 0 2/'')')
   write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
   write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
   write(nout,'(1p,3e12.3,''/'')') ymin,ymax,ystep
   write(nout,'(a,''<f>ission nubar'',a,''/'')') qu,qu
   write(nout,'(''/'')')
   write(nout,'(''/'')')
   write(nout,'(''0/'')')
   l=nu
   j=nint(xss(l))
   kf=j
   if (kf.lt.0) then
      l=l+iabs(kf)+1
      j=nint(xss(l))
   endif
   if (j.ne.2) then
      e=xss(esz)
      emax=xss(esz+nes-1)
      l=l+1
      n=nint(xss(l))
      do while (e.lt.emax)
         sum=xss(l+1)
         do i=2,n
            sum=sum+xss(l+i)*e**(i-1)
         enddo
         write(nout,'(1p,2e14.6,''/'')') e,sum
         e=e+step
      enddo
   else
      l=l+1
      nr1=nint(xss(l))
      if (nr1.gt.0) l=l+2*nr1
      l=l+1
      ne=nint(xss(l))
      do i=1,ne
         x=xss(l+i)
         y=xss(l+i+ne)
         write(nout,'(1p,2e14.6,''/'')') x,y
      enddo
   endif
   write(nout,'(''/'')')

   return
   end subroutine aplonu

   subroutine aplodd(nout,iwcol,hk)
   !-------------------------------------------------------------------
   ! Plot the distribution data in 3D form.
   !-------------------------------------------------------------------
   use acecm ! provides mtname
   ! externals
   integer::nout,iwcol
   character(70)::hk
   ! internals
   integer::it,i,n,mt,nlaw,l,j,m,law,ne,ie,loci,intt,nn
   integer::major,minor,iflag,nd,k
   real(kr)::test,e,ep,p,epl,ppl,thin,epnext,f0
   real(kr)::xstep,ystep,xlast,ylast
   real(kr)::xmin,xmax,ymin,ymax,zmin,zmax1,zmax2,zmax
   integer::locl(20)
   character(10)::name
   character(1)::qu=''''
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::small=1.e-12_kr
   real(kr),parameter::dp=.01e0_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1

   !--initialize
   it=1
   do i=1,70
      if (hk(i:i).ne.' ') it=i
   enddo

   !--loop over reactions
   do n=1,nr
      mt=iabs(nint(xss(mtr+n-1)))
      call mtname(mt,name,izai)
      if (name(1:1).eq.'(') then
         if (izai.eq.1) then
            name(2:2)='n'
         else if (izai.eq.1001) then
            name(2:2)='p'
         else if (izai.eq.1002) then
            name(2:2)='d'
         else if (izai.eq.1003) then
            name(2:2)='t'
         else if (izai.eq.2003) then
            name(2:2)='s'
         else if (izai.eq.2004) then
            name(2:2)='a'
         endif
      endif
      nlaw=1
      l=nint(xss(ldlw+n-1)+dlw-1)
      locl(1)=l
      do while (nint(xss(l)).ne.0)
         l=nint(xss(l))
         l=l+dlw-1
         nlaw=nlaw+1
         locl(nlaw)=l
      enddo
      l=iabs(nint(xss(tyr+n-1)))
      if (l.ge.100) then
         l=l-100+dlw-1
         j=nint(xss(l))
         if (j.ne.0) l=l+2*j
         l=l+1
         j=nint(xss(l))
         l=l+2*j
      endif

      !--loop over laws
      do m=1,nlaw
         l=locl(m)+1
         law=nint(xss(l))
         l=l+2
         j=nint(xss(l))
         if (j.ne.0) l=l+2*j
         l=l+1
         j=nint(xss(l))
         l=l+2*j+1

         !--law 18
         if (law.eq.18) then
            j=nint(xss(l))
            if (j.eq.2) then
               l=l+1
               j=nint(xss(l))
            endif

         !--law 4
         else if (law.eq.4) then
            zmin=1000
            zmax1=0
            zmax2=0
            j=nint(xss(l))
            if (j.ne.0) l=l+2*j
            l=l+1
            ne=nint(xss(l))
            ymin=xss(l+1)
            ymax=xss(l+ne)
            test=3
            if (ymax.gt.test*xss(l+ne-1)) ymax=xss(l+ne-1)
            do ie=2,ne
               e=xss(ie+l)
               test=ymax+ymax/1000
               if (e.le.test) then
                  loci=nint(xss(ie+ne+l)+dlw-1)
                  intt=nint(xss(loci))
                  nn=nint(xss(loci+1))
                  loci=loci+1
                  if (nn.gt.2) then
                     do j=1,nn
                        ep=xss(j+loci)
                        p=xss(j+nn+loci)
                        if (p.lt.zmin) zmin=p
                        if (p.gt.zmax2.and.p.lt.zmax1)zmax2=p
                        if (p.gt.zmax1) then
                           zmax2=zmax1
                           zmax1=p
                        endif
                     enddo
                  endif
               endif
            enddo
            if (zmax1.gt.zero) then
               if ((zmax2/zmax1.lt.1.e-4_kr).and.(zmax2.ne.zero)) then
                  zmax=zmax2
               else
                  zmax=zmax1
               endif
               zmin=zmax/10000
               call ascll(zmin,zmax)
               xmin=1000
               xmax=0
               do ie=2,ne
                  e=xss(ie+l)
                  test=test+test/1000
                  if (e.le.test) then
                     loci=nint(xss(ie+ne+l)+dlw-1)
                     intt=nint(xss(loci))
                     nn=nint(xss(loci+1))
                     loci=loci+1
                     if (nn.gt.2) then
                        epl=xss(1+loci)
                        ppl=xss(1+n+loci)
                        do j=1,nn
                           ep=xss(j+loci)
                           p=xss(j+nn+loci)
                           if (p.ge.zmin.and.epl.lt.xmin) xmin=epl
                           if (ppl.ge.zmin.and.ep.gt.xmax) xmax=ep
                           epl=ep
                           ppl=p
                        enddo
                     endif
                  endif
               enddo
               call ascle(2,xmin,xmax,major,minor)
               xstep=(xmax-xmin)/major
               call ascle(4,ymin,ymax,major,minor)
               ystep=(ymax-ymin)/major
               write(nout,'(''1'',i3,''/'')') iwcol
               write(nout,&
                  '(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
               if (izai.eq.1) write(nout,'(a,&
                 &''<n>eutron emission for '',a,a,''/'')') qu,name,qu
               if (izai.eq.1001) write(nout,'(a,&
                 &''<p>roton emission for '',a,a,''/'')') qu,name,qu
               if (izai.eq.1002) write(nout,'(a,&
                 &''<d>euteron emission for '',a,a,''/'')') qu,name,qu
               if (izai.eq.1003) write(nout,'(a,&
                 &''<t>riton emission for '',a,a,''/'')') qu,name,qu
               if (izai.eq.2003) write(nout,'(a,&
                 &''<H>e-3 emission for '',a,a,''/'')') qu,name,qu
               if (izai.eq.2004) write(nout,'(a,&
                 &''<a>lpha emission for '',a,a,''/'')') qu,name,qu
               write(nout,'(''-1 2/'')')
               write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
               write(nout,'(a,''<s>ec. <e>nergy'',a,''/'')') qu,qu
               write(nout,'(1p,3e12.3,''/'')') ymin,ymax,ystep
               write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
               write(nout,'(1p,3e12.3,''/'')') zmin,zmax,1.
               write(nout,'(a,''<p>rob/<m>e<v>'',a,''/'')') qu,qu
               write(nout,'(''/'')')
               write(nout,'(''/'')')
               write(nout,'(''1/'')')
               do ie=2,ne
                  e=xss(ie+l)
                  loci=nint(xss(ie+ne+l)+dlw-1)
                  intt=nint(xss(loci))
                  nn=nint(xss(loci+1))
                  if (nn.gt.2) then
                     loci=loci+1
                     write(nout,'(1p,e14.6,''/'')') e
                     xlast=small
                     ylast=zmin
                     do j=1,nn
                        ep=xss(j+loci)
                        p=xss(j+nn+loci)
                        if (p.lt.zmin) p=zmin
                        if (ep.le.xmax) then
                           if (intt.eq.1) write(nout,&
                              '(1p,2e14.6,''/'')') ep,ylast
                           write(nout,'(1p,2e14.6,''/'')') ep,p
                        endif
                        xlast=ep
                        ylast=p
                     enddo
                     write(nout,'(''/'')')
                  endif
               enddo
               write(nout,'(''/'')')
            endif

         !--law 44 or law 61
         else if (law.eq.44.or.law.eq.61) then
            j=nint(xss(l))
            if (j.ne.0) l=l+2*j
            l=l+1
            ne=nint(xss(l))
            ymin=xss(l+1)
            ymax=xss(l+ne)
            test=3
            if (ymax.gt.test*xss(l+ne-1)) ymax=xss(l+ne-1)
            zmin=1000
            zmax1=0
            zmax2=0
            do ie=2,ne
               e=xss(ie+l)
               test=ymax+ymax/1000
               if (e.le.test) then
                  loci=nint(xss(ie+ne+l)+dlw-1)
                  intt=mod(nint(xss(loci)),10)
                  nd=nint(xss(loci)/10)
                  nn=nint(xss(loci+1))
                  if (nn.gt.2) then
                     loci=loci+1
                     do j=1,nn
                        ep=xss(j+loci)
                        p=xss(j+nn+loci)
                        if (j.le.nd) p=100*p
                        if (p.lt.zmin) zmin=p
                        if (p.gt.zmax2.and.p.lt.zmax1)zmax2=p
                        if (p.gt.zmax1) then
                           zmax2=zmax1
                           zmax1=p
                        endif
                     enddo
                  endif
               endif
            enddo
            if (zmax1.gt.zero) then
               if ((zmax2/zmax1.lt.1.e-4_kr).and.(zmax2.ne.zero)) then
                  zmax=zmax2
               else
                  zmax=zmax1
               endif
               zmin=zmax/10000
               call ascll(zmin,zmax)
               xmin=1000
               xmax=0
               do ie=2,ne
                  e=xss(ie+l)
                  test=ymax+ymax/1000
                  if (e.le.test) then
                     loci=nint(xss(ie+ne+l)+dlw-1)
                     intt=mod(nint(xss(loci)),10)
                     nd=nint(xss(loci)/10)
                     nn=nint(xss(loci+1))
                     if (nn.gt.2) then
                        loci=loci+1
                        epl=xss(1+loci)
                        ppl=xss(1+nn+loci)
                        do j=1,nn
                           ep=xss(j+loci)
                           p=xss(j+nn+loci)
                           if (j.le.nd) p=100*p
                           if (p.ge.zmin.and.epl.lt.xmin) xmin=epl
                           if (ppl.ge.zmin.and.ep.gt.xmax) xmax=ep
                           epl=ep
                           ppl=p
                        enddo
                     endif
                  endif
               enddo
               call ascle(2,xmin,xmax,major,minor)
               xstep=(xmax-xmin)/major
               call ascle(4,ymin,ymax,major,minor)
               ystep=(ymax-ymin)/major
               zmin=zmax/10000
               call ascll(zmin,zmax)
               thin=(xmax-xmin)/500
               write(nout,'(''1'',i3,''/'')') iwcol
               write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
               if (izai.eq.1) write(nout,'(a,&
                 &''<n>eutron emission for '',a,a,''/'')') qu,name,qu
               if (izai.eq.1001) write(nout,'(a,&
                 &''<p>roton emission for '',a,a,''/'')') qu,name,qu
               if (izai.eq.1002) write(nout,'(a,&
                 &''<d>euteron emission for '',a,a,''/'')') qu,name,qu
               if (izai.eq.1003) write(nout,'(a,&
                 &''<t>riton emission for '',a,a,''/'')') qu,name,qu
               if (izai.eq.2003) write(nout,'(a,&
                 &''<H>e-3 emission for '',a,a,''/'')') qu,name,qu
               if (izai.eq.2004) write(nout,'(a,&
                 &''<a>lpha emission for '',a,a,''/'')') qu,name,qu
               write(nout,'(''-1 2/'')')
               write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
               write(nout,'(a,''<s>ec. <e>nergy'',a,''/'')') qu,qu
               write(nout,'(1p,3e12.3,''/'')') ymin,ymax,ystep
               write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
               write(nout,'(1p,3e12.3,''/'')') zmin,zmax,one
               write(nout,'(a,''<p>rob/<m>e<v>'',a,''/'')') qu,qu
               write(nout,'(''/'')')
               write(nout,'(''/'')')
               write(nout,'(''1/'')')
               do ie=2,ne
                  e=xss(ie+l)
                  test=ymax+ymax/1000
                  if (e.le.test) then
                     loci=nint(xss(ie+ne+l)+dlw-1)
                     intt=mod(nint(xss(loci)),10)
                     nd=nint(xss(loci)/10)
                     nn=nint(xss(loci+1))
                     loci=loci+1
                     write(nout,'(1p,e14.6,''/'')') e
                     xlast=small
                     ylast=zmin
                     k=0
                     do j=1,nn
                        ep=xss(j+loci)
                        iflag=0
                        if (ep.lt.xmin.or.ep.gt.xmax) iflag=1
                        if (nn.gt.500.and.&
                          ep.lt.xlast+thin.and.j.ne.1.and.j.ne.nn)&
                           iflag=1
                        if (iflag.eq.0) then
                           p=xss(j+nn+loci)
                           if (j.le.nd) p=100*p
                           if (p.lt.zmin) p=zmin
                           if (p.gt.zmax) p=zmax
                           if (j.le.nd) write(nout,&
                              '(1p,2e14.6,''/'')') ep-dp,zmin
                           if (intt.eq.1.and.j.gt.nd)&
                             write(nout,'(1p,2e14.6,''/'')') ep,ylast
                           write(nout,'(1p,2e14.6,''/'')') ep,p
                           if (j.le.nd) write(nout,&
                              '(1p,2e14.6,''/'')') ep+dp,zmin
                           k=k+1
                           xlast=ep
                           ylast=p
                        endif
                     enddo
                     if (k.lt.2) then
                        write(nout,'(1p,2e14.6,''/'')') xmin,zmin
                           write(nout,'(1p,2e14.6,''/'')') xmax,zmin
                     endif
                     write(nout,'(''/'')')
                  endif
               enddo
               write(nout,'(''/'')')
            endif

         !--law 67
         else if (law.eq.67) then
            xmin=1000
            xmax=0
            ymin=1000
            ymax=0
            zmin=1000
            zmax1=0
            zmax2=0
            j=nint(xss(l))
            if (j.ne.0) then
               l=l+2*j
            endif
            l=l+1
            ne=nint(xss(l))
            do ie=2,ne
               e=xss(ie+l)
               if (e.lt.ymin) ymin=e
               if (e.gt.ymax) ymax=e
               loci=nint(xss(ie+ne+l))
               ep=-1
               call getl7(ep,epnext,loci,xss(dlw),f0)
               do while (epnext.ne.emax)
                  ep=epnext
                  call getl7(ep,epnext,loci,xss(dlw),f0)
                  if (f0.gt.zmax2.and.f0.lt.zmax1)zmax2=f0
                  if (f0.gt.zmax1) then
                     zmax2=zmax1
                     zmax1=f0
                  endif
               enddo
            enddo
            if (zmax1.gt.zero) then
               if ((zmax2/zmax1.lt.1.e-4_kr).and.(zmax2.ne.zero)) then
                  zmax=zmax2
               else
                  zmax=zmax1
               endif
               zmin=zmax/10000
               call ascll(zmin,zmax)
               do ie=2,ne
                  e=xss(ie+l)
                  loci=nint(xss(ie+ne+l))
                  ep=-1
                  call getl7(ep,epnext,loci,xss(dlw),f0)
                  epl=ep
                  ppl=0
                  do while (epnext.ne.emax)
                     ep=epnext
                     call getl7(ep,epnext,loci,xss(dlw),f0)
                     if (f0.ge.zmin.and.epl.lt.xmin) xmin=epl
                     if (ppl.ge.zmin.and.ep.gt.xmax) xmax=ep
                     epl=ep
                     ppl=f0
                  enddo
               enddo
               call ascle(2,xmin,xmax,major,minor)
               xstep=(xmax-xmin)/major
               call ascle(4,ymin,ymax,major,minor)
               ystep=(ymax-ymin)/major
               write(nout,'(''1'',i3,''/'')') iwcol
               write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
               if (izai.eq.1) write(nout,'(a,&
                 &''<n>eutron emission for '',a,a,''/'')') qu,name,qu
               if (izai.eq.1001) write(nout,'(a,&
                 &''<p>roton emission for '',a,a,''/'')') qu,name,qu
               if (izai.eq.1002) write(nout,'(a,&
                 &''<d>euteron emission for '',a,a,''/'')') qu,name,qu
               if (izai.eq.1003) write(nout,'(a,&
                 &''<t>riton emission for '',a,a,''/'')') qu,name,qu
               if (izai.eq.2003) write(nout,'(a,&
                 &''<H>e-3 emission for '',a,a,''/'')') qu,name,qu
               if (izai.eq.2004) write(nout,'(a,&
                 &''<a>lpha emission for '',a,a,''/'')') qu,name,qu
               write(nout,'(''-1 2/'')')
               write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
               write(nout,'(a,''<s>ec. <e>nergy'',a,''/'')') qu,qu
               write(nout,'(1p,3e12.3,''/'')') ymin,ymax,ystep
               write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
               write(nout,'(1p,3e12.3,''/'')') zmin,zmax,one
               write(nout,'(a,''<p>rob/<m>e<v>'',a,''/'')') qu,qu
               write(nout,'(''/'')')
               write(nout,'(''/'')')
               write(nout,'(''1/'')')
               do ie=2,ne
                  e=xss(ie+l)
                  loci=nint(xss(ie+ne+l))
                  write(nout,'(1p,e14.6,''/'')') e
                  ep=-1
                  call getl7(ep,epnext,loci,xss(dlw),f0)
                  do while (epnext.ne.emax)
                     ep=epnext
                     call getl7(ep,epnext,loci,xss(dlw),f0)
                     if (f0.lt.zmin) f0=zmin
                     write(nout,'(1p,2e14.6,''/'')') ep,f0
                  enddo
                  write(nout,'(''/'')')
               enddo
               write(nout,'(''/'')')
            endif
         endif
      enddo
   enddo
   return
   end subroutine aplodd

   subroutine aplodn(nout,iwcol,hk)
   !-------------------------------------------------------------------
   ! Plot the delayed-neutron data.
   !-------------------------------------------------------------------
   !--externals
   integer::nout,iwcol
   character(70)::hk
   !--internals
   integer::l,j,ne,i,major,minor,loct,law,m,n,loci,intt,nn,it
   real(kr)::xmin,xmax,ymin,ymax,xstep,ystep,x,y,decay,frac,xtag,ytag
   character(1)::qu=''''
   real(kr),parameter::big=1.0e10_kr
   real(kr),parameter::small=1.e-12_kr
   real(kr),parameter::ten=10.e0_kr
   real(kr),parameter::step=1.2e0_kr
   real(kr),parameter::scale=1.e2_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1.e0_kr

   !--set up the page for the delayed nubar curve
   xmin=big
   xmax=0
   ymin=big
   ymax=-big
   l=nud
   j=nint(xss(l))
   l=l+1
   nr=nint(xss(l))
   if (nr.gt.0) l=l+2*nr
   l=l+1
   ne=nint(xss(l))
   do i=1,ne
      x=xss(l+i)
      y=xss(l+i+ne)
      if (x.lt.xmin) xmin=x
      if (x.gt.xmax) xmax=x
      if (y.lt.ymin) ymin=y
      if (y.gt.ymax) ymax=y
   enddo
   ymin=ymin/step
   ymax=ymax*step
   call ascle(4,ymin,ymax,major,minor)
   ystep=(ymax-ymin)/major
   xstep=(xmax-xmin)/4
   write(nout,'(''1'',i3,''/'')') iwcol
   it=1
   do i=1,70
      if (hk(i:i).ne.' ') it=i
   enddo
   write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
   write(nout,'(a,''<d>elayed nubar'',a,''/'')') qu,qu
   write(nout,'(''1 0 2/'')')
   write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
   write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
   write(nout,'(1p,3e12.3,''/'')') ymin,ymax,ystep
   write(nout,'(a,''<d>elayed nubar'',a,''/'')') qu,qu
   write(nout,'(''/'')')
   write(nout,'(''/'')')
   write(nout,'(''0/'')')
   l=nud
   j=nint(xss(l))
   l=l+1
   nr=nint(xss(l))
   if (nr.gt.0) l=l+2*nr
   l=l+1
   ne=nint(xss(l))
   do i=1,ne
      x=xss(l+i)
      y=xss(l+i+ne)
      write(nout,'(1p,2e14.6,''/'')') x,y
   enddo
   write(nout,'(''/'')')

   !--set up the page for the delayed spectra curves
   xmin=big
   xmax=0
   ymin=big
   ymax=-big
   do i=1,ndnf
      loct=nint(xss(i-1+ldnd)+dnd-1)
      law=nint(xss(loct+1))
      m=nint(xss(loct+3))
      loct=loct+3+2*m
      loct=loct+1
      n=nint(xss(loct))
      loct=loct+1+2*n
      m=nint(xss(loct))
      loct=loct+2*m
      loct=loct+1
      ne=nint(xss(loct))
      loci=nint(xss(1+ne+loct))+dnd-1
      intt=nint(xss(loci))
      n=nint(xss(loci+1))
      loci=loci+1
      l=dndat
      do j=1,ndnf
         if (j.eq.i) decay=xss(l)
         l=l+2
         nn=nint(xss(l))
         l=l+nn
         if (j.eq.i) frac=xss(l+1)
         l=l+nn+1
      enddo
      do j=1,n
         x=xss(j+loci)
         if (x.eq.zero) x=xss(j+1+loci)/10
         y=frac*xss(j+loci+n)
         if (x.lt.xmin) xmin=x
         if (x.gt.xmax) xmax=x
         if (y.lt.ymin) ymin=y
         if (y.gt.ymax) ymax=y
      enddo
   enddo
   call ascll(xmin,xmax)
   if (ymin.lt.ymax/scale) ymin=ymax/scale
   call ascll(ymin,ymax)
   write(nout,'(''1'',i3,''/'')') iwcol
   it=1
   do i=1,70
      if (hk(i:i).ne.' ') it=i
   enddo
   write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
   write(nout,'(a,''<d>elayed neutron spectra'',a,''/'')') qu,qu
   xtag=step*xmin
   ytag=ymax/30
   write(nout,'(''4 0 2 1'',2e12.4,''/'')') xtag,ytag
   write(nout,'(1p,3e12.3,''/'')') xmin,xmax,one
   write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
   write(nout,'(1p,3e12.3,''/'')') ymin,ymax,one
   write(nout,'(a,''<p>robability'',a,''/'')') qu,qu
   write(nout,'(''/'')')
   do i=1,ndnf
      if (i.gt.1) then
         write(nout,'(i2,''/'')') i
         write(nout,'(''/'')')
      endif
      if (iwcol.eq.0) then
         write(nout,'(''0 0'',i2,''/'')') i
      else
         write(nout,'(''0 0 0'',i2,''/'')') i
      endif
      l=dndat
      do j=1,ndnf
         if (j.eq.i) decay=xss(l)
         l=l+2
         n=nint(xss(l))
         l=l+n
         if (j.eq.i) frac=xss(l+1)
         l=l+n+1
      enddo
      write(nout,'(a,''group'',i2,'' frac'',f7.4,'' decay/shake'',&
        &1p,e10.3,a,''/'')') qu,i,frac,decay,qu
      write(nout,'(''0/'')')
      loct=nint(xss(i-1+ldnd)+dnd-1)
      law=nint(xss(loct+1))
      m=nint(xss(loct+3))
      loct=loct+3+2*m
      loct=loct+1
      n=nint(xss(loct))
      loct=loct+1+2*n
      m=nint(xss(loct))
      loct=loct+2*m
      loct=loct+1
      ne=nint(xss(loct))
      loci=nint(xss(1+ne+loct))+dnd-1
      intt=nint(xss(loci))
      n=nint(xss(loci+1))
      loci=loci+1
      m=n
      if (intt.eq.1) m=m-1
      do j=1,m
         x=xss(j+loci)
         if (x.eq.zero) x=xss(j+1+loci)/10
         y=frac*xss(j+loci+n)
         if (y.lt.ymin) y=ymin
         write(nout,'(1p,2e14.6,''/'')') x,y
         if (intt.eq.1) then
            x=xss(j+1+loci)
            y=frac*xss(j+loci+n)
            if (y.lt.ymin) y=ymin
            write(nout,'(1p,2e14.6,''/'')') x,y
         endif
      enddo
      write(nout,'(''/'')')
   enddo

   return
   end subroutine aplodn

   subroutine aplopp(nout,iwcol,hk)
   !-------------------------------------------------------------------
   ! Plot the detailed photon data.
   !-------------------------------------------------------------------
   use util ! provides sigfig
   use acecm ! provides mtname
   ! externals
   integer::nout,iwcol
   character(70)::hk
   ! internals
   integer::it,i,mti,loct,law,m,n,ne,ie1,ie,loci,intt,nd,j
   integer::major,minor,k,iflag,idone,iie,ii,mftype,mtmult
   integer::ir,mt,iaa,je,loc1,lp,iii,iimax
   real(kr)::egmth,egm14,test,e,ep,p,epl,ppl,xstep,ystep,thin
   real(kr)::xlast,ylast,delta,y,s,eg,f,eii,g
   real(kr)::xmin,xmax,ymin,ymax,zmin,zmax
   real(kr)::spect(2000)
   character(10)::name
   character(1)::qu=''''
   real(kr),parameter::small=1.e-12_kr
   real(kr),parameter::up=1.001e0_kr
   real(kr),parameter::delt=.01e0_kr
   real(kr),parameter::eth=2.53e-8_kr
   real(kr),parameter::e14=14.e0_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1

   !--initialize
   it=1
   do i=1,70
      if (hk(i:i).ne.' ') it=i
   enddo
   egmth=0
   egm14=0

   !--loop over reactions
   if (ntrp.ne.0) then
      do i=1,ntrp
         mti=nint(xss(i-1+mtrp)/1000)
         call mtname(mti,name,izai)
         loct=nint(xss(i-1+ldlwp)+dlwp-1)
         law=nint(xss(loct+1))
         if (law.eq.4) then
            loct=loct+3
            m=nint(xss(loct))
            if (m.ne.0) then
               loct=loct+2*m
            endif
            loct=loct+1
            n=nint(xss(loct))
            loct=loct+1+2*n
            m=nint(xss(loct))
            if (m.ne.0) then
               loct=loct+2*m
            endif
            loct=loct+1
            ne=nint(xss(loct))
            ymin=xss(loct+1)
            ymax=xss(loct+ne)
            test=3
            if (ne.gt.3.and.ymax.gt.test*xss(loct+ne-1))&
              ymax=xss(loct+ne-1)
            zmin=1000
            zmax=0
            ie1=1
            if (ne.gt.2.and.mti.ne.102.and.mti.ne.18) ie1=2
            do ie=ie1,ne
               e=xss(ie+loct)
               if (e.le.up*ymax) then
                  loci=nint(xss(ie+ne+loct)+dlwp-1)
                  intt=mod(nint(xss(loci)),10)
                  nd=nint(xss(loci))/10
                  n=nint(xss(loci+1))
                  loci=loci+1
                  if (nd.gt.0) then
                     do j=1,nd
                        ep=xss(j+loci)
                        p=xss(j+n+loci)*100
                        if (p.lt.zmin) zmin=p
                        if (p.gt.zmax) zmax=p
                     enddo
                  endif
                  if (n.gt.nd.and.n.gt.2) then
                     do j=nd+1,n
                        ep=xss(j+loci)
                        p=xss(j+n+loci)
                        if (intt.ne.1.or.j.ne.n) then
                           if (p.lt.zmin) zmin=p
                           if (p.gt.zmax) zmax=p
                        endif
                     enddo
                  endif
               endif
            enddo
            if (zmax.gt.zero) then
               zmin=zmax/10000
               call ascll(zmin,zmax)
               xmin=1000
               xmax=0
               do ie=ie1,ne
                  e=xss(ie+loct)
                  if (e.le.up*ymax) then
                     loci=nint(xss(ie+ne+loct)+dlwp-1)
                     intt=mod(nint(xss(loci)),10)
                     nd=nint(xss(loci))/10
                     n=nint(xss(loci+1))
                     loci=loci+1
                     if (nd.gt.0) then
                        do j=1,nd
                           ep=xss(j+loci)
                           if (ep.lt.zero) ep=-ep+e*aw0/(aw0+awi)
                           if (ep-delt.lt.xmin) xmin=ep-delt
                           if (ep+delt.gt.xmax) xmax=ep+delt
                        enddo
                     endif
                     if (n.gt.nd) then
                        epl=xss(1+nd+loci)
                        ppl=xss(1+nd+n+loci)
                        do j=nd+1,n
                           ep=xss(j+loci)
                           p=xss(j+n+loci)
                           if (p.ge.zmin.and.epl.lt.xmin) xmin=epl
                           if (ppl.ge.zmin.and.ep.gt.xmax) xmax=ep
                           epl=ep
                           ppl=p
                        enddo
                     endif
                     if (e.le.100*eth.and.xmax.gt.egmth)&
                        egmth=xmax
                     if (e.le.e14+e14/10.and.xmax.gt.egm14)&
                        egm14=xmax
                  endif
               enddo
               if (xmax-xmin.lt.xmax/10) then
                  xmin=xmin-xmin/10
                  xmax=xmax+xmax/10
               endif
               call ascle(1,xmin,xmax,major,minor)
               xstep=(xmax-xmin)/major
               call ascle(4,ymin,ymax,major,minor)
               ystep=(ymax-ymin)/major
               thin=(xmax-xmin)/500
               write(nout,'(''1'',i3,''/'')') iwcol
               write(nout,&
                 '(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
               write(nout,&
                 '(a,''<p>hoton emission for '',a,a,''/'')') qu,name,qu
               write(nout,'(''-1 2/'')')
               write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
               write(nout,&
                 '(a,''<e>#L.3H.8]g#HXLX> (<m>e<v>)'',a,''/'')') qu,qu
               write(nout,'(1p,3e12.3,''/'')') ymin,ymax,ystep
               write(nout,&
                 '(a,''<e>#L.3H.8>n#HXLX> (<m>e<v>)'',a,''/'')') qu,qu
               write(nout,'(1p,3e12.3,''/'')') zmin,zmax,one
               write(nout,'(a,''<p>rob/<m>e<v>'',a,''/'')') qu,qu
               write(nout,'(''/'')')
               write(nout,'(''/'')')
               write(nout,'(''1/'')')
               ie1=1
               if (ne.gt.2.and.mti.ne.102.and.mti.ne.18) ie1=2
               do ie=ie1,ne
                  e=xss(ie+loct)
                  if (e.le.up*ymax) then
                     loci=nint(xss(ie+ne+loct)+dlwp-1)
                     intt=mod(nint(xss(loci)),10)
                     nd=nint(xss(loci))/10
                     n=nint(xss(loci+1))
                     write(nout,'(1p,e14.6,''/'')') e
                     loci=loci+1
                     xlast=small
                     ylast=zmin
                     k=0
                     do j=1,n
                        ep=xss(j+loci)
                        if (j.le.nd.and.ep.lt.zero)&
                           ep=-ep+e*aw0/(aw0+awi)
                        iflag=0
                        if (ep.lt.xmin.or.ep.gt.xmax) iflag=1
                        if (n.gt.500.and.ep.lt.xlast+thin.and.&
                          j.ne.1.and.j.ne.n) iflag=1
                        if (iflag.eq.0) then
                           p=xss(j+n+loci)
                           if (j.lt.nd) p=100*p
                           if (p.lt.zmin) p=zmin
                           if (p.gt.zmax) p=zmax
                           if (j.le.nd) write(nout,&
                              '(1p,2e14.6,''/'')') ep-delt,zmin
                           if (intt.eq.1.and.j.gt.nd)&
                             write(nout,'(1p,2e14.6,''/'')')&
                             ep,ylast
                           write(nout,'(1p,2e14.6,''/'')') ep,p
                           if (j.le.nd) write(nout,&
                              '(1p,2e14.6,''/'')') ep+delt,zmin
                           k=k+1
                           xlast=ep
                           ylast=p
                        endif
                     enddo
                     if (k.lt.2) then
                        write(nout,'(1p,2e14.6,''/'')') xmin,zmin
                        write(nout,'(1p,2e14.6,''/'')') xmax,zmin
                     endif
                     write(nout,'(''/'')')
                  endif
               enddo
               write(nout,'(''/'')')
            endif
         endif
      enddo
   endif
   egmth=sigfig(egmth+egmth/20,2,0)
   egm14=sigfig(egm14+egm14/20,2,0)
   if (egmth.eq.zero) egmth=10
   if (egm14.eq.zero) egm14=25

   !--we only do the spectra for neutrons
   if (izai.ne.1) return

   !--plot the thermal capture photon spectrum
   delta=egmth/2000
   if (ntrp.ne.0) then
      idone=0
      iie=0
      do while (iie.lt.nes.and.idone.eq.0)
         iie=iie+1
         if (xss(esz+iie-1).gt.eth) idone=1
      enddo
      do ii=1,2000
         spect(ii)=0
      enddo
      ymin=0
      ymax=0
      do i=1,ntrp
         mti=nint(xss(i-1+mtrp)/1000)
         if (mti.eq.102) then
            loct=nint(xss(i-1+lsigp)+sigp-1)
            mftype=nint(xss(loct))
            loct=loct+1
            if (mftype.ne.13) then
               mtmult=nint(xss(loct))
               loct=loct+1
               m=nint(xss(loct))
               intt=2
               if (m.gt.0) intt=nint(xss(loct+2))
               loct=loct+1+2*m
               n=nint(xss(loct))
               idone=0
               ii=0
               do while (ii.lt.n-1.and.idone.eq.0)
                  ii=ii+1
                  if (xss(loct+ii+1).gt.eth) idone=1
               enddo
               if (intt.eq.1) then
                  y=xss(loct+ii+n)
               else
                  y=(eth-xss(loct+ii))*xss(loct+ii+1+n)&
                    /(xss(loct+ii+1)-xss(loct+ii))&
                    +(xss(loct+ii+1)-eth)*xss(loct+ii+n)&
                    /(xss(loct+ii+1)-xss(loct+ii))
               endif
               do ir=1,ntr
                  mt=nint(xss(mtr+ir-1))
                  if (mt.eq.mtmult) then
                     k=nint(xss(lsig+ir-1)+sig-1)
                     n=nint(xss(k+1))
                     iaa=nint(xss(k))
                     s=xss(k+2+iie-iaa)
                  endif
               enddo
            else
               je=nint(xss(loct))
               s=0
               if (iie.ge.je) s=xss(loct+2+iie-je)
               y=1
            endif
            loc1=nint(xss(i-1+ldlwp)+dlwp-1)
            law=nint(xss(loc1+1))
            if (law.eq.2) then
               lp=nint(xss(loc1+9))
               eg=xss(loc1+10)
               iii=max0(1,nint(eg/delta))
               if (iii.le.2000) then
                  spect(iii)=spect(iii)+y*s/delta
                  if (spect(iii).gt.ymax) ymax=spect(iii)
               endif
            endif
            if (law.eq.4) then
               loc1=loc1+3
               m=nint(xss(loc1))
               if (m.ne.0) then
                  loc1=loc1+2*m
               endif
               loc1=loc1+1
               n=nint(xss(loc1))
               loc1=loc1+1+2*n
               m=nint(xss(loc1))
               if (m.ne.0) then
                  loc1=loc1+2*m
               endif
               loc1=loc1+1
               ne=nint(xss(loc1))
               do ie=1,ne
                  iflag=0
                  if (ie.lt.ne.and.&
                    xss(ie+loc1).le.eth.and.xss(ie+1+loc1).gt.eth)&
                    iflag=1
                  if (ie.gt.1.and.&
                    xss(ie-1+loc1).le.eth.and.xss(ie+loc1).gt.eth)&
                    iflag=1
                  if (iflag.eq.1) then
                     if (xss(ie+loc1).le.eth) then
                        f=(xss(ie+1+loc1)-eth)&
                          /(xss(ie+1+loc1)-xss(ie+loc1))
                     else
                        f=(eth-xss(ie-1+loc1))&
                          /(xss(ie+loc1)-xss(ie-1+loc1))
                     endif
                     e=xss(ie+loc1)
                     loci=nint(xss(ie+ne+loc1)+dlwp-1)
                     intt=mod(nint(xss(loci)),10)
                     nd=nint(xss(loci))/10
                     n=nint(xss(loci+1))
                     loci=loci+1
                     if (nd.ne.0) then
                        do j=1,nd
                           ep=xss(j+loci)
                           p=xss(j+n+loci)
                           iii=max0(1,nint(ep/delta))
                           if (iii.le.2000) then
                              spect(iii)=spect(iii)+f*p/delta
                              if (spect(iii).gt.ymax) ymax=spect(iii)
                           endif
                        enddo
                     endif
                     do ii=1,2000
                        eii=delta*ii
                        do j=nd+1,n
                           iflag=0
                           if (intt.eq.1.and.j.lt.n.and.&
                             xss(j+loci).le.eii.and.&
                             xss(j+loci+1).gt.eii) iflag=1
                           if (intt.ne.1.and.j.lt.n.and.&
                             xss(j+loci).le.eii.and.&
                             xss(j+loci+1).gt.eii) iflag=1
                           if (intt.ne.1.and.j.gt.1.and.&
                             xss(j-1+loci).le.eii.and.&
                             xss(j+loci).gt.eii) iflag=1
                           if (iflag.eq.1) then
                              if (intt.eq.1) then
                                 g=1
                              else if (xss(j+loci).le.eii) then
                                 g=(xss(j+loci+1)-eii)&
                                    /(xss(j+loci+1)-xss(j+loci))
                              else
                                 g=(eii-xss(j+loci-1))&
                                    /(xss(j+loci)-xss(j+loci-1))
                              endif
                              ep=xss(j+loci)
                              p=xss(j+n+loci)
                              spect(ii)=spect(ii)+y*s*p*f*g
                              if (spect(ii).gt.ymax) ymax=spect(ii)
                           endif
                        enddo
                     enddo
                  endif
               enddo
            endif
         endif
      enddo
      if (ymax.gt.zero) then
         ymin=ymax/1000
         call ascll(ymin,ymax)
         do ii=1,2000
            if (spect(ii).lt.ymin) spect(ii)=ymin
            if (spect(ii).gt.ymin) iimax=ii
         enddo
         xmin=delta
         xmax=delta*(iimax+1)
         call ascle(2,xmin,xmax,major,minor)
         xstep=(xmax-xmin)/major
         write(nout,'(''1'',i3,''/'')') iwcol
         write(nout,&
           '(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
         write(nout,&
           '(a,''thermal capture photon spectrum'',a,''/'')') qu,qu
         write(nout,'(''2/'')')
         write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
         write(nout,&
           '(a,''<g>amma <e>nergy (<m>e<v>)'',a,''/'')') qu,qu
         write(nout,'(1p,3e12.3,''/'')') ymin,ymax,one
         write(nout,&
           '(a,''<g>amma <p>rod (barns/<m>e<v>)'',a,''/'')') qu,qu
         write(nout,'(''/'')')
         write(nout,'(''/'')')
         write(nout,'(''0/'')')
         do ie=1,iimax
            e=delta*ie
            p=spect(ie)
            if (p.lt.ymin) p=ymin
            write(nout,'(1p,2e14.6,''/'')') e,p
         enddo
         write(nout,'(''/'')')
      endif
   endif

   !--plot the net photon spectrum for 14 mev neutrons
   delta=egm14/2000
   if (ntrp.ne.0) then
      idone=0
      i=0
      do while (i.lt.nes.and.idone.eq.0)
         i=i+1
         iie=nes-i
         if (xss(esz+iie-1).le.e14) idone=1
      enddo
      do ii=1,2000
         spect(ii)=0
      enddo
      ymin=0
      ymax=0
      do i=1,ntrp
         mti=nint(xss(i-1+mtrp)/1000)
         loct=nint(xss(i-1+lsigp)+sigp-1)
         mftype=nint(xss(loct))
         loct=loct+1
         if (mftype.ne.13) then
            mtmult=nint(xss(loct))
            loct=loct+1
            m=nint(xss(loct))
            loct=loct+1
            n=nint(xss(loct))
            idone=0
            ii=0
            do while (idone.eq.0.and.ii.lt.n)
               ii=ii+1
               if (xss(loct+ii+1).gt.e14) idone=1
            enddo
            y=(e14-xss(loct+ii))*xss(loct+ii+1+n)&
              /(xss(loct+ii+1)-xss(loct+ii))&
              +(xss(loct+ii+1)-e14)*xss(loct+ii+n)&
              /(xss(loct+ii+1)-xss(loct+ii))
            do ir=1,ntr
               mt=nint(xss(mtr+ir-1))
               if (mt.eq.mtmult) then
                  k=nint(xss(lsig+ir-1)+sig-1)
                  n=nint(xss(k+1))
                  iaa=nint(xss(k))
                  s=xss(k+2+iie-iaa)
               endif
            enddo
         else
            je=nint(xss(loct))
            s=0
            if (iie.ge.je) s=xss(loct+2+iie-je)
            y=1
         endif
         loc1=nint(xss(i-1+ldlwp)+dlwp-1)
         law=nint(xss(loc1+1))
         if (law.eq.2) then
            lp=nint(xss(loc1+9))
            eg=xss(loc1+10)
            if (lp.eq.2) eg=eg+aw0*e14/(aw0+1)
            iii=max0(1,nint(eg/delta))
            if (iii.le.2000) then
               spect(iii)=spect(iii)+y*s/delta
               if (spect(iii).gt.ymax) ymax=spect(iii)
            endif
         endif
         if (law.eq.4) then
            loc1=loc1+3
            m=nint(xss(loc1))
            if (m.ne.0) then
               loc1=loc1+2*m
            endif
            loc1=loc1+1
            n=nint(xss(loc1))
            loc1=loc1+1+2*n
            m=nint(xss(loc1))
            if (m.ne.0) then
               loc1=loc1+2*m
            endif
            loc1=loc1+1
            ne=nint(xss(loc1))
            do ie=1,ne
               iflag=0
               if (ie.lt.ne.and.&
                 xss(ie+loc1).le.e14.and.&
                 xss(ie+1+loc1).gt.e14) iflag=1
               if (ie.gt.1.and.&
                 xss(ie-1+loc1).le.e14.and.&
                 xss(ie+loc1).gt.e14) iflag=1
               if (iflag.eq.1) then
                  if (xss(ie+loc1).le.e14) then
                     f=(xss(ie+1+loc1)-e14)&
                       /(xss(ie+1+loc1)-xss(ie+loc1))
                  else
                     f=(e14-xss(ie-1+loc1))&
                       /(xss(ie+loc1)-xss(ie-1+loc1))
                  endif
                  e=xss(ie+loc1)
                  loci=nint(xss(ie+ne+loc1)+dlwp-1)
                  intt=mod(nint(xss(loci)),10)
                  nd=nint(xss(loci))/10
                  n=nint(xss(loci+1))
                  loci=loci+1
                  if (nd.ne.0) then
                     do j=1,nd
                        ep=xss(j+loci)
                        if (ep.lt.zero) ep=-ep+e14*aw0/(aw0+awi)
                        p=xss(j+n+loci)
                        iii=max0(1,nint(ep/delta))
                        if (iii.le.2000) then
                           spect(iii)=spect(iii)+f*p/delta
                           if (spect(iii).gt.ymax) ymax=spect(iii)
                        endif
                     enddo
                  endif
                  if (n.gt.nd) then
                     do ii=1,2000
                        eii=delta*ii
                        do j=nd+1,n
                           iflag=0
                           if (intt.eq.1.and.j.lt.n.and.&
                             xss(j+loci).le.eii.and.&
                             xss(j+loci+1).gt.eii) iflag=1
                           if (intt.ne.1.and.j.lt.n.and.&
                             xss(j+loci).le.eii.and.&
                             xss(j+loci+1).gt.eii) iflag=1
                           if (intt.ne.1.and.j.gt.1.and.&
                             xss(j-1+loci).le.eii.and.&
                             xss(j+loci).gt.eii) iflag=1
                           if (iflag.eq.1) then
                              if (intt.eq.1) then
                                 g=1
                              else if (xss(j+loci).le.eii) then
                                 g=(xss(j+loci+1)-eii)&
                                    /(xss(j+loci+1)-xss(j+loci))
                              else
                                 g=(eii-xss(j+loci-1))&
                                    /(xss(j+loci)-xss(j+loci-1))
                              endif
                              ep=xss(j+loci)
                              p=xss(j+n+loci)
                              spect(ii)=spect(ii)+y*s*p*f*g
                              if (spect(ii).gt.ymax) ymax=spect(ii)
                           endif
                        enddo
                     enddo
                  endif
               endif
            enddo
         endif
      enddo
      if (ymax.gt.zero) then
         ymin=ymax
         do ii=1,2000
            if (spect(ii).gt.zero.and.spect(ii).lt.ymin) ymin=spect(ii)
         enddo
         ymin=ymin/2
         call ascll(ymin,ymax)
         do ii=1,2000
            if (spect(ii).lt.ymin) spect(ii)=ymin
            if (spect(ii).gt.ymin) iimax=ii
         enddo
         xmin=delta
         xmax=delta*(iimax+1)
         call ascle(2,xmin,xmax,major,minor)
         xstep=(xmax-xmin)/major
         write(nout,'(''1'',i3,''/'')') iwcol
         write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
         write(nout,'(a,''14 <m>e<v> photon spectrum'',a,''/'')') qu,qu
         write(nout,'(''2/'')')
         write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
         write(nout,'(a,''<g>amma <e>nergy (<m>e<v>)'',a,''/'')') qu,qu
         write(nout,'(1p,3e12.3,''/'')') ymin,ymax,one
         write(nout,'(a,''<g>amma <p>rod (barns/<m>e<v>)'',a,''/'')')&
           qu,qu
         write(nout,'(''/'')')
         write(nout,'(''/'')')
         write(nout,'(''0/'')')
         do ie=1,iimax
            e=delta*ie
            p=spect(ie)
            if (p.lt.ymin) p=ymin
            write(nout,'(1p,2e14.6,''/'')') e,p
         enddo
         write(nout,'(''/'')')
      endif
   endif

   !--finished
   return
   end subroutine aplopp

   subroutine aploxp(nout,iwcol,hk)
   !-------------------------------------------------------------------
   ! Plot the particle production data.
   !-------------------------------------------------------------------
   use acecm ! provides mtname
   ! externals
   integer::nout,iwcol
   character(70)::hk
   ! internals
   integer::i,iaa,naa,ie,major,minor,it,ii,j,imt,mt,l1,l2
   integer::ne,law,l3,loci,intt,nn,nd,k,na,nb,itwo,np
   integer::ic,l,iflag,i1
   integer::ipt,mtrh,nmtr
   integer hpd,lsigh,sigh,landh,andh,ldlwh,dlwh,yh
   real(kr)::test,e,xs,xstep,ystep,xtag,ytag,thin,xlast
   real(kr)::ep,pd,ylast,break,cc,pp,rat,stepm,elast
   real(kr)::xmin,xmax,ymin,ymax,zmin,zmax,zmax1,zmax2,heat
   character(10)::name
   character(1)::qu=''''
   integer,parameter::nden=4000
   real(kr),parameter::delta=.01e0_kr
   real(kr),parameter::up=1.01e0_kr
   real(kr),parameter::dn=0.99e0_kr
   real(kr),parameter::eps=1.e-5_kr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::small=1.e-12_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1

   !--lin-lin plot of particle heating values
   xmin=big
   xmax=0
   ymin=big
   ymax=-big
   do i=1,ntype
      ipt=nint(xss(ptype+i-1))
      hpd=nint(xss(ploct+10*(i-1)))
      iaa=nint(xss(hpd))
      naa=nint(xss(hpd+1))
      test=1
      test=test/5
      do ie=1,naa
         e=xss(esz+iaa+ie-2)
         if (e.gt.test) then
            xs=xss(hpd+1+naa+ie)
            if (e.lt.xmin) xmin=e
            if (e.gt.xmax) xmax=e
            if (xs.lt.ymin) ymin=xs
            if (xs.gt.ymax) ymax=xs
         endif
      enddo
   enddo
   if (ymax.ne.zero) then
      call ascle(4,xmin,xmax,major,minor)
      xstep=(xmax-xmin)/major
      call ascle(4,ymin,ymax,major,minor)
      ystep=(ymax-ymin)/major
      write(nout,'(''1'',i3,''/'')') iwcol
      it=1
      do i=1,70
         if (hk(i:i).ne.' ') it=i
      enddo
      write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
      write(nout,'(a,''<p>article heating contributions'',&
        &a,''/'')') qu,qu
      xtag=95*xmin/100+5*xmax/100
      ytag=ymin/10+9*ymax/10
      write(nout,'(''1 0 2 1'',2e12.4,''/'')') xtag,ytag
      write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
      write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
      write(nout,'(1p,3e12.3,''/'')') ymin,ymax,ystep
      write(nout,'(a,''<m>e<v>/collision'',a,''/'')') qu,qu
      write(nout,'(''/'')')
      ii=1
      ipt=nint(xss(ptype+ii-1))
      hpd=nint(xss(ploct+10*(ii-1)))
      iaa=nint(xss(hpd))
      naa=nint(xss(hpd+1))
      write(nout,'(''0 0 0'',i2,''/'')') ii-1
      if (ipt.eq.1) write(nout,'(a,''neutrons'',a,''/'')') qu,qu
      if (ipt.eq.2) write(nout,'(a,''photons'',a,''/'')') qu,qu
      if (ipt.eq.9) write(nout,'(a,''protons'',a,''/'')') qu,qu
      if (ipt.eq.31) write(nout,'(a,''deuterons'',a,''/'')') qu,qu
      if (ipt.eq.32) write(nout,'(a,''tritons'',a,''/'')') qu,qu
      if (ipt.eq.33) write(nout,'(a,''he-3'',a,''/'')') qu,qu
      if (ipt.eq.34) write(nout,'(a,''alphas'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      thin=(xmax-xmin)/nden
      xlast=small
      j=0
      test=1
      test=test/5
      do i=1,naa
         e=xss(esz+iaa-2+i)
         if (e.gt.test) then
            xs=xss(hpd+1+naa+i)
            if (naa.le.nden.or.e.ge.xlast+thin.or.i.eq.naa) then
               write(nout,'(1p,2e14.6,''/'')') e,xs
               xlast=e
               j=j+1
            endif
         endif
      enddo
      write(nout,'(''/'')')
      if (ntype.ne.1) then
         do ii=2,ntype
            ipt=nint(xss(ptype+ii-1))
            hpd=nint(xss(ploct+10*(ii-1)))
            iaa=nint(xss(hpd))
            naa=nint(xss(hpd+1))
            write(nout,'(''2/'')')
            write(nout,'(''/'')')
            write(nout,'(''0 0 0'',i2,''/'')') ii-1
            if (ipt.eq.1) write(nout,&
              '(a,''neutrons'',a,''/'')') qu,qu
            if (ipt.eq.2) write(nout,&
              '(a,''photons'',a,''/'')') qu,qu
            if (ipt.eq.9) write(nout,&
              '(a,''protons'',a,''/'')') qu,qu
            if (ipt.eq.31) write(nout,&
              '(a,''deuterons'',a,''/'')') qu,qu
            if (ipt.eq.32) write(nout,&
              '(a,''tritons'',a,''/'')') qu,qu
            if (ipt.eq.33) write(nout,&
              '(a,''he-3'',a,''/'')') qu,qu
            if (ipt.eq.34) write(nout,&
              '(a,''alphas'',a,''/'')') qu,qu
            write(nout,'(''0/'')')
            j=0
            xlast=small
            test=1
            test=test/5
            do i=1,naa
               e=xss(esz+iaa-2+i)
               if (e.gt.test) then
                  xs=xss(hpd+1+naa+i)
                  if (naa.le.nden.or.e.ge.xlast+thin.or.&
                    i.eq.naa) then
                     write(nout,'(1p,2e14.6,''/'')') e,xs
                     xlast=e
                     j=j+1
                  endif
               endif
            enddo
            write(nout,'(''/'')')
         enddo
      endif
   endif

   !-plot lin-lin recoil heating
   xmin=big
   xmax=0
   ymin=big
   ymax=-big
   do i=1,nes
      e=xss(esz-1+i)
      heat=xss(esz+4*nes-1+i)
      do j=1,ntype
         ipt=nint(xss(ptype+j-1))
         hpd=nint(xss(ploct+10*(j-1)))
         iaa=nint(xss(hpd))
         naa=nint(xss(hpd+1))
         ie=i-iaa-1
         if (ie.ge.1.and.ie.le.naa) heat=heat-xss(hpd+1+naa+ie)
      enddo
      if (e.lt.xmin) xmin=e
      if (e.gt.xmax) xmax=e
      if (heat.lt.ymin) ymin=heat
      if (heat.gt.ymax) ymax=heat
   enddo
   if (ymin.ne.zero.or.ymax.ne.zero) then
      if (ymin.lt.zero.and.ymax.lt.-ymin/2) ymax=-ymin/2
      call ascle(4,xmin,xmax,major,minor)
      xstep=(xmax-xmin)/major
      call ascle(4,ymin,ymax,major,minor)
      ystep=(ymax-ymin)/major
      write(nout,'(''1'',i3,''/'')') iwcol
      it=1
      do i=1,70
         if (hk(i:i).ne.' ') it=i
      enddo
      write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
      write(nout,'(a,''<r>ecoil <h>eating'',a,''/'')') qu,qu
      write(nout,'(''1 0 2 1/'')')
      write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
      write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
      write(nout,'(1p,3e12.3,''/'')') ymin,ymax,ystep
      write(nout,'(a,''<h>eating (<m>e<v>/reaction)'',a,''/'')') qu,qu
      write(nout,'(''/'')')
      write(nout,'(''/'')')
      write(nout,'(a,''recoil heating'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      thin=(xmax-xmin)/nden
      xlast=small
      j=0
      do i=1,nes
         e=xss(esz-1+i)
         test=1
         test=test/5
         if (e.ge.test) then
            if (nes.le.nden.or.e.ge.xlast+thin.or.i.eq.nes) then
               heat=xss(esz+4*nes-1+i)
               do k=1,ntype
                  ipt=nint(xss(ptype+k-1))
                  hpd=nint(xss(ploct+10*(k-1)))
                  iaa=nint(xss(hpd))
                  naa=nint(xss(hpd+1))
                  ie=i-iaa-1
                  if (ie.ge.1.and.ie.le.naa) heat=heat-xss(hpd+1+naa+ie)
               enddo
               j=j+1
               write(nout,'(1p,2e14.6,''/'')') e,heat
               xlast=e
            endif
         endif
      enddo
      write(nout,'(''/'')')
   endif

   !--lin-lin plot of particle production cross sections
   xmin=big
   xmax=0
   ymin=big
   ymax=-big
   do i=1,ntype
      ipt=nint(xss(ptype+i-1))
      hpd=nint(xss(ploct+10*(i-1)))
      iaa=nint(xss(hpd))
      naa=nint(xss(hpd+1))
      test=1
      test=test/5
      do ie=1,naa
         e=xss(esz+iaa+ie-2)
         if (e.gt.test) then
            xs=xss(hpd+1+ie)
            if (ipt.eq.2) xs=xs/5
            if (e.lt.xmin) xmin=e
            if (e.gt.xmax) xmax=e
            if (xs.lt.ymin) ymin=xs
            if (xs.gt.ymax) ymax=xs
         endif
      enddo
   enddo
   call ascle(4,xmin,xmax,major,minor)
   xstep=(xmax-xmin)/major
   call ascle(4,ymin,ymax,major,minor)
   ystep=(ymax-ymin)/major
   write(nout,'(''1'',i3,''/'')') iwcol
   it=1
   do i=1,70
      if (hk(i:i).ne.' ') it=i
   enddo
   write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
   write(nout,'(a,''<p>article production cross sections'',&
     &a,''/'')') qu,qu
   xtag=95*xmin/100+5*xmax/100
   test=30
   if (xmax.gt.test) xtag=35*xmin/100+65*xmax/100
   ytag=ymin/10+9*ymax/10
   write(nout,'(''1 0 2 1'',2e12.4,''/'')') xtag,ytag
   write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
   write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
   write(nout,'(1p,3e12.3,''/'')') ymin,ymax,ystep
   write(nout,'(a,''<c>ross section (barns)'',a,''/'')') qu,qu
   write(nout,'(''/'')')
   ii=1
   ipt=nint(xss(ptype+ii-1))
   hpd=nint(xss(ploct+10*(ii-1)))
   iaa=nint(xss(hpd))
   naa=nint(xss(hpd+1))
   write(nout,'(''0 0 0'',i2,''/'')') ii-1
   if (ipt.eq.1) write(nout,'(a,''neutrons'',a,''/'')') qu,qu
   if (ipt.eq.2) write(nout,'(a,''photons/5'',a,''/'')') qu,qu
   if (ipt.eq.9) write(nout,'(a,''protons'',a,''/'')') qu,qu
   if (ipt.eq.31) write(nout,'(a,''deuterons'',a,''/'')') qu,qu
   if (ipt.eq.32) write(nout,'(a,''tritons'',a,''/'')') qu,qu
   if (ipt.eq.33) write(nout,'(a,''he-3'',a,''/'')') qu,qu
   if (ipt.eq.34) write(nout,'(a,''alphas'',a,''/'')') qu,qu
   write(nout,'(''0/'')')
   j=0
   thin=(xmax-xmin)/nden
   xlast=small
   test=1
   test=test/5
   do i=1,naa
      e=xss(esz+iaa-2+i)
      if (e.gt.test) then
         xs=xss(hpd+1+i)
         if (ipt.eq.2) xs=xs/5
         if (naa.le.nden.or.e.ge.xlast+thin.or.i.eq.naa) then
            write(nout,'(1p,2e14.6,''/'')') e,xs
            xlast=e
            j=j+1
         endif
      endif
   enddo
   write(nout,'(''/'')')
   if (ntype.ne.1) then
      do ii=2,ntype
         ipt=nint(xss(ptype+ii-1))
         hpd=nint(xss(ploct+10*(ii-1)))
         iaa=nint(xss(hpd))
         naa=nint(xss(hpd+1))
         write(nout,'(''2/'')')
         write(nout,'(''/'')')
         write(nout,'(''0 0 0'',i2,''/'')') ii-1
         if (ipt.eq.1) write(nout,'(a,''neutrons'',a,''/'')') qu,qu
         if (ipt.eq.2) write(nout,'(a,''photons/5'',a,''/'')') qu,qu
         if (ipt.eq.9) write(nout,'(a,''protons'',a,''/'')') qu,qu
         if (ipt.eq.31) write(nout,'(a,''deuterons'',a,''/'')') qu,qu
         if (ipt.eq.32) write(nout,'(a,''tritons'',a,''/'')') qu,qu
         if (ipt.eq.33) write(nout,'(a,''he-3'',a,''/'')') qu,qu
         if (ipt.eq.34) write(nout,'(a,''alphas'',a,''/'')') qu,qu
         write(nout,'(''0/'')')
         j=0
         xlast=small
         test=1
         test=test/5
         do i=1,naa
            e=xss(esz+iaa-2+i)
            if (e.gt.test) then
               xs=xss(hpd+1+i)
               if (ipt.eq.2) xs=xs/5
               if (naa.le.nden.or.e.ge.xlast+thin.or.i.eq.naa) then
                  write(nout,'(1p,2e14.6,''/'')') e,xs
                  xlast=e
                  j=j+1
               endif
            endif
         enddo
         write(nout,'(''/'')')
      enddo
   endif

   !--loop over each production type
   do i=1,ntype
      ipt=nint(xss(ptype+i-1))
      hpd=nint(xss(ploct+10*(i-1)))
      mtrh=nint(xss(ploct+10*(i-1)+1))
      nmtr=nint(xss(ntro+i-1))
      lsigh=nint(xss(ploct+10*(i-1)+3))
      sigh=nint(xss(ploct+10*(i-1)+4))
      landh=nint(xss(ploct+10*(i-1)+5))
      andh=nint(xss(ploct+10*(i-1)+6))
      ldlwh=nint(xss(ploct+10*(i-1)+7))
      dlwh=nint(xss(ploct+10*(i-1)+8))
      yh=nint(xss(ploct+10*(i-1)+9))

      !--loop over the reactions contributing to this production
      do imt=1,nmtr
         mt=nint(xss(mtrh+imt-1))
         l1=nint(xss(lsigh+imt-1))
         l2=sigh+l1-1
         ne=nint(xss(l2+3))
         l1=nint(xss(ldlwh+imt-1))
         l2=dlwh+l1-1
         law=nint(xss(l2+1))
         if (law.eq.4.or.law.eq.44.or.law.eq.61) then

            !--make a 3d plot of the particle emission for reaction
            call mtname(mt,name,izai)
            if (name(1:1).eq.'(') then
               if (izai.eq.1) then
                  name(2:2)='n'
               else if (izai.eq.1001) then
                  name(2:2)='p'
               else if (izai.eq.1002) then
                  name(2:2)='d'
               else if (izai.eq.1003) then
                  name(2:2)='t'
               else if (izai.eq.2003) then
                  name(2:2)='s'
               else if (izai.eq.2004) then
                  name(2:2)='a'
               endif
            endif
            l3=dlwh+nint(xss(l2+2))-1
            ne=nint(xss(l3+1))
            zmin=1000
            zmax1=0
            zmax2=0
            ymin=xss(l3+2)
            ymax=xss(l3+2+ne-1)
            test=3
            i1=1
            if (ne.gt.2) then
               i1=2
               if (ymax.gt.test*xss(l3+2+ne-2)) ymax=xss(l3+2+ne-2)
            endif
            do ie=i1,ne
               e=xss(l3+2+ie-1)
               if (e.le.up*ymax) then
                  loci=nint(xss(l3+2+ne+ie-1))+dlwh-1
                  intt=nint(xss(loci))
                  nn=nint(xss(loci+1))
                  if (nn.ge.2) then
                     loci=loci+1
                     do j=1,nn
                        ep=xss(loci+j)
                        pd=xss(loci+nn+j)
                        if (pd.lt.zmin) zmin=pd
                        if (pd.gt.zmax2.and.pd.lt.zmax1)zmax2=pd
                        if (pd.gt.zmax1) then
                           zmax2=zmax1
                           zmax1=pd
                        endif
                     enddo
                  elseif(nn.eq.1) then
                     loci=loci+1
                     pd=xss(loci+2)
                     zmin=pd
                     zmax1=pd
                     zmax2=zmax1
                  endif
               endif
            enddo
            if (zmax2/zmax1.lt.1.e-4_kr.and.zmax2.ne.0) then
               zmax=zmax2
            else
               zmax=zmax1
            endif
            zmin=zmax/10000
            call ascll(zmin,zmax)
            xmin=1000
            xmax=0
            do ie=1,ne
               e=xss(l3+2+ie-1)
               if (e.le.up*ymax) then
                  loci=nint(xss(l3+2+ne+ie-1))+dlwh-1
                  intt=nint(xss(loci))
                  nn=nint(xss(loci+1))
                  loci=loci+1
                  do j=1,nn
                     ep=xss(loci+j)
                     pd=xss(loci+nn+j)
                     if (pd.ge.zmin.and.pd.le.zmax) then
                        if (ep.lt.xmin) xmin=ep
                        if (ep.gt.xmax) xmax=ep
                     endif
                  enddo
               endif
            enddo
            if (xmax-xmin.lt.xmax/10) then
               xmin=xmin-xmin/10
               xmax=xmax+xmax/10
            endif
            call ascle(1,xmin,xmax,major,minor)
            xstep=(xmax-xmin)/major
            call ascle(4,ymin,ymax,major,minor)
            ystep=(ymax-ymin)/major
            write(nout,'(''1'',i3,''/'')') iwcol
            it=1
            do j=1,70
               if (hk(j:j).ne.' ') it=j
            enddo
            write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
            if (ipt.eq.1) write(nout,&
               '(a,''neutrons from '',a,a,''/'')') qu,name,qu
            if (ipt.eq.9) write(nout,&
               '(a,''protons from '',a,a,''/'')') qu,name,qu
            if (ipt.eq.31) write(nout,&
               '(a,''deuterons from '',a,a,''/'')') qu,name,qu
            if (ipt.eq.32) write(nout,&
               '(a,''tritons from '',a,a,''/'')') qu,name,qu
            if (ipt.eq.33) write(nout,&
               '(a,''he3s from '',a,a,''/'')') qu,name,qu
            if (ipt.eq.34) write(nout,&
               '(a,''alphas from '',a,a,''/'')') qu,name,qu
            write(nout,'(''-1 2/'')')
            write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
            write(nout,'(a,''<s>ec. <e>nergy'',a,''/'')') qu,qu
            write(nout,'(1p,3e12.3,''/'')') ymin,ymax,ystep
            write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
            write(nout,'(1p,3e12.3,''/'')') zmin,zmax,one
            write(nout,'(a,''<p>rob/<m>e<v>'',a,''/'')') qu,qu
            write(nout,'(''/'')')
            write(nout,'(''/'')')
            write(nout,'(''1/'')')
            do ie=1,ne
               e=xss(l3+2+ie-1)
               if (e.le.up*ymax) then
                  loci=nint(xss(l3+2+ne+ie-1))+dlwh-1
                  intt=nint(xss(loci))
                  nd=nint(xss(loci))/10
                  nn=nint(xss(loci+1))
                  write(nout,'(1p,e14.6,''/'')') e
                  loci=loci+1
                  ylast=zmin
                  k=0
                  do j=1,nn
                     ep=xss(loci+j)
                     if (ep.ge.xmin.and.ep.le.xmax) then
                        pd=xss(loci+nn+j)
                        if (j.le.nd) pd=100*pd
                        if (pd.lt.zmin) pd=zmin
                        if (pd.gt.zmax) pd=zmax
                        if (j.le.nd) write(nout,&
                           '(1p,2e14.6,''/'')') ep-delta,zmin
                        if (intt.eq.1.and.j.gt.nd) then
                           write(nout,'(1p,2e14.6,''/'')') ep,ylast
                        endif
                        write(nout,'(1p,2e14.6,''/'')') ep,pd
                        if (j.le.nd) write(nout,&
                           '(1p,2e14.6,''/'')') ep+delta,zmin
                        k=k+1
                        ylast=pd
                     endif
                  enddo
                  if (k.lt.2) then
                     write(nout,'(1p,2e14.6,''/'')') xmin,zmin
                     write(nout,'(1p,2e14.6,''/'')') xmax,zmin
                  endif
                  write(nout,'(''/'')')
               endif
            enddo
            write(nout,'(''/'')')
         endif

         !--check for an angular distribution
         na=nint(xss(landh+imt-1))
         if (na.gt.0) then
            mt=iabs(nint(xss(mtrh+imt-1)))
            call mtname(mt,name,izai)
            if (name(1:1).eq.'(') then
               if (izai.eq.1) then
                  name(2:2)='n'
               else if (izai.eq.1001) then
                  name(2:2)='p'
               else if (izai.eq.1002) then
                  name(2:2)='d'
               else if (izai.eq.1003) then
                  name(2:2)='t'
               else if (izai.eq.2003) then
                  name(2:2)='s'
               else if (izai.eq.2004) then
                  name(2:2)='a'
               endif
            endif
            na=na+andh-1
            ne=nint(xss(na))
            nb=na+ne
            k=nint(xss(nb+1))
            xmin=-1
            xmax=1
            xstep=1
            xstep=xstep/2
            ymin=xss(na+1)
            ymax=xss(na+ne)
            itwo=0
            if (ne.gt.4) then
               test=3
               if (ymax.gt.test*xss(na+ne-1)) ymax=xss(na+ne-1)
               test=50
               if (ymax.ge.test) then
                  break=20
                  if (ymin.lt.dn*break.and.ymax.gt.up*break) itwo=1
               endif
            endif
            do while (itwo.ge.0.and.itwo.le.2)
               if (itwo.eq.1) then
                  ymax=break
               else if (itwo.eq.2) then
                  ymin=break
                  ymax=xss(na+ne)
               endif
               call ascle(4,ymin,ymax,major,minor)
               ystep=(ymax-ymin)/major
               zmin=big
               zmax=0
               do ie=1,ne
                  e=xss(na+ie)
                  if (e.le.(1+eps)*ymax.and.e.ge.ymin) then
                     k=nint(abs(xss(nb+ie)))+andh-1
                     np=nint(xss(k+1))
                     k=k+1
                     do j=1,np
                        cc=xss(k+j)
                        pp=xss(k+np+j)
                        if (pp.lt.zmin) zmin=pp
                        if (pp.gt.zmax) zmax=pp
                     enddo
                  endif
               enddo
               rat=100000
               if (zmin.le.zero) zmin=zmax/rat
               if (zmax/zmin.gt.rat)  zmin=zmax/rat
               call ascll(zmin,zmax)
               if (zmin.eq.zmax) then
                   zmax=2*zmax
                   zmin=zmax/10
               endif
               do ic=1,70
                  if (hk(ic:ic).ne.' ') it=ic
               enddo
               write(nout,'(''1'',i3,''/'')') iwcol
               write(nout,&
                  '(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
               l=len_trim(name)
               if (ipt.eq.1) then
                  write(nout,&
                    &'(a,''angular distribution for '',a,&
                    &'' neutron'',a,''/'')') qu,name(1:l),qu
               else if (ipt.eq.9) then
                  write(nout,&
                    &'(a,''angular distribution for '',a,&
                    &'' proton'',a,''/'')') qu,name(1:l),qu
               else if (ipt.eq.31) then
                  write(nout,&
                    &'(a,''angular distribution for '',a,&
                    &'' deuteron'',a,''/'')') qu,name(1:l),qu
               else if (ipt.eq.32) then
                  write(nout,&
                    &'(a,''angular distribution for '',a,&
                    &'' triton'',a,''/'')') qu,name(1:l),qu
               else if (ipt.eq.33) then
                  write(nout,&
                    &'(a,''angular distribution for '',a,&
                    &'' 3he'',a,''/'')') qu,name(1:l),qu
               else if (ipt.eq.34) then
                  write(nout,&
                    &'(a,''angular distribution for '',a,&
                    &'' alpha'',a,''/'')') qu,name(1:l),qu
               endif
               write(nout,'(''-1 2/'')')
               write(nout,'(1p,3e12.3,''/'')') xmin,xmax,xstep
               write(nout,'(a,''<c>osine'',a,''/'')') qu,qu
               write(nout,'(1p,3e12.3,''/'')') ymin,ymax,ystep
               write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
               write(nout,'(1p,3e12.3,''/'')') zmin,zmax,one
               write(nout,'(a,''<p>rob/<C>os'',a,''/'')') qu,qu
               write(nout,'(''/'')')
               write(nout,'('' 15. -15. 15. -2.5 6.5 2.5/'')')
               write(nout,'(''1/'')')
               stepm=(ymax-ymin)/150
               elast=0
               do ie=1,ne
                  e=xss(na+ie)
                  iflag=0
                  if (e.gt.(1+eps)*ymax) iflag=1
                  if (e.lt.ymin) iflag=1
                  if (ie.gt.1.and.ie.lt.ne.and.ne.gt.150.and.&
                    e-elast.lt.stepm) iflag=1
                  if (iflag.eq.0) then
                     write(nout,'(1p,e14.6,''/'')') e
                     k=nint(abs(xss(nb+ie)))+andh-1
                     intt=nint(xss(k))
                     np=nint(xss(k+1))
                     k=k+1
                     ylast=zmin
                     do j=1,np
                        cc=xss(k+j)
                        pp=xss(k+np+j)
                        if (pp.lt.zmin) pp=zmin
                        if (intt.eq.1) write(nout,&
                           '(1p,2e14.6,''/'')') cc,ylast
                        write(nout,'(1p,2e14.6,''/'')') cc,pp
                        ylast=pp
                     enddo
                     write(nout,'(''/'')')
                     elast=e
                  endif
               enddo
               write(nout,'(''/'')')
               if (itwo.eq.0) itwo=3
               if (itwo.eq.2) itwo=3
               if (itwo.eq.1) itwo=2
            enddo
         endif
      enddo
   enddo
   return
   end subroutine aploxp

   subroutine getl7(ep,epnext,loci,a7,f0)
   !-------------------------------------------------------------------
   ! Special routine for returning isotropic part of scattering
   ! for evaluations using law 7.  Call with ep.lt.0. to initialize.
   ! Thereafter, call with any ep and the routine will compute
   ! f0 and return epnext, the next grid value of ep.  a7 contains
   ! the dlw data for this reaction. loci is the index to the current
   ! incident energy.   Returns epnext=1.e10 when finished with the
   ! entire ep grid.
   !--***************************************************************
   use endf ! provides terp1
   ! externals
   integer::loci
   real(kr)::ep,epnext,a7(*),f0
   ! internals
   integer::intmu,nmu,locmu,imu,locimu,intep,npep,iep,inow,idone
   real(kr)::ulast,flast,u,ep1,ep2,epn,fact1,fact2,f
   real(kr),parameter::emax=1.e10_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1

   !--initialize
   if (ep.lt.zero) then
      intmu=nint(a7(loci))
      nmu=nint(a7(loci+1))
      locmu=loci+2
      epnext=emax
      do imu=1,nmu
         u=a7(locmu+imu-1)
         locimu=nint(a7(locmu+nmu+imu-1))
         intep=nint(a7(locimu))
         npep=nint(a7(locimu+1))
         ep=a7(locimu+2)
         if (ep.lt.epnext) epnext=ep
      enddo
      return
   endif

   !--normal entry
   intmu=nint(a7(loci))
   nmu=nint(a7(loci+1))
   locmu=loci+2
   epnext=emax
   f0=0
   ulast=0
   flast=0
   do imu=1,nmu
      u=a7(locmu+imu-1)
      locimu=nint(a7(locmu+nmu+imu-1))
      intep=nint(a7(locimu))
      npep=nint(a7(locimu+1))
      iep=0
      idone=0
      do while (iep.lt.npep.and.idone.eq.0)
         iep=iep+1
         inow=locimu+1+iep
         ep1=a7(inow)
         ep2=a7(inow+1)
         if (ep.ge.ep1.and.ep.lt.ep2) idone=1
      enddo
      if (ep.ge.ep2) then
         f=a7(inow+1+npep)
      else
         call terp1(ep1,a7(inow+npep),ep2,a7(inow+1+npep),ep,f,intep)
         epn=ep2
         fact1=one-one/10000
         fact2=one-2*one/10000
         if (intep.eq.1.and.ep.lt.fact2*ep2) epn=fact1*ep2
         if (epn.lt.epnext) epnext=epn
      endif
      if (intmu.eq.1.and.imu.gt.1) f0=f0+flast*(u-ulast)
      if (intmu.eq.2.and.imu.gt.1) f0=f0+(f+flast)*(u-ulast)/2
      ulast=u
      flast=f
   enddo
   return
   end subroutine getl7

   subroutine ascle(m,z1,z2,major,minor)
   !-------------------------------------------------------------------
   ! Automatic scaling routine for a linear axis.
   ! Borrowed from Los Alamos SC4020 library (with modifications).
   ! On input------------
   !   m              minimum number of major divisions desired
   !   z1,z2          min and max values of data to be plotted
   ! On output-----------
   !   z1,z2          min and max limits of axis
   !   major          number of major division on axis
   !   minor          number of minor divisions on axis
   !--***************************************************************
   ! externals
   integer::m,major,minor
   real(kr)::z1,z2
   ! internals
   integer::iflag,k,nm,n1,n2
   real(kr)::zmin,zmax,fm,zbar,z,p,tenk,test,dz,fn
   real(kr),parameter::zero=0
   real(kr),parameter::one=1
   real(kr),parameter::ten=10
   real(kr),parameter::tentho=10000

   !--check input
   !--exit if parameters are unreasonable
   zmin=z1
   zmax=z2
   if (zmax.le.zmin.or.m.le.0.or.m.gt.20) then
      major=1
      minor=0
      z1=0
      z2=2*z2
      return
   endif

   !--find minimum span per interval
   fm=m
   if (zmax.eq.zero.or.zmin.eq.zero) then
      zmax=zmax-abs(zmax)/1000000
      zmin=zmin+abs(zmin)/1000000
   else
      zbar=zmax/zmin
      if (abs(zbar)/1000.ge.one) zmin=0
      if (abs(zbar)*1000.le.one) then
         zmax=0
         zmax=zmax-abs(zmax)/1000000
         zmin=zmin+abs(zmin)/1000000
      else
         if ((abs(zbar-1)-5*fm/100000).gt.zero) then
            zmax=zmax-abs(zmax)/1000000
            zmin=zmin+abs(zmin)/1000000
         else
            zbar=(zmax+zmin)/2
            z=26*fm*abs(zbar)/1000000
            zmax=zbar+z
            zmin=zbar-z
         endif
      endif
   endif
   p=(zmax-zmin)/fm

   !--determine k such that
   !--    10**k.le.p.lt.10**(k+1)
   iflag=0
   tenk=1
   k=0
   if (p.lt.one) then
      iflag=1
      p=1/p
   endif
   do while (p.ge.tentho)
      p=p/10000
      tenk=tenk*10000
      k=k+4
   enddo
   do while (p.ge.ten)
      p=p/10
      tenk=tenk*10
      k=k+1
   enddo
   if (iflag.ne.0) then
      p=10/p
      tenk=one/10/tenk
      k=k-1
   endif

   !--determine dz such that
   !--    dz = 1*10**k, or
   !--    dz = 2*10**k, or
   !--    dz = 5*10**k.
   test=2+one/100
   if (p.le.test) then
      test=2-one/100
      if (p.le.test) then
         p=1
         nm=5
      else
         p=2
         nm=4
      endif
   else
      test=5-one/100
      if (p.lt.test) then
         p=2
         nm=4
      else
         p=5
         nm=5
      endif
   endif
   dz=p*tenk

   !--find integers n1 and n2 such that
   !--    n1*dz.le.zmin.lt.(n1+1)*dz, and
   !--    (n2-1)*dz.lt.zmax.le.n2*dz
   !--the new zmin and zmax become
   !--    zmin = n1*dz
   !--    zmax = n2*dz
   n1=int(zmin/dz)
   fn=n1
   z=fn*dz
   if (z.gt.zmin) then
      z=z-dz
      n1=n1-1
   endif
   zmin=z
   z1=zmin
   n2=int(zmax/dz)
   fn=n2
   z=fn*dz
   if (z.lt.zmax) then
      n2=n2+1
      z=z+dz
   endif
   zmax=z
   z2=zmax

   !--the values of major and minor can be set
   major=n2-n1
   minor=nm*major
   return
   end subroutine ascle

   subroutine ascll(amin,amax)
   !-------------------------------------------------------------------
   ! Adjust axis limits for a log scale.
   !-------------------------------------------------------------------
   ! externals
   real(kr)::amin,amax
   ! internals
   integer::itop,ibot
   real(kr)::top,bot
   real(kr),parameter::b1=0.e0_kr
   real(kr),parameter::b2=.301e0_kr
   real(kr),parameter::b3=.699e0_kr
   real(kr),parameter::b4=1.e0_kr
   real(kr),parameter::ohoh1=.001e0_kr

   real(kr),parameter::zero=0
   real(kr),parameter::ten=10

   top=log10(amax)
   itop=int(top)
   if (top.lt.zero) itop=itop-1
   amax=itop+1
   if (top.le.itop+b3+ohoh1) amax=itop+b3
   if (top.le.itop+b2+ohoh1) amax=itop+b2
   if (top.le.itop+b1+ohoh1) amax=itop+b1
   amax=ten**amax
   bot=log10(amin)
   ibot=int(bot)
   if (bot.lt.zero) ibot=ibot-1
   amin=ibot
   if (bot.ge.ibot+b2-ohoh1) amin=ibot+b2
   if (bot.ge.ibot+b3-ohoh1) amin=ibot+b3
   if (bot.ge.ibot+b4-ohoh1) amin=ibot+b4
   amin=ten**amin
   ! make sure amin and amax differ.
   if (amin.eq.amax) then
       amin=0.1*amin
       amax=10.*amax
   endif
   return
   end subroutine ascll

end module acefc

