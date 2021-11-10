module acepn
   ! provides photonuclear options for acer
   use locale
   use acecm, only: xss,nxss
   implicit none
   private

   !--Public routines
   public acephn,phnfix

   !--Private global variables

   ! ace file header information
   character(13)::hz
   character(10)::hd
   character(10)::hm
   real(kr)::aw0,tz

   ! parameters for photonuclear nxs block
   integer::lxs,za,nes,ntr,ntype,npixs,neixs,nxsd(8),tvn

   ! parameters for photonuclear jxs block
   integer::esz,tot,non,els,thn,mtr,lqr,lsig,sig,ixsa,ixs,jxsd(21)

contains

   subroutine acephn(nendf,npend,nace,ndir,matd,tempd,iprint,mcnpx,&
     ityp,suff,hk,izn,awn)
   !-------------------------------------------------------------------
   ! Prepare ACE photo-nuclear files.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use physics ! provides amassn,amu,ev,pi,clight
   use util  ! provides openz,closz,repoz,mess,error
   use mathm ! provides erfc
   use endf  ! provides endf routines and variables
   use acecm ! provides bachaa,eavl,ptleg2,pttab2
   ! externals
   integer::nendf,npend,nace,ndir,matd,iprint,mcnpx,ityp
   real(kr)::tempd,suff
   integer::izn(16)
   real(kr)::awn(16)
   character(70)::hk
   ! internals
   integer::nin,nb,nw,nwscr,nx,mtx,ielas,mf4,mf6,mt452,mt456,mtxnu
   integer::mt103,mt104,mt105,mt106,mt107
   integer::i,mfd,mtd,l,mttot,idis,nex,nexc,ir,j,idone,nnex,n
   integer::nneut,nphot,nprot,ndeut,ntrit,nhe3,nhe4
   integer::k,ia,iaa,nk,ik,lly,izai,izap,law,jscr,nrr,npp
   integer::ll,lll,lep,ne,llh,lld,ie,np,ip,mtt,lct,ii
   integer::icapt,jj,itype,it,jp,nr,il,llht,iie,lang,lleg,ileg
   integer::iint,nn,kk,m,intt,last,lf,jnt,ja,jb,ipp,irr
   integer::lee,lle,nd,na,ncyc,ng,ig,nnr,nnp,mf,mt
   integer::ipt,ntrp,pxs,phn,mtrp,tyrp,lsigp,sigp,landp,andp,ldlwp,dlwp
   integer::izarec,nl,iil,nexn,nexd,ki
   integer::nle
   integer::imu,intmu,nmu
   real(kr)::emc2,e,enext,s,y,ynext,heat,en,ep,g,h,epl
   real(kr)::tneut,tphot,tprot,tdeut,ttrit,the3,the4,thresh
   real(kr)::ss,tt,ubar,sum,renorm,ebar,hh,u,theta,x,anorm
   real(kr)::ee,amass,avadd,avlab,avll,test,rkal,akal
   real(kr)::eavi,avl,avcm,sign,dele,avav,zaid,gl,awp,awr,q
   real(kr)::av,del
   real(kr)::awprec,awpp
   integer,parameter::mmax=1000
   integer::mfm(mmax),mtm(mmax),nr6(mmax)
   real(kr)::fnubar(300)
   character(8)::hdt
   real(kr),dimension(:),allocatable::scr
   real(kr),parameter::emev=1.e6_kr
   real(kr),parameter::etop=1.e10_kr
   real(kr),parameter::small=1.e-12_kr
   real(kr),parameter::rmin=1.e-30_kr
   real(kr),parameter::eps=1.e-10_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1
   integer,parameter::ni=64
   character(66)::text
   emc2=amassn*amu*clight*clight/ev/emev
   tvn=1

   nxsd=0
   jxsd=0
   xss=0

   !--allocate scratch storage
   nwscr=50000
   allocate(scr(nwscr))

   !--assign input file
   nin=nendf
   call openz(nin,0)
   call openz(npend,0)

   !--determine what endf version is being used
   call repoz(nin)
   call tpidio(nin,0,0,scr,nb,nw)
   call contio(nin,0,0,scr,nb,nw)
   awr=c2h
   call contio(nin,0,0,scr,nb,nw)
   if (n1h.ne.0) then
      iverf=4
   else if (n2h.eq.0) then
      iverf=5
   else
      iverf=6
   endif
   write(nsyso,'(/'' using endf-'',i1,'' format'')') iverf

   !--process dictionary from endf tape.
   call findf(matd,1,451,nin)
   i=1
   call contio(nin,0,0,scr,nb,nw)
   nx=n2h
   izai=0
   if (iverf.ge.5) then
      i=i+6
      call contio(nin,0,0,scr(i),nb,nw)
   endif
   if (iverf.eq.6) then
      i=i+6
      call contio(nin,0,0,scr(i),nb,nw)
      izai=nint(scr(17))/10
      if (izai.ne.0) call error('acephn',&
         'endf input tape is not a photonuclear tape',' ')
   endif
   call hdatio(nin,0,0,scr,nb,nw)
   if (iverf.ge.5) nx=n2h
   do while (nb.ne.0)
      call moreio(nin,0,0,scr,nb,nw)
   enddo
   nw=nx
   ntr=0
   mtx=0
   ielas=0
   mf4=0
   mf6=0
   mt103=0
   mt104=0
   mt105=0
   mt106=0
   mt107=0
   mt452=0
   mt456=0
   call dictio(nin,0,0,scr,nb,nw)
   do i=1,nw,6
      mfd=nint(scr(i+2))
      mtd=nint(scr(i+3))
      if (mfd.eq.1.and.mtd.eq.452) mt452=1
      if (mfd.eq.1.and.mtd.eq.456) mt456=1
      if (mfd.ge.3.and.mfd.lt.30.and.(mtd.eq.2.or.mtd.gt.4)) then
         if (mfd.eq.3) ntr=ntr+1
         if (mtd.eq.2) ielas=1
         if (mfd.eq.3.and.(mtd.ge.600.and.mtd.le.649)) mt103=1
         if (mfd.eq.3.and.(mtd.ge.650.and.mtd.le.699)) mt104=1
         if (mfd.eq.3.and.(mtd.ge.700.and.mtd.le.749)) mt105=1
         if (mfd.eq.3.and.(mtd.ge.750.and.mtd.le.799)) mt106=1
         if (mfd.eq.3.and.(mtd.ge.800.and.mtd.le.849)) mt107=1
         mtx=mtx+1
         if (mtx.gt.mmax) call error('acephn',&
           'too many reactions in mtr list',' ')
         if (mfd.eq.6.and.mtd.ge.201.and.mtd.le.207)&
           call error('acephn','mf=6/mt=201-207 not supported.',&
           'does not conform to endf format.')
         mfm(mtx)=mfd
         mtm(mtx)=mtd
         nr6(mtx)=0
         if (mfd.eq.4) mf4=1
         if (mfd.eq.6) mf6=1
      endif
   enddo

   !--save the nubar tabulation, if present
   mtxnu=0
   if (mt452.eq.1) then
      call findf(matd,1,452,nin)
      call contio(nin,0,0,scr,nb,nw)
      ! if polynomial representation is used, we will need to linearise
      ! for now: error out and wait for this to come up to actually implement it
      if (scr(4).eq.1) then
        call error('acephn','mf=1/mt=452 uses polynomial representation.',&
        'this is currently unsupported for photonuclear ACE files.')
      endif
      call tab1io(nin,0,0,fnubar,nb,nw)
      mtxnu=452
   endif
   if (mt456.eq.1) then
      call findf(matd,1,456,nin)
      call contio(nin,0,0,scr,nb,nw)
      ! if polynomial representation is used, we will need to linearise
      ! for now: error out and wait for this to come up to actually implement it
      if (scr(4).eq.1) then
        call error('acephn','mf=1/mt=456 uses polynomial representation.',&
        'this is currently unsupported for photonuclear ACE files.')
      endif
      call tab1io(nin,0,0,fnubar,nb,nw)
      mtxnu=456
   endif

   !--locate and store energy grid of total cross section
   !--using the pendf input
   nin=npend
   l=0
   call findf(matd,3,0,nin)
   call contio(nin,0,0,scr,nb,nw)
   mttot=mth
   za=nint(scr(1))
   e=0
   call gety1(e,enext,idis,s,nin,scr)
   do while (enext.lt.etop)
      e=enext
      call gety1(e,enext,idis,s,nin,scr)
      if (l.eq.0.and.s.ne.zero) then
         l=l+1
         xss(l)=sigfig(e,7,-1)
      endif
      l=l+1
      xss(l)=e
   enddo
   call tosend(nin,0,0,scr)
   nes=l

   !--mttot is 5 for a lanl photonuclear files, which
   !--use only uses mt=5 in file 3 and file 6.
   !--mttot is 1 or 3 for a photonuclear file that starts
   !--with a normal total or nonelastic cross sections and
   !--that can contain partial reaction cross sections.

   !--define cross section locators
   esz=1
   tot=esz+nes
   if (ielas.eq.1) then
      non=tot+nes
      els=non+nes
      thn=els+nes
   else
      non=tot
      els=0
      thn=tot+nes
   endif
   do i=1,nes
      xss(thn+i)=0
   enddo
   nex=thn+nes

   !--assign locators for the partial cross sections
   mtr=nex
   lqr=mtr+ntr
   lsig=lqr+ntr
   sig=lsig+ntr
   nex=sig
   ir=0

   !--read the reaction cross sections from file 3
   if (mttot.eq.5) call findf(matd,3,5,nin)
   mfh=3
   do while (mfh.ne.0)
      call contio(nin,0,0,scr,nb,nw)
      if (mfh.ne.0) then
         if (mth.eq.2.or.mth.gt.4) then
            if (mth.eq.103.and.mt103.ne.0) go to 99
            if (mth.eq.104.and.mt104.ne.0) go to 99
            if (mth.eq.105.and.mt105.ne.0) go to 99
            if (mth.eq.106.and.mt106.ne.0) go to 99
            if (mth.eq.107.and.mt107.ne.0) go to 99
            e=0
            call gety1(e,enext,idis,s,nin,scr)
            e=enext
            call gety1(e,enext,idis,s,nin,scr)
            if (s.ne.zero) e=sigfig(e,7,-1)

            !--locate energy index for reaction threshold
            j=nes
            idone=0
            do while (idone.eq.0)
               if (j.ge.1) then
                  if (e.le.(1+eps)*xss(esz+j-1)) then
                     j=j-1
                  else
                     idone=1
                  endif
               else
                  idone=1
               endif
            enddo
            j=j+1

            !--store reaction parameters
            ir=ir+1
            xss(mtr+ir-1)=mth
            xss(lqr+ir-1)=sigfig(scr(2)/emev,7,0)
            xss(lsig+ir-1)=nex-sig+1
            xss(nex)=j
            nnex=nex+1
            nex=nex+2
            n=0

            !--store cross sections
            !--compute total from the partials
            do i=j,nes
               e=xss(i)
               call gety1(e,enext,idis,s,nin,scr)
               s=sigfig(s,7,0)
               xss(tot+i-1)=xss(tot+i-1)+s
               if (mth.eq.2) xss(els+i-1)=s
               xss(nex)=s
               n=n+1
               nex=nex+1
            enddo
            xss(nnex)=n
   99       continue
         endif
         call tosend(nin,0,0,scr)
      endif
   enddo

   !--make one pass through the distributions to
   !--count the different particles produced,
   !--to determine the production thresholds, and
   !--to accumulate the heating from mf=6 recoils.
   nin=nendf
   nneut=0
   nphot=0
   nprot=0
   ndeut=0
   ntrit=0
   nhe3=0
   nhe4=0
   tneut=etop
   tphot=etop
   tprot=etop
   tdeut=etop
   ttrit=etop
   the3=etop
   the4=etop
   do while (math.gt.0)
      call contio(nin,0,0,scr,nb,nw)
      if (math.gt.0.and.mfh.ne.0) then

         !--file 4
         if (mfh.eq.4) then
            ! file 4 is only to be used for secondary neutrons so if a reaction
            ! is present in mf4, it describes secondary neutrons so every
            ! reaction is counted
            nneut=nneut+1
            mtt=0
            ir=0
            do while (mtt.ne.mth)
               ir=ir+1
               mtt=nint(xss(mtr+ir-1))
               q=xss(lqr+ir-1)
               k=nint(xss(lsig+ir-1))+sig-1
               n=nint(xss(k+1))
               iaa=nint(xss(k))
            enddo
            thresh=xss(esz+iaa-1)
            if (thresh.lt.tneut) tneut=thresh

         !--file 6
         else if (mfh.eq.6) then
            lct=nint(scr(4))
            nk=nint(scr(5))
            ik=0
            ! as long as there are reaction products in the mf6 entry
            do while (ik.lt.nk)
               ik=ik+1
               lly=1
               ! read the multiplicity
               call tab1io(nin,0,0,scr,nb,nw)
               jscr=1+nw
               do while (nb.ne.0)
                  call moreio(nin,0,0,scr(jscr),nb,nw)
                  jscr=jscr+nw
               enddo
               ! retrieve izap and the law
               izap=nint(scr(1))
               law=nint(scr(4))
               ! count particle producing reactions
               if (izap.eq.1) nneut=nneut+1
               if (izap.eq.0) nphot=nphot+1
               if (izap.eq.1001) nprot=nprot+1
               if (izap.eq.1002) ndeut=ndeut+1
               if (izap.eq.1003) ntrit=ntrit+1
               if (izap.eq.2003) nhe3=nhe3+1
               if (izap.eq.2004) nhe4=nhe4+1
               ! if this is fission, check if the multiplicity is equal to nubar
               ! issue a warning if this is not the case and replace the yield
               if (mth.eq.18.and.izap.eq.1.and.mtxnu.gt.0) then
                  ! check yield != nubar
                  if (scr(6+2*nint(scr(5))+2).ne.fnubar(6+2*nint(fnubar(5))+2)) then
                    write(text,'(''the multiplicity will be replaced with nubar from mf=1/mt='',i3,''.'')')mtxnu
                    call mess('acephn','mf=6/mt=18 neutron multiplicity not consistent with nubar.',text)
                  endif
                  call copynubar(scr,fnubar,jscr)
               endif
               mtt=0
               ir=0
               ! look for the corresponding reaction in the XSS array
               do while (mtt.ne.mth)
                  ir=ir+1
                  mtt=nint(xss(mtr+ir-1))
                  q=xss(lqr+ir-1)
                  k=nint(xss(lsig+ir-1))+sig-1
                  n=nint(xss(k+1))
                  iaa=nint(xss(k))
               enddo
               thresh=xss(esz+iaa-1)

               !--for particles
               !--check for production thresholds
               if (izap.le.2004) then
                  nrr=nint(scr(5))
                  npp=nint(scr(6))
                  y=0
                  i=0
                  ynext=0
                  do while (ynext.le.zero)
                     i=i+1
                     e=scr(5+2*nrr+2*i)
                     ynext=scr(6+2*nrr+2*i+2)
                  enddo
                  if (e.lt.thresh) e=thresh
                  if (izap.eq.1.and.e.lt.tneut) tneut=e
                  if (izap.eq.0.and.e.lt.tphot) tphot=e
                  if (izap.eq.1001.and.e.lt.tprot) tprot=e
                  if (izap.eq.1002.and.e.lt.tdeut) tdeut=e
                  if (izap.eq.1003.and.e.lt.ttrit) ttrit=e
                  if (izap.eq.2003.and.e.lt.the3) the3=e
                  if (izap.eq.2004.and.e.lt.the4) the4=e

                  !--skip to the next subsection
                  call skip6(nin,0,0,scr,law)

               !--for recoil nuclides
               !--compute the heating contributions
               !--for law=1 continuous distributions
               else if (law.eq.1) then
                  ll=jscr
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
                     heat=0
                     na=nint(scr(lld+3))
                     np=nint(scr(lld+5))
                     call terpa(y,e,en,idis,scr(lly),npp,nrr)
                     do ip=1,np
                        ep=scr(lld+6+(na+2)*(ip-1))
                        g=scr(lld+7+(na+2)*(ip-1))
                        if (ip.gt.1) then
                           heat=heat+(ep-epl)*gl*(ep+epl)/2
                        endif
                        epl=ep
                        gl=g
                     enddo
                     scr(llh+6+2*ie)=e
                     scr(llh+7+2*ie)=y*heat
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
                     e=xss(esz+ie-1)
                     call terpa(h,e,en,idis,scr(llh),npp,nrr)
                     h=(h/emev)*xss(2+k+ie-iaa)
                     if (xss(tot+ie-1).ne.zero) h=h/xss(tot+ie-1)
                     xss(thn+ie-1)=sigfig(xss(thn+ie-1)+h,7,0)
                  enddo
                  do ii=1,mtx
                     if (mfm(ii).eq.6.and.mtm(ii).eq.mth)&
                       nr6(ii)=nr6(ii)+1
                  enddo

               !--heating from discrete recoils
               !--back up to corresponding law=2 distribution
               else if (law.eq.4) then
                  mf=mfh
                  mt=mth
                  call findf(matd,mf,mt,nin)
                  call contio(nin,0,0,scr,nb,nw)
                  lly=1
                  call tab1io(nin,0,0,scr,nb,nw)
                  izap=nint(scr(1))
                  awp=scr(2)
                  law=nint(scr(4))
                  jscr=1+nw
                  do while (nb.ne.0)
                     call moreio(nin,0,0,scr(jscr),nb,nw)
                     jscr=jscr+nw
                  enddo
                  ll=jscr
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
                     scr(llh+6+2*ie)=e
                     scr(llh+7+2*ie)=awp*(e+q*emev)/(awr-awp)
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
                     e=xss(esz+ie-1)
                     call terpa(h,e,en,idis,scr(llh),npp,nrr)
                     h=(h/emev)*xss(2+k+ie-iaa)
                     if (xss(tot+ie-1).ne.zero) h=h/xss(tot+ie-1)
                     xss(thn+ie-1)=sigfig(xss(thn+ie-1)+h,7,0)
                  enddo
                  do ii=1,mtx
                     if (mfm(ii).eq.6.and.mtm(ii).eq.mth)&
                       nr6(ii)=nr6(ii)+1
                  enddo

                  ! read subsections until we're back at the right one
                  ir=1
                  do while (ir.lt.ik)
                     ir=ir+1
                     call tab1io(nin,0,0,scr,nb,nw)
                     law=nint(scr(4))
                     jscr=1+nw
                     do while (nb.ne.0)
                        call moreio(nin,0,0,scr(jscr),nb,nw)
                        jscr=jscr+nw
                     enddo
                     call skip6(nin,0,0,scr,law)
                  enddo

               !--unknown distribution
               else if (law.eq.0) then
                  write(text,'(''recoil'',i6,'' in MT'',I4)')izap,mth
                  call mess('acephn','no heating info for ',text)

               !--this law is not currently handled
               else
                  write(text,'(''particle '',i5,'' law'',I4)')izap,law
                  call mess('acephn','file 6 law not coded for ',text)
               endif
            enddo
         endif
         call tosend(nin,0,0,scr)
      endif
   enddo

   !--add in e+q heating for capture when recoil is not given.
   icapt=0
   do ii=1,mtx
      if (mfm(ii).eq.3.and.mtm(ii).eq.102) icapt=1
      if (mfm(ii).eq.6.and.mtm(ii).eq.102.and.nr6(ii).gt.0) icapt=0
   enddo
   do ii=1,ntr
      mt=int(xss(mtr+ii-1))
      if (mt.eq.102.and.icapt.eq.1) then
         k=nint(xss(lsig+ii-1))+sig-1
         j=nint(xss(k))
         do ie=j,nes
            e=xss(esz+ie-1)
            if (xss(tot+ie-1).ne.zero)&
              xss(thn+ie-1)=xss(thn+ie-1)&
              +(e/emev)*xss(2+k+ie-j)/xss(tot+ie-1)
         enddo
      endif
   enddo

   !--add in e+q heating for mf=6 sections with no recoil
   do ii=1,ntr
      mt=int(xss(mtr+ii-1))
      do jj=1,mtx
         if (mfm(jj).eq.6.and.mtm(jj).eq.mt.and.nr6(jj).eq.0) then
            q=xss(lqr+ii-1)
            k=nint(xss(lsig+ii-1))+sig-1
            j=nint(xss(k))
            do ie=j,nes
               e=xss(esz+ie-1)
               if (xss(tot+ie-1).ne.zero)&
                 xss(thn+ie-1)=xss(thn+ie-1)&
                 +(e/emev+q)*xss(2+k+ie-j)/xss(tot+ie-1)
            enddo
         endif
      enddo
   enddo

   !--fill in the ixs array block
   ntype=0
   npixs=2
   neixs=12
   if (nneut.gt.0) ntype=ntype+1
   if (nphot.gt.0) ntype=ntype+1
   if (nprot.gt.0) ntype=ntype+1
   if (ndeut.gt.0) ntype=ntype+1
   if (ntrit.gt.0) ntype=ntype+1
   if (nhe3.gt.0) ntype=ntype+1
   if (nhe4.gt.0) ntype=ntype+1
   ixsa=nex
   itype=0
   if (nneut.gt.0) then
      itype=itype+1
      xss(ixsa+neixs*(itype-1))=1
      xss(ixsa+neixs*(itype-1)+1)=nneut
   endif
   if (nphot.gt.0) then
      itype=itype+1
      xss(ixsa+neixs*(itype-1))=2
      xss(ixsa+neixs*(itype-1)+1)=nphot
   endif
   if (nprot.gt.0) then
      itype=itype+1
      xss(ixsa+neixs*(itype-1))=9
      xss(ixsa+neixs*(itype-1)+1)=nprot
   endif
   if (ndeut.gt.0) then
      itype=itype+1
      xss(ixsa+neixs*(itype-1))=31
      xss(ixsa+neixs*(itype-1)+1)=ndeut
   endif
   if (ntrit.gt.0) then
      itype=itype+1
      xss(ixsa+neixs*(itype-1))=32
      xss(ixsa+neixs*(itype-1)+1)=ntrit
   endif
   if (nhe3.gt.0) then
      itype=itype+1
      xss(ixsa+neixs*(itype-1))=33
      xss(ixsa+neixs*(itype-1)+1)=nhe3
   endif
   if (nhe4.gt.0) then
      itype=itype+1
      xss(ixsa+neixs*(itype-1))=34
      xss(ixsa+neixs*(itype-1)+1)=nhe4
   endif
   nex=ixsa+neixs*ntype
   ixs=nex

   !--loop over each of the ntype productions
   !--to build the production data.
   do itype=1,ntype
      ipt=nint(xss(ixsa+neixs*(itype-1)))
      ntrp=nint(xss(ixsa+neixs*(itype-1)+1))
      if (ipt.eq.1) ip=1
      if (ipt.eq.2) ip=0
      if (ipt.eq.9) ip=1001
      if (ipt.eq.31) ip=1002
      if (ipt.eq.32) ip=1003
      if (ipt.eq.33) ip=2003
      if (ipt.eq.34) ip=2004

      !--determine the threshold
      !--and assign some locators
      if (ip.eq.1) thresh=tneut
      if (ip.eq.0) thresh=tphot
      if (ip.eq.1001) thresh=tprot
      if (ip.eq.1002) thresh=tdeut
      if (ip.eq.1003) thresh=ttrit
      if (ip.eq.2003) thresh=the3
      if (ip.eq.2004) thresh=the4
      it=1
      do while (xss(esz+it-1).lt.thresh)
         it=it+1
      enddo
      if (it.gt.1) it=it-1
      pxs=nex
      xss(ixsa+neixs*(itype-1)+2)=pxs
      xss(pxs)=it
      xss(pxs+1)=nes-it+1
      phn=pxs+2+nes-it+1
      xss(ixsa+neixs*(itype-1)+3)=phn
      xss(phn)=it
      xss(phn+1)=nes-it+1
      mtrp=phn+2+nes-it+1
      xss(ixsa+neixs*(itype-1)+4)=mtrp
      tyrp=mtrp+ntrp
      xss(ixsa+neixs*(itype-1)+5)=tyrp
      lsigp=tyrp+ntrp
      xss(ixsa+neixs*(itype-1)+6)=lsigp
      sigp=lsigp+ntrp
      xss(ixsa+neixs*(itype-1)+7)=sigp
      nex=sigp

      !--find each contribution to this production.
      call repoz(nin)
      jp=0

      !--here for mf4/5 representations - i.e. neutrons only
      if (mf4.eq.1) then
         do i=1,mtx
            if (mfm(i).eq.4.and.ip.eq.1) then
               jp=jp+1
               mt=mtm(i)
               mtt=0
               ir=0
               do while (mtt.ne.mt)
                  ir=ir+1
                  mtt=nint(xss(mtr+ir-1))
                  q=xss(lqr+ir-1)
                  k=nint(xss(lsig+ir-1))+sig-1
                  n=nint(xss(k+1))
                  iaa=nint(xss(k))
               enddo
               xss(mtrp+jp-1)=mt
               xss(tyrp+jp-1)=1
               xss(lsigp+jp-1)=nex-sigp+1
               if (mt.eq.18) then
                  nrr=1
                  npp=2
                  do j=iaa,nes
                     e=xss(esz+j-1)
                     call terpa(y,e,en,idis,fnubar,npp,nrr)
                     ss=xss(2+k+j-iaa)
                     tt=xss(pxs+2+j-it)+y*ss
                     xss(pxs+2+j-it)=sigfig(tt,7,0)
                     if (xss(tot+j-1).ne.zero)&
                       xss(thn+j-1)=xss(thn+j-1)&
                       +ss*(e/emev+q)/xss(tot+j-1)
                  enddo
                  nr=nint(fnubar(5))
                  ne=nint(fnubar(6))
                  xss(nex)=12
                  xss(nex+1)=mt
                  xss(nex+2)=0
                  xss(nex+3)=ne
                  do j=1,ne
                     xss(nex+3+j)=sigfig(fnubar(5+2*nr+2*j)/emev,7,0)
                     xss(nex+3+ne+j)=sigfig(fnubar(6+2*nr+2*j),7,0)
                  enddo
                  nex=nex+4+2*ne
               else
                  y=1
                  if (mth.eq.11.or.mth.eq.16.or.mth.eq.24.or.&
                      mth.eq.30.or.mth.eq.41.or.mth.eq.154.or.&
                      mth.eq.159.or.mth.eq.176.or.mth.eq.190.or.&
                      (mth.ge.875.and.mth.le.891)) then
                     y=2
                  elseif (mth.eq.17.or.mth.eq.25.or.mth.eq.42.or.&
                          mth.eq.157.or.mth.eq.172.or.mth.eq.177.or.&
                          mth.eq.179.or.(mth.ge.181.and.mth.le.199)) then
                     y=3
                  elseif (mth.eq.37.or.mth.eq.156.or.mth.eq.165.or.&
                          mth.eq.169.or.mth.eq.173.or.mth.eq.178.or.&
                          (mth.ge.194.and.mth.le.196)) then
                     y=4
                  elseif (mth.eq.152.or.mth.eq.162.or.mth.eq.166.or.&
                          mth.eq.170.or.mth.eq.174.or.mth.eq.200) then
                     y=5
                  elseif (mth.eq.153.or.mth.eq.163.or.mth.eq.167.or.&
                          mth.eq.171.or.mth.eq.175) then
                     y=6
                  elseif (mth.eq.160.or.mth.eq.164.or.mth.eq.168) then
                     y=7
                  elseif (mth.eq.161) then
                     y=8
                  endif
                  do j=iaa,nes
                     e=xss(esz+j-1)
                     ss=xss(2+k+j-iaa)
                     tt=xss(pxs+2+j-it)+y*ss
                     xss(pxs+2+j-it)=sigfig(tt,7,0)
                     if (xss(tot+j-1).ne.zero)&
                       xss(thn+j-1)=xss(thn+j-1)&
                       +ss*(e/emev+q)/xss(tot+j-1)
                  enddo
                  xss(nex)=12
                  xss(nex+1)=mt
                  xss(nex+2)=0
                  xss(nex+3)=2
                  xss(nex+4)=sigfig(xss(esz+iaa-1)/emev,7,0)
                  xss(nex+5)=sigfig(xss(esz+nes-1)/emev,7,0)
                  xss(nex+6)=y
                  xss(nex+7)=y
                  nex=nex+8
               endif
            endif
         enddo
      endif

      !--here for mf6 representations
      if (mf6.ne.0) then
         call findf(matd,6,0,nin)
         mfh=6
         do while (mfh.gt.0)
            call contio(nin,0,0,scr,nb,nw)
            if (mfh.gt.0) then
               lct=nint(scr(4))
               nk=nint(scr(5))
               mtt=0
               ir=0
               do while (mtt.ne.mth)
                  ir=ir+1
                  mtt=nint(xss(mtr+ir-1))
                  q=xss(lqr+ir-1)
                  k=nint(xss(lsig+ir-1))+sig-1
                  n=nint(xss(k+1))
                  iaa=nint(xss(k))
               enddo
               ik=0
               do while (ik.lt.nk)
                  ik=ik+1
                  ! read the multiplicity
                  call tab1io(nin,0,0,scr,nb,nw)
                  jscr=1+nw
                  do while (nb.ne.0)
                     call moreio(nin,0,0,scr(jscr),nb,nw)
                     jscr=jscr+nw
                  enddo
                  ! retrieve izap and the law
                  izap=nint(scr(1))
                  law=nint(scr(4))

                  !--find the desired particle
                  if (izap.eq.ip) then

                     ! if this is fission, replace the yield with the nubar - as before
                     if (mth.eq.18.and.izap.eq.1.and.mtxnu.gt.0) then
                        call copynubar(scr,fnubar,jscr)
                     endif

                     jp=jp+1
                     xss(mtrp+jp-1)=mth
                     xss(lsigp+jp-1)=nex-sigp+1
                     if (lct.eq.1) then
                        xss(tyrp+jp-1)=1
                     else
                        xss(tyrp+jp-1)=-1
                     endif

                     !--accumulate yield times cross section
                     nrr=1
                     npp=2
                     do i=iaa,nes
                        e=xss(esz+i-1)
                        call terpa(y,e,en,idis,scr,npp,nrr)
                        ss=xss(2+k+i-iaa)
                        tt=xss(pxs+2+i-it)+y*ss
                        xss(pxs+2+i-it)=sigfig(tt,7,0)
                     enddo

                     ! the next piece of code assumes the yield is given
                     ! using one lin-lin interpolation region
                     ! for now: error out and wait for this to come up to actually implement it
                     nr=nint(scr(5))
                     if (nr.gt.1.and.nint(scr(8)).ne.2) then
                        write(text,'(''no linearised multiplicity for izap='',i4,'' in mf=6/mt='',i3,''.'')')izap,mth
                        call mess('acephn',text,'this is currently unsupported for photonuclear ACE files.')
                     endif

                     !--store the yield
                     xss(nex)=6
                     xss(nex+1)=mth
                     xss(nex+2)=0
                     ne=nint(scr(6))
                     xss(nex+3)=ne
                     do i=1,ne
                        xss(nex+3+i)=sigfig(scr(7+2*i)/emev,7,0)
                        xss(nex+3+i+ne)=sigfig(scr(8+2*i),7,0)
                     enddo
                     nex=nex+4+2*ne
                  endif
                  call skip6(nin,0,0,scr,law)
               enddo
            endif
         enddo
      endif

      !--assign the pointers for the angular data
      landp=nex
      xss(ixsa+neixs*(itype-1)+8)=landp
      andp=landp+ntrp
      xss(ixsa+neixs*(itype-1)+9)=andp
      nex=andp

      !--angular distributions
      do i=1,ntrp
         mt=nint(xss(mtrp+i-1))
         mtt=0
         ir=0
         do while (mtt.ne.mt)
            ir=ir+1
            mtt=nint(xss(mtr+ir-1))
            q=xss(lqr+ir-1)
            k=nint(xss(lsig+ir-1))+sig-1
            n=nint(xss(k+1))
            iaa=nint(xss(k))
         enddo
         jp=0
         do j=1,mtx

            !--for file 4,
            !--only neutron producing reactions are
            !--recognized, because file 5 is only
            !--able to represent neutron emission.
            !--just assume isotropy for these.
            if (mfm(j).eq.4.and.mtm(j).eq.mt) then
               jp=i
               xss(landp+jp-1)=0

            !--for file 6,
            !--check for subsections for this production
            !--with law=2.
            else if (mfm(j).eq.6.and.mtm(j).eq.mt) then
               call findf(matd,6,mt,nin)
               call contio(nin,0,0,scr,nb,nw)
               nk=nint(scr(5))
               ik=0
               do while (ik.lt.nk)
                  ik=ik+1
                  call tab1io(nin,0,0,scr,nb,nw)
                  izap=nint(scr(1))
                  awp=scr(2)
                  law=nint(scr(4))
                  jscr=1+nw
                  do while (nb.ne.0)
                     call moreio(nin,0,0,scr(jscr),nb,nw)
                     jscr=jscr+nw
                  enddo

                  !--special steps for two-body recoil
                  !--back up to the corresponding law=2 distr.
                  izarec=-1
                  awprec=-1
                  if (izap.eq.ip.and.law.eq.4) then
                     izarec=izap
                     awprec=awp
                     mf=6
                     call findf(matd,mf,mt,nin)
                     call contio(nin,0,0,scr,nb,nw)
                     call tab1io(nin,0,0,scr,nb,nw)
                     izap=nint(c1h)
                     awp=c2h
                     law=l2h
                     jscr=1+nw
                     do while (nb.ne.0)
                        call moreio(nin,0,0,scr(jscr),nb,nw)
                        jscr=jscr+nw
                     enddo
                  endif

                  !--law2 angular distribution
                  !--also used for law 4 two-body recoils
                  if ((izap.eq.ip.and.law.eq.2).or.&
                    (izarec.eq.ip.and.law.eq.2)) then
                     lld=jscr
                     jp=i
                     ll=jscr
                     call tab2io(nin,0,0,scr(ll),nb,nw)
                     xss(landp+jp-1)=nex-andp+1
                     ne=nint(scr(ll+5))
                     xss(nex)=ne
                     ie=nex
                     il=ie+ne          ! index before the locators for the outgoing distributions
                     nex=il+ne+1       ! first outgoing distribution
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
                     if (izarec.eq.-1) then
                        awpp=awp
                     else
                        awpp=awprec
                     endif
                     do iie=1,ne
                        ll=lld
                        call listio(nin,0,0,scr(ll),nb,nw)
                        lang=nint(scr(lld+2))
                        if (lang.eq.0) then
                           if (izarec.ne.-1) then
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
                           scr(lld+6)=nn
                           scr(lld+7)=iint
                           call pttab2(scr(lld))
                        endif
                        xss(ie+iie)=sigfig(scr(lld+1)/emev,7,0) ! E(iie)
                        m=nint(scr(lld+4))
                        n=nint(scr(lld+5))
                        xss(il+iie)=nex-andp+1                  ! L(iie)
                        xss(il+iie)=-xss(il+iie)
                        intt=nint(scr(lld+7))
                        xss(nex)=intt                           ! intt
                        xss(nex+1)=n                            ! number outgoing energy values
                        if (nex+2+3*n.gt.nxss) call error('acephn',&
                          'insufficient storage for',&
                          ' angular distributions.')
                        do ii=1,n
                           xss(nex+1+ii)=&
                             sigfig(scr(lld+4+2*m+2*ii),7,0)    ! Eout(ii)
                           xss(nex+1+n+ii)=&
                             sigfig(scr(lld+5+2*m+2*ii),7,0)    ! PDF(ii)
                           if (xss(nex+1+n+ii).lt.rmin)&
                             xss(nex+1+n+ii)=0                  ! PDF(ii)
                           if (ii.eq.1) then
                              xss(nex+1+2*n+ii)=0               ! CDF(1)
                              ubar=0
                           endif
                           if (ii.gt.1.and.intt.eq.1) then
                              sum=xss(nex+1+2*n+ii-1)&
                                +xss(nex+1+n+ii-1)&
                                *(xss(nex+1+ii)-xss(nex+1+ii-1))
                              xss(nex+1+2*n+ii)=sum             ! CDF(ii)
                              ubar=ubar&
                                +xss(nex+1+n+ii-1)&
                                *(xss(nex+1+ii)-xss(nex+1+ii-1))&
                                +(xss(nex+1+ii)+xss(nex+1+ii-1))/2
                           endif
                           if (ii.gt.1.and.intt.eq.2) then
                              sum=xss(nex+1+2*n+ii-1)&
                               +(xss(nex+1+n+ii)+xss(nex+1+n+ii-1))&
                               *(xss(nex+1+ii)-xss(nex+1+ii-1))/2
                              xss(nex+1+2*n+ii)=sum             ! CDF(ii)
                              ubar=ubar&
                                +(xss(nex+1+n+ii)+xss(nex+1+n+ii-1))&
                                *(xss(nex+1+ii)-xss(nex+1+ii-1))&
                                *(xss(nex+1+ii)+xss(nex+1+ii-1))/4
                           endif
                        enddo
                        ! renormalize cummulative probabilities
                        renorm=one/xss(nex+1+3*n)
                        do ii=1,n
                           xss(nex+1+n+ii)=&
                             sigfig(renorm*xss(nex+1+n+ii),7,0)    ! PDF(ii)
                           xss(nex+1+2*n+ii)=&
                             sigfig(renorm*xss(nex+1+2*n+ii),9,0)  ! CDF(ii)
                        enddo
                        nex=nex+2+3*n       ! index for the next distribution
                        e=xss(ie+iie)
                        scr(llht+6+2*iie)=e
                        scr(llht+7+2*iie)=(awr-awpp)*(e+q)/awr
                     enddo
                     ! add in contribution to heating
                     nrr=1
                     npp=2
                     do ie=it,nes
                        e=xss(esz+ie-1)/emev
                        call terpa(h,e,en,idis,scr(llht),npp,nrr)
                        ss=0
                        if (ie.ge.iaa) ss=xss(2+k+ie-iaa)
                        xss(phn+2+ie-it)=xss(phn+2+ie-it)+h*ss
                        if (xss(tot+ie-1).ne.zero)&
                          xss(thn+ie-1)=xss(thn+ie-1)&
                          +h*ss/xss(tot+ie-1)
                     enddo
                  else
                     call skip6(nin,0,0,scr,law)
                  endif
               enddo
            endif
         enddo
      enddo

      !--assign the pointers for the distribution data
      ldlwp=nex
      xss(ixsa+neixs*(itype-1)+10)=ldlwp
      dlwp=ldlwp+ntrp
      xss(ixsa+neixs*(itype-1)+11)=dlwp
      nex=dlwp

      !--energy distributions from discrete levels in file 4.
      !--only neutron producing reactions are given here.
      do i=1,ntrp
         mt=nint(xss(mtrp+i-1))
         mtt=0
         ir=0
         do while (mtt.ne.mt)
            ir=ir+1
            mtt=nint(xss(mtr+ir-1))
            q=xss(lqr+ir-1)
            k=nint(xss(lsig+ir-1))+sig-1
            n=nint(xss(k+1))
            iaa=nint(xss(k))
         enddo
         jp=0
         do j=1,mtx
            if (mfm(j).eq.4.and.mtm(j).eq.mt&
              .and.mt.ge.50.and.mt.lt.91) jp=i
         enddo
         if (jp.gt.0) then
            ! distribution using law33
            xss(ldlwp+jp-1)=nex-dlwp+1
            last=nex
            xss(nex)=0
            xss(nex+1)=33
            nex=nex+3
            xss(nex)=0
            xss(nex+1)=2
            xss(nex+2)=sigfig(xss(esz+iaa-1)/emev,7,0)
            xss(nex+3)=sigfig(xss(esz+nes-1)/emev,7,0)
            xss(nex+4)=1
            xss(nex+5)=1
            nex=nex+2+2*2
            xss(last+2)=nex-dlwp+1
            ! amass=awr/awi
            ! aprime=awp/awi
            ! xss(nex)=sigfig((1+amass)*(-q)/amass,7,0)
            ! xss(nex+1)=&
            !  sigfig(amass*(amass+1-aprime)/(1+amass)**2,7,0)
            ! with awi = 0 for photons
            awp=1 ! MF4 is for outgoing neutrons
            xss(nex)=sigfig(-q,7,0)
            xss(nex+1)=sigfig((awr-awp)/awr,7,0)
            nex=nex+2
            ! neutron ebar for this reaction
            ! and local heating from recoil+photon
            ! neglecting photon momentum.
            do j=iaa,nes
               e=xss(esz+j-1)/emev
               ebar=awr*(e-abs(q))/(1+awr)
               hh=ebar*xss(2+k+j-iaa)
               xss(phn+2+j-it)=xss(phn+2+j-it)+hh
               hh=(e-abs(q))*xss(2+k+j-iaa)-hh
               if (xss(tot+j-1).ne.zero)&
                  xss(thn+j-1)=xss(thn+j-1)+hh/xss(tot+j-1)
            enddo
         endif
      enddo

      !--energy distributions from file 5.
      !--only neutron producing reactions are given here.
      do i=1,ntrp
         mt=nint(xss(mtrp+i-1))
         mtt=0
         ir=0
         do while (mtt.ne.mt)
            ir=ir+1
            mtt=nint(xss(mtr+ir-1))
            q=xss(lqr+ir-1)
            k=nint(xss(lsig+ir-1))+sig-1
            n=nint(xss(k+1))
            iaa=nint(xss(k))
         enddo
         jp=0
         do j=1,mtx
            if (mfm(j).eq.5.and.mtm(j).eq.mt) jp=i
         enddo
         if (jp.gt.0) then
            call findf(matd,5,mt,nin)
            call contio(nin,0,0,scr,nb,nw)
            nk=nint(scr(5))
            xss(ldlwp+jp-1)=nex-dlwp+1

            !--loop over subsections
            do ik=1,nk
               call tab1io(nin,0,0,scr,nb,nw)
               u=scr(1)
               lf=nint(scr(4))
               m=nint(scr(5))
               n=nint(scr(6))
               jnt=nint(scr(8))
               if (ik.gt.1) xss(last)=nex-dlwp+1
               last=nex
               xss(nex)=0
               xss(nex+1)=lf
               if (lf.eq.1) xss(nex+1)=4
               if (lf.eq.12) xss(nex+1)=4
               if (m.ne.1.or.jnt.ne.2) then
                  ja=nex+3
                  xss(ja)=m
                  jb=ja+m
                  do j=1,m
                     xss(ja+j)=scr(5+2*j)
                     xss(jb+j)=scr(6+2*j)
                  enddo
                  nex=jb+m+1
               else
                  xss(nex+3)=0
                  nex=nex+4
               endif
               xss(nex)=n
               l=5+2*m
               do j=1,n
                  xss(j+nex)=sigfig(scr(2*j+l)/emev,7,0)
                  xss(j+n+nex)=scr(2*j+1+l)
               enddo
               nex=nex+2*n+1
               xss(last+2)=nex-dlwp+1

               !--store data according to type
               if (lf.eq.1) then
                  call tab2io(nin,0,0,scr,nb,nw)
                  m=nint(scr(5))
                  n=nint(scr(6))
                  jnt=nint(scr(8))
                  jnt=mod(jnt,10)
                  if (jnt.gt.2) jnt=2
                  if (m.ne.1.or.jnt.ne.2) then
                     xss(nex)=m
                     do j=1,m
                        xss(j+nex)=scr(2*j+5)
                        jnt=nint(scr(2*j+6))
                        jnt=mod(jnt,10)
                        if (jnt.gt.2) jnt=2
                        xss(j+m+nex)=jnt
                     enddo
                     nex=nex+1+2*m
                  else
                     xss(nex)=0
                     nex=nex+1
                  endif
                  xss(nex)=n
                  nexn=nex+n
                  nexd=nexn+n+1
                  ne=n
                  do j=1,ne
                     call tab1io(nin,0,0,scr,nb,nw)
                     jscr=1
                     do while (nb.ne.0)
                        jscr=jscr+nw
                        call moreio(nin,0,0,scr(jscr),nb,nw)
                     enddo
                     e=c2h
                     xss(nex+j)=sigfig(e/emev,6,0)
                     xss(nexn+j)=nexd-dlwp+1
                     m=n1h
                     n=n2h
                     jnt=nint(scr(6+2*m))
                     xss(nexd)=jnt
                     xss(nexd+1)=n
                     nexd=nexd+1
                     xss(nexd+1+2*n)=0
                     do ki=1,n
                        ep=scr(5+2*m+2*ki)
                        ll=5+2*m+2*ki
                        xss(ki+nexd)=sigfig(scr(ll)/emev,7,0)
                        xss(ki+n+nexd)=sigfig(scr(ll+1)*emev,7,0)
                        if (xss(ki+n+nexd).lt.rmin) xss(ki+n+nexd)=0
                        if (ki.gt.1.and.jnt.eq.1) xss(ki+2*n+nexd)=&
                          xss(ki+2*n-1+nexd)+scr(ll-1)*(scr(ll)-scr(ll-2))
                        if (ki.gt.1.and.jnt.eq.2) xss(ki+2*n+nexd)=&
                          xss(ki+2*n-1+nexd)+((scr(ll-1)&
                          +scr(ll+1))/2)*(scr(ll)-scr(ll-2))
                     enddo
                     !renormalize
                     renorm=1
                     if (xss(3*n+nexd).ne.zero)&
                       renorm=1/xss(3*n+nexd)
                     do ki=1,n
                        xss(ki+n+nexd)=&
                          sigfig(xss(ki+n+nexd)*renorm,7,0)
                        xss(ki+2*n+nexd)=&
                          sigfig(xss(ki+2*n+nexd)*renorm,9,0)
                     enddo
                     nexd=nexd+3*n+1
                  enddo
                  nex=nexd
               else if (lf.eq.7.or.lf.eq.9) then
                  call tab1io(nin,0,0,scr,nb,nw)
                  m=nint(scr(5))
                  n=nint(scr(6))
                  jnt=nint(scr(8))
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
                  xss(nex)=n
                  do j=1,n
                     xss(j+nex)=sigfig(scr(5+2*j+2*m)/emev,7,0)
                     xss(j+n+nex)=sigfig(scr(6+2*j+2*m)/emev,7,0)
                  enddo
                  nex=nex+1+2*n
                  xss(nex)=u/emev
                  nex=nex+1
                  nrr=1
                  npp=2
                  do ie=iaa,nes
                     e=xss(esz+ie-1)
                     if (mth.eq.18) then
                        ipp=2
                        irr=1
                        call terpa(y,e,en,idis,fnubar,ipp,irr)
                     else
                        y=1
                        if (mth.eq.11.or.mth.eq.16.or.mth.eq.24.or.&
                            mth.eq.30.or.mth.eq.41.or.mth.eq.154.or.&
                            mth.eq.159.or.mth.eq.176.or.mth.eq.190.or.&
                            (mth.ge.875.and.mth.le.891)) then
                           y=2
                        elseif (mth.eq.17.or.mth.eq.25.or.mth.eq.42.or.&
                                mth.eq.157.or.mth.eq.172.or.mth.eq.177.or.&
                                mth.eq.179.or.(mth.ge.181.and.mth.le.199)) then
                           y=3
                        elseif (mth.eq.37.or.mth.eq.156.or.mth.eq.165.or.&
                                mth.eq.169.or.mth.eq.173.or.mth.eq.178.or.&
                                (mth.ge.194.and.mth.le.196)) then
                           y=4
                        elseif (mth.eq.152.or.mth.eq.162.or.mth.eq.166.or.&
                                mth.eq.170.or.mth.eq.174.or.mth.eq.200) then
                           y=5
                        elseif (mth.eq.153.or.mth.eq.163.or.mth.eq.167.or.&
                                mth.eq.171.or.mth.eq.175) then
                           y=6
                        elseif (mth.eq.160.or.mth.eq.164.or.mth.eq.168) then
                           y=7
                        elseif (mth.eq.161) then
                           y=8
                        endif
                     endif
                     call terpa(theta,e,en,idis,scr,npp,nrr)
                     x=0
                     if (theta.gt.zero) x=(e-u)/theta
                     ebar=0
                     if (lf.eq.7.and.x.gt.zero) then
                        anorm=sqrt(theta)**3*(sqrt(pi)*(1-erfc(x))/2&
                           -sqrt(x)*exp(-x))
                        ebar=3*theta/2-(sqrt(theta)**2&
                           *sqrt(x)**3*exp(-x))/anorm
                     else if (lf.eq.9.and.x.gt.zero) then
                        anorm=theta**2*(1-exp(-x)*(1+x))
                        ebar=2*theta-(theta**3*x**2*exp(-x))/anorm
                     endif
                     hh=y*(ebar/emev)*xss(2+k+ie-iaa)
                     xss(phn+2+ie-it)=xss(phn+2+ie-it)+hh
                     if (xss(tot+ie-1).ne.zero)&
                       xss(thn+ie-1)=xss(thn+ie-1)-hh/xss(tot+ie-1)
                  enddo
               else
                  call error('acephn','file 5 law not ready',' ')
               endif
            enddo
         endif
      enddo

      !--now check for file-6 distribution data
      !--find each subsection that contributes to this production
      if (mf6.ne.0) then
         jp=0
         do i=1,ntrp
            mt=nint(xss(mtrp+i-1))
            do j=1,mtx
               if (mfm(j).eq.6.and.mtm(j).eq.mt) jp=i
            enddo
            if (jp.ne.0) then
               call findf(matd,6,mt,nin)
               call contio(nin,0,0,scr,nb,nw)
               lct=nint(scr(4))
               nk=nint(scr(5))
               mtt=0
               ir=0
               do while (mtt.ne.mth)
                  ir=ir+1
                  mtt=nint(xss(mtr+ir-1))
                  q=xss(lqr+ir-1)
                  k=nint(xss(lsig+ir-1))+sig-1
                  n=nint(xss(k+1))
                  iaa=nint(xss(k))
               enddo
               ik=0
               do while (ik.lt.nk)
                  ik=ik+1
                  ! read the multiplicity
                  call tab1io(nin,0,0,scr,nb,nw)
                  jscr=1+nw
                  do while (nb.ne.0)
                     call moreio(nin,0,0,scr(jscr),nb,nw)
                     jscr=jscr+nw
                  enddo
                  ! retrieve izap, awp and the law
                  izap=nint(scr(1))
                  awp=scr(2)
                  law=nint(scr(4))

                  if (izap.ne.ip) then
                     call skip6(nin,0,0,scr,law)
                  else

                     ! if this is fission, replace the yield with the nubar - as before
                     if (mth.eq.18.and.izap.eq.1.and.mtxnu.gt.0) then
                        call copynubar(scr,fnubar,jscr)
                     endif

                     xss(ldlwp+jp-1)=nex-dlwp+1  ! locator, points to LNW
                     last=nex
                     xss(last)=0                 ! LNW
                     xss(last+1)=0               ! LAW set to 0
                     xss(last+2)=0               ! IDAT set to 0
                     nex=nex+3                   ! nex points to NR

                     !--we can only process law=1, 2, and 4 currently
                     if (law.ne.1.and.law.ne.2.and.law.ne.4) then
                        call mess('acephn','can only do law=1,',&
                          ' law=2, or law=4 currently')
                        call skip6(nin,0,0,scr,law)
                     else if (law.eq.1) then
                        ll=jscr
                        call tab2io(nin,0,0,scr(ll),nb,nw)
                        lang=nint(scr(ll+2))
                        lep=nint(scr(ll+3))
                        ne=nint(scr(ll+5))       ! number of incident energies
                        if (lang.eq.1) then      ! legendre polynomials to law=61
                           xss(last+1)=61        ! LAW
                        else if (lang.eq.2) then ! Kalbach-Mann to law=44
                           xss(last+1)=44        ! LAW
                        else
                           write(text,'(''lang='',i3,'' not supported for law='',i2)')lang,law
                           call error('acephn',text,'')
                        endif
                        xss(landp+jp-1)=-1     ! angular included in energy distribution
                        nr=0
                        xss(nex)=nr            ! NR set to 0
                        lee=nex                ! lee points to NR
                        nex=nex+2*nr+1
                        nle=2
                        xss(nex)=nle           ! number of energies, NE, default to 2
                        nex=nex+1+2*nle        ! leaving room for E(1:2), P(1:2).  nex points to LDAT(1)
                        xss(last+2)=nex-dlwp+1 ! IDAT
                        nr=0
                        xss(nex)=nr            ! LDAT(1) = NR set to 0
                        nex=nex+1+2*nr
                        xss(nex)=ne            ! LDAT(2) = NE set to number of incident energies
                        nex=nex+1
                        lle=nex                ! lle points to LDAT(3) = E(1)
                        nex=lle+2*ne           ! nex points to start of first distribution

                        ! scr(llh) up to scr(lld-1) is set up for heating
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

                        ! go over each incident energy value
                        do ie=1,ne

                           ! read distribution, store starting at scr(lld)
                           ll=lld
                           call listio(nin,0,0,scr(ll),nb,nw)
                           ll=ll+nw
                           if (ll.gt.nwscr) call error('acephn',&
                               'scr array overflow in file 6 list',' ')
                           do while (nb.ne.0)
                              call moreio(nin,0,0,scr(ll),nb,nw)
                              ll=ll+nw
                              if (ll.gt.nwscr) call error('acephn',&
                                  'scr array overflow in file 6 list',' ')
                           enddo
                           nd=nint(scr(lld+2))
                           na=nint(scr(lld+3))
                           ng=nint(scr(lld+5))
                           ncyc=na+2

                           ! set energy range and probability for this law
                           ! only first and last incident energy needed
                           ! probability set to 1 all energies
                           if (ie.eq.1) then
                              xss(lee+2)=sigfig(scr(lld+1)/emev,7,0)
                              xss(lee+4)=1
                           else if (ie.eq.ne) then
                              xss(lee+3)=sigfig(scr(lld+1)/emev,7,0)
                              xss(lee+5)=1
                           endif

                           ! set incident energy and locator for the current distribution
                           xss(lle+ie-1)=sigfig(scr(lld+1)/emev,7,0) ! Ein(ie)
                           ee=xss(lle+ie-1)
                           xss(lle+ne+ie-1)=nex-dlwp+1  ! locator for distribution
                           xss(nex)=lep+10*nd           ! INTT for this secondary energy distribution
                           xss(nex+1)=ng                 ! NP
                           if (lang.eq.1) nexc=nex+2+4*ng ! only needed if law=61

                           amass=awp*emc2
                           avadd=ee/(awr*emc2)
                           avlab=0
                           avll=0

                           ! go over the outgoing energies
                           do ig=1,ng
                              ! outgoing energy
                              xss(nex+1+ig)=&
                                sigfig(scr(lld+6+ncyc*(ig-1))/emev,7,0)
                              ! pdf
                              if (ig.le.nd) then
                                xss(nex+1+ig+ng)=&
                                  sigfig(scr(lld+7+ncyc*(ig-1)),7,0)
                              else
                                xss(nex+1+ig+ng)=&
                                  sigfig(scr(lld+7+ncyc*(ig-1))*emev,7,0)
                              endif
                              test=xss(nex+1+ig+ng)
                              if (test.gt.zero.and.test.lt.small)&
                                xss(nex+1+ig+ng)=small
                              ! cdf
                              if (ig.eq.1) then
                                 ! initial value
                                 if (nd.eq.0) then
                                    xss(nex+1+ig+2*ng)=0
                                 else
                                    xss(nex+1+ig+2*ng)=&
                                                  scr(lld+7+ncyc*(ig-1))
                                 endif
                              elseif (ig.le.nd) then
                                 ! discrete photon
                                 xss(nex+1+ig+2*ng)=xss(nex+ig+2*ng)+&
                                                  scr(lld+7+ncyc*(ig-1))
                              elseif (ig.eq.nd+1) then
                                 ! start of continuum
                                 xss(nex+1+ig+2*ng)=xss(nex+ig+2*ng)
                              endif
                              if (ig.gt.nd+1) then
                                 ! continuum
                                 if (lep.eq.1) then
                                    ! histogram
                                    xss(nex+1+ig+2*ng)=xss(nex+ig+2*ng)&
                                       +scr(lld+7+ncyc*(ig-2))&
                                       *(scr(lld+6+ncyc*(ig-1))&
                                       -scr(lld+6+ncyc*(ig-2)))
                                 elseif (lep.eq.2) then
                                    ! lin-lin
                                    xss(nex+1+ig+2*ng)=xss(nex+ig+2*ng)&
                                       +((scr(lld+7+ncyc*(ig-2))&
                                       +scr(lld+7+ncyc*(ig-1)))/2)&
                                       *(scr(lld+6+ncyc*(ig-1))&
                                       -scr(lld+6+ncyc*(ig-2)))
                                 endif
                              endif
                              if (lang.eq.2) then
                                 ! kalbach-mann
                                 rkal=scr(lld+8+ncyc*(ig-1))
                                 ep=xss(nex+1+ig)
                                 akal=bachaa(izai,izap,za,ee,ep)
                                 xss(nex+1+ig+3*ng)=sigfig(rkal,7,0) ! r
                                 xss(nex+1+ig+4*ng)=sigfig(akal,7,0) ! a
                              else

                                 if (lang.eq.1) then
                                    xss(nex+1+ig+3*ng)=nexc-dlwp+1  !pointer to angdist table
                                    ! convert lang=1 list in scr to a normalized P(1) to P(NA) list for ptleg2
                                    scr(ll)=0
                                    scr(ll+1)=scr(lld+6+ncyc*(ig-1))        !EOUT(ig)
                                    scr(ll+2)=0
                                    scr(ll+3)=0
                                    scr(ll+4)=na                            !P(l) order (zero is allowed)
                                    scr(ll+5)=0
                                    do ia=1,na
                                       lll=lld+7+ncyc*(ig-1)
                                       scr(ll+5+ia)=0
                                       if (scr(lll).ne.zero) then
                                          scr(ll+5+ia)=scr(lll+ia)/scr(lll) !P(n)/P(0)
                                       endif
                                    enddo

                                    call ptleg2(scr(ll))  !P(l) list in, tab1 (mu,pdf) out

                                    intmu=2
                                    xss(nexc)=intmu
                                    nmu=nint(scr(ll+5))
                                    xss(nexc+1)=nmu
                                    do imu=1,nmu
                                       xss(nexc+1+imu)=sigfig(scr(ll+6+2*imu),7,0)
                                       xss(nexc+1+nmu+imu)=sigfig(scr(ll+7+2*imu),7,0)
                                       if (imu.eq.1) then
                                           xss(nexc+1+2*nmu+imu)=0
                                       else if (imu.eq.nmu) then
                                           xss(nexc+1+2*nmu+imu)=1
                                       else
                                          del=scr(ll+6+2*imu)-scr(ll+4+2*imu)
                                          av=(scr(ll+7+2*imu)+scr(ll+5+2*imu))/2
                                          xss(nexc+1+2*nmu+imu)=&
                                                               xss(nexc+1+2*nmu+imu-1)+del*av
                                          xss(nexc+1+2*nmu+imu)=&
                                                        sigfig(xss(nexc+1+2*nmu+imu),7,0)
                                       endif
                                    enddo
                                    nexc=nexc+2+3*nmu
                                 endif

                              endif
                              ! average lab energy
                              if (ig.ne.1) then
                                 eavi=xss(nex+1+ig)
                                 if (lang.ne.2.or.na.eq.0) then
                                    avl=eavi
                                 else
                                    avcm=sqrt(2*eavi/amass)
                                    sign=1
                                    avl=eavl(akal,amass,avcm,avadd,&
                                      rkal,sign)
                                 endif
                                 dele=xss(nex+1+ig)-xss(nex+ig)
                                 if (lep.eq.1) then
                                    avav=xss(nex+ig+ng)*(avll+avl)/2
                                 else
                                    avav=(xss(nex+ig+ng)*avll&
                                      +xss(nex+ig+1+ng)*avl)/2
                                 endif
                                 avlab=avlab+avav*dele
                                 avll=avl
                              endif
                           enddo  !end of loop over secondary energies
                           ! renormalize cummulative probabilities
                           renorm=one/xss(nex+1+3*ng)
                           do ig=1,ng
                              xss(nex+1+ng+ig)=&
                                sigfig(renorm*xss(nex+1+ng+ig),7,0)
                              xss(nex+1+2*ng+ig)=&
                                sigfig(renorm*xss(nex+1+2*ng+ig),9,0)
                           enddo
                           scr(llh+6+2*ie)=ee
                           scr(llh+7+2*ie)=avlab
                           if(lang.eq.1)then
                              nex=nexc
                           else
                              nex=nex+2+(2*na+3)*ng
                           endif
                        enddo  !end of loop over incident energies

                        !add in contribution to heating
                        !for this subsection
                        nrr=1
                        npp=2
                        do ie=iaa,nes
                           e=xss(esz+ie-1)
                           call terpa(y,e,en,idis,scr(lly),npp,nrr)
                           scr(lld+ie-iaa)=y*xss(2+k+ie-iaa)
                        enddo
                        it=nint(xss(pxs))
                        nrr=1
                        npp=2
                        do ie=iaa,nes
                           e=xss(esz+ie-1)/emev
                           call terpa(h,e,en,idis,scr(llh),npp,nrr)
                           hh=h*scr(lld+ie-iaa)
                           xss(phn+2+ie-it)=xss(phn+2+ie-it)+hh
                           do ii=1,mtx
                              if (mfm(ii).eq.6.and.mtm(ii).eq.mth&
                                .and.nr6(ii).eq.0.and.&
                                xss(tot+ie-1).ne.0)&
                                xss(thn+ie-1)=xss(thn+ie-1)&
                                -hh/xss(tot+ie-1)
                           enddo
                        enddo
                     else if (law.eq.2) then
                        ll=jscr
                        xss(last+1)=33
                        xss(nex)=0
                        xss(nex+1)=2
                        nnr=nint(scr(5))
                        nnp=nint(scr(6))
                        xss(nex+2)=sigfig(scr(7+2*nnr)/emev,7,0)
                        xss(nex+3)=&
                           sigfig(scr(5+2*nnr+2*nnp)/emev,7,0)
                        xss(nex+4)=1
                        xss(nex+5)=1
                        nex=nex+2+2*2
                        xss(last+2)=nex-dlwp+1
                        ! amass=awr/awi
                        ! aprime=awp/awi
                        ! xss(nex)=sigfig((1+amass)*(-q)/amass,7,0)
                        ! xss(nex+1)=&
                        !  sigfig(amass*(amass+1-aprime)/(1+amass)**2,7,0)
                        ! with awi = 0 for photons
                        xss(nex)=sigfig(-q,7,0)
                        xss(nex+1)=sigfig((awr-awp)/awr,7,0)
                        nex=nex+2
                        call tab2io(nin,0,0,scr(ll),nb,nw)
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
                        enddo
                     else if (law.eq.4) then
                        xss(last+1)=33
                        xss(nex)=0
                        xss(nex+1)=2
                        nnr=nint(scr(5))
                        nnp=nint(scr(6))
                        xss(nex+2)=sigfig(scr(7+2*nnr)/emev,7,0)
                        xss(nex+3)=&
                           sigfig(scr(5+2*nnr+2*nnp)/emev,7,0)
                        xss(nex+4)=1
                        xss(nex+5)=1
                        nex=nex+2+2*2
                        xss(last+2)=nex-dlwp+1
                        ! amass=awr/awi
                        ! aprime=awp/awi
                        ! xss(nex)=sigfig((1+amass)*(-q)/amass,7,0)
                        ! xss(nex+1)=&
                        !  sigfig(amass*(amass+1-aprime)/(1+amass)**2,7,0)
                        ! with awi = 0 for photons
                        xss(nex)=sigfig(-q,7,0)
                        xss(nex+1)=sigfig((awr-awp)/awr,7,0)
                        nex=nex+2
                     endif
                  endif
               enddo
            endif
         enddo

         !--divide the heating contribution by the total xsec
         !--and add it into the total heating, if appropriate
         do ie=it,nes
            if (xss(tot+ie-1).ne.zero) then
               xss(phn+2+ie-it)=xss(phn+2+ie-it)/xss(tot+ie-1)
            endif
            if (ip.gt.1) then
               xss(thn+ie-1)=xss(thn+ie-1)+xss(phn+2+ie-it)
            endif
            xss(phn+2+ie-it)=sigfig(xss(phn+2+ie-it),7,0)
         enddo
      endif

   !--continue the loop over types
   enddo
   lxs=nex-1

   !--convert main energy grid to mev
   do i=1,nes
      xss(esz+i-1)=sigfig(xss(esz+i-1)/emev,7,0)
      xss(thn+i-1)=sigfig(xss(thn+i-1),7,0)
   enddo
   call closz(nin)

   !--print and write the photonuclear file
   zaid=za+suff
   if (mcnpx.eq.0) then
      write(hz,'(f9.2,''u'')') zaid
   else
      write(hz,'(f10.3,''pu '')') zaid
   endif
   aw0=awr
   tz=0
   call dater(hdt)
   hd='  '//hdt
   write(hm,'(''   mat'',i4)') matd
   if (iprint.gt.0) call phnprt(hk)
   call phnout(ityp,nace,ndir,mcnpx,hk,izn,awn)

   !--finished
   return
   end subroutine acephn

   subroutine phnfix(itype,nin,nout,ndir,iprint,nplot,mcnpx,suff,&
     nxtra,hk,izn,awn)
   !-------------------------------------------------------------------
   ! Print and or edit an ACE photonuclear file.
   !-------------------------------------------------------------------
   use util ! provides closz
   ! externals
   integer::itype,nin,nout,ndir,iprint,nplot,mcnpx,nxtra
   real(kr)::suff
   character(70)::hk
   integer::izn(16)
   real(kr)::awn(16)
   ! internals
   integer::l,j,n,max,iza,i
   integer::izo(16)
   real(kr)::awo(16)
   character(70)::hko
   real(kr)::zaid
   character(3)::ht
   character(9)::str
   real(kr),parameter::zero=0

   integer,parameter::ner=1
   integer,parameter::nbw=1

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
      read (nin,'(8i9)') lxs,za,nes,ntr,ntype,npixs,neixs,nxsd(1:8),tvn,&
        esz,tot,non,els,thn,mtr,lqr,lsig,sig,ixsa,ixs,jxsd(1:21)
      n=(lxs+3)/4
      l=0
      do i=1,n
         read (nin,'(4e20.0)') (xss(l+j),j=1,4)
         l=l+4
      enddo

   !--read type 2 ace format file
   else if (itype.eq.2) then
      if (mcnpx.eq.0) then
        read(nin) hz(1:10),aw0,tz,hd,hko,hm,(izo(i),awo(i),i=1,16),&
          lxs,za,nes,ntr,ntype,npixs,neixs,nxsd(1:8),tvn,&
          esz,tot,non,els,thn,mtr,lqr,lsig,sig,ixsa,ixs,jxsd(1:21)
      else
        read(nin) hz(1:13),aw0,tz,hd,hko,hm,(izo(i),awo(i),i=1,16),&
          lxs,za,nes,ntr,ntype,npixs,neixs,nxsd(1:8),tvn,&
          esz,tot,non,els,thn,mtr,lqr,lsig,sig,ixsa,ixs,jxsd(1:21)
      endif
      n=(lxs+ner-1)/ner
      l=0
      do i=1,n
         max=lxs-l
         if (max.gt.ner) max=ner
         read (nin) (xss(l+j),j=1,max)
         l=l+max
      enddo
      call closz(-nin)
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
      iza=nint(zaid)
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
   if (iprint.gt.0) call phnprt(hk)
   if (nout.gt.0) call phnout(itype,nout,ndir,mcnpx,hk,izn,awn)
   if (nplot.ne.0) call phnplo(nplot,hk)

   return
   end subroutine phnfix

   subroutine phnprt(hk)
   !-------------------------------------------------------------------
   ! Print ACE photo-nuclear data from memory.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util ! provides error
   use acecm ! provides mtname
   ! externals
   character(70)::hk
   ! internals
   integer::i,k,iaa,ib,ic,i1,i3,list,l,na,nb,nc,j,n,nc2,nkk
   integer::m,ii,imt,naa,ie,mt,l1,l2,ne,ll,nbina,nbin1
   integer::intt,np,ln,law,l3,loci,nn,nd
   integer::ipt,ntrp,pxs,phn,mtrp,tyrp,lsigp,sigp,landp,andp,ldlwp,dlwp
   integer::ip,locj,intmu,nmu,imu
   real(kr)::e,xs,heat,e2
   integer::imn(8),imx(8),loc(8)
   character(10)::name,title(16)
   character(15)::kk(40)
   character(6)::blank='      '
   character(6)::ek='energy'
   character(12)::dashes='------------'

   !--print information block
   write(nsyso,'(''1''////////&
     &38x,''zaid'',1x,a13/39x,''awr'',f10.3/&
     &6x,''***********************'',9x,''temp'',1p,e10.2/&
     &6x,''*                     *'',9x,''date'',a10/&
     &6x,''*    photo-nuclear    *'',10x,''mat'',a10/&
     &6x,''*                     *''/&
     &6x,''*   ace format file   *'',10x,''lxs'',i10/&
     &6x,''*                     *'',11x,''za'',i10/&
     &6x,''*     processed by    *'',10x,''nes'',i10/&
     &6x,''*                     *'',10x,''ntr'',i10/&
     &6x,''*        njoy         *'',8x,''ntype'',i10/&
     &6x,''*                     *'',8x,''npixs'',i10/&
     &6x,''***********************'',8x,''neixs'',i10/&
     &39x,''tvn'',i10/39x,''esz'',i10/39x,''mtr'',i10/&
     &39x,''lqr'',i10/38x,''lsig'',i10/39x,''sig'',i10/&
     &//6x,''hk--- '',a70///)')&
     hz,aw0,tz,hd,hm,lxs,za,nes,ntr,ntype,npixs,neixs,tvn,&
     esz,mtr,lqr,lsig,sig,hk

   !--photonuclear reaction descriptions
   write(nsyso,'(//&
     &7x,''photonuclear reactions''/7x,''----------------------'')')
   write(nsyso,'(/&
     &7x,''reaction'',8x,''mt'',6x,''lsig'',6x,&
     &''           emin           emax              q''/&
     &7x,''--------'',8x,''--'',6x,''----'',6x,&
     &''   ------------   ------------   ------------'')')
   do i=1,ntr
      k=nint(xss(lsig+i-1)+sig-1)
      iaa=nint(xss(k))
      ib=iaa+nint(xss(k+1))-1
      ic=nint(xss(mtr+i-1))
      call mtname(ic,name,0)
      if (name(1:1).eq.'(') then
         name(2:2)='g'
      endif
      i1=nint(xss(mtr+i-1))
      i3=nint(xss(lsig+i-1))
      write(nsyso,'(7x,a10,i8,i10,6x,1p,3e15.6)')&
        name,i1,i3,xss(esz+iaa-1),xss(esz+ib-1),xss(lqr+i-1)
   enddo

   !--print principal cross sections
   do i=1,nes
      if (mod(i,57).eq.1) then
         write(nsyso,'(''1''/&
           &5x,''i'',5x,''energy'',11x,'' total  '',7x,&
           &''nonelastic'',5x,''elastic'',8x,''heating''/&
           &3x,''---'',3x,''--------------'',&
           &4(3x,''------------''))')
      endif
      if (els.eq.0.and.thn.eq.0) then
         write(nsyso,'(i6,1p,e17.8,2e15.6)') i,&
           xss(esz+i-1),xss(tot+i-1),xss(non+i-1)
      else if (els.eq.0) then
         write(nsyso,'(i6,1p,e17.8,2e15.6,15x,e15.6)') i,&
           xss(esz+i-1),xss(tot+i-1),xss(non+i-1),xss(thn+i-1)
      else if (thn.eq.0) then
         write(nsyso,'(i6,1p,e17.8,3e15.6)') i,&
           xss(esz+i-1),xss(tot+i-1),xss(non+i-1),xss(els+i-1)
      else
         write(nsyso,'(i6,1p,e17.8,4e15.6)') i,xss(esz+i-1),&
           xss(tot+i-1),xss(non+i-1),xss(els+i-1),xss(thn+i-1)
      endif
   enddo

   !--print reaction cross sections
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
            k=nint(xss(lsig+n-1)+sig-1)
            imn(j)=nint(xss(k))
            iaa=min0(iaa,imn(j))
            k=k+1
            loc(j)=k
            imx(j)=imn(j)+nint(xss(k))-1
            ib=max0(ib,imx(j))
            k=iabs(nint(xss(mtr+n-1)))
            call mtname(k,title(j),0)
            if (title(j)(1:1).eq.'(') then
               title(j)(2:2)='g'
            endif
            j=j+1
         enddo
         nc2=nc
         nkk=nc
         do m=iaa,ib
            if (mod(m+1-iaa,57).eq.1) then
               write(nsyso,'(''1''/4x,''i'',5x,''energy'',11x,a10,&
                 &6(5x,a10))') (title(ii),ii=1,nc2)
               write(nsyso,'(3x,''---'',3x,''--------------'',&
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
            write(nsyso,'(1x,i5,1p,e17.8,7a15)')&
              m,xss(esz+m-1),(kk(i),i=1,nkk)
         enddo
      enddo
   endif

   !--print particle production data
   do i=1,ntype
      ipt=nint(xss(ixsa+neixs*(i-1)))
      ntrp=nint(xss(ixsa+neixs*(i-1)+1))
      if (ipt.eq.1) write(nsyso,'(''1''//&
        &'' ********** neutron production **********'')')
      if (ipt.eq.2) write(nsyso,'(''1''//&
        &'' ********** photon production **********'')')
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
      pxs=nint(xss(ixsa+neixs*(i-1)+2))
      phn=nint(xss(ixsa+neixs*(i-1)+3))
      mtrp=nint(xss(ixsa+neixs*(i-1)+4))
      tyrp=nint(xss(ixsa+neixs*(i-1)+5))
      lsigp=nint(xss(ixsa+neixs*(i-1)+6))
      sigp=nint(xss(ixsa+neixs*(i-1)+7))
      landp=nint(xss(ixsa+neixs*(i-1)+8))
      andp=nint(xss(ixsa+neixs*(i-1)+9))
      ldlwp=nint(xss(ixsa+neixs*(i-1)+10))
      dlwp=nint(xss(ixsa+neixs*(i-1)+11))
      write(nsyso,'(/&
        &''   mtrp   tyrp   lsigp   landp   ldlwp''/&
        &''   ----   ----   -----   -----   -----'')')
      do imt=1,ntrp
         write(nsyso,'(i7,i7,i8,i8,i8,i9)') nint(xss(mtrp+imt-1)),&
           nint(xss(tyrp+imt-1)),nint(xss(lsigp+imt-1)),&
           nint(xss(landp+imt-1)),nint(xss(ldlwp+imt-1))
      enddo

      !--print production cross section and heating
      iaa=nint(xss(pxs))
      naa=nint(xss(pxs+1))
      write(nsyso,'(/'' production cross section''/&
        &4x,''ie'',9x,''energy'',13x,''xs'',8x,''heating''/&
        &2x,''----'',3(2x,''-------------''))')
      do ie=1,naa
         e=xss(esz+iaa+ie-2)
         xs=xss(pxs+1+ie)
         heat=xss(phn+1+ie)
         write(nsyso,'(i6,1p,3e15.6)') iaa+ie-1,e,xs,heat
      enddo

      !--print production distributions
      if (ipt.eq.1)&
        write(nsyso,'(''1''/'' neutron production distributions'')')
      if (ipt.eq.2)&
        write(nsyso,'(''1''/'' photon production distributions'')')
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
      do imt=1,ntrp
         mt=nint(xss(mtrp+imt-1))
         write(nsyso,'(/'' contributions from mt='',i3,&
           &'' ------->'')') mt
         l1=nint(xss(lsigp+imt-1))
         l2=sigp+l1-1
         ne=nint(xss(l2+3))
         write(nsyso,'(/''   yields:''//&
           &''              energy         yield''/&
           &6x,2(2x,''------------''))')
         do ii=1,ne
            write(nsyso,'(6x,1p,2e14.6)')&
              xss(l2+4+ii-1),xss(l2+4+ne+ii-1)
         enddo
         ll=nint(xss(landp+imt-1))

         !--angular distributions
         if (ll.gt.0) then
            write(nsyso,'(/''   angular distribution:'')')
            ll=ll+andp-1
            ne=nint(xss(ll))
            na=ll
            nb=na+ne

            !--check on elastic format
            k=nint(xss(nb+1))
            if (k.ge.0) then

               !--equally-probable bins format
               list=(ne+7)/8
               nbina=32
               nbin1=nbina+1
               do l=1,list
                  iaa=(l-1)*8+1
                  ib=min0(ne,iaa+7)
                  ic=ib-iaa+1
                  j=1
                  do m=iaa,ib
                     k=nint(xss(m+nb))
                     if (k.gt.0) k=k+andp
                     loc(j)=k
                     j=j+1
                  enddo
                  if (l.eq.1) then
                     write(nsyso,'(6x,''ne ='',i4)') ne
                     write(nsyso,'(1x)')
                  endif
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
            else

               !--cummulative-distribution format
               write(nsyso,'(6x,''ne ='',i4)') ne
               do ii=1,ne
                  e=xss(na+ii)
                  k=nint(abs(xss(nb+ii)))+andp-1
                  intt=nint(xss(k))
                  np=nint(xss(k+1))
                  k=k+1
                  write(nsyso,'(/5x,'' incident particle energy ='',&
                    &1p,e14.6,''    int ='',i2,''    np ='',i4)')&
                    e,intt,np
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

         !--energy and energy-angle distributions
         l1=nint(xss(ldlwp+imt-1))
         l2=dlwp+l1-1
         ln=1
         do while (ln.ne.0)
            ln=nint(xss(l2))
            law=nint(xss(l2+1))
            write(nsyso,&
              '(/''   distributions (law'',i2,''):'')') law
            l3=dlwp+nint(xss(l2+2))-1

            !--fractional probability for this law
            if (ln.ne.0.or.l2.ne.dlwp+l1-1) then
               l=l2+3
               write(nsyso,'(8x,''probability of law'')')
               j=nint(xss(l))
               write(nsyso,'(12x,''nr ='',i4)') j
               if (j.gt.0) then
                  write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
                    (nint(xss(l+ii)),ii=1,j)
                  l=l+j
                  write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
                    (nint(xss(l+ii)),ii=1,j)
                  l=l+j
               endif
               l=l+1
               j=nint(xss(l))
               write(nsyso,'(12x,''ne ='',i4)') j
               write(nsyso,&
                 &'(12x,''e(i=1,ne) =   '',1p,6e14.6/(12x,7e14.6))')&
                 (xss(l+ii),ii=1,j)
               l=l+j
               write(nsyso,&
                 '(12x,''p(i=1,ne) =   '',1p,6e14.6/(12x,7e14.6))')&
                 (xss(l+ii),ii=1,j)
               l=l+j+1
               write(nsyso,'(8x,''data for law'')')
            endif

            !--law=4
            if (law.eq.4) then
               ne=nint(xss(l3+1))
               l=l3+2+2*ne
               do ie=1,ne
                  e2=xss(l3+2+ie-1)
                  loci=nint(xss(l3+2+ne+ie-1))+dlwp-1
                  intt=mod(nint(xss(loci)),10)
                  nd=nint(xss(loci)/10)
                  nn=nint(xss(loci+1))
                  loci=loci+1
                  write(nsyso,'(/6x,&
                    &'' incident energy = '',1p,e14.6,&
                    &''   intt ='',i2,''    nd = '',i4,&
                    &''    np = '',i4//&
                    &1x,8x,''energy'',11x,''pdf'',11x,''cdf'',&
                    &8x,''energy'',11x,''pdf'',11x,''cdf''/&
                    &1x,6(2x,''------------''))') e2,intt,nd,nn
                  do j=1,nn,2
                     if (j.lt.nn) then
                        write(nsyso,'(1x,1p,6e14.6)')&
                          xss(j+loci),xss(j+nn+loci),&
                          xss(j+2*nn+loci),xss(j+1+loci),&
                          xss(j+1+nn+loci),xss(j+1+2*nn+loci)
                     else
                        write(nsyso,'(1x,1p,3e14.6)')&
                          xss(j+loci),xss(j+nn+loci),&
                          xss(j+2*nn+loci)
                     endif
                  enddo
                  l=l+3*nn
               enddo
               l2=l

            !--law=44
            else if (law.eq.44) then
               ne=nint(xss(l3+1))
               l=l3+2+2*ne
               do ie=1,ne
                  e2=xss(l3+2+ie-1)
                  loci=nint(xss(l3+2+ne+ie-1))+dlwp-1
                  intt=mod(nint(xss(loci)),10)
                  nd=nint(xss(loci)/10)
                  nd=0
                  nn=nint(xss(loci+1))
                  loci=loci+1
                  write(nsyso,'(/6x,&
                    &'' incident energy = '',1p,e14.6,&
                    &''   intt ='',i2,''    nd = '',i4,''    np = '',&
                    &i4,//,1x,8x,''energy'',11x,''pdf'',11x,''cdf'',&
                    &13x,''r'',13x,''a''/&
                    &1x,5(2x,''------------'')/(1x,1p,5e14.6))')&
                    e2,intt,nd,nn,(xss(j+loci),xss(j+nn+loci),&
                    xss(j+2*nn+loci),xss(j+3*nn+loci),&
                    xss(j+4*nn+loci),j=1,nn)
                  l=l+4*nn
               enddo
               l2=l

            !--law 61
            else if (law.eq.61) then
               ne=nint(xss(l3+1))
               l=l3+2+2*ne
               do ie=1,ne
                  e2=xss(l3+2+ie-1)
                  loci=nint(xss(l3+2+ne+ie-1))+dlwp-1
                  intt=mod(nint(xss(loci)),10)
                  nd=nint(xss(loci)/10)
                  nn=nint(xss(loci+1))
                  loci=loci+1
                  write(nsyso,'(/6x,'' incident energy = '',1p,e14.6,&
                    &''   intt ='',i2,''    nd = '',i4,''    np = '',&
                    &i4)') e2,intt,nd,nn
                  do ip=1,nn
                     write(nsyso,'(/&
                       &6x,'' secondary energy = '',1p,e14.6/&
                       &6x,''              pdf = '',e14.6/&
                       &6x,''              cdf = '',e14.6)')&
                       xss(ip+loci),xss(ip+nn+loci),xss(ip+2*nn+loci)
                     locj=nint(xss(ip+3*nn+loci)+dlwp-1)
                     if (locj.ne.0) then
                        intmu=nint(xss(locj))
                        nmu=nint(xss(locj+1))
                        write(nsyso,'(&
                          &6x,''            intmu = '',i8/&
                          &6x,''              nmu = '',i8/&
                          &''         cosine           pdf           cdf'',&
                          &''        cosine           pdf           cdf''/&
                          &''   ------------  ------------  ------------'',&
                          &''  ------------  ------------  ------------'')')&
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
                      else
                        write(nsyso,'('' angular distribution is isotropic'')')
                      endif
                  enddo
                  l=l+4*nn
               enddo
               l2=l

            !--law=7 or 9
            else if (law.eq.7.or.law.eq.9) then
               l=l3
               j=nint(xss(l))
               write(nsyso,'(12x,''nr ='',i4)') j
               if (j.ne.0) then
                  write(nsyso,'(12x,''nbt(i=1,nr) = '',20i5)')&
                    (nint(xss(l+ii)),ii=1,j)
                  l=l+j
                  write(nsyso,'(12x,''int(i=1,nr) = '',20i5)')&
                    (nint(xss(l+ii)),ii=1,j)
                  l=l+j
               endif
               l=l+1
               j=nint(xss(l))
               write(nsyso,'(12x,''ne ='',i4)') j
               write(nsyso,'(12x,&
                 &''e(i=1,ne) = '',16x,1p,5e14.6/(12x,7e14.6))')&
                 (xss(l+ii),ii=1,j)
               l=l+j
               write(nsyso,'(12x,''theta(i=1,ne) = '',12x,1p,5e14.6/&
                 &(12x,7e14.6))') (xss(l+ii),ii=1,j)
               l=l+j+1
               write(nsyso,'(12x,''u =  '',9x,1p,e14.6)') xss(l)
               l2=l+1

            !--law=33
            else if (law.eq.33) then
               write(nsyso,'(12x,''eout = c*(e-ec)   ec ='',1p,&
                 &e14.6,5x,''c ='',e14.6)') xss(l3),xss(l3+1)

            !--law not installed
            else
               call error('phnprt','law not installed',' ')
            endif

         !--continue loop over partial distributions, if any
         enddo

      !--continue loop over mts for this particle type
      enddo

   !--continue loop over productions
   enddo
   return
   end subroutine phnprt

   subroutine phnout(itype,nout,ndir,mcnpx,hk,izn,awn)
   !-------------------------------------------------------------------
   ! Write photo-nuclear ACE data to output and directory files
   !-------------------------------------------------------------------
   use util  ! provides openz,closz,error
   use acecm ! provides write routines
   ! externals
   integer::itype,nout,ndir,mcnpx
   integer::izn(16)
   real(kr)::awn(16)
   character(70)::hk
   ! internals
   integer::l,n,ne,ip,mftype,nr,li,ir,nn,ll,k,np,nw,nmu,nrr
   integer::ii,lnw,law,kk,nern,lrec,j,i
   integer::ipt, ntrp, pxs, phn, mtrp, tyrp, lsigp, sigp, landp, andp, ldlwp, dlwp ! IXS
   integer::rlocator  ! locator index for reaction data
   integer::plocator  ! locator index for the particle IXS array
   integer::ielocator ! locator index for incident energy data
   integer::oelocator ! locator index for outgoing energy data
   character(66)::text

   integer::ner=1
   integer::nbw=1

   !--open the output file
   if (itype.eq.1) call openz(nout,1)
   if (itype.eq.2) call openz(-nout,1)

   !--write type-1 header block
   if (itype.eq.1) then
      if (mcnpx.eq.0) then
         write(nout,&
           '(a10,f12.6,1x,1p,e11.4,1x,a10/a70,a10)')&
           hz(1:10),aw0,tz,hd,hk,hm
      else
         write(nout,&
           '(a13,f12.6,1x,1p,e11.4,1x,a10/a70,a10)')&
           hz,aw0,tz,hd,hk,hm
      endif
      write(nout,'(4(i7,f11.0))') (izn(i),awn(i),i=1,16)
      write(nout,'(8i9)')lxs,za,nes,ntr,ntype,npixs,neixs,nxsd(1:8),tvn,&
        esz,tot,non,els,thn,mtr,lqr,lsig,sig,ixsa,ixs,jxsd(1:21)

      !--esz block
      l=1
      call advance_to_locator(nout,l,esz)
      call write_real_list(nout,l,2*nes)   ! energies and total cross section (2*nes values)

      !--els block
      if (els.ne.0) then
         call advance_to_locator(nout,l,els)
         call write_real_list(nout,l,nes)  ! elastic cross section (nes values)
      endif

      !--thn block
      if (thn.ne.0) then
         call advance_to_locator(nout,l,thn)
         call write_real_list(nout,l,nes)  ! non-elastic cross section (nes values)
      endif

      !--mtr block
      call advance_to_locator(nout,l,mtr)
      call write_integer_list(nout,l,ntr)

      !--lqr block
      call advance_to_locator(nout,l,lqr)
      call write_real_list(nout,l,ntr)

      !--lsig block
      call advance_to_locator(nout,l,lsig)
      rlocator=l
      call write_integer_list(nout,l,ntr)

      !--sig block
      call advance_to_locator(nout,l,sig)
      do i=1,ntr
         call advance_to_locator(nout,l,sig+nint(xss(rlocator))-1) ! sig=jxs(7)
         call write_integer(nout,l)
         ne=nint(xss(l))
         call write_integer(nout,l)
         call write_real_list(nout,l,ne)
         rlocator=rlocator+1
      enddo

      !--particle production blocks
      if (ntype.gt.0) then

         !--ixs arrays
         call advance_to_locator(nout,l,ixsa)
         plocator=l
         call write_integer_list(nout,l,neixs*ntype) ! IXS array (neixs values) for each IP

         !--loop over the ntype productions
         do ip=1,ntype

            ! IXS array entries
            ipt=nint(xss(plocator))
            ntrp=nint(xss(plocator+1))
            pxs=nint(xss(plocator+2))
            phn=nint(xss(plocator+3))
            mtrp=nint(xss(plocator+4))
            tyrp=nint(xss(plocator+5))
            lsigp=nint(xss(plocator+6))
            sigp=nint(xss(plocator+7))
            landp=nint(xss(plocator+8))
            andp=nint(xss(plocator+9))
            ldlwp=nint(xss(plocator+10))
            dlwp=nint(xss(plocator+11))

            !--pxs block
            call advance_to_locator(nout,l,pxs)
            call write_integer(nout,l)
            ne=nint(xss(l))
            call write_integer(nout,l)
            call write_real_list(nout,l,ne)

            !--phn block
            call advance_to_locator(nout,l,phn)
            call write_integer(nout,l)
            ne=nint(xss(l))
            call write_integer(nout,l)
            call write_real_list(nout,l,ne)

            !--mtrp block
            call advance_to_locator(nout,l,mtrp)
            call write_integer_list(nout,l,ntrp) ! MT (ntrp values)

            !--tyrp block
            call advance_to_locator(nout,l,tyrp)
            call write_integer_list(nout,l,ntrp) ! TYR (ntrp values)

            !--lsigp block
            call advance_to_locator(nout,l,lsigp)
            rlocator=l
            call write_integer_list(nout,l,ntrp) ! L (ntrp values)

            !--sigp block
            call advance_to_locator(nout,l,sigp)
            do i=1,ntrp
               call advance_to_locator(nout,l,sigp+nint(xss(rlocator))-1)

               mftype=nint(xss(l))
               call write_integer(nout,l)              ! MFTYPE
               if (mftype.eq.13) then                  ! MF=13
                  call write_integer(nout,l)           ! IE
                  ne=nint(xss(l))
                  call write_integer(nout,l)           ! NE
                  call write_real_list(nout,l,ne)      ! sigma (NE values)
               else
                  call write_integer(nout,l)           ! MTMULT
                  nr=nint(xss(l))
                  call write_integer(nout,l)           ! NR
                  if (nr.gt.0) then
                     call write_integer_list(nout,l,2*nr) ! NBT, INT (each NR values)
                  endif
                  ne=nint(xss(l))
                  call write_integer(nout,l)           ! NE
                  call write_real_list(nout,l,2*ne)    ! E, Y (each NE values)
               endif
               rlocator=rlocator+1
            enddo

            !--landp block
            call advance_to_locator(nout,l,landp)
            rlocator=l
            call write_integer_list(nout,l,ntrp)

            !--andp block
            call advance_to_locator(nout,l,andp)
            do i=1,ntrp
               nn=nint(xss(rlocator))                     ! relative locator position
               if (nn.gt.0) then
                  call advance_to_locator(nout,l,andp+nn-1)
                  ne=nint(xss(l))
                  call write_integer(nout,l)              ! NE
                  call write_real_list(nout,l,ne)         ! E (NE values)
                  ielocator=l
                  call write_integer_list(nout,l,ne)      ! L (NE values)
                  do j=1,ne
                     nn=nint(xss(ielocator))              ! relative locator position
                     call advance_to_locator(nout,l,andp+iabs(nn)-1)
                     if (nn.gt.0) then                    ! 32 equiprobable bins
                        call write_real_list(nout,l,33)
                     else if (nn.lt.0) then               ! tabulated angular
                        call write_integer(nout,l)        ! interpolation flag
                        np=nint(xss(l))
                        call write_integer(nout,l)        ! NP
                        call write_real_list(nout,l,3*np) ! CS, PDF, CDF (each NP values)
                     endif
                     ielocator=ielocator+1
                  enddo
               endif
               rlocator=rlocator+1
            enddo

            !--ldlwp block
            call advance_to_locator(nout,l,ldlwp)
            rlocator=l
            call write_integer_list(nout,l,ntrp)

            !--dlwp block
            call advance_to_locator(nout,l,dlwp)
            do i=1,ntrp
               nn=nint(xss(rlocator))                     ! relative locator position
               if (nn.gt.0) then
                  call advance_to_locator(nout,l,dlwp+nn-1)
                  lnw=1
                  do while (lnw.ne.0)
                     lnw=nint(xss(l))
                     call write_integer(nout,l)            ! LNW
                     law=nint(xss(l))
                     call write_integer(nout,l)            ! LAW
                     call write_integer(nout,l)            ! IDAT
                     nr=nint(xss(l))
                     call write_integer(nout,l)               ! NR
                     if (nr.gt.0) then
                        call write_integer_list(nout,l,2*nr)  ! NBT, INT (each NR values)
                     endif
                     ne=nint(xss(l))
                     call write_integer(nout,l)            ! NE
                     call write_real_list(nout,l,2*ne)     ! E and P (each NE values)

                     !--law 4
                     if (law.eq.4) then
                        nr=nint(xss(l))
                        call write_integer(nout,l)               ! NR
                        if (nr.gt.0) then
                           call write_integer_list(nout,l,2*nr)  ! NBT, INT (each NR values)
                        endif
                        ne=nint(xss(l))
                        call write_integer(nout,l)            ! NE
                        call write_real_list(nout,l,ne)       ! E (NE values)
                        ielocator=l
                        call write_integer_list(nout,l,ne)    ! L (NE values)
                        do j=1,ne
                           call advance_to_locator(nout,l,dlwp+nint(xss(ielocator))-1)
                           call write_integer(nout,l)         ! INTT
                           np=nint(xss(l))
                           call write_integer(nout,l)         ! NP
                           call write_real_list(nout,l,3*np)  ! Eout, PDF, CDF (each NP values)
                           ielocator=ielocator+1
                        enddo

                     !--law 44
                     else if (law.eq.44) then
                        nr=nint(xss(l))
                        call write_integer(nout,l)               ! NR
                        if (nr.gt.0) then
                           call write_integer_list(nout,l,2*nr) ! NBT, INT (each NR values)
                        endif
                        ne=nint(xss(l))
                        call write_integer(nout,l)            ! NE
                        call write_real_list(nout,l,ne)       ! E (NE values)
                        ielocator=l
                        call write_integer_list(nout,l,ne)    ! L (NE values)
                        do j=1,ne
                           call advance_to_locator(nout,l,dlwp+nint(xss(ielocator))-1)
                           call write_integer(nout,l)         ! INTT
                           np=nint(xss(l))
                           call write_integer(nout,l)         ! NP
                           call write_real_list(nout,l,5*np)  ! Eout, PDF, CDF, R, A (each NP values)
                           ielocator=ielocator+1
                        enddo

                     !--law 61
                     else if (law.eq.61) then
                        nrr=nint(xss(l))
                        call write_integer(nout,l)               ! NR
                        if (nrr.gt.0) then
                           call write_integer_list(nout,l,2*nrr) ! NBT, INT (each NR values)
                        endif
                        ne=nint(xss(l))
                        call write_integer(nout,l)            ! NE
                        call write_real_list(nout,l,ne)       ! E (NE values)
                        ielocator=l
                        call write_integer_list(nout,l,ne)    ! L (NE values)
                        do j=1,ne
                           call advance_to_locator(nout,l,dlwp+nint(xss(ielocator))-1)
                           call write_integer(nout,l)         ! INTT
                           np=nint(xss(l))
                           call write_integer(nout,l)         ! NP
                           call write_real_list(nout,l,3*np)  ! Eout, PDF, CDF (each NP values)
                           oelocator=l
                           call write_integer_list(nout,l,np) ! L (NP values)
                           do k=1,np
                              call advance_to_locator(nout,l,dlwp+nint(xss(oelocator))-1)
                              call write_integer(nout,l)         ! JJ
                              nmu=nint(xss(l))
                              call write_integer(nout,l)         ! NMU
                              call write_real_list(nout,l,3*nmu) ! Mu, PDF, CDF (each NMU values)
                              oelocator=oelocator+1
                           enddo
                           ielocator=ielocator+1
                        enddo

                     !--law 7 or 9
                     else if (law.eq.7.or.law.eq.9) then
                        nr=nint(xss(l))
                        call write_integer(nout,l)               ! NR
                        if (nr.gt.0) then
                           call write_integer_list(nout,l,2*nr) ! NBT, INT (each NR values)
                        endif
                        ne=nint(xss(l))
                        call write_integer(nout,l)            ! NE
                        call write_real_list(nout,l,2*ne)     ! E, theta (each NE values)
                        call write_real(nout,l)               ! U

                     !--law 33
                     else if (law.eq.33) then
                        call write_real(nout,l)
                        call write_real(nout,l)

                     ! unknown law
                     else
                        write(text,'(''Undefined law for dlwp block: '',i3)') law
                        call error('phnout',text,' ')
                     endif
                  enddo
               endif
               rlocator=rlocator+1
            enddo
            plocator=plocator+neixs
         !--continue loop over productions
         enddo
      endif
      call typen(0,nout,3)
      nern=0
      lrec=0
      call closz(nout)

   !--type-2 format
   else
      if (mcnpx.eq.0) then
         write(nout) hz(1:10),aw0,tz,hd,hk,hm,(izn(i),awn(i),i=1,16),&
           lxs,za,nes,ntr,ntype,npixs,neixs,nxsd(1:8),tvn,&
           esz,tot,non,els,thn,mtr,lqr,lsig,sig,ixsa,ixs,jxsd(1:21)
      else
         write(nout) hz(1:13),aw0,tz,hd,hk,hm,(izn(i),awn(i),i=1,16),&
           lxs,za,nes,ntr,ntype,npixs,neixs,nxsd(1:8),tvn,&
           esz,tot,non,els,thn,mtr,lqr,lsig,sig,ixsa,ixs,jxsd(1:21)
      endif
      ll=0
      nn=lxs
      do while (nn.gt.0)
         n=nn
         if (n.gt.ner) n=ner
         write(nout) (xss(ll+i),i=1,n)
         ll=ll+n
         nn=nn-n
      enddo
      nern=ner
      lrec=ner*nbw
      call closz(nout)
   endif

   !--output directory for mcnp
   call openz(ndir,1)
   if (mcnpx.eq.0) then
      write(ndir,&
        '(a10,f12.6,'' filename route'',i2,'' 1 '',i8,2i6,1p,e10.3)')&
        hz(1:10),aw0,itype,lxs,lrec,nern,tz
   else
      write(ndir,&
        '(a13,f12.6,'' filename route'',i2,'' 1 '',i8,2i6,1p,e10.3)')&
        hz(1:13),aw0,itype,lxs,lrec,nern,tz
   endif
   call closz(ndir)
   return
   end subroutine phnout

   subroutine phnplo(nout,hk)
   !-------------------------------------------------------------------
   ! Do standard plots of an ACE photonuclear file.
   !-------------------------------------------------------------------
   use util ! provides openz
   use acecm ! provides mtname
   ! externals
   integer::nout
   character(70)::hk
   ! internals
   integer::ipcol,iwcol,i,major,minor,it,j,mt,k,n,iaa,naa
   integer::ii,idash,icolor,ie,imt,l1,l2,ne,na,law,nb
   integer::itwo,np,ic,intt,iflag,l3,loci,nn,nd,l
   integer::ipt,ntrp,pxs,phn,mtrp,tyrp,lsigp,sigp,landp,andp,ldlwp,dlwp
   real(kr)::xmin,xmax,ymin,ymax,zmin,zmax
   real(kr)::selas,xstep,ystep,xtag,ytag,test,elas,x,y
   real(kr)::e,stot,heat,snon,xs,break,cc,pp,rat,stepm
   real(kr)::elast,ylast,ep,pd
   character(10)::name
   character(1)::qu=''''
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::delta=.01e0_kr
   real(kr),parameter::up=1.01e0_kr
   real(kr),parameter::dn=0.99e0_kr
   real(kr),parameter::eps=1.e-5_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1

   !--default colors
   ipcol=2
   iwcol=3

   !--start the viewr input text
   call openz(nout,1)
   write(nout,'(''1 2 .30'',i3,''/'')') ipcol

   !--plot lin-lin total, nonelastic, and elastic (if any)
   xmin=big
   xmax=0
   ymin=big
   ymax=-big
   do i=1,nes
      e=xss(esz-1+i)
      stot=xss(esz+nes-1+i)
      if (e.lt.xmin) xmin=e
      if (e.gt.xmax) xmax=e
      if (stot.lt.ymin) ymin=stot
      if (stot.gt.ymax) ymax=stot
      if (els.ne.0) then
         snon=xss(non-1+i)
         if (snon.lt.ymin) ymin=snon
         if (snon.gt.ymax) ymax=snon
         selas=xss(els-1+i)
         if (selas.lt.ymin) ymin=selas
         if (selas.gt.ymax) ymax=selas
      endif
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
   write(nout,'(a,''<p>rincipal cross sections'',a,''/'')') qu,qu
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
   write(nout,'(''/'')')
   write(nout,'(a,''total'',a,''/'')') qu,qu
   write(nout,'(''0/'')')
   j=0
   do i=1,nes
      e=xss(esz-1+i)
      stot=xss(esz+nes-1+i)
      write(nout,'(1p,2e14.6,''/'')') e,stot
   enddo
   write(nout,'(''/'')')
   if (els.gt.0) then
      write(nout,'(''2/'')')
      write(nout,'(''/'')')
      write(nout,'(''0 0 1/'')')
      write(nout,'(a,''nonelastic'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      j=0
      do i=1,nes
         e=xss(esz-1+i)
         snon=xss(non-1+i)
         write(nout,'(1p,2e14.6,''/'')') e,snon
      enddo
      write(nout,'(''/'')')
      write(nout,'(''3/'')')
      write(nout,'(''/'')')
      write(nout,'(''0 0 2/'')')
      write(nout,'(a,''elastic'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      j=0
      do i=1,nes
         e=xss(esz-1+i)
         elas=xss(els-1+i)
         write(nout,'(1p,2e14.6,''/'')') e,elas
      enddo
      write(nout,'(''/'')')
   endif

   !--plot lin-lin partial cross sections
   if (ntr.gt.1) then
      xmin=big
      xmax=0
      ymin=big
      ymax=-big
      do i=1,ntr
         mt=nint(xss(mtr+i-1))
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
      write(nout,'(a,''<p>artial cross sections'',a,''/'')') qu,qu
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
      mt=nint(xss(mtr+ii-1))
      call mtname(mt,name,0)
      if (name(1:1).eq.'(') then
         name(2:2)='g'
      endif
      k=nint(xss(lsig+ii-1)+sig-1)
      n=nint(xss(k+1))
      iaa=nint(xss(k))
      write(nout,'(''0 0 0'',i2,''/'')') ii-1
      write(nout,'(a,a,a,''/'')') qu,name,qu
      write(nout,'(''0/'')')
      do j=1,n
         x=xss(iaa+j-1)
         y=xss(k+2+j-1)
         write(nout,'(1p,2e14.6,''/'')') x,y
      enddo
      write(nout,'(''/'')')
      do ii=2,ntr
         mt=nint(xss(mtr+ii-1))
         call mtname(mt,name,0)
         if (name(1:1).eq.'(') then
            name(2:2)='g'
         endif
         k=nint(xss(lsig+ii-1)+sig-1)
         n=nint(xss(k+1))
         iaa=nint(xss(k))
         write(nout,'(''2/'')')
         write(nout,'(''/'')')
         idash=0
         icolor=ii-1
         do while (icolor.ge.8)
            idash=idash+1
            icolor=icolor-8
         enddo
         write(nout,'(''0 0'',2i2,''/'')') idash,icolor
         write(nout,'(a,a,a,''/'')') qu,name,qu
         write(nout,'(''0/'')')
         do j=1,n
            x=xss(iaa+j-1)
            y=xss(k+2+j-1)
            write(nout,'(1p,2e14.6,''/'')') x,y
         enddo
         write(nout,'(''/'')')
      enddo
   endif

   !--plot lin-lin heating per reaction
   xmin=big
   xmax=0
   ymin=big
   ymax=-big
   do i=1,nes
      e=xss(esz-1+i)
      heat=xss(thn-1+i)
      if (e.lt.xmin) xmin=e
      if (e.gt.xmax) xmax=e
      if (heat.lt.ymin) ymin=heat
      if (heat.gt.ymax) ymax=heat
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
   do i=1,nes
      e=xss(esz-1+i)
      heat=xss(thn-1+i)
      write(nout,'(1p,2e14.6,''/'')') e,heat
   enddo
   write(nout,'(''/'')')

   !--lin-lin plot of particle heating values
   xmin=big
   xmax=0
   ymin=big
   ymax=-big
   do i=1,ntype
      ipt=nint(xss(ixsa+12*(i-1)))
      ntrp=nint(xss(ixsa+12*(i-1)+1))
      pxs=nint(xss(ixsa+12*(i-1)+2))
      phn=nint(xss(ixsa+12*(i-1)+3))
      iaa=nint(xss(pxs))
      naa=nint(xss(pxs+1))
      do ie=1,naa
         e=xss(esz+iaa+ie-2)
         xs=xss(phn+1+ie)
         if (e.lt.xmin) xmin=e
         if (e.gt.xmax) xmax=e
         if (xs.lt.ymin) ymin=xs
         if (xs.gt.ymax) ymax=xs
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
      ipt=nint(xss(ixsa+12*(ii-1)))
      ntrp=nint(xss(ixsa+12*(ii-1)+1))
      pxs=nint(xss(ixsa+12*(ii-1)+2))
      phn=nint(xss(ixsa+12*(ii-1)+3))
      iaa=nint(xss(pxs))
      naa=nint(xss(pxs+1))
      write(nout,'(''0 0 0'',i2,''/'')') ii-1
      if (ipt.eq.1) write(nout,'(a,''neutrons'',a,''/'')') qu,qu
      if (ipt.eq.2) write(nout,'(a,''photons'',a,''/'')') qu,qu
      if (ipt.eq.9) write(nout,'(a,''protons'',a,''/'')') qu,qu
      if (ipt.eq.31) write(nout,'(a,''deuterons'',a,''/'')') qu,qu
      if (ipt.eq.32) write(nout,'(a,''tritons'',a,''/'')') qu,qu
      if (ipt.eq.33) write(nout,'(a,''he-3'',a,''/'')') qu,qu
      if (ipt.eq.34) write(nout,'(a,''alphas'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      j=0
      do i=1,naa
         e=xss(esz+iaa-2+i)
         xs=xss(phn+1+i)
         write(nout,'(1p,2e14.6,''/'')') e,xs
      enddo
      write(nout,'(''/'')')
      if (ntype.ne.1) then
         do ii=2,ntype
            ipt=nint(xss(ixsa+12*(ii-1)))
            ntrp=nint(xss(ixsa+12*(ii-1)+1))
            pxs=nint(xss(ixsa+12*(ii-1)+2))
            phn=nint(xss(ixsa+12*(ii-1)+3))
            iaa=nint(xss(pxs))
            naa=nint(xss(pxs+1))
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
            do i=1,naa
               e=xss(esz+iaa-2+i)
               xs=xss(phn+1+i)
               write(nout,'(1p,2e14.6,''/'')') e,xs
            enddo
         write(nout,'(''/'')')
         enddo
      endif
   endif

   !--lin-lin plot of particle production cross sections
   xmin=big
   xmax=0
   ymin=big
   ymax=-big
   do i=1,ntype
      ipt=nint(xss(ixsa+12*(i-1)))
      ntrp=nint(xss(ixsa+12*(i-1)+1))
      pxs=nint(xss(ixsa+12*(i-1)+2))
      phn=nint(xss(ixsa+12*(i-1)+3))
      iaa=nint(xss(pxs))
      naa=nint(xss(pxs+1))
      do ie=1,naa
         e=xss(esz+iaa+ie-2)
         xs=xss(pxs+1+ie)
         if (ipt.eq.2) xs=xs/5
         if (e.lt.xmin) xmin=e
         if (e.gt.xmax) xmax=e
         if (xs.lt.ymin) ymin=xs
         if (xs.gt.ymax) ymax=xs
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
   ipt=nint(xss(ixsa+12*(ii-1)))
   ntrp=nint(xss(ixsa+12*(ii-1)+1))
   pxs=nint(xss(ixsa+12*(ii-1)+2))
   iaa=nint(xss(pxs))
   naa=nint(xss(pxs+1))
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
   do i=1,naa
      e=xss(esz+iaa-2+i)
      xs=xss(pxs+1+i)
      if (ipt.eq.2) xs=xs/5
      write(nout,'(1p,2e14.6,''/'')') e,xs
   enddo
   write(nout,'(''/'')')
   if (ntype.ne.1) then
      do ii=2,ntype
         ipt=nint(xss(ixsa+12*(ii-1)))
         ntrp=nint(xss(ixsa+12*(ii-1)+1))
         pxs=nint(xss(ixsa+12*(ii-1)+2))
         iaa=nint(xss(pxs))
         naa=nint(xss(pxs+1))
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
         do i=1,naa
            e=xss(esz+iaa-2+i)
            xs=xss(pxs+1+i)
            if (ipt.eq.2) xs=xs/5
            write(nout,'(1p,2e14.6,''/'')') e,xs
         enddo
         write(nout,'(''/'')')
      enddo
   endif

   !--make 3-d plots of particle emission
   !--loop over each production type
   do i=1,ntype
      ipt=nint(xss(ixsa+12*(i-1)))
      ntrp=nint(xss(ixsa+12*(i-1)+1))
      pxs=nint(xss(ixsa+12*(i-1)+2))
      phn=nint(xss(ixsa+12*(i-1)+3))
      mtrp=nint(xss(ixsa+12*(i-1)+4))
      tyrp=nint(xss(ixsa+12*(i-1)+5))
      lsigp=nint(xss(ixsa+12*(i-1)+6))
      sigp=nint(xss(ixsa+12*(i-1)+7))
      landp=nint(xss(ixsa+12*(i-1)+8))
      andp=nint(xss(ixsa+12*(i-1)+9))
      ldlwp=nint(xss(ixsa+12*(i-1)+10))
      dlwp=nint(xss(ixsa+12*(i-1)+11))

      !--loop over the reactions contributing to this production
      do imt=1,ntrp
         mt=nint(xss(mtrp+imt-1))
         l1=nint(xss(lsigp+imt-1))
         l2=sigp+l1-1
         ne=nint(xss(l2+3))
         na=nint(xss(landp+imt-1))
         l1=nint(xss(ldlwp+imt-1))
         l2=dlwp+l1-1
         law=nint(xss(l2+1))
         call mtname(mt,name,0)
         if (name(1:1).eq.'(') then
            name(2:2)='g'
         endif

         !--check for an angular distribution
         if (na.gt.0) then
            na=na+andp-1
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
                     k=nint(abs(xss(nb+ie)))+andp-1
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
                     k=nint(abs(xss(nb+ie)))+andp-1
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

         !--make a 3d plot of the particle emission for this reaction
         if (law.eq.4.or.law.eq.44.or.law.eq.61) then
            l3=dlwp+nint(xss(l2+2))-1
            ne=nint(xss(l3+1))
            zmin=1000
            zmax=0
            ymin=xss(l3+2)
            ymax=xss(l3+2+ne-1)
            test=3
            if (ymax.gt.test*xss(l3+2+ne-2)) ymax=xss(l3+2+ne-2)
            do ie=2,ne
               e=xss(l3+2+ie-1)
               if (e.le.up*ymax) then
                  loci=nint(xss(l3+2+ne+ie-1))+dlwp-1
                  intt=nint(xss(loci))
                  nn=nint(xss(loci+1))
                  if (nn.gt.2) then
                     loci=loci+1
                     !--skip first two point that may be pseudo-threshold
                     do j=3,nn
                        ep=xss(loci+j)
                        pd=xss(loci+nn+j)
                        if (pd.lt.zmin) zmin=pd
                        if (pd.gt.zmax) zmax=pd
                     enddo
                  endif
               endif
            enddo
            if (zmax.eq.0) go to 100
            zmin=zmax/10000
            if (zmax.eq.0) cycle
            call ascll(zmin,zmax)
            xmin=1000
            xmax=0
            do ie=1,ne
               e=xss(l3+2+ie-1)
               if (e.le.up*ymax) then
                  loci=nint(xss(l3+2+ne+ie-1))+dlwp-1
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
            call ascle(2,xmin,xmax,major,minor)
            xstep=(xmax-xmin)/major
            call ascle(4,ymin,ymax,major,minor)
            ystep=(ymax-ymin)/major
            write(nout,'(''1'',i3,''/'')') iwcol
            it=1
            do j=1,70
               if (hk(j:j).ne.' ') it=j
            enddo
            l=len_trim(name)
            write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
            if (ipt.eq.1) write(nout,&
              '(a,''neutrons from '',a,a,''/'')') qu,name(1:l),qu
            if (ipt.eq.2) write(nout,&
              '(a,''photons from '',a,a,''/'')') qu,name(1:l),qu
            if (ipt.eq.9) write(nout,&
              '(a,''protons from '',a,a,''/'')') qu,name(1:l),qu
            if (ipt.eq.31) write(nout,&
              '(a,''deuterons from '',a,a,''/'')') qu,name(1:l),qu
            if (ipt.eq.32) write(nout,&
              '(a,''tritons from '',a,a,''/'')') qu,name(1:l),qu
            if (ipt.eq.33) write(nout,&
              '(a,''he3s from '',a,a,''/'')') qu,name(1:l),qu
            if (ipt.eq.34) write(nout,&
              '(a,''alphas from '',a,a,''/'')') qu,name(1:l),qu
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
                  loci=nint(xss(l3+2+ne+ie-1))+dlwp-1
                  intt=nint(xss(loci))
                  nd=nint(xss(loci))/10
                  nn=nint(xss(loci+1))
                  if (nn.gt.2) then
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
                              write(nout,'(1p,2e14.6,''/'')')&
                                ep,ylast
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
               endif
            enddo
            write(nout,'(''/'')')
         endif
  100    continue
      enddo
   enddo

   !--end the plotr input text
   write(nout,'(''99/'')')
   write(nout,'(''stop'')')
   return
   end subroutine phnplo

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
   return
   end subroutine ascll

   subroutine copynubar(scr,fnubar,jscr)
   !-------------------------------------------------------------------
   ! Copy the content of the nubar table to the scr array
   !-------------------------------------------------------------------
   ! externals
   real(kr)::scr(*)
   real(kr)::fnubar(*)
   integer::jscr
   ! internals
   integer::ii,nrr,npp

   ! replace the yield with the nubar
   nrr=nint(fnubar(5))
   npp=nint(fnubar(6))
   scr(5)=nrr
   scr(6)=npp
   do ii=1,nrr
      scr(5+2*ii)=fnubar(5+2*ii)
      scr(6+2*ii)=fnubar(6+2*ii)
   enddo
   do ii=1,npp
      scr(5+2*nrr+2*ii)=fnubar(5+2*nrr+2*ii)
      scr(6+2*nrr+2*ii)=fnubar(6+2*nrr+2*ii)
   enddo

   ! set the scr array index to the appropriate value
   jscr=6+2*nrr+2*npp+1

   return
   end subroutine copynubar

end module acepn
