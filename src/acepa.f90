module acepa
   ! provides ACE photo-atomic routines for acer
   use locale
   use acecm, only: xss,nxss
   implicit none
   private

   !--Public routines
   public acepho,phofix

   !--Private global variables

   ! ace header data for photoatomic format
   character(13)::hz
   character(10)::hd
   character(10)::hm
   real(kr)::aw0,tz

   ! parameters for photoatomic nsx block
   integer::len2,z,nes,nflo,nxsd(12)

   ! parameters for photoatomic jxs block
   integer::eszg,jinc,jcoh,jflo,lhnm,jxsd(27)

   ! parameter for scratch array
   integer,parameter::nwscr=50000

contains

   subroutine acepho(nin,nlax,nace,ndir,matd,mcnpx,iprint,itype,&
     suff,hk,izn,awn)
   !-------------------------------------------------------------------
   ! Prepare ACE photo-atomic files.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util ! provides openz,mess,closz,dater
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nlax,nace,ndir,matd,mcnpx,iprint,itype
   real(kr)::suff
   integer::izn(16)
   real(kr)::awn(16)
   character(70)::hk
   ! internals
   integer::nb,nw,l,iza,idis,iinc,icoh,iabs,ipair,next
   integer::i,ip,ir,nr,np,iz
   real(kr)::e,enext,s,v,vnext,v2,en,heat,siginc,zaid,tot
   character(8)::hdt
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(21),parameter::vi=(/0.e0_kr,.005e0_kr,&
     .01e0_kr,.05e0_kr,.1e0_kr,.15e0_kr,.2e0_kr,.3e0_kr,.4e0_kr,&
     .5e0_kr,.6e0_kr,.7e0_kr,.8e0_kr,.9e0_kr,1.e0_kr,1.5e0_kr,&
     2.e0_kr,3.e0_kr,4.e0_kr,5.e0_kr,8.e0_kr/)
   real(kr),dimension(55),parameter::vc=(/0.e0_kr,.01e0_kr,&
     .02e0_kr,.03e0_kr,.04e0_kr,.05e0_kr,.06e0_kr,.08e0_kr,.10e0_kr,&
     .12e0_kr,.15e0_kr,.18e0_kr,.20e0_kr,.25e0_kr,.30e0_kr,.35e0_kr,&
     .40e0_kr,.45e0_kr,.50e0_kr,.55e0_kr,.60e0_kr,.70e0_kr,.80e0_kr,&
     .90e0_kr,1.0e0_kr,1.1e0_kr,1.2e0_kr,1.3e0_kr,1.4e0_kr,1.5e0_kr,&
     1.6e0_kr,1.7e0_kr,1.8e0_kr,1.9e0_kr,2.0e0_kr,2.2e0_kr,2.4e0_kr,&
     2.6e0_kr,2.8e0_kr,3.0e0_kr,3.2e0_kr,3.4e0_kr,3.6e0_kr,3.8e0_kr,&
     4.0e0_kr,4.2e0_kr,4.4e0_kr,4.6e0_kr,4.8e0_kr,5.0e0_kr,5.2e0_kr,&
     5.4e0_kr,5.6e0_kr,5.8e0_kr,6.0e0_kr/)
   real(kr),parameter::emin=1.e3_kr
   real(kr),parameter::emev=1.e6_kr
   real(kr),parameter::emax=1.01e11_kr
   real(kr),parameter::epair=1.022e0_kr
   real(kr),parameter::zero=0

   nxsd=0
   jxsd=0
   xss=0

   !--allocate scratch storage
   allocate(scr(nwscr))

   !--assign input file
   call openz(nin,0)

   !--determine what endf version is being used
   call repoz(nin)
   call tpidio(nin,0,0,scr,nb,nw)
   call contio(nin,0,0,scr,nb,nw)
   call contio(nin,0,0,scr,nb,nw)
   if (n1h.ne.0) then
      iverf=4
   else if (n2h.eq.0) then
      iverf=5
   else
      iverf=6
   endif
   write(nsyso,'(/'' using endf-'',i1,'' format'')') iverf

   !--locate and store energy grid of total cross section
   l=0
   call findf(matd,23,501,nin)
   call contio(nin,0,0,scr,nb,nw)
   iza=nint(scr(1))
   aw0=scr(2)
   e=0
   call gety1(e,enext,idis,s,nin,scr)
   enext=emin
   do while (enext.lt.emax)
      e=sigfig(enext,9,0)
      if (idis.ne.0) then
         e=sigfig(e,9,-1)
         call gety1(e,enext,idis,s,nin,scr)
         l=l+1
         xss(l)=e
         e=sigfig(e,9,+2)
      endif
      call gety1(e,enext,idis,s,nin,scr)
      l=l+1
      xss(l)=e
   enddo
   call tosend(nin,0,0,scr)
   nes=l

   !--define cross section locators
   eszg=1
   iinc=eszg+nes
   icoh=iinc+nes
   iabs=icoh+nes
   ipair=iabs+nes
   next=ipair+nes

   !--loop over photon reactions on file 23
   mfh=23
   do while (mfh.ne.0)
      call contio(nin,0,0,scr,nb,nw)
      l=0
      if (mfh.ne.0) then
         if (mth.eq.502) l=icoh
         if (mth.eq.504) l=iinc
         if (mth.eq.516) l=ipair
         if (iverf.lt.6.and.mth.eq.602) l=iabs
         if (iverf.ge.6.and.mth.eq.522) l=iabs
         if (l.ne.0) then
            e=0
            call gety1(e,enext,idis,s,nin,scr)
            do i=1,nes
               e=xss(eszg-1+i)
               call gety1(e,enext,idis,s,nin,scr)
               xss(l-1+i)=s
            enddo
         endif
         call tosend(nin,0,0,scr)
      endif
   enddo

   !--fix up energy units
   do i=1,nes
      xss(eszg-1+i)=xss(eszg-1+i)/emev
   enddo

   !--set number of fluorescence lines
   iz=matd/100
   if (iz.lt.12) then
      nflo=0
   else if (iz.lt.20) then
      nflo=2
   else if (iz.lt.31) then
      nflo=4
   else if (iz.lt.37) then
      nflo=5
   else
      nflo=6
   endif

   !--assign locators for the other data blocks
   jinc=next
   jcoh=jinc+21
   jflo=jcoh+110
   lhnm=jflo+4*nflo
   len2=lhnm+nes-1

   !--process the coherent form factors
   call findf(matd,27,502,nin)
   call contio(nin,0,0,scr,nb,nw)
   z=nint(scr(1)/1000)
   call tab1io(nin,0,0,scr,nb,nw)
   l=1+nw
   do while (nb.ne.0)
      if (l.gt.nwscr) call error('acepho',&
              'storage exceeded for the coherent form factors',' ')
      call moreio(nin,0,0,scr(l),nb,nw)
      l=l+nw
   enddo
   ip=2
   ir=1
   do i=1,55
      v=vc(i)
      call terpa(s,v,vnext,idis,scr,ip,ir)
      xss(jcoh+54+i)=s
   enddo
   nr=nint(scr(5))
   np=nint(scr(6))
   l=1+6+2*nr
   do i=1,np
      scr(l)=scr(l)**2
      scr(l+1)=scr(l+1)**2/z**2
      l=l+2
   enddo
   do i=1,55
      ir=1
      ip=2
      v2=vc(i)**2
      call intega (xss(jcoh-1+i),zero,v2,scr,ip,ir)
   enddo

   !--process the incoherent scattering function
   call findf(matd,27,504,nin)
   call contio(nin,0,0,scr,nb,nw)
   call tab1io(nin,0,0,scr,nb,nw)
   l=1+nw
   do while (nb.ne.0)
      if (l.gt.nwscr) call error('acepho',&
              'storage exceeded for the incoherent scattering function',' ')
      call moreio(nin,0,0,scr(l),nb,nw)
      l=l+nw
   enddo
   ip=2
   ir=1
   do i=1,21
      v=vi(i)
      call terpa (s,v,vnext,idis,scr,ip,ir)
      xss(jinc-1+i)=s
   enddo
   do i=1,nes
      call iheat(xss(i)*emev,en,idis,scr,heat,siginc)
      xss(lhnm-1+i)=xss(iinc-1+i)*heat/emev
   enddo

   !--process photoelectric data
   !--for fluorescence photons
   !--and contribution to heating
   do i=1,nes
      xss(lhnm-1+i)=xss(lhnm-1+i)+xss(eszg-1+i)*xss(iabs-1+i)
   enddo
   if (nlax.eq.0) then
      call mess('acepho','no atomic relaxation data',&
        'fluorescence data not processed')
   else if (nflo.gt.0) then
      call alax(nin,nlax,matd,xss(jflo),xss(eszg),xss(iabs),xss(lhnm),scr,nes)
   endif

   !--add in heating contribution from pair production
   do i=1,nes
      xss(lhnm-1+i)=xss(lhnm-1+i)+xss(ipair-1+i)*(xss(eszg-1+i)-epair)
   enddo

   !--convert heating to a per collision basis
   do i=1,nes
      tot=xss(iinc-1+i)+xss(icoh-1+i)+xss(iabs-1+i)+xss(ipair-1+i)
      xss(lhnm-1+i)=xss(lhnm-1+i)/tot
   enddo

   !--finished with photon data
   call closz(nin)
   deallocate(scr)
   zaid=iza+suff
   if (mcnpx.eq.0) then
      write(hz,'(f9.2,''p'')') zaid
   else
      write(hz,'(f10.3,''pp '')') zaid
   endif
   tz=0
   call dater(hdt)
   hd='  '//hdt
   write(hm,'(''   mat'',i4)') matd
   if (iprint.gt.0) call phoprt(hk)
   call phoout(itype,nace,ndir,mcnpx,hk,izn,awn)
   return
   end subroutine acepho

   subroutine phofix(itype,nin,nout,ndir,iprint,nplot,mcnpx,&
     suff,nxtra,hk,izn,awn)
   !-------------------------------------------------------------------
   ! Print and or edit an ACE photoatomic file.
   !-------------------------------------------------------------------
   use util ! provides closz
   ! externals
   integer::itype,nin,nout,ndir,iprint,nplot,mcnpx,nxtra
   real(kr)::suff
   integer::izn(16)
   real(kr)::awn(16)
   character(70)::hk
   ! internals
   integer::l,j,n,max,iza,i
   real(kr)::zaid
   integer::izo(16)
   real(kr)::awo(16)
   character(70)::hko
   character(3)::ht
   character(9)::str
   real(kr),parameter::zero=0

   integer,parameter::ner=1

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
      read(nin,'(8i9)')&
        len2,z,nes,nflo,nxsd(1:12),&
        eszg,jinc,jcoh,jflo,lhnm,jxsd(1:27)
      n=(len2+3)/4
      l=0
      do i=1,n
         read (nin,'(4e20.0)') (xss(l+j),j=1,4)
         l=l+4
      enddo

   !--read type 2 ace format file
   else if (itype.eq.2) then
      if (mcnpx.eq.0) then
        read(nin) hz(1:10),aw0,tz,hd,hko,hm,(izo(i),awo(i),i=1,16),&
          len2,z,nes,nflo,nxsd(1:12),&
          eszg,jinc,jcoh,jflo,lhnm,jxsd(1:27)
      else
        read(nin) hz(1:13),aw0,tz,hd,hko,hm,(izo(i),awo(i),i=1,16),&
          len2,z,nes,nflo,nxsd(1:12),&
          eszg,jinc,jcoh,jflo,lhnm,jxsd(1:27)
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
   if (iprint.gt.0) call phoprt(hk)
   if (nout.gt.0) call phoout(itype,nout,ndir,mcnpx,hk,izn,awn)

   return
   end subroutine phofix

   subroutine iheat(e,enext,idis,a,heat,siginc)
   !-------------------------------------------------------------------
   ! Compute incoherent heating and cross section.
   !-------------------------------------------------------------------
   use endf ! provides terpa
   ! externals
   integer::idis
   real(kr)::e,enext,a(*),heat,siginc
   ! internals
   integer::ip,ir,idone,nq,iq
   real(kr)::zz,enow,enowi,enow2,pnow,xnow,xzz,q2m,snow,xnext
   real(kr)::arg,ebar,q2,unext,pnext,px,aq,bq,uq,wq,pnowi,unow
   real(kr)::rm2,rm,rt,dk,fact
   real(kr),dimension(6),parameter::qp6=(/-1.e0_kr,&
     -.76505532e0_kr,-.28523152e0_kr,.28523152e0_kr,.76505532e0_kr,&
     1.e0_kr/)
   real(kr),dimension(6),parameter::qw6=(/.06666667e0_kr,&
     .37847496e0_kr,.55485838e0_kr,.55485838e0_kr,.37847496e0_kr,&
     .06666667e0_kr/)
   real(kr),dimension(10),parameter::qp10=(/-1.e0_kr,&
      -.9195339082e0_kr,-.7387738651e0_kr,-.4779249498e0_kr,&
      -.1652789577e0_kr,.1652789577e0_kr,.4779249498e0_kr,&
      .7387738651e0_kr,.9195339082e0_kr,1.e0_kr/)
   real(kr),dimension(10),parameter::qw10=(/.0222222222e0_kr,&
     .1333059908e0_kr,.2248893420e0_kr,.2920426836e0_kr,&
     .3275397612e0_kr,.3275397612e0_kr,.2920426836e0_kr,&
     .2248893420e0_kr,.1333059908e0_kr,.0222222222e0_kr/)
   real(kr),parameter::c2=0.249467e0_kr
   real(kr),parameter::c3=1.95693e-6_kr
   real(kr),parameter::c4=0.0485262e0_kr
   real(kr),parameter::c5=20.60744e0_kr
   real(kr),parameter::rndoff=1.0000001e0_kr
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::one=1

   !--photon incoherent scattering and heating.
   zz=a(2)
   enow=c3*e
   enowi=1/enow
   enow2=enow*enow
   pnow=enow
   xnow=0
   xzz=c5*sqrt(enow/500)
   q2m=(2*enow*(1+enow)/(1+2*enow))**2
   if (xzz.le.zz) then
      ip=2
      ir=1
      call terpa(snow,xnow,xnext,idis,a(1),ip,ir)
   else
      snow=zz
      xnext=big
   endif
   arg=0
   siginc=0
   ebar=0

   !--loop over panels defined by form factor breaks.
   !--also limit panels by the fractional energy change
   idone=0
   do while (idone.eq.0)
      do while (idone.eq.0)
         q2=(c4*xnext)**2
         unext=-1
         if (q2.gt.q2m) then
            idone=1
         else
            unext=1-((1-q2*enowi)-sqrt(1+q2))/(q2-enow2-2*enow)
            if (unext.lt.-one) unext=-1
            if (unext.lt.0.99999) then
               idone=1
            else
               xnow=xnext*rndoff
               call terpa(snow,xnow,xnext,idis,a(1),ip,ir)
            endif
         endif
      enddo
      idone=0
      pnext=enow/(1+2*enow)
      px=enow/(1+enow*(1-unext))
      if (px.gt.pnext) pnext=px
      if (pnext.lt.pnow/2) pnext=pnow/2
      if (pnext.gt.pnow/rndoff) pnext=pnow/rndoff
      aq=(pnext+pnow)/2
      bq=(pnext-pnow)/2
      nq=10
      do iq=1,nq
         if (nq.eq.6) then
            uq=aq+bq*qp6(iq)
            wq=-c2*bq*qw6(iq)
         else if (nq.eq.10) then
            uq=aq+bq*qp10(iq)
            wq=-c2*bq*qw10(iq)
         endif
         pnow=uq
         if (iq.gt.1) then
            pnowi=1/pnow
            unow=1+enowi-pnowi
            if (xzz.le.zz) then
               rm2=(1-unow)/2
               rm=sqrt(rm2)
               rt=1+2*enow*rm2
               xnow=c5*2*enow*rm*sqrt(rt+enow2*rm2)/rt
               xnow=xnow*rndoff
               call terpa(snow,xnow,xnext,idis,a(1),ip,ir)
            endif
            dk=unow-1
            fact=snow*(enow*pnowi+pnow*enowi+dk*(2+dk))/enow2
            arg=fact
         endif
         siginc=siginc+wq*arg
         ebar=ebar+wq*arg*pnow/c3
      enddo
      if (unow.lt.-.9999) idone=1
   enddo

   !--compute final results
   ebar=ebar/siginc
   heat=e-ebar
   enext=big
   nq=nq+2
   return
   end subroutine iheat

   subroutine alax(nin,nlax,matd,fluor,ener,abs,heat,a,nes)
   !-------------------------------------------------------------------
   ! Generate fluorescence data in the Cashwell-Everett format.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util ! provides closz,error
   use endf ! provides endf parameters and routines
   ! externals
   integer::nin,nlax,matd,nes
   real(kr)::fluor(*)
   real(kr),dimension(nes)::ener,abs,heat
   real(kr)::a(*)
   ! internals
   integer::loc(50)
   integer::iz,nw,nb,nss,ll,iss,ntr,idis,i,jj,n,mm
   real(kr)::e,en,sig,slo,shi,ek,rhok,sum1,sum2
   real(kr)::el2,pl2,el3,pl3,tot,y,phi,rholt,elav,denom
   real(kr)::wt,ylt,flt,sum11,sum12,sum21,sum22,phik
   real(kr)::enl(3),rhol(3),wtl(3)
   real(kr),parameter::dn=.9999e0_kr
   real(kr),parameter::up=1.0001e0_kr
   real(kr),parameter::emev=1.e6_kr

   !--charge for desired material
   iz=matd/100

   !--read in the atomic relaxation file for the desired material
   call repoz(nin)
   call openz(nlax,0)
   call tpidio(nlax,0,0,a,nb,nw)
  110 call contio(nlax,0,0,a,nb,nw)
   if (math.gt.0) go to 120
   call error('alax','mat not found',' ')
  120 if (math.eq.matd) go to 130
   call tomend(nlax,0,0,a)
   go to 110
  130 call tofend(nlax,0,0,a)
   call contio(nlax,0,0,a,nb,nw)
   nss=n1h
   ll=1
   do iss=1,nss
      loc(iss)=ll
      if (ll.gt.nwscr) call error('alax',&
              'storage exceeded for the atomic relaxation data',' ')
      call listio(nlax,0,0,a(ll),nb,nw)
      ntr=n2h
      ll=ll+nw
      do while (nb.ne.0)
         if (ll.gt.nwscr) call error('alax',&
                 'storage exceeded for the atomic relaxation data',' ')
         call moreio(nlax,0,0,a(ll),nb,nw)
         ll=ll+nw
      enddo
   enddo

   !--read in the photoionization cross section for the material
   call openz(nin,0)
   call tpidio(nin,0,0,a(ll),nb,nw)
  210 call contio(nin,0,0,a(ll),nb,nw)
   if (math.gt.0) go to 220
   call error('spect','mat not found',' ')
  220 if (math.eq.matd) go to 230
   call tomend(nin,0,0,a(ll))
   go to 210
  230 call tofend(nin,0,0,a(ll))
   call findf(matd,23,522,nin)
   call contio(nin,0,0,a(ll),nb,nw)
   e=0
   call gety1(e,en,idis,sig,nin,a(ll))

   !--for z>30, get the l1, l2, and l3 edges and jumps
   if (iz.gt.30) then
      do i=1,3
         jj=loc(5-i)
         enl(4-i)=a(jj+6)
         e=dn*enl(4-i)
         call gety1(e,en,idis,slo,nin,a(ll))
         e=up*enl(4-i)
         call gety1(e,en,idis,shi,nin,a(ll))
         rhol(4-i)=slo/shi
      enddo
   endif

   !--get the energy and jump of the k edge
   ek=a(7)
   e=dn*ek
   call gety1(e,en,idis,slo,nin,a(ll))
   e=up*ek
   call gety1(e,en,idis,shi,nin,a(ll))
   rhok=slo/shi

   !--case of 11<z<20
   if (iz.gt.11.and.iz.lt.20) then

      !--average all the transitions to the k shell
      n=nint(a(6))
      sum1=0
      sum2=0
      do i=1,n
         jj=8+6*i
         if (nint(a(jj)).eq.0) then
            sum1=sum1+a(jj+2)
            sum2=sum2+a(jj+1)*a(jj+2)
         endif
      enddo
      sum2=sum2/sum1

      !--subtract the photon energy from the heating
      do i=1,nes
         if (ener(i).gt.ek/emev) heat(i)=heat(i)-abs(i)*sum1*sum2/emev
      enddo

      !--store the results
      fluor(1)=ek/emev
      fluor(2)=ek/emev
      fluor(3)=rhok
      fluor(4)=1
      fluor(5)=0
      fluor(6)=(1-rhok)*sum1
      fluor(7)=0
      fluor(8)=sum2/emev

   !--case of z>19 and z<31
   else if (iz.gt.19.and.iz.lt.31) then

      !--extract l2, l3, and total for higher shells
      n=nint(a(6))
      sum1=0
      sum2=0
      do i=1,n
         jj=8+6*i
         if (nint(a(jj)).eq.0.and.nint(a(jj-1)).eq.3) then
            el2=a(jj+1)
            pl2=a(jj+2)
         else if (nint(a(jj)).eq.0.and.nint(a(jj-1)).eq.4) then
            el3=a(jj+1)
            pl3=a(jj+2)
         else if (nint(a(jj)).eq.0.and.nint(a(jj-1)).gt.4) then
            sum1=sum1+a(jj+2)
            sum2=sum2+a(jj+1)*a(jj+2)
         endif
      enddo
      sum2=sum2/sum1

      !--subtract the photon energy from the heating
      do i=1,nes
         if (ener(i).gt.ek/emev) heat(i)=heat(i)-abs(i)&
           *(sum1*sum2+el2*pl2+el3*pl3)/emev
      enddo

      !--store the results
      tot=(pl2+pl3+sum1)/(1-rhok)
      y=0
      phi=rhok
      fluor(1)=ek/emev
      fluor(5)=phi
      fluor(9)=y
      fluor(13)=0
      phi=phi+pl3/tot
      y=y+(1-rhok)*pl3
      fluor(2)=ek/emev
      fluor(6)=phi
      fluor(10)=y
      fluor(14)=el3/emev
      phi=phi+pl2/tot
      y=y+(1-rhok)*pl2
      fluor(3)=ek/emev
      fluor(7)=phi
      fluor(11)=y
      fluor(15)=el2/emev
      phi=1
      y=y+(1-rhok)*sum1
      fluor(4)=ek/emev
      fluor(8)=phi
      fluor(12)=y
      fluor(16)=sum2/emev

   !--all other z values
   else
      rholt=rhol(1)*rhol(2)*rhol(3)
      elav=(enl(1)+enl(2)+enl(3))/3
      wtl(1)=1/rhol(1)
      wtl(2)=wtl(1)/rhol(2)
      wtl(3)=wtl(2)/rhol(3)
      denom=wtl(3)-1
      wtl(3)=(wtl(3)-wtl(2))/denom
      wtl(2)=(wtl(2)-wtl(1))/denom
      wtl(1)=(wtl(1)-1)/denom

      !--compute the average yield and energy for l fluorescence
      sum1=0
      sum2=0
      do iss=2,4
         jj=loc(iss)
         n=nint(a(jj+5))
         wt=wtl(iss-1)
         do i=1,n
            if (nint(a(jj+7+6*i)).eq.0) then
               sum1=sum1+a(jj+9+6*i)*wt
               sum2=sum2+a(jj+8+6*i)*a(jj+9+6*i)*wt
            endif
         enddo
      enddo
      sum2=sum2/sum1
      ylt=sum1
      flt=sum2
      if (flt.gt.enl(1)) then
         write(nsyso,'('' L edge problem'')')
         write(nsyso,'(1x,3f10.4)') flt,enl(1),elav
      endif

      !--extract kalpha1, kalpha2, kbeta1, and kbeta2
      n=nint(a(6))
      sum11=0
      sum12=0
      sum21=0
      sum22=0
      do i=1,n
         jj=8+6*i
         if (nint(a(jj)).eq.0) then
            mm=nint(a(jj-1))
            if (mm.eq.3) then
               el2=a(jj+1)
               pl2=a(jj+2)
            else if (mm.eq.4) then
               el3=a(jj+1)
               pl3=a(jj+2)
            else if (mm.ge.5.and.mm.le.9) then
               sum11=sum11+a(jj+2)
               sum12=sum12+a(jj+1)*a(jj+2)
            else if (mm.ge.10.and.mm.le.16) then
               sum21=sum21+a(jj+2)
               sum22=sum22+a(jj+1)*a(jj+2)
            endif
         endif
      enddo
      if (iz.ge.37) then
         sum22=sum22/sum21
      else
         sum11=sum11+sum21
         sum12=sum12+sum22
         sum21=0
      endif
      sum12=sum12/sum11

      !--subtract the photon energy from the heating
      do i=1,nes
         if (ener(i).gt.elav/emev) heat(i)=heat(i)-abs(i)*ylt*flt/emev
         if (ener(i).gt.ek/emev) heat(i)=heat(i)&
           -abs(i)*(sum11*sum12+sum21*sum22+el2*pl2+el3*pl3)/emev
      enddo

      !--store the results
      n=5
      if (iz.gt.36) n=6
      fluor(1)=elav/emev
      fluor(1+n)=rholt
      fluor(1+2*n)=0
      fluor(1+3*n)=0
      y=(1-rholt)*ylt
      fluor(2)=elav/emev
      fluor(2+n)=1
      fluor(2+2*n)=y
      fluor(2+3*n)=flt/emev
      phi=1/rhok
      phik=phi-1
      tot=(pl2+pl3+sum11+sum21)/phik
      phi=1
      phi=phi+pl3/tot
      y=y+phik*pl3
      fluor(3)=ek/emev
      fluor(3+n)=phi
      fluor(3+2*n)=y
      fluor(3+3*n)=el3/emev
      phi=phi+pl2/tot
      y=y+phik*pl2
      fluor(4)=ek/emev
      fluor(4+n)=phi
      fluor(4+2*n)=y
      fluor(4+3*n)=el2/emev
      phi=phi+sum11/tot
      y=y+phik*sum11
      fluor(5)=ek/emev
      fluor(5+n)=phi
      fluor(5+2*n)=y
      fluor(5+3*n)=sum12/emev
      if (iz.ge.37) then
         phi=phi+sum21/tot
         y=y+phik*sum21
         fluor(6)=ek/emev
         fluor(6+n)=phi
         fluor(6+2*n)=y
         fluor(6+3*n)=sum22/emev
      endif

   endif
   return
   end subroutine alax

   subroutine phoprt(hk)
   !-------------------------------------------------------------------
   ! Print ACE photon interaction data from memory.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   ! externals
   character(70)::hk
   ! internals
   integer::ieg,iinc,icoh,iabs,ipair,ihtng,i,j
   real(kr)::x
   character(14)::col(6)
   character(14)::blank='             '
   real(kr),parameter::zero=0

   !--print information block
   write(nsyso,'(''1''////////&
     &38x,''zaid'',1x,a13/39x,''awr'',f10.3/&
     &6x,''***********************'',9x,''temp'',1p,e10.2/&
     &6x,''*                     *'',9x,''date'',a10/&
     &6x,''*     photo-atomic    *'',10x,''mat'',a10/&
     &6x,''*                     *''/&
     &6x,''*   ace format file   *'',9x,''len2'',i10/&
     &6x,''*                     *'',12x,''z'',i10/&
     &6x,''*     processed by    *'',10x,''nes'',i10/&
     &6x,''*                     *'',9x,''nflo'',i10/&
     &6x,''*        njoy         *''/&
     &6x,''*                     *'',9x,''eszg'',i10/&
     &6x,''***********************'',9x,''jinc'',i10/&
     &38x,''jcoh'',i10/38x,''jflo'',i10/38x,''lhnm'',i10/&
     &//6x,''hk--- '',a70///)')&
     hz,aw0,tz,hd,hm,len2,z,nes,nflo,eszg,jinc,&
     jcoh,jflo,lhnm,hk

   !--print eszg block
   ieg=eszg-1
   iinc=ieg+nes
   icoh=iinc+nes
   iabs=icoh+nes
   ipair=iabs+nes
   ihtng=lhnm-1
   do i=1,nes
      if (mod(i,57).eq.1) write(nsyso,'(''1''/&
        &''     i'',8x,''energy'',4x,''incoherent'',&
        &6x,''coherent'',4x,''photoelect'',5x,''pair prod'',&
        &7x,''heating''/2x,''----'',4x,''----------'',4x,&
        &''----------'',4x,''----------'',4x,''----------'',&
        &4x,''----------'',4x,''----------'')')
      col(1)=blank
      x=xss(ieg+i)
      if (x.ne.zero) write(col(1),'(1p,e14.4)') x
      col(2)=blank
      x=xss(iinc+i)
      if (x.ne.zero) write(col(2),'(1p,e14.4)') x
      col(3)=blank
      x=xss(icoh+i)
      if (x.ne.zero) write(col(3),'(1p,e14.4)') x
      col(4)=blank
      x=xss(iabs+i)
      if (x.ne.zero) write(col(4),'(1p,e14.4)') x
      col(5)=blank
      x=xss(ipair+i)
      if (x.ne.zero) write(col(5),'(1p,e14.4)') x
      col(6)=blank
      x=xss(ihtng+i)
      if (x.ne.zero) write(col(6),'(1p,e14.4)') x
      write(nsyso,'(1x,i5,6a14)') i,col
   enddo
   write(nsyso,'(/&
     &'' the above numbers are stored as logs on the ace file'')')

   !--print the incoherent scattering function
   write(nsyso,'(//&
     &'' incoherent scattering function''/&
     &'' ------------------------------'')')
   write(nsyso,'(1x,1p,6e14.4)') (xss(jinc-1+i),i=1,21)

   !--print the coherent from factors
   write(nsyso,'(//&
     &'' coherent sampling integral''/&
     &'' --------------------------'')')
   write(nsyso,'(1x,1p,6e14.4)') (xss(jcoh-1+i),i=1,55)
   write(nsyso,'(//&
     &'' coherent scattering function''/&
     &'' ----------------------------'')')
   write(nsyso,'(1x,1p,6e14.4)') (xss(jcoh+54+i),i=1,55)

   !--print the fluorescence data
   if (nflo.gt.0) then
      write(nsyso,'(//&
        &'' fluorescence data''/&
        &'' -----------------'')')
      write(nsyso,'(/&
        &''          edge       phi         y           f''/&
        &''    ----------   -------   -------   ---------'')')
      do i=1,nflo
         write(nsyso,'(3x,f11.7,f10.4,f10.4,2x,f10.6)')&
           (xss(jflo+i-1+nflo*(j-1)),j=1,4)
      enddo
   endif
   return
   end subroutine phoprt

   subroutine phoout(itype,nout,ndir,mcnpx,hk,izn,awn)
   !-------------------------------------------------------------------
   ! Write photo-atomic ACE data to output and directory files.
   !-------------------------------------------------------------------
   use util  ! provides openz,closz,error
   use acecm ! provides write routines
   ! externals
   integer::itype,nout,ndir,mcnpx
   integer::izn(16)
   real(kr)::awn(16)
   character(70)::hk
   ! internals
   integer::l,n,i,ll,nn,nern,lrec

   integer::ner=1
   integer::nbw=1

   !--branch according to ace file type
   if (itype.eq.1) call openz(nout,1)
   if (itype.eq.2) call openz(-nout,1)

   !--type 1

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
      write(nout,'(8i9)')&
        len2,z,nes,nflo,nxsd(1:12),&
        eszg,jinc,jcoh,jflo,lhnm,jxsd(1:27)

      !--eszg block
      l=eszg
      n=5*nes
      do i=1,n
         if (xss(l).ne.0.) xss(l)=log(xss(l))
         call typen(l,nout,2)
         l=l+1
      enddo

      !--form factor blocks
      n=21+2*55
      l=jinc
      do i=1,n
         call typen(l,nout,2)
         l=l+1
      enddo

      !--fluorescence data block
      l=jflo
      if (nflo.ne.0) then
         do i=1,4*nflo
            call typen(l,nout,2)
            l=l+1
         enddo
      endif

      !--heating block
      l=lhnm
      do i=1,nes
         call typen(l,nout,2)
         l=l+1
      enddo

      !--finished
      call typen(l,nout,3)
      nern=0
      lrec=0

   !--write ace file in type 2 format
   else
      if (mcnpx.eq.0) then
         write(nout) hz(1:10),aw0,tz,hd,hk,hm,(izn(i),awn(i),i=1,16),&
           len2,z,nes,nflo,nxsd(1:12),&
           eszg,jinc,jcoh,jflo,lhnm,jxsd(1:27)
      else
         write(nout) hz(1:13),aw0,tz,hd,hk,hm,(izn(i),awn(i),i=1,16),&
           len2,z,nes,nflo,nxsd(1:12),&
           eszg,jinc,jcoh,jflo,lhnm,jxsd(1:27)
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
      nern=ner
      lrec=ner*nbw
   endif

   !--output directory for mcnp
   call openz(ndir,1)
   if (mcnpx.eq.0) then
      write(ndir,&
        '(a10,f12.6,'' filename route'',i2,'' 1 '',i8,2i6,1p,e10.3)')&
        hz(1:10),aw0,itype,len2,lrec,nern,tz
   else
      write(ndir,&
        '(a13,f12.6,'' filename route'',i2,'' 1 '',i8,2i6,1p,e10.3)')&
        hz(1:13),aw0,itype,len2,lrec,nern,tz
   endif
   call closz(ndir)
   return
   end subroutine phoout

end module acepa
