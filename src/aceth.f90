module aceth
   ! provides thermal ace stuff
   use locale
   implicit none
   private

   !--Public routines
   public acesix,thrfix

   !--Private global variables

   ! ace header parameters
   character(13)::hz
   character(10)::hd
   character(10)::hm
   real(kr)::aw0,tz

   ! ace nxs parameters for thermal data
   integer::len2,idpni,nil,nieb,idpnc,ncl,ifeng,nxsd(9)

   ! ace jxs parameters for thermal data
   integer::itie,itix,itxe,itce,itcx,itca,jxsd(26)

   ! body of the ace data
   integer,parameter::nxss=9000000
   real(kr)::xss(nxss)
   integer::nei
   real(kr),dimension(5000)::wt
   integer,dimension(5000)::ndp

contains

   subroutine acesix(nin,nace,ndir,matd,tempd,tname,suff,&
     hk,izn,awn,iza01,iza02,iza03,&
     mti,nbin,mte,ielas,nmix,emax,iwt,iprint,mcnpx)
   !-------------------------------------------------------------------
   ! Convert thermal matrices in njoy MF6 format to various ACE
   ! thermal formats.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util ! provides openz,repoz,mess,error,closz
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nace,ndir,matd,iza01,iza02,iza03
   integer::mti,nbin,mte,ielas,nmix,iwt,iprint,mcnpx
   real(kr)::tempd,suff,emax
   integer::izn(16)
   real(kr)::awn(16)
   character(70)::hk
   character(6)::tname
   ! internals
   integer::nw,nwscr,nb,indx,nwt,idis,law
   integer::ninmax,nea,j,idone,nee,i,neie,loc,iee,iea
   integer::ne,ie,nl,nep,is,nl1,l,isl,isn,k,itype,nang,nn
   real(kr)::fract,x
   real(kr)::temp,em,e,enext,clast,cnow,wtt,emin,sum
   real(kr)::grall,xl,yl,y,add,xn,yn,f,disc,sign,xbar,xlo,xhi
   real(kr)::area
   real(kr)::delta
   character(60)::strng
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::xs
   real(kr),dimension(:),allocatable::six
   real(kr),parameter::eps=1.e-5_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1
   real(kr),parameter::ttol=1

   integer::nout=10
   integer::nscr=11

   !--check for thermal data
   if (nace.eq.0.and.iprint.eq.0) then
      call closz(nscr)
      call closz(-nout)
      return
   endif

   !--process the thermal data
   nxsd=0
   jxsd=0
   ninmax=500000
   call openz(nin,0)
   nscr=iabs(nscr)
   if (nin.lt.0) nscr=-nscr
   call openz(nscr,1)
   call repoz(nscr)
   nout=iabs(nout)
   call openz(-nout,1)
   call repoz(-nout)
   call repoz(nin)

   !--allocate storage
   nw=npage+50
   allocate(xs(nw))
   allocate(six(ninmax))
   nwscr=500000
   allocate(scr(nwscr))

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

   !--find desired temperature on pendf tape
   call repoz(nin)
   call tpidio(nin,0,0,xs,nb,nw)
   temp=0
   delta=ttol
   do while (delta.ge.ttol)
      call contio(nin,0,0,xs,nb,nw)
      aw0=xs(2)
      if (iverf.ge.5) call contio(nin,0,0,xs,nb,nw)
      if (iverf.ge.6) call contio(nin,0,0,xs,nb,nw)
      call contio(nin,0,0,xs,nb,nw)
      temp=xs(1)
      delta=abs(temp-tempd)
      if (delta.ge.ttol) call tomend(nin,0,0,xs)
   enddo

   !--copy incoherent (and coherent elastic) cross section to nscr
   call findf(matd,3,mti,nin)
   call contio(nin,0,nscr,xs,nb,nw)
   call tosend(nin,0,nscr,xs)
   if (mte.ne.0.and.ielas.ne.0) then
      call findf(matd,3,mte,nin)
      call contio(nin,0,nscr,xs,nb,nw)
      call tosend(nin,0,nscr,xs)
   endif

   !--determine or check emax
   call repoz(nscr)
   call contio(nscr,0,0,scr,nb,nw)
   indx=1
   nwt=0
   call tab1io(nscr,0,0,scr(indx),nb,nw)
   indx=indx+nw
   nwt=nwt+nw
   do while (nb.ne.0)
      if (indx.gt.nwscr-npage) call error('aceth',&
        'exceeded scr storage',' ')
      call moreio(nscr,0,0,scr(indx),nb,nw)
      indx=indx+nw
      nwt=nwt+nw
   enddo
   em=scr(nwt-3)
   if (scr(nwt-2).lt.eps) em=scr(nwt-5)
   if (emax.gt.em) then
      write(strng,'(''emax reset to '',1p,e12.4)') em
      call mess('acesix',strng,' ')
   endif
   if (emax.gt.em) emax=em

   !--check for elastic data
   if (mte.eq.0) then
      write(nout) mte,mte
   else

      !--process coherent reaction if requested
      if (ielas.eq.0) then
         call findf(matd,3,mte,nin)
         call contio(nin,0,0,scr,nb,nw)
           e=0
            call gety1(e,enext,idis,x,nin,scr)
            nea=0

         !--locate bragg edges and compute cumulative probabilities
         j=0
         clast=0
         idone=0
         do while (idone.eq.0)
            e=enext
            call gety1(e,enext,idis,x,nin,scr)
            cnow=e*x
            if (enext.gt.emax.or.abs(enext-emax).le.eps) then
               idone=1
            else
               if ((cnow-clast).gt.(clast+1)/1000000) then
                  j=j+1
                  if (2*j.gt.ninmax) call error('acesix',&
                    'storage six exceeded for coherent reactions.',' ')
                  six(1+2*(j-1))=e
                  six(1+2*(j-1)+1)=cnow
                  clast=cnow
               endif
            endif
         enddo

         !--write out coherent reaction record
         nee=j
         nw=2*nee
         write(nout) nee,nea
         write(nout) (six(1-1+i),i=1,nw)

      !--incoherent elastic
      else
         call findf(matd,6,mte,nin)
         call contio(nin,0,0,scr,nb,nw)
         call tab1io(nin,0,0,scr,nb,nw)
         call tab2io(nin,0,0,scr,nb,nw)
         nee=nint(scr(6))
         neie=0
         call repoz(nscr)
         call findf(matd,3,mte,nscr)
         call contio(nscr,0,0,xs,nb,nw)
         e=0
         call gety1(e,enext,idis,x,nscr,xs)
         loc=1
         idone=0
         iee=0
         do while (iee.lt.nee.and.idone.eq.0)
            iee=iee+1
            indx=1
            call listio(nin,0,0,scr(indx),nb,nw)
            indx=indx+nw
            if (indx.gt.nwscr) call error('acesix',&
              'exceeded scr storage for incoherent elastic',' ')
            do while (nb.ne.0)
               if (indx.gt.nwscr-npage) call error('acesix',&
                 'exceeded scr storage for incoherent elastic',' ')
               call moreio(nin,0,0,scr(indx),nb,nw)
               indx=indx+nw
            enddo
            nea=nint(scr(6))-1
            e=scr(7)
            if (e.gt.emax.or.abs(e-emax).lt.eps) then
               idone=1
            else
               neie=neie+1
               call gety1(e,enext,idis,x,nscr,xs)
               if (loc.gt.ninmax-nea-1) call error('acesix',&
                 'exceeded six storage for incoherent elastic',' ')
               six(loc)=e
               six(loc+1)=x
               do iea=2,nea
                     six(iea+loc)=scr(iea+7)
               enddo
               loc=loc+1+nea
            endif
         enddo
         loc=loc-1
         nea=nea-1
         write(nout) neie,nea
         write(nout) (six(i),i=1,loc)
      endif
   endif

   !--find inelastic energy-angle distribution
   call findf(matd,6,mti,nin)
   call contio(nin,0,0,scr,nb,nw)
   call tab1io(nin,0,0,scr,nb,nw)
   law=l2h
   if (law.ne.1) call error('acesix',&
     'coded for E-E-mu ordering only.',' ')
   call tab2io(nin,0,0,scr,nb,nw)
   ne=nint(scr(6))
   call repoz(nscr)
   nei=0
   call contio(nscr,0,0,xs,nb,nw)
   e=0
   call gety1(e,enext,idis,x,nscr,xs)

   !--construct pattern of weights (constant or variable)
   !--if needed
   if (iwt.eq.0) then
      wtt=1
      wtt=wtt/10
      wtt=wtt/(nbin-3)
      do i=1,nbin
         wt(i)=10*wtt
         if (i.eq.1.or.i.eq.nbin) wt(i)=wtt
         if (i.eq.2.or.i.eq.nbin-1) wt(i)=4*wtt
      enddo
      write(nsyso,'(/&
        &'' relative weights for energy bins are  '',&
        &''1  4  10...10  4  1'')')
   else if (iwt.eq.1) then
      wtt=1
      wtt=wtt/nbin
      do i=1,nbin
         wt(i)=wtt
      enddo
      write(nsyso,'(/'' constant weighting for energy bins'')')
   else
      write(nsyso,'(/'' tabulated probability distribution''//&
        &'' original and modified number of secondary energies'')')
   endif

   !--loop over incident energies
   len2=0
   do 210 ie=1,ne
   call listio(nin,0,0,scr,nb,nw)
   e=scr(2)
   emin=1
   emin=emin/1000000
   emin=emin+emin/1000000
   if (ie.eq.1.and.e.lt.emin) e=emin
   if (e.gt.emax) go to 290
   if (abs(e-emax).lt.eps) go to 290
   nei=nei+1
   nl=nint(scr(6))
   nang=nl-2
   nep=nint(scr(5)/nl)
   indx=1+nw
   do while (nb.ne.0)
      if (indx.gt.nwscr-npage) call error('acesix',&
        &'exceeded scr storage for incoherent inelastic.',' ')
      call moreio(nin,0,0,scr(indx),nb,nw)
      indx=indx+nw
   enddo

   !--fix the angular distribution of the first and last points
   nn=2+nang
   isn=8
   do i=1,nang
      scr(isn+i)=scr(isn+i+nn)
   enddo
   isn=8+nl*(nep-1)
   do i=1,nang
      scr(isn+i)=scr(isn+i-nn)
   enddo

   !--get cross section
   call gety1(e,enext,idis,x,nscr,xs)
   six(1)=e
   six(2)=x
   loc=3
   is=1
   nl1=nl
   nn=0

   !--determine equal or variable probable energies
   if (iwt.le.1) then
      fract=wt(1)
   else
      fract=1
   endif
   sum=0
   grall=0
   xl=scr(is+6)
   yl=scr(is+7)
   i=0
   j=0
   l=2
  240 continue
   if (i.ge.nep) go to 285
   i=i+1
   x=scr(is+6+nl1*(i-1))
   y=scr(is+7+nl1*(i-1))
  250 continue
   add=(y+yl)*(x-xl)/2
   if (x.eq.xl) go to 280
   if (iwt.le.1.and.i.eq.nep.and.j.eq.nbin-1) then
      xn=x
      j=j+1
      go to 270
   endif
   if (iwt.gt.1.and.i.eq.nep) then
      xn=x
      fract=sum+add
      j=j+1
      go to 270
   endif
   if (iwt.gt.1) then
     if (sum+add.gt.eps/10) fract=sum+add
     !if (fract.gt.one/nbin) fract=one/nbin
   endif
   if (sum+add.ge.fract-fract/10000) go to 260
   sum=sum+add
   grall=grall+(yl-(y-yl)*xl/(x-xl))*(x*x-xl*xl)/2+&
     ((y-yl)/(x-xl))*(x**3-xl**3)/3
   go to 280
  260 continue
   j=j+1
   if (sum+add.lt.fract+fract/10000) then
      xn=x
   else if (abs(y-yl).gt.(y+yl)/100000) then
      f=(y-yl)/(x-xl)
      disc=(yl/f)**2+2*(fract-sum)/f
      if (disc.lt.zero) then
         write(strng,'(''disc='',1p,e12.4)') disc
         call mess('acesix',strng,'set to abs value and continue')
         disc=-disc
      endif
      sign=1
      if (f.lt.zero) sign=-1
      xn=xl-(yl/f)+sign*sqrt(disc)
      if (xn.le.xl.or.xn.gt.x) then
         call error('acesix','solution out of range1',' ')
      endif
   else
      xn=xl+(fract-sum)/yl
      if (xn.gt.x) then
         xn=x
      else
         if (xn.le.xl.or.xn.gt.x) then
            call error('acesix','solution out of range2',' ')
         endif
      endif
   endif
 270 continue
   yn=yl+(y-yl)*(xn-xl)/(x-xl)
   sum=sum+(yn+yl)*(xn-xl)/2
   grall=grall+(yl-(y-yl)*xl/(x-xl))*(xn*xn-xl*xl)/2+&
     ((y-yl)/(x-xl))*(xn**3-xl**3)/3
   xbar=grall/sum
   l=2
   do
      xlo=scr(7+nl1*(l-2))
      xhi=scr(7+nl1*(l-1))
      if (l.eq.nep) exit
      if (iwt.le.1) then
         if (xhi.ge.xbar) exit
      endif
      if (iwt.eq.2) then
         if (xhi.ge.xn) exit
      endif
      l=l+1
   enddo
   if (loc.gt.ninmax-nang) call error('acesix','exceeded storage in six',' ')
   nn=nn+1
   isl=is+7+nl1*(l-2)
   isn=is+7+nl1*(l-1)
   if (iwt.le.1) then
      six(loc)=xbar
   else
      six(loc)=xn
      loc=loc+1
      six(loc)=yn
      loc=loc+1
      six(loc)=fract
   endif
   do k=1,nang
      if (iwt.le.1) then
         six(k+loc)=scr(k+isl)&
           +(scr(k+isn)-scr(k+isl))*(xbar-xlo)/(xhi-xlo)
      else
         if (abs(xn-xhi).lt.eps) then
            six(k+loc)=scr(k+isn)
         else if (abs(xn-xlo).lt.eps) then
            six(k+loc)=scr(k+isl)
         else
            six(k+loc)=scr(k+isl)+&
             (scr(k+isn)-scr(k+isl))*(xn-xlo)/(xhi-xlo)
         endif
      endif
      if ((six(k+loc).lt.-1.0e0_kr).or.(six(k+loc).gt.1.0e0_kr)) then
         write(strng,'(''cosine '',f12.8,&
                     &'' outside [-1,1] range for e_in -> e_out '',&
                     &1p,e13.6,1x,e13.6,&
                     &'',  bin_mu= '',i4)') six(k+loc),e,xbar,k
         if (six(k+loc).lt.-1.0e0_kr) then
            six(k+loc) = -1.0e0_kr
            call mess('acesix',strng,'cosine set to -1.0')
         endif
         if (six(k+loc).gt.1.0e0_kr) then
            six(k+loc) = 1.0e0_kr
            call mess('acesix',strng,'cosine set to 1.0')
         endif
      endif
   enddo
   loc=loc+1+nang
   xl=xn
   yl=yn
   if (iwt.le.1) then
      fract=wt(j+1)
   else
      fract=1
   endif
   sum=0
   grall=0
   if (iwt.le.1.and.j.eq.nbin) go to 285
   if (xl.lt.x) go to 250
  280 continue
   xl=x
   yl=y
   go to 240
  285 continue

   !--check normalization for tabulated distribution
   if (iwt.gt.1) then
      area=0
      j=3
      do i=1,nn
         area=area+six(j+2)
         six(j+2)=area
         j=j+3+nang
      enddo
      j=3
      do i=1,nn
         six(j+1)=six(j+1)/area
         six(j+2)=six(j+2)/area
         j=j+3+nang
      enddo
      write(nsyso,'(6x,3i8)') ie,nep,nn
   endif

   !--write incoherent reaction record for this energy
   loc=loc-1
   len2=len2+loc-2
   ndp(ie)=loc
   write(nout) (six(i),i=1,loc)
  210 continue
  290 continue

   !--write ace format thermal tape
   call thrlod(nout,matd,tempd,tname,suff,izn,iza01,iza02,iza03,&
     nbin,nang,nmix,iwt,mcnpx)
   if (iprint.gt.0) call thrprt(hk)
   itype=1
   call throut(itype,nace,ndir,mcnpx,hk,izn,awn)

   !--finished with thermal data
   call closz(nscr)
   call closz(-nout)
   return
   end subroutine acesix

   subroutine thrfix(itype,nin,nout,ndir,iprint,nplot,mcnpx,suff,&
     nxtra,hk,izn,awn)
   !-------------------------------------------------------------------
   ! Print and or edit an ACE thermal file.
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
           '(a10,e12.0,e12.0,1x,a10/a70,a10)')&
           hz(1:10),aw0,tz,hd,hko,hm
      else
         read(nin,&
           '(a13,f12.0,e12.0,1x,a10/a70,a10)')&
           hz(1:13),aw0,tz,hd,hko,hm
      endif
      read (nin,'(4(i7,f11.0))') (izo(i),awo(i),i=1,16)
      read (nin,'(8i9)')&
        len2,idpni,nil,nieb,idpnc,ncl,ifeng,nxsd,&
        itie,itix,itxe,itce,itcx,itca,jxsd
      n=(len2+3)/4
      n=n-1
      l=0
      do i=1,n
         read (nin,'(4e20.0)') (xss(l+j),j=1,4)
         l=l+4
      enddo

   !--read type 2 ace format file
   else if (itype.eq.2) then
      if (mcnpx.eq.0) then
        read(nin) hz(1:10),aw0,tz,hd,hko,hm,(izo(i),awo(i),i=1,16),&
          len2,idpni,nil,nieb,idpnc,ncl,ifeng,nxsd,&
          itie,itix,itxe,itce,itcx,itca,jxsd
      else
        read(nin) hz(1:13),aw0,tz,hd,hko,hm,(izo(i),awo(i),i=1,16),&
          len2,idpni,nil,nieb,idpnc,ncl,ifeng,nxsd,&
          itie,itix,itxe,itce,itcx,itca,jxsd
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
   if (nxtra.eq.0) then
      do i=1,16
         izn(i)=izo(i)
         awn(i)=awo(i)
      enddo
   endif

   !--print, plot, and/or write the file.
   if (iprint.gt.0) call thrprt(hk)
   if (nout.gt.0) call throut(itype,nout,ndir,mcnpx,hk,izn,awn)
   if (nplot.ne.0) call tplots(nplot,hk)

   return
   end subroutine thrfix

   subroutine thrlod(nin,matd,tempd,tname,suff,izn,iza01,iza02,iza03,&
     nbini,nang,nmix,iwt,mcnpx)
   !-------------------------------------------------------------------
   ! write a thermal output tape in ace format.
   !-------------------------------------------------------------------
   use physics ! provides bk
   use util ! provides date,repoz
   use mainio ! provides nsyso
   ! externals
   integer::nin,matd,iza01,iza02,iza03,nbini,nang,nmix,iwt,mcnpx
   real(kr)::tempd,suff
   integer::izn(16)
   character(6)::tname
   ! internals
   integer::nwscr,nw,i,loc,indx,j,k,jang
   integer::nee,nea,nce,nie
   character(8)::hdt
   real(kr),dimension(:),allocatable::scr
   real(kr),parameter::emev=1.e6_kr

   !--initialize
   nwscr=500000
   allocate(scr(nwscr))
   write(hm,'(''   mat'',i4)') matd
   tz=tempd*bk/emev
   if (mcnpx.eq.0) then
      write(hz,'(a6,f3.2,''t'')') tname,suff
   else
      write(hz,'(a6,f4.3,''nt '')') tname,suff
   endif
   call dater(hdt)
   hd='  '//hdt
   izn(1)=iza01
   izn(2)=iza02
   izn(3)=iza03
   ncl=0
   idpni=3
   ifeng=0
   if (iwt.eq.0) ifeng=1
   if (iwt.eq.2) ifeng=2
   call repoz(-nin)

   !--read elastic pointers
   read(nin) nee,nea
   nce=nee
   nie=nei

   !--assign pointers in xss array
   itie=1
   itix=itie+1+nie
   itxe=itix+nie
   if (ifeng.gt.1) len2=len2+2*nie
   len2=len2+itxe-1
   if (nee.gt.0) then
      itce=len2+1
      itcx=itce+nee+1
      len2=itcx+nee-1
   endif
   itca=0
   if (nee.gt.0.and.nea.gt.0) then
      itca=itcx+nee
      len2=itca+nee*nea-1
   endif
   if (len2.gt.nxss) then
      write(nsyso,'(i10)') len2
      call error('thrlod','xss too small',' ')
   endif

   !--check for elastic data
   if (nee.ne.0) then
      xss(itce)=nce

      !--process coherent elastic data
      if (nea.le.0) then
         ncl=-1
         idpnc=4
         nw=2*nee
         if (nw.gt.nwscr) call error('thrlod','scr exceeded',' ')
         read(nin) (scr(i),i=1,nw)
         do i=1,nee
            xss(i+itce)=scr(1+2*(i-1))/emev
            xss(i+nee+itce)=scr(2+2*(i-1))/emev/nmix
         enddo

      !--process incoherent elastic data
      else
         ncl=nea-1
         idpnc=3
         nw=nee*(2+nea)
         if (nw.gt.nwscr) call error('thrlod','scr exceeded',' ')
         read(nin) (scr(i),i=1,nw)
         loc=itca-1
         indx=1
         do i=1,nee
            xss(itce+i)=scr(indx)/emev
            xss(itce+nee+i)=scr(indx+1)/nmix
            do j=1,nea
               xss(j+loc)=scr(indx+1+j)
            enddo
            loc=loc+nea
            indx=indx+nea+2
         enddo
      endif
   endif

   !--process inelastic data
   nie=nei
   nieb=nbini
   if (ifeng.le.1) then
      nil=nang-1
   else
      nil=nang+1
   endif
   xss(itie)=nei
   indx=itxe-1
   if (ifeng.gt.1) indx=indx+2*nei
   do i=1,nei
      nw=ndp(i)
      if (nw.gt.nwscr) call error('thrlod','scr exceeded',' ')
      read(nin) (scr(j),j=1,nw)
      nw=nw-2
      xss(itie+i)=scr(1)/emev
      xss(itie+nie+i)=scr(2)/nmix
      if (ifeng.gt.1) then
         xss(itxe-1+i)=indx
         xss(itxe-1+nei+i)=nw/(nang+3)
      endif
      k=0
      do while (k.lt.nw-1)
         k=k+1
         xss(indx+k)=scr(2+k)/emev
         if (ifeng.gt.1) then
            k=k+1
            xss(indx+k)=scr(2+k)*emev
            k=k+1
            xss(indx+k)=scr(2+k)
         endif
         do jang=1,nang
            k=k+1
            xss(indx+k)=scr(2+k)
         enddo
      enddo
      indx=indx+k
   enddo
   return
   end subroutine thrlod

   subroutine thrprt(hk)
   !-------------------------------------------------------------------
   ! Print thermal ACE data.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   ! externals
   character(70)::hk
   ! internals
   integer::nie,nee,loc,nang,nbini,nln,lim,lim1,i,lines
   integer::j,k,nea,ncol,npg,ipg,indx,mcol,ifini,ne
   integer::ind(4)
   real(kr)::b(4),c(4)
   character(10),parameter::labl1='bragg edge'
   character(7),parameter::labl2='e*cross'
   character(1),parameter::labl3='i'
   character(6),parameter::labl4='energy'
   character(7),parameter::labl5='section'
   character(10),parameter::dash4='    '
   character(10),parameter::dash10='          '

   !--print thermal header information
   nie=nint(xss(itie))
   nee=0
   if (itce.gt.0) nee=nint(xss(itce))
   write(nsyso,'(''1''///////&
     &38x,''zaid'',1x,a13/39x,''awr'',f10.3/&
     &38x,''temp'',1p,e10.2/38x,''date'',a10/39x,''mat'',a10/&
     &6x,''***********************''/&
     &6x,''*                     *'',9x,''len2'',i10/&
     &6x,''*       thermal       *'',8x,''idpni'',i10/&
     &6x,''*                     *'',10x,''nil'',i10/&
     &6x,''*   ace format file   *'',9x,''nieb'',i10/&
     &6x,''*                     *'',8x,''idpnc'',i10/&
     &6x,''*     processed by    *'',10x,''ncl'',i10/&
     &6x,''*                     *'',8x,''ifeng'',i10/&
     &6x,''*        njoy         *''/&
     &6x,''*                     *'',9x,''itie'',i10/&
     &6x,''***********************'',9x,''itix'',i10/&
     &38x,''itxe'',i10/38x,''itce'',i10/38x,''itcx'',i10/&
     &38x,''itca'',i10//39x,''nie'',i10/39x,''nee'',i10///&
     &6x,''hk---'',a70)')&
     hz,aw0,tz,hd,hm,len2,idpni,nil,nieb,idpnc,ncl,&
     ifeng,itie,itix,itxe,itce,itcx,itca,nie,nee,hk

   !--print inelastic data
   loc=itxe-1
   if (ifeng.gt.1) loc=loc+2*nie
   ne=nie
   do i=1,ne
      if (ifeng.le.1) then
         nang=nil+1
         nbini=nieb
      else
         nang=nil-1
         nbini=nint(xss(itxe-1+ne+i))
      endif
      nln=((nang+7)/8)*nbini
      if (ifeng.le.1) then
         lim=nang+1
         if (nang.gt.8) lim=9
      else
         lim=nang+3
        if (nang.gt.8) lim=11
      endif
      lim1=lim+1
      if (i.eq.1) then
         write(nsyso,'(/&
           &'' inelastic data - equally probable angles''/&
           &'' ----------------------------------------''/)')
         lines=4
      else
         if ((lines+nln+4).gt.58) then
            write(nsyso,'(''1'')')
            lines=1
         endif
      endif
       write(nsyso,'(/6x,''incident energy = '',1p,e13.6,8x,&
         &''cross section = '',e13.6)') xss(itie+i),xss(itie+ne+i)
       write(nsyso,'(5x,'' num of e-prime '',i5)') nbini
       if (ifeng.le.1) then
          write(nsyso,'(/&
            &9x,''exit energy'',5x,''cosines''/9x,''-----------'',&
            &4x,8(''------------''))')
       else
          write(nsyso,'(/&
            &9x,''exit energy'',8x,''pdf'',11x,''cdf'',10x,''cosines''/&
            &9x,''-----------'',4x,''-----------'',&
            &3x,''------------'',4x,8(''------------''))')
       endif
       lines=lines+4
       do j=1,nbini
         if (ifeng.le.1) then
            write(nsyso,'(7x,1p,e13.6,1x,0p,8f12.6)')&
              (xss(k+loc),k=1,lim)
            if (nang.gt.8) write(nsyso,'(21x,8f12.6)')&
              (xss(loc+k),k=lim1,nang+1)
         else
            write(nsyso,'(7x,1p,e13.6,e15.6,e15.6,1x,0p,8f12.6)')&
              (xss(k+loc),k=1,lim)
            if (nang.gt.8) write(nsyso,'(51x,8f12.6)')&
              (xss(loc+k),k=lim1,nang+3)
         endif
         if (ifeng.le.1) then
            loc=loc+nang+1
         else
            loc=loc+nang+3
         endif
      enddo
      lines=lines+nln
   enddo

   !--check for elastic data
   if (nee.eq.0) return

   !--print incoherent elastic data
   if (idpnc.ne.4) then
      nea=ncl+1
      write(nsyso,'(''1''/&
        &'' incoherent elastic data - equally probable angles''/&
        &'' -------------------------------------------------''/)')
      write(nsyso,'(/&
        &9x,''incident'',9x,''cross''/&
        &4x,''i'',5x,''energy'',9x,''section'',25x,''angles'')')
      lines=7
      nln=(nea+8)/9
      lim=nea
      if (nea.gt.8) lim=8
      lim1=lim+1
      loc=itca-1
      do i=1,nee
         if ((lines+nln).gt.58) then
            write(nsyso,'(''1'')')
            write(nsyso,'(/&
              &9x,''incident'',9x,''cross''/&
              &4x,''i'',5x,''energy'',9x,''section'',25x,''angles'')')
            lines=4
         endif
         write(nsyso,'(/2x,i3,1x,1p,e12.4,3x,e12.4,3x,0p,8f10.4)')&
           i,xss(itce+i),xss(itce+nee+i),(xss(loc+j),j=1,lim)
         if (nea.gt.8) write(nsyso,'(36x,8f10.4)')&
           (xss(loc+j),j=lim1,nea)
         loc=loc+nea
         lines=lines+nln
      enddo

   !--print coherent elastic data
   else
      write(nsyso,'(''1''/&
        &'' coherent elastic data - bragg edges and cumulative '',&
        &''intensity''/&
        &'' ---------------------------------------------------'',&
        &''---------'')')
      ncol=(nee+49)/50
      npg=(ncol+3)/4
      ipg=0
      indx=itce+1
      do while (ipg.lt.npg)
         ipg=ipg+1
         if (ipg.gt.1) write(nsyso,'(''1'')')
         mcol=ncol
         if (mcol.gt.4) mcol=4
         ncol=ncol-mcol
         write(nsyso,'(/4(9x,a10,5x,a7,1x))') (labl1,labl2,j=1,mcol)
         write(nsyso,'(5x,3(a1,5x,a6,7x,a7,6x),a1,5x,a6,7x,a7)')&
           (labl3,labl4,labl5,j=1,mcol)
            write(nsyso,'(4(2x,a4,3x,a10,3x,a10))')&
           (dash4,dash10,dash10,j=1,mcol)
         ifini=0
         i=0
         do while (i.lt.50.and.ifini.eq.0)
            i=i+1
            ind(1)=indx-itce
            b(1)=xss(indx)
            c(1)=xss(indx+nee)
            if (mcol.gt.1) then
               ind(2)=indx+50-itce
               b(2)=xss(indx+50)
               c(2)=xss(indx+nee+50)
               if (mcol.gt.2) then
                  ind(3)=indx+100-itce
                  b(3)=xss(indx+100)
                  c(3)=xss(indx+nee+100)
                  if (mcol.gt.3) then
                     ind(4)=indx+150-itce
                     b(4)=xss(indx+150)
                     c(4)=xss(indx+nee+150)
                  endif
               endif
            endif
            write(nsyso,&
              '(i6,1p,2e13.4,i6,2e13.4,i6,2e13.4,i6,2e13.4)')&
              (ind(j),b(j),c(j),j=1,mcol)
            indx=indx+1
            if (ipg.ge.npg) then
               if (indx-itce.gt.nee) ifini=1
               if (ifini.ne.1) then
                  if ((indx-itce+(mcol-1)*50).gt.nee) mcol=mcol-1
                  if (mcol.eq.0) ifini=1
               endif
            endif
         enddo
         if (ifini.eq.0) indx=itce+1+200*ipg
      enddo
   endif
   return
   end subroutine thrprt

   subroutine tplots(nout,hk)
   !-------------------------------------------------------------------
   ! Do standard plots of a thermal ACE file.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util   ! provides openz
   ! externals
   integer::nout
   character(70)::hk
   ! internals
   integer::ipcol,iwcol,nie,nee,i,it,idone,nang,nbini
   integer::loc,k,major,minor,j,ie,iskip
   real(kr)::xmin,xmax,ymin,ymax
   real(kr)::e,xnelas,xelas,e1,x1,e2,x2,tot
   real(kr)::xtag,ytag,xs,xn,ubar,ystep
   real(kr)::ell,bl,ei,ui,bi,ebar,eprime,sum,wt
   real(kr)::zmin,zmax,ep,epl,x,xl,p,pl,cdl,u,ul,un,skip
   character(1)::qu=''''
   real(kr),parameter::big=1.e10_kr
   real(kr),parameter::ten=10.e0_kr
   real(kr),parameter::scale=1.e6_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1

   !--black and white pages
   ipcol=0
   iwcol=0
   !--colored pages
   ipcol=2
   iwcol=3

   !--start the viewr input text
   call openz(nout,1)
   write(nout,'(''1 2 .30'',i3,''/'')') ipcol
   nie=nint(xss(itie))
   nee=0
   if (itce.gt.0) nee=nint(xss(itce))

   !--plot log-log total, inelastic, and elastic (if present)
   xmin=big
   xmax=0
   ymin=big
   ymax=-big
   do i=1,nie
      e=xss(itie+i)
      xnelas=xss(itie+nie+i)
      if (e.lt.xmin) xmin=e
      if (e.gt.xmax) xmax=e
      if (xnelas.lt.ymin) ymin=xnelas
      if (xnelas.gt.ymax) ymax=xnelas
      if (nee.ne.0) then
         j=0
         xelas=0
         do while (j.lt.nee)
            j=j+1
            e1=xss(itce+j)
            x1=xss(itce+nee+j)
            if (idpnc.eq.4) x1=x1/e1
            e2=xss(itce+j+1)
            x2=xss(itce+nee+j+1)
            if (idpnc.eq.4) x2=x2/e2
            if (e.ge.e1.and.e.le.e2) xelas=x1+(e-e1)*(x2-x1)/(e2-e1)
         enddo
         if (xelas.gt.zero) then
            tot=xnelas+xelas
            if (xelas.lt.ymin) ymin=xelas
            if (xelas.gt.ymax) ymax=xelas
            if (tot.lt.ymin) ymin=tot
            if (tot.gt.ymax) ymax=tot
         endif
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
   write(nout,'(a,''<t>hermal cross sections'',a,''/'')') qu,qu
   xtag=2*xmin
   ytag=7*log10(ymin)/10+3*log10(ymax)/10
   ytag=ten**ytag
   write(nout,'(''4 0 2 1'',2e12.4,''/'')') xtag,ytag
   write(nout,'(1p,3e12.3,''/'')') xmin,xmax,one
   write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
   write(nout,'(1p,3e12.3,''/'')') ymin,ymax,one
   write(nout,'(a,''<c>ross section (barns)'',a,''/'')') qu,qu
   write(nout,'(''/'')')
   if (iwcol.eq.0) then
      write(nout,'(''0 0 1/'')')
   else
      write(nout,'(''0 0 0 1/'')')
   endif
   write(nout,'(a,''inelastic'',a,''/'')') qu,qu
   write(nout,'(''0/'')')
   do i=1,nie
      e=xss(itie+i)
      xs=xss(itie+nie+i)
      write(nout,'(1p,2e14.6,''/'')') e,xs
   enddo
   write(nout,'(''/'')')
   if (nee.ne.0) then
      write(nout,'(''2/'')')
      write(nout,'(''/'')')
      if (iwcol.eq.0) then
         write(nout,'(''0 0 1/'')')
      else
         write(nout,'(''0 0 0 2/'')')
      endif
      write(nout,'(a,''elastic'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      if (idpnc.eq.4) then
         e=xss(itce+1)-xss(itce+1)/1000
         xs=xss(itce+nee+1)/e
         xs=xs/100
         write(nout,'(1p,2e14.6,''/'')') e,xs
      endif
      do i=1,nee
         e=xss(itce+i)
         xs=xss(itce+nee+i)
         if (idpnc.eq.4) xs=xs/e
         write(nout,'(1p,2e14.6,''/'')') e,xs
         if (idpnc.eq.4.and.i.lt.nee) then
            e=xss(itce+i+1)-xss(itce+i+1)/1000
            xs=xss(itce+nee+i)/e
            write(nout,'(1p,2e14.6,''/'')') e,xs
         endif
      enddo
      if (idpnc.eq.4) then
         do i=1,25
            e=e+e/10
            if (e.gt.xmax) exit
            xs=xss(itce+2*nee)/e
            write(nout,'(1p,2e14.6,''/'')') e,xs
         enddo
      endif
      write(nout,'(''/'')')
      write(nout,'(''3/'')')
      write(nout,'(''/'')')
      write(nout,'(''/'')')
      write(nout,'(a,''total'',a,''/'')') qu,qu
      write(nout,'(''0/'')')
      i=0
      idone=0
      e=0
      do while (idone.eq.0.and.i.lt.nee)
         i=i+1
         e=xss(itie+i)
         if (e.ge.xss(itce+1)) then
            idone=1
         else
            xs=xss(itie+nie+i)
            write(nout,'(1p,2e14.6,''/'')') e,xs
         endif
      enddo
      if (idpnc.eq.4) then
         e=xss(itce+1)-xss(itce+1)/1000
         j=0
         do while (j.lt.nie)
            j=j+1
            e1=xss(itie+j)
            x1=xss(itie+nie+j)
            e2=xss(itie+j+1)
            x2=xss(itie+nie+j+1)
            if (e.ge.e1.and.e.le.e2) xn=x1+(e-e1)*(x2-x1)/(e2-e1)
         enddo
         tot=xn
         write(nout,'(1p,2e14.6,''/'')') e,tot
      endif
      do i=1,nee
         e=xss(itce+i)
         xs=xss(itce+nee+i)
         if (idpnc.eq.4) xs=xs/e
         j=0
         do while (j.lt.nie)
            j=j+1
            e1=xss(itie+j)
            x1=xss(itie+nie+j)
            e2=xss(itie+j+1)
            x2=xss(itie+nie+j+1)
            if (e.ge.e1.and.e.le.e2) xn=x1+(e-e1)*(x2-x1)/(e2-e1)
         enddo
         tot=xs+xn
         write(nout,'(1p,2e14.6,''/'')') e,tot
         if (i.lt.nee.and.idpnc.eq.4) then
            e=xss(itce+i+1)-xss(itce+i+1)/1000
            xs=xss(itce+nee+i)/e
            j=0
            do while (j.lt.nie)
               j=j+1
               e1=xss(itie+j)
               x1=xss(itie+nie+j)
               e2=xss(itie+j+1)
               x2=xss(itie+nie+j+1)
               if (e.ge.e1.and.e.le.e2) xn=x1+(e-e1)*(x2-x1)/(e2-e1)
            enddo
            tot=xs+xn
            write(nout,'(1p,2e14.6,''/'')') e,tot
         endif
      enddo
      if (idpnc.eq.4) then
         do i=1,25
            e=e+e/10
            if (e.gt.xmax) exit
            xs=xss(itce+2*nee)/e
            j=0
            do while (j.lt.nie)
               j=j+1
               e1=xss(itie+j)
               x1=xss(itie+nie+j)
               e2=xss(itie+j+1)
               x2=xss(itie+nie+j+1)
               if (e.ge.e1.and.e.le.e2) xn=x1+(e-e1)*(x2-x1)/(e2-e1)
            enddo
            tot=xs+xn
            write(nout,'(1p,2e14.6,''/'')') e,tot
         enddo
      endif
      write(nout,'(''/'')')
   endif

   !--mubar plot
   nang=nil+1
   if (ifeng.gt.1) nang=nil-1
   nbini=nieb
   xmin=big
   xmax=0
   ymin=big
   ymax=-big
   loc=itxe-1
   if (ifeng.gt.1) loc=loc+2*nie
   do i=1,nie
      e=xss(itie+i)
      xnelas=xss(itie+nie+i)
      if (ifeng.gt.1) nbini=nint(xss(itxe-1+nie+i))
      ubar=0
      sum=0
      if (ifeng.le.1) then
         do j=1,nbini
            if (ifeng.eq.0) then
               wt=1
            else
               wt=10
               if (j.eq.1.or.j.eq.nbini) wt=1
               if (j.eq.2.or.j.eq.nbini-1) wt=4
            endif
            do k=2,nang+1
              ubar=ubar+wt*xss(loc+k)
              sum=sum+wt
            enddo
            loc=loc+nang+1
         enddo
      else
         cdl=0
         loc=loc+nang+3
         do j=2,nbini
            p=xss(loc+3)-cdl
            do k=1,nang
!              if ((xss(loc+3+k).lt.-1.0e0_kr).or.&
!                  (xss(loc+3+k).gt.1.0e0_kr)) then
!                 write(nsyso,'(/'' ---warning from tplots---'', &
!                       &'' cosine '',f12.8,'' outside [-1,1] range'')')&
!                       xss(loc+3+k)
!              endif
              ubar=ubar+xss(loc+3+k)*p/2
              ubar=ubar+xss(loc+3+k-nang-2)*p/2
              sum=sum+p
            enddo
            cdl=xss(loc+3)
            loc=loc+nang+3
         enddo
      endif
      ubar=ubar/sum
      if (e.lt.xmin) xmin=e
      if (e.gt.xmax) xmax=e
      if (ubar.lt.ymin) ymin=ubar
      if (ubar.gt.ymax) ymax=ubar
   enddo
   if (nee.gt.0) then
      ymin=-1
      ymax=1
   endif
   call ascll(xmin,xmax)
   call ascle(4,ymin,ymax,major,minor)
   ystep=(ymax-ymin)/major
   write(nout,'(''1'',i3,''/'')') iwcol
   it=1
   do i=1,70
      if (hk(i:i).ne.' ') it=i
   enddo
   write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
   write(nout,'(a,''<t>hermal mubar'',a,''/'')') qu,qu
   xtag=2*xmin
   ytag=7*ymin/10+3*ymax/10
   write(nout,'(''3 0 2 1'',2e12.4,''/'')') xtag,ytag
   write(nout,'(1p,3e12.3,''/'')') xmin,xmax,one
   write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
   write(nout,'(1p,3e12.3,''/'')') ymin,ymax,ystep
   write(nout,'(a,''<m>ubar'',a,''/'')') qu,qu
   write(nout,'(''/'')')
   if (iwcol.eq.0) then
      write(nout,'(''0 0 1/'')')
   else
      write(nout,'(''0 0 0 1/'')')
   endif
   write(nout,'(a,''inelastic'',a,''/'')') qu,qu
   write(nout,'(''0/'')')
   loc=itxe-1
   if (ifeng.gt.1) loc=loc+2*nie
   do i=1,nie
      e=xss(itie+i)
      xs=xss(itie+nie+i)
      if (ifeng.gt.1) nbini=nint(xss(itxe-1+nie+i))
      ubar=0
      sum=0
      if (ifeng.le.1) then
         do j=1,nbini
            if (ifeng.eq.0) then
               wt=1
            else
               wt=10
               if (j.eq.1.or.j.eq.nbini) wt=1
               if (j.eq.2.or.j.eq.nbini-1) wt=4
            endif
            do k=2,nang+1
              ubar=ubar+wt*xss(loc+k)
              sum=sum+wt
            enddo
            loc=loc+nang+1
         enddo
      else
         cdl=0
         loc=loc+nang+3
         do j=2,nbini
            p=xss(loc+3)-cdl
            do k=1,nang
              ubar=ubar+xss(loc+3+k)*p/2
              ubar=ubar+xss(loc+3+k-nang-2)*p/2
              sum=sum+p
            enddo
            cdl=xss(loc+3)
            loc=loc+nang+3
         enddo
      endif
      ubar=ubar/sum
      write(nout,'(1p,2e14.6,''/'')') e,ubar
   enddo
   write(nout,'(''/'')')
   if (nee.ne.0) then
      write(nout,'(''2/'')')
      write(nout,'(''/'')')
      if (iwcol.eq.0) then
         write(nout,'(''0 0 1/'')')
      else
         write(nout,'(''0 0 0 2/'')')
      endif
      if (idpnc.eq.4) then
         write(nout,'(a,''coherent elastic'',a,''/'')') qu,qu
         write(nout,'(''0/'')')
         e=xss(itce+1)
         ell=xss(itce+nee)
         do while (e.le.ell)
            idone=0
            i=0
            ubar=0
            sum=0
            bl=0
            do while (idone.eq.0.and.i.lt.nee)
               i=i+1
               ei=xss(itce+i)
               if (ei.gt.e) then
                  idone=1
               else
                  ui=1-2*ei/e
                  bi=(xss(itce+nee+i)-bl)
                  ubar=ubar+bi*ui/e
                  sum=sum+bi/e
               endif
               bl=xss(itce+nee+i)
            enddo
            ubar=ubar/sum
            write(nout,'(1p,2e14.6,''/'')') e,ubar
            e=e+e/50
         enddo
      else
         write(nout,'(a,''incoherent elastic'',a,''/'')') qu,qu
         write(nout,'(''0/'')')
         loc=itca-1
         do i=1,nee
            e=xss(itce+i)
            ubar=0
            do j=1,ncl+1
               ubar=ubar+xss(loc+j)/(ncl+1)
            enddo
            write(nout,'(1p,2e14.6,''/'')') e,ubar
            loc=loc+ncl+1
         enddo
      endif
      write(nout,'(''/'')')
   endif

   !--ebar plot
   nang=nil+1
   if (ifeng.gt.1) nang=nil-1
   nbini=nieb
   xmin=big
   xmax=0
   ymin=big
   ymax=-big
   loc=itxe-1
   if (ifeng.gt.1) loc=loc+2*nie
   do i=1,nie
      e=xss(itie+i)
      xnelas=xss(itie+nie+i)
      if (ifeng.gt.1) nbini=nint(xss(itxe-1+nie+i))
      ebar=0
      sum=0
      if (ifeng.le.1) then
         do j=1,nbini
            eprime=xss(loc+1)
            if (ifeng.eq.0) then
               wt=1
            else
               wt=10
               if (j.eq.1.or.j.eq.nbini) wt=1
               if (j.eq.2.or.j.eq.nbini-1) wt=4
            endif
            ebar=ebar+wt*eprime
            sum=sum+wt
            loc=loc+nang+1
         enddo
      else
         cdl=0
         xl=0
         pl=0
         loc=loc+nang+3
         do j=2,nbini
            x=xss(loc+1)
            p=xss(loc+2)
            u=(pl-(p-pl)*xl/(x-xl))*(x**2-xl**2)/2&
              +((p-pl)/(x-xl))*(x**3-xl**3)/3
            ul=(x-xl)*(p+pl)/2
            u=u/ul
            un=xss(loc+3)-cdl
            ebar=ebar+u*un
            sum=sum+un
            xl=x
            pl=p
            cdl=xss(loc+3)
            loc=loc+nang+3
         enddo
      endif
      ebar=ebar/sum
      if (e.lt.xmin) xmin=e
      if (e.gt.xmax) xmax=e
      if (ebar.lt.ymin) ymin=ebar
      if (ebar.gt.ymax) ymax=ebar
   enddo
   call ascll(xmin,xmax)
   call ascll(ymin,ymax)
   if (ymin.lt.ymax/scale) ymin=ymax/scale
   write(nout,'(''1'',i3,''/'')') iwcol
   it=1
   do i=1,70
      if (hk(i:i).ne.' ') it=i
   enddo
   write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
   write(nout,'(a,''<t>hermal ebar'',a,''/'')') qu,qu
   xtag=2*xmin
   ytag=2*log10(ymin)/10+8*log10(ymax)/10
   ytag=ten**ytag
   write(nout,'(''4 0 2 1'',2e12.4,''/'')') xtag,ytag
   write(nout,'(1p,3e12.3,''/'')') xmin,xmax,one
   write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
   write(nout,'(1p,3e12.3,''/'')') ymin,ymax,one
   write(nout,'(a,''<e>bar (<m>e<v>)'',a,''/'')') qu,qu
   write(nout,'(''/'')')
   if (iwcol.eq.0) then
      write(nout,'(''0 0 1/'')')
   else
      write(nout,'(''0 0 0 1/'')')
   endif
   write(nout,'(a,''inelastic'',a,''/'')') qu,qu
   write(nout,'(''0/'')')
   loc=itxe-1
   if (ifeng.gt.1) loc=loc+2*nie
   do i=1,nie
      e=xss(itie+i)
      xnelas=xss(itie+nie+i)
      if (ifeng.gt.1) nbini=nint(xss(itxe-1+nie+i))
      ebar=0
      sum=0
      if (ifeng.le.1) then
      do j=1,nbini
         eprime=xss(loc+1)
         if (ifeng.eq.0) then
            wt=1
         else
            wt=10
            if (j.eq.1.or.j.eq.nbini) wt=1
            if (j.eq.2.or.j.eq.nbini-1) wt=4
         endif
         ebar=ebar+wt*eprime
         sum=sum+wt
         loc=loc+nang+1
      enddo
      else
         cdl=0
         xl=0
         pl=0
         loc=loc+nang+3
         do j=2,nbini
            x=xss(loc+1)
            p=xss(loc+2)
            u=(pl-(p-pl)*xl/(x-xl))*(x**2-xl**2)/2&
              +((p-pl)/(x-xl))*(x**3-xl**3)/3
            ul=(x-xl)*(p+pl)/2
            u=u/ul
            un=xss(loc+3)-cdl
            ebar=ebar+u*un
            sum=sum+un
            xl=x
            pl=p
            cdl=xss(loc+3)
            loc=loc+nang+3
         enddo
      endif
      ebar=ebar/sum
      write(nout,'(1p,2e14.6,''/'')') e,ebar
   enddo
   write(nout,'(''/'')')
   if (ifeng.le.1) then
      write(nout,'(''99/'')')
      return
   endif

   !--3-d plots of thermal inelastic distributions
   !--when continuous probability distribution option is used

   ! plot energy distributions for low incident energies
   nang=nil-1
   loc=itxe-1
   loc=loc+2*nie
   zmin=1000
   zmax=0
   xmin=5/scale/100000
   xmax=.5/scale
   ymin=1/scale/100000
   ymax=2/scale/100
   do ie=1,nie
      e=xss(itie+ie)
      nbini=nint(xss(itxe-1+nie+ie))
      do j=1,nbini
         ep=xss(loc+1)
         if (e.ge.ymin.and.e.le.ymax.and.&
           ep.gt.xmin.and.ep.lt.xmax) then
            p=xss(loc+2)
            if (p.lt.zmin) zmin=p
            if (p.gt.zmax) zmax=p
         endif
         loc=loc+nang+3
      enddo
   enddo
   if (zmax.gt.zero) then
      zmin=zmax/1000
      call ascll(zmin,zmax)
      write(nout,'(''1'',i3,''/'')') iwcol
      write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
      write(nout,'(a,''thermal inelastic'',a,''/'')') qu,qu
      write(nout,'(''-4 2/'')')
      write(nout,'(1p,3e12.3,''/'')') xmin,xmax,1.
      write(nout,'(a,''<s>ec. <e>nergy'',a,''/'')') qu,qu
      write(nout,'(1p,3e12.3,''/'')') ymin,ymax,1.
      write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
      write(nout,'(1p,3e12.3,''/'')') zmin,zmax,1.
      write(nout,'(a,''<p>rob/<m>e<v>'',a,''/'')') qu,qu
      write(nout,'(''/'')')
      write(nout,'(''10 -15 15 3.5 5.5 2.5/'')')
      write(nout,'(''1/'')')
      loc=itxe-1
      loc=loc+2*nie
      do ie=1,nie
         e=xss(itie+ie)
         if (e.ge.ymin.and.e.le.ymax) write(nout,'(1p,e14.6,''/'')') e
         nbini=nint(xss(itxe-1+nie+ie))
         do j=1,nbini
            ep=xss(loc+1)
            p=xss(loc+2)
            if (p.lt.zmin) p=zmin
            if (e.ge.ymin.and.e.le.ymax.and.&
              ep.gt.xmin.and.ep.le.xmax) then
               write(nout,'(1p,2e14.6,''/'')') ep,p
            endif
            loc=loc+nang+3
         enddo
         if (e.ge.ymin.and.e.le.ymax) write(nout,'(''/'')')
      enddo
      write(nout,'(''/'')')
   endif

   ! plot energy distributions for middle incident energies
   nang=nil-1
   loc=itxe-1
   loc=loc+2*nie
   zmin=1000
   zmax=0
   xmin=2/scale/1000
   xmax=5/scale/10
   ymin=2/scale/100
   ymax=2/scale/10
   do ie=1,nie
      e=xss(itie+ie)
      nbini=nint(xss(itxe-1+nie+ie))
      epl=0
      do j=1,nbini
         ep=xss(loc+1)
         if (e.ge.ymin.and.e.le.ymax.and.&
           ep.gt.xmin.and.ep.lt.xmax) then
            p=xss(loc+2)
            if (p.lt.zmin) zmin=p
            if (p.gt.zmax) zmax=p
         endif
         loc=loc+nang+3
      enddo
   enddo
   if (zmax.gt.zero) then
      zmin=zmax/500
      call ascll(zmin,zmax)
      write(nout,'(''1'',i3,''/'')') iwcol
      write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
      write(nout,'(a,''thermal inelastic'',a,''/'')') qu,qu
      write(nout,'(''-4 2/'')')
      write(nout,'(1p,3e12.3,''/'')') xmin,xmax,1.
      write(nout,'(a,''<s>ec. <e>nergy'',a,''/'')') qu,qu
      write(nout,'(1p,3e12.3,''/'')') ymin,ymax,1.
      write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
      write(nout,'(1p,3e12.3,''/'')') zmin,zmax,1.
      write(nout,'(a,''<p>rob/<m>e<v>'',a,''/'')') qu,qu
      write(nout,'(''/'')')
      write(nout,'(''10 -15 15 3.5 5.5 2.5/'')')
      write(nout,'(''1/'')')
      loc=itxe-1
      loc=loc+2*nie
      do ie=1,nie
         e=xss(itie+ie)
         if (e.ge.ymin.and.e.le.ymax) write(nout,'(1p,e14.6,''/'')') e
         nbini=nint(xss(itxe-1+nie+ie))
         do j=1,nbini
            ep=xss(loc+1)
            p=xss(loc+2)
            if (p.lt.zmin) p=zmin
            if (e.ge.ymin.and.e.le.ymax.and.&
              ep.gt.xmin.and.ep.le.xmax) then
               write(nout,'(1p,2e14.6,''/'')') ep,p
            endif
            loc=loc+nang+3
         enddo
         if (e.ge.ymin.and.e.le.ymax) write(nout,'(''/'')')
      enddo
      write(nout,'(''/'')')
   endif

   ! plot energy distributions for higher incident energies
   nang=nil-1
   loc=itxe-1
   loc=loc+2*nie
   zmin=1000
   zmax=0
   xmin=2/scale/100
   xmax=2/scale
   ymin=2/scale/10
   ymax=2/scale
   do ie=1,nie
      e=xss(itie+ie)
      nbini=nint(xss(itxe-1+nie+ie))
      do j=1,nbini
         ep=xss(loc+1)
         if (e.ge.ymin.and.e.le.ymax.and.&
           ep.gt.xmin.and.ep.lt.xmax) then
            p=xss(loc+2)
            if (p.lt.zmin) zmin=p
            if (p.gt.zmax) zmax=p
         endif
         loc=loc+nang+3
      enddo
   enddo
   if (zmax.gt.zero) then
      zmin=zmax/500
      call ascll(zmin,zmax)
      write(nout,'(''1'',i3,''/'')') iwcol
      write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
      write(nout,'(a,''thermal inelastic'',a,''/'')') qu,qu
      write(nout,'(''-4 2/'')')
      write(nout,'(1p,3e12.3,''/'')') xmin,xmax,1.
      write(nout,'(a,''<s>ec. <e>nergy'',a,''/'')') qu,qu
      write(nout,'(1p,3e12.3,''/'')') ymin,ymax,1.
      write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
      write(nout,'(1p,3e12.3,''/'')') zmin,zmax,1.
      write(nout,'(a,''<p>rob/<m>e<v>'',a,''/'')') qu,qu
      write(nout,'(''/'')')
      write(nout,'(''10 -15 15 3.5 5.5 2.5/'')')
      write(nout,'(''1/'')')
      loc=itxe-1
      loc=loc+2*nie
      do ie=1,nie
         e=xss(itie+ie)
         if (e.ge.ymin.and.e.le.ymax) write(nout,'(1p,e14.6,''/'')') e
         nbini=nint(xss(itxe-1+nie+ie))
         do j=1,nbini
            ep=xss(loc+1)
            p=xss(loc+2)
            if (p.lt.zmin) p=zmin
            if (e.ge.ymin.and.e.le.ymax.and.&
              ep.gt.xmin.and.ep.le.xmax) then
               write(nout,'(1p,2e14.6,''/'')') ep,p
            endif
            loc=loc+nang+3
         enddo
         if (e.ge.ymin.and.e.le.ymax) write(nout,'(''/'')')
      enddo
      write(nout,'(''/'')')
   endif

   ! plot energy distributions for highest incident energies
   nang=nil-1
   loc=itxe-1
   loc=loc+2*nie
   zmin=1000
   zmax=0
   xmin=1/scale/100
   xmax=10/scale
   ymin=2/scale
   ymax=10/scale
   do ie=1,nie
      e=xss(itie+ie)
      nbini=nint(xss(itxe-1+nie+ie))
      do j=1,nbini
         ep=xss(loc+1)
         if (e.ge.ymin.and.e.le.ymax.and.&
           ep.gt.xmin.and.ep.lt.xmax) then
            p=xss(loc+2)
            if (p.lt.zmin) zmin=p
            if (p.gt.zmax) zmax=p
         endif
         loc=loc+nang+3
      enddo
   enddo
   if (zmax.gt.zero) then
      zmin=zmax/500
      call ascll(zmin,zmax)
      write(nout,'(''1'',i3,''/'')') iwcol
      write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
      write(nout,'(a,''thermal inelastic'',a,''/'')') qu,qu
      write(nout,'(''-4 2/'')')
      write(nout,'(1p,3e12.3,''/'')') xmin,xmax,1.
      write(nout,'(a,''<s>ec. <e>nergy'',a,''/'')') qu,qu
      write(nout,'(1p,3e12.3,''/'')') ymin,ymax,1.
      write(nout,'(a,''<e>nergy (<m>e<v>)'',a,''/'')') qu,qu
      write(nout,'(1p,3e12.3,''/'')') zmin,zmax,1.
      write(nout,'(a,''<p>rob/<m>e<v>'',a,''/'')') qu,qu
      write(nout,'(''/'')')
      write(nout,'(''10 -15 15 3.5 5.5 2.5/'')')
      write(nout,'(''1/'')')
      loc=itxe-1
      loc=loc+2*nie
      do ie=1,nie
         e=xss(itie+ie)
         if (e.ge.ymin.and.e.le.ymax) write(nout,'(1p,e14.6,''/'')') e
         nbini=nint(xss(itxe-1+nie+ie))
         do j=1,nbini
            ep=xss(loc+1)
            p=xss(loc+2)
            if (p.lt.zmin) p=zmin
            if (e.ge.ymin.and.e.le.ymax.and.&
              ep.gt.xmin.and.ep.le.xmax) then
               write(nout,'(1p,2e14.6,''/'')') ep,p
            endif
            loc=loc+nang+3
         enddo
         if (e.ge.ymin.and.e.le.ymax) write(nout,'(''/'')')
      enddo
      write(nout,'(''/'')')
   endif

   ! plot angle-energy distribution for several incident energies
   ie=19
   iskip=(nie-ie)/4
   do while (ie.lt.nie)
      nang=nil-1
      loc=nint(xss(itxe-1+ie))
      nbini=nint(xss(itxe-1+nie+ie))
      zmin=1
      zmin=zmin/100
      zmax=10
      ymin=10
      ymax=1/scale/10000000
      xmin=-1
      xmax=+1
      e=xss(itie+ie)
      epl=0
      loc=loc+nang+3
      do j=2,nbini-1
         ep=(xss(loc+1)+epl)/2
         if (xss(loc+1).gt.epl) then
            if (ep.lt.ymin) ymin=ep
            if (ep.gt.ymax) ymax=ep
            epl=xss(loc+1)
            ul=-1
            do k=1,nang
               u=xss(loc+3+k)
               if (k.lt.nang) then
                  un=(u+xss(loc+3+k+1))/2
               else
                  un=1
               endif
               if (k.eq.1.and.u-ul.gt.5*(un-u)) ul=u-3*(un-u)
               if (k.eq.nang.and.un-u.gt.5*(u-ul)) un=u+3*(u-ul)
               p=1
               p=p/nang
               p=p/(un-ul)
               ul=un
            enddo
         endif
         loc=loc+nang+3
      enddo
      call ascll(ymin,ymax)
      write(nout,'(''1'',i3,''/'')') iwcol
      write(nout,'(a,''<'',a,''>'',a,''/'')') qu,hk(1:it),qu
      write(nout,'(a,''thermal inelastic for e='',1p,e10.3,'' MeV'',a,''/'')')&
        qu,xss(itie+ie),qu
      write(nout,'(''-2 2/'')')
      write(nout,'(1p,3e12.3,''/'')') xmin,xmax,.5
      write(nout,'(a,''<c>osine'',a,''/'')') qu,qu
      write(nout,'(1p,3e12.3,''/'')') ymin,ymax,1.
      write(nout,'(a,''<s>ec. <e>nergy'',a,''/'')') qu,qu
      write(nout,'(1p,3e12.3,''/'')') zmin,zmax,1.
      write(nout,'(a,''<p>rob/cosine'',a,''/'')') qu,qu
      write(nout,'(''/'')')
      write(nout,'(''10 -15 15 3.5 5.5 2.5/'')')
      write(nout,'(''1/'')')
      loc=nint(xss(itxe-1+ie))
      nbini=nint(xss(itxe-1+nie+ie))
      skip=1+log10(ymax/ymin)/60
      epl=0
      loc=loc+nang+3
      do j=2,nbini-1
         ep=(xss(loc+1)+epl)/2
         if (xss(loc+1).gt.epl) then
            if (j.eq.2.or.ep/epl.gt.skip) then
               write(nout,'(1p,e14.6,''/'')') ep
               ul=-1
               write(nout,'(1p,2e14.6,''/'')') -1.,zmin
               do k=1,nang
                  u=xss(loc+3+k)
                  if (k.lt.nang) then
                     un=(u+xss(loc+3+k+1))/2
                  else
                     un=1
                  endif
                  if (k.eq.1.and.u-ul.gt.5*(un-u)) ul=u-3*(un-u)
                  if (k.eq.nang.and.un-u.gt.5*(u-ul)) un=u+3*(u-ul)
                  p=1
                  p=p/nang
                  p=p/(un-ul)
                  if (k.eq.1) write(nout,'(1p,2e14.6,''/'')') ul,zmin
                  write(nout,'(1p,2e14.6,''/'')') u,p
                  if (k.eq.nang) write(nout,'(1p,2e14.6,''/'')') un,zmin
                  ul=un
               enddo
               write(nout,'(1p,2e14.6,''/'')') 1.,zmin
               write(nout,'(''/'')')
               epl=xss(loc+1)
            endif
         endif
         loc=loc+nang+3
      enddo
      write(nout,'(''/'')')
      ie=ie+iskip
   enddo

   write(nout,'(''99/'')')

   return
   end subroutine tplots

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

   subroutine throut(itype,nout,ndir,mcnpx,hk,izn,awn)
   !-------------------------------------------------------------------
   ! Write thermal ACE data to output and directory files.
   !-------------------------------------------------------------------
   use util ! provides openz,closz
   ! externals
   integer::itype,nout,ndir,mcnpx
   integer::izn(16)
   real(kr)::awn(16)
   character(70)::hk
   ! internals
   integer::l,ne,n,nexe,nern,lrec,ll,nn,i
   integer::ner=1
   integer::nbw=1

   !--open ace file
   if (itype.eq.1) call openz(nout,1)
   if (itype.eq.2) call openz(-nout,1)

   !--branch according to ace file type
   if (itype.eq.1) then

      !--write type-1 header block
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
      write(nout,'(8i9)') &
        len2,idpni,nil,nieb,idpnc,ncl,ifeng,nxsd,&
        itie,itix,itxe,itce,itcx,itca,jxsd

      !--itie block
      l=itie
      ne=nint(xss(l))
      call typen(l,nout,1)
      l=l+1
      n=2*ne
      do i=1,n
         call typen(l,nout,2)
         l=l+1
      enddo

      !--itxe block
      l=itxe
      if (ifeng.le.1) then
         n=ne*nieb*(nil+2)
         do i=1,n
            call typen(l,nout,2)
            l=l+1
         enddo
      else
         n=2*ne
         do i=1,ne
            n=n+nint(xss(l+ne+i-1))*(nil+2)
         enddo
         do i=1,n
            call typen(l,nout,2)
            l=l+1
         enddo
      endif

      !--itce block
      if (itce.ne.0) then
         l=itce
         ne=nint(xss(l))
         nexe=ne
         call typen(l,nout,1)
         l=l+1
         n=2*ne
         do i=1,n
            call typen(l,nout,2)
            l=l+1
         enddo
      endif

      !--itca block
      if (itce.ne.0.and.ncl.ne.-1) then
         n=nexe*(ncl+1)
         do i=1,n
            call typen(l,nout,2)
            l=l+1
         enddo
      endif

      !--finished
      call typen(l,nout,3)
      nern=0
      lrec=0
      call closz(nout)

   !--write ace file in type 2 format
   else if (itype.eq.2) then
      if (mcnpx.eq.0) then
         write(nout) hz(1:10),aw0,tz,hd,hk,hm,&
            (izn(i),awn(i),i=1,16),&
           len2,idpni,nil,nieb,idpnc,ncl,ifeng,nxsd,&
           itie,itix,itxe,itce,itcx,itca,jxsd
      else
         write(nout) hz(1:13),aw0,tz,hd,hk,hm,&
            (izn(i),awn(i),i=1,16),&
           len2,idpni,nil,nieb,idpnc,ncl,ifeng,nxsd,&
           itie,itix,itxe,itce,itcx,itca,jxsd
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
      call closz(-nout)
   endif
   call openz(ndir,1)
   if (mcnpx.eq.0) then
      write(ndir,&
        '(a10,f12.6,'' filename route'',i2,'' 1 '',i9,2i6,1p,e10.3)')&
        hz(1:10),aw0,itype,len2,lrec,nern,tz
   else
      write(ndir,&
        '(a13,f12.6,'' filename route'',i2,'' 1 '',i9,2i6,1p,e10.3)')&
        hz(1:13),aw0,itype,len2,lrec,nern,tz
   endif
   call closz(ndir)
   return
   end subroutine throut

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

end module aceth

