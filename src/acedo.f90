module acedo
   ! provides ace dosimetry formats for acer
   use locale
   use acecm, only: xss,nxss
   implicit none
   private

   !--Public routines
   public acedos,dosfix

   !--Private global variables

   ! incident particle
   real(kr)::awi
   integer::izai

   ! header info for dosimetry format
   character(13)::hz
   character(10)::hd
   character(10)::hm
   real(kr)::aw0,tz

   ! parameters for dosimetry nxs block
   integer::len2,za,nxs3,ntr,nxsd(12)

   ! parameters for dosimetry jxs block
   integer::lone,jxs2,mtr,jxs4,jxs5,lsig,sigd,jxsd(14),end,jxsd2(10)

contains

   subroutine acedos(matd,tempd,nin,nace,ndir,itype,iprint,mcnpx,&
     suff,hk,izn,awn)
   !-------------------------------------------------------------------
   ! Load dosimetry data in ACE format.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use util ! provides openz,closz,dater
   use endf ! provides endf routines and variables
   use physics ! provides bk
   ! externals
   integer::matd,nin,nace,ndir,itype,iprint,mcnpx
   real(kr)::tempd,suff
   integer::izn(16)
   real(kr)::awn(16)
   character(70)::hk
   ! internals
   integer::nwscr,nb,nw,j,l,nr,i,ne,k,nmax,jscr,intr
   integer::l1,nn,n,is,ns,lfs
   real(kr)::temp,awr,zaid
   character(8)::hdt
   real(kr),dimension(:),allocatable::scr,tmpscr
   real(kr),parameter::emev=1.e6_kr

   !--allocate main container array
   nxs3=0
   nxsd=0
   jxs2=0
   jxs4=0
   jxs5=0
   jxsd2=0
   xss=0

   !--allocate scratch storage
   nwscr=250000
   allocate(scr(nwscr))

   !--determine what endf version is being used
   call openz(nin,0)
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

   !--find desired temperature on endf tape
   call repoz(nin)
   call tpidio(nin,0,0,scr,nb,nw)
   temp=0
   do while (abs(temp-tempd).ge.tempd/100+1)
      call contio(nin,0,0,scr,nb,nw)
      if (math.lt.0)&
        call error('acedos','desired mat and temp not found.',' ')
      za=nint(scr(1))
      awr=scr(2)
      if (iverf.ge.5) call contio(nin,0,0,scr,nb,nw)
      if (iverf.eq.6) then
         call contio(nin,0,0,scr,nb,nw)
         awi=scr(1)
         izai=nint(scr(5)/10)
      else
         awi=0
         izai=1 !assume incident neutron if not a version 6 file
      endif

      call contio(nin,0,0,scr,nb,nw)
      temp=scr(1)
      if (abs(temp-tempd).ge.tempd/100+1) then
         call tomend(nin,0,0,scr)
      endif
   enddo

   !--locate first reaction in file 3 or file 10
   do while (mfh.ne.3 .and. mfh.ne.10)
      call tofend(nin,0,0,scr)
      call contio(nin,0,0,scr,nb,nw)
      if (mfh.gt.10 .or. math.le.0)&
         call error('acedos','no xs data for desired mat.',' ')
   enddo

   !--define locators for dosimetry data
   nmax=350
   lone=1
   mtr=lone
   lsig=mtr+nmax
   sigd=lsig+nmax
   l=sigd
   j=1
   if (mfh.gt.3) go to 110

   !--loop over file 3 reactions on nin
   do while (mfh.ne.0)
      xss(lsig-1+j)=l
      if (mfh.ne.0) then
         if (mth.ne.1) then
            xss(mtr-1+j)=mth
            if(allocated(tmpscr)) deallocate(tmpscr)
            allocate(tmpscr(npage+46))
            jscr=1
            call tab1io(nin,0,0,tmpscr(jscr),nb,nw)
            if (nb.gt.nwscr) then
               nwscr=nb+nw+1
               if(allocated(scr)) deallocate(scr)
               allocate(scr(nwscr))
            endif
            scr(1:nw)=tmpscr(1:nw)
            nr=nint(scr(5))
            ne=nint(scr(6))
            intr=nint(scr(8))
            jscr=jscr+nw
            if (jscr.gt.nwscr)&
               call error('acedos','scr array storage exceeded',' ')
            do while (nb.ne.0)
               call moreio(nin,0,0,scr(jscr),nb,nw)
               jscr=jscr+nw
               if (jscr.gt.nwscr)&
                  call error('acedos','scr array storage exceeded',' ')
            enddo
            if (nr.ne.1.or.intr.ne.2) then
               xss(l)=nr
               l=l+1
               do i=1,nr
                  xss(l+i-1)=scr(5+2*i)
                  xss(l+nr+i-1)=scr(6+2*i)
               enddo
               l=l+2*nr
            else
               xss(l)=0
               l=l+1
            endif
            xss(l)=ne
            k=7+2*nr
            l=l+1
            do i=1,ne
               xss(l)=scr(k)/emev
               xss(l+ne)=scr(k+1)
               l=l+1
               k=k+2
            enddo
            l=l+ne
            j=j+1
            if (j.gt.nmax) call error('acedos','too many reactions',&
                                      'need larger nmax')
         endif
         call tosend(nin,0,0,scr)
         call contio(nin,0,0,scr,nb,nw)
      endif
   enddo

   !--locate first reaction in file 10
   call contio(nin,0,0,scr,nb,nw)
   if (mfh.gt.10 .or. math.le.0) go to 120
   do while (mfh.ne.10)
      call tofend(nin,0,0,scr)
      call contio(nin,0,0,scr,nb,nw)
      if (mfh.gt.10 .or. math.le.0) go to 120
   end do
  110 continue
   !--loop over reactions on nin
      do while (mfh.ne.0)
         ns=max(1,n1h)
         xss(lsig-1+j)=l
         if (mfh.ne.0) then
            if (mth.ne.1) then
               do is=1,ns
                  if(allocated(tmpscr)) deallocate(tmpscr)
                  allocate(tmpscr(npage+46))
                  jscr=1
                  call tab1io(nin,0,0,tmpscr(jscr),nb,nw)
                  if (nb.gt.nwscr) then
                     nwscr=nb+nw+1
                     if(allocated(scr)) deallocate(scr)
                     allocate(scr(nwscr))
                  endif
                  scr(1:nw)=tmpscr(1:nw)
                  lfs=l2h
                  xss(mtr-1+j)=mth+1000*(10+lfs)
                  nr=nint(scr(5))
                  ne=nint(scr(6))
                  intr=nint(scr(8))
                  jscr=jscr+nw
                  if (jscr.gt.nwscr)&
                     call error('acedos','array storage exceeded',' ')
                  do while (nb.ne.0)
                     call moreio(nin,0,0,scr(jscr),nb,nw)
                     jscr=jscr+nw
                     if (jscr.gt.nwscr)&
                        call error('acedos','array storage exceeded',&
                                   ' ')
                  enddo
                  if (nr.ne.1.or.intr.ne.2) then
                     xss(l)=nr
                     l=l+1
                     do i=1,nr
                        xss(l+2*i-2)=scr(5+2*i)
                        xss(l+2*i-1)=scr(6+2*i)
                     enddo
                     l=l+2*nr
                  else
                     xss(l)=0
                     l=l+1
                  endif
                  xss(l)=ne
                  k=7+2*nr
                  l=l+1
                  do i=1,ne
                     xss(l)=scr(k)/emev
                     xss(l+ne)=scr(k+1)
                     l=l+1
                     k=k+2
                  enddo
                  l=l+ne
                  j=j+1
                  if (j.gt.nmax)&
                     call error('acedos','too many reactions',&
                                         'need larger nmax')
                  if (is.ne.ns) xss(lsig-1+j)=l
               enddo
            endif
            call tosend(nin,0,0,scr)
            call contio(nin,0,0,scr,nb,nw)
         endif
      enddo
  120 continue
   ntr=j-1

   !--squeeze out space in mtr and lsig
   l1=nint(xss(lsig))
   nn=l1+mtr-2
   do i=1,ntr
      xss(mtr+ntr-1+i)=xss(lsig-1+i)-nn
   enddo
   lsig=mtr+ntr
   n=l-sigd
   do i=1,n
      xss(mtr+2*ntr-1+i)=xss(sigd-1+i)
   enddo
   sigd=mtr+2*ntr
   len2=sigd+n-1
   end=len2

   !--dosimetry load is finished
   call closz(nin)
   deallocate(scr)
   zaid=za+suff
   if (mcnpx.eq.0) then
      write(hz,'(f9.2,''y'')') zaid
   else
      write(hz,'(f10.3,''ny '')') zaid
   endif
   tz=temp*bk/emev
   call dater(hdt)
   hd='  '//hdt
   write(hm,'(''   mat'',i4)') matd
   aw0=awr
   if (iprint.gt.0) call dosprt(hk)
   call dosout(itype,nace,ndir,mcnpx,hk,izn,awn)
   return
   end subroutine acedos

   subroutine dosfix(itype,nin,nout,ndir,iprint,nplot,mcnpx,&
     suff,nxtra,hk,izn,awn)
   !-------------------------------------------------------------------
   ! Print and or edit an ACE dosimetry file.
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
        len2,za,nxs3,ntr,nxsd(1:12),&
        lone,jxs2,mtr,jxs4,jxs5,lsig,sigd,jxsd(1:14),end,jxsd2(1:10)
      n=(lone+3)/4
      l=0
      do i=1,n
         read (nin,'(4e20.0)') (xss(l+j),j=1,4)
         l=l+4
      enddo

   !--read type 2 ace format file
   else if (itype.eq.2) then
      if (mcnpx.eq.0) then
        read(nin) hz(1:10),aw0,tz,hd,hko,hm,(izo(i),awo(i),i=1,16),&
          len2,za,nxs3,ntr,nxsd(1:12),&
          lone,jxs2,mtr,jxs4,jxs5,lsig,sigd,jxsd(14),end,jxsd2(10)
      else
        read(nin) hz(1:13),aw0,tz,hd,hko,hm,(izo(i),awo(i),i=1,16),&
          len2,za,nxs3,ntr,nxsd(1:12),&
          lone,jxs2,mtr,jxs4,jxs5,lsig,sigd,jxsd(14),end,jxsd2(10)
      endif
      n=(lone+ner-1)/ner
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
   if (iprint.gt.0) call dosprt(hk)
   if (nout.gt.0) call dosout(itype,nout,ndir,mcnpx,hk,izn,awn)

   return
   end subroutine dosfix

   subroutine dosprt(hk)
   !-------------------------------------------------------------------
   ! Print dosimetry cross sections from data in memory.
   !-------------------------------------------------------------------
   use mainio ! provides nsyso
   use acecm  ! provides mtname
   ! externals
   character(70)::hk
   ! internals
   integer::j,mt,l,nr,lim,i,ne,k,line
   real(kr)::e,s
   character(10)::name
   character(14)::col(6)

   !--print information block
   write(nsyso,'(''1''/////&
     &38x,''zaid'',1x,a13/&
     &6x,''***********************'',10x,''awr'',f10.3/&
     &6x,''*                     *'',9x,''temp'',1p,e10.2/&
     &6x,''*      dosimetry      *'',9x,''date'',a10/&
     &6x,''*                     *'',10x,''mat'',a10/&
     &6x,''*   ace format file   *''/&
     &6x,''*                     *'',9x,''len2'',i10/&
     &6x,''*     processed by    *'',11x,''za'',i10/&
     &6x,''*                     *''/&
     &6x,''*        njoy         *'',10x,''ntr'',i10/&
     &6x,''*                     *''/&
     &6x,''***********************'',9x,''lone'',i10/&
     &39x,''mtr'',i10/38x,''lsig'',i10/38x,''sigd'',i10///&
     &6x,''hk>>> '',a70)')&
     hz,aw0,tz,hd,hm,len2,za,ntr,lone,mtr,lsig,sigd,hk

   !--print sigd block
   do j=1,ntr
      mt=nint(xss(mtr-1+j))
      call mtname(mt,name,izai)
      l=nint(xss(lsig-1+j)+sigd-1)
      nr=nint(xss(l))
      l=l+1
      if (nr.ne.0) then
         lim=2*nr
         write(nsyso,'(''1''///&
           &'' reaction mt = '',i6,3x,a10/'' interpolation: '',12i6)')&
           mt,name,(nint(xss(l-1+i)),i=1,lim)
         l=l+2*nr
      else
         write(nsyso,'(''1''///&
           &'' reaction mt = '',i6,3x,a10/'' linear interpolation'')')&
           mt,name
      endif
      ne=nint(xss(l))
      l=l+1
      k=1
      line=1
      do i=1,ne
         e=xss(l)
         s=xss(l+ne)
         l=l+1
         write(col(k),'(1p,e14.4)') e
         write(col(k+1),'(1p,e14.4)') s
         k=k+2
         if (k.ge.7) then
            if (line.eq.1) write(nsyso,'(/&
              &1x,3(8x,''energy'',10x,''xsec'')/&
              &1x,3(8x,''------'',10x,''----''))')
            write(nsyso,'(1x,6a14)') col
            k=1
            line=line+1
         endif
      enddo
      k=k-1
      if (k.ne.0) then
         write(nsyso,'(1x,6a14)') (col(i),i=1,k)
      endif
   enddo
   return
   end subroutine dosprt

   subroutine dosout(itype,nout,ndir,mcnpx,hk,izn,awn)
   !-------------------------------------------------------------------
   ! Write out the dosimetry file.
   !-------------------------------------------------------------------
   use util  ! provides openz,closz,error
   use acecm ! provides write routines
   ! externals
   integer::itype,nout,ndir,mcnpx
   integer::izn(16)
   real(kr)::awn(16)
   character(70)::hk
   ! internals
   integer::l,nr,n,j,ne,ll,nn,i,nern,lrec

   integer,parameter::ner=1
   integer,parameter::nbw=1

   !--process according to ace file type
   if (itype.eq.1) call openz(nout,1)
   if (itype.eq.2) call openz(-nout,1)

   !--type 1
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
      write(nout,'(8i9)')&
        len2,za,nxs3,ntr,nxsd(1:12),&
        lone,jxs2,mtr,jxs4,jxs5,lsig,sigd,jxsd(1:14),end,jxsd2(1:10)

      !--mtr block
      l=mtr
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

      !--sigd block
      l=sigd
      do i=1,ntr
         nr=nint(xss(l))
         call typen(l,nout,1)
         l=l+1
         if (nr.ne.0) then
            n=2*nr
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
      call typen(0,nout,3)
      call closz(nout)
      nern=0
      lrec=0

   !--write ace file in type 2 format
   else
      if (mcnpx.eq.0) then
         write(nout) hz(1:10),aw0,tz,hd,hk,hm,(izn(i),awn(i),i=1,16),&
           len2,za,nxs3,ntr,nxsd(1:12),&
           lone,jxs2,mtr,jxs4,jxs5,lsig,sigd,jxsd(1:14),end,jxsd2(1:10)
      else
         write(nout) hz(1:13),aw0,tz,hd,hk,hm,(izn(i),awn(i),i=1,16),&
           len2,za,nxs3,ntr,nxsd(1:12),&
           lone,jxs2,mtr,jxs4,jxs5,lsig,sigd,jxsd(1:14),end,jxsd2(1:10)
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
   end subroutine dosout

end module acedo
