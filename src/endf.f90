module endf
   ! Routines for handling ENDF formats.
   use locale
   implicit none
   private

   !--Public routines
   public contio,listio,tab1io,tab2io,moreio
   public tpidio,hdatio,dictio,intgio
   public tosend,tofend,tomend,totend
   public asend,afend,amend,atend
   public skip6,findf
   public terp1,terpa,intega
   public gety1,gety2
   public gral
   public a11

   !--Public variables
   integer,public::npage=306
   integer,public::iverf
   real(kr),public::c1h,c2h
   integer,public::l1h,l2h,n1h,n2h,math,mfh,mth,nsh,nsp,nsc
   real(kr),public::thr6

   !--Private global variables
   integer::nah
   real(kr)::ah(6)
   integer::ihol
   character(4)::hol(17)
   real(kr)::rhol(17)
   equivalence(hol(1),rhol(1))
   integer::nr
   integer::nbt(20),jnt(20)

contains

   subroutine contio(nin,nout,nscr,a,nb,nw)
   !--------------------------------------------------------------------
   ! Utility routine for ENDF coded and blocked binary tapes.
   ! Read, write, and/or convert one control record.
   ! Positive units are coded, negative ones are blocked binary.
   ! If any unit is zero, it is not used.
   !--------------------------------------------------------------------
   ! externals
   integer::nin,nout,nscr,nb,nw
   real(kr)::a(*)
   ! internals
   integer::inin,inout,inscr,i
   character(11)::field(2)

   !--input
   if (nin.lt.0) then
      inin=iabs(nin)
      read(inin) math,mfh,mth,nb,nw,(a(i),i=1,6)
   else if (nin.gt.0) then
      read(nin,'(6e11.0,i4,i2,i3,i5)') (a(i),i=1,6),math,mfh,mth,nsp
   endif

   !--output
   nb=0
   nw=6
   c1h=a(1)
   c2h=a(2)
   l1h=nint(a(3))
   l2h=nint(a(4))
   n1h=nint(a(5))
   n2h=nint(a(6))
   a(3)=l1h
   a(4)=l2h
   a(5)=n1h
   a(6)=n2h
   if (nout.eq.0.and.nscr.eq.0) return
   inout=iabs(nout)
   if (nout.lt.0) then
      write(inout) math,mfh,mth,nb,nw,(a(i),i=1,6)
      inout=0
   endif
   inscr=iabs(nscr)
   if (nscr.lt.0) then
      write(inscr) math,mfh,mth,nb,nw,(a(i),i=1,6)
      inscr=0
   endif
   if (nout.le.0.and.nscr.le.0) return

   !--special path to write blanks on *end* cards
   if (math.eq.-1) call atend(inout,inscr)
   if (math.eq.-1) return
   if (math.eq.0) call amend(inout,inscr)
   if (math.eq.0) return
   if (mfh.eq.0) call afend(inout,inscr)
   if (mfh.eq.0) return
   if (mth.eq.0) call asend(inout,inscr)
   if (mth.eq.0) return

   !--format the output
   call a11(c1h,field(1))
   call a11(c2h,field(2))
   if (nscr.gt.0) then
      write(nscr,'(2a11,4i11,i4,i2,i3,i5)') field(1),field(2),&
        l1h,l2h,n1h,n2h,math,mfh,mth,nsc
      nsc=nsc+1
      if (nsc.gt.99999) nsc=1
   endif
   if (nout.gt.0) then
      write(nout,'(2a11,4i11,i4,i2,i3,i5)') field(1),field(2),&
        l1h,l2h,n1h,n2h,math,mfh,mth,nsh
      nsh=nsh+1
      if (nsh.gt.99999) nsh=1
   endif
   return
   end subroutine contio

   subroutine listio(nin,nout,nscr,a,nb,nw)
   !--------------------------------------------------------------------
   ! Utility routine for ENDF coded and blocked binary tapes.
   ! Read, write, and/or convert the first page or record of a
   ! LIST structure.  If any unit is zero, it is not used.
   ! Positive units are coded, negative ones are blocked binary.
   !--------------------------------------------------------------------
   ! externals
   integer::nin,nout,nscr,nb,nw
   real(kr)::a(*)
   ! internals
   integer::inin,inout,inscr
   integer::i,j,k,l,iend
   integer::nbx,nwx

   ihol=0

   !--input
   if (nin.lt.0) then
      inin=iabs(nin)
      read(inin) math,mfh,mth,nb,nw,(a(i),i=1,nw)
      c1h=a(1)
      c2h=a(2)
      l1h=nint(a(3))
      l2h=nint(a(4))
      n1h=nint(a(5))
      n2h=nint(a(6))
      a(3)=l1h
      a(4)=l2h
      a(5)=n1h
      a(6)=n2h
   else if (nin.gt.0) then
      call contio(nin,0,0,a,nb,nw)
      nb=6+n1h
      iend=npage
      if (n1h.lt.npage) iend=n1h
      nw=iend+nw
      k=6
      do i=1,iend,6
         call lineio(nin,0,0)
         l=n1h-k+6
         if (l.ne.0) then
            if (l.gt.6) l=6
            do j=1,l
               k=k+1
               a(k)=ah(j)
            enddo
         endif
      enddo
      a(1)=c1h
      a(2)=c2h
      a(3)=l1h
      a(4)=l2h
      a(5)=n1h
      a(6)=n2h
   endif

   !--output
   n1h=nint(a(5))
   if (nin.eq.0.and.n1h.le.npage) nw=6+n1h
   if (nin.eq.0.and.n1h.gt.npage) nw=6+npage
   nb=6+n1h-nw
   if (nout.eq.0.and.nscr.eq.0) return
   c1h=a(1)
   c2h=a(2)
   l1h=nint(a(3))
   l2h=nint(a(4))
   n2h=nint(a(6))
   inout=iabs(nout)
   if (nout.lt.0) then
      write(inout) math,mfh,mth,nb,nw,(a(i),i=1,nw)
      inout=0
   endif
   inscr=iabs(nscr)
   if (nscr.lt.0) then
      write(inscr) math,mfh,mth,nb,nw,(a(i),i=1,nw)
      inscr=0
   endif
   call contio(0,inout,inscr,a,nbx,nwx)
   k=6
   iend=nw-6
   do i=1,iend,6
      l=iend-k+6
      if (l.ne.0) then
         if (l.gt.6) l=6
         do j=1,l
            k=k+1
            ah(j)=a(k)
         enddo
         nah=l
         call lineio(0,inout,inscr)
      endif
   enddo
   return
   end subroutine listio

   subroutine tab1io(nin,nout,nscr,a,nb,nw)
   !--------------------------------------------------------------------
   ! Utility routine for ENDF bcd and blocked binary tapes.
   ! Read, write, and/or convert the first page or record of a
   ! TAB1 structure.  If any unit is zero, it is not used.
   ! Positive units are bcd, negative ones are blocked binary.
   !--------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::nin,nout,nscr,nb,nw
   real(kr)::a(*)
   ! internals
   character(len=60)::strng
   integer::inin,inout,inscr
   integer::i,j,k,l,m,lim,kp1,nbx,nwx

   ihol=0

   !--input
   if (nin.lt.0) then
      inin=iabs(nin)
      read(inin)math,mfh,mth,nb,nw,(a(i),i=1,nw)
      c1h=a(1)
      c2h=a(2)
      l1h=nint(a(3))
      l2h=nint(a(4))
      n1h=nint(a(5))
      n2h=nint(a(6))
      if (n1h.le.0) then
         write(strng,'(" illegal TAB1, nr<=0 for mf/mt = ",i2,"/",i3)')&
                      mfh,mth
         call error('endf',strng,'')
      endif
      a(3)=l1h
      a(4)=l2h
      a(5)=n1h
      a(6)=n2h
      do i=1,n1h
         nbt(i)=nint(a(2*i+5))
         jnt(i)=nint(a(2*i+6))
      enddo
   else if (nin.gt.0) then
      call contio(nin,0,0,a,nb,nw)
      if (n1h.le.0) then
         write(strng,'(" illegal TAB1, nr<=0 for mf/mt = ",i2,"/",i3)')&
                      mfh,mth
         call error('endf',strng,'')
      endif
      nb=6+2*n1h+2*n2h
      call tablio(nin,0,0)
      k=6
      do i=1,nr
         k=k+1
         a(k)=nbt(i)
         k=k+1
         a(k)=jnt(i)
      enddo
      lim=npage
      if ((2*n2h).lt.npage) lim=2*n2h
      j=0
      do i=1,lim,6
         call lineio(nin,0,0)
         l=6
         if ((j+6).gt.lim) l=lim-j
         do m=1,l
            j=j+1
            a(k+j)=ah(m)
         enddo
      enddo
      a(1)=c1h
      a(2)=c2h
      a(3)=l1h
      a(4)=l2h
      a(5)=n1h
      a(6)=n2h
      nw=k+j
   endif
   if (nin.ne.0) then
      if (nbt(n1h).ne.n2h) then
         write(strng,'(" illegal TAB1, nbt(nr)/=np for mf/mt = ",&
                     &i2,"/",i3)')mfh,mth
         call error('endf',strng,'')
      endif
   endif

   !--output
   n1h=nint(a(5))
   n2h=nint(a(6))
   if (nin.eq.0.and.2*n2h.le.npage) nw=6+2*n1h+2*n2h
   if (nin.eq.0.and.2*n2h.gt.npage) nw=6+2*n1h+npage
   nb=6+2*n1h+2*n2h-nw
   if (nout.eq.0.and.nscr.eq.0) return
   c1h=a(1)
   c2h=a(2)
   l1h=nint(a(3))
   l2h=nint(a(4))
   do i=1,n1h
      nbt(i)=nint(a(2*i+5))
      jnt(i)=nint(a(2*i+6))
   enddo
   inout=iabs(nout)
   if (nout.lt.0) then
      write(inout) math,mfh,mth,nb,nw,(a(i),i=1,nw)
      inout=0
   endif
   inscr=iabs(nscr)
   if (nscr.lt.0) then
      write(inscr) math,mfh,mth,nb,nw,(a(i),i=1,nw)
      inscr=0
   endif
   if (nout.le.0.and.nscr.le.0) return
   call contio(0,inout,inscr,a,nbx,nwx)
   do i=1,n1h
      nbt(i)=nint(a(2*i+5))
      jnt(i)=nint(a(2*i+6))
   enddo
   nr=n1h
   call tablio(0,inout,inscr)
   k=6+2*nr
   kp1=k+1
   do i=kp1,nw,6
      l=6
      if ((k+6).gt.nw) l=nw-k
      do j=1,l
         k=k+1
         ah(j)=a(k)
      enddo
      nah=l
      call lineio(0,inout,inscr)
   enddo
   return
   end subroutine tab1io

   subroutine tab2io(nin,nout,nscr,a,nb,nw)
   !--------------------------------------------------------------------
   ! Utility routine for ENDF bcd and blocked binary tapes.
   ! Read, write, and/or convert a TAB2 structure.
   ! If any unit is zero, it is not used.
   ! Positive units are bcd, negative ones are blocked binary.
   !--------------------------------------------------------------------
   ! externals
   integer::nin,nout,nscr,nb,nw
   real(kr)::a(*)
   ! internals
   integer::inin,inout,inscr
   integer::i,k,nbx,nwx

   !--input
   if (nin.lt.0) then
      inin=iabs(nin)
      read(inin) math,mfh,mth,nb,nw,(a(i),i=1,nw)
      c1h=a(1)
      c2h=a(2)
      l1h=nint(a(3))
      l2h=nint(a(4))
      n1h=nint(a(5))
      n2h=nint(a(6))
      a(3)=l1h
      a(4)=l2h
      a(5)=n1h
      a(6)=n2h
   else if (nin.gt.0) then
      call contio(nin,0,0,a,nb,nw)
      call tablio(nin,0,0)
      a(1)=c1h
      a(2)=c2h
      a(3)=l1h
      a(4)=l2h
      a(5)=n1h
      a(6)=n2h
      k=6
      do i=1,n1h
         k=k+1
         a(k)=nbt(i)
         k=k+1
         a(k)=jnt(i)
      enddo
      nw=nw+2*n1h
   endif

   !--output
   nb=0
   if (nout.eq.0.and.nscr.eq.0) return
   c1h=a(1)
   c2h=a(2)
   l1h=nint(a(3))
   l2h=nint(a(4))
   n1h=nint(a(5))
   n2h=nint(a(6))
   do i=1,n1h
      nbt(i)=nint(a(2*i+5))
      jnt(i)=nint(a(2*i+6))
   enddo
   inout=iabs(nout)
   if (nout.lt.0) then
      write(inout) math,mfh,mth,nb,nw,(a(i),i=1,nw)
      inout=0
   endif
   inscr=iabs(nscr)
   if (nscr.lt.0) then
      write(inscr) math,mfh,mth,nb,nw,(a(i),i=1,nw)
      inscr=0
   endif
   if (nout.le.0.and.nscr.le.0) return
   call contio(0,inout,inscr,a,nbx,nwx)
   do i=1,n1h
      nbt(i)=nint(a(2*i+5))
      jnt(i)=nint(a(2*i+6))
   enddo
   nr=n1h
   call tablio(0,inout,inscr)
   return
   end subroutine tab2io

   subroutine moreio(nin,nout,nscr,a,nb,nw)
   !--------------------------------------------------------------------
   ! Utility routine for ENDF bcd and blocked binary tapes.
   ! Read, write, and/or convert the next page or record from a
   ! TAB1 or LIST structure.  If any unit is zero, it is not used.
   ! Positive units are bcd, negative ones are blocked binary.
   !--------------------------------------------------------------------
   ! externals
   integer::nin,nout,nscr,nb,nw
   real(kr)::a(*)
   ! internals
   integer::inin,inout,inscr
   integer::i,j,k,l

   !--input
   if (nin.lt.0) then
      inin=iabs(nin)
      read(inin) math,mfh,mth,nb,nw,(a(i),i=1,nw)
   else if (nin.gt.0) then
      nw=nb
      if (nw.gt.npage) nw=npage
      nb=nb-nw
      k=0
      if (ihol.ne.1) then
         ! normal data
         do i=1,nw,6
            call lineio(nin,0,0)
            l=6
            if ((k+6).gt.nw) l=nw-k
            do j=1,l
               k=k+1
               a(k)=ah(j)
            enddo
         enddo
      else
         ! hollerith descriptive data
         do i=1,nw,17
            call hollio(nin,0,0)
            do j=1,17
               k=k+1
               a(k)=rhol(j)
            enddo
         enddo
      endif
   endif

   !--output
   inout=iabs(nout)
   inscr=iabs(nscr)
   if (nin.eq.0) nw=nb
   if (nw.gt.npage) nw=npage
   if (nin.eq.0) nb=nb-nw
   if (nout.lt.0) then
      write(inout) math,mfh,mth,nb,nw,(a(i),i=1,nw)
      inout=0
   endif
   if (nscr.lt.0) then
      write(inscr) math,mfh,mth,nb,nw,(a(i),i=1,nw)
      inscr=0
   endif
   if (nout.le.0.and.nscr.le.0) return
   if (ihol.ne.1) then
      ! normal data
      k=0
      do i=1,nw,6
         l=6
         if ((l+k).gt.nw) l=nw-k
         do j=1,l
            k=k+1
            ah(j)=a(k)
         enddo
         nah=l
         call lineio(0,inout,inscr)
      enddo
   else
      ! hollerith descriptive data
      k=0
      do i=1,nw,17
         do j=1,17
            k=k+1
            rhol(j)=a(k)
         enddo
         call hollio(0,inout,inscr)
      enddo
   endif
   return
   end subroutine moreio

   subroutine tpidio(nin,nout,nscr,a,nb,nw)
   !--------------------------------------------------------------------
   ! Utility routine for ENDF bcd and binary tapes.  Read, write,
   ! and/or convert the tape identification record to/from a.  If any
   ! unit is zero, it is not used.  Positive units are bcd, and
   ! negative units are binary.
   !--------------------------------------------------------------------
   ! externals
   integer::nin,nout,nscr,nb,nw
   real(kr)::a(*)
   ! internals
   real(kr)::rb(17)
   character(4)::hb(17)
   equivalence (rb(1),hb(1))
   integer::inin,inout,inscr,i

   !--input.
   if (nin.lt.0) then
      inin=iabs(nin)
      read(inin) math,mfh,mth,nb,nw,(a(i),i=1,17)
   else if (nin.gt.0) then
      read(nin,'(16a4,a2,i4,i2,i3,i5)') (hb(i),i=1,17),math,mfh,mth
      nw=17
      do i=1,17
         a(i)=rb(i)
      enddo
   endif

   !--output
   nb=0
   inout=iabs(nout)
   inscr=iabs(nscr)
   if (nout.lt.0) then
      write(inout) math,mfh,mth,nb,nw,(a(i),i=1,17)
   endif
   if (nscr.lt.0) then
      write(inscr) math,mfh,mth,nb,nw,(a(i),i=1,17)
   endif
   if (nout.le.0.and.nscr.le.0) return
   if (nout.gt.0) then
      do i=1,17
         rb(i)=a(i)
      enddo
      write(nout,'(16a4,a2,i4,i2,i3,i5)') (hb(i),i=1,17),math,mfh,mth,nsh
      nsh=nsh+1
   endif
   if (nscr.gt.0) then
      do i=1,17
         rb(i)=a(i)
      enddo
      write(nscr,'(16a4,a2,i4,i2,i3,i5)') (hb(i),i=1,17),math,mfh,mth,nsc
      nsc=nsc+1
   endif
   return
   end subroutine tpidio

   subroutine hdatio(nin,nout,nscr,a,nb,nw)
   !--------------------------------------------------------------------
   ! Utility routine for ENDF bcd and binary tapes.  Read, write,
   ! and/or convert the Hollerith descriptive data to/from a.
   ! If any unit is zero, it is not used.  Positive units are bcd,
   ! and negative units are binary.
   !--------------------------------------------------------------------
   ! externals
   integer::nin,nout,nscr,nb,nw
   real(kr)::a(*)
   ! internals
   integer::inin,inout,inscr
   integer::i,j,k,iend
   integer::ncds,nbx,nwx

   ihol=1

   !--input
   if (nin.lt.0) then
      inin=iabs(nin)
      read(inin) math,mfh,mth,nb,nw,(a(i),i=1,nw)
      c1h=a(1)
      c2h=a(2)
      l1h=nint(a(3))
      l2h=nint(a(4))
      n1h=nint(a(5))/17
      n2h=nint(a(6))
      a(3)=l1h
      a(4)=l2h
      a(5)=n1h*17
      a(6)=n2h
   else if (nin.gt.0) then
      call contio(nin,0,0,a,nb,nw)
      a(5)=a(5)*17
      nb=6+nint(a(5))
      iend=nint(a(5))
      if (iend.gt.npage) iend=npage
      k=6
      do i=1,iend,17
         call hollio(nin,0,0)
         do j=1,17
            k=k+1
            a(k)=rhol(j)
         enddo
      enddo
      nw=k
   endif

   !--output
   n1h=nint(a(5))/17
   if (nin.eq.0.and.a(5).le.npage) nw=6+nint(a(5))
   if (nin.eq.0.and.a(5).gt.npage) nw=6+npage
   nb=6+nint(a(5))-nw
   inout=iabs(nout)
   inscr=iabs(nscr)
   a(5)=n1h*17
   if (nout.lt.0) then
      write(inout) math,mfh,mth,nb,nw,(a(i),i=1,nw)
      inout=0
   endif
   if (nscr.lt.0) then
      write(inscr) math,mfh,mth,nb,nw,(a(i),i=1,nw)
      inscr=0
   endif
   if (nout.le.0.and.nscr.le.0) return
   a(5)=n1h
   ncds=(nw-6)/17
   call contio(0,inout,inscr,a,nbx,nwx)
   j=6
   do i=1,ncds
      do k=1,17
         j=j+1
         rhol(k)=a(j)
      enddo
      call hollio(0,inout,inscr)
   enddo
   return
   end subroutine hdatio

   subroutine dictio(nin,nout,nscr,a,nb,nw)
   !--------------------------------------------------------------------
   ! Utility routine for ENDF bcd and binary tapes.  Read, write,
   ! and/or convert the material dictionary to/from a.  Any unit
   ! which is zero is not used.  Positive units are bcd and
   ! negative units are binary.  On entry, nw is the number of records
   ! in the dictionary.
   !--------------------------------------------------------------------
   ! externals
   integer::nin,nout,nscr,nb,nw
   real(kr)::a(*)
   ! internals
   integer::inin,inout,inscr
   integer::nx,i,j,k,n
   nx=nw

   !--input
   if (nin.lt.0) then
      inin=iabs(nin)
      k=0
      do j=1,nx
         read(inin) math,mfh,mth,nb,nw,(a(i+k),i=1,6)
         k=k+6
      enddo
      nw=k
   else if (nin.gt.0) then
      k=0
      do j=1,nx
         read(nin,'(6e11.0,i4,i2,i3,i5)')&
           (a(i+k),i=1,6),math,mfh,mth,nsp
         k=k+6
      enddo
      nw=k
   endif

   !--output
   nb=0
   if (nout.eq.0.and.nscr.eq.0) return
   inout=iabs(nout)
   inscr=iabs(nscr)
   n=6
   k=0
   do j=1,nx
      l1h=nint(a(k+3))
      l2h=nint(a(k+4))
      n1h=nint(a(k+5))
      n2h=nint(a(k+6))
      a(k+3)=l1h
      a(k+4)=l2h
      a(k+5)=n1h
      a(k+6)=n2h
      if (nout.lt.0) then
         write(inout) math,mfh,mth,nb,n,(a(i+k),i=1,n)
      endif
      if (nout.gt.0) then
         write(nout,'(22x,4i11,i4,i2,i3,i5)')&
           l1h,l2h,n1h,n2h,math,mfh,mth,nsh
         nsh=nsh+1
         if (nsh.gt.99999) nsh=1
      endif
      if (nscr.lt.0) then
         write(inscr) math,mfh,mth,nb,n,(a(i+k),i=1,n)
      endif
      if (nscr.gt.0) then
         write(nscr,'(22x,4i11,i4,i2,i3,i5)')&
           l1h,l2h,n1h,n2h,math,mfh,mth,nsc
         nsc=nsc+1
         if (nsc.gt.99999) nsc=1
      endif
      k=k+6
   enddo
   nw=k
   return
   end subroutine dictio

   subroutine intgio(nin,nout,nscr,a,nb,nw)
   !--------------------------------------------------------------------
   ! utility routine for endf/b coded and blocked binary tapes.
   ! read, write, and/or convert one intg record.
   ! positive units are coded, negative ones are blocked binary.
   ! if any unit is zero, it is not used. Parameter nw determines
   ! the number of entries read or written, as well as the format
   ! for formatted files; nw is an input quantity.  If 0<nw<7 then
   ! nw is interpreted as ndigit and the corresponding array length
   ! is set internally by redefining nw. Parameter nb is always zero.
   !--------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::nin,nout,nscr,nb,nw
   real(kr)::a(*)
   ! internals
   integer::inin,inout,inscr
   integer::j
   integer::ia(20)
   integer,dimension(6),parameter::nwdig=(/0,20,15,13,11,10/)

   !--check number of digits
   nb=0
   if (nw.gt.0.and.nw.le.6) nw=nwdig(nw)

   !--input
   if (nin.lt.0) then
      inin=iabs(nin)
      read(inin) math,mfh,mth,nb,nw,(a(j),j=1,nw)
   else if (nin.gt.0) then
      if (nw.eq.20) then
        read(nin,'(2i5,1x,18i3,1x,i4,i2,i3,i5)')&
          (ia(j),j=1,nw),math,mfh,mth,nsp
      else if (nw.eq.15) then
        read(nin,'(2i5,1x,13i4,3x,i4,i2,i3,i5)')&
          (ia(j),j=1,nw),math,mfh,mth,nsp
      else if (nw.eq.13) then
        read(nin,'(2i5,1x,11i5,i4,i2,i3,i5)')&
          (ia(j),j=1,nw),math,mfh,mth,nsp
      else if (nw.eq.11) then
        read(nin,'(2i5,1x,9i6,1x,i4,i2,i3,i5)')&
          (ia(j),j=1,nw),math,mfh,mth,nsp
      else if (nw.eq.10) then
        read(nin,'(2i5,1x,8i7,i4,i2,i3,i5)')&
          (ia(j),j=1,nw),math,mfh,mth,nsp
      else
        call error('intgio','1invalid nw',' ')
      end if
      do j=1,nw
         a(j)=ia(j)
      enddo
   endif

   !--output
   if (nout.eq.0.and.nscr.eq.0) return
   inout=iabs(nout)
   if (nout.lt.0) then
      write(inout) math,mfh,mth,nb,nw,(a(j),j=1,nw)
      inout=0
   endif
   inscr=iabs(nscr)
   if (nscr.lt.0) then
      write(inscr) math,mfh,mth,nb,nw,(a(j),j=1,nw)
      inscr=0
   endif

   !--format the output
   do j=1,nw
     ia(j)=nint(a(j))
   end do
   if (nscr.gt.0) then
      if (nw.eq.20) then
        write(nscr,'(2i5,1x,18i3,1x,i4,i2,i3,i5)')&
          (ia(j),j=1,nw),math,mfh,mth,nsp
      else if (nw.eq.15) then
        write(nscr,'(2i5,1x,13i4,3x,i4,i2,i3,i5)')&
          (ia(j),j=1,nw),math,mfh,mth,nsp
      else if (nw.eq.13) then
        write(nscr,'(2i5,1x,11i5,i4,i2,i3,i5)')&
          (ia(j),j=1,nw),math,mfh,mth,nsp
      else if (nw.eq.11) then
        write(nscr,'(2i5,1x,9i6,1x,i4,i2,i3,i5)')&
          (ia(j),j=1,nw),math,mfh,mth,nsp
      else if (nw.eq.10) then
        write(nscr,'(2i5,1x,8i7,i4,i2,i3,i5)')&
          (ia(j),j=1,nw),math,mfh,mth,nsp
      else
        call error('intgio','2invalid nw',' ')
      end if
      nsc=nsc+1
      if (nsc.gt.99999) nsc=0
   endif
   if (nout.gt.0) then
      if (nw.eq.20) then
        write(nout,'(2i5,1x,18i3,1x,i4,i2,i3,i5)')&
          (ia(j),j=1,nw),math,mfh,mth,nsp
      else if (nw.eq.15) then
        write(nout,'(2i5,1x,13i4,3x,i4,i2,i3,i5)')&
          (ia(j),j=1,nw),math,mfh,mth,nsp
      else if (nw.eq.13) then
        write(nout,'(2i5,1x,11i5,i4,i2,i3,i5)')&
          (ia(j),j=1,nw),math,mfh,mth,nsp
      else if (nw.eq.11) then
        write(nout,'(2i5,1x,9i6,1x,i4,i2,i3,i5)')&
          (ia(j),j=1,nw),math,mfh,mth,nsp
      else if (nw.eq.10) then
        write(nout,'(2i5,1x,8i7,i4,i2,i3,i5)')&
          (ia(j),j=1,nw),math,mfh,mth,nsp
      else
        call error('intgio','3invalid nw',' ')
      end if
      nsh=nsh+1
      if (nsh.gt.99999) nsh=0
   endif
   return
   end subroutine intgio

   subroutine lineio(nin,nout,nscr)
   !--------------------------------------------------------------------
   ! Utility routine for bcd ENDF tapes.
   ! Read and/or write one line of floating point data from nin
   ! to nout and nscr.  If any unit is zero, it is not used.
   !--------------------------------------------------------------------
   ! externals
   integer::nin,nout,nscr
   ! internals
   integer::k
   character(11)::field(6)

   !--input
   if (nin.gt.0) then
      read(nin,'(6e11.0,i4,i2,i3,i5)') ah,math,mfh,mth,nsp
      nah=6
   endif

   !--output
   if (nout.ne.0) then
      do k=1,6
         field(k)='           '
         if (k.le.nah) then
            call a11(ah(k),field(k))
         endif
      enddo
      write(nout,'(6a11,i4,i2,i3,i5)') (field(k),k=1,6),math,mfh,mth,nsh
      nsh=nsh+1
      if (nsh.gt.99999) nsh=1
   endif
   if (nscr.ne.0) then
      do k=1,6
         field(k)='           '
         if (k.le.nah) then
            call a11(ah(k),field(k))
         endif
      enddo
      write(nscr,'(6a11,i4,i2,i3,i5)') (field(k),k=1,6),math,mfh,mth,nsc
      nsc=nsc+1
      if (nsc.gt.99999) nsc=1
   endif
   return
   end subroutine lineio

   subroutine a11(x,hx)
   !--------------------------------------------------------------------
   ! Convert x to ENDF 11-column format as a string in hx.  Normally,
   ! hx has the forms +1.234567+6, +1.23456-38, or -1.2345+308.
   ! If ix=1 and x is between 1.0 and 1.e7, the 9 significant-figure
   ! extended precision forms +1.23456789, +12.3456789, ...,
   ! +123456.789 in f format are allowed as sometimes needed for
   ! resonance shapes.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::x
   character(11)::hx
   ! internals
   character(1)::s
   real(kr)::ff,xx,f
   integer::n
   real(kr),parameter::zero=0
   real(kr),parameter::tenth=0.1e0_kr
   real(kr),parameter::onem=.999999999e0_kr
   real(kr),parameter::one=1
   real(kr),parameter::ten=10
   real(kr),parameter::tmil=.999999999e7_kr
   real(kr),parameter::up=.000000001e0_kr
   real(kr),parameter::top7=9.9999995e0_kr
   real(kr),parameter::top6=9.999995e0_kr
   real(kr),parameter::top5=9.99995e0_kr
   real(kr),parameter::top9=9.999999995e0_kr
   real(kr),parameter::bot9=9.99999995e0_kr

   !--check for zero as a special case
   if (x.ne.zero) go to 100
      f=0
      s='+'
      n=0
      go to 160
   100 if (x.gt.tenth.and.x.lt.tmil) go to 130

   !--normal seven, six, and five sig-fig forms
   xx=x
   n=int(log10(abs(xx)))
   if (abs(xx).lt.one) go to 110
      ff=xx/ten**n
      s='+'
      if (iabs(n).lt.10.and.abs(ff).lt.top7) go to 120
      if (iabs(n).lt.100.and.iabs(n).ge.10.and.&
        abs(ff).lt.top6) go to 120
      if (iabs(n).ge.100.and.abs(ff).lt.top5) go to 120
      ff=ff/ten
      ff=ff+up
      n=n+1
      go to 120
   110 n=1-n
      ff=xx*ten**n
      s='-'
      if (iabs(n).lt.10.and.abs(ff).lt.top7) go to 120
      if (iabs(n).lt.100.and.iabs(n).ge.10.and.&
        abs(ff).lt.top6) go to 120
      if (iabs(n).ge.100.and.abs(ff).lt.top5) go to 120
      ff=ff/ten
      ff=ff+up
      n=n-1
      if (n.gt.0) go to 120
      s='+'
   120 f=ff
   160 if (iabs(n).lt.10) then
         write(hx,'(f9.6,a,i1)') f,s,n
      else if (iabs(n).lt.100) then
         write(hx,'(f8.5,a,i2)') f,s,n
      else
         write(hx,'(f7.4,a,i3)') f,s,n
      endif
      return

   !--extended nine sig-fig form
  130 n=int(log10(abs(x)))
      if (abs(x).lt.onem) go to 140
      f=x/ten**n
      s='+'
      if (iabs(n).lt.10.and.abs(f).lt.top9) go to 150
      if (iabs(n).ge.10.and.abs(f).lt.bot9) go to 150
      f=f/ten
      n=n+1
      go to 150
  140 f=x
      s='-'
      n=-1
  150 continue
   if (n.le.0) write(hx,'(f11.8)') f
   if (n.eq.1) write(hx,'(1p,f11.7)') f
   if (n.eq.2) write(hx,'(2p,f11.6)') f
   if (n.eq.3) write(hx,'(3p,f11.5)') f
   if (n.eq.4) write(hx,'(4p,f11.4)') f
   if (n.eq.5) write(hx,'(5p,f11.3)') f
   if (n.eq.6) write(hx,'(6p,f11.2)') f
   if (n.eq.-1) n=1
   if (f.gt.onem.and.hx(10:11).eq.'00') write(hx,'(f9.6,a,i1)') f,s,n
   if (f.gt.tenth.and.f.lt.onem.and.hx(11:11).eq.'0')&
     write(hx,'(1p,f9.6,a,i1)') f,s,n
   return
   end subroutine a11

   subroutine tablio(nin,nout,nscr)
   !--------------------------------------------------------------------
   ! Utility routine for bcd ENDF tapes.
   ! Read and/or write the interpolation table for a TAB1 or TAB2
   ! record.  If any unit is zero, it is not used.
   !--------------------------------------------------------------------
   ! externals
   integer::nin,nout,nscr
   ! internals
   integer::i,j,ip2,ipend,nterm

   if (nin.ne.0) nr=n1h
   do i=1,nr,3
      ip2=i+2
      ipend=ip2
      if (ip2.gt.nr) ipend=nr
      nterm=ipend-i+1
      if (nin.ne.0) then
         read(nin,'(6i11,i4,i2,i3,i5)')&
           (nbt(j),jnt(j),j=i,ip2),math,mfh,mth,nsp
      endif
      if (nout.ne.0) then
         if (nterm.eq.1) then
            write(nout,'(2i11,44x,i4,i2,i3,i5)')&
              (nbt(j),jnt(j),j=i,ipend),math,mfh,mth,nsh
         else if (nterm.eq.2) then
            write(nout,'(4i11,22x,i4,i2,i3,i5)')&
              (nbt(j),jnt(j),j=i,ipend),math,mfh,mth,nsh
         else if (nterm.eq.3) then
            write(nout,'(6i11,i4,i2,i3,i5)')&
              (nbt(j),jnt(j),j=i,ipend),math,mfh,mth,nsh
         endif
         nsh=nsh+1
         if (nsh.gt.99999) nsh=1
      endif
      if (nscr.ne.0) then
         if (nterm.eq.1) then
            write(nscr,'(2i11,44x,i4,i2,i3,i5)')&
              (nbt(j),jnt(j),j=i,ipend),math,mfh,mth,nsc
         else if (nterm.eq.2) then
            write(nscr,'(4i11,22x,i4,i2,i3,i5)')&
              (nbt(j),jnt(j),j=i,ipend),math,mfh,mth,nsc
         else if (nterm.eq.3) then
            write(nscr,'(6i11,i4,i2,i3,i5)')&
              (nbt(j),jnt(j),j=i,ipend),math,mfh,mth,nsc
         endif
         nsc=nsc+1
         if (nsc.gt.99999) nsc=1
      endif
   enddo
   return
   end subroutine tablio

   subroutine hollio(nin,nout,nscr)
   !--------------------------------------------------------------------
   ! Utility routine for bcd ENDF tapes.
   ! Read and/or write one line of Hollerith information from
   ! nin to nout and nscr.  If any unit is zero, it is not used.
   !--------------------------------------------------------------------
   ! externals
   integer::nin,nout,nscr

   !--input
   if (nin.ne.0) then
      read(nin,'(16a4,a2,i4,i2,i3,i5)') hol,math,mfh,mth,nsp
   endif

   !--output
   if (nout.ne.0) then
      write(nout,'(16a4,a2,i4,i2,i3,i5)') hol,math,mfh,mth,nsh
      if (nsh.gt.99999) nsh=1
      nsh=nsh+1
   endif
   if (nscr.ne.0) then
      write(nscr,'(16a4,a2,i4,i2,i3,i5)') hol,math,mfh,mth,nsc
      if (nsc.gt.99999) nsc=1
      nsc=nsc+1
   endif
   return
   end subroutine hollio

   subroutine tosend(nin,nout,nscr,a)
   !--------------------------------------------------------------------
   ! Utility routine for ENDF bcd or blocked binary tapes.
   ! Skip or copy through the SEND card.
   ! Positive units are bcd, negative ones are blocked binary.
   !--------------------------------------------------------------------
   use util ! provide error
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::inin,inout,inscr
   integer::mat,mf,mt,ns,i,nb,nw
   character(4)::t(17)

   !--mode conversion not allowed
   if (nin.gt.0.and.(nout.lt.0.or.nscr.lt.0))&
     call error('tosend','mode conversion not allowed.',' ')
   if (nin.lt.0.and.(nout.gt.0.or.nscr.gt.0))&
     call error('tosend','mode conversion not allowed.',' ')

   !--bcd/binary input/output until send is read
   mt=1
   do while (mt.gt.0)
      if (nin.gt.0) then
         read(nin,'(16a4,a2,i4,i2,i3,i5)') (t(i),i=1,17),mat,mf,mt,ns
         if (nscr.ne.0) then
            write(nscr,'(16a4,a2,i4,i2,i3,i5)')&
              (t(i),i=1,17),mat,mf,mt,nsc
            nsc=nsc+1
            if (nsc.gt.99999) nsc=1
         endif
         if (nout.ne.0) then
            write(nout,'(16a4,a2,i4,i2,i3,i5)')&
              (t(i),i=1,17),mat,mf,mt,nsh
            nsh=nsh+1
            if (nsh.gt.99999) nsh=1
         endif
      else
         inin=iabs(nin)
         read(inin)mat,mf,mt,nb,nw,(a(i),i=1,nw)
         if (nscr.ne.0) then
            inscr=iabs(nscr)
            write(inscr)mat,mf,mt,nb,nw,(a(i),i=1,nw)
         endif
         if (nout.ne.0) then
            inout=iabs(nout)
            write(inout)mat,mf,mt,nb,nw,(a(i),i=1,nw)
         endif
      endif
   enddo
   return
   end subroutine tosend

   subroutine tofend(nin,nout,nscr,a)
   !--------------------------------------------------------------------
   ! Utility routine for ENDF bcd or blocked binary tapes.
   ! Skip or copy through the FEND card.
   ! Positive units are bcd, negative ones are blocked binary.
   !--------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::inin,inout,inscr
   integer::mat,mf,mt,ns,i,nb,nw
   character(4)::t(17)

   !--mode conversion not allowed
   if (nin.gt.0.and.(nout.lt.0.or.nscr.lt.0))&
     call error('tofend','mode conversion not allowed.',' ')
   if (nin.lt.0.and.(nout.gt.0.or.nscr.gt.0))&
     call error('tofend','mode conversion not allowed.',' ')

   !--bcd/binary input/output until fend is read
   mf=1
   do while (mf.gt.0)
      if (nin.gt.0) then
         read(nin,'(16a4,a2,i4,i2,i3,i5)') (t(i),i=1,17),mat,mf,mt,ns
         if (nscr.ne.0) then
            write(nscr,'(16a4,a2,i4,i2,i3,i5)')&
             (t(i),i=1,17),mat,mf,mt,nsc
            nsc=nsc+1
            if (nsc.gt.99999) nsc=1
         endif
         if (nout.ne.0) then
            write(nout,'(16a4,a2,i4,i2,i3,i5)')&
              (t(i),i=1,17),mat,mf,mt,nsh
            nsh=nsh+1
            if (nsh.gt.99999) nsh=1
         endif
      else
         inin=iabs(nin)
         read(inin) mat,mf,mt,nb,nw,(a(i),i=1,nw)
         if (nscr.ne.0) then
            inscr=iabs(nscr)
            write(inscr) mat,mf,mt,nb,nw,(a(i),i=1,nw)
         endif
         if (nout.ne.0) then
            inout=iabs(nout)
            write(inout) mat,mf,mt,nb,nw,(a(i),i=1,nw)
         endif
      endif
   enddo
   return
   end subroutine tofend

   subroutine tomend(nin,nout,nscr,a)
   !--------------------------------------------------------------------
   ! Utility routine for ENDF bcd or blocked binary tapes.
   ! Skip or copy through the MEND card.
   ! Positive units are bcd, negative ones are blocked binary.
   !--------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::inin,inout,inscr
   integer::mat,mf,mt,ns,i,nb,nw
   character(4)::t(17)

   !--mode conversion not allowed
   if (nin.gt.0.and.(nout.lt.0.or.nscr.lt.0))&
     call error('tomend','mode conversion not allowed.',' ')
   if (nin.lt.0.and.(nout.gt.0.or.nscr.gt.0))&
     call error('tomend','mode conversion not allowed.',' ')

   !--bcd/binary input/output until mend is read
   mat=1
   do while (mat.gt.0)
      if (nin.gt.0) then
         read(nin,'(16a4,a2,i4,i2,i3,i5)') (t(i),i=1,17),mat,mf,mt,ns
         if (nscr.ne.0) then
            write(nscr,'(16a4,a2,i4,i2,i3,i5)')&
              (t(i),i=1,17),mat,mf,mt,nsc
            nsc=nsc+1
            if (nsc.gt.99999) nsc=1
         endif
         if (nout.ne.0) then
            write(nout,'(16a4,a2,i4,i2,i3,i5)')&
              (t(i),i=1,17),mat,mf,mt,nsh
            nsh=nsh+1
            if (nsh.gt.99999) nsh=1
         endif
      else
         inin=iabs(nin)
         read(inin) mat,mf,mt,nb,nw,(a(i),i=1,nw)
         if (nscr.ne.0) then
            inscr=iabs(nscr)
            write(inscr) mat,mf,mt,nb,nw,(a(i),i=1,nw)
         endif
         if (nout.ne.0) then
            inout=iabs(nout)
            write(inout) mat,mf,mt,nb,nw,(a(i),i=1,nw)
         endif
      endif
   enddo
   return
   end subroutine tomend

   subroutine totend(nin,nout,nscr,a)
   !--------------------------------------------------------------------
   ! Utility routine for ENDF bcd or blocked binary tapes.
   ! Skip or copy through the TEND card.
   ! Positive units are bcd, negative ones are blocked binary.
   !--------------------------------------------------------------------
   use util ! provides error
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::inin,inout,inscr
   integer::mat,mf,mt,ns,i,nb,nw
   character(4)::t(17)

   !--mode conversion not allowed
   if (nin.gt.0.and.(nout.lt.0.or.nscr.lt.0))&
     call error('totend','mode conversion not allowed.',' ')
   if (nin.lt.0.and.(nout.gt.0.or.nscr.gt.0))&
     call error('totend','mode conversion not allowed.',' ')

   !--bcd/binary input/output until tend is read
   mat=1
   do while (mat.gt.-1)
      if (nin.gt.0) then
         read(nin,'(16a4,a2,i4,i2,i3,i5)') (t(i),i=1,17),mat,mf,mt,ns
         if (nscr.ne.0) then
            write(nscr,'(16a4,a2,i4,i2,i3,i5)')&
              (t(i),i=1,17),mat,mf,mt,nsc
            nsc=nsc+1
            if (nsc.gt.99999) nsc=1
         endif
         if (nout.ne.0) then
            write(nout,'(16a4,a2,i4,i2,i3,i5)')&
              (t(i),i=1,17),mat,mf,mt,nsh
            nsh=nsh+1
            if (nsh.gt.99999) nsh=1
         endif
      else
         inin=iabs(nin)
         read(inin) mat,mf,mt,nb,nw,(a(i),i=1,nw)
         if (nscr.ne.0) then
            inscr=iabs(nscr)
            write(inscr) mat,mf,mt,nb,nw,(a(i),i=1,nw)
         endif
         if (nout.ne.0) then
            inout=iabs(nout)
            write(inout) mat,mf,mt,nb,nw,(a(i),i=1,nw)
         endif
      endif
   enddo
   return
   end subroutine totend

   subroutine asend(nout,nscr)
   !--------------------------------------------------------------------
   ! Utility routine for ENDF bcd or blocked binary tapes.
   ! Write a SEND card with blanks.
   ! If a unit is zero, it is not used.
   !--------------------------------------------------------------------
   ! externals
   integer::nout,nscr
   ! internals
   integer::inout,inscr
   integer::nb,nw
   real(kr)::z

   nb=0
   nw=6
   z=0
   mth=0
   if (nout.gt.0) then
      nsh=99999
      write(nout,'(66x,i4,i2,i3,i5)') math,mfh,mth,nsh
      nsh=1
   endif
   if (nout.lt.0) then
      inout=iabs(nout)
      write(inout)math,mfh,mth,nb,nw,z,z,z,z,z,z
   endif
   if (nscr.gt.0) then
      nsc=99999
      write(nscr,'(66x,i4,i2,i3,i5)') math,mfh,mth,nsc
      nsc=1
   endif
   if (nscr.lt.0) then
      inscr=iabs(nscr)
      write(inscr)math,mfh,mth,nb,nw,z,z,z,z,z,z
   endif
   return
   end subroutine asend

   subroutine afend(nout,nscr)
   !--------------------------------------------------------------------
   ! Utility routine for ENDF bcd or block binary tapes.
   ! Write a FEND card with blanks.
   ! If a unit is zero, it is not used.
   !--------------------------------------------------------------------
   ! externals
   integer::nout,nscr
   ! internals
   integer::inout,inscr
   integer::nb,nw
   real(kr)::z

   nb=0
   nw=6
   z=0
   mfh=0
   mth=0
   if (nout.gt.0) then
      nsh=0
      write(nout,'(66x,i4,i2,i3,i5)') math,mfh,mth,nsh
      nsh=nsh+1
   endif
   if (nout.lt.0) then
      inout=iabs(nout)
      write(inout) math,mfh,mth,nb,nw,z,z,z,z,z,z
   endif
   if (nscr.gt.0) then
      nsc=0
      write(nscr,'(66x,i4,i2,i3,i5)') math,mfh,mth,nsc
      nsc=nsc+1
   endif
   if (nscr.lt.0) then
      inscr=iabs(nscr)
      write(inscr) math,mfh,mth,nb,nw,z,z,z,z,z,z
   endif
   return
   end subroutine afend

   subroutine amend(nout,nscr)
   !--------------------------------------------------------------------
   ! Utility routine for ENDF bcd or block binary tapes.
   ! Write a MEND card with blanks.
   ! If a unit is zero, it is not used.
   !--------------------------------------------------------------------
   ! externals
   integer::nout,nscr
   ! internals
   integer::inout,inscr
   integer::nb,nw
   real(kr)::z

   nb=0
   nw=6
   z=0
   math=0
   mfh=0
   mth=0
   if (nout.gt.0) then
      nsh=0
      write(nout,'(66x,i4,i2,i3,i5)') math,mfh,mth,nsh
      nsh=nsh+1
   endif
   if (nout.lt.0) then
      inout=iabs(nout)
      write(inout) math,mfh,mth,nb,nw,z,z,z,z,z,z
   endif
   if (nscr.gt.0) then
      nsc=0
      write(nscr,'(66x,i4,i2,i3,i5)') math,mfh,mth,nsc
      nsc=nsc+1
   endif
   if (nscr.lt.0) then
      inscr=iabs(nscr)
      write(inscr) math,mfh,mth,nb,nw,z,z,z,z,z,z
   endif
   return
   end subroutine amend

   subroutine atend(nout,nscr)
   !--------------------------------------------------------------------
   ! Utility routine for ENDF bcd or block binary tapes.
   ! Write a TEND card with blanks.
   ! If a unit is zero, it is not used.
   !--------------------------------------------------------------------
   ! externals
   integer::nout,nscr
   ! internals
   integer::inout,inscr
   integer::nb,nw
   real(kr)::z

   nb=0
   nw=6
   z=0
   math=-1
   mfh=0
   mth=0
   if (nout.gt.0) then
      nsh=0
      write(nout,'(66x,i4,i2,i3,i5)') math,mfh,mth,nsh
      nsh=nsh+1
   endif
   if (nout.lt.0) then
      inout=iabs(nout)
      write(inout) math,mfh,mth,nb,nw,z,z,z,z,z,z
   endif
   if (nscr.gt.0) then
      nsc=0
      write(nscr,'(66x,i4,i2,i3,i5)') math,mfh,mth,nsc
      nsc=nsc+1
   endif
   if (nscr.lt.0) then
      inscr=iabs(nscr)
      write(inscr) math,mfh,mth,nb,nw,z,z,z,z,z,z
   endif
   return
   end subroutine atend

   subroutine skip6(nin,nout,nscr,a,law)
   !--------------------------------------------------------------------
   ! Utility routine for ENDF File 6.
   ! Skip the next subsection in the current section (MT).
   !--------------------------------------------------------------------
   ! externals
   integer::nin,nout,nscr,law
   real(kr)::a(*)
   ! internals
   integer::ne,ie,nmu,imu,nb,nw

   ! negative laws, law 0, 3 and 4 have no law dependent structure so no
   ! need to skip over it
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
         call tab2io(nin,nout,nscr,a(1),nb,nw)
         nmu=n2h
         do imu=1,nmu
            call tab1io(nin,nout,nscr,a(1),nb,nw)
            do while (nb.ne.0)
               call moreio(nin,nout,nscr,a(1),nb,nw)
            enddo
         enddo
      enddo
   endif
   return
   end subroutine skip6

   subroutine findf(mat,mf,mt,ntape)
   !--------------------------------------------------------------------
   ! Find specified section on an ENDF format tape.
   ! If ntape lt 0, it is assumed to be in binary mode.
   ! If mt=0, find the first section in MAT,MF.
   ! If mf=0, find the first file in MAT.
   ! Routine searches up or down and leaves tape positioned
   ! to read the *head* card for the section requested.
   !--------------------------------------------------------------------
   use util ! provides skiprz,error
   ! externals
   integer::mat,mf,mt,ntape
   ! internals
   integer::math,mfh,mth,matl,mfl,itape,igo
   character strng*60

   !---read first card and determine whether to read up or down.
   itape=iabs(ntape)
   100 continue
   if (ntape.gt.0) then
      read(itape,'(66x,i4,i2,i3)') math,mfh,mth
   else
      read(itape)math,mfh,mth
   endif

   !--test for mat
   if (math.eq.0) go to 100
   if (math.eq.mat) go to 120
   igo=+1
   if (math.eq.-1) igo=-1
   if (math.gt.mat) igo=-1
   go to 200

   !--test for file
   120 continue
   if (mf.eq.0) go to 300
   if (mfh.eq.mf) go to 130
   if (mfh.eq.0) go to 100
   igo=+1
   if (mfh.gt.mf) igo=-1
   go to 200

   !--test for section
   130 continue
   if (mt.eq.0) go to 300
   if (mth.eq.mt) go to 300
   if (mth.eq.0) call skiprz(ntape,-2)
   if (mth.eq.0) go to 100
   igo=+1
   if (mth.gt.mt) igo=-1
   go to 200

   !--search up or down for section
   200 continue
   if (igo.lt.0) call skiprz(ntape,-2)
   if (igo.lt.0) matl=math
   if (igo.lt.0) mfl=mfh
   if (ntape.gt.0) then
      read(itape,'(66x,i4,i2,i3)') math,mfh,mth
   else
      read(itape)math,mfh,mth
   endif

   !--test for mat
   if (math.eq.mat) go to 230
   if (igo.lt.0) go to 220
   if (math.gt.mat.or.math.eq.-1) go to 400
   go to 200
   220 continue
   if (math.lt.mat.and.math.ne.0) go to 400
   if (math.eq.matl.and.mfl.eq.0.and.mfh.eq.0) igo=+1
   go to 200

   !--test for file
   230 continue
   if (mf.eq.0) go to 300
   if (mfh.eq.mf) go to 240
   if (mfh.eq.0) go to 200
   if (igo.lt.0.and.mfh.lt.mf) go to 400
   if (igo.gt.0.and.mfh.gt.mf) go to 400
   go to 200

   !--test for section
   240 continue
   if (mt.eq.0) go to 300
   if (mth.eq.mt) go to 300
   if (mth.eq.0) go to 200
   go to 200

   !--desired section has been found
   !--backspace to card before head card
   300 continue
   call skiprz(ntape,-2)
   if (ntape.gt.0) then
      read(itape,'(66x,i4,i2,i3)') math,mfh,mth
   else
      read(itape)math,mfh,mth
   endif
   if (mth.gt.0) go to 300
   if (mt.eq.0.and.mfh.gt.0) go to 300
   if (mf.eq.0.and.math.eq.mat) go to 300
   return

   !--desired section is not on itape.  write message.
   400 continue
   write(strng,'('' mat'',i4,'' mf'',i2,'' mt'',i3,&
     &'' not on tape '',i3)') mat,mf,mt,itape
   call error('findf',strng,' ')
   return
   end subroutine findf

   subroutine terp1(x1,y1,x2,y2,x,y,i)
   !--------------------------------------------------------------------
   ! Interpolate one point, where
   ! (x1,y1) and (x2,y2) are the end points of the line,
   ! (x,y) is the interpolated point,
   ! i is the interpolation code, and
   ! thr6 is the kinematic threshold for i=6 (thr6.ge.0.).
   !--------------------------------------------------------------------
   ! externals
   real(kr)::x1,y1,x2,y2,x,y
   integer::i
   ! internals
   real(kr)::a,b,t
   real(kr),parameter::zero=0

   !--make sure x2 .ne. x1
   if (x2.eq.x1) then
      y=y1
      return
   endif

   !--y is constant
   if (i.eq.1.or.y2.eq.y1.or.x.eq.x1) then
      y=y1

   !--y is linear in x
   else if (i.eq.2) then
      y=y1+(x-x1)*(y2-y1)/(x2-x1)

   !--y is linear in ln(x)
   else if (i.eq.3) then
      y=y1+log(x/x1)*(y2-y1)/log(x2/x1)

   !--ln(y) is linear in x
   else if (i.eq.4) then
      y=y1*exp((x-x1)*log(y2/y1)/(x2-x1))

   !--ln(y) is linear in ln(x)
   else if (i.eq.5) then
      if (y1.eq.zero) then
         y=y1
      else
         y=y1*exp(log(x/x1)*log(y2/y1)/log(x2/x1))
      endif

   !--coulomb penetrability law (charged particles only)
   else if (i.eq.6) then
      if (y1.eq.zero) then
         y=y1
      else
         t=sqrt(x1-thr6)
         b=log((x2*y2)/(x1*y1))
         b=b/(1/t-1/sqrt(x2-thr6))
         a=exp(b/t)*x1*y1
         y=(a/x)*exp(-b/sqrt(x-thr6))
      endif
   endif
   return
   end subroutine terp1

   subroutine intega(f,x1,x2,a,ip,ir)
   !--------------------------------------------------------------------
   ! Integrate from x1 to x2 in the TAB1 record packed in a.
   ! Assume function is zero outside the range of the table.
   ! The contribution of each panel in the data is computed
   ! analytically in function gral.  On entry, ip and ir are
   ! starting estimates for the first data point greater than
   ! x1 and the corresponding interpolation range.  Initialize
   ! them to 2 and 1 before first call.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::f,x1,x2
   real(kr)::a(*)
   integer::ip,ir
   ! internals
   integer::nr,np,jr,jp,int,it
   real(kr)::xlo,xhi

   !--initialize integral.  set up limits and pointers
   f=0
   nr=nint(a(5))
   np=nint(a(6))
   if (ip.gt.np) ip=np
   if (ir.gt.nr) ir=nr
   jr=5+2*ir
   jp=5+2*nr+2*ip

   !--locate first data point greater than x1
   110 continue
   if (x1.lt.a(jp)) go to 120
   if (ip.eq.np) go to 170
   ! move up
   jp=jp+2
   ip=ip+1
   it=nint(a(jr))
   if (ip.le.it) go to 110
   jr=jr+2
   ir=ir+1
   go to 110
   120 continue
   if (x1.ge.a(jp-2)) go to 130
   if (ip.gt.2) go to 125
   if (x2.le.a(jp-2)) go to 170
   xlo=a(jp-2)
   go to 150
   125 continue
   ! move down
   jp=jp-2
   ip=ip-1
   if (ir.eq.1) go to 110
   it=nint(a(jr-2))
   if (ip.gt.it) go to 110
   jr=jr-2
   ir=ir-1
   go to 110

   !--accumulate contributions to integral
   130 continue
   xlo=x1
   150 continue
   xhi=a(jp)
   if (xhi.gt.x2) xhi=x2
   int=nint(a(jr+1))
   f=f+gral(a(jp-2),a(jp-1),a(jp),a(jp+1),xlo,xhi,int)
   if (xhi.eq.x2) go to 170
   if (ip.eq.np) go to 170
   xlo=xhi
   jp=jp+2
   ip=ip+1
   it=nint(a(jr))
   if (ip.le.it) go to 150
   jr=jr+2
   ir=ir+1
   go to 150

   !--integral is complete
   170 continue
   return
   end subroutine intega

   subroutine terpa(y,x,xnext,idis,a,ip,ir)
   !--------------------------------------------------------------------
   ! Interpolate for y(x) in the TAB1 record packed in a.  Return zero
   ! if x is outside the range of the table.   Here xnext is the next
   ! data grid point greater than x.  On entry, ip and ir are starting
   ! estimates for the first data point greater than x and the
   ! corresponding interpolation range.  Initialize them to 2 and 1
   ! before first call to routine.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::y,x,xnext
   real(kr)::a(*)
   integer::idis,ip,ir
   integer::nr,np,jr,jp,int,it
   real(kr),parameter::shade=1.00001e0_kr
   real(kr),parameter::xbig=1.e12_kr
   real(kr),parameter::zero=0

   !--set up limits and pointers
   nr=nint(a(5))
   np=nint(a(6))
   if (ir.gt.nr) ir=nr
   if (ip.gt.np) ip=np
   jr=5+2*ir
   jp=5+2*nr+2*ip
   idis=0

   !--locate interpolation interval and law for x
   110 continue
   if (x.lt.a(jp)) go to 120
   if (ip.eq.np) go to 150
   ! move up
   jp=jp+2
   ip=ip+1
   it=nint(a(jr))
   if (ip.le.it) go to 110
   jr=jr+2
   ir=ir+1
   go to 110
   120 continue
   if (x.gt.a(jp-2)) go to 130
   if (x.eq.a(jp-2)) go to 140
   if (ip.eq.2) go to 170
   ! move down
   jp=jp-2
   ip=ip-1
   if (ir.eq.1) go to 110
   it=nint(a(jr-2))
   if (ip.gt.it) go to 110
   jr=jr-2
   ir=ir-1
   go to 110

   !--interpolate for y in this interval
   130 continue
   int=nint(a(jr+1))
   call terp1(a(jp-2),a(jp-1),a(jp),a(jp+1),x,y,int)
   xnext=a(jp)
   if (int.eq.1) idis=1
   if (ip.eq.np) return
   if (a(jp+2).eq.xnext) idis=1
   return
   140 continue
   y=a(jp-1)
   int=nint(a(jr+1))
   xnext=a(jp)
   if (int.eq.1) idis=1
   if (ip.eq.np) return
   if (a(jp+2).eq.xnext) idis=1
   return

   !--special branch for last point and above
   150 continue
   if (x.lt.shade*a(jp)) go to 160
   y=0
   xnext=xbig
   return
   160 continue
   y=a(jp+1)
   xnext=xbig
   if (y.gt.zero) xnext=shade*shade*a(jp)
   return

   !--special branch for x below first point
   170 continue
   y=0
   xnext=a(jp-2)
   idis=1
   return
   end subroutine terpa

   real(kr) function gral(xl,yl,xh,yh,x1,x2,int)
   !--------------------------------------------------------------------
   ! Compute the integral from x1 to x2 of the function described
   ! by the endpoints (xl,yl) and (xh,yh) and the ENDF interpolation \
   ! scheme int.
   !--------------------------------------------------------------------
   use util ! provides error
   ! externals
   real(kr)::xl,yl,xh,yh,x1,x2
   integer::int
   ! internals
   real(kr)::a,b,z
   real(kr),parameter::break=0.1e0_kr
   real(kr),parameter::zero=0
   real(kr),parameter::one=1

   !--trivial case.
   gral=0
   if (x2.eq.x1) return

   !--process according to interpolation scheme.
   if (x2.lt.x1) call error('gral','x2 lt x1.',' ')

   !--y is constant
   if (int.eq.1) then
      gral=(x2-x1)*yl

   !--y is linear in x.
   else if (int.eq.2) then
      b=(yh-yl)/(xh-xl)
      a=yl-b*xl
      gral=(x2-x1)*(a+b*(x2+x1)/2)

   !--y is linear in ln(x)
   else if (int.eq.3) then
      if ((xl.le.zero).or.(xh.le.zero)) then
         ! default to linear-linear
         b=(yh-yl)/(xh-xl)
         a=yl-b*xl
         gral=(x2-x1)*(a+b*(x2+x1)/2)
      else
         b=(yh-yl)/log(xh/xl)
         z=(x2-x1)/x1
         if (abs(z).le.break) then
            gral=(x2-x1)*(yl+b*log(x1/xl))+b*x1*z*z*(1+&
              z*(-one/3+z*(one/6-z/10)))/2
         else
            gral=(x2-x1)*(yl+b*log(x1/xl))+b*x1*(1+&
              (x2/x1)*(log(x2/x1)-1))
         endif
      endif

   !--ln(y) linear in x.
   else if (int.eq.4) then
      if ((yl.lt.zero).or.(yh.lt.zero)) then
         ! default to linear-linear
         b=(yh-yl)/(xh-xl)
         a=yl-b*xl
         gral=(x2-x1)*(a+b*(x2+x1)/2)
      else
         b=log(yh/yl)/(xh-xl)
         a=log(yl)-b*xl
         z=(x2-x1)*b
         if (abs(z).le.break) then
            gral=exp(a+b*x1)*(x2-x1)*(1+z*(one/2+z/6))
         else
            gral=exp(a+b*x1)*(exp(z)-1)/b
         endif
      endif

   !--ln(y) is linear in ln(x).
   else if (int.eq.5) then
      if ((xl.le.zero).or.(xh.le.zero)) then
         ! default to ln(y) linear in x
         if ((yl.lt.zero).or.(yh.lt.zero)) then
            ! default to linear-linear
            b=(yh-yl)/(xh-xl)
            a=yl-b*xl
            gral=(x2-x1)*(a+b*(x2+x1)/2)
         else
            b=log(yh/yl)/(xh-xl)
            a=log(yl)-b*xl
            z=(x2-x1)*b
            if (abs(z).le.break) then
               gral=exp(a+b*x1)*(x2-x1)*(1+z*(one/2+z/6))
            else
               gral=exp(a+b*x1)*(exp(z)-1)/b
            endif
         endif
      else if ((yl.lt.zero).or.(yh.lt.zero)) then
         ! default to y is linear in ln(x)
         if ((xl.le.zero).or.(xh.le.zero)) then
            ! default to linear-linear
            b=(yh-yl)/(xh-xl)
            a=yl-b*xl
            gral=(x2-x1)*(a+b*(x2+x1)/2)
         else
            b=(yh-yl)/log(xh/xl)
            z=(x2-x1)/x1
            if (abs(z).le.break) then
               gral=(x2-x1)*(yl+b*log(x1/xl))+b*x1*z*z*(1+&
                 z*(-one/3+z*(one/6-z/10)))/2
            else
               gral=(x2-x1)*(yl+b*log(x1/xl))+b*x1*(1+&
                 (x2/x1)*(log(x2/x1)-1))
            endif
         endif
      else
         b=log(yh/yl)/log(xh/xl)
         z=(b+1)*log(x2/x1)
         if (abs(z).le.break) then
            gral=yl*x1*((x1/xl)**b)*log(x2/x1)*(1+z*(one/2+z/6))
         else
            gral=yl*x1*((x1/xl)**b)*(((x2/x1)**(b+1))-1)/(b+1)
         endif
      endif
   endif
   return
   end function gral

   subroutine gety1(x,xnext,idis,y1,itape,a)
   !--------------------------------------------------------------------
   ! Retrieve y1(x) from an ENDF TAB1 structure using paged bcd or
   ! blocked binary formats.  Call with x=0 to read in first page
   ! or block of data and initialize pointers.  Routine assumes values
   ! will be called in ascending order.  On return, xnext is the
   ! first data grid point greater than x unless x is the last point.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::x,xnext,y1
   real(kr)::a(*)
   integer::idis,itape
   ! internals
   integer::nw,loc,ln,int
   integer::nwtot,nr,np,lt,ip1,ip2,ip,ir,nb,nbt
   real(kr)::xlast,ylast,xn
   save nwtot,nr,np,lt,ip1,ip2,ip,ir,nb,nbt,xlast,ylast
   real(kr),parameter::down=.999999e0_kr
   real(kr),parameter::xbig=1.e12_kr
   real(kr),parameter::zero=0

   !--read first page or block of data and initialize
   idis=0
   if (x.gt.zero) go to 110
   call tab1io(itape,0,0,a,nb,nw)
   nwtot=nw
   nr=nint(a(5))
   np=nint(a(6))
   lt=6+2*nr
   ip1=1
   ip2=(nw-lt)/2
   if (nb.eq.0) ip2=ip2+2
   ir=1
   ip=1
   ! check for zero extension as in mf13.
   100 continue
   loc=lt+2*(ip-ip1)+1
   xlast=a(loc)
   ylast=a(loc+1)
   xnext=xlast
   y1=ylast
   if (ylast.ne.zero) xlast=down*xlast
   if (ylast.ne.zero) return
   if (ip.ge.np-1) return
   if (a(loc+3).ne.zero.and.ip.gt.1) xnext=down*xnext
   if (a(loc+3).ne.zero) return
   ip=ip+1
   nbt=nint(a(5+2*ir))
   if (ip.gt.nbt) ir=ir+1
   if ((ip+2).le.ip2) go to 100
   a(lt+1)=a(nwtot-3)
   a(lt+2)=a(nwtot-2)
   a(lt+3)=a(nwtot-1)
   a(lt+4)=a(nwtot)
   call moreio(itape,0,0,a(lt+5),nb,nw)
   nwtot=nw+lt+4
   ip1=ip
   ip2=ip1+nw/2+1
   if (nb.eq.0) ip2=ip2+2
   go to 100

   !--is x in this panel
   110 continue
   if (x.lt.xlast) go to 200
   ln=2*(ip-ip1)+lt
   if (x.lt.a(ln+1)) go to 120
   if (ip.eq.np) go to 300

   !--no.  move up to next range.
   !--read in new page of data if needed.
   xlast=a(ln+1)
   ylast=a(ln+2)
   ip=ip+1
   nbt=nint(a(5+2*ir))
   if (ip.gt.nbt) ir=ir+1
   if ((ip+2).le.ip2) go to 110
   if (nb.eq.0) go to 120
   a(lt+1)=a(nwtot-3)
   a(lt+2)=a(nwtot-2)
   a(lt+3)=a(nwtot-1)
   a(lt+4)=a(nwtot)
   call moreio(itape,0,0,a(lt+5),nb,nw)
   nwtot=nw+lt+4
   ip1=ip
   ip2=ip1+nw/2+1
   if (nb.eq.0) ip2=ip2+2
   go to 110

   !--yes.  interpolate for desired value
   120 continue
   int=nint(a(2*ir+6))
   if (int.eq.1) idis=1
   call terp1(xlast,ylast,a(ln+1),a(ln+2),x,y1,int)
   xnext=a(ln+1)
   if ((ln+3).gt.nwtot.and.nb.eq.0) return
   xn=a(ln+3)
   if (xn.eq.xnext) idis=1
   return

   !--special branch for x outside range of table
   200 continue
   y1=0
   xnext=xlast
   return

   !--special branch for last point
   300 continue
   y1=a(ln+2)
   xlast=xbig
   xnext=xbig
   return
   end subroutine gety1

   subroutine gety2(x,xnext,idis,y2,itape,a)
   !--------------------------------------------------------------------
   ! Retrieve y2(x) from an ENDF TAB1 structure using paged bcd or
   ! blocked binary formats.  Call with x=0 to read in first page
   ! or block of data and initialize pointers.  Routine assumes
   ! values will be called in ascending order.  On return, xnext is
   ! the first data grid point greater than x unless x is the last
   ! point.  This subroutine is similar to gety1, for use when
   ! two y quantities must be retrieved simultaneously.
   !--------------------------------------------------------------------
   ! externals
   real(kr)::x,xnext,y2
   real(kr)::a(*)
   integer::idis,itape
   ! internals
   integer::nw,ln,int
   integer::nwtot,nr,np,lt,ip1,ip2,ip,ir,nb,nbt
   real(kr)::xlast,ylast,xn
   save nwtot,nr,np,lt,ip1,ip2,ip,ir,nb,nbt,xlast,ylast
   real(kr),parameter::down=.999999e0_kr
   real(kr),parameter::zero=0
   real(kr),parameter::xbig=1.e12_kr

   !--read first page or block of data and initialize
   idis=0
   if (x.gt.zero) go to 110
   call tab1io(itape,0,0,a,nb,nw)
   nwtot=nw
   nr=nint(a(5))
   np=nint(a(6))
   lt=6+2*nr
   ip1=1
   ip2=(nw-lt)/2
   if (nb.eq.0) ip2=ip2+2
   ir=1
   ip=1
   xlast=a(lt+1)
   ylast=a(lt+2)
   xnext=xlast
   y2=ylast
   if (ylast.ne.zero) xlast=down*xlast
   return

   !--is x in this panel
   110 continue
   if (x.lt.xlast) go to 200
   ln=2*(ip-ip1)+lt
   if (x.lt.a(ln+1)) go to 120
   if (ip.eq.np) go to 300

   !--no.  move up to next range.
   !--read in new page of data if needed.
   xlast=a(ln+1)
   ylast=a(ln+2)
   ip=ip+1
   nbt=nint(a(5+2*ir))
   if (ip.gt.nbt) ir=ir+1
   if ((ip+2).le.ip2) go to 110
   if (nb.eq.0) go to 120
   a(lt+1)=a(nwtot-3)
   a(lt+2)=a(nwtot-2)
   a(lt+3)=a(nwtot-1)
   a(lt+4)=a(nwtot)
   call moreio(itape,0,0,a(lt+5),nb,nw)
   nwtot=nw+lt+4
   ip1=ip
   ip2=ip1+nw/2+1
   if (nb.eq.0) ip2=ip2+2
   go to 110

   !--yes.  interpolate for desired value
   120 continue
   int=nint(a(2*ir+6))
   if (int.eq.1) idis=1
   call terp1(xlast,ylast,a(ln+1),a(ln+2),x,y2,int)
   xnext=a(ln+1)
   if ((ln+3).gt.nwtot.and.nb.eq.0) return
   xn=a(ln+3)
   if (xn.eq.xnext) idis=1
   return

   !--special branch for x outside range of table
   200 continue
   y2=0
   xnext=xlast
   return

   !--special branch for last point
   300 continue
   y2=a(ln+2)
   xlast=xbig
   xnext=xbig
   return
   end subroutine gety2

end module endf
