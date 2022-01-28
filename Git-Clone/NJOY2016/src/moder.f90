module modem
   ! provides moder for NJOY2016
   use locale
   implicit none
   private
   public moder

   ! global variables
   integer::ntw,ng

contains

   subroutine moder
   !-------------------------------------------------------------------
   !
   !  Change the mode of an ENDF tape.
   !
   !  Converts ENDF-type files between binary and text formats.
   !  Also works for PENDF, GENDF and covariance tapes.
   !  Moder can also be used to select materials from an ENDF-type
   !  tape, or to merge several materials into a new tape.
   !
   !---input specifications (free format)---------------------------
   !
   ! card 1       unit numbers
   !      nin     input unit
   !      nout    output unit
   !
   !  a positive unit is coded (mode 3).
   !  a negative unit is blocked binary (njoy mode).
   !
   !  note: abs(nin) ge 1 and le 19 is a flag to select various
   !        materials from one or more input tapes, with or
   !        without mode conversion.  the kind of data to be
   !        processed is keyed to nin as follows:
   !             nin=1, for endf or pendf input and output,
   !                 2, for gendf input and output,
   !                 3, for errorr-format input and output.
   !
   !      cards 2 and 3 for abs (nin) ge 1 and le 19 only.
   !
   ! card 2
   !      tpid    tapeid for nout.  66 characters allowed
   !              (delimited with *, ended with /)
   ! card 3
   !      nin     input unit
   !              terminate moder by setting nin=0
   !      matd    material on this tape to add to nout
   !-------------------------------------------------------------------
   use mainio ! provides nysi,nsyso
   use util   ! provides timer,openz,repoz,error,mess,skiprz,closz
   use endf   ! provides endf routines and variables
   ! internals
   integer::nin,nout,no,loop,i,nz,inout
   integer::matd,mend,ig,nk,ik,nb,nw
   integer::nscr,ninl,matl
   real(kr)::time
   character(105)::strng
   character(4)::hb(17)
   real(kr)::rb(17)
   equivalence(hb(1),rb(1))
   real(kr)::a(5200)
   real(kr)::z(20)
   character(4)::zc(20)
   equivalence(zc(1),z(1))

   nscr=0
   ninl=0
   matl=-2

   !--read user input, initialize,
   !--and write output header.
   call timer(time)
   write(nsyso,'(/&
     &'' moder...change the mode of an endf tape or njoy '',&
     &''output tape'',9x,f8.1,''s'')') time
   write(nsyse,'(/'' moder...'',60x,f8.1,''s'')') time
   read(nsysi,*) nin,nout
   call openz(nout,1)
   if (iabs(nin).ge.20) call openz(nin,0)
   nsh=0
   nsp=0
   nsc=0
   call repoz(nout)
   no=nout
   loop=0
   if (iabs(nin).ge.1.and.iabs(nin).le.19) loop=1
   if (loop.eq.1) go to 110
   write(nsyso,'(/&
     &'' input unit (+ for coded, - for bb) ... '',i10/&
     &'' output unit (+ for coded, - for bb) .. '',i10)')&
     nin,nout
   go to 140

  110 continue
   write(nsyso,'(/&
     &'' put materials from various tapes on output tape '',i3)') nout
   strng=' '
   read(nsysi,*) strng
   read(strng,'(16a4,a2)') (zc(i),i=1,17)
   nz=17
   write(nsyso,'(/&
     &'' tape id for nout''/'' ----------------''/1x,16a4,a2)')&
     (zc(i),i=1,nz)
   math=1
   mfh=0
   mth=0
   inout=1
   if (nin.eq.2) inout=2
   if (nin.eq.3) inout=3
   if (inout.eq.1) then
      write(nsyso,'(/'' processing endf or pendf tape.'')')
      write(nsyse,'(/'' processing endf or pendf tape.'')')
   else if (inout.eq.2) then
      write(nsyso,'(/'' processing gendf tape'')')
      write(nsyse,'(/'' processing gendf tape'')')
   else if (inout.eq.3) then
      write(nsyso,'(/'' processing covariance tape'')')
      write(nsyse,'(/'' processing covariance tape'')')
   endif
   call tpidio(0,nout,nscr,z,nb,nz)
   no=0
   write(nsyso,'(/10x,''nin'',10x,''matd''/&
     &10x,''---'',10x,''----'')')

  130 continue
   nz=2
   nin=0
   read(nsysi,*) nin,matd
   if (nin.eq.0) go to 1000
   write(nsyso,'(10x,i3,10x,i4)') nin,matd
   if (inout.eq.1.and.matd.le.matl) call error('moder',&
     'endf materials must be in ascending order.',' ')
   if (nin.ne.ninl) call closz(ninl)
   if (nin.ne.ninl) call openz(nin,0)
   if (nin.eq.ninl.and.inout.ge.2) go to 140
   if (nin.eq.ninl.and.matd.gt.matl) go to 180

  140 continue
   call repoz(nin)
   call repoz(nscr)
   call tpidio(nin,no,nscr,a,nb,nw)
   if (loop.gt.0) go to 180
   do i=1,17
      rb(i)=a(i)
   enddo
   write(strng,'(16a4,a2)') (hb(i),i=1,17)
   if (strng.ne.' ') then
      write(nsyso,&
        '(/'' tape label''/1x,9(''----''),''--''/2x,17a4)')&
        (hb(i),i=1,17)
   else
      write(nsyso,'('' tape id is blank'')')
   endif
   go to 205
  180 continue
   call contio(nin,0,0,a,nb,nw)
   if (n1h.eq.-1) go to 190
   call findf(matd,0,0,nin)
   ninl=nin
   matl=matd
   go to 205
  190 continue
   if (math.ne.matd) go to 200
   ninl=nin
   matl=matd
   go to 700
  200 continue
   call tomend(nin,0,0,a)
   call contio(nin,0,0,a,nb,nw)
   if (math.ne.-1) go to 190
   write(strng,'(''mat='',i4,'' not found on gendf tape'')') matd
   call mess('moder',strng,' ')
   go to 130
  205 continue
   iverf=0

   !--process endf or pendf files
  210 continue
   call contio(nin,0,0,a,nb,nw)
   if (mfh.ne.1.or.mth.ne.451) go to 215
   nsh=1
   if (n1h.eq.-1) go to 700
   if (iverf.ne.0) go to 215
   call contio(nin,0,0,a(7),nb,nw)
   if (n1h.ne.0) then
      iverf=4
   else if (n2h.eq.0) then
      iverf=5
   else
      iverf=6
   endif
   write(nsyso,'(/'' using endf-'',i1,'' format'')') iverf
   call skiprz(nin,-1)
  215 continue
   if (loop.eq.0.or.math.eq.0) go to 220
   if (math.ne.matd) go to 130
  220 continue
   call contio(0,nout,nscr,a,nb,nw)
   if (loop.gt.0.and.math.eq.0) go to 130
   if (math.eq.-1) go to 1010
   if (mfh.gt.0) go to 230
   if (math.gt.0) go to 210
   nsh=1
   nsp=1
   nsc=1
   go to 210
  230 continue
   if (mfh.eq.1) go to 270
   if (mfh.eq.2) then
      call file2(nin,nout,nscr,a)
   else if (mfh.eq.3) then
      call file3(nin,nout,nscr,a)
   else if (mfh.eq.4) then
      call file4(nin,nout,nscr,a)
   else if (mfh.eq.5) then
      call file5(nin,nout,nscr,a)
   else if (mfh.eq.6) then
      call file6(nin,nout,nscr,a)
   else if (mfh.eq.7) then
      call file7(nin,nout,nscr,a)
   else if (mfh.eq.8) then
      call file8(nin,nout,nscr,a)
   else if (mfh.eq.9.or.mfh.eq.10) then
      call file9(nin,nout,nscr,a)
   else if (mfh.eq.12) then
      call file1x(nin,nout,nscr,a)
   else if (mfh.eq.13) then
      call file1x(nin,nout,nscr,a)
   else if (mfh.eq.14) then
      call file14(nin,nout,nscr,a)
   else if (mfh.eq.15) then
      call file15(nin,nout,nscr,a)
   else if (mfh.eq.23) then
      call file3(nin,nout,nscr,a)
   else if (mfh.eq.24) then
      call file4(nin,nout,nscr,a)
   else if (mfh.eq.25) then
      call file5(nin,nout,nscr,a)
   else if (mfh.eq.26) then
      call file6(nin,nout,nscr,a)
   else if (mfh.eq.27) then
      call file3(nin,nout,nscr,a)
   else if (mfh.eq.30) then
      call file1x(nin,nout,nscr,a)
   else if (mfh.eq.31) then
      call file3x(nin,nout,nscr,a)
   else if (mfh.eq.32) then
      call file32(nin,nout,nscr,a)
   else if (mfh.eq.33) then
      call file3x(nin,nout,nscr,a)
   else if (mfh.eq.34) then
      call file34(nin,nout,nscr,a)
   else if (mfh.eq.35) then
      call file35(nin,nout,nscr,a)
   else if (mfh.eq.40) then
      call file40(nin,nout,nscr,a)
   else
      write(strng,'(''conversion not coded for mf='',i3)') mfh
      call error('moder',strng,' ')
   endif
   go to 410
  270 continue
   if (mth.eq.451.and.n1h.eq.-11) go to 800
   if (mth.eq.451.and.n1h.eq.-12) go to 800
   if (mth.eq.451.and.n1h.eq.-14) go to 800
   if (loop.ne.0) then
      if (inout.eq.2)&
        call error('moder','this material is not a gendf material.',' ')
      if (inout.eq.3)&
        call error('moder','input is not an errorr output tape.',' ')
   endif
   call file1(nin,nout,nscr,a)
   go to 410
  410 continue
   call contio(nin,nout,nscr,a,nb,nw)
   if (mth.ne.0) then
      write(strng,'(1x,5i6,1p,6e12.5)') math,mfh,mth,nb,nw,(a(i),i=1,6)
      call error('moder','should have found send card',strng)
   endif
   go to 210

   !--process gendf tapes
  700 continue
   if (loop.ne.0) then
      if (inout.eq.1)&
        call error('moder','input is not an endf or pendf tape.',' ')
      if (inout.eq.3)&
        call error('moder','input is not an errorr output tape.',' ')
   endif
   ntw=n2h
  705 continue
   if (n1h.ne.-1) then
      write(strng,'(''mat='',i4,'' is not a gendf material'')') math
      call error('moder',strng,' ')
   endif
   call glstio(nin,nout,nscr,a,nb,nw)
   mend=0
   go to 740
  710 continue
   call listio(nin,nout,nscr,a,nb,nw)
   ig=n2h
   do while (nb.ne.0)
      call moreio(nin,nout,nscr,a,nb,nw)
   enddo
   if (mfh.ne.1.and.ig.ne.ng) go to 710
  740 continue
   call contio(nin,0,0,a,nb,nw)
   if (loop.gt.0) then
       if (math.ne.matd.and.math.ne.0) go to 130
   endif
   if (math.eq.-1) go to 1000
   if (math.gt.0) go to 750
   mend=1
   call amend(nout,nscr)
   nsp=1
   nsc=1
   go to 740
  750 continue
   if (mend.gt.0) go to 705
   call contio(0,nout,nscr,a,nb,nw)
   if (math.le.0) go to 740
   if (mfh.eq.0) go to 740
   if (mth.eq.0) go to 740
   go to 710

   !--process covariance tape
  800 continue
   if (loop.ne.0) then
      if (inout.eq.1)&
        call error('moder','input is not an endf tape.',' ')
      if (inout.eq.2)&
        call error('moder','input is not a gendf tape.',' ')
   endif
  805 call listio(nin,nout,nscr,a,nb,nw)
   ng=l1h
   do while (nb.ne.0)
      call moreio(nin,nout,nscr,a,nb,nw)
   enddo
  820 continue
   call contio(nin,nout,nscr,a,nb,nw)
   if (loop.gt.0.and.math.eq.0) go to 130
   if (math.eq.-1) go to 1010
   if (mth.eq.0) go to 820
   if (mfh.eq.1) go to 805

   !--covariance tape mf3
   if (mfh.eq.3.or.mfh.eq.5) then
      ! files 3 & 5 have no head cards.  backspace and read list record
      call skiprz(nin,-1)
      call skiprz(nout,-1)
      if (nout.gt.0) nsh=nsh-1
      call listio(nin,nout,nscr,a,nb,nw)
      do while (nb.ne.0)
         call moreio(nin,nout,nscr,a,nb,nw)
      enddo

   !--covariance tapes mf33 or mf40
   else if (mfh.eq.33.or.mfh.eq.40) then
      nk=n2h
      ! loop over subsections
      do ik=1,nk
         call contio(nin,nout,nscr,a,nb,nw)
         ! loop over groups
         ig=0
         do while (ig.lt.ng)
            call listio(nin,nout,nscr,a,nb,nw)
            ig=n2h
            do while (nb.ne.0)
               call moreio(nin,nout,nscr,a,nb,nw)
            enddo
         enddo
      enddo

   !--covariance tape mf34 or mf35
   else if (mfh.eq.34.or.mfh.eq.35) then
      call contio(nin,nout,nscr,a,nb,nw)
      ng=n2h
      ig=0
      ! loop over groups
      do while (ig.lt.ng)
         call listio(nin,nout,nscr,a,nb,nw)
         ig=n2h
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
      enddo

   !--only allow mf1, mf3 and mf33, 34, 35 & 40 on covariance tapes
   else
      write(strng,'(''illegal covariance mf='',i3)') mfh
      call error('moder',strng,' ')
   endif
   go to 820

   !--moder complete
 1000 continue
   call atend(nout,nscr)
 1010 continue
   call repoz(nout)
   call repoz(nin)
   call closz(nin)
   call closz(nout)
   call timer(time)
   write(nsyso,'(69x,f8.1,''s''/&
     &1x,7(''**********''),''*******'')') time
   return
   end subroutine moder

   subroutine file1(nin,nout,nscr,a)
   !-------------------------------------------------------------------
   ! Convert mode of File 1.
   ! General information.
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::nx,ns,i,npr,j,lep1,lnd,lnp,ns457,nb,nw,ng,ng460,lfc,nfc,lo,nply

   !--hollerith descriptive data and tape dictionary.
   if (mth.eq.451) then
      if (iverf.le.4) nx=n2h
      if (iverf.ge.5) call contio(nin,nout,nscr,a,nb,nw)
      if (iverf.ge.6) call contio(nin,nout,nscr,a,nb,nw)
      call hdatio(nin,nout,nscr,a,nb,nw)
      if (iverf.ge.5) nx=n2h
      do while (nb.ne.0)
         call moreio(nin,nout,nscr,a,nb,nw)
      enddo
      nw=nx
      call dictio(nin,nout,nscr,a,nb,nw)

   !--total fission yield (nu-bar).
   else if (mth.eq.452) then
      if (l2h.ne.2) then
         call listio(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
      else
         call tab1io(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
      endif

   !--induced reaction branching ratios.
   else if (mth.eq.453) then
      ns=n1h
      do i=1,ns
         call listio(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
         npr=n2h
         do j=1,npr
            call listio(nin,nout,nscr,a,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,nout,nscr,a,nb,nw)
            enddo
         enddo
      enddo

   !--fission product yield data.
   else if (mth.eq.454) then
      lep1=l1h
      do i=1,lep1
         call listio(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
      enddo

   !--delayed fission yield (nu-bar).
   else if (mth.eq.455) then
      lnd=l2h
      call listio(nin,nout,nscr,a,nb,nw)
      do while (nb.ne.0)
         call moreio(nin,nout,nscr,a,nb,nw)
      enddo
      if (lnd.ne.2) then
         call listio(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
      else
         call tab1io(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
      endif

   !--prompt fission yield (nu-bar).
   else if (mth.eq.456) then
      lnp=l2h
      if (lnp.ne.2) then
         call listio(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
      else
         call tab1io(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
      endif

   !--radioactive decay data.
   else if (mth.eq.457) then
      ns457=n2h
      do i=1,2
         call listio(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
      enddo
      ns=0
      do while (ns.lt.ns457)
         call listio(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
         ns=ns+1
      enddo

   !--components of energy release due to fission.
   else if (mth.eq.458) then
      lfc=l2h
      nfc=n2h
      call listio(nin,nout,nscr,a,nb,nw)
      nply=l2h
      do while (nb.ne.0)
         call moreio(nin,nout,nscr,a,nb,nw)
      enddo
      if (lfc.eq.1) then
         if (nply.ne.0) then
            call error('file1','bad NPLY in mt=458.',' ')
         endif
         do i=1,nfc
            call tab1io(nin,nout,nscr,a,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,nout,nscr,a,nb,nw)
            enddo
         enddo
      else if (lfc.ne.0) then
         call error('file1','bad LFC in mt=458.',' ')
      endif

   !--beta-delayed photon spectra
   else if (mth.eq.460) then
      lo=l1h
      if (lo.eq.1) then
         ng460=n1h
         do ng=1,ng460
            call tab1io(nin,nout,nscr,a,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,nout,nscr,a,nb,nw)
            enddo
         enddo
      else if (lo.eq.2) then
         call listio(nin,nout,nscr,a,nb,nw)
      else
         call error('file1','bad LO in mt=460.',' ')
      endif

   !--illegal mt for file 1
   else
      call error('file1','illegal mt.',' ')
   endif
   return
   end subroutine file1

   subroutine file2(nin,nout,nscr,a)
   !-------------------------------------------------------------------
   ! Convert mode of File 2.
   ! Resonance parameters.
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)

   !--normal endf resonance parameters
   if (mth.eq.151) then
      call file2a(nin,nout,nscr,a)

   !--special point-unresolved and probability-table formats
   else if (mth.eq.152.or.mth.eq.153) then
      call file2b(nin,nout,nscr,a)

   !--illegal mt for file 2.
   else
      call error('file2','illegal mt.',' ')
   endif
   return
   end subroutine file2

   subroutine file2a(nin,nout,nscr,a)
   !-------------------------------------------------------------------
   ! Convert mode of File 2.
   ! Resonance parameters.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::nis,i,ner,lfw,n,lru,lrf,nro,nls,l,njs,j,nb,nw
   integer::kbk,kps,lbk,lps

   nis=n1h
   do i=1,nis
      call contio(nin,nout,nscr,a,nb,nw)
      ner=n1h
      lfw=l2h
      do n=1,ner
         call contio(nin,nout,nscr,a,nb,nw)
         lru=l1h
         lrf=l2h
         nro=n1h
         if (nro.ne.0) call tab1io(nin,nout,nscr,a,nb,nw)
         if (lru.eq.2.and.lrf.eq.1.and.lfw.eq.1) then
            call listio(nin,nout,nscr,a,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,nout,nscr,a,nb,nw)
            enddo
            nls=n2h
            do l=1,nls
               call contio(nin,nout,nscr,a,nb,nw)
               njs=n1h
               do j=1,njs
                  call listio(nin,nout,nscr,a,nb,nw)
                  do while (nb.ne.0)
                     call moreio(nin,nout,nscr,a,nb,nw)
                  enddo
               enddo
            enddo
         else
            call contio(nin,nout,nscr,a,nb,nw)
            if (lrf.ne.7) then
               nls=n1h
               if (nls.ne.0) then
                  do l=1,nls
                     if (lru.eq.2.and.lrf.eq.2) then
                        call contio(nin,nout,nscr,a,nb,nw)
                        njs=n1h
                        do j=1,njs
                           call listio(nin,nout,nscr,a,nb,nw)
                           do while (nb.ne.0)
                              call moreio(nin,nout,nscr,a,nb,nw)
                           enddo
                        enddo
                     else
                        call listio(nin,nout,nscr,a,nb,nw)
                        do while (nb.ne.0)
                           call moreio(nin,nout,nscr,a,nb,nw)
                        enddo
                        if ((lru.eq.1.and.lrf.eq.4)&
                          .or.(lru.eq.2.and.lfw.ne.0)) then
                           call contio(nin,nout,nscr,a,nb,nw)
                           njs=n1h
                           do j=1,njs
                              call listio(nin,nout,nscr,a,nb,nw)
                              do while (nb.ne.0)
                                 call moreio(nin,nout,nscr,a,nb,nw)
                              enddo
                           enddo
                        endif
                     endif
                  enddo
               endif
            else
               njs=n1h
               call listio(nin,nout,nscr,a,nb,nw)
               do while (nb.ne.0)
                  call moreio(nin,nout,nscr,a,nb,nw)
               enddo
               do j=1,njs
                  call listio(nin,nout,nscr,a,nb,nw)
                  do while (nb.ne.0)
                     call moreio(nin,nout,nscr,a,nb,nw)
                  enddo
                  kbk=l1h
                  kps=l2h
                  call listio(nin,nout,nscr,a,nb,nw)
                  do while (nb.ne.0)
                     call moreio(nin,nout,nscr,a,nb,nw)
                  enddo
                  if (kbk.gt.0) then
                     call listio(nin,nout,nscr,a,nb,nw)
                     lbk=n1h
                     if (lbk.eq.1) then
                        call tab1io(nin,nout,nscr,a,nb,nw)
                        do while (nb.ne.0)
                           call moreio(nin,nout,nscr,a,nb,nw)
                        enddo
                        call tab1io(nin,nout,nscr,a,nb,nw)
                        do while (nb.ne.0)
                           call moreio(nin,nout,nscr,a,nb,nw)
                        enddo
                     endif
                  endif
                  if (kps.eq.1)then
                     call listio(nin,nout,nscr,a,nb,nw)
                     lps=n1h
                     if (lps.eq.1) then
                        call tab1io(nin,nout,nscr,a,nb,nw)
                        do while (nb.ne.0)
                           call moreio(nin,nout,nscr,a,nb,nw)
                        enddo
                     endif
                  endif
               enddo
            endif
         endif
      enddo
   enddo
   return
   end subroutine file2a

   subroutine file2b(nin,nout,nscr,a)
   !-------------------------------------------------------------------
   ! Convert mode of File 2 for
   ! special NJOY point-unresolved and probability-table formats.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::nb,nw

   call  listio(nin,nout,nscr,a,nb,nw)
   do while (nb.ne.0)
      call moreio(nin,nout,nscr,a,nb,nw)
   enddo
   return
   end subroutine file2b

   subroutine file3(nin,nout,nscr,a)
   !-------------------------------------------------------------------
   ! Convert mode of File 3.
   ! Cross sections.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::lt,nb,nw

   call tab1io(nin,nout,nscr,a,nb,nw)
   lt=l1h
   do while (nb.ne.0)
      call moreio(nin,nout,nscr,a,nb,nw)
   enddo
   do while (lt.gt.0)
      call listio(nin,nout,nscr,a,nb,nw)
      do while (nb.ne.0)
         call moreio(nin,nout,nscr,a,nb,nw)
      enddo
      lt=lt-1
   enddo
   end subroutine file3

   subroutine file4(nin,nout,nscr,a)
   !-------------------------------------------------------------------
   ! Convert mode of File 4.
   ! Angular distributions.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::lvt,ltt,nk,ik,ne,ie,lt,nb,nw

   lvt=l1h
   ltt=l2h
   nk=n1h
   if (nk.eq.0) nk=1
   do ik=1,nk
      if (lvt.eq.0) call contio(nin,nout,nscr,a,nb,nw)
      if (lvt.ne.0.or.l1h.ne.1) then
         if (lvt.eq.1) then
            call listio(nin,nout,nscr,a,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,nout,nscr,a,nb,nw)
            enddo
         endif
         call tab2io(nin,nout,nscr,a,nb,nw)
         ne=n2h
         if (ltt.ne.2) then
            do ie=1,ne
               call listio(nin,nout,nscr,a,nb,nw)
               lt=l1h
               do while (nb.ne.0)
                  call moreio(nin,nout,nscr,a,nb,nw)
               enddo
               do while (lt.gt.0)
                  call listio(nin,nout,nscr,a,nb,nw)
                  do while (nb.ne.0)
                     call moreio(nin,nout,nscr,a,nb,nw)
                  enddo
                  lt=lt-1
               enddo
            enddo
         else
            do ie=1,ne
               call tab1io(nin,nout,nscr,a,nb,nw)
               lt=l1h
               do while (nb.ne.0)
                  call moreio(nin,nout,nscr,a,nb,nw)
               enddo
               do while (lt.gt.0)
                  call listio(nin,nout,nscr,a,nb,nw)
                  do while (nb.ne.0)
                     call moreio(nin,nout,nscr,a,nb,nw)
                  enddo
                  lt=lt-1
               enddo
            enddo
         endif
      endif
      if (ltt.ge.3) then
         call tab2io(nin,nout,nscr,a,nb,nw)
         ne=n2h
         do ie=1,ne
            call tab1io(nin,nout,nscr,a,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,nout,nscr,a,nb,nw)
            enddo
         enddo
      endif
   enddo
   return
   end subroutine file4

   subroutine file5(nin,nout,nscr,a)
   !-------------------------------------------------------------------
   ! Convert mode of File 5.
   ! Secondary energy distributions.
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::nk,k,lf,ne,ie,nb,nw

   nk=n1h
   do k=1,nk
      call tab1io(nin,nout,nscr,a,nb,nw)
      do while (nb.ne.0)
         call moreio(nin,nout,nscr,a,nb,nw)
      enddo
      lf=l2h
      if (lf.eq.1) then
         call tab2io(nin,nout,nscr,a,nb,nw)
         ne=n2h
         do ie=1,ne
            call tab1io(nin,nout,nscr,a,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,nout,nscr,a,nb,nw)
            enddo
         enddo
      else if (lf.eq.5.or.lf.eq.11) then
         call tab1io(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
         call tab1io(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
      else if (lf.eq.7.or.lf.eq.9.or.lf.eq.12) then
         call tab1io(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
      else if (lf.eq.10) then
         call listio(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
      else if (lf.ne.3) then
         call error('file5','illegal lf.',' ')
      endif
   enddo
   return
   end subroutine file5

   subroutine file6(nin,nout,nscr,a)
   !-------------------------------------------------------------------
   ! Convert mode of File 6.
   ! Particle energy-angle distributions
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::nk,ik,law,ne,ie,nmu,imu,nb,nw
   character(60)::strng

   nk=n1h
   do ik=1,nk
      call tab1io(nin,nout,nscr,a,nb,nw)
      law=l2h
      do while (nb.ne.0)
         call moreio(nin,nout,nscr,a,nb,nw)
      enddo
      if (law.eq.6) then
         call contio(nin,nout,nscr,a,nb,nw)
      else if (law.eq.1.or.law.eq.2.or.law.eq.5) then
         call tab2io(nin,nout,nscr,a,nb,nw)
         ne=n2h
         do ie=1,ne
            call listio(nin,nout,nscr,a,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,nout,nscr,a,nb,nw)
            enddo
         enddo
      else if (law.eq.7) then
         call tab2io(nin,nout,nscr,a,nb,nw)
         ne=n2h
         do ie=1,ne
            call tab2io(nin,nout,nscr,a,nb,nw)
            nmu=n2h
            do imu=1,nmu
               call tab1io(nin,nout,nscr,a,nb,nw)
               do while (nb.ne.0)
                  call moreio(nin,nout,nscr,a,nb,nw)
               enddo
            enddo
         enddo
      else if (law.lt.0) then
         if (mth.ne.18) then
            write(strng,'(''illegal endf6 law for mt='',i3)') mth
            call error('file6',strng,' ')
         endif
      else if (law.ne.0.and.law.ne.3.and.law.ne.4) then
         call error('file6','illegal endf6 law.',' ')
      endif
   enddo
   return
   end subroutine file6

   subroutine file7(nin,nout,nscr,a)
   !-------------------------------------------------------------------
   ! Convert mode of File 7.
   ! Thermal neutron scattering law.
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::it2,n,i,ii,lt,lthr,nb,nw
   real(kr)::b
   integer::ia1(3)
   character(60)::strng
   real(kr),parameter::zero=0

   !--incoherent inelastic (mt=4)
   if (mth.eq.4) then
      call listio(nin,nout,nscr,a,nb,nw)
      b=a(7)
      it2=0
      ia1(1)=0
      ia1(2)=0
      ia1(3)=0
      if (iverf.ge.6) then
         it2=n2h
         if (it2.gt.0.and.it2.le.3) then
            do i=1,it2
               ii=6*i+7
               ia1(i)=nint(a(ii))
            enddo
         elseif (it2.gt.3) then
            call error('file7','bad NS in mt=4.',' ')
         endif
      endif
      do while (nb.ne.0)
         call moreio(nin,nout,nscr,a,nb,nw)
      enddo
      if (b.eq.zero) return
      call tab2io(nin,nout,nscr,a,nb,nw)
      n=n2h
      do i=1,n
         call tab1io(nin,nout,nscr,a,nb,nw)
         lt=l1h
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
         do while (lt.gt.0)
            call listio(nin,nout,nscr,a,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,nout,nscr,a,nb,nw)
            enddo
            lt=lt-1
         enddo
      enddo

      !--effective temperatures
      if (iverf.ge.6) then
         call tab1io(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
         if (it2.ne.0) then
            do i=1,it2
               if (ia1(i).eq.0) then
                  call tab1io(nin,nout,nscr,a,nb,nw)
                  do while (nb.ne.0)
                     call moreio(nin,nout,nscr,a,nb,nw)
                  enddo
               endif
            enddo
         endif
      endif

   !--elastic (mt=2)
   else if (mth.eq.2.and.iverf.ge.6) then
      lthr=l1h

      !--coherent (lthr=1)
      if (lthr.eq.1) then
         call tab1io(nin,nout,nscr,a,nb,nw)
         lt=l1h
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
         do while (lt.gt.0)
            call listio(nin,nout,nscr,a,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,nout,nscr,a,nb,nw)
            enddo
            lt=lt-1
         enddo

      !--incoherent (lthr=2)
      else if (lthr.eq.2) then
         call tab1io(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo

      !--illegal lthr for mt=2
      else
         write(strng,'(''illegal value of lthr='',i4)') lthr
         call error('file7',strng,' ')
      endif

   !--illegal mt
   else
      write(strng,'(''illegal mt='',i3)') mth
      call error('file7',strng,' ')
   endif
   return
   end subroutine file7

   subroutine file8(nin,nout,nscr,a)
   !-------------------------------------------------------------------
   ! Convert mode of File 8.
   ! Radioactive decay and fission product yield data.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::ne,ie,ns,is,lcon,no,lfs,nb,nw

   !--fission product yield data.
   if (mth.eq.454.or.mth.eq.459) then
      ne=l1h
      do ie=1,ne
         call listio(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
      enddo

   !--radioactive decay data.
   else if (mth.eq.457) then
      ns=n2h
      call listio(nin,nout,nscr,a,nb,nw)
      call listio(nin,nout,nscr,a,nb,nw)
      do while (nb.ne.0)
         call moreio(nin,nout,nscr,a,nb,nw)
      enddo
      if (ns.gt.0) then
         do is=1,ns
            call listio(nin,nout,nscr,a,nb,nw)
            lcon=l1h
            if (lcon.ne.1) then
               ne=n2h
               do ie=1,ne
                  call listio(nin,nout,nscr,a,nb,nw)
               enddo
            endif
            if (lcon.ne.0) then
               call tab1io(nin,nout,nscr,a,nb,nw)
               do while (nb.ne.0)
                  call moreio(nin,nout,nscr,a,nb,nw)
               enddo
            endif
         enddo
      endif

   !--radioactive nuclide production.
   else
      ns=n1h
      no=n2h
      if (no.ne.1) then
         do lfs=1,ns
            call listio(nin,nout,nscr,a,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,nout,nscr,a,nb,nw)
            enddo
         enddo
      else
         do lfs=1,ns
            call contio(nin,nout,nscr,a,nb,nw)
         enddo
      endif
   endif
   return
   end subroutine file8

   subroutine file9(nin,nout,nscr,a)
   !-------------------------------------------------------------------
   ! Convert mode of File 9 or File 10.
   ! Yields or cross sections for radioactive nuclide production.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::ns,lfs,nb,nw

   ns=n1h
   do lfs=1,ns
      call tab1io(nin,nout,nscr,a,nb,nw)
      do while (nb.ne.0)
         call moreio(nin,nout,nscr,a,nb,nw)
      enddo
   enddo
   return
   end subroutine file9

   subroutine file1x(nin,nout,nscr,a)
   !-------------------------------------------------------------------
   ! Convert mode of File 12 or File 13.
   ! Photon yields (MF12) or photon cross sections (MF13).
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides error
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::lo,nk,nkp,k,nb,nw,ng,ng460

   !--special path for beta delayed gammas (mt=460)
   lo=l1h
   if (mfh.eq.12.and.mth.eq.460) then
      if (lo.eq.1) then
         ng460=n1h
         do ng=1,ng460+1
            call tab1io(nin,nout,nscr,a,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,nout,nscr,a,nb,nw)
            enddo
         enddo
      else if (lo.eq.2) then
         call listio(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
      else
         call error('file12','bad LO in mt=460.',' ')
      endif

   !--normal path for other mf1x formats
   else
      if (lo.ne.2) then
         nk=n1h
         nkp=nk+1
         if (nk.eq.1) nkp=nkp-1
         do k=1,nkp
            call tab1io(nin,nout,nscr,a,nb,nw)
            do while (nb.ne.0)
              call moreio(nin,nout,nscr,a,nb,nw)
            enddo
         enddo
      else
         call listio(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
      endif
   endif
   return
   end subroutine file1x

   subroutine file14(nin,nout,nscr,a)
   !-------------------------------------------------------------------
   ! Convert mode of File 14.
   ! Photon angular distributions.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::li,ltt,nk,ni,i,ne,j,nb,nw

   li=l1h
   if (li.eq.1) return
   ltt=l2h
   nk=n1h
   ni=n2h
   do i=1,nk
      if (i.le.ni) then
         call contio(nin,nout,nscr,a,nb,nw)
      else
         call tab2io(nin,nout,nscr,a,nb,nw)
         ne=n2h
         if (ltt.ne.2) then
            do j=1,ne
               call listio(nin,nout,nscr,a,nb,nw)
               do while (nb.ne.0)
                  call moreio(nin,nout,nscr,a,nb,nw)
               enddo
            enddo
         else
            do j=1,ne
               call tab1io(nin,nout,nscr,a,nb,nw)
                  do while (nb.ne.0)
                  call moreio(nin,nout,nscr,a,nb,nw)
               enddo
            enddo
         endif
      endif
   enddo
   return
   end subroutine file14

   subroutine file15(nin,nout,nscr,a)
   !-------------------------------------------------------------------
   ! Convert mode of File 15.
   ! Photon secondary energy distributions.
   !-------------------------------------------------------------------
   use util ! provides error
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::lf,ne,i,nb,nw

   call tab1io(nin,nout,nscr,a,nb,nw)
   lf=l2h
   do while (nb.ne.0)
      call moreio(nin,nout,nscr,a,nb,nw)
   enddo

   !--tabulated distribution.
   if (lf.eq.1) then
      call tab2io(nin,nout,nscr,a,nb,nw)
      ne=n2h
      do i=1,ne
         call tab1io(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
      enddo

   !--special lanl format.
   else if (lf.eq.2) then
      call tab2io(nin,nout,nscr,a,nb,nw)
      ne=l2h
      do i=1,ne
         call listio(nin,nout,nscr,a,nb,nw)
         do while (nb.ne.0)
            call moreio(nin,nout,nscr,a,nb,nw)
         enddo
      enddo

   !--illegal lf
   else
      call error('file15','illegal lf.',' ')
   endif
   return
   end subroutine file15

   subroutine file3x(nin,nout,nscr,a)
   !-------------------------------------------------------------------
   ! Convert mode of File 31 or File 33.
   ! Cross section covariances.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::nl,l,nc,ni,n,i,nb,nw

   nl=l2h
   if (iverf.ge.5) nl=n2h
   if (nl.eq.0) return
   do l=1,nl
      call contio(nin,nout,nscr,a,nb,nw)
      nc=n1h
      ni=n2h
      if (nc.ne.0) then
         do n=1,nc
            if (iverf.ge.5) call contio(nin,nout,nscr,a,nb,nw)
            call listio(nin,nout,nscr,a,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,nout,nscr,a,nb,nw)
            enddo
         enddo
      endif
      if (ni.ne.0) then
         do i=1,ni
            call listio(nin,nout,nscr,a,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,nout,nscr,a,nb,nw)
            enddo
         enddo
      endif
   enddo
   return
   end subroutine file3x

   subroutine file32(nin,nout,nscr,a)
   !-------------------------------------------------------------------
   ! Convert mode of File 32.
   ! Resonance-parameter covariances.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides error
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::nis,i,nls,l,nb,nw,ner,j,lru,nro,ni,k,lcomp,nsrs,nlrs,nn,ndigit
   integer::isr,lrf,njsx,njs,nsg

   nis=n1h
   do i=1,nis
      call contio(nin,nout,nscr,a,nb,nw)
      ner=n1h
      do j=1,ner
         call contio(nin,nout,nscr,a,nb,nw)
         lru=l1h
         lrf=l2h
         nro=n1h
         if (nro.gt.0) then
            call contio(nin,nout,nscr,a,nb,nw)
            ni=n2h
            do k=1,ni
               call listio(nin,nout,nscr,a,nb,nw)
               do while (nb.ne.0)
                  call moreio(nin,nout,nscr,a,nb,nw)
               enddo
            enddo
         endif
         call contio(nin,nout,nscr,a,nb,nw)
         lcomp=l2h
         nls=n1h
         isr=n2h
         if (lrf.ne.7) then
            if (isr.ne.0) then
               if (isr.eq.1.and.lrf.le.2) then
                  call contio(nin,nout,nscr,a,nb,nw)
               else if (isr.eq.1.and.lrf.eq.3) then
                  call listio(nin,nout,nscr,a,nb,nw)
               else
                  call error('file32','unknown isr',' ')
               endif
            endif
            if (lru.eq.2) then
               do l=1,nls
                  call listio(nin,nout,nscr,a,nb,nw)
                  do while (nb.ne.0)
                     call moreio(nin,nout,nscr,a,nb,nw)
                  enddo
               enddo
               call listio(nin,nout,nscr,a,nb,nw)
               do while (nb.ne.0)
                  call moreio(nin,nout,nscr,a,nb,nw)
               enddo
            else if (lcomp.eq.0) then
               do l=1,nls
                  call listio(nin,nout,nscr,a,nb,nw)
                  do while (nb.ne.0)
                     call moreio(nin,nout,nscr,a,nb,nw)
                  enddo
               enddo
            else if (lcomp.eq.1) then
               call contio(nin,nout,nscr,a,nb,nw)
               nsrs=n1h
               nlrs=n2h
               do k=1,nsrs
                  call listio(nin,nout,nscr,a,nb,nw)
                  do while (nb.ne.0)
                     call moreio(nin,nout,nscr,a,nb,nw)
                  enddo
               enddo
               if (nlrs.gt.0) then
                  do k=1,nsrs
                     call listio(nin,nout,nscr,a,nb,nw)
                     do while (nb.ne.0)
                        call moreio(nin,nout,nscr,a,nb,nw)
                     enddo
                  enddo
               endif
            else if (lcomp.eq.2) then
               call listio(nin,nout,nscr,a,nb,nw)
               do while (nb.ne.0)
                  call moreio(nin,nout,nscr,a,nb,nw)
               enddo
               call contio(nin,nout,nscr,a,nb,nw)
               nn=n1h
               ndigit=l1h
               if (ndigit.lt.2 .or. ndigit.gt.6) then
                  if(ndigit.eq.0) then
                    ndigit=2
                    call mess('file32','1illegal value of ndigit',&
                              'set default, ndigit=2')
                  else
                    call error('file32','illegal value of ndigit',' ')
                  endif
               endif
               do k=1,nn
                  nw=ndigit
                  call intgio(nin,nout,nscr,a,nb,nw)
               end do
            else
               call error('file32','illegal value of lcomp',' ')
            endif
         else
            if (isr.ne.0) call error('file32',&
              'isr>0 not supported for lrf=7',' ')
            if (lcomp.eq.1) then
               call contio(nin,nout,nscr,a,nb,nw)
               nsg=n1h
               do k=1,nsg
                  call contio(nin,nout,nscr,a,nb,nw)
                  njsx=l1h
                  do l=1,njsx
                     call listio(nin,nout,nscr,a,nb,nw)
                     do while (nb.ne.0)
                        call moreio(nin,nout,nscr,a,nb,nw)
                     enddo
                  enddo
               enddo
               call listio(nin,nout,nscr,a,nb,nw)
               do while (nb.ne.0)
                  call moreio(nin,nout,nscr,a,nb,nw)
               enddo
            else if (lcomp.eq.2) then
               njs=n1h
               call listio(nin,nout,nscr,a,nb,nw)
               do while (nb.ne.0)
                  call moreio(nin,nout,nscr,a,nb,nw)
               enddo
               do k=1,njs
                  call listio(nin,nout,nscr,a,nb,nw)
                  do while (nb.ne.0)
                     call moreio(nin,nout,nscr,a,nb,nw)
                  enddo
                  call listio(nin,nout,nscr,a,nb,nw)
                  do while (nb.ne.0)
                     call moreio(nin,nout,nscr,a,nb,nw)
                  enddo
               enddo
               call contio(nin,nout,nscr,a,nb,nw)
               nn=n1h
               ndigit=l1h
               if (ndigit.lt.2.or.ndigit.gt.6) then
                  if (ndigit.eq.0) then
                    ndigit=2
                    call mess('file32','2illegal value of ndigit',&
                              'set default, ndigit=2')
                  else
                    call error('file32','illegal value of ndigit',' ')
                  endif
               endif
               do k=1,nn
                  nw=ndigit
                  call intgio(nin,nout,nscr,a,nb,nw)
               enddo
            else
               call error('file32','illegal value of lcomp',' ')
            endif
         endif
      enddo
   enddo
   return
   end subroutine file32

   subroutine file34(nin,nout,nscr,a)
   !-------------------------------------------------------------------
   ! Convert mode of File 34
   ! Covariances for angular distributions.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::nmt1,i,mt1,nl,nl1,l,l1,ni,j,nb,nw

   nmt1=n2h
   do i=1,nmt1
      call contio(nin,nout,nscr,a,nb,nw)
      mt1=l2h
      nl=n1h
      nl1=n2h
      do l=1,nl
         do l1=1,nl1
            if (mt1.ne.mth.or.l1.ge.l) then
               call contio(nin,nout,nscr,a,nb,nw)
               ni=n2h
               do j=1,ni
                  call listio(nin,nout,nscr,a,nb,nw)
                  do while (nb.ne.0)
                     call moreio(nin,nout,nscr,a,nb,nw)
                  enddo
               enddo
            endif
         enddo
      enddo
   enddo
   return
   end subroutine file34

   subroutine file35(nin,nout,nscr,a)
   !-------------------------------------------------------------------
   ! Convert mode of File 35.
   ! Covariances for secondary energy distributions.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::nk,i,nb,nw

   nk=n1h
   do i=1,nk
      call listio(nin,nout,nscr,a,nb,nw)
      do while (nb.ne.0)
         call moreio(nin,nout,nscr,a,nb,nw)
      enddo
   enddo
   return
   end subroutine file35

   subroutine file40(nin,nout,nscr,a)
   !-------------------------------------------------------------------
   ! Convert mode of File 40.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   ! external
   integer::nin,nout,nscr
   real(kr)::a(*)
   ! internals
   integer::ns,i,nl,j,nc,ni,k,l,nb,nw

   ns=n1h
   do i=1,ns
      call contio(nin,nout,nscr,a,nb,nw)
      nl=n2h
      do j=1,nl
         call contio(nin,nout,nscr,a,nb,nw)
         nc=n1h
         ni=n2h
         do k=1,nc
            call contio(nin,nout,nscr,a,nb,nw)
            call listio(nin,nout,nscr,a,nb,nw)
            do while (nb.ne.0)
               call moreio(nin,nout,nscr,a,nb,nw)
            enddo
         enddo
         if (ni.ne.0) then
            do l=1,ni
               call listio(nin,nout,nscr,a,nb,nw)
               do while (nb.ne.0)
                  call moreio(nin,nout,nscr,a,nb,nw)
               enddo
            enddo
         endif
      enddo
   enddo
   return
   end subroutine file40

   subroutine glstio(nin,nout,nscr,a,nb,nw)
   !-------------------------------------------------------------------
   ! Convert a GENDF File 1 MT451 list record to opposite format.
   ! glstio processes the entire record.
   ! Text GOUT tapes do not have titles.
   !-------------------------------------------------------------------
   use endf ! provides endf routines and variables
   ! externals
   integer::nin,nout,nscr,nb,nw,nbt,nwt
   real(kr)::a(*)
   ! internals
   integer::n1,n2,l,noff,i

   !--read in rest of header record.
   !--binary outputs can also be written.
   n1=nout
   if (n1.gt.0) n1=0
   n2=nscr
   if (n2.gt.0) n2=0
   call contio(0,n1,n2,a(1),nb,nw)
   ntw=n2h
   l=7
   call listio(nin,n1,n2,a(l),nb,nw)
   ng=l1h
   do while (nb.ne.0)
      l=l+nw
      call moreio(nin,n1,n2,a(l),nb,nw)
   enddo
   nw=l+nw-7
   if (nout.le.0.and.nscr.le.0) return

   !--remove title from binary input and write bcd outputs, if any.
   n1=nout
   if (n1.lt.0) n1=0
   n2=nscr
   if (n2.lt.0) n2=0
   a(6)=1
   call contio(0,n1,n2,a(1),nbt,nwt)
   noff=ntw-1
   nw=nw-noff
   do i=2,nw
      l=12+i
      a(l)=a(l+noff)
   enddo
   a(11)=a(11)-noff
   a(13)=0
   l=7
   call listio(0,n1,n2,a(l),nb,nw)
   do while (nb.ne.0)
      l=l+nw
      call moreio(0,n1,n2,a(l),nb,nw)
   enddo
   return
   end subroutine glstio

end module modem

